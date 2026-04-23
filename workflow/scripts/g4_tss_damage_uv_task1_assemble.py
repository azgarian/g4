#!/usr/bin/env python3
"""Task 1: Build promoter-G4 locus and gene lesion tables.

Intersects canonical TSS windows with the merged G4 peak set to define
promoter-G4 loci, then quantifies UV lesion burden (real and simulated DS
reads) at each locus. Collapses to gene-level burden metrics.

DS bed files (plus and minus strand) are strand-collapsed by concatenation
before intersection. Simulated reads represent the random-placement null
for the same library; log2(real_rpkm / sim_rpkm) is the primary burden metric.

Primary predictor: damage_ds0 = damage_max at DS time_after_exposure == 0.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


PSEUDO = 1e-3
DS0_TIME = 0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tss-windows", required=True,
                   help="canonical_tss_windows_1kb.bed (col4 = gene_id)")
    p.add_argument("--g4-merged", required=True,
                   help="g4_hela_chip_cuttag_merged.bed")
    p.add_argument("--annotation", required=True,
                   help="tss_group_annotation.tsv")
    p.add_argument("--chip-peaks", required=True,
                   help="g4_hela_peaks_prepared.tsv (ChIP)")
    p.add_argument("--cuttag-peaks", required=True,
                   help="g4_hela_peaks_prepared.tsv (CUT&Tag)")
    p.add_argument("--gc-rich-bg", required=True,
                   help="gc_rich_bg_prepared.bed")
    p.add_argument("--sample-manifest", required=True,
                   help="JSON file with DS sample metadata")
    p.add_argument("--out-locus-map", required=True)
    p.add_argument("--out-locus-table", required=True)
    p.add_argument("--out-gene-table", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def bedtools_intersect(a: str, b: str, flags: list[str] | None = None) -> str:
    """Run bedtools intersect and return stdout as string."""
    cmd = ["bedtools", "intersect", "-a", a, "-b", b] + (flags or [])
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return result.stdout


def count_reads_in_bed(path: str) -> int:
    """Count non-header lines in a BED file."""
    n = 0
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#") and line.strip():
                n += 1
    return n


def bedtools_count_per_locus(loci_bed: str, reads_bed: str) -> pd.Series:
    """
    Use bedtools intersect -c to count reads per locus.
    Returns Series indexed by locus_id (col 4 of loci_bed) with count values.
    loci_bed col 4 must be the locus_id.
    """
    cmd = ["bedtools", "intersect", "-a", loci_bed, "-b", reads_bed, "-c"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or not result.stdout.strip():
        return pd.Series(dtype=int)
    rows = []
    for line in result.stdout.strip().split("\n"):
        parts = line.split("\t")
        locus_id = parts[3]
        count = int(parts[-1])
        rows.append((locus_id, count))
    df = pd.DataFrame(rows, columns=["locus_id", "count"])
    return df.groupby("locus_id")["count"].sum()


def merge_beds_to_temp(paths: list[str]) -> str:
    """Concatenate multiple BED files into a single temp file, return path."""
    tf = tempfile.NamedTemporaryFile(suffix=".bed", mode="w", delete=False)
    for p in paths:
        if Path(p).exists():
            with open(p) as fh:
                for line in fh:
                    if line.strip():
                        tf.write(line)
    tf.close()
    return tf.name


def rank_normalize(series: pd.Series) -> pd.Series:
    ranked = series.rank(method="average", na_option="keep")
    n = ranked.notna().sum()
    if n <= 1:
        return ranked
    return (ranked - 1) / (n - 1)


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else sys.stdout

    def log(msg: str) -> None:
        print(msg, file=log_fh, flush=True)

    for path in [args.out_locus_map, args.out_locus_table, args.out_gene_table]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    log("=== Task 1: Build promoter-G4 locus and gene lesion tables ===")
    log(f"Primary predictor: DS time_after_exposure == {DS0_TIME} (DS0)")

    # ── Load inputs ───────────────────────────────────────────────────────────
    manifest = json.loads(Path(args.sample_manifest).read_text())
    log(f"Sample manifest: {len(manifest)} DS samples")

    anno = pd.read_csv(args.annotation, sep="\t")
    anno["gene_id_nv"] = anno["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    log(f"Annotation: {len(anno):,} genes; groups: {anno['group'].value_counts().to_dict()}")

    tss_cols = ["chrom", "start", "end", "gene_id"]
    tss = pd.read_csv(args.tss_windows, sep="\t", header=None, usecols=[0, 1, 2, 3],
                      names=tss_cols)
    log(f"TSS windows: {len(tss):,}")

    # ── Step 1: Intersect TSS windows with merged G4 peaks ───────────────────
    log("\n--- Step 1: Define promoter-G4 loci ---")
    with tempfile.NamedTemporaryFile(suffix=".bed", mode="w", delete=False) as tf:
        loci_raw_path = tf.name
        # Write TSS windows with gene_id in col 4 (required by bedtools)
        for _, row in tss.iterrows():
            tf.write(f"{row.chrom}\t{row.start}\t{row.end}\t{row.gene_id}\n")

    intersect_out = bedtools_intersect(
        loci_raw_path, args.g4_merged, ["-wa", "-wb"]
    )
    Path(loci_raw_path).unlink(missing_ok=True)

    if not intersect_out.strip():
        log("ERROR: no overlap between TSS windows and G4 peaks")
        sys.exit(1)

    loci_rows = []
    seen_loci = set()
    for line in intersect_out.strip().split("\n"):
        parts = line.split("\t")
        tss_chrom, tss_start, tss_end, gene_id = parts[0], int(parts[1]), int(parts[2]), parts[3]
        g4_chrom, g4_start, g4_end = parts[4], int(parts[5]), int(parts[6])
        # Locus = intersection of TSS window and G4 peak
        loc_start = max(tss_start, g4_start)
        loc_end = min(tss_end, g4_end)
        if loc_end <= loc_start:
            continue
        tss_window_id = f"{tss_chrom}:{tss_start}-{tss_end}:{gene_id}"
        promoter_g4_locus_id = f"{g4_chrom}:{g4_start}-{g4_end}:{gene_id}"
        key = (gene_id, promoter_g4_locus_id)
        if key in seen_loci:
            continue
        seen_loci.add(key)
        gene_id_nv = str(gene_id).split(".")[0]
        name_row = anno[anno["gene_id_nv"] == gene_id_nv]
        gene_name = name_row["gene_name"].iloc[0] if len(name_row) and "gene_name" in anno.columns else ""
        loci_rows.append({
            "gene_id": gene_id,
            "gene_id_nv": gene_id_nv,
            "gene_name": gene_name,
            "chrom": g4_chrom,
            "start": g4_start,
            "end": g4_end,
            "locus_start": loc_start,
            "locus_end": loc_end,
            "tss_window_id": tss_window_id,
            "promoter_g4_locus_id": promoter_g4_locus_id,
            "locus_length": loc_end - loc_start,
        })

    locus_map = pd.DataFrame(loci_rows)
    log(f"Promoter-G4 loci: {len(locus_map):,} (gene×locus pairs)")
    log(f"Unique genes with promoter-G4: {locus_map['gene_id'].nunique():,}")

    loci_per_gene = locus_map.groupby("gene_id").size()
    for n in [1, 2, 3]:
        label = f"{n}+" if n == 3 else str(n)
        count = (loci_per_gene >= n).sum() if n == 3 else (loci_per_gene == n).sum()
        log(f"  Genes with {label} promoter-G4 loci: {count:,}")

    # ── Write locus map ───────────────────────────────────────────────────────
    locus_map.drop(columns=["gene_id_nv"]).to_csv(args.out_locus_map, sep="\t", index=False)
    log(f"\nWritten locus map: {args.out_locus_map}")

    # ── Step 2: Build locus BED for intersection ──────────────────────────────
    with tempfile.NamedTemporaryFile(suffix=".bed", mode="w", delete=False) as tf:
        loci_bed_path = tf.name
        for _, row in locus_map.drop_duplicates("promoter_g4_locus_id").iterrows():
            tf.write(
                f"{row.chrom}\t{row.locus_start}\t{row.locus_end}"
                f"\t{row.promoter_g4_locus_id}\t.\t.\n"
            )
    log(f"Loci BED written: {locus_map['promoter_g4_locus_id'].nunique():,} unique loci")

    # ── Step 3: Quantify lesion burden per sample ─────────────────────────────
    log("\n--- Step 3: Quantify lesion burden per sample ---")
    locus_records = []

    for sample in manifest:
        sid = sample["sample_id"]
        product = sample["product"]
        time_min = sample["time_after_exposure"]
        name = sample["name"]

        real_beds = [sample["real_plus_bed"], sample["real_minus_bed"]]
        sim_beds = [sample["sim_plus_bed"], sample["sim_minus_bed"]]

        existing_real = [b for b in real_beds if Path(b).exists()]
        existing_sim = [b for b in sim_beds if Path(b).exists()]

        if not existing_real:
            log(f"  SKIP {name}: no real bed files found")
            continue

        # Total read counts (for RPKM normalization)
        total_real = sum(count_reads_in_bed(b) for b in existing_real)
        total_sim = sum(count_reads_in_bed(b) for b in existing_sim) if existing_sim else 0

        # Merge strands into temp files
        real_merged = merge_beds_to_temp(existing_real)
        sim_merged = merge_beds_to_temp(existing_sim) if existing_sim else None

        # Count real reads per locus
        real_counts = bedtools_count_per_locus(loci_bed_path, real_merged)
        sim_counts = bedtools_count_per_locus(loci_bed_path, sim_merged) if sim_merged else pd.Series(dtype=int)

        Path(real_merged).unlink(missing_ok=True)
        if sim_merged:
            Path(sim_merged).unlink(missing_ok=True)

        log(f"  {name} ({product} t={time_min}): total_real={total_real:,}, "
            f"loci_with_real_reads={real_counts.astype(bool).sum()}")

        # Compute RPKM per locus
        unique_loci = locus_map.drop_duplicates("promoter_g4_locus_id")[
            ["promoter_g4_locus_id", "locus_length"]
        ]
        for _, row in unique_loci.iterrows():
            lid = row["promoter_g4_locus_id"]
            llen = row["locus_length"]
            real_c = real_counts.get(lid, 0)
            sim_c = sim_counts.get(lid, 0)

            real_rpkm = (real_c * 1e3 * 1e6) / (llen * total_real) if total_real > 0 else np.nan
            sim_rpkm = (sim_c * 1e3 * 1e6) / (llen * total_sim) if total_sim > 0 else np.nan

            if not np.isnan(real_rpkm) and not np.isnan(sim_rpkm):
                log2_ratio = np.log2((real_rpkm + PSEUDO) / (sim_rpkm + PSEUDO))
            elif not np.isnan(real_rpkm):
                log2_ratio = np.log2(real_rpkm + PSEUDO)
            else:
                log2_ratio = np.nan

            locus_records.append({
                "sample_id": sid,
                "product": product,
                "lesion_time_min": time_min,
                "replicate": name,
                "promoter_g4_locus_id": lid,
                "real_count": real_c,
                "sim_count": sim_c,
                "real_rpkm": real_rpkm,
                "sim_rpkm": sim_rpkm,
                "log2_real_over_sim": log2_ratio,
            })

    Path(loci_bed_path).unlink(missing_ok=True)

    locus_df = pd.DataFrame(locus_records)
    log(f"\nLocus-level records: {len(locus_df):,}")

    # ── Step 4: Aggregate replicates to locus-level summaries ────────────────
    log("\n--- Step 4: Aggregate replicates ---")
    grp_cols = ["product", "lesion_time_min", "promoter_g4_locus_id"]
    locus_agg = locus_df.groupby(grp_cols).agg(
        n_replicates=("replicate", "count"),
        real_rpkm_mean=("real_rpkm", "mean"),
        real_rpkm_max=("real_rpkm", "max"),
        sim_rpkm_mean=("sim_rpkm", "mean"),
        damage_mean=("log2_real_over_sim", "mean"),
        damage_max=("log2_real_over_sim", "max"),
    ).reset_index()
    log(f"Locus-level summary rows: {len(locus_agg):,}")

    # Join gene_id back
    locus_agg = locus_agg.merge(
        locus_map[["promoter_g4_locus_id", "gene_id", "gene_name"]].drop_duplicates(),
        on="promoter_g4_locus_id",
        how="left",
    )

    # Write locus table
    locus_agg.to_csv(args.out_locus_table, sep="\t", index=False)
    log(f"Written locus table: {args.out_locus_table}")

    # ── Step 5: Collapse locus-level to gene-level ────────────────────────────
    log("\n--- Step 5: Collapse to gene level ---")
    gene_grp = ["gene_id", "product", "lesion_time_min"]
    gene_agg = locus_agg.groupby(gene_grp).agg(
        n_promoter_g4_loci=("promoter_g4_locus_id", "nunique"),
        damage_max=("damage_max", "max"),
        damage_mean=("damage_mean", "mean"),
        damage_sum=("damage_mean", "sum"),
        real_rpkm_max=("real_rpkm_max", "max"),
        real_rpkm_mean=("real_rpkm_mean", "mean"),
    ).reset_index()

    # Add gene_name
    gene_agg = gene_agg.merge(
        locus_map[["gene_id", "gene_name"]].drop_duplicates("gene_id"),
        on="gene_id",
        how="left",
    )

    log(f"Gene-level rows: {len(gene_agg):,}")

    # ── Step 6 & 7: Primary burden predictor (DS0) and z-score ───────────────
    log("\n--- Steps 6-7: DS0 primary predictor and z-score ---")
    ds0 = gene_agg[gene_agg["lesion_time_min"] == DS0_TIME].copy()
    ds0 = ds0.rename(columns={"damage_max": "damage_ds0"})[
        ["gene_id", "product", "damage_ds0"]
    ]

    for prod in ds0["product"].unique():
        mask = ds0["product"] == prod
        vals = ds0.loc[mask, "damage_ds0"]
        mu, sd = vals.mean(), vals.std()
        ds0.loc[mask, "damage_ds0_z"] = (vals - mu) / sd if sd > 0 else 0.0
        log(f"  {prod} DS0 damage_ds0: mean={mu:.3f}, sd={sd:.3f}, n={mask.sum():,}")

    gene_agg = gene_agg.merge(ds0, on=["gene_id", "product"], how="left")

    # ── Step 8: Join annotation and filter to G4_TSS ─────────────────────────
    log("\n--- Step 8: Join annotation ---")
    gene_id_nv_map = anno[["gene_id", "gene_id_nv", "group"]].copy()
    gene_agg["gene_id_nv"] = gene_agg["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    gene_agg = gene_agg.merge(gene_id_nv_map, on="gene_id_nv", how="left", suffixes=("", "_anno"))
    gene_agg["group"] = gene_agg["group"].fillna("unknown")

    g4_tss_mask = gene_agg["group"] == "G4_TSS"
    log(f"G4_TSS genes with DS0 burden (any product): "
        f"{gene_agg[g4_tss_mask & (gene_agg['lesion_time_min'] == DS0_TIME)]['gene_id'].nunique():,}")

    gene_agg = gene_agg.drop(columns=["gene_id_nv"], errors="ignore")
    gene_agg.to_csv(args.out_gene_table, sep="\t", index=False)
    log(f"Written gene table: {args.out_gene_table}")

    # ── Step 9: Burden distribution and mapping checks ────────────────────────
    log("\n--- Step 9: Burden distribution and mapping checks ---")
    log(f"Total unique promoter-G4 loci: {locus_map['promoter_g4_locus_id'].nunique():,}")
    g4_tss_genes_with_burden = gene_agg[
        (gene_agg["group"] == "G4_TSS") & (gene_agg["lesion_time_min"] == DS0_TIME)
    ]["gene_id"].nunique()
    log(f"G4_TSS genes with >=1 quantified locus: {g4_tss_genes_with_burden:,}")

    loci_per_gene = locus_map.groupby("gene_id").size()
    for n in [1, 2]:
        log(f"  Genes with exactly {n} promoter-G4 loci: {(loci_per_gene == n).sum():,}")
    log(f"  Genes with 3+ promoter-G4 loci: {(loci_per_gene >= 3).sum():,}")

    # Correlation with G4 signal
    try:
        from g4_tss_damage_uv_task2_master_join import bedtools_max_g4_signal
    except ImportError:
        bedtools_max_g4_signal = None

    # G4 signal correlation (re-implement here for Task 1 report)
    try:
        chip = pd.read_csv(args.chip_peaks, sep="\t")
        cuttag = pd.read_csv(args.cuttag_peaks, sep="\t")

        def max_signal_per_gene(peaks_df: pd.DataFrame, loci: pd.DataFrame) -> pd.Series:
            """Max peak signal for peaks overlapping any promoter-G4 locus of each gene."""
            result = {}
            for _, locus in loci.iterrows():
                ovlp = peaks_df[
                    (peaks_df["chrom"] == locus["chrom"]) &
                    (peaks_df["start"] < locus["locus_end"]) &
                    (peaks_df["end"] > locus["locus_start"])
                ]
                if not ovlp.empty and "source_signal" in ovlp.columns:
                    gid = locus["gene_id"]
                    result[gid] = max(result.get(gid, 0), ovlp["source_signal"].max())
            return pd.Series(result)

        chip_sig = max_signal_per_gene(chip, locus_map)
        cuttag_sig = max_signal_per_gene(cuttag, locus_map)
        g4_sig = pd.DataFrame({"chip": chip_sig, "cuttag": cuttag_sig}).fillna(0)
        g4_sig["norm"] = (
            rank_normalize(g4_sig["chip"]) + rank_normalize(g4_sig["cuttag"])
        ) / 2

        ds0_gene = gene_agg[
            (gene_agg["lesion_time_min"] == DS0_TIME)
        ].drop_duplicates("gene_id").set_index("gene_id")["damage_ds0"]

        for prod in gene_agg["product"].unique():
            ds0_prod = gene_agg[
                (gene_agg["lesion_time_min"] == DS0_TIME) & (gene_agg["product"] == prod)
            ].set_index("gene_id")["damage_ds0"]
            common = ds0_prod.index.intersection(g4_sig.index)
            if len(common) < 5:
                continue
            pear = stats.pearsonr(ds0_prod[common], g4_sig.loc[common, "norm"])
            spear = stats.spearmanr(ds0_prod[common], g4_sig.loc[common, "norm"])
            log(f"\nCorrelation damage_ds0 ({prod}) vs G4 signal (n={len(common)}):")
            log(f"  Pearson  r={pear[0]:.3f}  p={pear[1]:.3e}")
            log(f"  Spearman r={spear[0]:.3f}  p={spear[1]:.3e}")
    except Exception as exc:
        log(f"  WARNING: G4 signal correlation skipped: {exc}")

    log("\n=== Task 1 complete ===")
    if args.log:
        log_fh.close()


if __name__ == "__main__":
    main()
