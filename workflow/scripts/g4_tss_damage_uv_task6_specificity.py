#!/usr/bin/env python3
"""Task 6: Specificity and sensitivity checks.

Tests whether the lesion-burden effect is promoter-G4-specific and robust to
alternate burden definitions.

Steps:
1. Build promoter lesion burden for No_overlap genes (no G4 filter, raw TSS windows)
2. Match G4_TSS and No_overlap by baseline-expression decile
3. Fit specificity interaction model: lfc_t ~ damage_ds0_z * promoter_group + log2_tpm_t0
4. Alternate burden metrics (damage_mean, damage_sum, real_rpkm_max)
5. GC-rich background control
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf


TIMEPOINTS = ["12", "30", "60"]
PSEUDO = 1e-3
DS0_TIME = 0
MIN_GENES_GC = 10


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cpd-table", required=True)
    p.add_argument("--pp64-table", required=True)
    p.add_argument("--gene-table", required=True,
                   help="promoter_g4_damage_gene_table.tsv")
    p.add_argument("--annotation", required=True)
    p.add_argument("--tss-windows", required=True)
    p.add_argument("--gc-rich-bg", required=True)
    p.add_argument("--uv-master", required=True)
    p.add_argument("--baseline-tpm", required=True)
    p.add_argument("--sample-manifest", required=True)
    p.add_argument("--out-interaction", required=True)
    p.add_argument("--out-matched-summary", required=True)
    p.add_argument("--out-cpd-sensitivity", required=True)
    p.add_argument("--out-64pp-sensitivity", required=True)
    p.add_argument("--out-interaction-fig", required=True)
    p.add_argument("--out-gc-bg-summary", required=True)
    p.add_argument("--out-gc-bg-lfc", required=True)
    p.add_argument("--out-gc-bg-vs-g4", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def strip_version(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace(r"\.\d+$", "", regex=True)


def count_reads_in_bed(path: str) -> int:
    n = 0
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#") and line.strip():
                n += 1
    return n


def bedtools_count_per_locus(loci_bed: str, reads_bed: str) -> pd.Series:
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
    tf = tempfile.NamedTemporaryFile(suffix=".bed", mode="w", delete=False)
    for p in paths:
        if Path(p).exists():
            with open(p) as fh:
                for line in fh:
                    if line.strip():
                        tf.write(line)
    tf.close()
    return tf.name


def compute_burden_for_windows(
    windows_df: pd.DataFrame,
    manifest: list[dict],
    product_filter: str,
    time_filter: int,
    log_fn,
) -> pd.DataFrame:
    """Compute log2(real/sim) burden in arbitrary genomic windows.

    windows_df must have columns: gene_id, chrom, start, end, locus_id.
    Returns gene-level damage_max for DS0.
    """
    # Write loci BED
    with tempfile.NamedTemporaryFile(suffix=".bed", mode="w", delete=False) as tf:
        loci_bed_path = tf.name
        for _, row in windows_df.drop_duplicates("locus_id").iterrows():
            tf.write(f"{row.chrom}\t{row.start}\t{row.end}\t{row.locus_id}\t.\t.\n")

    samples = [s for s in manifest
               if s["product"] == product_filter and s["time_after_exposure"] == time_filter]
    log_fn(f"  Computing burden for {len(samples)} {product_filter} t={time_filter} samples")

    locus_records = []
    for sample in samples:
        real_beds = [b for b in [sample["real_plus_bed"], sample["real_minus_bed"]]
                     if Path(b).exists()]
        sim_beds = [b for b in [sample["sim_plus_bed"], sample["sim_minus_bed"]]
                    if Path(b).exists()]
        if not real_beds:
            continue
        total_real = sum(count_reads_in_bed(b) for b in real_beds)
        total_sim = sum(count_reads_in_bed(b) for b in sim_beds) if sim_beds else 0

        real_merged = merge_beds_to_temp(real_beds)
        sim_merged = merge_beds_to_temp(sim_beds) if sim_beds else None
        real_counts = bedtools_count_per_locus(loci_bed_path, real_merged)
        sim_counts = bedtools_count_per_locus(loci_bed_path, sim_merged) if sim_merged else pd.Series(dtype=int)
        Path(real_merged).unlink(missing_ok=True)
        if sim_merged:
            Path(sim_merged).unlink(missing_ok=True)

        for _, row in windows_df.drop_duplicates("locus_id").iterrows():
            lid = row["locus_id"]
            llen = row["end"] - row["start"]
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
                "gene_id": row["gene_id"],
                "locus_id": lid,
                "replicate": sample["name"],
                "log2_real_over_sim": log2_ratio,
            })

    Path(loci_bed_path).unlink(missing_ok=True)

    if not locus_records:
        return pd.DataFrame(columns=["gene_id", "damage_max"])

    locus_df = pd.DataFrame(locus_records)
    locus_agg = locus_df.groupby(["gene_id", "locus_id"]).agg(
        damage_max=("log2_real_over_sim", "max"),
        damage_mean=("log2_real_over_sim", "mean"),
    ).reset_index()
    gene_agg = locus_agg.groupby("gene_id").agg(
        damage_max=("damage_max", "max"),
        damage_mean=("damage_mean", "mean"),
    ).reset_index()
    return gene_agg


def standardize(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    df = df.copy()
    for col in cols:
        if col in df.columns:
            mu, sd = df[col].mean(), df[col].std()
            df[col] = (df[col] - mu) / sd if sd > 0 else 0.0
    return df


def sensitivity_analysis(df: pd.DataFrame, product: str, log_fn) -> pd.DataFrame:
    """Test alternate burden metrics via Spearman correlation with lfc_60."""
    if "lfc_60" not in df.columns:
        return pd.DataFrame()
    g4 = df[df["group"] == "G4_TSS"].copy() if "group" in df.columns else df.copy()
    rows = []
    alt_metrics = ["damage_ds0", "damage_mean_ds0", "damage_sum_ds0", "real_rpkm_max_ds0"]
    for metric in alt_metrics:
        if metric not in g4.columns:
            continue
        sub = g4[[metric, "lfc_60"]].dropna()
        if len(sub) < 5:
            continue
        r, p = stats.spearmanr(sub[metric], sub["lfc_60"])
        log_fn(f"  {product} {metric} vs lfc_60: r={r:.3f}, p={p:.3e}, n={len(sub)}")
        rows.append({
            "product": product, "metric": metric,
            "n": len(sub), "spearman_r": r, "pvalue": p,
        })
    return pd.DataFrame(rows)


def plot_interaction(results: pd.DataFrame, out_path: str) -> None:
    if results.empty:
        plt.figure()
        plt.text(0.5, 0.5, "No interaction results", ha="center", va="center")
        plt.savefig(out_path)
        plt.close()
        return

    fig, ax = plt.subplots(figsize=(7, 5))
    colors = {"G4_TSS": "#1f77b4", "No_overlap": "#7f7f7f"}
    tps = results["timepoint"].unique() if "timepoint" in results.columns else []
    for group, color in colors.items():
        g = results[results.get("group", pd.Series()) == group] if "group" in results.columns else pd.DataFrame()
        if g.empty:
            continue
        ax.plot(g["timepoint"], g.get("coef_damage", pd.Series()), marker="o",
                label=group, color=color, linewidth=1.5)
    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
    ax.set_xlabel("Timepoint (min)")
    ax.set_ylabel("Damage coefficient (lfc ~ damage)")
    ax.set_title("Promoter-group specificity: damage effect by group")
    ax.legend()
    plt.tight_layout()
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else sys.stdout

    def log(msg: str) -> None:
        print(msg, file=log_fh, flush=True)

    for path in [args.out_interaction]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    log("=== Task 6: Specificity and sensitivity checks ===")

    manifest = json.loads(Path(args.sample_manifest).read_text())
    anno = pd.read_csv(args.annotation, sep="\t")
    anno["gene_id_nv"] = strip_version(anno["gene_id"])

    uv = pd.read_csv(args.uv_master, sep="\t")
    uv["gene_id_nv"] = strip_version(uv["gene_id"])

    tpm = pd.read_csv(args.baseline_tpm, sep="\t")
    tpm_id_col = "gene_id" if "gene_id" in tpm.columns else tpm.columns[0]
    tpm_val_col = "mean_tpm" if "mean_tpm" in tpm.columns else tpm.columns[1]
    tpm = tpm.rename(columns={tpm_id_col: "gene_id_nv", tpm_val_col: "mean_tpm_t0"})
    tpm["gene_id_nv"] = strip_version(tpm["gene_id_nv"])
    tpm["log2_tpm_t0"] = np.log2(tpm["mean_tpm_t0"].fillna(0) + 1)

    gene_table = pd.read_csv(args.gene_table, sep="\t")

    # ── Step 1: Build No_overlap promoter lesion table ────────────────────────
    log("\n--- Step 1: No_overlap promoter lesion burden ---")
    no_overlap_genes = anno[anno["group"] == "No_overlap"].copy()
    log(f"No_overlap genes: {len(no_overlap_genes):,}")

    tss = pd.read_csv(args.tss_windows, sep="\t", header=None, usecols=[0, 1, 2, 3],
                      names=["chrom", "start", "end", "gene_id"])
    tss["gene_id_nv"] = strip_version(tss["gene_id"])

    no_overlap_tss = tss[tss["gene_id_nv"].isin(no_overlap_genes["gene_id_nv"])].copy()
    no_overlap_tss["locus_id"] = (
        no_overlap_tss["chrom"] + ":" +
        no_overlap_tss["start"].astype(str) + "-" +
        no_overlap_tss["end"].astype(str) + ":" +
        no_overlap_tss["gene_id"]
    )
    log(f"No_overlap TSS windows: {len(no_overlap_tss):,}")

    interaction_rows = []
    matched_rows = []
    cpd = pd.read_csv(args.cpd_table, sep="\t")
    pp64 = pd.read_csv(args.pp64_table, sep="\t")

    for product, master_df in [("CPD", cpd), ("64-PP", pp64)]:
        log(f"\n  Computing No_overlap burden for {product}...")
        no_burden = compute_burden_for_windows(
            no_overlap_tss, manifest, product, DS0_TIME, log
        )
        no_burden["gene_id_nv"] = strip_version(no_burden["gene_id"])
        no_burden = no_burden.merge(tpm[["gene_id_nv", "log2_tpm_t0"]], on="gene_id_nv", how="left")
        no_burden = no_burden.merge(uv[["gene_id_nv", "lfc_60", "uv_trajectory"]], on="gene_id_nv", how="left")
        no_burden["promoter_group"] = "No_overlap"
        no_burden = no_burden.rename(columns={"damage_max": "damage_ds0"})
        log(f"  No_overlap genes with burden: {no_burden['damage_ds0'].notna().sum():,}")

        # G4_TSS from master table
        g4 = master_df[master_df.get("group", pd.Series("G4_TSS", index=master_df.index)) == "G4_TSS"].copy() \
            if "group" in master_df.columns else master_df.copy()
        g4 = g4.dropna(subset=["damage_ds0"])
        g4["promoter_group"] = "G4_TSS"

        # ── Step 2: Match by expression decile ───────────────────────────────
        log(f"\n  Matching by expression decile ({product})...")
        g4["expr_decile"] = pd.qcut(g4["log2_tpm_t0"].fillna(0), q=10, labels=False, duplicates="drop")
        no_burden["expr_decile"] = pd.qcut(no_burden["log2_tpm_t0"].fillna(0), q=10, labels=False, duplicates="drop")

        matched_rows_prod = []
        for decile in range(10):
            g4_d = g4[g4["expr_decile"] == decile]["damage_ds0"].dropna()
            no_d = no_burden[no_burden["expr_decile"] == decile]["damage_ds0"].dropna()
            if len(g4_d) >= 3 and len(no_d) >= 3:
                matched_rows_prod.append({
                    "product": product, "expr_decile": decile,
                    "n_g4_tss": len(g4_d), "n_no_overlap": len(no_d),
                    "median_burden_g4": g4_d.median(),
                    "median_burden_no": no_d.median(),
                })
        matched_rows.extend(matched_rows_prod)

        # ── Step 3: Interaction model ─────────────────────────────────────────
        log(f"\n  Interaction model ({product})...")
        combined = pd.concat([
            g4[["damage_ds0", "lfc_60", "log2_tpm_t0", "promoter_group"]].dropna(),
            no_burden[["damage_ds0", "lfc_60", "log2_tpm_t0", "promoter_group"]].dropna(),
        ], ignore_index=True)
        combined["damage_ds0_z"] = (combined["damage_ds0"] - combined["damage_ds0"].mean()) / \
                                    combined["damage_ds0"].std()
        combined["promoter_group"] = combined["promoter_group"].astype("category")

        if len(combined) >= 30 and "lfc_60" in combined.columns:
            for tp in TIMEPOINTS:
                lfc_col = f"lfc_{tp}"
                # For this analysis, lfc_60 is the available column from merged data
                if lfc_col not in combined.columns:
                    # We only have lfc_60 in the combined frame
                    if tp != "60" or "lfc_60" not in combined.columns:
                        continue
                    lfc_col = "lfc_60"

                sub = combined[["damage_ds0_z", lfc_col, "log2_tpm_t0", "promoter_group"]].dropna()
                if len(sub) < 30:
                    continue
                formula = f"{lfc_col} ~ damage_ds0_z * promoter_group + log2_tpm_t0"
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        model = smf.ols(formula, data=sub).fit()
                    # Find interaction term
                    inter_terms = [p for p in model.params.index
                                   if "damage_ds0_z" in p and "promoter_group" in p]
                    for term in inter_terms:
                        coef = model.params[term]
                        pval = model.pvalues[term]
                        log(f"    {product} t={tp}: interaction {term}: coef={coef:.3f}, p={pval:.3e}")
                        interaction_rows.append({
                            "product": product, "timepoint": tp,
                            "interaction_term": term, "coef": coef,
                            "se": model.bse.get(term, np.nan),
                            "pvalue": pval, "n": len(sub),
                        })
                except Exception as exc:
                    log(f"    Interaction model failed: {exc}")

    pd.DataFrame(interaction_rows).to_csv(args.out_interaction, sep="\t", index=False)
    pd.DataFrame(matched_rows).to_csv(args.out_matched_summary, sep="\t", index=False)
    log(f"\nWritten: {args.out_interaction}")
    log(f"Written: {args.out_matched_summary}")

    # Interaction plot
    plot_interaction(pd.DataFrame(interaction_rows), args.out_interaction_fig)

    # ── Steps 4-5: Sensitivity analyses ──────────────────────────────────────
    log("\n--- Step 4: Alternate burden metrics ---")
    cpd_sens = sensitivity_analysis(cpd, "CPD", log)
    pp64_sens = sensitivity_analysis(pp64, "64-PP", log)
    cpd_sens.to_csv(args.out_cpd_sensitivity, sep="\t", index=False)
    pp64_sens.to_csv(args.out_64pp_sensitivity, sep="\t", index=False)

    # ── Step 6: GC-rich background control ───────────────────────────────────
    log("\n--- Step 6: GC-rich background control ---")
    gc_bed = pd.read_csv(args.gc_rich_bg, sep="\t", header=None,
                         names=["chrom", "start", "end"] + [f"c{i}" for i in range(10)])
    log(f"GC-rich background loci: {len(gc_bed):,}")

    # Use locus_id as position identifier
    gc_bed["locus_id"] = (
        gc_bed["chrom"] + ":" + gc_bed["start"].astype(str) + "-" + gc_bed["end"].astype(str)
    )
    gc_bed["gene_id"] = gc_bed["locus_id"]  # treat each locus as its own "gene"

    gc_summary_rows = []
    gc_lfc_rows = []
    gc_vs_g4_rows = []

    for product, master_df in [("CPD", cpd), ("64-PP", pp64)]:
        log(f"\n  GC-rich control: {product}")
        gc_burden = compute_burden_for_windows(gc_bed, manifest, product, DS0_TIME, log)
        n_gc = gc_burden["damage_max"].notna().sum()
        log(f"  GC loci with quantified burden: {n_gc}")

        gc_summary_rows.append({
            "product": product,
            "n_loci_total": len(gc_bed),
            "n_loci_quantified": n_gc,
            "median_damage_max": gc_burden["damage_max"].median() if n_gc > 0 else np.nan,
            "mean_damage_max": gc_burden["damage_max"].mean() if n_gc > 0 else np.nan,
        })

        if n_gc < MIN_GENES_GC:
            log(f"  Skipping correlation: too few GC loci (n={n_gc}, min={MIN_GENES_GC})")
            gc_lfc_rows.append({"product": product, "note": f"underpowered (n={n_gc})"})
        else:
            gc_lfc_rows.append({
                "product": product, "n": n_gc,
                "note": "GC-rich loci have no RNA data; correlation not applicable",
            })

        # Compare GC-rich burden vs promoter-G4 burden
        g4_burden = master_df[master_df.get("group", pd.Series()) == "G4_TSS"][
            "damage_ds0"
        ].dropna() if "group" in master_df.columns else master_df["damage_ds0"].dropna()

        if n_gc >= MIN_GENES_GC and len(g4_burden) >= 5:
            u_stat, p_val = stats.mannwhitneyu(g4_burden, gc_burden["damage_max"].dropna())
            gc_vs_g4_rows.append({
                "product": product,
                "n_g4_tss": len(g4_burden),
                "n_gc_rich": n_gc,
                "median_burden_g4": g4_burden.median(),
                "median_burden_gc": gc_burden["damage_max"].median(),
                "mann_whitney_p": p_val,
            })

    pd.DataFrame(gc_summary_rows).to_csv(args.out_gc_bg_summary, sep="\t", index=False)
    pd.DataFrame(gc_lfc_rows).to_csv(args.out_gc_bg_lfc, sep="\t", index=False)
    pd.DataFrame(gc_vs_g4_rows).to_csv(args.out_gc_bg_vs_g4, sep="\t", index=False)
    log(f"\nWritten GC-rich background summaries.")

    log("\n=== Task 6 complete ===")
    if args.log:
        log_fh.close()


if __name__ == "__main__":
    main()
