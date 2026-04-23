#!/usr/bin/env python3
"""Task 2: Assemble the lesion × RNA master table.

Joins promoter-G4 lesion burden (DS0) to the existing UV-response master table,
baseline expression, and G4 strength covariates. Produces separate tables for
CPD and 64-PP.

Sign convention: positive lfc = higher at t=0 than post-UV = repression.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


LFC_THRESH = 0.5
PADJ_THRESH = 0.05
DS0_TIME = 0
TIMEPOINTS = ["12", "30", "60"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gene-table", required=True,
                   help="promoter_g4_damage_gene_table.tsv")
    p.add_argument("--uv-master", required=True,
                   help="results/g4_tss_uv/uv_master_table.tsv")
    p.add_argument("--baseline-tpm", required=True,
                   help="results/g4_tss/baseline_tpm.tsv")
    p.add_argument("--tss-windows", required=True,
                   help="canonical_tss_windows_1kb.bed")
    p.add_argument("--g4chip-peaks", required=True,
                   help="g4_hela_peaks_prepared.tsv (ChIP)")
    p.add_argument("--g4cuttag-peaks", required=True,
                   help="g4_hela_peaks_prepared.tsv (CUT&Tag)")
    p.add_argument("--out-cpd", required=True)
    p.add_argument("--out-64pp", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def strip_version(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace(r"\.\d+$", "", regex=True)


def rank_normalize(series: pd.Series) -> pd.Series:
    ranked = series.rank(method="average", na_option="keep")
    n = ranked.notna().sum()
    if n <= 1:
        return ranked
    return (ranked - 1) / (n - 1)


def bedtools_max_g4_signal(
    windows_bed: str,
    peaks_tsv: str,
    assay_name: str,
    log_fn,
) -> pd.Series:
    """Return max source_signal per gene_id from bedtools intersect."""
    peaks = pd.read_csv(peaks_tsv, sep="\t")
    if "source_signal" not in peaks.columns:
        log_fn(f"  WARNING: source_signal missing from {assay_name} peaks")
        return pd.Series(dtype=float, name=f"max_g4_signal_{assay_name}")

    with tempfile.NamedTemporaryFile(suffix=".bed", mode="w", delete=False) as tf:
        peaks_bed_path = tf.name
        for _, row in peaks.iterrows():
            tf.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t"
                     f"{row.get('name', '.')}\t{row['source_signal']}\n")

    cmd = ["bedtools", "intersect", "-a", windows_bed, "-b", peaks_bed_path, "-wo"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        lines = result.stdout.strip().split("\n") if result.stdout.strip() else []
        log_fn(f"  {assay_name}: {len(lines)} overlap records")
    except subprocess.CalledProcessError as e:
        log_fn(f"  WARNING: bedtools failed for {assay_name}: {e.stderr[:200]}")
        Path(peaks_bed_path).unlink(missing_ok=True)
        return pd.Series(dtype=float, name=f"max_g4_signal_{assay_name}")

    Path(peaks_bed_path).unlink(missing_ok=True)

    if not lines or lines == [""]:
        return pd.Series(dtype=float, name=f"max_g4_signal_{assay_name}")

    rows = []
    for line in lines:
        parts = line.split("\t")
        gene_id = parts[3]
        signal = float(parts[-2])
        rows.append((gene_id, signal))

    df = pd.DataFrame(rows, columns=["gene_id", "signal"])
    df["gene_id"] = strip_version(df["gene_id"])
    return df.groupby("gene_id")["signal"].max().rename(f"max_g4_signal_{assay_name}")


def add_response_indicators(df: pd.DataFrame) -> pd.DataFrame:
    for tp in TIMEPOINTS:
        lfc_col = f"lfc_{tp}"
        padj_col = f"padj_{tp}"
        if lfc_col not in df.columns or padj_col not in df.columns:
            continue
        df[f"repressed_{tp}"] = (
            (df[padj_col] < PADJ_THRESH) & (df[lfc_col] > LFC_THRESH)
        )
        df[f"induced_{tp}"] = (
            (df[padj_col] < PADJ_THRESH) & (df[lfc_col] < -LFC_THRESH)
        )
    if "uv_trajectory" in df.columns:
        df["repressed_sustained"] = df["uv_trajectory"] == "repressed_sustained"
    return df


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else sys.stdout

    def log(msg: str) -> None:
        print(msg, file=log_fh, flush=True)

    for path in [args.out_cpd, args.out_64pp]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    log("=== Task 2: Assemble lesion × RNA master table ===")

    # ── Load gene-level lesion burden (DS0 only) ──────────────────────────────
    gene_table = pd.read_csv(args.gene_table, sep="\t")
    ds0 = gene_table[gene_table["lesion_time_min"] == DS0_TIME].copy()
    ds0["gene_id_nv"] = strip_version(ds0["gene_id"])
    log(f"Gene table (DS0): {len(ds0):,} rows, "
        f"products: {ds0['product'].unique().tolist()}")

    # ── Load UV RNA master table ───────────────────────────────────────────────
    uv = pd.read_csv(args.uv_master, sep="\t")
    uv["gene_id_nv"] = strip_version(uv["gene_id"])
    required_rna = [
        "lfc_12", "lfc_30", "lfc_60",
        "padj_12", "padj_30", "padj_60",
        "lrt_padj", "uv_trajectory", "group",
    ]
    missing = [c for c in required_rna if c not in uv.columns]
    if missing:
        log(f"WARNING: UV master missing columns: {missing}")
    log(f"UV master table: {len(uv):,} genes")

    # ── Load baseline TPM ──────────────────────────────────────────────────────
    tpm = pd.read_csv(args.baseline_tpm, sep="\t")
    tpm["gene_id_nv"] = strip_version(tpm.get("gene_id", tpm.iloc[:, 0]))
    tpm_col = "mean_tpm" if "mean_tpm" in tpm.columns else tpm.columns[1]
    tpm = tpm.rename(columns={tpm_col: "mean_tpm_t0"})[["gene_id_nv", "mean_tpm_t0"]]
    log(f"Baseline TPM: {len(tpm):,} genes")

    # ── Compute G4 signal covariates ──────────────────────────────────────────
    log("\nComputing G4 signal covariates via bedtools intersect...")
    chip_sig = bedtools_max_g4_signal(args.tss_windows, args.g4chip_peaks, "chip", log)
    cuttag_sig = bedtools_max_g4_signal(args.tss_windows, args.g4cuttag_peaks, "cuttag", log)

    g4_cov = pd.DataFrame({
        "gene_id_nv": chip_sig.index.union(cuttag_sig.index),
    }).set_index("gene_id_nv")
    g4_cov["max_g4_signal_chip"] = chip_sig
    g4_cov["max_g4_signal_cuttag"] = cuttag_sig
    g4_cov = g4_cov.fillna(0).reset_index()
    g4_cov["max_g4_signal_chip_norm"] = rank_normalize(g4_cov["max_g4_signal_chip"])
    g4_cov["max_g4_signal_cuttag_norm"] = rank_normalize(g4_cov["max_g4_signal_cuttag"])
    g4_cov["max_g4_signal_norm"] = (
        g4_cov["max_g4_signal_chip_norm"] + g4_cov["max_g4_signal_cuttag_norm"]
    ) / 2
    log(f"G4 signal covariates: {len(g4_cov):,} genes with signal")

    # ── Build product-specific master tables ───────────────────────────────────
    log("\n--- Building product-specific master tables ---")

    product_map = {"CPD": args.out_cpd, "64-PP": args.out_64pp}

    for product, out_path in product_map.items():
        log(f"\nProduct: {product}")
        burden = ds0[ds0["product"] == product].copy()
        burden = burden.rename(columns={
            "damage_ds0": "damage_ds0",
            "damage_ds0_z": "damage_ds0_z",
            "damage_mean": "damage_mean_ds0",
            "damage_sum": "damage_sum_ds0",
            "real_rpkm_max": "real_rpkm_max_ds0",
            "real_rpkm_mean": "real_rpkm_mean_ds0",
        })
        n_burden = len(burden)
        log(f"  Burden rows (DS0, {product}): {n_burden:,}")

        # Join RNA master
        uv_join_cols = ["gene_id_nv"]
        for col in uv.columns:
            if col in {"gene_id", "gene_id_nv", "mean_tpm_t0", "log2_tpm_t0"}:
                continue
            if col in burden.columns:
                continue
            uv_join_cols.append(col)
        merged = burden.merge(
            uv[uv_join_cols],
            on="gene_id_nv",
            how="inner",
        )
        log(f"  After RNA join: {len(merged):,} (lost {n_burden - len(merged):,})")

        # Join baseline TPM
        n_pre = len(merged)
        merged = merged.merge(tpm, on="gene_id_nv", how="left")
        merged["log2_tpm_t0"] = np.log2(merged["mean_tpm_t0"].fillna(0) + 1)
        log(f"  After TPM join: {len(merged):,} (lost {n_pre - len(merged):,})")

        # Join G4 covariates
        merged = merged.merge(g4_cov, on="gene_id_nv", how="left")
        merged["max_g4_signal_norm"] = merged["max_g4_signal_norm"].fillna(0)

        # Add response indicators
        merged = add_response_indicators(merged)

        # Filter to G4_TSS as primary analysis set
        g4_tss = merged[merged["group"] == "G4_TSS"].copy()
        log(f"  G4_TSS genes in master table: {len(g4_tss):,}")

        # Keep all groups in the output but flag primary set
        merged["is_g4_tss"] = merged["group"] == "G4_TSS"

        # Compute abs_lfc_60
        if "lfc_60" in merged.columns:
            merged["abs_lfc_60"] = merged["lfc_60"].abs()

        # Drop gene_id_nv helper
        merged = merged.drop(columns=["gene_id_nv"], errors="ignore")

        merged.to_csv(out_path, sep="\t", index=False)
        log(f"  Written: {out_path}")

        # Join loss summary
        log(f"  Join loss summary:")
        log(f"    Burden (DS0): {n_burden:,}")
        log(f"    After RNA join: {len(merged):,} ({len(merged) / n_burden * 100:.1f}%)")

    log("\n=== Task 2 complete ===")
    if args.log:
        log_fh.close()


if __name__ == "__main__":
    main()
