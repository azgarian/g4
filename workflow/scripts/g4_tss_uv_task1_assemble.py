#!/usr/bin/env python3
"""Task 1: Assemble the master UV-response table.

Sign convention for all pairwise DESeq2 contrasts (0_vs_X):
  positive log2FoldChange -> higher at t=0 than post-UV timepoint -> repression
  negative log2FoldChange -> higher at post-UV timepoint than t=0 -> induction
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


GROUP_ORDER = ["G4_TSS", "GC_bg_TSS", "No_overlap"]
LFC_THRESH = 0.5
PADJ_THRESH = 0.05
LRT_PADJ_THRESH = 0.05


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tss-annotation", required=True)
    p.add_argument("--lrt-results", required=True)
    p.add_argument("--pairwise-12", required=True)
    p.add_argument("--pairwise-30", required=True)
    p.add_argument("--pairwise-60", required=True)
    p.add_argument("--baseline-tpm", required=True)
    p.add_argument("--out-table", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def strip_version(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace(r"\.\d+$", "", regex=True)


def assign_direction(lfc: pd.Series, padj: pd.Series) -> pd.Series:
    direction = pd.Series("not_sig", index=lfc.index)
    repressed = (padj < PADJ_THRESH) & (lfc > LFC_THRESH)
    induced = (padj < PADJ_THRESH) & (lfc < -LFC_THRESH)
    direction[repressed] = "repressed"
    direction[induced] = "induced"
    return direction


def first_sig(row: pd.Series, cols: list[str], directions: list[pd.Series]) -> str:
    for tp, d in zip(["12", "30", "60"], directions):
        if row[d] in ("repressed", "induced"):
            return tp
    return "NA"


def last_sig(row: pd.Series, cols: list[str], directions: list[pd.Series]) -> str:
    result = "NA"
    for tp, d in zip(["12", "30", "60"], directions):
        if row[d] in ("repressed", "induced"):
            result = tp
    return result


def assign_trajectory(row: pd.Series) -> str:
    d12, d30, d60 = row["dir_12"], row["dir_30"], row["dir_60"]
    has_rep = any(d == "repressed" for d in [d12, d30, d60])
    has_ind = any(d == "induced" for d in [d12, d30, d60])
    lrt_sig = (not pd.isna(row["lrt_padj"])) and (row["lrt_padj"] < LRT_PADJ_THRESH)

    if has_rep and has_ind:
        return "complex"
    if has_rep:
        if d60 == "repressed":
            return "repressed_sustained"
        else:
            return "repressed_transient"
    if has_ind:
        if d60 == "induced":
            return "induced_sustained"
        else:
            return "induced_transient"
    if lrt_sig:
        return "lrt_only"
    return "not_responsive"


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else sys.stdout

    def log(msg: str) -> None:
        print(msg, file=log_fh, flush=True)

    for outpath in [args.out_table]:
        Path(outpath).parent.mkdir(parents=True, exist_ok=True)
    Path(args.log).parent.mkdir(parents=True, exist_ok=True) if args.log else None

    log("=== Task 1: Assemble master UV-response table ===")
    log(
        "Sign convention: positive log2FoldChange = higher at t=0 than post-UV = repression; "
        "negative = induction"
    )

    # --- Load annotation ---
    anno = pd.read_csv(args.tss_annotation, sep="\t")
    anno["gene_id_nv"] = strip_version(anno["gene_id"])
    log(f"Annotation loaded: {len(anno):,} genes")
    for g in GROUP_ORDER:
        log(f"  {g}: {(anno['group'] == g).sum():,}")

    # --- Load LRT results ---
    lrt = pd.read_csv(args.lrt_results, sep="\t", compression="gzip")
    lrt_cols = {"gene_id_stripped": "gene_id_nv", "padj": "lrt_padj", "stat": "lrt_stat"}
    lrt = lrt.rename(columns=lrt_cols)[["gene_id_nv", "lrt_padj", "lrt_stat"]]
    log(f"LRT results loaded: {len(lrt):,} genes")

    # --- Load pairwise results ---
    def load_pairwise(path: str, suffix: str) -> pd.DataFrame:
        df = pd.read_csv(path, sep="\t", compression="gzip")
        df = df.rename(columns={"gene_id_stripped": "gene_id_nv"})
        return df[["gene_id_nv", "log2FoldChange", "padj"]].rename(
            columns={
                "log2FoldChange": f"lfc_{suffix}",
                "padj": f"padj_{suffix}",
            }
        )

    pw12 = load_pairwise(args.pairwise_12, "12")
    pw30 = load_pairwise(args.pairwise_30, "30")
    pw60 = load_pairwise(args.pairwise_60, "60")
    log(f"Pairwise 12-min: {len(pw12):,} genes")
    log(f"Pairwise 30-min: {len(pw30):,} genes")
    log(f"Pairwise 60-min: {len(pw60):,} genes")

    # --- Load baseline TPM ---
    tpm = pd.read_csv(args.baseline_tpm, sep="\t")
    tpm = tpm.rename(columns={"gene_id": "gene_id_nv", "mean_tpm": "mean_tpm_t0"})
    tpm["gene_id_nv"] = strip_version(tpm["gene_id_nv"])
    tpm = tpm[["gene_id_nv", "mean_tpm_t0"]]
    log(f"Baseline TPM loaded: {len(tpm):,} genes")

    # --- Sequential inner joins with exclusion tracking ---
    merged = anno[["gene_id", "gene_id_nv", "gene_name", "group"]].copy()
    n0 = len(merged)

    merged = merged.merge(lrt, on="gene_id_nv", how="inner")
    log(f"After join with LRT: {len(merged):,} (excluded {n0 - len(merged):,})")

    n1 = len(merged)
    merged = merged.merge(pw12, on="gene_id_nv", how="inner")
    log(f"After join with 0_vs_12: {len(merged):,} (excluded {n1 - len(merged):,})")

    n2 = len(merged)
    merged = merged.merge(pw30, on="gene_id_nv", how="inner")
    log(f"After join with 0_vs_30: {len(merged):,} (excluded {n2 - len(merged):,})")

    n3 = len(merged)
    merged = merged.merge(pw60, on="gene_id_nv", how="inner")
    log(f"After join with 0_vs_60: {len(merged):,} (excluded {n3 - len(merged):,})")

    n4 = len(merged)
    merged = merged.merge(tpm, on="gene_id_nv", how="inner")
    log(f"After join with baseline TPM: {len(merged):,} (excluded {n4 - len(merged):,})")
    log(f"Total excluded from annotation: {n0 - len(merged):,}")

    # --- Baseline expression ---
    merged["log2_tpm_t0"] = np.log2(merged["mean_tpm_t0"] + 1)

    # --- Magnitude metrics ---
    merged["abs_lfc_12"] = merged["lfc_12"].abs()
    merged["abs_lfc_30"] = merged["lfc_30"].abs()
    merged["abs_lfc_60"] = merged["lfc_60"].abs()
    merged["abs_lfc_max"] = merged[["abs_lfc_12", "abs_lfc_30", "abs_lfc_60"]].max(axis=1)

    # --- Direction labels ---
    merged["dir_12"] = assign_direction(merged["lfc_12"], merged["padj_12"])
    merged["dir_30"] = assign_direction(merged["lfc_30"], merged["padj_30"])
    merged["dir_60"] = assign_direction(merged["lfc_60"], merged["padj_60"])

    # --- Trajectory label ---
    merged["uv_trajectory"] = merged.apply(assign_trajectory, axis=1)

    # --- Onset timing ---
    dir_cols = ["dir_12", "dir_30", "dir_60"]

    def _first(row):
        for tp, col in zip(["12", "30", "60"], dir_cols):
            if row[col] in ("repressed", "induced"):
                return tp
        return "NA"

    def _last(row):
        result = "NA"
        for tp, col in zip(["12", "30", "60"], dir_cols):
            if row[col] in ("repressed", "induced"):
                result = tp
        return result

    merged["first_sig_timepoint"] = merged.apply(_first, axis=1)
    merged["last_sig_timepoint"] = merged.apply(_last, axis=1)

    # Drop helper column
    merged = merged.drop(columns=["gene_id_nv"])

    # --- Group size check ---
    log("\nGroup sizes after assembly:")
    for g in GROUP_ORDER:
        log(f"  {g}: {(merged['group'] == g).sum():,}")

    # --- NA rates ---
    log("\nNA rates per column:")
    for col in ["lfc_12", "lfc_30", "lfc_60", "lrt_padj", "mean_tpm_t0"]:
        na_rate = merged[col].isna().mean()
        log(f"  {col}: {na_rate:.3f}")

    # --- Trajectory distribution ---
    log("\nTrajectory distribution:")
    for traj, cnt in merged["uv_trajectory"].value_counts().items():
        log(f"  {traj}: {cnt:,}")

    # --- Write output ---
    merged.to_csv(args.out_table, sep="\t", index=False)
    log(f"\nWritten: {args.out_table}")

    if args.log:
        log_fh.close()


if __name__ == "__main__":
    main()
