#!/usr/bin/env python3
"""Task 2: Compile baseline (t=0) expression summaries for descriptive and inferential use.

Produces:
  - baseline_tpm.tsv            per-gene per-replicate TPM at t00
  - baseline_tpm_summary.tsv    per-gene mean/SD/log2(mean+1) and expression class
  - baseline_normalized_counts.tsv  per-gene normalized counts at t00
  - gene_expression_by_group.tsv    normalized counts joined with TSS group annotation
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tss-annotation", required=True, help="tss_group_annotation.tsv from Task 1")
    p.add_argument("--tpm-matrix", required=True, help="gene_tpm.tsv.gz")
    p.add_argument("--norm-counts-matrix", required=True, help="normalized_counts.tsv.gz")
    p.add_argument("--sample-manifest", default=None,
                   help="Optional sample manifest used to select an explicit RNA-seq timepoint")
    p.add_argument("--rna-timepoint", type=int, default=None,
                   help="RNA-seq timepoint to select from the sample manifest")
    p.add_argument("--analysis-label", default=None,
                   help="Optional label used only for logging context")
    p.add_argument("--out-baseline-tpm", required=True)
    p.add_argument("--out-baseline-tpm-summary", required=True)
    p.add_argument("--out-baseline-norm-counts", required=True)
    p.add_argument("--out-gene-expression-by-group", required=True)
    p.add_argument("--log", default=None)
    args = p.parse_args()
    if (args.sample_manifest is None) != (args.rna_timepoint is None):
        p.error("--sample-manifest and --rna-timepoint must be provided together")
    return args


def t00_columns(df: pd.DataFrame) -> list[str]:
    """Return columns that correspond to t00 samples (contain 't00' or 'time0' etc.)."""
    candidates = [c for c in df.columns if c not in ("gene_id", "gene_name")]
    t00 = [c for c in candidates if "t00" in c or "_00" in c or "00_" in c]
    if not t00:
        # fall back: look for columns with '00' substring in name
        t00 = [c for c in candidates if "00" in c]
    return t00


def manifest_columns(
    df: pd.DataFrame,
    manifest: pd.DataFrame,
    rna_timepoint: int,
) -> tuple[list[str], list[str], list[str]]:
    """Return requested, present, and missing sample columns for an RNA-seq timepoint."""
    time_col = "timepoint" if "timepoint" in manifest.columns else "timepoint_label"
    sample_ids = manifest.loc[
        manifest[time_col].astype(str) == str(rna_timepoint), "sample_id"
    ].astype(str).tolist()
    present = [sample_id for sample_id in sample_ids if sample_id in df.columns]
    missing = [sample_id for sample_id in sample_ids if sample_id not in df.columns]
    return sample_ids, present, missing


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else sys.stdout

    def log(msg: str) -> None:
        print(msg, file=log_fh, flush=True)

    if args.analysis_label:
        log(f"Analysis label: {args.analysis_label}")

    anno = pd.read_csv(args.tss_annotation, sep="\t")
    log(f"TSS annotation: {len(anno):,} genes")

    manifest = None
    if args.sample_manifest:
        manifest = pd.read_csv(args.sample_manifest, sep="\t")
        log(
            f"Selecting RNA-seq samples for timepoint {args.rna_timepoint} using "
            f"{args.sample_manifest}"
        )

    # --- TPM matrix ---
    tpm = pd.read_csv(args.tpm_matrix, sep="\t", compression="gzip")
    log(f"TPM matrix columns: {list(tpm.columns[:8])}{'...' if len(tpm.columns) > 8 else ''}")

    if manifest is not None:
        requested_tpm, t00_cols, missing_tpm = manifest_columns(tpm, manifest, args.rna_timepoint)
        log(f"Requested TPM samples for RNA timepoint {args.rna_timepoint}: {requested_tpm}")
        if missing_tpm:
            log(f"WARNING: TPM samples missing from matrix: {missing_tpm}")
    else:
        t00_cols = t00_columns(tpm)
        log(f"t00 TPM columns identified: {t00_cols}")
        if not t00_cols:
            log("WARNING: no t00 columns found in TPM matrix; using all non-ID columns")
            t00_cols = [c for c in tpm.columns if c not in ("gene_id", "gene_name")]
    if not t00_cols:
        raise ValueError("No TPM sample columns were selected for expression summarization")

    id_col = "gene_id" if "gene_id" in tpm.columns else tpm.columns[0]
    baseline_tpm = tpm[[id_col] + t00_cols].copy()
    baseline_tpm = baseline_tpm.rename(columns={id_col: "gene_id"})

    # strip version from gene_id if present
    baseline_tpm["gene_id"] = baseline_tpm["gene_id"].str.replace(r"\.\d+$", "", regex=True)
    anno["gene_id_noversion"] = anno["gene_id"].str.replace(r"\.\d+$", "", regex=True)

    tpm_vals = baseline_tpm[t00_cols].values.astype(float)
    baseline_tpm["mean_tpm"] = tpm_vals.mean(axis=1)
    baseline_tpm["sd_tpm"] = tpm_vals.std(axis=1, ddof=1)
    baseline_tpm["log2_mean_tpm_p1"] = np.log2(baseline_tpm["mean_tpm"] + 1)

    # expression class: silent if all t00 TPM = 0, else quartile of log2(mean+1)
    all_zero = (tpm_vals == 0).all(axis=1)
    expressed = baseline_tpm[~all_zero].copy()
    q_labels = ["Q1", "Q2", "Q3", "Q4"]
    quartiles = pd.qcut(expressed["log2_mean_tpm_p1"], q=4, labels=q_labels, duplicates="drop")
    baseline_tpm["expression_class"] = "silent"
    baseline_tpm.loc[~all_zero, "expression_class"] = quartiles.values

    # join with annotation
    merged_tpm = anno[["gene_id", "gene_id_noversion", "gene_name", "group"]].merge(
        baseline_tpm, left_on="gene_id_noversion", right_on="gene_id", how="inner",
        suffixes=("_anno", "")
    )
    # keep annotation gene_id
    if "gene_id_anno" in merged_tpm.columns:
        merged_tpm = merged_tpm.drop(columns=["gene_id_anno"])
    merged_tpm = merged_tpm.drop(columns=["gene_id_noversion"], errors="ignore")

    for outdir in [args.out_baseline_tpm, args.out_baseline_tpm_summary,
                   args.out_baseline_norm_counts, args.out_gene_expression_by_group]:
        Path(outdir).parent.mkdir(parents=True, exist_ok=True)

    merged_tpm.to_csv(args.out_baseline_tpm, sep="\t", index=False)
    log(f"baseline_tpm.tsv: {len(merged_tpm):,} genes")

    summary_cols = ["gene_id", "gene_name", "group", "mean_tpm", "sd_tpm",
                    "log2_mean_tpm_p1", "expression_class"]
    merged_tpm[summary_cols].to_csv(args.out_baseline_tpm_summary, sep="\t", index=False)

    # --- Normalized counts matrix ---
    nc = pd.read_csv(args.norm_counts_matrix, sep="\t", compression="gzip")
    log(f"Normalized counts matrix columns: {list(nc.columns[:8])}{'...' if len(nc.columns) > 8 else ''}")

    if manifest is not None:
        requested_nc, t00_nc, missing_nc = manifest_columns(nc, manifest, args.rna_timepoint)
        log(f"Requested normalized-count samples for RNA timepoint {args.rna_timepoint}: {requested_nc}")
        if missing_nc:
            log(f"WARNING: normalized-count samples missing from matrix: {missing_nc}")
    else:
        t00_nc = t00_columns(nc)
        log(f"t00 normalized count columns identified: {t00_nc}")
        if not t00_nc:
            t00_nc = [c for c in nc.columns if c not in ("gene_id", "gene_name")]
    if not t00_nc:
        raise ValueError("No normalized-count sample columns were selected for expression summarization")

    id_col_nc = "gene_id" if "gene_id" in nc.columns else nc.columns[0]
    baseline_nc = nc[[id_col_nc] + t00_nc].copy()
    baseline_nc = baseline_nc.rename(columns={id_col_nc: "gene_id"})
    baseline_nc["gene_id"] = baseline_nc["gene_id"].str.replace(r"\.\d+$", "", regex=True)

    nc_vals = baseline_nc[t00_nc].values.astype(float)
    baseline_nc["mean_normalized_count"] = nc_vals.mean(axis=1)
    baseline_nc["mean_log2"] = np.log2(baseline_nc["mean_normalized_count"] + 1)

    merged_nc = anno[["gene_id", "gene_id_noversion", "gene_name", "group"]].merge(
        baseline_nc, left_on="gene_id_noversion", right_on="gene_id", how="inner",
        suffixes=("_anno", "")
    )
    if "gene_id_anno" in merged_nc.columns:
        merged_nc = merged_nc.drop(columns=["gene_id_anno"])
    merged_nc = merged_nc.drop(columns=["gene_id_noversion"], errors="ignore")

    merged_nc.to_csv(args.out_baseline_norm_counts, sep="\t", index=False)
    log(f"baseline_normalized_counts.tsv: {len(merged_nc):,} genes")

    # gene_expression_by_group: merge tpm summary + normalized counts
    tpm_sub = merged_tpm[["gene_id", "gene_name", "group", "mean_tpm", "log2_mean_tpm_p1",
                           "expression_class"]].copy()
    nc_sub = merged_nc[["gene_id", "mean_normalized_count", "mean_log2"]].copy()
    by_group = tpm_sub.merge(nc_sub, on="gene_id", how="inner")
    by_group.to_csv(args.out_gene_expression_by_group, sep="\t", index=False)
    log(f"gene_expression_by_group.tsv: {len(by_group):,} genes")

    # summary statistics
    timepoint_label = (
        f"RNA timepoint {args.rna_timepoint}" if args.rna_timepoint is not None else "baseline"
    )
    log("\nSummary notes:")
    log("  TPM-derived quartiles and 'silent' class used for descriptive stratification and Task 7.")
    log("  mean_log2 from normalized counts is the inferential metric used in Tasks 3-5.")
    log(f"\n  Expressed genes for {timepoint_label}: {(~all_zero).sum():,}")
    log(f"  Silent genes for {timepoint_label}: {all_zero.sum():,}")

    ec = merged_tpm["expression_class"].value_counts()
    for cls in ["silent", "Q1", "Q2", "Q3", "Q4"]:
        log(f"  {cls}: {ec.get(cls, 0):,}")

    log("\n  Group-wise medians (mean_log2 from normalized counts):")
    for grp in ["G4_TSS", "GC_bg_TSS", "No_overlap"]:
        sub = by_group[by_group["group"] == grp]["mean_log2"]
        log(f"    {grp}: median={sub.median():.3f}, n={len(sub):,}")

    log("\n  Group-wise medians (mean_tpm):")
    for grp in ["G4_TSS", "GC_bg_TSS", "No_overlap"]:
        sub = by_group[by_group["group"] == grp]["mean_tpm"]
        log(f"    {grp}: median={sub.median():.3f}, n={len(sub):,}")

    if args.log:
        log_fh.close()


if __name__ == "__main__":
    main()
