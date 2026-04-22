#!/usr/bin/env python3
"""Task 9: Build per-gene promoter-G4 × LRT intersection table (pG4_DEG_intersect.tsv).

Joins tss_group_annotation.tsv with the omnibus DESeq2 LRT results and the three
pairwise contrasts (0 vs 12, 0 vs 30, 0 vs 60 min), adds normalized counts across all
time points, and assigns a simple response label to each LRT-significant gene.

Sign convention inherited from DESeq2 (0 as numerator, X min as denominator):
  positive log2FC → higher at t=0 → repressed after UV
  negative log2FC → higher at t=X → induced after UV
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tss-annotation", required=True)
    p.add_argument("--lrt-results", required=True)
    p.add_argument("--norm-counts", required=True)
    p.add_argument("--baseline-tpm", required=True)
    p.add_argument("--pairwise-0-vs-12", required=True)
    p.add_argument("--pairwise-0-vs-30", required=True)
    p.add_argument("--pairwise-0-vs-60", required=True)
    p.add_argument("--lrt-padj-threshold", type=float, default=0.05)
    p.add_argument("--fc-threshold", type=float, default=0.5,
                   help="|log2FC| threshold for response label assignment")
    p.add_argument("--out-intersect", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def strip_version(s: pd.Series) -> pd.Series:
    return s.astype(str).str.replace(r"\.\d+$", "", regex=True)


def assign_response_label(
    row: pd.Series,
    padj_thresh: float = 0.05,
    lfc_thresh: float = 0.5,
) -> str:
    """Classify the UV transcriptional response of one LRT-significant gene.

    Categories (only assigned when lrt_sig is True):
      induced_sustained / repressed_sustained  — sig at 60 min
      induced_transient / repressed_transient  — sig at 12 or 30 min but not 60 min
      induced / repressed                      — sig at any timepoint, consistent direction
      ambiguous                                — mixed directions across timepoints
      sig_lrt_no_pairwise                      — LRT-sig but no pairwise contrast passes
      not_sig                                  — LRT padj >= threshold
    """
    if not row.get("lrt_sig", False):
        return "not_sig"

    sig_tps: dict[str, str] = {}
    for tp in ("12", "30", "60"):
        lfc = row.get(f"lfc_0_vs_{tp}", np.nan)
        padj = row.get(f"padj_0_vs_{tp}", np.nan)
        if pd.notna(lfc) and pd.notna(padj) and padj < padj_thresh:
            if lfc < -lfc_thresh:
                sig_tps[tp] = "induced"
            elif lfc > lfc_thresh:
                sig_tps[tp] = "repressed"

    if not sig_tps:
        return "sig_lrt_no_pairwise"

    directions = set(sig_tps.values())
    if len(directions) > 1:
        return "ambiguous"

    direction = directions.pop()
    if "60" in sig_tps:
        return f"{direction}_sustained"
    return f"{direction}_transient"


def load_pairwise(path: str, label: str) -> pd.DataFrame:
    pw = pd.read_csv(path, sep="\t", compression="gzip")
    pw["gene_id_nv"] = strip_version(pw["gene_id"])
    rename = {
        "log2FoldChange": f"lfc_0_vs_{label}",
        "padj": f"padj_0_vs_{label}",
    }
    pw = pw.rename(columns=rename)
    keep = ["gene_id_nv", f"lfc_0_vs_{label}", f"padj_0_vs_{label}"]
    if f"lfcSE" in pw.columns:
        pw = pw.rename(columns={"lfcSE": f"lfcSE_0_vs_{label}"})
        keep.append(f"lfcSE_0_vs_{label}")
    return pw[keep].copy()


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else None

    def log(msg: str) -> None:
        if log_fh:
            print(msg, file=log_fh, flush=True)

    # ── annotation ────────────────────────────────────────────────────────────
    anno = pd.read_csv(args.tss_annotation, sep="\t")
    anno["gene_id_nv"] = strip_version(anno["gene_id"])
    log(f"Annotation loaded: {len(anno):,} genes")
    for grp in ("G4_TSS", "GC_bg_TSS", "No_overlap"):
        log(f"  {grp}: {(anno['group'] == grp).sum():,}")

    # ── LRT results ───────────────────────────────────────────────────────────
    lrt = pd.read_csv(args.lrt_results, sep="\t", compression="gzip")
    lrt["gene_id_nv"] = strip_version(lrt["gene_id"])
    lrt = lrt.rename(columns={
        "baseMean": "lrt_baseMean",
        "stat":     "lrt_stat",
        "pvalue":   "lrt_pvalue",
        "padj":     "lrt_padj",
    })
    lrt_cols = ["gene_id_nv", "lrt_baseMean", "lrt_stat", "lrt_pvalue", "lrt_padj"]
    lrt_sub = lrt[[c for c in lrt_cols if c in lrt.columns]].copy()
    log(f"LRT results loaded: {len(lrt_sub):,} genes tested")

    # ── pairwise contrasts ────────────────────────────────────────────────────
    pw12 = load_pairwise(args.pairwise_0_vs_12, "12")
    pw30 = load_pairwise(args.pairwise_0_vs_30, "30")
    pw60 = load_pairwise(args.pairwise_0_vs_60, "60")
    log(f"Pairwise contrasts loaded: 0v12={len(pw12):,}  0v30={len(pw30):,}  0v60={len(pw60):,}")

    # ── baseline TPM ─────────────────────────────────────────────────────────
    btpm = pd.read_csv(args.baseline_tpm, sep="\t")
    btpm["gene_id_nv"] = strip_version(btpm["gene_id"])
    btpm_sub = btpm[["gene_id_nv", "mean_tpm"]].copy()
    log(f"Baseline TPM loaded: {len(btpm_sub):,} genes")

    # ── normalized counts (all time points) ──────────────────────────────────
    nc = pd.read_csv(args.norm_counts, sep="\t", compression="gzip")
    nc_id = "gene_id" if "gene_id" in nc.columns else nc.columns[0]
    nc = nc.rename(columns={nc_id: "gene_id"})
    nc["gene_id_nv"] = strip_version(nc["gene_id"])
    nc_sample_cols = [c for c in nc.columns if c not in ("gene_id", "gene_id_nv", "gene_name")]
    nc_sub = nc[["gene_id_nv"] + nc_sample_cols].copy()
    log(f"Normalized counts loaded: {len(nc_sub):,} genes, {len(nc_sample_cols)} samples")

    # ── build merged table ────────────────────────────────────────────────────
    meta_cols = ["gene_id", "gene_id_nv", "gene_name", "chrom", "tss", "strand",
                 "group", "has_g4_tss", "has_gc_bg_tss"]
    merged = anno[[c for c in meta_cols if c in anno.columns]].copy()

    merged = (
        merged
        .merge(lrt_sub,  on="gene_id_nv", how="left")
        .merge(pw12,     on="gene_id_nv", how="left")
        .merge(pw30,     on="gene_id_nv", how="left")
        .merge(pw60,     on="gene_id_nv", how="left")
        .merge(btpm_sub, on="gene_id_nv", how="left")
        .merge(nc_sub,   on="gene_id_nv", how="left")
    )
    merged = merged.drop(columns=["gene_id_nv"])

    # ── LRT significance and response labels ──────────────────────────────────
    merged["lrt_sig"] = merged["lrt_padj"] < args.lrt_padj_threshold
    merged["response_label"] = merged.apply(
        assign_response_label,
        axis=1,
        padj_thresh=args.lrt_padj_threshold,
        lfc_thresh=args.fc_threshold,
    )

    Path(args.out_intersect).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.out_intersect, sep="\t", index=False)
    log(f"Written: {args.out_intersect} ({len(merged):,} genes total)")

    # ── summary ───────────────────────────────────────────────────────────────
    log(f"\nLRT padj threshold: {args.lrt_padj_threshold}")
    for grp in ("G4_TSS", "GC_bg_TSS", "No_overlap"):
        sub = merged[merged["group"] == grp]
        n_sig = sub["lrt_sig"].sum()
        n_tested = sub["lrt_padj"].notna().sum()
        pct = 100 * n_sig / n_tested if n_tested > 0 else float("nan")
        log(f"\n{grp}: {len(sub):,} genes, {n_sig:,}/{n_tested:,} LRT-sig ({pct:.1f}%)")
        lbl_counts = sub["response_label"].value_counts()
        for lbl, cnt in lbl_counts.items():
            log(f"  {lbl}: {cnt:,}")

    if log_fh:
        log_fh.close()


if __name__ == "__main__":
    main()
