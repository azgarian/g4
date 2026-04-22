#!/usr/bin/env python3
"""Task 6: Measure UV-induced transcriptional response across the three promoter groups.

Sign convention for 0_vs_60:
  positive log2FoldChange  → higher at t=0 than t=60  → repression after UV
  negative log2FoldChange  → higher at t=60 than t=0  → induction after UV
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


GROUP_ORDER = ["G4_TSS", "GC_bg_TSS", "No_overlap"]
GROUP_PALETTE = {
    "G4_TSS": "#d62728",
    "GC_bg_TSS": "#ff7f0e",
    "No_overlap": "#1f77b4",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tss-annotation", required=True)
    p.add_argument("--lrt-results", required=True, help="time_lrt_results.tsv.gz")
    p.add_argument("--pairwise-results", dest="pairwise_results", default=None,
                   help="Pairwise DESeq2 results for the matched RNA-seq comparison")
    p.add_argument("--pairwise-0-vs-60", dest="pairwise_results", default=None,
                   help="Backward-compatible alias for --pairwise-results")
    p.add_argument("--analysis-label", default=None)
    p.add_argument("--contrast-label", default=None,
                   help="Display label for the DESeq2 pairwise comparison, e.g. '0 vs 12 min'")
    p.add_argument("--out-lrt-summary", required=True)
    p.add_argument("--out-fc-stats", required=True)
    p.add_argument("--out-fc-plot", required=True)
    p.add_argument("--out-volcano-plot", required=True)
    p.add_argument("--log", default=None)
    args = p.parse_args()
    if args.pairwise_results is None:
        p.error("one of --pairwise-results or --pairwise-0-vs-60 is required")
    return args


def strip_version(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace(r"\.\d+$", "", regex=True)


def infer_contrast_parts(pairwise_df: pd.DataFrame) -> tuple[str | None, str | None]:
    if {"numerator_timepoint", "denominator_timepoint"}.issubset(pairwise_df.columns):
        numerator = pairwise_df["numerator_timepoint"].dropna()
        denominator = pairwise_df["denominator_timepoint"].dropna()
        if not numerator.empty and not denominator.empty:
            return str(int(numerator.iloc[0])), str(int(denominator.iloc[0]))
    return None, None


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else None

    def log(msg: str) -> None:
        if log_fh:
            print(msg, file=log_fh, flush=True)

    if args.analysis_label:
        log(f"Analysis label: {args.analysis_label}")

    anno = pd.read_csv(args.tss_annotation, sep="\t")
    anno["gene_id_nv"] = strip_version(anno["gene_id"])

    lrt = pd.read_csv(args.lrt_results, sep="\t", compression="gzip")
    pw = pd.read_csv(args.pairwise_results, sep="\t", compression="gzip")

    lrt["gene_id_nv"] = strip_version(lrt["gene_id"])
    pw["gene_id_nv"] = strip_version(pw["gene_id"])
    numerator_tp, denominator_tp = infer_contrast_parts(pw)
    contrast_label = args.contrast_label
    if contrast_label is None:
        if numerator_tp is not None and denominator_tp is not None:
            contrast_label = f"{numerator_tp} vs {denominator_tp} min"
        else:
            contrast_label = "0 vs 60 min"
    sign_label = (
        f"positive log2FC = higher at {numerator_tp} min than {denominator_tp} min"
        if numerator_tp is not None and denominator_tp is not None
        else "positive log2FC = higher in the numerator condition than the denominator condition"
    )

    merged = anno[["gene_id", "gene_id_nv", "gene_name", "group"]].merge(
        lrt[["gene_id_nv", "padj"]].rename(columns={"padj": "lrt_padj"}),
        on="gene_id_nv", how="left"
    ).merge(
        pw[["gene_id_nv", "log2FoldChange", "padj"]].rename(
            columns={"padj": "pw_padj"}
        ),
        on="gene_id_nv", how="left"
    )
    merged = merged.drop(columns=["gene_id_nv"])

    for outpath in [args.out_lrt_summary, args.out_fc_stats,
                    args.out_fc_plot, args.out_volcano_plot]:
        Path(outpath).parent.mkdir(parents=True, exist_ok=True)

    log(f"Merged annotation with DESeq2 results: {len(merged):,} genes")
    log(f"Sign convention: {sign_label}")

    # --- LRT summary per group ---
    lrt_rows = []
    for grp in GROUP_ORDER:
        sub = merged[merged["group"] == grp]
        n = len(sub)
        n_lrt = sub["lrt_padj"].notna().sum()
        n_sig = (sub["lrt_padj"] < 0.05).sum()
        lrt_rows.append({
            "group": grp,
            "n_genes": n,
            "n_tested_lrt": int(n_lrt),
            "n_significant_lrt": int(n_sig),
            "fraction_significant_lrt": n_sig / n_lrt if n_lrt > 0 else np.nan,
        })
    lrt_df = pd.DataFrame(lrt_rows)
    lrt_df.to_csv(args.out_lrt_summary, sep="\t", index=False)
    log(f"Written: {args.out_lrt_summary}")

    # --- Fold-change stats per group ---
    fc_rows = []
    for grp in GROUP_ORDER:
        sub = merged[(merged["group"] == grp) & merged["log2FoldChange"].notna()].copy()
        sub_sig = sub[sub["pw_padj"] < 0.05]
        n = len(sub)
        n_repressed = ((sub_sig["log2FoldChange"] > 0.5)).sum()
        n_induced = ((sub_sig["log2FoldChange"] < -0.5)).sum()

        # compare FC distributions
        fc_vals = sub["log2FoldChange"].values
        fc_rows.append({
            "group": grp,
            "n_genes": n,
            "mean_log2FC": fc_vals.mean() if n > 0 else np.nan,
            "median_log2FC": np.median(fc_vals) if n > 0 else np.nan,
            "n_repressed_UV": int(n_repressed),
            "n_induced_UV": int(n_induced),
            "fraction_repressed": n_repressed / n if n > 0 else np.nan,
            "fraction_induced": n_induced / n if n > 0 else np.nan,
        })
    fc_df = pd.DataFrame(fc_rows)
    fc_df.to_csv(args.out_fc_stats, sep="\t", index=False)
    log(f"Written: {args.out_fc_stats}")

    # Kruskal-Wallis on fold-changes
    fc_by_grp = [
        merged[(merged["group"] == g) & merged["log2FoldChange"].notna()]["log2FoldChange"].values
        for g in GROUP_ORDER
    ]
    if all(len(x) > 0 for x in fc_by_grp):
        kw_stat, kw_p = stats.kruskal(*fc_by_grp)
        log(f"Kruskal-Wallis on {contrast_label} log2FC: H={kw_stat:.4f}, p={kw_p:.4e}")

    # --- Violin + boxplot of fold changes ---
    plot_df = merged[merged["log2FoldChange"].notna()].copy()
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.violinplot(data=plot_df, x="group", y="log2FoldChange",
                   order=GROUP_ORDER, palette=GROUP_PALETTE,
                   inner=None, linewidth=0.8, alpha=0.7, ax=ax)
    sns.boxplot(data=plot_df, x="group", y="log2FoldChange",
                order=GROUP_ORDER, width=0.15, fliersize=1.5,
                linewidth=0.8, color="white", ax=ax)
    ax.axhline(0, color="black", lw=0.8, ls="--")
    ax.set_xlabel("")
    ax.set_ylabel(f"log₂ fold change ({contrast_label})")
    title_prefix = f"{args.analysis_label}: " if args.analysis_label else ""
    ax.set_title(
        f"{title_prefix}UV-induced fold change by promoter group\n"
        f"(positive = higher at {numerator_tp} min)"
        if numerator_tp is not None
        else f"{title_prefix}UV-induced fold change by promoter group"
    )
    plt.tight_layout()
    fig.savefig(args.out_fc_plot, dpi=150)
    plt.close(fig)

    # --- Volcano plot ---
    vdf = merged[merged["log2FoldChange"].notna() & merged["pw_padj"].notna()].copy()
    vdf["neg_log10_padj"] = -np.log10(vdf["pw_padj"].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(7, 5))
    for grp in GROUP_ORDER:
        sub = vdf[vdf["group"] == grp]
        ax.scatter(sub["log2FoldChange"], sub["neg_log10_padj"],
                   s=3, alpha=0.4, color=GROUP_PALETTE[grp], label=grp, rasterized=True)
    ax.axhline(-np.log10(0.05), color="gray", lw=0.8, ls="--")
    ax.axvline(0.5, color="gray", lw=0.5, ls=":")
    ax.axvline(-0.5, color="gray", lw=0.5, ls=":")
    ax.set_xlabel(f"log₂FC ({contrast_label})")
    ax.set_ylabel("-log₁₀(adj. p-value)")
    ax.set_title(f"{title_prefix}UV-response volcano ({contrast_label})")
    ax.legend(markerscale=3, fontsize=8)
    plt.tight_layout()
    fig.savefig(args.out_volcano_plot, dpi=150)
    plt.close(fig)

    if log_fh:
        log_fh.close()


if __name__ == "__main__":
    main()
