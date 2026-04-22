#!/usr/bin/env python3
"""Task 4: Violin + boxplot of baseline expression and group-specific TSS BED files.

Produces:
  - expression_violin_by_group.pdf
  - tss_G4_TSS.bed, tss_GC_bg_TSS.bed, tss_No_overlap.bed  (for deepTools)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


GROUP_ORDER = ["G4_TSS", "GC_bg_TSS", "No_overlap"]
GROUP_PALETTE = {
    "G4_TSS": "#d62728",
    "GC_bg_TSS": "#ff7f0e",
    "No_overlap": "#1f77b4",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gene-expression-by-group", required=True)
    p.add_argument("--tss-annotation", required=True)
    p.add_argument("--expression-statistics", required=True,
                   help="expression_group_statistics.tsv for p-value annotations")
    p.add_argument("--out-violin-pdf", required=True)
    p.add_argument("--out-bed-g4-tss", required=True)
    p.add_argument("--out-bed-gc-bg-tss", required=True)
    p.add_argument("--out-bed-no-overlap", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def format_pval(p: float) -> str:
    if p < 0.001:
        return f"p={p:.2e}"
    return f"p={p:.3f}"


def sample_points_per_group(expr: pd.DataFrame, max_points: int = 800) -> pd.DataFrame:
    """Downsample large groups so the point overlay stays readable."""
    sampled = []
    for group in GROUP_ORDER:
        sub = expr[expr["group"] == group]
        if len(sub) > max_points:
            sub = sub.sample(n=max_points, random_state=7)
        sampled.append(sub)
    return pd.concat(sampled, ignore_index=True)


def main() -> None:
    args = parse_args()

    expr = pd.read_csv(args.gene_expression_by_group, sep="\t")
    stats_df = pd.read_csv(args.expression_statistics, sep="\t")
    anno = pd.read_csv(args.tss_annotation, sep="\t")
    expr_points = sample_points_per_group(expr)
    group_counts = expr["group"].value_counts().reindex(GROUP_ORDER).fillna(0).astype(int)

    for outpath in [args.out_violin_pdf, args.out_bed_g4_tss,
                    args.out_bed_gc_bg_tss, args.out_bed_no_overlap]:
        Path(outpath).parent.mkdir(parents=True, exist_ok=True)

    # --- Violin + boxplot ---
    fig, ax = plt.subplots(figsize=(7.2, 5.4))

    sns.violinplot(
        data=expr,
        x="group",
        y="mean_log2",
        order=GROUP_ORDER,
        palette=GROUP_PALETTE,
        inner=None,
        cut=0,
        scale="width",
        linewidth=0.8,
        ax=ax,
    )
    sns.stripplot(
        data=expr_points,
        x="group",
        y="mean_log2",
        order=GROUP_ORDER,
        jitter=0.22,
        size=1.6,
        alpha=0.18,
        color="black",
        ax=ax,
        zorder=1,
    )
    sns.boxplot(
        data=expr,
        x="group",
        y="mean_log2",
        order=GROUP_ORDER,
        width=0.15,
        fliersize=1.5,
        linewidth=0.8,
        color="white",
        boxprops=dict(zorder=2),
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel(r"Baseline expression (log$_2$(mean norm. count + 1))")
    ax.set_title("Baseline expression by TSS promoter group")
    ax.set_xticklabels([f"{group}\n(n={group_counts[group]:,})" for group in GROUP_ORDER])
    ax.grid(axis="y", color="#d9d9d9", lw=0.6)
    ax.set_axisbelow(True)

    # annotate with BH-adjusted p-values from Mann-Whitney
    mw_rows = stats_df[stats_df["test"].str.startswith("Mann")]
    y_max = expr["mean_log2"].max()
    offsets = [y_max + 0.6, y_max + 1.4]
    comparisons = [
        ("G4_TSS_vs_GC_bg_TSS", 0, 1),
        ("G4_TSS_vs_No_overlap", 0, 2),
    ]
    for i, (comp_name, xi, xj) in enumerate(comparisons):
        row = mw_rows[mw_rows["comparison"] == comp_name]
        if row.empty:
            continue
        p_adj = row["p_adj_BH"].values[0]
        y_line = offsets[i]
        ax.plot([xi, xi, xj, xj], [y_line, y_line + 0.1, y_line + 0.1, y_line],
                lw=0.8, color="black")
        ax.text((xi + xj) / 2, y_line + 0.12, format_pval(p_adj),
                ha="center", va="bottom", fontsize=7)
    ax.set_ylim(top=y_max + 2.2)

    plt.tight_layout()
    fig.savefig(args.out_violin_pdf, dpi=150)
    plt.close(fig)

    # --- Group TSS BED files ---
    # anno has gene_id, chrom, tss, strand + group via merge
    anno_expr = anno.merge(expr[["gene_id", "group"]], on="gene_id", how="left",
                           suffixes=("", "_expr"))
    group_col = "group" if "group" in anno.columns else "group_expr"

    group_bed_map = {
        "G4_TSS": args.out_bed_g4_tss,
        "GC_bg_TSS": args.out_bed_gc_bg_tss,
        "No_overlap": args.out_bed_no_overlap,
    }

    for grp, bed_path in group_bed_map.items():
        sub = anno[anno["group"] == grp].copy() if "group" in anno.columns else \
              anno_expr[anno_expr["group_expr"] == grp].copy()
        sub["start"] = sub["tss"]
        sub["end"] = sub["tss"] + 1
        sub["score"] = 0
        sub[["chrom", "start", "end", "gene_id", "score", "strand"]].to_csv(
            bed_path, sep="\t", header=False, index=False
        )

    if args.log:
        import sys
        with open(args.log, "w") as fh:
            print(f"Violin plot: {args.out_violin_pdf}", file=fh)
            for grp, bed_path in group_bed_map.items():
                n = len(anno[anno["group"] == grp]) if "group" in anno.columns else "?"
                print(f"  {grp} BED: {bed_path} ({n} genes)", file=fh)


if __name__ == "__main__":
    main()
