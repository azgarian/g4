#!/usr/bin/env python3
"""Summarize kept vs filtered genes per TSS group after RNA-seq prefiltering."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


GROUP_ORDER = ["G4_TSS", "GC_bg_TSS", "No_overlap"]
STATUS_ORDER = ["kept", "filtered"]
STATUS_COLORS = {
    "kept": "#4c78a8",
    "filtered": "#d9d9d9",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tss-annotation", required=True)
    p.add_argument("--gene-expression-by-group", required=True)
    p.add_argument("--analysis-label", default=None)
    p.add_argument("--out-summary-tsv", required=True)
    p.add_argument("--out-barplot-pdf", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def strip_version(values: pd.Series) -> pd.Series:
    return values.astype(str).str.replace(r"\.\d+$", "", regex=True)


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else None

    def log(msg: str) -> None:
        if log_fh:
            print(msg, file=log_fh, flush=True)

    anno = pd.read_csv(args.tss_annotation, sep="\t")
    expr = pd.read_csv(args.gene_expression_by_group, sep="\t")

    anno["gene_id_nv"] = strip_version(anno["gene_id"])
    expr["gene_id_nv"] = strip_version(expr["gene_id"])

    for outpath in [args.out_summary_tsv, args.out_barplot_pdf]:
        Path(outpath).parent.mkdir(parents=True, exist_ok=True)

    rows = []
    expr_ids = set(expr["gene_id_nv"])
    for grp in GROUP_ORDER:
        anno_sub = anno[anno["group"] == grp].copy()
        anno_ids = set(anno_sub["gene_id_nv"])
        kept_n = len(anno_ids & expr_ids)
        annotated_n = len(anno_ids)
        filtered_n = annotated_n - kept_n
        kept_pct = 100 * kept_n / annotated_n if annotated_n else np.nan
        filtered_pct = 100 * filtered_n / annotated_n if annotated_n else np.nan
        rows.append({
            "group": grp,
            "annotated_n": annotated_n,
            "kept_n": kept_n,
            "filtered_n": filtered_n,
            "kept_pct": kept_pct,
            "filtered_pct": filtered_pct,
        })
        log(
            f"{grp}: annotated={annotated_n:,}, kept={kept_n:,} "
            f"({kept_pct:.1f}%), filtered={filtered_n:,} ({filtered_pct:.1f}%)"
        )

    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(args.out_summary_tsv, sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    x = np.arange(len(GROUP_ORDER))
    bottoms = np.zeros(len(GROUP_ORDER))

    for status in STATUS_ORDER:
        values = summary_df[f"{status}_pct"].values
        ax.bar(
            x,
            values,
            bottom=bottoms,
            width=0.7,
            color=STATUS_COLORS[status],
            edgecolor="white",
            linewidth=0.8,
            label=status.capitalize(),
        )
        bottoms += values

    xticklabels = []
    for _, row in summary_df.iterrows():
        xticklabels.append(f"{row['group']}\n(n={int(row['annotated_n']):,})")
    ax.set_xticks(x)
    ax.set_xticklabels(xticklabels)
    ax.set_ylabel("Genes in annotated group (%)")
    title_prefix = f"{args.analysis_label}: " if args.analysis_label else ""
    ax.set_title(f"{title_prefix}RNA-seq prefilter outcome by TSS group")
    ax.set_ylim(0, 100)
    ax.grid(axis="y", color="#d9d9d9", lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False)

    for i, row in summary_df.iterrows():
        for status in STATUS_ORDER:
            pct = float(row[f"{status}_pct"])
            count = int(row[f"{status}_n"])
            label = f"{count:,} ({pct:.1f}%)"
            segment_bottom = 0.0 if status == "kept" else float(row["kept_pct"])
            segment_mid = segment_bottom + pct / 2

            if pct >= 8:
                ax.text(
                    i,
                    segment_mid,
                    label,
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="white" if status == "kept" else "#444444",
                )
            else:
                y = min(99.3, segment_bottom + pct + 1.0)
                va = "bottom" if y < 99.3 else "top"
                ax.text(
                    i,
                    y,
                    label,
                    ha="center",
                    va=va,
                    fontsize=8,
                    color="#444444",
                )

    plt.tight_layout()
    fig.savefig(args.out_barplot_pdf, dpi=150, bbox_inches="tight")
    plt.close(fig)

    if log_fh:
        log_fh.close()


if __name__ == "__main__":
    main()
