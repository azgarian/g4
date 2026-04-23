#!/usr/bin/env python3
"""Task 2: Compare UV fold-change magnitude across promoter groups at each timepoint.

Uses the full distribution of |lfc| (not just significant genes).
Primary contrast: G4_TSS vs No_overlap.
GC_bg_TSS treated as descriptive/sensitivity comparator.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests


GROUP_ORDER = ["G4_TSS", "GC_bg_TSS", "No_overlap"]
GROUP_PALETTE = {
    "G4_TSS": "#1f77b4",
    "GC_bg_TSS": "#ff7f0e",
    "No_overlap": "#7f7f7f",
}
TIMEPOINTS = ["12", "30", "60"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--master-table", required=True)
    p.add_argument("--out-stats", required=True)
    p.add_argument("--out-summary", required=True)
    p.add_argument("--out-figure", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def rank_biserial(x: np.ndarray, y: np.ndarray) -> float:
    """Rank-biserial correlation from Mann-Whitney U."""
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return np.nan
    u_stat, _ = stats.mannwhitneyu(x, y, alternative="two-sided")
    return 1 - (2 * u_stat) / (n1 * n2)


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else sys.stdout

    def log(msg: str) -> None:
        print(msg, file=log_fh, flush=True)

    for p in [args.out_stats, args.out_summary, args.out_figure]:
        Path(p).parent.mkdir(parents=True, exist_ok=True)

    log("=== Task 2: UV fold-change magnitude comparison ===")

    df = pd.read_csv(args.master_table, sep="\t")
    # Restrict to expressed genes
    df = df[df["mean_tpm_t0"] > 0].copy()
    log(f"Expressed genes (TPM > 0): {len(df):,}")

    # --- Kruskal-Wallis and pairwise Mann-Whitney ---
    pairs = [("G4_TSS", "No_overlap"), ("G4_TSS", "GC_bg_TSS"), ("GC_bg_TSS", "No_overlap")]
    pairwise_label = ["G4_TSS vs No_overlap", "G4_TSS vs GC_bg_TSS", "GC_bg_TSS vs No_overlap"]

    stats_rows = []
    pvals_all = []

    for tp in TIMEPOINTS:
        col = f"abs_lfc_{tp}"
        grp_data = {g: df.loc[df["group"] == g, col].dropna().values for g in GROUP_ORDER}

        # Kruskal-Wallis
        kw_arrays = [grp_data[g] for g in GROUP_ORDER if len(grp_data[g]) > 0]
        if len(kw_arrays) >= 2:
            kw_stat, kw_p = stats.kruskal(*kw_arrays)
        else:
            kw_stat, kw_p = np.nan, np.nan
        log(f"Kruskal-Wallis |lfc_{tp}|: H={kw_stat:.4f}, p={kw_p:.4e}")

        for (g1, g2), label in zip(pairs, pairwise_label):
            x, y = grp_data[g1], grp_data[g2]
            if len(x) > 0 and len(y) > 0:
                u_stat, mw_p = stats.mannwhitneyu(x, y, alternative="two-sided")
                rb = rank_biserial(x, y)
            else:
                u_stat, mw_p, rb = np.nan, np.nan, np.nan
            pvals_all.append(mw_p)
            stats_rows.append({
                "timepoint": tp,
                "group1": g1,
                "group2": g2,
                "comparison": label,
                "n_group1": len(x),
                "n_group2": len(y),
                "kw_stat": kw_stat,
                "kw_pval": kw_p,
                "mw_u_stat": u_stat,
                "mw_pval": mw_p,
                "rank_biserial_r": rb,
            })

    # BH correction across 9 pairwise tests
    reject, padj, _, _ = multipletests(
        [r["mw_pval"] if not np.isnan(r["mw_pval"]) else 1.0 for r in stats_rows],
        method="fdr_bh",
    )
    for i, row in enumerate(stats_rows):
        row["mw_padj_bh"] = padj[i]

    stats_df = pd.DataFrame(stats_rows)
    stats_df.to_csv(args.out_stats, sep="\t", index=False)
    log(f"Written: {args.out_stats}")

    # --- Descriptive summary ---
    summary_rows = []
    for tp in TIMEPOINTS:
        col = f"abs_lfc_{tp}"
        for g in GROUP_ORDER:
            sub = df.loc[df["group"] == g, col].dropna()
            n = len(sub)
            summary_rows.append({
                "timepoint": tp,
                "group": g,
                "n": n,
                "mean_abs_lfc": sub.mean() if n > 0 else np.nan,
                "median_abs_lfc": sub.median() if n > 0 else np.nan,
                "Q1_abs_lfc": sub.quantile(0.25) if n > 0 else np.nan,
                "Q3_abs_lfc": sub.quantile(0.75) if n > 0 else np.nan,
                "pct_above_1": (sub > 1).mean() if n > 0 else np.nan,
                "pct_above_2": (sub > 2).mean() if n > 0 else np.nan,
            })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(args.out_summary, sep="\t", index=False)
    log(f"Written: {args.out_summary}")

    # --- Multi-panel figure ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=False)
    for ax, tp in zip(axes, TIMEPOINTS):
        col = f"abs_lfc_{tp}"
        plot_df = df[[col, "group"]].dropna().copy()
        plot_df.columns = ["abs_lfc", "group"]

        sns.violinplot(
            data=plot_df, x="group", y="abs_lfc",
            order=GROUP_ORDER,
            palette=GROUP_PALETTE,
            inner=None, linewidth=0.8, alpha=0.7, ax=ax,
        )
        sns.boxplot(
            data=plot_df, x="group", y="abs_lfc",
            order=GROUP_ORDER,
            width=0.15, fliersize=1.5, linewidth=0.8,
            color="white", ax=ax,
        )

        # Annotate primary contrast p-value
        primary = stats_df[(stats_df["timepoint"] == tp) & (stats_df["group1"] == "G4_TSS") & (stats_df["group2"] == "No_overlap")]
        if not primary.empty:
            padj_val = primary["mw_padj_bh"].iloc[0]
            label = f"BH adj. p = {padj_val:.3e}" if not np.isnan(padj_val) else ""
            ax.text(0.5, 0.97, label, transform=ax.transAxes,
                    ha="center", va="top", fontsize=8)

        ax.set_xlabel("")
        ax.set_ylabel("|log₂FC|")
        ax.set_title(f"{tp} min post-UV")
        ax.set_xticklabels(GROUP_ORDER, rotation=20, ha="right", fontsize=8)

    fig.suptitle("UV fold-change magnitude by promoter group", fontsize=12, y=1.01)
    plt.tight_layout()
    fig.savefig(args.out_figure, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"Written: {args.out_figure}")

    if args.log:
        log_fh.close()


if __name__ == "__main__":
    main()
