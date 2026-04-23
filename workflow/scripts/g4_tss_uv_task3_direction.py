#!/usr/bin/env python3
"""Task 3: Quantify direction of UV response — repression vs induction asymmetry.

Restricted to LRT-significant genes (lrt_padj < 0.05).
Primary contrast: G4_TSS vs No_overlap. GC_bg_TSS treated descriptively.

Sign convention:
  positive lfc -> higher at t=0 -> repression
  negative lfc -> higher post-UV -> induction
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
LRT_PADJ_THRESH = 0.05
DIRECTION_COLORS = {"repressed": "#d62728", "induced": "#2ca02c", "not_sig": "#aec7e8"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--master-table", required=True)
    p.add_argument("--out-direction-stats", required=True)
    p.add_argument("--out-fisher", required=True)
    p.add_argument("--out-stacked-bar", required=True)
    p.add_argument("--out-density", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def safe_fisher(a: int, b: int, c: int, d: int):
    """2x2 Fisher test: [[a, b], [c, d]]."""
    try:
        table = np.array([[a, b], [c, d]])
        result = stats.fisher_exact(table)
        return result[0], result[1]
    except Exception:
        return np.nan, np.nan


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else sys.stdout

    def log(msg: str) -> None:
        print(msg, file=log_fh, flush=True)

    for p in [args.out_direction_stats, args.out_fisher, args.out_stacked_bar, args.out_density]:
        Path(p).parent.mkdir(parents=True, exist_ok=True)

    log("=== Task 3: Direction asymmetry analysis ===")

    df = pd.read_csv(args.master_table, sep="\t")
    lrt_sig = df[df["lrt_padj"] < LRT_PADJ_THRESH].copy()
    log(f"LRT-significant genes: {len(lrt_sig):,}")
    for g in GROUP_ORDER:
        log(f"  {g}: {(lrt_sig['group'] == g).sum():,}")

    # --- Direction counts per group × timepoint ---
    dir_rows = []
    for tp in TIMEPOINTS:
        dir_col = f"dir_{tp}"
        lfc_col = f"lfc_{tp}"
        for g in GROUP_ORDER:
            sub = lrt_sig[lrt_sig["group"] == g]
            n_sig = len(sub)
            n_rep = (sub[dir_col] == "repressed").sum()
            n_ind = (sub[dir_col] == "induced").sum()
            n_other = n_sig - n_rep - n_ind
            signed_mean = sub[lfc_col].mean()
            dir_rows.append({
                "timepoint": tp,
                "group": g,
                "n_sig": n_sig,
                "n_repressed": int(n_rep),
                "n_induced": int(n_ind),
                "n_other": int(n_other),
                "pct_repressed": n_rep / n_sig if n_sig > 0 else np.nan,
                "pct_induced": n_ind / n_sig if n_sig > 0 else np.nan,
                "repression_ratio": n_rep / n_ind if n_ind > 0 else np.nan,
                "signed_mean_lfc": signed_mean,
            })
    dir_df = pd.DataFrame(dir_rows)
    dir_df.to_csv(args.out_direction_stats, sep="\t", index=False)
    log(f"Written: {args.out_direction_stats}")

    # --- Fisher exact tests ---
    # Per timepoint: repression and induction enrichment in G4_TSS vs No_overlap
    fisher_rows = []
    for tp in TIMEPOINTS:
        dir_col = f"dir_{tp}"
        g4 = lrt_sig[lrt_sig["group"] == "G4_TSS"]
        no = lrt_sig[lrt_sig["group"] == "No_overlap"]
        g4_rep = (g4[dir_col] == "repressed").sum()
        g4_not_rep = len(g4) - g4_rep
        no_rep = (no[dir_col] == "repressed").sum()
        no_not_rep = len(no) - no_rep
        or_rep, p_rep = safe_fisher(g4_rep, g4_not_rep, no_rep, no_not_rep)

        g4_ind = (g4[dir_col] == "induced").sum()
        g4_not_ind = len(g4) - g4_ind
        no_ind = (no[dir_col] == "induced").sum()
        no_not_ind = len(no) - no_ind
        or_ind, p_ind = safe_fisher(g4_ind, g4_not_ind, no_ind, no_not_ind)

        fisher_rows.append({
            "timepoint": tp, "test": "repression",
            "g4_tss_focal": int(g4_rep), "g4_tss_other": int(g4_not_rep),
            "no_overlap_focal": int(no_rep), "no_overlap_other": int(no_not_rep),
            "odds_ratio": or_rep, "pval": p_rep,
        })
        fisher_rows.append({
            "timepoint": tp, "test": "induction",
            "g4_tss_focal": int(g4_ind), "g4_tss_other": int(g4_not_ind),
            "no_overlap_focal": int(no_ind), "no_overlap_other": int(no_not_ind),
            "odds_ratio": or_ind, "pval": p_ind,
        })

    pvals = [r["pval"] if not np.isnan(r["pval"]) else 1.0 for r in fisher_rows]
    _, padj, _, _ = multipletests(pvals, method="fdr_bh")
    for i, row in enumerate(fisher_rows):
        row["padj_bh"] = padj[i]

    fisher_df = pd.DataFrame(fisher_rows)
    fisher_df.to_csv(args.out_fisher, sep="\t", index=False)
    log(f"Written: {args.out_fisher}")
    for _, row in fisher_df.iterrows():
        log(f"  {row['timepoint']} {row['test']}: OR={row['odds_ratio']:.3f}, p={row['pval']:.4e}, BH={row['padj_bh']:.4e}")

    # --- Stacked bar chart ---
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)
    for ax, tp in zip(axes, TIMEPOINTS):
        dir_col = f"dir_{tp}"
        plot_data = []
        for g in GROUP_ORDER:
            sub = lrt_sig[lrt_sig["group"] == g]
            n = len(sub)
            if n == 0:
                plot_data.append({"group": g, "repressed": 0, "induced": 0, "not_sig": 0})
                continue
            plot_data.append({
                "group": g,
                "repressed": (sub[dir_col] == "repressed").sum() / n,
                "induced": (sub[dir_col] == "induced").sum() / n,
                "not_sig": (sub[dir_col] == "not_sig").sum() / n,
            })
        pdata = pd.DataFrame(plot_data).set_index("group")
        bottom_rep = np.zeros(len(GROUP_ORDER))
        bottom_ind = pdata["repressed"].values
        bottom_ns = bottom_ind + pdata["induced"].values
        ax.bar(GROUP_ORDER, pdata["repressed"], color=DIRECTION_COLORS["repressed"],
               label="repressed")
        ax.bar(GROUP_ORDER, pdata["induced"], bottom=bottom_ind,
               color=DIRECTION_COLORS["induced"], label="induced")
        ax.bar(GROUP_ORDER, pdata["not_sig"], bottom=bottom_ns,
               color=DIRECTION_COLORS["not_sig"], label="not_sig at tp")
        ax.set_title(f"{tp} min post-UV")
        ax.set_ylabel("Fraction of LRT-sig genes")
        ax.set_xticklabels(GROUP_ORDER, rotation=20, ha="right", fontsize=8)
        ax.set_ylim(0, 1)
        if ax == axes[0]:
            ax.legend(fontsize=7)
    fig.suptitle("Direction of UV response among LRT-significant genes", fontsize=11)
    plt.tight_layout()
    fig.savefig(args.out_stacked_bar, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"Written: {args.out_stacked_bar}")

    # --- Signed LFC density plots ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=False)
    for ax, tp in zip(axes, TIMEPOINTS):
        lfc_col = f"lfc_{tp}"
        plotted = False
        for g in GROUP_ORDER:
            sub = lrt_sig.loc[lrt_sig["group"] == g, lfc_col].dropna()
            if len(sub) < 5:
                continue
            sub.plot.kde(ax=ax, label=g, color=GROUP_PALETTE[g], alpha=0.7, linewidth=1.5)
            ax.fill_between(
                *_kde_xy(sub, ax), alpha=0.15, color=GROUP_PALETTE[g]
            )
            plotted = True
        for xval in [0, 0.5, -0.5, 1, -1]:
            ax.axvline(xval, color="gray", lw=0.7, ls="--" if xval != 0 else "-", alpha=0.6)
        ax.set_xlabel("log₂FC")
        ax.set_ylabel("Density")
        ax.set_title(f"{tp} min post-UV")
        ax.legend(fontsize=7)
    fig.suptitle("Signed log₂FC distribution (LRT-significant genes)", fontsize=11)
    plt.tight_layout()
    fig.savefig(args.out_density, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"Written: {args.out_density}")

    if args.log:
        log_fh.close()


def _kde_xy(series: pd.Series, ax) -> tuple:
    """Return (x, y) arrays for KDE fill. Returns empty arrays on failure."""
    try:
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(series.dropna())
        xmin, xmax = ax.get_xlim() if ax.get_xlim() != (0, 1) else (series.min() - 1, series.max() + 1)
        x = np.linspace(xmin, xmax, 300)
        return x, kde(x)
    except Exception:
        return np.array([]), np.array([])


if __name__ == "__main__":
    main()
