#!/usr/bin/env python3
"""Task 3: Compare baseline expression across G4_TSS, GC_bg_TSS, and No_overlap.

Statistical tests:
  - One-sided Mann-Whitney U (G4_TSS > controls), BH-corrected, rank-biserial r
  - Kruskal-Wallis across all three groups
  - Fisher exact tests for expression prevalence (mean_log2 > log2(2))
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gene-expression-by-group", required=True)
    p.add_argument("--out-statistics", required=True)
    p.add_argument("--out-fisher", required=True)
    p.add_argument("--out-summary", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def rank_biserial(x: np.ndarray, y: np.ndarray) -> float:
    """Rank-biserial correlation for Mann-Whitney U."""
    n1, n2 = len(x), len(y)
    u_stat, _ = stats.mannwhitneyu(x, y, alternative="two-sided")
    return 1 - (2 * u_stat) / (n1 * n2)


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else sys.stdout

    def log(msg: str) -> None:
        print(msg, file=log_fh, flush=True)

    df = pd.read_csv(args.gene_expression_by_group, sep="\t")
    for outpath in [args.out_statistics, args.out_fisher, args.out_summary]:
        Path(outpath).parent.mkdir(parents=True, exist_ok=True)

    # split into groups; exclude mean_log2 == 0 for hypothesis tests
    g4 = df[df["group"] == "G4_TSS"]["mean_log2"].values
    gc = df[df["group"] == "GC_bg_TSS"]["mean_log2"].values
    no = df[df["group"] == "No_overlap"]["mean_log2"].values

    g4_expr = g4[g4 > 0]
    gc_expr = gc[gc > 0]
    no_expr = no[no > 0]

    log(f"Expressed genes (mean_log2 > 0): G4_TSS={len(g4_expr)}, "
        f"GC_bg_TSS={len(gc_expr)}, No_overlap={len(no_expr)}")

    # --- Mann-Whitney U (one-sided: G4_TSS > control) ---
    u1, p1 = stats.mannwhitneyu(g4_expr, gc_expr, alternative="greater")
    u2, p2 = stats.mannwhitneyu(g4_expr, no_expr, alternative="greater")
    rr1 = rank_biserial(g4_expr, gc_expr)
    rr2 = rank_biserial(g4_expr, no_expr)

    _, p_adj, _, _ = multipletests([p1, p2], method="fdr_bh")

    # Kruskal-Wallis
    kw_stat, kw_p = stats.kruskal(g4_expr, gc_expr, no_expr)
    log(f"Kruskal-Wallis: H={kw_stat:.4f}, p={kw_p:.4e}")

    stat_rows = [
        {
            "comparison": "G4_TSS_vs_GC_bg_TSS",
            "test": "Mann-Whitney U (one-sided)",
            "U_stat": u1,
            "p_value": p1,
            "p_adj_BH": p_adj[0],
            "rank_biserial_r": rr1,
            "n_g4": len(g4_expr),
            "n_ctrl": len(gc_expr),
        },
        {
            "comparison": "G4_TSS_vs_No_overlap",
            "test": "Mann-Whitney U (one-sided)",
            "U_stat": u2,
            "p_value": p2,
            "p_adj_BH": p_adj[1],
            "rank_biserial_r": rr2,
            "n_g4": len(g4_expr),
            "n_ctrl": len(no_expr),
        },
        {
            "comparison": "all_groups",
            "test": "Kruskal-Wallis",
            "U_stat": kw_stat,
            "p_value": kw_p,
            "p_adj_BH": np.nan,
            "rank_biserial_r": np.nan,
            "n_g4": len(g4_expr),
            "n_ctrl": len(gc_expr) + len(no_expr),
        },
    ]
    pd.DataFrame(stat_rows).to_csv(args.out_statistics, sep="\t", index=False)
    log(f"Written: {args.out_statistics}")

    # --- Fisher exact: expressed if mean_log2 > log2(2) ---
    threshold = np.log2(2)

    def fisher_row(grp_a_vals, grp_b_vals, grp_a_name, grp_b_name, all_a, all_b):
        exp_a = (grp_a_vals > threshold).sum()
        not_a = len(all_a) - exp_a
        exp_b = (grp_b_vals > threshold).sum()
        not_b = len(all_b) - exp_b
        oddsratio, p_fish = stats.fisher_exact([[exp_a, not_a], [exp_b, not_b]],
                                               alternative="greater")
        return {
            "comparison": f"{grp_a_name}_vs_{grp_b_name}",
            "expressed_a": int(exp_a),
            "total_a": int(len(all_a)),
            "expressed_b": int(exp_b),
            "total_b": int(len(all_b)),
            "odds_ratio": oddsratio,
            "p_value": p_fish,
        }

    g4_all = df[df["group"] == "G4_TSS"]["mean_log2"].values
    gc_all = df[df["group"] == "GC_bg_TSS"]["mean_log2"].values
    no_all = df[df["group"] == "No_overlap"]["mean_log2"].values

    fish_rows = [
        fisher_row(g4_all, gc_all, "G4_TSS", "GC_bg_TSS", g4_all, gc_all),
        fisher_row(g4_all, no_all, "G4_TSS", "No_overlap", g4_all, no_all),
    ]
    fish_df = pd.DataFrame(fish_rows)
    _, fish_padj, _, _ = multipletests(fish_df["p_value"].values, method="fdr_bh")
    fish_df["p_adj_BH"] = fish_padj
    fish_df.to_csv(args.out_fisher, sep="\t", index=False)
    log(f"Written: {args.out_fisher}")

    # --- Descriptive summary ---
    summary_rows = []
    for grp, vals in [("G4_TSS", g4_all), ("GC_bg_TSS", gc_all), ("No_overlap", no_all)]:
        summary_rows.append({
            "group": grp,
            "n": len(vals),
            "mean": vals.mean(),
            "median": np.median(vals),
            "Q1": np.percentile(vals, 25),
            "Q3": np.percentile(vals, 75),
            "n_expressed_log2gt1": int((vals > threshold).sum()),
        })
    pd.DataFrame(summary_rows).to_csv(args.out_summary, sep="\t", index=False)
    log(f"Written: {args.out_summary}")

    if args.log:
        log_fh.close()


if __name__ == "__main__":
    main()
