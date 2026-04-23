#!/usr/bin/env python3
"""Task 4: Time-course trajectory analysis.

Tests whether temporal dynamics of UV response differ between G4-promoter and
non-G4 genes: early onset, transient vs sustained, trajectory class enrichment.
Primary contrast: G4_TSS vs No_overlap. GC_bg_TSS reported descriptively.
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
TRAJECTORY_CLASSES = [
    "repressed_sustained", "repressed_transient",
    "induced_sustained", "induced_transient",
    "complex", "lrt_only", "not_responsive",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--master-table", required=True)
    p.add_argument("--out-trajectory-counts", required=True)
    p.add_argument("--out-trajectory-fisher", required=True)
    p.add_argument("--out-trajectory-plot", required=True)
    p.add_argument("--out-baseline-expr", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def safe_fisher(a: int, b: int, c: int, d: int):
    try:
        return stats.fisher_exact([[a, b], [c, d]])
    except Exception:
        return np.nan, np.nan


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else sys.stdout

    def log(msg: str) -> None:
        print(msg, file=log_fh, flush=True)

    for p in [args.out_trajectory_counts, args.out_trajectory_fisher,
              args.out_trajectory_plot, args.out_baseline_expr]:
        Path(p).parent.mkdir(parents=True, exist_ok=True)

    log("=== Task 4: Trajectory analysis ===")

    df = pd.read_csv(args.master_table, sep="\t")
    lrt_sig = df[df["lrt_padj"] < LRT_PADJ_THRESH].copy()
    log(f"LRT-significant genes: {len(lrt_sig):,}")

    # --- Trajectory class × group contingency ---
    count_rows = []
    for g in GROUP_ORDER:
        sub = lrt_sig[lrt_sig["group"] == g]
        n_total = len(sub)
        for traj in TRAJECTORY_CLASSES:
            n = (sub["uv_trajectory"] == traj).sum()
            count_rows.append({
                "group": g,
                "trajectory": traj,
                "n": int(n),
                "n_total_lrt_sig": n_total,
                "fraction": n / n_total if n_total > 0 else np.nan,
            })
    counts_df = pd.DataFrame(count_rows)
    counts_df.to_csv(args.out_trajectory_counts, sep="\t", index=False)
    log(f"Written: {args.out_trajectory_counts}")

    # Log contingency
    log("\nTrajectory distribution among LRT-sig genes:")
    pivot = counts_df.pivot(index="trajectory", columns="group", values="n").fillna(0).astype(int)
    log(str(pivot))

    # --- Fisher exact tests (G4_TSS vs No_overlap) ---
    g4 = lrt_sig[lrt_sig["group"] == "G4_TSS"]
    no = lrt_sig[lrt_sig["group"] == "No_overlap"]

    fisher_rows = []

    # Test 1: Among repressed genes, early onset (12) vs later
    rep_g4 = g4[g4["uv_trajectory"].isin(["repressed_sustained", "repressed_transient"])]
    rep_no = no[no["uv_trajectory"].isin(["repressed_sustained", "repressed_transient"])]
    g4_early = (rep_g4["first_sig_timepoint"] == "12").sum()
    g4_late = (rep_g4["first_sig_timepoint"].isin(["30", "60"])).sum()
    no_early = (rep_no["first_sig_timepoint"] == "12").sum()
    no_late = (rep_no["first_sig_timepoint"].isin(["30", "60"])).sum()
    or1, p1 = safe_fisher(g4_early, g4_late, no_early, no_late)
    fisher_rows.append({
        "test": "repressed_early_onset",
        "description": "Among repressed genes: early (12min) vs later onset in G4_TSS vs No_overlap",
        "g4_tss_focal": int(g4_early), "g4_tss_other": int(g4_late),
        "no_overlap_focal": int(no_early), "no_overlap_other": int(no_late),
        "odds_ratio": or1, "pval": p1,
    })

    # Test 2: Among induced genes, early onset (12) vs later
    ind_g4 = g4[g4["uv_trajectory"].isin(["induced_sustained", "induced_transient"])]
    ind_no = no[no["uv_trajectory"].isin(["induced_sustained", "induced_transient"])]
    g4_early2 = (ind_g4["first_sig_timepoint"] == "12").sum()
    g4_late2 = (ind_g4["first_sig_timepoint"].isin(["30", "60"])).sum()
    no_early2 = (ind_no["first_sig_timepoint"] == "12").sum()
    no_late2 = (ind_no["first_sig_timepoint"].isin(["30", "60"])).sum()
    or2, p2 = safe_fisher(g4_early2, g4_late2, no_early2, no_late2)
    fisher_rows.append({
        "test": "induced_early_onset",
        "description": "Among induced genes: early (12min) vs later onset in G4_TSS vs No_overlap",
        "g4_tss_focal": int(g4_early2), "g4_tss_other": int(g4_late2),
        "no_overlap_focal": int(no_early2), "no_overlap_other": int(no_late2),
        "odds_ratio": or2, "pval": p2,
    })

    # Test 3: complex trajectory enrichment
    g4_complex = (g4["uv_trajectory"] == "complex").sum()
    g4_not_complex = len(g4) - g4_complex
    no_complex = (no["uv_trajectory"] == "complex").sum()
    no_not_complex = len(no) - no_complex
    or3, p3 = safe_fisher(g4_complex, g4_not_complex, no_complex, no_not_complex)
    fisher_rows.append({
        "test": "complex_trajectory",
        "description": "Enrichment of complex trajectory in G4_TSS vs No_overlap",
        "g4_tss_focal": int(g4_complex), "g4_tss_other": int(g4_not_complex),
        "no_overlap_focal": int(no_complex), "no_overlap_other": int(no_not_complex),
        "odds_ratio": or3, "pval": p3,
    })

    pvals = [r["pval"] if not np.isnan(r["pval"]) else 1.0 for r in fisher_rows]
    _, padj, _, _ = multipletests(pvals, method="fdr_bh")
    for i, row in enumerate(fisher_rows):
        row["padj_bh"] = padj[i]

    fisher_df = pd.DataFrame(fisher_rows)
    fisher_df.to_csv(args.out_trajectory_fisher, sep="\t", index=False)
    log(f"Written: {args.out_trajectory_fisher}")
    for _, row in fisher_df.iterrows():
        log(f"  {row['test']}: OR={row['odds_ratio']:.3f}, p={row['pval']:.4e}, BH={row['padj_bh']:.4e}")

    # --- Trajectory line plot: median signed lfc per group × timepoint ---
    fig, ax = plt.subplots(figsize=(7, 4))
    tp_vals = [12, 30, 60]
    for g in GROUP_ORDER:
        sub = lrt_sig[lrt_sig["group"] == g]
        medians = [sub[f"lfc_{tp}"].median() for tp in TIMEPOINTS]
        ax.plot(tp_vals, medians, marker="o", label=g,
                color=GROUP_PALETTE[g], linewidth=2)
    ax.axhline(0, color="black", lw=0.8, ls="--")
    ax.set_xlabel("Time post-UV (min)")
    ax.set_ylabel("Median signed log₂FC")
    ax.set_title("UV-response trajectory (LRT-significant genes)")
    ax.set_xticks(tp_vals)
    ax.legend()
    plt.tight_layout()
    fig.savefig(args.out_trajectory_plot, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"Written: {args.out_trajectory_plot}")

    # --- Baseline expression within sustained-repression class ---
    sust_rep = df[df["uv_trajectory"] == "repressed_sustained"].copy()
    log(f"\nSustained-repression class: {len(sust_rep):,} genes")
    expr_rows = []
    g4_rep_expr = sust_rep.loc[sust_rep["group"] == "G4_TSS", "log2_tpm_t0"].dropna()
    no_rep_expr = sust_rep.loc[sust_rep["group"] == "No_overlap", "log2_tpm_t0"].dropna()
    gc_rep_expr = sust_rep.loc[sust_rep["group"] == "GC_bg_TSS", "log2_tpm_t0"].dropna()

    for g, vals in zip(GROUP_ORDER, [g4_rep_expr, gc_rep_expr, no_rep_expr]):
        expr_rows.append({
            "group": g,
            "n": len(vals),
            "mean_log2_tpm_t0": vals.mean() if len(vals) > 0 else np.nan,
            "median_log2_tpm_t0": vals.median() if len(vals) > 0 else np.nan,
        })

    # One-sided Mann-Whitney: G4_TSS more highly expressed than No_overlap?
    if len(g4_rep_expr) > 0 and len(no_rep_expr) > 0:
        mw_stat, mw_p = stats.mannwhitneyu(
            g4_rep_expr, no_rep_expr, alternative="greater"
        )
        log(f"MWU (G4_TSS > No_overlap log2_tpm_t0, sustained-rep): U={mw_stat:.1f}, p={mw_p:.4e}")
        for row in expr_rows:
            if row["group"] == "G4_TSS":
                row["mwu_vs_no_overlap_stat"] = mw_stat
                row["mwu_vs_no_overlap_pval_one_sided"] = mw_p

    expr_df = pd.DataFrame(expr_rows)
    expr_df.to_csv(args.out_baseline_expr, sep="\t", index=False)
    log(f"Written: {args.out_baseline_expr}")

    if args.log:
        log_fh.close()


if __name__ == "__main__":
    main()
