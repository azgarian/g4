#!/usr/bin/env python3
"""Task 6: Baseline-expression-matched comparison of UV fold-change magnitude.

Tests whether G4-TSS vs No_overlap magnitude differences persist after
controlling for baseline expression level (log2_tpm_t0).
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
import statsmodels.formula.api as smf


GROUP_PALETTE = {
    "G4_TSS": "#1f77b4",
    "GC_bg_TSS": "#ff7f0e",
    "No_overlap": "#7f7f7f",
}
TIMEPOINTS = ["12", "30", "60"]
N_BINS = 10


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--master-table", required=True)
    p.add_argument("--out-wilcoxon", required=True)
    p.add_argument("--out-lm", required=True)
    p.add_argument("--out-bin-plot", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else sys.stdout

    def log(msg: str) -> None:
        print(msg, file=log_fh, flush=True)

    for p in [args.out_wilcoxon, args.out_lm, args.out_bin_plot]:
        Path(p).parent.mkdir(parents=True, exist_ok=True)

    log("=== Task 6: Baseline-expression-matched comparison ===")

    df = pd.read_csv(args.master_table, sep="\t")
    expressed = df[df["mean_tpm_t0"] > 0].copy()
    log(f"Expressed genes: {len(expressed):,}")

    # Bin into 10 equal-size expression bins on log2_tpm_t0
    expressed["expr_bin"] = pd.qcut(
        expressed["log2_tpm_t0"], q=N_BINS, labels=False, duplicates="drop"
    )
    n_bins_actual = expressed["expr_bin"].nunique()
    log(f"Expression bins created: {n_bins_actual}")

    # Primary dataset: G4_TSS and No_overlap only
    primary = expressed[expressed["group"].isin(["G4_TSS", "No_overlap"])].copy()

    # --- Within-bin Wilcoxon on |lfc_60| ---
    wil_rows = []
    for bin_idx in sorted(expressed["expr_bin"].dropna().unique()):
        bin_data = primary[primary["expr_bin"] == bin_idx]
        g4_vals = bin_data.loc[bin_data["group"] == "G4_TSS", "abs_lfc_60"].dropna().values
        no_vals = bin_data.loc[bin_data["group"] == "No_overlap", "abs_lfc_60"].dropna().values
        if len(g4_vals) > 0 and len(no_vals) > 0:
            u_stat, mw_p = stats.mannwhitneyu(g4_vals, no_vals, alternative="two-sided")
            direction = "G4_TSS_higher" if g4_vals.mean() > no_vals.mean() else "No_overlap_higher"
        else:
            u_stat, mw_p, direction = np.nan, np.nan, "insufficient_data"

        bin_expr_range = expressed.loc[expressed["expr_bin"] == bin_idx, "log2_tpm_t0"]
        wil_rows.append({
            "expr_bin": int(bin_idx),
            "log2_tpm_range_lo": bin_expr_range.min(),
            "log2_tpm_range_hi": bin_expr_range.max(),
            "n_g4_tss": int(len(g4_vals)),
            "n_no_overlap": int(len(no_vals)),
            "mean_abs_lfc60_g4_tss": float(g4_vals.mean()) if len(g4_vals) > 0 else np.nan,
            "mean_abs_lfc60_no_overlap": float(no_vals.mean()) if len(no_vals) > 0 else np.nan,
            "mwu_stat": u_stat,
            "mwu_pval": mw_p,
            "direction": direction,
        })

    pvals = [r["mwu_pval"] if not np.isnan(r["mwu_pval"]) else 1.0 for r in wil_rows]
    _, padj, _, _ = multipletests(pvals, method="fdr_bh")
    for i, row in enumerate(wil_rows):
        row["mwu_padj_bh"] = padj[i]

    wil_df = pd.DataFrame(wil_rows)
    wil_df.to_csv(args.out_wilcoxon, sep="\t", index=False)
    log(f"Written: {args.out_wilcoxon}")

    n_sig_g4_higher = ((wil_df["mwu_padj_bh"] < 0.05) & (wil_df["direction"] == "G4_TSS_higher")).sum()
    n_sig_no_higher = ((wil_df["mwu_padj_bh"] < 0.05) & (wil_df["direction"] == "No_overlap_higher")).sum()
    log(f"Bins with G4_TSS significantly higher |lfc_60|: {n_sig_g4_higher}/{n_bins_actual}")
    log(f"Bins with No_overlap significantly higher |lfc_60|: {n_sig_no_higher}/{n_bins_actual}")

    # --- Linear models ---
    lm_rows = []

    for tp in TIMEPOINTS:
        abs_col = f"abs_lfc_{tp}"
        model_df = primary[[abs_col, "log2_tpm_t0", "group"]].dropna().copy()
        model_df["group"] = pd.Categorical(model_df["group"], categories=["No_overlap", "G4_TSS"])

        # Baseline-only model
        try:
            lm_base = smf.ols(f"{abs_col} ~ log2_tpm_t0", data=model_df).fit()
            r2_base = lm_base.rsquared_adj
        except Exception as e:
            log(f"  Baseline model failed for {tp}: {e}")
            r2_base = np.nan

        # Extended model with group
        try:
            lm_ext = smf.ols(f"{abs_col} ~ log2_tpm_t0 + C(group, Treatment('No_overlap'))",
                             data=model_df).fit()
            r2_ext = lm_ext.rsquared_adj
            g4_coef = lm_ext.params.get("C(group, Treatment('No_overlap'))[T.G4_TSS]", np.nan)
            g4_se = lm_ext.bse.get("C(group, Treatment('No_overlap'))[T.G4_TSS]", np.nan)
            g4_pval = lm_ext.pvalues.get("C(group, Treatment('No_overlap'))[T.G4_TSS]", np.nan)
            log(f"  {tp} min: G4_TSS coef={g4_coef:.4f}, SE={g4_se:.4f}, p={g4_pval:.4e}, R2_base={r2_base:.4f}, R2_ext={r2_ext:.4f}")
        except Exception as e:
            log(f"  Extended model failed for {tp}: {e}")
            r2_ext, g4_coef, g4_se, g4_pval = np.nan, np.nan, np.nan, np.nan

        lm_rows.append({
            "timepoint": tp,
            "n_genes": len(model_df),
            "g4_tss_coefficient": g4_coef,
            "g4_tss_se": g4_se,
            "g4_tss_pval": g4_pval,
            "r2_adj_baseline_only": r2_base,
            "r2_adj_extended": r2_ext,
            "delta_r2_adj": (r2_ext - r2_base) if not (np.isnan(r2_ext) or np.isnan(r2_base)) else np.nan,
        })

    lm_df = pd.DataFrame(lm_rows)
    lm_df.to_csv(args.out_lm, sep="\t", index=False)
    log(f"Written: {args.out_lm}")

    # --- Bin-level plot: median |lfc_60| per expression bin by group ---
    fig, ax = plt.subplots(figsize=(9, 5))
    for g in ["G4_TSS", "No_overlap"]:
        sub = primary[primary["group"] == g].copy()
        bin_medians = sub.groupby("expr_bin")["abs_lfc_60"].median()
        bin_centers = sub.groupby("expr_bin")["log2_tpm_t0"].median()
        ax.plot(
            bin_centers.values, bin_medians.values,
            marker="o", label=g, color=GROUP_PALETTE[g], linewidth=2, markersize=5,
        )

    ax.set_xlabel("Median log₂(TPM+1) at t=0 (expression bin)")
    ax.set_ylabel("Median |log₂FC| at 60 min")
    ax.set_title("UV magnitude vs baseline expression, by promoter group")
    ax.legend()
    plt.tight_layout()
    fig.savefig(args.out_bin_plot, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"Written: {args.out_bin_plot}")

    if args.log:
        log_fh.close()


if __name__ == "__main__":
    main()
