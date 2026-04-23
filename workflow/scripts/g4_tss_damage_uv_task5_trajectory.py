#!/usr/bin/env python3
"""Task 5: Trajectory-level and sustained-repression analysis.

Tests whether higher promoter-G4 lesion burden is specifically enriched in
genes with sustained repression trajectories.

Kruskal-Wallis across all trajectory classes, pairwise Mann-Whitney tests,
violin plots, expression-stratified comparison, and logistic model restricted
to repressed_sustained vs not_responsive.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf


TRAJECTORY_ORDER = [
    "repressed_sustained", "repressed_transient",
    "induced_sustained", "induced_transient",
    "lrt_only", "not_responsive", "complex",
]
TRAJ_PALETTE = {
    "repressed_sustained": "#d62728",
    "repressed_transient": "#ff9896",
    "induced_sustained": "#2ca02c",
    "induced_transient": "#98df8a",
    "lrt_only": "#ff7f0e",
    "not_responsive": "#aec7e8",
    "complex": "#9467bd",
}
PAIRWISE_TARGETS = [
    ("repressed_sustained", "not_responsive"),
    ("repressed_sustained", "induced_sustained"),
    ("repressed_sustained", "lrt_only"),
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cpd-table", required=True)
    p.add_argument("--pp64-table", required=True)
    p.add_argument("--out-cpd-kruskal", required=True)
    p.add_argument("--out-64pp-kruskal", required=True)
    p.add_argument("--out-cpd-pairwise", required=True)
    p.add_argument("--out-64pp-pairwise", required=True)
    p.add_argument("--out-cpd-violin", required=True)
    p.add_argument("--out-64pp-violin", required=True)
    p.add_argument("--out-cpd-logistic", required=True)
    p.add_argument("--out-64pp-logistic", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def standardize(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    df = df.copy()
    for col in cols:
        if col in df.columns:
            mu, sd = df[col].mean(), df[col].std()
            df[col] = (df[col] - mu) / sd if sd > 0 else 0.0
    return df


def kruskal_trajectory(df: pd.DataFrame, burden_col: str = "damage_ds0") -> pd.DataFrame:
    groups = {}
    for traj in TRAJECTORY_ORDER:
        sub = df[df["uv_trajectory"] == traj][burden_col].dropna().values
        if len(sub) >= 3:
            groups[traj] = sub

    if len(groups) < 2:
        return pd.DataFrame()

    kw_stat, kw_p = stats.kruskal(*groups.values())
    rows = []
    for traj, vals in groups.items():
        rows.append({
            "trajectory": traj,
            "n": len(vals),
            "median_burden": np.median(vals),
            "mean_burden": np.mean(vals),
        })
    result = pd.DataFrame(rows).sort_values("median_burden", ascending=False)
    result["kruskal_stat"] = kw_stat
    result["kruskal_pvalue"] = kw_p
    return result


def pairwise_tests(df: pd.DataFrame, burden_col: str = "damage_ds0") -> pd.DataFrame:
    rows = []
    for g1, g2 in PAIRWISE_TARGETS:
        x = df[df["uv_trajectory"] == g1][burden_col].dropna().values
        y = df[df["uv_trajectory"] == g2][burden_col].dropna().values
        if len(x) < 3 or len(y) < 3:
            continue
        # One-sided: repressed_sustained has higher burden
        u_stat, mw_p = stats.mannwhitneyu(x, y, alternative="greater")
        rb = 1 - (2 * u_stat) / (len(x) * len(y))
        rows.append({
            "group1": g1, "group2": g2,
            "n1": len(x), "n2": len(y),
            "median_burden_g1": np.median(x),
            "median_burden_g2": np.median(y),
            "mann_whitney_u": u_stat,
            "pvalue_raw": mw_p,
            "rank_biserial": rb,
        })
    df_out = pd.DataFrame(rows)
    if not df_out.empty and "pvalue_raw" in df_out.columns:
        valid = df_out["pvalue_raw"].notna()
        if valid.sum() > 0:
            _, padj, _, _ = multipletests(df_out.loc[valid, "pvalue_raw"], method="fdr_bh")
            df_out.loc[valid, "pvalue_bh"] = padj
    return df_out


def plot_violin(df: pd.DataFrame, product: str, out_path: str) -> None:
    if "uv_trajectory" not in df.columns or "damage_ds0" not in df.columns:
        return
    present = [t for t in TRAJECTORY_ORDER if t in df["uv_trajectory"].values]
    plot_df = df[df["uv_trajectory"].isin(present)][["uv_trajectory", "damage_ds0"]].dropna()
    if plot_df.empty:
        return

    medians = plot_df.groupby("uv_trajectory")["damage_ds0"].median()
    order = [t for t in TRAJECTORY_ORDER if t in medians.index]
    order.sort(key=lambda t: -medians.get(t, 0))

    fig, ax = plt.subplots(figsize=(max(8, len(order) * 1.5), 5))
    palette = [TRAJ_PALETTE.get(t, "grey") for t in order]
    sns.violinplot(
        data=plot_df, x="uv_trajectory", y="damage_ds0", order=order,
        palette=palette, inner="box", ax=ax, linewidth=0.8,
    )
    counts = plot_df.groupby("uv_trajectory").size()
    ax.set_xticklabels(
        [f"{t}\n(n={counts.get(t, 0)})" for t in order], rotation=30, ha="right", fontsize=8
    )
    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
    ax.set_xlabel("")
    ax.set_ylabel("damage_ds0 (log2 real/sim)")
    ax.set_title(f"{product}: Lesion burden by UV trajectory")
    plt.tight_layout()
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()


def expression_stratified_test(df: pd.DataFrame, burden_col: str, log_fn) -> None:
    if "log2_tpm_t0" not in df.columns or "uv_trajectory" not in df.columns:
        return
    med_tpm = df["log2_tpm_t0"].median()
    for label, mask in [("low_expr", df["log2_tpm_t0"] <= med_tpm),
                         ("high_expr", df["log2_tpm_t0"] > med_tpm)]:
        sub = df[mask]
        x = sub[sub["uv_trajectory"] == "repressed_sustained"][burden_col].dropna().values
        y = sub[sub["uv_trajectory"] == "not_responsive"][burden_col].dropna().values
        if len(x) < 3 or len(y) < 3:
            log_fn(f"  {label}: too few genes (repressed_sustained={len(x)}, not_responsive={len(y)})")
            continue
        _, p = stats.mannwhitneyu(x, y, alternative="greater")
        log_fn(f"  {label}: repressed_sustained (n={len(x)}) vs not_responsive (n={len(y)}) "
               f"burden p={p:.3e}")


def fit_logistic_trajectory(df: pd.DataFrame) -> pd.DataFrame:
    sub = df[df["uv_trajectory"].isin(["repressed_sustained", "not_responsive"])].copy()
    sub["repressed_sustained"] = (sub["uv_trajectory"] == "repressed_sustained").astype(int)

    preds = ["damage_ds0_z", "log2_tpm_t0", "max_g4_signal_norm"]
    available = [p for p in preds if p in sub.columns]
    sub = sub[["repressed_sustained"] + available].dropna()

    if len(sub) < 20 or sub["repressed_sustained"].sum() < 5:
        return pd.DataFrame()

    # Standardize if not already
    sub = standardize(sub, available)

    formula = "repressed_sustained ~ " + " + ".join(available)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = smf.logit(formula, data=sub).fit(disp=0)

        rows = []
        for param in model.params.index:
            coef = model.params[param]
            se = model.bse[param]
            pval = model.pvalues[param]
            rows.append({
                "predictor": param,
                "coef": coef,
                "se": se,
                "odds_ratio": np.exp(coef),
                "or_ci_lo": np.exp(coef - 1.96 * se),
                "or_ci_hi": np.exp(coef + 1.96 * se),
                "pvalue": pval,
                "n": len(sub),
                "n_repressed_sustained": int(sub["repressed_sustained"].sum()),
                "aic": model.aic,
                "pseudo_r2": model.prsquared,
            })
        return pd.DataFrame(rows)
    except Exception as exc:
        return pd.DataFrame([{"error": str(exc)}])


def run_product(
    df: pd.DataFrame, product: str,
    out_kruskal: str, out_pairwise: str, out_violin: str, out_logistic: str,
    log_fn,
) -> None:
    log_fn(f"\n=== Product: {product} ===")

    g4 = df[df["group"] == "G4_TSS"].copy() if "group" in df.columns else df.copy()
    g4 = g4.dropna(subset=["damage_ds0"])
    log_fn(f"G4_TSS genes with damage_ds0: {len(g4):,}")

    if "uv_trajectory" not in g4.columns:
        log_fn("  WARNING: no uv_trajectory column")
        for out in [out_kruskal, out_pairwise, out_logistic]:
            pd.DataFrame().to_csv(out, sep="\t", index=False)
        return

    # Trajectory distribution
    traj_counts = g4["uv_trajectory"].value_counts()
    log_fn(f"Trajectory counts:\n{traj_counts.to_string()}")

    # Kruskal-Wallis
    kw_df = kruskal_trajectory(g4)
    if not kw_df.empty:
        kw_stat = kw_df["kruskal_stat"].iloc[0]
        kw_p = kw_df["kruskal_pvalue"].iloc[0]
        log_fn(f"\nKruskal-Wallis: stat={kw_stat:.3f}, p={kw_p:.3e}")
    kw_df.to_csv(out_kruskal, sep="\t", index=False)

    # Pairwise tests
    pw_df = pairwise_tests(g4)
    if not pw_df.empty:
        for _, row in pw_df.iterrows():
            log_fn(f"  {row['group1']} vs {row['group2']}: "
                   f"p_raw={row['pvalue_raw']:.3e}, rb={row.get('rank_biserial', np.nan):.3f}")
    pw_df.to_csv(out_pairwise, sep="\t", index=False)

    # Expression-stratified tests
    log_fn("\nExpression-stratified repressed_sustained vs not_responsive:")
    expression_stratified_test(g4, "damage_ds0", log_fn)

    # Violin plot
    plot_violin(g4, product, out_violin)
    log_fn(f"Violin figure written: {out_violin}")

    # Logistic model
    logistic_df = fit_logistic_trajectory(g4)
    if not logistic_df.empty:
        dam_row = logistic_df[logistic_df.get("predictor", pd.Series()) == "damage_ds0_z"]
        if not dam_row.empty:
            r = dam_row.iloc[0]
            log_fn(f"\nLogistic (repressed_sustained vs not_responsive):")
            log_fn(f"  damage_ds0_z: OR={r.get('odds_ratio', np.nan):.3f}, "
                   f"p={r.get('pvalue', np.nan):.3e}")
    logistic_df.to_csv(out_logistic, sep="\t", index=False)
    log_fn(f"Logistic model written: {out_logistic}")


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else sys.stdout

    def log(msg: str) -> None:
        print(msg, file=log_fh, flush=True)

    for path in [args.out_cpd_kruskal]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    log("=== Task 5: Trajectory-level and sustained-repression analysis ===")

    cpd = pd.read_csv(args.cpd_table, sep="\t")
    pp64 = pd.read_csv(args.pp64_table, sep="\t")

    run_product(
        cpd, "CPD",
        args.out_cpd_kruskal, args.out_cpd_pairwise,
        args.out_cpd_violin, args.out_cpd_logistic,
        log,
    )
    run_product(
        pp64, "64-PP",
        args.out_64pp_kruskal, args.out_64pp_pairwise,
        args.out_64pp_violin, args.out_64pp_logistic,
        log,
    )

    log("\n=== Task 5 complete ===")
    if args.log:
        log_fh.close()


if __name__ == "__main__":
    main()
