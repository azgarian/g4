#!/usr/bin/env python3
"""Task 3: Primary burden-stratified test.

Bins G4_TSS genes into equal-count lesion-burden tertiles and tests whether
high-burden genes show stronger UV repression than low-burden genes.

Sign convention: positive lfc = repression after UV.
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


TIMEPOINTS = ["12", "30", "60"]
TERTILE_ORDER = ["damage_low", "damage_mid", "damage_high"]
TERTILE_PALETTE = {
    "damage_low": "#2196F3",
    "damage_mid": "#FF9800",
    "damage_high": "#F44336",
}
FIGSIZE_VIOLIN = (10, 6)
FIGSIZE_LINE = (6, 5)
FIGSIZE_STACK = (8, 5)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cpd-table", required=True)
    p.add_argument("--pp64-table", required=True)
    p.add_argument("--out-cpd-stats", required=True)
    p.add_argument("--out-64pp-stats", required=True)
    p.add_argument("--out-cpd-summary", required=True)
    p.add_argument("--out-64pp-summary", required=True)
    p.add_argument("--out-cpd-lfc-fig", required=True)
    p.add_argument("--out-64pp-lfc-fig", required=True)
    p.add_argument("--out-cpd-traj-fig", required=True)
    p.add_argument("--out-64pp-traj-fig", required=True)
    p.add_argument("--out-cpd-comp-fig", required=True)
    p.add_argument("--out-64pp-comp-fig", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def assign_tertiles(df: pd.DataFrame, col: str = "damage_ds0") -> pd.DataFrame:
    df = df.copy()
    df["damage_tertile"] = pd.qcut(
        df[col],
        q=3,
        labels=["damage_low", "damage_mid", "damage_high"],
        duplicates="drop",
    )
    return df


def rank_biserial(x: np.ndarray, y: np.ndarray) -> float:
    """Rank-biserial correlation as effect size for Mann-Whitney U."""
    nx, ny = len(x), len(y)
    u, _ = stats.mannwhitneyu(x, y, alternative="two-sided")
    return 1 - (2 * u) / (nx * ny)


def tertile_summary(df: pd.DataFrame, tp: str) -> pd.DataFrame:
    lfc_col = f"lfc_{tp}"
    rep_col = f"repressed_{tp}"
    ind_col = f"induced_{tp}"
    if lfc_col not in df.columns:
        return pd.DataFrame()

    rows = []
    for tert in TERTILE_ORDER:
        sub = df[df["damage_tertile"] == tert][lfc_col].dropna()
        if sub.empty:
            continue
        n = len(sub)
        n_rep = int(df[df["damage_tertile"] == tert][rep_col].sum()) if rep_col in df.columns else np.nan
        n_ind = int(df[df["damage_tertile"] == tert][ind_col].sum()) if ind_col in df.columns else np.nan
        rows.append({
            "timepoint": tp,
            "tertile": tert,
            "n": n,
            "median_lfc": sub.median(),
            "mean_lfc": sub.mean(),
            "median_abs_lfc": sub.abs().median(),
            "pct_repressed": n_rep / n * 100 if n > 0 else np.nan,
            "pct_induced": n_ind / n * 100 if n > 0 else np.nan,
            "pct_lfc_gt_1": (sub > 1).mean() * 100,
            "pct_lfc_gt_2": (sub > 2).mean() * 100,
        })
    return pd.DataFrame(rows)


def run_tests(df: pd.DataFrame, tp: str, log_fn) -> pd.DataFrame:
    lfc_col = f"lfc_{tp}"
    if lfc_col not in df.columns:
        return pd.DataFrame()

    groups = {t: df[df["damage_tertile"] == t][lfc_col].dropna().values
              for t in TERTILE_ORDER}
    if any(len(v) < 3 for v in groups.values()):
        log_fn(f"  t={tp}: skipping (insufficient data per tertile)")
        return pd.DataFrame()

    rows = []

    # Kruskal-Wallis
    kw_stat, kw_p = stats.kruskal(*groups.values())
    rows.append({
        "timepoint": tp, "test": "kruskal_wallis", "comparison": "all_tertiles",
        "statistic": kw_stat, "pvalue_raw": kw_p, "effect_size": np.nan,
    })

    # Mann-Whitney pairwise (one-sided: high > low → repression direction)
    pairs = [
        ("damage_high", "damage_low"),
        ("damage_high", "damage_mid"),
    ]
    for g1, g2 in pairs:
        x, y = groups[g1], groups[g2]
        # One-sided: does high-burden have greater (more positive) lfc?
        u_stat, mw_p = stats.mannwhitneyu(x, y, alternative="greater")
        rb = rank_biserial(x, y)
        rows.append({
            "timepoint": tp,
            "test": "mann_whitney_onesided",
            "comparison": f"{g1}_vs_{g2}",
            "statistic": u_stat,
            "pvalue_raw": mw_p,
            "effect_size": rb,
        })
        # Unsigned analysis: abs(lfc)
        u_abs, mw_abs_p = stats.mannwhitneyu(
            np.abs(x), np.abs(y), alternative="greater"
        )
        rows.append({
            "timepoint": tp,
            "test": "mann_whitney_abs_lfc",
            "comparison": f"{g1}_vs_{g2}",
            "statistic": u_abs,
            "pvalue_raw": mw_abs_p,
            "effect_size": rank_biserial(np.abs(x), np.abs(y)),
        })

    return pd.DataFrame(rows)


def categorical_tests(df: pd.DataFrame, tp: str, log_fn) -> pd.DataFrame:
    rep_col = f"repressed_{tp}"
    ind_col = f"induced_{tp}"
    if rep_col not in df.columns:
        return pd.DataFrame()

    rows = []

    def _response_cat(row):
        if row.get(rep_col, False):
            return "repressed"
        if row.get(ind_col, False):
            return "induced"
        return "not_sig"

    df = df.copy()
    df["_resp"] = df.apply(_response_cat, axis=1)

    # Chi-squared: tertile × response category
    ct = pd.crosstab(df["damage_tertile"], df["_resp"])
    if ct.shape == (3, 3) or ct.shape[0] >= 2:
        chi2, chi_p, dof, _ = stats.chi2_contingency(ct)
        rows.append({
            "timepoint": tp, "test": "chi2_tertile_vs_response",
            "comparison": "all_tertiles_x_direction",
            "statistic": chi2, "pvalue_raw": chi_p, "effect_size": np.nan,
        })

    # Fisher: damage_high enriched for repressed_60 / repressed_sustained
    high_mask = df["damage_tertile"] == "damage_high"
    rest_mask = ~high_mask
    for outcome_col, label in [(rep_col, f"repressed_{tp}"),
                                ("repressed_sustained", "repressed_sustained")]:
        if outcome_col not in df.columns:
            continue
        a = int(df[high_mask][outcome_col].sum())
        b = int(high_mask.sum() - a)
        c = int(df[rest_mask][outcome_col].sum())
        d = int(rest_mask.sum() - c)
        ct2 = [[a, b], [c, d]]
        if min(a, b, c, d) >= 0:
            or_val, fp = stats.fisher_exact(ct2, alternative="greater")
            rows.append({
                "timepoint": tp, "test": "fisher_high_enrichment",
                "comparison": f"damage_high_vs_rest_{label}",
                "statistic": or_val, "pvalue_raw": fp, "effect_size": or_val,
            })
    return pd.DataFrame(rows)


def plot_lfc_by_tertile(df: pd.DataFrame, product: str, out_path: str) -> None:
    fig, axes = plt.subplots(1, 3, figsize=FIGSIZE_VIOLIN, sharey=True)
    for ax, tp in zip(axes, TIMEPOINTS):
        lfc_col = f"lfc_{tp}"
        if lfc_col not in df.columns:
            continue
        plot_df = df[["damage_tertile", lfc_col]].dropna()
        plot_df["damage_tertile"] = pd.Categorical(
            plot_df["damage_tertile"], categories=TERTILE_ORDER
        )
        sns.violinplot(
            data=plot_df, x="damage_tertile", y=lfc_col, order=TERTILE_ORDER,
            palette=TERTILE_PALETTE, inner="box", ax=ax, linewidth=0.8,
        )
        ax.axhline(0, color="black", linewidth=0.6, linestyle="--")
        ax.set_title(f"t = {tp} min")
        ax.set_xlabel("")
        ax.set_ylabel("signed LFC (repression > 0)" if tp == "12" else "")
        ax.set_xticklabels(["Low", "Mid", "High"], fontsize=8)
    fig.suptitle(f"{product}: LFC by lesion-burden tertile", fontsize=11)
    plt.tight_layout()
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()


def plot_median_trajectory(df: pd.DataFrame, product: str, out_path: str) -> None:
    fig, ax = plt.subplots(figsize=FIGSIZE_LINE)
    for tert, color in TERTILE_PALETTE.items():
        sub = df[df["damage_tertile"] == tert]
        medians = []
        for tp in TIMEPOINTS:
            col = f"lfc_{tp}"
            if col in sub.columns:
                medians.append(sub[col].dropna().median())
            else:
                medians.append(np.nan)
        ax.plot(TIMEPOINTS, medians, marker="o", label=tert.replace("damage_", ""),
                color=color, linewidth=1.5)
    ax.axhline(0, color="black", linewidth=0.5, linestyle="--")
    ax.set_xlabel("Time post-UV (min)")
    ax.set_ylabel("Median signed LFC (repression > 0)")
    ax.set_title(f"{product}: Median LFC trajectory by burden tertile")
    ax.legend(title="Burden", fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()


def plot_trajectory_composition(df: pd.DataFrame, product: str, out_path: str) -> None:
    if "uv_trajectory" not in df.columns:
        plt.figure()
        plt.text(0.5, 0.5, "No trajectory column", ha="center", va="center")
        plt.savefig(out_path)
        plt.close()
        return

    traj_cats = [
        "repressed_sustained", "repressed_transient", "not_responsive",
        "induced_sustained", "induced_transient", "lrt_only", "complex",
    ]
    palette = {
        "repressed_sustained": "#d62728", "repressed_transient": "#ff9896",
        "not_responsive": "#aec7e8", "induced_sustained": "#2ca02c",
        "induced_transient": "#98df8a", "lrt_only": "#ff7f0e", "complex": "#9467bd",
    }

    ct = pd.crosstab(df["damage_tertile"], df["uv_trajectory"], normalize="index") * 100
    ct = ct.reindex(TERTILE_ORDER).fillna(0)
    cats_present = [c for c in traj_cats if c in ct.columns]

    fig, ax = plt.subplots(figsize=FIGSIZE_STACK)
    bottom = np.zeros(len(ct))
    for cat in cats_present:
        vals = ct[cat].values
        ax.bar(range(len(ct)), vals, bottom=bottom,
               label=cat, color=palette.get(cat, "grey"), edgecolor="white", linewidth=0.3)
        bottom += vals

    ax.set_xticks(range(len(ct)))
    ax.set_xticklabels(["Low", "Mid", "High"])
    ax.set_xlabel("Burden tertile")
    ax.set_ylabel("% genes")
    ax.set_title(f"{product}: UV trajectory composition by burden tertile")
    ax.legend(fontsize=7, loc="upper right", bbox_to_anchor=(1.35, 1))
    plt.tight_layout()
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()


def run_product(df: pd.DataFrame, product: str, out_stats: str, out_summary: str,
                out_lfc_fig: str, out_traj_fig: str, out_comp_fig: str, log_fn) -> None:
    log_fn(f"\n=== Product: {product} ===")

    g4 = df[df["group"] == "G4_TSS"].copy() if "group" in df.columns else df.copy()
    g4 = g4.dropna(subset=["damage_ds0"])
    log_fn(f"G4_TSS genes with damage_ds0: {len(g4):,}")

    if len(g4) < 9:
        log_fn("  WARN: too few genes for tertile analysis")
        pd.DataFrame().to_csv(out_stats, sep="\t", index=False)
        pd.DataFrame().to_csv(out_summary, sep="\t", index=False)
        return

    g4 = assign_tertiles(g4)
    for tert in TERTILE_ORDER:
        log_fn(f"  {tert}: n={int((g4['damage_tertile'] == tert).sum())}")

    # Summary stats
    all_summary = []
    for tp in TIMEPOINTS:
        s = tertile_summary(g4, tp)
        if not s.empty:
            all_summary.append(s)
    summary_df = pd.concat(all_summary, ignore_index=True) if all_summary else pd.DataFrame()

    # Statistical tests
    all_stats = []
    for tp in TIMEPOINTS:
        tests = run_tests(g4, tp, log_fn)
        cat_tests = categorical_tests(g4, tp, log_fn)
        for t in [tests, cat_tests]:
            if not t.empty:
                all_stats.append(t)

    stats_df = pd.concat(all_stats, ignore_index=True) if all_stats else pd.DataFrame()

    # BH correction across all planned tests per product
    if not stats_df.empty and "pvalue_raw" in stats_df.columns:
        valid = stats_df["pvalue_raw"].notna()
        if valid.sum() > 0:
            _, padj, _, _ = multipletests(stats_df.loc[valid, "pvalue_raw"], method="fdr_bh")
            stats_df.loc[valid, "pvalue_bh"] = padj
            log_fn(f"\nBH-corrected tests: {valid.sum()}")
            sig = stats_df[valid & (stats_df["pvalue_bh"] < 0.05)]
            log_fn(f"  Significant (BH < 0.05): {len(sig)}")

    stats_df.to_csv(out_stats, sep="\t", index=False)
    summary_df.to_csv(out_summary, sep="\t", index=False)
    log_fn(f"Written stats: {out_stats}")
    log_fn(f"Written summary: {out_summary}")

    # Figures
    plot_lfc_by_tertile(g4, product, out_lfc_fig)
    plot_median_trajectory(g4, product, out_traj_fig)
    plot_trajectory_composition(g4, product, out_comp_fig)
    log_fn(f"Figures written.")


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else sys.stdout

    def log(msg: str) -> None:
        print(msg, file=log_fh, flush=True)

    for path in [
        args.out_cpd_stats, args.out_64pp_stats,
        args.out_cpd_summary, args.out_64pp_summary,
        args.out_cpd_lfc_fig, args.out_64pp_lfc_fig,
    ]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    log("=== Task 3: Primary burden-stratified test ===")

    cpd = pd.read_csv(args.cpd_table, sep="\t")
    pp64 = pd.read_csv(args.pp64_table, sep="\t")
    log(f"CPD master table: {len(cpd):,} rows")
    log(f"64-PP master table: {len(pp64):,} rows")

    run_product(cpd, "CPD", args.out_cpd_stats, args.out_cpd_summary,
                args.out_cpd_lfc_fig, args.out_cpd_traj_fig, args.out_cpd_comp_fig, log)

    run_product(pp64, "64-PP", args.out_64pp_stats, args.out_64pp_summary,
                args.out_64pp_lfc_fig, args.out_64pp_traj_fig, args.out_64pp_comp_fig, log)

    log("\n=== Task 3 complete ===")
    if args.log:
        log_fh.close()


if __name__ == "__main__":
    main()
