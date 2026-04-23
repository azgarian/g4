#!/usr/bin/env python3
"""Task 4: Continuous association and covariate-adjusted modeling.

Tests whether lesion burden is continuously associated with UV repression, and
whether the association persists after accounting for baseline expression and
G4 occupancy strength.

Models:
  linear:   lfc_t ~ damage_ds0_z + log2_tpm_t0 + max_g4_signal_norm
  logistic: repressed_t ~ damage_ds0_z + log2_tpm_t0 + max_g4_signal_norm
  logistic: repressed_sustained ~ damage_ds0_z + log2_tpm_t0 + max_g4_signal_norm

Sign convention: positive lfc = repression after UV.
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
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
import statsmodels.api as sm


TIMEPOINTS = ["12", "30", "60"]
PREDICTORS = ["damage_ds0_z", "log2_tpm_t0", "max_g4_signal_norm"]
FIGSIZE = (6, 5)
COEF_FIGSIZE = (7, 5)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cpd-table", required=True)
    p.add_argument("--pp64-table", required=True)
    p.add_argument("--out-cpd-corr", required=True)
    p.add_argument("--out-64pp-corr", required=True)
    p.add_argument("--out-cpd-partial", required=True)
    p.add_argument("--out-64pp-partial", required=True)
    p.add_argument("--out-cpd-models", required=True)
    p.add_argument("--out-64pp-models", required=True)
    p.add_argument("--out-cpd-model-cmp", required=True)
    p.add_argument("--out-64pp-model-cmp", required=True)
    p.add_argument("--out-cpd-scatter", required=True)
    p.add_argument("--out-64pp-scatter", required=True)
    p.add_argument("--out-cpd-partial-fig", required=True)
    p.add_argument("--out-64pp-partial-fig", required=True)
    p.add_argument("--out-cpd-coef-fig", required=True)
    p.add_argument("--out-64pp-coef-fig", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[float, float]:
    """Partial Spearman correlation of x and y controlling for z via rank residuals."""
    valid = ~(np.isnan(x) | np.isnan(y) | np.isnan(z))
    x, y, z = x[valid], y[valid], z[valid]
    if len(x) < 5:
        return np.nan, np.nan
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    rz = stats.rankdata(z)
    # Residuals of rank regression
    def _resid(r_dep, r_ind):
        slope, intercept, _, _, _ = stats.linregress(r_ind, r_dep)
        return r_dep - (slope * r_ind + intercept)
    ex = _resid(rx, rz)
    ey = _resid(ry, rz)
    return stats.pearsonr(ex, ey)


def standardize(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    df = df.copy()
    for col in cols:
        if col in df.columns:
            mu, sd = df[col].mean(), df[col].std()
            df[col] = (df[col] - mu) / sd if sd > 0 else 0.0
    return df


def fit_linear(df: pd.DataFrame, outcome: str, predictors: list[str]) -> dict:
    sub = df[[outcome] + predictors].dropna()
    if len(sub) < 20:
        return {}
    formula = f"{outcome} ~ " + " + ".join(predictors)
    try:
        model = smf.ols(formula, data=sub).fit()
        return {
            "outcome": outcome, "n": len(sub),
            "r2_adj": model.rsquared_adj,
            "aic": model.aic,
            "predictor": predictors[0],
            "coef_damage": model.params.get("damage_ds0_z", np.nan),
            "se_damage": model.bse.get("damage_ds0_z", np.nan),
            "pvalue_damage": model.pvalues.get("damage_ds0_z", np.nan),
            "coef_tpm": model.params.get("log2_tpm_t0", np.nan),
            "coef_g4": model.params.get("max_g4_signal_norm", np.nan),
        }
    except Exception:
        return {}


def fit_logistic(df: pd.DataFrame, outcome: str, predictors: list[str]) -> dict:
    sub = df[[outcome] + predictors].dropna()
    if len(sub) < 20 or sub[outcome].sum() < 5:
        return {}
    formula = f"{outcome} ~ " + " + ".join(predictors)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = smf.logit(formula, data=sub).fit(disp=0)
        coef = model.params.get("damage_ds0_z", np.nan)
        se = model.bse.get("damage_ds0_z", np.nan)
        pval = model.pvalues.get("damage_ds0_z", np.nan)
        return {
            "outcome": outcome, "n": len(sub), "n_positive": int(sub[outcome].sum()),
            "aic": model.aic, "pseudo_r2": model.prsquared,
            "or_damage": np.exp(coef),
            "or_ci_lo": np.exp(coef - 1.96 * se),
            "or_ci_hi": np.exp(coef + 1.96 * se),
            "coef_damage": coef,
            "se_damage": se,
            "pvalue_damage": pval,
        }
    except Exception:
        return {}


def compare_models(df: pd.DataFrame, outcome: str, full_preds: list[str]) -> dict:
    """Compare nested models with and without damage_ds0_z."""
    sub = df[[outcome] + full_preds].dropna()
    if len(sub) < 20:
        return {}
    base_preds = [p for p in full_preds if p != "damage_ds0_z"]
    formula_full = f"{outcome} ~ " + " + ".join(full_preds)
    formula_base = f"{outcome} ~ " + " + ".join(base_preds) if base_preds else f"{outcome} ~ 1"
    try:
        m_full = smf.ols(formula_full, data=sub).fit()
        m_base = smf.ols(formula_base, data=sub).fit()
        delta_r2 = m_full.rsquared_adj - m_base.rsquared_adj
        delta_aic = m_base.aic - m_full.aic
        lr_stat = 2 * (m_full.llf - m_base.llf)
        lr_p = stats.chi2.sf(lr_stat, df=1)
        return {
            "outcome": outcome, "n": len(sub),
            "r2_adj_full": m_full.rsquared_adj,
            "r2_adj_base": m_base.rsquared_adj,
            "delta_r2_adj": delta_r2,
            "aic_full": m_full.aic,
            "aic_base": m_base.aic,
            "delta_aic": delta_aic,
            "lr_stat": lr_stat,
            "lr_pvalue": lr_p,
        }
    except Exception:
        return {}


def plot_scatter(df: pd.DataFrame, product: str, out_path: str) -> None:
    if "damage_ds0" not in df.columns or "lfc_60" not in df.columns:
        return
    sub = df[["damage_ds0", "lfc_60"]].dropna()
    if len(sub) < 5:
        return
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.scatter(sub["damage_ds0"], sub["lfc_60"], alpha=0.3, s=12, color="#1f77b4")
    # LOESS-like smooth using polynomial
    try:
        from numpy.polynomial import polynomial as Poly
        x = sub["damage_ds0"].values
        y = sub["lfc_60"].values
        sort_idx = np.argsort(x)
        xs = x[sort_idx]
        c = Poly.polyfit(xs, y[sort_idx], deg=2)
        ax.plot(xs, Poly.polyval(xs, c), color="red", linewidth=1.5, label="quadratic fit")
        ax.legend(fontsize=8)
    except Exception:
        pass
    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
    ax.set_xlabel(f"damage_ds0 ({product})")
    ax.set_ylabel("LFC at 60 min (repression > 0)")
    ax.set_title(f"{product}: Lesion burden vs LFC at 60 min")
    plt.tight_layout()
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()


def plot_partial_effect(df: pd.DataFrame, product: str, out_path: str) -> None:
    """Added-variable plot for damage_ds0_z on lfc_60 after removing confounders."""
    cols = ["lfc_60", "damage_ds0_z", "log2_tpm_t0", "max_g4_signal_norm"]
    sub = df[cols].dropna()
    if len(sub) < 20:
        plt.figure()
        plt.text(0.5, 0.5, "Insufficient data", ha="center", va="center")
        plt.savefig(out_path)
        plt.close()
        return

    confounders = ["log2_tpm_t0", "max_g4_signal_norm"]
    try:
        def _resid(dep_col, ind_cols):
            X = sm.add_constant(sub[ind_cols])
            return smf.ols(f"{dep_col} ~ " + " + ".join(ind_cols), data=sub).fit().resid

        e_y = _resid("lfc_60", confounders)
        e_x = _resid("damage_ds0_z", confounders)

        fig, ax = plt.subplots(figsize=FIGSIZE)
        ax.scatter(e_x, e_y, alpha=0.3, s=12, color="#1f77b4")
        m, b, r, p, se = stats.linregress(e_x, e_y)
        xs = np.linspace(e_x.min(), e_x.max(), 100)
        ax.plot(xs, m * xs + b, color="red", linewidth=1.5)
        ax.set_xlabel(f"Residual damage_ds0_z | confounders ({product})")
        ax.set_ylabel("Residual LFC_60 | confounders")
        ax.set_title(f"{product}: Partial effect of lesion burden on LFC_60\n"
                     f"Slope={m:.3f}, r={r:.3f}, p={p:.3e}")
        ax.axhline(0, color="grey", linewidth=0.4, linestyle="--")
        ax.axvline(0, color="grey", linewidth=0.4, linestyle="--")
        plt.tight_layout()
        plt.savefig(out_path, bbox_inches="tight")
        plt.close()
    except Exception as exc:
        plt.figure()
        plt.text(0.5, 0.5, f"Plot failed: {exc}", ha="center", va="center", fontsize=8)
        plt.savefig(out_path)
        plt.close()


def plot_coefficients(model_rows: list[dict], product: str, out_path: str) -> None:
    if not model_rows:
        plt.figure()
        plt.text(0.5, 0.5, "No model results", ha="center", va="center")
        plt.savefig(out_path)
        plt.close()
        return

    df_plot = pd.DataFrame(model_rows)
    df_plot = df_plot.dropna(subset=["coef_damage", "se_damage"])
    if df_plot.empty:
        plt.figure()
        plt.text(0.5, 0.5, "No coefficient data", ha="center", va="center")
        plt.savefig(out_path)
        plt.close()
        return

    fig, ax = plt.subplots(figsize=COEF_FIGSIZE)
    y_pos = range(len(df_plot))
    ax.errorbar(
        df_plot["coef_damage"], y_pos,
        xerr=1.96 * df_plot["se_damage"],
        fmt="o", color="#1f77b4", capsize=4, markersize=5,
    )
    ax.axvline(0, color="grey", linewidth=0.8, linestyle="--")
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(df_plot["outcome"].tolist(), fontsize=8)
    ax.set_xlabel("Coefficient (damage_ds0_z) ± 95% CI")
    ax.set_title(f"{product}: Lesion burden effect across models")
    plt.tight_layout()
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()


def run_product(
    df: pd.DataFrame, product: str,
    out_corr: str, out_partial: str, out_models: str, out_cmp: str,
    out_scatter: str, out_partial_fig: str, out_coef_fig: str,
    log_fn,
) -> None:
    log_fn(f"\n=== Product: {product} ===")

    g4 = df[df.get("group", df.get("is_g4_tss", pd.Series(True, index=df.index))).apply(
        lambda v: v == "G4_TSS" if isinstance(v, str) else bool(v)
    )].copy() if "group" in df.columns else df.copy()
    g4 = g4.dropna(subset=["damage_ds0"])
    log_fn(f"G4_TSS genes: {len(g4):,}")

    # ── Spearman correlations ─────────────────────────────────────────────────
    log_fn("\n--- Spearman correlations ---")
    corr_rows = []
    targets = [f"lfc_{tp}" for tp in TIMEPOINTS] + ["abs_lfc_60"]
    for target in targets:
        if target not in g4.columns:
            continue
        sub = g4[["damage_ds0", target]].dropna()
        if len(sub) < 5:
            continue
        r, p = spearmanr(sub["damage_ds0"], sub[target])
        log_fn(f"  damage_ds0 vs {target}: r={r:.3f}, p={p:.3e}, n={len(sub)}")
        corr_rows.append({
            "comparison": f"damage_ds0 vs {target}",
            "target": target, "n": len(sub),
            "spearman_r": r, "pvalue_raw": p,
        })

    # Lesion-vs-G4 collinearity check
    if "max_g4_signal_norm" in g4.columns:
        sub = g4[["damage_ds0", "max_g4_signal_norm"]].dropna()
        if len(sub) >= 5:
            r, p = spearmanr(sub["damage_ds0"], sub["max_g4_signal_norm"])
            log_fn(f"  damage_ds0 vs max_g4_signal_norm (collinearity): r={r:.3f}, p={p:.3e}")
            corr_rows.append({
                "comparison": "damage_ds0 vs max_g4_signal_norm (collinearity)",
                "target": "max_g4_signal_norm", "n": len(sub),
                "spearman_r": r, "pvalue_raw": p,
            })

    # Correlation with baseline expression
    if "log2_tpm_t0" in g4.columns:
        sub = g4[["damage_ds0", "log2_tpm_t0"]].dropna()
        if len(sub) >= 5:
            r, p = spearmanr(sub["damage_ds0"], sub["log2_tpm_t0"])
            log_fn(f"  damage_ds0 vs log2_tpm_t0: r={r:.3f}, p={p:.3e}")
            corr_rows.append({
                "comparison": "damage_ds0 vs log2_tpm_t0",
                "target": "log2_tpm_t0", "n": len(sub),
                "spearman_r": r, "pvalue_raw": p,
            })

    corr_df = pd.DataFrame(corr_rows)
    if not corr_df.empty and "pvalue_raw" in corr_df.columns:
        valid = corr_df["pvalue_raw"].notna()
        if valid.sum() > 0:
            _, padj, _, _ = multipletests(corr_df.loc[valid, "pvalue_raw"], method="fdr_bh")
            corr_df.loc[valid, "pvalue_bh"] = padj
    corr_df.to_csv(out_corr, sep="\t", index=False)
    log_fn(f"Written correlations: {out_corr}")

    # ── Partial Spearman ──────────────────────────────────────────────────────
    log_fn("\n--- Partial Spearman correlation ---")
    partial_rows = []
    if all(c in g4.columns for c in ["damage_ds0", "lfc_60", "log2_tpm_t0"]):
        sub = g4[["damage_ds0", "lfc_60", "log2_tpm_t0"]].dropna()
        r_part, p_part = partial_spearman(
            sub["damage_ds0"].values, sub["lfc_60"].values, sub["log2_tpm_t0"].values
        )
        log_fn(f"  Partial r (damage_ds0 ~ lfc_60 | log2_tpm_t0): r={r_part:.3f}, p={p_part:.3e}")
        partial_rows.append({
            "comparison": "damage_ds0 vs lfc_60 | log2_tpm_t0",
            "n": len(sub), "partial_spearman_r": r_part, "pvalue": p_part,
        })
    pd.DataFrame(partial_rows).to_csv(out_partial, sep="\t", index=False)

    # ── Linear and logistic models (standardized predictors) ──────────────────
    log_fn("\n--- Covariate-adjusted models ---")
    g4_std = standardize(g4, PREDICTORS)

    model_rows = []
    cmp_rows = []

    # Linear models for each timepoint
    for tp in TIMEPOINTS:
        outcome = f"lfc_{tp}"
        if outcome not in g4_std.columns:
            continue
        result = fit_linear(g4_std, outcome, PREDICTORS)
        if result:
            result["model_type"] = "linear"
            model_rows.append(result)
            log_fn(f"  Linear lfc_{tp}: coef_damage={result.get('coef_damage', np.nan):.3f}, "
                   f"p={result.get('pvalue_damage', np.nan):.3e}, R2_adj={result.get('r2_adj', np.nan):.3f}")
            cmp = compare_models(g4_std, outcome, PREDICTORS)
            if cmp:
                cmp["model_type"] = "linear"
                cmp_rows.append(cmp)

    # Logistic models
    for outcome in [f"repressed_{tp}" for tp in TIMEPOINTS] + ["repressed_sustained"]:
        if outcome not in g4_std.columns:
            continue
        result = fit_logistic(g4_std, outcome, PREDICTORS)
        if result:
            result["model_type"] = "logistic"
            model_rows.append(result)
            log_fn(f"  Logistic {outcome}: OR={result.get('or_damage', np.nan):.3f}, "
                   f"p={result.get('pvalue_damage', np.nan):.3e}")

    pd.DataFrame(model_rows).to_csv(out_models, sep="\t", index=False)
    pd.DataFrame(cmp_rows).to_csv(out_cmp, sep="\t", index=False)
    log_fn(f"Written models: {out_models}")
    log_fn(f"Written model comparison: {out_cmp}")

    # ── Figures ───────────────────────────────────────────────────────────────
    plot_scatter(g4, product, out_scatter)
    plot_partial_effect(g4_std, product, out_partial_fig)
    plot_coefficients(model_rows, product, out_coef_fig)
    log_fn("Figures written.")


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else sys.stdout

    def log(msg: str) -> None:
        print(msg, file=log_fh, flush=True)

    for path in [args.out_cpd_corr, args.out_64pp_corr]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    log("=== Task 4: Continuous association and covariate-adjusted modeling ===")

    cpd = pd.read_csv(args.cpd_table, sep="\t")
    pp64 = pd.read_csv(args.pp64_table, sep="\t")

    run_product(
        cpd, "CPD",
        args.out_cpd_corr, args.out_cpd_partial,
        args.out_cpd_models, args.out_cpd_model_cmp,
        args.out_cpd_scatter, args.out_cpd_partial_fig, args.out_cpd_coef_fig,
        log,
    )
    run_product(
        pp64, "64-PP",
        args.out_64pp_corr, args.out_64pp_partial,
        args.out_64pp_models, args.out_64pp_model_cmp,
        args.out_64pp_scatter, args.out_64pp_partial_fig, args.out_64pp_coef_fig,
        log,
    )

    log("\n=== Task 4 complete ===")
    if args.log:
        log_fh.close()


if __name__ == "__main__":
    main()
