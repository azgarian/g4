#!/usr/bin/env python3
"""Correlate promoter ATAC signal with RNA-seq expression by promoter group.

Steps:
  1. Extract per-gene mean ATAC signal from bigWig files at ±500 bp TSS windows.
  2. Extract matched RNA-seq normalised counts per timepoint.
  3. Build merged analysis table (group annotation + ATAC + RNA).
  4. Compute within-group Spearman correlations at each timepoint.
  5. Permutation test for correlation differences between groups.
  6. Plot results (scatter panel, rho heatmap, rho-over-time line plot).
  7. Write summary table with bootstrap confidence intervals.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
import seaborn as sns
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


TIMEPOINTS = ["t00", "t12", "t30", "t60"]
GROUPS = ["G4_TSS", "GC_bg_TSS", "No_overlap"]
GROUP_PALETTE = {
    "G4_TSS": "#1f77b4",
    "GC_bg_TSS": "#ff7f0e",
    "No_overlap": "#7f7f7f",
}
SMALL_GROUP_THRESHOLD = 50


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tss-annotation", required=True)
    p.add_argument("--promoter-windows-bed", required=True)
    p.add_argument("--norm-counts", required=True)
    p.add_argument("--bw-t00", required=True)
    p.add_argument("--bw-t12", required=True)
    p.add_argument("--bw-t30", required=True)
    p.add_argument("--bw-t60", required=True)
    p.add_argument("--out-atac-signal", required=True)
    p.add_argument("--out-rna-expression", required=True)
    p.add_argument("--out-merged", required=True)
    p.add_argument("--out-silent-genes", required=True)
    p.add_argument("--out-spearman", required=True)
    p.add_argument("--out-permutation", required=True)
    p.add_argument("--out-summary", required=True)
    p.add_argument("--out-scatter-pdf", required=True)
    p.add_argument("--out-heatmap-pdf", required=True)
    p.add_argument("--out-lineplot-pdf", required=True)
    p.add_argument("--n-permutations", type=int, default=10000)
    p.add_argument("--n-bootstrap", type=int, default=1000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--log", default=None)
    return p.parse_args()


def strip_version(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace(r"\.\d+$", "", regex=True)


# ─── Step 1: ATAC signal ──────────────────────────────────────────────────────

def extract_atac_signal(windows_bed: str, bw_paths: dict[str, str], log) -> pd.DataFrame:
    log("Step 1: Extracting promoter ATAC signal from bigWig files...")
    windows = pd.read_csv(
        windows_bed, sep="\t", header=None,
        names=["chrom", "start", "end", "gene_id", "score", "strand"],
    )
    log(f"  Loaded {len(windows)} promoter windows.")

    records: dict[str, list] = {"gene_id": list(windows["gene_id"])}
    for tp in TIMEPOINTS:
        records[f"atac_mean_{tp}"] = []

    for tp, bw_path in bw_paths.items():
        log(f"  Reading {bw_path}...")
        with pyBigWig.open(bw_path) as bw:
            vals = []
            for _, row in windows.iterrows():
                chrom = row["chrom"]
                start = int(row["start"])
                end = int(row["end"])
                try:
                    v = bw.stats(chrom, start, end, type="mean")[0]
                    vals.append(v if v is not None else 0.0)
                except RuntimeError:
                    vals.append(0.0)
        records[f"atac_mean_{tp}"] = vals

    df = pd.DataFrame(records)
    df["gene_id_nv"] = strip_version(df["gene_id"])
    for tp in TIMEPOINTS:
        col = f"atac_mean_{tp}"
        df[col] = df[col].fillna(0.0)
        df[f"atac_log2_{tp}"] = np.log2(df[col] + 1)

    log(f"  ATAC signal extracted for {len(df)} genes.")
    return df


# ─── Step 2: RNA-seq expression ───────────────────────────────────────────────

SAMPLE_GROUPS = {
    "t00": ["SU_100", "SU_200", "SU_300"],
    "t12": ["SU_112", "SU_212", "SU_312"],
    "t30": ["SU_130", "SU_230", "SU_330"],
    "t60": ["SU_160", "SU_260", "SU_360"],
}


def extract_rna_expression(norm_counts_path: str, log) -> pd.DataFrame:
    log("Step 2: Extracting RNA-seq expression per timepoint...")
    counts = pd.read_csv(norm_counts_path, sep="\t", compression="gzip", index_col=0)
    log(f"  Loaded counts matrix: {counts.shape[0]} genes x {counts.shape[1]} samples.")
    log(f"  Available samples: {list(counts.columns)}")

    df = pd.DataFrame({"gene_id": counts.index})
    df["gene_id_nv"] = strip_version(df["gene_id"])

    for tp, sample_cols in SAMPLE_GROUPS.items():
        present = [c for c in sample_cols if c in counts.columns]
        if not present:
            raise ValueError(
                f"None of the expected samples for {tp} found in counts matrix. "
                f"Expected: {sample_cols}. Available: {list(counts.columns)}"
            )
        missing = set(sample_cols) - set(present)
        if missing:
            log(f"  WARNING: samples missing for {tp}: {missing}")
        df[f"rna_mean_{tp}"] = counts[present].mean(axis=1).values
        df[f"rna_log2_{tp}"] = np.log2(df[f"rna_mean_{tp}"] + 1)

    log(f"  RNA expression computed for {len(df)} genes.")
    return df


# ─── Step 3: Merge ────────────────────────────────────────────────────────────

def build_merged(
    annotation_path: str, atac_df: pd.DataFrame, rna_df: pd.DataFrame, log
) -> tuple[pd.DataFrame, pd.DataFrame]:
    log("Step 3: Building merged analysis table...")

    anno = pd.read_csv(annotation_path, sep="\t")
    anno["gene_id_nv"] = strip_version(anno["gene_id"])
    n_anno = len(anno)

    atac_cols = ["gene_id_nv"] + [f"atac_log2_{tp}" for tp in TIMEPOINTS]
    merged = anno[["gene_id_nv", "gene_id", "gene_name", "chrom", "tss", "strand", "group"]].merge(
        atac_df[atac_cols], on="gene_id_nv", how="inner",
    )
    n_after_atac = len(merged)
    log(f"  After joining ATAC: {n_anno} -> {n_after_atac} genes (dropped {n_anno - n_after_atac})")

    rna_cols = ["gene_id_nv"] + [f"rna_log2_{tp}" for tp in TIMEPOINTS]
    merged = merged.merge(rna_df[rna_cols], on="gene_id_nv", how="inner")
    n_after_rna = len(merged)
    log(f"  After joining RNA: {n_after_atac} -> {n_after_rna} genes (dropped {n_after_atac - n_after_rna})")

    silent_mask = (merged["rna_log2_t00"] == 0) & (merged["atac_log2_t00"] == 0)
    silent = merged[silent_mask].copy()
    active = merged[~silent_mask].copy()
    log(f"  Silent genes (rna_log2_t00==0 AND atac_log2_t00==0): {len(silent)}")
    log(f"  Active genes retained: {len(active)}")

    keep_cols = (
        ["gene_id", "gene_name", "chrom", "tss", "strand", "group"]
        + [f"atac_log2_{tp}" for tp in TIMEPOINTS]
        + [f"rna_log2_{tp}" for tp in TIMEPOINTS]
    )
    return active[keep_cols].reset_index(drop=True), silent[keep_cols].reset_index(drop=True)


# ─── Step 4: Spearman correlations ────────────────────────────────────────────

def compute_spearman(merged: pd.DataFrame, log) -> pd.DataFrame:
    log("Step 4: Computing within-group Spearman correlations...")
    rows = []
    for tp in TIMEPOINTS:
        for grp in GROUPS:
            sub = merged[merged["group"] == grp]
            n = len(sub)
            x = sub[f"atac_log2_{tp}"].values
            y = sub[f"rna_log2_{tp}"].values
            if n >= 3:
                result = spearmanr(x, y)
                rho, pval = float(result.statistic), float(result.pvalue)
            else:
                rho, pval = np.nan, np.nan
            small_flag = (grp == "GC_bg_TSS") and (n < SMALL_GROUP_THRESHOLD)
            rows.append({
                "group": grp, "timepoint": tp, "n": n,
                "rho": rho, "p_value": pval, "small_n_flag": small_flag,
            })

    df = pd.DataFrame(rows)
    valid = df["p_value"].notna()
    df["p_adj"] = np.nan
    if valid.any():
        _, padj, _, _ = multipletests(df.loc[valid, "p_value"], method="fdr_bh")
        df.loc[valid, "p_adj"] = padj

    for _, r in df.iterrows():
        flag = " [small n — interpret cautiously]" if r["small_n_flag"] else ""
        log(f"  {r['group']} {r['timepoint']}: n={r['n']}, rho={r['rho']:.3f}, "
            f"p_adj={r['p_adj']:.3e}{flag}")
    return df


# ─── Step 5: Permutation tests ────────────────────────────────────────────────

def permutation_test(
    merged: pd.DataFrame, n_perm: int, seed: int, log
) -> pd.DataFrame:
    log(f"Step 5: Running permutation tests ({n_perm} permutations)...")
    rng = np.random.default_rng(seed)
    comparisons = [("G4_TSS", "No_overlap"), ("G4_TSS", "GC_bg_TSS")]
    rows = []

    for tp in TIMEPOINTS:
        for grp_a, grp_b in comparisons:
            sub_a = merged[merged["group"] == grp_a]
            sub_b = merged[merged["group"] == grp_b]
            na, nb = len(sub_a), len(sub_b)

            if na < 3 or nb < 3:
                rows.append({
                    "timepoint": tp, "group_a": grp_a, "group_b": grp_b,
                    "rho_a": np.nan, "rho_b": np.nan,
                    "observed_diff": np.nan, "empirical_p": np.nan,
                    "n_a": na, "n_b": nb,
                })
                continue

            combined = pd.concat([sub_a, sub_b], ignore_index=True)
            atac_all = combined[f"atac_log2_{tp}"].values
            rna_all = combined[f"rna_log2_{tp}"].values
            labels = combined["group"].values

            rho_a = float(spearmanr(sub_a[f"atac_log2_{tp}"].values,
                                    sub_a[f"rna_log2_{tp}"].values).statistic)
            rho_b = float(spearmanr(sub_b[f"atac_log2_{tp}"].values,
                                    sub_b[f"rna_log2_{tp}"].values).statistic)
            obs_diff = abs(rho_a - rho_b)

            perm_diffs = np.empty(n_perm)
            for i in range(n_perm):
                perm = rng.permutation(labels)
                ia = perm == grp_a
                ib = perm == grp_b
                if ia.sum() < 3 or ib.sum() < 3:
                    perm_diffs[i] = 0.0
                    continue
                ra = float(spearmanr(atac_all[ia], rna_all[ia]).statistic)
                rb = float(spearmanr(atac_all[ib], rna_all[ib]).statistic)
                perm_diffs[i] = abs(ra - rb)

            emp_p = float((perm_diffs >= obs_diff).mean())
            log(f"  {tp} {grp_a} vs {grp_b}: |rho_diff|={obs_diff:.3f}, empirical_p={emp_p:.4f}")
            rows.append({
                "timepoint": tp, "group_a": grp_a, "group_b": grp_b,
                "rho_a": rho_a, "rho_b": rho_b,
                "observed_diff": obs_diff, "empirical_p": emp_p,
                "n_a": na, "n_b": nb,
            })

    return pd.DataFrame(rows)


# ─── Step 6: Plots ────────────────────────────────────────────────────────────

def _lowess(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    from statsmodels.nonparametric.smoothers_lowess import lowess
    order = np.argsort(x)
    smoothed = lowess(y[order], x[order], frac=0.5, return_sorted=True)
    return smoothed[:, 0], smoothed[:, 1]


def plot_scatter_panel(
    merged: pd.DataFrame, spearman_df: pd.DataFrame, out_path: str, log
) -> None:
    log("Step 6a: Plotting scatter panel...")
    fig, axes = plt.subplots(len(GROUPS), len(TIMEPOINTS),
                             figsize=(4 * len(TIMEPOINTS), 4 * len(GROUPS)))

    for gi, grp in enumerate(GROUPS):
        sub_grp = merged[merged["group"] == grp]
        for ti, tp in enumerate(TIMEPOINTS):
            ax = axes[gi, ti]
            x = sub_grp[f"atac_log2_{tp}"].values
            y = sub_grp[f"rna_log2_{tp}"].values

            ax.scatter(x, y, s=4, alpha=0.3, color=GROUP_PALETTE[grp], rasterized=True)

            if len(x) >= 10:
                try:
                    xs, ys = _lowess(x, y)
                    ax.plot(xs, ys, color="black", lw=1.2)
                except Exception:
                    pass

            row = spearman_df[
                (spearman_df["group"] == grp) & (spearman_df["timepoint"] == tp)
            ]
            if not row.empty:
                rho = row["rho"].values[0]
                padj = row["p_adj"].values[0]
                if not np.isnan(rho):
                    sig = "**" if padj < 0.01 else ("*" if padj < 0.05 else "ns")
                    label = f"rho={rho:.2f}\np_adj={padj:.2e} {sig}"
                else:
                    label = "n<3"
                ax.text(0.05, 0.95, label, transform=ax.transAxes, fontsize=7,
                        va="top", ha="left",
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7))

            ax.set_xlabel("ATAC log2(mean+1)" if gi == len(GROUPS) - 1 else "")
            ax.set_ylabel("RNA log2(count+1)" if ti == 0 else "")
            title = f"{grp} | {tp}" if ti == 0 else tp
            ax.set_title(f"{grp}\n{tp}" if gi == 0 else (grp if ti == 0 else tp))

    plt.suptitle("Promoter ATAC vs RNA-seq expression by group and timepoint",
                 fontsize=11, y=1.01)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"  -> {out_path}")


def plot_rho_heatmap(spearman_df: pd.DataFrame, out_path: str, log) -> None:
    log("Step 6b: Plotting rho heatmap...")
    pivot_rho = spearman_df.pivot(index="group", columns="timepoint", values="rho")
    pivot_rho = pivot_rho.reindex(index=GROUPS, columns=TIMEPOINTS)
    pivot_padj = spearman_df.pivot(index="group", columns="timepoint", values="p_adj")
    pivot_padj = pivot_padj.reindex(index=GROUPS, columns=TIMEPOINTS)

    annot = np.empty(pivot_rho.shape, dtype=object)
    for i, grp in enumerate(GROUPS):
        for j, tp in enumerate(TIMEPOINTS):
            r = pivot_rho.loc[grp, tp]
            p = pivot_padj.loc[grp, tp]
            if pd.isna(r):
                annot[i, j] = "n/a"
            elif pd.isna(p):
                annot[i, j] = f"{r:.2f}"
            elif p < 0.01:
                annot[i, j] = f"{r:.2f}**"
            elif p < 0.05:
                annot[i, j] = f"{r:.2f}*"
            else:
                annot[i, j] = f"{r:.2f} ns"

    vmax = max(0.5, float(pivot_rho.abs().max().max()))
    fig, ax = plt.subplots(figsize=(7, 3.5))
    sns.heatmap(
        pivot_rho.astype(float), ax=ax,
        cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax,
        annot=annot, fmt="", linewidths=0.5,
        cbar_kws={"label": "Spearman rho"},
    )
    ax.set_title("Spearman rho: promoter ATAC vs RNA-seq expression")
    ax.set_xlabel("Timepoint")
    ax.set_ylabel("Promoter group")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"  -> {out_path}")


def _bootstrap_rho(
    x: np.ndarray, y: np.ndarray, n_boot: int, rng: np.random.Generator
) -> tuple[float, float]:
    n = len(x)
    boot = np.array([
        spearmanr(x[idx := rng.integers(0, n, size=n)],
                  y[idx]).statistic
        for _ in range(n_boot)
    ])
    return float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


def plot_rho_over_time(
    merged: pd.DataFrame, spearman_df: pd.DataFrame,
    n_boot: int, seed: int, out_path: str, log
) -> pd.DataFrame:
    log("Step 6c: Plotting rho over time with bootstrap CI...")
    rng = np.random.default_rng(seed)
    tp_nums = [0, 12, 30, 60]
    summary_rows = []

    fig, ax = plt.subplots(figsize=(7, 5))
    for grp in GROUPS:
        rhos, ci_los, ci_his = [], [], []
        for tp, tp_n in zip(TIMEPOINTS, tp_nums):
            sub = merged[merged["group"] == grp]
            n = len(sub)
            x = sub[f"atac_log2_{tp}"].values
            y = sub[f"rna_log2_{tp}"].values

            row = spearman_df[
                (spearman_df["group"] == grp) & (spearman_df["timepoint"] == tp)
            ]
            rho = float(row["rho"].values[0]) if not row.empty else np.nan
            pval = float(row["p_value"].values[0]) if not row.empty else np.nan
            padj = float(row["p_adj"].values[0]) if not row.empty else np.nan

            if n >= 3 and not np.isnan(rho):
                ci_lo, ci_hi = _bootstrap_rho(x, y, n_boot, rng)
            else:
                ci_lo, ci_hi = np.nan, np.nan

            rhos.append(rho)
            ci_los.append(ci_lo)
            ci_his.append(ci_hi)

            sig = "**" if (not np.isnan(padj) and padj < 0.01) else (
                "*" if (not np.isnan(padj) and padj < 0.05) else "ns"
            )
            summary_rows.append({
                "group": grp, "timepoint": tp, "n": n,
                "rho": rho, "p_value": pval, "p_adj": padj,
                "ci_lo": ci_lo, "ci_hi": ci_hi, "sig": sig,
            })

        color = GROUP_PALETTE[grp]
        xs = np.array(tp_nums, dtype=float)
        rhos_arr = np.array(rhos, dtype=float)
        ci_lo_arr = np.array(ci_los, dtype=float)
        ci_hi_arr = np.array(ci_his, dtype=float)
        valid = ~np.isnan(rhos_arr)

        ax.plot(xs[valid], rhos_arr[valid], marker="o", color=color, label=grp)
        if valid.sum() >= 2:
            ax.fill_between(xs[valid], ci_lo_arr[valid], ci_hi_arr[valid],
                            color=color, alpha=0.2)

    ax.axhline(0, color="gray", lw=0.8, ls="--")
    ax.set_xlabel("Time post-UV (min)")
    ax.set_ylabel("Spearman rho (ATAC vs RNA)")
    ax.set_xticks(tp_nums)
    ax.set_xticklabels(["0", "12", "30", "60"])
    ax.set_title("Promoter ATAC-RNA correlation over UV time-course")
    ax.legend(title="Group")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"  -> {out_path}")

    return pd.DataFrame(summary_rows)


# ─── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else None

    def log(msg: str) -> None:
        print(msg, flush=True)
        if log_fh:
            print(msg, file=log_fh, flush=True)

    bw_paths = {
        "t00": args.bw_t00,
        "t12": args.bw_t12,
        "t30": args.bw_t30,
        "t60": args.bw_t60,
    }

    for out in [args.out_atac_signal, args.out_rna_expression, args.out_merged,
                args.out_silent_genes, args.out_spearman, args.out_permutation,
                args.out_summary, args.out_scatter_pdf, args.out_heatmap_pdf,
                args.out_lineplot_pdf]:
        Path(out).parent.mkdir(parents=True, exist_ok=True)

    # Step 1
    atac_df = extract_atac_signal(args.promoter_windows_bed, bw_paths, log)
    atac_df.to_csv(args.out_atac_signal, sep="\t", index=False)
    log(f"  -> {args.out_atac_signal}")

    # Step 2
    rna_df = extract_rna_expression(args.norm_counts, log)
    rna_df.to_csv(args.out_rna_expression, sep="\t", index=False)
    log(f"  -> {args.out_rna_expression}")

    # Step 3
    merged, silent = build_merged(args.tss_annotation, atac_df, rna_df, log)
    merged.to_csv(args.out_merged, sep="\t", index=False)
    silent.to_csv(args.out_silent_genes, sep="\t", index=False)
    log(f"  -> {args.out_merged}")
    log(f"  -> {args.out_silent_genes}")

    # Step 4
    spearman_df = compute_spearman(merged, log)
    spearman_df.to_csv(args.out_spearman, sep="\t", index=False)
    log(f"  -> {args.out_spearman}")

    # Step 5
    perm_df = permutation_test(merged, n_perm=args.n_permutations,
                               seed=args.seed, log=log)
    perm_df.to_csv(args.out_permutation, sep="\t", index=False)
    log(f"  -> {args.out_permutation}")

    # Step 6
    plot_scatter_panel(merged, spearman_df, args.out_scatter_pdf, log)
    plot_rho_heatmap(spearman_df, args.out_heatmap_pdf, log)

    # Steps 6c + 7
    summary_df = plot_rho_over_time(
        merged, spearman_df,
        n_boot=args.n_bootstrap, seed=args.seed,
        out_path=args.out_lineplot_pdf, log=log,
    )
    summary_df.to_csv(args.out_summary, sep="\t", index=False)
    log(f"  -> {args.out_summary}")

    log("Done.")
    if log_fh:
        log_fh.close()


if __name__ == "__main__":
    main()
