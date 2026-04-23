#!/usr/bin/env python3
"""Task 5: G4 strength as a continuous predictor of UV-response magnitude and direction.

Assigns max G4 peak signal per G4_TSS gene (from ChIP and CUT&Tag separately),
normalizes within assay, then tests Spearman correlations and tertile comparisons.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests


GROUP_PALETTE = {
    "G4_TSS": "#1f77b4",
    "GC_bg_TSS": "#ff7f0e",
    "No_overlap": "#7f7f7f",
}
DIR_PALETTE = {"repressed": "#d62728", "induced": "#2ca02c", "not_sig": "#aec7e8"}
TIMEPOINTS = ["12", "30", "60"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--master-table", required=True)
    p.add_argument("--tss-windows", required=True,
                   help="canonical_tss_windows_1kb.bed")
    p.add_argument("--g4chip-peaks", required=True,
                   help="g4_hela_peaks_prepared.tsv (ChIP)")
    p.add_argument("--g4cuttag-peaks", required=True,
                   help="g4_hela_peaks_prepared.tsv (CUT&Tag)")
    p.add_argument("--out-correlations", required=True)
    p.add_argument("--out-direction-by-tertile", required=True)
    p.add_argument("--out-scatter", required=True)
    p.add_argument("--out-violin", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def strip_version(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace(r"\.\d+$", "", regex=True)


def rank_normalize(series: pd.Series) -> pd.Series:
    """Rank-normalize to [0, 1] interval."""
    ranked = series.rank(method="average", na_option="keep")
    return (ranked - 1) / (ranked.notna().sum() - 1)


def bedtools_intersect_max_signal(
    windows_bed: str, peaks_tsv: str, assay_name: str, log_fn
) -> pd.Series:
    """
    Run bedtools intersect -a windows -b peaks and return max source_signal per gene_id.
    Windows BED must have gene_id in column 4 (name field).
    Returns a Series indexed by gene_id (no version suffix).
    """
    peaks = pd.read_csv(peaks_tsv, sep="\t")
    # Write peaks as BED with source_signal in a usable column
    with tempfile.NamedTemporaryFile(suffix=".bed", mode="w", delete=False) as tf:
        peaks_bed_path = tf.name
        for _, row in peaks.iterrows():
            tf.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['name']}\t{row['source_signal']}\n")

    with tempfile.NamedTemporaryFile(suffix=".tsv", mode="w", delete=False) as tf:
        out_path = tf.name

    cmd = [
        "bedtools", "intersect",
        "-a", windows_bed,
        "-b", peaks_bed_path,
        "-wo",
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        lines = result.stdout.strip().split("\n") if result.stdout.strip() else []
        log_fn(f"  {assay_name}: {len(lines)} overlap records from bedtools intersect")
    except subprocess.CalledProcessError as e:
        log_fn(f"  WARNING: bedtools failed for {assay_name}: {e.stderr}")
        return pd.Series(dtype=float, name=f"max_g4_signal_{assay_name}")

    Path(peaks_bed_path).unlink(missing_ok=True)

    if not lines:
        return pd.Series(dtype=float, name=f"max_g4_signal_{assay_name}")

    # Parse overlap output
    # Columns: window_chrom, window_start, window_end, gene_id, [window extra cols...], peak_chrom, peak_start, peak_end, peak_name, source_signal, overlap_bp
    records = []
    for line in lines:
        parts = line.split("\t")
        # windows BED has variable columns; gene_id is always col 3 (0-indexed)
        # peaks BED has 5 columns, so:
        # a_cols + b_cols + 1 (overlap)
        # a has however many cols the window BED has; b has 5
        gene_id = parts[3]
        # source_signal is 5th col of b (index = n_a_cols + 4)
        # Count a columns: everything before last 6 (5 b cols + 1 overlap)
        n_b_plus_overlap = 6
        n_a = len(parts) - n_b_plus_overlap
        signal_idx = n_a + 4
        try:
            signal = float(parts[signal_idx])
        except (IndexError, ValueError):
            continue
        records.append({"gene_id": gene_id, "source_signal": signal})

    if not records:
        return pd.Series(dtype=float, name=f"max_g4_signal_{assay_name}")

    overlap_df = pd.DataFrame(records)
    overlap_df["gene_id_nv"] = strip_version(overlap_df["gene_id"])
    max_signal = overlap_df.groupby("gene_id_nv")["source_signal"].max()
    max_signal.name = f"max_g4_signal_{assay_name}"
    log_fn(f"  {assay_name}: {len(max_signal):,} genes with overlapping G4 peak")
    return max_signal


def bootstrap_spearman_ci(
    x: np.ndarray, y: np.ndarray, n_boot: int = 1000, ci: float = 0.95
) -> tuple[float, float, float, float]:
    """Returns rho, pval, ci_lo, ci_hi."""
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    if len(x) < 5:
        return np.nan, np.nan, np.nan, np.nan
    rho, pval = stats.spearmanr(x, y)
    boot_rhos = []
    rng = np.random.default_rng(42)
    for _ in range(n_boot):
        idx = rng.integers(0, len(x), size=len(x))
        r, _ = stats.spearmanr(x[idx], y[idx])
        boot_rhos.append(r)
    lo = np.percentile(boot_rhos, (1 - ci) / 2 * 100)
    hi = np.percentile(boot_rhos, (1 + ci) / 2 * 100)
    return rho, pval, lo, hi


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else sys.stdout

    def log(msg: str) -> None:
        print(msg, file=log_fh, flush=True)

    for p in [args.out_correlations, args.out_direction_by_tertile,
              args.out_scatter, args.out_violin]:
        Path(p).parent.mkdir(parents=True, exist_ok=True)

    log("=== Task 5: G4 strength as continuous predictor ===")

    df = pd.read_csv(args.master_table, sep="\t")
    g4_df = df[df["group"] == "G4_TSS"].copy()
    g4_df["gene_id_nv"] = strip_version(g4_df["gene_id"])
    log(f"G4_TSS genes: {len(g4_df):,}")

    # --- Assign max G4 signal per gene from each assay ---
    log("Running bedtools intersect for ChIP...")
    chip_signal = bedtools_intersect_max_signal(
        args.tss_windows, args.g4chip_peaks, "chip", log
    )
    log("Running bedtools intersect for CUT&Tag...")
    cuttag_signal = bedtools_intersect_max_signal(
        args.tss_windows, args.g4cuttag_peaks, "cuttag", log
    )

    g4_df = g4_df.join(chip_signal, on="gene_id_nv")
    g4_df = g4_df.join(cuttag_signal, on="gene_id_nv")

    log(f"G4_TSS genes with ChIP signal: {g4_df['max_g4_signal_chip'].notna().sum():,}")
    log(f"G4_TSS genes with CUT&Tag signal: {g4_df['max_g4_signal_cuttag'].notna().sum():,}")

    # --- Normalize within assay ---
    # Log-transform then rank-normalize to [0,1]
    for assay in ["chip", "cuttag"]:
        col = f"max_g4_signal_{assay}"
        log_col = f"log_g4_signal_{assay}"
        norm_col = f"norm_g4_signal_{assay}"
        g4_df[log_col] = np.log1p(g4_df[col])
        g4_df[norm_col] = rank_normalize(g4_df[log_col])

    # Combined predictor: max of the two normalized signals
    g4_df["max_g4_signal_norm"] = g4_df[["norm_g4_signal_chip", "norm_g4_signal_cuttag"]].max(axis=1)
    log(f"Genes with combined norm signal: {g4_df['max_g4_signal_norm'].notna().sum():,}")

    # --- Spearman correlations with bootstrap CIs ---
    targets = {
        "abs_lfc_60": "abs_lfc_60",
        "lfc_60": "lfc_60 (signed)",
        "abs_lfc_max": "abs_lfc_max",
        "lrt_stat": "lrt_stat (omnibus)",
    }
    predictors = {
        "max_g4_signal_norm": "combined_norm",
        "max_g4_signal_chip": "chip_raw",
        "max_g4_signal_cuttag": "cuttag_raw",
    }

    corr_rows = []
    for pred_col, pred_label in predictors.items():
        x = g4_df[pred_col].values
        for tgt_col, tgt_label in targets.items():
            y = g4_df[tgt_col].values
            rho, pval, ci_lo, ci_hi = bootstrap_spearman_ci(x, y)
            corr_rows.append({
                "predictor": pred_label,
                "target": tgt_col,
                "spearman_rho": rho,
                "pval": pval,
                "ci_lo_95": ci_lo,
                "ci_hi_95": ci_hi,
                "n": int((~(np.isnan(x) | np.isnan(y))).sum()),
            })
            log(f"  {pred_label} ~ {tgt_col}: rho={rho:.3f}, p={pval:.3e}, CI=[{ci_lo:.3f},{ci_hi:.3f}]")

    corr_df = pd.DataFrame(corr_rows)
    corr_df.to_csv(args.out_correlations, sep="\t", index=False)
    log(f"Written: {args.out_correlations}")

    # --- Tertile direction statistics ---
    valid = g4_df[g4_df["max_g4_signal_norm"].notna()].copy()
    valid["g4_tertile"] = pd.qcut(
        valid["max_g4_signal_norm"], q=3, labels=["G4_low", "G4_mid", "G4_high"]
    )
    log(f"\nTertile sizes: {valid['g4_tertile'].value_counts().to_dict()}")

    dir_rows = []
    for tp in TIMEPOINTS:
        dir_col = f"dir_{tp}"
        lfc_col = f"lfc_{tp}"
        for tertile in ["G4_low", "G4_mid", "G4_high"]:
            sub = valid[valid["g4_tertile"] == tertile]
            n = len(sub)
            n_rep = (sub[dir_col] == "repressed").sum()
            n_ind = (sub[dir_col] == "induced").sum()
            dir_rows.append({
                "g4_tertile": tertile,
                "timepoint": tp,
                "n": n,
                "n_repressed": int(n_rep),
                "n_induced": int(n_ind),
                "pct_repressed": n_rep / n if n > 0 else np.nan,
                "pct_induced": n_ind / n if n > 0 else np.nan,
                "signed_mean_lfc": sub[lfc_col].mean(),
            })

    dir_tertile_df = pd.DataFrame(dir_rows)
    dir_tertile_df.to_csv(args.out_direction_by_tertile, sep="\t", index=False)
    log(f"Written: {args.out_direction_by_tertile}")

    # --- Scatter: max_g4_signal_norm vs lfc_60 ---
    scatter_df = g4_df[["max_g4_signal_norm", "lfc_60", "dir_60"]].dropna()
    fig, ax = plt.subplots(figsize=(7, 5))
    for direction, color in DIR_PALETTE.items():
        sub = scatter_df[scatter_df["dir_60"] == direction]
        ax.scatter(sub["max_g4_signal_norm"], sub["lfc_60"],
                   s=8, alpha=0.5, color=color, label=direction, rasterized=True)
    # LOWESS smoothing
    if len(scatter_df) > 20:
        from statsmodels.nonparametric.smoothers_lowess import lowess
        smoothed = lowess(
            scatter_df["lfc_60"].values,
            scatter_df["max_g4_signal_norm"].values,
            frac=0.3, return_sorted=True,
        )
        ax.plot(smoothed[:, 0], smoothed[:, 1], color="black", lw=2, label="LOWESS")
    ax.axhline(0, color="gray", lw=0.8, ls="--")
    ax.set_xlabel("G4 signal (combined, normalized)")
    ax.set_ylabel("log₂FC at 60 min")
    ax.set_title("G4 strength vs UV response (60 min)")
    ax.legend(fontsize=8, markerscale=2)
    plt.tight_layout()
    fig.savefig(args.out_scatter, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"Written: {args.out_scatter}")

    # --- Violin/boxplot: abs_lfc_60 by tertile ---
    violin_df = valid[["g4_tertile", "abs_lfc_60"]].dropna()
    tertile_order = ["G4_low", "G4_mid", "G4_high"]
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.violinplot(data=violin_df, x="g4_tertile", y="abs_lfc_60",
                   order=tertile_order, color=GROUP_PALETTE["G4_TSS"],
                   inner=None, linewidth=0.8, alpha=0.7, ax=ax)
    sns.boxplot(data=violin_df, x="g4_tertile", y="abs_lfc_60",
                order=tertile_order, width=0.15, fliersize=1.5,
                linewidth=0.8, color="white", ax=ax)

    # Pairwise Wilcoxon BH-corrected annotations
    pairs = [("G4_high", "G4_low"), ("G4_high", "G4_mid"), ("G4_mid", "G4_low")]
    pvals_w = []
    for g1, g2 in pairs:
        x1 = violin_df.loc[violin_df["g4_tertile"] == g1, "abs_lfc_60"].values
        x2 = violin_df.loc[violin_df["g4_tertile"] == g2, "abs_lfc_60"].values
        if len(x1) > 0 and len(x2) > 0:
            _, p = stats.mannwhitneyu(x1, x2, alternative="two-sided")
        else:
            p = 1.0
        pvals_w.append(p)
    _, padj_w, _, _ = multipletests(pvals_w, method="fdr_bh")
    y_max = violin_df["abs_lfc_60"].max()
    for i, ((g1, g2), padj_val) in enumerate(zip(pairs, padj_w)):
        y_pos = y_max * (1.05 + i * 0.07)
        x1_pos = tertile_order.index(g1)
        x2_pos = tertile_order.index(g2)
        ax.plot([x1_pos, x2_pos], [y_pos, y_pos], color="black", lw=0.8)
        ax.text((x1_pos + x2_pos) / 2, y_pos + y_max * 0.01,
                f"p={padj_val:.3e}", ha="center", fontsize=7)

    ax.set_xlabel("G4 signal tertile")
    ax.set_ylabel("|log₂FC| at 60 min")
    ax.set_title("UV magnitude by G4 signal strength (G4_TSS genes)")
    plt.tight_layout()
    fig.savefig(args.out_violin, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"Written: {args.out_violin}")

    if args.log:
        log_fh.close()


if __name__ == "__main__":
    main()
