#!/usr/bin/env python3
"""Task 10: Functional interpretation and G4-strength stratification.

Two analyses:
  1. Over-representation analysis (ORA) on pG4 LRT-significant genes vs pG4 background.
     Gene sets are loaded from a GMT file or per-file manifest (TSV with set_name and
     gene list columns). If no gene sets are provided the enrichment table is empty.

  2. G4-strength stratification: bin G4_TSS genes by promoter-level G4 signal
     aggregated from overlapping peaks, then test whether LRT significance and
     fold-change magnitude differ across bins using chi-squared and Kruskal-Wallis tests.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


# ── helpers ────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pg4-deg-intersect", required=True,
                   help="pG4_DEG_intersect.tsv from Task 9")
    p.add_argument("--tss-windows-bed", required=True,
                   help="canonical_tss_windows_1kb.bed used for promoter overlaps")
    p.add_argument("--g4-intersect-bed", required=True,
                   help="windows_g4_intersect.bed (bedtools -c output) for G4 peak counts")
    p.add_argument("--g4chip-source-bed", required=True,
                   help="Prepared BED6 with source_signal for G4 ChIP peaks")
    p.add_argument("--g4cuttag-source-bed", required=True,
                   help="Prepared BED6 with source_signal for G4 CUT&Tag peaks")
    p.add_argument("--gene-set-gmt", default=None,
                   help="GMT-format gene set file for ORA (optional)")
    p.add_argument("--gene-set-manifest", default=None,
                   help="TSV manifest with columns set_name, genes_file (optional)")
    p.add_argument("--gene-name-col", default="gene_name",
                   help="Column used to match gene sets (default: gene_name)")
    p.add_argument("--lrt-padj-threshold", type=float, default=0.05)
    p.add_argument("--ora-padj-threshold", type=float, default=0.1,
                   help="BH-adjusted p-value cutoff for ORA reporting")
    p.add_argument("--n-strength-bins", type=int, default=4,
                   help="Number of G4-strength quantile bins (default: 4)")
    p.add_argument("--strength-metric", default="max_signal",
                   choices=["max_signal", "sum_signal", "mean_signal", "peak_count"],
                   help="Promoter-level G4 strength metric used for stratification")
    p.add_argument("--out-enrichment-table", required=True)
    p.add_argument("--out-g4-strength-table", required=True)
    p.add_argument("--out-g4-strength-plot", required=True)
    p.add_argument("--out-g4-strength-direction-plot", required=True)
    p.add_argument("--out-summary", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def strip_version(s: pd.Series) -> pd.Series:
    return s.astype(str).str.replace(r"\.\d+$", "", regex=True)


# ── gene set loading ───────────────────────────────────────────────────────────

def load_gmt(gmt_path: str) -> dict[str, set[str]]:
    """Load gene sets from GMT format (name \\t description \\t gene1 \\t gene2 ...)."""
    gene_sets: dict[str, set[str]] = {}
    with open(gmt_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            genes = {g.strip() for g in parts[2:] if g.strip()}
            if genes:
                gene_sets[name] = genes
    return gene_sets


def load_manifest_gene_sets(manifest_path: str) -> dict[str, set[str]]:
    """Load gene sets from a TSV manifest pointing to per-set files."""
    manifest = pd.read_csv(manifest_path, sep="\t")
    manifest_dir = Path(manifest_path).parent
    gene_sets: dict[str, set[str]] = {}

    name_col = next((c for c in ("set_name", "name", "geneset") if c in manifest.columns), None)
    file_col = next((c for c in ("genes_file", "file", "filename", "path") if c in manifest.columns), None)
    gene_col = next((c for c in ("genes", "gene_list", "gene_ids") if c in manifest.columns), None)

    if name_col is None:
        return gene_sets

    for _, row in manifest.iterrows():
        sname = str(row[name_col])
        if file_col is not None and pd.notna(row.get(file_col)):
            fpath = Path(row[file_col])
            if not fpath.is_absolute():
                fpath = manifest_dir / fpath
            if fpath.exists():
                gene_sets[sname] = set(pd.read_csv(fpath, header=None)[0].astype(str))
        elif gene_col is not None and pd.notna(row.get(gene_col)):
            gene_sets[sname] = {g.strip() for g in str(row[gene_col]).split(",") if g.strip()}

    return gene_sets


# ── ORA ───────────────────────────────────────────────────────────────────────

def run_ora(
    query_genes: set[str],
    background_genes: set[str],
    gene_sets: dict[str, set[str]],
    min_set_size: int = 5,
    max_set_size: int = 500,
) -> pd.DataFrame:
    """Fisher's exact test ORA, one-sided (greater).

    For each gene set, compares the fraction of query genes in the set against
    the background fraction.  Returns a DataFrame sorted by raw p-value.
    """
    rows = []
    n_query = len(query_genes)
    n_bg = len(background_genes)

    for set_name, set_genes in gene_sets.items():
        set_genes_in_bg = set_genes & background_genes
        n_set = len(set_genes_in_bg)
        if not (min_set_size <= n_set <= max_set_size):
            continue
        n_query_in_set = len(query_genes & set_genes_in_bg)
        # contingency table:
        # [query & set ,  query & ~set]
        # [bg & set    ,  bg & ~set   ]  (bg excludes query to be conservative)
        bg_only = background_genes - query_genes
        n_bg_only = len(bg_only)
        n_bg_in_set = len(bg_only & set_genes_in_bg)

        table = [
            [n_query_in_set,       n_query - n_query_in_set],
            [n_bg_in_set,          n_bg_only - n_bg_in_set],
        ]
        if any(v < 0 for row in table for v in row):
            continue
        try:
            or_val, pval = stats.fisher_exact(table, alternative="greater")
        except Exception:
            or_val, pval = np.nan, np.nan

        rows.append({
            "set_name":          set_name,
            "n_set_in_bg":       n_set,
            "n_query_in_set":    n_query_in_set,
            "n_query":           n_query,
            "odds_ratio":        or_val,
            "p_value":           pval,
            "overlap_genes":     ",".join(sorted(query_genes & set_genes_in_bg)),
        })

    if not rows:
        return pd.DataFrame(columns=[
            "set_name", "n_set_in_bg", "n_query_in_set", "n_query",
            "odds_ratio", "p_value", "p_adj_BH", "overlap_genes",
        ])

    df = pd.DataFrame(rows).sort_values("p_value")
    _, p_adj, _, _ = multipletests(df["p_value"].fillna(1).values, method="fdr_bh")
    df["p_adj_BH"] = p_adj
    return df


# ── G4 strength ───────────────────────────────────────────────────────────────

def load_g4_counts(bed_path: str) -> pd.DataFrame:
    """Parse bedtools intersect -c output (BED6 + count) into gene_id → g4_peak_count."""
    df = pd.read_csv(bed_path, sep="\t", header=None)
    # col 3 = gene_id (name field), col 6 = overlap count
    result = pd.DataFrame({
        "gene_id": strip_version(df.iloc[:, 3].astype(str)),
        "g4_peak_count": df.iloc[:, 6].astype(int),
    })
    return result.drop_duplicates(subset=["gene_id"])


def _read_windows_bed(windows_bed: str) -> pd.DataFrame:
    df = pd.read_csv(windows_bed, sep="\t", header=None)
    return pd.DataFrame({
        "chrom": df.iloc[:, 0].astype(str),
        "start": df.iloc[:, 1].astype(int),
        "end": df.iloc[:, 2].astype(int),
        "gene_id": strip_version(df.iloc[:, 3].astype(str)),
    })


def _read_signal_bed6(bed_path: str) -> pd.DataFrame:
    df = pd.read_csv(bed_path, sep="\t", header=None)
    return pd.DataFrame({
        "chrom": df.iloc[:, 0].astype(str),
        "start": df.iloc[:, 1].astype(int),
        "end": df.iloc[:, 2].astype(int),
        "source_signal": pd.to_numeric(df.iloc[:, 4], errors="coerce"),
    }).dropna(subset=["source_signal"])


def _summarize_signal_hits(hit_rows: list[dict[str, Any]]) -> pd.DataFrame:
    if not hit_rows:
        return pd.DataFrame(columns=[
            "gene_id",
            "g4_signal_sum",
            "g4_signal_mean",
            "g4_signal_max",
            "g4_signal_peak_count",
        ])

    hits = pd.DataFrame(hit_rows)
    return (
        hits.groupby("gene_id", as_index=False)
        .agg(
            g4_signal_sum=("source_signal", "sum"),
            g4_signal_mean=("source_signal", "mean"),
            g4_signal_max=("source_signal", "max"),
            g4_signal_peak_count=("source_signal", "size"),
        )
    )


def _load_g4_signal_summary_python(windows_bed: str, source_beds: list[str]) -> pd.DataFrame:
    windows = _read_windows_bed(windows_bed)
    peaks = pd.concat([_read_signal_bed6(path) for path in source_beds], ignore_index=True)

    hit_rows: list[dict[str, Any]] = []
    for chrom, win_sub in windows.groupby("chrom", sort=False):
        peak_sub = peaks[peaks["chrom"] == chrom].sort_values(["start", "end"], kind="stable")
        if peak_sub.empty:
            continue

        peak_rows = list(peak_sub[["start", "end", "source_signal"]].itertuples(index=False, name=None))
        active: list[tuple[int, int, float]] = []
        next_peak_idx = 0

        for win in win_sub.sort_values(["start", "end"], kind="stable").itertuples(index=False):
            while next_peak_idx < len(peak_rows) and peak_rows[next_peak_idx][0] < win.end:
                active.append(peak_rows[next_peak_idx])
                next_peak_idx += 1

            active = [peak for peak in active if peak[1] > win.start]
            for peak_start, peak_end, peak_signal in active:
                if peak_start < win.end and peak_end > win.start:
                    hit_rows.append({
                        "gene_id": win.gene_id,
                        "source_signal": peak_signal,
                    })

    return _summarize_signal_hits(hit_rows)


def load_g4_signal_summary(windows_bed: str, source_beds: list[str]) -> pd.DataFrame:
    """Aggregate overlapping peak signals per promoter from one or more BED6 files."""
    hit_rows: list[dict[str, Any]] = []

    try:
        for bed_path in source_beds:
            result = subprocess.run(
                ["bedtools", "intersect", "-nonamecheck", "-wa", "-wb",
                 "-a", windows_bed, "-b", bed_path],
                capture_output=True,
                text=True,
                check=True,
            )
            for line in result.stdout.strip().split("\n"):
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) < 11:
                    continue
                hit_rows.append({
                    "gene_id": str(fields[3]).split(".", 1)[0],
                    "source_signal": float(fields[10]),
                })
    except FileNotFoundError:
        return _load_g4_signal_summary_python(windows_bed, source_beds)

    return _summarize_signal_hits(hit_rows)


def assign_quantile_bins(values: pd.Series, n_bins: int) -> pd.Series:
    """Assign Q1..Qn labels while gracefully handling tied distributions."""
    out = pd.Series(index=values.index, dtype="object")
    valid = values.dropna()
    if valid.empty:
        return out

    raw_bins = pd.qcut(valid, q=n_bins, duplicates="drop")
    n_actual = len(raw_bins.cat.categories)
    if n_actual == 0:
        return out

    labels = [f"Q{i+1}" for i in range(n_actual)]
    raw_bins = raw_bins.cat.rename_categories(labels)
    out.loc[valid.index] = raw_bins.astype(str).to_numpy()
    return out


def strength_metric_label(strength_metric: str) -> str:
    labels = {
        "peak_count": "promoter peak count",
        "max_signal": "max overlapping peak signal",
        "sum_signal": "sum of overlapping peak signals",
        "mean_signal": "mean overlapping peak signal",
    }
    return labels.get(strength_metric, strength_metric)


def attach_g4_strength_metrics(
    df: pd.DataFrame,
    g4_counts: pd.DataFrame,
    g4_signal: pd.DataFrame,
) -> pd.DataFrame:
    out = df.copy()
    out["gene_id_nv"] = strip_version(out["gene_id"])
    out = out.merge(
        g4_counts.rename(columns={"gene_id": "gene_id_nv"}),
        on="gene_id_nv", how="left",
    )
    out = out.merge(
        g4_signal.rename(columns={"gene_id": "gene_id_nv"}),
        on="gene_id_nv", how="left",
    )
    out["g4_peak_count"] = out["g4_peak_count"].fillna(0).astype(int)
    for col in ("g4_signal_sum", "g4_signal_mean", "g4_signal_max", "g4_signal_peak_count"):
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def stratify_g4_strength(
    g4_sub: pd.DataFrame,
    n_bins: int,
    lrt_padj_thresh: float,
    strength_metric: str,
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    """Bin G4_TSS genes by a chosen G4-strength metric and test associations.

    Returns the per-gene binned table and a list of test-result dicts.
    """
    g4_sub = g4_sub.copy()
    strength_col = {
        "peak_count": "g4_peak_count",
        "max_signal": "g4_signal_max",
        "sum_signal": "g4_signal_sum",
        "mean_signal": "g4_signal_mean",
    }[strength_metric]
    g4_sub["g4_strength_metric"] = strength_metric

    # Genes with count 0 shouldn't be here (all G4_TSS genes have >=1 overlap),
    # but guard against edge cases.
    has_peaks = g4_sub["g4_peak_count"] > 0
    if not has_peaks.any():
        g4_sub["g4_strength_bin"] = "no_peaks"
        return g4_sub, []

    # Quantile bins on genes that have at least one peak
    g4_sub.loc[has_peaks, "g4_strength_bin"] = assign_quantile_bins(
        g4_sub.loc[has_peaks, strength_col],
        n_bins,
    )
    g4_sub.loc[has_peaks & g4_sub["g4_strength_bin"].isna(), "g4_strength_bin"] = "Q1"
    g4_sub.loc[~has_peaks, "g4_strength_bin"] = "Q0_noPeaks"

    actual_bins = [b for b in g4_sub["g4_strength_bin"].unique() if b != "Q0_noPeaks"]

    # ── Chi-squared: LRT significance vs strength bin ─────────────────────────
    test_results: list[dict[str, Any]] = []
    sub_binned = g4_sub[g4_sub["g4_strength_bin"].isin(actual_bins)].copy()
    ct = pd.crosstab(sub_binned["g4_strength_bin"], sub_binned["lrt_sig"])
    if ct.shape[1] == 2 and ct.shape[0] > 1:
        chi2, chi_p, chi_df, _ = stats.chi2_contingency(ct.values)
        test_results.append({
            "test": f"chi2_lrt_sig_vs_g4_strength_bin_{strength_metric}",
            "statistic": chi2,
            "df": chi_df,
            "p_value": chi_p,
            "n_genes": len(sub_binned),
        })

    # ── Kruskal-Wallis: |log2FC at 60 min| vs strength bin ───────────────────
    fc_col = "lfc_0_vs_60"
    if fc_col in g4_sub.columns:
        fc_groups = [
            sub_binned.loc[sub_binned["g4_strength_bin"] == b, fc_col]
            .abs().dropna().values
            for b in actual_bins
        ]
        fc_groups = [g for g in fc_groups if len(g) > 0]
        if len(fc_groups) > 1:
            kw_stat, kw_p = stats.kruskal(*fc_groups)
            test_results.append({
                "test": f"kruskal_abs_lfc60_vs_g4_strength_bin_{strength_metric}",
                "statistic": kw_stat,
                "df": len(fc_groups) - 1,
                "p_value": kw_p,
                "n_genes": sum(len(g) for g in fc_groups),
            })

    # ── Spearman: selected strength metric vs lrt_padj (-log10) ──────────────
    valid = g4_sub[
        g4_sub["lrt_padj"].notna() & has_peaks & g4_sub[strength_col].notna()
    ].copy()
    if len(valid) > 10:
        neg_log_padj = -np.log10(valid["lrt_padj"].clip(lower=1e-300))
        rho, sp_p = stats.spearmanr(valid[strength_col], neg_log_padj)
        test_results.append({
            "test": f"spearman_{strength_metric}_vs_neglog10_lrt_padj",
            "statistic": rho,
            "df": len(valid) - 2,
            "p_value": sp_p,
            "n_genes": len(valid),
        })

    return g4_sub, test_results


def plot_g4_strength(
    g4_sub: pd.DataFrame,
    out_path: str,
    strength_metric: str,
) -> None:
    bins = [b for b in g4_sub["g4_strength_bin"].unique() if b != "Q0_noPeaks"]
    bins = sorted(bins)
    metric_label = strength_metric_label(strength_metric)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Panel A: fraction LRT-significant per bin
    frac_data = []
    for b in bins:
        sub = g4_sub[g4_sub["g4_strength_bin"] == b]
        tested = sub["lrt_padj"].notna()
        n_sig = sub.loc[tested, "lrt_sig"].sum()
        n_test = tested.sum()
        frac_data.append({
            "bin": b,
            "fraction": n_sig / n_test if n_test > 0 else np.nan,
            "n_label": n_test,
        })
    frac_df = pd.DataFrame(frac_data)

    ax = axes[0]
    ax.bar(frac_df["bin"], frac_df["fraction"], color="#d62728", edgecolor="black", linewidth=0.7)
    ax.set_xlabel(f"G4-strength bin (by {metric_label})")
    ax.set_ylabel("Fraction LRT-significant")
    ax.set_title("LRT significance vs G4 strength (G4_TSS genes)")
    for _, row in frac_df.iterrows():
        ax.text(row["bin"], (row["fraction"] or 0) + 0.005,
                f"n={int(row['n_label'])}", ha="center", fontsize=7)

    # Panel B: |log2FC at 60 min| per bin (violin)
    fc_col = "lfc_0_vs_60"
    ax2 = axes[1]
    if fc_col in g4_sub.columns:
        plot_data = [
            g4_sub.loc[g4_sub["g4_strength_bin"] == b, fc_col].abs().dropna().values
            for b in bins
        ]
        positions = list(range(1, len(bins) + 1))
        valid_pos = [(pos, dat) for pos, dat in zip(positions, plot_data) if len(dat) > 0]
        if valid_pos:
            ax2.violinplot(
                [d for _, d in valid_pos],
                positions=[p for p, _ in valid_pos],
                showmedians=True,
            )
        ax2.set_xticks(positions)
        ax2.set_xticklabels(bins)
        ax2.set_xlabel("G4-strength bin")
        ax2.set_ylabel("|log₂FC| (0 vs 60 min)")
        ax2.set_title("Fold-change magnitude vs G4 strength")
    else:
        ax2.text(0.5, 0.5, "lfc_0_vs_60 not available",
                 ha="center", va="center", transform=ax2.transAxes)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_g4_strength_direction(
    g4_sub: pd.DataFrame,
    out_path: str,
    strength_metric: str,
) -> None:
    """Violin plots of log2FC stratified by G4 strength, split by direction.

    Left panel: increased genes (lfc_0_vs_60 > 0), actual positive values.
    Right panel: decreased genes (lfc_0_vs_60 < 0), actual negative values.
    """
    bins = sorted(b for b in g4_sub["g4_strength_bin"].unique() if b != "Q0_noPeaks")
    fc_col = "lfc_0_vs_60"
    metric_label = strength_metric_label(strength_metric)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    if fc_col not in g4_sub.columns:
        for ax in axes:
            ax.text(0.5, 0.5, f"{fc_col} not available",
                    ha="center", va="center", transform=ax.transAxes)
        plt.tight_layout()
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        return

    panels = [
        (axes[0], g4_sub[g4_sub[fc_col] > 0], "Increased", "#d62728"),
        (axes[1], g4_sub[g4_sub[fc_col] < 0], "Decreased", "#1f77b4"),
    ]

    positions = list(range(1, len(bins) + 1))
    for ax, subset, direction, color in panels:
        plot_data = [
            subset.loc[subset["g4_strength_bin"] == b, fc_col].dropna().values
            for b in bins
        ]
        valid = [(p, d) for p, d in zip(positions, plot_data) if len(d) > 0]
        if valid:
            parts = ax.violinplot(
                [d for _, d in valid],
                positions=[p for p, _ in valid],
                showmedians=True,
            )
            for pc in parts.get("bodies", []):
                pc.set_facecolor(color)
                pc.set_alpha(0.7)
            for partname in ("cmedians", "cbars", "cmins", "cmaxes"):
                if partname in parts:
                    parts[partname].set_edgecolor("black")
                    parts[partname].set_linewidth(0.8)
        ax.axhline(0, color="black", linewidth=0.5, linestyle="--")
        ax.set_xticks(positions)
        ax.set_xticklabels(bins)
        ax.set_xlabel(f"G4-strength bin (by {metric_label})")
        ax.set_ylabel("log₂FC (0 vs 60 min)")
        ax.set_title(f"{direction} genes: log₂FC vs G4 strength")
        for pos, dat in zip(positions, plot_data):
            ymin = ax.get_ylim()[0]
            ax.text(pos, ymin, f"n={len(dat)}", ha="center", va="bottom", fontsize=7)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else None

    def log(msg: str) -> None:
        if log_fh:
            print(msg, file=log_fh, flush=True)

    for outpath in (
        args.out_enrichment_table,
        args.out_g4_strength_table,
        args.out_g4_strength_plot,
        args.out_g4_strength_direction_plot,
        args.out_summary,
    ):
        Path(outpath).parent.mkdir(parents=True, exist_ok=True)

    # ── load intersect table ──────────────────────────────────────────────────
    tbl = pd.read_csv(args.pg4_deg_intersect, sep="\t")
    gene_col = args.gene_name_col
    log(f"Loaded pG4_DEG_intersect: {len(tbl):,} genes")

    g4 = tbl[tbl["group"] == "G4_TSS"].copy()
    log(f"G4_TSS genes: {len(g4):,}")

    # ── ORA ───────────────────────────────────────────────────────────────────
    gene_sets: dict[str, set[str]] = {}
    if args.gene_set_gmt:
        gene_sets = load_gmt(args.gene_set_gmt)
        log(f"Loaded {len(gene_sets):,} gene sets from GMT: {args.gene_set_gmt}")
    elif args.gene_set_manifest:
        gene_sets = load_manifest_gene_sets(args.gene_set_manifest)
        log(f"Loaded {len(gene_sets):,} gene sets from manifest: {args.gene_set_manifest}")
    else:
        log("No gene sets provided — ORA skipped. Supply --gene-set-gmt or --gene-set-manifest.")

    if gene_sets:
        # background = all G4_TSS genes with a valid identifier
        background = set(g4[gene_col].dropna().astype(str))
        sig_genes = set(
            g4.loc[g4["lrt_padj"] < args.lrt_padj_threshold, gene_col]
            .dropna().astype(str)
        )
        log(f"ORA: {len(sig_genes):,} query genes, {len(background):,} background genes")
        ora_df = run_ora(sig_genes, background, gene_sets)
        log(f"ORA results: {len(ora_df):,} gene sets tested")
        sig_ora = ora_df[ora_df["p_adj_BH"] < args.ora_padj_threshold] if not ora_df.empty else ora_df
        log(f"  Significant (BH padj < {args.ora_padj_threshold}): {len(sig_ora):,}")
    else:
        ora_df = pd.DataFrame(columns=[
            "set_name", "n_set_in_bg", "n_query_in_set", "n_query",
            "odds_ratio", "p_value", "p_adj_BH", "overlap_genes",
        ])

    ora_df.to_csv(args.out_enrichment_table, sep="\t", index=False)
    log(f"Written: {args.out_enrichment_table}")

    # ── G4-strength stratification ─────────────────────────────────────────────
    g4_counts = load_g4_counts(args.g4_intersect_bed)
    log(f"G4 peak counts loaded: {len(g4_counts):,} entries")
    g4_signal = load_g4_signal_summary(
        args.tss_windows_bed,
        [args.g4chip_source_bed, args.g4cuttag_source_bed],
    )
    log(f"G4 promoter signal summaries loaded: {len(g4_signal):,} entries")

    g4_with_counts = attach_g4_strength_metrics(g4, g4_counts, g4_signal)

    log(f"G4_TSS genes matched with peak counts: "
        f"{g4_with_counts['g4_peak_count'].gt(0).sum():,} with ≥1 peak")
    log(f"  Peak count distribution:\n"
        f"  {g4_with_counts['g4_peak_count'].describe().to_string()}")
    for col in ("g4_signal_max", "g4_signal_sum", "g4_signal_mean"):
        if col in g4_with_counts.columns and g4_with_counts[col].notna().any():
            log(f"  {col} distribution:\n"
                f"  {g4_with_counts[col].dropna().describe().to_string()}")

    g4_binned, test_results = stratify_g4_strength(
        g4_with_counts, n_bins=args.n_strength_bins,
        lrt_padj_thresh=args.lrt_padj_threshold,
        strength_metric=args.strength_metric,
    )
    g4_binned = g4_binned.drop(columns=["gene_id_nv"], errors="ignore")
    g4_binned.to_csv(args.out_g4_strength_table, sep="\t", index=False)
    log(f"Written: {args.out_g4_strength_table}")

    for t in test_results:
        log(f"  {t['test']}: stat={t['statistic']:.4f}, p={t['p_value']:.4e}, n={t['n_genes']}")

    plot_g4_strength(
        g4_binned,
        args.out_g4_strength_plot,
        args.strength_metric,
    )
    log(f"Written: {args.out_g4_strength_plot}")

    plot_g4_strength_direction(
        g4_binned,
        args.out_g4_strength_direction_plot,
        args.strength_metric,
    )
    log(f"Written: {args.out_g4_strength_direction_plot}")

    # ── summary table ──────────────────────────────────────────────────────────
    summary_rows = []

    n_g4_total = len(g4)
    n_g4_lrt_sig = (g4["lrt_padj"] < args.lrt_padj_threshold).sum()
    summary_rows.append({
        "section": "ORA_background",
        "metric": "n_G4_TSS_genes_total",
        "value": str(n_g4_total),
    })
    summary_rows.append({
        "section": "ORA_background",
        "metric": "n_G4_TSS_LRT_significant",
        "value": str(int(n_g4_lrt_sig)),
    })
    summary_rows.append({
        "section": "G4_strength",
        "metric": "strength_metric",
        "value": args.strength_metric,
    })
    summary_rows.append({
        "section": "ORA_background",
        "metric": "n_gene_sets_tested",
        "value": str(len(ora_df)),
    })
    if not ora_df.empty and "p_adj_BH" in ora_df.columns:
        summary_rows.append({
            "section": "ORA_background",
            "metric": f"n_significant_padj_{args.ora_padj_threshold}",
            "value": str((ora_df["p_adj_BH"] < args.ora_padj_threshold).sum()),
        })
    for t in test_results:
        summary_rows.append({
            "section": "G4_strength",
            "metric": t["test"] + "_pvalue",
            "value": f"{t['p_value']:.4e}",
        })

    pd.DataFrame(summary_rows).to_csv(args.out_summary, sep="\t", index=False)
    log(f"Written: {args.out_summary}")

    if log_fh:
        log_fh.close()


if __name__ == "__main__":
    main()
