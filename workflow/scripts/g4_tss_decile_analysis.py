#!/usr/bin/env python3
"""Task 5: Quantify G4 and GC-background occupancy across baseline expression deciles.

Assigns genes to expression deciles (0 = silent, 1-10 = expressed low to high),
computes overlap fractions with merged G4 peaks and GC-rich background per decile,
calculates Spearman correlations, and plots the trend.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gene-expression-by-group", required=True)
    p.add_argument("--tss-windows-bed", required=True,
                   help="canonical_tss_windows_1kb.bed")
    p.add_argument("--g4-merged-bed", required=True)
    p.add_argument("--gc-bg-bed", required=True)
    p.add_argument("--analysis-label", default=None)
    p.add_argument("--g4-label", default="G4 merged")
    p.add_argument("--out-deciles", required=True)
    p.add_argument("--out-overlap-fractions", required=True)
    p.add_argument("--out-correlation-stats", required=True)
    p.add_argument("--out-enrichment-plot", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def normalize_gene_id(gene_id: str) -> str:
    """Drop Ensembl version suffixes to match matrix gene IDs."""
    return str(gene_id).split(".", 1)[0]


def bedtools_intersect_fraction(windows_bed: str, peaks_bed: str) -> dict[str, bool]:
    """Return dict gene_id -> has_overlap (True/False)."""
    result = subprocess.run(
        ["bedtools", "intersect", "-nonamecheck", "-c", "-a", windows_bed, "-b", peaks_bed],
        capture_output=True, text=True, check=True
    )
    overlaps = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        gene_id = normalize_gene_id(fields[3])
        count = int(fields[-1])
        overlaps[gene_id] = count > 0
    return overlaps


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else None

    def log(msg: str) -> None:
        if log_fh:
            print(msg, file=log_fh, flush=True)

    if args.analysis_label:
        log(f"Analysis label: {args.analysis_label}")

    expr = pd.read_csv(args.gene_expression_by_group, sep="\t")
    expr["gene_id_norm"] = expr["gene_id"].map(normalize_gene_id)

    # assign deciles
    silent_mask = expr["mean_log2"] == 0
    expr["decile"] = 0
    expressed = expr[~silent_mask].copy()
    if len(expressed) > 0:
        decile_labels = list(range(1, 11))
        expressed["decile"] = pd.qcut(
            expressed["mean_log2"], q=10, labels=decile_labels, duplicates="drop"
        ).astype(int)
        expr.loc[~silent_mask, "decile"] = expressed["decile"].values

    log(f"Silent genes (decile 0): {silent_mask.sum():,}")
    log(f"Expressed genes (deciles 1-10): {(~silent_mask).sum():,}")

    for outpath in [args.out_deciles, args.out_overlap_fractions,
                    args.out_correlation_stats, args.out_enrichment_plot]:
        Path(outpath).parent.mkdir(parents=True, exist_ok=True)

    # compute overlaps
    log("Computing G4 overlaps...")
    g4_overlaps = bedtools_intersect_fraction(args.tss_windows_bed, args.g4_merged_bed)
    log("Computing GC-bg overlaps...")
    gc_overlaps = bedtools_intersect_fraction(args.tss_windows_bed, args.gc_bg_bed)

    expr["has_g4"] = expr["gene_id_norm"].map(g4_overlaps).fillna(False)
    expr["has_gc_bg"] = expr["gene_id_norm"].map(gc_overlaps).fillna(False)
    expr.drop(columns=["gene_id_norm"]).to_csv(args.out_deciles, sep="\t", index=False)

    frac_rows = []
    for decile in range(0, 11):
        sub = expr[expr["decile"] == decile]
        n = len(sub)
        if n == 0:
            continue
        frac_rows.append({
            "decile": decile,
            "n_genes": n,
            "g4_overlap_fraction": sub["has_g4"].mean(),
            "gc_bg_overlap_fraction": sub["has_gc_bg"].mean(),
        })
    frac_df = pd.DataFrame(frac_rows)
    frac_df.to_csv(args.out_overlap_fractions, sep="\t", index=False)

    # Spearman on deciles 1-10
    expressed_deciles = frac_df[frac_df["decile"] > 0]
    rho_g4, p_g4 = stats.spearmanr(expressed_deciles["decile"],
                                    expressed_deciles["g4_overlap_fraction"])
    rho_gc, p_gc = stats.spearmanr(expressed_deciles["decile"],
                                    expressed_deciles["gc_bg_overlap_fraction"])

    log(f"Spearman G4: rho={rho_g4:.4f}, p={p_g4:.4e}")
    log(f"Spearman GC-bg: rho={rho_gc:.4f}, p={p_gc:.4e}")

    corr_rows = [
        {"track": "G4_merged", "spearman_rho": rho_g4, "p_value": p_g4},
        {"track": "GC_bg", "spearman_rho": rho_gc, "p_value": p_gc},
    ]
    pd.DataFrame(corr_rows).to_csv(args.out_correlation_stats, sep="\t", index=False)

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(7, 4.5))
    x = frac_df["decile"].values
    ax.plot(x, frac_df["g4_overlap_fraction"].values * 100, "o-",
            color="#d62728", label=f"{args.g4_label} (ρ={rho_g4:.2f})")
    ax.plot(x, frac_df["gc_bg_overlap_fraction"].values * 100, "s--",
            color="#ff7f0e", label=f"GC-rich BG (ρ={rho_gc:.2f})")

    ax.axvline(0.5, color="gray", lw=0.8, ls=":")
    ax.set_xticks(range(0, 11))
    ax.set_xticklabels(["Silent"] + [str(i) for i in range(1, 11)])
    ax.set_xlabel("Expression decile (0 = silent)")
    ax.set_ylabel("Genes with TSS overlap (%)")
    title_prefix = f"{args.analysis_label}: " if args.analysis_label else ""
    ax.set_title(f"{title_prefix}G4 and GC-bg overlap fraction across expression deciles")
    ax.grid(axis="y", color="#d9d9d9", lw=0.6)
    ax.legend(fontsize=9)
    plt.tight_layout()
    fig.savefig(args.out_enrichment_plot, dpi=150)
    plt.close(fig)

    if log_fh:
        log_fh.close()


if __name__ == "__main__":
    main()
