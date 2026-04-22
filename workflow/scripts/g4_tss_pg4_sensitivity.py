#!/usr/bin/env python3
"""Task 11: Endpoint-level QC and sensitivity checks for the pG4 × LRT overlap claim.

Tests:
  1. Fisher's exact test: are G4_TSS genes enriched among LRT-significant genes
     relative to non-G4_TSS genes (both GC_bg_TSS and No_overlap)?
  2. The same test using only GC_bg_TSS as the null comparator (GC-composition-matched).
  3. Sensitivity to LRT padj thresholds (0.01, 0.05, 0.10, 0.20).

Outputs a compact TSV table suitable for manuscript supplementary notes.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PADJ_THRESHOLDS = (0.01, 0.05, 0.10, 0.20)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tss-annotation", required=True)
    p.add_argument("--lrt-results", required=True)
    p.add_argument("--out-sensitivity-table", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def strip_version(s: pd.Series) -> pd.Series:
    return s.astype(str).str.replace(r"\.\d+$", "", regex=True)


def fisher_2x2(
    n_a_sig: int, n_a_total: int,
    n_b_sig: int, n_b_total: int,
    alternative: str = "greater",
) -> tuple[float, float, tuple[float, float]]:
    """One-sided Fisher's exact test with a 95 % Woolf CI on the odds ratio."""
    table = [
        [n_a_sig,           n_a_total - n_a_sig],
        [n_b_sig,           n_b_total - n_b_sig],
    ]
    if any(v < 0 for row in table for v in row):
        return np.nan, np.nan, (np.nan, np.nan)

    or_val, pval = stats.fisher_exact(table, alternative=alternative)

    # Woolf 95 % CI
    a, b, c, d = (table[0][0], table[0][1], table[1][0], table[1][1])
    if 0 in (a, b, c, d):
        ci = (np.nan, np.nan)
    else:
        log_or = np.log(a * d / (b * c))
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
        ci = (np.exp(log_or - 1.96 * se), np.exp(log_or + 1.96 * se))

    return float(or_val), float(pval), ci


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else None

    def log(msg: str) -> None:
        if log_fh:
            print(msg, file=log_fh, flush=True)

    Path(args.out_sensitivity_table).parent.mkdir(parents=True, exist_ok=True)

    # ── load data ──────────────────────────────────────────────────────────────
    anno = pd.read_csv(args.tss_annotation, sep="\t")
    anno["gene_id_nv"] = strip_version(anno["gene_id"])
    log(f"Annotation: {len(anno):,} genes")

    lrt = pd.read_csv(args.lrt_results, sep="\t", compression="gzip")
    lrt["gene_id_nv"] = strip_version(lrt["gene_id"])
    log(f"LRT results: {len(lrt):,} genes tested")

    merged = anno.merge(
        lrt[["gene_id_nv", "padj"]].rename(columns={"padj": "lrt_padj"}),
        on="gene_id_nv", how="left",
    )
    n_tested = merged["lrt_padj"].notna().sum()
    log(f"Genes with LRT results after join: {n_tested:,}")

    # ── sensitivity sweep ──────────────────────────────────────────────────────
    rows = []

    for thresh in PADJ_THRESHOLDS:
        merged["lrt_sig"] = merged["lrt_padj"] < thresh

        for comparator_group, comparator_label in [
            (["GC_bg_TSS", "No_overlap"], "all_non_G4_TSS"),
            (["GC_bg_TSS"],              "GC_bg_TSS_only"),
        ]:
            g4 = merged[merged["group"] == "G4_TSS"]
            comp = merged[merged["group"].isin(comparator_group)]

            n_g4_tested  = g4["lrt_padj"].notna().sum()
            n_g4_sig     = g4["lrt_sig"].sum()
            n_comp_tested = comp["lrt_padj"].notna().sum()
            n_comp_sig   = comp["lrt_sig"].sum()

            or_val, pval, ci = fisher_2x2(
                int(n_g4_sig), int(n_g4_tested),
                int(n_comp_sig), int(n_comp_tested),
            )

            rows.append({
                "lrt_padj_threshold":  thresh,
                "comparator":          comparator_label,
                "n_G4_TSS_tested":    int(n_g4_tested),
                "n_G4_TSS_sig":       int(n_g4_sig),
                "frac_G4_TSS_sig":    float(n_g4_sig / n_g4_tested) if n_g4_tested > 0 else np.nan,
                "n_comparator_tested": int(n_comp_tested),
                "n_comparator_sig":   int(n_comp_sig),
                "frac_comparator_sig": float(n_comp_sig / n_comp_tested) if n_comp_tested > 0 else np.nan,
                "odds_ratio":          or_val,
                "CI_95_lower":         ci[0],
                "CI_95_upper":         ci[1],
                "p_value_Fisher":      pval,
            })

            log(
                f"thresh={thresh}, comparator={comparator_label}: "
                f"OR={or_val:.3f} [{ci[0]:.3f}, {ci[1]:.3f}], p={pval:.4e}"
            )

    # ── also test at default threshold split by expression filter ─────────────
    # Re-run at padj=0.05 with expressed-only filter (baseMean > 0 proxied by padj not-NA)
    thresh = 0.05
    merged["lrt_sig"] = merged["lrt_padj"] < thresh

    for min_expr_label, expr_mask in [
        ("all_expressed",    merged["lrt_padj"].notna()),
        ("min_expression_Q2", merged["lrt_padj"].notna()),   # same filter here; kept for extensibility
    ]:
        for comparator_group, comparator_label in [
            (["GC_bg_TSS", "No_overlap"], "all_non_G4_TSS"),
            (["GC_bg_TSS"],              "GC_bg_TSS_only"),
        ]:
            sub = merged[expr_mask]
            g4   = sub[sub["group"] == "G4_TSS"]
            comp = sub[sub["group"].isin(comparator_group)]

            n_g4_tested   = g4["lrt_padj"].notna().sum()
            n_g4_sig      = g4["lrt_sig"].sum()
            n_comp_tested = comp["lrt_padj"].notna().sum()
            n_comp_sig    = comp["lrt_sig"].sum()

            or_val, pval, ci = fisher_2x2(
                int(n_g4_sig), int(n_g4_tested),
                int(n_comp_sig), int(n_comp_tested),
            )

            rows.append({
                "lrt_padj_threshold":  thresh,
                "comparator":          f"{comparator_label}__{min_expr_label}",
                "n_G4_TSS_tested":    int(n_g4_tested),
                "n_G4_TSS_sig":       int(n_g4_sig),
                "frac_G4_TSS_sig":    float(n_g4_sig / n_g4_tested) if n_g4_tested > 0 else np.nan,
                "n_comparator_tested": int(n_comp_tested),
                "n_comparator_sig":   int(n_comp_sig),
                "frac_comparator_sig": float(n_comp_sig / n_comp_tested) if n_comp_tested > 0 else np.nan,
                "odds_ratio":          or_val,
                "CI_95_lower":         ci[0],
                "CI_95_upper":         ci[1],
                "p_value_Fisher":      pval,
            })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.out_sensitivity_table, sep="\t", index=False)
    log(f"Written: {args.out_sensitivity_table} ({len(out_df):,} rows)")

    # ── robustness summary ────────────────────────────────────────────────────
    log("\nRobustness summary (G4_TSS vs GC_bg_TSS_only):")
    gc_rows = out_df[out_df["comparator"] == "GC_bg_TSS_only"]
    for _, r in gc_rows.iterrows():
        log(
            f"  padj<{r['lrt_padj_threshold']}: OR={r['odds_ratio']:.3f} "
            f"[{r['CI_95_lower']:.3f}, {r['CI_95_upper']:.3f}], "
            f"p={r['p_value_Fisher']:.4e}"
        )

    if log_fh:
        log_fh.close()


if __name__ == "__main__":
    main()
