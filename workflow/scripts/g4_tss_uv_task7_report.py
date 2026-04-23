#!/usr/bin/env python3
"""Task 7: Synthesis HTML report — magnitude and direction of UV response in G4-promoter genes."""

from __future__ import annotations

import argparse
import base64
import importlib.metadata
import sys
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    # Task 2
    p.add_argument("--magnitude-stats", required=True)
    p.add_argument("--magnitude-summary", required=True)
    p.add_argument("--magnitude-figure", required=True)
    # Task 3
    p.add_argument("--direction-stats", required=True)
    p.add_argument("--direction-fisher", required=True)
    p.add_argument("--stacked-bar", required=True)
    p.add_argument("--density-plot", required=True)
    # Task 4
    p.add_argument("--trajectory-counts", required=True)
    p.add_argument("--trajectory-fisher", required=True)
    p.add_argument("--trajectory-plot", required=True)
    # Task 5
    p.add_argument("--g4-correlations", required=True)
    p.add_argument("--g4-direction-tertile", required=True)
    p.add_argument("--g4-scatter", required=True)
    p.add_argument("--g4-violin", required=True)
    # Task 6
    p.add_argument("--expr-matched-wilcoxon", required=True)
    p.add_argument("--expr-matched-lm", required=True)
    p.add_argument("--expr-matched-plot", required=True)
    # Outputs
    p.add_argument("--out-html", required=True)
    p.add_argument("--out-versions", required=True)
    return p.parse_args()


def img_tag(path: str, alt: str = "", width: str = "700px") -> str:
    """Embed a PDF or PNG as a base64-encoded PNG in an <img> tag.
    Falls back to a placeholder if the file cannot be read."""
    p = Path(path)
    if not p.exists():
        return f'<p class="missing">[Figure not found: {p.name}]</p>'
    try:
        if p.suffix.lower() == ".pdf":
            # Render first page of PDF to PNG via matplotlib
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_pdf import PdfPages
            # Read as binary and embed as PDF object if rendering unavailable
            data = p.read_bytes()
            b64 = base64.b64encode(data).decode()
            return (
                f'<embed src="data:application/pdf;base64,{b64}" '
                f'width="{width}" height="500px" type="application/pdf"/>'
            )
        else:
            data = p.read_bytes()
            b64 = base64.b64encode(data).decode()
            ext = p.suffix.lstrip(".").lower()
            mime = {"png": "image/png", "jpg": "image/jpeg", "svg": "image/svg+xml"}.get(ext, "image/png")
            return f'<img src="data:{mime};base64,{b64}" alt="{alt}" style="max-width:{width};"/>'
    except Exception as e:
        return f'<p class="missing">[Could not embed {p.name}: {e}]</p>'


def df_to_html(df: pd.DataFrame, float_fmt: str = ".4f") -> str:
    return df.to_html(
        index=False, border=0, classes="data-table",
        float_format=lambda x: f"{x:{float_fmt}}" if pd.notna(x) else "NA",
    )


def software_versions() -> dict[str, str]:
    pkgs = ["pandas", "numpy", "scipy", "statsmodels", "matplotlib", "seaborn"]
    versions = {}
    for pkg in pkgs:
        try:
            versions[pkg] = importlib.metadata.version(pkg)
        except Exception:
            versions[pkg] = "unknown"
    import platform
    versions["python"] = platform.python_version()
    versions["report_date"] = str(date.today())
    return versions


def build_executive_summary(
    mag_stats: pd.DataFrame,
    dir_fisher: pd.DataFrame,
    traj_fisher: pd.DataFrame,
    g4_corr: pd.DataFrame,
    lm_df: pd.DataFrame,
) -> str:
    # Magnitude at 60 min
    row60 = mag_stats[
        (mag_stats["timepoint"] == "60")
        & (mag_stats["group1"] == "G4_TSS")
        & (mag_stats["group2"] == "No_overlap")
    ]
    mag_sig = ""
    if not row60.empty:
        padj = row60["mw_padj_bh"].iloc[0]
        rb = row60["rank_biserial_r"].iloc[0]
        sig_str = "significant" if (not np.isnan(padj) and padj < 0.05) else "not significant"
        mag_sig = (
            f"At 60 min post-UV, the |log₂FC| magnitude difference between G4_TSS and "
            f"No_overlap genes is <strong>{sig_str}</strong> (BH-adj. p = "
            f"{padj:.3e}, rank-biserial r = {rb:.3f})."
        )

    # Direction at 60 min
    dir60 = dir_fisher[(dir_fisher["timepoint"] == "60")]
    dir_rep = dir60[dir60["test"] == "repression"]
    dir_ind = dir60[dir60["test"] == "induction"]
    dir_str = ""
    if not dir_rep.empty:
        or_r = dir_rep["odds_ratio"].iloc[0]
        p_r = dir_rep["padj_bh"].iloc[0]
        dir_str = (
            f"Repression enrichment in G4_TSS vs No_overlap at 60 min: "
            f"OR = {or_r:.3f}, BH-adj. p = {p_r:.3e}. "
        )
    if not dir_ind.empty:
        or_i = dir_ind["odds_ratio"].iloc[0]
        p_i = dir_ind["padj_bh"].iloc[0]
        dir_str += (
            f"Induction enrichment: OR = {or_i:.3f}, BH-adj. p = {p_i:.3e}."
        )

    # G4 strength correlation
    norm_corr = g4_corr[
        (g4_corr["predictor"] == "combined_norm") & (g4_corr["target"] == "abs_lfc_60")
    ]
    g4_str = ""
    if not norm_corr.empty:
        rho = norm_corr["spearman_rho"].iloc[0]
        pv = norm_corr["pval"].iloc[0]
        g4_str = (
            f"Within G4_TSS genes, G4 signal strength correlates with |lfc_60| "
            f"(Spearman ρ = {rho:.3f}, p = {pv:.3e})."
        )

    # Expression-matched
    lm60 = lm_df[lm_df["timepoint"] == "60"]
    lm_str = ""
    if not lm60.empty:
        coef = lm60["g4_tss_coefficient"].iloc[0]
        pv = lm60["g4_tss_pval"].iloc[0]
        dr2 = lm60["delta_r2_adj"].iloc[0]
        lm_str = (
            f"After controlling for baseline expression, the G4_TSS group coefficient "
            f"at 60 min is {coef:.4f} (p = {pv:.3e}); adding group membership "
            f"increases adjusted R² by {dr2:.4f}."
        )

    para1 = " ".join(filter(None, [mag_sig, g4_str]))
    para2 = " ".join(filter(None, [dir_str, lm_str]))
    return f"<p>{para1}</p>\n<p>{para2}</p>"


CSS = """
body { font-family: Arial, sans-serif; max-width: 1100px; margin: 40px auto; font-size: 14px; line-height: 1.5; }
h1 { color: #1f3a5f; border-bottom: 2px solid #1f3a5f; }
h2 { color: #2c5f8a; border-bottom: 1px solid #aaa; margin-top: 2em; }
h3 { color: #3a7abf; }
.data-table { border-collapse: collapse; width: 100%; margin-bottom: 1em; font-size: 12px; }
.data-table th { background: #e6f0ff; border: 1px solid #ccc; padding: 4px 8px; text-align: left; }
.data-table td { border: 1px solid #eee; padding: 3px 8px; }
.data-table tr:nth-child(even) { background: #f9f9f9; }
.missing { color: #c00; font-style: italic; }
.section { margin-bottom: 2em; }
.figure-block { margin: 1em 0; }
.methods { background: #f5f5f5; border-left: 4px solid #2c5f8a; padding: 10px 16px; font-size: 13px; }
.data-avail { font-size: 12px; }
"""


def main() -> None:
    args = parse_args()
    for p in [args.out_html, args.out_versions]:
        Path(p).parent.mkdir(parents=True, exist_ok=True)

    # Load key result tables
    mag_stats = pd.read_csv(args.magnitude_stats, sep="\t")
    mag_summary = pd.read_csv(args.magnitude_summary, sep="\t")
    dir_stats = pd.read_csv(args.direction_stats, sep="\t")
    dir_fisher = pd.read_csv(args.direction_fisher, sep="\t")
    traj_counts = pd.read_csv(args.trajectory_counts, sep="\t")
    traj_fisher = pd.read_csv(args.trajectory_fisher, sep="\t")
    g4_corr = pd.read_csv(args.g4_correlations, sep="\t")
    g4_tertile = pd.read_csv(args.g4_direction_tertile, sep="\t")
    wil_df = pd.read_csv(args.expr_matched_wilcoxon, sep="\t")
    lm_df = pd.read_csv(args.expr_matched_lm, sep="\t")

    exec_summary = build_executive_summary(mag_stats, dir_fisher, traj_fisher, g4_corr, lm_df)

    # Build output file list
    output_files = [
        ("uv_master_table.tsv", "Task 1 — master table (group, DESeq2 contrasts, trajectories)"),
        ("uv_magnitude_stats.tsv", "Task 2 — Kruskal-Wallis and Mann-Whitney test statistics"),
        ("uv_magnitude_summary.tsv", "Task 2 — descriptive summary of |lfc| per group × timepoint"),
        ("uv_magnitude_by_group.pdf", "Task 2 — violin plots of magnitude"),
        ("uv_direction_stats.tsv", "Task 3 — direction counts and fractions per group × timepoint"),
        ("uv_direction_fisher.tsv", "Task 3 — Fisher exact tests (repression/induction enrichment)"),
        ("uv_direction_stacked_bar.pdf", "Task 3 — stacked bar chart"),
        ("uv_signed_lfc_density.pdf", "Task 3 — signed lfc density plots"),
        ("uv_trajectory_counts.tsv", "Task 4 — trajectory class × group contingency"),
        ("uv_trajectory_fisher.tsv", "Task 4 — Fisher tests (early onset, complex trajectory)"),
        ("uv_trajectory_median_lfc.pdf", "Task 4 — trajectory line plot"),
        ("uv_trajectory_baseline_expr.tsv", "Task 4 — baseline expression within sustained-repression"),
        ("g4_strength_uv_correlations.tsv", "Task 5 — Spearman correlations with bootstrap CIs"),
        ("g4_strength_direction_by_tertile.tsv", "Task 5 — direction stats by G4 signal tertile"),
        ("g4_strength_lfc_scatter.pdf", "Task 5 — scatter plot"),
        ("g4_strength_lfc_violin.pdf", "Task 5 — violin by tertile"),
        ("uv_magnitude_expr_matched.tsv", "Task 6 — within-bin Wilcoxon results"),
        ("uv_magnitude_lm.tsv", "Task 6 — linear model coefficients"),
        ("uv_magnitude_by_bin.pdf", "Task 6 — median |lfc| per expression bin"),
    ]

    html_parts = [
        "<!DOCTYPE html>",
        "<html lang='en'><head><meta charset='UTF-8'>",
        f"<title>G4-UV Magnitude and Direction Report</title>",
        f"<style>{CSS}</style></head><body>",
        "<h1>G4-Promoter Genes: Magnitude and Direction of UV-Induced Transcriptional Change</h1>",
        f"<p><em>Generated: {date.today()}</em></p>",

        "<h2>Executive Summary</h2>",
        exec_summary,

        "<h2>Task 2 — Fold-Change Magnitude Across Promoter Groups</h2>",
        '<div class="section">',
        "<h3>Statistical tests (Kruskal-Wallis + pairwise Mann-Whitney, BH-corrected)</h3>",
        df_to_html(mag_stats),
        "<h3>Descriptive summary</h3>",
        df_to_html(mag_summary, float_fmt=".4f"),
        '<div class="figure-block">',
        img_tag(args.magnitude_figure, "Magnitude by group"),
        "</div></div>",

        "<h2>Task 3 — Repression vs Induction Asymmetry</h2>",
        '<div class="section">',
        "<h3>Direction counts and fractions</h3>",
        df_to_html(dir_stats),
        "<h3>Fisher exact tests (G4_TSS vs No_overlap)</h3>",
        df_to_html(dir_fisher),
        '<div class="figure-block">',
        img_tag(args.stacked_bar, "Stacked bar chart"),
        img_tag(args.density_plot, "Signed LFC density"),
        "</div></div>",

        "<h2>Task 4 — Trajectory Analysis</h2>",
        '<div class="section">',
        "<h3>Trajectory class counts</h3>",
        df_to_html(traj_counts),
        "<h3>Fisher tests (onset timing, complex trajectory)</h3>",
        df_to_html(traj_fisher),
        '<div class="figure-block">',
        img_tag(args.trajectory_plot, "Trajectory line plot"),
        "</div></div>",

        "<h2>Task 5 — G4 Strength as Continuous Predictor</h2>",
        '<div class="section">',
        "<h3>Spearman correlations with bootstrap 95% CIs</h3>",
        df_to_html(g4_corr),
        "<h3>Direction statistics by G4 signal tertile</h3>",
        df_to_html(g4_tertile),
        '<div class="figure-block">',
        img_tag(args.g4_scatter, "G4 strength scatter"),
        img_tag(args.g4_violin, "G4 tertile violin"),
        "</div></div>",

        "<h2>Task 6 — Baseline-Expression-Matched Comparison</h2>",
        '<div class="section">',
        "<h3>Within-bin Wilcoxon results (|lfc_60|, BH-corrected)</h3>",
        df_to_html(wil_df),
        "<h3>Linear model: |lfc| ~ log2_tpm_t0 + group</h3>",
        df_to_html(lm_df),
        '<div class="figure-block">',
        img_tag(args.expr_matched_plot, "Expression-matched bin plot"),
        "</div></div>",

        "<h2>Methods</h2>",
        '<div class="methods">',
        "<p><strong>Sign convention:</strong> All pairwise DESeq2 contrasts are coded as "
        "<em>0_vs_X</em> (numerator = time 0, denominator = post-UV timepoint). "
        "A positive log₂FoldChange indicates higher expression at t=0 than post-UV "
        "(<strong>repression</strong>); negative indicates "
        "<strong>induction</strong>.</p>",
        "<p><strong>Significance thresholds:</strong> LRT FDR &lt; 0.05 for omnibus "
        "time-series significance; pairwise |lfc| &gt; 0.5 with padj &lt; 0.05 for "
        "directional calls.</p>",
        "<p><strong>Multiple testing:</strong> Benjamini-Hochberg FDR correction applied "
        "within each task: 9 pairwise tests (Task 2), 6 Fisher tests (Task 3), "
        "3 Fisher tests (Task 4), 10 within-bin tests (Task 6).</p>",
        "<p><strong>Promoter groups:</strong> G4_TSS = gene has G4 peak overlapping "
        "±500 bp TSS window; GC_bg_TSS = sampled GC-rich canonical promoter control "
        "with no merged-G4 or OQS overlap; No_overlap = neither.</p>",
        "<p><strong>Color scheme:</strong> G4_TSS = blue (#1f77b4), "
        "GC_bg_TSS = orange (#ff7f0e), No_overlap = gray (#7f7f7f).</p>",
        "</div>",

        "<h2>Data Availability</h2>",
        '<div class="data-avail"><ul>',
        *[f"<li><code>results/g4_tss_uv/{f}</code> — {desc}</li>" for f, desc in output_files],
        "</ul></div>",

        "</body></html>",
    ]

    html_str = "\n".join(html_parts)
    Path(args.out_html).write_text(html_str, encoding="utf-8")
    print(f"Written: {args.out_html}")

    # Software versions
    vers = software_versions()
    lines = [f"{k}: {v}" for k, v in vers.items()]
    Path(args.out_versions).write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Written: {args.out_versions}")


if __name__ == "__main__":
    main()
