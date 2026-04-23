#!/usr/bin/env python3
"""Task 7: Synthesis HTML report.

Produces a self-contained HTML report answering the central question:
Do genes with high UV lesion burden at their promoter G4 show stronger
transcriptional suppression after UV irradiation?

Embeds all key figures as base64 and includes summary tables from all tasks.
"""

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
    # Task 1
    p.add_argument("--locus-map", required=True)
    p.add_argument("--gene-table", required=True)
    # Task 2
    p.add_argument("--cpd-table", required=True)
    p.add_argument("--pp64-table", required=True)
    # Task 3
    p.add_argument("--cpd-tertile-stats", required=True)
    p.add_argument("--pp64-tertile-stats", required=True)
    p.add_argument("--cpd-tertile-summary", required=True)
    p.add_argument("--pp64-tertile-summary", required=True)
    p.add_argument("--cpd-lfc-fig", required=True)
    p.add_argument("--pp64-lfc-fig", required=True)
    p.add_argument("--cpd-traj-fig", required=True)
    p.add_argument("--pp64-traj-fig", required=True)
    p.add_argument("--cpd-comp-fig", required=True)
    p.add_argument("--pp64-comp-fig", required=True)
    # Task 4
    p.add_argument("--cpd-corr", required=True)
    p.add_argument("--pp64-corr", required=True)
    p.add_argument("--cpd-partial", required=True)
    p.add_argument("--pp64-partial", required=True)
    p.add_argument("--cpd-models", required=True)
    p.add_argument("--pp64-models", required=True)
    p.add_argument("--cpd-model-cmp", required=True)
    p.add_argument("--pp64-model-cmp", required=True)
    p.add_argument("--cpd-scatter", required=True)
    p.add_argument("--pp64-scatter", required=True)
    p.add_argument("--cpd-partial-fig", required=True)
    p.add_argument("--pp64-partial-fig", required=True)
    p.add_argument("--cpd-coef-fig", required=True)
    p.add_argument("--pp64-coef-fig", required=True)
    # Task 5
    p.add_argument("--cpd-traj-kruskal", required=True)
    p.add_argument("--pp64-traj-kruskal", required=True)
    p.add_argument("--cpd-traj-pairwise", required=True)
    p.add_argument("--pp64-traj-pairwise", required=True)
    p.add_argument("--cpd-traj-violin", required=True)
    p.add_argument("--pp64-traj-violin", required=True)
    p.add_argument("--cpd-traj-logistic", required=True)
    p.add_argument("--pp64-traj-logistic", required=True)
    # Task 6
    p.add_argument("--interaction", required=True)
    p.add_argument("--interaction-fig", required=True)
    p.add_argument("--gc-bg-summary", required=True)
    p.add_argument("--gc-bg-vs-g4", required=True)
    # Output
    p.add_argument("--out-html", required=True)
    p.add_argument("--out-versions", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def embed_figure(path: str, alt: str = "", width: str = "680px") -> str:
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return f'<p class="missing">[Figure not available: {p.name}]</p>'
    try:
        data = p.read_bytes()
        b64 = base64.b64encode(data).decode()
        suffix = p.suffix.lower()
        if suffix == ".pdf":
            mime = "application/pdf"
            return (
                f'<embed src="data:{mime};base64,{b64}" '
                f'type="{mime}" width="{width}" height="500px" '
                f'alt="{alt}">'
            )
        elif suffix in (".png", ".jpg", ".jpeg"):
            mime = "image/png" if suffix == ".png" else "image/jpeg"
            return f'<img src="data:{mime};base64,{b64}" width="{width}" alt="{alt}">'
        else:
            return f'<p class="missing">[Unsupported figure format: {p.suffix}]</p>'
    except Exception as exc:
        return f'<p class="missing">[Figure embed failed: {exc}]</p>'


def df_to_html(path: str, max_rows: int = 20) -> str:
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return "<p><em>Table not available.</em></p>"
    try:
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            return "<p><em>Table is empty.</em></p>"
        df = df.head(max_rows)
        # Round floats
        df = df.applymap(lambda v: f"{v:.4g}" if isinstance(v, float) else v)
        return df.to_html(index=False, border=0, classes="data-table")
    except Exception as exc:
        return f"<p><em>Could not load table: {exc}</em></p>"


def safe_read(path: str) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.DataFrame()


def executive_summary(
    cpd_corr: pd.DataFrame,
    pp64_corr: pd.DataFrame,
    cpd_models: pd.DataFrame,
    pp64_models: pd.DataFrame,
    cpd_kruskal: pd.DataFrame,
    pp64_kruskal: pd.DataFrame,
    interaction: pd.DataFrame,
) -> str:
    lines = []

    def _spearman_lfc60(df: pd.DataFrame) -> str:
        if df.empty or "spearman_r" not in df.columns:
            return "N/A"
        row = df[df.get("target", pd.Series(dtype=object)) == "lfc_60"]
        if row.empty:
            row = df[df.get("comparison", pd.Series(dtype=object)).str.contains("lfc_60", na=False)]
        if row.empty:
            return "N/A"
        r = row.iloc[0].get("spearman_r", np.nan)
        p = row.iloc[0].get("pvalue_bh", row.iloc[0].get("pvalue_raw", np.nan))
        return f"r={r:.3f}, BH-p={p:.2e}" if not np.isnan(r) else "N/A"

    def _damage_coef(df: pd.DataFrame, outcome: str = "lfc_60") -> str:
        if df.empty or "coef_damage" not in df.columns:
            return "N/A"
        row = df[df.get("outcome", pd.Series(dtype=object)) == outcome]
        if row.empty:
            return "N/A"
        coef = row.iloc[0].get("coef_damage", np.nan)
        pval = row.iloc[0].get("pvalue_damage", np.nan)
        return f"β={coef:.3f}, p={pval:.2e}" if not np.isnan(coef) else "N/A"

    def _kruskal_p(df: pd.DataFrame) -> str:
        if df.empty or "kruskal_pvalue" not in df.columns:
            return "N/A"
        p = df["kruskal_pvalue"].iloc[0]
        return f"p={p:.2e}"

    lines.append("<ul>")
    lines.append(
        f"<li><strong>Continuous association (damage_ds0 vs LFC_60):</strong> "
        f"CPD: {_spearman_lfc60(cpd_corr)}; "
        f"64-PP: {_spearman_lfc60(pp64_corr)}</li>"
    )
    lines.append(
        f"<li><strong>Covariate-adjusted linear model (lfc_60):</strong> "
        f"CPD: {_damage_coef(cpd_models, 'lfc_60')}; "
        f"64-PP: {_damage_coef(pp64_models, 'lfc_60')}</li>"
    )
    lines.append(
        f"<li><strong>Trajectory Kruskal-Wallis:</strong> "
        f"CPD: {_kruskal_p(cpd_kruskal)}; "
        f"64-PP: {_kruskal_p(pp64_kruskal)}</li>"
    )
    if not interaction.empty and "pvalue" in interaction.columns:
        sig_inter = interaction[interaction["pvalue"] < 0.05]
        lines.append(
            f"<li><strong>Promoter-group specificity (interaction term):</strong> "
            f"{len(sig_inter)} of {len(interaction)} tests significant at p&lt;0.05</li>"
        )
    lines.append("</ul>")
    return "\n".join(lines)


def software_versions() -> dict[str, str]:
    pkgs = ["pandas", "numpy", "scipy", "statsmodels", "matplotlib", "seaborn"]
    versions = {}
    for pkg in pkgs:
        try:
            versions[pkg] = importlib.metadata.version(pkg)
        except Exception:
            versions[pkg] = "unknown"
    import subprocess
    try:
        r = subprocess.run(["bedtools", "--version"], capture_output=True, text=True)
        versions["bedtools"] = r.stdout.strip().split("\n")[0]
    except Exception:
        versions["bedtools"] = "unknown"
    return versions


CSS = """
body { font-family: Arial, sans-serif; max-width: 1100px; margin: 40px auto; color: #222; }
h1 { color: #1565C0; }
h2 { color: #1976D2; border-bottom: 1px solid #ccc; padding-bottom: 4px; margin-top: 36px; }
h3 { color: #333; margin-top: 24px; }
.data-table { border-collapse: collapse; font-size: 12px; width: 100%; margin: 12px 0; }
.data-table th { background: #E3F2FD; border: 1px solid #bbb; padding: 4px 8px; text-align: left; }
.data-table td { border: 1px solid #ddd; padding: 3px 8px; }
.data-table tr:nth-child(even) { background: #fafafa; }
.missing { color: #999; font-style: italic; font-size: 12px; }
.sign-convention { background: #FFF9C4; border-left: 4px solid #F9A825; padding: 8px 12px;
                   margin: 12px 0; font-size: 13px; }
.product-section { border: 1px solid #E0E0E0; border-radius: 6px; padding: 16px;
                   margin: 16px 0; background: #FAFAFA; }
.fig-row { display: flex; gap: 20px; flex-wrap: wrap; margin: 12px 0; }
.fig-cell { flex: 1; min-width: 300px; }
.summary-box { background: #E8F5E9; border-left: 4px solid #388E3C; padding: 10px 14px;
               margin: 12px 0; }
"""


def build_html(args: argparse.Namespace) -> str:
    # Load tables
    locus_map = safe_read(args.locus_map)
    gene_table = safe_read(args.gene_table)
    cpd_table = safe_read(args.cpd_table)
    pp64_table = safe_read(args.pp64_table)
    cpd_corr = safe_read(args.cpd_corr)
    pp64_corr = safe_read(args.pp64_corr)
    cpd_models = safe_read(args.cpd_models)
    pp64_models = safe_read(args.pp64_models)
    cpd_kruskal = safe_read(args.cpd_traj_kruskal)
    pp64_kruskal = safe_read(args.pp64_traj_kruskal)
    interaction = safe_read(args.interaction)
    gc_summary = safe_read(args.gc_bg_summary)
    gc_vs_g4 = safe_read(args.gc_bg_vs_g4)

    exec_sum = executive_summary(
        cpd_corr, pp64_corr, cpd_models, pp64_models,
        cpd_kruskal, pp64_kruskal, interaction,
    )

    # Compute simple header stats
    n_loci = locus_map["promoter_g4_locus_id"].nunique() if not locus_map.empty else "N/A"
    n_g4_tss_genes = (
        gene_table[gene_table["group"] == "G4_TSS"]["gene_id"].nunique()
        if not gene_table.empty and "group" in gene_table.columns else "N/A"
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>G4 Promoter UV Lesion Burden and Transcriptional Suppression</title>
<style>{CSS}</style>
</head>
<body>

<h1>G4 Promoter UV Lesion Burden and Transcriptional Suppression</h1>
<p><strong>Date:</strong> {date.today().isoformat()}</p>
<p><strong>Pipeline:</strong> g4_tss_damage_uv Tasks 1–7</p>

<div class="sign-convention">
<strong>Sign convention:</strong> positive log2FoldChange in 0_vs_T contrasts means expression
is higher at time 0 than at time T &rarr; positive LFC = <em>repression after UV</em>.
Primary lesion predictor: <code>damage_ds0</code> = max log2(real/sim) RPKM
at DS time_after_exposure = 0. Products analyzed separately: CPD and 64-PP.
</div>

<h2>Executive Summary</h2>
<div class="summary-box">
{exec_sum}
</div>

<!-- ═══════════════════════════════════════════════════════════════════ -->
<h2>Task 1 — Promoter-G4 Locus Definition and Lesion Burden</h2>

<p>Promoter-G4 loci were defined by intersecting canonical 1-kb TSS windows with the merged
ChIP + CUT&amp;Tag G4 peak set. Lesion burden was quantified as
log2(real_RPKM / sim_RPKM) per locus per DS sample, then collapsed to gene-level
<code>damage_max</code>.</p>

<table class="data-table">
  <tr><th>Metric</th><th>Value</th></tr>
  <tr><td>Total promoter-G4 loci (gene × locus pairs)</td><td>{n_loci}</td></tr>
  <tr><td>G4_TSS genes with ≥1 quantified locus</td><td>{n_g4_tss_genes}</td></tr>
</table>

<h3>Promoter-G4 Gene–Locus Map (first 20 rows)</h3>
{df_to_html(args.locus_map)}

<h3>Gene-level Burden Table — DS0 (first 20 rows)</h3>
{df_to_html(args.gene_table)}

<!-- ═══════════════════════════════════════════════════════════════════ -->
<h2>Task 2 — Lesion × RNA Master Tables</h2>

<p>DS0 burden joined to the UV-response master table (LFC at 12, 30, 60 min post-UV),
baseline expression (log2 TPM), and G4 signal covariates.</p>

<h3>CPD Master Table (first 20 rows)</h3>
{df_to_html(args.cpd_table)}

<h3>64-PP Master Table (first 20 rows)</h3>
{df_to_html(args.pp64_table)}

<!-- ═══════════════════════════════════════════════════════════════════ -->
<h2>Task 3 — Primary Burden-Stratified Test</h2>

<p>G4_TSS genes ranked by <code>damage_ds0</code> and binned into equal-count tertiles
(low / mid / high). Kruskal-Wallis and one-sided Mann-Whitney tests on signed LFC;
chi-squared and Fisher tests on response categories.</p>

<div class="product-section">
<h3>CPD</h3>
<h4>Tertile Summary</h4>
{df_to_html(args.cpd_tertile_summary)}
<h4>Statistical Tests</h4>
{df_to_html(args.cpd_tertile_stats)}
<div class="fig-row">
  <div class="fig-cell">
    <p><em>LFC by tertile (violin)</em></p>
    {embed_figure(args.cpd_lfc_fig, "CPD LFC by tertile")}
  </div>
  <div class="fig-cell">
    <p><em>Median LFC trajectory</em></p>
    {embed_figure(args.cpd_traj_fig, "CPD trajectory")}
  </div>
</div>
<p><em>Trajectory composition by tertile</em></p>
{embed_figure(args.cpd_comp_fig, "CPD composition")}
</div>

<div class="product-section">
<h3>64-PP</h3>
<h4>Tertile Summary</h4>
{df_to_html(args.pp64_tertile_summary)}
<h4>Statistical Tests</h4>
{df_to_html(args.pp64_tertile_stats)}
<div class="fig-row">
  <div class="fig-cell">
    <p><em>LFC by tertile (violin)</em></p>
    {embed_figure(args.pp64_lfc_fig, "64-PP LFC by tertile")}
  </div>
  <div class="fig-cell">
    <p><em>Median LFC trajectory</em></p>
    {embed_figure(args.pp64_traj_fig, "64-PP trajectory")}
  </div>
</div>
<p><em>Trajectory composition by tertile</em></p>
{embed_figure(args.pp64_comp_fig, "64-PP composition")}
</div>

<!-- ═══════════════════════════════════════════════════════════════════ -->
<h2>Task 4 — Continuous Association and Covariate-Adjusted Modeling</h2>

<p>Spearman correlations between <code>damage_ds0</code> and LFC at each timepoint.
Linear models: <code>lfc_t ~ damage_ds0_z + log2_tpm_t0 + max_g4_signal_norm</code>.
Logistic models for repression indicators. All predictors standardized before fitting.</p>

<div class="product-section">
<h3>CPD</h3>
<h4>Spearman Correlations</h4>
{df_to_html(args.cpd_corr)}
<h4>Partial Spearman (damage vs LFC_60 | log2_tpm_t0)</h4>
{df_to_html(args.cpd_partial)}
<h4>Covariate Models</h4>
{df_to_html(args.cpd_models)}
<h4>Model Comparison (with vs without damage)</h4>
{df_to_html(args.cpd_model_cmp)}
<div class="fig-row">
  <div class="fig-cell">
    <p><em>damage_ds0 vs LFC_60 scatter</em></p>
    {embed_figure(args.cpd_scatter, "CPD scatter")}
  </div>
  <div class="fig-cell">
    <p><em>Partial effect of burden on LFC_60</em></p>
    {embed_figure(args.cpd_partial_fig, "CPD partial effect")}
  </div>
</div>
<p><em>Coefficient plot (damage_ds0_z across models)</em></p>
{embed_figure(args.cpd_coef_fig, "CPD coefficients")}
</div>

<div class="product-section">
<h3>64-PP</h3>
<h4>Spearman Correlations</h4>
{df_to_html(args.pp64_corr)}
<h4>Partial Spearman (damage vs LFC_60 | log2_tpm_t0)</h4>
{df_to_html(args.pp64_partial)}
<h4>Covariate Models</h4>
{df_to_html(args.pp64_models)}
<h4>Model Comparison (with vs without damage)</h4>
{df_to_html(args.pp64_model_cmp)}
<div class="fig-row">
  <div class="fig-cell">
    <p><em>damage_ds0 vs LFC_60 scatter</em></p>
    {embed_figure(args.pp64_scatter, "64-PP scatter")}
  </div>
  <div class="fig-cell">
    <p><em>Partial effect of burden on LFC_60</em></p>
    {embed_figure(args.pp64_partial_fig, "64-PP partial effect")}
  </div>
</div>
<p><em>Coefficient plot (damage_ds0_z across models)</em></p>
{embed_figure(args.pp64_coef_fig, "64-PP coefficients")}
</div>

<!-- ═══════════════════════════════════════════════════════════════════ -->
<h2>Task 5 — Trajectory-Level and Sustained-Repression Analysis</h2>

<p>Kruskal-Wallis test across all UV trajectory classes; pairwise Mann-Whitney tests
comparing repressed_sustained vs not_responsive, induced_sustained, and lrt_only.
Logistic model restricted to {{repressed_sustained, not_responsive}} genes.</p>

<div class="product-section">
<h3>CPD</h3>
<h4>Kruskal-Wallis</h4>
{df_to_html(args.cpd_traj_kruskal)}
<h4>Pairwise Tests</h4>
{df_to_html(args.cpd_traj_pairwise)}
<h4>Logistic Model (repressed_sustained vs not_responsive)</h4>
{df_to_html(args.cpd_traj_logistic)}
<p><em>Burden by trajectory (violin)</em></p>
{embed_figure(args.cpd_traj_violin, "CPD trajectory violin")}
</div>

<div class="product-section">
<h3>64-PP</h3>
<h4>Kruskal-Wallis</h4>
{df_to_html(args.pp64_traj_kruskal)}
<h4>Pairwise Tests</h4>
{df_to_html(args.pp64_traj_pairwise)}
<h4>Logistic Model (repressed_sustained vs not_responsive)</h4>
{df_to_html(args.pp64_traj_logistic)}
<p><em>Burden by trajectory (violin)</em></p>
{embed_figure(args.pp64_traj_violin, "64-PP trajectory violin")}
</div>

<!-- ═══════════════════════════════════════════════════════════════════ -->
<h2>Task 6 — Specificity and Sensitivity Checks</h2>

<p>Specificity: interaction model testing whether the damage effect is stronger in G4_TSS
than in No_overlap genes matched by expression decile. Sensitivity: alternate burden metrics
(damage_mean, damage_sum, real_rpkm_max). GC-rich background: descriptive control.</p>

<h3>Promoter-Group Interaction Results</h3>
{df_to_html(args.interaction)}
<p><em>Interaction plot</em></p>
{embed_figure(args.interaction_fig, "Interaction plot")}

<h3>GC-rich Background Control</h3>
{df_to_html(args.gc_bg_summary)}
<h4>GC-rich vs G4_TSS Burden Comparison</h4>
{df_to_html(args.gc_bg_vs_g4)}

<!-- ═══════════════════════════════════════════════════════════════════ -->
<h2>Methods Notes</h2>
<ul>
  <li><strong>Lesion input:</strong> DS (damage-seq) BED files (plus and minus strand),
      strand-collapsed by concatenation. Simulated reads = random placements of the same
      library, used as the null for RPKM normalization.</li>
  <li><strong>Burden metric:</strong> log2(real_RPKM / sim_RPKM) per locus,
      pseudocount 1e-3 added to each RPKM. Gene-level: damage_max across loci.</li>
  <li><strong>Primary predictor:</strong> <code>damage_ds0</code> = damage_max at
      DS time_after_exposure = 0, capturing lesion formation before RNA response onset.</li>
  <li><strong>Sign convention:</strong> positive LFC = repression after UV.</li>
  <li><strong>Multiple testing:</strong> BH (Benjamini-Hochberg) correction applied
      within each product and test family.</li>
  <li><strong>Caveats:</strong> DS0 and t=0 RNA are from separate assays; their
      relationship is correlational. Baseline expression and G4 occupancy strength are
      potential confounders and are included in adjusted models. The GC-rich background
      control may be underpowered if the group is small.</li>
</ul>

<!-- ═══════════════════════════════════════════════════════════════════ -->
<h2>Software Versions</h2>
<p id="versions"><em>See versions file.</em></p>

</body>
</html>
"""
    return html


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else sys.stdout

    def log(msg: str) -> None:
        print(msg, file=log_fh, flush=True)

    Path(args.out_html).parent.mkdir(parents=True, exist_ok=True)

    log("=== Task 7: Synthesis HTML report ===")

    html = build_html(args)
    Path(args.out_html).write_text(html, encoding="utf-8")
    log(f"Written: {args.out_html}")

    # Software versions
    vers = software_versions()
    lines = [f"{pkg}: {ver}" for pkg, ver in sorted(vers.items())]
    Path(args.out_versions).write_text("\n".join(lines) + "\n", encoding="utf-8")
    log(f"Written: {args.out_versions}")

    log("\n=== Task 7 complete ===")
    if args.log:
        log_fh.close()


if __name__ == "__main__":
    main()
