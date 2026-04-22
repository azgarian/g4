#!/usr/bin/env python3
"""Task 8: Generate a self-contained HTML report summarising all tasks."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd


TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>G4s at TSSs and Active Transcription</title>
<style>
  body {{ font-family: Arial, sans-serif; margin: 40px; max-width: 1100px; }}
  h1 {{ color: #333; }}
  h2 {{ color: #555; border-bottom: 1px solid #ccc; }}
  table {{ border-collapse: collapse; width: 100%; margin-bottom: 1em; }}
  th, td {{ border: 1px solid #ccc; padding: 6px 10px; text-align: left; }}
  th {{ background: #f5f5f5; }}
  .note {{ background: #fffbe6; padding: 8px 14px; border-left: 4px solid #f0c040; }}
  img {{ max-width: 100%; }}
  .figure {{ margin: 1em 0; text-align: center; }}
  .figure img {{ border: 1px solid #ddd; }}
</style>
</head>
<body>
<h1>Are G4s at TSSs Linked to Active Transcription?</h1>
<p class="note"><strong>Summary:</strong> {executive_summary}</p>

<h2>1. Group sizes (TSS annotation)</h2>
{group_size_table}

<h2>2. Baseline expression by group (Task 3)</h2>
{expr_summary_table}
{expr_stats_table}

<h2>3. Fisher exact — expression prevalence (Task 3)</h2>
{fisher_table}

<h2>4. Expression violin (Task 4)</h2>
<div class="figure"><img src="{violin_plot}"></div>

<h2>5. RNA-seq TSS metaprofile (Task 4)</h2>
<div class="figure"><img src="{metaprofile_plot}"></div>

<h2>6. G4 overlap fraction across expression deciles (Task 5)</h2>
{decile_corr_table}
<div class="figure"><img src="{decile_plot}"></div>

<h2>7. UV response — omnibus LRT (Task 6)</h2>
{lrt_table}

<h2>8. UV response — fold change by group (Task 6)</h2>
{fc_stats_table}
<div class="figure"><img src="{fc_plot}"></div>

<h2>9. UV volcano (Task 6)</h2>
<div class="figure"><img src="{volcano_plot}"></div>

<h2>10. G4 structure class enrichment (Task 7)</h2>
{enrichment_table}
<div class="figure"><img src="{structure_plot}"></div>

<h2>Software versions</h2>
<pre>{software_versions}</pre>

</body>
</html>
"""


def df_to_html(df: pd.DataFrame) -> str:
    return df.round(4).fillna("NA").to_html(index=False, border=0)


def read_safe(path: str) -> pd.DataFrame | None:
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return None


def get_software_versions() -> str:
    lines = []
    for tool, cmd in [
        ("python", ["python3", "--version"]),
        ("pandas", ["python3", "-c", "import pandas; print('pandas', pandas.__version__)"]),
        ("bedtools", ["bedtools", "--version"]),
        ("deeptools", ["deeptools", "--version"]),
    ]:
        try:
            out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True).strip()
            lines.append(out)
        except Exception:
            lines.append(f"{tool}: not found")
    return "\n".join(lines)


def build_executive_summary(expr_summary, expr_stats, decile_corr, lrt_summary, fc_stats) -> str:
    parts = []

    if expr_summary is not None and expr_stats is not None:
        g4_median = expr_summary[expr_summary["group"] == "G4_TSS"]["median"].values
        gc_median = expr_summary[expr_summary["group"] == "GC_bg_TSS"]["median"].values
        no_median = expr_summary[expr_summary["group"] == "No_overlap"]["median"].values
        mw_row = expr_stats[expr_stats["comparison"] == "G4_TSS_vs_GC_bg_TSS"]
        if len(g4_median) and len(gc_median) and len(no_median) and not mw_row.empty:
            sig = "significantly" if mw_row["p_adj_BH"].values[0] < 0.05 else "not significantly"
            parts.append(
                f"G4_TSS genes are {sig} more highly expressed at baseline than GC_bg_TSS "
                f"(median log₂: {g4_median[0]:.2f} vs {gc_median[0]:.2f} vs {no_median[0]:.2f} for No_overlap; "
                f"BH-adjusted p = {mw_row['p_adj_BH'].values[0]:.2e})."
            )

    if decile_corr is not None:
        g4_rho = decile_corr[decile_corr["track"] == "G4_merged"]["spearman_rho"].values
        if len(g4_rho):
            direction = "increases" if g4_rho[0] > 0 else "does not increase"
            parts.append(
                f"G4 occupancy at promoters {direction} across expression deciles "
                f"(Spearman ρ = {g4_rho[0]:.3f})."
            )

    if lrt_summary is not None:
        g4_frac = lrt_summary[lrt_summary["group"] == "G4_TSS"]["fraction_significant_lrt"].values
        no_frac = lrt_summary[lrt_summary["group"] == "No_overlap"]["fraction_significant_lrt"].values
        if len(g4_frac) and len(no_frac):
            parts.append(
                f"G4_TSS genes show a {g4_frac[0]:.1%} LRT-significant UV response rate "
                f"vs {no_frac[0]:.1%} for No_overlap genes."
            )

    if not parts:
        parts = ["Results summary pending full run completion."]

    return " ".join(parts)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tss-annotation", required=True)
    p.add_argument("--expr-summary", required=True)
    p.add_argument("--expr-statistics", required=True)
    p.add_argument("--expr-fisher", required=True)
    p.add_argument("--violin-plot", required=True)
    p.add_argument("--metaprofile-plot", required=True)
    p.add_argument("--decile-overlap", required=True)
    p.add_argument("--decile-corr", required=True)
    p.add_argument("--decile-plot", required=True)
    p.add_argument("--lrt-summary", required=True)
    p.add_argument("--fc-stats", required=True)
    p.add_argument("--fc-plot", required=True)
    p.add_argument("--volcano-plot", required=True)
    p.add_argument("--enrichment-stats", required=True)
    p.add_argument("--structure-plot", required=True)
    p.add_argument("--out-html", required=True)
    p.add_argument("--out-versions", required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    for outpath in [args.out_html, args.out_versions]:
        Path(outpath).parent.mkdir(parents=True, exist_ok=True)

    versions = get_software_versions()
    Path(args.out_versions).write_text(versions)

    anno = read_safe(args.tss_annotation)
    expr_summary = read_safe(args.expr_summary)
    expr_stats = read_safe(args.expr_statistics)
    expr_fisher = read_safe(args.expr_fisher)
    decile_overlap = read_safe(args.decile_overlap)
    decile_corr = read_safe(args.decile_corr)
    lrt_summary = read_safe(args.lrt_summary)
    fc_stats = read_safe(args.fc_stats)
    enrichment = read_safe(args.enrichment_stats)

    group_size = anno["group"].value_counts().reset_index() if anno is not None else pd.DataFrame()
    group_size.columns = ["group", "n_genes"]

    exec_summary = build_executive_summary(expr_summary, expr_stats, decile_corr, lrt_summary, fc_stats)

    # Use relative paths for images so HTML is portable within results/g4_tss/
    def rel(path: str) -> str:
        try:
            return Path(path).name
        except Exception:
            return path

    html = TEMPLATE.format(
        executive_summary=exec_summary,
        group_size_table=df_to_html(group_size) if not group_size.empty else "<p>N/A</p>",
        expr_summary_table=df_to_html(expr_summary) if expr_summary is not None else "<p>N/A</p>",
        expr_stats_table=df_to_html(expr_stats) if expr_stats is not None else "<p>N/A</p>",
        fisher_table=df_to_html(expr_fisher) if expr_fisher is not None else "<p>N/A</p>",
        violin_plot=rel(args.violin_plot),
        metaprofile_plot=rel(args.metaprofile_plot),
        decile_corr_table=df_to_html(decile_corr) if decile_corr is not None else "<p>N/A</p>",
        decile_plot=rel(args.decile_plot),
        lrt_table=df_to_html(lrt_summary) if lrt_summary is not None else "<p>N/A</p>",
        fc_stats_table=df_to_html(fc_stats) if fc_stats is not None else "<p>N/A</p>",
        fc_plot=rel(args.fc_plot),
        volcano_plot=rel(args.volcano_plot),
        enrichment_table=df_to_html(enrichment.head(20)) if enrichment is not None else "<p>N/A</p>",
        structure_plot=rel(args.structure_plot),
        software_versions=versions,
    )

    Path(args.out_html).write_text(html)
    print(f"Written: {args.out_html}", file=sys.stdout)


if __name__ == "__main__":
    main()
