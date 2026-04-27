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
<title>{page_title}</title>
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
<h1>{page_title}</h1>
<p class="note"><strong>Summary:</strong> {executive_summary}</p>
{sections}

<h2>Software versions</h2>
<pre>{software_versions}</pre>

</body>
</html>
"""


def df_to_html(df: pd.DataFrame) -> str:
    return df.round(4).fillna("NA").to_html(index=False, border=0)


def read_safe(path: str) -> pd.DataFrame | None:
    if path is None:
        return None
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return None


def render_section(title: str, body: str) -> str:
    return f"<h2>{title}</h2>\n{body}"


def render_figure(path: str | None) -> str:
    if path is None:
        return "<p>N/A</p>"
    return f'<div class="figure"><img src="{Path(path).name}"></div>'


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


def build_executive_summary(
    expr_summary: pd.DataFrame | None,
    expr_stats: pd.DataFrame | None,
    decile_corr: pd.DataFrame | None,
    lrt_summary: pd.DataFrame | None,
    expression_label: str,
) -> str:
    parts = []
    expression_context = (
        "at baseline"
        if expression_label.lower().startswith("baseline")
        else "at the matched RNA-seq timepoint"
    )

    if expr_summary is not None and expr_stats is not None:
        g4_median = expr_summary[expr_summary["group"] == "G4_TSS"]["median"].values
        gc_median = expr_summary[expr_summary["group"] == "GC_bg_TSS"]["median"].values
        no_median = expr_summary[expr_summary["group"] == "No_overlap"]["median"].values
        mw_row = expr_stats[expr_stats["comparison"] == "G4_TSS_vs_GC_bg_TSS"]
        if len(g4_median) and len(gc_median) and len(no_median) and not mw_row.empty:
            sig = "significantly" if mw_row["p_adj_BH"].values[0] < 0.05 else "not significantly"
            parts.append(
                f"G4_TSS genes are {sig} more highly expressed {expression_context} than GC_bg_TSS "
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
    p.add_argument("--lrt-summary", default=None)
    p.add_argument("--fc-stats", default=None)
    p.add_argument("--fc-plot", default=None)
    p.add_argument("--volcano-plot", default=None)
    p.add_argument("--lrt-sig-fc-stats", default=None)
    p.add_argument("--lrt-sig-fc-plot", default=None)
    p.add_argument("--lrt-sig-volcano-plot", default=None)
    p.add_argument("--enrichment-stats", required=True)
    p.add_argument("--structure-plot", required=True)
    p.add_argument("--analysis-label", default=None)
    p.add_argument("--expression-label", default="Baseline expression")
    p.add_argument("--filter-summary", default=None)
    p.add_argument("--uv-label", default="0 vs 60 min")
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
    filter_summary = read_safe(args.filter_summary)
    decile_overlap = read_safe(args.decile_overlap)
    decile_corr = read_safe(args.decile_corr)
    lrt_summary = read_safe(args.lrt_summary)
    fc_stats = read_safe(args.fc_stats)
    lrt_sig_fc_stats = read_safe(args.lrt_sig_fc_stats)
    enrichment = read_safe(args.enrichment_stats)

    group_size = anno["group"].value_counts().reset_index() if anno is not None else pd.DataFrame()
    group_size.columns = ["group", "n_genes"]

    exec_summary = build_executive_summary(
        expr_summary, expr_stats, decile_corr, lrt_summary, args.expression_label
    )
    page_title = "Are G4s at TSSs Linked to Active Transcription?"
    if args.analysis_label:
        page_title = f"{page_title} ({args.analysis_label})"

    sections = []
    section_index = 1

    if filter_summary is not None:
        sections.append(
            render_section(
                f"{section_index}. ATAC-filtered G4 summary",
                df_to_html(filter_summary),
            )
        )
        section_index += 1

    sections.append(
        render_section(
            f"{section_index}. Group sizes (TSS annotation)",
            df_to_html(group_size) if not group_size.empty else "<p>N/A</p>",
        )
    )
    section_index += 1

    expr_tables = [
        df_to_html(expr_summary) if expr_summary is not None else "<p>N/A</p>",
        df_to_html(expr_stats) if expr_stats is not None else "<p>N/A</p>",
    ]
    sections.append(
        render_section(
            f"{section_index}. {args.expression_label} by group",
            "\n".join(expr_tables),
        )
    )
    section_index += 1

    sections.append(
        render_section(
            f"{section_index}. Fisher exact — expression prevalence",
            df_to_html(expr_fisher) if expr_fisher is not None else "<p>N/A</p>",
        )
    )
    section_index += 1

    sections.append(render_section(f"{section_index}. Expression violin", render_figure(args.violin_plot)))
    section_index += 1

    sections.append(
        render_section(f"{section_index}. RNA-seq TSS metaprofile", render_figure(args.metaprofile_plot))
    )
    section_index += 1

    sections.append(
        render_section(
            f"{section_index}. G4 overlap fraction across expression deciles",
            (
                f"{df_to_html(decile_corr) if decile_corr is not None else '<p>N/A</p>'}\n"
                f"{render_figure(args.decile_plot)}"
            ),
        )
    )
    section_index += 1

    if lrt_summary is not None and fc_stats is not None and args.fc_plot and args.volcano_plot:
        sections.append(
            render_section(
                f"{section_index}. UV response — omnibus LRT ({args.uv_label})",
                df_to_html(lrt_summary),
            )
        )
        section_index += 1

        sections.append(
            render_section(
                f"{section_index}. UV response — fold change by group ({args.uv_label})",
                f"{df_to_html(fc_stats)}\n{render_figure(args.fc_plot)}",
            )
        )
        section_index += 1

        sections.append(
            render_section(
                f"{section_index}. UV volcano ({args.uv_label})",
                render_figure(args.volcano_plot),
            )
        )
        section_index += 1

        lrt_sig_available = (
            lrt_sig_fc_stats is not None
            and args.lrt_sig_fc_plot
            and args.lrt_sig_volcano_plot
        )
        if lrt_sig_available:
            sections.append(
                render_section(
                    f"{section_index}. UV fold change — LRT-significant genes ({args.uv_label})",
                    f"{df_to_html(lrt_sig_fc_stats)}\n{render_figure(args.lrt_sig_fc_plot)}",
                )
            )
            section_index += 1

            sections.append(
                render_section(
                    f"{section_index}. UV volcano — LRT-significant genes ({args.uv_label})",
                    render_figure(args.lrt_sig_volcano_plot),
                )
            )
            section_index += 1
    else:
        sections.append(
            render_section(
                f"{section_index}. UV response",
                "<p>Not applicable for this branch.</p>",
            )
        )
        section_index += 1

    sections.append(
        render_section(
            f"{section_index}. G4 structure class enrichment",
            (
                f"{df_to_html(enrichment.head(20)) if enrichment is not None else '<p>N/A</p>'}\n"
                f"{render_figure(args.structure_plot)}"
            ),
        )
    )

    html = TEMPLATE.format(
        page_title=page_title,
        executive_summary=exec_summary,
        sections="\n".join(sections),
        software_versions=versions,
    )

    Path(args.out_html).write_text(html)
    print(f"Written: {args.out_html}", file=sys.stdout)


if __name__ == "__main__":
    main()
