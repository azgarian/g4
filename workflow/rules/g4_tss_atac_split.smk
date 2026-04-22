"""
ATAC-filtered G4-TSS transcription analysis.

Reuses shared TSS annotation prep from g4_tss, filters merged G4 peaks by
timepoint-specific ATAC boundaries, and reruns the same Tasks 1-8 style
analysis using matched RNA-seq timepoints.
"""

from pathlib import Path


G4_TSS_ATAC_SPLIT_CFG = config["g4_tss_atac_split"]
G4_TSS_ATAC_SPLIT_TIMEPOINTS = [str(tp) for tp in G4_TSS_ATAC_SPLIT_CFG["timepoints"]]
G4_TSS_ATAC_SPLIT_ATAC_BEDS = {
    str(tp): path for tp, path in G4_TSS_ATAC_SPLIT_CFG["atac_peak_beds"].items()
}
G4_TSS_ATAC_SPLIT_RNASEQ_TIMEPOINTS = {
    str(tp): int(rna_tp)
    for tp, rna_tp in G4_TSS_ATAC_SPLIT_CFG["rnaseq_timepoint_map"].items()
}
G4_TSS_ATAC_SPLIT_BIGWIGS = {
    str(tp): path for tp, path in G4_TSS_ATAC_SPLIT_CFG["rnaseq_bigwig_map"].items()
}
G4_TSS_ATAC_SPLIT_PAIRWISE = {
    str(tp): path for tp, path in G4_TSS_ATAC_SPLIT_CFG["deseq_pairwise_map"].items()
}
G4_TSS_ATAC_SPLIT_UV_TIMEPOINTS = [
    tp for tp in G4_TSS_ATAC_SPLIT_TIMEPOINTS if tp in G4_TSS_ATAC_SPLIT_PAIRWISE
]
G4_TSS_ATAC_SPLIT_UV_PATTERN = "|".join(G4_TSS_ATAC_SPLIT_UV_TIMEPOINTS)


def g4_tss_atac_split_atac_bed(wildcards):
    return G4_TSS_ATAC_SPLIT_ATAC_BEDS[wildcards.timepoint]


def g4_tss_atac_split_rna_timepoint(wildcards):
    return str(G4_TSS_ATAC_SPLIT_RNASEQ_TIMEPOINTS[wildcards.timepoint])


def g4_tss_atac_split_bigwig(wildcards):
    return G4_TSS_ATAC_SPLIT_BIGWIGS[wildcards.timepoint]


def g4_tss_atac_split_pairwise(wildcards):
    return G4_TSS_ATAC_SPLIT_PAIRWISE[wildcards.timepoint]


def g4_tss_atac_split_analysis_label(wildcards):
    return wildcards.timepoint


def g4_tss_atac_split_expression_label(wildcards):
    rna_timepoint = G4_TSS_ATAC_SPLIT_RNASEQ_TIMEPOINTS[wildcards.timepoint]
    if rna_timepoint == 0:
        return "Baseline expression"
    return f"{rna_timepoint} min RNA-seq expression"


def g4_tss_atac_split_uv_label(wildcards):
    rna_timepoint = G4_TSS_ATAC_SPLIT_RNASEQ_TIMEPOINTS[wildcards.timepoint]
    return f"0 vs {rna_timepoint} min"


def g4_tss_atac_split_bigwig_label(wildcards):
    return Path(G4_TSS_ATAC_SPLIT_BIGWIGS[wildcards.timepoint]).stem


def g4_tss_atac_split_metaprofile_title(wildcards):
    rna_timepoint = G4_TSS_ATAC_SPLIT_RNASEQ_TIMEPOINTS[wildcards.timepoint]
    if rna_timepoint == 0:
        return "RNA-seq metaprofile centred on TSS (baseline RNA-seq)"
    return f"RNA-seq metaprofile centred on TSS ({rna_timepoint} min RNA-seq)"


rule g4_tss_atac_split_filter_g4:
    input:
        g4="results/g4chip_g4cuttag/g4_hela_chip_cuttag_merged.bed",
        atac=g4_tss_atac_split_atac_bed,
    output:
        bed="results/g4_tss_atac_split/{timepoint}/g4_atac_center_filtered.bed",
        summary="results/g4_tss_atac_split/{timepoint}/g4_atac_center_filter_summary.tsv",
    log:
        "logs/g4_tss_atac_split/{timepoint}/filter_g4.log",
    conda:
        "../envs/g4_tss.yaml"
    params:
        analysis_label=g4_tss_atac_split_analysis_label,
    shell:
        """
        (echo "`date -R`: Filtering merged G4 loci by ATAC center overlap..." &&
        mkdir -p "$(dirname {output.bed})" "$(dirname {log})" &&
        python3 workflow/scripts/g4_tss_filter_g4_by_atac.py \
          --g4-bed {input.g4} \
          --atac-bed {input.atac} \
          --analysis-label "{params.analysis_label}" \
          --out-bed {output.bed} \
          --out-summary {output.summary} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule g4_tss_atac_split_intersect_g4:
    input:
        windows=rules.g4_tss_slop_windows.output.windows,
        peaks=rules.g4_tss_atac_split_filter_g4.output.bed,
    output:
        "results/g4_tss_atac_split/{timepoint}/windows_g4_intersect.bed",
    log:
        "logs/g4_tss_atac_split/{timepoint}/intersect_g4.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Intersecting TSS windows with ATAC-filtered G4 peaks (-c)..." &&
        bedtools intersect \
          -nonamecheck \
          -c \
          -a {input.windows} \
          -b {input.peaks} \
          > {output} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule g4_tss_atac_split_build_annotation:
    input:
        gene_table=rules.g4_tss_extract_canonical_tss.output.gene_table,
        g4_intersect=rules.g4_tss_atac_split_intersect_g4.output,
        gc_bg_intersect=rules.g4_tss_intersect_gc_bg.output,
    output:
        tsv="results/g4_tss_atac_split/{timepoint}/tss_group_annotation.tsv",
    log:
        "logs/g4_tss_atac_split/{timepoint}/build_annotation.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Building ATAC-split TSS group annotation table..." &&
        python3 workflow/scripts/g4_tss_build_annotation.py \
          --gene-table {input.gene_table} \
          --g4-intersect {input.g4_intersect} \
          --gc-bg-intersect {input.gc_bg_intersect} \
          --out-tsv {output.tsv} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1
        """


rule g4_tss_atac_split_baseline_expression:
    input:
        annotation=rules.g4_tss_atac_split_build_annotation.output.tsv,
        tpm="results/rnaseq/matrices/gene_tpm.tsv.gz",
        norm_counts="results/rnaseq/matrices/normalized_counts.tsv.gz",
        sample_manifest="results/rnaseq/metadata/sample_manifest.tsv",
    output:
        baseline_tpm="results/g4_tss_atac_split/{timepoint}/baseline_tpm.tsv",
        tpm_summary="results/g4_tss_atac_split/{timepoint}/baseline_tpm_summary.tsv",
        norm_counts="results/g4_tss_atac_split/{timepoint}/baseline_normalized_counts.tsv",
        by_group="results/g4_tss_atac_split/{timepoint}/gene_expression_by_group.tsv",
    log:
        "logs/g4_tss_atac_split/{timepoint}/baseline_expression.log",
    conda:
        "../envs/g4_tss.yaml"
    params:
        rna_timepoint=g4_tss_atac_split_rna_timepoint,
        analysis_label=g4_tss_atac_split_analysis_label,
    shell:
        """
        (echo "`date -R`: Compiling matched-time expression summaries..." &&
        python3 workflow/scripts/g4_tss_baseline_expression.py \
          --tss-annotation {input.annotation} \
          --tpm-matrix {input.tpm} \
          --norm-counts-matrix {input.norm_counts} \
          --sample-manifest {input.sample_manifest} \
          --rna-timepoint {params.rna_timepoint} \
          --analysis-label "{params.analysis_label}" \
          --out-baseline-tpm {output.baseline_tpm} \
          --out-baseline-tpm-summary {output.tpm_summary} \
          --out-baseline-norm-counts {output.norm_counts} \
          --out-gene-expression-by-group {output.by_group} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule g4_tss_atac_split_expression_stats:
    input:
        by_group=rules.g4_tss_atac_split_baseline_expression.output.by_group,
    output:
        statistics="results/g4_tss_atac_split/{timepoint}/expression_group_statistics.tsv",
        fisher="results/g4_tss_atac_split/{timepoint}/expression_group_fisher.tsv",
        summary="results/g4_tss_atac_split/{timepoint}/expression_group_summary.tsv",
    log:
        "logs/g4_tss_atac_split/{timepoint}/expression_stats.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Running matched-time expression group statistics..." &&
        python3 workflow/scripts/g4_tss_expression_stats.py \
          --gene-expression-by-group {input.by_group} \
          --out-statistics {output.statistics} \
          --out-fisher {output.fisher} \
          --out-summary {output.summary} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule g4_tss_atac_split_visualize:
    input:
        by_group=rules.g4_tss_atac_split_baseline_expression.output.by_group,
        annotation=rules.g4_tss_atac_split_build_annotation.output.tsv,
        statistics=rules.g4_tss_atac_split_expression_stats.output.statistics,
    output:
        violin="results/g4_tss_atac_split/{timepoint}/expression_violin_by_group.pdf",
        bed_g4="results/g4_tss_atac_split/{timepoint}/tss_G4_TSS.bed",
        bed_gc="results/g4_tss_atac_split/{timepoint}/tss_GC_bg_TSS.bed",
        bed_no="results/g4_tss_atac_split/{timepoint}/tss_No_overlap.bed",
    log:
        "logs/g4_tss_atac_split/{timepoint}/visualize.log",
    conda:
        "../envs/g4_tss.yaml"
    params:
        analysis_label=g4_tss_atac_split_analysis_label,
        expression_label=g4_tss_atac_split_expression_label,
    shell:
        """
        (echo "`date -R`: Generating matched-time violin plot and group BED files..." &&
        python3 workflow/scripts/g4_tss_visualize.py \
          --gene-expression-by-group {input.by_group} \
          --tss-annotation {input.annotation} \
          --expression-statistics {input.statistics} \
          --analysis-label "{params.analysis_label}" \
          --expression-label "{params.expression_label}" \
          --out-violin-pdf {output.violin} \
          --out-bed-g4-tss {output.bed_g4} \
          --out-bed-gc-bg-tss {output.bed_gc} \
          --out-bed-no-overlap {output.bed_no} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule g4_tss_atac_split_compute_matrix:
    input:
        bw=g4_tss_atac_split_bigwig,
        bed_g4=rules.g4_tss_atac_split_visualize.output.bed_g4,
        bed_gc=rules.g4_tss_atac_split_visualize.output.bed_gc,
        bed_no=rules.g4_tss_atac_split_visualize.output.bed_no,
    output:
        matrix="results/g4_tss_atac_split/{timepoint}/rnaseq_tss_matrix.gz",
    log:
        "logs/g4_tss_atac_split/{timepoint}/compute_matrix.log",
    threads: 8
    conda:
        "../envs/deeptools.yaml"
    params:
        sample_label=g4_tss_atac_split_bigwig_label,
    shell:
        """
        (echo "`date -R`: Running deepTools computeMatrix for matched RNA-seq signal..." &&
        computeMatrix reference-point \
          --referencePoint TSS \
          --scoreFileName {input.bw} \
          --regionsFileName {input.bed_g4} {input.bed_gc} {input.bed_no} \
          --beforeRegionStartLength 2000 \
          --afterRegionStartLength 2000 \
          --binSize 50 \
          --missingDataAsZero \
          --skipZeros \
          --samplesLabel {params.sample_label} \
          --numberOfProcessors {threads} \
          --outFileName {output.matrix} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule g4_tss_atac_split_plot_profile:
    input:
        matrix=rules.g4_tss_atac_split_compute_matrix.output.matrix,
    output:
        profile="results/g4_tss_atac_split/{timepoint}/rnaseq_tss_metaprofile.pdf",
    log:
        "logs/g4_tss_atac_split/{timepoint}/plot_profile.log",
    conda:
        "../envs/deeptools.yaml"
    params:
        plot_title=g4_tss_atac_split_metaprofile_title,
    shell:
        """
        (echo "`date -R`: Running deepTools plotProfile..." &&
        plotProfile \
          --matrixFile {input.matrix} \
          --outFileName {output.profile} \
          --regionsLabel "G4_TSS" "GC_bg_TSS" "No_overlap" \
          --plotTitle "{params.plot_title}" \
          --yAxisLabel "Mean RPKM" \
          --refPointLabel "TSS" \
          --plotType lines \
          --perGroup &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule g4_tss_atac_split_decile_analysis:
    input:
        by_group=rules.g4_tss_atac_split_baseline_expression.output.by_group,
        windows=rules.g4_tss_slop_windows.output.windows,
        g4_merged=rules.g4_tss_atac_split_filter_g4.output.bed,
        gc_bg="results/gc_rich_bg/gc_rich_bg_prepared.bed",
    output:
        deciles="results/g4_tss_atac_split/{timepoint}/expression_deciles.tsv",
        overlap_fractions="results/g4_tss_atac_split/{timepoint}/decile_overlap_fractions.tsv",
        correlation_stats="results/g4_tss_atac_split/{timepoint}/decile_correlation_stats.tsv",
        enrichment_plot="results/g4_tss_atac_split/{timepoint}/decile_enrichment_plot.pdf",
    log:
        "logs/g4_tss_atac_split/{timepoint}/decile_analysis.log",
    conda:
        "../envs/g4_tss.yaml"
    params:
        analysis_label=g4_tss_atac_split_analysis_label,
    shell:
        """
        (echo "`date -R`: Running matched-time expression decile analysis..." &&
        python3 workflow/scripts/g4_tss_decile_analysis.py \
          --gene-expression-by-group {input.by_group} \
          --tss-windows-bed {input.windows} \
          --g4-merged-bed {input.g4_merged} \
          --gc-bg-bed {input.gc_bg} \
          --analysis-label "{params.analysis_label}" \
          --g4-label "ATAC-filtered G4" \
          --out-deciles {output.deciles} \
          --out-overlap-fractions {output.overlap_fractions} \
          --out-correlation-stats {output.correlation_stats} \
          --out-enrichment-plot {output.enrichment_plot} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule g4_tss_atac_split_uv_response:
    wildcard_constraints:
        timepoint=G4_TSS_ATAC_SPLIT_UV_PATTERN,
    input:
        annotation=rules.g4_tss_atac_split_build_annotation.output.tsv,
        lrt="results/rnaseq/deseq2/time_lrt_results.tsv.gz",
        pairwise=g4_tss_atac_split_pairwise,
    output:
        lrt_summary="results/g4_tss_atac_split/{timepoint}/uv_group_lrt_summary.tsv",
        fc_stats="results/g4_tss_atac_split/{timepoint}/uv_group_fold_change_stats.tsv",
        fc_plot="results/g4_tss_atac_split/{timepoint}/uv_fold_change_by_group.pdf",
        volcano_plot="results/g4_tss_atac_split/{timepoint}/uv_volcano_by_group.pdf",
    log:
        "logs/g4_tss_atac_split/{timepoint}/uv_response.log",
    conda:
        "../envs/g4_tss.yaml"
    params:
        analysis_label=g4_tss_atac_split_analysis_label,
        contrast_label=g4_tss_atac_split_uv_label,
    shell:
        """
        (echo "`date -R`: Analysing UV-induced response by ATAC-split promoter group..." &&
        python3 workflow/scripts/g4_tss_uv_response.py \
          --tss-annotation {input.annotation} \
          --lrt-results {input.lrt} \
          --pairwise-results {input.pairwise} \
          --analysis-label "{params.analysis_label}" \
          --contrast-label "{params.contrast_label}" \
          --out-lrt-summary {output.lrt_summary} \
          --out-fc-stats {output.fc_stats} \
          --out-fc-plot {output.fc_plot} \
          --out-volcano-plot {output.volcano_plot} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule g4_tss_atac_split_structure_enrichment:
    input:
        annotation=rules.g4_tss_atac_split_build_annotation.output.tsv,
        windows=rules.g4_tss_slop_windows.output.windows,
        baseline_tpm=rules.g4_tss_atac_split_baseline_expression.output.baseline_tpm,
        g4chip_tsv="results/g4chip/g4_hela_peaks_prepared.tsv",
        g4cuttag_tsv="results/g4cuttag/g4_hela_peaks_prepared.tsv",
    output:
        struct_by_expr="results/g4_tss_atac_split/{timepoint}/g4_structure_by_expression_class.tsv",
        struct_plot="results/g4_tss_atac_split/{timepoint}/structure_class_by_expression_class.pdf",
        enrichment_stats="results/g4_tss_atac_split/{timepoint}/structure_class_enrichment_stats.tsv",
    log:
        "logs/g4_tss_atac_split/{timepoint}/structure_enrichment.log",
    conda:
        "../envs/g4_tss.yaml"
    params:
        analysis_label=g4_tss_atac_split_analysis_label,
    shell:
        """
        (echo "`date -R`: Analysing G4 structure class enrichment..." &&
        python3 workflow/scripts/g4_tss_structure_enrichment.py \
          --tss-annotation {input.annotation} \
          --tss-windows-bed {input.windows} \
          --baseline-tpm {input.baseline_tpm} \
          --g4chip-tsv {input.g4chip_tsv} \
          --g4cuttag-tsv {input.g4cuttag_tsv} \
          --analysis-label "{params.analysis_label}" \
          --out-structure-by-expression {output.struct_by_expr} \
          --out-plot {output.struct_plot} \
          --out-enrichment-stats {output.enrichment_stats} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule g4_tss_atac_split_report:
    wildcard_constraints:
        timepoint=G4_TSS_ATAC_SPLIT_UV_PATTERN,
    input:
        filter_summary=rules.g4_tss_atac_split_filter_g4.output.summary,
        annotation=rules.g4_tss_atac_split_build_annotation.output.tsv,
        expr_summary=rules.g4_tss_atac_split_expression_stats.output.summary,
        expr_statistics=rules.g4_tss_atac_split_expression_stats.output.statistics,
        expr_fisher=rules.g4_tss_atac_split_expression_stats.output.fisher,
        violin=rules.g4_tss_atac_split_visualize.output.violin,
        metaprofile=rules.g4_tss_atac_split_plot_profile.output.profile,
        decile_overlap=rules.g4_tss_atac_split_decile_analysis.output.overlap_fractions,
        decile_corr=rules.g4_tss_atac_split_decile_analysis.output.correlation_stats,
        decile_plot=rules.g4_tss_atac_split_decile_analysis.output.enrichment_plot,
        lrt_summary=rules.g4_tss_atac_split_uv_response.output.lrt_summary,
        fc_stats=rules.g4_tss_atac_split_uv_response.output.fc_stats,
        fc_plot=rules.g4_tss_atac_split_uv_response.output.fc_plot,
        volcano=rules.g4_tss_atac_split_uv_response.output.volcano_plot,
        enrichment_stats=rules.g4_tss_atac_split_structure_enrichment.output.enrichment_stats,
        struct_plot=rules.g4_tss_atac_split_structure_enrichment.output.struct_plot,
    output:
        html="results/g4_tss_atac_split/{timepoint}/g4_tss_atac_split_report.html",
        versions="results/g4_tss_atac_split/{timepoint}/software_versions.txt",
    log:
        "logs/g4_tss_atac_split/{timepoint}/report.log",
    conda:
        "../envs/g4_tss.yaml"
    params:
        analysis_label=g4_tss_atac_split_analysis_label,
        expression_label=g4_tss_atac_split_expression_label,
        uv_label=g4_tss_atac_split_uv_label,
    shell:
        """
        (echo "`date -R`: Generating ATAC-split HTML report..." &&
        python3 workflow/scripts/g4_tss_report.py \
          --tss-annotation {input.annotation} \
          --expr-summary {input.expr_summary} \
          --expr-statistics {input.expr_statistics} \
          --expr-fisher {input.expr_fisher} \
          --violin-plot {input.violin} \
          --metaprofile-plot {input.metaprofile} \
          --decile-overlap {input.decile_overlap} \
          --decile-corr {input.decile_corr} \
          --decile-plot {input.decile_plot} \
          --lrt-summary {input.lrt_summary} \
          --fc-stats {input.fc_stats} \
          --fc-plot {input.fc_plot} \
          --volcano-plot {input.volcano} \
          --enrichment-stats {input.enrichment_stats} \
          --structure-plot {input.struct_plot} \
          --analysis-label "{params.analysis_label}" \
          --expression-label "{params.expression_label}" \
          --filter-summary {input.filter_summary} \
          --uv-label "{params.uv_label}" \
          --out-html {output.html} \
          --out-versions {output.versions} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule g4_tss_atac_split_report_nouv:
    wildcard_constraints:
        timepoint="noUV",
    input:
        filter_summary=rules.g4_tss_atac_split_filter_g4.output.summary,
        annotation=rules.g4_tss_atac_split_build_annotation.output.tsv,
        expr_summary=rules.g4_tss_atac_split_expression_stats.output.summary,
        expr_statistics=rules.g4_tss_atac_split_expression_stats.output.statistics,
        expr_fisher=rules.g4_tss_atac_split_expression_stats.output.fisher,
        violin=rules.g4_tss_atac_split_visualize.output.violin,
        metaprofile=rules.g4_tss_atac_split_plot_profile.output.profile,
        decile_overlap=rules.g4_tss_atac_split_decile_analysis.output.overlap_fractions,
        decile_corr=rules.g4_tss_atac_split_decile_analysis.output.correlation_stats,
        decile_plot=rules.g4_tss_atac_split_decile_analysis.output.enrichment_plot,
        enrichment_stats=rules.g4_tss_atac_split_structure_enrichment.output.enrichment_stats,
        struct_plot=rules.g4_tss_atac_split_structure_enrichment.output.struct_plot,
    output:
        html="results/g4_tss_atac_split/{timepoint}/g4_tss_atac_split_report.html",
        versions="results/g4_tss_atac_split/{timepoint}/software_versions.txt",
    log:
        "logs/g4_tss_atac_split/{timepoint}/report.log",
    conda:
        "../envs/g4_tss.yaml"
    params:
        analysis_label=g4_tss_atac_split_analysis_label,
        expression_label=g4_tss_atac_split_expression_label,
    shell:
        """
        (echo "`date -R`: Generating ATAC-split HTML report..." &&
        python3 workflow/scripts/g4_tss_report.py \
          --tss-annotation {input.annotation} \
          --expr-summary {input.expr_summary} \
          --expr-statistics {input.expr_statistics} \
          --expr-fisher {input.expr_fisher} \
          --violin-plot {input.violin} \
          --metaprofile-plot {input.metaprofile} \
          --decile-overlap {input.decile_overlap} \
          --decile-corr {input.decile_corr} \
          --decile-plot {input.decile_plot} \
          --enrichment-stats {input.enrichment_stats} \
          --structure-plot {input.struct_plot} \
          --analysis-label "{params.analysis_label}" \
          --expression-label "{params.expression_label}" \
          --filter-summary {input.filter_summary} \
          --out-html {output.html} \
          --out-versions {output.versions} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """
