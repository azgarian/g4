"""
G4-TSS transcription analysis module — Tasks 1-8.

Answers the question: Are G4s at TSSs linked to active transcription?

Prerequisites (must already exist):
  results/g4chip_g4cuttag/g4_hela_chip_cuttag_merged.bed
  results/gc_rich_bg_promoter/gc_rich_bg_promoter_sampled.bed
  results/g4chip/g4_hela_peaks_prepared.tsv
  results/g4cuttag/g4_hela_peaks_prepared.tsv
  resources/ref_genomes/hg38_rnaseq_v35/gencode.v35.annotation.gtf
  resources/ref_genomes/hg38/genome_hg38.fa.fai
  results/rnaseq/matrices/gene_tpm.tsv.gz
  results/rnaseq/matrices/normalized_counts.tsv.gz
  results/rnaseq/deseq2/time_lrt_results.tsv.gz
  results/rnaseq/deseq2/pairwise/0_vs_60_results.tsv.gz
  resources/rna-seq/merged_t00.bw
"""

# ─── Task 1a: Extract canonical TSS from GENCODE v35 GTF ──────────────────────

rule g4_tss_extract_canonical_tss:
    input:
        gtf="resources/ref_genomes/hg38_rnaseq_v35/gencode.v35.annotation.gtf",
    output:
        bed="results/g4_tss/canonical_tss.bed",
        gene_table="results/g4_tss/gene_name_table.tsv",
    log:
        "logs/g4_tss/annotation.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Extracting canonical TSS from GENCODE v35 GTF..." &&
        mkdir -p results/g4_tss logs/g4_tss &&
        python3 workflow/scripts/g4_tss_canonical_tss.py \
          --gtf {input.gtf} \
          --out-bed {output.bed} \
          --out-gene-table {output.gene_table} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1
        """


# ─── Task 1b: Slop TSS to ±500 bp promoter windows ────────────────────────────

rule g4_tss_slop_windows:
    input:
        bed=rules.g4_tss_extract_canonical_tss.output.bed,
        fai="resources/ref_genomes/hg38/genome_hg38.fa.fai",
    output:
        windows="results/g4_tss/canonical_tss_windows_1kb.bed",
    log:
        "logs/g4_tss/slop_windows.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Expanding TSS to ±500 bp windows..." &&
        bedtools slop \
          -i {input.bed} \
          -g {input.fai} \
          -b 500 \
          > {output.windows} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 1c: Intersect promoter windows with G4 merged peaks ─────────────────

rule g4_tss_intersect_g4:
    input:
        windows=rules.g4_tss_slop_windows.output.windows,
        peaks="results/g4chip_g4cuttag/g4_hela_chip_cuttag_merged.bed",
    output:
        "results/g4_tss/windows_g4_intersect.bed",
    log:
        "logs/g4_tss/intersect_g4.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Intersecting TSS windows with G4 merged peaks (-c)..." &&
        bedtools intersect \
          -c \
          -a {input.windows} \
          -b {input.peaks} \
          > {output} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 1d: Build final annotation table with group assignments ──────────────

rule g4_tss_build_annotation:
    input:
        gene_table=rules.g4_tss_extract_canonical_tss.output.gene_table,
        g4_intersect=rules.g4_tss_intersect_g4.output,
        gc_bg_bed="results/gc_rich_bg_promoter/gc_rich_bg_promoter_sampled.bed",
    output:
        tsv="results/g4_tss/tss_group_annotation.tsv",
    log:
        "logs/g4_tss/build_annotation.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Building TSS group annotation table..." &&
        python3 workflow/scripts/g4_tss_build_annotation.py \
          --gene-table {input.gene_table} \
          --g4-intersect {input.g4_intersect} \
          --gc-bg-bed {input.gc_bg_bed} \
          --out-tsv {output.tsv} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1
        """


# ─── Task 2: Baseline expression summaries ────────────────────────────────────

rule g4_tss_baseline_expression:
    input:
        annotation=rules.g4_tss_build_annotation.output.tsv,
        tpm="results/rnaseq/matrices/gene_tpm.tsv.gz",
        norm_counts="results/rnaseq/matrices/normalized_counts.tsv.gz",
    output:
        baseline_tpm="results/g4_tss/baseline_tpm.tsv",
        tpm_summary="results/g4_tss/baseline_tpm_summary.tsv",
        norm_counts="results/g4_tss/baseline_normalized_counts.tsv",
        by_group="results/g4_tss/gene_expression_by_group.tsv",
    log:
        "logs/g4_tss/baseline_expression.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Compiling baseline expression summaries..." &&
        python3 workflow/scripts/g4_tss_baseline_expression.py \
          --tss-annotation {input.annotation} \
          --tpm-matrix {input.tpm} \
          --norm-counts-matrix {input.norm_counts} \
          --out-baseline-tpm {output.baseline_tpm} \
          --out-baseline-tpm-summary {output.tpm_summary} \
          --out-baseline-norm-counts {output.norm_counts} \
          --out-gene-expression-by-group {output.by_group} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 2b: Retained vs removed genes after RNA-seq prefilter ──────────────

rule g4_tss_group_retention:
    input:
        annotation=rules.g4_tss_build_annotation.output.tsv,
        by_group=rules.g4_tss_baseline_expression.output.by_group,
    output:
        summary="results/g4_tss/group_gene_retention_summary.tsv",
        barplot="results/g4_tss/group_gene_retention_barplot.pdf",
    log:
        "logs/g4_tss/group_retention.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Summarizing retained vs removed genes by TSS group..." &&
        python3 workflow/scripts/g4_tss_group_retention.py \
          --tss-annotation {input.annotation} \
          --gene-expression-by-group {input.by_group} \
          --out-summary-tsv {output.summary} \
          --out-barplot-pdf {output.barplot} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 3: Compare baseline expression across groups ────────────────────────

rule g4_tss_expression_stats:
    input:
        by_group=rules.g4_tss_baseline_expression.output.by_group,
    output:
        statistics="results/g4_tss/expression_group_statistics.tsv",
        fisher="results/g4_tss/expression_group_fisher.tsv",
        summary="results/g4_tss/expression_group_summary.tsv",
    log:
        "logs/g4_tss/expression_stats.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Running expression group statistics..." &&
        python3 workflow/scripts/g4_tss_expression_stats.py \
          --gene-expression-by-group {input.by_group} \
          --out-statistics {output.statistics} \
          --out-fisher {output.fisher} \
          --out-summary {output.summary} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 4a: Violin plot and group TSS BED files ─────────────────────────────

rule g4_tss_visualize:
    input:
        by_group=rules.g4_tss_baseline_expression.output.by_group,
        annotation=rules.g4_tss_build_annotation.output.tsv,
        statistics=rules.g4_tss_expression_stats.output.statistics,
    output:
        violin="results/g4_tss/expression_violin_by_group.pdf",
        bed_g4="results/g4_tss/tss_G4_TSS.bed",
        bed_gc="results/g4_tss/tss_GC_bg_TSS.bed",
        bed_no="results/g4_tss/tss_No_overlap.bed",
    log:
        "logs/g4_tss/visualize.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Generating violin plot and group BED files..." &&
        python3 workflow/scripts/g4_tss_visualize.py \
          --gene-expression-by-group {input.by_group} \
          --tss-annotation {input.annotation} \
          --expression-statistics {input.statistics} \
          --out-violin-pdf {output.violin} \
          --out-bed-g4-tss {output.bed_g4} \
          --out-bed-gc-bg-tss {output.bed_gc} \
          --out-bed-no-overlap {output.bed_no} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 4b: deepTools computeMatrix ─────────────────────────────────────────

rule g4_tss_compute_matrix:
    input:
        bw="resources/rna-seq/merged_t00.bw",
        bed_g4=rules.g4_tss_visualize.output.bed_g4,
        bed_gc=rules.g4_tss_visualize.output.bed_gc,
        bed_no=rules.g4_tss_visualize.output.bed_no,
    output:
        matrix="results/g4_tss/rnaseq_tss_matrix.gz",
    log:
        "logs/g4_tss/compute_matrix.log",
    threads: 8
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        (echo "`date -R`: Running deepTools computeMatrix..." &&
        computeMatrix reference-point \
          --referencePoint TSS \
          --scoreFileName {input.bw} \
          --regionsFileName {input.bed_g4} {input.bed_gc} {input.bed_no} \
          --beforeRegionStartLength 2000 \
          --afterRegionStartLength 2000 \
          --binSize 50 \
          --missingDataAsZero \
          --skipZeros \
          --samplesLabel merged_t00 \
          --numberOfProcessors {threads} \
          --outFileName {output.matrix} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 4c: deepTools plotProfile ───────────────────────────────────────────

rule g4_tss_plot_profile:
    input:
        matrix=rules.g4_tss_compute_matrix.output.matrix,
    output:
        profile="results/g4_tss/rnaseq_tss_metaprofile.pdf",
    log:
        "logs/g4_tss/plot_profile.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        (echo "`date -R`: Running deepTools plotProfile..." &&
        plotProfile \
          --matrixFile {input.matrix} \
          --outFileName {output.profile} \
          --regionsLabel "G4_TSS" "GC_bg_TSS" "No_overlap" \
          --plotTitle "RNA-seq metaprofile centred on TSS" \
          --yAxisLabel "Mean RPKM" \
          --refPointLabel "TSS" \
          --plotType lines \
          --colors "#d62728" "#ff7f0e" "#1f77b4" &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 5: Decile analysis ───────────────────────────────────────────────────

rule g4_tss_decile_analysis:
    input:
        by_group=rules.g4_tss_baseline_expression.output.by_group,
        windows=rules.g4_tss_slop_windows.output.windows,
        g4_merged="results/g4chip_g4cuttag/g4_hela_chip_cuttag_merged.bed",
        gc_bg="results/gc_rich_bg_promoter/gc_rich_bg_promoter_sampled.bed",
    output:
        deciles="results/g4_tss/expression_deciles.tsv",
        overlap_fractions="results/g4_tss/decile_overlap_fractions.tsv",
        correlation_stats="results/g4_tss/decile_correlation_stats.tsv",
        enrichment_plot="results/g4_tss/decile_enrichment_plot.pdf",
    log:
        "logs/g4_tss/decile_analysis.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Running expression decile analysis..." &&
        python3 workflow/scripts/g4_tss_decile_analysis.py \
          --gene-expression-by-group {input.by_group} \
          --tss-windows-bed {input.windows} \
          --g4-merged-bed {input.g4_merged} \
          --gc-bg-bed {input.gc_bg} \
          --out-deciles {output.deciles} \
          --out-overlap-fractions {output.overlap_fractions} \
          --out-correlation-stats {output.correlation_stats} \
          --out-enrichment-plot {output.enrichment_plot} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 6: UV-induced transcriptional response ──────────────────────────────

rule g4_tss_uv_response:
    input:
        annotation=rules.g4_tss_build_annotation.output.tsv,
        lrt="results/rnaseq/deseq2/time_lrt_results.tsv.gz",
        pairwise_0_60="results/rnaseq/deseq2/pairwise/0_vs_60_results.tsv.gz",
    output:
        lrt_summary="results/g4_tss/uv_group_lrt_summary.tsv",
        fc_stats="results/g4_tss/uv_group_fold_change_stats.tsv",
        fc_plot="results/g4_tss/uv_fold_change_by_group.pdf",
        volcano_plot="results/g4_tss/uv_volcano_by_group.pdf",
        lrt_sig_fc_stats="results/g4_tss/uv_group_fold_change_stats_lrt_sig.tsv",
        lrt_sig_fc_plot="results/g4_tss/uv_fold_change_by_group_lrt_sig.pdf",
        lrt_sig_volcano_plot="results/g4_tss/uv_volcano_by_group_lrt_sig.pdf",
    log:
        "logs/g4_tss/uv_response.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Analysing UV-induced response by TSS group..." &&
        python3 workflow/scripts/g4_tss_uv_response.py \
          --tss-annotation {input.annotation} \
          --lrt-results {input.lrt} \
          --pairwise-0-vs-60 {input.pairwise_0_60} \
          --out-lrt-summary {output.lrt_summary} \
          --out-fc-stats {output.fc_stats} \
          --out-fc-plot {output.fc_plot} \
          --out-volcano-plot {output.volcano_plot} \
          --out-lrt-sig-fc-stats {output.lrt_sig_fc_stats} \
          --out-lrt-sig-fc-plot {output.lrt_sig_fc_plot} \
          --out-lrt-sig-volcano-plot {output.lrt_sig_volcano_plot} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 7: G4 structure class enrichment across expression classes ───────────

rule g4_tss_structure_enrichment:
    input:
        annotation=rules.g4_tss_build_annotation.output.tsv,
        windows=rules.g4_tss_slop_windows.output.windows,
        baseline_tpm=rules.g4_tss_baseline_expression.output.baseline_tpm,
        g4chip_tsv="results/g4chip/g4_hela_peaks_prepared.tsv",
        g4cuttag_tsv="results/g4cuttag/g4_hela_peaks_prepared.tsv",
    output:
        struct_by_expr="results/g4_tss/g4_structure_by_expression_class.tsv",
        struct_plot="results/g4_tss/structure_class_by_expression_class.pdf",
        enrichment_stats="results/g4_tss/structure_class_enrichment_stats.tsv",
    log:
        "logs/g4_tss/structure_enrichment.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Analysing G4 structure class enrichment..." &&
        python3 workflow/scripts/g4_tss_structure_enrichment.py \
          --tss-annotation {input.annotation} \
          --tss-windows-bed {input.windows} \
          --baseline-tpm {input.baseline_tpm} \
          --g4chip-tsv {input.g4chip_tsv} \
          --g4cuttag-tsv {input.g4cuttag_tsv} \
          --out-structure-by-expression {output.struct_by_expr} \
          --out-plot {output.struct_plot} \
          --out-enrichment-stats {output.enrichment_stats} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 9: Per-gene promoter-G4 × LRT intersection table ───────────────────

rule g4_tss_pg4_deg_intersect:
    input:
        annotation=rules.g4_tss_build_annotation.output.tsv,
        lrt="results/rnaseq/deseq2/time_lrt_results.tsv.gz",
        norm_counts="results/rnaseq/matrices/normalized_counts.tsv.gz",
        baseline_tpm=rules.g4_tss_baseline_expression.output.baseline_tpm,
        pairwise_0_12="results/rnaseq/deseq2/pairwise/0_vs_12_results.tsv.gz",
        pairwise_0_30="results/rnaseq/deseq2/pairwise/0_vs_30_results.tsv.gz",
        pairwise_0_60="results/rnaseq/deseq2/pairwise/0_vs_60_results.tsv.gz",
    output:
        intersect="results/g4_tss/pG4_DEG_intersect.tsv",
    log:
        "logs/g4_tss/pg4_deg_intersect.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Building promoter-G4 × LRT intersection table..." &&
        python3 workflow/scripts/g4_tss_pg4_deg_intersect.py \
          --tss-annotation {input.annotation} \
          --lrt-results {input.lrt} \
          --norm-counts {input.norm_counts} \
          --baseline-tpm {input.baseline_tpm} \
          --pairwise-0-vs-12 {input.pairwise_0_12} \
          --pairwise-0-vs-30 {input.pairwise_0_30} \
          --pairwise-0-vs-60 {input.pairwise_0_60} \
          --out-intersect {output.intersect} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 10: Functional enrichment and G4-strength stratification ─────────────

rule g4_tss_pg4_enrichment:
    input:
        pg4_deg=rules.g4_tss_pg4_deg_intersect.output.intersect,
        windows=rules.g4_tss_slop_windows.output.windows,
        g4_intersect=rules.g4_tss_intersect_g4.output,
        g4chip_source_bed="results/g4chip/g4_hela_peaks_prepared.bed",
        g4cuttag_source_bed="results/g4cuttag/g4_hela_peaks_prepared.bed",
    params:
        gene_set_gmt=lambda w: (
            config.get("g4_tss_pg4", {}).get("gene_set_gmt", "")
        ),
        gene_set_manifest=lambda w: (
            config.get("rna_seq_g4_context", {}).get("geneset_manifest", "")
        ),
        strength_metric=lambda w: (
            config.get("g4_tss_pg4", {}).get("strength_metric", "max_signal")
        ),
    output:
        enrichment_table="results/g4_tss/pG4_pathway_enrichment.tsv",
        g4_strength_table="results/g4_tss/pG4_strength_stratification.tsv",
        g4_strength_plot="results/g4_tss/pG4_strength_stratification.pdf",
        g4_strength_direction_plot="results/g4_tss/pG4_strength_stratification_direction.pdf",
        summary="results/g4_tss/pG4_enrichment_summary.tsv",
    log:
        "logs/g4_tss/pg4_enrichment.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Running G4-strength stratification and pathway enrichment..." &&
        GENE_SET_GMT_ARG="" &&
        GENE_SET_MANIFEST_ARG="" &&
        [ -n "{params.gene_set_gmt}" ] && [ -f "{params.gene_set_gmt}" ] && \
          GENE_SET_GMT_ARG="--gene-set-gmt {params.gene_set_gmt}" || true &&
        [ -n "{params.gene_set_manifest}" ] && [ -f "{params.gene_set_manifest}" ] && \
          GENE_SET_MANIFEST_ARG="--gene-set-manifest {params.gene_set_manifest}" || true &&
        python3 workflow/scripts/g4_tss_pg4_enrichment.py \
          --pg4-deg-intersect {input.pg4_deg} \
          --tss-windows-bed {input.windows} \
          --g4-intersect-bed {input.g4_intersect} \
          --g4chip-source-bed {input.g4chip_source_bed} \
          --g4cuttag-source-bed {input.g4cuttag_source_bed} \
          --strength-metric {params.strength_metric} \
          $GENE_SET_GMT_ARG \
          $GENE_SET_MANIFEST_ARG \
          --out-enrichment-table {output.enrichment_table} \
          --out-g4-strength-table {output.g4_strength_table} \
          --out-g4-strength-plot {output.g4_strength_plot} \
          --out-g4-strength-direction-plot {output.g4_strength_direction_plot} \
          --out-summary {output.summary} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 11: Endpoint QC and sensitivity checks ───────────────────────────────

rule g4_tss_pg4_sensitivity:
    input:
        annotation=rules.g4_tss_build_annotation.output.tsv,
        lrt="results/rnaseq/deseq2/time_lrt_results.tsv.gz",
    output:
        sensitivity="results/g4_tss/pG4_sensitivity_qc.tsv",
    log:
        "logs/g4_tss/pg4_sensitivity.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Running endpoint QC and sensitivity checks..." &&
        python3 workflow/scripts/g4_tss_pg4_sensitivity.py \
          --tss-annotation {input.annotation} \
          --lrt-results {input.lrt} \
          --out-sensitivity-table {output.sensitivity} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 8: HTML report ───────────────────────────────────────────────────────

rule g4_tss_report:
    input:
        annotation=rules.g4_tss_build_annotation.output.tsv,
        expr_summary=rules.g4_tss_expression_stats.output.summary,
        expr_statistics=rules.g4_tss_expression_stats.output.statistics,
        expr_fisher=rules.g4_tss_expression_stats.output.fisher,
        violin=rules.g4_tss_visualize.output.violin,
        metaprofile=rules.g4_tss_plot_profile.output.profile,
        decile_overlap=rules.g4_tss_decile_analysis.output.overlap_fractions,
        decile_corr=rules.g4_tss_decile_analysis.output.correlation_stats,
        decile_plot=rules.g4_tss_decile_analysis.output.enrichment_plot,
        lrt_summary=rules.g4_tss_uv_response.output.lrt_summary,
        fc_stats=rules.g4_tss_uv_response.output.fc_stats,
        fc_plot=rules.g4_tss_uv_response.output.fc_plot,
        volcano=rules.g4_tss_uv_response.output.volcano_plot,
        lrt_sig_fc_stats=rules.g4_tss_uv_response.output.lrt_sig_fc_stats,
        lrt_sig_fc_plot=rules.g4_tss_uv_response.output.lrt_sig_fc_plot,
        lrt_sig_volcano=rules.g4_tss_uv_response.output.lrt_sig_volcano_plot,
        enrichment_stats=rules.g4_tss_structure_enrichment.output.enrichment_stats,
        struct_plot=rules.g4_tss_structure_enrichment.output.struct_plot,
    output:
        html="results/g4_tss/g4_tss_transcription_report.html",
        versions="results/g4_tss/software_versions.txt",
    log:
        "logs/g4_tss/report.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Generating HTML report..." &&
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
          --lrt-sig-fc-stats {input.lrt_sig_fc_stats} \
          --lrt-sig-fc-plot {input.lrt_sig_fc_plot} \
          --lrt-sig-volcano-plot {input.lrt_sig_volcano} \
          --enrichment-stats {input.enrichment_stats} \
          --structure-plot {input.struct_plot} \
          --out-html {output.html} \
          --out-versions {output.versions} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 12: ATAC–RNA correlation by promoter group across UV time-course ────

rule g4_tss_atac_rna_corr:
    input:
        annotation=rules.g4_tss_build_annotation.output.tsv,
        windows=rules.g4_tss_slop_windows.output.windows,
        norm_counts="results/rnaseq/matrices/normalized_counts.tsv.gz",
        bw_t00="resources/rna-seq/merged_t00.bw",
        bw_t12="resources/rna-seq/merged_t12.bw",
        bw_t30="resources/rna-seq/merged_t30.bw",
        bw_t60="resources/rna-seq/merged_t60.bw",
    output:
        atac_signal="results/g4_tss_atac_rna_corr/atac_promoter_signal.tsv",
        rna_expression="results/g4_tss_atac_rna_corr/rna_expression_by_timepoint.tsv",
        merged="results/g4_tss_atac_rna_corr/merged_atac_rna_by_group.tsv",
        silent_genes="results/g4_tss_atac_rna_corr/silent_genes.tsv",
        spearman="results/g4_tss_atac_rna_corr/atac_rna_spearman_by_group_timepoint.tsv",
        permutation="results/g4_tss_atac_rna_corr/group_correlation_permutation_tests.tsv",
        summary="results/g4_tss_atac_rna_corr/atac_rna_corr_summary.tsv",
        scatter_pdf="results/g4_tss_atac_rna_corr/scatter_atac_vs_rna_by_group_timepoint.pdf",
        heatmap_pdf="results/g4_tss_atac_rna_corr/rho_heatmap_group_timepoint.pdf",
        lineplot_pdf="results/g4_tss_atac_rna_corr/rho_over_time_by_group.pdf",
    log:
        "logs/g4_tss_atac_rna_corr/atac_rna_corr.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Correlating promoter ATAC signal with RNA-seq expression..." &&
        mkdir -p results/g4_tss_atac_rna_corr logs/g4_tss_atac_rna_corr &&
        python3 workflow/scripts/g4_tss_atac_rna_corr.py \
          --tss-annotation {input.annotation} \
          --promoter-windows-bed {input.windows} \
          --norm-counts {input.norm_counts} \
          --bw-t00 {input.bw_t00} \
          --bw-t12 {input.bw_t12} \
          --bw-t30 {input.bw_t30} \
          --bw-t60 {input.bw_t60} \
          --out-atac-signal {output.atac_signal} \
          --out-rna-expression {output.rna_expression} \
          --out-merged {output.merged} \
          --out-silent-genes {output.silent_genes} \
          --out-spearman {output.spearman} \
          --out-permutation {output.permutation} \
          --out-summary {output.summary} \
          --out-scatter-pdf {output.scatter_pdf} \
          --out-heatmap-pdf {output.heatmap_pdf} \
          --out-lineplot-pdf {output.lineplot_pdf} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """
