"""
G4-TSS UV magnitude and direction analysis module — Tasks 1-7.

Answers: Do genes with a G4 at their promoter change expression more strongly
(magnitude) or in a different direction (repression vs induction) after UV,
compared with non-G4 genes?

Prerequisites (must already exist):
  results/g4_tss/tss_group_annotation.tsv
  results/g4_tss/baseline_tpm.tsv
  results/g4_tss/canonical_tss_windows_1kb.bed
  results/rnaseq/deseq2/time_lrt_results.tsv.gz
  results/rnaseq/deseq2/pairwise/0_vs_12_results.tsv.gz
  results/rnaseq/deseq2/pairwise/0_vs_30_results.tsv.gz
  results/rnaseq/deseq2/pairwise/0_vs_60_results.tsv.gz
  results/g4chip/g4_hela_peaks_prepared.tsv
  results/g4cuttag/g4_hela_peaks_prepared.tsv

Sign convention: positive log2FoldChange = higher at t=0 than post-UV = repression.
"""

# ─── Task 1: Assemble master UV-response table ─────────────────────────────────

rule g4_tss_uv_assemble:
    input:
        annotation="results/g4_tss/tss_group_annotation.tsv",
        lrt="results/rnaseq/deseq2/time_lrt_results.tsv.gz",
        pw12="results/rnaseq/deseq2/pairwise/0_vs_12_results.tsv.gz",
        pw30="results/rnaseq/deseq2/pairwise/0_vs_30_results.tsv.gz",
        pw60="results/rnaseq/deseq2/pairwise/0_vs_60_results.tsv.gz",
        tpm="results/g4_tss/baseline_tpm.tsv",
    output:
        table="results/g4_tss_uv/uv_master_table.tsv",
    log:
        "logs/g4_tss_uv/task1_assembly.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Assembling master UV-response table..." &&
        mkdir -p results/g4_tss_uv logs/g4_tss_uv &&
        python3 workflow/scripts/g4_tss_uv_task1_assemble.py \
          --tss-annotation {input.annotation} \
          --lrt-results {input.lrt} \
          --pairwise-12 {input.pw12} \
          --pairwise-30 {input.pw30} \
          --pairwise-60 {input.pw60} \
          --baseline-tpm {input.tpm} \
          --out-table {output.table} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1
        """


# ─── Task 2: Fold-change magnitude comparison ──────────────────────────────────

rule g4_tss_uv_magnitude:
    input:
        master=rules.g4_tss_uv_assemble.output.table,
    output:
        stats="results/g4_tss_uv/uv_magnitude_stats.tsv",
        summary="results/g4_tss_uv/uv_magnitude_summary.tsv",
        figure="results/g4_tss_uv/uv_magnitude_by_group.pdf",
    log:
        "logs/g4_tss_uv/task2_magnitude.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Comparing UV fold-change magnitude across groups..." &&
        python3 workflow/scripts/g4_tss_uv_task2_magnitude.py \
          --master-table {input.master} \
          --out-stats {output.stats} \
          --out-summary {output.summary} \
          --out-figure {output.figure} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 3: Direction asymmetry ───────────────────────────────────────────────

rule g4_tss_uv_direction:
    input:
        master=rules.g4_tss_uv_assemble.output.table,
    output:
        stats="results/g4_tss_uv/uv_direction_stats.tsv",
        fisher="results/g4_tss_uv/uv_direction_fisher.tsv",
        stacked_bar="results/g4_tss_uv/uv_direction_stacked_bar.pdf",
        density="results/g4_tss_uv/uv_signed_lfc_density.pdf",
    log:
        "logs/g4_tss_uv/task3_direction.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Quantifying UV response direction asymmetry..." &&
        python3 workflow/scripts/g4_tss_uv_task3_direction.py \
          --master-table {input.master} \
          --out-direction-stats {output.stats} \
          --out-fisher {output.fisher} \
          --out-stacked-bar {output.stacked_bar} \
          --out-density {output.density} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 4: Trajectory analysis ──────────────────────────────────────────────

rule g4_tss_uv_trajectory:
    input:
        master=rules.g4_tss_uv_assemble.output.table,
    output:
        counts="results/g4_tss_uv/uv_trajectory_counts.tsv",
        fisher="results/g4_tss_uv/uv_trajectory_fisher.tsv",
        plot="results/g4_tss_uv/uv_trajectory_median_lfc.pdf",
        baseline_expr="results/g4_tss_uv/uv_trajectory_baseline_expr.tsv",
    log:
        "logs/g4_tss_uv/task4_trajectory.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Running trajectory analysis..." &&
        python3 workflow/scripts/g4_tss_uv_task4_trajectory.py \
          --master-table {input.master} \
          --out-trajectory-counts {output.counts} \
          --out-trajectory-fisher {output.fisher} \
          --out-trajectory-plot {output.plot} \
          --out-baseline-expr {output.baseline_expr} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 5: G4 strength as continuous predictor ──────────────────────────────

rule g4_tss_uv_g4strength:
    input:
        master=rules.g4_tss_uv_assemble.output.table,
        tss_windows="results/g4_tss/canonical_tss_windows_1kb.bed",
        chip_peaks="results/g4chip/g4_hela_peaks_prepared.tsv",
        cuttag_peaks="results/g4cuttag/g4_hela_peaks_prepared.tsv",
    output:
        correlations="results/g4_tss_uv/g4_strength_uv_correlations.tsv",
        direction_tertile="results/g4_tss_uv/g4_strength_direction_by_tertile.tsv",
        scatter="results/g4_tss_uv/g4_strength_lfc_scatter.pdf",
        violin="results/g4_tss_uv/g4_strength_lfc_violin.pdf",
    log:
        "logs/g4_tss_uv/task5_g4strength.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Testing G4 strength as continuous UV-response predictor..." &&
        python3 workflow/scripts/g4_tss_uv_task5_g4strength.py \
          --master-table {input.master} \
          --tss-windows {input.tss_windows} \
          --g4chip-peaks {input.chip_peaks} \
          --g4cuttag-peaks {input.cuttag_peaks} \
          --out-correlations {output.correlations} \
          --out-direction-by-tertile {output.direction_tertile} \
          --out-scatter {output.scatter} \
          --out-violin {output.violin} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 6: Baseline-expression-matched comparison ───────────────────────────

rule g4_tss_uv_expr_matched:
    input:
        master=rules.g4_tss_uv_assemble.output.table,
    output:
        wilcoxon="results/g4_tss_uv/uv_magnitude_expr_matched.tsv",
        lm="results/g4_tss_uv/uv_magnitude_lm.tsv",
        plot="results/g4_tss_uv/uv_magnitude_by_bin.pdf",
    log:
        "logs/g4_tss_uv/task6_expr_matched.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Running baseline-expression-matched comparison..." &&
        python3 workflow/scripts/g4_tss_uv_task6_expr_matched.py \
          --master-table {input.master} \
          --out-wilcoxon {output.wilcoxon} \
          --out-lm {output.lm} \
          --out-bin-plot {output.plot} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 7: Synthesis HTML report ────────────────────────────────────────────

rule g4_tss_uv_report:
    input:
        magnitude_stats=rules.g4_tss_uv_magnitude.output.stats,
        magnitude_summary=rules.g4_tss_uv_magnitude.output.summary,
        magnitude_figure=rules.g4_tss_uv_magnitude.output.figure,
        direction_stats=rules.g4_tss_uv_direction.output.stats,
        direction_fisher=rules.g4_tss_uv_direction.output.fisher,
        stacked_bar=rules.g4_tss_uv_direction.output.stacked_bar,
        density=rules.g4_tss_uv_direction.output.density,
        trajectory_counts=rules.g4_tss_uv_trajectory.output.counts,
        trajectory_fisher=rules.g4_tss_uv_trajectory.output.fisher,
        trajectory_plot=rules.g4_tss_uv_trajectory.output.plot,
        g4_correlations=rules.g4_tss_uv_g4strength.output.correlations,
        g4_direction_tertile=rules.g4_tss_uv_g4strength.output.direction_tertile,
        g4_scatter=rules.g4_tss_uv_g4strength.output.scatter,
        g4_violin=rules.g4_tss_uv_g4strength.output.violin,
        expr_wilcoxon=rules.g4_tss_uv_expr_matched.output.wilcoxon,
        expr_lm=rules.g4_tss_uv_expr_matched.output.lm,
        expr_plot=rules.g4_tss_uv_expr_matched.output.plot,
    output:
        html="results/g4_tss_uv/g4_uv_magnitude_direction_report.html",
        versions="results/g4_tss_uv/software_versions.txt",
    log:
        "logs/g4_tss_uv/task7_report.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Generating synthesis HTML report..." &&
        python3 workflow/scripts/g4_tss_uv_task7_report.py \
          --magnitude-stats {input.magnitude_stats} \
          --magnitude-summary {input.magnitude_summary} \
          --magnitude-figure {input.magnitude_figure} \
          --direction-stats {input.direction_stats} \
          --direction-fisher {input.direction_fisher} \
          --stacked-bar {input.stacked_bar} \
          --density-plot {input.density} \
          --trajectory-counts {input.trajectory_counts} \
          --trajectory-fisher {input.trajectory_fisher} \
          --trajectory-plot {input.trajectory_plot} \
          --g4-correlations {input.g4_correlations} \
          --g4-direction-tertile {input.g4_direction_tertile} \
          --g4-scatter {input.g4_scatter} \
          --g4-violin {input.g4_violin} \
          --expr-matched-wilcoxon {input.expr_wilcoxon} \
          --expr-matched-lm {input.expr_lm} \
          --expr-matched-plot {input.expr_plot} \
          --out-html {output.html} \
          --out-versions {output.versions} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """
