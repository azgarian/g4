"""
G4-TSS UV lesion burden and transcriptional suppression analysis — Tasks 1-7.

Answers: Do genes with higher UV lesion burden at their promoter G4 loci show
stronger transcriptional suppression after UV irradiation?

Sign convention: positive log2FoldChange = higher at t=0 than post-UV = repression.

Prerequisites:
  results/g4_tss/tss_group_annotation.tsv
  results/g4_tss/canonical_tss_windows_1kb.bed
  results/g4_tss/baseline_tpm.tsv
  results/g4chip/g4_hela_peaks_prepared.tsv
  results/g4cuttag/g4_hela_peaks_prepared.tsv
  results/g4chip_g4cuttag/g4_hela_chip_cuttag_merged.bed
  results/g4_tss_uv/uv_master_table.tsv
  results/rnaseq/deseq2/pairwise/*.tsv.gz
  results/gc_rich_bg/gc_rich_bg_prepared.bed
  resources/samples/xr_ds/*.bed
"""

import json
from pathlib import Path

XR_DS_DIR = "resources/samples/xr_ds"


def _ds_sample_manifest(cfg):
    """Build sample manifest for all DS-method samples with existing bed files."""
    manifest = []
    for sample_id, meta in cfg.get("samples", {}).items():
        if meta.get("method") != "DS":
            continue
        base = f"{XR_DS_DIR}/{sample_id}"
        entry = {
            "sample_id": sample_id,
            "product": meta["product"],
            "time_after_exposure": int(meta["time_after_exposure"]),
            "name": meta["name"],
            "real_plus_bed": f"{base}_plus.bed",
            "real_minus_bed": f"{base}_minus.bed",
            "sim_plus_bed": f"{base}_sim_plus.bed",
            "sim_minus_bed": f"{base}_sim_minus.bed",
        }
        manifest.append(entry)
    return manifest


def _ds_bed_inputs(cfg):
    """Return all existing DS bed file paths for Snakemake input declarations."""
    beds = []
    for entry in _ds_sample_manifest(cfg):
        for key in ["real_plus_bed", "real_minus_bed", "sim_plus_bed", "sim_minus_bed"]:
            p = Path(entry[key])
            if p.exists():
                beds.append(str(p))
    return beds


# ─── Manifest: write sample manifest JSON ────────────────────────────────────

rule g4_tss_damage_uv_manifest:
    output:
        manifest="results/g4_tss_damage_uv/ds_sample_manifest.json",
    run:
        Path(output.manifest).parent.mkdir(parents=True, exist_ok=True)
        manifest = _ds_sample_manifest(config)
        Path(output.manifest).write_text(json.dumps(manifest, indent=2))


# ─── Task 1: Build promoter-G4 locus and gene lesion tables ──────────────────

rule g4_tss_damage_uv_assemble:
    input:
        tss_windows="results/g4_tss/canonical_tss_windows_1kb.bed",
        g4_merged="results/g4chip_g4cuttag/g4_hela_chip_cuttag_merged.bed",
        annotation="results/g4_tss/tss_group_annotation.tsv",
        chip_peaks="results/g4chip/g4_hela_peaks_prepared.tsv",
        cuttag_peaks="results/g4cuttag/g4_hela_peaks_prepared.tsv",
        gc_rich_bg="results/gc_rich_bg/gc_rich_bg_prepared.bed",
        manifest=rules.g4_tss_damage_uv_manifest.output.manifest,
        ds_beds=lambda wc: _ds_bed_inputs(config),
    output:
        locus_map="results/g4_tss_damage_uv/promoter_g4_gene_locus_map.tsv",
        locus_table="results/g4_tss_damage_uv/promoter_g4_damage_locus_table.tsv",
        gene_table="results/g4_tss_damage_uv/promoter_g4_damage_gene_table.tsv",
    log:
        "logs/g4_tss_damage_uv/task1_damage_assembly.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Building promoter-G4 locus and gene lesion tables..." &&
        mkdir -p results/g4_tss_damage_uv logs/g4_tss_damage_uv &&
        python3 workflow/scripts/g4_tss_damage_uv_task1_assemble.py \
          --tss-windows {input.tss_windows} \
          --g4-merged {input.g4_merged} \
          --annotation {input.annotation} \
          --chip-peaks {input.chip_peaks} \
          --cuttag-peaks {input.cuttag_peaks} \
          --gc-rich-bg {input.gc_rich_bg} \
          --sample-manifest {input.manifest} \
          --out-locus-map {output.locus_map} \
          --out-locus-table {output.locus_table} \
          --out-gene-table {output.gene_table} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1
        """


# ─── Task 2: Assemble lesion × RNA master table ───────────────────────────────

rule g4_tss_damage_uv_master_join:
    input:
        gene_table=rules.g4_tss_damage_uv_assemble.output.gene_table,
        uv_master="results/g4_tss_uv/uv_master_table.tsv",
        baseline_tpm="results/g4_tss/baseline_tpm.tsv",
        tss_windows="results/g4_tss/canonical_tss_windows_1kb.bed",
        chip_peaks="results/g4chip/g4_hela_peaks_prepared.tsv",
        cuttag_peaks="results/g4cuttag/g4_hela_peaks_prepared.tsv",
    output:
        cpd_table="results/g4_tss_damage_uv/cpd_damage_suppression_master_table.tsv",
        pp64_table="results/g4_tss_damage_uv/64pp_damage_suppression_master_table.tsv",
    log:
        "logs/g4_tss_damage_uv/task2_master_join.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Assembling lesion x RNA master table..." &&
        python3 workflow/scripts/g4_tss_damage_uv_task2_master_join.py \
          --gene-table {input.gene_table} \
          --uv-master {input.uv_master} \
          --baseline-tpm {input.baseline_tpm} \
          --tss-windows {input.tss_windows} \
          --g4chip-peaks {input.chip_peaks} \
          --g4cuttag-peaks {input.cuttag_peaks} \
          --out-cpd {output.cpd_table} \
          --out-64pp {output.pp64_table} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 3: Primary burden-stratified test ───────────────────────────────────

rule g4_tss_damage_uv_burden_tertiles:
    input:
        cpd_table=rules.g4_tss_damage_uv_master_join.output.cpd_table,
        pp64_table=rules.g4_tss_damage_uv_master_join.output.pp64_table,
    output:
        cpd_stats="results/g4_tss_damage_uv/cpd_burden_repression_stats.tsv",
        pp64_stats="results/g4_tss_damage_uv/64pp_burden_repression_stats.tsv",
        cpd_summary="results/g4_tss_damage_uv/cpd_burden_repression_summary.tsv",
        pp64_summary="results/g4_tss_damage_uv/64pp_burden_repression_summary.tsv",
        cpd_lfc_fig="results/g4_tss_damage_uv/cpd_burden_lfc_by_tertile.pdf",
        pp64_lfc_fig="results/g4_tss_damage_uv/64pp_burden_lfc_by_tertile.pdf",
        cpd_traj_fig="results/g4_tss_damage_uv/cpd_burden_median_trajectory.pdf",
        pp64_traj_fig="results/g4_tss_damage_uv/64pp_burden_median_trajectory.pdf",
        cpd_comp_fig="results/g4_tss_damage_uv/cpd_burden_trajectory_composition.pdf",
        pp64_comp_fig="results/g4_tss_damage_uv/64pp_burden_trajectory_composition.pdf",
    log:
        "logs/g4_tss_damage_uv/task3_burden_tertiles.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Running burden-stratified tests..." &&
        python3 workflow/scripts/g4_tss_damage_uv_task3_burden_tertiles.py \
          --cpd-table {input.cpd_table} \
          --pp64-table {input.pp64_table} \
          --out-cpd-stats {output.cpd_stats} \
          --out-64pp-stats {output.pp64_stats} \
          --out-cpd-summary {output.cpd_summary} \
          --out-64pp-summary {output.pp64_summary} \
          --out-cpd-lfc-fig {output.cpd_lfc_fig} \
          --out-64pp-lfc-fig {output.pp64_lfc_fig} \
          --out-cpd-traj-fig {output.cpd_traj_fig} \
          --out-64pp-traj-fig {output.pp64_traj_fig} \
          --out-cpd-comp-fig {output.cpd_comp_fig} \
          --out-64pp-comp-fig {output.pp64_comp_fig} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 4: Continuous association and covariate-adjusted modeling ───────────

rule g4_tss_damage_uv_covariate_models:
    input:
        cpd_table=rules.g4_tss_damage_uv_master_join.output.cpd_table,
        pp64_table=rules.g4_tss_damage_uv_master_join.output.pp64_table,
    output:
        cpd_corr="results/g4_tss_damage_uv/cpd_burden_correlations.tsv",
        pp64_corr="results/g4_tss_damage_uv/64pp_burden_correlations.tsv",
        cpd_partial="results/g4_tss_damage_uv/cpd_partial_spearman.tsv",
        pp64_partial="results/g4_tss_damage_uv/64pp_partial_spearman.tsv",
        cpd_models="results/g4_tss_damage_uv/cpd_covariate_models.tsv",
        pp64_models="results/g4_tss_damage_uv/64pp_covariate_models.tsv",
        cpd_model_cmp="results/g4_tss_damage_uv/cpd_covariate_model_comparison.tsv",
        pp64_model_cmp="results/g4_tss_damage_uv/64pp_covariate_model_comparison.tsv",
        cpd_scatter="results/g4_tss_damage_uv/cpd_burden_lfc60_scatter.pdf",
        pp64_scatter="results/g4_tss_damage_uv/64pp_burden_lfc60_scatter.pdf",
        cpd_partial_fig="results/g4_tss_damage_uv/cpd_partial_effect_lfc60.pdf",
        pp64_partial_fig="results/g4_tss_damage_uv/64pp_partial_effect_lfc60.pdf",
        cpd_coef_fig="results/g4_tss_damage_uv/cpd_covariate_coefficients.pdf",
        pp64_coef_fig="results/g4_tss_damage_uv/64pp_covariate_coefficients.pdf",
    log:
        "logs/g4_tss_damage_uv/task4_covariate_models.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Fitting continuous association and covariate models..." &&
        python3 workflow/scripts/g4_tss_damage_uv_task4_covariate_models.py \
          --cpd-table {input.cpd_table} \
          --pp64-table {input.pp64_table} \
          --out-cpd-corr {output.cpd_corr} \
          --out-64pp-corr {output.pp64_corr} \
          --out-cpd-partial {output.cpd_partial} \
          --out-64pp-partial {output.pp64_partial} \
          --out-cpd-models {output.cpd_models} \
          --out-64pp-models {output.pp64_models} \
          --out-cpd-model-cmp {output.cpd_model_cmp} \
          --out-64pp-model-cmp {output.pp64_model_cmp} \
          --out-cpd-scatter {output.cpd_scatter} \
          --out-64pp-scatter {output.pp64_scatter} \
          --out-cpd-partial-fig {output.cpd_partial_fig} \
          --out-64pp-partial-fig {output.pp64_partial_fig} \
          --out-cpd-coef-fig {output.cpd_coef_fig} \
          --out-64pp-coef-fig {output.pp64_coef_fig} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 5: Trajectory-level and sustained-repression analysis ───────────────

rule g4_tss_damage_uv_trajectory:
    input:
        cpd_table=rules.g4_tss_damage_uv_master_join.output.cpd_table,
        pp64_table=rules.g4_tss_damage_uv_master_join.output.pp64_table,
    output:
        cpd_kruskal="results/g4_tss_damage_uv/cpd_trajectory_burden_kruskal.tsv",
        pp64_kruskal="results/g4_tss_damage_uv/64pp_trajectory_burden_kruskal.tsv",
        cpd_pairwise="results/g4_tss_damage_uv/cpd_trajectory_burden_pairwise.tsv",
        pp64_pairwise="results/g4_tss_damage_uv/64pp_trajectory_burden_pairwise.tsv",
        cpd_violin="results/g4_tss_damage_uv/cpd_trajectory_burden_violin.pdf",
        pp64_violin="results/g4_tss_damage_uv/64pp_trajectory_burden_violin.pdf",
        cpd_logistic="results/g4_tss_damage_uv/cpd_trajectory_burden_logistic.tsv",
        pp64_logistic="results/g4_tss_damage_uv/64pp_trajectory_burden_logistic.tsv",
    log:
        "logs/g4_tss_damage_uv/task5_trajectory.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Running trajectory-level burden analysis..." &&
        python3 workflow/scripts/g4_tss_damage_uv_task5_trajectory.py \
          --cpd-table {input.cpd_table} \
          --pp64-table {input.pp64_table} \
          --out-cpd-kruskal {output.cpd_kruskal} \
          --out-64pp-kruskal {output.pp64_kruskal} \
          --out-cpd-pairwise {output.cpd_pairwise} \
          --out-64pp-pairwise {output.pp64_pairwise} \
          --out-cpd-violin {output.cpd_violin} \
          --out-64pp-violin {output.pp64_violin} \
          --out-cpd-logistic {output.cpd_logistic} \
          --out-64pp-logistic {output.pp64_logistic} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 6: Specificity and sensitivity checks ───────────────────────────────

rule g4_tss_damage_uv_specificity:
    input:
        cpd_table=rules.g4_tss_damage_uv_master_join.output.cpd_table,
        pp64_table=rules.g4_tss_damage_uv_master_join.output.pp64_table,
        gene_table=rules.g4_tss_damage_uv_assemble.output.gene_table,
        annotation="results/g4_tss/tss_group_annotation.tsv",
        tss_windows="results/g4_tss/canonical_tss_windows_1kb.bed",
        gc_rich_bg="results/gc_rich_bg/gc_rich_bg_prepared.bed",
        uv_master="results/g4_tss_uv/uv_master_table.tsv",
        baseline_tpm="results/g4_tss/baseline_tpm.tsv",
        manifest=rules.g4_tss_damage_uv_manifest.output.manifest,
        ds_beds=lambda wc: _ds_bed_inputs(config),
    output:
        interaction="results/g4_tss_damage_uv/promoter_group_damage_interaction.tsv",
        matched_summary="results/g4_tss_damage_uv/promoter_group_damage_matched_summary.tsv",
        cpd_sensitivity="results/g4_tss_damage_uv/cpd_burden_sensitivity.tsv",
        pp64_sensitivity="results/g4_tss_damage_uv/64pp_burden_sensitivity.tsv",
        interaction_fig="results/g4_tss_damage_uv/promoter_group_damage_interaction.pdf",
        gc_bg_summary="results/g4_tss_damage_uv/gc_bg_damage_summary.tsv",
        gc_bg_lfc="results/g4_tss_damage_uv/gc_bg_damage_lfc_correlation.tsv",
        gc_bg_vs_g4="results/g4_tss_damage_uv/gc_bg_vs_g4_damage_comparison.tsv",
    log:
        "logs/g4_tss_damage_uv/task6_specificity_sensitivity.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Running specificity and sensitivity checks..." &&
        python3 workflow/scripts/g4_tss_damage_uv_task6_specificity.py \
          --cpd-table {input.cpd_table} \
          --pp64-table {input.pp64_table} \
          --gene-table {input.gene_table} \
          --annotation {input.annotation} \
          --tss-windows {input.tss_windows} \
          --gc-rich-bg {input.gc_rich_bg} \
          --uv-master {input.uv_master} \
          --baseline-tpm {input.baseline_tpm} \
          --sample-manifest {input.manifest} \
          --out-interaction {output.interaction} \
          --out-matched-summary {output.matched_summary} \
          --out-cpd-sensitivity {output.cpd_sensitivity} \
          --out-64pp-sensitivity {output.pp64_sensitivity} \
          --out-interaction-fig {output.interaction_fig} \
          --out-gc-bg-summary {output.gc_bg_summary} \
          --out-gc-bg-lfc {output.gc_bg_lfc} \
          --out-gc-bg-vs-g4 {output.gc_bg_vs_g4} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


# ─── Task 7: Synthesis HTML report ───────────────────────────────────────────

rule g4_tss_damage_uv_report:
    input:
        locus_map=rules.g4_tss_damage_uv_assemble.output.locus_map,
        gene_table=rules.g4_tss_damage_uv_assemble.output.gene_table,
        cpd_table=rules.g4_tss_damage_uv_master_join.output.cpd_table,
        pp64_table=rules.g4_tss_damage_uv_master_join.output.pp64_table,
        cpd_stats=rules.g4_tss_damage_uv_burden_tertiles.output.cpd_stats,
        pp64_stats=rules.g4_tss_damage_uv_burden_tertiles.output.pp64_stats,
        cpd_summary=rules.g4_tss_damage_uv_burden_tertiles.output.cpd_summary,
        pp64_summary=rules.g4_tss_damage_uv_burden_tertiles.output.pp64_summary,
        cpd_lfc_fig=rules.g4_tss_damage_uv_burden_tertiles.output.cpd_lfc_fig,
        pp64_lfc_fig=rules.g4_tss_damage_uv_burden_tertiles.output.pp64_lfc_fig,
        cpd_traj_fig=rules.g4_tss_damage_uv_burden_tertiles.output.cpd_traj_fig,
        pp64_traj_fig=rules.g4_tss_damage_uv_burden_tertiles.output.pp64_traj_fig,
        cpd_comp_fig=rules.g4_tss_damage_uv_burden_tertiles.output.cpd_comp_fig,
        pp64_comp_fig=rules.g4_tss_damage_uv_burden_tertiles.output.pp64_comp_fig,
        cpd_corr=rules.g4_tss_damage_uv_covariate_models.output.cpd_corr,
        pp64_corr=rules.g4_tss_damage_uv_covariate_models.output.pp64_corr,
        cpd_partial=rules.g4_tss_damage_uv_covariate_models.output.cpd_partial,
        pp64_partial=rules.g4_tss_damage_uv_covariate_models.output.pp64_partial,
        cpd_models=rules.g4_tss_damage_uv_covariate_models.output.cpd_models,
        pp64_models=rules.g4_tss_damage_uv_covariate_models.output.pp64_models,
        cpd_model_cmp=rules.g4_tss_damage_uv_covariate_models.output.cpd_model_cmp,
        pp64_model_cmp=rules.g4_tss_damage_uv_covariate_models.output.pp64_model_cmp,
        cpd_scatter=rules.g4_tss_damage_uv_covariate_models.output.cpd_scatter,
        pp64_scatter=rules.g4_tss_damage_uv_covariate_models.output.pp64_scatter,
        cpd_partial_fig=rules.g4_tss_damage_uv_covariate_models.output.cpd_partial_fig,
        pp64_partial_fig=rules.g4_tss_damage_uv_covariate_models.output.pp64_partial_fig,
        cpd_coef_fig=rules.g4_tss_damage_uv_covariate_models.output.cpd_coef_fig,
        pp64_coef_fig=rules.g4_tss_damage_uv_covariate_models.output.pp64_coef_fig,
        cpd_kruskal=rules.g4_tss_damage_uv_trajectory.output.cpd_kruskal,
        pp64_kruskal=rules.g4_tss_damage_uv_trajectory.output.pp64_kruskal,
        cpd_traj_pairwise=rules.g4_tss_damage_uv_trajectory.output.cpd_pairwise,
        pp64_traj_pairwise=rules.g4_tss_damage_uv_trajectory.output.pp64_pairwise,
        cpd_traj_violin=rules.g4_tss_damage_uv_trajectory.output.cpd_violin,
        pp64_traj_violin=rules.g4_tss_damage_uv_trajectory.output.pp64_violin,
        cpd_traj_logistic=rules.g4_tss_damage_uv_trajectory.output.cpd_logistic,
        pp64_traj_logistic=rules.g4_tss_damage_uv_trajectory.output.pp64_logistic,
        interaction=rules.g4_tss_damage_uv_specificity.output.interaction,
        interaction_fig=rules.g4_tss_damage_uv_specificity.output.interaction_fig,
        gc_bg_summary=rules.g4_tss_damage_uv_specificity.output.gc_bg_summary,
        gc_bg_vs_g4=rules.g4_tss_damage_uv_specificity.output.gc_bg_vs_g4,
    output:
        html="results/g4_tss_damage_uv/g4_promoter_damage_suppression_report.html",
        versions="results/g4_tss_damage_uv/software_versions.txt",
    log:
        "logs/g4_tss_damage_uv/task7_report.log",
    conda:
        "../envs/g4_tss.yaml"
    shell:
        """
        (echo "`date -R`: Generating synthesis HTML report..." &&
        python3 workflow/scripts/g4_tss_damage_uv_task7_report.py \
          --locus-map {input.locus_map} \
          --gene-table {input.gene_table} \
          --cpd-table {input.cpd_table} \
          --pp64-table {input.pp64_table} \
          --cpd-tertile-stats {input.cpd_stats} \
          --pp64-tertile-stats {input.pp64_stats} \
          --cpd-tertile-summary {input.cpd_summary} \
          --pp64-tertile-summary {input.pp64_summary} \
          --cpd-lfc-fig {input.cpd_lfc_fig} \
          --pp64-lfc-fig {input.pp64_lfc_fig} \
          --cpd-traj-fig {input.cpd_traj_fig} \
          --pp64-traj-fig {input.pp64_traj_fig} \
          --cpd-comp-fig {input.cpd_comp_fig} \
          --pp64-comp-fig {input.pp64_comp_fig} \
          --cpd-corr {input.cpd_corr} \
          --pp64-corr {input.pp64_corr} \
          --cpd-partial {input.cpd_partial} \
          --pp64-partial {input.pp64_partial} \
          --cpd-models {input.cpd_models} \
          --pp64-models {input.pp64_models} \
          --cpd-model-cmp {input.cpd_model_cmp} \
          --pp64-model-cmp {input.pp64_model_cmp} \
          --cpd-scatter {input.cpd_scatter} \
          --pp64-scatter {input.pp64_scatter} \
          --cpd-partial-fig {input.cpd_partial_fig} \
          --pp64-partial-fig {input.pp64_partial_fig} \
          --cpd-coef-fig {input.cpd_coef_fig} \
          --pp64-coef-fig {input.pp64_coef_fig} \
          --cpd-traj-kruskal {input.cpd_kruskal} \
          --pp64-traj-kruskal {input.pp64_kruskal} \
          --cpd-traj-pairwise {input.cpd_traj_pairwise} \
          --pp64-traj-pairwise {input.pp64_traj_pairwise} \
          --cpd-traj-violin {input.cpd_traj_violin} \
          --pp64-traj-violin {input.pp64_traj_violin} \
          --cpd-traj-logistic {input.cpd_traj_logistic} \
          --pp64-traj-logistic {input.pp64_traj_logistic} \
          --interaction {input.interaction} \
          --interaction-fig {input.interaction_fig} \
          --gc-bg-summary {input.gc_bg_summary} \
          --gc-bg-vs-g4 {input.gc_bg_vs_g4} \
          --out-html {output.html} \
          --out-versions {output.versions} \
          --log {log} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """
