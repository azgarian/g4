#!/bin/env python

import os
import warnings

################### Helper Functions ###########################################

def allInput():

    inputList = []

    # smk: g4p_cut_tag
    inputList.append("results/g4cuttag/g4_hela_peaks_prepared.bed")

    # smk: g4p_chip_seq
    inputList.append("results/g4chip/g4_hela_peaks_prepared.bed")

    # smk: g4chip_g4cuttag_merge
    inputList.append("results/g4chip_g4cuttag/g4_hela_chip_cuttag_merged.bed")

    # smk: prepare_oqs
    inputList.append("results/oqs/oqs_K_prepared.bed")
    inputList.append("results/oqs/oqs_PDS_prepared.bed")

    # smk: gc_rich_bg
    inputList.append("results/gc_rich_bg/gc_rich_bg_prepared.bed")
    inputList.append("results/gc_rich_bg/gc_rich_bg_selection_summary.tsv")

    # smk: g4_tss (Tasks 1-8)
    inputList.append("results/g4_tss/tss_group_annotation.tsv")
    inputList.append("results/g4_tss/gene_expression_by_group.tsv")
    inputList.append("results/g4_tss/expression_group_statistics.tsv")
    inputList.append("results/g4_tss/expression_violin_by_group.pdf")
    inputList.append("results/g4_tss/rnaseq_tss_metaprofile.pdf")
    inputList.append("results/g4_tss/decile_enrichment_plot.pdf")
    inputList.append("results/g4_tss/uv_volcano_by_group.pdf")
    inputList.append("results/g4_tss/structure_class_by_expression_class.pdf")
    inputList.append("results/g4_tss/g4_tss_transcription_report.html")

    # smk: g4_tss (Tasks 9-11 — pG4 × LRT intersection, enrichment, sensitivity)
    inputList.append("results/g4_tss/pG4_DEG_intersect.tsv")
    inputList.append("results/g4_tss/pG4_pathway_enrichment.tsv")
    inputList.append("results/g4_tss/pG4_strength_stratification.tsv")
    inputList.append("results/g4_tss/pG4_strength_stratification.pdf")
    inputList.append("results/g4_tss/pG4_enrichment_summary.tsv")
    inputList.append("results/g4_tss/pG4_sensitivity_qc.tsv")

    # smk: g4_tss_atac_split
    for timepoint in config.get("g4_tss_atac_split", {}).get("timepoints", []):
        inputList.append(
            f"results/g4_tss_atac_split/{timepoint}/g4_tss_atac_split_report.html"
        )

    return inputList
