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

    return inputList
