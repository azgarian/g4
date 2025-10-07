#!/bin/env python

#### configuration file ####
configfile: "config/config.yaml"

# Optional: singularity image to use
# containerized: "docker://your-image:tag"

# Include common functions and utilities
# Shared helpers (e.g., small Python utilities) available to all rules
include: "workflow/rules/common.smk"

# Wildcard constraints
# Keep patterns broad but unambiguous; tailor to your sample naming scheme
wildcard_constraints:
    sample="[^/]+",

# Final targets for this workflow; dependencies are inferred from inputs/outputs
rule all:
    input:
        expand("results/master/{region}/mapped_all.csv", region=config["region_names"]),
        expand("results/g4_repair/{region}/{region}_64_normalized.png",
            region=[r for r in config["region_names"] if "oqs" in r]),
        # expand("results/deeptools/{region}_atac/profile.pdf", region=["oqs_0", "oqs_12", "oqs_30", "oqs_60", "oqs_lost_after_uv", "oqs_formed_after_uv", "oqs_persistent", "oqs_never_formed", "oqs_cur_alg"]),
        # expand("results/g4_miner/{region}.promoter_info.tsv", region=["oqs_0", "oqs_12", "oqs_30", "oqs_60"]),
        "results/g4_txn/g4_promoter_counts.tsv",
        # "resources/Hela_15m_R3_ATAC_AGGCAGAA-TATCCTCT_S2.shifted.bw",
        # "resources/Hela_1h_R3_ATAC.shifted.bw",
        # "resources/Hela_1m_R3_ATAC_AGGCAGAA-CTCTCTAT_S3.shifted.bw",
        # "resources/Hela_30m_R3_ATAC.shifted.bw",
        # "resources/Hela_4h_R3_ATAC_TCCTGAGC-CTCTCTAT_S4.shifted.bw",
        # "resources/Hela_8h_R3_ATAC_TCCTGAGC-TATCCTCT_S6.shifted.bw",
        # "resources/HelanoUV_R1_ATAC_TAAGGCGA-TAGATCGC_S7.shifted.bw",


include: "workflow/rules/mapping.smk"

# deepTools rules
include: "workflow/rules/deeptools.smk"

# G4 transcription analysis rules
include: "workflow/rules/g4_txn.smk"
include: "workflow/rules/g4_promoters.smk"
include: "workflow/rules/g4_repair.smk"
