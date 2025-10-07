########################################################

# deepTools: map bigWig signal on regions and plot
########################################################

import os

rule index_bam:
    input:
        bam="resources/atac/{sample}.bam"
    output:
        bai="resources/atac/{sample}.bam.bai"
    log:
        "logs/deeptools/index_bam_{sample}.log"
    benchmark:
        "logs/deeptools/index_bam_{sample}.benchmark.log"
    conda:
        "../../workflow/envs/samtools.yaml"
    threads: 16
    shell:
        r"""
        samtools index -@ {threads} {input.bam} > {log} 2>&1
        """

rule bam_to_bigwig:
    input:
        bam="resources/atac/{sample}.bam",
        bai="resources/atac/{sample}.bam.bai",
    output:
        bw="resources/atac/{sample}.bw",
    log:
        "logs/deeptools/bam_to_bigwig_{sample}.log"
    benchmark:
        "logs/deeptools/bam_to_bigwig_{sample}.benchmark.log"
    conda:
        "../../workflow/envs/deeptools.yaml"
    threads: 8
    resources:
        mem_mb=16000
    shell:
        r"""
        bamCoverage \
            -b {input.bam} \
            -o {output.bw} \
            --normalizeUsing CPM \
            --numberOfProcessors {threads} > {log} 2>&1
        """

rule compute_matrix_regions:
    input:
        bw=[
            "resources/atac/HelanoUV_R1_ATAC_TAAGGCGA-TAGATCGC_S7.shifted.bw",
            "resources/atac/Hela_1m_R3_ATAC_AGGCAGAA-CTCTCTAT_S3.shifted.bw",
            "resources/atac/Hela_15m_R3_ATAC_AGGCAGAA-TATCCTCT_S2.shifted.bw",
            "resources/atac/Hela_30m_R3_ATAC.shifted.bw",
            "resources/atac/Hela_1h_R3_ATAC.shifted.bw",
            "resources/atac/Hela_4h_R3_ATAC_TCCTGAGC-CTCTCTAT_S4.shifted.bw",
            "resources/atac/Hela_8h_R3_ATAC_TCCTGAGC-TATCCTCT_S6.shifted.bw",
        ],
        regions="results/g4_miner/{region}.bed",
    output:
        matrix="results/deeptools/{region}_atac/matrix.gz",
    log:
        "logs/deeptools/compute_matrix_{region}.log"
    benchmark:
        "logs/deeptools/compute_matrix_{region}.benchmark.log"
    conda:
        "../../workflow/envs/deeptools.yaml"
    threads: 16
    shell:
        r"""
        computeMatrix reference-point \
            --referencePoint center \
            --scoreFileName {input.bw} \
            --regionsFileName {input.regions} \
            --beforeRegionStartLength 2500 \
            --afterRegionStartLength 2500 \
            --skipZeros \
            --numberOfProcessors {threads} \
            --outFileName {output.matrix} > {log} 2>&1
        """

rule plot_profile_regions:
    input:
        matrix="results/deeptools/{region}_atac/matrix.gz"
    output:
        profile="results/deeptools/{region}_atac/profile.pdf",
        heatmap="results/deeptools/{region}_atac/heatmap.pdf"
    log:
        "logs/deeptools/plot_{region}_atac.log"
    benchmark:
        "logs/deeptools/plot_{region}_atac.benchmark.log"
    conda:
        "../../workflow/envs/deeptools.yaml"
    threads: 16
    shell:
        r"""
        plotProfile -m {input.matrix} \
        --samplesLabel "no UV" "1m" "15m" "30m" "1h" "4h" "8h" \
        --yAxisLabel "Normalized signal" \
        -out {output.profile} > {log} 2>&1

        plotHeatmap -m {input.matrix} \
        --samplesLabel "no UV" "1m" "15m" "30m" "1h" "4h" "8h" \
        --xAxisLabel "Distance from G4 center (bp)" \
        --yAxisLabel "Normalized signal" \
        -out {output.heatmap} >> {log} 2>&1
        """


