########################################################

# G4–Transcription (bioframe) analysis rule
########################################################

rule g4_transcription_bioframe:
    input:
        gtf="resources/ref_genomes/hg38/gencode.v49.annotation.gtf.gz",
        g4=[
            "results/g4_miner/oqs_0.bed",
            "results/g4_miner/oqs_12.bed",
            "results/g4_miner/oqs_30.bed",
            "results/g4_miner/oqs_60.bed",
        ],
        rna=[
            "resources/rna-seq/SU_100/quant.sf",
            "resources/rna-seq/SU_112/quant.sf",
            "resources/rna-seq/SU_130/quant.sf",
            "resources/rna-seq/SU_160/quant.sf",
        ]
    output:
        outdir=directory("results/g4_txn"),
        g4_promoter_counts="results/g4_txn/g4_promoter_counts.tsv"
    params:
        assembly="hg38",
        rna_names=["T0", "T12", "T30", "T60"],
        rna_names_str="T0 T12 T30 T60",
        prom_up=2000,
        prom_down=500,
        gene_id_col="gene_id",
        gene_name_col="gene_name"
    log:
        "logs/g4_txn/g4_txn.log"
    benchmark:
        "logs/g4_txn/g4_txn.benchmark.log"
    conda:
        "../../workflow/envs/g4_txn.yaml"
    threads: 4
    resources:
        mem_mb=16000
    shell:
        r"""
        mkdir -p {output.outdir} logs/g4_txn
        python3 workflow/scripts/g4_transcription_bioframe.py \
          --gtf {input.gtf} \
          --assembly {params.assembly} \
          --g4 {input.g4} \
          --g4-names {params.rna_names_str} \
          --rna {input.rna} \
          --rna-sample-names {params.rna_names} \
          --tmm --prefer-counts \
          --prom-up {params.prom_up} --prom-down {params.prom_down} \
          --gene-id-col {params.gene_id_col} --gene-name-col {params.gene_name_col} \
          --outdir {output.outdir} > {log} 2>&1
        """


