
rule g4_promoter_annotate:
    input:
        g4="results/g4_miner/{region}.bed",
        promoters="results/g4_txn/promoters.bed.gz"
    output:
        info="results/g4_miner/{region}.promoter_info.tsv"
    conda:
        "../../workflow/envs/g4_txn.yaml"
    threads: 2
    resources:
        mem_mb=4000
    log:
        "logs/g4_txn/g4_promoter_annotate_{region}.log"
    benchmark:
        "logs/g4_txn/g4_promoter_annotate_{region}.benchmark.log"
    shell:
        r"""
        python3 workflow/scripts/g4_promoter_annotate.py \
          --g4 {input.g4} \
          --promoters {input.promoters} \
          --out {output.info} > {log} 2>&1
        """


