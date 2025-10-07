rule g4_repair_aggregate:
    input:
        csv="results/master/{region}/mapped_all.csv"
    output:
        "results/g4_repair/{region}/{region}_64_normalized.png"
    conda:
        "../../workflow/envs/g4_repair.yaml"
    threads: 24
    resources:
        mem_mb=250000
    log:
        "logs/g4_repair/{region}.log"
    benchmark:
        "logs/g4_repair/{region}.benchmark.log"
    shell:
        r"""
        python3 workflow/scripts/g4_repair.py \
          --name {wildcards.region} \
          --input {input.csv} \
          --outdir results/g4_repair/{wildcards.region} > {log} 2>&1
        """


