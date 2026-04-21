rule merge_g4chip_g4cuttag:
    input:
        g4chip=rules.prepare_g4chip_structured.output.source_bed,
        g4cuttag=rules.prepare_g4cuttag_structured.output.source_bed,
    output:
        bed="results/g4chip_g4cuttag/g4_hela_chip_cuttag_merged.bed",
    log:
        "logs/g4chip_g4cuttag/merge_g4chip_g4cuttag.log",
    benchmark:
        "logs/g4chip_g4cuttag/merge_g4chip_g4cuttag.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        mkdir -p results/g4chip_g4cuttag logs/g4chip_g4cuttag

        (echo "`date -R`: Merging g4chip and g4cuttag intervals..." &&
        cat {input.g4chip} {input.g4cuttag} |
        sort -k1,1 -k2,2n |
        bedtools merge -i - > {output.bed} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """
