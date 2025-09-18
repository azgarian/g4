rule example_rule:
    output:
        "results/{sample}_raw.txt"
    log:
        "logs/example_rule_{sample}.log"
    benchmark:
        "logs/example_rule_{sample}.benchmark.log"
    shell:
        """
        (echo 'Hello from Snakemake on HPC!\nCurrent sample is {wildcards.sample}' \
        > {output}) &> {log}
        """

rule example_rule2:
    input:
        "results/{sample}_raw.txt"
    output:
        "results/{sample}_number_of_lines.txt"
    log:
        "logs/example_rule2_{sample}.log"
    benchmark:
        "logs/example_rule2_{sample}.benchmark.log"
    shell:
        "(wc -l {input} > {output}) &> {log}"

rule example_rule3:
    input:
        "results/{sample}_raw.txt"
    output:
        "results/{sample}_length.txt"
    log:
        "logs/example_rule3_{sample}.log"
    benchmark:
        "logs/example_rule3_{sample}.benchmark.log"
    params:
        is_odd = lambda wildcards, input: is_odd(str(input[0]))
    shell:
        "(echo 'Length of {input} is {params.is_odd}' > {output}) &> {log}"
