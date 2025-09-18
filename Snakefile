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
        # Use wildcards to reference all per-sample final outputs
        expand("results/{sample}_number_of_lines.txt", sample=config["samples"]),
        expand("results/{sample}_length.txt", sample=config["samples"])

# You can either include a rule file/s or define the rule/s directly here
# For larger projects, prefer modular rule files under workflow/rules/

## Include your rule files here
## include: "workflow/rules/example_rules.smk"

## Example rule - replace with your actual rules
# Produce a simple per-sample artifact (placeholder for your first step)
rule example_rule:
    output:
        "results/{sample}_raw.txt"
    log:
        "logs/example_rule_{sample}.log"
    benchmark:
        "logs/example_rule_{sample}.benchmark.log"
    shell:
        """
        (echo 'Hello from Snakemake on HPC! \nCurrent sample is {wildcards.sample}' \
        > {output}) &> {log}
        """

# Derive a basic metric from the raw artifact (illustrative second step)
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

# Compute a simple property using a helper in params (demonstrates params usage)
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
