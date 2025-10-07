########################################################

# mapping to regions of interest
########################################################
rule map_xr_ds:
    input:
        bed=lambda wildcards: f"results/g4_miner/oqs_{"_".join(wildcards.region.split('_')[1:-2])}.bed",
        plus="resources/samples/{sample}_plus.bed",
        minus="resources/samples/{sample}_minus.bed",
        plus_sim="resources/samples/{sample}_sim_plus.bed",
        minus_sim="resources/samples/{sample}_sim_minus.bed"
    output:
        mapped="results/master/{region}/{sample}_mapped.csv",
        mapped_sim="results/master/{region}/{sample}_sim_mapped.csv"
    log:
        "logs/mapping/mapping_{sample}_{region}.log"
    benchmark:
        "logs/mapping/mapping_{sample}_{region}.benchmark.log"
    conda:
        "../../workflow/envs/mapping.yaml"
    threads: 8
    resources:
        mem_mb=64000
    wildcard_constraints:
        region="(?!tss).*",
    shell:
        """
        python3 workflow/scripts/mapping.py -i resources/samples/{wildcards.sample} \
        -r {wildcards.region} > {log} 2>&1
        """

rule map_xr_ds_tss:
    input:
        bed="resources/ref_genomes/hg38/hg38_proteinCoding_genes.bed",
        plus="resources/samples/{sample}_plus.bed",
        minus="resources/samples/{sample}_minus.bed",
        plus_sim="resources/samples/{sample}_sim_plus.bed",
        minus_sim="resources/samples/{sample}_sim_minus.bed"
    output:
        mapped="results/master/{region}/{sample}_mapped.csv",
        mapped_sim="results/master/{region}/{sample}_sim_mapped.csv"
    log:
        "logs/mapping/mapping_{sample}_{region}.log"
    benchmark:
        "logs/mapping/mapping_{sample}_{region}.benchmark.log"
    conda:
        "../../workflow/envs/mapping.yaml"
    threads: 8
    resources:
        mem_mb=64000
    wildcard_constraints:
        region=".*tss.*"
    shell:
        """
        python3 workflow/scripts/mapping.py -i resources/samples/{wildcards.sample} \
        -r {wildcards.region} -t True > {log} 2>&1
        """

########################################################

# combine mapped results
########################################################
rule map_xr_ds_combine:
    input:
        mapped_paths=lambda wildcards: expand(
            "results/master/{region}/{sample}_mapped.csv",
            region=[wildcards.region],
            sample=list(config["samples"].keys()),
        ),
        sim_paths=lambda wildcards: expand(
            "results/master/{region}/{sample}_sim_mapped.csv",
            region=[wildcards.region],
            sample=list(config["samples"].keys()),
        ),
    output:
        mapped_all="results/master/{region}/mapped_all.csv"
    params:
        config_file="config/config.yaml"
    log:
        "logs/mapping_combine/mapping_combine_{region}.log"
    benchmark:
        "logs/mapping_combine/mapping_combine_{region}.benchmark.log"
    conda:
        "../../workflow/envs/mapping.yaml"
    threads: 8
    shell:
        """
        python3 workflow/scripts/mapping_combine.py -r {wildcards.region} -p {input.mapped_paths} \
        -s {input.sim_paths} -c {params.config_file} > {log} 2>&1
        """