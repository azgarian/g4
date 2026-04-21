

rule download_bed:
    output:
        bed="results/oqs/{oqs_name}.bed.gz",
        bed_unzip=temp("results/oqs/{oqs_name}.bed"),
    params:
        url=lambda wc: config["oqs"][wc.oqs_name]["bed_url"],
    log:
        "logs/prepare_oqs/download_bed/{oqs_name}.log"
    benchmark:
        "logs/prepare_oqs/download_bed/{oqs_name}.benchmark.log"
    shell:
        "(mkdir -p results/oqs && wget -O {output.bed} {params.url} && gunzip -k {output.bed}) > {log} 2>&1"

rule prep_oqs_beds:
    input:
        bed="results/oqs/{oqs_name}.bed",
    output:
        bed="results/oqs/{oqs_name}_prep.bed",
    params:
        strand=lambda wc: "+" if "minus" in wc.oqs_name else "-",
    log:
        "logs/prepare_oqs/prep_oqs_beds/{oqs_name}.log"
    benchmark:
        "logs/prepare_oqs/prep_oqs_beds/{oqs_name}.benchmark.log"
    shell:
        """
        (echo "`date -R`: Preparing OQS beds..." &&
        awk 'BEGIN{{OFS="\\t"}} {{print $1, $2, $3, "peak_" NR, $4, "{params.strand}"}}' {input.bed} > {output.bed} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule liftOver_prep_oqs_bed:
    input:
        bed="results/oqs/{oqs_name}_prep.bed",
        chain=rules.download_liftOver_chain.output.chain_unzip,
    output:
        lifted_bed="results/oqs/{oqs_name}_prep_lifted2hg38.bed",
        unmapped_bed="results/oqs/{oqs_name}_prep_lifted_unmapped.bed",
    log:
        "logs/prepare_oqs/liftOver_prep_oqs_bed/{oqs_name}.log"
    benchmark:
        "logs/prepare_oqs/liftOver_prep_oqs_bed/{oqs_name}.benchmark.log"
    conda:
        "../envs/liftOver.yaml"
    shell:
        "liftOver {input.bed} {input.chain} {output.lifted_bed} {output.unmapped_bed} > {log} 2>&1"


rule liftOver_prep_oqs_bed_blacklist_filtered:
    input:
        bed="results/oqs/{oqs_name}_prep_lifted2hg38.bed",
        blacklist=rules.hg38_blacklist_download.output.bed,
    output:
        bed="results/oqs/{oqs_name}_prep_lifted2hg38_filtered.bed",
    log:
        "logs/prepare_oqs/liftOver_prep_oqs_bed_blacklist_filtered/{oqs_name}.log"
    benchmark:
        "logs/prepare_oqs/liftOver_prep_oqs_bed_blacklist_filtered/{oqs_name}.benchmark.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        (echo "`date -R`: Removing hg38 ENCODE blacklist overlaps from lifted prepared OQS BED for {wildcards.oqs_name}..." &&
        bedtools intersect -nonamecheck -v -a {input.bed} -b {input.blacklist} | sort -k1,1 -k2,2n > {output.bed} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; rm -f {output.bed}; exit 1; }} ) > {log} 2>&1
        """


rule oqs_strand_combine:
    input:
        beds=lambda wc: [
            f"results/oqs/minus_{wc.condition}_prep_lifted2hg38_filtered.bed",
            f"results/oqs/plus_{wc.condition}_prep_lifted2hg38_filtered.bed",
        ],
    output:
        bed="results/oqs/oqs_{condition}.bed",
    log:
        "logs/prepare_oqs/oqs_strand_combine/{condition}.log"
    benchmark:
        "logs/prepare_oqs/oqs_strand_combine/{condition}.benchmark.log"
    shell:
        """
        (echo "`date -R`: Combining lifted prepared OQS tracks for condition={wildcards.condition}..." &&
        cat {input.beds} | sort -k1,1 -k2,2n > {output.bed} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; rm -f {output.bed}; exit 1; }} ) > {log} 2>&1
        """


rule prepare_oqs_condition_structured:
    input:
        bed="results/oqs/oqs_{condition}.bed",
        ref=rules.genome_download.output.fa,
    output:
        tsv="results/oqs/oqs_{condition}_prepared.tsv",
        source_bed="results/oqs/oqs_{condition}_prepared.bed",
        center_bed="results/oqs/oqs_{condition}_prepared_center.bed",
    wildcard_constraints:
        condition="K|PDS",
    log:
        "logs/prepare_oqs/prepare_oqs_condition_structured/{condition}.log"
    benchmark:
        "logs/prepare_oqs/prepare_oqs_condition_structured/{condition}.benchmark.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        (echo "`date -R`: Preparing structure-aware OQS loci for condition={wildcards.condition}..." &&
        python3 workflow/scripts/prepare_external_g4_dataset.py \
          --input-bed {input.bed} \
          --ref {input.ref} \
          --source-dataset oqs_{wildcards.condition} \
          --mode fixed \
          --signal-column-index 4 \
          --name-column-index 3 \
          --strand-column-index 5 \
          --out-tsv {output.tsv} \
          --out-source-bed {output.source_bed} \
          --out-center-bed {output.center_bed} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """
