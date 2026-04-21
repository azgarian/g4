rule genome_download:
    output:
        fa="resources/ref_genomes/hg38/genome_hg38.fa",
        info="resources/ref_genomes/hg38/genome_hg38.info",
    params:
        link=config["genome"]["link"],
    log:
        "logs/prepare_genome/genome_download/hg38.log",
    benchmark:
        "logs/prepare_genome/genome_download/hg38.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Downloading hg38 genome..." &&
        wget {params.link} -O {output.fa}.gz &&
        gunzip {output.fa}.gz &&
        echo {params.link} > {output.info} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

rule genome_tar_zip:
    input:
        rules.genome_download.output.fa,
    output:
        "resources/ref_genomes/hg38/GRCh38.tar.gz"
    log:
        "logs/prepare_genome/genome_tar_zip/hg38.log",
    benchmark:
        "logs/prepare_genome/genome_tar_zip/hg38.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Zipping genome..." &&
        tar -czvf {output} {input} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

rule genome_split2chr:
    input:
        rules.genome_download.output.fa,
    output:
        expand("resources/ref_genomes/hg38/chr/chr{chr}.fa", 
            chr=[*range(1,23)]+["X"]),
    log:
        "logs/prepare_genome/genome_split2chr/hg38.log",
    benchmark:
        "logs/prepare_genome/genome_split2chr/hg38.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Splitting chromosomes into separate files..." &&
        mkdir -p resources/ref_genomes/hg38/chr/ &&
        csplit -s -z {input} '/>/' '{{*}}' \
        --prefix resources/ref_genomes/hg38/chr/ &&
        for i in resources/ref_genomes/hg38/chr/* ; do \
        n=$(sed 's/>// ; s/ .*// ; 1q' "$i") ; \
        mv "$i" resources/ref_genomes/hg38/chr/"$n.fa" ; \
        done &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

rule genome_build:
    input:
        rules.genome_download.output.fa,
    output:
        multiext(
        "resources/ref_genomes/hg38/Bowtie2/genome_hg38",
        ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    log: 
        "logs/prepare_genome/genome_build/hg38.log",
    benchmark:
        "logs/prepare_genome/genome_build/hg38.benchmark.txt",
    params:
        extra="", 
        name="resources/ref_genomes/hg38/Bowtie2/genome_hg38",
    conda:
        "../envs/bowtie2.yaml" 
    threads: 
        16
    shell: 
        """
        (echo "`date -R`: Building indexes..." &&
        bowtie2-build --threads {threads} \
        {params.extra} \
        {input} \
        {params.name} &&
        echo "`date -R`: Success! Indexes are build." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """

rule bwa_index:
    input:
        rules.genome_download.output.fa,
    output:
        multiext("resources/ref_genomes/hg38/genome_hg38.fa", 
        ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log: 
        "logs/prepare_genome/bwa_index/hg38.log",
    benchmark:
        "logs/prepare_genome/bwa_index/hg38.benchmark.txt",
    conda:
        "../envs/bwa.yaml"
    threads: 
        16
    shell: 
        """    
        (echo "`date -R`: Building indexes..." &&
        bwa index -a bwtsw {input} &&
        echo "`date -R`: Success! Indexes are build." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """

rule genome_indexing:
    input:
        rules.genome_download.output.fa,
    output:
        "resources/ref_genomes/hg38/genome_hg38.fa.fai",
    log: 
        "logs/prepare_genome/genome_indexing/hg38.log",
    benchmark:
        "logs/prepare_genome/genome_indexing/hg38.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        (echo "`date -R`: Creating fai file..." &&
        samtools faidx {input} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """

rule genome_bed:
    input:
        rules.genome_indexing.output,
    output:
        "resources/ref_genomes/hg38/genome_hg38.bed",
    log: 
        "logs/prepare_genome/genome_bed/hg38.log",
    benchmark:
        "logs/prepare_genome/genome_bed/hg38.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Creating bed file..." &&
        awk '{{print $1"\\t"0"\\t"$2"\\t"".""\\t"".""\\t""."}}' {input} \
        > {output} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """


rule hg38_blacklist_download:
    output:
        bed="resources/ref_genomes/hg38/hg38_encode_blacklist.bed",
    params:
        url=config["genome"]["blacklist_url"],
    log:
        "logs/prepare_genome/hg38_blacklist_download.log",
    benchmark:
        "logs/prepare_genome/hg38_blacklist_download.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Downloading hg38 ENCODE blacklist..." &&
        mkdir -p resources/ref_genomes/hg38 &&
        curl -L {params.url} -o {output.bed}.gz &&
        gunzip -f {output.bed}.gz &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule repeatmasker_download:
    output:
        bed="resources/ref_genomes/hg38/hg38_repeatmasker.bed",
    params:
        url=config["genome"]["repeatmasker_url"],
    log:
        "logs/prepare_genome/repeatmasker_download/hg38.log",
    benchmark:
        "logs/prepare_genome/repeatmasker_download/hg38.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Downloading hg38 RepeatMasker annotations..." &&
        mkdir -p resources/ref_genomes/hg38 &&
        curl -L {params.url} -o {output.bed}.txt.gz &&
        gunzip -c {output.bed}.txt.gz | \
          awk 'BEGIN {{OFS="\\t"}} {{print $6, $7, $8, $11 "|" $12 "|" $13, $2, $10}}' \
          > {output.bed} &&
        rm -f {output.bed}.txt.gz &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; rm -f {output.bed}; exit 1; }} ) > {log} 2>&1
        """


rule mappability_download:
    output:
        bed="resources/ref_genomes/hg38/hg38_k50_umap.bed",
    params:
        url=config["genome"]["mappability_url"],
    log:
        "logs/prepare_genome/mappability_download/hg38.log",
    benchmark:
        "logs/prepare_genome/mappability_download/hg38.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Downloading hg38 Umap k50 mappability regions..." &&
        curl -L {params.url} -o {output.bed}.gz &&
        gunzip {output.bed}.gz &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; rm -f {output.bed}; exit 1; }} ) > {log} 2>&1
        """


rule gencode_gtf_download:
    output:
        gtf = "resources/ref_genomes/hg38/gencode.gtf",
    params:
        url = config["genome"]["gtf_url"],
        outdir = "resources/ref_genomes/hg38",
    log:
        "logs/prepare_genome/gencode_gtf_download/hg38.log",
    benchmark:
        "logs/prepare_genome/gencode_gtf_download/hg38.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Downloading Gencode GTF (hg38)..." &&
        mkdir -p {params.outdir} &&
        wget -q {params.url} -O {output.gtf}.gz &&
        gunzip -f {output.gtf}.gz &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule download_liftOver_chain:
    output:
        chain="resources/ref_genomes/hg38/hg19ToHg38.over.chain.gz",
        chain_unzip=temp("resources/ref_genomes/hg38/hg19ToHg38.over.chain"),
    params:
        url=config["genome"]["liftOver_chain_url"],
    log:
        "logs/prepare_genome/download_liftOver_chain.log"
    benchmark:
        "logs/prepare_genome/download_liftOver_chain.benchmark.log"
    shell:
        """
        (echo "`date -R`: Downloading UCSC chain hg19->hg38..." &&
        mkdir -p resources logs/prepare_genome &&
        curl -L -o {output.chain} {params.url} &&
        gunzip -k {output.chain} > {log} 2>&1
        echo "`date -R`: done!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

# rule g4predict:
#     output:
#         "results/g4predict_downladed.txt",
#     log:
#         "logs/prepare_genome/g4predict.log",
#     benchmark:
#         "logs/prepare_genome/g4predict.benchmark.txt",
#     conda:
#         "../envs/g4_miner.yaml"
#     shell:
#         """
#         (echo "`date -R`: build..." &&
#         mydir=($(pwd)/workflow/scripts/g4predict/) &&
#         if [ -n "${{PYTHONPATH:-}}" ]; then
#             export PYTHONPATH="${{PYTHONPATH}}:$mydir"
#         else
#             export PYTHONPATH="$mydir"
#         fi &&
#         python workflow/scripts/g4predict/setup.py build &&
#         echo "`date -R`: Success!" ||
#         {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1

#         (echo "`date -R`: install..." &&
#         python workflow/scripts/g4predict/setup.py install &&
#         echo "g4predict is ready!" > {output} &&
#         echo "`date -R`: Success!" ||
#         {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1
#         """