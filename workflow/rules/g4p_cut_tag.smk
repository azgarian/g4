
rule align_cut_tag:
    input:
        r1 = "resources/samples/g4cuttag/g4_hela_cut_tag_1.fastq",
        r2 = "resources/samples/g4cuttag/g4_hela_cut_tag_2.fastq",
        bowtie2 = "resources/ref_genomes/hg38/Bowtie2/genome_hg38.1.bt2",
    output:
        temp("results/g4cuttag/g4_hela_cut_tag_hg38.sam"),
    params:
        ref_genome = "resources/ref_genomes/hg38/Bowtie2/genome_hg38",
        extra = "--no-unal --no-discordant --no-mixed --seed 1 --reorder",
    threads: 32 
    log:
        "logs/g4p_cut_tag/align_cut_tag.log",
    benchmark:
        "logs/g4p_cut_tag/align_cut_tag.benchmark.txt",
    conda:
        "../envs/bowtie2.yaml"
    shell:  
        """
        mkdir -p results/g4cuttag/

        (echo "`date -R`: Aligning..." &&
        bowtie2 --threads {threads} \
        -x {params.ref_genome} \
        -1 {input.r1} -2 {input.r2} -S {output} \
        {params.extra} && 
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

rule samtools:
    input:
        rules.align_cut_tag.output,
    output:
        bam = temp("results/g4cuttag/g4_hela_cut_tag_hg38.bam"),
        sort = temp("results/g4cuttag/g4_hela_cut_tag_hg38_sorted.bam"),
    params:
        tmpdir = "results/g4cuttag/",
    threads: 16
    log:
        "logs/g4p_cut_tag/samtools.log",
    benchmark:
        "logs/g4p_cut_tag/samtools.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Removing Low Quality..." &&
        samtools view -bh -q 20 -o {output.bam} {input} && 
        samtools sort {output.bam} -o {output.sort} \
        -@ {threads} -T {params.tmpdir} &&
        samtools index {output.sort} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

rule filter_blacklist_g4cuttag:
    input:
        bam=rules.samtools.output.sort,
        blacklist=rules.hg38_blacklist_download.output.bed,
    output:
        tmp_bam=temp("results/g4cuttag/g4_hela_cut_tag_hg38_blacklist_filtered.tmp.bam"),
        bam="results/g4cuttag/g4_hela_cut_tag_hg38_blacklist_filtered.bam",
    log:
        "logs/g4p_cut_tag/filter_blacklist_g4cuttag.log",
    benchmark:
        "logs/g4p_cut_tag/filter_blacklist_g4cuttag.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    threads: 16
    shell:
        """
        (echo "date -R: Filtering blacklist regions..." &&
        bedtools intersect -nonamecheck -v -abam {input.bam} -b {input.blacklist} > {output.tmp_bam} &&
        samtools sort -@ {threads} -O bam -o {output.bam} {output.tmp_bam} &&
        samtools index -@ {threads} {output.bam} &&
        echo "date -R: Success! Filtering blacklist regions is done." ||
        {{ echo "date -R: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

rule macs2_hela_cut_tag:
    input:
        bam = rules.filter_blacklist_g4cuttag.output.bam,
    output:
        np = "results/g4cuttag/g4_hela_peaks.narrowPeak",
        summits = "results/g4cuttag/g4_hela_summits.bed",
    params:
        name = "g4_hela", 
        outdir = "results/g4cuttag",
    log:
        "logs/g4p_cut_tag/macs2_hela_cut_tag.log",
    benchmark:
        "logs/g4p_cut_tag/macs2_hela_cut_tag.benchmark.txt",
    conda:
        "../envs/macs2.yaml"
    shell:  
        """
        (echo "`date -R`: Peak calling (narrowPeak)..." && 
        macs2 callpeak \
        -q 0.0001 \
        -f BAMPE \
        -t {input.bam} \
        -n {params.name} \
        --outdir {params.outdir} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

rule prepare_g4cuttag_structured:
    input:
        bed=rules.macs2_hela_cut_tag.output.np,
        ref=rules.genome_download.output.fa,
    output:
        tsv="results/g4cuttag/g4_hela_peaks_prepared.tsv",
        source_bed="results/g4cuttag/g4_hela_peaks_prepared.bed",
        center_bed="results/g4cuttag/g4_hela_peaks_prepared_center.bed",
    log:
        "logs/g4p_cut_tag/prepare_g4cuttag_structured.log",
    benchmark:
        "logs/g4p_cut_tag/prepare_g4cuttag_structured.benchmark.txt",
    conda:
        "../envs/python.yaml"
    shell:
        """
        (echo "`date -R`: Preparing structure-aware g4cuttag loci from narrowPeak..." &&
        python3 workflow/scripts/prepare_external_g4_dataset.py \
          --input-bed {input.bed} \
          --ref {input.ref} \
          --source-dataset g4cuttag \
          --mode strandless \
          --signal-column-index 6 \
          --name-column-index 3 \
          --out-tsv {output.tsv} \
          --out-source-bed {output.source_bed} \
          --out-center-bed {output.center_bed} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """
