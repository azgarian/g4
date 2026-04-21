
rule align:
    input:
        r1 = "resources/samples/g4chip/g4_hela_1.fastq",
        r2 = "resources/samples/g4chip/g4_hela_2.fastq",
        bowtie2 = "resources/ref_genomes/hg38/Bowtie2/genome_hg38.1.bt2",
    output:
        sam = temp("results/g4chip/g4_hela_hg38.sam"),
    params:
        ref_genome = "resources/ref_genomes/hg38/Bowtie2/genome_hg38",
        extra = "--sensitive-local --no-unal --no-discordant --no-mixed \
            --seed 1 --reorder",
    threads: 32 
    log:
        "logs/g4p_chip_seq/align.log",
    benchmark:
        "logs/g4p_chip_seq/align.benchmark.txt",
    conda:
        "../envs/bowtie2.yaml"
    shell:  
        """
        mkdir -p results/g4chip/

        (echo "`date -R`: Aligning..." &&
        bowtie2 --threads {threads} \
        -x {params.ref_genome} \
        -1 {input.r1} -2 {input.r2} -S {output} \
        {params.extra} && 
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

rule rmdup:
    input:
        rules.align.output,
    output:
        bam = temp("results/g4chip/g4_hela_hg38.bam"),
        sort = temp("results/g4chip/g4_hela_hg38_sorted.bam"),
        rmdup = temp("results/g4chip/g4_hela_hg38_rmdup.bam"),
    params:
        tmpdir = "results/g4chip/",
    threads: 16
    log:
        "logs/g4p_chip_seq/rmdup.log",
    benchmark:
        "logs/g4p_chip_seq/rmdup.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Removing Duplicates..." &&
        samtools view -bh -q 20 -o {output.bam} {input} && 
        samtools sort {output.bam} -o {output.sort} \
        -@ {threads} -T {params.tmpdir} &&
        samtools index {output.sort} &&
        samtools rmdup {output.sort} {output.rmdup} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

rule filter_blacklist_g4chip:
    input:
        bam=rules.rmdup.output.rmdup,
        blacklist=rules.hg38_blacklist_download.output.bed,
    output:
        tmp_bam=temp("results/g4chip/g4_hela_hg38_blacklist_filtered.tmp.bam"),
        bam="results/g4chip/g4_hela_hg38_blacklist_filtered.bam",
    log:
        "logs/g4p_chip_seq/filter_blacklist_g4chip.log",
    benchmark:
        "logs/g4p_chip_seq/filter_blacklist_g4chip.benchmark.txt",
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

rule macs2_hela:
    input:
        bam = rules.filter_blacklist_g4chip.output.bam,
    output:
        np = "results/g4chip/g4_hela_peaks.narrowPeak",
        summits = "results/g4chip/g4_hela_summits.bed",
    params:
        name = "g4_hela", 
        outdir = "results/g4chip",
    log:
        "logs/g4p_chip_seq/macs2.log",
    benchmark:
        "logs/g4p_chip_seq/macs2.benchmark.txt",
    conda:
        "../envs/macs2.yaml"
    shell:  
        """
        (echo "`date -R`: Peak calling (narrowPeak)..." && 
        macs2 callpeak \
        -q 0.001 \
        --keep-dup 1 \
        -f BAMPE \
        -t {input.bam} \
        -n {params.name} \
        --outdir {params.outdir} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

rule prepare_g4chip_structured:
    input:
        bed=rules.macs2_hela.output.np,
        ref=rules.genome_download.output.fa,
    output:
        tsv="results/g4chip/g4_hela_peaks_prepared.tsv",
        source_bed="results/g4chip/g4_hela_peaks_prepared.bed",
        center_bed="results/g4chip/g4_hela_peaks_prepared_center.bed",
    log:
        "logs/g4p_chip_seq/prepare_g4chip_structured.log",
    benchmark:
        "logs/g4p_chip_seq/prepare_g4chip_structured.benchmark.txt",
    conda:
        "../envs/python.yaml"
    shell:
        """
        (echo "`date -R`: Preparing structure-aware g4chip loci from narrowPeak..." &&
        python3 workflow/scripts/prepare_external_g4_dataset.py \
          --input-bed {input.bed} \
          --ref {input.ref} \
          --source-dataset g4chip \
          --mode strandless \
          --signal-column-index 6 \
          --name-column-index 3 \
          --out-tsv {output.tsv} \
          --out-source-bed {output.source_bed} \
          --out-center-bed {output.center_bed} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """
