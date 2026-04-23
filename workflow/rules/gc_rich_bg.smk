GC_RICH_BG_TILE_SIZE = 150
GC_RICH_BG_GC_THRESHOLD = 0.28
GC_RICH_BG_SAMPLE_SIZE = 40000
GC_RICH_BG_SEED = 42
GC_RICH_BG_MIN_PREPARED = 7000
GC_RICH_BG_PROMOTER_SAMPLE_SIZE = 1000


rule gc_rich_bg_select:
    input:
        ref=rules.genome_download.output.fa,
        fai=rules.genome_indexing.output,
        blacklist=rules.hg38_blacklist_download.output.bed,
        oqs_k_bed="results/oqs/oqs_K_prepared.bed",
        oqs_pds_bed="results/oqs/oqs_PDS_prepared.bed",
        g4_bed="results/g4chip_g4cuttag/g4_hela_chip_cuttag_merged.bed",
    output:
        sampled_bed="results/gc_rich_bg/gc_rich_bg_sampled.bed",
        filtered_bed="results/gc_rich_bg/gc_rich_bg_sampled_filtered.bed",
        summary_tsv="results/gc_rich_bg/gc_rich_bg_selection_summary.tsv",
    log:
        "logs/gc_rich_bg/gc_rich_bg_select.log",
    benchmark:
        "logs/gc_rich_bg/gc_rich_bg_select.benchmark.txt",
    conda:
        "../envs/python.yaml"
    shell:
        """
        (echo "`date -R`: Sampling GC-rich background loci and filtering OQS/G4 overlaps..." &&
        python3 workflow/scripts/sample_gc_rich_background.py \
          --ref {input.ref} \
          --genome-fai {input.fai} \
          --blacklist-bed {input.blacklist} \
          --oqs-bed {input.oqs_k_bed} \
          --oqs-bed {input.oqs_pds_bed} \
          --g4-bed {input.g4_bed} \
          --tile-size {GC_RICH_BG_TILE_SIZE} \
          --gc-threshold {GC_RICH_BG_GC_THRESHOLD} \
          --sample-size {GC_RICH_BG_SAMPLE_SIZE} \
          --seed {GC_RICH_BG_SEED} \
          --out-sampled-bed {output.sampled_bed} \
          --out-filtered-bed {output.filtered_bed} \
          --out-summary-tsv {output.summary_tsv} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule gc_rich_bg_prepare_structured:
    input:
        bed=rules.gc_rich_bg_select.output.filtered_bed,
        ref=rules.genome_download.output.fa,
    output:
        tsv="results/gc_rich_bg/gc_rich_bg_prepared.tsv",
        source_bed="results/gc_rich_bg/gc_rich_bg_prepared.bed",
        center_bed="results/gc_rich_bg/gc_rich_bg_prepared_center.bed",
    log:
        "logs/gc_rich_bg/gc_rich_bg_prepare_structured.log",
    benchmark:
        "logs/gc_rich_bg/gc_rich_bg_prepare_structured.benchmark.txt",
    conda:
        "../envs/python.yaml"
    shell:
        """
        (echo "`date -R`: Preparing structure-aware GC-rich background loci..." &&
        python3 workflow/scripts/prepare_external_g4_dataset.py \
          --input-bed {input.bed} \
          --ref {input.ref} \
          --source-dataset gc_rich_bg \
          --mode strandless \
          --signal-column-index 4 \
          --name-column-index 3 \
          --out-tsv {output.tsv} \
          --out-source-bed {output.source_bed} \
          --out-center-bed {output.center_bed} &&
        python3 workflow/scripts/assert_min_tsv_rows.py \
          --tsv {output.tsv} \
          --min-rows {GC_RICH_BG_MIN_PREPARED} \
          --label gc_rich_bg &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule gc_rich_bg_promoter_select:
    input:
        promoter_bed="results/g4_tss/canonical_tss_windows_1kb.bed",
        ref=rules.genome_download.output.fa,
        oqs_k_bed="results/oqs/oqs_K_prepared.bed",
        oqs_pds_bed="results/oqs/oqs_PDS_prepared.bed",
        g4_bed="results/g4chip_g4cuttag/g4_hela_chip_cuttag_merged.bed",
    output:
        sampled_bed="results/gc_rich_bg_promoter/gc_rich_bg_promoter_sampled.bed",
        summary_tsv="results/gc_rich_bg_promoter/gc_rich_bg_promoter_selection_summary.tsv",
    log:
        "logs/gc_rich_bg_promoter/gc_rich_bg_promoter_select.log",
    benchmark:
        "logs/gc_rich_bg_promoter/gc_rich_bg_promoter_select.benchmark.txt",
    conda:
        "../envs/python.yaml"
    shell:
        """
        (echo "`date -R`: Sampling promoter-specific GC-rich background loci..." &&
        python3 workflow/scripts/sample_gc_rich_promoter_background.py \
          --promoter-bed {input.promoter_bed} \
          --ref {input.ref} \
          --oqs-bed {input.oqs_k_bed} \
          --oqs-bed {input.oqs_pds_bed} \
          --g4-bed {input.g4_bed} \
          --gc-threshold {GC_RICH_BG_GC_THRESHOLD} \
          --sample-size {GC_RICH_BG_PROMOTER_SAMPLE_SIZE} \
          --seed {GC_RICH_BG_SEED} \
          --out-sampled-bed {output.sampled_bed} \
          --out-summary-tsv {output.summary_tsv} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """
