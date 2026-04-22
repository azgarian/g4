import re
from itertools import combinations


RNASEQ_CFG = config["rna_seq"]
RNASEQ_SALMON_DIR = str(RNASEQ_CFG["salmon_dir"]).rstrip("/")
RNASEQ_SAMPLE_REGEX = re.compile(RNASEQ_CFG["sample_regex"])
RNASEQ_TIME_ORDER = [str(tp).strip() for tp in RNASEQ_CFG["time_order"]]
RNASEQ_ANNOTATION_DIR = str(RNASEQ_CFG["annotation_dir"]).rstrip("/")
RNASEQ_V35_GTF = f"{RNASEQ_ANNOTATION_DIR}/gencode.v35.annotation.gtf"
RNASEQ_RESULTS_DIR = "results/rnaseq"
RNASEQ_METADATA_DIR = f"{RNASEQ_RESULTS_DIR}/metadata"
RNASEQ_MATRICES_DIR = f"{RNASEQ_RESULTS_DIR}/matrices"
RNASEQ_DESEQ2_DIR = f"{RNASEQ_RESULTS_DIR}/deseq2"
RNASEQ_QC_DIR = f"{RNASEQ_RESULTS_DIR}/qc"
RNASEQ_PAIRWISE_CONTRASTS = list(combinations(RNASEQ_TIME_ORDER, 2))
RNASEQ_PAIRWISE_LABELS = [
    f"{numerator}_vs_{denominator}" for numerator, denominator in RNASEQ_PAIRWISE_CONTRASTS
]
RNASEQ_PAIRWISE_DIR = f"{RNASEQ_DESEQ2_DIR}/pairwise"


def discover_rnaseq_samples():
    samples = []
    root = Path(RNASEQ_SALMON_DIR)
    for entry in root.glob("SU_*"):
        if not entry.is_dir():
            continue
        match = RNASEQ_SAMPLE_REGEX.match(entry.name)
        if match is None:
            continue
        rep = int(match.group("rep"))
        time_label = str(int(match.group("time")))
        if time_label not in RNASEQ_TIME_ORDER:
            continue
        samples.append((RNASEQ_TIME_ORDER.index(time_label), rep, entry.name))
    samples.sort()
    return [sample for _, _, sample in samples]


RNASEQ_SAMPLES = discover_rnaseq_samples()


rule rnaseq_v35_gtf_download:
    output:
        gtf=RNASEQ_V35_GTF,
    params:
        url=RNASEQ_CFG["gencode_v35_gtf_url"],
        outdir=RNASEQ_ANNOTATION_DIR,
    log:
        "logs/rnaseq/rnaseq_v35_gtf_download.log",
    benchmark:
        "logs/rnaseq/rnaseq_v35_gtf_download.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Downloading Gencode v35 GTF for RNA-seq..." &&
        mkdir -p {params.outdir} "$(dirname {log})" &&
        wget -q {params.url} -O {output.gtf}.gz &&
        gunzip -f {output.gtf}.gz &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; rm -f {output.gtf} {output.gtf}.gz; exit 1; }} ) > {log} 2>&1
        """


rule rnaseq_sample_manifest:
    input:
        quant=expand(f"{RNASEQ_SALMON_DIR}/{{sample}}/quant.sf", sample=RNASEQ_SAMPLES),
    output:
        tsv=f"{RNASEQ_METADATA_DIR}/sample_manifest.tsv",
    params:
        salmon_dir=RNASEQ_SALMON_DIR,
        sample_regex=RNASEQ_CFG["sample_regex"],
        time_order=" ".join(RNASEQ_TIME_ORDER),
    log:
        "logs/rnaseq/rnaseq_sample_manifest.log",
    benchmark:
        "logs/rnaseq/rnaseq_sample_manifest.benchmark.txt",
    conda:
        "../envs/python.yaml"
    shell:
        """
        (echo "`date -R`: Building RNA-seq sample manifest from Salmon outputs..." &&
        mkdir -p "$(dirname {output.tsv})" "$(dirname {log})" &&
        python3 workflow/scripts/rnaseq_sample_manifest.py \
          --salmon-dir {params.salmon_dir} \
          --sample-regex '{params.sample_regex}' \
          --time-order {params.time_order} \
          --out-tsv {output.tsv} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; rm -f {output.tsv}; exit 1; }} ) > {log} 2>&1
        """


rule rnaseq_tx2gene:
    input:
        gtf=rules.rnaseq_v35_gtf_download.output.gtf,
    output:
        tx2gene=f"{RNASEQ_METADATA_DIR}/tx2gene.tsv",
        gene_annotation=f"{RNASEQ_METADATA_DIR}/gene_annotation.tsv",
    log:
        "logs/rnaseq/rnaseq_tx2gene.log",
    benchmark:
        "logs/rnaseq/rnaseq_tx2gene.benchmark.txt",
    conda:
        "../envs/python.yaml"
    shell:
        """
        (echo "`date -R`: Parsing RNA-seq GTF into tx2gene and gene annotation tables..." &&
        mkdir -p "$(dirname {output.tx2gene})" "$(dirname {log})" &&
        python3 workflow/scripts/rnaseq_tx2gene_from_gtf.py \
          --gtf {input.gtf} \
          --out-tx2gene {output.tx2gene} \
          --out-gene-annotation {output.gene_annotation} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; rm -f {output.tx2gene} {output.gene_annotation}; exit 1; }} ) > {log} 2>&1
        """


rule rnaseq_time_series_deseq2:
    input:
        sample_manifest=rules.rnaseq_sample_manifest.output.tsv,
        tx2gene=rules.rnaseq_tx2gene.output.tx2gene,
        gene_annotation=rules.rnaseq_tx2gene.output.gene_annotation,
    output:
        mapping_summary=f"{RNASEQ_METADATA_DIR}/tximport_mapping_summary.tsv",
        gene_counts=f"{RNASEQ_MATRICES_DIR}/gene_counts.tsv.gz",
        gene_tpm=f"{RNASEQ_MATRICES_DIR}/gene_tpm.tsv.gz",
        gene_lengths=f"{RNASEQ_MATRICES_DIR}/gene_lengths.tsv.gz",
        normalized_counts=f"{RNASEQ_MATRICES_DIR}/normalized_counts.tsv.gz",
        lrt_results=f"{RNASEQ_DESEQ2_DIR}/time_lrt_results.tsv.gz",
        lrt_significant=f"{RNASEQ_DESEQ2_DIR}/time_lrt_significant.tsv.gz",
        pairwise_summary=f"{RNASEQ_PAIRWISE_DIR}/pairwise_significant_gene_counts.tsv",
        pairwise_results=expand(
            f"{RNASEQ_PAIRWISE_DIR}" + "/{contrast}_results.tsv.gz",
            contrast=RNASEQ_PAIRWISE_LABELS,
        ),
        pairwise_significant=expand(
            f"{RNASEQ_PAIRWISE_DIR}" + "/{contrast}_significant.tsv.gz",
            contrast=RNASEQ_PAIRWISE_LABELS,
        ),
        pca_tsv=f"{RNASEQ_QC_DIR}/vst_pca.tsv",
        pca_png=f"{RNASEQ_QC_DIR}/vst_pca.png",
        pca_pdf=f"{RNASEQ_QC_DIR}/vst_pca.pdf",
        sample_distance_tsv=f"{RNASEQ_QC_DIR}/sample_distance.tsv",
        sample_distance_png=f"{RNASEQ_QC_DIR}/sample_distance.png",
        sample_distance_pdf=f"{RNASEQ_QC_DIR}/sample_distance.pdf",
        summary_md=f"{RNASEQ_RESULTS_DIR}/summary.md",
    params:
        alpha=RNASEQ_CFG["alpha"],
        prefilter_min_count=RNASEQ_CFG["prefilter_min_count"],
        prefilter_min_samples=RNASEQ_CFG["prefilter_min_samples"],
        time_order=",".join(RNASEQ_TIME_ORDER),
        pairwise_dir=RNASEQ_PAIRWISE_DIR,
    log:
        "logs/rnaseq/rnaseq_time_series_deseq2.log",
    benchmark:
        "logs/rnaseq/rnaseq_time_series_deseq2.benchmark.txt",
    conda:
        "../envs/rnaseq_deseq2.yaml"
    threads: 8
    shell:
        """
        (echo "`date -R`: Running RNA-seq tximport and DESeq2 time-series analysis..." &&
        mkdir -p "$(dirname {output.summary_md})" "$(dirname {log})" &&
        Rscript workflow/scripts/rnaseq_time_series_deseq2.R \
          --sample-manifest {input.sample_manifest} \
          --tx2gene {input.tx2gene} \
          --gene-annotation {input.gene_annotation} \
          --prefilter-min-count {params.prefilter_min_count} \
          --prefilter-min-samples {params.prefilter_min_samples} \
          --alpha {params.alpha} \
          --time-order {params.time_order} \
          --threads {threads} \
          --out-mapping-summary {output.mapping_summary} \
          --out-gene-counts {output.gene_counts} \
          --out-gene-tpm {output.gene_tpm} \
          --out-gene-lengths {output.gene_lengths} \
          --out-normalized-counts {output.normalized_counts} \
          --out-lrt-results {output.lrt_results} \
          --out-lrt-significant {output.lrt_significant} \
          --out-pairwise-dir {params.pairwise_dir} \
          --out-pairwise-summary {output.pairwise_summary} \
          --out-pca-tsv {output.pca_tsv} \
          --out-pca-png {output.pca_png} \
          --out-pca-pdf {output.pca_pdf} \
          --out-sample-distance-tsv {output.sample_distance_tsv} \
          --out-sample-distance-png {output.sample_distance_png} \
          --out-sample-distance-pdf {output.sample_distance_pdf} \
          --out-summary-md {output.summary_md} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """
