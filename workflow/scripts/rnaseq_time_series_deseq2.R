#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(BiocParallel)
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(readr)
  library(tibble)
  library(tximport)
})

parser <- ArgumentParser(description = "Import Salmon RNA-seq quantifications and run DESeq2 LRT over time.")
parser$add_argument("--sample-manifest", required = TRUE, help = "Sample manifest TSV.")
parser$add_argument("--tx2gene", required = TRUE, help = "tx2gene TSV.")
parser$add_argument("--gene-annotation", required = TRUE, help = "Gene annotation TSV.")
parser$add_argument("--prefilter-min-count", required = TRUE, type = "integer", help = "Minimum count per sample for prefilter.")
parser$add_argument("--prefilter-min-samples", required = TRUE, type = "integer", help = "Minimum number of samples meeting the count threshold.")
parser$add_argument("--alpha", required = TRUE, type = "double", help = "Adjusted p-value cutoff.")
parser$add_argument("--time-order", required = TRUE, help = "Comma-separated ordered timepoint labels.")
parser$add_argument("--threads", required = TRUE, type = "integer", help = "Thread count for DESeq2.")
parser$add_argument("--out-mapping-summary", required = TRUE, help = "Output mapping summary TSV.")
parser$add_argument("--out-gene-counts", required = TRUE, help = "Output gene counts TSV.gz.")
parser$add_argument("--out-gene-tpm", required = TRUE, help = "Output gene TPM TSV.gz.")
parser$add_argument("--out-gene-lengths", required = TRUE, help = "Output gene effective length TSV.gz.")
parser$add_argument("--out-normalized-counts", required = TRUE, help = "Output normalized counts TSV.gz.")
parser$add_argument("--out-lrt-results", required = TRUE, help = "Output LRT results TSV.gz.")
parser$add_argument("--out-lrt-significant", required = TRUE, help = "Output significant LRT results TSV.gz.")
parser$add_argument("--out-pairwise-dir", required = TRUE, help = "Output directory for pairwise DESeq2 tables.")
parser$add_argument("--out-pairwise-summary", required = TRUE, help = "Output pairwise summary TSV.")
parser$add_argument("--out-pca-tsv", required = TRUE, help = "Output PCA TSV.")
parser$add_argument("--out-pca-png", required = TRUE, help = "Output PCA PNG.")
parser$add_argument("--out-pca-pdf", required = TRUE, help = "Output PCA PDF.")
parser$add_argument("--out-sample-distance-tsv", required = TRUE, help = "Output sample distance TSV.")
parser$add_argument("--out-sample-distance-png", required = TRUE, help = "Output sample distance heatmap PNG.")
parser$add_argument("--out-sample-distance-pdf", required = TRUE, help = "Output sample distance heatmap PDF.")
parser$add_argument("--out-summary-md", required = TRUE, help = "Output markdown summary.")
args <- parser$parse_args()

strip_version <- function(values) {
  sub("\\.[^.]+$", "", values)
}

ensure_parent <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

write_tsv_with_dirs <- function(df, path) {
  ensure_parent(path)
  readr::write_tsv(df, path)
}

manifest <- readr::read_tsv(args$sample_manifest, show_col_types = FALSE, progress = FALSE)
tx2gene <- readr::read_tsv(args$tx2gene, show_col_types = FALSE, progress = FALSE)
gene_annotation <- readr::read_tsv(args$gene_annotation, show_col_types = FALSE, progress = FALSE)

required_manifest_cols <- c("sample_id", "quant_sf", "replicate", "timepoint", "timepoint_label")
missing_manifest_cols <- setdiff(required_manifest_cols, colnames(manifest))
if (length(missing_manifest_cols) > 0) {
  stop("Sample manifest is missing required columns: ", paste(missing_manifest_cols, collapse = ", "))
}

required_tx2gene_cols <- c("transcript_id", "transcript_id_stripped", "gene_id", "gene_id_stripped")
missing_tx2gene_cols <- setdiff(required_tx2gene_cols, colnames(tx2gene))
if (length(missing_tx2gene_cols) > 0) {
  stop("tx2gene table is missing required columns: ", paste(missing_tx2gene_cols, collapse = ", "))
}

required_gene_cols <- c("gene_id", "gene_id_stripped", "gene_name", "gene_type")
missing_gene_cols <- setdiff(required_gene_cols, colnames(gene_annotation))
if (length(missing_gene_cols) > 0) {
  stop("Gene annotation table is missing required columns: ", paste(missing_gene_cols, collapse = ", "))
}

manifest <- manifest %>%
  mutate(
    replicate = as.integer(replicate),
    timepoint = as.character(timepoint),
    timepoint_label = as.character(timepoint_label)
  )

time_levels <- trimws(strsplit(args$time_order, ",", fixed = TRUE)[[1]])
pairwise_defs <- combn(time_levels, 2, simplify = FALSE)
if (length(pairwise_defs) == 0) {
  stop("At least two timepoints are required for pairwise DE analysis.")
}
manifest$timepoint <- factor(manifest$timepoint, levels = time_levels, ordered = FALSE)
if (any(is.na(manifest$timepoint))) {
  stop("Sample manifest contains timepoints not present in --time-order.")
}

files <- manifest$quant_sf
names(files) <- manifest$sample_id
if (!all(file.exists(files))) {
  missing_files <- files[!file.exists(files)]
  stop("Missing quant.sf files: ", paste(missing_files, collapse = ", "))
}

manifest_df <- as.data.frame(manifest)
rownames(manifest_df) <- manifest_df$sample_id

tx2gene_import <- tx2gene %>%
  select(transcript_id_stripped, gene_id_stripped) %>%
  distinct()

tx2gene_lookup <- unique(tx2gene_import$transcript_id_stripped)
mapping_summary <- lapply(seq_along(files), function(i) {
  quant_names <- readr::read_tsv(
    files[[i]],
    show_col_types = FALSE,
    progress = FALSE,
    col_types = cols_only(Name = col_character())
  )$Name
  quant_ids <- unique(strip_version(quant_names))
  mapped_n <- sum(quant_ids %in% tx2gene_lookup)
  tibble(
    sample_id = names(files)[[i]],
    quant_transcripts = length(quant_ids),
    mapped_transcripts = mapped_n,
    unmapped_transcripts = length(quant_ids) - mapped_n,
    mapped_fraction = ifelse(length(quant_ids) > 0, mapped_n / length(quant_ids), NA_real_)
  )
}) %>% bind_rows()

write_tsv_with_dirs(mapping_summary, args$out_mapping_summary)

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = as.data.frame(tx2gene_import),
  txOut = FALSE,
  ignoreTxVersion = TRUE,
  countsFromAbundance = "no"
)

genes_before_prefilter <- nrow(txi$counts)
prefilter_keep <- rowSums(txi$counts >= args$prefilter_min_count) >= args$prefilter_min_samples
prefilter_n <- sum(prefilter_keep)
if (prefilter_n == 0) {
  stop("No genes passed the RNA-seq prefilter.")
}

txi$counts <- txi$counts[prefilter_keep, , drop = FALSE]
txi$abundance <- txi$abundance[prefilter_keep, , drop = FALSE]
txi$length <- txi$length[prefilter_keep, , drop = FALSE]

build_dds <- function() {
  DESeqDataSetFromTximport(
    txi = txi,
    colData = manifest_df,
    design = ~ timepoint
  )
}

if (args$threads > 1) {
  BiocParallel::register(BiocParallel::MulticoreParam(workers = args$threads))
}

dds_lrt <- DESeq(build_dds(), test = "LRT", reduced = ~ 1, parallel = args$threads > 1)
lrt_res <- results(dds_lrt, alpha = args$alpha)
dds_wald <- DESeq(build_dds(), test = "Wald", parallel = args$threads > 1)

gene_annotation_unique <- gene_annotation %>%
  distinct(gene_id_stripped, .keep_all = TRUE)

annotated_matrix_df <- function(mat) {
  tibble(gene_id_stripped = rownames(mat)) %>%
    left_join(gene_annotation_unique, by = "gene_id_stripped") %>%
    bind_cols(as_tibble(mat, .name_repair = "minimal"))
}

annotate_deseq_results_df <- function(res_obj) {
  as.data.frame(res_obj) %>%
    rownames_to_column("gene_id_stripped") %>%
    as_tibble() %>%
    left_join(gene_annotation_unique, by = "gene_id_stripped") %>%
    select(gene_id, gene_id_stripped, gene_name, gene_type, everything()) %>%
    arrange(padj, pvalue)
}

annotate_pairwise_results_df <- function(res_obj, contrast_id, numerator, denominator) {
  annotate_deseq_results_df(res_obj) %>%
    mutate(
      contrast_id = contrast_id,
      numerator_timepoint = numerator,
      denominator_timepoint = denominator,
      .before = gene_id
    )
}

write_tsv_with_dirs(annotated_matrix_df(txi$counts), args$out_gene_counts)
write_tsv_with_dirs(annotated_matrix_df(txi$abundance), args$out_gene_tpm)
write_tsv_with_dirs(annotated_matrix_df(txi$length), args$out_gene_lengths)
write_tsv_with_dirs(annotated_matrix_df(counts(dds_wald, normalized = TRUE)), args$out_normalized_counts)

res_tbl <- annotate_deseq_results_df(lrt_res)

sig_tbl <- res_tbl %>%
  filter(!is.na(padj) & padj <= args$alpha)

write_tsv_with_dirs(res_tbl, args$out_lrt_results)
write_tsv_with_dirs(sig_tbl, args$out_lrt_significant)

pairwise_summary_tbl <- lapply(pairwise_defs, function(pair) {
  numerator <- pair[[1]]
  denominator <- pair[[2]]
  contrast_id <- paste0(numerator, "_vs_", denominator)
  pairwise_res <- results(
    dds_wald,
    contrast = c("timepoint", numerator, denominator),
    alpha = args$alpha
  )
  pairwise_tbl <- annotate_pairwise_results_df(pairwise_res, contrast_id, numerator, denominator)
  pairwise_sig_tbl <- pairwise_tbl %>%
    filter(!is.na(padj) & padj <= args$alpha)

  write_tsv_with_dirs(pairwise_tbl, file.path(args$out_pairwise_dir, paste0(contrast_id, "_results.tsv.gz")))
  write_tsv_with_dirs(
    pairwise_sig_tbl,
    file.path(args$out_pairwise_dir, paste0(contrast_id, "_significant.tsv.gz"))
  )

  tibble(
    contrast_id = contrast_id,
    numerator_timepoint = numerator,
    denominator_timepoint = denominator,
    tested_gene_count = sum(!is.na(pairwise_tbl$padj)),
    significant_gene_count = nrow(pairwise_sig_tbl)
  )
}) %>% bind_rows()

write_tsv_with_dirs(pairwise_summary_tbl, args$out_pairwise_summary)

vst_obj <- vst(dds_wald, blind = TRUE)
pca_data <- plotPCA(vst_obj, intgroup = c("timepoint", "replicate"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
pca_data$sample_id <- rownames(pca_data)
pca_tbl <- as_tibble(pca_data) %>%
  select(sample_id, timepoint, replicate, PC1, PC2)
write_tsv_with_dirs(pca_tbl, args$out_pca_tsv)

pca_plot <- ggplot(pca_tbl, aes(x = PC1, y = PC2, color = timepoint, shape = factor(replicate))) +
  geom_point(size = 3) +
  xlab(sprintf("PC1: %s%% variance", percent_var[1])) +
  ylab(sprintf("PC2: %s%% variance", percent_var[2])) +
  labs(color = "Timepoint", shape = "Replicate", title = "RNA-seq VST PCA") +
  theme_bw()

ensure_parent(args$out_pca_png)
ggsave(args$out_pca_png, plot = pca_plot, width = 7, height = 5, dpi = 300)
ensure_parent(args$out_pca_pdf)
ggsave(args$out_pca_pdf, plot = pca_plot, width = 7, height = 5)

sample_dist <- dist(t(assay(vst_obj)))
sample_dist_mat <- as.matrix(sample_dist)
ensure_parent(args$out_sample_distance_tsv)
write.table(
  sample_dist_mat,
  file = args$out_sample_distance_tsv,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

sample_labels <- manifest %>%
  transmute(label = paste0(sample_id, " | t", timepoint_label, " | rep", replicate)) %>%
  pull(label)
names(sample_labels) <- manifest$sample_id
rownames(sample_dist_mat) <- sample_labels[rownames(sample_dist_mat)]
colnames(sample_dist_mat) <- sample_labels[colnames(sample_dist_mat)]

ensure_parent(args$out_sample_distance_png)
png(args$out_sample_distance_png, width = 1800, height = 1600, res = 200)
pheatmap(sample_dist_mat, main = "RNA-seq sample distances (VST)")
dev.off()

ensure_parent(args$out_sample_distance_pdf)
pdf(args$out_sample_distance_pdf, width = 9, height = 8)
pheatmap(sample_dist_mat, main = "RNA-seq sample distances (VST)")
dev.off()

summary_lines <- c(
  "# RNA-seq Time-Series Summary",
  "",
  sprintf("- Samples analyzed: `%d`.", nrow(manifest)),
  sprintf("- Timepoints: `%s`.", paste(time_levels, collapse = ", ")),
  sprintf("- Genes before prefilter: `%d`.", genes_before_prefilter),
  sprintf("- Genes after prefilter: `%d`.", prefilter_n),
  sprintf("- Omnibus time-course LRT significant genes at adjusted p-value <= %.2f: `%d`.", args$alpha, nrow(sig_tbl)),
  sprintf(
    "- Mean transcript mapping fraction across samples: `%.4f`.",
    mean(mapping_summary$mapped_fraction, na.rm = TRUE)
  ),
  "",
  "## Pairwise DESeq2 contrasts",
  "",
  vapply(
    seq_len(nrow(pairwise_summary_tbl)),
    function(i) {
      sprintf(
        "- `%s`: `%d` significant genes.",
        pairwise_summary_tbl$contrast_id[[i]],
        pairwise_summary_tbl$significant_gene_count[[i]]
      )
    },
    character(1)
  ),
  "",
  "Outputs include gene-level counts, TPM, effective lengths, normalized counts, DESeq2 LRT results, pairwise Wald contrasts, PCA, and sample-distance QC."
)
ensure_parent(args$out_summary_md)
writeLines(summary_lines, con = args$out_summary_md)
