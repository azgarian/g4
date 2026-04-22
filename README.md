# G-Quadruplex UV-Response Pipeline

This repository contains the active Snakemake workflow for preparing HeLa G4 reference sets, quantifying RNA-seq time-series data, and testing whether promoter-proximal G4 loci are associated with transcriptional state and UV-response behavior.

The current workflow entry point is [`workflow/Snakefile`](/cta/users/cazgari/pipelines/g4/workflow/Snakefile:1). Run Snakemake with `-s workflow/Snakefile` from the repository root.

## Active workflow modules

The active DAG includes all of the modules below:

1. Reference-genome preparation for `hg38`
2. OQS import, liftOver, filtering, and structure annotation
3. GC-rich genomic background sampling and structure annotation
4. G4 ChIP-seq peak calling and structure-aware peak preparation
5. G4 CUT&Tag peak calling and structure-aware peak preparation
6. Merging G4 ChIP-seq and CUT&Tag loci into a union G4 set
7. RNA-seq time-series import from Salmon, gene-level matrices, DESeq2, and QC
8. `g4_tss`: promoter/TSS transcription analysis (Tasks 1-8)
9. `g4_tss_atac_split`: ATAC-filtered versions of the TSS analysis for `noUV`, `15m`, `30m`, and `60m`

## Inputs and configuration used by the active DAG

The active workflow consumes the following configuration blocks in [`config/config.yaml`](/cta/users/cazgari/pipelines/g4/config/config.yaml:1):

- `genome`
- `oqs`
- `rna_seq`
- `g4_tss_atac_split`

Key configured inputs currently used by the workflow:

- `hg38` primary assembly, GENCODE v48 annotation, ENCODE blacklist, RepeatMasker, Umap/Bismap, and `hg19` to `hg38` liftOver chain
- strand-specific OQS BED files for `K` and `PDS` conditions
- local paired-end FASTQs for G4 ChIP-seq and G4 CUT&Tag
- Salmon quantifications in `resources/rna-seq/SU_*`
- RNA-seq bigWigs `merged_t00.bw`, `merged_t12.bw`, `merged_t30.bw`, `merged_t60.bw`
- ATAC peak BED files for `noUV`, `15m`, `30m`, and `60m`

## Workflow summary

### 1. Reference preparation

[`workflow/rules/prepare_genome.smk`](/cta/users/cazgari/pipelines/g4/workflow/rules/prepare_genome.smk:1) downloads and prepares shared `hg38` resources:

- genome FASTA
- Bowtie2 index
- FASTA index
- genome BED
- ENCODE blacklist BED
- RepeatMasker BED
- mappability BED
- GENCODE v48 GTF
- `hg19` to `hg38` liftOver chain

### 2. Experimental and reference G4 sets

The workflow builds four complementary G4-related datasets:

- [`workflow/rules/g4p_chip_seq.smk`](/cta/users/cazgari/pipelines/g4/workflow/rules/g4p_chip_seq.smk:1): aligns G4 ChIP-seq FASTQs, filters BAMs, calls MACS2 peaks, and annotates structure-supported loci
- [`workflow/rules/g4p_cut_tag.smk`](/cta/users/cazgari/pipelines/g4/workflow/rules/g4p_cut_tag.smk:1): does the same for G4 CUT&Tag
- [`workflow/rules/g4chip_g4cuttag_merge.smk`](/cta/users/cazgari/pipelines/g4/workflow/rules/g4chip_g4cuttag_merge.smk:1): merges prepared ChIP-seq and CUT&Tag intervals into a non-redundant union BED
- [`workflow/rules/prepare_oqs.smk`](/cta/users/cazgari/pipelines/g4/workflow/rules/prepare_oqs.smk:1): downloads OQS BEDs, standardizes fields, lifts to `hg38`, filters blacklist overlaps, combines strands, and annotates structure class

### 3. GC-rich background set

[`workflow/rules/gc_rich_bg.smk`](/cta/users/cazgari/pipelines/g4/workflow/rules/gc_rich_bg.smk:1) samples GC-rich windows, removes blacklist and OQS overlaps, and applies the same structure-aware annotation used for G4 datasets.

Current parameters in the rule file:

- tile size: `150 bp`
- GC threshold: `0.28`
- requested sample size: `20000`
- random seed: `42`

### 4. RNA-seq time-series branch

[`workflow/rules/rnaseq.smk`](/cta/users/cazgari/pipelines/g4/workflow/rules/rnaseq.smk:1) discovers Salmon quantifications matching `^SU_(?P<rep>[123])(?P<time>00|12|30|60)$`, then:

1. builds a sample manifest
2. downloads GENCODE v35 for RNA-seq gene models
3. derives `tx2gene` and gene-annotation tables
4. imports transcript quantification to gene level
5. writes counts, TPM, lengths, and normalized-count matrices
6. runs omnibus DESeq2 likelihood-ratio testing across time
7. runs all pairwise DESeq2 contrasts across `0`, `12`, `30`, `60`
8. generates PCA and sample-distance QC outputs
9. writes a markdown summary of the RNA-seq branch

The current RNA-seq summary in [`results/rnaseq/summary.md`](/cta/users/cazgari/pipelines/g4/results/rnaseq/summary.md:1) reports:

- `12` samples across `0`, `12`, `30`, `60`
- `17064` genes retained after prefiltering
- `5741` significant genes in the omnibus time-course LRT at adjusted `p <= 0.05`

### 5. `g4_tss`: promoter/TSS transcription analysis

[`workflow/rules/g4_tss.smk`](/cta/users/cazgari/pipelines/g4/workflow/rules/g4_tss.smk:1) answers the question: are promoter-proximal G4s linked to active transcription and UV response?

It performs eight linked analyses:

1. extracts canonical TSSs from the RNA-seq GTF and expands them to `±500 bp` windows
2. classifies promoter windows into `G4_TSS`, `GC_bg_TSS`, or `No_overlap`
3. summarizes baseline expression by group from TPM and normalized counts
4. tests expression differences across groups and generates violin plots plus group BED files
5. generates RNA-seq TSS metaprofiles with deepTools
6. tests G4 enrichment across expression deciles
7. tests UV-response behavior by promoter group using DESeq2 LRT and pairwise contrasts
8. tests G4 structure-class enrichment across expression classes and assembles an HTML report

### 6. `g4_tss_atac_split`: ATAC-filtered promoter analyses

[`workflow/rules/g4_tss_atac_split.smk`](/cta/users/cazgari/pipelines/g4/workflow/rules/g4_tss_atac_split.smk:1) reuses the shared TSS annotation and reruns the promoter analysis after filtering merged G4 loci by ATAC peak overlap.

Configured analyses:

- `noUV` uses `merged_t00.bw` and baseline expression only
- `15m` uses ATAC peaks at `15m`, RNA-seq expression at `12 min`, and the `0_vs_12` DESeq2 contrast
- `30m` uses ATAC peaks at `30m`, RNA-seq expression at `30 min`, and the `0_vs_30` contrast
- `60m` uses ATAC peaks at `60m`, RNA-seq expression at `60 min`, and the `0_vs_60` contrast

Each timepoint produces a matched set of annotation tables, expression summaries, decile analysis outputs, structure-enrichment outputs, metaprofiles, and an HTML report. UV-response plots are produced for `15m`, `30m`, and `60m`; the `noUV` report omits UV-response sections.

## Current analysis highlights

The current results support a strong link between promoter-proximal G4 overlap and active transcription in this HeLa UV-response dataset.

- The RNA-seq branch analyzes `12` samples across `0`, `12`, `30`, and `60` minutes, retains `17064` genes after prefiltering, and detects `5741` significant time-responsive genes in the omnibus DESeq2 LRT.
- In the main `g4_tss` analysis, `9637` genes are classified as `G4_TSS`, compared with `3530` `No_overlap` genes and only `15` `GC_bg_TSS` genes. Baseline expression is highest in the `G4_TSS` group (mean `10.04`, median `10.22`) versus `No_overlap` (mean `7.52`, median `7.37`) and `GC_bg_TSS` (mean `8.22`, median `8.58`).
- Baseline expression differences are statistically strong: `G4_TSS` exceeds `No_overlap` in the one-sided Mann-Whitney test (`p_adj = 0`, rank-biserial `r = -0.49`) and also exceeds `GC_bg_TSS` despite the very small GC-rich control set (`p_adj = 0.00446`).
- G4 promoter overlap rises sharply with expression rank. The fraction of genes with promoter G4 overlap increases from `24.3%` in decile `1` to `86.0%` in decile `10`, while GC-rich control overlap stays near background (`~0.2-0.7%`). The decile trend is strong for merged G4s (`Spearman rho = 0.903`, `p = 3.44e-4`) and absent for GC-rich controls.
- UV responsiveness is also enriched among promoter-G4 genes. At `lrt_padj <= 0.05`, `3887/9637` (`40.3%`) `G4_TSS` genes are significant versus `919/3543` (`25.9%`) in all non-`G4_TSS` genes (`odds ratio = 1.93`, Fisher `p = 2.94e-54`).
- Directionality is mixed rather than uniformly induced: among `G4_TSS` genes in `pG4_DEG_intersect.tsv`, the largest response classes are `repressed_sustained` (`1035` genes), `induced_sustained` (`461`), `repressed_transient` (`416`), and `induced_transient` (`214`), with an additional `2677` genes significant only in the omnibus LRT.
- Promoter G4 strength is associated with UV-response significance. The `max_signal` stratification shows a significant association between strength bin and LRT significance (`chi-square p = 1.50e-5`), and continuous G4 signal is associated with stronger LRT significance (`Spearman p = 1.34e-6`).
- Structure-class composition is broadly stable across expression classes. No structure reaches `FDR < 0.05` in the `Q4_vs_Q1` enrichment tests, although `twoTetrads`, `loop1_3`, and `loop4_5` show nominal upward trends in the highest-expression class.
- The ATAC-restricted follow-up preserves the same baseline pattern. Filtering by G4-center overlap with ATAC peaks retains `54.96%` of merged G4 intervals in `noUV`, `36.92%` at `15m`, `57.61%` at `30m`, and `68.66%` at `60m`, and the `G4_TSS` group still has the highest baseline expression in every split (mean `10.03-10.10` vs `8.51-9.02` for `No_overlap`).
- Expression-decile enrichment becomes clearest in the later accessible-G4 analyses: the merged-G4 decile correlation is weak in `noUV` and `15m`, but significant at `30m` (`rho = 0.733`, `p = 0.0158`) and `60m` (`rho = 0.818`, `p = 0.00381`).
- UV-response enrichment also persists after ATAC filtering. The fraction of significant LRT genes in accessible `G4_TSS` promoters is `40.4%` at `15m`, `40.2%` at `30m`, and `40.7%` at `60m`, compared with `34.7%`, `33.3%`, and `31.5%` in the corresponding `No_overlap` sets.

## Output layout

### G4 preparation outputs

`results/g4chip/`

- `g4_hela_peaks.narrowPeak`, `g4_hela_peaks.xls`, `g4_hela_summits.bed`
- `g4_hela_peaks.bed`
- `g4_hela_peaks_prepared.bed`
- `g4_hela_peaks_prepared.tsv`
- `g4_hela_peaks_prepared_center.bed`

`results/g4cuttag/`

- `g4_hela_peaks.narrowPeak`, `g4_hela_peaks.xls`, `g4_hela_summits.bed`
- `g4_hela_peaks_cuttag.bed`
- `g4_hela_peaks_prepared.bed`
- `g4_hela_peaks_prepared.tsv`
- `g4_hela_peaks_prepared_center.bed`

`results/g4chip_g4cuttag/`

- `g4_hela_chip_cuttag_merged.bed`

### OQS and GC-rich background outputs

`results/oqs/`

- strand-specific downloaded and lifted intermediates for `plus/minus` and `K/PDS`
- final merged sets:
  - `oqs_K_prepared.bed`
  - `oqs_K_prepared.tsv`
  - `oqs_K_prepared_center.bed`
  - `oqs_PDS_prepared.bed`
  - `oqs_PDS_prepared.tsv`
  - `oqs_PDS_prepared_center.bed`

`results/gc_rich_bg/`

- `gc_rich_bg_sampled.bed`
- `gc_rich_bg_sampled_no_oqs.bed`
- `gc_rich_bg_prepared.bed`
- `gc_rich_bg_prepared.tsv`
- `gc_rich_bg_prepared_center.bed`
- `gc_rich_bg_selection_summary.tsv`

### RNA-seq outputs

`results/rnaseq/metadata/`

- `sample_manifest.tsv`
- `tx2gene.tsv`
- `gene_annotation.tsv`
- `tximport_mapping_summary.tsv`

`results/rnaseq/matrices/`

- `gene_counts.tsv.gz`
- `gene_tpm.tsv.gz`
- `gene_lengths.tsv.gz`
- `normalized_counts.tsv.gz`

`results/rnaseq/deseq2/`

- `time_lrt_results.tsv.gz`
- `time_lrt_significant.tsv.gz`
- `pairwise/0_vs_12_results.tsv.gz`
- `pairwise/0_vs_12_significant.tsv.gz`
- `pairwise/0_vs_30_results.tsv.gz`
- `pairwise/0_vs_30_significant.tsv.gz`
- `pairwise/0_vs_60_results.tsv.gz`
- `pairwise/0_vs_60_significant.tsv.gz`
- `pairwise/12_vs_30_results.tsv.gz`
- `pairwise/12_vs_30_significant.tsv.gz`
- `pairwise/12_vs_60_results.tsv.gz`
- `pairwise/12_vs_60_significant.tsv.gz`
- `pairwise/30_vs_60_results.tsv.gz`
- `pairwise/30_vs_60_significant.tsv.gz`
- `pairwise/pairwise_significant_gene_counts.tsv`

`results/rnaseq/qc/`

- `vst_pca.tsv`, `vst_pca.png`, `vst_pca.pdf`
- `sample_distance.tsv`, `sample_distance.png`, `sample_distance.pdf`

Top-level RNA-seq outputs:

- `results/rnaseq/summary.md`

### `g4_tss` outputs

`results/g4_tss/`

- annotation and promoter windows:
  - `canonical_tss.bed`
  - `canonical_tss_windows_1kb.bed`
  - `gene_name_table.tsv`
  - `windows_g4_intersect.bed`
  - `windows_gc_bg_intersect.bed`
  - `tss_group_annotation.tsv`
- baseline expression outputs:
  - `baseline_tpm.tsv`
  - `baseline_tpm_summary.tsv`
  - `baseline_normalized_counts.tsv`
  - `gene_expression_by_group.tsv`
- statistics and visualizations:
  - `expression_group_statistics.tsv`
  - `expression_group_fisher.tsv`
  - `expression_group_summary.tsv`
  - `expression_violin_by_group.pdf`
  - `tss_G4_TSS.bed`
  - `tss_GC_bg_TSS.bed`
  - `tss_No_overlap.bed`
  - `rnaseq_tss_matrix.gz`
  - `rnaseq_tss_metaprofile.pdf`
- decile analysis:
  - `expression_deciles.tsv`
  - `decile_overlap_fractions.tsv`
  - `decile_correlation_stats.tsv`
  - `decile_enrichment_plot.pdf`
- UV-response analysis:
  - `uv_group_lrt_summary.tsv`
  - `uv_group_fold_change_stats.tsv`
  - `uv_fold_change_by_group.pdf`
  - `uv_volcano_by_group.pdf`
- structure enrichment:
  - `g4_structure_by_expression_class.tsv`
  - `structure_class_enrichment_stats.tsv`
  - `structure_class_by_expression_class.pdf`
- final report:
  - `g4_tss_transcription_report.html`
  - `software_versions.txt`

### `g4_tss_atac_split` outputs

`results/g4_tss_atac_split/{timepoint}/` for `noUV`, `15m`, `30m`, `60m`

Shared outputs at every timepoint:

- `g4_atac_center_filtered.bed`
- `g4_atac_center_filter_summary.tsv`
- `windows_g4_intersect.bed`
- `tss_group_annotation.tsv`
- `baseline_tpm.tsv`
- `baseline_tpm_summary.tsv`
- `baseline_normalized_counts.tsv`
- `gene_expression_by_group.tsv`
- `expression_group_statistics.tsv`
- `expression_group_fisher.tsv`
- `expression_group_summary.tsv`
- `expression_violin_by_group.pdf`
- `tss_G4_TSS.bed`
- `tss_GC_bg_TSS.bed`
- `tss_No_overlap.bed`
- `rnaseq_tss_matrix.gz`
- `rnaseq_tss_metaprofile.pdf`
- `expression_deciles.tsv`
- `decile_overlap_fractions.tsv`
- `decile_correlation_stats.tsv`
- `decile_enrichment_plot.pdf`
- `g4_structure_by_expression_class.tsv`
- `structure_class_enrichment_stats.tsv`
- `structure_class_by_expression_class.pdf`
- `g4_tss_atac_split_report.html`
- `software_versions.txt`

Additional UV-response outputs for `15m`, `30m`, and `60m`:

- `uv_group_lrt_summary.tsv`
- `uv_group_fold_change_stats.tsv`
- `uv_fold_change_by_group.pdf`
- `uv_volcano_by_group.pdf`

## Structure annotation schema

Prepared G4-like datasets such as `*_prepared.tsv` use the structure annotation generated by [`workflow/scripts/prepare_external_g4_dataset.py`](/cta/users/cazgari/pipelines/g4/workflow/scripts/prepare_external_g4_dataset.py:1) and [`workflow/scripts/regex_structures_common.py`](/cta/users/cazgari/pipelines/g4/workflow/scripts/regex_structures_common.py:1).

The TSV schema is:

```text
chrom  start  end  name  source_signal  strand  structure_start  structure_end  structure_center  structure  source_dataset
```

Current structure labels include:

- `loop1_3`
- `loop4_5`
- `loop6_7`
- `longLoop`
- `simpleBulge`
- `complexBulge`
- `twoTetrads`

Rows without a supported structure match are dropped during preparation.

## Final targets requested by `rule all`

`rule all` is defined through [`workflow/rules/common.smk`](/cta/users/cazgari/pipelines/g4/workflow/rules/common.smk:1). It currently requests:

- core prepared G4 and background outputs
- the main `g4_tss` report and representative intermediate deliverables
- one `g4_tss_atac_split_report.html` per configured ATAC timepoint

This means a default run builds the full RNA-seq branch plus both promoter-analysis branches, not just the original G4 preparation steps.

## Repository layout

```text
config/
  config.yaml
  slurm/config.yaml

logs/
  cluster/
  g4_tss/
  g4_tss_atac_split/
  rnaseq/

resources/
  ref_genomes/
  rna-seq/
  samples/

results/
  g4chip/
  g4cuttag/
  g4chip_g4cuttag/
  gc_rich_bg/
  g4_tss/
  g4_tss_atac_split/
  oqs/
  rnaseq/

workflow/
  Snakefile
  envs/
  rules/
  scripts/
```

## Requirements

Software used by the active workflow:

- Snakemake
- Conda or Mamba
- Bowtie2
- samtools
- bedtools
- MACS2
- UCSC `liftOver`
- deepTools
- `wget`, `curl`, `gunzip`, `tar`, `awk`, `csplit`
- R with the packages required by [`workflow/envs/rnaseq_deseq2.yaml`](/cta/users/cazgari/pipelines/g4/workflow/envs/rnaseq_deseq2.yaml:1)

If you use the bundled SLURM profile, Snakemake also needs the `cluster-generic` executor plugin.

## Running the pipeline

Dry-run:

```bash
snakemake -s workflow/Snakefile -n -p
```

Local execution:

```bash
snakemake -s workflow/Snakefile --cores 8 --use-conda -p
```

SLURM execution with the bundled profile:

```bash
snakemake -s workflow/Snakefile --profile config/slurm
```

Before cluster execution, review [`config/slurm/config.yaml`](/cta/users/cazgari/pipelines/g4/config/slurm/config.yaml:1) and update site-specific defaults such as `account`, `partition`, `qos`, and memory settings.

## Useful commands

Run a single branch target:

```bash
snakemake -s workflow/Snakefile --use-conda results/g4_tss/g4_tss_transcription_report.html
```

Run one ATAC-split report:

```bash
snakemake -s workflow/Snakefile --use-conda results/g4_tss_atac_split/30m/g4_tss_atac_split_report.html
```

Unlock after an interrupted run:

```bash
snakemake -s workflow/Snakefile --unlock
```

Show a summary of tracked outputs:

```bash
snakemake -s workflow/Snakefile --summary
```
