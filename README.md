# G4 Pipeline (Reduced Spin-Off)

This repository contains a smaller, restart-from-scratch version of the original g-quadruplex pipeline. The overall goal is the same: collect G4-related loci from several sources, place everything on a common `hg38` coordinate system, and generate structure-aware outputs that can be compared downstream. The difference is that this branch currently works with a much smaller set of inputs and a narrower active workflow.

The resource layout is intentionally kept close to the original project. Reference files live under `resources/ref_genomes/`, sample-level inputs live under `resources/samples/`, workflow rules live under `workflow/rules/`, and final outputs are written to `results/`.

## What is active in this branch

The current `workflow/Snakefile` builds five data products plus a GC-rich background set:

1. `g4chip`: process one paired-end G4 ChIP-seq FASTQ pair into MACS2 peaks and structure-aware loci
2. `g4cuttag`: process one paired-end G4 CUT&Tag FASTQ pair into MACS2 peaks and structure-aware loci
3. `g4chip_g4cuttag`: merge the prepared `g4chip` and `g4cuttag` BED intervals into one union BED with `bedtools merge`
4. `oqs_K` and `oqs_PDS`: download published OQS BED tracks, lift them from `hg19` to `hg38`, merge strands, and annotate structure
5. `gc_rich_bg`: sample GC-rich genomic windows that avoid blacklist and OQS regions, then annotate structure

Compared with the original pipeline, this spin-off is deliberately limited. It is not yet driving a large sample sheet or the broader multi-omics analyses that still appear elsewhere in the repository.

## Entry point

The active Snakemake entry point is:

```bash
workflow/Snakefile
```

There is currently no top-level `Snakefile`, so run Snakemake with `-s workflow/Snakefile`.

## Workflow summary

### 1. Reference preparation

Rules in [`workflow/rules/prepare_genome.smk`](workflow/rules/prepare_genome.smk) download and prepare shared `hg38` resources:

- genome FASTA
- Bowtie2 index
- FASTA index (`.fai`)
- genome BED
- ENCODE blacklist BED
- RepeatMasker BED
- Umap mappability BED
- GENCODE GTF
- `hg19` to `hg38` liftOver chain

These are pulled from URLs defined in `config/config.yaml`.

### 2. G4 ChIP-seq branch

Rules in [`workflow/rules/g4p_chip_seq.smk`](workflow/rules/g4p_chip_seq.smk) expect the following inputs:

- `resources/samples/g4chip/g4_hela_1.fastq`
- `resources/samples/g4chip/g4_hela_2.fastq`

Processing steps:

1. align reads to `hg38` with Bowtie2
2. filter low-quality alignments and remove duplicates
3. remove ENCODE blacklist overlaps
4. call peaks with MACS2
5. convert peaks into structure-aware loci with `workflow/scripts/prepare_external_g4_dataset.py`

### 3. G4 CUT&Tag branch

Rules in [`workflow/rules/g4p_cut_tag.smk`](workflow/rules/g4p_cut_tag.smk) expect:

- `resources/samples/g4cuttag/g4_hela_cut_tag_1.fastq`
- `resources/samples/g4cuttag/g4_hela_cut_tag_2.fastq`

Processing steps are similar to the ChIP-seq branch:

1. Bowtie2 alignment
2. BAM filtering and sorting
3. ENCODE blacklist removal
4. MACS2 peak calling
5. structure-aware peak preparation

### 4. Combined G4 peak set

Rules in [`workflow/rules/g4chip_g4cuttag_merge.smk`](workflow/rules/g4chip_g4cuttag_merge.smk):

1. take the prepared BED outputs from the `g4chip` and `g4cuttag` branches
2. concatenate and coordinate-sort the intervals
3. run `bedtools merge` to create a non-redundant union BED

This merged file is a simple interval-level union. It is not currently passed back through the structure-annotation script, so it contains merged genomic spans only.

### 5. OQS import and harmonization

Rules in [`workflow/rules/prepare_oqs.smk`](workflow/rules/prepare_oqs.smk):

1. download four strand-specific OQS BED files from GEO
2. standardize the BED fields
3. lift coordinates from `hg19` to `hg38`
4. remove ENCODE blacklist overlaps
5. combine strands into `oqs_K` and `oqs_PDS`
6. annotate each locus with a regex-based G4 structure class

The active OQS inputs are controlled by the `oqs` section of `config/config.yaml`.

### 6. GC-rich background generation

Rules in [`workflow/rules/gc_rich_bg.smk`](workflow/rules/gc_rich_bg.smk):

1. scan `chr1`-`chr22` and `chrX` in the `hg38` FASTA
2. sample fixed-width GC-rich windows
3. drop windows overlapping the ENCODE blacklist
4. drop windows overlapping prepared OQS loci
5. annotate the remaining windows with the same structure-aware procedure

Current rule parameters:

- tile size: `150 bp`
- GC threshold: `0.28`
- requested sample size: `20000`
- random seed: `42`

Note: the output filenames still contain the older `300k` stem even though the current sample size is `20000`.

## Final targets

`rule all` currently requests these files:

```text
results/g4cuttag/g4_hela_peaks_prepared.bed
results/g4chip/g4_hela_peaks_prepared.bed
results/g4chip_g4cuttag/g4_hela_chip_cuttag_merged.bed
results/oqs/oqs_K_prepared.bed
results/oqs/oqs_PDS_prepared.bed
results/gc_rich_bg/gc_rich_bg_prepared.bed
results/gc_rich_bg/gc_rich_bg_selection_summary.tsv
```

Each prepared dataset also produces:

- a tabular file with structure annotations: `*.tsv`
- a source-locus BED: `*.bed`
- a 1 bp center BED: `*_center.bed`

The prepared TSV files have the schema:

```text
chrom  start  end  name  source_signal  strand  structure_start  structure_end  structure_center  structure  source_dataset
```

The structure labels come from [`workflow/scripts/regex_structures_common.py`](workflow/scripts/regex_structures_common.py) and currently include:

- `loop1_3`
- `loop4_5`
- `loop6_7`
- `longLoop`
- `simpleBulge`
- `complexBulge`
- `twoTetrads`

Rows without a supported structure match are dropped during preparation.

## Repository layout

```text
config/
  config.yaml              Project configuration
  slurm/config.yaml        SLURM profile

resources/
  ref_genomes/             Shared genome resources
  samples/                 Input FASTQs and other source datasets

results/
  g4chip/
  g4cuttag/
  g4chip_g4cuttag/
  oqs/
  gc_rich_bg/

workflow/
  Snakefile                Active workflow entry point
  rules/                   Modular Snakemake rules
  envs/                    Per-rule conda environments
  scripts/                 Helper scripts for structure annotation and sampling
```

## Requirements

### Software

- Snakemake
- Conda or Mamba
- standard command-line tools such as `wget`, `curl`, `gunzip`, `tar`, `awk`, and `csplit`

If you want to use the bundled SLURM profile, your Snakemake installation also needs the `cluster-generic` executor plugin.

### Conda environments

Per-rule environments are defined in `workflow/envs/`. The active workflow uses environments for:

- Bowtie2
- bedtools and samtools
- liftOver
- MACS2
- Python-based helper scripts

## Running the pipeline

### Dry-run

From the repository root:

```bash
snakemake -s workflow/Snakefile -n -p
```

### Local execution

```bash
snakemake -s workflow/Snakefile --cores 8 --use-conda -p
```

Adjust `--cores` to match your machine.

### SLURM execution

The bundled profile lives in `config/slurm/config.yaml`.

```bash
snakemake -s workflow/Snakefile --profile config/slurm
```

Before using it, update the cluster-specific defaults in [`config/slurm/config.yaml`](config/slurm/config.yaml), especially:

- `account`
- `partition`
- `qos`
- memory defaults

## Configuration notes

The active workflow currently reads only a small part of `config/config.yaml`:

- `genome`
- `oqs`

Other sections are retained from the original project layout for future expansion, but they are not used by the current `workflow/Snakefile`.

## Common outputs and logs

- rule logs: `logs/<module>/...`
- cluster stderr logs: `logs/cluster/{rule}_{jobid}.err`
- prepared BED and TSV outputs: `results/`

## Practical notes for this branch

- This spin-off assumes the reduced input set hard-coded in the current rule files.
- The G4 ChIP-seq and G4 CUT&Tag branches are not yet parameterized over multiple samples.
- The repository still contains resources and configuration carried over from the original pipeline, but the active workflow is intentionally much smaller.
- If you already have intermediate files in `results/`, Snakemake will only rebuild the outdated portions of the DAG.

## Useful commands

Run one branch only:

```bash
snakemake -s workflow/Snakefile --use-conda results/oqs/oqs_K_prepared.bed
```

Unlock a working directory after an interrupted run:

```bash
snakemake -s workflow/Snakefile --unlock
```

Show a summary of tracked outputs:

```bash
snakemake -s workflow/Snakefile --summary
```

This README describes the workflow as it is currently wired in `workflow/Snakefile`. If more of the original g-quadruplex project is brought back into scope later, this document should be extended module by module rather than reverting to the old generic template.
