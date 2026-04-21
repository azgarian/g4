# G-Quadruplex UV-Response Pipeline

## Active Datasets

The tables below list the dataset families that are actively consumed by the current workflow. They are limited to live inputs used by uncommented modules in [`workflow/Snakefile`](workflow/Snakefile).

### Reference and genome resources

| Dataset family | Assay / source | Active conditions | Role in pipeline | Where defined |
| --- | --- | --- | --- | --- |
| Reference genome | hg38 primary assembly | hg38 only | Core reference for alignment, sequence extraction, structure annotation, and interval generation | [`config/config.yaml`](config/config.yaml) |
| Gene annotation | Gencode v48 GTF | hg38 only | Promoter/TSS annotation and gene-context analyses | [`config/config.yaml`](config/config.yaml) |
| RNA-seq gene annotation | Gencode v35 GTF | hg38 only | Transcript-to-gene mapping and gene-level annotation for the active RNA-seq time-series branch | [`config/config.yaml`](config/config.yaml), [`workflow/rules/rnaseq.smk`](workflow/rules/rnaseq.smk) |
| ENCODE blacklist | ENCODE blacklist BED | hg38 only | Removes problematic genomic regions before downstream analyses | [`config/config.yaml`](config/config.yaml) |
| Repeat annotations | UCSC RepeatMasker table | hg38 only | Reference annotation used by active genome-prep resources | [`config/config.yaml`](config/config.yaml) |
| Mappability track | Umap/Bismap k50 BED | hg38 only | Reference resource prepared for active genome context support | [`config/config.yaml`](config/config.yaml) |
| LiftOver chain | hg19-to-hg38 chain | hg19 to hg38 | Required to lift observed quadruplex tracks into the active genome build | [`config/config.yaml`](config/config.yaml) |

### G4 and PQS reference datasets

| Dataset family | Assay / source | Active conditions | Role in pipeline | Where defined |
| --- | --- | --- | --- | --- |
| OQS, K-stabilized | GEO observed quadruplex BEDs | plus/minus strand, K condition | Orthogonal G4 support set used to judge whether inferred UV-responsive loci overlap experimentally observed quadruplex-rich regions | [`config/config.yaml`](config/config.yaml) |
| OQS, PDS-stabilized | GEO observed quadruplex BEDs | plus/minus strand, PDS condition | Orthogonal G4 support set used to test whether inferred UV-responsive loci remain consistent with an independent observed quadruplex condition | [`config/config.yaml`](config/config.yaml) |
| GC-rich background | Genome-derived sampled BED | structure-supported GC-rich non-OQS loci | Baseline/background comparison set used to test whether inferred mG4 behavior exceeds generic GC-rich sequence context | [`workflow/rules/gc_rich_bg.smk`](workflow/rules/gc_rich_bg.smk) |

### Experimental G4 discovery datasets

| Dataset family | Assay / source | Active conditions | Role in pipeline | Where defined |
| --- | --- | --- | --- | --- |
| G4 ChIP-seq | local paired FASTQs | HeLa, hg38 | Orthogonal experimental G4 support set used to evaluate robustness of inferred loci against assay-based G4 mapping | [`workflow/rules/g4p_chip_seq.smk`](workflow/rules/g4p_chip_seq.smk) |
| G4 CUT&Tag | local paired FASTQs | HeLa, hg38 | Orthogonal experimental G4 support set used to evaluate robustness of inferred loci against an independent assay chemistry | [`workflow/rules/g4p_cut_tag.smk`](workflow/rules/g4p_cut_tag.smk) |

### UV damage and repair datasets

| Dataset family | Assay / source | Active conditions | Role in pipeline | Where defined |
| --- | --- | --- | --- | --- |
| XR-seq style repair inputs | HeLa XR samples | `64-PP` and `CPD`, active time `12` min | Repair layer used to test how inferred UV-responsive G4 states relate to repair signal after UV | [`config/config.yaml`](config/config.yaml) |
| Damage-seq style lesion inputs | HeLa DS samples | `64-PP` and `CPD`, active times `0`, `1`, `15`, `30`, `60`, `240`, `480` min | Damage layer used to test how inferred UV-responsive G4 states relate to lesion accumulation through time | [`config/config.yaml`](config/config.yaml) |
| Simulated UV controls | simulated BED/BW tracks paired with XR/DS inputs | matched to active XR/DS samples | Null or background reference for real-vs-sim damage profile comparisons | generated in active UV damage workflows |

### Chromatin and expression context datasets

| Dataset family | Assay / source | Active conditions | Role in pipeline | Where defined |
| --- | --- | --- | --- | --- |
| ATAC-seq bigWigs | local pooled coverage bigWigs | noUV, `15m`, `30m`, `1h` | Regulatory-context layer used to test whether UV-responsive G4 groups occupy different accessibility states across the UV response |
| RNA-seq bigWig | baseline merged RNA-seq signal | `merged_t00.bw` | Regulatory-context layer used to place inferred G4 loci relative to baseline transcriptional activity |
| RNA-seq Salmon quantifications | local Salmon `quant.sf` directories | `0`, `12`, `30`, `60` minutes with discovered matching replicates | Transcriptome-wide expression time series used for tximport, DESeq2 likelihood-ratio testing, and RNA-seq QC across the UV response | [`config/config.yaml`](config/config.yaml), [`workflow/rules/rnaseq.smk`](workflow/rules/rnaseq.smk) |
| Histone ChIP-seq bigWigs | ENCODE HeLa-S3 released fold-change bigWigs | resolved active manifest tracks | Regulatory chromatin layer used to ask whether UV-responsive or damage-sensitive G4 groups occupy distinct chromatin environments | [`resources/HELA_S3/HeLa_S3_histones.tsv`](resources/HELA_S3/HeLa_S3_histones.tsv) |
| TF ChIP-seq peak BEDs | ENCODE HeLa-S3 released IDR thresholded peaks | resolved active manifest tracks | Regulatory occupancy layer used to ask whether UV-responsive or external-support G4 groups align with distinct TF-associated contexts | [`resources/HELA_S3/HeLa_S3_tfs.tsv`](resources/HELA_S3/HeLa_S3_tfs.tsv) |
| HeLa-S3 ChIP peak catalogs | ENCODE HeLa-S3 released TF and histone IDR thresholded peak BEDs | resolved active manifest tracks | [`config/config.yaml`](config/config.yaml) |

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
