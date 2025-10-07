#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G4–Transcription Analysis with Bioframe

Pipeline:
1) Load genome chromsizes (from file or UCSC via bioframe).
2) Build gene models and strand-aware promoters from GTF/BED12.
3) Read multiple G4 BEDs (e.g., g4_0.bed, g4_12.bed, g4_30.bed, g4_60.bed).
4) Count G4s in promoters using bioframe overlap ops.
5) Load RNA-seq quant (Salmon/Kallisto): TPMs and counts/NumReads.
6) Normalize expression:
   - Default: TPM
   - Optional: TMM via R edgeR using rpy2 (if --tmm and edgeR available)
7) Correlate promoter G4 burden vs expression per timepoint; also Δ vs Δ.
8) Export tidy CSVs and plots (PNG/PDF).

Examples
--------
python g4_transcription_bioframe.py \
  --gtf resources/ref_genomes/hg38/gencode.v49.annotation.gtf.gz \
  --assembly hg38 \
  --g4 results/g4_miner/oqs_0.bed results/g4_miner/oqs_12.bed results/g4_miner/oqs_30.bed results/g4_miner/oqs_60.bed \
  --rna resources/rna_seq/SU_10/quant.sf resources/rna_seq/SU_112/quant.sf resources/rna_seq/SU_130/quant.sf resources/rna_seq/SU_160/quant.sf \
  --rna-sample-names T0 T12 T30 T60 \
  --prom-up 2000 --prom-down 500 \
  --gene-id-col gene_id --gene-name-col gene_name \
  --outdir results/g4_txn --tmm

Notes
-----
- If you don’t have a chromsizes file, omit --chromsizes and set --assembly hg38 (or hg19, mm10, etc.)
- For GTFs, provide the attributes to extract gene IDs/names (defaults sensible for GENCODE).
- RNA quant assumed Salmon quant.sf (columns: Name, TPM, NumReads, etc.). Kallisto tsv (abundance.tsv) is also supported.
- TMM requires R, edgeR, and rpy2. Fallback is TPM if TMM not requested or fails.
"""

from __future__ import annotations
import argparse
import gzip
import os
import scipy
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import bioframe as bf

# ----------------------------- I/O utils -----------------------------

def parse_gtf_attributes(attribute_series: pd.Series, target_keys: List[str]) -> pd.DataFrame:
    """
    Parse the GTF attribute field into a DataFrame with requested keys.
    Example attribute string: key "value"; key2 "value2";
    """
    parsed_rows: List[Dict[str, Optional[str]]] = []
    for raw in attribute_series.fillna("").astype(str):
        row: Dict[str, Optional[str]] = {k: None for k in target_keys}
        # Split on semicolons, which delimit key-value pairs in GTF
        for field in raw.strip().split(";"):
            field = field.strip()
            if not field:
                continue
            # Split only on the first space to keep quoted value intact
            parts = field.split(" ", 1)
            if len(parts) != 2:
                continue
            key, value = parts[0], parts[1].strip()
            # Remove surrounding quotes if present
            if value.startswith('"') and value.endswith('"'):
                value = value[1:-1]
            # Some GTFs may have trailing semicolons without space
            value = value.strip().strip('"')
            if key in row:
                row[key] = value
        parsed_rows.append(row)
    return pd.DataFrame(parsed_rows)

def read_chromsizes(assembly: Optional[str], chromsizes_path: Optional[str]) -> pd.Series:
    if chromsizes_path:
        cs = pd.read_csv(chromsizes_path, sep="\t", header=None, names=["chrom", "length"])
        return pd.Series(cs.length.values, index=cs.chrom)
    if assembly:
        return bf.fetch_chromsizes(assembly, filter_chroms=True)
    raise ValueError("Provide either --chromsizes or --assembly.")

def open_any(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")

def ensure_outdir(p: str | Path) -> Path:
    p = Path(p)
    p.mkdir(parents=True, exist_ok=True)
    return p

# ------------------------ Gene & promoter models ---------------------

def load_genes_from_gtf(
    gtf_path: str,
    gene_id_attr: str = "gene_id",
    gene_name_attr: str = "gene_name",
) -> pd.DataFrame:
    """
    Returns BED-like gene table: chrom, start, end, strand, gene_id, gene_name, feature
    """
    # Standard GTF columns
    gtf_cols = [
        "seqname", "source", "feature", "start", "end",
        "score", "strand", "frame", "attribute"
    ]
    df = bf.read_table(gtf_path, names=gtf_cols, sep="\t", comment="#")
    genes = df[df["feature"] == "gene"].copy()
    # Parse attributes locally to avoid dependency on bioframe internals
    attrs = parse_gtf_attributes(genes["attribute"], [gene_id_attr, gene_name_attr])
    for col in [gene_id_attr, gene_name_attr]:
        if col not in attrs.columns:
            raise ValueError(f"Attribute '{col}' not found in GTF attributes.")
    genes = genes.rename(columns={"seqname": "chrom", "start": "start", "end": "end"})
    genes = genes[["chrom", "start", "end", "strand"]].join(
        attrs[[gene_id_attr, gene_name_attr]]
    )
    genes = genes.rename(columns={gene_id_attr: "gene_id", gene_name_attr: "gene_name"})
    genes["feature"] = "gene"
    # Enforce ints
    genes["start"] = genes["start"].astype(int)
    genes["end"] = genes["end"].astype(int)
    return bf.sort_bedframe(genes)

def load_tx2gene_from_gtf(gtf_path: str) -> pd.DataFrame:
    """
    Build transcript->gene mapping from GTF transcript records.
    Returns DataFrame with columns: transcript_id, gene_id
    """
    gtf_cols = [
        "seqname", "source", "feature", "start", "end",
        "score", "strand", "frame", "attribute"
    ]
    df = bf.read_table(gtf_path, names=gtf_cols, sep="\t", comment="#")
    tx = df[df["feature"] == "transcript"]["attribute"].copy()
    attrs = parse_gtf_attributes(tx, ["transcript_id", "gene_id"])
    attrs = attrs.dropna(subset=["transcript_id", "gene_id"]).drop_duplicates()
    return attrs.rename(columns={"transcript_id": "transcript_id", "gene_id": "gene_id"})

def load_genes_from_bed12(bed12_path: str) -> pd.DataFrame:
    """
    Minimal BED12 gene model: chrom, start, end, name, score, strand
    """
    cols = ["chrom","start","end","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"]
    df = bf.read_table(bed12_path, names=cols)
    df = df[["chrom","start","end","strand","name"]].rename(columns={"name":"gene_id"})
    df["gene_name"] = df["gene_id"]
    df["feature"] = "gene"
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return bf.sort_bedframe(df)

def promoters_from_genes(
    genes: pd.DataFrame,
    upstream: int = 1000,
    downstream: int = 100,
    chromsizes: Optional[pd.Series] = None,
) -> pd.DataFrame:
    """
    Strand-aware promoter windows around TSS:
    - '+' strand: [TSS - upstream, TSS + downstream)
    - '-' strand: [TSS - downstream, TSS + upstream)
    """
    g = genes.copy()
    tss = np.where(g["strand"] == "+", g["start"], g["end"])
    prom_start = np.where(g["strand"] == "+", tss - upstream, tss - downstream)
    prom_end   = np.where(g["strand"] == "+", tss + downstream, tss + upstream)

    promoters = pd.DataFrame({
        "chrom": g["chrom"].values,
        "start": prom_start.astype(int),
        "end":   prom_end.astype(int),
        "strand": g["strand"].values,
        "gene_id": g["gene_id"].values,
        "gene_name": g["gene_name"].values
    })

    if chromsizes is not None:
        # Manually clamp intervals to chromosome bounds using provided chromsizes
        promoters = promoters.merge(
            chromsizes.rename("chrom_length"),
            left_on="chrom",
            right_index=True,
            how="left",
        )
        promoters["start"] = promoters["start"].clip(lower=0)
        promoters["end"] = promoters[["end", "chrom_length"]].min(axis=1).astype(int)
        promoters["start"] = promoters["start"].astype(int)
        promoters = promoters[promoters["end"] > promoters["start"]]
        promoters = promoters.drop(columns=["chrom_length"])
    promoters["feature"] = "promoter"
    return bf.sort_bedframe(promoters)

# --------------------------- G4 processing ---------------------------

def read_g4_bed(path: str) -> pd.DataFrame:
    """
    Read BED3+ G4 intervals.
    """
    # Read as generic BED, enforce required column names
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype={0: str})
    # Assign up to 6 standard BED columns if present
    bed_cols = ["chrom", "start", "end", "name", "score", "strand"]
    ncols = min(len(df.columns), len(bed_cols))
    df = df.iloc[:, :ncols]
    df.columns = bed_cols[:ncols]
    # Ensure integer coordinates
    if "start" in df.columns:
        df["start"] = df["start"].astype(int)
    if "end" in df.columns:
        df["end"] = df["end"].astype(int)
    return df

def count_g4_in_promoters(
    promoters: pd.DataFrame,
    g4_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Count overlaps (interval hits) per promoter.
    """
    ovl = bf.overlap(promoters, g4_df, suffixes=("_prom", "_g4"))
    # Each row in ovl is one overlap; count per promoter (gene_id)
    gid_col = "gene_id" if "gene_id" in ovl.columns else ("gene_id_prom" if "gene_id_prom" in ovl.columns else None)
    if gid_col is None:
        raise KeyError("gene_id column not found in overlap result")
    counts = ovl.groupby(gid_col).size().rename("g4_count").to_frame()
    out = promoters[["gene_id", "gene_name"]].drop_duplicates().merge(
        counts, left_on="gene_id", right_index=True, how="left"
    )
    out["g4_count"] = out["g4_count"].fillna(0).astype(int)
    return out

# --------------------------- RNA quant I/O ---------------------------

def read_salmon_quant(path: str) -> pd.DataFrame:
    """
    Expects Salmon quant.sf with columns: Name, TPM, NumReads, Length, EffectiveLength
    """
    df = pd.read_csv(path, sep="\t")
    if {"Name","TPM"}.issubset(df.columns):
        # Keep transcript-level; aggregate later via tx2gene
        df = df.rename(columns={"Name":"transcript_id"})
        return df
    raise ValueError(f"{path} does not look like Salmon quant.sf")

def read_kallisto_abundance(path: str) -> pd.DataFrame:
    """
    Kallisto abundance.tsv: target_id length eff_length est_counts tpm
    """
    df = pd.read_csv(path, sep="\t")
    if {"target_id","tpm","est_counts"}.issubset(df.columns):
        df = df.rename(columns={"target_id":"transcript_id","tpm":"TPM","est_counts":"NumReads"})
        return df
    raise ValueError(f"{path} does not look like Kallisto abundance.tsv")

def load_rna_quant(paths: List[str]) -> List[pd.DataFrame]:
    out = []
    for p in paths:
        p = str(p)
        if p.endswith("quant.sf"):
            out.append(read_salmon_quant(p))
        elif p.endswith("abundance.tsv"):
            out.append(read_kallisto_abundance(p))
        else:
            # try generic TSV with gene_id, TPM, NumReads
            df = pd.read_csv(p, sep="\t")
            if "gene_id" in df.columns and ("TPM" in df.columns or "tpm" in df.columns):
                if "tpm" in df.columns:
                    df["TPM"] = df["tpm"]
                out.append(df)
            else:
                raise ValueError(f"Unsupported RNA quant format: {p}")
    return out

def aggregate_tx_to_gene(df: pd.DataFrame, tx2gene: pd.DataFrame) -> pd.DataFrame:
    """
    Map transcript-level quant to gene-level by summing TPM and NumReads per gene_id.
    Expects columns: transcript_id, TPM, optionally NumReads.
    Returns DataFrame with columns: gene_id, TPM, optionally NumReads.
    """
    if "transcript_id" not in df.columns:
        # Already gene-level
        return df
    merged = df.merge(tx2gene, on="transcript_id", how="left")
    # Drop transcripts without gene mapping
    merged = merged.dropna(subset=["gene_id"]) 
    agg_cols = {"TPM": "sum"}
    if "NumReads" in merged.columns:
        agg_cols["NumReads"] = "sum"
    gene_df = merged.groupby("gene_id", as_index=False).agg(agg_cols)
    return gene_df

# --------------------------- Normalization --------------------------

def tmm_normalize_with_rpy2(count_mtx: pd.DataFrame) -> Optional[pd.DataFrame]:
    """
    Returns CPM normalized via edgeR::calcNormFactors + cpm() if possible.
    Requires rpy2 + edgeR. If unavailable, returns None.
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()
        ro.r("suppressPackageStartupMessages(library(edgeR))")
        r_counts = pandas2ri.py2rpy(count_mtx)
        ro.globalenv["count_mtx"] = r_counts
        ro.r("""
        d <- DGEList(counts = count_mtx)
        d <- calcNormFactors(d, method="TMM")
        cpm_mat <- cpm(d, normalized.lib.sizes=TRUE, log=FALSE)
        """)
        cpm = ro.r("as.data.frame(cpm_mat)")
        cpm = pandas2ri.rpy2py(cpm)
        cpm.index = count_mtx.index
        cpm.columns = count_mtx.columns
        return cpm
    except Exception as e:
        print(f"[WARN] TMM normalization failed or unavailable: {e}")
        return None

def prepare_expression_matrix(
    rna_list: List[pd.DataFrame],
    sample_names: List[str],
    gene_id_col: str = "gene_id",
    prefer_counts: bool = True,
    use_tmm: bool = False
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns (TPM_matrix, COUNT_or_TMM_matrix) indexed by gene_id with samples as columns.
    """
    if len(rna_list) != len(sample_names):
        raise ValueError("Number of RNA files must match number of sample names.")

    # TPM
    tpms = []
    for df, s in zip(rna_list, sample_names):
        col = "TPM" if "TPM" in df.columns else "tpm"
        t = df[[gene_id_col, col]].rename(columns={col: s})
        tpms.append(t)
    tpm_mat = tpms[0]
    for t in tpms[1:]:
        tpm_mat = tpm_mat.merge(t, on=gene_id_col, how="outer")
    tpm_mat = tpm_mat.fillna(0.).set_index(gene_id_col)

    # Counts (if present)
    count_present = all(("NumReads" in df.columns) for df in rna_list)
    count_mat = None
    if count_present:
        counts = []
        for df, s in zip(rna_list, sample_names):
            c = df[[gene_id_col, "NumReads"]].rename(columns={"NumReads": s})
            counts.append(c)
        count_mat = counts[0]
        for c in counts[1:]:
            count_mat = count_mat.merge(c, on=gene_id_col, how="outer")
        count_mat = count_mat.fillna(0.).set_index(gene_id_col)

    # TMM via edgeR if requested and counts available
    tmm_mat = None
    if use_tmm and count_mat is not None:
        tmm_mat = tmm_normalize_with_rpy2(count_mat.astype(float))

    # Choose expression matrix for correlation: TMM > TPM
    expr_for_corr = tmm_mat if (use_tmm and tmm_mat is not None) else tpm_mat
    return tpm_mat, (expr_for_corr if prefer_counts else tpm_mat)

# ------------------------------ Analysis ----------------------------

def collect_g4_promoter_counts(
    promoters: pd.DataFrame,
    g4_paths: List[str],
    g4_names: List[str],
) -> pd.DataFrame:
    matrices = []
    for pth, name in zip(g4_paths, g4_names):
        g4 = read_g4_bed(pth)
        cnt = count_g4_in_promoters(promoters, g4)
        cnt = cnt.rename(columns={"g4_count": name})
        matrices.append(cnt[["gene_id", name]])
    out = matrices[0]
    for m in matrices[1:]:
        out = out.merge(m, on="gene_id", how="outer")
    out = out.fillna(0).merge(promoters[["gene_id","gene_name"]].drop_duplicates(), on="gene_id", how="left")
    return out.set_index("gene_id")

def correlate_g4_expression(
    g4_counts: pd.DataFrame, expr: pd.DataFrame, sample_names: List[str], gene_name_col: str = "gene_name"
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    For each sample/timepoint name present in both (g4_counts columns & expr columns),
    compute Pearson and Spearman correlations across genes.
    Returns (per_timepoint_stats_df, per_gene_tidy_df).
    """
    common = [s for s in sample_names if s in g4_counts.columns and s in expr.columns]
    stats = []
    tidy_rows = []
    for s in common:
        x = g4_counts[s].astype(float)
        y = expr[s].astype(float).reindex(x.index)
        ok = (x.notna() & y.notna())
        xv = x[ok]; yv = y[ok]
        pear = xv.corr(yv, method="pearson")
        spear = xv.corr(yv, method="spearman")
        stats.append({"sample": s, "pearson_r": pear, "spearman_rho": spear, "n_genes": ok.sum()})
        # Tidy for scatter plotting
        tmp = pd.DataFrame({"gene_id": xv.index, "sample": s, "g4_promoter": xv.values, "expr": yv.values})
        tmp[gene_name_col] = g4_counts.loc[xv.index, gene_name_col].values
        tidy_rows.append(tmp)
    stats_df = pd.DataFrame(stats)
    tidy_df = pd.concat(tidy_rows, ignore_index=True) if tidy_rows else pd.DataFrame()
    return stats_df, tidy_df

def delta_vs_delta(
    g4_counts: pd.DataFrame, expr: pd.DataFrame, baseline: str, compare_to: List[str]
) -> pd.DataFrame:
    """
    ΔG4 (count) vs ΔExpression (expr) relative to baseline for each target sample.
    Returns tidy long DF.
    """
    rows = []
    for s in compare_to:
        if baseline not in g4_counts.columns or s not in g4_counts.columns:
            continue
        if baseline not in expr.columns or s not in expr.columns:
            continue
        d_g4 = g4_counts[s] - g4_counts[baseline]
        d_expr = expr[s] - expr[baseline]
        ok = d_g4.notna() & d_expr.notna()
        rows.append(pd.DataFrame({
            "gene_id": g4_counts.index[ok],
            "sample": s,
            "delta_g4": d_g4[ok].values,
            "delta_expr": d_expr[ok].values
        }))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()

# ------------------------------- Plots -------------------------------

def scatter_plot(
    df: pd.DataFrame,
    xcol: str, ycol: str,
    title: str, outpath: Path,
    alpha: float = 0.25
):
    plt.figure(figsize=(5,5))
    plt.scatter(df[xcol], df[ycol], s=6, alpha=alpha)
    plt.xlabel(xcol)
    plt.ylabel(ycol)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath.with_suffix(".png"), dpi=300)
    plt.savefig(outpath.with_suffix(".pdf"))
    plt.close()

# ------------------------------- CLI --------------------------------

def main():
    ap = argparse.ArgumentParser(description="Correlate promoter G4 burden with transcription using bioframe.")
    # References
    ref = ap.add_argument_group("References")
    ref.add_argument("--gtf", type=str, help="GTF annotation (preferred).")
    ref.add_argument("--bed12", type=str, help="BED12 gene annotation (alternative to GTF).")
    ref.add_argument("--chromsizes", type=str, help="Tab-separated chrom sizes: chrom<TAB>length.")
    ref.add_argument("--assembly", type=str, help="Assembly name for chromsizes via bioframe (e.g., hg38).")
    ref.add_argument("--gene-id-col", type=str, default="gene_id", help="GTF attribute for gene id.")
    ref.add_argument("--gene-name-col", type=str, default="gene_name", help="GTF attribute for gene name.")

    # Promoters
    prom = ap.add_argument_group("Promoters")
    prom.add_argument("--prom-up", type=int, default=1000, help="Upstream bp for promoter.")
    prom.add_argument("--prom-down", type=int, default=100, help="Downstream bp for promoter.")

    # Inputs
    inp = ap.add_argument_group("Inputs")
    inp.add_argument("--g4", nargs="+", required=True, help="G4 BED files (order = timepoints).")
    inp.add_argument("--g4-names", nargs="+", help="Names for G4 files (e.g., T0 T12 T30 T60).")
    inp.add_argument("--rna", nargs="+", required=True, help="RNA quant files (Salmon quant.sf or Kallisto abundance.tsv).")
    inp.add_argument("--rna-sample-names", nargs="+", required=True, help="Names for RNA samples matching order of --rna.")

    # Normalization
    norm = ap.add_argument_group("Normalization")
    norm.add_argument("--tmm", action="store_true", help="Use TMM via rpy2+edgeR if available.")
    norm.add_argument("--prefer-counts", action="store_true", help="Prefer counts/TMM over TPM for correlation if available.")

    # Δ analysis
    delta = ap.add_argument_group("Delta Analysis")
    delta.add_argument("--delta-baseline", type=str, default=None, help="Baseline sample name for Δ (e.g., T0).")
    delta.add_argument("--delta-compare", nargs="+", default=None, help="Samples to compare vs baseline (e.g., T12 T30 T60).")

    # Output
    out = ap.add_argument_group("Output")
    out.add_argument("--outdir", type=str, required=True, help="Output directory.")

    args = ap.parse_args()

    outdir = ensure_outdir(args.outdir)

    # Chromsizes
    chromsizes = read_chromsizes(args.assembly, args.chromsizes)

    # Genes
    if args.gtf:
        genes = load_genes_from_gtf(args.gtf, gene_id_attr=args.gene_id_col, gene_name_attr=args.gene_name_col)
    elif args.bed12:
        genes = load_genes_from_bed12(args.bed12)
    else:
        raise ValueError("Provide --gtf or --bed12.")

    # Promoters
    promoters = promoters_from_genes(
        genes, upstream=args.prom_up, downstream=args.prom_down, chromsizes=chromsizes
    )
    promoters.to_csv(outdir / "promoters.bed.gz", sep="\t", index=False, compression="gzip")

    # G4 names
    if args.g4_names:
        if len(args.g4_names) != len(args.g4):
            raise ValueError("--g4-names length must match --g4")
        g4_names = args.g4_names
    else:
        # Derive names from filenames
        g4_names = [Path(p).stem for p in args.g4]

    # G4 promoter counts
    g4_counts = collect_g4_promoter_counts(promoters, args.g4, g4_names)
    g4_counts.to_csv(outdir / "g4_promoter_counts.tsv", sep="\t")

    # RNA quant
    rna_list = load_rna_quant(args.rna)
    # Aggregate transcript-level to gene-level using GTF mapping
    tx2gene = load_tx2gene_from_gtf(args.gtf)
    rna_list_gene = [aggregate_tx_to_gene(df, tx2gene) for df in rna_list]
    tpm_mat, expr_mat = prepare_expression_matrix(
        rna_list_gene, args.rna_sample_names,
        gene_id_col="gene_id",
        prefer_counts=args.prefer_counts,
        use_tmm=args.tmm
    )
    # Log2-transform TMM (or counts-for-corr) if TMM requested
    if args.tmm:
        expr_mat = np.log2(expr_mat.astype(float) + 1.0)
    tpm_mat.to_csv(outdir / "rna_TPM_matrix.tsv", sep="\t")
    expr_mat.to_csv(outdir / ("rna_TMM_or_TPM_matrix.tsv" if args.tmm else "rna_expr_for_corr.tsv"), sep="\t")

    # Align gene universe
    # Merge via gene_id; use gene_name from g4_counts (if present)
    all_genes = g4_counts.index.union(expr_mat.index)
    g4_aln = g4_counts.reindex(all_genes).fillna(0)
    expr_aln = expr_mat.reindex(all_genes).fillna(0)

    # Correlations per timepoint
    stats_df, tidy_df = correlate_g4_expression(
        g4_aln, expr_aln, sample_names=args.rna_sample_names, gene_name_col="gene_name"
    )
    stats_df.to_csv(outdir / "correlations_per_timepoint.tsv", sep="\t", index=False)
    tidy_df.to_csv(outdir / "scatter_tidy.tsv", sep="\t", index=False)

    # Plots per timepoint (only if correlations were computed)
    if not stats_df.empty and "sample" in stats_df.columns:
        for s in stats_df["sample"]:
            sdf = tidy_df[tidy_df["sample"] == s]
            if not sdf.empty:
                scatter_plot(
                    sdf, "g4_promoter", "expr",
                    f"G4 promoter burden vs expression ({s})",
                    outdir / f"scatter_g4_vs_expr_{s}"
                )

    # Δ analysis (optional)
    if args.delta_baseline and args.delta_compare:
        ddf = delta_vs_delta(g4_aln, expr_aln, baseline=args.delta_baseline, compare_to=args.delta_compare)
        ddf.to_csv(outdir / "delta_vs_delta.tsv", sep="\t", index=False)
        for s in ddf["sample"].unique():
            sdf = ddf[ddf["sample"] == s]
            scatter_plot(
                sdf, "delta_g4", "delta_expr",
                f"ΔG4 vs ΔExpression ({s} - {args.delta_baseline})",
                outdir / f"scatter_delta_{s}_vs_{args.delta_baseline}"
            )

    # Small QC summary
    with open(outdir / "SUMMARY.txt", "w") as fh:
        fh.write("G4–Transcriptome Bioframe Analysis\n")
        fh.write(f"Genes: {genes.shape[0]}\n")
        fh.write(f"Promoters: {promoters.shape[0]} (up={args.prom_up}, down={args.prom_down})\n")
        fh.write(f"G4 files: {len(args.g4)} -> columns: {', '.join(g4_names)}\n")
        fh.write(f"RNA samples: {len(args.rna)} -> {', '.join(args.rna_sample_names)}\n")
        if args.tmm:
            fh.write("Normalization: TMM requested (fallback to TPM if TMM failed)\n")
        else:
            fh.write("Normalization: TPM (default)\n")
        fh.write(f"Correlation rows: {len(tidy_df)}\n")
        if args.delta_baseline and args.delta_compare:
            fh.write(f"Delta baseline: {args.delta_baseline}; compare: {', '.join(args.delta_compare)}\n")

if __name__ == "__main__":
    main()
