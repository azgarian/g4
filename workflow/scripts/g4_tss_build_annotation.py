#!/usr/bin/env python3
"""Task 1: Build TSS group annotation table after bedtools intersect.

Reads the gene_name side table and two bedtools-intersect -c outputs
(G4 overlap counts and GC-background overlap counts), then assigns
each gene to G4_TSS / GC_bg_TSS / No_overlap using the precedence rule.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gene-table", required=True, help="TSV from g4_tss_canonical_tss.py")
    p.add_argument("--g4-intersect", required=True,
                   help="Output of bedtools intersect -c with G4 peaks (BED6+count)")
    p.add_argument("--gc-bg-intersect", required=True,
                   help="Output of bedtools intersect -c with GC-rich background (BED6+count)")
    p.add_argument("--out-tsv", required=True, help="Output annotation TSV")
    p.add_argument("--log", default=None)
    return p.parse_args()


def read_intersect_counts(path: str) -> pd.Series:
    """Return a Series: gene_id -> overlap_count (bedtools -c output col 4 and col 7)."""
    df = pd.read_csv(path, sep="\t", header=None)
    # BED6+count: cols 0-5 are BED fields, col 6 is count
    # col 3 is gene_id (name field from our BED6)
    gene_id = df.iloc[:, 3].astype(str)
    count = df.iloc[:, 6].astype(int)
    return pd.Series(count.values, index=gene_id.values)


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else sys.stdout

    def log(msg: str) -> None:
        print(msg, file=log_fh, flush=True)

    gene_df = pd.read_csv(args.gene_table, sep="\t")

    g4_counts = read_intersect_counts(args.g4_intersect)
    gc_counts = read_intersect_counts(args.gc_bg_intersect)

    gene_df = gene_df.set_index("gene_id")
    gene_df["has_g4_tss"] = (g4_counts.reindex(gene_df.index).fillna(0) > 0).astype(int)
    gene_df["has_gc_bg_tss"] = (gc_counts.reindex(gene_df.index).fillna(0) > 0).astype(int)
    gene_df = gene_df.reset_index()

    both = (gene_df["has_g4_tss"] == 1) & (gene_df["has_gc_bg_tss"] == 1)
    log(f"Genes with BOTH G4 and GC-bg overlap (before precedence): {both.sum():,}")

    def assign_group(row) -> str:
        if row["has_g4_tss"] == 1:
            return "G4_TSS"
        if row["has_gc_bg_tss"] == 1:
            return "GC_bg_TSS"
        return "No_overlap"

    gene_df["group"] = gene_df.apply(assign_group, axis=1)

    counts = gene_df["group"].value_counts()
    log(f"Group counts after precedence:")
    for grp, n in counts.items():
        log(f"  {grp}: {n:,}")
    log(f"  Total genes: {len(gene_df):,}")

    out_cols = ["gene_id", "gene_name", "chrom", "tss", "strand",
                "group", "has_g4_tss", "has_gc_bg_tss"]
    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    gene_df[out_cols].to_csv(args.out_tsv, sep="\t", index=False)
    log(f"Written: {args.out_tsv}")

    if args.log:
        log_fh.close()


if __name__ == "__main__":
    main()
