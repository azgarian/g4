#!/usr/bin/env python3
"""Task 1: Extract one canonical TSS per protein-coding gene from GENCODE v35 GTF.

Selects the longest protein-coding transcript per gene and records its TSS.
Writes a BED6 file (canonical_tss.bed) and a gene_name side table.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gtf", required=True, help="GENCODE GTF (plain or gzipped)")
    p.add_argument("--out-bed", required=True, help="Output BED6 file of canonical TSS loci")
    p.add_argument("--out-gene-table", required=True, help="TSV with gene_id, gene_name, chrom, tss, strand")
    p.add_argument("--log", default=None, help="Log file path")
    return p.parse_args()


_ATTR_RE = re.compile(r'(\w+) "([^"]+)"')


def parse_attributes(attr_str: str) -> dict[str, str]:
    return dict(_ATTR_RE.findall(attr_str))


def open_gtf(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def main() -> None:
    args = parse_args()

    log_fh = open(args.log, "w") if args.log else sys.stdout

    def log(msg: str) -> None:
        print(msg, file=log_fh, flush=True)

    log("Parsing GENCODE GTF for protein-coding transcripts...")

    # transcript_id -> {gene_id, gene_name, chrom, strand, tx_start, tx_end, exon_length}
    transcripts: dict[str, dict] = {}
    # exon lengths accumulate per transcript
    exon_lengths: dict[str, int] = defaultdict(int)

    with open_gtf(args.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature not in ("transcript", "exon"):
                continue

            attrs = parse_attributes(fields[8])
            gene_type = attrs.get("gene_type", "")
            transcript_type = attrs.get("transcript_type", "")

            if gene_type != "protein_coding" or transcript_type != "protein_coding":
                continue

            gene_id = attrs.get("gene_id", "")
            transcript_id = attrs.get("transcript_id", "")
            gene_name = attrs.get("gene_name", gene_id)
            chrom = fields[0]
            start = int(fields[3])  # 1-based
            end = int(fields[4])    # 1-based inclusive
            strand = fields[6]

            if feature == "transcript":
                transcripts[transcript_id] = {
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "chrom": chrom,
                    "strand": strand,
                    "tx_start": start,
                    "tx_end": end,
                }
            elif feature == "exon":
                exon_lengths[transcript_id] += (end - start + 1)

    log(f"  Parsed {len(transcripts):,} protein-coding transcripts")

    # select longest transcript per gene
    gene_best: dict[str, tuple[int, str]] = {}  # gene_id -> (exon_length, transcript_id)
    for tx_id, info in transcripts.items():
        length = exon_lengths.get(tx_id, 0)
        gid = info["gene_id"]
        if gid not in gene_best or length > gene_best[gid][0]:
            gene_best[gid] = (length, tx_id)

    log(f"  Selected canonical transcripts for {len(gene_best):,} protein-coding genes")

    # build output rows
    rows = []
    for gid, (_, tx_id) in gene_best.items():
        info = transcripts[tx_id]
        chrom = info["chrom"]
        strand = info["strand"]
        # TSS: 0-based position
        if strand == "+":
            tss_0based = info["tx_start"] - 1
        else:
            tss_0based = info["tx_end"] - 1  # tx_end is 1-based inclusive

        rows.append({
            "chrom": chrom,
            "start": tss_0based,
            "end": tss_0based + 1,
            "gene_id": gid,
            "score": 0,
            "strand": strand,
            "gene_name": info["gene_name"],
            "tss": tss_0based,
        })

    df = pd.DataFrame(rows)

    # sort by chrom, start
    chrom_order = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y", "M"]]
    df["chrom_cat"] = pd.Categorical(df["chrom"], categories=chrom_order, ordered=True)
    df = df.sort_values(["chrom_cat", "start"]).drop(columns=["chrom_cat"])

    # write BED6
    Path(args.out_bed).parent.mkdir(parents=True, exist_ok=True)
    bed_cols = ["chrom", "start", "end", "gene_id", "score", "strand"]
    df[bed_cols].to_csv(args.out_bed, sep="\t", header=False, index=False)
    log(f"  Written BED6: {args.out_bed} ({len(df):,} rows)")

    # write gene_name side table
    Path(args.out_gene_table).parent.mkdir(parents=True, exist_ok=True)
    table_cols = ["gene_id", "gene_name", "chrom", "tss", "strand"]
    df[table_cols].to_csv(args.out_gene_table, sep="\t", index=False)
    log(f"  Written gene table: {args.out_gene_table}")

    if args.log:
        log_fh.close()


if __name__ == "__main__":
    main()
