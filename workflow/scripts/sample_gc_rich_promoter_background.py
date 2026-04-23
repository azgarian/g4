#!/usr/bin/env python3
"""Sample GC-rich canonical promoters that avoid G4 and OQS overlaps."""

from __future__ import annotations

import argparse
import csv
import gzip
import random
import re
from pathlib import Path


TARGET_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
TARGET_CHROM_SET = set(TARGET_CHROMS)
CHROM_ORDER = {chrom: idx for idx, chrom in enumerate(TARGET_CHROMS)}
NON_ACGT_RE = re.compile(r"[^ACGT]")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--promoter-bed", required=True, help="Canonical promoter BED6.")
    p.add_argument("--ref", required=True, help="Reference FASTA.")
    p.add_argument(
        "--oqs-bed",
        dest="oqs_beds",
        action="append",
        required=True,
        help="Prepared OQS BED to exclude. Repeat for multiple files.",
    )
    p.add_argument("--g4-bed", required=True, help="Merged G4 BED to exclude.")
    p.add_argument("--gc-threshold", type=float, default=0.28)
    p.add_argument("--sample-size", type=int, default=1000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out-sampled-bed", required=True)
    p.add_argument("--out-summary-tsv", required=True)
    return p.parse_args()


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals.sort()
    merged: list[tuple[int, int]] = []
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= cur_end:
            if end > cur_end:
                cur_end = end
            continue
        merged.append((cur_start, cur_end))
        cur_start, cur_end = start, end
    merged.append((cur_start, cur_end))
    return merged


def read_merged_bed(paths: list[Path]) -> dict[str, list[tuple[int, int]]]:
    interval_map: dict[str, list[tuple[int, int]]] = {chrom: [] for chrom in TARGET_CHROMS}
    for path in paths:
        if not path.exists():
            raise FileNotFoundError(f"Missing BED file: {path}")
        with open_text(path) as handle:
            for line_no, line in enumerate(handle, start=1):
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                fields = s.split("\t")
                if len(fields) < 3:
                    raise ValueError(f"{path}:{line_no}: expected at least 3 BED columns.")
                chrom = str(fields[0]).strip()
                if chrom not in TARGET_CHROM_SET:
                    continue
                try:
                    start = int(fields[1])
                    end = int(fields[2])
                except ValueError as exc:
                    raise ValueError(f"{path}:{line_no}: invalid BED coordinates.") from exc
                if end <= start:
                    continue
                interval_map[chrom].append((start, end))
    return {chrom: merge_intervals(intervals) for chrom, intervals in interval_map.items()}


def read_promoter_bed(
    path: Path,
) -> tuple[dict[str, list[tuple[int, int, str, str]]], int]:
    promoters_by_chrom: dict[str, list[tuple[int, int, str, str]]] = {chrom: [] for chrom in TARGET_CHROMS}
    n_total = 0
    with path.open("r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            fields = s.split("\t")
            if len(fields) < 6:
                raise ValueError(f"{path}:{line_no}: expected BED6 promoter rows.")
            chrom = str(fields[0]).strip()
            if chrom not in TARGET_CHROM_SET:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_no}: invalid promoter coordinates.") from exc
            if end <= start:
                continue
            gene_id = str(fields[3]).strip()
            strand = str(fields[5]).strip()
            promoters_by_chrom[chrom].append((start, end, gene_id, strand))
            n_total += 1
    for chrom in TARGET_CHROMS:
        promoters_by_chrom[chrom].sort(key=lambda row: (row[0], row[1], row[2]))
    return promoters_by_chrom, n_total


def has_interval_overlap(
    intervals: list[tuple[int, int]],
    start: int,
    end: int,
    interval_idx: int,
) -> tuple[bool, int]:
    while interval_idx < len(intervals) and intervals[interval_idx][1] <= start:
        interval_idx += 1
    if interval_idx < len(intervals):
        cur_start, cur_end = intervals[interval_idx]
        if cur_start < end and cur_end > start:
            return True, interval_idx
    return False, interval_idx


def iter_fasta_chromosomes(fasta_path: Path):
    chrom: str | None = None
    chunks: list[str] = []
    with open_text(fasta_path) as handle:
        for line in handle:
            if line.startswith(">"):
                if chrom in TARGET_CHROM_SET:
                    yield chrom, "".join(chunks).upper()
                chrom = line[1:].strip().split()[0]
                chunks = []
                continue
            if chrom in TARGET_CHROM_SET:
                chunks.append(line.strip())
    if chrom in TARGET_CHROM_SET:
        yield chrom, "".join(chunks).upper()


def select_gc_rich_promoters(
    *,
    promoters_by_chrom: dict[str, list[tuple[int, int, str, str]]],
    fasta_path: Path,
    oqs_map: dict[str, list[tuple[int, int]]],
    g4_map: dict[str, list[tuple[int, int]]],
    gc_threshold: float,
) -> tuple[list[tuple[str, int, int, str, float, str]], dict[str, int]]:
    eligible: list[tuple[str, int, int, str, float, str]] = []
    counts = {
        "n_promoters_oqs_only_dropped": 0,
        "n_promoters_g4_only_dropped": 0,
        "n_promoters_oqs_and_g4_dropped": 0,
        "n_promoters_non_acgt_dropped": 0,
        "n_promoters_exclusion_pass": 0,
        "n_promoters_gc_pass": 0,
        "n_promoters_sequence_checked": 0,
    }

    for chrom, chrom_seq in iter_fasta_chromosomes(fasta_path):
        promoters = promoters_by_chrom.get(chrom, [])
        if not promoters:
            continue
        oqs_intervals = oqs_map.get(chrom, [])
        g4_intervals = g4_map.get(chrom, [])
        oqs_idx = 0
        g4_idx = 0
        chrom_len = len(chrom_seq)
        for start, end, gene_id, strand in promoters:
            hit_oqs, oqs_idx = has_interval_overlap(oqs_intervals, start, end, oqs_idx)
            hit_g4, g4_idx = has_interval_overlap(g4_intervals, start, end, g4_idx)
            if hit_oqs or hit_g4:
                if hit_oqs and hit_g4:
                    counts["n_promoters_oqs_and_g4_dropped"] += 1
                elif hit_oqs:
                    counts["n_promoters_oqs_only_dropped"] += 1
                else:
                    counts["n_promoters_g4_only_dropped"] += 1
                continue

            counts["n_promoters_exclusion_pass"] += 1
            counts["n_promoters_sequence_checked"] += 1
            if end > chrom_len:
                raise ValueError(
                    f"Promoter interval {chrom}:{start}-{end} extends past FASTA chromosome length."
                )
            seq = chrom_seq[start:end]
            if NON_ACGT_RE.search(seq) is not None:
                counts["n_promoters_non_acgt_dropped"] += 1
                continue
            gc_frac = float(seq.count("G") + seq.count("C")) / float(len(seq))
            if gc_frac <= gc_threshold:
                continue
            counts["n_promoters_gc_pass"] += 1
            eligible.append((chrom, start, end, gene_id, gc_frac, strand))

    counts["n_promoters_oqs_dropped"] = (
        counts["n_promoters_oqs_only_dropped"] + counts["n_promoters_oqs_and_g4_dropped"]
    )
    counts["n_promoters_g4_dropped"] = (
        counts["n_promoters_g4_only_dropped"] + counts["n_promoters_oqs_and_g4_dropped"]
    )
    counts["n_promoters_excluded_any"] = (
        counts["n_promoters_oqs_only_dropped"]
        + counts["n_promoters_g4_only_dropped"]
        + counts["n_promoters_oqs_and_g4_dropped"]
    )
    return eligible, counts


def write_bed(rows: list[tuple[str, int, int, str, float, str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for chrom, start, end, name, score, strand in rows:
            writer.writerow([chrom, start, end, name, f"{score:.6f}", strand])


def write_summary(out_path: Path, summary: dict[str, int | float | str]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = list(summary.keys())
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(summary)


def main() -> None:
    args = parse_args()
    if not (0.0 <= float(args.gc_threshold) <= 1.0):
        raise ValueError("--gc-threshold must be between 0 and 1.")
    if int(args.sample_size) <= 0:
        raise ValueError("--sample-size must be positive.")

    promoter_path = Path(args.promoter_bed)
    ref_path = Path(args.ref)
    oqs_paths = [Path(path) for path in args.oqs_beds]
    g4_path = Path(args.g4_bed)

    if not promoter_path.exists():
        raise FileNotFoundError(f"Missing promoter BED: {promoter_path}")
    if not ref_path.exists():
        raise FileNotFoundError(f"Missing reference FASTA: {ref_path}")

    promoters_by_chrom, n_promoters_total = read_promoter_bed(promoter_path)
    oqs_map = read_merged_bed(oqs_paths)
    g4_map = read_merged_bed([g4_path])

    eligible_promoters, counts = select_gc_rich_promoters(
        promoters_by_chrom=promoters_by_chrom,
        fasta_path=ref_path,
        oqs_map=oqs_map,
        g4_map=g4_map,
        gc_threshold=float(args.gc_threshold),
    )

    requested = int(args.sample_size)
    if len(eligible_promoters) < requested:
        raise ValueError(
            "GC-rich promoter pool is smaller than the requested sample size: "
            f"{len(eligible_promoters):,} < {requested:,}."
        )

    rng = random.Random(int(args.seed))
    sampled = rng.sample(eligible_promoters, requested)
    sampled.sort(key=lambda row: (CHROM_ORDER[row[0]], row[1], row[2], row[3]))

    write_bed(sampled, Path(args.out_sampled_bed))
    summary = {
        "promoter_bed": str(promoter_path),
        "gc_threshold": float(args.gc_threshold),
        "sample_size_requested": requested,
        "seed": int(args.seed),
        "n_promoters_total": int(n_promoters_total),
        "n_promoters_oqs_only_dropped": int(counts["n_promoters_oqs_only_dropped"]),
        "n_promoters_g4_only_dropped": int(counts["n_promoters_g4_only_dropped"]),
        "n_promoters_oqs_and_g4_dropped": int(counts["n_promoters_oqs_and_g4_dropped"]),
        "n_promoters_oqs_dropped": int(counts["n_promoters_oqs_dropped"]),
        "n_promoters_g4_dropped": int(counts["n_promoters_g4_dropped"]),
        "n_promoters_excluded_any": int(counts["n_promoters_excluded_any"]),
        "n_promoters_exclusion_pass": int(counts["n_promoters_exclusion_pass"]),
        "n_promoters_non_acgt_dropped": int(counts["n_promoters_non_acgt_dropped"]),
        "n_promoters_gc_pass": int(counts["n_promoters_gc_pass"]),
        "n_promoters_sampled": len(sampled),
    }
    write_summary(Path(args.out_summary_tsv), summary)


if __name__ == "__main__":
    main()
