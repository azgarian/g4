#!/usr/bin/env python3
"""Sample a GC-rich genomic background and remove OQS-overlapping loci."""

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
    p.add_argument("--ref", required=True, help="Reference FASTA.")
    p.add_argument("--genome-fai", required=True, help="Genome FASTA index.")
    p.add_argument("--blacklist-bed", required=True, help="Blacklist BED.")
    p.add_argument(
        "--oqs-bed",
        dest="oqs_beds",
        action="append",
        required=True,
        help="Prepared OQS BED to exclude after sampling. Repeat for multiple files.",
    )
    p.add_argument("--tile-size", type=int, default=150)
    p.add_argument("--gc-threshold", type=float, default=0.28)
    p.add_argument("--sample-size", type=int, default=300000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out-sampled-bed", required=True)
    p.add_argument("--out-no-oqs-bed", required=True)
    p.add_argument("--out-summary-tsv", required=True)
    return p.parse_args()


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_fai(path: Path) -> dict[str, int]:
    chrom_sizes: dict[str, int] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_no}: expected at least 2 columns.")
            chrom = str(fields[0]).strip()
            if chrom not in TARGET_CHROM_SET:
                continue
            try:
                chrom_sizes[chrom] = int(fields[1])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_no}: invalid chromosome size '{fields[1]}'.") from exc
    missing = [chrom for chrom in TARGET_CHROMS if chrom not in chrom_sizes]
    if missing:
        raise ValueError(f"{path} is missing required chromosomes: {', '.join(missing)}")
    return chrom_sizes


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


def iter_target_tiles(fasta_path: Path, tile_size: int):
    chrom: str | None = None
    tile_start = 0
    tile_chunks: list[str] = []
    tile_len = 0

    with open_text(fasta_path) as handle:
        for line in handle:
            if line.startswith(">"):
                chrom = line[1:].strip().split()[0]
                if chrom not in TARGET_CHROM_SET:
                    chrom = None
                tile_start = 0
                tile_chunks = []
                tile_len = 0
                continue

            if chrom is None:
                continue

            seq = line.strip().upper()
            if not seq:
                continue

            idx = 0
            while idx < len(seq):
                need = tile_size - tile_len
                take = min(need, len(seq) - idx)
                tile_chunks.append(seq[idx : idx + take])
                tile_len += take
                idx += take
                if tile_len == tile_size:
                    tile_seq = "".join(tile_chunks)
                    yield chrom, tile_start, tile_start + tile_size, tile_seq
                    tile_start += tile_size
                    tile_chunks = []
                    tile_len = 0


def reservoir_sample_tiles(
    *,
    fasta_path: Path,
    blacklist_map: dict[str, list[tuple[int, int]]],
    chrom_sizes: dict[str, int],
    tile_size: int,
    gc_threshold: float,
    sample_size: int,
    seed: int,
) -> tuple[list[tuple[str, int, int, float]], dict[str, int]]:
    rng = random.Random(int(seed))
    reservoir: list[tuple[str, int, int, float]] = []
    counts = {
        "n_tiles_scanned": 0,
        "n_tiles_blacklist_dropped": 0,
        "n_tiles_non_acgt_dropped": 0,
        "n_tiles_gc_pass": 0,
    }
    interval_idx_by_chrom = {chrom: 0 for chrom in TARGET_CHROMS}

    for chrom, start, end, tile_seq in iter_target_tiles(fasta_path, tile_size):
        chrom_len = int(chrom_sizes[chrom])
        if end > chrom_len:
            continue

        counts["n_tiles_scanned"] += 1

        intervals = blacklist_map.get(chrom, [])
        interval_idx = interval_idx_by_chrom[chrom]
        while interval_idx < len(intervals) and intervals[interval_idx][1] <= start:
            interval_idx += 1
        interval_idx_by_chrom[chrom] = interval_idx
        if interval_idx < len(intervals):
            blk_start, blk_end = intervals[interval_idx]
            if blk_start < end and blk_end > start:
                counts["n_tiles_blacklist_dropped"] += 1
                continue

        if NON_ACGT_RE.search(tile_seq) is not None:
            counts["n_tiles_non_acgt_dropped"] += 1
            continue

        gc_frac = float(tile_seq.count("G") + tile_seq.count("C")) / float(tile_size)
        if gc_frac <= gc_threshold:
            continue

        counts["n_tiles_gc_pass"] += 1
        seen = counts["n_tiles_gc_pass"]
        record = (chrom, start, end, gc_frac)
        if len(reservoir) < sample_size:
            reservoir.append(record)
            continue

        replace_idx = rng.randrange(seen)
        if replace_idx < sample_size:
            reservoir[replace_idx] = record

    if counts["n_tiles_gc_pass"] < sample_size:
        raise ValueError(
            "GC-rich tile pool is smaller than the requested sample size: "
            f"{counts['n_tiles_gc_pass']:,} < {sample_size:,}."
        )

    return reservoir, counts


def filter_sampled_against_oqs(
    sampled_tiles: list[tuple[str, int, int, float]],
    oqs_map: dict[str, list[tuple[int, int]]],
) -> tuple[list[tuple[str, int, int, float]], int]:
    retained: list[tuple[str, int, int, float]] = []
    dropped = 0
    tiles_by_chrom: dict[str, list[tuple[str, int, int, float]]] = {chrom: [] for chrom in TARGET_CHROMS}
    for tile in sampled_tiles:
        tiles_by_chrom[tile[0]].append(tile)

    for chrom in TARGET_CHROMS:
        tiles = sorted(tiles_by_chrom[chrom], key=lambda row: (row[1], row[2]))
        intervals = oqs_map.get(chrom, [])
        interval_idx = 0
        for record in tiles:
            _, start, end, _ = record
            while interval_idx < len(intervals) and intervals[interval_idx][1] <= start:
                interval_idx += 1
            if interval_idx < len(intervals):
                oqs_start, oqs_end = intervals[interval_idx]
                if oqs_start < end and oqs_end > start:
                    dropped += 1
                    continue
            retained.append(record)

    retained.sort(key=lambda row: (CHROM_ORDER[row[0]], row[1], row[2]))
    return retained, dropped


def assign_names(records: list[tuple[str, int, int, float]]) -> list[tuple[str, int, int, str, float, str]]:
    named_rows: list[tuple[str, int, int, str, float, str]] = []
    for idx, (chrom, start, end, gc_frac) in enumerate(records, start=1):
        named_rows.append((chrom, start, end, f"gc_rich_bg_{idx}", gc_frac, "."))
    return named_rows


def write_bed(rows: list[tuple[str, int, int, str, float, str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for chrom, start, end, name, score, strand in rows:
            writer.writerow([chrom, start, end, name, f"{score:.6f}", strand])


def write_summary(out_path: Path, summary: dict[str, int | float]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = list(summary.keys())
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(summary)


def main() -> None:
    args = parse_args()
    if int(args.tile_size) <= 0:
        raise ValueError("--tile-size must be positive.")
    if not (0.0 <= float(args.gc_threshold) <= 1.0):
        raise ValueError("--gc-threshold must be between 0 and 1.")
    if int(args.sample_size) <= 0:
        raise ValueError("--sample-size must be positive.")

    ref_path = Path(args.ref)
    fai_path = Path(args.genome_fai)
    blacklist_path = Path(args.blacklist_bed)
    oqs_paths = [Path(path) for path in args.oqs_beds]

    if not ref_path.exists():
        raise FileNotFoundError(f"Missing reference FASTA: {ref_path}")
    chrom_sizes = read_fai(fai_path)
    blacklist_map = read_merged_bed([blacklist_path])
    oqs_map = read_merged_bed(oqs_paths)

    sampled_tiles, counts = reservoir_sample_tiles(
        fasta_path=ref_path,
        blacklist_map=blacklist_map,
        chrom_sizes=chrom_sizes,
        tile_size=int(args.tile_size),
        gc_threshold=float(args.gc_threshold),
        sample_size=int(args.sample_size),
        seed=int(args.seed),
    )
    sampled_tiles.sort(key=lambda row: (CHROM_ORDER[row[0]], row[1], row[2]))

    retained_tiles, oqs_dropped = filter_sampled_against_oqs(sampled_tiles, oqs_map)

    sampled_named = assign_names(sampled_tiles)
    sampled_no_oqs_named = assign_names(retained_tiles)

    write_bed(sampled_named, Path(args.out_sampled_bed))
    write_bed(sampled_no_oqs_named, Path(args.out_no_oqs_bed))

    summary = {
        "tile_size_bp": int(args.tile_size),
        "gc_threshold": float(args.gc_threshold),
        "sample_size_requested": int(args.sample_size),
        "seed": int(args.seed),
        "n_tiles_scanned": int(counts["n_tiles_scanned"]),
        "n_tiles_blacklist_dropped": int(counts["n_tiles_blacklist_dropped"]),
        "n_tiles_non_acgt_dropped": int(counts["n_tiles_non_acgt_dropped"]),
        "n_tiles_gc_pass": int(counts["n_tiles_gc_pass"]),
        "n_tiles_sampled": len(sampled_named),
        "n_tiles_oqs_dropped": int(oqs_dropped),
        "n_tiles_oqs_retained": len(sampled_no_oqs_named),
    }
    write_summary(Path(args.out_summary_tsv), summary)


if __name__ == "__main__":
    main()
