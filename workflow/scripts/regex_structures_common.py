#!/usr/bin/env python3

"""Shared helpers for regex-based G4 structure assignment."""

import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Optional, Tuple

import pandas as pd

STRUCTURE_ORDER = [
    "loop1_3",
    "loop4_5",
    "loop6_7",
    "longLoop",
    "simpleBulge",
    "complexBulge",
    "twoTetrads",
    "other",
]

STRUCTURE_LABELS = [
    "Loop 1-3",
    "Loop 4-5",
    "Loop 6-7",
    "Long loops",
    "Simple bulges",
    "Complex bulges",
    "2-Tetrads",
    "Others",
]

STRUCTURE_REGEX = {
    "loop1_3": [
        # PQS with short loops 1-3
        r"([gG]{3,}\w{1,3}){3,}[gG]{3,}",
    ],
    "loop4_5": [
        # PQS with short loops 4-5
        r"([gG]{3,}\w{4,5}){3,}[gG]{3,}",
    ],
    "loop6_7": [
        # PQS with short loops 6-7
        r"([gG]{3,}\w{6,7}){3,}[gG]{3,}",
    ],
    "longLoop": [
        # PQS with longer loops 1-12
        r"([gG]{3,}\w{1,12}){3,}[gG]{3,}",
        # PQS with very long central loop loops 8-21
        r"[gG]{3,}\w{1,7}[gG]{3,}\w{13,21}[gG]{3,}\w{1,7}[gG]{3,}",
    ],
    "simpleBulge": [
        # PQ with simple buldge first tetrad
        r"[gG]{1,}\w{1,7}[gG]{2,}\w{1,7}([gG]{3,}\w{1,7}){2,}[gG]{3,}",
        # PQ with simple buldge second tetrad
        r"[gG]{3,}\w{1,7}[gG]{1,}\w{1,7}[gG]{2,}\w{1,7}[gG]{3,}\w{1,7}[gG]{3,}",
        # PQ with simple buldge third tetrad
        r"([gG]{3,}\w{1,7}){2,}[gG]{1,}\w{1,7}[gG]{2,}\w{1,7}[gG]{3,}",
        # PQ with simple buldge fourth tetrad
        r"([gG]{3,}\w{1,7}){3,}[gG]{1,}\w{1,7}[gG]{2,}",
    ],
    "twoTetrads": [
        # PQS with 2 tetrads loops 1-12
        r'([gG]{2,}\w{1,7}){3,}[gG]{2,}',
    ],
    "complexBulge": [
        # PQ with complex buldge first and third tetrad
        r"[gG]{1,}\w{1,5}[gG]{2,}\w{1,7}[gG]{3,}\w{1,7}[gG]{1,}\w{1,5}[gG]{2,}\w{1,7}[gG]{3,}",
        # PQ with complex buldge second and fourth tetrad
        r"[gG]{3,}\w{1,7}[gG]{1,}\w{1,5}[gG]{2,}\w{1,7}[gG]{3,}\w{1,7}[gG]{1,}\w{1,5}[gG]{2,}",
    ],
}

COMPILED_STRUCTURE_REGEX = {
    name: [re.compile(pattern) for pattern in pattern_list]
    for name, pattern_list in STRUCTURE_REGEX.items()
}


def classify_sequence(
    sequence: str,
) -> Tuple[str, Optional[Tuple[int, int]]]:
    """Assign sequence class and return matched span in sequence coordinates."""

    for structure in STRUCTURE_ORDER[:-1]:
        for regex in COMPILED_STRUCTURE_REGEX[structure]:
            match = regex.search(sequence)
            if match is not None:
                return structure, match.span()
    return "other", None


def normalize_structure(value: str) -> str:
    """Normalize arbitrary structure labels into known classes."""

    return value if value in STRUCTURE_ORDER[:-1] else "other"


def bed_row_coords(row: pd.Series) -> Tuple[int, int, str]:
    """Return (start, end, strand) from a BED row."""

    start = int(float(row.iloc[1]))
    end = int(float(row.iloc[2]))
    strand = "."
    if len(row) > 5:
        value = str(row.iloc[5]).strip()
        if value:
            strand = value
    return start, end, strand


def bed_row_to_fasta_header(row: pd.Series) -> str:
    """Map a BED row to the default bedtools -s FASTA header format."""

    chrom = str(row.iloc[0])
    start, end, strand = bed_row_coords(row)

    return f"{chrom}:{start}-{end}({strand})"


def df_to_fasta_dict(df: pd.DataFrame, genome: str) -> Dict[str, str]:
    """Extract FASTA for BED coordinates and return header->sequence map."""

    if shutil.which("bedtools") is None:
        raise RuntimeError("bedtools not found on PATH; required for FASTA extraction.")

    tmp_root = Path("tmp/bedtools")
    tmp_root.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".bed",
        prefix="regex_",
        dir=str(tmp_root),
        delete=False,
    ) as bed_f, tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".fa",
        prefix="regex_",
        dir=str(tmp_root),
        delete=False,
    ) as fa_f:
        bed_path = Path(bed_f.name)
        fa_path = Path(fa_f.name)

    try:
        bed_df = pd.DataFrame()
        bed_df[0] = df.iloc[:, 0]
        bed_df[1] = pd.to_numeric(df.iloc[:, 1], errors="coerce")
        bed_df[2] = pd.to_numeric(df.iloc[:, 2], errors="coerce")
        # Keep strand-aware extraction; default to "." when missing.
        bed_df[3] = "."
        bed_df[4] = "."
        bed_df[5] = df.iloc[:, 5].astype(str) if df.shape[1] > 5 else "."
        bed_df = bed_df.dropna(subset=[1, 2]).copy()
        bed_df[1] = bed_df[1].astype(int)
        bed_df[2] = bed_df[2].astype(int)
        bed_df = bed_df[bed_df[2] > bed_df[1]]
        bed_df.to_csv(bed_path, sep="\t", header=False, index=False)

        cmd = ["bedtools", "getfasta", "-fi", str(genome), "-bed", str(bed_path), "-s", "-fo", str(fa_path)]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            detail = proc.stderr.strip() or proc.stdout.strip() or "no output"
            raise RuntimeError(f"bedtools getfasta failed: {detail}")

        fasta_sequences = fa_path.read_text(encoding="utf-8").split(">")
        fasta_dict: Dict[str, str] = {}
        for fasta in fasta_sequences:
            lines = fasta.strip().split("\n")
            if not lines or not lines[0]:
                continue
            header = lines[0]
            sequence = "".join(lines[1:])
            fasta_dict[header] = sequence
    finally:
        try:
            bed_path.unlink(missing_ok=True)
        except OSError:
            pass
        try:
            fa_path.unlink(missing_ok=True)
        except OSError:
            pass

    return fasta_dict


def annotate_bed_with_structure(
    bed_df: pd.DataFrame,
    genome_file: str,
    midpoint_col_name: str = "structure_midpoint",
    structure_col_name: str = "structure",
    match_start_col_name: Optional[str] = None,
    match_end_col_name: Optional[str] = None,
) -> pd.DataFrame:
    """Append structure metadata using regex-based sequence calls."""

    fasta_dict = df_to_fasta_dict(bed_df, genome_file)
    call_by_header = {
        header: classify_sequence(sequence)
        for header, sequence in fasta_dict.items()
    }

    assigned_midpoints = []
    assigned_match_starts = []
    assigned_match_ends = []
    assigned_structures = []
    for _, row in bed_df.iterrows():
        try:
            header = bed_row_to_fasta_header(row)
            structure, span = call_by_header.get(header, ("other", None))
            assigned_structures.append(structure)

            if span is None:
                assigned_midpoints.append(pd.NA)
                assigned_match_starts.append(pd.NA)
                assigned_match_ends.append(pd.NA)
                continue

            start, end, strand = bed_row_coords(row)
            span_start, span_end = span
            if strand == "-":
                match_start = end - span_end
                match_end = end - span_start
            else:
                match_start = start + span_start
                match_end = start + span_end

            assigned_midpoints.append((match_start + match_end) // 2)
            assigned_match_starts.append(match_start)
            assigned_match_ends.append(match_end)
        except Exception:
            assigned_midpoints.append(pd.NA)
            assigned_match_starts.append(pd.NA)
            assigned_match_ends.append(pd.NA)
            assigned_structures.append("other")

    annotated_df = bed_df.copy()
    annotated_df[midpoint_col_name] = assigned_midpoints
    if match_start_col_name:
        annotated_df[match_start_col_name] = assigned_match_starts
    if match_end_col_name:
        annotated_df[match_end_col_name] = assigned_match_ends
    annotated_df[structure_col_name] = assigned_structures
    return annotated_df
