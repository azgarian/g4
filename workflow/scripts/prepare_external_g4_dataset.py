#!/usr/bin/env python3
"""Prepare external G4 datasets with structure-aware metadata."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from regex_structures_common import annotate_bed_with_structure, normalize_structure


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-bed", required=True)
    p.add_argument("--ref", required=True)
    p.add_argument("--source-dataset", required=True)
    p.add_argument("--mode", choices=["strandless", "fixed"], required=True)
    p.add_argument("--signal-column-index", type=int, required=True)
    p.add_argument("--name-column-index", type=int, default=3)
    p.add_argument("--strand-column-index", type=int, default=-1)
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--out-source-bed", required=True)
    p.add_argument("--out-center-bed", required=True)
    return p.parse_args()


def read_bedlike(
    path: Path,
    *,
    source_dataset: str,
    signal_idx: int,
    name_idx: int,
    strand_idx: int,
    mode: str,
) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing BED-like input: {path}")

    df = pd.read_csv(path, sep="\t", header=None, comment="#", low_memory=False)
    needed = max(2, signal_idx, name_idx, strand_idx)
    if df.shape[1] <= needed:
        raise ValueError(f"{path} is missing required columns up to index {needed}.")

    out = pd.DataFrame(
        {
            "chrom": df.iloc[:, 0].astype(str),
            "start": pd.to_numeric(df.iloc[:, 1], errors="raise").astype(np.int64),
            "end": pd.to_numeric(df.iloc[:, 2], errors="raise").astype(np.int64),
            "name": df.iloc[:, name_idx].astype(str).str.strip(),
            "source_signal": pd.to_numeric(df.iloc[:, signal_idx], errors="raise"),
        }
    )
    out = out[out["end"] > out["start"]].copy()
    if mode == "fixed":
        out["strand"] = df.iloc[:, strand_idx].astype(str).str.strip()
        bad = out[~out["strand"].isin({"+", "-"})]
        if not bad.empty:
            raise ValueError(f"{path} contains invalid fixed-strand values.")
    else:
        out["strand"] = "."

    empty_name = out["name"].isin({"", ".", "nan", "None"})
    if empty_name.any():
        out.loc[empty_name, "name"] = [f"{source_dataset}_{i}" for i in out.index[empty_name]]
    out["source_dataset"] = source_dataset
    return out.reset_index(drop=True)


def _finalize_structure_rows(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["structure"] = out["structure"].astype(str).map(normalize_structure)
    for col in ["structure_start", "structure_end", "structure_center"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    out = out[out["structure"] != "other"].copy()
    out = out.dropna(subset=["structure_start", "structure_end", "structure_center"]).copy()
    for col in ["start", "end", "structure_start", "structure_end", "structure_center"]:
        out[col] = out[col].astype(np.int64)
    out = out[out["structure_end"] > out["structure_start"]].copy()
    out = out[(out["structure_start"] >= out["start"]) & (out["structure_end"] <= out["end"])].copy()
    out = out[(out["structure_center"] >= out["structure_start"]) & (out["structure_center"] < out["structure_end"])].copy()
    out = out[out["strand"].isin({"+", "-"})].copy()
    cols = [
        "chrom",
        "start",
        "end",
        "name",
        "source_signal",
        "strand",
        "structure_start",
        "structure_end",
        "structure_center",
        "structure",
        "source_dataset",
    ]
    missing = [col for col in cols if col not in out.columns]
    if missing:
        raise ValueError(
            "Annotated structure rows are missing required columns: "
            f"{', '.join(missing)}."
        )
    return out[cols].sort_values(["chrom", "start", "end", "strand", "name"], kind="stable").reset_index(drop=True)


def annotate_fixed_strand(df: pd.DataFrame, ref_fasta: str) -> pd.DataFrame:
    annotated = annotate_bed_with_structure(
        bed_df=df[["chrom", "start", "end", "name", "source_signal", "strand", "source_dataset"]],
        genome_file=ref_fasta,
        midpoint_col_name="structure_center",
        structure_col_name="structure",
        match_start_col_name="structure_start",
        match_end_col_name="structure_end",
    )
    return _finalize_structure_rows(annotated)


def annotate_strandless(df: pd.DataFrame, ref_fasta: str) -> pd.DataFrame:
    plus_df = df.copy()
    plus_df["strand"] = "+"
    minus_df = df.copy()
    minus_df["strand"] = "-"

    plus_annot = annotate_bed_with_structure(
        bed_df=plus_df[["chrom", "start", "end", "name", "source_signal", "strand"]],
        genome_file=ref_fasta,
        midpoint_col_name="structure_center",
        structure_col_name="structure",
        match_start_col_name="structure_start",
        match_end_col_name="structure_end",
    )
    minus_annot = annotate_bed_with_structure(
        bed_df=minus_df[["chrom", "start", "end", "name", "source_signal", "strand"]],
        genome_file=ref_fasta,
        midpoint_col_name="structure_center",
        structure_col_name="structure",
        match_start_col_name="structure_start",
        match_end_col_name="structure_end",
    )

    plus_struct = plus_annot["structure"].astype(str).map(normalize_structure)
    minus_struct = minus_annot["structure"].astype(str).map(normalize_structure)
    choose_plus = plus_struct != "other"
    choose_minus = (~choose_plus) & (minus_struct != "other")

    out = df.copy()
    out["strand"] = "."
    out["structure"] = "other"
    for col in ["structure_start", "structure_end", "structure_center"]:
        out[col] = pd.NA

    for mask, strand, annot in [(choose_plus, "+", plus_annot), (choose_minus, "-", minus_annot)]:
        out.loc[mask, "strand"] = strand
        out.loc[mask, "structure"] = annot.loc[mask, "structure"].astype(str).to_numpy()
        out.loc[mask, "structure_start"] = annot.loc[mask, "structure_start"].to_numpy()
        out.loc[mask, "structure_end"] = annot.loc[mask, "structure_end"].to_numpy()
        out.loc[mask, "structure_center"] = annot.loc[mask, "structure_center"].to_numpy()

    return _finalize_structure_rows(out)


def write_bed_outputs(prepared_df: pd.DataFrame, source_bed: Path, center_bed: Path) -> None:
    source_bed.parent.mkdir(parents=True, exist_ok=True)
    center_bed.parent.mkdir(parents=True, exist_ok=True)

    source_df = prepared_df[["chrom", "start", "end", "name", "source_signal", "strand"]].copy()
    source_df.to_csv(source_bed, sep="\t", header=False, index=False)

    center_df = pd.DataFrame(
        {
            "chrom": prepared_df["chrom"],
            "start": prepared_df["structure_center"].astype(np.int64),
            "end": prepared_df["structure_center"].astype(np.int64) + 1,
            "name": prepared_df["name"],
            "source_signal": prepared_df["source_signal"],
            "strand": prepared_df["strand"],
        }
    )
    center_df.to_csv(center_bed, sep="\t", header=False, index=False)


def main() -> None:
    args = parse_args()
    input_df = read_bedlike(
        Path(args.input_bed),
        source_dataset=str(args.source_dataset).strip(),
        signal_idx=int(args.signal_column_index),
        name_idx=int(args.name_column_index),
        strand_idx=int(args.strand_column_index),
        mode=str(args.mode).strip(),
    )
    if args.mode == "strandless":
        prepared_df = annotate_strandless(input_df, str(args.ref))
    else:
        prepared_df = annotate_fixed_strand(input_df, str(args.ref))
    if prepared_df.empty:
        raise ValueError(f"No structure-supported loci remain for dataset '{args.source_dataset}'.")

    out_tsv = Path(args.out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    prepared_df.to_csv(out_tsv, sep="\t", index=False)
    write_bed_outputs(prepared_df, Path(args.out_source_bed), Path(args.out_center_bed))


if __name__ == "__main__":
    main()
