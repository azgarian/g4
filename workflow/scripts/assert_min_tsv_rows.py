#!/usr/bin/env python3
"""Fail if a TSV has fewer than the required number of data rows."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tsv", required=True, help="Input TSV file with a header row.")
    parser.add_argument("--min-rows", type=int, required=True, help="Minimum number of data rows required.")
    parser.add_argument("--label", default="table", help="Short label used in error messages.")
    return parser.parse_args()


def count_data_rows(path: Path) -> int:
    if not path.exists():
        raise FileNotFoundError(f"Missing TSV file: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            next(reader)
        except StopIteration as exc:
            raise ValueError(f"{path} is empty.") from exc

        return sum(1 for row in reader if row and any(field.strip() for field in row))


def main() -> None:
    args = parse_args()
    min_rows = int(args.min_rows)
    if min_rows < 0:
        raise ValueError("--min-rows must be non-negative.")

    tsv_path = Path(args.tsv)
    n_rows = count_data_rows(tsv_path)
    if n_rows < min_rows:
        raise ValueError(
            f"{args.label} has too few rows: observed {n_rows}, required at least {min_rows}."
        )

    print(f"{args.label}: {n_rows} rows >= {min_rows}")


if __name__ == "__main__":
    main()
