#!/usr/bin/env python3
"""Filter merged G4 intervals by whether their center falls inside ATAC peaks."""

from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--g4-bed", required=True)
    p.add_argument("--atac-bed", required=True)
    p.add_argument("--analysis-label", default=None)
    p.add_argument("--out-bed", required=True)
    p.add_argument("--out-summary", required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    Path(args.out_bed).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_summary).parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as centers_fh:
        centers_path = Path(centers_fh.name)
        with open(args.g4_bed, "r", encoding="utf-8") as fh:
            for idx, line in enumerate(fh):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                chrom, start_str, end_str = fields[:3]
                start = int(start_str)
                end = int(end_str)
                center = (start + end) // 2
                centers_fh.write(f"{chrom}\t{center}\t{center + 1}\tg4_{idx}\n")

    try:
        result = subprocess.run(
            [
                "bedtools",
                "intersect",
                "-nonamecheck",
                "-u",
                "-a",
                str(centers_path),
                "-b",
                args.atac_bed,
            ],
            capture_output=True,
            text=True,
            check=True,
        )
    finally:
        centers_path.unlink(missing_ok=True)

    keep_ids = {
        line.split("\t")[3]
        for line in result.stdout.strip().splitlines()
        if line.strip()
    }

    total = 0
    retained = 0
    with open(args.out_bed, "w", encoding="utf-8") as out_fh, open(
        args.g4_bed, "r", encoding="utf-8"
    ) as in_fh:
        for idx, line in enumerate(in_fh):
            total += 1
            if f"g4_{idx}" in keep_ids:
                out_fh.write(line)
                retained += 1

    with open(args.out_summary, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "analysis_label",
                "g4_bed",
                "atac_bed",
                "filter_rule",
                "total_g4_intervals",
                "retained_g4_intervals",
                "retained_fraction",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "analysis_label": args.analysis_label or "",
                "g4_bed": args.g4_bed,
                "atac_bed": args.atac_bed,
                "filter_rule": "g4_center_in_atac_peak",
                "total_g4_intervals": total,
                "retained_g4_intervals": retained,
                "retained_fraction": f"{(retained / total) if total else 0:.6f}",
            }
        )

    print(f"Filtered {retained:,} of {total:,} G4 intervals using ATAC center overlap.")


if __name__ == "__main__":
    main()
