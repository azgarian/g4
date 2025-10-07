import argparse
from pathlib import Path
import pandas as pd

try:
    import bioframe as bf
except ImportError:
    raise SystemExit("bioframe is required; run within env with bioframe installed")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--g4", required=True, help="Input G4 BED file (BED3+) to annotate")
    ap.add_argument("--promoters", required=True, help="Promoters BED file (e.g., results/g4_txn/promoters.bed.gz)")
    ap.add_argument("--out", required=True, help="Output TSV with G4 and on_promoter flag")
    args = ap.parse_args()

    g4 = bf.read_table(args.g4, schema="bed")
    prom = bf.read_table(args.promoters, schema="bed")
    g4 = bf.sort_bedframe(g4)
    prom = bf.sort_bedframe(prom)

    # Compute overlaps; label whether each G4 overlaps any promoter
    ov = bf.overlap(g4, prom, how="left", suffixes=("", "_prom"))
    # Rows with NA in promoter coordinates mean no overlap
    on_prom = ~ov["start_prom"].isna()

    # Build output table; include identifying columns from input
    id_cols = [c for c in ["chrom", "start", "end", "name", "score", "strand"] if c in g4.columns]
    out_df = ov[id_cols].copy()
    out_df["on_promoter"] = on_prom.values

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()


