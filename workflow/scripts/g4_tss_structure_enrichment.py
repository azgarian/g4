#!/usr/bin/env python3
"""Task 7: Test G4 structure-class enrichment across baseline expression classes.

For G4_TSS genes, intersects canonical TSS windows with G4 ChIP and CUT&Tag structure
annotations, deduplicates at gene level, and tests enrichment of each structure class
in Q4 vs Q1 and silent genes.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tss-annotation", required=True)
    p.add_argument("--tss-windows-bed", required=True)
    p.add_argument("--baseline-tpm", required=True)
    p.add_argument("--g4chip-tsv", required=True)
    p.add_argument("--g4cuttag-tsv", required=True)
    p.add_argument("--analysis-label", default=None)
    p.add_argument("--out-structure-by-expression", required=True)
    p.add_argument("--out-plot", required=True)
    p.add_argument("--out-enrichment-stats", required=True)
    p.add_argument("--log", default=None)
    return p.parse_args()


def normalize_gene_id(gene_id: str) -> str:
    """Drop Ensembl version suffixes so tables can be joined safely."""
    return str(gene_id).split(".", 1)[0]


def tsv_to_bed6(tsv_path: str, tmp_bed: str) -> None:
    """Write a temporary BED6 from the prepared TSV structure file."""
    df = pd.read_csv(tsv_path, sep="\t")
    bed = df[["chrom", "start", "end", "name", "source_signal", "strand"]].copy()
    bed.to_csv(tmp_bed, sep="\t", header=False, index=False)


def intersect_structures(windows_bed: str, struct_tsv: str) -> pd.DataFrame:
    """Return DataFrame: gene_id, structure — one row per (gene, structure) pair."""
    import tempfile
    with tempfile.NamedTemporaryFile(suffix=".bed", mode="w", delete=False) as tmp:
        tmp_path = tmp.name
    tsv_to_bed6(struct_tsv, tmp_path)

    result = subprocess.run(
        ["bedtools", "intersect", "-nonamecheck", "-wa", "-wb", "-a", windows_bed, "-b", tmp_path],
        capture_output=True, text=True, check=True
    )
    Path(tmp_path).unlink(missing_ok=True)

    rows = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        gene_id = fields[3]  # col 4 of -a = gene_id
        # -b BED6: col 7 is name, but we need structure from TSV
        # bedtools -wb appends all -b columns; structure col offset depends on TSV
        # TSV structure-annotated file has: chrom,start,end,name,signal,strand = 6 cols
        # plus structure_start,structure_end,structure_center,structure,source_dataset
        # After bedtools intersect -wa -wb the -b cols start at index 6 (0-based)
        # name is at 6+3=9, but we wrote BED6 to tmp so only 6 -b cols
        # we need structure from the original TSV — approach: use name to map back
        rows.append({"gene_id": gene_id, "peak_name": fields[9] if len(fields) > 9 else ""})

    if not rows:
        return pd.DataFrame(columns=["gene_id", "peak_name"])
    return pd.DataFrame(rows)


def intersect_structures_full(windows_bed: str, struct_tsv: str) -> pd.DataFrame:
    """Return DataFrame: gene_id, structure using full TSV with all columns."""
    import tempfile

    df_tsv = pd.read_csv(struct_tsv, sep="\t")

    # write full TSV as tab-separated for bedtools (use first 6 cols for BED, keep structure)
    with tempfile.NamedTemporaryFile(suffix=".bed", mode="w", delete=False) as tmp:
        tmp_path = tmp.name
    # write extended BED: chrom start end name signal strand structure
    ext = df_tsv[["chrom", "start", "end", "name", "source_signal",
                   "strand", "structure"]].copy()
    ext.to_csv(tmp_path, sep="\t", header=False, index=False)

    result = subprocess.run(
        ["bedtools", "intersect", "-nonamecheck", "-wa", "-wb", "-a", windows_bed, "-b", tmp_path],
        capture_output=True, text=True, check=True
    )
    Path(tmp_path).unlink(missing_ok=True)

    rows = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        gene_id = fields[3]
        # -b has 7 cols (extended BED), structure is at index 6+6=12
        structure = fields[12] if len(fields) > 12 else "other"
        rows.append({"gene_id": gene_id, "structure": structure})

    if not rows:
        return pd.DataFrame(columns=["gene_id", "structure"])
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    log_fh = open(args.log, "w") if args.log else None

    def log(msg: str) -> None:
        if log_fh:
            print(msg, file=log_fh, flush=True)

    if args.analysis_label:
        log(f"Analysis label: {args.analysis_label}")

    anno = pd.read_csv(args.tss_annotation, sep="\t")
    tpm = pd.read_csv(args.baseline_tpm, sep="\t")

    # keep only G4_TSS genes
    g4_genes = set(
        anno.loc[anno["group"] == "G4_TSS", "gene_id"].map(normalize_gene_id).tolist()
    )
    log(f"G4_TSS genes: {len(g4_genes):,}")

    for outpath in [args.out_structure_by_expression, args.out_plot,
                    args.out_enrichment_stats]:
        Path(outpath).parent.mkdir(parents=True, exist_ok=True)

    # intersect TSS windows with structure TSVs
    chip_hits = intersect_structures_full(args.tss_windows_bed, args.g4chip_tsv)
    ct_hits = intersect_structures_full(args.tss_windows_bed, args.g4cuttag_tsv)
    all_hits = pd.concat([chip_hits, ct_hits], ignore_index=True)
    all_hits["gene_id_norm"] = all_hits["gene_id"].map(normalize_gene_id)
    all_hits = all_hits[all_hits["gene_id_norm"].isin(g4_genes)]

    # deduplicate: one occurrence per (gene_id, structure)
    all_hits = all_hits.drop_duplicates(subset=["gene_id_norm", "structure"])
    log(f"Unique (gene, structure) pairs: {len(all_hits):,}")

    # join with expression class
    tpm_sub = tpm[["gene_id", "expression_class"]].copy()
    tpm_sub["gene_id"] = tpm_sub["gene_id"].map(normalize_gene_id)

    merged = all_hits.merge(
        tpm_sub,
        left_on="gene_id_norm",
        right_on="gene_id",
        how="left",
        suffixes=("", "_tpm"),
    )

    # per structure class × expression class counts
    CLASS_ORDER = ["silent", "Q1", "Q2", "Q3", "Q4"]
    structure_classes = [s for s in merged["structure"].dropna().unique() if s != "other"]
    g4_tpm = tpm_sub[tpm_sub["gene_id"].isin(g4_genes)].copy()
    totals_by_class = g4_tpm["expression_class"].value_counts().to_dict()

    rows = []
    for struct in sorted(structure_classes):
        sub = merged[merged["structure"] == struct]
        for cls in CLASS_ORDER:
            n_genes = (sub["expression_class"] == cls).sum()
            # fraction among all G4_TSS genes in that expression class
            total_cls = totals_by_class.get(cls, 0)
            fraction = n_genes / total_cls if total_cls > 0 else np.nan
            rows.append({
                "structure": struct,
                "expression_class": cls,
                "n_genes": int(n_genes),
                "fraction_in_class": fraction,
            })

    struct_df = pd.DataFrame(rows)
    struct_df.to_csv(args.out_structure_by_expression, sep="\t", index=False)

    # --- Fisher exact: Q4 vs Q1 and Q4 vs silent ---
    enrich_rows = []
    all_g4_tss = g4_tpm
    n_q4 = (all_g4_tss["expression_class"] == "Q4").sum()
    n_q1 = (all_g4_tss["expression_class"] == "Q1").sum()
    n_sil = (all_g4_tss["expression_class"] == "silent").sum()

    for struct in sorted(structure_classes):
        sub = merged[merged["structure"] == struct]
        genes_with_struct = set(sub["gene_id_norm"].dropna().unique())

        for ref_cls, n_ref in [("Q1", n_q1), ("silent", n_sil)]:
            n_q4_struct = sum(
                1 for g in genes_with_struct
                if (tpm_sub[tpm_sub["gene_id"] == g]["expression_class"] == "Q4").any()
            )
            n_ref_struct = sum(
                1 for g in genes_with_struct
                if (tpm_sub[tpm_sub["gene_id"] == g]["expression_class"] == ref_cls).any()
            )
            table = [
                [n_q4_struct, n_q4 - n_q4_struct],
                [n_ref_struct, n_ref - n_ref_struct],
            ]
            if table[0][1] < 0 or table[1][1] < 0:
                continue
            try:
                or_val, p_fish = stats.fisher_exact(table, alternative="greater")
            except Exception:
                or_val, p_fish = np.nan, np.nan
            enrich_rows.append({
                "structure": struct,
                "comparison": f"Q4_vs_{ref_cls}",
                "n_q4_with_struct": n_q4_struct,
                "n_ref_with_struct": n_ref_struct,
                "odds_ratio": or_val,
                "p_value": p_fish,
            })

    if enrich_rows:
        enrich_df = pd.DataFrame(enrich_rows)
        p_vals = enrich_df["p_value"].fillna(1).values
        _, p_adj, _, _ = multipletests(p_vals, method="fdr_bh")
        enrich_df["p_adj_BH"] = p_adj
        enrich_df.to_csv(args.out_enrichment_stats, sep="\t", index=False)
        log(f"Written: {args.out_enrichment_stats}")

    # --- Stacked bar chart ---
    pivot = struct_df.pivot(index="structure", columns="expression_class",
                            values="fraction_in_class").reindex(columns=CLASS_ORDER)
    if pivot.empty:
        log("No structure data to plot.")
    else:
        fig, ax = plt.subplots(figsize=(max(6, len(structure_classes) * 0.8 + 2), 4.5))
        pivot.plot(kind="bar", ax=ax, colormap="RdYlBu_r", width=0.75)
        ax.set_xlabel("G4 structure class")
        ax.set_ylabel("Fraction of G4_TSS genes in class")
        title_prefix = f"{args.analysis_label}: " if args.analysis_label else ""
        ax.set_title(f"{title_prefix}G4 structure class distribution across expression classes")
        ax.legend(title="Expression class", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
        plt.xticks(rotation=45, ha="right", fontsize=8)
        plt.tight_layout()
        fig.savefig(args.out_plot, dpi=150)
        plt.close(fig)

    if log_fh:
        log_fh.close()


if __name__ == "__main__":
    main()
