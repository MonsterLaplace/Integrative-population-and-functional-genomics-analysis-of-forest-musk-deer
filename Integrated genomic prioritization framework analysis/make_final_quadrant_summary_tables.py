#!/usr/bin/env python3

"""
Create final all-quadrant summary tables:

1. final_quadrant_actual_windows_summary.tsv
   - actual window counts for all four quadrants

2. final_quadrant_unique_genes_summary.tsv
   - unique gene counts for all four quadrants

These are all-window summaries, not top1% summaries.
"""

import argparse
import os
import re
import pandas as pd


QUADRANT_ORDER = [
    "conservation_priority",
    "breeding_function_priority",
    "conflict_quadrant",
    "background",
]

QUADRANT_LABELS = {
    "conservation_priority": "Conservation priority",
    "breeding_function_priority": "Breeding/function priority",
    "conflict_quadrant": "Conflict quadrant",
    "background": "Background",
}

BACKGROUND_ALIASES = {
    "background_low_low",
    "background_quadrant",
    "background_region",
    "background",
    "other",
    "other_windows",
    "Other windows",
}


def normalize_quadrant(x):
    x = "" if pd.isna(x) else str(x).strip()
    if x in BACKGROUND_ALIASES:
        return "background"
    return x


def valid_gene(x):
    if pd.isna(x):
        return False
    x = str(x).strip()
    return x not in {"", ".", "NA", "NaN", "NULL", "None", "none"}


def split_genes(x):
    if not valid_gene(x):
        return []
    genes = []
    for g in re.split(r"[,;|]", str(x)):
        g = g.strip()
        if valid_gene(g):
            genes.append(g)
    return genes


def main():
    parser = argparse.ArgumentParser(description="Make all-quadrant actual-window and unique-gene summary tables.")
    parser.add_argument("--outdir", default="postprocess_results")
    parser.add_argument("--quadrant-table", default=None)
    parser.add_argument("--expanded-all", default=None)
    parser.add_argument("--window-out", default=None)
    parser.add_argument("--gene-out", default=None)
    args = parser.parse_args()

    quadrant_file = args.quadrant_table or os.path.join(args.outdir, "conflict_quadrant_table.tsv")
    expanded_file = args.expanded_all or os.path.join(args.outdir, "conflict_window_gene_expanded.tsv")
    window_out = args.window_out or os.path.join(args.outdir, "final_quadrant_actual_windows_summary.tsv")
    gene_out = args.gene_out or os.path.join(args.outdir, "final_quadrant_unique_genes_summary.tsv")

    qtab = pd.read_csv(
        quadrant_file,
        sep="\t",
        usecols=lambda c: c in {"CHROM", "start", "end", "quadrant"},
        low_memory=False,
    )
    qtab["quadrant"] = qtab["quadrant"].map(normalize_quadrant)
    qtab["start"] = pd.to_numeric(qtab["start"], errors="coerce")
    qtab["end"] = pd.to_numeric(qtab["end"], errors="coerce")
    qtab = qtab.dropna(subset=["CHROM", "start", "end"])
    qtab = qtab.drop_duplicates(["CHROM", "start", "end", "quadrant"])

    total_windows = qtab.drop_duplicates(["CHROM", "start", "end"]).shape[0]
    window_rows = []
    for q in QUADRANT_ORDER:
        n = qtab.loc[qtab["quadrant"] == q].drop_duplicates(["CHROM", "start", "end"]).shape[0]
        window_rows.append({
            "quadrant": q,
            "quadrant_label": QUADRANT_LABELS.get(q, q),
            "actual_windows": int(n),
            "percent_of_all_windows": (n / total_windows * 100) if total_windows else 0,
        })
    window_rows.append({
        "quadrant": "TOTAL",
        "quadrant_label": "TOTAL",
        "actual_windows": int(total_windows),
        "percent_of_all_windows": 100.0 if total_windows else 0,
    })
    window_df = pd.DataFrame(window_rows)
    window_df.to_csv(window_out, sep="\t", index=False)

    # Map each window to final quadrant, then join with all-window gene-expanded table.
    key_map = qtab[["CHROM", "start", "end", "quadrant"]].rename(columns={
        "CHROM": "window_CHROM",
        "start": "window_start",
        "end": "window_end",
    })
    key_map["window_start"] = pd.to_numeric(key_map["window_start"], errors="coerce")
    key_map["window_end"] = pd.to_numeric(key_map["window_end"], errors="coerce")
    key_map["window_CHROM"] = key_map["window_CHROM"].astype(str)

    expanded = pd.read_csv(
        expanded_file,
        sep="\t",
        usecols=lambda c: c in {"window_CHROM", "window_start", "window_end", "gene_id", "gene_std", "gene_input"},
        low_memory=False,
    )
    expanded["window_start"] = pd.to_numeric(expanded["window_start"], errors="coerce")
    expanded["window_end"] = pd.to_numeric(expanded["window_end"], errors="coerce")
    expanded["window_CHROM"] = expanded["window_CHROM"].astype(str)
    expanded = expanded.merge(key_map, on=["window_CHROM", "window_start", "window_end"], how="inner")

    gene_col = "gene_id" if "gene_id" in expanded.columns else ("gene_std" if "gene_std" in expanded.columns else "gene_input")
    expanded[gene_col] = expanded[gene_col].astype(str).str.strip()
    expanded = expanded[expanded[gene_col].map(valid_gene)].copy()

    total_unique_genes_global = expanded[gene_col].nunique()
    gene_rows = []
    for q in QUADRANT_ORDER:
        sub = expanded[expanded["quadrant"] == q]
        unique_genes = sub[gene_col].dropna().astype(str).str.strip()
        unique_genes = unique_genes[unique_genes.map(valid_gene)].unique()
        n = len(unique_genes)
        gene_rows.append({
            "quadrant": q,
            "quadrant_label": QUADRANT_LABELS.get(q, q),
            "unique_genes_within_quadrant": int(n),
            "percent_of_quadrant_level_unique_gene_sum": 0.0,
        })

    gene_df = pd.DataFrame(gene_rows)
    quadrant_gene_sum = int(gene_df["unique_genes_within_quadrant"].sum())
    if quadrant_gene_sum:
        gene_df["percent_of_quadrant_level_unique_gene_sum"] = (
            gene_df["unique_genes_within_quadrant"] / quadrant_gene_sum * 100
        )
    gene_df.loc[len(gene_df)] = {
        "quadrant": "TOTAL_QUADRANT_SUM",
        "quadrant_label": "TOTAL_QUADRANT_SUM",
        "unique_genes_within_quadrant": quadrant_gene_sum,
        "percent_of_quadrant_level_unique_gene_sum": 100.0 if quadrant_gene_sum else 0.0,
    }
    gene_df.loc[len(gene_df)] = {
        "quadrant": "GLOBAL_UNIQUE",
        "quadrant_label": "GLOBAL_UNIQUE",
        "unique_genes_within_quadrant": int(total_unique_genes_global),
        "percent_of_quadrant_level_unique_gene_sum": pd.NA,
    }
    gene_df.to_csv(gene_out, sep="\t", index=False)

    print("[OK] final all-quadrant summary tables written")
    print(f"  windows: {window_out}")
    print(f"  genes  : {gene_out}")
    print("\nWindow summary:")
    print(window_df.to_string(index=False))
    print("\nGene summary:")
    print(gene_df.to_string(index=False))


if __name__ == "__main__":
    main()
