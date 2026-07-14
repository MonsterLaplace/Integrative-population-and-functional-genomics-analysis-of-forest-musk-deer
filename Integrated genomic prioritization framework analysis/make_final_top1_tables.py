#!/usr/bin/env python3

"""
Create final top1% tables:

1. final_top1_actual_windows.tsv
   - final non-background quadrant top1% windows
   - class is synchronised to quadrant

2. final_top1_unique_genes.tsv
   - globally unique genes appearing in the final top1% windows
   - if one gene appears in multiple top1 windows/quadrants, it is represented
     once with summarised quadrant/window information
"""

import argparse
import os
import pandas as pd


VALID_QUADRANTS = [
    "conservation_priority",
    "breeding_function_priority",
    "conflict_quadrant",
]


def numeric_window_coords(df, chrom_col, start_col, end_col):
    df = df.copy()
    df[chrom_col] = df[chrom_col].astype(str)
    df[start_col] = pd.to_numeric(df[start_col], errors="coerce")
    df[end_col] = pd.to_numeric(df[end_col], errors="coerce")
    return df


def first_nonempty(series):
    for x in series:
        if pd.notna(x) and str(x).strip() not in {"", ".", "NA", "NaN", "None", "none", "NULL"}:
            return x
    return ""


def join_unique(series):
    vals = []
    seen = set()
    for x in series:
        if pd.isna(x):
            continue
        for y in str(x).split(","):
            y = y.strip()
            if y and y not in {".", "NA", "NaN", "None", "none", "NULL"} and y not in seen:
                vals.append(y)
                seen.add(y)
    return ",".join(vals)


def main():
    parser = argparse.ArgumentParser(description="Make final top1 actual-window and unique-gene tables.")
    parser.add_argument("--outdir", default="postprocess_results")
    parser.add_argument("--top1-windows", default=None)
    parser.add_argument("--top1-expanded", default=None)
    parser.add_argument("--windows-out", default=None)
    parser.add_argument("--genes-out", default=None)
    args = parser.parse_args()

    top1_windows_file = args.top1_windows or os.path.join(args.outdir, "conflict_quadrant_top1_windows.tsv")
    top1_expanded_file = args.top1_expanded or os.path.join(args.outdir, "conflict_quadrant_top1_window_gene_expanded.tsv")
    windows_out = args.windows_out or os.path.join(args.outdir, "final_top1_actual_windows.tsv")
    genes_out = args.genes_out or os.path.join(args.outdir, "final_top1_unique_genes.tsv")
    audit_out = os.path.join(args.outdir, "final_top1_tables_audit.tsv")

    windows = pd.read_csv(top1_windows_file, sep="\t", low_memory=False)
    required = {"CHROM", "start", "end", "quadrant"}
    if not required.issubset(windows.columns):
        raise ValueError(f"{top1_windows_file} must contain columns: {sorted(required)}")

    windows = windows[windows["quadrant"].isin(VALID_QUADRANTS)].copy()
    windows = numeric_window_coords(windows, "CHROM", "start", "end")
    windows = windows.dropna(subset=["CHROM", "start", "end"])
    windows = windows.drop_duplicates(["CHROM", "start", "end", "quadrant"]).copy()

    if "class" in windows.columns and "legacy_class" not in windows.columns:
        windows["legacy_class"] = windows["class"]
    windows["class"] = windows["quadrant"]

    # Put the most useful identity columns first, then keep all original columns.
    preferred = [
        "CHROM", "start", "end", "mid", "quadrant", "class", "legacy_class",
        "quadrant_top1_rank", "quadrant_top1_score", "quadrant_top1_metric",
        "quadrant_top1_cutoff", "C_score", "S_score", "M_score", "T_score",
        "selection_function_score", "genes_in_window", "gene_input_in_window",
    ]
    ordered_cols = [c for c in preferred if c in windows.columns] + [c for c in windows.columns if c not in preferred]
    windows = windows[ordered_cols]
    windows.to_csv(windows_out, sep="\t", index=False)

    expanded = pd.read_csv(top1_expanded_file, sep="\t", low_memory=False)
    required_expanded = {"window_CHROM", "window_start", "window_end", "top1_quadrant"}
    if not required_expanded.issubset(expanded.columns):
        raise ValueError(f"{top1_expanded_file} must contain columns: {sorted(required_expanded)}")

    expanded = expanded[expanded["top1_quadrant"].isin(VALID_QUADRANTS)].copy()
    expanded = numeric_window_coords(expanded, "window_CHROM", "window_start", "window_end")
    if "class" in expanded.columns and "legacy_class" not in expanded.columns:
        expanded["legacy_class"] = expanded["class"]
    expanded["class"] = expanded["top1_quadrant"]

    gene_key_candidates = ["gene_id", "gene_std", "gene_input", "gene_name", "display_name"]
    gene_key = next((c for c in gene_key_candidates if c in expanded.columns), None)
    if gene_key is None:
        raise ValueError("Could not find a gene identifier column in expanded table.")

    expanded[gene_key] = expanded[gene_key].astype(str).str.strip()
    expanded = expanded[
        expanded[gene_key].notna()
        & ~expanded[gene_key].isin(["", ".", "NA", "NaN", "None", "none", "NULL"])
    ].copy()

    expanded["window_id"] = (
        expanded["window_CHROM"].astype(str)
        + ":"
        + expanded["window_start"].astype("Int64").astype(str)
        + "-"
        + expanded["window_end"].astype("Int64").astype(str)
    )

    agg_spec = {
        "top1_quadrant": join_unique,
        "window_id": join_unique,
        "window_CHROM": join_unique,
        "window_start": "min",
        "window_end": "max",
    }
    optional_first = [
        "gene_std", "gene_id", "gene_id_dash", "gene_name", "display_name", "gene_label",
        "gene_biotype", "gene_CHROM", "gene_start", "gene_end", "candidate_class",
        "display_description",
    ]
    optional_max = [
        "C_score", "S_score", "M_score", "T_score", "selection_function_score",
        "musk_score_raw", "celltype_specificity", "peak2gene_max_score",
        "mature_bulk_score", "bulkDE_log2FC", "bulk_age_sig",
        "bulk_tissue_log2FC", "bulk_tissue_sig", "tissue_bonus",
        "bulk_interaction_log2FC", "bulk_interaction_sig", "interaction_bonus",
    ]
    for c in optional_first:
        if c in expanded.columns and c != gene_key:
            agg_spec[c] = first_nonempty
    for c in optional_max:
        if c in expanded.columns:
            expanded[c] = pd.to_numeric(expanded[c], errors="coerce")
            agg_spec[c] = "max"

    genes = expanded.groupby(gene_key, dropna=False).agg(agg_spec).reset_index()
    genes = genes.rename(columns={
        gene_key: "gene_key",
        "top1_quadrant": "top1_quadrants",
        "window_id": "top1_window_ids",
        "window_CHROM": "top1_window_chromosomes",
        "window_start": "min_top1_window_start",
        "window_end": "max_top1_window_end",
    })

    # Counts per unique gene.
    window_count = expanded.groupby(gene_key)["window_id"].nunique().rename("n_top1_windows").reset_index()
    quadrant_count = expanded.groupby(gene_key)["top1_quadrant"].nunique().rename("n_top1_quadrants").reset_index()
    window_count = window_count.rename(columns={gene_key: "gene_key"})
    quadrant_count = quadrant_count.rename(columns={gene_key: "gene_key"})
    genes = genes.merge(window_count, on="gene_key", how="left").merge(quadrant_count, on="gene_key", how="left")

    preferred_gene_cols = [
        "gene_key", "gene_std", "gene_id", "gene_id_dash", "gene_name", "display_name",
        "gene_label", "gene_biotype", "gene_CHROM", "gene_start", "gene_end",
        "top1_quadrants", "n_top1_quadrants", "n_top1_windows", "top1_window_ids",
        "top1_window_chromosomes", "min_top1_window_start", "max_top1_window_end",
        "C_score", "S_score", "M_score", "T_score", "selection_function_score",
        "musk_score_raw", "candidate_class", "display_description",
    ]
    gene_cols = [c for c in preferred_gene_cols if c in genes.columns] + [c for c in genes.columns if c not in preferred_gene_cols]
    genes = genes[gene_cols]
    sort_cols = [c for c in ["n_top1_quadrants", "n_top1_windows", "T_score", "C_score"] if c in genes.columns]
    if sort_cols:
        genes = genes.sort_values(sort_cols, ascending=[False] * len(sort_cols), na_position="last")
    genes.to_csv(genes_out, sep="\t", index=False)

    audit_rows = []
    for q in VALID_QUADRANTS:
        qwin = windows[windows["quadrant"] == q]
        qexp = expanded[expanded["top1_quadrant"] == q]
        audit_rows.append({
            "quadrant": q,
            "actual_top1_windows": qwin.drop_duplicates(["CHROM", "start", "end"]).shape[0],
            "expanded_distinct_windows": qexp.drop_duplicates(["window_CHROM", "window_start", "window_end"]).shape[0],
            "expanded_gene_rows": qexp.shape[0],
            "unique_genes_within_quadrant": qexp[gene_key].nunique(),
        })
    audit = pd.DataFrame(audit_rows)
    audit.loc[len(audit)] = {
        "quadrant": "TOTAL",
        "actual_top1_windows": windows.drop_duplicates(["CHROM", "start", "end"]).shape[0],
        "expanded_distinct_windows": expanded.drop_duplicates(["window_CHROM", "window_start", "window_end"]).shape[0],
        "expanded_gene_rows": expanded.shape[0],
        "unique_genes_within_quadrant": expanded[gene_key].nunique(),
    }
    audit.to_csv(audit_out, sep="\t", index=False)

    print("[OK] final top1 tables written")
    print(f"  windows: {windows_out}")
    print(f"  genes  : {genes_out}")
    print(f"  audit  : {audit_out}")
    print(audit.to_string(index=False))


if __name__ == "__main__":
    main()
