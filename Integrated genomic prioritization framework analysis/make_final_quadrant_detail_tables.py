#!/usr/bin/env python3

"""
Create final all-quadrant detail tables:

1. final_quadrant_actual_windows.tsv
   One row per valid actual window across the four quadrants.

2. final_quadrant_unique_genes.tsv
   One row per globally unique gene appearing in the four-quadrant windows.
   If a gene appears in multiple quadrants/windows, those memberships are
   summarised in the same row.
"""

import argparse
import os
import pandas as pd


QUADRANT_ORDER = [
    "conservation_priority",
    "breeding_function_priority",
    "conflict_quadrant",
    "background",
]

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


def first_nonempty(series):
    for x in series:
        if valid_gene(x):
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
            if valid_gene(y) and y not in seen:
                vals.append(y)
                seen.add(y)
    return ",".join(vals)


def numeric_window_coords(df, chrom_col, start_col, end_col):
    df = df.copy()
    df[chrom_col] = df[chrom_col].astype(str)
    df[start_col] = pd.to_numeric(df[start_col], errors="coerce")
    df[end_col] = pd.to_numeric(df[end_col], errors="coerce")
    return df


def main():
    parser = argparse.ArgumentParser(description="Make final all-quadrant actual-window and unique-gene detail tables.")
    parser.add_argument("--outdir", default="postprocess_results")
    parser.add_argument("--quadrant-table", default=None)
    parser.add_argument("--expanded-all", default=None)
    parser.add_argument("--windows-out", default=None)
    parser.add_argument("--genes-out", default=None)
    args = parser.parse_args()

    quadrant_file = args.quadrant_table or os.path.join(args.outdir, "conflict_quadrant_table.tsv")
    expanded_file = args.expanded_all or os.path.join(args.outdir, "conflict_window_gene_expanded.tsv")
    windows_out = args.windows_out or os.path.join(args.outdir, "final_quadrant_actual_windows.tsv")
    genes_out = args.genes_out or os.path.join(args.outdir, "final_quadrant_unique_genes.tsv")
    audit_out = os.path.join(args.outdir, "final_quadrant_detail_tables_audit.tsv")

    windows = pd.read_csv(quadrant_file, sep="\t", low_memory=False)
    required = {"CHROM", "start", "end", "quadrant"}
    if not required.issubset(windows.columns):
        raise ValueError(f"{quadrant_file} must contain columns: {sorted(required)}")

    windows["quadrant"] = windows["quadrant"].map(normalize_quadrant)
    windows = numeric_window_coords(windows, "CHROM", "start", "end")
    windows = windows.dropna(subset=["CHROM", "start", "end"])
    windows = windows[windows["quadrant"].isin(QUADRANT_ORDER)].copy()
    windows = windows.drop_duplicates(["CHROM", "start", "end", "quadrant"]).copy()

    if "class" in windows.columns and "legacy_class" not in windows.columns:
        windows["legacy_class"] = windows["class"]
    windows["class"] = windows["quadrant"]
    windows["quadrant_order"] = windows["quadrant"].map({q: i + 1 for i, q in enumerate(QUADRANT_ORDER)})

    preferred_window_cols = [
        "CHROM", "start", "end", "mid", "quadrant", "class", "legacy_class", "quadrant_order",
        "C_score", "S_score", "M_score", "T_score", "selection_function_score",
        "genes_in_window", "gene_input_in_window",
        "quadrant_x_cutoff", "quadrant_y_cutoff",
        "quadrant_top1_flag", "quadrant_top1_rank", "quadrant_top1_score",
        "quadrant_top1_metric", "quadrant_top1_cutoff",
    ]
    window_cols = [c for c in preferred_window_cols if c in windows.columns] + [c for c in windows.columns if c not in preferred_window_cols]
    windows = windows[window_cols].sort_values(["quadrant_order", "CHROM", "start", "end"], na_position="last")
    windows.to_csv(windows_out, sep="\t", index=False)

    # Build unique-gene detail table from all-window expanded table plus final quadrant assignment.
    key_map = windows[["CHROM", "start", "end", "quadrant"]].rename(columns={
        "CHROM": "window_CHROM",
        "start": "window_start",
        "end": "window_end",
    })
    key_map = numeric_window_coords(key_map, "window_CHROM", "window_start", "window_end")

    expanded = pd.read_csv(expanded_file, sep="\t", low_memory=False)
    required_expanded = {"window_CHROM", "window_start", "window_end"}
    if not required_expanded.issubset(expanded.columns):
        raise ValueError(f"{expanded_file} must contain columns: {sorted(required_expanded)}")

    expanded = numeric_window_coords(expanded, "window_CHROM", "window_start", "window_end")
    expanded = expanded.merge(key_map, on=["window_CHROM", "window_start", "window_end"], how="inner")
    if "class" in expanded.columns and "legacy_class" not in expanded.columns:
        expanded["legacy_class"] = expanded["class"]
    expanded["class"] = expanded["quadrant"]

    gene_key_candidates = ["gene_id", "gene_std", "gene_input", "gene_name", "display_name"]
    gene_key = next((c for c in gene_key_candidates if c in expanded.columns), None)
    if gene_key is None:
        raise ValueError("Could not find a gene identifier column in expanded table.")

    expanded[gene_key] = expanded[gene_key].astype(str).str.strip()
    expanded = expanded[expanded[gene_key].map(valid_gene)].copy()
    expanded["window_id"] = (
        expanded["window_CHROM"].astype(str)
        + ":"
        + expanded["window_start"].astype("Int64").astype(str)
        + "-"
        + expanded["window_end"].astype("Int64").astype(str)
    )

    agg_spec = {
        "quadrant": join_unique,
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
        "quadrant": "quadrants",
        "window_id": "window_ids",
        "window_CHROM": "window_chromosomes",
        "window_start": "min_window_start",
        "window_end": "max_window_end",
    })

    window_count = expanded.groupby(gene_key)["window_id"].nunique().rename("n_windows").reset_index()
    quadrant_count = expanded.groupby(gene_key)["quadrant"].nunique().rename("n_quadrants").reset_index()
    window_count = window_count.rename(columns={gene_key: "gene_key"})
    quadrant_count = quadrant_count.rename(columns={gene_key: "gene_key"})
    genes = genes.merge(window_count, on="gene_key", how="left").merge(quadrant_count, on="gene_key", how="left")

    preferred_gene_cols = [
        "gene_key", "gene_std", "gene_id", "gene_id_dash", "gene_name", "display_name",
        "gene_label", "gene_biotype", "gene_CHROM", "gene_start", "gene_end",
        "quadrants", "n_quadrants", "n_windows", "window_ids", "window_chromosomes",
        "min_window_start", "max_window_end",
        "C_score", "S_score", "M_score", "T_score", "selection_function_score",
        "musk_score_raw", "candidate_class", "display_description",
    ]
    gene_cols = [c for c in preferred_gene_cols if c in genes.columns] + [c for c in genes.columns if c not in preferred_gene_cols]
    genes = genes[gene_cols]
    sort_cols = [c for c in ["n_quadrants", "n_windows", "T_score", "C_score"] if c in genes.columns]
    if sort_cols:
        genes = genes.sort_values(sort_cols, ascending=[False] * len(sort_cols), na_position="last")
    genes.to_csv(genes_out, sep="\t", index=False)

    audit_rows = []
    for q in QUADRANT_ORDER:
        qwin = windows[windows["quadrant"] == q]
        qexp = expanded[expanded["quadrant"] == q]
        audit_rows.append({
            "quadrant": q,
            "actual_windows": qwin.drop_duplicates(["CHROM", "start", "end"]).shape[0],
            "expanded_distinct_windows": qexp.drop_duplicates(["window_CHROM", "window_start", "window_end"]).shape[0],
            "expanded_gene_rows": qexp.shape[0],
            "unique_genes_within_quadrant": qexp[gene_key].nunique(),
        })
    audit = pd.DataFrame(audit_rows)
    audit.loc[len(audit)] = {
        "quadrant": "TOTAL",
        "actual_windows": windows.drop_duplicates(["CHROM", "start", "end"]).shape[0],
        "expanded_distinct_windows": expanded.drop_duplicates(["window_CHROM", "window_start", "window_end"]).shape[0],
        "expanded_gene_rows": expanded.shape[0],
        "unique_genes_within_quadrant": expanded[gene_key].nunique(),
    }
    audit.to_csv(audit_out, sep="\t", index=False)

    print("[OK] final all-quadrant detail tables written")
    print(f"  windows: {windows_out}")
    print(f"  genes  : {genes_out}")
    print(f"  audit  : {audit_out}")
    print(audit.to_string(index=False))


if __name__ == "__main__":
    main()
