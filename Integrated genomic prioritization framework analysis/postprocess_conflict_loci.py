#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import numpy as np
import pandas as pd

# optional plotting
try:
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except Exception:
    HAVE_MPL = False


# =========================================================
# Utility
# =========================================================

def zscore(series):
    s = pd.to_numeric(series, errors="coerce")
    mu = s.mean(skipna=True)
    sd = s.std(skipna=True)
    if pd.isna(sd) or sd == 0:
        return pd.Series([0] * len(s), index=s.index)
    return (s - mu) / sd

def safe_numeric(df, cols):
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def neglog10_p(series, min_p=1e-300):
    s = pd.to_numeric(series, errors="coerce")
    s = s.clip(lower=min_p)
    return -np.log10(s)

def first_nonempty(*vals):
    for v in vals:
        if pd.notna(v) and str(v).strip() != "":
            return v
    return np.nan


# =========================================================
# Annotation helpers
# =========================================================

def build_annotation_lookup(anno):
    anno = anno.copy()

    keep_cols = [c for c in [
        "gene_id", "gene_id_dash", "gene_name", "display_name",
        "display_description", "gene_biotype"
    ] if c in anno.columns]
    anno = anno[keep_cols].copy()

    lookup_tables = {}

    if "gene_id_dash" in anno.columns:
        tmp = anno.dropna(subset=["gene_id_dash"]).drop_duplicates(subset=["gene_id_dash"])
        lookup_tables["gene_id_dash"] = tmp.set_index("gene_id_dash")

    if "gene_id" in anno.columns:
        tmp = anno.dropna(subset=["gene_id"]).drop_duplicates(subset=["gene_id"])
        lookup_tables["gene_id"] = tmp.set_index("gene_id")

    if "gene_name" in anno.columns:
        tmp = anno.dropna(subset=["gene_name"]).drop_duplicates(subset=["gene_name"])
        lookup_tables["gene_name"] = tmp.set_index("gene_name")

    if "display_name" in anno.columns:
        tmp = anno.dropna(subset=["display_name"]).drop_duplicates(subset=["display_name"])
        lookup_tables["display_name"] = tmp.set_index("display_name")

    return lookup_tables

def resolve_gene_identifier(gene_value, lookup_tables):
    if pd.isna(gene_value):
        return None

    g = str(gene_value).strip()
    if g == "":
        return None

    if "gene_id_dash" in lookup_tables and g in lookup_tables["gene_id_dash"].index:
        row = lookup_tables["gene_id_dash"].loc[g]
        return {
            "match_strategy": "gene_id_dash",
            "gene_std": row.get("gene_id", g),
            "gene_id": row.get("gene_id", np.nan),
            "gene_id_dash": row.get("gene_id_dash", np.nan),
            "gene_name": row.get("gene_name", np.nan),
            "display_name": row.get("display_name", np.nan),
            "display_description": row.get("display_description", np.nan),
            "gene_biotype": row.get("gene_biotype", np.nan)
        }

    if "gene_id" in lookup_tables and g in lookup_tables["gene_id"].index:
        row = lookup_tables["gene_id"].loc[g]
        return {
            "match_strategy": "gene_id",
            "gene_std": row.get("gene_id", g),
            "gene_id": row.get("gene_id", np.nan),
            "gene_id_dash": row.get("gene_id_dash", np.nan),
            "gene_name": row.get("gene_name", np.nan),
            "display_name": row.get("display_name", np.nan),
            "display_description": row.get("display_description", np.nan),
            "gene_biotype": row.get("gene_biotype", np.nan)
        }

    if "gene_name" in lookup_tables and g in lookup_tables["gene_name"].index:
        row = lookup_tables["gene_name"].loc[g]
        return {
            "match_strategy": "gene_name",
            "gene_std": row.get("gene_id", g),
            "gene_id": row.get("gene_id", np.nan),
            "gene_id_dash": row.get("gene_id_dash", np.nan),
            "gene_name": row.get("gene_name", np.nan),
            "display_name": row.get("display_name", np.nan),
            "display_description": row.get("display_description", np.nan),
            "gene_biotype": row.get("gene_biotype", np.nan)
        }

    if "display_name" in lookup_tables and g in lookup_tables["display_name"].index:
        row = lookup_tables["display_name"].loc[g]
        return {
            "match_strategy": "display_name",
            "gene_std": row.get("gene_id", g),
            "gene_id": row.get("gene_id", np.nan),
            "gene_id_dash": row.get("gene_id_dash", np.nan),
            "gene_name": row.get("gene_name", np.nan),
            "display_name": row.get("display_name", np.nan),
            "display_description": row.get("display_description", np.nan),
            "gene_biotype": row.get("gene_biotype", np.nan)
        }

    return {
        "match_strategy": "unresolved",
        "gene_std": g,
        "gene_id": np.nan,
        "gene_id_dash": np.nan,
        "gene_name": np.nan,
        "display_name": np.nan,
        "display_description": np.nan,
        "gene_biotype": np.nan
    }


# =========================================================
# Load gene table
# =========================================================

def load_and_unify_musk_genes(musk_file, anno_file, coord_file):
    musk = pd.read_csv(musk_file, sep="\t", low_memory=False)
    anno = pd.read_csv(anno_file, sep="\t", low_memory=False)
    coords = pd.read_csv(coord_file, sep="\t", low_memory=False)

    lookup_tables = build_annotation_lookup(anno)

    musk = musk.rename(columns={"gene": "gene_input"}).copy()

    resolved = musk["gene_input"].apply(lambda x: resolve_gene_identifier(x, lookup_tables))
    resolved_df = pd.DataFrame(list(resolved))

    mg = pd.concat([musk.reset_index(drop=True), resolved_df.reset_index(drop=True)], axis=1)

    coords = coords.rename(columns={"gene": "gene_std"}).copy()
    coords["start"] = pd.to_numeric(coords["start"], errors="coerce")
    coords["end"] = pd.to_numeric(coords["end"], errors="coerce")

    mg = pd.merge(mg, coords, on="gene_std", how="left")

    numeric_cols = [
        "scRNA_marker","celltype_specificity","ATAC_peak_support","motif_support","musk_score_raw",
        "bulk_support_score","peak2gene_support","peak2gene_n_links","peak2gene_max_score",
        "mature_bulk_score","tissue_bonus","interaction_bonus","adult_up",
        "bulkDE_log2FC","bulkDE_padj",
        "bulk_tissue_log2FC","bulk_tissue_padj",
        "bulk_interaction_log2FC","bulk_interaction_padj"
    ]
    mg = safe_numeric(mg, [c for c in numeric_cols if c in mg.columns])

    mg["bulk_age_sig"] = neglog10_p(mg["bulkDE_padj"]) if "bulkDE_padj" in mg.columns else np.nan
    mg["bulk_tissue_sig"] = neglog10_p(mg["bulk_tissue_padj"]) if "bulk_tissue_padj" in mg.columns else np.nan
    mg["bulk_interaction_sig"] = neglog10_p(mg["bulk_interaction_padj"]) if "bulk_interaction_padj" in mg.columns else np.nan

    return mg


# =========================================================
# Expand window -> gene rows
# =========================================================

def build_window_gene_expanded(conflict_df, mg):
    records = []

    mg2 = mg.dropna(subset=["CHROM", "start", "end"]).copy()

    for _, w in conflict_df.iterrows():
        chrom = w["CHROM"]
        wstart = pd.to_numeric(w["start"], errors="coerce")
        wend = pd.to_numeric(w["end"], errors="coerce")

        if pd.isna(chrom) or pd.isna(wstart) or pd.isna(wend):
            continue

        hits = mg2[
            (mg2["CHROM"] == chrom) &
            (pd.to_numeric(mg2["end"], errors="coerce") >= wstart) &
            (pd.to_numeric(mg2["start"], errors="coerce") <= wend)
        ].copy()

        if hits.empty:
            continue

        for _, g in hits.iterrows():
            gene_label = first_nonempty(
                g.get("display_name", np.nan),
                g.get("gene_name", np.nan),
                g.get("gene_std", np.nan),
                g.get("gene_input", np.nan)
            )

            records.append({
                "window_CHROM": w["CHROM"],
                "window_start": w["start"],
                "window_end": w["end"],
                "class": w.get("class", np.nan),
                "high_conflict_flag": w.get("high_conflict_flag", 0),
                "C_score": w.get("C_score", np.nan),
                "S_score": w.get("S_score", np.nan),
                "M_score": w.get("M_score", np.nan),
                "T_score": w.get("T_score", np.nan),
                "selection_function_score": 0.4 * pd.to_numeric(w.get("S_score", np.nan), errors="coerce") +
                                            0.6 * pd.to_numeric(w.get("M_score", np.nan), errors="coerce"),

                "gene_input": g.get("gene_input", np.nan),
                "gene_std": g.get("gene_std", np.nan),
                "gene_id": g.get("gene_id", np.nan),
                "gene_id_dash": g.get("gene_id_dash", np.nan),
                "gene_name": g.get("gene_name", np.nan),
                "display_name": g.get("display_name", np.nan),
                "gene_label": gene_label,
                "match_strategy": g.get("match_strategy", np.nan),
                "display_description": g.get("display_description", np.nan),
                "gene_biotype": g.get("gene_biotype", np.nan),

                "gene_CHROM": g.get("CHROM", np.nan),
                "gene_start": g.get("start", np.nan),
                "gene_end": g.get("end", np.nan),

                "musk_score_raw": g.get("musk_score_raw", np.nan),
                "candidate_class": g.get("candidate_class", np.nan),

                "scRNA_marker": g.get("scRNA_marker", np.nan),
                "adult_up": g.get("adult_up", np.nan),
                "celltype_specificity": g.get("celltype_specificity", np.nan),
                "ATAC_peak_support": g.get("ATAC_peak_support", np.nan),
                "peak2gene_support": g.get("peak2gene_support", np.nan),
                "motif_support": g.get("motif_support", np.nan),
                "peak2gene_max_score": g.get("peak2gene_max_score", np.nan),

                "mature_bulk_score": g.get("mature_bulk_score", np.nan),
                "bulkDE_log2FC": g.get("bulkDE_log2FC", np.nan),
                "bulkDE_padj": g.get("bulkDE_padj", np.nan),
                "bulk_age_sig": g.get("bulk_age_sig", np.nan),

                "bulk_tissue_log2FC": g.get("bulk_tissue_log2FC", np.nan),
                "bulk_tissue_padj": g.get("bulk_tissue_padj", np.nan),
                "bulk_tissue_sig": g.get("bulk_tissue_sig", np.nan),
                "tissue_bonus": g.get("tissue_bonus", np.nan),

                "bulk_interaction_log2FC": g.get("bulk_interaction_log2FC", np.nan),
                "bulk_interaction_padj": g.get("bulk_interaction_padj", np.nan),
                "bulk_interaction_sig": g.get("bulk_interaction_sig", np.nan),
                "interaction_bonus": g.get("interaction_bonus", np.nan)
            })

    expanded = pd.DataFrame(records)

    if expanded.empty:
        return expanded

    numeric_cols = [
        "high_conflict_flag", "C_score", "S_score", "M_score", "T_score", "selection_function_score",
        "musk_score_raw", "scRNA_marker", "adult_up", "celltype_specificity", "ATAC_peak_support",
        "peak2gene_support", "motif_support", "peak2gene_max_score", "mature_bulk_score",
        "bulkDE_log2FC", "bulkDE_padj", "bulk_age_sig",
        "bulk_tissue_log2FC", "bulk_tissue_padj", "bulk_tissue_sig", "tissue_bonus",
        "bulk_interaction_log2FC", "bulk_interaction_padj", "bulk_interaction_sig", "interaction_bonus"
    ]
    expanded = safe_numeric(expanded, [c for c in numeric_cols if c in expanded.columns])

    return expanded


# =========================================================
# Gene-level ranking
# =========================================================

def build_top_conflict_genes(expanded):
    if expanded.empty:
        return expanded

    tmp = expanded.copy()

    tmp["z_window_T"] = zscore(tmp["T_score"])
    tmp["z_musk_raw"] = zscore(tmp["musk_score_raw"])
    tmp["z_celltype_spec"] = zscore(tmp["celltype_specificity"])
    tmp["z_peak2gene_score"] = zscore(tmp["peak2gene_max_score"])

    tmp["gene_conflict_priority_row"] = (
        0.60 * tmp["z_window_T"].fillna(0) +
        0.25 * tmp["z_musk_raw"].fillna(0) +
        0.05 * tmp["z_celltype_spec"].fillna(0) +
        0.05 * tmp["z_peak2gene_score"].fillna(0) +
        0.03 * tmp["adult_up"].fillna(0) +
        0.01 * tmp["scRNA_marker"].fillna(0) +
        0.01 * tmp["peak2gene_support"].fillna(0)
    )

    grouped = tmp.groupby("gene_std", dropna=False)

    rows = []
    for gene_std, gdf in grouped:
        best_idx = gdf["gene_conflict_priority_row"].idxmax()
        best_row = gdf.loc[best_idx]

        rows.append({
            "gene_std": gene_std,
            "gene_input": best_row.get("gene_input", np.nan),
            "gene_id": best_row.get("gene_id", np.nan),
            "gene_id_dash": best_row.get("gene_id_dash", np.nan),
            "gene_name": best_row.get("gene_name", np.nan),
            "display_name": best_row.get("display_name", np.nan),
            "gene_label": best_row.get("gene_label", np.nan),
            "match_strategy": best_row.get("match_strategy", np.nan),
            "display_description": best_row.get("display_description", np.nan),
            "gene_biotype": best_row.get("gene_biotype", np.nan),

            "gene_CHROM": best_row.get("gene_CHROM", np.nan),
            "gene_start": best_row.get("gene_start", np.nan),
            "gene_end": best_row.get("gene_end", np.nan),

            "n_conflict_windows": int((gdf["class"] == "conflict").sum()),
            "n_all_windows": int(gdf.shape[0]),
            "any_high_conflict": int((gdf["high_conflict_flag"].fillna(0) > 0).any()),

            "best_window_CHROM": best_row.get("window_CHROM", np.nan),
            "best_window_start": best_row.get("window_start", np.nan),
            "best_window_end": best_row.get("window_end", np.nan),
            "best_window_class": best_row.get("class", np.nan),

            "best_C_score": gdf["C_score"].max(skipna=True),
            "best_S_score": gdf["S_score"].max(skipna=True),
            "best_M_score": gdf["M_score"].max(skipna=True),
            "best_T_score": gdf["T_score"].max(skipna=True),
            "best_selection_function_score": gdf["selection_function_score"].max(skipna=True),

            "musk_score_raw": gdf["musk_score_raw"].max(skipna=True),
            "mean_musk_score_raw": gdf["musk_score_raw"].mean(skipna=True),

            "scRNA_marker": gdf["scRNA_marker"].max(skipna=True),
            "adult_up": gdf["adult_up"].max(skipna=True),
            "celltype_specificity": gdf["celltype_specificity"].max(skipna=True),
            "ATAC_peak_support": gdf["ATAC_peak_support"].max(skipna=True),
            "peak2gene_support": gdf["peak2gene_support"].max(skipna=True),
            "motif_support": gdf["motif_support"].max(skipna=True),
            "peak2gene_max_score": gdf["peak2gene_max_score"].max(skipna=True),

            "mature_bulk_score": gdf["mature_bulk_score"].max(skipna=True),
            "bulkDE_log2FC": gdf["bulkDE_log2FC"].max(skipna=True),
            "bulkDE_padj": gdf["bulkDE_padj"].min(skipna=True),
            "bulk_age_sig": gdf["bulk_age_sig"].max(skipna=True),

            "bulk_tissue_log2FC": gdf["bulk_tissue_log2FC"].max(skipna=True),
            "bulk_tissue_padj": gdf["bulk_tissue_padj"].min(skipna=True),
            "bulk_tissue_sig": gdf["bulk_tissue_sig"].max(skipna=True),
            "tissue_bonus": gdf["tissue_bonus"].max(skipna=True),

            "bulk_interaction_log2FC": gdf["bulk_interaction_log2FC"].max(skipna=True),
            "bulk_interaction_padj": gdf["bulk_interaction_padj"].min(skipna=True),
            "bulk_interaction_sig": gdf["bulk_interaction_sig"].max(skipna=True),
            "interaction_bonus": gdf["interaction_bonus"].max(skipna=True),

            "candidate_class_best": best_row.get("candidate_class", np.nan),
            "candidate_classes_all": ",".join(sorted(set(gdf["candidate_class"].dropna().astype(str).tolist()))),

            "genes_priority_row_max": gdf["gene_conflict_priority_row"].max(skipna=True)
        })

    gene_df = pd.DataFrame(rows)

    gene_df["z_best_T_score"] = zscore(gene_df["best_T_score"])
    gene_df["z_musk_score_raw"] = zscore(gene_df["musk_score_raw"])
    gene_df["z_celltype_specificity"] = zscore(gene_df["celltype_specificity"])
    gene_df["z_peak2gene_max_score"] = zscore(gene_df["peak2gene_max_score"])
    gene_df["z_n_conflict_windows"] = zscore(gene_df["n_conflict_windows"])

    gene_df["gene_conflict_priority_score"] = (
        0.50 * gene_df["z_best_T_score"].fillna(0) +
        0.20 * gene_df["z_musk_score_raw"].fillna(0) +
        0.10 * gene_df["z_celltype_specificity"].fillna(0) +
        0.05 * gene_df["z_peak2gene_max_score"].fillna(0) +
        0.05 * gene_df["z_n_conflict_windows"].fillna(0) +
        0.05 * gene_df["adult_up"].fillna(0) +
        0.03 * gene_df["scRNA_marker"].fillna(0) +
        0.02 * gene_df["peak2gene_support"].fillna(0)
    )

    gene_df = gene_df.sort_values(
        by=["any_high_conflict", "gene_conflict_priority_score", "best_T_score", "musk_score_raw"],
        ascending=[False, False, False, False],
        na_position="last"
    ).reset_index(drop=True)

    gene_df["gene_rank"] = np.arange(1, gene_df.shape[0] + 1)

    return gene_df


# =========================================================
# Quadrant table / plot
# =========================================================

def build_quadrant_table(conflict_df):
    df = conflict_df.copy()

    df["C_score"] = pd.to_numeric(df["C_score"], errors="coerce")
    df["S_score"] = pd.to_numeric(df["S_score"], errors="coerce")
    df["M_score"] = pd.to_numeric(df["M_score"], errors="coerce")
    df["T_score"] = pd.to_numeric(df["T_score"], errors="coerce")

    df["selection_function_score"] = (
        0.4 * df["S_score"] +
        0.6 * df["M_score"]
    )

    qC = df["C_score"].quantile(0.90)
    qY = df["selection_function_score"].quantile(0.90)

    df["quadrant"] = "background_low_low"
    df.loc[(df["C_score"] >= qC) & (df["selection_function_score"] < qY), "quadrant"] = "conservation_priority"
    df.loc[(df["C_score"] < qC) & (df["selection_function_score"] >= qY), "quadrant"] = "breeding_function_priority"
    df.loc[(df["C_score"] >= qC) & (df["selection_function_score"] >= qY), "quadrant"] = "conflict_quadrant"

    df["quadrant_x_cutoff"] = qC
    df["quadrant_y_cutoff"] = qY

    # Within the three priority quadrants, keep only the top 1% as the main
    # discussion set.  Use a quadrant-aware ranking metric:
    #   - conservation priority: conservation risk (C_score)
    #   - breeding/function priority: selection/function score
    #   - conflict quadrant: final interaction score (T_score)
    # This avoids forcing breeding-priority windows to sit near the C cutoff
    # merely because T_score multiplies by C_score.
    df["quadrant_top1_flag"] = 0
    df["quadrant_top1_rank"] = np.nan
    df["quadrant_top1_cutoff"] = np.nan
    df["quadrant_top1_metric"] = ""
    df["quadrant_top1_score"] = np.nan
    df.loc[df["quadrant"] == "conservation_priority", "quadrant_top1_score"] = df["C_score"]
    df.loc[df["quadrant"] == "conservation_priority", "quadrant_top1_metric"] = "C_score"
    df.loc[df["quadrant"] == "breeding_function_priority", "quadrant_top1_score"] = df["selection_function_score"]
    df.loc[df["quadrant"] == "breeding_function_priority", "quadrant_top1_metric"] = "selection_function_score"
    df.loc[df["quadrant"] == "conflict_quadrant", "quadrant_top1_score"] = df["T_score"]
    df.loc[df["quadrant"] == "conflict_quadrant", "quadrant_top1_metric"] = "T_score"

    for quad, idx in df.groupby("quadrant").groups.items():
        if quad == "background_low_low":
            continue
        sub = df.loc[idx].copy()
        valid = sub.dropna(subset=["quadrant_top1_score"]).copy()
        if valid.empty:
            continue

        n_top = max(1, int(np.ceil(valid.shape[0] * 0.01)))
        valid = valid.sort_values("quadrant_top1_score", ascending=False)
        top_idx = valid.head(n_top).index
        cutoff = valid.iloc[n_top - 1]["quadrant_top1_score"]

        df.loc[top_idx, "quadrant_top1_flag"] = 1
        df.loc[valid.index, "quadrant_top1_rank"] = np.arange(1, valid.shape[0] + 1)
        df.loc[idx, "quadrant_top1_cutoff"] = cutoff

    return df


def plot_quadrant(df, gene_df, out_pdf, top_n_labels=5):
    if not HAVE_MPL:
        print("[WARN] matplotlib not available, skip plot.")
        return

    try:
        from adjustText import adjust_text
        HAVE_ADJUSTTEXT = True
    except Exception:
        HAVE_ADJUSTTEXT = False

    plt.rcParams.update({
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "legend.fontsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.8
    })

    plot_df = df.copy()
    plot_df["C_score"] = pd.to_numeric(plot_df["C_score"], errors="coerce")
    plot_df["selection_function_score"] = pd.to_numeric(plot_df["selection_function_score"], errors="coerce")
    plot_df = plot_df.dropna(subset=["C_score", "selection_function_score"]).copy()

    if plot_df.empty:
        print("[WARN] No valid points for plotting.")
        return

    qC = plot_df["quadrant_x_cutoff"].iloc[0]
    qY = plot_df["quadrant_y_cutoff"].iloc[0]

    color_map = {
        "background_low_low": "#6F6F6F",
        "conservation_priority": "#2C7FB8",
        "breeding_function_priority": "#4DAF6B",
        "conflict_quadrant": "#D84B3A"
    }
    label_map = {
        "conservation_priority": "Conservation priority top 1%",
        "breeding_function_priority": "Breeding/function priority top 1%",
        "conflict_quadrant": "Conflict quadrant top 1%"
    }

    x_global_max = plot_df["C_score"].max()
    x_min = min(plot_df["C_score"].min(), -0.5) - 0.3

    top1_priority = plot_df[
        (plot_df["quadrant"] != "background_low_low") &
        (plot_df["quadrant_top1_flag"].fillna(0) > 0)
    ].copy()
    top1_x_max = top1_priority["C_score"].max() if not top1_priority.empty else np.nan
    x_main_max = max(plot_df["C_score"].quantile(0.98), qC + 2.5)
    # No inset is drawn now, so the main panel must include all highlighted
    # priority top1% windows.  Otherwise conflict top1% points can be clipped.
    if pd.notna(top1_x_max):
        x_main_max = max(x_main_max, top1_x_max)
    x_main_max = x_main_max + max((x_main_max - x_min) * 0.03, 0.2)

    y_min = plot_df["selection_function_score"].min()
    y_max = plot_df["selection_function_score"].max()
    y_pad = max((y_max - y_min) * 0.05, 0.2)
    y_main_min = y_min - y_pad
    y_main_max = y_max + y_pad

    fig, ax = plt.subplots(figsize=(8.2, 6.6))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    order = [
        "conservation_priority",
        "breeding_function_priority",
        "conflict_quadrant"
    ]

    # Draw all non-top1 windows first as neutral background context.
    non_top = plot_df[(plot_df["quadrant_top1_flag"].fillna(0) <= 0) &
                      (plot_df["C_score"] <= x_main_max)].copy()
    if not non_top.empty:
        ax.scatter(
            non_top["C_score"], non_top["selection_function_score"],
            s=13, alpha=0.18, color="#C9C9C9",
            edgecolors="none", label="Other windows", zorder=1,
            rasterized=True
        )

    # Then draw the top 1% from each quadrant using quadrant-specific colors.
    for quad in order:
        sub_main = plot_df[
            (plot_df["quadrant"] == quad) &
            (plot_df["quadrant_top1_flag"].fillna(0) > 0) &
            (plot_df["C_score"] <= x_main_max)
        ].copy()
        if sub_main.empty:
            continue

        edge = "white" if quad == "conflict_quadrant" else "none"
        lw = 0.25 if quad == "conflict_quadrant" else 0
        ax.scatter(
            sub_main["C_score"], sub_main["selection_function_score"],
            s=26 if quad == "conflict_quadrant" else 22,
            alpha=0.86,
            color=color_map[quad],
            edgecolors=edge,
            linewidths=lw,
            label=label_map[quad],
            zorder=3,
            rasterized=True
        )

    # Quadrant thresholds: grey dashed lines.  The upper-right inset panel is
    # intentionally removed, so the figure contains only this main quadrant.
    ax.axvline(qC, color="#8C8C8C", linestyle="--", linewidth=1.0, zorder=0, label="90% quadrant thresholds")
    ax.axhline(qY, color="#8C8C8C", linestyle="--", linewidth=1.0, zorder=0)

    ax.set_xlim(x_min, x_main_max)
    ax.set_ylim(y_main_min, y_main_max)

    # Top 1% thresholds within each priority quadrant.  The threshold shape
    # follows the metric used for that quadrant:
    #   conservation: vertical C-score cutoff
    #   breeding/function: horizontal selection/function cutoff
    #   conflict: iso-T-score curve, y = T_cutoff / x
    top1_label_used = False
    x_grid = np.linspace(x_min, x_main_max, 1200)
    x_grid[np.abs(x_grid) < 1e-8] = np.nan
    for quad in order:
        qsub = plot_df[plot_df["quadrant"] == quad].copy()
        if qsub.empty or "quadrant_top1_cutoff" not in qsub.columns:
            continue
        t_cut = pd.to_numeric(qsub["quadrant_top1_cutoff"], errors="coerce").dropna()
        if t_cut.empty:
            continue
        t_cut = float(t_cut.iloc[0])

        if quad == "conservation_priority":
            if x_min <= t_cut <= x_main_max:
                ax.vlines(
                    t_cut, y_main_min, min(qY, y_main_max),
                    color="#6F6F6F",
                    linestyle=(0, (4, 3)),
                    linewidth=1.15,
                    alpha=0.85,
                    zorder=2,
                    label="Top 1% quadrant cutoff" if not top1_label_used else None
                )
                top1_label_used = True
            continue

        if quad == "breeding_function_priority":
            if y_main_min <= t_cut <= y_main_max:
                ax.hlines(
                    t_cut, x_min, min(qC, x_main_max),
                    color="#6F6F6F",
                    linestyle=(0, (4, 3)),
                    linewidth=1.15,
                    alpha=0.85,
                    zorder=2,
                    label="Top 1% quadrant cutoff" if not top1_label_used else None
                )
                top1_label_used = True
            continue

        if quad == "conflict_quadrant":
            y_curve = t_cut / x_grid
            mask = np.isfinite(x_grid) & np.isfinite(y_curve)
            mask &= (y_curve >= y_main_min) & (y_curve <= y_main_max)
            mask &= (x_grid >= qC) & (y_curve >= qY)
            y_plot = np.where(mask, y_curve, np.nan)
        else:
            continue

        if np.isfinite(y_plot).sum() >= 2:
            ax.plot(
                x_grid, y_plot,
                color="#6F6F6F",
                linestyle=(0, (4, 3)),
                linewidth=1.15,
                alpha=0.85,
                zorder=2,
                label="Top 1% quadrant cutoff" if not top1_label_used else None
            )
            top1_label_used = True

    ax.set_xlabel("Conservation-risk score (C score)")
    ax.set_ylabel("Selection/function composite score [0.4×S + 0.6×M]")

    # top genes
    texts = []
    if gene_df is not None and not gene_df.empty:
        top = gene_df.head(top_n_labels).copy()
        top["best_C_score"] = pd.to_numeric(top["best_C_score"], errors="coerce")
        top["best_selection_function_score"] = pd.to_numeric(top["best_selection_function_score"], errors="coerce")

        top_main = top[top["best_C_score"] <= x_main_max].copy()

        if not top_main.empty:
            ax.scatter(
                top_main["best_C_score"],
                top_main["best_selection_function_score"],
                s=44, facecolors="none", edgecolors="black", linewidths=0.65, zorder=5
            )

            for _, r in top_main.iterrows():
                gx = r.get("best_C_score", np.nan)
                gy = r.get("best_selection_function_score", np.nan)
                glabel = first_nonempty(
                    r.get("display_name", np.nan),
                    r.get("gene_name", np.nan),
                    r.get("gene_std", np.nan)
                )

                if pd.notna(gx) and pd.notna(gy) and pd.notna(glabel):
                    txt = ax.text(
                        gx, gy, str(glabel),
                        fontsize=8.5,
                        color="black",
                        zorder=6,
                        bbox=dict(
                            boxstyle="round,pad=0.18",
                            fc="white",
                            ec="none",
                            alpha=0.80
                        )
                    )
                    texts.append(txt)

            if HAVE_ADJUSTTEXT and len(texts) > 0:
                adjust_text(
                    texts,
                    ax=ax,
                    expand_points=(1.12, 1.20),
                    expand_text=(1.08, 1.16),
                    force_points=0.22,
                    force_text=0.25,
                    arrowprops=dict(
                        arrowstyle="-",
                        color="#666666",
                        lw=0.55,
                        alpha=0.75
                    )
                )

    ax.legend(
        frameon=True,
        loc="upper right",
        bbox_to_anchor=(0.985, 0.985),
        facecolor="white",
        edgecolor="none",
        framealpha=0.78,
        handletextpad=0.3,
        borderpad=0.12,
        labelspacing=0.35,
        scatterpoints=1
    )

    plt.tight_layout()
    # Scatter layers are rasterized to keep PDF/CorelDRAW/Illustrator/PS
    # responsive, while axes, text, labels, legends and threshold lines remain
    # vector objects.  dpi controls only the rasterized point layers.
    plt.savefig(out_pdf, bbox_inches="tight", dpi=500)
    png_file = os.path.splitext(out_pdf)[0] + ".png"
    plt.savefig(png_file, dpi=500, bbox_inches="tight")
    plt.close()


# =========================================================
# Main
# =========================================================

def main():
    parser = argparse.ArgumentParser(
        description="Postprocess conflict loci results into gene-level candidate table and final publication-style quadrant plot."
    )
    parser.add_argument("--windows", required=True, help="conflict_loci_integrated_table.tsv")
    parser.add_argument("--musk-genes", required=True, help="musk_function_genes.tsv")
    parser.add_argument("--annotation", required=True, help="FMdeer_unified_gene_annotation.unique.tsv")
    parser.add_argument("--gene-coords", required=True, help="gene_coordinates.tsv")
    parser.add_argument("--outdir", required=True, help="output directory")
    parser.add_argument("--top-labels", type=int, default=5, help="number of top genes to label in main panel")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("[1/5] Loading integrated windows...")
    win = pd.read_csv(args.windows, sep="\t", low_memory=False)

    needed_num = [
        "start", "end", "C_score", "S_score", "M_score", "T_score", "high_conflict_flag"
    ]
    win = safe_numeric(win, [c for c in needed_num if c in win.columns])

    print("[2/5] Loading and unifying musk genes...")
    mg = load_and_unify_musk_genes(args.musk_genes, args.annotation, args.gene_coords)

    print("[3/5] Expanding windows to gene rows...")
    expanded = build_window_gene_expanded(win, mg)
    expanded_file = os.path.join(args.outdir, "conflict_window_gene_expanded.tsv")
    expanded.to_csv(expanded_file, sep="\t", index=False)

    print("[4/5] Building gene-level candidate table...")
    gene_df = build_top_conflict_genes(expanded)
    gene_file = os.path.join(args.outdir, "top_conflict_candidate_genes.tsv")
    gene_df.to_csv(gene_file, sep="\t", index=False)

    print("[5/5] Building quadrant table and plot...")
    quad = build_quadrant_table(win)
    quad_file = os.path.join(args.outdir, "conflict_quadrant_table.tsv")
    quad.to_csv(quad_file, sep="\t", index=False)

    quad_top1 = quad[quad["quadrant_top1_flag"].fillna(0) > 0].copy()
    quad_top1_file = os.path.join(args.outdir, "conflict_quadrant_top1_windows.tsv")
    quad_top1.to_csv(quad_top1_file, sep="\t", index=False)

    if not expanded.empty and not quad_top1.empty:
        top1_keys = quad_top1[["CHROM", "start", "end", "quadrant", "quadrant_top1_rank", "quadrant_top1_cutoff"]].copy()
        top1_keys = top1_keys.rename(columns={
            "CHROM": "window_CHROM",
            "start": "window_start",
            "end": "window_end",
            "quadrant": "top1_quadrant"
        })
        top1_keys["window_start"] = pd.to_numeric(top1_keys["window_start"], errors="coerce")
        top1_keys["window_end"] = pd.to_numeric(top1_keys["window_end"], errors="coerce")
        expanded_top1 = expanded.copy()
        expanded_top1["window_start"] = pd.to_numeric(expanded_top1["window_start"], errors="coerce")
        expanded_top1["window_end"] = pd.to_numeric(expanded_top1["window_end"], errors="coerce")
        expanded_top1 = pd.merge(
            expanded_top1,
            top1_keys,
            on=["window_CHROM", "window_start", "window_end"],
            how="inner"
        )
    else:
        expanded_top1 = expanded.iloc[0:0].copy()

    top1_gene_expanded_file = os.path.join(args.outdir, "conflict_quadrant_top1_window_gene_expanded.tsv")
    expanded_top1.to_csv(top1_gene_expanded_file, sep="\t", index=False)

    top1_gene_df = build_top_conflict_genes(expanded_top1) if not expanded_top1.empty else expanded_top1
    top1_gene_file = os.path.join(args.outdir, "conflict_quadrant_top1_candidate_genes.tsv")
    top1_gene_df.to_csv(top1_gene_file, sep="\t", index=False)

    plot_file = os.path.join(args.outdir, "conflict_quadrant_plot.pdf")
    plot_quadrant(quad, top1_gene_df, plot_file, top_n_labels=args.top_labels)

    print("[OK] Done.")
    print(f"  Expanded table : {expanded_file}")
    print(f"  Gene table     : {gene_file}")
    print(f"  Quadrant table : {quad_file}")
    print(f"  Top1 windows   : {quad_top1_file}")
    print(f"  Top1 gene rows : {top1_gene_expanded_file}")
    print(f"  Top1 genes     : {top1_gene_file}")
    if HAVE_MPL:
        print(f"  Quadrant plot  : {plot_file}")
        print(f"  Quadrant plot  : {os.path.splitext(plot_file)[0] + '.png'}")
    else:
        print("  Quadrant plot  : skipped (matplotlib not available)")


if __name__ == "__main__":
    main()
