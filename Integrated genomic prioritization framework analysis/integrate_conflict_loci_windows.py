#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np

# =========================================================
# Utility functions
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

def row_nanmean(df_or_array):
    """
    Safe row-wise mean ignoring NaN without RuntimeWarning.
    Returns NaN when all values in a row are NaN.
    This is intentionally used for score integration: NA means the window
    lacks that evidence type, not that the whole window should be discarded.
    """
    if isinstance(df_or_array, pd.DataFrame):
        return df_or_array.mean(axis=1, skipna=True)
    return pd.DataFrame(df_or_array).mean(axis=1, skipna=True)

def normalize_window_coords(df, chrom_col="CHROM", start_col="start", end_col="end", start_is_1based=True):
    df = df.copy()
    df = df.rename(columns={chrom_col: "CHROM", start_col: "start", end_col: "end"})
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    if start_is_1based:
        df["start"] = df["start"] - 1
    return df

# =========================================================
# Load window-level tables
# =========================================================

def load_selection_table(selection_file):
    df = pd.read_csv(selection_file, sep="\t", low_memory=False)
    df = normalize_window_coords(df, chrom_col="CHROM", start_col="start", end_col="end", start_is_1based=True)

    numeric_cols = [
        "pi_domestic","pi_wild","pi_ratio","log2_pi_ratio","Fst","XPCLR","XPCLR_NORM",
        "XPEHH_mean","XPEHH_mean_abs","XPEHH_max","XPEHH_min","XPEHH_max_abs",
        "iHS_mean","iHS_mean_abs","iHS_max","iHS_min","iHS_max_abs",
        "z_FstPi","z_Fst","z_pi_ratio","z_XPCLR","z_XPEHH_mean_abs","z_iHS_mean_abs",
        "integrated_score","top1_integrated_score",
        "hit_FstPi","hit_Fst","hit_piRatio","hit_XPCLR","hit_XPEHH","hit_iHS","hit_sum","selected_candidate"
    ]
    df = safe_numeric(df, [c for c in numeric_cols if c in df.columns])
    return df

def load_deleterious_table(del_file):
    df = pd.read_csv(del_file, sep="\t", low_memory=False)
    df = normalize_window_coords(df, chrom_col="CHROM", start_col="start", end_col="end", start_is_1based=True)

    numeric_cols = [
        "window_size","LoF_count","severe_missense_count","missense_all_count","synonymous_count",
        "total_del_count","LoF_density","severe_density","missense_all_density","synonymous_density","del_to_syn_ratio"
    ]
    df = safe_numeric(df, [c for c in numeric_cols if c in df.columns])
    return df

def load_roh_table(roh_file):
    df = pd.read_csv(roh_file, sep="\t", low_memory=False)
    df = normalize_window_coords(df, chrom_col="CHROM", start_col="start", end_col="end", start_is_1based=False)

    numeric_cols = ["roh_freq", "mean_roh_count"]
    df = safe_numeric(df, [c for c in numeric_cols if c in df.columns])
    return df

def load_tajima_table(taj_file):
    df = pd.read_csv(taj_file, sep="\t", low_memory=False)
    # TajimaD windows are already aligned with master coordinates
    df = normalize_window_coords(df, chrom_col="chrom", start_col="start", end_col="end", start_is_1based=False)

    numeric_cols = ["tajimaD_domestic", "tajimaD_wild"]
    df = safe_numeric(df, [c for c in numeric_cols if c in df.columns])
    return df

def load_pixy_pi_table(pi_file):
    df = pd.read_csv(pi_file, sep="\t", low_memory=False)
    df = df.rename(columns={"chromosome":"CHROM", "window_pos_1":"start", "window_pos_2":"end"})
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df["avg_pi"] = pd.to_numeric(df["avg_pi"], errors="coerce")

    dom = df[df["pop"] == "domestic"][["CHROM","start","end","avg_pi"]].rename(columns={"avg_pi":"pixy_pi_domestic"})
    wild = df[df["pop"] == "wild"][["CHROM","start","end","avg_pi"]].rename(columns={"avg_pi":"pixy_pi_wild"})
    out = pd.merge(dom, wild, on=["CHROM","start","end"], how="outer")
    out["pixy_pi_ratio"] = out["pixy_pi_wild"] / out["pixy_pi_domestic"]
    return out

# =========================================================
# Load and unify musk genes with multiple ID strategies
# =========================================================

def build_annotation_lookup(anno):
    """
    Build multiple lookup tables from unified annotation:
    - gene_id_dash
    - gene_id
    - gene_name
    - display_name
    """
    anno = anno.copy()

    keep_cols = [c for c in [
        "gene_id", "gene_id_dash", "gene_name", "display_name", "display_description", "gene_biotype"
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
    """
    Try to resolve a gene token by multiple strategies:
    1) gene_id_dash
    2) gene_id
    3) gene_name
    4) display_name
    """
    if pd.isna(gene_value):
        return None

    g = str(gene_value).strip()
    if g == "":
        return None

    # 1. gene_id_dash
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

    # 2. gene_id
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

    # 3. gene_name
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

    # 4. display_name
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

def load_and_unify_musk_genes(musk_file, anno_file, coord_file):
    musk = pd.read_csv(musk_file, sep="\t", low_memory=False)
    anno = pd.read_csv(anno_file, sep="\t", low_memory=False)
    coords = pd.read_csv(coord_file, sep="\t", low_memory=False)

    lookup_tables = build_annotation_lookup(anno)

    musk = musk.rename(columns={"gene":"gene_input"}).copy()

    resolved = musk["gene_input"].apply(lambda x: resolve_gene_identifier(x, lookup_tables))
    resolved_df = pd.DataFrame(list(resolved))

    mg = pd.concat([musk.reset_index(drop=True), resolved_df.reset_index(drop=True)], axis=1)

    # coordinates: use standardized gene_id
    coords = coords.rename(columns={"gene":"gene_std"}).copy()
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

    # compute significance metrics
    mg["bulk_age_sig"] = neglog10_p(mg["bulkDE_padj"]) if "bulkDE_padj" in mg.columns else np.nan
    mg["bulk_tissue_sig"] = neglog10_p(mg["bulk_tissue_padj"]) if "bulk_tissue_padj" in mg.columns else np.nan
    mg["bulk_interaction_sig"] = neglog10_p(mg["bulk_interaction_padj"]) if "bulk_interaction_padj" in mg.columns else np.nan

    return mg

# =========================================================
# Map functional genes to windows
# =========================================================

def map_musk_genes_to_windows(master, mg):
    master = master.copy()
    records = []

    master["start"] = pd.to_numeric(master["start"], errors="coerce")
    master["end"] = pd.to_numeric(master["end"], errors="coerce")
    n_bad_windows = master[["CHROM", "start", "end"]].isna().any(axis=1).sum()
    if n_bad_windows > 0:
        print(f"[WARN] Skip {n_bad_windows} windows with missing CHROM/start/end only during musk-gene mapping.")
    master = master.dropna(subset=["CHROM", "start", "end"]).copy()
    master = master[master["start"].le(master["end"])].copy()

    mg2 = mg.dropna(subset=["CHROM","start","end"]).copy()
    mg2["start"] = pd.to_numeric(mg2["start"], errors="coerce")
    mg2["end"] = pd.to_numeric(mg2["end"], errors="coerce")
    mg2 = mg2.dropna(subset=["CHROM", "start", "end"]).copy()
    mg2 = mg2[mg2["start"].le(mg2["end"])].copy()

    for chrom in master["CHROM"].dropna().unique():
        win_chr = master[master["CHROM"] == chrom].copy()
        gene_chr = mg2[mg2["CHROM"] == chrom].copy()

        for _, w in win_chr.iterrows():
            if pd.isna(w["start"]) or pd.isna(w["end"]):
                continue
            wstart = int(w["start"])
            wend = int(w["end"])
            hits = gene_chr[(gene_chr["end"] >= wstart) & (gene_chr["start"] <= wend)].copy()

            if hits.empty:
                records.append({
                    "CHROM": w["CHROM"],
                    "start": w["start"],
                    "end": w["end"],
                    "musk_gene_count": 0,
                    "matched_gene_count": 0,
                    "unresolved_gene_count": 0,

                    "scRNA_support_count": 0,
                    "ATAC_support_count": 0,
                    "motif_support_count": 0,
                    "peak2gene_support_count": 0,
                    "adult_up_count": 0,

                    "max_musk_score_raw": np.nan,
                    "mean_musk_score_raw": np.nan,
                    "max_celltype_specificity": np.nan,
                    "max_peak2gene_score": np.nan,

                    "max_bulk_age_log2FC": np.nan,
                    "min_bulk_age_padj": np.nan,
                    "max_bulk_age_sig": np.nan,
                    "max_mature_bulk_score": np.nan,

                    "max_bulk_tissue_log2FC": np.nan,
                    "min_bulk_tissue_padj": np.nan,
                    "max_bulk_tissue_sig": np.nan,
                    "max_tissue_bonus": np.nan,

                    "max_bulk_interaction_log2FC": np.nan,
                    "min_bulk_interaction_padj": np.nan,
                    "max_bulk_interaction_sig": np.nan,
                    "max_interaction_bonus": np.nan,

                    "genes_in_window": "",
                    "gene_input_in_window": "",
                    "match_strategies_in_window": "",
                    "candidate_classes_in_window": ""
                })
            else:
                label_candidates = []
                for _, r in hits.iterrows():
                    label = r.get("display_name", np.nan)
                    if pd.isna(label) or str(label).strip() == "":
                        label = r.get("gene_name", np.nan)
                    if pd.isna(label) or str(label).strip() == "":
                        label = r.get("gene_std", np.nan)
                    if pd.notna(label):
                        label_candidates.append(str(label))

                cand_classes = hits["candidate_class"].dropna().astype(str).unique().tolist() if "candidate_class" in hits.columns else []
                match_strategies = hits["match_strategy"].dropna().astype(str).unique().tolist() if "match_strategy" in hits.columns else []

                records.append({
                    "CHROM": w["CHROM"],
                    "start": w["start"],
                    "end": w["end"],

                    "musk_gene_count": hits["gene_std"].nunique(),
                    "matched_gene_count": hits["match_strategy"].ne("unresolved").sum() if "match_strategy" in hits.columns else 0,
                    "unresolved_gene_count": hits["match_strategy"].eq("unresolved").sum() if "match_strategy" in hits.columns else 0,

                    "scRNA_support_count": hits["scRNA_marker"].fillna(0).gt(0).sum() if "scRNA_marker" in hits.columns else 0,
                    "ATAC_support_count": hits["ATAC_peak_support"].fillna(0).gt(0).sum() if "ATAC_peak_support" in hits.columns else 0,
                    "motif_support_count": hits["motif_support"].fillna(0).gt(0).sum() if "motif_support" in hits.columns else 0,
                    "peak2gene_support_count": hits["peak2gene_support"].fillna(0).gt(0).sum() if "peak2gene_support" in hits.columns else 0,
                    "adult_up_count": hits["adult_up"].fillna(0).gt(0).sum() if "adult_up" in hits.columns else 0,

                    "max_musk_score_raw": pd.to_numeric(hits["musk_score_raw"], errors="coerce").max() if "musk_score_raw" in hits.columns else np.nan,
                    "mean_musk_score_raw": pd.to_numeric(hits["musk_score_raw"], errors="coerce").mean() if "musk_score_raw" in hits.columns else np.nan,
                    "max_celltype_specificity": pd.to_numeric(hits["celltype_specificity"], errors="coerce").max() if "celltype_specificity" in hits.columns else np.nan,
                    "max_peak2gene_score": pd.to_numeric(hits["peak2gene_max_score"], errors="coerce").max() if "peak2gene_max_score" in hits.columns else np.nan,

                    "max_bulk_age_log2FC": pd.to_numeric(hits["bulkDE_log2FC"], errors="coerce").max() if "bulkDE_log2FC" in hits.columns else np.nan,
                    "min_bulk_age_padj": pd.to_numeric(hits["bulkDE_padj"], errors="coerce").min() if "bulkDE_padj" in hits.columns else np.nan,
                    "max_bulk_age_sig": pd.to_numeric(hits["bulk_age_sig"], errors="coerce").max() if "bulk_age_sig" in hits.columns else np.nan,
                    "max_mature_bulk_score": pd.to_numeric(hits["mature_bulk_score"], errors="coerce").max() if "mature_bulk_score" in hits.columns else np.nan,

                    "max_bulk_tissue_log2FC": pd.to_numeric(hits["bulk_tissue_log2FC"], errors="coerce").max() if "bulk_tissue_log2FC" in hits.columns else np.nan,
                    "min_bulk_tissue_padj": pd.to_numeric(hits["bulk_tissue_padj"], errors="coerce").min() if "bulk_tissue_padj" in hits.columns else np.nan,
                    "max_bulk_tissue_sig": pd.to_numeric(hits["bulk_tissue_sig"], errors="coerce").max() if "bulk_tissue_sig" in hits.columns else np.nan,
                    "max_tissue_bonus": pd.to_numeric(hits["tissue_bonus"], errors="coerce").max() if "tissue_bonus" in hits.columns else np.nan,

                    "max_bulk_interaction_log2FC": pd.to_numeric(hits["bulk_interaction_log2FC"], errors="coerce").max() if "bulk_interaction_log2FC" in hits.columns else np.nan,
                    "min_bulk_interaction_padj": pd.to_numeric(hits["bulk_interaction_padj"], errors="coerce").min() if "bulk_interaction_padj" in hits.columns else np.nan,
                    "max_bulk_interaction_sig": pd.to_numeric(hits["bulk_interaction_sig"], errors="coerce").max() if "bulk_interaction_sig" in hits.columns else np.nan,
                    "max_interaction_bonus": pd.to_numeric(hits["interaction_bonus"], errors="coerce").max() if "interaction_bonus" in hits.columns else np.nan,

                    "genes_in_window": ",".join(sorted(set(label_candidates))),
                    "gene_input_in_window": ",".join(sorted(set(hits["gene_input"].astype(str).tolist()))) if "gene_input" in hits.columns else "",
                    "match_strategies_in_window": ",".join(sorted(set(match_strategies))),
                    "candidate_classes_in_window": ",".join(sorted(set(cand_classes)))
                })

    return pd.DataFrame(records)

# =========================================================
# Classification
# =========================================================

def classify_windows(df):
    qC = df["C_score"].quantile(0.90)
    qS = df["S_score"].quantile(0.90)
    qM = df["M_score"].quantile(0.90)

    df["class"] = "background"

    df.loc[
        (df["C_score"] >= qC) & ((df["S_score"] >= qS) | (df["M_score"] >= qM)),
        "class"
    ] = "conflict"

    df.loc[
        (df["C_score"] >= qC) & (df["S_score"] < qS) & (df["M_score"] < qM),
        "class"
    ] = "conservation_priority"

    df.loc[
        (df["C_score"] < qC) & ((df["S_score"] >= qS) | (df["M_score"] >= qM)),
        "class"
    ] = "breeding_priority"

    return df

# =========================================================
# Main
# =========================================================

def main():
    parser = argparse.ArgumentParser(description="Integrate protection, selection and musk functional genomics evidence into conflict loci windows.")
    parser.add_argument("--deleterious", required=True, help="deleterious_window_density.tsv")
    parser.add_argument("--selection", required=True, help="selection_integrated_windows.tsv")
    parser.add_argument("--musk-genes", required=True, help="musk_function_genes.tsv")
    parser.add_argument("--roh", required=True, help="roh_window_frequency.tsv")
    parser.add_argument("--tajima", required=True, help="tajimaD_50kb.tsv")
    parser.add_argument("--pixy-pi", required=False, help="FMdeer_pixy.merged_pi.txt (optional)")
    parser.add_argument("--annotation", required=True, help="FMdeer_unified_gene_annotation.unique.tsv")
    parser.add_argument("--gene-coords", required=True, help="gene_coordinates.tsv")
    parser.add_argument("-o", "--output", required=True, help="output integrated table")
    args = parser.parse_args()

    # 1. master table
    master = load_selection_table(args.selection)

    # 2. merge protection-related data
    dele = load_deleterious_table(args.deleterious)
    master = pd.merge(master, dele, on=["CHROM","start","end"], how="left")

    roh = load_roh_table(args.roh)
    master = pd.merge(master, roh, on=["CHROM","start","end"], how="left")

    taj = load_tajima_table(args.tajima)
    master = pd.merge(master, taj, on=["CHROM","start","end"], how="left")

    if args.pixy_pi:
        pixy = load_pixy_pi_table(args.pixy_pi)
        master = pd.merge(master, pixy, on=["CHROM","start","end"], how="left")

    # 3. functional genes
    mg = load_and_unify_musk_genes(args.musk_genes, args.annotation, args.gene_coords)
    funcwin = map_musk_genes_to_windows(master[["CHROM","start","end"]].drop_duplicates().copy(), mg)
    master = pd.merge(master, funcwin, on=["CHROM","start","end"], how="left")

    # -----------------------------------------------------
    # C-score: conservation risk
    # -----------------------------------------------------
    master["C_pi_low"] = zscore(-pd.to_numeric(master["pi_domestic"], errors="coerce"))
    master["C_roh"] = zscore(pd.to_numeric(master.get("roh_freq", np.nan), errors="coerce"))
    master["C_lof"] = zscore(pd.to_numeric(master.get("LoF_density", np.nan), errors="coerce"))
    master["C_severe"] = zscore(pd.to_numeric(master.get("severe_density", np.nan), errors="coerce"))
    master["C_del_ratio"] = zscore(pd.to_numeric(master.get("del_to_syn_ratio", np.nan), errors="coerce"))

    c_cols = ["C_pi_low","C_roh","C_lof","C_severe","C_del_ratio"]
    master["C_score"] = row_nanmean(master[c_cols])

    # -----------------------------------------------------
    # S-score: selection
    # -----------------------------------------------------
    # New selection tables merge Fst and pi into ONE joint signal:
    #   z_FstPi = mean(z(Fst), z(abs(log2_pi_ratio)))
    # Use z_FstPi as a single component so Fst and pi are not double counted.
    # If an older selection table is supplied, compute the same joint component
    # on the fly as a backward-compatible fallback.
    if "z_FstPi" not in master.columns:
        if "Fst" in master.columns and "log2_pi_ratio" in master.columns:
            master["z_FstPi"] = (
                zscore(master["Fst"]) +
                zscore(pd.to_numeric(master["log2_pi_ratio"], errors="coerce").abs())
            ) / 2
        elif "Fst" in master.columns and "pi_ratio" in master.columns:
            master["z_FstPi"] = (
                zscore(master["Fst"]) +
                zscore(np.log2(pd.to_numeric(master["pi_ratio"], errors="coerce")).abs())
            ) / 2
        else:
            master["z_FstPi"] = np.nan

    if "z_XPCLR" not in master.columns:
        if "XPCLR_NORM" in master.columns:
            master["z_XPCLR"] = zscore(master["XPCLR_NORM"])
        elif "XPCLR" in master.columns:
            master["z_XPCLR"] = zscore(master["XPCLR"])
        else:
            master["z_XPCLR"] = np.nan
    if "z_XPEHH_mean_abs" not in master.columns:
        master["z_XPEHH_mean_abs"] = zscore(master["XPEHH_mean_abs"]) if "XPEHH_mean_abs" in master.columns else np.nan
    if "z_iHS_mean_abs" not in master.columns:
        master["z_iHS_mean_abs"] = zscore(master["iHS_mean_abs"]) if "iHS_mean_abs" in master.columns else np.nan

    s_cols = [c for c in ["z_FstPi","z_XPCLR","z_XPEHH_mean_abs","z_iHS_mean_abs"] if c in master.columns]
    # NA in one selection statistic only means this evidence is unavailable
    # for that window.  Other available statistics are still used.
    master["S_n_components_used"] = master[s_cols].notna().sum(axis=1)
    master["S_score"] = row_nanmean(master[s_cols])

    # -----------------------------------------------------
    # M-score: musk-related function (deduplicated version)
    # Keep M_musk_score only for inspection, but DO NOT use it
    # in final M_score to avoid double counting.
    # -----------------------------------------------------
    master["M_gene_count"] = zscore(pd.to_numeric(master.get("musk_gene_count", np.nan), errors="coerce"))

    # inspection only
    master["M_musk_score"] = zscore(pd.to_numeric(master.get("max_musk_score_raw", np.nan), errors="coerce"))

    # non-bulk functional evidence
    master["M_scRNA"] = zscore(pd.to_numeric(master.get("scRNA_support_count", np.nan), errors="coerce"))
    master["M_ATAC"] = zscore(pd.to_numeric(master.get("ATAC_support_count", np.nan), errors="coerce"))
    master["M_motif"] = zscore(pd.to_numeric(master.get("motif_support_count", np.nan), errors="coerce"))
    master["M_peak2gene"] = zscore(pd.to_numeric(master.get("peak2gene_support_count", np.nan), errors="coerce"))
    master["M_celltype_spec"] = zscore(pd.to_numeric(master.get("max_celltype_specificity", np.nan), errors="coerce"))
    master["M_peak2gene_score"] = zscore(pd.to_numeric(master.get("max_peak2gene_score", np.nan), errors="coerce"))

    master["M_nonbulk"] = row_nanmean(
        master[[
            "M_scRNA",
            "M_ATAC",
            "M_motif",
            "M_peak2gene",
            "M_celltype_spec",
            "M_peak2gene_score"
        ]]
    )

    # bulk age evidence
    master["M_adult_up"] = zscore(pd.to_numeric(master.get("adult_up_count", np.nan), errors="coerce"))
    master["M_bulk_age_fc"] = zscore(pd.to_numeric(master.get("max_bulk_age_log2FC", np.nan), errors="coerce"))
    master["M_bulk_age_sig"] = zscore(pd.to_numeric(master.get("max_bulk_age_sig", np.nan), errors="coerce"))
    master["M_bulk_age_score"] = zscore(pd.to_numeric(master.get("max_mature_bulk_score", np.nan), errors="coerce"))

    master["M_bulk_age"] = row_nanmean(
        master[[
            "M_adult_up",
            "M_bulk_age_fc",
            "M_bulk_age_sig",
            "M_bulk_age_score"
        ]]
    )

    # bulk tissue + interaction context evidence
    master["M_bulk_tissue_fc"] = zscore(pd.to_numeric(master.get("max_bulk_tissue_log2FC", np.nan), errors="coerce"))
    master["M_bulk_tissue_sig"] = zscore(pd.to_numeric(master.get("max_bulk_tissue_sig", np.nan), errors="coerce"))
    master["M_bulk_tissue_bonus"] = zscore(pd.to_numeric(master.get("max_tissue_bonus", np.nan), errors="coerce"))

    master["M_bulk_inter_fc"] = zscore(pd.to_numeric(master.get("max_bulk_interaction_log2FC", np.nan), errors="coerce"))
    master["M_bulk_inter_sig"] = zscore(pd.to_numeric(master.get("max_bulk_interaction_sig", np.nan), errors="coerce"))
    master["M_bulk_inter_bonus"] = zscore(pd.to_numeric(master.get("max_interaction_bonus", np.nan), errors="coerce"))

    master["M_bulk_context"] = row_nanmean(
        master[[
            "M_bulk_tissue_fc",
            "M_bulk_tissue_sig",
            "M_bulk_tissue_bonus",
            "M_bulk_inter_fc",
            "M_bulk_inter_sig",
            "M_bulk_inter_bonus"
        ]]
    )

    master["M_score"] = (
        0.15 * master["M_gene_count"] +
        0.45 * master["M_nonbulk"] +
        0.25 * master["M_bulk_age"] +
        0.15 * master["M_bulk_context"]
    )

    # -----------------------------------------------------
    # T-score
    # -----------------------------------------------------
    master["T_score"] = master["C_score"] * (0.4 * master["S_score"] + 0.6 * master["M_score"])

    # class
    master = classify_windows(master)

    # support flags
    master["high_conflict_flag"] = 0
    if "top1_integrated_score" in master.columns:
        master.loc[
            (master["class"] == "conflict") &
            (pd.to_numeric(master["top1_integrated_score"], errors="coerce") > 0),
            "high_conflict_flag"
        ] = 1

    master["has_musk_function_gene"] = master["musk_gene_count"].fillna(0).gt(0).astype(int)

    # sort
    master = master.sort_values(
        by=["high_conflict_flag","T_score","C_score","S_score","M_score"],
        ascending=[False, False, False, False, False],
        na_position="last"
    )

    master.to_csv(args.output, sep="\t", index=False)
    print(f"[OK] Integrated conflict loci table written to: {args.output}")

if __name__ == "__main__":
    main()
