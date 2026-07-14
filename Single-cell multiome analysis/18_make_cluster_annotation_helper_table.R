suppressPackageStartupMessages({
  library(data.table)
})

outDir <- "06.multiome/04.seurat_signac"
markerDir <- file.path(outDir, "markers")

rna_file <- file.path(markerDir, "RNA_markers_top20.unified.tsv")
act_file <- file.path(markerDir, "ACTIVITY_markers_top20.unified.tsv")
sample_file <- file.path(markerDir, "cluster_sample_distribution.tsv")

if (!file.exists(rna_file)) stop("Missing file: ", rna_file)
if (!file.exists(act_file)) stop("Missing file: ", act_file)

rna <- fread(rna_file, sep = "\t", header = TRUE, fill = TRUE)
act <- fread(act_file, sep = "\t", header = TRUE, fill = TRUE)

get_fc_col <- function(dt) {
  for (cc in c("avg_log2FC", "avg_logFC")) {
    if (cc %in% colnames(dt)) return(cc)
  }
  stop("No logFC column found.")
}

get_label_col <- function(dt) {
  for (cc in c("display_label", "display_name", "gene_name", "gene")) {
    if (cc %in% colnames(dt)) return(cc)
  }
  stop("No suitable label column found.")
}

get_desc_col <- function(dt) {
  for (cc in c("display_description", "description", "product", "eggNOG_Description")) {
    if (cc %in% colnames(dt)) return(cc)
  }
  return(NULL)
}

make_top_summary <- function(dt, prefix, topn = 5) {
  fc_col <- get_fc_col(dt)
  label_col <- get_label_col(dt)
  desc_col <- get_desc_col(dt)

  if (!"cluster" %in% colnames(dt)) {
    stop("Missing 'cluster' column.")
  }

  dt <- dt[order(cluster, -get(fc_col))]

  out <- dt[, .(
    tmp_gene  = paste(head(gene, topn), collapse = ", "),
    tmp_label = paste(head(get(label_col), topn), collapse = ", "),
    tmp_desc  = if (!is.null(desc_col)) paste(head(get(desc_col), topn), collapse = " | ") else NA_character_
  ), by = cluster]

  setnames(
    out,
    old = c("tmp_gene", "tmp_label", "tmp_desc"),
    new = c(
      paste0(prefix, "_genes"),
      paste0(prefix, "_labels"),
      paste0(prefix, "_descriptions")
    )
  )

  out
}

rna_sum <- make_top_summary(rna, prefix = "rna", topn = 5)
act_sum <- make_top_summary(act, prefix = "activity", topn = 5)

helper <- merge(rna_sum, act_sum, by = "cluster", all = TRUE)
helper[, cluster := as.character(cluster)]

if (file.exists(sample_file)) {
  samp <- fread(sample_file, sep = "\t", header = TRUE)
  samp[, cluster := as.character(cluster)]
  samp <- samp[n_cells > 0]

  if (all(c("cluster", "sample", "n_cells") %in% colnames(samp))) {
    samp_sum <- samp[, .(
      sample_distribution = paste0(sample, ":", n_cells, collapse = "; ")
    ), by = cluster]

    helper <- merge(helper, samp_sum, by = "cluster", all.x = TRUE)
  }
}

# Keep this helper table in sync with 13_featureplot_dotplot_composition.R.
# "detailed" is the cluster-level call, "minor" merges the closest epithelial
# programs, and "major" is the robust lineage used for composition analysis.
detailed_map <- c(
  "0" = "Neuronal_like_cells", "1" = "Ciliated_neuroepithelial_cells",
  "2" = "Epithelial_KRT79_positive", "3" = "Fibroblasts",
  "4" = "Stress_response_cells", "5" = "Myeloid_cells",
  "6" = "Epithelial_differentiated_WNT", "7" = "Stromal_fibroblasts",
  "8" = "Activated_stromal_cells", "9" = "Epithelial_KRT79_positive",
  "10" = "Epithelial_differentiated_KLK_SPINK5", "11" = "Lymphatic_endothelial_cells",
  "12" = "Endothelial_cells", "13" = "Lymphatic_endothelial_cells",
  "14" = "Epithelial_KRT14_basal", "15" = "Smooth_muscle_cells",
  "16" = "B_cells", "17" = "Secretory_epithelial_cells", "18" = "T_cells",
  "19" = "ACKR1_positive_endothelial_cells", "20" = "DC_like_myeloid_cells",
  "21" = "SOX10_positive_glial_like_cells", "22" = "SOX10_positive_specialized_cells",
  "23" = "Epithelial_specialized_metabolic", "24" = "Epithelial_specialized_lipid_metabolic",
  "25" = "Epithelial_specialized", "26" = "P2RY12_positive_myeloid_cells",
  "27" = "Neurotrophic_stromal_cells", "28" = "C1QB_positive_myeloid_cells",
  "29" = "Contractile_smooth_muscle_cells", "30" = "Activated_B_cells",
  "31" = "Cycling_cells", "32" = "Rare_neurovascular_like_cells",
  "33" = "Rare_mesenchymal_cells", "34" = "Epithelial_specialized_2"
)

major_map <- c(
  "0" = "Neuronal_like_cells", "1" = "Epithelial_cells", "2" = "Epithelial_cells",
  "3" = "Stromal_fibroblast_cells", "4" = "Unresolved_lineage", "5" = "Myeloid_cells",
  "6" = "Epithelial_cells", "7" = "Stromal_fibroblast_cells", "8" = "Stromal_fibroblast_cells",
  "9" = "Epithelial_cells", "10" = "Epithelial_cells", "11" = "Lymphatic_endothelial_cells",
  "12" = "Endothelial_cells", "13" = "Lymphatic_endothelial_cells", "14" = "Epithelial_cells",
  "15" = "Smooth_muscle_cells", "16" = "B_cells", "17" = "Epithelial_cells", "18" = "T_cells",
  "19" = "Endothelial_cells", "20" = "Myeloid_cells", "21" = "SOX10_positive_glial_cells",
  "22" = "SOX10_positive_glial_cells", "23" = "Epithelial_cells", "24" = "Epithelial_cells",
  "25" = "Epithelial_cells", "26" = "Myeloid_cells", "27" = "Stromal_fibroblast_cells",
  "28" = "Myeloid_cells", "29" = "Smooth_muscle_cells", "30" = "B_cells",
  "31" = "Unresolved_lineage", "32" = "Rare_neurovascular_like_cells",
  "33" = "Rare_mesenchymal_cells", "34" = "Epithelial_cells"
)

minor_map <- c(
  "0" = "Neuronal_like_cells", "1" = "Ciliated_neuroepithelial",
  "2" = "Epithelial_KRT79_positive", "9" = "Epithelial_KRT79_positive",
  "6" = "Epithelial_differentiated_keratinizing", "10" = "Epithelial_differentiated_keratinizing",
  "14" = "Epithelial_KRT14_positive_basal", "17" = "Epithelial_secretory",
  "23" = "Epithelial_specialized_or_metabolic", "24" = "Epithelial_specialized_or_metabolic",
  "25" = "Epithelial_specialized_or_metabolic", "34" = "Epithelial_specialized_or_metabolic"
)
state_map <- c("4" = "Stress_response", "31" = "Cycling")

helper[, celltype_detailed := unname(detailed_map[cluster])]
helper[, celltype_minor := unname(minor_map[cluster])]
helper[is.na(celltype_minor), celltype_minor := celltype_detailed]
helper[, celltype_major := unname(major_map[cluster])]
helper[, cell_state := unname(state_map[cluster])]
helper[is.na(cell_state), cell_state := "None_detected"]

# The proposed call is pre-filled; final_celltype and notes remain editable.
helper[, suggested_celltype := celltype_detailed]
helper[, final_celltype := celltype_major]
helper[, notes := NA_character_]

setcolorder(helper, c(
  "cluster",
  "rna_labels", "rna_genes", "rna_descriptions",
  "activity_labels", "activity_genes", "activity_descriptions",
  "sample_distribution",
  "celltype_detailed", "celltype_minor", "celltype_major", "cell_state",
  "suggested_celltype", "final_celltype", "notes"
))

outfile <- file.path(outDir, "cluster_annotation_helper_table.tsv")
fwrite(helper, outfile, sep = "\t", quote = FALSE)

message("Saved: ", outfile)
