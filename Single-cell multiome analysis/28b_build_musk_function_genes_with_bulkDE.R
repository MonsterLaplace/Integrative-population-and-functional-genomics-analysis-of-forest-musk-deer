#!/usr/bin/env Rscript
# Integrate single-cell, chromatin, and pseudobulk evidence into candidate genes.

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
})

outDir <- "06.multiome/04.seurat_signac"
inferredFile <- file.path(outDir, "FMdeer_multiome.wnn.final_annotated.state_lineage_inferred.rds")
objFile <- if (file.exists(inferredFile)) inferredFile else file.path(outDir, "FMdeer_multiome.wnn.final_annotated.rds")
args <- commandArgs(trailingOnly = TRUE)
analysis_level <- if (length(args) >= 1) args[[1]] else "major"
specs <- list(
  major = list(group_col = "celltype_major_inferred", marker_label = "major_celltype", de_label = "major_celltype"),
  epithelial_minor = list(group_col = "celltype_minor", marker_label = "epithelial_minor", de_label = "epithelial_minor"),
  cell_state = list(group_col = "cell_state", marker_label = "cell_state", de_label = "cell_state")
)
if (!analysis_level %in% names(specs)) stop("analysis_level must be one of: ", paste(names(specs), collapse = ", "))
spec <- specs[[analysis_level]]

rnaMarkerFile <- file.path(outDir, "final_celltype_markers", paste0("RNA_", spec$marker_label, "_markers_all.tsv"))
actMarkerFile <- file.path(outDir, "final_celltype_markers", paste0("ACTIVITY_", spec$marker_label, "_markers_all.tsv"))
bulkSummaryFile <- file.path(outDir, "pseudobulk", paste0("edgeR_by_", spec$de_label), "bulkDE_gene_level_summary.adult_vs_young.tsv")
outFile <- file.path(outDir, if (analysis_level == "major") "musk_function_genes.tsv" else paste0("musk_function_genes_", analysis_level, ".tsv"))

# Optional external evidence tables; their paths may be changed if available.
peak2geneFile <- "peak2gene_links.tsv"
motifSupportFile <- "motif_support.tsv"
if (!file.exists(objFile)) stop("Missing object: ", objFile)
if (!file.exists(rnaMarkerFile)) stop("Missing RNA marker file: ", rnaMarkerFile, ". Run script 22 first.")

obj <- readRDS(objFile)
if (!spec$group_col %in% colnames(obj@meta.data)) stop("Object metadata missing: ", spec$group_col)
if (!"RNA" %in% Assays(obj)) stop("RNA assay not found.")
Idents(obj) <- spec$group_col

get_fc_col <- function(dt) {
  out <- intersect(c("avg_log2FC", "avg_logFC", "log2FC"), names(dt))[1]
  if (is.na(out)) stop("Marker file has no logFC column.")
  out
}
marker_flag <- function(marker_file, flag_name) {
  if (!file.exists(marker_file)) {
    out <- data.table(gene = character())
    out[, (flag_name) := integer()]
    return(out)
  }
  dt <- fread(marker_file)
  if (!"gene" %in% names(dt)) stop("Marker file missing gene: ", marker_file)
  fc_col <- get_fc_col(dt)
  padj_col <- intersect(c("p_val_adj", "padj", "FDR"), names(dt))[1]
  dt <- dt[get(fc_col) >= 0.25]
  if (!is.na(padj_col)) dt <- dt[get(padj_col) <= 0.05]
  out <- unique(dt[, .(gene)])
  out[, (flag_name) := 1L]
  out
}

rna_markers <- fread(rnaMarkerFile)
rna_marker_flag <- marker_flag(rnaMarkerFile, "scRNA_marker")
atac_support <- marker_flag(actMarkerFile, "ATAC_peak_support")

# Cell-type specificity is calculated within the selected analysis level.
DefaultAssay(obj) <- "RNA"
expr_mat <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "data"),
                     error = function(e) GetAssayData(obj, assay = "RNA", slot = "data"))
labels <- as.character(Idents(obj))
label_levels <- unique(labels[!is.na(labels)])
avg_expr_list <- lapply(label_levels, function(label) Matrix::rowMeans(expr_mat[, labels == label, drop = FALSE]))
names(avg_expr_list) <- label_levels
avg_expr_dt <- as.data.table(avg_expr_list)
avg_expr_dt[, gene := rownames(expr_mat)]
expr_cols <- setdiff(names(avg_expr_dt), "gene")
avg_expr_dt[, max_expr := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = expr_cols]
avg_expr_dt[, sum_expr := rowSums(.SD), .SDcols = expr_cols]
specificity_dt <- avg_expr_dt[, .(gene, celltype_specificity = fifelse(sum_expr > 0, max_expr / sum_expr, 0))]

if (file.exists(bulkSummaryFile)) {
  bulk_dt <- fread(bulkSummaryFile)
  required_bulk <- c("gene", "bulkDE_log2FC", "bulkDE_padj")
  if (!all(required_bulk %in% names(bulk_dt))) stop("Bulk summary missing required columns.")
  if (!"bulkDE_annotation" %in% names(bulk_dt)) bulk_dt[, bulkDE_annotation := NA_character_]
  if (!"bulkDE_annotation_level" %in% names(bulk_dt)) bulk_dt[, bulkDE_annotation_level := spec$group_col]
} else {
  bulk_dt <- data.table(gene = character(), bulkDE_log2FC = numeric(), bulkDE_padj = numeric(),
                        bulkDE_annotation = character(), bulkDE_annotation_level = character())
}

load_binary_support <- function(file, column) {
  if (!file.exists(file)) {
    out <- data.table(gene = character())
    out[, (column) := integer()]
    return(out)
  }
  dt <- fread(file)
  if (!"gene" %in% names(dt)) stop("Evidence table missing gene: ", file)
  if (!column %in% names(dt)) dt <- unique(dt[, .(gene)])[, (column) := 1L] else dt <- unique(dt[, .(gene, get(column))])
  setnames(dt, names(dt)[2], column)
  dt
}
p2g_dt <- load_binary_support(peak2geneFile, "peak2gene_support")
motif_dt <- load_binary_support(motifSupportFile, "motif_support")

gene_universe <- unique(c(rownames(expr_mat), rna_markers$gene, atac_support$gene, bulk_dt$gene, p2g_dt$gene, motif_dt$gene))
res <- data.table(gene = gene_universe)
for (dt in list(rna_marker_flag, bulk_dt, specificity_dt, atac_support, p2g_dt, motif_dt)) {
  res <- merge(res, dt, by = "gene", all.x = TRUE, sort = FALSE)
}
for (cc in c("scRNA_marker", "ATAC_peak_support", "peak2gene_support", "motif_support")) res[is.na(get(cc)), (cc) := 0L]
res[is.na(celltype_specificity), celltype_specificity := 0]
res[, bulk_support_score := 0]
res[!is.na(bulkDE_log2FC) & !is.na(bulkDE_padj) & bulkDE_padj < 0.05 & bulkDE_log2FC > 0,
    bulk_support_score := pmin(bulkDE_log2FC, 3)]
res[, selected_analysis_level := analysis_level]
res[, musk_score_raw := 1.0 * scRNA_marker + 1.0 * bulk_support_score + 1.5 * celltype_specificity +
      0.5 * ATAC_peak_support + 0.8 * peak2gene_support + 0.5 * motif_support]
setorder(res, -musk_score_raw, -bulk_support_score, -celltype_specificity)

wanted <- c("gene", "selected_analysis_level", "scRNA_marker", "bulkDE_log2FC", "bulkDE_padj",
            "bulkDE_annotation", "bulkDE_annotation_level", "celltype_specificity", "ATAC_peak_support",
            "peak2gene_support", "motif_support", "musk_score_raw")
setcolorder(res, c(wanted, setdiff(names(res), wanted)))
fwrite(res, outFile, sep = "\t", quote = FALSE)
message("Done. Saved: ", outFile)
message("bulkDE_log2FC > 0 means higher expression in adult relative to young.")
