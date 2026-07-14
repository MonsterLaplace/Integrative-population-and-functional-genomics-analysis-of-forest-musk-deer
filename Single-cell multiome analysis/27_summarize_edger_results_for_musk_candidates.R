#!/usr/bin/env Rscript
# Summarize adult-versus-young edgeR results at one selected annotation level.

suppressPackageStartupMessages({
  library(data.table)
})

outDir <- "06.multiome/04.seurat_signac"
args <- commandArgs(trailingOnly = TRUE)
analysis_level <- if (length(args) >= 1) args[[1]] else "major"
dir_labels <- c(major = "major_celltype", epithelial_minor = "epithelial_minor", cell_state = "cell_state")
if (!analysis_level %in% names(dir_labels)) {
  stop("analysis_level must be one of: ", paste(names(dir_labels), collapse = ", "))
}
deDir <- file.path(outDir, "pseudobulk", paste0("edgeR_by_", dir_labels[[analysis_level]]))
if (!dir.exists(deDir)) stop("edgeR result directory not found: ", deDir)

files <- list.files(deDir, pattern = "^edgeR_.*\\.tsv$", full.names = TRUE)
files <- files[!grepl("summary|top100", basename(files))]
if (length(files) == 0) stop("No edgeR result files found in: ", deDir)

all_list <- lapply(files, function(f) {
  dt <- fread(f)
  # Accept old result files while preferring the new standardized names.
  if ("celltype" %in% names(dt) && !"annotation" %in% names(dt)) setnames(dt, "celltype", "annotation")
  if (!"annotation_level" %in% names(dt)) dt[, annotation_level := dir_labels[[analysis_level]]]
  required <- c("gene", "log2FC", "padj", "annotation", "annotation_level", "comparison")
  if (!all(required %in% names(dt))) {
    message("Skipping file with missing required columns: ", basename(f))
    return(NULL)
  }
  dt <- dt[comparison == "adult_vs_young"]
  if (nrow(dt) == 0) return(NULL)
  dt[, source_file := basename(f)]
  dt
})
all_dt <- rbindlist(all_list, fill = TRUE)
if (nrow(all_dt) == 0) stop("No usable adult_vs_young edgeR records found.")

all_dt[, selected_analysis_level := analysis_level]
fwrite(all_dt, file.path(deDir, "edgeR_all_annotations_merged.adult_vs_young.tsv"), sep = "\t")

# For each gene, retain the strongest significant adult-up signal across
# annotations in the selected analysis level.
adult_up <- all_dt[!is.na(padj) & !is.na(log2FC) & padj < 0.05 & log2FC > 0]
if (nrow(adult_up) > 0) {
  setorder(adult_up, gene, padj, -log2FC)
  best_dt <- adult_up[, .SD[1], by = gene][, .(
    gene,
    bulkDE_log2FC = log2FC,
    bulkDE_padj = padj,
    bulkDE_annotation = annotation,
    bulkDE_annotation_level = annotation_level,
    bulkDE_comparison = comparison
  )]
} else {
  best_dt <- data.table(
    gene = character(), bulkDE_log2FC = numeric(), bulkDE_padj = numeric(),
    bulkDE_annotation = character(), bulkDE_annotation_level = character(), bulkDE_comparison = character()
  )
}

bulk_summary <- merge(data.table(gene = unique(all_dt$gene)), best_dt, by = "gene", all.x = TRUE, sort = FALSE)
bulk_summary[, selected_analysis_level := analysis_level]
fwrite(bulk_summary, file.path(deDir, "bulkDE_gene_level_summary.adult_vs_young.tsv"), sep = "\t")

message("Saved results in: ", deDir)
message("Interpretation: bulkDE_log2FC > 0 means higher expression in adult relative to young.")
