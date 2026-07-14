#!/usr/bin/env Rscript
# Assign candidate genes to the strongest marker-supported annotation at one level.

suppressPackageStartupMessages({
  library(data.table)
})

outDir <- "06.multiome/04.seurat_signac"
args <- commandArgs(trailingOnly = TRUE)
analysis_level <- if (length(args) >= 1) args[[1]] else "major"
specs <- list(
  major = "major_celltype",
  epithelial_minor = "epithelial_minor",
  cell_state = "cell_state"
)
if (!analysis_level %in% names(specs)) stop("analysis_level must be one of: ", paste(names(specs), collapse = ", "))
marker_label <- specs[[analysis_level]]

muskFile <- file.path(outDir, if (analysis_level == "major") "musk_function_genes.tsv" else paste0("musk_function_genes_", analysis_level, ".tsv"))
rnaMarkerFile <- file.path(outDir, "final_celltype_markers", paste0("RNA_", marker_label, "_markers_all.tsv"))
actMarkerFile <- file.path(outDir, "final_celltype_markers", paste0("ACTIVITY_", marker_label, "_markers_all.tsv"))
assignDir <- file.path(outDir, if (analysis_level == "major") "musk_function_gene_assignment" else paste0("musk_function_gene_assignment_", analysis_level))
dir.create(assignDir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(muskFile)) stop("Missing file: ", muskFile)
if (!file.exists(rnaMarkerFile)) stop("Missing RNA marker file: ", rnaMarkerFile)
musk <- fread(muskFile)
rna <- fread(rnaMarkerFile)
if (!"gene" %in% names(musk)) stop("Musk file missing gene column")
if (!all(c("gene", "cluster") %in% names(rna))) stop("RNA marker file missing gene/cluster columns")

get_fc_col <- function(dt) {
  x <- intersect(c("avg_log2FC", "avg_logFC", "log2FC", "avg_diff"), names(dt))[1]
  if (is.na(x)) stop("No logFC column found.")
  x
}
best_annotation <- function(dt, assay_label) {
  if (is.null(dt) || nrow(dt) == 0 || !all(c("gene", "cluster") %in% names(dt))) {
    return(data.table(gene = character(), annotation = character(), logFC = numeric()))
  }
  fc_col <- get_fc_col(dt)
  x <- copy(dt)
  x[, tmp_fc := get(fc_col)]
  setorder(x, gene, -tmp_fc)
  x <- x[, .SD[1], by = gene]
  x[, .(gene, annotation = as.character(cluster), logFC = tmp_fc, assay = assay_label)]
}

rna_best <- best_annotation(rna, "RNA_marker")
act_best <- if (file.exists(actMarkerFile)) best_annotation(fread(actMarkerFile), "ACTIVITY_marker") else
  data.table(gene = character(), annotation = character(), logFC = numeric(), assay = character())
setnames(rna_best, c("annotation", "logFC", "assay"), c("best_annotation_rna", "best_annotation_rna_logFC", "rna_source"))
setnames(act_best, c("annotation", "logFC", "assay"), c("best_annotation_activity", "best_annotation_activity_logFC", "activity_source"))

res <- merge(musk, rna_best, by = "gene", all.x = TRUE)
res <- merge(res, act_best, by = "gene", all.x = TRUE)
res[, best_annotation := fifelse(!is.na(best_annotation_rna) & best_annotation_rna != "",
                                 best_annotation_rna, best_annotation_activity)]
res[, assignment_source := fifelse(!is.na(best_annotation_rna) & best_annotation_rna != "", "RNA_marker",
                                   fifelse(!is.na(best_annotation_activity) & best_annotation_activity != "", "ACTIVITY_marker", NA_character_))]
res[, selected_analysis_level := analysis_level]

for (cc in c("bulkDE_log2FC", "bulkDE_padj", "celltype_specificity", "musk_score_raw")) {
  if (!cc %in% names(res)) res[, (cc) := NA_real_]
}
for (cc in c("scRNA_marker", "ATAC_peak_support")) if (!cc %in% names(res)) res[, (cc) := 0L]
res[, adult_up := as.integer(!is.na(bulkDE_log2FC) & !is.na(bulkDE_padj) & bulkDE_log2FC > 0 & bulkDE_padj < 0.05)]
res[, candidate_class := fifelse(adult_up == 1 & scRNA_marker == 1 & ATAC_peak_support == 1, "adult_up_scRNA_ATAC",
                          fifelse(adult_up == 1 & scRNA_marker == 1, "adult_up_scRNA",
                          fifelse(adult_up == 1, "adult_up_only", "other")))]
setorder(res, -musk_score_raw, -celltype_specificity, -bulkDE_log2FC)
fwrite(res, file.path(assignDir, "musk_function_genes.with_annotation_assignment.tsv"), sep = "\t", quote = FALSE)

assigned <- res[!is.na(best_annotation) & best_annotation != ""]
setorder(assigned, best_annotation, -musk_score_raw, -bulkDE_log2FC)
top_by_annotation <- assigned[, head(.SD, 20), by = best_annotation]
fwrite(top_by_annotation, file.path(assignDir, "top20_musk_function_genes_by_annotation.tsv"), sep = "\t", quote = FALSE)
top_wide <- top_by_annotation[, .(gene, rank = seq_len(.N)), by = best_annotation]
fwrite(dcast(top_wide, rank ~ best_annotation, value.var = "gene"),
       file.path(assignDir, "top20_musk_function_genes_by_annotation.wide.tsv"), sep = "\t", quote = FALSE)

summary_annotation <- assigned[, .(
  n_assigned_genes = .N, n_adult_up = sum(adult_up == 1, na.rm = TRUE),
  n_high_score = sum(musk_score_raw >= 2, na.rm = TRUE), median_score = median(musk_score_raw, na.rm = TRUE)
), by = best_annotation]
setorder(summary_annotation, -n_assigned_genes)
fwrite(summary_annotation, file.path(assignDir, "musk_function_gene_assignment_summary_by_annotation.tsv"), sep = "\t", quote = FALSE)

strict_dt <- res[adult_up == 1 & scRNA_marker == 1 & celltype_specificity >= 0.3]
fwrite(strict_dt, file.path(assignDir, "musk_function_genes.strict_with_annotation.tsv"), sep = "\t", quote = FALSE)
strict_top <- strict_dt[!is.na(best_annotation) & best_annotation != ""]
setorder(strict_top, best_annotation, -musk_score_raw, -bulkDE_log2FC)
fwrite(strict_top[, head(.SD, 20), by = best_annotation],
       file.path(assignDir, "top20_strict_musk_function_genes_by_annotation.tsv"), sep = "\t", quote = FALSE)

message("Done. Outputs saved to: ", assignDir)
