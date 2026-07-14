#!/usr/bin/env Rscript
# Pseudobulk edgeR: adult versus young within major lineages, epithelial subtypes, or states.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
})

outDir <- "06.multiome/04.seurat_signac"
pbDir <- file.path(outDir, "pseudobulk")

# Optional command-line argument: major (default), epithelial_minor, or cell_state.
args <- commandArgs(trailingOnly = TRUE)
analysis_level <- if (length(args) >= 1) args[[1]] else "major"
specs <- list(
  major = list(count = "RNA_pseudobulk_counts.rds", meta = "pseudobulk_metadata_major_celltype.tsv", label = "major_celltype"),
  epithelial_minor = list(count = "RNA_pseudobulk_epithelial_minor_counts.rds", meta = "pseudobulk_metadata_epithelial_minor.tsv", label = "epithelial_minor"),
  cell_state = list(count = "RNA_pseudobulk_cell_state_counts.rds", meta = "pseudobulk_metadata_cell_state.tsv", label = "cell_state")
)
if (!analysis_level %in% names(specs)) {
  stop("analysis_level must be one of: ", paste(names(specs), collapse = ", "))
}
spec <- specs[[analysis_level]]
countFile <- file.path(pbDir, spec$count)
pbMetaFile <- file.path(pbDir, spec$meta)
deDir <- file.path(pbDir, paste0("edgeR_by_", spec$label))
dir.create(deDir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(countFile)) stop("Missing count file: ", countFile, ". Run scripts 24 and 25 first.")
if (!file.exists(pbMetaFile)) stop("Missing pseudobulk design metadata: ", pbMetaFile, ". Run script 25 first.")

pb_counts <- readRDS(countFile)
pb_meta <- fread(pbMetaFile, sep = "\t")
required_cols <- c("pb_group", "sample", "annotation", "annotation_level", "n_cells", "group")
if (!all(required_cols %in% names(pb_meta))) {
  stop("Pseudobulk design metadata is missing: ", paste(setdiff(required_cols, names(pb_meta)), collapse = ", "))
}
if (any(is.na(pb_meta$group) | pb_meta$group == "")) stop("Some samples have missing age group.")

pb_meta[, group := factor(group, levels = c("young", "adult"))]
if (any(is.na(pb_meta$group))) stop("Found groups outside expected levels: young/adult")
pb_meta <- pb_meta[match(colnames(pb_counts), pb_group)]
if (anyNA(pb_meta$pb_group) || !all(pb_meta$pb_group == colnames(pb_counts))) {
  stop("Pseudobulk metadata and count-matrix columns do not match.")
}

summary_list <- list()
for (annotation_value in unique(pb_meta$annotation)) {
  message("Processing ", spec$label, ": ", annotation_value)
  meta_one <- pb_meta[annotation == annotation_value]
  grp_tab <- table(meta_one$group)
  summary_row <- data.table(
    annotation = annotation_value, annotation_level = unique(meta_one$annotation_level),
    comparison = "adult_vs_young", n_samples = nrow(meta_one),
    n_young = if ("young" %in% names(grp_tab)) unname(grp_tab[["young"]]) else 0,
    n_adult = if ("adult" %in% names(grp_tab)) unname(grp_tab[["adult"]]) else 0
  )
  if (length(grp_tab) < 2 || any(grp_tab < 2)) {
    summary_list[[annotation_value]] <- cbind(summary_row, status = "skipped_insufficient_replicates")
    next
  }

  counts_one <- pb_counts[, meta_one$pb_group, drop = FALSE]
  y <- DGEList(counts = counts_one)
  keep <- filterByExpr(y, group = meta_one$group)
  y <- y[keep, , keep.lib.sizes = FALSE]
  if (nrow(y) == 0) {
    summary_list[[annotation_value]] <- cbind(summary_row, status = "skipped_no_genes_after_filtering")
    next
  }

  y <- calcNormFactors(y)
  design <- model.matrix(~ group, data = meta_one)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = 2) # adult relative to young
  tt <- as.data.table(topTags(qlf, n = Inf)$table, keep.rownames = "gene")
  if ("logFC" %in% names(tt)) setnames(tt, "logFC", "log2FC")
  if ("FDR" %in% names(tt)) setnames(tt, "FDR", "padj")
  tt[, `:=`(annotation = annotation_value, annotation_level = unique(meta_one$annotation_level), comparison = "adult_vs_young")]

  safe_annotation <- gsub("[^A-Za-z0-9_.-]", "_", annotation_value)
  fwrite(tt, file.path(deDir, paste0("edgeR_", safe_annotation, ".tsv")), sep = "\t")
  fwrite(tt[order(padj, -abs(log2FC))][1:min(.N, 100)],
         file.path(deDir, paste0("edgeR_", safe_annotation, "_top100.tsv")), sep = "\t")
  summary_list[[annotation_value]] <- cbind(
    summary_row, status = "ok", n_genes_tested = nrow(tt),
    n_sig_fdr_0_05 = sum(tt$padj < 0.05, na.rm = TRUE),
    n_sig_fdr_0_05_logFC_1 = sum(tt$padj < 0.05 & abs(tt$log2FC) >= 1, na.rm = TRUE)
  )
}

summary_dt <- rbindlist(summary_list, fill = TRUE)
fwrite(summary_dt, file.path(deDir, paste0("edgeR_by_", spec$label, "_summary.tsv")), sep = "\t")
message("Done. Results saved to: ", deDir)
message("log2FC > 0 means higher expression in adult relative to young.")
