#!/usr/bin/env Rscript
# Build RNA pseudobulk matrices by sample at major, epithelial-minor, and state levels.

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(data.table)
  library(Matrix)
})

outDir <- "06.multiome/04.seurat_signac"
inferredFile <- file.path(outDir, "FMdeer_multiome.wnn.final_annotated.state_lineage_inferred.rds")
objFile <- if (file.exists(inferredFile)) inferredFile else file.path(outDir, "FMdeer_multiome.wnn.final_annotated.rds")
if (!file.exists(objFile)) stop("Annotated object not found: ", objFile)
obj <- readRDS(objFile)

required_columns <- c("sample", "celltype_major_inferred", "celltype_minor", "cell_state")
if (!all(required_columns %in% colnames(obj@meta.data))) {
  stop("Required annotation metadata missing. Run 19_apply_cluster_annotation_from_table.R first.")
}
if (!"RNA" %in% Assays(obj)) stop("RNA assay not found.")

pbDir <- file.path(outDir, "pseudobulk")
dir.create(pbDir, recursive = TRUE, showWarnings = FALSE)
DefaultAssay(obj) <- "RNA"
rna_counts <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "counts"),
                       error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts"))
if (is.null(rna_counts) || nrow(rna_counts) == 0 || ncol(rna_counts) == 0) stop("RNA counts matrix is empty.")

meta <- as.data.table(obj@meta.data, keep.rownames = "cell")
min_cells_per_group <- 20

aggregate_pseudobulk <- function(dt, annotation_col, prefix) {
  if (nrow(dt) == 0) {
    message("Skipping ", prefix, ": no cells available.")
    return(invisible(NULL))
  }
  dt <- copy(dt)
  dt <- dt[!is.na(get(annotation_col)) & get(annotation_col) != ""]
  dt[, annotation := as.character(get(annotation_col))]
  dt[, pb_group := paste(sample, annotation, sep = "__")]
  pb_meta <- dt[, .(n_cells = .N), by = .(pb_group, sample, annotation)]
  setorder(pb_meta, sample, annotation)
  fwrite(pb_meta, file.path(pbDir, paste0(prefix, "_group_metadata.tsv")), sep = "\t")

  pb_keep <- pb_meta[n_cells >= min_cells_per_group]
  fwrite(pb_keep, file.path(pbDir, paste0(prefix, "_group_metadata.filtered.tsv")), sep = "\t")
  if (nrow(pb_keep) == 0) {
    message("Skipping ", prefix, ": no groups have >= ", min_cells_per_group, " cells.")
    return(invisible(NULL))
  }

  dt_keep <- dt[pb_group %in% pb_keep$pb_group]
  common_cells <- intersect(colnames(rna_counts), dt_keep$cell)
  if (length(common_cells) == 0) stop("No pseudobulk cells found in RNA count matrix for ", prefix)
  group_levels <- pb_keep$pb_group
  cell_to_group <- setNames(dt_keep$pb_group, dt_keep$cell)
  group_factor <- factor(cell_to_group[common_cells], levels = group_levels)
  design <- Matrix::sparse.model.matrix(~ 0 + group_factor)
  colnames(design) <- sub("^group_factor", "", colnames(design))
  pb_counts <- as(rna_counts[, common_cells, drop = FALSE] %*% design, "dgCMatrix")

  saveRDS(pb_counts, file.path(pbDir, paste0(prefix, "_counts.rds")))
  Matrix::writeMM(pb_counts, file.path(pbDir, paste0(prefix, "_counts.mtx")))
  fwrite(data.table(gene = rownames(pb_counts)), file.path(pbDir, paste0(prefix, "_genes.tsv")), sep = "\t")
  fwrite(pb_keep[, .(pb_group, sample, annotation, n_cells)],
         file.path(pbDir, paste0(prefix, "_samples.tsv")), sep = "\t")
  fwrite(data.table(gene = rownames(pb_counts), as.data.frame(as.matrix(pb_counts))),
         file.path(pbDir, paste0(prefix, "_counts.tsv")), sep = "\t")

  fwrite(data.table(gene = rownames(pb_counts), total_counts = as.numeric(Matrix::rowSums(pb_counts))),
         file.path(pbDir, paste0(prefix, "_gene_total_counts.tsv")), sep = "\t")
  fwrite(data.table(pb_group = colnames(pb_counts), total_counts = as.numeric(Matrix::colSums(pb_counts))),
         file.path(pbDir, paste0(prefix, "_sample_total_counts.tsv")), sep = "\t")
  message("Saved ", ncol(pb_counts), " pseudobulk groups for ", prefix)
}

# Primary analysis: merged lineages. State-only clusters are not assigned a
# lineage and therefore are excluded here.
aggregate_pseudobulk(
  meta[celltype_major_inferred != "Unresolved_lineage" & !is.na(celltype_major_inferred)],
  "celltype_major_inferred", "RNA_pseudobulk"
)

# Secondary analysis: epithelial subtypes, normalized within epithelial cells.
aggregate_pseudobulk(
  meta[celltype_major_inferred == "Epithelial_cells" & !is.na(celltype_minor)],
  "celltype_minor", "RNA_pseudobulk_epithelial_minor"
)

# State analysis: stress versus cycling only, without heterogeneous non-state cells.
aggregate_pseudobulk(
  meta[cell_state %in% c("Stress_response", "Cycling")],
  "cell_state", "RNA_pseudobulk_cell_state"
)

message("Done. Pseudobulk outputs saved to: ", pbDir)
