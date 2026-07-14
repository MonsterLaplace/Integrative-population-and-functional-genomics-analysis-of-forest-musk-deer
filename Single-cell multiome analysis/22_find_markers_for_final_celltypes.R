#!/usr/bin/env Rscript
# Find markers at the merged lineage, epithelial subtype, and state levels.

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(data.table)
  library(dplyr)
})

outDir <- "06.multiome/04.seurat_signac"
inferredFile <- file.path(outDir, "FMdeer_multiome.wnn.final_annotated.state_lineage_inferred.rds")
objFile <- if (file.exists(inferredFile)) inferredFile else file.path(outDir, "FMdeer_multiome.wnn.final_annotated.rds")
if (!file.exists(objFile)) stop("Annotated object not found: ", objFile)
obj <- readRDS(objFile)

required_metadata <- c("celltype_minor", "celltype_major_inferred", "cell_state")
if (!all(required_metadata %in% colnames(obj@meta.data))) {
  stop("Annotation-level metadata missing. Run 19_apply_cluster_annotation_from_table.R first.")
}
markerDir <- file.path(outDir, "final_celltype_markers")
dir.create(markerDir, recursive = TRUE, showWarnings = FALSE)

safe_prepare_assay <- function(x, assay_name) {
  if (!assay_name %in% Assays(x)) return(NULL)
  DefaultAssay(x) <- assay_name
  x <- tryCatch({ x[[assay_name]] <- JoinLayers(x[[assay_name]]); x }, error = function(e) x)
  has_data <- tryCatch({
    dat <- GetAssayData(x, assay = assay_name, layer = "data")
    nrow(dat) > 0 && ncol(dat) > 0
  }, error = function(e) FALSE)
  if (!has_data) x <- NormalizeData(x, assay = assay_name, verbose = FALSE)
  x
}

write_marker_outputs <- function(markers, prefix, topn = 20) {
  markers <- as.data.frame(markers)
  fwrite(markers, file.path(markerDir, paste0(prefix, "_all.tsv")), sep = "\t")
  if (nrow(markers) == 0) return(invisible(NULL))
  fc_col <- intersect(c("avg_log2FC", "avg_logFC", "log2FC", "avg_diff"), colnames(markers))[1]
  if (is.na(fc_col) || !"cluster" %in% colnames(markers)) return(invisible(NULL))

  top <- markers %>% group_by(cluster) %>%
    slice_max(order_by = .data[[fc_col]], n = topn, with_ties = FALSE) %>% ungroup()
  fwrite(as.data.table(top), file.path(markerDir, paste0(prefix, "_top", topn, ".tsv")), sep = "\t")

  top10 <- markers %>% group_by(cluster) %>%
    slice_max(order_by = .data[[fc_col]], n = 10, with_ties = FALSE) %>%
    mutate(rank = row_number()) %>% ungroup()
  gene_col <- if ("gene" %in% colnames(top10)) "gene" else colnames(top10)[1]
  fwrite(dcast(as.data.table(top10), rank ~ cluster, value.var = gene_col),
         file.path(markerDir, paste0(prefix, "_top10_wide.tsv")), sep = "\t")
}

run_marker_analysis <- function(x, assay_name, group_col, prefix, cells = NULL) {
  if (!assay_name %in% Assays(x)) return(invisible(NULL))
  if (!is.null(cells)) x <- subset(x, cells = cells)
  groups <- droplevels(factor(x@meta.data[[group_col]]))
  if (length(unique(groups[!is.na(groups)])) < 2) {
    message("Skipping ", prefix, ": fewer than two groups.")
    return(invisible(NULL))
  }
  x@meta.data[[group_col]] <- groups
  Idents(x) <- group_col
  x <- safe_prepare_assay(x, assay_name)
  if (is.null(x)) return(invisible(NULL))
  DefaultAssay(x) <- assay_name
  message("Finding ", assay_name, " markers: ", prefix)
  markers <- tryCatch(FindAllMarkers(x, only.pos = TRUE, min.pct = 0.1,
                                     logfc.threshold = 0.25, test.use = "wilcox", verbose = TRUE),
                      error = function(e) { message("FindAllMarkers failed: ", e$message); data.frame() })
  write_marker_outputs(markers, prefix)
}

# Primary marker set: merged lineages. Unresolved stress/cycling-only clusters
# are excluded, since a state is not a terminal lineage.
major_cells <- rownames(obj@meta.data)[obj@meta.data$celltype_major_inferred != "Unresolved_lineage" &
                                       !is.na(obj@meta.data$celltype_major_inferred)]
for (assay_name in c("RNA", "ACTIVITY")) {
  run_marker_analysis(obj, assay_name, "celltype_major_inferred", paste0(assay_name, "_major_celltype_markers"), major_cells)
}

# Secondary marker set: epithelial variation only.
epi_cells <- rownames(obj@meta.data)[obj@meta.data$celltype_major_inferred == "Epithelial_cells" &
                                     !is.na(obj@meta.data$celltype_minor)]
for (assay_name in c("RNA", "ACTIVITY")) {
  run_marker_analysis(obj, assay_name, "celltype_minor", paste0(assay_name, "_epithelial_minor_markers"), epi_cells)
}

# State markers compare stress versus cycling cells only; "None_detected" is
# intentionally excluded because it is a heterogeneous background group.
state_cells <- rownames(obj@meta.data)[obj@meta.data$cell_state %in% c("Stress_response", "Cycling")]
for (assay_name in c("RNA", "ACTIVITY")) {
  run_marker_analysis(obj, assay_name, "cell_state", paste0(assay_name, "_cell_state_markers"), state_cells)
}

for (cc in required_metadata) {
  counts <- as.data.table(table(obj@meta.data[[cc]]))
  setnames(counts, c(cc, "n_cells"))
  setorder(counts, -n_cells)
  fwrite(counts, file.path(markerDir, paste0(cc, "_counts.tsv")), sep = "\t")
}

saveRDS(obj, file.path(outDir, "FMdeer_multiome.wnn.final_annotated.markers.rds"))
message("Done. Results saved to: ", markerDir)
