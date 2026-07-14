#!/usr/bin/env Rscript
# Plot marker panels created by script 22 at major, epithelial-minor, and state levels.

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

outDir <- "06.multiome/04.seurat_signac"
objFile <- file.path(outDir, "FMdeer_multiome.wnn.final_annotated.markers.rds")
markerDir <- file.path(outDir, "final_celltype_markers")
if (!file.exists(objFile)) stop("Object not found: ", objFile)
obj <- readRDS(objFile)

required_metadata <- c("celltype_minor", "celltype_major_inferred", "cell_state")
if (!all(required_metadata %in% colnames(obj@meta.data))) {
  stop("Annotation-level metadata missing. Run 19_apply_cluster_annotation_from_table.R first.")
}
plotDir <- file.path(outDir, "final_celltype_marker_plots")
dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)

get_existing_genes <- function(x, genes, assay) {
  if (!assay %in% Assays(x)) return(character(0))
  genes[genes %in% rownames(x[[assay]])]
}
get_top_markers <- function(marker_file, topn = 3) {
  if (!file.exists(marker_file)) return(data.table(cluster = character(), gene = character()))
  dt <- fread(marker_file, sep = "\t")
  if (nrow(dt) == 0 || !"cluster" %in% colnames(dt)) return(data.table(cluster = character(), gene = character()))
  fc_col <- intersect(c("avg_log2FC", "avg_logFC", "log2FC", "avg_diff"), names(dt))[1]
  gene_col <- if ("gene" %in% names(dt)) "gene" else names(dt)[1]
  if (is.na(fc_col)) return(data.table(cluster = character(), gene = character()))
  dt <- dt[order(cluster, -get(fc_col))]
  out <- dt[, head(.SD, topn), by = cluster]
  setnames(out, gene_col, "gene")
  out[, .(cluster, gene)]
}

panel_rows <- list()
plot_marker_panel <- function(assay, level, group_col, cells = NULL) {
  marker_file <- file.path(markerDir, paste0(assay, "_", level, "_markers_top20.tsv"))
  panel <- get_top_markers(marker_file, topn = 3)
  if (nrow(panel) == 0 || !assay %in% Assays(obj)) return(invisible(NULL))

  x <- obj
  if (!is.null(cells)) x <- subset(x, cells = cells)
  genes <- get_existing_genes(x, unique(panel$gene), assay)
  if (length(genes) == 0) return(invisible(NULL))
  DefaultAssay(x) <- assay
  prefix <- paste0(assay, "_", level)

  p_dot <- DotPlot(x, features = genes, group.by = group_col, assay = assay) +
    RotatedAxis() + ggtitle(paste(assay, "top markers by", gsub("_", " ", level)))
  ggsave(file.path(plotDir, paste0(prefix, "_DotPlot.pdf")), p_dot, width = 18, height = 8, limitsize = FALSE)

  if ("wnn.umap" %in% Reductions(x)) {
    pdf(file.path(plotDir, paste0(prefix, "_FeaturePlots.pdf")), width = 12, height = 8)
    for (g in genes) print(FeaturePlot(x, features = g, reduction = "wnn.umap", order = TRUE) + ggtitle(paste0(assay, ": ", g)))
    dev.off()
  }
  if (length(genes) > 1) {
    pdf(file.path(plotDir, paste0(prefix, "_Heatmap.pdf")), width = 12, height = 10)
    print(DoHeatmap(x, features = genes, group.by = group_col, assay = assay) +
            ggtitle(paste(assay, "marker heatmap by", gsub("_", " ", level))))
    dev.off()
  }
  panel_rows[[length(panel_rows) + 1]] <<- cbind(assay = assay, annotation_level = level, panel)
}

major_cells <- rownames(obj@meta.data)[obj@meta.data$celltype_major_inferred != "Unresolved_lineage" & !is.na(obj@meta.data$celltype_major_inferred)]
epi_cells <- rownames(obj@meta.data)[obj@meta.data$celltype_major_inferred == "Epithelial_cells" & !is.na(obj@meta.data$celltype_minor)]
state_cells <- rownames(obj@meta.data)[obj@meta.data$cell_state %in% c("Stress_response", "Cycling")]

for (assay in c("RNA", "ACTIVITY")) {
  plot_marker_panel(assay, "major_celltype", "celltype_major_inferred", major_cells)
  plot_marker_panel(assay, "epithelial_minor", "celltype_minor", epi_cells)
  plot_marker_panel(assay, "cell_state", "cell_state", state_cells)
}

panel_dt <- if (length(panel_rows) > 0) rbindlist(panel_rows, fill = TRUE) else data.table(
  assay = character(), annotation_level = character(), cluster = character(), gene = character()
)
fwrite(panel_dt, file.path(plotDir, "marker_panel_gene_list.tsv"), sep = "\t")
message("Done. Marker panel plots saved to: ", plotDir)
