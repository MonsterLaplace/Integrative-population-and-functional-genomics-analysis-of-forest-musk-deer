#!/usr/bin/env Rscript
# Plot UMAPs, marker panels, and sample composition from the state-lineage
# inferred Seurat object without rebuilding cluster annotations.

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(data.table)
  library(ggplot2)
  library(scales)
})

outDir <- "06.multiome/04.seurat_signac"
objFile <- file.path(outDir, "FMdeer_multiome.wnn.final_annotated.state_lineage_inferred.rds")
plotDir <- file.path(outDir, "13b_featureplot_dotplot_composition_inferred")
if (!file.exists(objFile)) stop("Inferred object not found: ", objFile)
dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(objFile)
required_metadata <- c("sample", "celltype_detailed", "celltype_minor", "celltype_major_inferred", "cell_state")
if (!all(required_metadata %in% colnames(obj@meta.data))) {
  stop("Object is missing: ", paste(setdiff(required_metadata, colnames(obj@meta.data)), collapse = ", "))
}
if (!"wnn.umap" %in% Reductions(obj)) stop("wnn.umap reduction not found.")

save_plot <- function(p, filename, width = 12, height = 8) {
  ggsave(file.path(plotDir, filename), p, width = width, height = height, limitsize = FALSE)
}

meta <- as.data.table(obj@meta.data, keep.rownames = "cell")
fwrite(meta, file.path(plotDir, "cell_metadata_with_inferred_lineages.tsv"), sep = "\t")

# UMAPs at the three complementary annotation layers.
save_plot(
  DimPlot(obj, reduction = "wnn.umap", group.by = "celltype_major_inferred", label = TRUE, repel = TRUE) +
    ggtitle("WNN UMAP: inferred major cell types"),
  "UMAP_major_celltypes_inferred.pdf", 14, 9
)
save_plot(
  DimPlot(obj, reduction = "wnn.umap", group.by = "celltype_minor", label = TRUE, repel = TRUE) +
    ggtitle("WNN UMAP: minor cell types"),
  "UMAP_minor_celltypes.pdf", 14, 9
)
save_plot(
  DimPlot(obj, reduction = "wnn.umap", group.by = "cell_state", label = TRUE, repel = TRUE) +
    ggtitle("WNN UMAP: stress and cycling states"),
  "UMAP_cell_states.pdf", 12, 8
)

# Broad marker panel; features absent from the deer RNA annotation are skipped.
if ("RNA" %in% Assays(obj)) {
  DefaultAssay(obj) <- "RNA"
  markers <- c("KRT14", "KRT79", "KRT77", "KRT5", "DMKN", "RBP2", "KLK5", "SPINK5",
               "BPIFB1", "CSF1R", "C1QB", "CD3D", "IL7R", "POU2AF1", "COL1A1", "DCN",
               "RGS5", "VWF", "FLT4", "SOX10", "ACTG2", "ATF3", "HJURP", "DEPDC1B")
  markers <- intersect(markers, rownames(obj[["RNA"]]))
  if (length(markers) > 0) {
    save_plot(FeaturePlot(obj, features = markers, reduction = "wnn.umap", ncol = 4),
              "RNA_FeaturePlot_markers.pdf", 14, max(6, ceiling(length(markers) / 4) * 3.2))
    p_major <- DotPlot(obj, features = markers, group.by = "celltype_major_inferred") +
      RotatedAxis() + ggtitle("RNA markers by inferred major cell type")
    save_plot(p_major, "RNA_DotPlot_major_celltypes_inferred.pdf", 14, 8)

    epithelial_cells <- WhichCells(obj, expression = celltype_major_inferred == "Epithelial_cells")
    if (length(epithelial_cells) > 0) {
      p_epi <- DotPlot(subset(obj, cells = epithelial_cells), features = markers, group.by = "celltype_minor") +
        RotatedAxis() + ggtitle("RNA markers across epithelial minor cell types")
      save_plot(p_epi, "RNA_DotPlot_epithelial_minor_celltypes.pdf", 14, 8)
    }
  }
}

# Major-lineage composition includes accepted inferred state cells in their
# inferred source lineage. Low-confidence cells remain unresolved and are shown
# only in the state composition table.
major_comp <- meta[!is.na(celltype_major_inferred) & celltype_major_inferred != "Unresolved_lineage",
                   .(n_cells = .N), by = .(sample, celltype_major_inferred)]
major_comp[, fraction := n_cells / sum(n_cells), by = sample]
fwrite(major_comp, file.path(plotDir, "major_celltype_inferred_composition_by_sample.tsv"), sep = "\t")
save_plot(
  ggplot(major_comp, aes(x = sample, y = fraction, fill = celltype_major_inferred)) +
    geom_col(width = 0.8) + scale_y_continuous(labels = percent_format()) +
    labs(x = "Sample", y = "Fraction of lineage-assigned cells", fill = "Major cell type",
         title = "Inferred major cell-type composition by sample") +
    theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 30, hjust = 1)),
  "major_celltype_inferred_composition_by_sample.pdf", 12, 7
)

state_comp <- meta[, .(n_cells = .N), by = .(sample, cell_state)]
state_comp[, fraction := n_cells / sum(n_cells), by = sample]
fwrite(state_comp, file.path(plotDir, "cell_state_composition_by_sample.tsv"), sep = "\t")

message("Done. Plots and tables saved to: ", normalizePath(plotDir))
