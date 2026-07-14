#!/usr/bin/env Rscript
# Plot sample composition at major lineage, epithelial minor-subtype, and state levels.

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(scales)
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

plotDir <- file.path(outDir, "composition_plots")
dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)
meta <- as.data.table(obj@meta.data, keep.rownames = "cell")

make_composition <- function(dt, group_col, denominator_label) {
  out <- dt[, .N, by = c("sample", group_col)]
  setnames(out, "N", "n_cells")
  setnames(out, group_col, "annotation")
  out[, frac_in_sample := n_cells / sum(n_cells), by = sample]
  out[, frac_in_annotation := n_cells / sum(n_cells), by = annotation]
  out[, denominator := denominator_label]
  out
}

save_composition_plots <- function(comp, prefix, title_prefix, width = 11, height = 7) {
  p_frac <- ggplot(comp, aes(x = sample, y = frac_in_sample, fill = annotation)) +
    geom_col(width = 0.8) +
    scale_y_continuous(labels = percent_format()) +
    labs(x = "Sample", y = "Fraction of cells", fill = NULL,
         title = paste(title_prefix, "composition by sample")) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
  ggsave(file.path(plotDir, paste0(prefix, "_fraction_by_sample.pdf")), p_frac, width = width, height = height)

  p_count <- ggplot(comp, aes(x = sample, y = n_cells, fill = annotation)) +
    geom_col(width = 0.8) +
    labs(x = "Sample", y = "Number of cells", fill = NULL,
         title = paste(title_prefix, "counts by sample")) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
  ggsave(file.path(plotDir, paste0(prefix, "_counts_by_sample.pdf")), p_count, width = width, height = height)

  p_heat <- ggplot(comp, aes(x = sample, y = annotation, fill = frac_in_sample)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "firebrick", labels = percent_format()) +
    labs(x = "Sample", y = NULL, fill = "Fraction", title = paste(title_prefix, "fraction heatmap")) +
    theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(plotDir, paste0(prefix, "_fraction_heatmap.pdf")), p_heat, width = 9, height = height)
}

# Primary composition: use merged lineages. Stress/cycling-only clusters are
# excluded because they are cell states and do not have a reliable lineage call.
major_meta <- meta[celltype_major_inferred != "Unresolved_lineage" & !is.na(celltype_major_inferred)]
major_comp <- make_composition(major_meta, "celltype_major_inferred", "lineage_assigned_cells")
fwrite(major_comp, file.path(plotDir, "major_celltype_composition_by_sample.tsv"), sep = "\t")
save_composition_plots(major_comp, "01_major_celltype", "Major cell-type")

major_count_wide <- dcast(major_comp, annotation ~ sample, value.var = "n_cells", fill = 0)
major_frac_wide <- dcast(major_comp, annotation ~ sample, value.var = "frac_in_sample", fill = 0)
fwrite(major_count_wide, file.path(plotDir, "major_celltype_by_sample_counts.wide.tsv"), sep = "\t")
fwrite(major_frac_wide, file.path(plotDir, "major_celltype_by_sample_fractions.wide.tsv"), sep = "\t")

# Secondary composition: inspect only epithelial cells at the minor subtype level.
epi_meta <- meta[celltype_major_inferred == "Epithelial_cells" & !is.na(celltype_minor)]
if (nrow(epi_meta) > 0) {
  epi_comp <- make_composition(epi_meta, "celltype_minor", "epithelial_cells")
  fwrite(epi_comp, file.path(plotDir, "epithelial_minor_celltype_composition_by_sample.tsv"), sep = "\t")
  save_composition_plots(epi_comp, "02_epithelial_minor_celltype", "Epithelial minor cell-type")
}

# State composition is reported separately and never mixed with the lineage plot.
state_comp <- make_composition(meta, "cell_state", "all_cells")
fwrite(state_comp, file.path(plotDir, "cell_state_composition_by_sample.tsv"), sep = "\t")
save_composition_plots(state_comp, "03_cell_state", "Cell-state")

sample_counts <- meta[, .N, by = sample]
setnames(sample_counts, "N", "n_cells")
fwrite(sample_counts, file.path(plotDir, "sample_total_cell_counts.tsv"), sep = "\t")

excluded <- meta[celltype_major_inferred == "Unresolved_lineage", .N, by = .(sample, cell_state)]
setnames(excluded, "N", "n_cells")
fwrite(excluded, file.path(plotDir, "unresolved_lineage_state_counts_by_sample.tsv"), sep = "\t")

message("Done. Composition plots saved to: ", plotDir)
