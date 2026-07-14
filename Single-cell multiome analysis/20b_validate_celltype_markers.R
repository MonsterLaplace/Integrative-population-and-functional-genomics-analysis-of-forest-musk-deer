#!/usr/bin/env Rscript
# Validate detailed/minor/major/states annotations after applying script 19.

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

outDir <- "06.multiome/04.seurat_signac"
inferredFile <- file.path(outDir, "FMdeer_multiome.wnn.final_annotated.state_lineage_inferred.rds")
objFile <- if (file.exists(inferredFile)) inferredFile else file.path(outDir, "FMdeer_multiome.wnn.final_annotated.rds")
helperFile <- file.path(outDir, "cluster_annotation_helper_table.tsv")
if (!file.exists(objFile)) stop("Annotated object not found: ", objFile)
if (!file.exists(helperFile)) stop("Helper table not found: ", helperFile)

obj <- readRDS(objFile)
helper <- fread(helperFile, sep = "\t", header = TRUE)
required_metadata <- c("celltype_detailed", "celltype_minor", "celltype_major", "celltype_major_inferred", "cell_state")
if (!all(required_metadata %in% colnames(obj@meta.data))) {
  stop("Annotation-level metadata missing. Run 19_apply_cluster_annotation_from_table.R first.")
}
if (!"RNA" %in% Assays(obj)) stop("RNA assay not found in the annotated object.")
if (!"rna_genes" %in% colnames(helper)) stop("rna_genes column not found in helper table.")

plotDir <- file.path(outDir, "validation_plots_helper")
dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)
Idents(obj) <- "celltype_major_inferred"

split_marker_string <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  y <- trimws(unlist(strsplit(x, ",\\s*")))
  y[y != ""]
}
get_existing_genes <- function(obj, genes, assay = "RNA") {
  if (!assay %in% Assays(obj)) return(character(0))
  genes[genes %in% rownames(obj[[assay]])]
}
save_dotplot <- function(plot, filename, width = 20, height = 9) {
  ggsave(file.path(plotDir, filename), plot, width = width, height = height, limitsize = FALSE)
}

# Use the top helper-table markers as a data-driven validation panel.
helper[, rna_gene_vec := lapply(rna_genes, split_marker_string)]
cluster2genes <- lapply(helper$rna_gene_vec, function(x) unique(head(x, 5)))
all_rna_genes <- unique(unlist(cluster2genes))
all_rna_genes_exist <- get_existing_genes(obj, all_rna_genes)
fwrite(data.table(gene = all_rna_genes, exists_in_RNA = all_rna_genes %in% rownames(obj[["RNA"]])),
       file.path(plotDir, "RNA_marker_existence_check.tsv"), sep = "\t")

if (length(all_rna_genes_exist) > 0) {
  DefaultAssay(obj) <- "RNA"
  p_major <- DotPlot(obj, features = all_rna_genes_exist, group.by = "celltype_major_inferred", assay = "RNA") +
    RotatedAxis() + ggtitle("RNA markers by merged major cell type")
  save_dotplot(p_major, "01_RNA_DotPlot_major_celltypes.pdf")

  # Broad Epithelial_cells are deliberately split here only for subtype validation.
  epithelial_cells <- WhichCells(obj, expression = celltype_major_inferred == "Epithelial_cells")
  if (length(epithelial_cells) > 0) {
    p_epi <- DotPlot(subset(obj, cells = epithelial_cells), features = all_rna_genes_exist,
                     group.by = "celltype_minor", assay = "RNA") +
      RotatedAxis() + ggtitle("RNA markers across epithelial minor cell types")
    save_dotplot(p_epi, "02_RNA_DotPlot_epithelial_minor_celltypes.pdf")
  }
}

# Feature plots retain cluster-derived markers so their spatial patterns remain visible.
feature_genes <- unique(unlist(lapply(helper$rna_gene_vec, function(x) {
  head(get_existing_genes(obj, x, assay = "RNA"), 2)
})))
if (length(feature_genes) > 0 && "wnn.umap" %in% Reductions(obj)) {
  DefaultAssay(obj) <- "RNA"
  pdf(file.path(plotDir, "03_RNA_FeaturePlot_helper_markers.pdf"), width = 12, height = 8)
  for (g in feature_genes) {
    print(FeaturePlot(obj, features = g, reduction = "wnn.umap", order = TRUE) + ggtitle(paste0("RNA: ", g)))
  }
  dev.off()
}

vln_genes <- head(all_rna_genes_exist, 20)
if (length(vln_genes) > 0) {
  DefaultAssay(obj) <- "RNA"
  p_vln <- VlnPlot(obj, features = vln_genes, group.by = "celltype_major_inferred", assay = "RNA", pt.size = 0, ncol = 4)
  ggsave(file.path(plotDir, "04_RNA_VlnPlot_major_celltypes.pdf"), p_vln, width = 16, height = 10, limitsize = FALSE)
}

if ("ACTIVITY" %in% Assays(obj) && "activity_genes" %in% colnames(helper)) {
  helper[, activity_gene_vec := lapply(activity_genes, split_marker_string)]
  all_act_genes <- unique(unlist(helper$activity_gene_vec))
  all_act_genes_exist <- get_existing_genes(obj, all_act_genes, assay = "ACTIVITY")
  fwrite(data.table(gene = all_act_genes, exists_in_ACTIVITY = all_act_genes %in% rownames(obj[["ACTIVITY"]])),
         file.path(plotDir, "ACTIVITY_marker_existence_check.tsv"), sep = "\t")
  if (length(all_act_genes_exist) > 0) {
    DefaultAssay(obj) <- "ACTIVITY"
    p_act <- DotPlot(obj, features = all_act_genes_exist, group.by = "celltype_major_inferred", assay = "ACTIVITY") +
      RotatedAxis() + ggtitle("ACTIVITY markers by merged major cell type")
    save_dotplot(p_act, "05_ACTIVITY_DotPlot_major_celltypes.pdf")
  }
}

# Validate stress and cycling separately because they are states, not final lineages.
state_markers <- get_existing_genes(obj, c("ATF3", "HJURP", "DEPDC1B", "NEIL3", "IQGAP3"), assay = "RNA")
if (length(state_markers) > 0) {
  DefaultAssay(obj) <- "RNA"
  p_state <- DotPlot(obj, features = state_markers, group.by = "cell_state", assay = "RNA") +
    RotatedAxis() + ggtitle("RNA markers for stress and cycling states")
  save_dotplot(p_state, "06_RNA_DotPlot_cell_states.pdf", width = 10, height = 5)
}

for (cc in required_metadata) {
  counts <- as.data.table(table(obj@meta.data[[cc]]))
  setnames(counts, c(cc, "n_cells"))
  fwrite(counts, file.path(plotDir, paste0(cc, "_counts.validation.tsv")), sep = "\t")
}

message("Done. Outputs saved to: ", plotDir)
