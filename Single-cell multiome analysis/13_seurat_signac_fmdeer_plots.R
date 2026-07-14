suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(patchwork)
})

outDir <- "06.multiome/04.seurat_signac"
obj <- readRDS(file.path(outDir, "FMdeer_multiome.wnn.robust.markers.rds"))

pdf(file.path(outDir, "05_featureplots_examples.pdf"), width = 12, height = 10)

# 这里换成你实际存在的基因名
genes_to_plot <- c("CD3D", "MS4A1", "LYZ")

DefaultAssay(obj) <- "RNA"
for (g in genes_to_plot) {
  if (g %in% rownames(obj[["RNA"]])) {
    print(FeaturePlot(obj, features = g, reduction = "wnn.umap"))
  }
}

DefaultAssay(obj) <- "ACTIVITY"
for (g in genes_to_plot) {
  if (g %in% rownames(obj[["ACTIVITY"]])) {
    print(FeaturePlot(obj, features = g, reduction = "wnn.umap"))
  }
}

dev.off()
