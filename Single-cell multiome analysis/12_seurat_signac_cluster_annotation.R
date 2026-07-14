suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

outDir <- "06.multiome/04.seurat_signac"
obj <- readRDS(file.path(outDir, "FMdeer_multiome.wnn.robust.markers.rds"))

if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
  stop("seurat_clusters not found in metadata.")
}

Idents(obj) <- "seurat_clusters"

# ---------------------------------------------------------
# cluster ids
# ---------------------------------------------------------
cluster_ids <- levels(Idents(obj))

# 默认占位注释，后面你可手工改
annotation_map <- setNames(
  paste0("Cluster_", cluster_ids),
  cluster_ids
)

# 示例：
# annotation_map["0"] <- "T_cell"
# annotation_map["1"] <- "B_cell"
# annotation_map["2"] <- "Myeloid"

# ---------------------------------------------------------
# 生成按细胞顺序排列的 celltype 向量
# 注意去掉 names，避免 Seurat 按名字误对齐
# ---------------------------------------------------------
celltype_vec <- annotation_map[as.character(Idents(obj))]
celltype_vec <- unname(as.character(celltype_vec))

# 直接写入 meta.data
obj@meta.data$celltype <- celltype_vec
obj@meta.data$celltype <- factor(obj@meta.data$celltype, levels = unique(unname(annotation_map)))

# =========================================================
# 导出 cell metadata
# =========================================================
meta <- obj@meta.data
meta$cell <- rownames(meta)
fwrite(as.data.table(meta), file.path(outDir, "cell_metadata_with_celltype.tsv"), sep = "\t")

# =========================================================
# cluster -> celltype 对照表
# =========================================================
anno_dt <- data.table(
  cluster = names(annotation_map),
  celltype = unname(annotation_map)
)
fwrite(anno_dt, file.path(outDir, "cluster_to_celltype.tsv"), sep = "\t")

# =========================================================
# 每个 sample 的 celltype 统计
# =========================================================
if ("sample" %in% colnames(obj@meta.data)) {
  ct_sample <- as.data.table(table(obj@meta.data$celltype, obj@meta.data$sample))
  colnames(ct_sample) <- c("celltype", "sample", "n_cells")
  fwrite(ct_sample, file.path(outDir, "celltype_by_sample_counts.tsv"), sep = "\t")
}

# =========================================================
# 画图
# =========================================================
pdf(file.path(outDir, "06_annotated_UMAPs.pdf"), width = 14, height = 10)

print(DimPlot(
  obj,
  reduction = "wnn.umap",
  group.by = "seurat_clusters",
  label = TRUE
) + ggtitle("WNN UMAP by cluster"))

print(DimPlot(
  obj,
  reduction = "wnn.umap",
  group.by = "celltype",
  label = TRUE
) + ggtitle("WNN UMAP by celltype"))

if ("sample" %in% colnames(obj@meta.data)) {
  print(DimPlot(
    obj,
    reduction = "wnn.umap",
    group.by = "sample"
  ) + ggtitle("WNN UMAP by sample"))
}

dev.off()

saveRDS(obj, file.path(outDir, "FMdeer_multiome.wnn.annotated.rds"))
message("Done.")
