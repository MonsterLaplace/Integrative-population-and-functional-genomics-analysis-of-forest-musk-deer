suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(BSgenome)
  library(BSgenome.FMdeer.Custom.v1)
})

outDir <- "06.multiome/04.seurat_signac"
obj <- readRDS(file.path(outDir, "FMdeer_multiome.wnn.robust.markers.rds"))

DefaultAssay(obj) <- "ATAC"

# 如果此前没算过 RegionStats，可再跑一次
obj <- RegionStats(obj, genome = BSgenome.FMdeer.Custom.v1)

obj <- LinkPeaks(
  object = obj,
  peak.assay = "ATAC",
  expression.assay = "RNA"
)

saveRDS(obj, file.path(outDir, "FMdeer_multiome.wnn.linkpeaks.rds"))
message("Done.")
