suppressPackageStartupMessages({
  library(ArchR)
})

addArchRThreads(threads = 16)

proj <- readRDS("06.multiome/03.archr/ArchR_proj/ArchRProject.afterPeaks.rds")

# Marker genes from GeneScoreMatrix
markersGS <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

saveRDS(markersGS, "06.multiome/03.archr/ArchR_proj/markersGeneScore.rds")

markerGenes <- getMarkers(
  markersGS,
  cutOff = "FDR <= 0.01 & Log2FC >= 1"
)

saveRDS(markerGenes, "06.multiome/03.archr/ArchR_proj/markerGeneList.rds")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  transpose = TRUE
)

plotPDF(
  heatmapGS,
  name = "04_markerGeneScore_heatmap.pdf",
  ArchRProj = proj
)

saveRDS(proj, file = "06.multiome/03.archr/ArchR_proj/ArchRProject.afterGeneScore.rds")
