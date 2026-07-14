suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome)
  library(BSgenome.FMdeer.Custom.v1)
})

addArchRThreads(threads = 16)

proj <- readRDS("06.multiome/03.archr/ArchR_proj/ArchRProject.afterGeneScore.rds")

# --------------------------------------------------
# 1. Browser track examples
# geneSymbol 需要替换成你的注释中真实存在的基因名
# --------------------------------------------------
candidateGenes <- c("CD3D", "MS4A1", "LYZ")

for (g in candidateGenes) {
  try({
    p <- plotBrowserTrack(
      ArchRProj = proj,
      groupBy = "Clusters",
      geneSymbol = g
    )

    plotPDF(
      p,
      name = paste0("06_browser_", g, ".pdf"),
      ArchRProj = proj
    )

    message("Browser track plotted for: ", g)
  }, silent = TRUE)
}

# --------------------------------------------------
# 2. Check available matrices
# --------------------------------------------------
availableMatrices <- getAvailableMatrices(proj)
message("Available matrices:")
print(availableMatrices)

# --------------------------------------------------
# 3. Peak2Gene only if expression/integration matrix exists
# --------------------------------------------------
if ("GeneIntegrationMatrix" %in% availableMatrices) {
  message("Using GeneIntegrationMatrix for peak2gene links...")
  proj <- addPeak2GeneLinks(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    useMatrix = "GeneIntegrationMatrix"
  )
  saveRDS(proj, file = "06.multiome/03.archr/ArchR_proj/ArchRProject.afterPeak2Gene.rds")
  message("Saved project with peak2gene links.")
} else if ("GeneExpressionMatrix" %in% availableMatrices) {
  message("Using GeneExpressionMatrix for peak2gene links...")
  proj <- addPeak2GeneLinks(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    useMatrix = "GeneExpressionMatrix"
  )
  saveRDS(proj, file = "06.multiome/03.archr/ArchR_proj/ArchRProject.afterPeak2Gene.rds")
  message("Saved project with peak2gene links.")
} else {
  message("No GeneIntegrationMatrix or GeneExpressionMatrix found.")
  message("Skipping addPeak2GeneLinks().")
}
