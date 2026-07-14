suppressPackageStartupMessages({
  library(ArchR)
})

addArchRThreads(threads = 16)

proj <- readRDS("06.multiome/03.archr/ArchR_proj/ArchRProject.rds")

# 1. Doublet score
proj <- addDoubletScores(
  input = proj,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1
)

# 2. Filter doublets
proj <- filterDoublets(proj)

# 3. Dimensional reduction
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list(
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30
)

# 4. Clustering
proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)

# 5. UMAP
proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine"
)

# 6. Plot
p1 <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP"
)

p2 <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Clusters",
  embedding = "UMAP"
)

plotPDF(
  p1, p2,
  name = "UMAP_by_sample_cluster.pdf",
  ArchRProj = proj
)

saveRDS(proj, file = "06.multiome/03.archr/ArchR_proj/ArchRProject.afterQC.rds")
