suppressPackageStartupMessages({
  library(ArchR)
})

addArchRThreads(threads = 16)

geneAnnotation <- readRDS("06.multiome/03.archr/FMdeer_geneAnnotation.rds")
genomeAnnotation <- readRDS("06.multiome/03.archr/FMdeer_genomeAnnotation.rds")

inputFiles <- c(
  FMD3      = "06.multiome/cellranger_arc/FMD3/outs/atac_fragments.tsv.gz",
  FMD4      = "06.multiome/cellranger_arc/FMD4/outs/atac_fragments.tsv.gz",
  J21010063 = "06.multiome/cellranger_arc/J21010063/outs/atac_fragments.tsv.gz",
  J21090141 = "06.multiome/cellranger_arc/J21090141/outs/atac_fragments.tsv.gz"
)

missingFiles <- inputFiles[!file.exists(inputFiles)]
if (length(missingFiles) > 0) {
  stop("Missing fragment files:\n", paste0(names(missingFiles), " : ", missingFiles, collapse = "\n"))
}

message("Input fragment files:")
print(inputFiles)

arrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  minTSS = 0,
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

print(arrowFiles)

proj <- ArchRProject(
  ArrowFiles = arrowFiles,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  outputDirectory = "06.multiome/03.archr/ArchR_proj",
  copyArrows = TRUE
)

saveRDS(proj, file = "06.multiome/03.archr/ArchR_proj/ArchRProject.rds")
print(proj)
