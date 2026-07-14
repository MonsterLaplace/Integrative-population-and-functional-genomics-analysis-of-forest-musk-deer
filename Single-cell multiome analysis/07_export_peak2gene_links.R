suppressPackageStartupMessages({
  library(ArchR)
  library(data.table)
})

proj <- readRDS("06.multiome/03.archr/ArchR_proj/ArchRProject.afterPeak2Gene.rds")

p2g <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = 0.45,
  resolution = 1
)

# p2g 是 GRanges / DataFrame 风格对象
dt <- as.data.table(as.data.frame(p2g))
fwrite(dt, "06.multiome/03.archr/ArchR_proj/peak2gene_links.tsv", sep = "\t")

message("Saved: 06.multiome/03.archr/ArchR_proj/peak2gene_links.tsv")
