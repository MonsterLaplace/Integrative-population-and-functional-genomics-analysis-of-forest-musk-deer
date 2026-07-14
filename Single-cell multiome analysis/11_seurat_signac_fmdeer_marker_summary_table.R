suppressPackageStartupMessages({
  library(data.table)
})

outDir <- "06.multiome/04.seurat_signac"
markerDir <- file.path(outDir, "markers")

make_wide_top_table <- function(file, outfile, topn = 10, gene_col = "gene") {
  dt <- fread(file)
  if (!all(c("cluster", gene_col) %in% colnames(dt))) {
    stop("Missing required columns in: ", file)
  }

  dt_top <- dt[, head(.SD, topn), by = cluster]

  dt_top[, rank := sequence(.N), by = cluster]
  wide <- dcast(dt_top, rank ~ cluster, value.var = gene_col)

  fwrite(wide, outfile, sep = "\t")
}

rna_file <- file.path(markerDir, "RNA_markers_top20.tsv")
if (file.exists(rna_file)) {
  make_wide_top_table(
    file = rna_file,
    outfile = file.path(markerDir, "RNA_markers_top10_wide.tsv"),
    topn = 10,
    gene_col = "gene"
  )
}

activity_file <- file.path(markerDir, "ACTIVITY_markers_top20.tsv")
if (file.exists(activity_file)) {
  make_wide_top_table(
    file = activity_file,
    outfile = file.path(markerDir, "ACTIVITY_markers_top10_wide.tsv"),
    topn = 10,
    gene_col = "gene"
  )
}

message("Done.")
