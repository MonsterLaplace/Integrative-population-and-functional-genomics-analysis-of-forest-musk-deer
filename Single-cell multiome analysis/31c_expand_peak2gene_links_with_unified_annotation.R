#!/usr/bin/env Rscript
# Expand global LinkPeaks evidence to gene-name and gene-ID identifier variants.

suppressPackageStartupMessages({
  library(data.table)
})

outDir <- "06.multiome/04.seurat_signac"
realP2GFile <- file.path(outDir, "peak2gene_links.tsv")
annoFile <- file.path(outDir, "FMdeer_unified_gene_annotation.unique.tsv")
if (!file.exists(realP2GFile)) stop("Missing file: ", realP2GFile)
if (!file.exists(annoFile)) stop("Missing file: ", annoFile)

p2g <- fread(realP2GFile)
anno <- fread(annoFile)
if (!"gene" %in% names(p2g)) stop("peak2gene_links.tsv missing gene column")
if (!all(c("gene_id", "gene_id_dash", "gene_name") %in% names(anno))) {
  stop("Unified annotation must contain gene_id, gene_id_dash, gene_name")
}

# LinkPeaks evidence was computed across all cells. Preserve that provenance
# during ID expansion rather than assigning an artificial cell-type specificity.
if (!"link_scope" %in% names(p2g)) p2g[, link_scope := "all_cells"]
if (!"annotation_level" %in% names(p2g)) p2g[, annotation_level := "global"]
for (cc in c("peak2gene_support", "peak2gene_n_links", "peak2gene_max_score")) {
  if (!cc %in% names(p2g)) p2g[, (cc) := NA_real_]
}

anno2 <- unique(anno[, .(gene_id, gene_id_dash, gene_name)])
p2g_expand <- merge(p2g, anno2, by.x = "gene", by.y = "gene_name", all.x = TRUE, sort = FALSE)

make_keys <- function(dt, key_col, match_type) {
  if (!key_col %in% names(dt)) return(data.table())
  x <- dt[!is.na(get(key_col)) & get(key_col) != "", .(
    match_key = as.character(get(key_col)),
    match_type = match_type,
    peak2gene_support, peak2gene_n_links, peak2gene_max_score,
    link_scope, annotation_level
  )]
  unique(x)
}

expanded <- rbindlist(list(
  make_keys(p2g_expand, "gene", "gene_name"),
  make_keys(p2g_expand, "gene_id", "gene_id"),
  make_keys(p2g_expand, "gene_id_dash", "gene_id_dash")
), fill = TRUE)

# If a key has several forms, keep the most strongly supported link while
# retaining the all-cell/global provenance columns.
setorder(expanded, match_key, -peak2gene_max_score, na.last = TRUE)
expanded_best <- expanded[, .SD[1], by = match_key]
fwrite(expanded_best, file.path(outDir, "peak2gene_links.expanded.tsv"), sep = "\t", quote = FALSE)
message("Saved expanded peak2gene file: ", file.path(outDir, "peak2gene_links.expanded.tsv"))
