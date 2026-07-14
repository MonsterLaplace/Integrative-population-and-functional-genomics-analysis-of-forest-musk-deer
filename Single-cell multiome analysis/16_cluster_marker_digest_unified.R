suppressPackageStartupMessages({
  library(data.table)
})

outDir <- "06.multiome/04.seurat_signac"
markerDir <- file.path(outDir, "markers")

make_digest <- function(infile, outfile, topn = 10) {
  if (!file.exists(infile)) {
    stop("File not found: ", infile)
  }

  dt <- fread(infile)

  # 优先显示统一标签
  label_col <- NULL
  for (cc in c("display_label", "display_name", "gene_name", "gene")) {
    if (cc %in% colnames(dt)) {
      label_col <- cc
      break
    }
  }
  if (is.null(label_col)) {
    stop("No suitable label column found in: ", infile)
  }

  # 功能描述列
  desc_col <- NULL
  for (cc in c("display_description", "description", "product", "eggNOG_Description")) {
    if (cc %in% colnames(dt)) {
      desc_col <- cc
      break
    }
  }

  if (!"cluster" %in% colnames(dt)) {
    stop("Missing 'cluster' column in: ", infile)
  }

  fc_col <- NULL
  for (cc in c("avg_log2FC", "avg_logFC")) {
    if (cc %in% colnames(dt)) {
      fc_col <- cc
      break
    }
  }
  if (is.null(fc_col)) {
    stop("Missing logFC column in: ", infile)
  }

  dt <- dt[order(cluster, -get(fc_col))]

  digest_dt <- dt[, .(
    top_markers = paste(head(get(label_col), topn), collapse = ", "),
    top_descriptions = if (!is.null(desc_col)) {
      paste(head(get(desc_col), topn), collapse = " | ")
    } else {
      NA_character_
    }
  ), by = cluster]

  fwrite(digest_dt, outfile, sep = "\t")
  message("Saved: ", outfile)
}

make_digest(
  infile = file.path(markerDir, "RNA_markers_top20.unified.tsv"),
  outfile = file.path(markerDir, "RNA_cluster_marker_digest.unified.tsv"),
  topn = 10
)

make_digest(
  infile = file.path(markerDir, "ACTIVITY_markers_top20.unified.tsv"),
  outfile = file.path(markerDir, "ACTIVITY_cluster_marker_digest.unified.tsv"),
  topn = 10
)

message("Done.")
