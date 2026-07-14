suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(data.table)
  library(dplyr)
})

outDir <- "06.multiome/04.seurat_signac"
obj <- readRDS(file.path(outDir, "FMdeer_multiome.wnn.robust.rds"))

markerDir <- file.path(outDir, "markers")
dir.create(markerDir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------
# helper: safely join layers for Seurat v5
# --------------------------------------------
safe_join_layers <- function(obj, assay_name) {
  if (!assay_name %in% Assays(obj)) {
    message("Assay ", assay_name, " not found, skipping JoinLayers.")
    return(obj)
  }

  message("Checking layers for assay: ", assay_name)
  DefaultAssay(obj) <- assay_name

  tryCatch({
    obj[[assay_name]] <- JoinLayers(obj[[assay_name]])
    message("JoinLayers done for assay: ", assay_name)
  }, error = function(e) {
    message("JoinLayers skipped/failed for assay ", assay_name, ": ", e$message)
  })

  obj
}

# --------------------------------------------
# helper: ensure assay data layer exists
# --------------------------------------------
safe_prepare_assay <- function(obj, assay_name) {
  if (!assay_name %in% Assays(obj)) {
    message("Assay ", assay_name, " not found, skipping preparation.")
    return(obj)
  }

  DefaultAssay(obj) <- assay_name
  obj <- safe_join_layers(obj, assay_name)

  need_norm <- FALSE
  tryCatch({
    x <- GetAssayData(obj, assay = assay_name, layer = "data")
    if (is.null(x) || nrow(x) == 0 || ncol(x) == 0) need_norm <- TRUE
  }, error = function(e) {
    need_norm <<- TRUE
  })

  if (need_norm) {
    message(assay_name, " data layer missing/empty, running NormalizeData() ...")
    obj <- tryCatch({
      NormalizeData(obj, assay = assay_name, verbose = FALSE)
    }, error = function(e) {
      message("NormalizeData failed for assay ", assay_name, ": ", e$message)
      obj
    })
  }

  obj
}

# --------------------------------------------
# helper: write markers robustly
# --------------------------------------------
write_top_markers <- function(markers, prefix, outdir, topn = 20) {
  if (is.null(markers) || nrow(markers) == 0) {
    message("No markers found for ", prefix, ". Writing empty file.")
    fwrite(data.frame(), file.path(outdir, paste0(prefix, "_all.tsv")), sep = "\t")
    return(invisible(NULL))
  }

  markers <- as.data.frame(markers)
  fwrite(markers, file.path(outdir, paste0(prefix, "_all.tsv")), sep = "\t")

  fc_col <- NULL
  for (x in c("avg_log2FC", "avg_logFC", "log2FC", "avg_diff")) {
    if (x %in% colnames(markers)) {
      fc_col <- x
      break
    }
  }

  if (is.null(fc_col) || !"cluster" %in% colnames(markers)) {
    message("Cannot find cluster/logFC columns for ", prefix, ", skipping top tables.")
    return(invisible(NULL))
  }

  top_dt <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = .data[[fc_col]], n = topn, with_ties = FALSE) %>%
    ungroup()

  fwrite(as.data.table(top_dt), file.path(outdir, paste0(prefix, "_top", topn, ".tsv")), sep = "\t")
}

# --------------------------------------------
# helper: run FindAllMarkers robustly
# --------------------------------------------
run_markers <- function(obj, assay_name, prefix, outdir, only.pos = TRUE,
                        min.pct = 0.1, logfc.threshold = 0.25, topn = 20) {

  if (!assay_name %in% Assays(obj)) {
    message("Assay ", assay_name, " not found, skipping ", prefix)
    return(obj)
  }

  DefaultAssay(obj) <- assay_name

  if (assay_name %in% c("RNA", "ACTIVITY")) {
    obj <- safe_prepare_assay(obj, assay_name)
  }

  message("Running FindAllMarkers for assay: ", assay_name)

  markers <- tryCatch({
    FindAllMarkers(
      object = obj,
      only.pos = only.pos,
      min.pct = min.pct,
      logfc.threshold = logfc.threshold,
      test.use = "wilcox",
      verbose = TRUE
    )
  }, error = function(e) {
    message("FindAllMarkers failed for ", assay_name, ": ", e$message)
    data.frame()
  })

  write_top_markers(markers, prefix, outdir, topn = topn)
  obj
}

# --------------------------------------------
# ensure cluster identities
# --------------------------------------------
if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
  stop("Cannot find seurat_clusters in metadata.")
}
Idents(obj) <- "seurat_clusters"

# --------------------------------------------
# 1. RNA markers
# --------------------------------------------
obj <- run_markers(
  obj = obj,
  assay_name = "RNA",
  prefix = "RNA_markers",
  outdir = markerDir,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  topn = 20
)

# --------------------------------------------
# 2. ACTIVITY markers
# --------------------------------------------
if ("ACTIVITY" %in% Assays(obj)) {
  obj <- run_markers(
    obj = obj,
    assay_name = "ACTIVITY",
    prefix = "ACTIVITY_markers",
    outdir = markerDir,
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.25,
    topn = 20
  )
} else {
  message("ACTIVITY assay not found, skipping.")
}

# --------------------------------------------
# 3. ATAC peak markers
# --------------------------------------------
if ("ATAC" %in% Assays(obj)) {
  obj <- run_markers(
    obj = obj,
    assay_name = "ATAC",
    prefix = "ATAC_peak_markers",
    outdir = markerDir,
    only.pos = TRUE,
    min.pct = 0.05,
    logfc.threshold = 0.25,
    topn = 20
  )
} else {
  message("ATAC assay not found, skipping.")
}

# --------------------------------------------
# 4. cluster summary report
# --------------------------------------------
if ("sample" %in% colnames(obj@meta.data)) {
  cluster_summary <- as.data.table(table(obj$seurat_clusters, obj$sample))
  colnames(cluster_summary) <- c("cluster", "sample", "n_cells")
  fwrite(cluster_summary, file.path(markerDir, "cluster_sample_distribution.tsv"), sep = "\t")
} else {
  message("sample column not found in metadata, skipping cluster_sample_distribution.tsv")
}

saveRDS(obj, file.path(outDir, "FMdeer_multiome.wnn.robust.markers.rds"))
message("Done.")
