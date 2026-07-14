suppressPackageStartupMessages({
  library(ArchR)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(IRanges)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# =========================================================
# Standalone marker peak heatmap
# ---------------------------------------------------------
# This script only reads the marker peak result saved by
# 03_peak_calling_marker_analysis.R and makes a publication-style heatmap.
# It does NOT re-run group coverage, peak calling, PeakMatrix generation, or
# marker testing, so it is much faster and safer to re-run for figure tuning.
# =========================================================

addArchRThreads(threads = 16)

archr_proj_dir <- "06.multiome/03.archr/ArchR_proj"
markers_file <- file.path(archr_proj_dir, "markersPeaks.rds")
archr_project_file <- file.path(archr_proj_dir, "ArchRProject.afterPeaks.rds")
plot_dir <- file.path(archr_proj_dir, "Plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(plot_dir, "03_markerPeaks_heatmap_custom_publication.pdf")
out_png <- file.path(plot_dir, "03_markerPeaks_heatmap_custom_publication.png")
out_label_table <- file.path(plot_dir, "03_markerPeaks_heatmap_custom_publication_peak_labels.tsv")

# Tune these if needed.
fdr_cutoff <- 0.01
log2fc_cutoff <- 1
top_n_per_cluster <- 4L
max_total_peaks <- 140L
show_peak_labels_each_side <- 4L
show_peak_labels_total <- 22L

if (!file.exists(markers_file)) {
  stop(
    "Cannot find marker peak result: ", markers_file, "\n",
    "Please run 03_peak_calling_marker_analysis.R first. It should save markersPeaks.rds before plotting."
  )
}

markersPeaks <- readRDS(markers_file)

pick_assay <- function(se, candidates) {
  available <- SummarizedExperiment::assayNames(se)
  hit <- candidates[candidates %in% available]
  if (length(hit) > 0) return(hit[[1]])
  available[[1]]
}

log2fc_assay <- pick_assay(markersPeaks, c("Log2FC", "log2FC", "Log2FCMatrix"))
fdr_assay <- pick_assay(markersPeaks, c("FDR", "fdr", "padj", "qval"))
signal_assay <- log2fc_assay

log2fc <- as.matrix(SummarizedExperiment::assay(markersPeaks, log2fc_assay))
fdr <- as.matrix(SummarizedExperiment::assay(markersPeaks, fdr_assay))
signal <- as.matrix(SummarizedExperiment::assay(markersPeaks, signal_assay))

if (nrow(signal) == 0 || ncol(signal) == 0) {
  stop("The marker peak assay is empty: ", markers_file)
}

cluster_names <- colnames(signal)
if (is.null(cluster_names)) cluster_names <- paste0("C", seq_len(ncol(signal)))
cluster_names_clean <- sub("^C(?=\\d+$)", "", cluster_names, perl = TRUE)

format_peak_midpoint_label <- function(seqname, start, end) {
  mid <- (as.numeric(start) + as.numeric(end)) / 2
  ifelse(
    mid >= 1e6,
    sprintf("%s:%.2f Mb", seqname, mid / 1e6),
    ifelse(
      mid >= 1e3,
      sprintf("%s:%.1f kb", seqname, mid / 1e3),
      sprintf("%s:%d bp", seqname, round(mid))
    )
  )
}

format_peak_full_coordinate <- function(seqname, start, end) {
  paste0(seqname, ":", as.character(start), "-", as.character(end))
}

is_valid_peak_coordinate <- function(seqname, start, end) {
  !is.na(seqname) &&
    grepl("^(chr|Chr|CHR)[A-Za-z0-9._-]+$", as.character(seqname)) &&
    is.finite(as.numeric(start)) &&
    is.finite(as.numeric(end)) &&
    as.numeric(end) >= as.numeric(start)
}

get_peak_coordinates_from_project <- function(project_file, peak_ids, n_signal_rows) {
  if (!file.exists(project_file)) return(NULL)
  message("Trying to recover genomic coordinates from ArchR project peak set: ", project_file)
  proj <- tryCatch(readRDS(project_file), error = function(e) {
    message("Could not read ArchR project for peak-coordinate recovery: ", e$message)
    NULL
  })
  if (is.null(proj)) return(NULL)

  peak_set <- tryCatch(ArchR::getPeakSet(proj), error = function(e) {
    message("Could not read peak set from ArchR project: ", e$message)
    NULL
  })
  if (is.null(peak_set) || length(peak_set) == 0) return(NULL)

  peak_ids_numeric <- suppressWarnings(as.integer(peak_ids))
  if (all(is.finite(peak_ids_numeric)) &&
      min(peak_ids_numeric, na.rm = TRUE) >= 1L &&
      max(peak_ids_numeric, na.rm = TRUE) <= length(peak_set)) {
    idx <- peak_ids_numeric
  } else if (length(peak_set) == n_signal_rows) {
    idx <- seq_len(n_signal_rows)
  } else {
    return(NULL)
  }

  data.frame(
    chr = as.character(GenomicRanges::seqnames(peak_set))[idx],
    start = IRanges::start(peak_set)[idx],
    end = IRanges::end(peak_set)[idx],
    stringsAsFactors = FALSE
  )
}

# Prefer genomic coordinates from rowRanges. In some ArchR objects, rownames are
# only numeric peak IDs, which are ambiguous in a publication figure.
rr <- tryCatch(SummarizedExperiment::rowRanges(markersPeaks), error = function(e) NULL)
peak_ids <- rownames(signal)
if (is.null(peak_ids)) peak_ids <- as.character(seq_len(nrow(signal)))

if (!is.null(rr) &&
    length(rr) == nrow(signal) &&
    any(mapply(
      is_valid_peak_coordinate,
      as.character(GenomicRanges::seqnames(rr)),
      IRanges::start(rr),
      IRanges::end(rr)
    ))) {
  peak_chr <- as.character(GenomicRanges::seqnames(rr))
  peak_start <- IRanges::start(rr)
  peak_end <- IRanges::end(rr)
  peak_names <- format_peak_midpoint_label(peak_chr, peak_start, peak_end)
  peak_full_coordinates <- format_peak_full_coordinate(peak_chr, peak_start, peak_end)
  peak_label_source <- rep("genomic_coordinate", length(peak_names))
} else {
  recovered <- get_peak_coordinates_from_project(
    project_file = archr_project_file,
    peak_ids = peak_ids,
    n_signal_rows = nrow(signal)
  )
  if (!is.null(recovered) &&
      nrow(recovered) == nrow(signal) &&
      any(mapply(is_valid_peak_coordinate, recovered$chr, recovered$start, recovered$end))) {
    peak_names <- format_peak_midpoint_label(recovered$chr, recovered$start, recovered$end)
    peak_full_coordinates <- format_peak_full_coordinate(recovered$chr, recovered$start, recovered$end)
    peak_label_source <- rep("archr_project_peakset", length(peak_names))
  } else {
    # Final fallback: keep the original peak IDs in the table only. The plot
    # itself will use compact P1/P2 labels to avoid ambiguous numeric labels.
    peak_names <- paste0("PeakID_", peak_ids)
    peak_full_coordinates <- rep(NA_character_, length(peak_names))
    peak_label_source <- rep("peak_id_fallback", length(peak_names))
    message(
      "Warning: genomic peak coordinates were not found in rowRanges(markersPeaks) ",
      "or in the ArchR project peak set. The plot will use compact P labels, ",
      "and peak IDs will be reported in the TSV table."
    )
  }
}

selected_by_cluster <- lapply(seq_along(cluster_names), function(j) {
  ok <- which(is.finite(fdr[, j]) & is.finite(log2fc[, j]) &
                fdr[, j] <= fdr_cutoff & log2fc[, j] >= log2fc_cutoff)
  if (length(ok) == 0) return(integer(0))
  ok <- ok[order(log2fc[ok, j], decreasing = TRUE)]
  head(ok, top_n_per_cluster)
})
names(selected_by_cluster) <- cluster_names

selected <- unique(unlist(selected_by_cluster, use.names = FALSE))
if (length(selected) == 0) {
  stop(
    "No marker peaks passed cutoffs: FDR <= ", fdr_cutoff,
    " and Log2FC >= ", log2fc_cutoff,
    ". Try lowering cutoffs at the top of this script."
  )
}

if (length(selected) > max_total_peaks) {
  selected <- selected[seq_len(max_total_peaks)]
}

owner <- rep(NA_character_, length(selected))
names(owner) <- selected
for (cl in names(selected_by_cluster)) {
  idx <- intersect(selected_by_cluster[[cl]], selected)
  owner[as.character(idx)] <- cl
}
owner[is.na(owner)] <- "Marker peak"

mat <- t(signal[selected, , drop = FALSE])
rownames(mat) <- cluster_names_clean
plot_peak_ids <- paste0("P", seq_along(selected))
rownames(mat) <- cluster_names_clean
colnames(mat) <- plot_peak_ids

label_table <- data.frame(
  column_index = seq_along(selected),
  plot_label = plot_peak_ids,
  marker_cluster = owner[as.character(selected)],
  peak_label = peak_names[selected],
  full_coordinate = peak_full_coordinates[selected],
  label_source = peak_label_source[selected],
  peak_source_index = selected,
  original_peak_id = peak_ids[selected],
  stringsAsFactors = FALSE
)
utils::write.table(
  label_table,
  file = out_label_table,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Column z-score, matching the usual ArchR marker heatmap interpretation while
# keeping the display robust and independent from plotMarkerHeatmap internals.
z <- t(scale(t(mat)))
z[!is.finite(z)] <- 0
z <- pmax(pmin(z, 2), -2)

make_peak_coordinate_labels <- function(labels,
                                        keep_each_side = 4L,
                                        max_labels = 22L,
                                        add_ellipsis = TRUE) {
  labels <- as.character(labels)
  n <- length(labels)
  if (n <= max_labels) return(labels)

  out <- rep("", n)

  left_idx <- seq_len(min(keep_each_side, n))
  right_idx <- seq(max(1, n - keep_each_side + 1), n)

  interior_capacity <- max(0L, max_labels - length(unique(c(left_idx, right_idx))))
  interior_pool <- setdiff(seq_len(n), c(left_idx, right_idx))
  interior_idx <- integer(0)
  if (interior_capacity > 0L && length(interior_pool) > 0L) {
    interior_idx <- unique(round(seq(
      from = min(interior_pool),
      to = max(interior_pool),
      length.out = min(interior_capacity, length(interior_pool))
    )))
  }

  keep_idx <- sort(unique(c(left_idx, interior_idx, right_idx)))
  out[keep_idx] <- labels[keep_idx]

  # Add a few small ellipses in long blank stretches to make it visually clear
  # that peak coordinates are intentionally sampled rather than missing.
  if (isTRUE(add_ellipsis)) {
    gaps <- split(seq_len(n), cumsum(seq_len(n) %in% keep_idx))
    gaps <- gaps[lengths(gaps) >= 10L]
    if (length(gaps) > 0) {
      ellipsis_idx <- vapply(gaps, function(x) x[ceiling(length(x) / 2)], integer(1))
      out[ellipsis_idx] <- "..."
    }
  }

  out
}

group_palette <- setNames(
  grDevices::colorRampPalette(c(
    "#F7B7A3", "#F3DFA2", "#B8E0D2", "#A8DADC", "#BFD7EA",
    "#CDB4DB", "#FFC8DD", "#D9D9D9"
  ))(length(unique(owner))),
  unique(owner)
)

peak_labels <- make_peak_coordinate_labels(
  plot_peak_ids,
  keep_each_side = show_peak_labels_each_side,
  max_labels = show_peak_labels_total,
  add_ellipsis = TRUE
)

top_anno <- HeatmapAnnotation(
  cluster = owner[as.character(selected)],
  col = list(cluster = group_palette),
  annotation_name_gp = gpar(fontsize = 7),
  annotation_legend_param = list(
    title = "Marker cluster",
    title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6),
    nrow = 2
  ),
  simple_anno_size = unit(2.2, "mm")
)

ht <- Heatmap(
  z,
  name = "Column Z-score",
  col = colorRamp2(c(-2, 0, 2), c("#C9DFF2", "#FAFAFA", "#D99A9A")),
  top_annotation = top_anno,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 7),
  column_names_side = "top",
  column_labels = peak_labels,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 5.2),
  column_title = "Marker peaks (P labels; see peak_labels.tsv)",
  column_title_gp = gpar(fontsize = 7, fontface = "bold"),
  border = TRUE,
  rect_gp = gpar(col = NA),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6),
    legend_height = unit(22, "mm")
  )
)

pdf(out_pdf, width = 7.2, height = 4.8, useDingbats = FALSE)
draw(
  ht,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  padding = unit(c(3, 3, 3, 3), "mm")
)
dev.off()

png(out_png, width = 7.2, height = 4.8, units = "in", res = 600, type = "cairo")
draw(
  ht,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  padding = unit(c(3, 3, 3, 3), "mm")
)
dev.off()

message("Done. Marker peak heatmap saved to:")
message("  ", out_pdf)
message("  ", out_png)
message("  ", out_label_table)
