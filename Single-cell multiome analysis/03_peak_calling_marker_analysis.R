suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome)
  library(BSgenome.FMdeer.Custom.v1)
  library(data.table)
})

addArchRThreads(threads = 16)

proj <- readRDS("06.multiome/03.archr/ArchR_proj/ArchRProject.afterQC.rds")

# =========================================================
# 1. Calculate genome size from FAI
# =========================================================
faiFile <- "01.reference/FMdeer.primary.fa.fai"
if (!file.exists(faiFile)) {
  stop("FAI file not found: ", faiFile)
}

fai <- fread(faiFile, header = FALSE)
genomeSize <- sum(as.numeric(fai[[2]]), na.rm = TRUE)

message("Genome size calculated from FAI: ", genomeSize)

# =========================================================
# 2. Add group coverages
# =========================================================
proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = "Clusters"
)

# =========================================================
# 3. Call reproducible peaks
# =========================================================
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "Clusters",
  pathToMacs2 = "/data/xb/miniconda3/envs/macs3_env/bin/macs3",
  genomeSize = genomeSize
)

# =========================================================
# 4. Add peak matrix
# =========================================================
proj <- addPeakMatrix(proj)

# =========================================================
# 5. Marker peaks
# =========================================================
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

saveRDS(markersPeaks, "06.multiome/03.archr/ArchR_proj/markersPeaks.rds")

markerList <- getMarkers(
  markersPeaks,
  cutOff = "FDR <= 0.01 & Log2FC >= 1"
)

saveRDS(markerList, "06.multiome/03.archr/ArchR_proj/markerPeakList.rds")

# =========================================================
# 6. Heatmap
# =========================================================
peak_heatmap_palette <- grDevices::colorRampPalette(c(
  "#C9DFF2", # soft blue
  "#FAFAFA", # near-white midpoint
  "#D99A9A"  # soft rose
))(100)

compress_peak_coordinate_labels <- function(labels, keep_each_side = 3L) {
  if (is.null(labels)) return(labels)
  labels <- as.character(labels)
  n <- length(labels)
  if (n <= keep_each_side * 2 + 1) return(labels)

  # Only compress coordinate-like labels; leave simple cluster labels alone.
  coord_like <- grepl("^(chr|Chr|CHR)[A-Za-z0-9._-]+[:-][0-9,]+[-:][0-9,]+", labels)
  if (mean(coord_like, na.rm = TRUE) < 0.25) return(labels)

  out <- rep("", n)
  left_idx <- seq_len(min(keep_each_side, n))
  right_idx <- seq(max(1, n - keep_each_side + 1), n)
  mid_idx <- ceiling(n / 2)
  out[left_idx] <- labels[left_idx]
  out[right_idx] <- labels[right_idx]
  out[mid_idx] <- "..."
  out
}

clean_archr_marker_peak_heatmap <- function(ht) {
  clean_one_heatmap <- function(h) {
    if (is.null(h) || !inherits(h, "Heatmap")) return(h)

    # Compress very dense peak-coordinate labels: show a few on the left and
    # right edges, and use a single ellipsis in the middle.
    h <- tryCatch({
      if ("column_names_param" %in% slotNames(h) && !is.null(h@column_names_param)) {
        h@column_names_param$show <- TRUE
        if (!is.null(h@column_names_param$labels)) {
          h@column_names_param$labels <- compress_peak_coordinate_labels(
            h@column_names_param$labels,
            keep_each_side = 3L
          )
        }
        h@column_names_param$rot <- 60
        if (is.null(h@column_names_param$gp)) {
          h@column_names_param$gp <- grid::gpar(fontsize = 5)
        } else {
          h@column_names_param$gp$fontsize <- 5
        }
      }
      h
    }, error = function(e) {
      message("Skipping peak-coordinate label compression for one heatmap component: ", e$message)
      h
    })

    # Simplify ArchR cluster labels, e.g. C1 -> 1, to reduce right-side crowding.
    h <- tryCatch({
      if ("row_names_param" %in% slotNames(h) &&
          !is.null(h@row_names_param) &&
          !is.null(h@row_names_param$labels)) {
        h@row_names_param$labels <- sub("^C(?=\\d+$)", "", h@row_names_param$labels, perl = TRUE)
      }
      h
    }, error = function(e) {
      message("Skipping row-label simplification for one heatmap component: ", e$message)
      h
    })
    h
  }

  if (is.null(ht)) return(ht)
  if (inherits(ht, "HeatmapList")) {
    ht <- tryCatch({
      for (i in seq_along(ht@ht_list)) {
        if (!is.null(ht@ht_list[[i]]) && inherits(ht@ht_list[[i]], "Heatmap")) {
          ht@ht_list[[i]] <- clean_one_heatmap(ht@ht_list[[i]])
        }
      }
      ht
    }, error = function(e) {
      message("Skipping HeatmapList-level cleanup: ", e$message)
      ht
    })
  } else if (inherits(ht, "Heatmap")) {
    ht <- clean_one_heatmap(ht)
  }
  ht
}

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks,
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  transpose = TRUE,
  labelMarkers = NULL,
  pal = peak_heatmap_palette
)
heatmapPeaks <- clean_archr_marker_peak_heatmap(heatmapPeaks)

plotPDF(
  heatmapPeaks,
  name = "03_markerPeaks_heatmap.pdf",
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  ArchRProj = proj
)

# =========================================================
# 7. Save project
# =========================================================
saveRDS(proj, file = "06.multiome/03.archr/ArchR_proj/ArchRProject.afterPeaks.rds")
