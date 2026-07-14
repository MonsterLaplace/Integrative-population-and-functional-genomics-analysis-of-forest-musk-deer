suppressPackageStartupMessages({
  library(data.table)
  library(S4Vectors)
  library(GenomicRanges)
})

outDir <- "06.multiome/03.archr/ArchR_proj/MarkerTables"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

safe_as_dt <- function(x) {
  # 尽量把不同对象转成 data.table
  if (is.null(x)) return(NULL)

  if (inherits(x, "data.table")) {
    return(copy(x))
  }

  if (inherits(x, "data.frame")) {
    return(as.data.table(x))
  }

  if (inherits(x, "DataFrame")) {
    return(as.data.table(as.data.frame(x)))
  }

  if (inherits(x, "GRanges")) {
    return(as.data.table(as.data.frame(x)))
  }

  # 兜底
  return(as.data.table(as.data.frame(x)))
}

export_marker_object <- function(markerObj, prefix, outDir) {
  message("Exporting: ", prefix)
  message("Class: ", paste(class(markerObj), collapse = ", "))

  # --------------------------------------------------
  # 情况1：named list
  # --------------------------------------------------
  if (is.list(markerObj) && !is.data.frame(markerObj) && !inherits(markerObj, "DataFrame")) {
    nms <- names(markerObj)

    # 1a. 有 names
    if (!is.null(nms) && any(nzchar(nms))) {
      summary_list <- list()

      for (nm in nms) {
        x <- markerObj[[nm]]
        dt <- safe_as_dt(x)
        if (is.null(dt) || nrow(dt) == 0) next

        dt[, cluster := nm]

        fwrite(dt, file.path(outDir, paste0(prefix, "_", nm, ".tsv")), sep = "\t")
        fwrite(dt, file.path(outDir, paste0(prefix, "_", nm, ".csv")))

        summary_list[[nm]] <- data.table(cluster = nm, n_markers = nrow(dt))
      }

      if (length(summary_list) > 0) {
        summary_dt <- rbindlist(summary_list, fill = TRUE)
        fwrite(summary_dt, file.path(outDir, paste0(prefix, "_summary.tsv")), sep = "\t")
        fwrite(summary_dt, file.path(outDir, paste0(prefix, "_summary.csv")))
      }

      return(invisible(TRUE))
    }

    # 1b. list 但没 names，用序号导出
    summary_list <- list()
    for (i in seq_along(markerObj)) {
      x <- markerObj[[i]]
      dt <- safe_as_dt(x)
      if (is.null(dt) || nrow(dt) == 0) next

      cl <- paste0("cluster_", i)
      dt[, cluster := cl]

      fwrite(dt, file.path(outDir, paste0(prefix, "_", cl, ".tsv")), sep = "\t")
      fwrite(dt, file.path(outDir, paste0(prefix, "_", cl, ".csv")))

      summary_list[[cl]] <- data.table(cluster = cl, n_markers = nrow(dt))
    }

    if (length(summary_list) > 0) {
      summary_dt <- rbindlist(summary_list, fill = TRUE)
      fwrite(summary_dt, file.path(outDir, paste0(prefix, "_summary.tsv")), sep = "\t")
      fwrite(summary_dt, file.path(outDir, paste0(prefix, "_summary.csv")))
    }

    return(invisible(TRUE))
  }

  # --------------------------------------------------
  # 情况2：单个 data.frame / DataFrame / GRanges
  # 尝试按常见分组列拆分
  # --------------------------------------------------
  dt <- safe_as_dt(markerObj)
  if (is.null(dt) || nrow(dt) == 0) {
    warning(prefix, ": object is empty.")
    return(invisible(FALSE))
  }

  fwrite(dt, file.path(outDir, paste0(prefix, "_all.tsv")), sep = "\t")
  fwrite(dt, file.path(outDir, paste0(prefix, "_all.csv")))

  possible_group_cols <- c("cluster", "Cluster", "group", "Group", "name", "Name")
  group_col <- possible_group_cols[possible_group_cols %in% colnames(dt)]

  if (length(group_col) > 0) {
    gc <- group_col[1]
    split_list <- split(dt, by = gc, keep.by = TRUE)

    summary_list <- list()
    for (nm in names(split_list)) {
      subdt <- split_list[[nm]]
      safe_nm <- gsub("[^A-Za-z0-9_.-]", "_", nm)

      fwrite(subdt, file.path(outDir, paste0(prefix, "_", safe_nm, ".tsv")), sep = "\t")
      fwrite(subdt, file.path(outDir, paste0(prefix, "_", safe_nm, ".csv")))

      summary_list[[nm]] <- data.table(cluster = nm, n_markers = nrow(subdt))
    }

    if (length(summary_list) > 0) {
      summary_dt <- rbindlist(summary_list, fill = TRUE)
      fwrite(summary_dt, file.path(outDir, paste0(prefix, "_summary.tsv")), sep = "\t")
      fwrite(summary_dt, file.path(outDir, paste0(prefix, "_summary.csv")))
    }
  } else {
    message(prefix, ": no obvious grouping column found, exported only *_all.tsv/csv")
  }

  invisible(TRUE)
}

# --------------------------------------------------
# marker peaks
# --------------------------------------------------
peakFile <- "06.multiome/03.archr/ArchR_proj/markerPeakList.rds"
if (file.exists(peakFile)) {
  markerPeakList <- readRDS(peakFile)
  export_marker_object(markerPeakList, prefix = "markerPeaks", outDir = outDir)
} else {
  message("Not found: ", peakFile)
}

# --------------------------------------------------
# marker genes
# --------------------------------------------------
geneFile <- "06.multiome/03.archr/ArchR_proj/markerGeneList.rds"
if (file.exists(geneFile)) {
  markerGeneList <- readRDS(geneFile)
  export_marker_object(markerGeneList, prefix = "markerGenes", outDir = outDir)
} else {
  message("Not found: ", geneFile)
}

message("Done. Output dir: ", outDir)
