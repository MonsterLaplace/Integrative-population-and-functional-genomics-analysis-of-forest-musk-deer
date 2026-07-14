suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

outDir <- "06.multiome/04.seurat_signac"
objFile <- file.path(outDir, "FMdeer_multiome.wnn.final_annotated.linkpeaks.rds")

candidateFile1 <- file.path(outDir, "musk_function_genes.strict_confidence.re_scored.tsv")
candidateFile2 <- file.path(outDir, "musk_function_genes.high_confidence.re_scored.tsv")
candidateFile3 <- file.path(outDir, "musk_function_genes.tsv")

plotDir <- file.path(outDir, "peak2gene_coverage_plots_sci")
dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)

# Keep all auxiliary outputs inside the current analysis output tree.
# Also write a compatibility status table to the legacy fallback directory,
# because 36c_plot_peak2gene_and_coverage_per_gene.R reads that location.
localOutputDir <- plotDir
compatOutputDir <- file.path(outDir, "peak2gene_coverage_plots_fallback")
dir.create(localOutputDir, recursive = TRUE, showWarnings = FALSE)
dir.create(compatOutputDir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------
# user-tunable display settings
# --------------------------------------------------
# Plot all annotated cell types. Keep labels as numeric IDs in the coverage
# tracks and export a separate ID legend to avoid unreadable left-side labels.
use_all_celltypes <- TRUE

# Publication-style export settings.
plot_width <- 14
plot_height_min <- 11
plot_height_max <- 24
png_dpi <- 600
coverage_window <- 18
coverage_ymax <- "q90"
downsample_rate <- 0.35
max_downsample <- 12000

# Region selection. Long windows make peak signals look needle-like. When the
# gene plus linked peaks span more than max_focus_width, automatically zoom to
# the densest linked-peak neighborhood.
target_focus_width <- 45000
max_focus_width <- 65000
min_focus_width <- 22000
dense_peak_top_n <- 6L
promoter_focus_flank <- 25000

# Highlight policy. Automatic red coverage backgrounds are disabled by default:
# LinkPeaks-derived ranges can be link spans rather than exact peak intervals,
# which makes the red bars look misplaced in genes such as SMPDL3A and NRCAM.
# The Peaks track still shows peak positions. Add exact coordinates to
# manual_peak_highlight_overrides if a specific peak needs to be emphasized.
highlight_gene_body <- FALSE
highlight_promoter <- FALSE
highlight_linked_peaks <- FALSE
max_highlight_peak_width <- 2500
linked_peak_highlight_color <- "#D73027"

manual_peak_highlight_overrides <- list(
  # Example:
  # SMPDL3A = data.table(chr = "Chr10", start = 33409500, end = 33410500)
)

# For genes that need a publication-style zoomed-in view,
# define manual plotting ranges here.
# Example requested by user: KCNC2 only show 8,000,000-8,050,000 bp
manual_region_overrides <- list(
  KCNC2 = list(chr = "Chr06", start = 8000000, end = 8050000)
)

# Optional legacy filter. It is disabled by default because the current figure
# should show all cell types.
secretion_related_patterns <- c(
  "secretory",
  "epithelial",
  "krt14",
  "krt79",
  "klk",
  "spink",
  "metabolic"
)

# Preferred ordering for the remaining cell types in the main figure.
preferred_celltype_order <- c(
  "Secretory_epithelial_cells",
  "Specialized_epithelial_cells",
  "Specialized_epithelial_cells_2",
  "Epithelial_cells_lipid_metabolic",
  "Metabolic_epithelial_cells",
  "Epithelial_cells_KRT14",
  "Epithelial_cells_KRT79",
  "Epithelial_cells_KRT79_like",
  "Epithelial_cells_KLK_SPINK5",
  "Epithelial_cells_WNT"
)

normalize_label <- function(x) {
  tolower(gsub("[^A-Za-z0-9]+", "_", x))
}

if (!file.exists(objFile)) stop("Object not found: ", objFile)

candidateFile <- NULL
for (f in c(candidateFile1, candidateFile2, candidateFile3)) {
  if (file.exists(f)) {
    candidateFile <- f
    break
  }
}
if (is.null(candidateFile)) stop("No candidate gene file found.")

message("Using candidate file: ", candidateFile)

obj <- readRDS(objFile)
cand <- fread(candidateFile)

if (!"gene" %in% colnames(cand)) stop("Candidate file missing gene column.")
if (!"celltype" %in% colnames(obj@meta.data)) stop("Metadata missing 'celltype' column.")

celltype_raw <- as.character(obj@meta.data$celltype)
celltype_raw <- celltype_raw[!is.na(celltype_raw) & nzchar(celltype_raw)]
if (is.factor(obj@meta.data$celltype)) {
  all_celltypes <- levels(obj@meta.data$celltype)
  all_celltypes <- all_celltypes[all_celltypes %in% celltype_raw]
} else {
  all_celltypes <- unique(celltype_raw)
}

if (use_all_celltypes) {
  selected_celltypes <- all_celltypes
} else {
  all_celltypes_norm <- normalize_label(all_celltypes)
  keep_idx <- Reduce(
    `|`,
    lapply(secretion_related_patterns, function(pat) {
      grepl(normalize_label(pat), all_celltypes_norm, fixed = TRUE)
    })
  )

  selected_celltypes <- all_celltypes[keep_idx]
  if (length(selected_celltypes) == 0) {
    warning("No secretion-related cell types matched the current pattern list; using all cell types.")
    selected_celltypes <- all_celltypes
  }
}

preferred_norm <- normalize_label(preferred_celltype_order)
selected_norm <- normalize_label(selected_celltypes)
ordered_hits <- selected_celltypes[match(preferred_norm, selected_norm, nomatch = 0)]
selected_celltypes <- unique(c(ordered_hits, selected_celltypes[!selected_norm %in% preferred_norm]))

obj <- subset(
  x = obj,
  cells = rownames(obj@meta.data)[as.character(obj@meta.data$celltype) %in% selected_celltypes]
)

celltype_levels <- unique(as.character(obj@meta.data$celltype))
celltype_map <- data.table(
  celltype = celltype_levels,
  celltype_id = as.character(seq_along(celltype_levels))
)
celltype_map[, n_cells := as.integer(table(as.character(obj@meta.data$celltype))[celltype])]
obj@meta.data$celltype_id <- factor(
  celltype_map$celltype_id[match(as.character(obj@meta.data$celltype), celltype_map$celltype)],
  levels = celltype_map$celltype_id
)
fwrite(
  celltype_map,
  file.path(localOutputDir, "36b_celltype_id_legend.tsv"),
  sep = "\t"
)
fwrite(
  data.table(selected_celltype = selected_celltypes),
  file.path(localOutputDir, "36b_selected_musk_secretion_related_celltypes.tsv"),
  sep = "\t"
)

fwrite(
  celltype_map,
  file.path(localOutputDir, "36b_all_plotted_celltypes.tsv"),
  sep = "\t"
)

message("Cell types included in coverage plots: ", length(celltype_levels))

anno <- Annotation(obj[["ATAC"]])
if (is.null(anno) || length(anno) == 0) {
  stop("ATAC annotation missing.")
}

anno_df <- as.data.table(as.data.frame(anno))

sanitize_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

coalesce_chr <- function(x, y) {
  if (!is.null(x) && length(x) > 0 && !is.na(x[[1]]) && nzchar(x[[1]])) x[[1]] else y
}

first_nonempty <- function(...) {
  vals <- unlist(list(...), use.names = FALSE)
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (length(vals) == 0) return(NA_character_)
  vals[[1]]
}

find_first_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

match_gene_to_annotation <- function(g, anno_df) {
  g2 <- gsub("-", "_", g)

  if ("gene_name" %in% colnames(anno_df)) {
    hit <- anno_df[gene_name == g]
    if (nrow(hit) > 0) return(hit[1])
  }

  if ("gene_id" %in% colnames(anno_df)) {
    hit <- anno_df[gene_id == g]
    if (nrow(hit) > 0) return(hit[1])

    hit <- anno_df[gene_id == g2]
    if (nrow(hit) > 0) return(hit[1])
  }

  if ("gene_name" %in% colnames(anno_df)) {
    hit <- anno_df[gene_name == g2]
    if (nrow(hit) > 0) return(hit[1])
  }

  NULL
}

calc_flank_size <- function(start, end,
                            min_flank = 8000,
                            max_flank = 20000,
                            target_total_width = 40000) {
  gene_width <- max(1, as.numeric(end) - as.numeric(start) + 1)
  flank <- ceiling((target_total_width - gene_width) / 2)
  flank <- max(min_flank, flank)
  flank <- min(max_flank, flank)
  flank
}

parse_peak_strings <- function(peak_strings) {
  peak_strings <- as.character(peak_strings)
  peak_strings <- peak_strings[!is.na(peak_strings) & nzchar(peak_strings)]
  if (length(peak_strings) == 0) return(NULL)

  peak_strings <- gsub(":", "-", peak_strings)
  peak_strings <- gsub("_", "-", peak_strings)

  parsed <- rbindlist(lapply(peak_strings, function(x) {
    parts <- strsplit(x, "-", fixed = TRUE)[[1]]
    if (length(parts) < 3) return(NULL)
    data.table(
      seqnames = parts[[1]],
      start = suppressWarnings(as.numeric(parts[[2]])),
      end = suppressWarnings(as.numeric(parts[[3]]))
    )
  }), fill = TRUE)

  if (nrow(parsed) == 0) return(NULL)
  parsed <- parsed[!is.na(start) & !is.na(end) & end >= start]
  if (nrow(parsed) == 0) return(NULL)

  GRanges(
    seqnames = parsed$seqnames,
    ranges = IRanges(start = parsed$start, end = parsed$end)
  )
}

get_gene_specific_links <- function(obj, gene_candidates, seqname_keep = NULL) {
  lk <- tryCatch(Links(obj[["ATAC"]]), error = function(e) NULL)
  if (is.null(lk) || length(lk) == 0) return(NULL)

  lk_df <- as.data.table(as.data.frame(lk))
  if (nrow(lk_df) == 0) return(NULL)

  gene_col <- find_first_col(
    lk_df,
    c("gene", "gene_name", "gene.id", "gene_id", "target", "feature", "name")
  )
  if (is.null(gene_col)) return(NULL)

  keep_genes <- unique(gene_candidates[!is.na(gene_candidates) & nzchar(gene_candidates)])
  hit <- lk_df[get(gene_col) %in% keep_genes]
  if (!is.null(seqname_keep) && "seqnames" %in% colnames(hit)) {
    hit <- hit[as.character(seqnames) == as.character(seqname_keep)]
  }
  if (nrow(hit) == 0) return(NULL)

  if (!all(c("start", "end") %in% colnames(hit))) return(NULL)

  peak_col <- find_first_col(
    hit,
    c("peak", "peak_name", "peak.id", "peak_id", "query_region", "query", "region")
  )
  gr <- NULL
  if (!is.null(peak_col)) {
    gr <- parse_peak_strings(hit[[peak_col]])
    if (!is.null(seqname_keep) && !is.null(gr) && length(gr) > 0) {
      gr <- gr[as.character(seqnames(gr)) == as.character(seqname_keep)]
    }
  }

  if (is.null(gr) || length(gr) == 0) {
    gr <- GRanges(
      seqnames = hit$seqnames,
      ranges = IRanges(start = hit$start, end = hit$end)
    )
  }

  meta_cols <- setdiff(colnames(hit), c("seqnames", "start", "end", "width", "strand"))
  if (length(meta_cols) > 0 && length(gr) == nrow(hit)) {
    mcols(gr) <- S4Vectors::DataFrame(hit[, ..meta_cols])
  }

  gr
}

keep_narrow_peak_ranges <- function(gr, max_width = max_highlight_peak_width) {
  if (is.null(gr) || length(gr) == 0) return(NULL)
  gr2 <- gr[width(gr) <= max_width]
  if (length(gr2) == 0) return(NULL)
  gr2
}

rank_linked_peaks <- function(linked_peak_gr) {
  if (is.null(linked_peak_gr) || length(linked_peak_gr) == 0) return(integer(0))

  md <- as.data.table(as.data.frame(mcols(linked_peak_gr)))
  candidate_score_cols <- intersect(
    c("score", "zscore", "z_score", "correlation", "Correlation", "avg_log2FC", "pvalue", "p_val"),
    colnames(md)
  )

  if (length(candidate_score_cols) > 0) {
    score_col <- candidate_score_cols[[1]]
    score <- suppressWarnings(as.numeric(md[[score_col]]))
    if (score_col %in% c("pvalue", "p_val")) {
      score <- -log10(pmax(score, .Machine$double.xmin))
    } else {
      score <- abs(score)
    }
    score[is.na(score)] <- 0
    return(order(score, decreasing = TRUE))
  }

  order(width(linked_peak_gr), decreasing = TRUE)
}

make_peak_focused_region <- function(gene_gr, linked_peak_gr,
                                     target_width = target_focus_width,
                                     max_width = max_focus_width,
                                     top_n = dense_peak_top_n) {
  if (is.null(linked_peak_gr) || length(linked_peak_gr) == 0) {
    gene_center <- round((start(gene_gr)[1] + end(gene_gr)[1]) / 2)
    half_width <- round(target_width / 2)
    return(GRanges(
      seqnames = as.character(seqnames(gene_gr))[1],
      ranges = IRanges(start = max(1, gene_center - half_width), end = gene_center + half_width)
    ))
  }

  ranked_idx <- rank_linked_peaks(linked_peak_gr)
  ranked_idx <- ranked_idx[seq_len(min(length(ranked_idx), top_n))]
  peak_mid <- round((start(linked_peak_gr) + end(linked_peak_gr)) / 2)
  candidate_centers <- unique(c(
    peak_mid[ranked_idx],
    round((start(gene_gr)[1] + end(gene_gr)[1]) / 2),
    start(gene_gr)[1],
    end(gene_gr)[1]
  ))

  best <- NULL
  best_score <- -Inf
  half_width <- round(target_width / 2)
  for (center in candidate_centers) {
    win_start <- max(1, center - half_width)
    win_end <- win_start + target_width - 1
    inside <- start(linked_peak_gr) <= win_end & end(linked_peak_gr) >= win_start
    n_inside <- sum(inside)
    gene_overlap <- start(gene_gr)[1] <= win_end && end(gene_gr)[1] >= win_start
    ranked_inside_bonus <- sum(ranked_idx %in% which(inside))
    score <- n_inside * 10 + ranked_inside_bonus * 3 + as.integer(gene_overlap)
    if (score > best_score) {
      best_score <- score
      best <- c(win_start, win_end)
    }
  }

  chosen <- linked_peak_gr[start(linked_peak_gr) <= best[[2]] & end(linked_peak_gr) >= best[[1]]]
  if (length(chosen) > 0) {
    region_start <- min(start(chosen), best[[1]])
    region_end <- max(end(chosen), best[[2]])
  } else {
    region_start <- best[[1]]
    region_end <- best[[2]]
  }

  center <- round((region_start + region_end) / 2)
  final_width <- min(max_width, max(min_focus_width, region_end - region_start + 1))
  half_width <- round(final_width / 2)
  GRanges(
    seqnames = as.character(seqnames(gene_gr))[1],
    ranges = IRanges(start = max(1, center - half_width), end = center + half_width)
  )
}

make_focus_region <- function(gene_gr, linked_peak_gr = NULL,
                              pad_min = 2500,
                              pad_max = 8000,
                              max_width = max_focus_width) {
  base_gr <- gene_gr
  if (!is.null(linked_peak_gr) && length(linked_peak_gr) > 0) {
    base_gr <- suppressWarnings(c(gene_gr, linked_peak_gr))
  }

  region_start <- min(start(base_gr))
  region_end <- max(end(base_gr))
  span <- max(1, region_end - region_start + 1)
  if (span > max_width) {
    return(make_peak_focused_region(gene_gr, linked_peak_gr))
  }

  pad <- max(pad_min, ceiling(span * 0.18))
  pad <- min(pad, pad_max)

  region_gr <- GRanges(
    seqnames = as.character(seqnames(gene_gr))[1],
    ranges = IRanges(
      start = max(1, region_start - pad),
      end = region_end + pad
    )
  )

  if (width(region_gr)[1] > max_width) {
    center <- round((start(region_gr)[1] + end(region_gr)[1]) / 2)
    half_width <- round(max_width / 2)
    region_gr <- GRanges(
      seqnames = as.character(seqnames(gene_gr))[1],
      ranges = IRanges(start = max(1, center - half_width), end = center + half_width)
    )
  }

  region_gr
}

get_manual_region_override <- function(gene_candidates, override_list) {
  keep_genes <- unique(gene_candidates[!is.na(gene_candidates) & nzchar(gene_candidates)])
  for (nm in names(override_list)) {
    if (nm %in% keep_genes) return(override_list[[nm]])
  }
  NULL
}

get_manual_peak_highlight_override <- function(gene_candidates, override_list) {
  keep_genes <- unique(gene_candidates[!is.na(gene_candidates) & nzchar(gene_candidates)])
  for (nm in names(override_list)) {
    if (!nm %in% keep_genes) next
    dt <- as.data.table(override_list[[nm]])
    if (!all(c("chr", "start", "end") %in% colnames(dt))) {
      stop("Manual peak highlight for ", nm, " must contain chr/start/end columns.")
    }
    dt <- dt[!is.na(chr) & !is.na(start) & !is.na(end) & end >= start]
    if (nrow(dt) == 0) return(NULL)
    gr <- GRanges(
      seqnames = dt$chr,
      ranges = IRanges(start = as.numeric(dt$start), end = as.numeric(dt$end))
    )
    gr$color <- linked_peak_highlight_color
    return(gr)
  }
  NULL
}

make_gene_body_highlight <- function(hit_row) {
  gr <- GRanges(
    seqnames = hit_row$seqnames[[1]],
    ranges = IRanges(start = hit_row$start[[1]], end = hit_row$end[[1]])
  )
  gr$color <- "#E6F2FF"
  gr
}

make_promoter_highlight <- function(hit_row, promoter_up = 1500, promoter_down = 500) {
  strand_use <- if ("strand" %in% colnames(hit_row)) as.character(hit_row$strand[[1]]) else "*"
  start_use <- as.numeric(hit_row$start[[1]])
  end_use <- as.numeric(hit_row$end[[1]])

  if (strand_use == "-") {
    prom_start <- max(1, end_use - promoter_down)
    prom_end <- end_use + promoter_up
  } else {
    prom_start <- max(1, start_use - promoter_up)
    prom_end <- start_use + promoter_down
  }

  gr <- GRanges(
    seqnames = hit_row$seqnames[[1]],
    ranges = IRanges(start = prom_start, end = prom_end)
  )
  gr$color <- "#FEE8C8"
  gr
}

trim_ranges_to_region <- function(gr, region_gr) {
  if (is.null(gr) || length(gr) == 0) return(gr)

  seq_keep <- as.character(seqnames(region_gr))[1]
  region_start <- start(region_gr)[1]
  region_end <- end(region_gr)[1]
  keep <- as.character(seqnames(gr)) == seq_keep & end(gr) >= region_start & start(gr) <= region_end
  if (!any(keep)) return(gr[FALSE])

  gr2 <- gr[keep]
  start(gr2) <- pmax(start(gr2), region_start)
  end(gr2) <- pmin(end(gr2), region_end)
  gr2
}

calc_plot_height <- function(n_celltypes) {
  max(plot_height_min, min(plot_height_max, 6.5 + 0.30 * n_celltypes))
}

style_coverage_plot <- function(p, title_text, n_celltypes, highlight_note = NULL) {
  subtitle_text <- if (is.null(highlight_note) || !nzchar(highlight_note)) {
    paste0(n_celltypes, " cell types shown as numeric IDs; peak positions are shown in the Peaks track.")
  } else {
    paste0(n_celltypes, " cell types shown as numeric IDs; ", highlight_note)
  }

  p +
    plot_annotation(
      title = title_text,
      subtitle = subtitle_text,
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0, margin = margin(b = 3)),
        plot.subtitle = element_text(size = 9.5, colour = "grey25", hjust = 0, margin = margin(b = 8))
      )
    ) &
    theme_classic(base_size = 9, base_family = "sans") &
    theme(
      axis.title.x = element_text(size = 9.5, face = "bold", colour = "black"),
      axis.title.y = element_text(size = 9.5, face = "bold", colour = "black"),
      axis.text.x = element_text(size = 8.5, colour = "black"),
      axis.text.y = element_text(size = 7.5, colour = "black"),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.30, colour = "black"),
      strip.text = element_text(size = 7.8, face = "bold", colour = "black"),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.35),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(6, 10, 6, 6)
    )
}

save_plot_dual <- function(p, out_prefix, width = plot_width, height = calc_plot_height(length(celltype_levels))) {
  ggsave(
    filename = paste0(out_prefix, ".pdf"),
    plot = p,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )

  ggsave(
    filename = paste0(out_prefix, ".png"),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = png_dpi,
    bg = "white"
  )

  ggsave(
    filename = paste0(out_prefix, ".tiff"),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = png_dpi,
    compression = "lzw",
    bg = "white"
  )
}

save_celltype_legend_plot <- function(celltype_map, out_prefix, ncol = NULL) {
  legend_dt <- copy(celltype_map)
  if (is.null(ncol)) {
    ncol <- if (nrow(legend_dt) <= 18) 2 else if (nrow(legend_dt) <= 36) 3 else 4
  }
  legend_dt[, label := paste0(celltype_id, "  =  ", celltype, "  (n=", n_cells, ")")]
  legend_dt[, idx := .I]
  nrow_block <- ceiling(nrow(legend_dt) / ncol)
  legend_dt[, col_id := ((idx - 1) %/% nrow_block) + 1]
  legend_dt[, row_id := ((idx - 1) %% nrow_block) + 1]

  p <- ggplot(legend_dt, aes(x = col_id, y = -row_id, label = label)) +
    geom_text(hjust = 0, size = 3.4, family = "sans", colour = "black") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.15))) +
    theme_void(base_size = 10, base_family = "sans") +
    labs(title = "Cell type ID legend") +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.margin = margin(12, 20, 12, 20)
    )

  ggsave(
    filename = paste0(out_prefix, ".pdf"),
    plot = p,
    width = 12,
    height = max(4.5, 0.24 * nrow_block + 1.4),
    units = "in",
    bg = "white"
  )
  ggsave(
    filename = paste0(out_prefix, ".png"),
    plot = p,
    width = 12,
    height = max(4.5, 0.24 * nrow_block + 1.4),
    units = "in",
    dpi = png_dpi,
    bg = "white"
  )
}

cand2 <- copy(cand)
if (!"peak2gene_support" %in% colnames(cand2)) cand2[, peak2gene_support := 0L]
if (!"musk_score_raw" %in% colnames(cand2)) cand2[, musk_score_raw := NA_real_]

cand2 <- cand2[order(-peak2gene_support, -musk_score_raw)]

manual_genes <- c("KCNC2", "IL7R", "SMPDL3A", "NRCAM", "COL1A1", "ACTG2", "VWF", "GABRR2")
auto_genes <- unique(cand2$gene[1:min(12, nrow(cand2))])
plot_genes <- unique(c(manual_genes, auto_genes))

message("Genes selected for plotting:")
print(plot_genes)

save_celltype_legend_plot(
  celltype_map = celltype_map,
  out_prefix = file.path(plotDir, "36b_celltype_id_legend")
)

status_dt <- data.table(
  input_gene = character(),
  matched_gene_name = character(),
  matched_gene_id = character(),
  success = logical(),
  mode = character(),
  flank_bp = numeric(),
  region_strategy = character(),
  region_width_bp = numeric(),
  n_linked_peaks_total = integer(),
  n_linked_peaks_narrow = integer(),
  n_linked_peaks_plotted = integer(),
  highlight_strategy = character(),
  note = character()
)
linked_peak_candidates_dt <- data.table(
  input_gene = character(),
  matched_gene_name = character(),
  seqnames = character(),
  start = numeric(),
  end = numeric(),
  width = numeric(),
  plotted_in_current_region = logical()
)

combined_pdf <- file.path(plotDir, "36b_peak2gene_coverage_candidate_genes.sci.pdf")
pdf(
  combined_pdf,
  width = plot_width,
  height = calc_plot_height(length(celltype_levels)),
  onefile = TRUE
)

for (g in plot_genes) {
  message("Processing gene: ", g)

  hit <- match_gene_to_annotation(g, anno_df)

  if (is.null(hit) || nrow(hit) == 0) {
    status_dt <- rbind(
      status_dt,
      data.table(
        input_gene = g,
        matched_gene_name = NA_character_,
        matched_gene_id = NA_character_,
        success = FALSE,
        mode = "no_match",
        flank_bp = NA_real_,
        region_strategy = "no_match",
        region_width_bp = NA_real_,
        n_linked_peaks_total = NA_integer_,
        n_linked_peaks_narrow = NA_integer_,
        n_linked_peaks_plotted = NA_integer_,
        highlight_strategy = "no_match",
        note = "No matching gene found in ATAC annotation"
      ),
      fill = TRUE
    )
    next
  }

  gene_name_use <- if ("gene_name" %in% colnames(hit)) as.character(hit$gene_name[[1]]) else NA_character_
  gene_id_use <- if ("gene_id" %in% colnames(hit)) as.character(hit$gene_id[[1]]) else NA_character_
  gene_label <- coalesce_chr(gene_name_use, g)
  feature_use <- first_nonempty(gene_name_use, gene_id_use, g)

  gene_gr <- GRanges(
    seqnames = hit$seqnames[[1]],
    ranges = IRanges(start = hit$start[[1]], end = hit$end[[1]])
  )

  linked_peak_gr <- get_gene_specific_links(
    obj = obj,
    gene_candidates = unique(c(g, gene_name_use, gene_id_use)),
    seqname_keep = hit$seqnames[[1]]
  )
  linked_peak_gr_narrow <- keep_narrow_peak_ranges(linked_peak_gr)

  manual_override <- get_manual_region_override(
    gene_candidates = unique(c(g, gene_name_use, gene_id_use)),
    override_list = manual_region_overrides
  )
  manual_peak_highlight <- get_manual_peak_highlight_override(
    gene_candidates = unique(c(g, gene_name_use, gene_id_use)),
    override_list = manual_peak_highlight_overrides
  )

  if (!is.null(manual_override)) {
    region_gr <- GRanges(
      seqnames = manual_override$chr,
      ranges = IRanges(start = manual_override$start, end = manual_override$end)
    )
    flank_bp <- NA_real_
    region_strategy <- "manual_override"
  } else if (!is.null(linked_peak_gr_narrow) && length(linked_peak_gr_narrow) > 0) {
    full_span <- max(end(suppressWarnings(c(gene_gr, linked_peak_gr_narrow)))) -
      min(start(suppressWarnings(c(gene_gr, linked_peak_gr_narrow)))) + 1
    region_gr <- make_focus_region(gene_gr, linked_peak_gr_narrow)
    flank_bp <- NA_real_
    region_strategy <- if (full_span > max_focus_width) "peak_focused_zoom" else "gene_plus_linked_peaks"
  } else {
    gene_width <- end(gene_gr)[1] - start(gene_gr)[1] + 1
    if (gene_width > max_focus_width) {
      gene_center <- round((start(gene_gr)[1] + end(gene_gr)[1]) / 2)
      region_gr <- GRanges(
        seqnames = hit$seqnames[[1]],
        ranges = IRanges(
          start = max(1, gene_center - promoter_focus_flank),
          end = gene_center + promoter_focus_flank
        )
      )
      flank_bp <- NA_real_
      region_strategy <- "long_gene_center_zoom"
    } else {
      flank_bp <- calc_flank_size(hit$start[[1]], hit$end[[1]])
      region_gr <- GRanges(
        seqnames = hit$seqnames[[1]],
        ranges = IRanges(
          start = max(1, hit$start[[1]] - flank_bp),
          end = hit$end[[1]] + flank_bp
        )
      )
      region_strategy <- "gene_flank"
    }
  }

  region_highlight <- GRanges()
  if (highlight_gene_body) {
    region_highlight <- suppressWarnings(c(region_highlight, make_gene_body_highlight(hit)))
  }
  if (highlight_promoter) {
    region_highlight <- suppressWarnings(c(region_highlight, make_promoter_highlight(hit)))
  }
  region_highlight <- trim_ranges_to_region(region_highlight, region_gr)

  linked_peak_gr_plot <- trim_ranges_to_region(linked_peak_gr_narrow, region_gr)
  manual_peak_highlight_plot <- trim_ranges_to_region(manual_peak_highlight, region_gr)
  n_linked_peaks_total <- if (!is.null(linked_peak_gr)) length(linked_peak_gr) else 0L
  n_linked_peaks_narrow <- if (!is.null(linked_peak_gr_narrow)) length(linked_peak_gr_narrow) else 0L
  n_linked_peaks_plotted <- if (!is.null(linked_peak_gr_plot)) length(linked_peak_gr_plot) else 0L
  if (!is.null(linked_peak_gr_narrow) && length(linked_peak_gr_narrow) > 0) {
    narrow_dt <- as.data.table(as.data.frame(linked_peak_gr_narrow))
    narrow_dt[, input_gene := g]
    narrow_dt[, matched_gene_name := gene_label]
    narrow_dt[, plotted_in_current_region := seqnames == as.character(seqnames(region_gr))[1] &
      start <= end(region_gr)[1] & end >= start(region_gr)[1]]
    linked_peak_candidates_dt <- rbind(
      linked_peak_candidates_dt,
      narrow_dt[, .(
        input_gene,
        matched_gene_name,
        seqnames = as.character(seqnames),
        start = as.numeric(start),
        end = as.numeric(end),
        width = as.numeric(width),
        plotted_in_current_region
      )],
      fill = TRUE
    )
  }

  highlight_strategy <- "disabled"
  highlight_note <- NULL
  if (!is.null(manual_peak_highlight_plot) && length(manual_peak_highlight_plot) > 0) {
    region_highlight <- suppressWarnings(c(region_highlight, manual_peak_highlight_plot))
    highlight_strategy <- "manual_exact_peak"
    highlight_note <- "red shading marks manually curated peak coordinates."
  } else if (highlight_linked_peaks && !is.null(linked_peak_gr_plot) && length(linked_peak_gr_plot) > 0) {
    linked_peak_gr_plot$color <- linked_peak_highlight_color
    region_highlight <- suppressWarnings(c(region_highlight, linked_peak_gr_plot))
    highlight_strategy <- "automatic_narrow_linked_peaks"
    highlight_note <- "red shading marks narrow linked peaks inferred from Links metadata."
  }

  has_links <- !is.null(linked_peak_gr_plot) && length(linked_peak_gr_plot) > 0
  region_highlight_use <- if (length(region_highlight) > 0) region_highlight else NULL
  title_use <- if (has_links) {
    paste0(gene_label, " | peak-to-gene accessibility")
  } else {
    paste0(gene_label, " | chromatin accessibility")
  }
  if (!is.null(manual_override)) {
    title_use <- paste0(
      gene_label,
      " | zoomed chromatin accessibility"
    )
  }
  file_stub <- file.path(
    plotDir,
    paste0("36b_", sanitize_filename(gene_label), "_coverage_sci")
  )

  ok <- FALSE
  mode <- NA_character_
  note <- ""

  tryCatch({
    p <- CoveragePlot(
      object = obj,
      region = region_gr,
      features = feature_use,
      expression.assay = "RNA",
      peaks = TRUE,
      links = if (has_links) gene_label else FALSE,
      show.bulk = TRUE,
      region.highlight = region_highlight_use,
      group.by = "celltype_id",
      window = coverage_window,
      ymax = coverage_ymax,
      downsample.rate = downsample_rate,
      max.downsample = max_downsample
    )
    p2 <- style_coverage_plot(p, title_use, length(celltype_levels), highlight_note)
    print(p2)
    save_plot_dual(p2, file_stub)
    ok <- TRUE
    mode <- "region+feature+gene_specific_links"
    note <- "success"
  }, error = function(e) {
    note <<- e$message
  })

  if (!ok) {
    tryCatch({
      p <- CoveragePlot(
        object = obj,
        region = region_gr,
        peaks = TRUE,
        links = if (has_links) gene_label else FALSE,
        show.bulk = TRUE,
        region.highlight = region_highlight_use,
        group.by = "celltype_id",
        window = coverage_window,
        ymax = coverage_ymax,
        downsample.rate = downsample_rate,
        max.downsample = max_downsample
      )
      p2 <- style_coverage_plot(p, paste0(gene_label, " | peak-to-gene accessibility"), length(celltype_levels), highlight_note)
      print(p2)
      save_plot_dual(p2, file_stub)
      ok <- TRUE
      mode <- "region+gene_specific_links"
      note <- "success_after_feature_fallback"
    }, error = function(e) {
      note <<- paste(note, " | fallback1:", e$message)
    })
  }

  if (!ok) {
    tryCatch({
      p <- CoveragePlot(
        object = obj,
        region = region_gr,
        peaks = TRUE,
        show.bulk = TRUE,
        region.highlight = region_highlight_use,
        group.by = "celltype_id",
        window = coverage_window,
        ymax = coverage_ymax,
        downsample.rate = downsample_rate,
        max.downsample = max_downsample
      )
      p2 <- style_coverage_plot(p, paste0(gene_label, " | chromatin accessibility"), length(celltype_levels), highlight_note)
      print(p2)
      save_plot_dual(p2, file_stub)
      ok <- TRUE
      mode <- "region_only"
      note <- "success_accessibility_only"
    }, error = function(e) {
      note <<- paste(note, " | fallback2:", e$message)
    })
  }

  status_dt <- rbind(
    status_dt,
    data.table(
      input_gene = g,
      matched_gene_name = gene_name_use,
      matched_gene_id = gene_id_use,
      success = ok,
      mode = mode,
      flank_bp = flank_bp,
      region_strategy = region_strategy,
      region_width_bp = as.numeric(width(region_gr)[1]),
      n_linked_peaks_total = n_linked_peaks_total,
      n_linked_peaks_narrow = n_linked_peaks_narrow,
      n_linked_peaks_plotted = n_linked_peaks_plotted,
      highlight_strategy = highlight_strategy,
      note = note
    ),
    fill = TRUE
  )
}

dev.off()

fwrite(
  status_dt,
  file.path(localOutputDir, "36b_peak2gene_coverage_candidate_genes.sci.status.tsv"),
  sep = "\t"
)

fwrite(
  status_dt,
  file.path(plotDir, "36b_peak2gene_coverage_candidate_genes.status.tsv"),
  sep = "\t"
)

fwrite(
  status_dt,
  file.path(compatOutputDir, "36b_peak2gene_coverage_candidate_genes.status.tsv"),
  sep = "\t"
)

fwrite(
  linked_peak_candidates_dt,
  file.path(plotDir, "36b_linked_peak_candidates_for_manual_highlight.tsv"),
  sep = "\t"
)

message("Done. SCI-style outputs saved to: ", plotDir)
