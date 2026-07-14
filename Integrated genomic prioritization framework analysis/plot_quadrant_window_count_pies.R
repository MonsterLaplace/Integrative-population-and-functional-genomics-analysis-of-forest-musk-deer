#!/usr/bin/env Rscript

# ================================================================
# Pie charts for quadrant window counts
# ------------------------------------------------
# Input:
#   postprocess_results/conflict_quadrant_table.tsv
#   postprocess_results/conflict_quadrant_top1_windows.tsv
#
# Output:
#   postprocess_results/quadrant_all_windows_pie.pdf/png
#   postprocess_results/quadrant_top1_windows_pie.pdf/png
#   postprocess_results/quadrant_window_count_pies.pdf/png
#   postprocess_results/quadrant_all_windows_unique_genes_pie.pdf/png
#   postprocess_results/quadrant_top1_windows_unique_genes_pie.pdf/png
#   postprocess_results/quadrant_unique_gene_count_pies.pdf/png
#   postprocess_results/quadrant_window_counts_summary.tsv
#   postprocess_results/quadrant_unique_gene_counts_summary.tsv
# ================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0 || hit == length(args)) return(default)
  args[hit + 1]
}

all_file <- get_arg("--all", file.path("postprocess_results", "conflict_quadrant_table.tsv"))
top1_file <- get_arg("--top1", file.path("postprocess_results", "conflict_quadrant_top1_windows.tsv"))
outdir <- get_arg("--outdir", "postprocess_results")
prefix <- get_arg("--prefix", "quadrant")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

quadrant_levels <- c(
  "conservation_priority",
  "breeding_function_priority",
  "conflict_quadrant",
  "background"
)

quadrant_labels <- c(
  conservation_priority = "Conservation priority",
  breeding_function_priority = "Breeding/function priority",
  conflict_quadrant = "Conflict quadrant",
  background = "Background"
)

quadrant_colors <- c(
  conservation_priority = "#B9DDF3",
  breeding_function_priority = "#BFE7CD",
  conflict_quadrant = "#F4B6B0",
  background = "#E6E6E6"
)
quadrant_border_colors <- c(
  conservation_priority = "#5DA5D6",
  breeding_function_priority = "#63B97A",
  conflict_quadrant = "#D84B3A",
  background = "#B5B5B5"
)

normalize_quadrant <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c(
    "background_low_low",
    "background_quadrant",
    "background_region",
    "background",
    "other",
    "other_windows",
    "Other windows"
  )] <- "background"
  x
}

read_quadrant_col <- function(path) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path)
  }

  if (requireNamespace("data.table", quietly = TRUE)) {
    x <- data.table::fread(
      path,
      select = "quadrant",
      data.table = FALSE,
      showProgress = FALSE
    )
    return(normalize_quadrant(x$quadrant))
  }

  if (requireNamespace("readr", quietly = TRUE)) {
    x <- readr::read_tsv(
      path,
      col_select = "quadrant",
      show_col_types = FALSE,
      progress = FALSE
    )
    return(normalize_quadrant(x$quadrant))
  }

  message("Neither data.table nor readr is installed; falling back to read.delim(), which may be slow for large files.")
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"quadrant" %in% names(x)) {
    stop("Column 'quadrant' not found in: ", path)
  }
  normalize_quadrant(x$quadrant)
}

read_quadrant_gene_cols <- function(path) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path)
  }

  if (requireNamespace("data.table", quietly = TRUE)) {
    x <- data.table::fread(
      path,
      select = c("quadrant", "genes_in_window"),
      data.table = FALSE,
      showProgress = FALSE
    )
  } else if (requireNamespace("readr", quietly = TRUE)) {
    x <- readr::read_tsv(
      path,
      col_select = c("quadrant", "genes_in_window"),
      show_col_types = FALSE,
      progress = FALSE
    )
    x <- as.data.frame(x)
  } else {
    message("Neither data.table nor readr is installed; falling back to read.delim(), which may be slow for large files.")
    x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    x <- x[, intersect(c("quadrant", "genes_in_window"), names(x)), drop = FALSE]
  }

  if (!all(c("quadrant", "genes_in_window") %in% names(x))) {
    stop("Columns 'quadrant' and/or 'genes_in_window' not found in: ", path)
  }

  x$quadrant <- normalize_quadrant(x$quadrant)
  x$genes_in_window <- as.character(x$genes_in_window)
  x
}

split_gene_list <- function(x) {
  x <- x[!is.na(x)]
  x <- x[nzchar(trimws(x))]
  if (length(x) == 0) return(character())
  genes <- unlist(strsplit(x, "[,;|]", perl = TRUE), use.names = FALSE)
  genes <- trimws(genes)
  genes <- genes[nzchar(genes)]
  genes <- genes[!genes %in% c(".", "NA", "NaN", "NULL", "None", "none")]
  unique(genes)
}

count_quadrants <- function(path, panel_label) {
  q <- read_quadrant_col(path)
  tab <- as.data.frame(table(factor(q, levels = quadrant_levels)), stringsAsFactors = FALSE)
  names(tab) <- c("quadrant", "n")

  # Keep any unexpected quadrant names instead of silently dropping them.
  extra <- setdiff(unique(q), quadrant_levels)
  if (length(extra) > 0) {
    extra_tab <- as.data.frame(table(q[q %in% extra]), stringsAsFactors = FALSE)
    names(extra_tab) <- c("quadrant", "n")
    tab <- rbind(tab, extra_tab)
    quadrant_labels[extra] <<- extra
    quadrant_colors[extra] <<- "#999999"
    quadrant_border_colors[extra] <<- "#666666"
  }

  tab$panel <- panel_label
  tab$n <- as.integer(tab$n)
  tab$total <- sum(tab$n, na.rm = TRUE)
  tab$percent <- ifelse(tab$total > 0, tab$n / tab$total * 100, 0)
  tab$quadrant_label <- ifelse(
    tab$quadrant %in% names(quadrant_labels),
    unname(quadrant_labels[tab$quadrant]),
    tab$quadrant
  )
  tab$label <- sprintf("%s\n%s (%.1f%%)", tab$quadrant_label, format(tab$n, big.mark = ","), tab$percent)
  tab
}

count_unique_genes_by_quadrant <- function(path, panel_label) {
  x <- read_quadrant_gene_cols(path)
  tab <- data.frame(
    quadrant = quadrant_levels,
    n = integer(length(quadrant_levels)),
    stringsAsFactors = FALSE
  )

  for (q in quadrant_levels) {
    genes <- split_gene_list(x$genes_in_window[x$quadrant == q])
    tab$n[tab$quadrant == q] <- length(genes)
  }

  extra <- setdiff(unique(x$quadrant), quadrant_levels)
  if (length(extra) > 0) {
    extra_tab <- data.frame(
      quadrant = extra,
      n = vapply(extra, function(q) length(split_gene_list(x$genes_in_window[x$quadrant == q])), integer(1)),
      stringsAsFactors = FALSE
    )
    tab <- rbind(tab, extra_tab)
    quadrant_labels[extra] <<- extra
    quadrant_colors[extra] <<- "#999999"
    quadrant_border_colors[extra] <<- "#666666"
  }

  tab$panel <- panel_label
  tab$total <- sum(tab$n, na.rm = TRUE)
  tab$percent <- ifelse(tab$total > 0, tab$n / tab$total * 100, 0)
  tab$quadrant_label <- ifelse(
    tab$quadrant %in% names(quadrant_labels),
    unname(quadrant_labels[tab$quadrant]),
    tab$quadrant
  )
  tab$label <- sprintf("%s\n%s (%.1f%%)", tab$quadrant_label, format(tab$n, big.mark = ","), tab$percent)
  tab
}

adjust_label_y <- function(y, min_gap = 0.17, lower = -1.20, upper = 1.20) {
  if (length(y) <= 1) return(y)
  ord <- order(y)
  yy <- y[ord]
  for (i in seq_along(yy)[-1]) {
    if (yy[i] - yy[i - 1] < min_gap) {
      yy[i] <- yy[i - 1] + min_gap
    }
  }
  overflow <- max(yy, na.rm = TRUE) - upper
  if (is.finite(overflow) && overflow > 0) yy <- yy - overflow
  underflow <- lower - min(yy, na.rm = TRUE)
  if (is.finite(underflow) && underflow > 0) yy <- yy + underflow
  out <- y
  out[ord] <- yy
  out
}

prepare_donut_layers <- function(df, n_angle = 120) {
  panels <- unique(as.character(df$panel))
  if (length(panels) > 1) {
    parts <- lapply(panels, function(pan) {
      prepare_donut_layers(df[df$panel == pan, , drop = FALSE], n_angle = n_angle)
    })
    return(list(
      poly = do.call(rbind, lapply(parts, `[[`, "poly")),
      labels = do.call(rbind, lapply(parts, `[[`, "labels"))
    ))
  }

  slice_df <- df[df$n > 0, , drop = FALSE]
  slice_df$quadrant <- factor(slice_df$quadrant, levels = quadrant_levels)
  slice_df <- slice_df[order(slice_df$quadrant), , drop = FALSE]

  total <- sum(slice_df$n, na.rm = TRUE)
  slice_df$fraction <- slice_df$n / total
  slice_df$cum0 <- c(0, head(cumsum(slice_df$fraction), -1))
  slice_df$cum1 <- cumsum(slice_df$fraction)
  slice_df$theta_start <- pi / 2 - 2 * pi * slice_df$cum0
  slice_df$theta_end <- pi / 2 - 2 * pi * slice_df$cum1
  slice_df$theta_mid <- (slice_df$theta_start + slice_df$theta_end) / 2

  outer_r <- 1.00
  inner_r <- 0.43
  label_r0 <- 1.03
  label_r1 <- 1.17
  label_x_abs <- 1.42

  polygon_df <- do.call(rbind, lapply(seq_len(nrow(slice_df)), function(i) {
    row <- slice_df[i, ]
    theta <- seq(row$theta_start, row$theta_end, length.out = n_angle)
    outer <- data.frame(
      quadrant = row$quadrant,
      quadrant_label = row$quadrant_label,
      panel = row$panel,
      x = outer_r * cos(theta),
      y = outer_r * sin(theta)
    )
    inner <- data.frame(
      quadrant = row$quadrant,
      quadrant_label = row$quadrant_label,
      panel = row$panel,
      x = inner_r * cos(rev(theta)),
      y = inner_r * sin(rev(theta))
    )
    poly <- rbind(outer, inner)
    poly$piece <- paste(row$panel, row$quadrant, sep = "__")
    poly
  }))

  label_df <- slice_df
  label_df$x0 <- label_r0 * cos(label_df$theta_mid)
  label_df$y0 <- label_r0 * sin(label_df$theta_mid)
  label_df$x1 <- label_r1 * cos(label_df$theta_mid)
  label_df$y1 <- label_r1 * sin(label_df$theta_mid)
  label_df$side <- ifelse(cos(label_df$theta_mid) >= 0, "right", "left")
  label_df$x_label <- ifelse(label_df$side == "right", label_x_abs, -label_x_abs)
  label_df$hjust <- ifelse(label_df$side == "right", 0, 1)
  label_df$y_label <- label_df$y1

  label_df <- do.call(rbind, lapply(split(label_df, list(label_df$panel, label_df$side), drop = TRUE), function(z) {
    z$y_label <- adjust_label_y(z$y_label)
    z
  }))
  label_df$x1 <- ifelse(label_df$side == "right", label_df$x_label - 0.05, label_df$x_label + 0.05)

  list(poly = polygon_df, labels = label_df)
}

all_counts <- count_quadrants(all_file, "All windows")
top1_counts <- count_quadrants(top1_file, "Top 1% windows")
all_gene_counts <- count_unique_genes_by_quadrant(all_file, "All-window genes")
top1_gene_counts <- count_unique_genes_by_quadrant(top1_file, "Top 1% window genes")

summary_df <- rbind(all_counts, top1_counts)
write.table(
  summary_df[, c("panel", "quadrant", "quadrant_label", "n", "percent")],
  file.path(outdir, paste0(prefix, "_window_counts_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

gene_summary_df <- rbind(all_gene_counts, top1_gene_counts)
write.table(
  gene_summary_df[, c("panel", "quadrant", "quadrant_label", "n", "percent")],
  file.path(outdir, paste0(prefix, "_unique_gene_counts_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

make_pie <- function(df, title, subtitle = NULL) {
  layers <- prepare_donut_layers(df)

  ggplot() +
    geom_polygon(
      data = layers$poly,
      aes(x = x, y = y, group = piece, fill = quadrant),
      color = "white",
      linewidth = 0.55
    ) +
    geom_segment(
      data = layers$labels,
      aes(x = x0, xend = x1, y = y0, yend = y_label, color = quadrant),
      linewidth = 0.35
    ) +
    geom_text(
      data = layers$labels,
      aes(x = x_label, y = y_label, label = label, hjust = hjust),
      vjust = 0.5,
      size = 3.0,
      lineheight = 0.95,
      color = "black"
    ) +
    coord_fixed(xlim = c(-1.72, 1.72), ylim = c(-1.32, 1.32), clip = "off") +
    scale_fill_manual(
      values = quadrant_colors,
      breaks = names(quadrant_colors),
      labels = quadrant_labels,
      drop = FALSE,
      name = NULL
    ) +
    scale_color_manual(values = quadrant_border_colors, guide = "none") +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_size = 11) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"),
      legend.position = "right",
      legend.text = element_text(size = 9),
      plot.margin = margin(8, 8, 8, 8)
    )
}

all_pie <- make_pie(
  all_counts,
  "Window counts across four quadrants",
  sprintf("Total windows = %s", format(sum(all_counts$n), big.mark = ","))
)

top1_pie <- make_pie(
  top1_counts,
  "Top 1% window counts across quadrants",
  sprintf("Total top 1%% windows = %s; background is shown as 0 if absent", format(sum(top1_counts$n), big.mark = ","))
)

all_gene_pie <- make_pie(
  all_gene_counts,
  "Unique gene counts across four quadrants",
  sprintf("Total quadrant-level unique genes = %s", format(sum(all_gene_counts$n), big.mark = ","))
)

top1_gene_pie <- make_pie(
  top1_gene_counts,
  "Unique gene counts in top 1% windows",
  sprintf("Total quadrant-level unique genes = %s", format(sum(top1_gene_counts$n), big.mark = ","))
)

combined_layers <- prepare_donut_layers(rbind(all_counts, top1_counts))
combined_layers$poly$panel <- factor(combined_layers$poly$panel, levels = c("All windows", "Top 1% windows"))
combined_layers$labels$panel <- factor(combined_layers$labels$panel, levels = c("All windows", "Top 1% windows"))

combined_pie <- ggplot() +
  geom_polygon(
    data = combined_layers$poly,
    aes(x = x, y = y, group = piece, fill = quadrant),
    color = "white",
    linewidth = 0.55
  ) +
  geom_segment(
    data = combined_layers$labels,
    aes(x = x0, xend = x1, y = y0, yend = y_label, color = quadrant),
    linewidth = 0.32
  ) +
  geom_text(
    data = combined_layers$labels,
    aes(x = x_label, y = y_label, label = label, hjust = hjust),
    vjust = 0.5,
    size = 2.65,
    lineheight = 0.93,
    color = "black"
  ) +
  coord_fixed(xlim = c(-1.78, 1.78), ylim = c(-1.35, 1.35), clip = "off") +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_manual(
    values = quadrant_colors,
    breaks = names(quadrant_colors),
    labels = quadrant_labels,
    drop = FALSE,
    name = NULL
  ) +
  scale_color_manual(values = quadrant_border_colors, guide = "none") +
  labs(title = "Quadrant window count distributions") +
  theme_void(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    plot.margin = margin(8, 8, 8, 8)
  )

gene_combined_layers <- prepare_donut_layers(rbind(all_gene_counts, top1_gene_counts))
gene_combined_layers$poly$panel <- factor(gene_combined_layers$poly$panel, levels = c("All-window genes", "Top 1% window genes"))
gene_combined_layers$labels$panel <- factor(gene_combined_layers$labels$panel, levels = c("All-window genes", "Top 1% window genes"))

gene_combined_pie <- ggplot() +
  geom_polygon(
    data = gene_combined_layers$poly,
    aes(x = x, y = y, group = piece, fill = quadrant),
    color = "white",
    linewidth = 0.55
  ) +
  geom_segment(
    data = gene_combined_layers$labels,
    aes(x = x0, xend = x1, y = y0, yend = y_label, color = quadrant),
    linewidth = 0.32
  ) +
  geom_text(
    data = gene_combined_layers$labels,
    aes(x = x_label, y = y_label, label = label, hjust = hjust),
    vjust = 0.5,
    size = 2.65,
    lineheight = 0.93,
    color = "black"
  ) +
  coord_fixed(xlim = c(-1.78, 1.78), ylim = c(-1.35, 1.35), clip = "off") +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_manual(
    values = quadrant_colors,
    breaks = names(quadrant_colors),
    labels = quadrant_labels,
    drop = FALSE,
    name = NULL
  ) +
  scale_color_manual(values = quadrant_border_colors, guide = "none") +
  labs(title = "Unique gene count distributions across quadrants") +
  theme_void(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    plot.margin = margin(8, 8, 8, 8)
  )

ggsave(file.path(outdir, paste0(prefix, "_all_windows_pie.pdf")), all_pie, width = 7.6, height = 5.2, device = cairo_pdf)
ggsave(file.path(outdir, paste0(prefix, "_all_windows_pie.png")), all_pie, width = 7.6, height = 5.2, dpi = 400, bg = "white")

ggsave(file.path(outdir, paste0(prefix, "_top1_windows_pie.pdf")), top1_pie, width = 7.6, height = 5.2, device = cairo_pdf)
ggsave(file.path(outdir, paste0(prefix, "_top1_windows_pie.png")), top1_pie, width = 7.6, height = 5.2, dpi = 400, bg = "white")

ggsave(file.path(outdir, paste0(prefix, "_window_count_pies.pdf")), combined_pie, width = 13.5, height = 5.6, device = cairo_pdf)
ggsave(file.path(outdir, paste0(prefix, "_window_count_pies.png")), combined_pie, width = 13.5, height = 5.6, dpi = 400, bg = "white")

ggsave(file.path(outdir, paste0(prefix, "_all_windows_unique_genes_pie.pdf")), all_gene_pie, width = 7.6, height = 5.2, device = cairo_pdf)
ggsave(file.path(outdir, paste0(prefix, "_all_windows_unique_genes_pie.png")), all_gene_pie, width = 7.6, height = 5.2, dpi = 400, bg = "white")

ggsave(file.path(outdir, paste0(prefix, "_top1_windows_unique_genes_pie.pdf")), top1_gene_pie, width = 7.6, height = 5.2, device = cairo_pdf)
ggsave(file.path(outdir, paste0(prefix, "_top1_windows_unique_genes_pie.png")), top1_gene_pie, width = 7.6, height = 5.2, dpi = 400, bg = "white")

ggsave(file.path(outdir, paste0(prefix, "_unique_gene_count_pies.pdf")), gene_combined_pie, width = 13.5, height = 5.6, device = cairo_pdf)
ggsave(file.path(outdir, paste0(prefix, "_unique_gene_count_pies.png")), gene_combined_pie, width = 13.5, height = 5.6, dpi = 400, bg = "white")

cat("\nSaved files:\n")
cat("  ", file.path(outdir, paste0(prefix, "_window_counts_summary.tsv")), "\n", sep = "")
cat("  ", file.path(outdir, paste0(prefix, "_unique_gene_counts_summary.tsv")), "\n", sep = "")
cat("  ", file.path(outdir, paste0(prefix, "_all_windows_pie.pdf")), "\n", sep = "")
cat("  ", file.path(outdir, paste0(prefix, "_all_windows_pie.png")), "\n", sep = "")
cat("  ", file.path(outdir, paste0(prefix, "_top1_windows_pie.pdf")), "\n", sep = "")
cat("  ", file.path(outdir, paste0(prefix, "_top1_windows_pie.png")), "\n", sep = "")
cat("  ", file.path(outdir, paste0(prefix, "_window_count_pies.pdf")), "\n", sep = "")
cat("  ", file.path(outdir, paste0(prefix, "_window_count_pies.png")), "\n", sep = "")
cat("  ", file.path(outdir, paste0(prefix, "_all_windows_unique_genes_pie.pdf")), "\n", sep = "")
cat("  ", file.path(outdir, paste0(prefix, "_all_windows_unique_genes_pie.png")), "\n", sep = "")
cat("  ", file.path(outdir, paste0(prefix, "_top1_windows_unique_genes_pie.pdf")), "\n", sep = "")
cat("  ", file.path(outdir, paste0(prefix, "_top1_windows_unique_genes_pie.png")), "\n", sep = "")
cat("  ", file.path(outdir, paste0(prefix, "_unique_gene_count_pies.pdf")), "\n", sep = "")
cat("  ", file.path(outdir, paste0(prefix, "_unique_gene_count_pies.png")), "\n", sep = "")
