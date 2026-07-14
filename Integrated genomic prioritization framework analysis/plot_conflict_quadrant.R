#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# Optional packages.  If unavailable, the script still runs with normal points.
has_ggrastr <- requireNamespace("ggrastr", quietly = TRUE)
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

# ================================================================
# Parameters
# ================================================================
result_dir <- "postprocess_results"
quadrant_file <- file.path(result_dir, "conflict_quadrant_table.tsv")
gene_file <- file.path(result_dir, "conflict_quadrant_top1_candidate_genes.tsv")
fallback_gene_file <- file.path(result_dir, "top_conflict_candidate_genes.tsv")

out_prefix <- file.path(result_dir, "conflict_quadrant_plot_R")

top_label_n <- 5
raster_dpi <- 500

# If TRUE, dense scatter points are rasterized while axes, labels, legend,
# text, and threshold lines remain vector objects in the PDF.
rasterize_points <- TRUE

# ================================================================
# Helpers
# ================================================================
num <- function(x) suppressWarnings(as.numeric(x))
first_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else x[1]
}

add_raster_points <- function(p, data, mapping, ..., raster = rasterize_points) {
  if (raster && has_ggrastr) {
    p + ggrastr::geom_point_rast(data = data, mapping = mapping, ..., raster.dpi = raster_dpi)
  } else {
    p + geom_point(data = data, mapping = mapping, ...)
  }
}

rebuild_quadrants_and_top1 <- function(df) {
  df <- df %>%
    mutate(
      C_score = num(C_score),
      S_score = num(S_score),
      M_score = num(M_score),
      T_score = num(T_score),
      selection_function_score = if ("selection_function_score" %in% names(.)) {
        num(selection_function_score)
      } else {
        0.4 * S_score + 0.6 * M_score
      }
    )

  qC <- quantile(df$C_score, 0.90, na.rm = TRUE)
  qY <- quantile(df$selection_function_score, 0.90, na.rm = TRUE)

  df <- df %>%
    mutate(
      quadrant = case_when(
        C_score >= qC & selection_function_score >= qY ~ "conflict_quadrant",
        C_score >= qC & selection_function_score <  qY ~ "conservation_priority",
        C_score <  qC & selection_function_score >= qY ~ "breeding_function_priority",
        TRUE ~ "background_low_low"
      ),
      quadrant_x_cutoff = qC,
      quadrant_y_cutoff = qY,
      quadrant_top1_flag = 0L,
      quadrant_top1_rank = NA_real_,
      quadrant_top1_cutoff = NA_real_,
      quadrant_top1_metric = "",
      quadrant_top1_score = case_when(
        quadrant == "conservation_priority" ~ C_score,
        quadrant == "breeding_function_priority" ~ selection_function_score,
        quadrant == "conflict_quadrant" ~ T_score,
        TRUE ~ NA_real_
      ),
      quadrant_top1_metric = case_when(
        quadrant == "conservation_priority" ~ "C_score",
        quadrant == "breeding_function_priority" ~ "selection_function_score",
        quadrant == "conflict_quadrant" ~ "T_score",
        TRUE ~ ""
      )
    )

  for (quad in c("conservation_priority", "breeding_function_priority", "conflict_quadrant")) {
    idx <- which(df$quadrant == quad & !is.na(df$quadrant_top1_score))
    if (length(idx) == 0) next

    ord <- idx[order(df$quadrant_top1_score[idx], decreasing = TRUE)]
    n_top <- max(1, ceiling(length(ord) * 0.01))
    top_idx <- ord[seq_len(n_top)]
    cutoff <- df$quadrant_top1_score[ord[n_top]]

    df$quadrant_top1_flag[top_idx] <- 1L
    df$quadrant_top1_rank[ord] <- seq_along(ord)
    df$quadrant_top1_cutoff[idx] <- cutoff
  }

  df
}

# ================================================================
# Load and prepare data
# ================================================================
if (!file.exists(quadrant_file)) {
  stop("Cannot find quadrant table: ", quadrant_file)
}

quad <- read.delim(quadrant_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
  rebuild_quadrants_and_top1()

write.table(
  quad,
  file.path(result_dir, "conflict_quadrant_table.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

quad_top1 <- quad %>% filter(quadrant_top1_flag == 1)
write.table(
  quad_top1,
  file.path(result_dir, "conflict_quadrant_top1_windows.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

qC <- unique(quad$quadrant_x_cutoff)[1]
qY <- unique(quad$quadrant_y_cutoff)[1]

plot_df <- quad %>%
  filter(!is.na(C_score), !is.na(selection_function_score))

priority_top1 <- plot_df %>%
  filter(quadrant != "background_low_low", quadrant_top1_flag == 1)

x_min <- min(min(plot_df$C_score, na.rm = TRUE), -0.5) - 0.3
x_main_max <- max(
  quantile(plot_df$C_score, 0.98, na.rm = TRUE),
  qC + 2.5,
  ifelse(nrow(priority_top1) > 0, max(priority_top1$C_score, na.rm = TRUE), NA_real_),
  na.rm = TRUE
)
x_main_max <- x_main_max + max((x_main_max - x_min) * 0.03, 0.2)

y_min <- min(plot_df$selection_function_score, na.rm = TRUE)
y_max <- max(plot_df$selection_function_score, na.rm = TRUE)
y_pad <- max((y_max - y_min) * 0.05, 0.2)
y_main_min <- y_min - y_pad
y_main_max <- y_max + y_pad

plot_df <- plot_df %>%
  filter(C_score <= x_main_max) %>%
  mutate(
    plot_group = case_when(
      quadrant_top1_flag == 1 & quadrant == "conservation_priority" ~ "Conservation priority top 1%",
      quadrant_top1_flag == 1 & quadrant == "breeding_function_priority" ~ "Breeding/function priority top 1%",
      quadrant_top1_flag == 1 & quadrant == "conflict_quadrant" ~ "Conflict quadrant top 1%",
      TRUE ~ "Other windows"
    ),
    plot_group = factor(
      plot_group,
      levels = c(
        "Other windows",
        "Conservation priority top 1%",
        "Breeding/function priority top 1%",
        "Conflict quadrant top 1%"
      )
    )
  )

other_df <- plot_df %>% filter(plot_group == "Other windows")
top_df <- plot_df %>% filter(plot_group != "Other windows")

# Threshold data
conservation_cut <- quad %>%
  filter(quadrant == "conservation_priority") %>%
  pull(quadrant_top1_cutoff) %>%
  first_or_na()

breeding_cut <- quad %>%
  filter(quadrant == "breeding_function_priority") %>%
  pull(quadrant_top1_cutoff) %>%
  first_or_na()

conflict_cut <- quad %>%
  filter(quadrant == "conflict_quadrant") %>%
  pull(quadrant_top1_cutoff) %>%
  first_or_na()

curve_df <- tibble(x = seq(x_min, x_main_max, length.out = 1600)) %>%
  mutate(y = conflict_cut / x) %>%
  filter(is.finite(y), x >= qC, y >= qY, y >= y_main_min, y <= y_main_max)

quadrant_cut_v <- tibble(
  threshold_type = "90% quadrant cutoff",
  xintercept = qC
)
quadrant_cut_h <- tibble(
  threshold_type = "90% quadrant cutoff",
  yintercept = qY
)
top1_conservation_line <- tibble(
  threshold_type = "Conservation top 1% cutoff",
  x = conservation_cut, xend = conservation_cut,
  y = y_main_min, yend = min(qY, y_main_max)
) %>% filter(!is.na(x))
top1_breeding_line <- tibble(
  threshold_type = "Breeding/function top 1% cutoff",
  x = x_min, xend = min(qC, x_main_max),
  y = breeding_cut, yend = breeding_cut
) %>% filter(!is.na(y))
curve_df <- curve_df %>%
  mutate(threshold_type = "Conflict top 1% cutoff")

# Label top genes if available.
gene_path <- if (file.exists(gene_file)) gene_file else fallback_gene_file
label_df <- tibble()
if (file.exists(gene_path)) {
  genes <- read.delim(gene_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (all(c("best_C_score", "best_selection_function_score") %in% names(genes))) {
    for (cc in c("display_name", "gene_name", "gene_std")) {
      if (!cc %in% names(genes)) genes[[cc]] <- NA_character_
    }
    label_df <- genes %>%
      mutate(
        best_C_score = num(best_C_score),
        best_selection_function_score = num(best_selection_function_score),
        label = coalesce(
          na_if(as.character(display_name), ""),
          na_if(as.character(gene_name), ""),
          na_if(as.character(gene_std), "")
        )
      ) %>%
      filter(!is.na(best_C_score), !is.na(best_selection_function_score), best_C_score <= x_main_max) %>%
      slice_head(n = top_label_n)
  }
}

# ================================================================
# Plot
# ================================================================
cols <- c(
  "Other windows" = "#C9C9C9",
  "Conservation priority top 1%" = "#2C7FB8",
  "Breeding/function priority top 1%" = "#4DAF6B",
  "Conflict quadrant top 1%" = "#D84B3A"
)

threshold_cols <- c(
  "90% quadrant cutoff" = "#8C8C8C",
  "Conservation top 1% cutoff" = "#2C7FB8",
  "Breeding/function top 1% cutoff" = "#4DAF6B",
  "Conflict top 1% cutoff" = "#D84B3A"
)
threshold_ltys <- c(
  "90% quadrant cutoff" = "dashed",
  "Conservation top 1% cutoff" = "dashed",
  "Breeding/function top 1% cutoff" = "dashed",
  "Conflict top 1% cutoff" = "dashed"
)

p <- ggplot()

p <- add_raster_points(
  p, other_df,
  aes(x = C_score, y = selection_function_score, color = plot_group),
  size = 1.15, alpha = 0.18, show.legend = FALSE
)

p <- add_raster_points(
  p, top_df,
  aes(x = C_score, y = selection_function_score, color = plot_group),
  size = 1.75, alpha = 0.88, show.legend = FALSE
)

# Vector-only legend keys for point classes.  The actual dense point layers
# above can be rasterized, but these dummy points keep the legend editable as
# vector objects in PDF/CorelDRAW/AI.
point_legend_df <- tibble(
  C_score = x_min,
  selection_function_score = y_main_min,
  plot_group = factor(names(cols), levels = names(cols))
)

p <- p +
  geom_point(
    data = point_legend_df,
    aes(x = C_score, y = selection_function_score, color = plot_group),
    size = 1.9,
    alpha = 0,
    show.legend = TRUE
  )

p <- p +
  geom_vline(data = quadrant_cut_v,
             aes(xintercept = xintercept, linetype = threshold_type),
             color = threshold_cols["90% quadrant cutoff"],
             linewidth = 0.45) +
  geom_hline(data = quadrant_cut_h,
             aes(yintercept = yintercept, linetype = threshold_type),
             color = threshold_cols["90% quadrant cutoff"],
             linewidth = 0.45) +
  scale_color_manual(
    values = cols,
    breaks = names(cols),
    name = NULL,
    drop = FALSE,
    guide = guide_legend(
      order = 1,
      override.aes = list(
        alpha = c(0.28, 0.95, 0.95, 0.95),
        size = c(1.6, 2.1, 2.1, 2.1)
      )
    )
  ) +
  scale_linetype_manual(
    values = threshold_ltys,
    breaks = names(threshold_cols),
    name = "Thresholds",
    drop = FALSE,
    guide = guide_legend(
      order = 2,
      override.aes = list(
        color = unname(threshold_cols),
        linewidth = rep(0.7, length(threshold_cols))
      )
    )
  ) +
  coord_cartesian(xlim = c(x_min, x_main_max), ylim = c(y_main_min, y_main_max), clip = "off") +
  labs(
    x = "Conservation-risk score (C score)",
    y = "Selection/function composite score [0.4×S + 0.6×M]"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.line = element_line(linewidth = 0.65, colour = "black"),
    axis.ticks = element_line(linewidth = 0.55, colour = "black"),
    axis.text = element_text(colour = "black"),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.box = "horizontal",
    legend.box.just = "right",
    legend.spacing.x = unit(10, "pt"),
    legend.background = element_rect(fill = scales::alpha("white", 0.78), colour = NA),
    legend.key = element_blank(),
    legend.text = element_text(size = 10.5),
    plot.margin = margin(8, 12, 8, 8)
  )

if (nrow(top1_conservation_line) > 0) {
  p <- p +
    geom_segment(
      data = top1_conservation_line,
      aes(x = x, xend = xend, y = y, yend = yend,
          linetype = threshold_type),
      color = threshold_cols["Conservation top 1% cutoff"],
      linewidth = 0.55
    )
}

if (nrow(top1_breeding_line) > 0) {
  p <- p +
    geom_segment(
      data = top1_breeding_line,
      aes(x = x, xend = xend, y = y, yend = yend,
          linetype = threshold_type),
      color = threshold_cols["Breeding/function top 1% cutoff"],
      linewidth = 0.55
    )
}

if (nrow(curve_df) > 0) {
  p <- p +
    geom_line(
      data = curve_df,
      aes(x = x, y = y, linetype = threshold_type),
      color = threshold_cols["Conflict top 1% cutoff"],
      linewidth = 0.55
    )
}

if (nrow(label_df) > 0) {
  p <- p +
    geom_point(
      data = label_df,
      aes(x = best_C_score, y = best_selection_function_score),
      shape = 21, fill = NA, color = "black", size = 2.9, linewidth = 0.45
    )

  if (has_ggrepel) {
    p <- p +
      ggrepel::geom_text_repel(
        data = label_df,
        aes(x = best_C_score, y = best_selection_function_score, label = label),
        size = 3.1,
        color = "black",
        min.segment.length = 0,
        segment.color = "#666666",
        segment.linewidth = 0.25,
        box.padding = 0.25,
        point.padding = 0.15,
        seed = 1
      )
  } else {
    p <- p +
      geom_text(
        data = label_df,
        aes(x = best_C_score, y = best_selection_function_score, label = label),
        size = 3.1,
        hjust = -0.05,
        vjust = -0.3,
        color = "black"
      )
  }
}

# Standard R/Cairo PDF tends to be friendlier to CorelDRAW than matplotlib PDF.
# Points are rasterized only if ggrastr is installed; all non-point layers remain vector.
ggsave(paste0(out_prefix, ".pdf"), p, width = 8.2, height = 6.6, device = cairo_pdf)
ggsave(paste0(out_prefix, ".png"), p, width = 8.2, height = 6.6, dpi = 500)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(paste0(out_prefix, ".svg"), p, width = 8.2, height = 6.6, device = svglite::svglite)
}

cat("Saved:\n")
cat("  ", paste0(out_prefix, ".pdf"), "\n", sep = "")
cat("  ", paste0(out_prefix, ".png"), "\n", sep = "")
if (file.exists(paste0(out_prefix, ".svg"))) {
  cat("  ", paste0(out_prefix, ".svg"), "\n", sep = "")
}
