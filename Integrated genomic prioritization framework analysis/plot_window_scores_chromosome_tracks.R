#!/usr/bin/env Rscript

# ================================================================
# Chromosome-by-chromosome window score tracks
# ------------------------------------------------
# A GenomeSyn-like single-genome track figure:
#   - one chromosome per row
#   - x-axis in Mb for each chromosome
#   - C/S/M scores as line tracks above the chromosome body
#   - T score as heatmap blocks inside the chromosome body
#
# Input:
#   conflict_loci_integrated_table.tsv
#   or postprocess_results/conflict_quadrant_table.tsv
#
# Output:
#   FMdeer_chromosome_score_tracks.pdf/png/svg
# ================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

has_ggrastr <- requireNamespace("ggrastr", quietly = TRUE)
if (!has_ggrastr) {
  message("Package 'ggrastr' is not installed; T score heatmap will remain vector. Install ggrastr to rasterize only this layer.")
}

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0 || hit == length(args)) return(default)
  args[hit + 1]
}

input_file <- get_arg("--input", "conflict_loci_integrated_table.tsv")
out_prefix <- get_arg("--out", "FMdeer_chromosome_score_tracks")

if (!file.exists(input_file) && file.exists(file.path("postprocess_results", "conflict_quadrant_table.tsv"))) {
  input_file <- file.path("postprocess_results", "conflict_quadrant_table.tsv")
}

chrom_levels <- c(sprintf("Chr%02d", 1:28), "ChrX", "ChrY")
score_cols <- c("C_score", "S_score", "M_score", "T_score")
line_scores <- c("C_score", "S_score", "M_score")

# Soft pastel palette
line_cols <- c(
  C_score = "#8DBFE3",
  S_score = "#E9A59A",
  M_score = "#9FD6B2"
)

score_labels <- c(
  C_score = "C",
  S_score = "S",
  M_score = "M",
  T_score = "T"
)

# Display settings
chrom_body_ymin <- 0.00
chrom_body_ymax <- 0.22
line_base_y <- c(C_score = 0.36, S_score = 0.61, M_score = 0.86)
line_track_height <- 0.18
row_gap <- 1.05

tick_step_mb <- 10
show_top_axis <- TRUE

rasterize_tscore <- TRUE
raster_dpi <- 500

pdf_width <- 16
pdf_height <- 16
png_dpi <- 500

# Cap line tracks by quantile for readability.
trim_for_plot <- TRUE
trim_lower_quantile <- 0.01
trim_upper_quantile <- 0.99

# Priority quadrant top1% windows from postprocess_conflict_loci.py.
# These will be highlighted as translucent blocks over the chromosome body.
highlight_priority_top1 <- TRUE
priority_top1_file <- file.path("postprocess_results", "conflict_quadrant_top1_windows.tsv")

priority_highlight_cols <- c(
  conservation_priority = "#B9DDF3",
  breeding_function_priority = "#BFE7CD",
  conflict_quadrant = "#F4B6B0"
)
priority_border_cols <- c(
  conservation_priority = "#5DA5D6",
  breeding_function_priority = "#63B97A",
  conflict_quadrant = "#D84B3A"
)
priority_labels <- c(
  conservation_priority = "Conservation priority top 1%",
  breeding_function_priority = "Breeding/function priority top 1%",
  conflict_quadrant = "Conflict quadrant top 1%"
)

# ================================================================
# Helpers
# ================================================================

num <- function(x) suppressWarnings(as.numeric(x))

normalize_chr <- function(x) {
  x <- trimws(as.character(x))
  x <- stringr::str_replace(x, stringr::regex("^chromosome", ignore_case = TRUE), "")
  x <- stringr::str_replace(x, stringr::regex("^chr", ignore_case = TRUE), "")
  x <- trimws(x)
  dplyr::case_when(
    toupper(x) == "X" ~ "ChrX",
    toupper(x) == "Y" ~ "ChrY",
    grepl("^[0-9]+$", x) ~ sprintf("Chr%02d", as.integer(x)),
    TRUE ~ paste0("Chr", x)
  )
}

rescale_01 <- function(x) {
  x <- num(x)
  if (all(is.na(x))) return(x)
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
    return(rep(0.5, length(x)))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}

cap_quantile <- function(x, lo = 0.01, hi = 0.99) {
  x <- num(x)
  if (all(is.na(x))) return(x)
  qlo <- quantile(x, lo, na.rm = TRUE)
  qhi <- quantile(x, hi, na.rm = TRUE)
  pmin(pmax(x, qlo), qhi)
}

geom_tscore_rects <- function(data, mapping, ...) {
  layer <- geom_rect(data = data, mapping = mapping, ...)
  if (rasterize_tscore && has_ggrastr) {
    ggrastr::rasterise(layer, dpi = raster_dpi)
  } else {
    layer
  }
}

# ================================================================
# Read data
# ================================================================

if (!file.exists(input_file)) {
  stop("Cannot find input file: ", input_file)
}

dat <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)

if (!all(c("CHROM", "start", "end") %in% names(dat))) {
  stop("Input file must contain CHROM, start and end columns.")
}
missing_scores <- setdiff(score_cols, names(dat))
if (length(missing_scores) > 0) {
  stop("Missing score columns: ", paste(missing_scores, collapse = ", "))
}

dat <- dat %>%
  mutate(
    CHROM = normalize_chr(CHROM),
    CHROM = factor(CHROM, levels = chrom_levels, ordered = TRUE),
    start = num(start),
    end = num(end),
    mid = (start + end) / 2,
    mid_mb = mid / 1e6,
    start_mb = start / 1e6,
    end_mb = end / 1e6
  ) %>%
  filter(!is.na(CHROM), !is.na(start), !is.na(end), start <= end) %>%
  arrange(CHROM, start, end)

if (nrow(dat) == 0) stop("No valid windows after filtering.")

chr_info <- dat %>%
  group_by(CHROM) %>%
  summarise(chr_len = max(end, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHROM) %>%
  mutate(
    chr_len_mb = chr_len / 1e6,
    row_id = row_number(),
    y0 = (n() - row_id) * row_gap,
    chr_label = as.character(CHROM)
  )

dat <- dat %>%
  left_join(select(chr_info, CHROM, chr_len_mb, y0), by = "CHROM")

cat("Input: ", input_file, "\n", sep = "")
cat("Windows: ", nrow(dat), "\n", sep = "")
cat("Chromosomes: ", nrow(chr_info), "\n", sep = "")

# ================================================================
# Prepare line tracks and heatmap blocks
# ================================================================

line_long <- dat %>%
  select(CHROM, mid_mb, y0, all_of(line_scores)) %>%
  pivot_longer(all_of(line_scores), names_to = "score_type", values_to = "score_value") %>%
  mutate(score_value = num(score_value)) %>%
  filter(!is.na(score_value)) %>%
  group_by(score_type) %>%
  mutate(
    score_plot = if (trim_for_plot) cap_quantile(score_value, trim_lower_quantile, trim_upper_quantile) else score_value,
    score_scaled = rescale_01(score_plot),
    y = y0 + line_base_y[score_type] + score_scaled * line_track_height
  ) %>%
  ungroup()

heat_df <- dat %>%
  transmute(
    CHROM, start_mb, end_mb, y0,
    T_score = num(T_score)
  ) %>%
  filter(!is.na(T_score)) %>%
  group_by(CHROM) %>%
  mutate(
    T_scaled = rescale_01(T_score)
  ) %>%
  ungroup()

priority_top1 <- tibble()
if (highlight_priority_top1 && file.exists(priority_top1_file)) {
  priority_top1_raw <- read.delim(priority_top1_file, check.names = FALSE, stringsAsFactors = FALSE)
  quadrant_col <- intersect(c("quadrant", "priority_group", "region", "class"), names(priority_top1_raw))[1]
  if (all(c("CHROM", "start", "end") %in% names(priority_top1_raw)) && !is.na(quadrant_col)) {
    priority_top1 <- priority_top1_raw %>%
      mutate(
        CHROM = normalize_chr(CHROM),
        CHROM = factor(CHROM, levels = chrom_levels, ordered = TRUE),
        start = num(start),
        end = num(end),
        quadrant = as.character(.data[[quadrant_col]])
      ) %>%
      filter(
        !is.na(CHROM), !is.na(start), !is.na(end), start <= end,
        quadrant %in% names(priority_highlight_cols)
      ) %>%
      left_join(select(chr_info, CHROM, y0), by = "CHROM") %>%
      mutate(
        start_mb = start / 1e6,
        end_mb = end / 1e6,
        priority_label = priority_labels[quadrant]
      ) %>%
      filter(!is.na(y0))
  } else {
    warning("Priority top1 file exists but lacks CHROM/start/end plus a quadrant-like column: ", priority_top1_file)
  }
}

if (nrow(priority_top1) > 0) {
  cat("Priority top1 windows highlighted: ", nrow(priority_top1), "\n", sep = "")
  print(priority_top1 %>% count(quadrant, name = "n"))
} else {
  cat("No priority top1 window highlight file found or no priority top1 rows used.\n")
}

chrom_body <- chr_info %>%
  transmute(
    CHROM,
    xmin = 0,
    xmax = chr_len_mb,
    ymin = y0 + chrom_body_ymin,
    ymax = y0 + chrom_body_ymax,
    ymid = (ymin + ymax) / 2,
    label_x = -max(chr_info$chr_len_mb, na.rm = TRUE) * 0.018,
    label = chr_label
  )

line_baselines <- tidyr::crossing(
  chr_info %>% select(CHROM, chr_len_mb, y0),
  score_type = line_scores
) %>%
  mutate(
    y = y0 + line_base_y[score_type],
    score_label = score_labels[score_type]
  )

max_mb <- max(chr_info$chr_len_mb, na.rm = TRUE)
axis_ticks <- tibble(
  x = seq(0, ceiling(max_mb / tick_step_mb) * tick_step_mb, by = tick_step_mb)
)

top_y <- max(chr_info$y0, na.rm = TRUE) + row_gap * 1.05

# ================================================================
# Plot
# ================================================================

p <- ggplot()

# top Mb axis
if (show_top_axis) {
  p <- p +
    geom_segment(aes(x = 0, xend = max_mb, y = top_y, yend = top_y),
                 linewidth = 0.35, color = "black") +
    geom_segment(data = axis_ticks,
                 aes(x = x, xend = x, y = top_y - 0.025, yend = top_y + 0.025),
                 linewidth = 0.25, color = "black") +
    geom_text(data = axis_ticks,
              aes(x = x, y = top_y + 0.055, label = x),
              size = 2.6, vjust = 0) +
    annotate("text", x = max_mb + max_mb * 0.025, y = top_y + 0.055,
             label = "Mb", hjust = 0, vjust = 0, size = 2.8)
}

# score line baselines
p <- p +
  geom_segment(
    data = line_baselines,
    aes(x = 0, xend = chr_len_mb, y = y, yend = y, color = score_type),
    linewidth = 0.22,
    alpha = 0.45,
    show.legend = FALSE
  )

# score line tracks
p <- p +
  geom_line(
    data = line_long,
    aes(x = mid_mb, y = y, color = score_type, group = interaction(CHROM, score_type)),
    linewidth = 0.18,
    alpha = 0.95,
    lineend = "round",
    show.legend = FALSE
  )

# chromosome body background; the visible border is redrawn after the
# T-score and priority-highlight layers so it will not be covered.
p <- p +
  geom_rect(
    data = chrom_body,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "white",
    color = NA,
    linewidth = 0
  )

# T-score heatmap inside chromosome body
p <- p +
  geom_tscore_rects(
    data = heat_df,
    aes(
      xmin = start_mb,
      xmax = end_mb,
      ymin = y0 + chrom_body_ymin,
      ymax = y0 + chrom_body_ymax,
      fill = T_score
    ),
    color = NA,
    alpha = 0.95
  ) +
  scale_fill_gradientn(
    colors = c("#FAFBFC", "#DDEEF6", "#AED4E6", "#C7B7E2", "#9E8AC8"),
    name = "T score"
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      barheight = grid::unit(26, "mm"),
      barwidth = grid::unit(4, "mm")
    )
  )

# Priority top1% windows, highlighted on top of T-score heatmap.
# Use one constant-color layer per quadrant to avoid interfering with
# the T-score heatmap fill scale.
if (nrow(priority_top1) > 0) {
  for (q in names(priority_highlight_cols)) {
    q_df <- priority_top1 %>% filter(quadrant == q)
    if (nrow(q_df) == 0) next
    p <- p +
      geom_rect(
        data = q_df,
        aes(
          xmin = start_mb,
          xmax = end_mb,
          ymin = y0 + chrom_body_ymin - 0.014,
          ymax = y0 + chrom_body_ymax + 0.014
        ),
        inherit.aes = FALSE,
        fill = priority_highlight_cols[[q]],
        color = priority_border_cols[[q]],
        linewidth = 0.16,
        alpha = 0.62
      )
  }
}

# chromosome body border, drawn above T-score and priority-highlight layers
p <- p +
  geom_rect(
    data = chrom_body,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = NA,
    color = "grey15",
    linewidth = 0.38
  )

# chromosome labels
p <- p +
  geom_text(
    data = chrom_body,
    aes(x = label_x, y = ymid, label = label),
    hjust = 1,
    vjust = 0.5,
    size = 2.5,
    fontface = "bold"
  )

# score legend as vector dummy segments, placed on right
legend_x0 <- max_mb * 1.08
legend_x1 <- max_mb * 1.20
legend_text_x <- max_mb * 1.23

legend_df <- tibble(
  score_type = line_scores,
  score_label = c("C score", "S score", "M score"),
  x = legend_x0,
  xend = max_mb * 1.16,
  y = top_y - c(0.25, 0.43, 0.61),
  yend = y
)

p <- p +
  annotate(
    "text",
    x = legend_x0,
    y = top_y - 0.08,
    label = "Line tracks",
    hjust = 0,
    vjust = 0.5,
    size = 3.0,
    fontface = "bold"
  ) +
  geom_segment(
    data = legend_df,
    aes(x = x, xend = xend, y = y, yend = yend, color = score_type),
    linewidth = 0.55
  ) +
  geom_text(
    data = legend_df,
    aes(x = legend_text_x, y = y, label = score_label),
    hjust = 0,
    vjust = 0.5,
    size = 2.8
  )

if (nrow(priority_top1) > 0) {
  present_priority <- names(priority_highlight_cols)[names(priority_highlight_cols) %in% unique(priority_top1$quadrant)]
  p <- p +
    annotate(
      "text",
      x = legend_x0,
      y = top_y - 1.05,
      label = "Top 1% windows",
      hjust = 0,
      vjust = 0.5,
      size = 3.0,
      fontface = "bold"
    )
  for (i in seq_along(present_priority)) {
    q <- present_priority[[i]]
    y_mid <- top_y - 1.32 - (i - 1) * 0.32
    p <- p +
      geom_rect(
        aes(
          xmin = legend_x0,
          xmax = legend_x1,
          ymin = y_mid - 0.085,
          ymax = y_mid + 0.085
        ),
        inherit.aes = FALSE,
        fill = priority_highlight_cols[[q]],
        color = priority_border_cols[[q]],
        linewidth = 0.35,
        alpha = 0.85
      ) +
      annotate(
        "text",
        x = legend_text_x,
        y = y_mid,
        label = priority_labels[[q]],
        hjust = 0,
        vjust = 0.5,
        size = 2.8
      )
  }
}

p <- p +
  scale_color_manual(values = line_cols, guide = "none")

p <- p +
  coord_cartesian(
    xlim = c(-max_mb * 0.04, max_mb * 1.82),
    ylim = c(-0.25, top_y + 0.18),
    clip = "off"
  ) +
  labs(x = NULL, y = NULL, title = "Window scores across chromosomes") +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 6)),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(8, 70, 8, 18)
  )

ggsave(paste0(out_prefix, ".pdf"), p, width = pdf_width, height = pdf_height, device = cairo_pdf)
ggsave(paste0(out_prefix, ".png"), p, width = pdf_width, height = pdf_height, dpi = png_dpi)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(paste0(out_prefix, ".svg"), p, width = pdf_width, height = pdf_height, device = svglite::svglite)
}

write.table(line_long, paste0(out_prefix, "_line_tracks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(heat_df, paste0(out_prefix, "_T_heatmap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nSaved files:\n")
cat("  ", paste0(out_prefix, ".pdf"), "\n", sep = "")
cat("  ", paste0(out_prefix, ".png"), "\n", sep = "")
if (file.exists(paste0(out_prefix, ".svg"))) {
  cat("  ", paste0(out_prefix, ".svg"), "\n", sep = "")
}
cat("  ", paste0(out_prefix, "_line_tracks.tsv"), "\n", sep = "")
cat("  ", paste0(out_prefix, "_T_heatmap.tsv"), "\n", sep = "")
