rm(list = ls())

library(tidyverse)
library(UpSetR)
library(ggrastr)
library(patchwork)

# ================================================================
# Separate autosome / sex-chromosome selection scans
# ================================================================
file_fstpi <- "selection_Fst_piRatio_windows.tsv"
file_xpclr <- "selection_XPCLR_windows.tsv"
file_xpehh <- "selection_XPEHH_windows.tsv"
file_ihs   <- "selection_iHS_windows.tsv"

autosome_chr <- sprintf("Chr%02d", 1:28)
sex_chr <- c("ChrX", "ChrY")
cutoff_quantile <- 0.99
trim_for_plot <- TRUE
trim_quantile <- 0.995
# Only the dense, non-significant point layer is rasterised.  Larger values
# look smoother but make PDF files heavier; 600 is a good publication balance.
raster_dpi <- 600
min_hit_for_candidate <- 2

# Publication-friendly muted palette. Avoid red/blue and avoid over-pale
# "foggy" backgrounds: use sage/sand chromosome stripes and warm amber hits.
palette_chr <- c("#9FB79D", "#D8C79D")      # alternating chromosome windows
palette_highlight <- "#D39B2A"              # warm amber for top windows
palette_threshold <- "#8A6F32"              # muted umber threshold line
palette_upset_main <- "#7E7E73"             # soft olive-grey
palette_upset_sets <- "#A8B98E"             # muted sage
point_alpha_bg <- 0.62
point_alpha_hit <- 0.92

normalize_chr <- function(x) {
  x <- trimws(as.character(x))
  x <- str_replace(x, regex("^chromosome", ignore_case = TRUE), "")
  x <- str_replace(x, regex("^chr", ignore_case = TRUE), "")
  x <- trimws(x)
  case_when(
    toupper(x) == "X" ~ "ChrX",
    toupper(x) == "Y" ~ "ChrY",
    grepl("^[0-9]+$", x) ~ sprintf("Chr%02d", as.integer(x)),
    TRUE ~ paste0("Chr", x)
  )
}

zscore_na <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- sd(x, na.rm = TRUE); m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

cap_outlier <- function(x, q = 0.995) {
  if (all(is.na(x))) return(x)
  pmin(x, quantile(x, q, na.rm = TRUE))
}

# ================================================================
# 1. Read, harmonise, and merge input tables (done only once)
# ================================================================
fstpi <- read.delim(file_fstpi, check.names = FALSE) %>% mutate(start = start + 1)
xpclr <- read.delim(file_xpclr, check.names = FALSE)
xpehh <- read.delim(file_xpehh, check.names = FALSE)
ihs <- read.delim(file_ihs, check.names = FALSE)

fstpi2 <- fstpi %>% select(CHROM, start, end, pi_domestic, pi_wild, pi_ratio, log2_pi_ratio, Fst)
xpclr2 <- xpclr %>% select(CHROM, start, end, XPCLR, XPCLR_NORM, POS, nSNPs, nSNPs_avail)
xpehh2 <- xpehh %>%
  select(CHROM, start, end, XPEHH_mean, XPEHH_mean_abs, XPEHH_max, XPEHH_min,
         XPEHH_max_abs, n_snp, n_crit, crit_frac) %>%
  rename(XPEHH_n_snp = n_snp, XPEHH_n_crit = n_crit, XPEHH_crit_frac = crit_frac)
ihs2 <- ihs %>%
  select(CHROM, start, end, iHS_mean, iHS_mean_abs, iHS_max, iHS_min,
         iHS_max_abs, n_snp, n_crit, crit_frac) %>%
  rename(iHS_n_snp = n_snp, iHS_n_crit = n_crit, iHS_crit_frac = crit_frac)

for (nm in c("fstpi2", "xpclr2", "xpehh2", "ihs2")) {
  x <- get(nm); x$CHROM <- normalize_chr(x$CHROM); assign(nm, x)
}

sel_all <- fstpi2 %>%
  full_join(xpclr2, by = c("CHROM", "start", "end")) %>%
  full_join(xpehh2, by = c("CHROM", "start", "end")) %>%
  full_join(ihs2, by = c("CHROM", "start", "end")) %>%
  mutate(start = as.numeric(start), end = as.numeric(end), mid = (start + end) / 2)

cat("Merged chromosome labels:\n")
print(sort(unique(sel_all$CHROM)))

# ================================================================
# 2. Analyse one chromosome class independently
# Each invocation calculates its OWN quantiles and z scores.
# ================================================================
run_selection_scan <- function(sel_input, chr_order, group_name) {
  out_prefix <- paste0("selection_", group_name)
  sel <- sel_input %>%
    filter(CHROM %in% chr_order) %>%
    mutate(CHROM = factor(CHROM, levels = chr_order, ordered = TRUE)) %>%
    arrange(CHROM, start, end)

  if (nrow(sel) == 0) stop("No windows found for ", group_name)
  cat("\n", group_name, " chromosomes: ", paste(unique(as.character(sel$CHROM)), collapse = ", "), "\n", sep = "")
  cat(group_name, " windows: ", nrow(sel), "\n", sep = "")

  # Thresholds are computed only within this chromosome class.
  fst_cut <- quantile(sel$Fst, cutoff_quantile, na.rm = TRUE)
  # Fst/Pi is one joint signal: high Fst plus either tail of log2(pi ratio).
  pi_left_cut <- quantile(sel$log2_pi_ratio, 1 - cutoff_quantile, na.rm = TRUE)
  pi_right_cut <- quantile(sel$log2_pi_ratio, cutoff_quantile, na.rm = TRUE)
  xpclr_cut <- quantile(sel$XPCLR, cutoff_quantile, na.rm = TRUE)
  xpehh_cut <- quantile(sel$XPEHH_mean_abs, cutoff_quantile, na.rm = TRUE)
  ihs_cut <- quantile(sel$iHS_mean_abs, cutoff_quantile, na.rm = TRUE)

  sel <- sel %>%
    mutate(
      hit_FstPi = if_else(
        !is.na(Fst) & !is.na(log2_pi_ratio) & Fst >= fst_cut &
          (log2_pi_ratio <= pi_left_cut | log2_pi_ratio >= pi_right_cut), 1L, 0L
      ),
      hit_XPCLR = if_else(!is.na(XPCLR) & XPCLR >= xpclr_cut, 1L, 0L),
      hit_XPEHH = if_else(!is.na(XPEHH_mean_abs) & XPEHH_mean_abs >= xpehh_cut, 1L, 0L),
      hit_iHS = if_else(!is.na(iHS_mean_abs) & iHS_mean_abs >= ihs_cut, 1L, 0L),
      hit_sum = hit_FstPi + hit_XPCLR + hit_XPEHH + hit_iHS,
      selected_candidate = if_else(hit_sum >= min_hit_for_candidate, 1L, 0L),
      # One equally weighted Fst/Pi component: high Fst plus extreme Pi in
      # either direction.  This avoids counting Fst and Pi as two signals.
      z_FstPi = (zscore_na(Fst) + zscore_na(abs(log2_pi_ratio))) / 2,
      z_XPCLR = zscore_na(XPCLR),
      z_XPEHH_mean_abs = zscore_na(XPEHH_mean_abs),
      z_iHS_mean_abs = zscore_na(iHS_mean_abs)
    ) %>%
    mutate(integrated_score = rowSums(cbind(z_FstPi, z_XPCLR,
                                             z_XPEHH_mean_abs, z_iHS_mean_abs), na.rm = TRUE))

  score_cut <- quantile(sel$integrated_score, cutoff_quantile, na.rm = TRUE)
  sel <- sel %>% mutate(top1_integrated_score = if_else(integrated_score >= score_cut, 1L, 0L))

  threshold_df <- tibble(
    chromosome_group = group_name,
    statistic = c("Fst", "log2_pi_ratio_left", "log2_pi_ratio_right", "XPCLR", "XPEHH_mean_abs", "iHS_mean_abs", "integrated_score"),
    cutoff = c(fst_cut, pi_left_cut, pi_right_cut, xpclr_cut, xpehh_cut, ihs_cut, score_cut)
  )

  # Export separate results.
  write.table(sel, paste0(out_prefix, "_integrated_windows.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(filter(sel, selected_candidate == 1), paste0(out_prefix, "_candidates_by_hits.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(filter(sel, top1_integrated_score == 1), paste0(out_prefix, "_candidates_by_score.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(filter(sel, selected_candidate == 1 & top1_integrated_score == 1),
              paste0(out_prefix, "_candidates_highconf.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(threshold_df, paste0(out_prefix, "_thresholds.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  # Cumulative coordinates are rebuilt for this class only.
  chr_sizes <- sel %>% group_by(CHROM) %>%
    summarise(chr_len = as.numeric(max(end, na.rm = TRUE)), .groups = "drop") %>%
    arrange(CHROM) %>%
    mutate(tot = cumsum(chr_len) - chr_len, center = tot + chr_len / 2)
  sel_plot <- sel %>% left_join(select(chr_sizes, CHROM, tot), by = "CHROM") %>%
    mutate(pos_cum = mid + tot)
  axisdf <- select(chr_sizes, CHROM, center)
  chr_colours <- setNames(rep(palette_chr, length.out = length(chr_order)), chr_order)

  sel_plot2 <- sel_plot %>% mutate(
    integrated_score_plot = if (trim_for_plot) cap_outlier(integrated_score, trim_quantile) else integrated_score,
    highlight_score = if_else(integrated_score >= score_cut, "Top1%", "Background")
  )
  # Background points are rasterised; highlighted hits, threshold line, labels and
  # axes remain vector objects in PDF/SVG-aware editors.
  p1 <- ggplot(sel_plot2, aes(pos_cum, integrated_score_plot)) +
    ggrastr::geom_point_rast(
      data = filter(sel_plot2, highlight_score == "Background"),
      aes(color = CHROM), size = 0.7, alpha = point_alpha_bg, raster.dpi = raster_dpi
    ) +
    geom_point(data = filter(sel_plot2, highlight_score == "Top1%"), color = palette_highlight, size = 0.9, alpha = point_alpha_hit) +
    geom_hline(yintercept = score_cut, linetype = "dashed", color = palette_threshold, linewidth = 0.5) +
    scale_x_continuous(breaks = axisdf$center, labels = as.character(axisdf$CHROM), expand = expansion(mult = c(0.01, 0.01))) +
    scale_color_manual(values = chr_colours, drop = FALSE) +
    labs(x = "Chromosome", y = "Integrated score", title = paste0("Integrated selection scan: ", group_name)) +
    theme_bw() + theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
                       axis.text.x = element_text(angle = 45, hjust = 1, size = 8), plot.title = element_text(hjust = 0.5, face = "bold"))

  # Pi remains part of the joint Fst/Pi scoring above, but is not displayed in
  # the Manhattan panel to keep the final figure compact and easier to read.
  cut_df <- tibble(stat = c("Fst", "XPCLR", "XPEHH_mean_abs", "iHS_mean_abs"),
                   cutoff = c(fst_cut, xpclr_cut, xpehh_cut, ihs_cut))
  stat_levels <- c("Fst", "XPCLR", "XPEHH_mean_abs", "iHS_mean_abs")
  stat_labels <- c("Fst", "XPCLR", "XPEHH", "IHS")
  plot_long <- sel_plot %>%
    select(CHROM, start, end, mid, pos_cum, Fst, XPCLR, XPEHH_mean_abs, iHS_mean_abs) %>%
    pivot_longer(c(Fst, XPCLR, XPEHH_mean_abs, iHS_mean_abs), names_to = "stat", values_to = "value") %>%
    group_by(stat) %>% mutate(value_plot = if (trim_for_plot) cap_outlier(value, trim_quantile) else value) %>%
    ungroup() %>%
    mutate(
      highlight = case_when(
        stat == "Fst" & !is.na(value) & value >= fst_cut ~ "Top1%",
        stat == "XPCLR" & !is.na(value) & value >= xpclr_cut ~ "Top1%",
        stat == "XPEHH_mean_abs" & !is.na(value) & value >= xpehh_cut ~ "Top1%",
        stat == "iHS_mean_abs" & !is.na(value) & value >= ihs_cut ~ "Top1%",
        TRUE ~ "Background"
      ),
      stat = factor(stat,
                    levels = stat_levels,
                    labels = stat_labels)
    )
  # geom_hline() also participates in faceting.  Give it exactly the same
  # factor labels as the point data to prevent duplicate/empty facets.
  cut_df_plot <- cut_df %>%
    mutate(stat = factor(stat, levels = stat_levels, labels = stat_labels))
  write.table(plot_long, paste0(out_prefix, "_long_for_plot.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  p2 <- ggplot(plot_long, aes(pos_cum, value_plot)) +
    ggrastr::geom_point_rast(
      data = filter(plot_long, highlight == "Background"),
      aes(color = CHROM), size = 0.45, alpha = point_alpha_bg, raster.dpi = raster_dpi
    ) +
    geom_point(data = filter(plot_long, highlight == "Top1%"), color = palette_highlight, size = 0.5, alpha = point_alpha_hit) +
    geom_hline(data = cut_df_plot, aes(yintercept = cutoff), linetype = "dashed", color = palette_threshold, linewidth = 0.4) +
    scale_x_continuous(breaks = axisdf$center, labels = as.character(axisdf$CHROM), expand = expansion(mult = c(0.01, 0.01))) +
    scale_color_manual(values = chr_colours, drop = FALSE) + facet_wrap(~stat, scales = "free_y", ncol = 1) +
    labs(x = "Chromosome", y = "Statistic value", title = paste0("Selection statistics: ", group_name)) +
    theme_bw() + theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
                       strip.background = element_rect(fill = "white"), strip.text = element_text(face = "bold"),
                       axis.text.x = element_text(angle = 45, hjust = 1, size = 8), plot.title = element_text(hjust = 0.5, face = "bold"))

  summary_df <- tibble(chromosome_group = group_name, chromosomes = paste(as.character(axisdf$CHROM), collapse = ","),
                       total_windows = nrow(sel), candidate_by_hits = sum(sel$selected_candidate == 1),
                       candidate_by_score = sum(sel$top1_integrated_score == 1),
                       highconf_candidate = sum(sel$selected_candidate == 1 & sel$top1_integrated_score == 1))
  write.table(summary_df, paste0(out_prefix, "_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  # The hit table is returned for one combined UpSet plot after BOTH groups
  # have been scored using their own group-specific thresholds.
  upset_df <- sel %>%
    transmute(chromosome_group = group_name,
              hit_FstPi, hit_XPCLR, hit_XPEHH, hit_iHS)
  return(list(thresholds = threshold_df, summary = summary_df, upset = upset_df,
              integrated_plot = p1, each_stat_plot = p2))
}

# ================================================================
# 3. Run the two independent analyses
# ================================================================
autosome_result <- run_selection_scan(sel_all, autosome_chr, "autosome")
sexchr_result <- run_selection_scan(sel_all, sex_chr, "sexchr")

# ================================================================
# 4. One combined Manhattan figure
# Left: autosomes; right: sex chromosomes.  The six-to-one column ratio
# preserves the requested narrow sex-chromosome display without stretching it.
# ================================================================
# Put Fst, XPCLR, XPEHH and IHS (top to bottom) above the integrated selection
# scan, which is always the bottom panel. Pi is retained in the Fst/Pi scoring
# but omitted from the Manhattan display.
autosome_column <- (autosome_result$each_stat_plot / autosome_result$integrated_plot) +
  plot_layout(heights = c(10.5, 5))
sexchr_column <- (sexchr_result$each_stat_plot / sexchr_result$integrated_plot) +
  plot_layout(heights = c(10.5, 5))
combined_manhattan <- (autosome_column | sexchr_column) +
  plot_layout(widths = c(6, 1))

ggsave("selection_autosome_sexchr_manhattan_combined.pdf",
       combined_manhattan, width = 18.67, height = 15.5, device = cairo_pdf)
ggsave("selection_autosome_sexchr_manhattan_combined.png",
       combined_manhattan, width = 18.67, height = 15.5, dpi = 300)

# ================================================================
# 5. One combined UpSet plot
# The rows from both chromosome classes are combined here only after the
# threshold calls above; each row therefore retains its own class-specific
# threshold decision.
# ================================================================
upset_combined <- bind_rows(autosome_result$upset, sexchr_result$upset)
upset_input <- upset_combined %>%
  select(hit_FstPi, hit_XPCLR, hit_XPEHH, hit_iHS) %>%
  as.data.frame()

pdf("selection_upset.pdf", width = 9, height = 6)
upset(upset_input,
      sets = names(upset_input),
      keep.order = TRUE,
      order.by = "freq",
      main.bar.color = palette_upset_main,
      sets.bar.color = palette_upset_sets,
      text.scale = 1.2)
dev.off()

png("selection_upset.png", width = 2700, height = 1800, res = 300)
upset(upset_input,
      sets = names(upset_input),
      keep.order = TRUE,
      order.by = "freq",
      main.bar.color = palette_upset_main,
      sets.bar.color = palette_upset_sets,
      text.scale = 1.2)
dev.off()

# Retain the chromosome-group column for downstream inspection of which
# windows contributed to the combined plot.
write.table(upset_combined, "selection_upset_combined_input.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(bind_rows(autosome_result$thresholds, sexchr_result$thresholds),
            "selection_thresholds_autosome_vs_sexchr.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bind_rows(autosome_result$summary, sexchr_result$summary),
            "selection_summary_autosome_vs_sexchr.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("All separate analyses completed.\n")
