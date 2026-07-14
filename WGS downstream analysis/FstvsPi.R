library(tidyverse)
library(patchwork)
library(scales)
library(ggrastr)

#========================
# 1. 输入文件
#========================
fst_file <- "domestic_vs_wild_fst.txt"
pi_file  <- "domestic_vs_wild_pi.txt"

#========================
# 2. 读取数据
#========================
fst_raw <- read.delim(fst_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pi_raw  <- read.delim(pi_file,  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#========================
# 3. 整理 FST
#========================
fst_dat <- fst_raw %>%
  transmute(
    CHROM = as.character(chromosome),
    start = as.numeric(window_pos_1),
    end   = as.numeric(window_pos_2),
    Fst   = as.numeric(avg_wc_fst),
    no_snps = as.numeric(no_snps)
  )

#========================
# 4. 整理 PI，长表转宽表
#========================
pi_dat <- pi_raw %>%
  transmute(
    pop   = as.character(pop),
    CHROM = as.character(chromosome),
    start = as.numeric(window_pos_1),
    end   = as.numeric(window_pos_2),
    avg_pi = as.numeric(avg_pi)
  ) %>%
  pivot_wider(
    names_from = pop,
    values_from = avg_pi
  ) %>%
  rename(
    pi_domestic = domestic,
    pi_wild = wild
  )

#========================
# 5. 合并并计算 ratio
#========================
dat_all <- fst_dat %>%
  inner_join(pi_dat, by = c("CHROM", "start", "end")) %>%
  mutate(
    pi_ratio = pi_wild / pi_domestic,
    log2_pi_ratio = log2(pi_ratio)
  ) %>%
  dplyr::filter(
    !is.na(Fst),
    !is.na(pi_domestic),
    !is.na(pi_wild),
    !is.na(pi_ratio),
    !is.na(log2_pi_ratio),
    is.finite(Fst),
    is.finite(pi_ratio),
    is.finite(log2_pi_ratio),
    pi_domestic > 0,
    pi_wild > 0
  )

#========================
# 6. 区分常染色体和性染色体
#========================
dat_all <- dat_all %>%
  mutate(
    chrom_type = case_when(
      CHROM %in% c("ChrX", "ChrY", "X", "Y", "chrX", "chrY", "Chrx", "Chry") ~ "sexchr",
      TRUE ~ "autosome"
    )
  )

#========================
# 7. 颜色与主题
#========================
col_whole <- "grey78"
col_A <- "#3f60ad"
col_B <- "#58b45e"
col_fst_high <- "#ef6c3b"

cols_main <- c(
  "Whole genome" = col_whole,
  "Selected region (A region)" = col_A,
  "Selected region (B region)" = col_B
)

theme_paper <- theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(linewidth = 0.7, colour = "black"),
    axis.ticks = element_line(linewidth = 0.6, colour = "black"),
    axis.text = element_text(colour = "black"),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

#========================
# 8. 直方图辅助函数
#========================
make_hist_df <- function(x, bins = 60) {
  h <- hist(x, breaks = bins, plot = FALSE)
  tibble(
    xmid = h$mids,
    xmin = h$breaks[-length(h$breaks)],
    xmax = h$breaks[-1],
    count = h$counts,
    freq = h$counts / sum(h$counts) * 100
  )
}

#========================
# 9. 多格式保存函数
#========================
save_plot_multi <- function(plot_obj, filename_base, width = 10.5, height = 8.2, dpi = 400) {
  ggsave(
    paste0(filename_base, ".pdf"),
    plot_obj,
    width = width,
    height = height,
    device = cairo_pdf
  )

  ggsave(
    paste0(filename_base, ".png"),
    plot_obj,
    width = width,
    height = height,
    dpi = dpi
  )

  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(
      paste0(filename_base, ".svg"),
      plot_obj,
      width = width,
      height = height,
      device = svglite::svglite
    )
  } else {
    message("Package 'svglite' is not available, skip SVG output: ", filename_base)
  }
}

#========================
# 10. 作图函数：对一个数据子集单独计算阈值、作图、导表
#========================
plot_joint_by_subset <- function(dat_sub, prefix, title_text = NULL, raster_dpi = 600) {
  
  if (nrow(dat_sub) < 10) {
    message(prefix, ": 数据行太少，跳过。")
    return(NULL)
  }
  
  # 重新计算阈值：每个子集单独算
  fst_cut <- quantile(dat_sub$Fst, 0.99, na.rm = TRUE)
  x_left  <- quantile(dat_sub$log2_pi_ratio, 0.01, na.rm = TRUE)
  x_right <- quantile(dat_sub$log2_pi_ratio, 0.99, na.rm = TRUE)
  
  # 分类
  dat_sub <- dat_sub %>%
    mutate(
      group = case_when(
        Fst >= fst_cut & log2_pi_ratio <= x_left  ~ "Selected region (A region)",
        Fst >= fst_cut & log2_pi_ratio >= x_right ~ "Selected region (B region)",
        TRUE ~ "Whole genome"
      )
    )
  
  dat_sub$group <- factor(
    dat_sub$group,
    levels = c("Whole genome", "Selected region (A region)", "Selected region (B region)")
  )
  
  # 导出候选区域
  candidate_A <- dat_sub %>% dplyr::filter(group == "Selected region (A region)")
  candidate_B <- dat_sub %>% dplyr::filter(group == "Selected region (B region)")
  
  write.table(dat_sub, paste0(prefix, "_merged.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(candidate_A, paste0(prefix, "_candidate_A_regions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(candidate_B, paste0(prefix, "_candidate_B_regions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  # 顶部边际图数据
  top_hist <- make_hist_df(dat_sub$log2_pi_ratio, bins = 80) %>%
    mutate(
      zone = case_when(
        xmid <= x_left  ~ "left",
        xmid >= x_right ~ "right",
        TRUE ~ "middle"
      ),
      cum = cumsum(freq)
    )
  
  # 右侧边际图数据
  right_hist <- make_hist_df(dat_sub$Fst, bins = 60) %>%
    mutate(
      zone = if_else(xmid >= fst_cut, "high", "low"),
      cum = cumsum(freq)
    )
  
  top_freq_max <- max(top_hist$freq, na.rm = TRUE)
  if (top_freq_max == 0) top_freq_max <- 1
  
  right_freq_max <- max(right_hist$freq, na.rm = TRUE)
  if (right_freq_max == 0) right_freq_max <- 1
  
  binw_top <- ifelse(nrow(top_hist) > 1, median(diff(top_hist$xmid)), 0.1)
  binw_right <- ifelse(nrow(right_hist) > 1, median(diff(right_hist$xmid)), 0.01)
  
  # 主图：仅背景灰点栅格化；候选点保持矢量
  p_main <- ggplot() +
    ggrastr::geom_point_rast(
      data = dat_sub %>% dplyr::filter(group == "Whole genome"),
      aes(x = log2_pi_ratio, y = Fst),
      color = col_whole,
      size = 2.0,
      alpha = 0.85,
      raster.dpi = raster_dpi
    ) +
    geom_point(
      data = dat_sub %>% dplyr::filter(group == "Selected region (A region)"),
      aes(x = log2_pi_ratio, y = Fst, color = group),
      size = 2.7,
      alpha = 0.95
    ) +
    geom_point(
      data = dat_sub %>% dplyr::filter(group == "Selected region (B region)"),
      aes(x = log2_pi_ratio, y = Fst, color = group),
      size = 2.7,
      alpha = 0.95
    ) +
    geom_vline(xintercept = x_left, linetype = "dashed", color = "grey60", linewidth = 0.9) +
    geom_vline(xintercept = x_right, linetype = "dashed", color = "grey60", linewidth = 0.9) +
    geom_hline(yintercept = fst_cut, linetype = "dashed", color = "grey60", linewidth = 0.9) +
    scale_color_manual(values = cols_main) +
    labs(
      x = expression(log[2](pi[wild] / pi[domestic])),
      y = "FST",
      title = title_text
    ) +
    annotate(
      "text",
      x = x_left,
      y = max(dat_sub$Fst, na.rm = TRUE) * 0.985,
      label = paste0("Left cutoff = ", round(x_left, 3)),
      hjust = 1.05,
      vjust = 1,
      size = 4,
      color = "grey30"
    ) +
    annotate(
      "text",
      x = x_right,
      y = max(dat_sub$Fst, na.rm = TRUE) * 0.985,
      label = paste0("Right cutoff = ", round(x_right, 3)),
      hjust = -0.05,
      vjust = 1,
      size = 4,
      color = "grey30"
    ) +
    annotate(
      "text",
      x = min(dat_sub$log2_pi_ratio, na.rm = TRUE),
      y = fst_cut,
      label = paste0("FST cutoff = ", round(fst_cut, 3)),
      hjust = 0,
      vjust = -0.6,
      size = 4,
      color = "grey30"
    ) +
    coord_cartesian(clip = "off") +
    theme_paper +
    theme(
      legend.position = c(0.64, 0.18),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # 顶部边际图
  p_top <- ggplot(top_hist, aes(x = xmid, y = freq)) +
    geom_col(
      data = top_hist %>% dplyr::filter(zone == "middle"),
      fill = "grey82",
      width = binw_top * 0.95
    ) +
    geom_col(
      data = top_hist %>% dplyr::filter(zone == "left"),
      fill = col_A,
      width = binw_top * 0.95
    ) +
    geom_col(
      data = top_hist %>% dplyr::filter(zone == "right"),
      fill = col_B,
      width = binw_top * 0.95
    ) +
    geom_line(
      aes(y = cum / 100 * top_freq_max),
      color = "black",
      linewidth = 1
    ) +
    geom_vline(xintercept = x_left, linetype = "dashed", color = "grey60", linewidth = 0.9) +
    geom_vline(xintercept = x_right, linetype = "dashed", color = "grey60", linewidth = 0.9) +
    scale_y_continuous(
      name = "Frequency (%)",
      sec.axis = sec_axis(~ . / top_freq_max * 100, name = "Cumulative (%)")
    ) +
    theme_paper +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "none"
    )
  
  # 右侧边际图
  p_right <- ggplot(right_hist, aes(x = freq, y = xmid)) +
    geom_col(
      data = right_hist %>% dplyr::filter(zone == "low"),
      fill = "grey82",
      width = binw_right * 0.95
    ) +
    geom_col(
      data = right_hist %>% dplyr::filter(zone == "high"),
      fill = col_fst_high,
      width = binw_right * 0.95
    ) +
    geom_line(
      aes(x = cum / 100 * right_freq_max),
      color = "black",
      linewidth = 1
    ) +
    geom_hline(yintercept = fst_cut, linetype = "dashed", color = "grey60", linewidth = 0.9) +
    scale_x_continuous(
      name = "Frequency (%)",
      sec.axis = sec_axis(~ . / right_freq_max * 100, name = "Cumulative (%)")
    ) +
    theme_paper +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none"
    )
  
  p_blank <- ggplot() + theme_void()
  
  final_plot <- (p_top + p_blank) / (p_main + p_right) +
    plot_layout(
      widths = c(4.6, 1.15),
      heights = c(1.25, 4.3)
    )
  
  save_plot_multi(
    plot_obj = final_plot,
    filename_base = paste0("Fig_", prefix, "_Fst_log2piRatio_joint_paper"),
    width = 10.5,
    height = 8.2,
    dpi = 400
  )
  
  cat("\n============================\n")
  cat("Group:", prefix, "\n")
  cat("Total windows:", nrow(dat_sub), "\n")
  cat("A region:", nrow(candidate_A), "\n")
  cat("B region:", nrow(candidate_B), "\n")
  cat("FST cutoff =", round(fst_cut, 6), "\n")
  cat("Left cutoff =", round(x_left, 6), "\n")
  cat("Right cutoff =", round(x_right, 6), "\n")
  cat("Raster DPI =", raster_dpi, "\n")
  cat("============================\n")
  
  return(list(
    data = dat_sub,
    candidate_A = candidate_A,
    candidate_B = candidate_B,
    fst_cut = fst_cut,
    x_left = x_left,
    x_right = x_right
  ))
}

#========================
# 11. 分别处理常染色体和性染色体
#========================
dat_autosome <- dat_all %>% dplyr::filter(chrom_type == "autosome")
dat_sexchr   <- dat_all %>% dplyr::filter(chrom_type == "sexchr")

res_autosome <- plot_joint_by_subset(
  dat_sub = dat_autosome,
  prefix = "autosome",
  title_text = "Autosomes",
  raster_dpi = 600
)

res_sexchr <- plot_joint_by_subset(
  dat_sub = dat_sexchr,
  prefix = "sexchr",
  title_text = "Sex chromosomes (ChrX + ChrY)",
  raster_dpi = 600
)

cat("\n全部分析完成。\n")
