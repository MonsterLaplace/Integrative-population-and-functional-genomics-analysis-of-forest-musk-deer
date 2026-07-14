#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(gridExtra)
  library(CMplot)
})

# =========================================================
# 1. 输入文件
# =========================================================
roh_file    <- "all_samples.roh.hom"
pop_file    <- "pop.info"
genome_file <- "genome.bed"

# 基因组总大小（bp）
genome_size_bp <- 2847697564

# 滑窗参数
window_size <- 50000
step_size   <- 10000

# 输出前缀
prefix <- "FMdeer_ROH"

# =========================================================
# 2. 参数设置
# =========================================================
# 只保留两类，移除 >2.5 Mb
breaks_bp <- c(0, 200e3, 2.5e6)
labels_bp <- c("<200 kb", "200 kb–2.5 Mb")
max_roh_bp <- 2.5e6

fill_cols  <- c("Domestic" = "#E9786B", "Wild" = "#18B7C4")
point_cols <- c("Domestic" = "#E9786B", "Wild" = "#18B7C4")

point_alpha <- 0.75
point_size  <- 2.5

# CMplot 染色体交替颜色
manhattan_cols <- c(
  "#4E79A7", "#A0CBE8",
  "#F28E2B", "#FFBE7D",
  "#59A14F", "#8CD17D",
  "#E15759", "#FF9D9A",
  "#B07AA1", "#D4A6C8",
  "#76B7B2", "#9CDED9"
)

cmplot_cex <- 0.35

# =========================================================
# 3. 读取数据
# =========================================================
cat("[INFO] Reading files...\n")

pop <- read.table(pop_file, header = FALSE, sep = "", stringsAsFactors = FALSE)
if (ncol(pop) < 2) stop("pop.info.tsv 至少需要两列：IID 和 pop")
pop <- pop[, 1:2]
colnames(pop) <- c("IID", "pop")
pop$IID <- as.character(pop$IID)
pop$pop[pop$pop == "domestic"] <- "Domestic"
pop$pop[pop$pop == "wild"] <- "Wild"

genome <- read.table(genome_file, header = FALSE, sep = "", stringsAsFactors = FALSE)
if (ncol(genome) < 2) stop("genome.bed.tsv 至少需要两列：CHROM 和 LEN")
genome <- genome[, 1:2]
colnames(genome) <- c("CHROM", "LEN")
genome$CHROM <- as.character(genome$CHROM)
genome$LEN <- as.numeric(genome$LEN)

roh <- read.table(roh_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
colnames(roh) <- trimws(colnames(roh))

required_roh <- c("IID", "CHR", "POS1", "POS2", "KB")
missing_roh <- setdiff(required_roh, colnames(roh))
if (length(missing_roh) > 0) {
  stop(paste("ROH文件缺少必要列：", paste(missing_roh, collapse = ", ")))
}

# =========================================================
# 4. 整理 ROH 数据
# =========================================================
cat("[INFO] Processing ROH data...\n")

convert_chr <- function(x) {
  x <- as.character(x)
  out <- x
  is_chr_style <- grepl("^Chr", x)
  out[is_chr_style] <- x[is_chr_style]

  is_num <- suppressWarnings(!is.na(as.numeric(x)))
  out[!is_chr_style & is_num] <- paste0("Chr", sprintf("%02d", as.numeric(x[!is_chr_style & is_num])))
  out
}

roh2 <- roh %>%
  mutate(
    IID = as.character(IID),
    CHR = convert_chr(CHR),
    POS1 = as.numeric(POS1),
    POS2 = as.numeric(POS2),
    KB = as.numeric(KB)
  ) %>%
  mutate(
    ROH_length_bp = POS2 - POS1 + 1,
    ROH_length_Mb = ROH_length_bp / 1e6
  ) %>%
  left_join(pop, by = "IID") %>%
  filter(!is.na(pop))

if (nrow(roh2) == 0) {
  stop("ROH 与 pop.info.tsv 合并后没有数据，请检查 IID 是否一致。")
}

# 过滤掉 >2.5 Mb
roh2 <- roh2 %>% filter(ROH_length_bp < max_roh_bp)

if (nrow(roh2) == 0) {
  stop("过滤 >2.5 Mb ROH 后没有剩余数据。")
}

pop$pop <- factor(pop$pop, levels = c("Domestic", "Wild"))
roh2$pop <- factor(roh2$pop, levels = c("Domestic", "Wild"))

# =========================================================
# 5. 染色体排序与映射
# =========================================================
extract_chr_num <- function(x) {
  x2 <- gsub("^Chr", "", x)
  suppressWarnings(as.numeric(x2))
}

genome <- genome %>%
  mutate(chr_num = extract_chr_num(CHROM)) %>%
  arrange(is.na(chr_num), chr_num, CHROM)

chrom_levels <- as.character(genome$CHROM)

# 给CMplot用数值型染色体编号
genome$CHR_NUM <- seq_len(nrow(genome))

chr_map <- genome %>%
  select(CHROM, CHR_NUM)

# =========================================================
# 6. Panel C：ROH长度类别分布
# =========================================================
cat("[INFO] Preparing panel C...\n")

roh2 <- roh2 %>%
  mutate(
    roh_class = cut(
      ROH_length_bp,
      breaks = c(0, 200e3, 2.5e6),
      labels = labels_bp,
      include.lowest = TRUE,
      right = FALSE
    )
  )

all_ids <- pop %>% distinct(IID, pop)
all_classes <- data.frame(roh_class = factor(labels_bp, levels = labels_bp))

# 每个个体各类别 ROH 数量
count_ind <- roh2 %>%
  group_by(IID, pop, roh_class) %>%
  summarise(roh_count = n(), .groups = "drop")

count_ind_full <- merge(all_ids, all_classes, all = TRUE) %>%
  left_join(count_ind, by = c("IID", "pop", "roh_class")) %>%
  mutate(roh_count = ifelse(is.na(roh_count), 0, roh_count))

count_sum <- count_ind_full %>%
  group_by(pop, roh_class) %>%
  summarise(
    mean_count = mean(roh_count),
    se_count = sd(roh_count) / sqrt(n()),
    .groups = "drop"
  )

# 每个个体各类别 ROH 总长度
length_ind <- roh2 %>%
  group_by(IID, pop, roh_class) %>%
  summarise(total_roh_Mb = sum(ROH_length_Mb), .groups = "drop")

length_ind_full <- merge(all_ids, all_classes, all = TRUE) %>%
  left_join(length_ind, by = c("IID", "pop", "roh_class")) %>%
  mutate(total_roh_Mb = ifelse(is.na(total_roh_Mb), 0, total_roh_Mb))

length_sum <- length_ind_full %>%
  group_by(pop, roh_class) %>%
  summarise(
    mean_length = mean(total_roh_Mb),
    se_length = sd(total_roh_Mb) / sqrt(n()),
    .groups = "drop"
  )

pC_count <- ggplot(count_sum, aes(x = roh_class, y = mean_count, fill = pop)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72, color = "black", linewidth = 0.4) +
  geom_errorbar(
    aes(ymin = mean_count - se_count, ymax = mean_count + se_count),
    position = position_dodge(width = 0.78),
    width = 0.18,
    linewidth = 0.4
  ) +
  scale_fill_manual(values = fill_cols) +
  labs(
    x = "ROH length category",
    y = "Mean number of ROH segments per individual",
    fill = NULL
  ) +
  theme_gray(base_size = 15) +
  theme(
    panel.grid.major = element_line(color = "#B0B0B0", linewidth = 0.65),
    panel.grid.minor = element_line(color = "#C9C9C9", linewidth = 0.4),
    panel.border = element_rect(color = "#4A4A4A", fill = NA, linewidth = 0.9),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(angle = 15, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "top"
  )

pC_length <- ggplot(length_sum, aes(x = roh_class, y = mean_length, fill = pop)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72, color = "black", linewidth = 0.4) +
  geom_errorbar(
    aes(ymin = mean_length - se_length, ymax = mean_length + se_length),
    position = position_dodge(width = 0.78),
    width = 0.18,
    linewidth = 0.4
  ) +
  scale_fill_manual(values = fill_cols) +
  labs(
    x = "ROH length category",
    y = "Mean total ROH length per individual (Mb)",
    fill = NULL
  ) +
  theme_gray(base_size = 15) +
  theme(
    panel.grid.major = element_line(color = "#B0B0B0", linewidth = 0.65),
    panel.grid.minor = element_line(color = "#C9C9C9", linewidth = 0.4),
    panel.border = element_rect(color = "#4A4A4A", fill = NA, linewidth = 0.9),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(angle = 15, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "top"
  )

# =========================================================
# 7. Panel D：FROH 箱线图
# =========================================================
cat("[INFO] Preparing panel D...\n")

froh_df <- roh2 %>%
  group_by(IID, pop) %>%
  summarise(
    total_roh_bp = sum(ROH_length_bp),
    total_roh_Mb = sum(ROH_length_Mb),
    FROH = total_roh_bp / genome_size_bp,
    .groups = "drop"
  )

missing_ids <- setdiff(pop$IID, froh_df$IID)
if (length(missing_ids) > 0) {
  froh_missing <- pop %>%
    filter(IID %in% missing_ids) %>%
    mutate(
      total_roh_bp = 0,
      total_roh_Mb = 0,
      FROH = 0
    )
  froh_df <- bind_rows(froh_df, froh_missing)
}

froh_df$pop <- factor(froh_df$pop, levels = c("Domestic", "Wild"))

wt_froh <- wilcox.test(FROH ~ pop, data = froh_df, exact = FALSE)

if (wt_froh$p.value < 2.2e-16) {
  p_label_froh <- "P < 2.2 × 10^-16"
} else {
  p_label_froh <- paste0("P = ", signif(wt_froh$p.value, 3))
}

ymax_froh <- max(froh_df$FROH, na.rm = TRUE)
ymin_froh <- min(froh_df$FROH, na.rm = TRUE)
yrange_froh <- ymax_froh - ymin_froh
if (yrange_froh <= 0) yrange_froh <- 0.02

y_bracket_froh <- ymax_froh + yrange_froh * 0.08
y_text_froh    <- ymax_froh + yrange_froh * 0.15
y_upper_froh   <- ymax_froh + yrange_froh * 0.24

pD <- ggplot() +
  geom_jitter(
    data = froh_df,
    aes(x = pop, y = FROH, color = pop),
    width = 0.13,
    size = point_size,
    alpha = point_alpha
  ) +
  geom_boxplot(
    data = froh_df,
    aes(x = pop, y = FROH, fill = pop),
    width = 0.38,
    alpha = 0.95,
    outlier.shape = NA,
    color = "black",
    linewidth = 1.0
  ) +
  stat_summary(
    data = froh_df,
    aes(x = pop, y = FROH),
    fun = median,
    geom = "crossbar",
    width = 0.38,
    fatten = 0,
    color = "black",
    linewidth = 0.7
  ) +
  scale_fill_manual(values = fill_cols) +
  scale_color_manual(values = point_cols) +
  labs(
    x = NULL,
    y = expression("Genomic inbreeding coefficient (" * F[ROH] * ")")
  ) +
  annotate("segment", x = 1, xend = 2, y = y_bracket_froh, yend = y_bracket_froh, linewidth = 0.9) +
  annotate("segment", x = 1, xend = 1, y = y_bracket_froh - yrange_froh * 0.015, yend = y_bracket_froh, linewidth = 0.9) +
  annotate("segment", x = 2, xend = 2, y = y_bracket_froh - yrange_froh * 0.015, yend = y_bracket_froh, linewidth = 0.9) +
  annotate("text", x = 1.5, y = y_text_froh, label = p_label_froh, size = 7) +
  coord_cartesian(ylim = c(0, y_upper_froh)) +
  theme_gray(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "#B0B0B0", linewidth = 0.65),
    panel.grid.minor = element_line(color = "#C9C9C9", linewidth = 0.4),
    panel.border = element_rect(color = "#4A4A4A", fill = NA, linewidth = 0.9),
    axis.title.y = element_text(size = 18, face = "bold", color = "black"),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  )

# =========================================================
# 8. 加速计算 ROH frequency
# =========================================================
cat("[INFO] Calculating ROH frequency...\n")

calc_roh_freq_fast <- function(roh_sub, sample_ids, genome_df, window_size, step_size) {
  sample_ids <- unique(sample_ids)
  n_ind <- length(sample_ids)

  out_list <- vector("list", nrow(genome_df))

  for (k in seq_len(nrow(genome_df))) {
    chr <- as.character(genome_df$CHROM[k])
    chr_len <- as.numeric(genome_df$LEN[k])

    starts <- seq(0, chr_len - 1, by = step_size)
    ends <- starts + window_size
    ends[ends > chr_len] <- chr_len
    n_win <- length(starts)

    coverage_count <- integer(n_win)

    roh_chr <- roh_sub[roh_sub$CHR == chr, , drop = FALSE]

    if (nrow(roh_chr) > 0) {
      split_by_id <- split(roh_chr, roh_chr$IID)

      for (sid in names(split_by_id)) {
        dfi <- split_by_id[[sid]]
        diff_vec <- integer(n_win + 1)

        for (i in seq_len(nrow(dfi))) {
          s <- dfi$POS1[i]
          e <- dfi$POS2[i]

          j_start <- ceiling((s - window_size) / step_size)
          j_end   <- floor(e / step_size)

          if (is.na(j_start) || is.na(j_end)) next

          j_start <- max(0, j_start)
          j_end   <- min(n_win - 1, j_end)

          if (j_start <= j_end) {
            diff_vec[j_start + 1] <- diff_vec[j_start + 1] + 1
            diff_vec[j_end + 2]   <- diff_vec[j_end + 2] - 1
          }
        }

        cov_id <- cumsum(diff_vec)[1:n_win]
        coverage_count <- coverage_count + as.integer(cov_id > 0)
      }
    }

    out_list[[k]] <- data.frame(
      CHROM = chr,
      start = starts,
      end = ends,
      n_overlap = coverage_count,
      roh_freq = coverage_count / n_ind,
      stringsAsFactors = FALSE
    )
  }

  bind_rows(out_list)
}

dom_ids <- pop %>% filter(pop == "Domestic") %>% pull(IID)
wild_ids <- pop %>% filter(pop == "Wild") %>% pull(IID)

roh_dom <- roh2 %>% filter(IID %in% dom_ids)
roh_wild <- roh2 %>% filter(IID %in% wild_ids)

freq_dom <- calc_roh_freq_fast(roh_dom, dom_ids, genome, window_size, step_size) %>%
  mutate(pop = "Domestic")

freq_wild <- calc_roh_freq_fast(roh_wild, wild_ids, genome, window_size, step_size) %>%
  mutate(pop = "Wild")

freq_by_pop <- bind_rows(freq_dom, freq_wild)

# 映射CMplot所需数值染色体
freq_by_pop <- freq_by_pop %>%
  left_join(chr_map, by = "CHROM") %>%
  mutate(
    BIN_POS = floor((start + end) / 2)
  )

freq_wide <- freq_by_pop %>%
  select(CHROM, CHR_NUM, start, end, BIN_POS, pop, roh_freq) %>%
  pivot_wider(names_from = pop, values_from = roh_freq)

freq_wide <- freq_wide %>%
  mutate(
    Delta = Domestic - Wild,
    Domestic_biased = ifelse(Delta > 0, Delta, 0)
  )

# =========================================================
# 9. 输出CMplot输入文件
# =========================================================
cat("[INFO] Writing CMplot input tables...\n")

cm_dom <- freq_wide %>%
  mutate(Marker = paste0("win_", CHROM, "_", start, "_", end)) %>%
  select(Marker, CHR_NUM, BIN_POS, Domestic)

cm_delta <- freq_wide %>%
  mutate(Marker = paste0("win_", CHROM, "_", start, "_", end)) %>%
  select(Marker, CHR_NUM, BIN_POS, Domestic_biased)

write.table(
  cm_dom,
  file = paste0(prefix, "_CMplot_domestic_frequency.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  cm_delta,
  file = paste0(prefix, "_CMplot_domestic_biased_delta.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# =========================================================
# 10. 用CMplot作图（兼容旧版 + 彩色）
# =========================================================
cat("[INFO] Plotting Manhattan with CMplot...\n")

# 10.1 Domestic frequency Manhattan
pdf(file = paste0(prefix, "_CMplot_domestic_frequency.pdf"), width = 13, height = 4.5)
CMplot(
  cm_dom,
  plot.type = "m",
  LOG10 = FALSE,
  col = manhattan_cols,
  cex = cmplot_cex,
  ylab = "Domestic ROH frequency",
  amplify = FALSE,
  threshold = NULL,
  file.output = FALSE,
  verbose = FALSE
)
dev.off()

png(filename = paste0(prefix, "_CMplot_domestic_frequency.png"),
    width = 13, height = 4.5, units = "in", res = 600)
CMplot(
  cm_dom,
  plot.type = "m",
  LOG10 = FALSE,
  col = manhattan_cols,
  cex = cmplot_cex,
  ylab = "Domestic ROH frequency",
  amplify = FALSE,
  threshold = NULL,
  file.output = FALSE,
  verbose = FALSE
)
dev.off()

# 10.2 Domestic-biased delta Manhattan
pdf(file = paste0(prefix, "_CMplot_domestic_biased_delta.pdf"), width = 13, height = 4.5)
CMplot(
  cm_delta,
  plot.type = "m",
  LOG10 = FALSE,
  col = manhattan_cols,
  cex = cmplot_cex,
  ylab = "Domestic-biased ΔROH frequency",
  amplify = FALSE,
  threshold = NULL,
  file.output = FALSE,
  verbose = FALSE
)
dev.off()

png(filename = paste0(prefix, "_CMplot_domestic_biased_delta.png"),
    width = 13, height = 4.5, units = "in", res = 600)
CMplot(
  cm_delta,
  plot.type = "m",
  LOG10 = FALSE,
  col = manhattan_cols,
  cex = cmplot_cex,
  ylab = "Domestic-biased ΔROH frequency",
  amplify = FALSE,
  threshold = NULL,
  file.output = FALSE,
  verbose = FALSE
)
dev.off()

# =========================================================
# 11. 输出其它结果
# =========================================================
write.table(
  freq_by_pop,
  file = paste0(prefix, "_recalculated_roh_frequency_by_pop.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  freq_wide,
  file = paste0(prefix, "_delta_roh_frequency.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  froh_df,
  file = paste0(prefix, "_individual_FROH.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

sink(paste0(prefix, "_statistics.txt"))
cat("ROH analyses summary\n")
cat("====================\n\n")
cat("Genome size (bp): ", genome_size_bp, "\n", sep = "")
cat("Window size (bp): ", window_size, "\n", sep = "")
cat("Step size (bp): ", step_size, "\n", sep = "")
cat("ROH used: < 2.5 Mb only\n\n")

cat("ROH length classes:\n")
cat("  - <200 kb\n")
cat("  - 200 kb–2.5 Mb\n\n")

cat("FROH Wilcoxon rank-sum test\n")
cat("---------------------------\n")
print(wt_froh)
cat("\n")

cat("FROH summary by population\n")
cat("--------------------------\n")
print(
  froh_df %>%
    group_by(pop) %>%
    summarise(
      n = n(),
      mean_FROH = mean(FROH, na.rm = TRUE),
      median_FROH = median(FROH, na.rm = TRUE),
      sd_FROH = sd(FROH, na.rm = TRUE),
      min_FROH = min(FROH, na.rm = TRUE),
      max_FROH = max(FROH, na.rm = TRUE),
      .groups = "drop"
    )
)
cat("\n")
sink()

# =========================================================
# 12. 保存 Panel C / D
# =========================================================
ggsave(paste0(prefix, "_panelC_count.pdf"),  pC_count,  width = 6.2, height = 5.5)
ggsave(paste0(prefix, "_panelC_count.png"),  pC_count,  width = 6.2, height = 5.5, dpi = 600)

ggsave(paste0(prefix, "_panelC_length.pdf"), pC_length, width = 6.2, height = 5.5)
ggsave(paste0(prefix, "_panelC_length.png"), pC_length, width = 6.2, height = 5.5, dpi = 600)

ggsave(paste0(prefix, "_panelD_FROH_boxplot.pdf"), pD, width = 5.2, height = 6.0)
ggsave(paste0(prefix, "_panelD_FROH_boxplot.png"), pD, width = 5.2, height = 6.0, dpi = 600)

combined_CD <- gridExtra::arrangeGrob(pC_count, pC_length, pD, ncol = 3)
ggsave(paste0(prefix, "_panelC_D_combined.pdf"), combined_CD, width = 16.5, height = 5.8)
ggsave(paste0(prefix, "_panelC_D_combined.png"), combined_CD, width = 16.5, height = 5.8, dpi = 600)

cat("[INFO] Done!\n")
cat("[INFO] Output files generated.\n")
