#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# ==============================
# 1. 参数设置
# ==============================
input_dir <- "XPCLR"
output_prefix <- file.path(input_dir, "all_chr")

# ==============================
# 2. 查找结果文件
# ==============================
files <- list.files(
  path = input_dir,
  pattern = "\\.xpclr\\.tsv\\..*\\.xpclr$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No xpclrs result files found in: ", input_dir)
}

cat("Found", length(files), "xpclrs result files\n")

extract_chr_from_filename <- function(f) {
  x <- basename(f)
  x <- sub("\\.xpclr\\.tsv\\..*$", "", x)
  return(x)
}

# ==============================
# 3. 读取并合并文件
# ==============================
read_one_file <- function(f) {
  message("Reading: ", f)
  df <- readr::read_tsv(f, show_col_types = FALSE)
  df$SOURCE_FILE <- basename(f)
  df$FILE_CHROM <- extract_chr_from_filename(f)
  return(df)
}

lst <- lapply(files, read_one_file)
xpclr_raw <- bind_rows(lst)

write_tsv(xpclr_raw, paste0(output_prefix, ".xpclrs.raw.tsv"))

cat("Merged rows:", nrow(xpclr_raw), "\n")
cat("Merged columns:\n")
print(colnames(xpclr_raw))

# ==============================
# 4. 检查必要列
# ==============================
required_cols <- c("chrom", "start", "stop", "xpclr")
missing_cols <- setdiff(required_cols, colnames(xpclr_raw))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ==============================
# 5. 标准化结果表
# ==============================
xpclr_std <- xpclr_raw %>%
  mutate(
    CHROM = ifelse(is.na(chrom) | chrom == "", FILE_CHROM, chrom),
    start = as.numeric(start),
    end   = as.numeric(stop),
    POS   = floor((start + end) / 2),
    XPCLR = as.numeric(xpclr),
    XPCLR_NORM = if ("xpclr_norm" %in% colnames(.)) as.numeric(xpclr_norm) else NA_real_
  ) %>%
  select(
    CHROM, POS, start, end, XPCLR, XPCLR_NORM,
    pos_start, pos_stop, modelL, nullL, sel_coef, nSNPs, nSNPs_avail,
    SOURCE_FILE, FILE_CHROM, everything()
  ) %>%
  arrange(CHROM, start, end)

write_tsv(xpclr_std, paste0(output_prefix, ".xpclrs.standardized.tsv"))

# ==============================
# 6. 基本统计
# ==============================
stats_df <- xpclr_std %>%
  summarise(
    n_windows = n(),
    min_xpclr = min(XPCLR, na.rm = TRUE),
    max_xpclr = max(XPCLR, na.rm = TRUE),
    mean_xpclr = mean(XPCLR, na.rm = TRUE),
    median_xpclr = median(XPCLR, na.rm = TRUE),
    sd_xpclr = sd(XPCLR, na.rm = TRUE),
    min_xpclr_norm = min(XPCLR_NORM, na.rm = TRUE),
    max_xpclr_norm = max(XPCLR_NORM, na.rm = TRUE),
    mean_xpclr_norm = mean(XPCLR_NORM, na.rm = TRUE),
    median_xpclr_norm = median(XPCLR_NORM, na.rm = TRUE),
    sd_xpclr_norm = sd(XPCLR_NORM, na.rm = TRUE)
  )

write_tsv(stats_df, paste0(output_prefix, ".xpclrs.summary.tsv"))
print(stats_df)

# ==============================
# 7. 提取候选窗口（基于 XPCLR）
# ==============================
cutoff_1pct <- quantile(xpclr_std$XPCLR, 0.99, na.rm = TRUE)
cutoff_0.5pct <- quantile(xpclr_std$XPCLR, 0.995, na.rm = TRUE)

candidate_1pct <- xpclr_std %>%
  filter(!is.na(XPCLR), XPCLR >= cutoff_1pct) %>%
  arrange(desc(XPCLR))

candidate_0.5pct <- xpclr_std %>%
  filter(!is.na(XPCLR), XPCLR >= cutoff_0.5pct) %>%
  arrange(desc(XPCLR))

write_tsv(candidate_1pct, paste0(output_prefix, ".top1pct.tsv"))
write_tsv(candidate_0.5pct, paste0(output_prefix, ".top0.5pct.tsv"))

cat("Top 1% cutoff   =", cutoff_1pct, "\n")
cat("Top 0.5% cutoff =", cutoff_0.5pct, "\n")
cat("Top 1% windows   =", nrow(candidate_1pct), "\n")
cat("Top 0.5% windows =", nrow(candidate_0.5pct), "\n")

# ==============================
# 8. 提取候选窗口（基于 XPCLR_NORM）
# ==============================
if (all(is.na(xpclr_std$XPCLR_NORM))) {
  cat("XPCLR_NORM column not available or all NA, skip normalized-score candidate extraction.\n")
  candidate_norm_1pct <- tibble()
  candidate_norm_0.5pct <- tibble()
} else {
  cutoff_norm_1pct <- quantile(xpclr_std$XPCLR_NORM, 0.99, na.rm = TRUE)
  cutoff_norm_0.5pct <- quantile(xpclr_std$XPCLR_NORM, 0.995, na.rm = TRUE)

  candidate_norm_1pct <- xpclr_std %>%
    filter(!is.na(XPCLR_NORM), XPCLR_NORM >= cutoff_norm_1pct) %>%
    arrange(desc(XPCLR_NORM))

  candidate_norm_0.5pct <- xpclr_std %>%
    filter(!is.na(XPCLR_NORM), XPCLR_NORM >= cutoff_norm_0.5pct) %>%
    arrange(desc(XPCLR_NORM))

  write_tsv(candidate_norm_1pct, paste0(output_prefix, ".norm.top1pct.tsv"))
  write_tsv(candidate_norm_0.5pct, paste0(output_prefix, ".norm.top0.5pct.tsv"))

  cat("Top 1% cutoff (XPCLR_NORM)   =", cutoff_norm_1pct, "\n")
  cat("Top 0.5% cutoff (XPCLR_NORM) =", cutoff_norm_0.5pct, "\n")
  cat("Top 1% windows (XPCLR_NORM)   =", nrow(candidate_norm_1pct), "\n")
  cat("Top 0.5% windows (XPCLR_NORM) =", nrow(candidate_norm_0.5pct), "\n")
}

# ==============================
# 9. Z-score 标准化（基于 XPCLR）
# ==============================
xpclr_z <- xpclr_std %>%
  mutate(
    XPCLR_Z = (XPCLR - mean(XPCLR, na.rm = TRUE)) / sd(XPCLR, na.rm = TRUE)
  )

write_tsv(xpclr_z, paste0(output_prefix, ".zscore.tsv"))

candidate_z3 <- xpclr_z %>%
  filter(!is.na(XPCLR_Z), XPCLR_Z >= 3) %>%
  arrange(desc(XPCLR_Z))

write_tsv(candidate_z3, paste0(output_prefix, ".zscore_ge3.tsv"))

cat("Z >= 3 windows =", nrow(candidate_z3), "\n")

# ==============================
# 10. 合并相邻候选窗口
# ==============================
merge_adjacent <- function(df, gap = 1, score_col = "XPCLR") {
  if (nrow(df) == 0) return(df)

  df <- df %>%
    arrange(CHROM, start, end)

  merged <- list()
  current <- df[1, , drop = FALSE]

  if (nrow(df) == 1) return(df)

  for (i in 2:nrow(df)) {
    x <- df[i, , drop = FALSE]

    same_chr <- current$CHROM == x$CHROM
    overlap_or_close <- x$start <= (current$end + gap)

    if (same_chr && overlap_or_close) {
      current$end <- max(current$end, x$end)
      current[[score_col]] <- max(current[[score_col]], x[[score_col]], na.rm = TRUE)

      if ("XPCLR_NORM" %in% colnames(df)) {
        current$XPCLR_NORM <- max(current$XPCLR_NORM, x$XPCLR_NORM, na.rm = TRUE)
      }
      if ("XPCLR_Z" %in% colnames(df)) {
        current$XPCLR_Z <- max(current$XPCLR_Z, x$XPCLR_Z, na.rm = TRUE)
      }
    } else {
      merged[[length(merged) + 1]] <- current
      current <- x
    }
  }

  merged[[length(merged) + 1]] <- current
  bind_rows(merged)
}

candidate_1pct_merged <- merge_adjacent(candidate_1pct, gap = 1, score_col = "XPCLR")
candidate_0.5pct_merged <- merge_adjacent(candidate_0.5pct, gap = 1, score_col = "XPCLR")
candidate_z3_merged <- merge_adjacent(candidate_z3, gap = 1, score_col = "XPCLR_Z")

write_tsv(candidate_1pct_merged, paste0(output_prefix, ".top1pct.merged.tsv"))
write_tsv(candidate_0.5pct_merged, paste0(output_prefix, ".top0.5pct.merged.tsv"))
write_tsv(candidate_z3_merged, paste0(output_prefix, ".zscore_ge3.merged.tsv"))

cat("Merged top1% regions   =", nrow(candidate_1pct_merged), "\n")
cat("Merged top0.5% regions =", nrow(candidate_0.5pct_merged), "\n")
cat("Merged Z>=3 regions    =", nrow(candidate_z3_merged), "\n")

# 如果有 XPCLR_NORM 候选，也做合并
if (exists("candidate_norm_1pct") && nrow(candidate_norm_1pct) > 0) {
  candidate_norm_1pct_merged <- merge_adjacent(candidate_norm_1pct, gap = 1, score_col = "XPCLR_NORM")
  write_tsv(candidate_norm_1pct_merged, paste0(output_prefix, ".norm.top1pct.merged.tsv"))
}

if (exists("candidate_norm_0.5pct") && nrow(candidate_norm_0.5pct) > 0) {
  candidate_norm_0.5pct_merged <- merge_adjacent(candidate_norm_0.5pct, gap = 1, score_col = "XPCLR_NORM")
  write_tsv(candidate_norm_0.5pct_merged, paste0(output_prefix, ".norm.top0.5pct.merged.tsv"))
}

# ==============================
# 11. 导出 BED
# ==============================
candidate_1pct_merged %>%
  transmute(CHROM, start = pmax(0, start - 1), end) %>%
  write_tsv(paste0(output_prefix, ".top1pct.merged.bed"), col_names = FALSE)

candidate_0.5pct_merged %>%
  transmute(CHROM, start = pmax(0, start - 1), end) %>%
  write_tsv(paste0(output_prefix, ".top0.5pct.merged.bed"), col_names = FALSE)

candidate_z3_merged %>%
  transmute(CHROM, start = pmax(0, start - 1), end) %>%
  write_tsv(paste0(output_prefix, ".zscore_ge3.merged.bed"), col_names = FALSE)

cat("BED files exported.\n")
cat("Done.\n")
