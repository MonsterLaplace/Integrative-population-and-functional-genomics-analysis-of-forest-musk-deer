library(tidyverse)
library(rstatix)
library(ggpubr)

#========================
# 1. 读取数据
#========================
load_df <- read.delim("genetic_load_per_sample.tsv",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)

meta <- read.delim("wgs_samples.tsv",
                   sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
  rename(sample = sample_id)

dat <- left_join(load_df, meta, by = "sample")

dat <- dat %>%
  filter(class %in% c("LoF", "missense_severe", "missense_all")) %>%
  mutate(
    class = factor(class, levels = c("LoF", "missense_severe", "missense_all")),
    group = factor(group, levels = c("wild", "domestic"))
  )

#========================
# 2. 汇总统计
#========================
summary_stats <- dat %>%
  group_by(group, class) %>%
  summarise(
    n = n(),
    mean_load = mean(genetic_load, na.rm = TRUE),
    median_load = median(genetic_load, na.rm = TRUE),
    sd_load = sd(genetic_load, na.rm = TRUE),
    min_load = min(genetic_load, na.rm = TRUE),
    max_load = max(genetic_load, na.rm = TRUE),
    .groups = "drop"
  )

write.table(summary_stats,
            "genetic_load_group_summary.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

#========================
# 3. domestic vs wild 统计检验
#========================
stat_test <- dat %>%
  group_by(class) %>%
  wilcox_test(genetic_load ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

write.table(stat_test,
            "genetic_load_wilcox_domestic_vs_wild.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

#========================
# 4. 配色
#========================
col_group_fill <- c(
  "wild"     = "#18B7C9",
  "domestic" = "#F07F72"
)

col_group_point <- c(
  "wild"     = "#18B7C9",
  "domestic" = "#F07F72"
)


#========================
# 5. 统一主题
#========================
theme_fm <- function(base_size = 15) {
  theme_gray(base_size = base_size) +
    theme(
      panel.background = element_rect(fill = "#EBEBEB", color = NA),
      panel.grid.major = element_line(color = "#A8A8A8", linewidth = 0.7),
      panel.grid.minor = element_line(color = "#D0D0D0", linewidth = 0.45),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),

      axis.title = element_text(color = "black", size = base_size + 1),
      axis.text = element_text(color = "black", size = base_size),
      axis.text.x = element_text(size = base_size + 1),
      axis.ticks = element_line(color = "black", linewidth = 0.7),

      strip.background = element_rect(fill = "#D9D9D9", color = "black", linewidth = 0.8),
      strip.text = element_text(size = base_size, color = "black"),

      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),

      plot.margin = margin(10, 12, 10, 10),
      panel.spacing = unit(1.0, "lines")
    )
}

#========================
# 6. P值格式和标注位置
#========================
format_p <- function(p) {
  ifelse(p < 2.2e-16,
         "P < 2.2e-16",
         paste0("P = ", format(p, scientific = TRUE, digits = 3)))
}

stat_label <- dat %>%
  group_by(class) %>%
  summarise(y.position = max(genetic_load, na.rm = TRUE) * 1.12, .groups = "drop") %>%
  left_join(
    stat_test %>%
      mutate(
        xmin = "wild",
        xmax = "domestic",
        label = format_p(p.adj)
      ),
    by = "class"
  )

#========================
# 7. 作图
#========================
p <- ggplot(dat, aes(x = group, y = genetic_load)) +
  geom_boxplot(aes(fill = group),
               outlier.shape = NA,
               alpha = 0.88,
               width = 0.55,
               color = "black",
               linewidth = 0.8,
               median.linewidth = 0.9) +
  geom_point(aes(color = group),
             position = position_jitter(width = 0.06, height = 0),
             size = 2.6,
             alpha = 0.8) +
  stat_pvalue_manual(
    stat_label,
    label = "label",
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    tip.length = 0.012,
    bracket.size = 0.7,
    size = 4.8
  ) +
  facet_wrap(~class, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = col_group_fill) +
  scale_color_manual(values = col_group_point) +
  labs(x = "", y = "Genetic load") +
  theme_fm(base_size = 15) +
  theme(
    legend.position = "none"
  )

ggsave("Fig_genetic_load_domestic_vs_wild.pdf", p, width = 10.5, height = 4.8)
ggsave("Fig_genetic_load_domestic_vs_wild.jpg", p, width = 10.5, height = 4.8, dpi = 600)

#========================
# 8. 输出整合表
#========================
write.table(dat,
            "genetic_load_merged_for_plot.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
