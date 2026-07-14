library(tidyverse)
library(rstatix)
library(ggpubr)

#========================
# 1. 读取数据
#========================
dat <- read.delim("genetic_load_inside_outside_ROH.tsv",
                  sep = "\t", header = TRUE, stringsAsFactors = FALSE)

meta <- read.delim("wgs_samples.tsv",
                   sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
  rename(sample = sample_id)

dat <- left_join(dat, meta, by = "sample")

# 强制统一 group 命名和顺序
dat <- dat %>%
  filter(class %in% c("LoF", "missense_severe", "missense_all")) %>%
  mutate(
    class = factor(class, levels = c("LoF", "missense_severe", "missense_all")),
    region = factor(region, levels = c("outside_ROH", "inside_ROH")),
    group = case_when(
      group %in% c("wild", "Wild") ~ "Wild",
      group %in% c("domestic", "Domestic") ~ "Domestic",
      TRUE ~ group
    ),
    group = factor(group, levels = c("Domestic", "Wild"))
  )

#========================
# 2. 统计检验
#========================
paired_test <- dat %>%
  group_by(class) %>%
  wilcox_test(genetic_load ~ region, paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

write.table(paired_test,
            "genetic_load_inside_vs_outside_ROH_paired.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

paired_group_test <- dat %>%
  group_by(group, class) %>%
  wilcox_test(genetic_load ~ region, paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

write.table(paired_group_test,
            "genetic_load_inside_vs_outside_ROH_by_group.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

#========================
# 3. 配色
#========================
# 统一成参考图配色
col_group_fill <- c(
  "Domestic" = "#F07F72",
  "Wild"     = "#18B7C9"
)

col_group_point <- c(
  "Domestic" = "#F07F72",
  "Wild"     = "#18B7C9"
)

col_region_fill <- c(
  "outside_ROH" = "#D9D9D9",
  "inside_ROH"  = "#E88A80"
)

#========================
# 4. 统一主题
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

      plot.margin = margin(10, 14, 10, 10),
      panel.spacing = unit(1.0, "lines")
    )
}

#========================
# 5. P值标签
#========================
format_p <- function(p) {
  ifelse(p < 2.2e-16,
         "P < 2.2e-16",
         paste0("P = ", format(p, scientific = TRUE, digits = 3)))
}

p1_label <- dat %>%
  group_by(class) %>%
  summarise(y.position = max(genetic_load, na.rm = TRUE) * 1.14, .groups = "drop") %>%
  left_join(
    paired_test %>%
      mutate(
        xmin = "outside_ROH",
        xmax = "inside_ROH",
        label = format_p(p.adj)
      ),
    by = "class"
  )

p2_label <- dat %>%
  group_by(group, class) %>%
  summarise(y.position = max(genetic_load, na.rm = TRUE) * 1.16, .groups = "drop") %>%
  left_join(
    paired_group_test %>%
      mutate(
        xmin = "outside_ROH",
        xmax = "inside_ROH",
        label = format_p(p.adj)
      ),
    by = c("group", "class")
  )

#========================
# 6. 图1：总体 inside vs outside
#========================
p1 <- ggplot(dat, aes(x = region, y = genetic_load)) +
  geom_line(aes(group = sample),
            color = "grey65", alpha = 0.35, linewidth = 0.8) +
  geom_boxplot(aes(fill = region),
               width = 0.50,
               outlier.shape = NA,
               alpha = 0.92,
               color = "black",
               linewidth = 0.8,
               median.linewidth = 0.9) +
  geom_point(aes(color = group),
             position = position_jitter(width = 0.05, height = 0),
             size = 2.6,
             alpha = 0.60) +
  stat_pvalue_manual(
    p1_label,
    label = "label",
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    tip.length = 0.012,
    bracket.size = 0.7,
    size = 4.8
  ) +
  facet_wrap(~class, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = col_region_fill) +
  scale_color_manual(values = col_group_point) +
  labs(x = "", y = "Genetic load", color = "group", fill = "region") +
  theme_fm(base_size = 15) +
  theme(
    legend.position = "right"
  )

ggsave("Fig_genetic_load_inside_vs_outside_ROH.pdf", p1, width = 11, height = 4.8)
ggsave("Fig_genetic_load_inside_vs_outside_ROH.jpg", p1, width = 11, height = 4.8, dpi = 600)

#========================
# 7. 图2：按 group 分层 inside vs outside
#========================
p2 <- ggplot(dat, aes(x = region, y = genetic_load)) +
  geom_boxplot(aes(fill = group),
               width = 0.45,
               outlier.shape = NA,
               alpha = 0.78,
               color = "black",
               linewidth = 0.75,
               median.linewidth = 0.85) +
  geom_point(aes(color = group),
             position = position_jitter(width = 0.045, height = 0),
             size = 2.5,
             alpha = 0.60) +
  stat_pvalue_manual(
    p2_label,
    label = "label",
    xmin = "xmin",
    xmax = "xmax",
    y.position = "y.position",
    tip.length = 0.012,
    bracket.size = 0.65,
    size = 3.7
  ) +
  facet_grid(group ~ class, scales = "free_y") +
  scale_fill_manual(values = col_group_fill) +
  scale_color_manual(values = col_group_point) +
  labs(x = "", y = "Genetic load", fill = "group", color = "group") +
  theme_fm(base_size = 13.5) +
  theme(
    legend.position = "right",
    panel.spacing.x = unit(1.1, "lines"),
    panel.spacing.y = unit(1.2, "lines"),
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 11),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 12)
  )

ggsave("Fig_genetic_load_inside_vs_outside_ROH_by_group.pdf", p2, width = 10.5, height = 8.6)
ggsave("Fig_genetic_load_inside_vs_outside_ROH_by_group.jpg", p2, width = 10.5, height = 8.6, dpi = 600)

#========================
# 8. 输出合并数据
#========================
write.table(dat,
            "genetic_load_inside_outside_ROH_merged.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
