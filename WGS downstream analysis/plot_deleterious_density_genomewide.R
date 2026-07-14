library(tidyverse)

if (!requireNamespace("ggrastr", quietly = TRUE)) {
  stop("Package 'ggrastr' is required to rasterize point layers while keeping other plot elements vector. Install it with install.packages('ggrastr').")
}
if (!requireNamespace("ragg", quietly = TRUE)) {
  stop("Package 'ragg' is required for transparent point-layer rasterization. Install it with install.packages('ragg').")
}

point_raster_dpi <- 600
point_raster_dev <- "ragg_png"

# =========================
# 1. Read data
# =========================
dat <- read.delim("deleterious_window_density.tsv",
                  sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# 清理染色体名称
dat$CHROM <- trimws(as.character(dat$CHROM))

# 固定染色体顺序
chrom_levels <- c(sprintf("Chr%02d", 1:28), "ChrX", "ChrY")
dat$CHROM <- factor(dat$CHROM, levels = chrom_levels, ordered = TRUE)

# =========================
# 2. Build cumulative genome coordinates
# =========================
chr_sizes <- dat %>%
  filter(!is.na(CHROM)) %>%
  group_by(CHROM) %>%
  summarise(chr_len = as.numeric(max(end, na.rm = TRUE)), .groups = "drop") %>%
  arrange(CHROM) %>%
  mutate(
    chr_len = as.numeric(chr_len),
    tot = cumsum(chr_len) - chr_len
  )

dat <- dat %>%
  left_join(chr_sizes, by = "CHROM") %>%
  mutate(
    start = as.numeric(start),
    end   = as.numeric(end),
    BPcum = start + tot
  )

axisdf <- chr_sizes %>%
  mutate(center = tot + chr_len / 2)

bgdf <- chr_sizes %>%
  mutate(
    xmin = tot,
    xmax = tot + chr_len,
    ymin = -Inf,
    ymax = Inf,
    bg_group = rep(c("odd", "even"), length.out = n())
  )

# =========================
# 3. Deterministic grouped downsampling
# =========================
step_single <- 3
step_combined <- 4

dat_plot1 <- dat %>%
  arrange(CHROM, start) %>%
  group_by(CHROM) %>%
  slice(seq(1, n(), by = step_single)) %>%
  ungroup()

dat_plot2 <- dat %>%
  arrange(CHROM, start) %>%
  group_by(CHROM) %>%
  slice(seq(1, n(), by = step_single)) %>%
  ungroup()

dat_long <- dat %>%
  select(CHROM, start, end, BPcum, LoF_density, severe_density) %>%
  pivot_longer(
    cols = c(LoF_density, severe_density),
    names_to = "class",
    values_to = "density"
  )

# panel 标题改成正式名称
dat_long$class <- factor(dat_long$class,
                         levels = c("LoF_density", "severe_density"),
                         labels = c("Loss-of-function", "Severe missense"))

dat_long_plot <- dat_long %>%
  arrange(class, CHROM, start) %>%
  group_by(class, CHROM) %>%
  slice(seq(1, n(), by = step_combined)) %>%
  ungroup()

# 为 combined 图创建联合颜色变量
dat_long_plot <- dat_long_plot %>%
  mutate(class_chr = paste(class, CHROM, sep = "__"))

# =========================
# 4. Colors
# =========================
chr_colors_blue <- rep(c("#4E79A7", "#A0CBE8"), length.out = length(chrom_levels))
names(chr_colors_blue) <- chrom_levels
chr_colors_blue["ChrX"] <- "#7B6FD0"
chr_colors_blue["ChrY"] <- "#B07AA1"

chr_colors_red <- rep(c("#E15759", "#FF9D9A"), length.out = length(chrom_levels))
names(chr_colors_red) <- chrom_levels
chr_colors_red["ChrX"] <- "#C44E52"
chr_colors_red["ChrY"] <- "#8172B2"

combined_colors <- c(
  setNames(chr_colors_blue, paste("Loss-of-function", names(chr_colors_blue), sep = "__")),
  setNames(chr_colors_red,  paste("Severe missense", names(chr_colors_red), sep = "__"))
)

# =========================
# 5. Themes
# =========================
theme_genome <- theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# 组合图专用主题：去掉黑框，只保留坐标轴
theme_combined <- theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.border = element_blank(),              # 去掉 panel 黑色边框
    strip.background = element_blank(),          # 去掉 strip 背景框
    strip.text = element_text(face = "bold"),    # 保留 panel 标题文字
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(colour = "black"),  # 保留 X/Y 坐标轴线
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# =========================
# 6. LoF density plot
# =========================
p1 <- ggplot() +
  geom_rect(data = bgdf,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = bg_group),
            inherit.aes = FALSE, alpha = 0.18, colour = NA) +
  ggrastr::geom_point_rast(data = dat_plot1,
                           aes(x = BPcum, y = LoF_density, color = CHROM),
                           size = 0.5, alpha = 0.65,
                           raster.dpi = point_raster_dpi,
                           dev = point_raster_dev) +
  scale_fill_manual(values = c("odd" = "grey92", "even" = "white"), guide = "none") +
  scale_color_manual(values = chr_colors_blue, drop = FALSE) +
  scale_x_continuous(
    breaks = axisdf$center,
    labels = as.character(axisdf$CHROM),
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  labs(x = "Chromosome", y = "LoF density", title = "Genome-wide LoF density") +
  theme_genome

# =========================
# 7. Severe missense density plot
# =========================
p2 <- ggplot() +
  geom_rect(data = bgdf,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = bg_group),
            inherit.aes = FALSE, alpha = 0.18, colour = NA) +
  ggrastr::geom_point_rast(data = dat_plot2,
                           aes(x = BPcum, y = severe_density, color = CHROM),
                           size = 0.5, alpha = 0.65,
                           raster.dpi = point_raster_dpi,
                           dev = point_raster_dev) +
  scale_fill_manual(values = c("odd" = "grey92", "even" = "white"), guide = "none") +
  scale_color_manual(values = chr_colors_red, drop = FALSE) +
  scale_x_continuous(
    breaks = axisdf$center,
    labels = as.character(axisdf$CHROM),
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  labs(x = "Chromosome", y = "Severe missense density", title = "Genome-wide severe missense density") +
  theme_genome

# =========================
# 8. Combined plot
# =========================
p3 <- ggplot() +
  geom_rect(data = bgdf,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = bg_group),
            inherit.aes = FALSE, alpha = 0.18, colour = NA) +
  ggrastr::geom_point_rast(data = dat_long_plot,
                           aes(x = BPcum, y = density, color = class_chr),
                           size = 0.5, alpha = 0.65,
                           raster.dpi = point_raster_dpi,
                           dev = point_raster_dev) +
  facet_wrap(~class, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("odd" = "grey92", "even" = "white"), guide = "none") +
  scale_color_manual(values = combined_colors, drop = FALSE) +
  scale_x_continuous(
    breaks = axisdf$center,
    labels = as.character(axisdf$CHROM),
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  labs(x = "Chromosome", y = "Density", title = "Genome-wide deleterious variant density") +
  theme_combined

# =========================
# 9. Save plots
# =========================
ggsave("Fig3D_LoF_density_genomewide.pdf", p1, width = 18, height = 4.5)
ggsave("Fig3D_severe_density_genomewide.pdf", p2, width = 18, height = 4.5)
ggsave("Fig3D_deleterious_density_combined.pdf", p3, width = 18, height = 6.5)

ggsave("Fig3D_LoF_density_genomewide.png", p1, width = 18, height = 4.5, dpi = 300)
ggsave("Fig3D_severe_density_genomewide.png", p2, width = 18, height = 4.5, dpi = 300)
ggsave("Fig3D_deleterious_density_combined.png", p3, width = 18, height = 6.5, dpi = 300)

# =========================
# 10. Export plotting table
# =========================
write.table(dat, "deleterious_density_for_plot.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
