library(DESeq2)
library(tidyverse)
library(pheatmap)

if (!requireNamespace("ggrastr", quietly = TRUE)) {
  stop(
    "Package 'ggrastr' is required to rasterize volcano plot points. ",
    "Install it first with: install.packages('ggrastr')"
  )
}

dir.create("05.bulkRNA/05.deseq2", showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1. Read count matrix
# -----------------------------
counts <- read.delim(
  "05.bulkRNA/04.counts/gene_counts.txt",
  comment.char = "#",
  check.names = FALSE
)

count_mat <- as.matrix(counts[, 7:ncol(counts), drop = FALSE])
rownames(count_mat) <- counts$Geneid
colnames(count_mat) <- gsub("05.bulkRNA/03.bam/|\\.bam$", "", colnames(count_mat))
storage.mode(count_mat) <- "integer"

# -----------------------------
# 2. Read metadata
# -----------------------------
meta <- read.delim(
  "00.metadata/bulkRNA_samples.tsv",
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rownames(meta) <- meta$sample_id
meta <- meta[colnames(count_mat), , drop = FALSE]

stopifnot(identical(colnames(count_mat), rownames(meta)))

meta$tissue <- factor(meta$tissue, levels = c("muscle", "musk_gland"))
meta$age_group <- factor(meta$age_group, levels = c("9m", "3y"))

# create group label for plotting
meta$group <- paste(meta$tissue, meta$age_group, sep = "_")

# -----------------------------
# 3. Build DESeq2 object
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = meta,
  design = ~ age_group + tissue + age_group:tissue
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

vsd <- vst(dds)

# -----------------------------
# 4. PCA
# -----------------------------
pdf("05.bulkRNA/05.deseq2/PCA.pdf", width = 6, height = 5)
plotPCA(vsd, intgroup = c("tissue", "age_group"))
dev.off()

print(resultsNames(dds))

# -----------------------------
# 5. Differential expression results
# -----------------------------

# tissue effect at age_group = 9m
res_tissue <- results(dds, name = "tissue_musk_gland_vs_muscle")
res_tissue_df <- as.data.frame(res_tissue) %>%
  rownames_to_column("gene")

# interaction
res_inter <- results(dds, name = "age_group3y.tissuemusk_gland")
res_inter_df <- as.data.frame(res_inter) %>%
  rownames_to_column("gene")

# age effect in musk gland
dds_gland <- dds[, dds$tissue == "musk_gland"]
design(dds_gland) <- ~ age_group
dds_gland <- DESeq(dds_gland)

res_gland_age <- results(dds_gland, contrast = c("age_group", "3y", "9m"))
res_gland_age_df <- as.data.frame(res_gland_age) %>%
  rownames_to_column("gene")

# age effect in muscle
dds_muscle <- dds[, dds$tissue == "muscle"]
design(dds_muscle) <- ~ age_group
dds_muscle <- DESeq(dds_muscle)

res_muscle_age <- results(dds_muscle, contrast = c("age_group", "3y", "9m"))
res_muscle_age_df <- as.data.frame(res_muscle_age) %>%
  rownames_to_column("gene")

# -----------------------------
# 6. Add DEG labels
# -----------------------------
add_deg_label <- function(df, padj_cutoff = 0.05, lfc_cutoff = 1) {
  df %>%
    mutate(DEG = case_when(
      !is.na(padj) & padj < padj_cutoff & log2FoldChange >= lfc_cutoff  ~ "UP",
      !is.na(padj) & padj < padj_cutoff & log2FoldChange <= -lfc_cutoff ~ "DOWN",
      TRUE ~ "NS"
    ))
}

res_tissue_df <- add_deg_label(res_tissue_df)
res_inter_df <- add_deg_label(res_inter_df)
res_gland_age_df <- add_deg_label(res_gland_age_df)
res_muscle_age_df <- add_deg_label(res_muscle_age_df)

# -----------------------------
# 7. Write all result tables
# -----------------------------
write.table(
  res_tissue_df,
  "05.bulkRNA/05.deseq2/res_tissue.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  res_inter_df,
  "05.bulkRNA/05.deseq2/res_interaction.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  res_gland_age_df,
  "05.bulkRNA/05.deseq2/res_gland_age.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  res_muscle_age_df,
  "05.bulkRNA/05.deseq2/res_muscle_age.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# DEG-only tables
write.table(
  res_tissue_df %>% filter(DEG != "NS"),
  "05.bulkRNA/05.deseq2/deg_tissue.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  res_inter_df %>% filter(DEG != "NS"),
  "05.bulkRNA/05.deseq2/deg_interaction.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  res_gland_age_df %>% filter(DEG != "NS"),
  "05.bulkRNA/05.deseq2/deg_gland_age.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  res_muscle_age_df %>% filter(DEG != "NS"),
  "05.bulkRNA/05.deseq2/deg_muscle_age.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# -----------------------------
# 8. Volcano plots
# -----------------------------
plot_volcano <- function(
  df,
  outfile,
  padj_cutoff = 0.05,
  lfc_cutoff = 1,
  raster_dpi = 300,
  width = 8,
  height = 8
) {
  plot_df <- df %>%
    mutate(
      padj_plot = case_when(
        is.na(padj) ~ NA_real_,
        padj <= 0 ~ .Machine$double.xmin,
        TRUE ~ padj
      ),
      neg_log10_padj = -log10(padj_plot),
      Regulate = case_when(
        !is.na(padj) & padj < padj_cutoff & log2FoldChange <= -lfc_cutoff ~ "down-regulated",
        !is.na(padj) & padj < padj_cutoff & log2FoldChange >= lfc_cutoff ~ "up-regulated",
        TRUE ~ "unchanged"
      ),
      Regulate = factor(
        Regulate,
        levels = c("down-regulated", "unchanged", "up-regulated")
      )
    ) %>%
    filter(is.finite(log2FoldChange), is.finite(neg_log10_padj))

  x_limit <- ceiling(max(abs(plot_df$log2FoldChange), lfc_cutoff, na.rm = TRUE))
  y_limit <- ceiling(max(plot_df$neg_log10_padj, -log10(padj_cutoff), na.rm = TRUE))
  p_threshold <- -log10(padj_cutoff)

  p <- ggplot(plot_df, aes(x = log2FoldChange, y = neg_log10_padj, color = Regulate)) +
    ggrastr::geom_point_rast(
      size = 1.45,
      alpha = 0.95,
      stroke = 0,
      raster.dpi = raster_dpi
    ) +
    geom_vline(
      xintercept = c(-lfc_cutoff, lfc_cutoff),
      linetype = "dashed",
      size = 0.5,
      color = "black"
    ) +
    geom_hline(
      yintercept = p_threshold,
      linetype = "dashed",
      size = 0.5,
      color = "black"
    ) +
    annotate(
      "label",
      x = -lfc_cutoff,
      y = y_limit * 0.985,
      label = paste0("log2FC = -", lfc_cutoff),
      hjust = 1.04,
      vjust = 1,
      size = 4.2,
      label.size = 0,
      fill = "white",
      color = "black"
    ) +
    annotate(
      "label",
      x = lfc_cutoff,
      y = y_limit * 0.985,
      label = paste0("log2FC = ", lfc_cutoff),
      hjust = -0.04,
      vjust = 1,
      size = 4.2,
      label.size = 0,
      fill = "white",
      color = "black"
    ) +
    annotate(
      "label",
      x = -x_limit * 0.985,
      y = p_threshold,
      label = paste0("padj = ", padj_cutoff),
      hjust = 0,
      vjust = -0.45,
      size = 4.2,
      label.size = 0,
      fill = "white",
      color = "black"
    ) +
    scale_color_manual(
      name = "Regulate",
      values = c(
        "down-regulated" = "#00A6B2",
        "unchanged" = "#9E9E9E",
        "up-regulated" = "#FB4B0B"
      ),
      drop = FALSE
    ) +
    scale_x_continuous(
      limits = c(-x_limit, x_limit),
      breaks = pretty(c(-x_limit, x_limit), n = 5),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      limits = c(0, y_limit * 1.04),
      breaks = pretty(c(0, y_limit), n = 5),
      expand = expansion(mult = c(0, 0.02))
    ) +
    guides(color = guide_legend(override.aes = list(size = 4.2, alpha = 1))) +
    labs(x = "log2FC", y = "-Log10P.Value") +
    theme_bw(base_size = 18) +
    theme(
      panel.grid.major = element_line(color = "#E9E9E9", size = 0.45),
      panel.grid.minor = element_line(color = "#F2F2F2", size = 0.35),
      panel.border = element_rect(color = "#333333", fill = NA, size = 0.55),
      axis.title = element_text(color = "black", size = 18),
      axis.text = element_text(color = "#3A3A3A", size = 14),
      legend.position = c(0.81, 0.96),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.title = element_text(color = "black", size = 17),
      legend.text = element_text(color = "black", size = 14),
      plot.margin = margin(6, 8, 6, 6)
    )

  ggsave(outfile, p, width = width, height = height, device = "pdf")
}

plot_volcano(
  res_tissue_df,
  "05.bulkRNA/05.deseq2/volcano_tissue.pdf"
)

plot_volcano(
  res_gland_age_df,
  "05.bulkRNA/05.deseq2/volcano_gland_age.pdf"
)

plot_volcano(
  res_inter_df,
  "05.bulkRNA/05.deseq2/volcano_interaction.pdf"
)

plot_volcano(
  res_inter_df,
  "05.bulkRNA/05.deseq2/Figure6C_volcano_interaction.pdf"
)

plot_volcano(
  res_muscle_age_df,
  "05.bulkRNA/05.deseq2/volcano_muscle_age.pdf"
)

# -----------------------------
# 9. Heatmap: top50 genes from tissue result
# -----------------------------
top50_tissue <- res_tissue_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice(1:min(50, n())) %>%
  pull(gene)

if (length(top50_tissue) >= 2) {
  mat_tissue <- assay(vsd)[top50_tissue, , drop = FALSE]
  mat_tissue <- t(scale(t(mat_tissue)))

  pdf("05.bulkRNA/05.deseq2/top50_tissue_heatmap.pdf", width = 8, height = 8)
  pheatmap(
    mat_tissue,
    annotation_col = meta[, c("tissue", "age_group", "group"), drop = FALSE],
    show_rownames = FALSE,
    cluster_cols = TRUE,
    scale = "none"
  )
  dev.off()
}

# -----------------------------
# 10. Heatmap: top50 genes from gland age result
# -----------------------------
top50_gland_age <- res_gland_age_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice(1:min(50, n())) %>%
  pull(gene)

if (length(top50_gland_age) >= 2) {
  gland_samples <- rownames(meta)[meta$tissue == "musk_gland"]
  mat_gland <- assay(vsd)[top50_gland_age, gland_samples, drop = FALSE]
  ann_gland <- meta[gland_samples, c("tissue", "age_group", "group"), drop = FALSE]
  mat_gland <- t(scale(t(mat_gland)))

  pdf("05.bulkRNA/05.deseq2/top50_gland_age_heatmap.pdf", width = 7, height = 8)
  pheatmap(
    mat_gland,
    annotation_col = ann_gland,
    show_rownames = FALSE,
    cluster_cols = TRUE,
    scale = "none"
  )
  dev.off()
}

# -----------------------------
# 11. Heatmap: top50 genes from muscle age result
# -----------------------------
top50_muscle_age <- res_muscle_age_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice(1:min(50, n())) %>%
  pull(gene)

if (length(top50_muscle_age) >= 2) {
  muscle_samples <- rownames(meta)[meta$tissue == "muscle"]
  mat_muscle <- assay(vsd)[top50_muscle_age, muscle_samples, drop = FALSE]
  ann_muscle <- meta[muscle_samples, c("tissue", "age_group", "group"), drop = FALSE]
  mat_muscle <- t(scale(t(mat_muscle)))

  pdf("05.bulkRNA/05.deseq2/top50_muscle_age_heatmap.pdf", width = 7, height = 8)
  pheatmap(
    mat_muscle,
    annotation_col = ann_muscle,
    show_rownames = FALSE,
    cluster_cols = TRUE,
    scale = "none"
  )
  dev.off()
}

# -----------------------------
# 12. DEG summary
# -----------------------------
deg_summary <- tibble(
  comparison = c("tissue", "interaction", "gland_age", "muscle_age"),
  up = c(
    sum(res_tissue_df$DEG == "UP"),
    sum(res_inter_df$DEG == "UP"),
    sum(res_gland_age_df$DEG == "UP"),
    sum(res_muscle_age_df$DEG == "UP")
  ),
  down = c(
    sum(res_tissue_df$DEG == "DOWN"),
    sum(res_inter_df$DEG == "DOWN"),
    sum(res_gland_age_df$DEG == "DOWN"),
    sum(res_muscle_age_df$DEG == "DOWN")
  )
)

write.table(
  deg_summary,
  "05.bulkRNA/05.deseq2/DEG_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(deg_summary)
