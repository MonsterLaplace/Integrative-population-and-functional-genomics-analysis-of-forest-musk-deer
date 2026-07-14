suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

outDir <- "06.multiome/04.seurat_signac"
args <- commandArgs(trailingOnly = TRUE)
analysis_level <- if (length(args) >= 1) args[[1]] else "major"
if (!analysis_level %in% c("major", "epithelial_minor", "cell_state")) {
  stop("analysis_level must be one of: major, epithelial_minor, cell_state")
}
inFile <- file.path(
  outDir,
  if (analysis_level == "major") "musk_function_genes.tsv" else paste0("musk_function_genes_", analysis_level, ".tsv")
)

if (!file.exists(inFile)) {
  stop("Input file not found: ", inFile)
}

plotDir <- file.path(
  outDir,
  if (analysis_level == "major") "musk_function_gene_plots" else paste0("musk_function_gene_plots_", analysis_level)
)
dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(inFile, sep = "\t", header = TRUE)

if (!all(c("gene", "musk_score_raw") %in% colnames(dt))) {
  stop("musk_function_genes.tsv must contain at least gene and musk_score_raw")
}

# =========================================================
# 1. add helper columns
# =========================================================
if (!"scRNA_marker" %in% colnames(dt)) dt[, scRNA_marker := 0]
if (!"bulkDE_log2FC" %in% colnames(dt)) dt[, bulkDE_log2FC := NA_real_]
if (!"bulkDE_padj" %in% colnames(dt)) dt[, bulkDE_padj := NA_real_]
if (!"celltype_specificity" %in% colnames(dt)) dt[, celltype_specificity := 0]
if (!"ATAC_peak_support" %in% colnames(dt)) dt[, ATAC_peak_support := 0]
if (!"peak2gene_support" %in% colnames(dt)) dt[, peak2gene_support := 0]
if (!"motif_support" %in% colnames(dt)) dt[, motif_support := 0]
if (!"selected_analysis_level" %in% colnames(dt)) dt[, selected_analysis_level := analysis_level]
if (!"bulkDE_annotation" %in% colnames(dt)) dt[, bulkDE_annotation := NA_character_]
if (!"bulkDE_annotation_level" %in% colnames(dt)) dt[, bulkDE_annotation_level := NA_character_]

# adult-up support
dt[, adult_up := ifelse(!is.na(bulkDE_log2FC) & !is.na(bulkDE_padj) &
                          bulkDE_log2FC > 0 & bulkDE_padj < 0.05, 1L, 0L)]

# candidate class
dt[, candidate_class := fifelse(
  adult_up == 1 & scRNA_marker == 1 & ATAC_peak_support == 1,
  "adult_up_scRNA_ATAC",
  fifelse(
    adult_up == 1 & scRNA_marker == 1,
    "adult_up_scRNA",
    fifelse(
      adult_up == 1,
      "adult_up_only",
      "other"
    )
  )
)]

# sort
setorder(dt, -musk_score_raw, -celltype_specificity, -bulkDE_log2FC)

# =========================================================
# 2. save ranked full table
# =========================================================
fwrite(dt, file.path(plotDir, "musk_function_genes.ranked.tsv"), sep = "\t")

# high-confidence subset
high_conf <- dt[
  musk_score_raw >= 2 &
    adult_up == 1 &
    scRNA_marker == 1
]

fwrite(high_conf, file.path(plotDir, "musk_function_genes.high_confidence.tsv"), sep = "\t")

# top 100
top100 <- dt[1:min(.N, 100)]
fwrite(top100, file.path(plotDir, "musk_function_genes.top100.tsv"), sep = "\t")

# =========================================================
# 3. Top 30 barplot
# =========================================================
topN <- min(30, nrow(dt))
top30 <- copy(dt[1:topN])
top30[, gene := factor(gene, levels = rev(gene))]

p_bar <- ggplot(top30, aes(x = gene, y = musk_score_raw, fill = candidate_class)) +
  geom_col(width = 0.8) +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    x = "Gene",
    y = "musk_score_raw",
    title = paste("Top candidate musk gland function genes:", analysis_level)
  ) +
  theme(
    panel.grid = element_blank()
  )

ggsave(
  filename = file.path(plotDir, "01_top30_musk_function_genes_barplot.pdf"),
  plot = p_bar,
  width = 9,
  height = 8
)

# =========================================================
# 4. bulkDE vs specificity scatter
# only genes with bulk DE available
# =========================================================
scatter_dt <- dt[!is.na(bulkDE_log2FC)]

if (nrow(scatter_dt) > 0) {
  scatter_dt[, label_flag := musk_score_raw >= quantile(musk_score_raw, 0.99, na.rm = TRUE)]

  p_scatter <- ggplot(scatter_dt, aes(
    x = bulkDE_log2FC,
    y = celltype_specificity,
    color = musk_score_raw,
    shape = factor(scRNA_marker)
  )) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_gradient(low = "grey70", high = "firebrick") +
    theme_bw(base_size = 12) +
    labs(
      x = "bulkDE_log2FC (adult vs young)",
      y = "celltype specificity",
      color = "musk score",
      shape = "scRNA marker",
      title = paste("Bulk DE vs annotation specificity:", analysis_level)
    )

  ggsave(
    filename = file.path(plotDir, "02_bulkDE_vs_specificity_scatter.pdf"),
    plot = p_scatter,
    width = 8,
    height = 6
  )
}

# =========================================================
# 5. evidence heatmap-like plot for top 50 genes
# =========================================================
top50 <- copy(dt[1:min(.N, 50)])

heat_dt <- rbindlist(list(
  data.table(gene = top50$gene, evidence = "scRNA_marker", value = top50$scRNA_marker),
  data.table(gene = top50$gene, evidence = "adult_up_bulkDE", value = top50$adult_up),
  data.table(gene = top50$gene, evidence = "celltype_specificity", value = top50$celltype_specificity),
  data.table(gene = top50$gene, evidence = "ATAC_peak_support", value = top50$ATAC_peak_support),
  data.table(gene = top50$gene, evidence = "peak2gene_support", value = top50$peak2gene_support),
  data.table(gene = top50$gene, evidence = "motif_support", value = top50$motif_support)
))

heat_dt[, gene := factor(gene, levels = rev(top50$gene))]
heat_dt[, evidence := factor(evidence, levels = c(
  "scRNA_marker",
  "adult_up_bulkDE",
  "celltype_specificity",
  "ATAC_peak_support",
  "peak2gene_support",
  "motif_support"
))]

p_heat <- ggplot(heat_dt, aes(x = evidence, y = gene, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_bw(base_size = 12) +
  labs(
    x = "Evidence",
    y = "Gene",
    fill = "Value",
    title = paste("Integrated evidence for top 50 candidate genes:", analysis_level)
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(plotDir, "03_top50_evidence_heatmap.pdf"),
  plot = p_heat,
  width = 10,
  height = 12
)

# =========================================================
# 6. adult-up candidate summary
# =========================================================
adult_up_dt <- dt[adult_up == 1]
fwrite(adult_up_dt, file.path(plotDir, "musk_function_genes.adult_up.tsv"), sep = "\t")

adult_up_strict <- dt[
  adult_up == 1 &
    scRNA_marker == 1 &
    celltype_specificity >= 0.3
]
fwrite(adult_up_strict, file.path(plotDir, "musk_function_genes.adult_up_strict.tsv"), sep = "\t")

# =========================================================
# 7. summary counts
# =========================================================
summary_dt <- data.table(
  metric = c(
    "n_total_genes",
    "n_scRNA_marker",
    "n_adult_up_bulkDE",
    "n_ATAC_supported",
    "n_peak2gene_supported",
    "n_motif_supported",
    "n_high_confidence",
    "n_adult_up_strict"
  ),
  value = c(
    nrow(dt),
    sum(dt$scRNA_marker == 1, na.rm = TRUE),
    sum(dt$adult_up == 1, na.rm = TRUE),
    sum(dt$ATAC_peak_support == 1, na.rm = TRUE),
    sum(dt$peak2gene_support == 1, na.rm = TRUE),
    sum(dt$motif_support == 1, na.rm = TRUE),
    nrow(high_conf),
    nrow(adult_up_strict)
  )
)
summary_dt[, selected_analysis_level := analysis_level]

fwrite(summary_dt, file.path(plotDir, "musk_function_genes.summary.tsv"), sep = "\t")

message("Done. Outputs saved to: ", plotDir)
