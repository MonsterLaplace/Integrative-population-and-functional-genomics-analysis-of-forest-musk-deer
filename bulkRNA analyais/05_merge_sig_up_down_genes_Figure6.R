library(tidyverse)

# Merge significant UP/DOWN genes from:
# 1) musk gland vs muscle tissue DE
# 2) musk gland age-related DE
# 3) tissue-by-age interaction DE

padj_cutoff <- 0.05
lfc_cutoff <- 1

find_deseq2_dir <- function() {
  candidates <- c(
    "05.bulkRNA/05.deseq2",
    "outputs/05.deseq2",
    "05.deseq2"
  )

  existing <- candidates[dir.exists(candidates)]
  if (length(existing) == 0) {
    stop("Cannot find DESeq2 result directory.")
  }

  existing[[1]]
}

classify_deg <- function(df, comparison, padj_cutoff = 0.05, lfc_cutoff = 1) {
  df %>%
    mutate(
      comparison = comparison,
      direction = case_when(
        !is.na(padj) & padj < padj_cutoff & log2FoldChange >= lfc_cutoff ~ "UP",
        !is.na(padj) & padj < padj_cutoff & log2FoldChange <= -lfc_cutoff ~ "DOWN",
        TRUE ~ "NS"
      )
    ) %>%
    filter(direction != "NS") %>%
    select(gene, comparison, direction, log2FoldChange, padj, baseMean)
}

deseq2_dir <- find_deseq2_dir()
out_dir <- file.path(deseq2_dir, "merged_sig_genes")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

result_files <- tribble(
  ~comparison, ~file,
  "tissue_musk_gland_vs_muscle", "res_tissue.tsv",
  "musk_gland_age_3y_vs_9m", "res_gland_age.tsv",
  "tissue_by_age_interaction", "res_interaction.tsv"
)

deg_long <- result_files %>%
  mutate(path = file.path(deseq2_dir, file)) %>%
  pmap_dfr(function(comparison, file, path) {
    if (!file.exists(path)) {
      stop("Cannot find required result file: ", path)
    }
    read.delim(path, sep = "\t", stringsAsFactors = FALSE) %>%
      classify_deg(comparison, padj_cutoff = padj_cutoff, lfc_cutoff = lfc_cutoff)
  })

up_detail <- deg_long %>%
  filter(direction == "UP") %>%
  arrange(gene, comparison)

down_detail <- deg_long %>%
  filter(direction == "DOWN") %>%
  arrange(gene, comparison)

up_ids <- up_detail %>%
  distinct(gene) %>%
  arrange(gene)

down_ids <- down_detail %>%
  distinct(gene) %>%
  arrange(gene)

up_summary <- up_detail %>%
  group_by(gene) %>%
  summarise(
    n_comparisons = n_distinct(comparison),
    comparisons = paste(sort(unique(comparison)), collapse = ";"),
    max_log2FoldChange = max(log2FoldChange, na.rm = TRUE),
    min_padj = min(padj, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_comparisons), min_padj, gene)

down_summary <- down_detail %>%
  group_by(gene) %>%
  summarise(
    n_comparisons = n_distinct(comparison),
    comparisons = paste(sort(unique(comparison)), collapse = ";"),
    min_log2FoldChange = min(log2FoldChange, na.rm = TRUE),
    min_padj = min(padj, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_comparisons), min_padj, gene)

comparison_counts <- deg_long %>%
  count(comparison, direction) %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  arrange(comparison)

merged_counts <- tibble(
  set = c("all_significant_UP_gene_union", "all_significant_DOWN_gene_union"),
  n_genes = c(nrow(up_ids), nrow(down_ids))
)

write.table(
  up_ids,
  file.path(out_dir, "all_significant_UP_gene_ids.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  down_ids,
  file.path(out_dir, "all_significant_DOWN_gene_ids.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  up_summary,
  file.path(out_dir, "all_significant_UP_gene_ids_with_sources.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  down_summary,
  file.path(out_dir, "all_significant_DOWN_gene_ids_with_sources.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  deg_long,
  file.path(out_dir, "all_significant_UP_DOWN_gene_long_table.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  comparison_counts,
  file.path(out_dir, "significant_gene_counts_by_comparison.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  merged_counts,
  file.path(out_dir, "merged_significant_gene_union_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("DESeq2 directory: ", deseq2_dir)
message("Output directory: ", out_dir)
message("UP gene union: ", nrow(up_ids))
message("DOWN gene union: ", nrow(down_ids))
