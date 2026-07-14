suppressPackageStartupMessages({
  library(data.table)
})

outDir <- "06.multiome/04.seurat_signac"
inFile <- file.path(outDir, "musk_function_genes.tsv")

if (!file.exists(inFile)) {
  stop("Input file not found: ", inFile)
}

dt <- fread(inFile, sep = "\t", header = TRUE)

if (!"gene" %in% colnames(dt)) {
  stop("Input file must contain gene column")
}

# =========================================================
# ensure required columns
# =========================================================
ensure_col <- function(dt, col, default) {
  if (!col %in% colnames(dt)) {
    dt[, (col) := default]
  }
  dt
}

dt <- ensure_col(dt, "scRNA_marker", 0L)
dt <- ensure_col(dt, "bulkDE_log2FC", NA_real_)
dt <- ensure_col(dt, "bulkDE_padj", NA_real_)
dt <- ensure_col(dt, "celltype_specificity", 0)
dt <- ensure_col(dt, "ATAC_peak_support", 0L)
dt <- ensure_col(dt, "peak2gene_support", 0L)
dt <- ensure_col(dt, "motif_support", 0L)

dt <- ensure_col(dt, "bulk_tissue_log2FC", NA_real_)
dt <- ensure_col(dt, "bulk_tissue_padj", NA_real_)
dt <- ensure_col(dt, "bulk_interaction_log2FC", NA_real_)
dt <- ensure_col(dt, "bulk_interaction_padj", NA_real_)

# type coercion
int_cols <- c("scRNA_marker", "ATAC_peak_support", "peak2gene_support", "motif_support")
for (cc in intersect(int_cols, colnames(dt))) {
  dt[is.na(get(cc)), (cc) := 0L]
  dt[, (cc) := as.integer(get(cc))]
}

num_cols <- c(
  "bulkDE_log2FC", "bulkDE_padj", "celltype_specificity",
  "bulk_tissue_log2FC", "bulk_tissue_padj",
  "bulk_interaction_log2FC", "bulk_interaction_padj"
)
for (cc in intersect(num_cols, colnames(dt))) {
  dt[, (cc) := as.numeric(get(cc))]
}

dt[is.na(celltype_specificity), celltype_specificity := 0]

# =========================================================
# 1. mature/adult bulk score
# only adult-up contributes positively
# =========================================================
dt[, mature_bulk_score := 0]

dt[
  !is.na(bulkDE_log2FC) & !is.na(bulkDE_padj) &
    bulkDE_log2FC > 0 & bulkDE_padj < 0.05,
  mature_bulk_score := pmin(bulkDE_log2FC, 3)
]

# =========================================================
# 2. optional tissue bonus
# if tissue effect is also positive and significant
# =========================================================
dt[, tissue_bonus := 0]

dt[
  !is.na(bulk_tissue_log2FC) & !is.na(bulk_tissue_padj) &
    bulk_tissue_log2FC > 0 & bulk_tissue_padj < 0.05,
  tissue_bonus := pmin(bulk_tissue_log2FC, 2)
]

# =========================================================
# 3. optional interaction bonus
# use absolute significant interaction as small bonus
# =========================================================
dt[, interaction_bonus := 0]

dt[
  !is.na(bulk_interaction_log2FC) & !is.na(bulk_interaction_padj) &
    bulk_interaction_padj < 0.05,
  interaction_bonus := pmin(abs(bulk_interaction_log2FC), 2)
]

# =========================================================
# 4. recompute raw score
# weights tuned for mature musk gland candidate prioritization
# =========================================================
dt[, musk_score_raw :=
      1.2 * scRNA_marker +
      1.5 * mature_bulk_score +
      1.5 * celltype_specificity +
      0.5 * ATAC_peak_support +
      1.0 * peak2gene_support +
      0.5 * motif_support +
      0.3 * tissue_bonus +
      0.2 * interaction_bonus
]

# =========================================================
# 5. helper labels
# =========================================================
dt[, adult_up := fifelse(
  !is.na(bulkDE_log2FC) & !is.na(bulkDE_padj) &
    bulkDE_log2FC > 0 & bulkDE_padj < 0.05,
  1L, 0L
)]

dt[, candidate_class := fifelse(
  adult_up == 1 & scRNA_marker == 1 & peak2gene_support == 1 & ATAC_peak_support == 1,
  "adult_up_scRNA_peak2gene_ATAC",
  fifelse(
    adult_up == 1 & scRNA_marker == 1 & peak2gene_support == 1,
    "adult_up_scRNA_peak2gene",
    fifelse(
      adult_up == 1 & scRNA_marker == 1,
      "adult_up_scRNA",
      fifelse(
        adult_up == 1,
        "adult_up_only",
        "other"
      )
    )
  )
)]

# =========================================================
# 6. sort
# =========================================================
setorder(dt, -musk_score_raw, -mature_bulk_score, -celltype_specificity)

# =========================================================
# 7. save updated main table
# =========================================================
fwrite(dt, inFile, sep = "\t", quote = FALSE)

backupFile <- file.path(outDir, "musk_function_genes.re_scored.tsv")
fwrite(dt, backupFile, sep = "\t", quote = FALSE)

# =========================================================
# 8. save useful subsets
# =========================================================
top100 <- dt[1:min(.N, 100)]
fwrite(top100, file.path(outDir, "musk_function_genes.top100.re_scored.tsv"), sep = "\t", quote = FALSE)

high_conf <- dt[
  adult_up == 1 &
    scRNA_marker == 1 &
    celltype_specificity >= 0.3
]
fwrite(high_conf, file.path(outDir, "musk_function_genes.high_confidence.re_scored.tsv"), sep = "\t", quote = FALSE)

strict_conf <- dt[
  adult_up == 1 &
    scRNA_marker == 1 &
    celltype_specificity >= 0.3 &
    peak2gene_support == 1
]
fwrite(strict_conf, file.path(outDir, "musk_function_genes.strict_confidence.re_scored.tsv"), sep = "\t", quote = FALSE)

# =========================================================
# 9. summary
# =========================================================
summary_dt <- data.table(
  metric = c(
    "n_total_genes",
    "n_adult_up",
    "n_scRNA_marker",
    "n_peak2gene_support",
    "n_ATAC_support",
    "n_high_confidence",
    "n_strict_confidence"
  ),
  value = c(
    nrow(dt),
    sum(dt$adult_up == 1, na.rm = TRUE),
    sum(dt$scRNA_marker == 1, na.rm = TRUE),
    sum(dt$peak2gene_support == 1, na.rm = TRUE),
    sum(dt$ATAC_peak_support == 1, na.rm = TRUE),
    nrow(high_conf),
    nrow(strict_conf)
  )
)

fwrite(summary_dt, file.path(outDir, "musk_function_genes.re_scored.summary.tsv"), sep = "\t", quote = FALSE)

message("Done.")
message("Updated: ", inFile)
message("Backup:  ", backupFile)
