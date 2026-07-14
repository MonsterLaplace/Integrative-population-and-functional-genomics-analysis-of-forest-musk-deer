suppressPackageStartupMessages({
  library(data.table)
})

outDir <- "06.multiome/04.seurat_signac"

expandedP2GFile <- file.path(outDir, "peak2gene_links.expanded.tsv")
muskFile <- file.path(outDir, "musk_function_genes.tsv")

assignDir <- file.path(outDir, "musk_function_gene_assignment")
assignFile <- file.path(assignDir, "musk_function_genes.with_celltype_assignment.tsv")
top20File <- file.path(assignDir, "top20_musk_function_genes_by_celltype.tsv")

reportDir <- file.path(outDir, "peak2gene_qc")
dir.create(reportDir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(expandedP2GFile)) stop("Missing file: ", expandedP2GFile)
if (!file.exists(muskFile)) stop("Missing file: ", muskFile)
if (!file.exists(assignFile)) stop("Missing file: ", assignFile)
if (!file.exists(top20File)) stop("Missing file: ", top20File)

p2g <- fread(expandedP2GFile)
musk <- fread(muskFile)
assign_dt <- fread(assignFile)
top20 <- fread(top20File)

# --------------------------------------------------
# helper
# --------------------------------------------------
calc_recovery <- function(dt, table_name) {
  if (!"gene" %in% colnames(dt)) stop(table_name, " missing gene column")
  if (!"peak2gene_support" %in% colnames(dt)) stop(table_name, " missing peak2gene_support column")

  data.table(
    table_name = table_name,
    n_rows = nrow(dt),
    n_unique_genes = uniqueN(dt$gene),
    n_supported_rows = sum(dt$peak2gene_support == 1, na.rm = TRUE),
    frac_supported_rows = sum(dt$peak2gene_support == 1, na.rm = TRUE) / nrow(dt),
    n_supported_genes = uniqueN(dt[peak2gene_support == 1, gene]),
    frac_supported_genes = uniqueN(dt[peak2gene_support == 1, gene]) / uniqueN(dt$gene)
  )
}

# =========================================================
# 1. basic summary of expanded peak2gene
# =========================================================
p2g_summary <- data.table(
  metric = c(
    "n_rows_expanded_peak2gene",
    "n_unique_match_keys",
    "n_supported_match_keys"
  ),
  value = c(
    nrow(p2g),
    if ("match_key" %in% colnames(p2g)) uniqueN(p2g$match_key) else NA_integer_,
    if ("match_key" %in% colnames(p2g) && "peak2gene_support" %in% colnames(p2g)) {
      uniqueN(p2g[peak2gene_support == 1, match_key])
    } else {
      NA_integer_
    }
  )
)

fwrite(p2g_summary, file.path(reportDir, "peak2gene_expanded_summary.tsv"), sep = "\t")

# =========================================================
# 2. overall recovery in target tables
# =========================================================
rec_musk <- calc_recovery(musk, "musk_function_genes.tsv")
rec_assign <- calc_recovery(assign_dt, "musk_function_genes.with_celltype_assignment.tsv")
rec_top20 <- calc_recovery(top20, "top20_musk_function_genes_by_celltype.tsv")

recovery_summary <- rbindlist(list(rec_musk, rec_assign, rec_top20), fill = TRUE)
fwrite(recovery_summary, file.path(reportDir, "peak2gene_recovery_summary.tsv"), sep = "\t")

# =========================================================
# 3. recovery by celltype in top20 table
# =========================================================
if ("best_celltype" %in% colnames(top20)) {
  by_ct <- top20[, .(
    n_rows = .N,
    n_supported_rows = sum(peak2gene_support == 1, na.rm = TRUE),
    frac_supported_rows = sum(peak2gene_support == 1, na.rm = TRUE) / .N,
    n_unique_genes = uniqueN(gene),
    n_supported_genes = uniqueN(gene[peak2gene_support == 1]),
    frac_supported_genes = uniqueN(gene[peak2gene_support == 1]) / uniqueN(gene)
  ), by = best_celltype]

  setorder(by_ct, -frac_supported_genes, -n_supported_genes)
  fwrite(by_ct, file.path(reportDir, "peak2gene_recovery_by_celltype_top20.tsv"), sep = "\t")
}

# =========================================================
# 4. high-confidence subset recovery
# =========================================================
if ("musk_score_raw" %in% colnames(assign_dt)) {
  high_conf <- assign_dt[musk_score_raw >= 2]
  if (nrow(high_conf) > 0) {
    rec_high <- calc_recovery(high_conf, "high_confidence_musk_score_ge_2")
    fwrite(rec_high, file.path(reportDir, "peak2gene_recovery_high_confidence.tsv"), sep = "\t")
  }
}

# =========================================================
# 5. adult-up subset recovery
# =========================================================
if (all(c("bulkDE_log2FC", "bulkDE_padj") %in% colnames(assign_dt))) {
  adult_up <- assign_dt[
    !is.na(bulkDE_log2FC) & !is.na(bulkDE_padj) &
      bulkDE_log2FC > 0 & bulkDE_padj < 0.05
  ]

  if (nrow(adult_up) > 0) {
    rec_adult <- calc_recovery(adult_up, "adult_up_genes")
    fwrite(rec_adult, file.path(reportDir, "peak2gene_recovery_adult_up.tsv"), sep = "\t")
  }
}

# =========================================================
# 6. strict candidate subset recovery
# =========================================================
if (all(c("bulkDE_log2FC", "bulkDE_padj", "scRNA_marker", "celltype_specificity") %in% colnames(assign_dt))) {
  strict_dt <- assign_dt[
    !is.na(bulkDE_log2FC) & !is.na(bulkDE_padj) &
      bulkDE_log2FC > 0 & bulkDE_padj < 0.05 &
      scRNA_marker == 1 &
      celltype_specificity >= 0.3
  ]

  if (nrow(strict_dt) > 0) {
    rec_strict <- calc_recovery(strict_dt, "strict_candidates")
    fwrite(rec_strict, file.path(reportDir, "peak2gene_recovery_strict_candidates.tsv"), sep = "\t")
  }
}

message("Done. QC reports saved to: ", reportDir)
