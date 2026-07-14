suppressPackageStartupMessages({
  library(data.table)
})

outDir <- "06.multiome/04.seurat_signac"

expandedP2GFile <- file.path(outDir, "peak2gene_links.expanded.tsv")
muskFile <- file.path(outDir, "musk_function_genes.tsv")

assignDir <- file.path(outDir, "musk_function_gene_assignment")
assignFile <- file.path(assignDir, "musk_function_genes.with_celltype_assignment.tsv")
top20File <- file.path(assignDir, "top20_musk_function_genes_by_celltype.tsv")

if (!file.exists(expandedP2GFile)) stop("Missing file: ", expandedP2GFile)
if (!file.exists(muskFile)) stop("Missing file: ", muskFile)
if (!file.exists(assignFile)) stop("Missing file: ", assignFile)
if (!file.exists(top20File)) stop("Missing file: ", top20File)

p2g <- fread(expandedP2GFile)

if (!all(c("match_key", "peak2gene_support") %in% colnames(p2g))) {
  stop("Expanded peak2gene file missing required columns.")
}

# --------------------------------------------------
# helper: drop old peak2gene-related columns
# --------------------------------------------------
drop_old_p2g_cols <- function(dt) {
  old_cols <- c(
    "peak2gene_support",
    "peak2gene_source",
    "peak2gene_n_links",
    "peak2gene_max_score",
    "match_type"
  )

  old_cols_present <- intersect(old_cols, colnames(dt))
  if (length(old_cols_present) > 0) {
    dt[, (old_cols_present) := NULL]
  }

  # 也清理可能残留的 .x/.y 列
  extra_cols <- grep("^peak2gene_.*\\.[xy]$|^match_type\\.[xy]$", colnames(dt), value = TRUE)
  if (length(extra_cols) > 0) {
    dt[, (extra_cols) := NULL]
  }

  dt
}

# --------------------------------------------------
# helper: merge expanded peak2gene support
# --------------------------------------------------
merge_real_p2g <- function(dt, p2g) {
  if (!"gene" %in% colnames(dt)) stop("Input table missing gene column")

  dt <- copy(dt)
  dt <- drop_old_p2g_cols(dt)

  out <- merge(
    dt,
    p2g,
    by.x = "gene",
    by.y = "match_key",
    all.x = TRUE,
    sort = FALSE
  )

  if (!"peak2gene_support" %in% colnames(out)) {
    out[, peak2gene_support := 0L]
  } else {
    out[is.na(peak2gene_support), peak2gene_support := 0L]
  }

  # 可选：其余信息允许为空
  if ("peak2gene_n_links" %in% colnames(out)) {
    out[is.na(peak2gene_n_links), peak2gene_n_links := 0L]
  }
  if ("peak2gene_max_score" %in% colnames(out)) {
    out[is.na(peak2gene_max_score), peak2gene_max_score := NA_real_]
  }

  out
}

# =========================================================
# 1. musk_function_genes.tsv
# =========================================================
musk <- fread(muskFile)
musk2 <- merge_real_p2g(musk, p2g)
fwrite(musk2, muskFile, sep = "\t", quote = FALSE)
message("Updated: ", muskFile)

# =========================================================
# 2. with_celltype_assignment
# =========================================================
assign_dt <- fread(assignFile)
assign_dt2 <- merge_real_p2g(assign_dt, p2g)
fwrite(assign_dt2, assignFile, sep = "\t", quote = FALSE)
message("Updated: ", assignFile)

# =========================================================
# 3. top20 by celltype
# =========================================================
top20 <- fread(top20File)
top20_2 <- merge_real_p2g(top20, p2g)
fwrite(top20_2, top20File, sep = "\t", quote = FALSE)
message("Updated: ", top20File)

message("Done.")
