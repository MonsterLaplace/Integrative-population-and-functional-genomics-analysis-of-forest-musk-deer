suppressPackageStartupMessages({
  library(data.table)
})

outDir <- "06.multiome/04.seurat_signac"
muskFile <- file.path(outDir, "musk_function_genes.tsv")
annoFile <- file.path(outDir, "FMdeer_unified_gene_annotation.unique.tsv")

bulkDir <- "05.bulkRNA/05.deseq2"
glandFile <- file.path(bulkDir, "res_gland_age.tsv")
tissueFile <- file.path(bulkDir, "res_tissue.tsv")
interFile <- file.path(bulkDir, "res_interaction.tsv")

if (!file.exists(muskFile)) stop("Missing: ", muskFile)
if (!file.exists(annoFile)) stop("Missing: ", annoFile)
if (!file.exists(glandFile)) stop("Missing: ", glandFile)
if (!file.exists(tissueFile)) stop("Missing: ", tissueFile)
if (!file.exists(interFile)) stop("Missing: ", interFile)

musk <- fread(muskFile)
anno <- fread(annoFile)

# --------------------------------------------------
# check annotation columns
# --------------------------------------------------
need_anno <- c("gene_id", "gene_id_dash", "gene_name")
miss_anno <- setdiff(need_anno, colnames(anno))
if (length(miss_anno) > 0) {
  stop("Annotation missing columns: ", paste(miss_anno, collapse = ", "))
}

anno_map <- unique(anno[, .(gene_id, gene_id_dash, gene_name)])
anno_map <- anno_map[!is.na(gene_id) & gene_id != ""]
anno_map <- anno_map[order(gene_id), .SD[1], by = gene_id]

# --------------------------------------------------
# helper: read and standardize DESeq2 result
# --------------------------------------------------
read_bulk_res <- function(file, prefix) {
  dt <- fread(file)

  if (!all(c("gene", "log2FoldChange", "padj") %in% colnames(dt))) {
    stop("Bulk file missing required columns: ", file)
  }

  dt <- unique(dt[, .(
    gene,
    log2FoldChange,
    padj
  )])

  setnames(
    dt,
    c("log2FoldChange", "padj"),
    c(paste0(prefix, "_log2FC"), paste0(prefix, "_padj"))
  )

  dt
}

# --------------------------------------------------
# helper: expand bulk genes to multiple match keys
# KEEP ONLY:
#   match_key + value_cols
# do NOT keep original gene column, to avoid merge duplication
# --------------------------------------------------
expand_bulk_keys <- function(bulk_dt, anno_map, value_cols) {
  stopifnot("gene" %in% colnames(bulk_dt))

  x <- merge(
    bulk_dt,
    anno_map,
    by.x = "gene",
    by.y = "gene_id",
    all.x = TRUE,
    sort = FALSE
  )

  # route 1: original bulk gene_id
  key1 <- x[, c("gene", value_cols), with = FALSE]
  setnames(key1, "gene", "match_key")

  # route 2: dashed gene_id
  key2 <- x[!is.na(gene_id_dash) & gene_id_dash != "", c("gene_id_dash", value_cols), with = FALSE]
  setnames(key2, "gene_id_dash", "match_key")

  # route 3: gene_name
  key3 <- x[!is.na(gene_name) & gene_name != "", c("gene_name", value_cols), with = FALSE]
  setnames(key3, "gene_name", "match_key")

  expanded <- rbindlist(list(key1, key2, key3), use.names = TRUE, fill = TRUE)

  expanded <- expanded[!is.na(match_key) & match_key != ""]
  setorder(expanded, match_key)
  expanded <- expanded[, .SD[1], by = match_key]

  expanded
}

# --------------------------------------------------
# helper: merge expanded bulk into musk
# --------------------------------------------------
merge_expanded_bulk <- function(musk, expanded_bulk, value_cols, source_name = NULL) {
  res <- copy(musk)

  # 删除旧列，避免重复
  old_cols <- intersect(c(value_cols, if (!is.null(source_name)) "bulkDE_source"), colnames(res))
  if (length(old_cols) > 0) {
    res[, (old_cols) := NULL]
  }

  # expanded_bulk 应该只有 match_key + value_cols
  keep_cols <- unique(c("match_key", value_cols))
  expanded_bulk2 <- copy(expanded_bulk[, ..keep_cols])

  out <- merge(
    res,
    expanded_bulk2,
    by.x = "gene",
    by.y = "match_key",
    all.x = TRUE,
    sort = FALSE
  )

  if (!is.null(source_name)) {
    out[, bulkDE_source := source_name]
  }

  out
}

# --------------------------------------------------
# read bulk result tables
# --------------------------------------------------
bulk_gland <- read_bulk_res(glandFile, "bulkDE")
bulk_tissue <- read_bulk_res(tissueFile, "bulk_tissue")
bulk_inter <- read_bulk_res(interFile, "bulk_interaction")

# --------------------------------------------------
# expand keys
# --------------------------------------------------
bulk_gland_exp <- expand_bulk_keys(
  bulk_dt = bulk_gland,
  anno_map = anno_map,
  value_cols = c("bulkDE_log2FC", "bulkDE_padj")
)

bulk_tissue_exp <- expand_bulk_keys(
  bulk_dt = bulk_tissue,
  anno_map = anno_map,
  value_cols = c("bulk_tissue_log2FC", "bulk_tissue_padj")
)

bulk_inter_exp <- expand_bulk_keys(
  bulk_dt = bulk_inter,
  anno_map = anno_map,
  value_cols = c("bulk_interaction_log2FC", "bulk_interaction_padj")
)

# 保存展开后的表，方便检查
fwrite(bulk_gland_exp, file.path(outDir, "bulk_gland_age.expanded_keys.tsv"), sep = "\t")
fwrite(bulk_tissue_exp, file.path(outDir, "bulk_tissue.expanded_keys.tsv"), sep = "\t")
fwrite(bulk_inter_exp, file.path(outDir, "bulk_interaction.expanded_keys.tsv"), sep = "\t")

# --------------------------------------------------
# merge into musk table
# --------------------------------------------------
musk2 <- merge_expanded_bulk(
  musk = musk,
  expanded_bulk = bulk_gland_exp,
  value_cols = c("bulkDE_log2FC", "bulkDE_padj"),
  source_name = "res_gland_age"
)

musk2 <- merge_expanded_bulk(
  musk = musk2,
  expanded_bulk = bulk_tissue_exp,
  value_cols = c("bulk_tissue_log2FC", "bulk_tissue_padj")
)

musk2 <- merge_expanded_bulk(
  musk = musk2,
  expanded_bulk = bulk_inter_exp,
  value_cols = c("bulk_interaction_log2FC", "bulk_interaction_padj")
)

# --------------------------------------------------
# save updated table
# --------------------------------------------------
fwrite(musk2, muskFile, sep = "\t", quote = FALSE)

outfile2 <- file.path(outDir, "musk_function_genes.with_bulk_deseq2.tsv")
fwrite(musk2, outfile2, sep = "\t", quote = FALSE)

message("Done.")
message("Updated main file: ", muskFile)
message("Saved backup file: ", outfile2)

# --------------------------------------------------
# integration summary
# --------------------------------------------------
summary_dt <- data.table(
  source = c("gland_age", "tissue", "interaction"),
  n_nonNA_log2FC = c(
    sum(!is.na(musk2$bulkDE_log2FC)),
    sum(!is.na(musk2$bulk_tissue_log2FC)),
    sum(!is.na(musk2$bulk_interaction_log2FC))
  ),
  n_nonNA_padj = c(
    sum(!is.na(musk2$bulkDE_padj)),
    sum(!is.na(musk2$bulk_tissue_padj)),
    sum(!is.na(musk2$bulk_interaction_padj))
  )
)

print(summary_dt)
fwrite(summary_dt, file.path(outDir, "bulk_deseq2_integration_summary.tsv"), sep = "\t")
