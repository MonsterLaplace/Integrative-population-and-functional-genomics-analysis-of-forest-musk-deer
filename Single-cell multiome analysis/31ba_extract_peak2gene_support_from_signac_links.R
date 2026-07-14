suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(data.table)
})

outDir <- "06.multiome/04.seurat_signac"
objFile <- file.path(outDir, "FMdeer_multiome.wnn.final_annotated.linkpeaks.rds")

if (!file.exists(objFile)) {
  stop("Object with links not found: ", objFile)
}

obj <- readRDS(objFile)

# Links were computed across all cells for statistical power. They are global
# regulatory evidence, not celltype-specific links, but retain the annotation
# hierarchy in an audit table for downstream interpretation.
required_metadata <- c("celltype_detailed", "celltype_minor", "celltype_major", "celltype_major_inferred", "cell_state")
if (!all(required_metadata %in% colnames(obj@meta.data))) {
  stop("Annotation-level metadata missing. Run 19_apply_cluster_annotation_from_table.R first.")
}
annotation_context <- rbindlist(lapply(required_metadata, function(cc) {
  x <- as.data.table(table(obj@meta.data[[cc]]))
  setnames(x, c("annotation", "n_cells"))
  x[, annotation_level := cc]
  x
}))
annotation_context[, link_scope := "all_cells"]
fwrite(annotation_context, file.path(outDir, "peak2gene_link_annotation_context.tsv"), sep = "\t")

# --------------------------------------------------
# 1. extract links
# --------------------------------------------------
lnk <- Links(obj)

if (is.null(lnk) || length(lnk) == 0) {
  stop("No peak-gene links found in object.")
}

lnk_df <- as.data.frame(lnk)
lnk_dt <- as.data.table(lnk_df)

# 看看有哪些列
message("Columns in Links(object):")
print(colnames(lnk_dt))

# --------------------------------------------------
# 2. 尽量识别 gene 列和 score 列
# 不同版本 Signac 字段名可能略有差异
# --------------------------------------------------
gene_col <- NULL
for (cc in c("gene", "gene_name", "name", "target")) {
  if (cc %in% colnames(lnk_dt)) {
    gene_col <- cc
    break
  }
}

if (is.null(gene_col)) {
  # 尝试从 metadata 列中找
  stop("Cannot identify gene column in Links(obj). Please inspect column names.")
}

score_col <- NULL
for (cc in c("score", "cor", "correlation")) {
  if (cc %in% colnames(lnk_dt)) {
    score_col <- cc
    break
  }
}

# --------------------------------------------------
# 3. gene-level summary
# --------------------------------------------------
if (!is.null(score_col)) {
  gene_support <- lnk_dt[, .(
    peak2gene_support = 1L,
    peak2gene_n_links = .N,
    peak2gene_max_score = max(get(score_col), na.rm = TRUE)
  ), by = c(gene_col)]

  setnames(gene_support, gene_col, "gene")
} else {
  gene_support <- unique(lnk_dt[, .(tmp_gene = get(gene_col))])
  setnames(gene_support, "tmp_gene", "gene")
  gene_support[, peak2gene_support := 1L]
  gene_support[, peak2gene_n_links := 1L]
  gene_support[, peak2gene_max_score := NA_real_]
}
gene_support[, `:=`(link_scope = "all_cells", annotation_level = "global")]

# --------------------------------------------------
# 4. save
# --------------------------------------------------
fwrite(lnk_dt, file.path(outDir, "peak2gene_links.raw.tsv"), sep = "\t")
fwrite(gene_support, file.path(outDir, "peak2gene_links.tsv"), sep = "\t")

message("Saved:")
message("  ", file.path(outDir, "peak2gene_links.raw.tsv"))
message("  ", file.path(outDir, "peak2gene_links.tsv"))
