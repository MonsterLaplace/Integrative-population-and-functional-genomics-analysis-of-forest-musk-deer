suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(data.table)
  library(BSgenome)
  library(BSgenome.FMdeer.Custom.v1)
})

outDir <- "06.multiome/04.seurat_signac"
inferredFile <- file.path(outDir, "FMdeer_multiome.wnn.final_annotated.state_lineage_inferred.rds")
objFile <- if (file.exists(inferredFile)) inferredFile else file.path(outDir, "FMdeer_multiome.wnn.final_annotated.rds")
annoFile <- file.path(outDir, "FMdeer_unified_gene_annotation.unique.tsv")

if (!file.exists(objFile)) stop("Object not found: ", objFile)
if (!file.exists(annoFile)) stop("Annotation file not found: ", annoFile)

obj <- readRDS(objFile)
anno_dt <- fread(annoFile)

# LinkPeaks is intentionally computed on all cells to preserve power.  The
# object must nevertheless carry the merged annotation hierarchy so downstream
# peak-gene links can be interpreted at major, minor, or state level.
required_metadata <- c("celltype_detailed", "celltype_minor", "celltype_major", "celltype_major_inferred", "cell_state")
if (!all(required_metadata %in% colnames(obj@meta.data))) {
  stop("Annotation-level metadata missing. Run 19_apply_cluster_annotation_from_table.R first.")
}
annotation_counts <- rbindlist(lapply(required_metadata, function(cc) {
  x <- as.data.table(table(obj@meta.data[[cc]]))
  setnames(x, c("annotation", "n_cells"))
  x[, annotation_level := cc]
  x
}))
fwrite(annotation_counts, file.path(outDir, "linkpeaks_annotation_level_counts.tsv"), sep = "\t")

if (!all(c("RNA", "ATAC") %in% Assays(obj))) {
  stop("Object must contain both RNA and ATAC assays.")
}

if (!all(c("gene_id", "gene_id_dash", "gene_name") %in% colnames(anno_dt))) {
  stop("Annotation file must contain gene_id, gene_id_dash, gene_name")
}

# --------------------------------------------------
# 1. get RNA counts
# --------------------------------------------------
rna_counts <- tryCatch({
  GetAssayData(obj, assay = "RNA", layer = "counts")
}, error = function(e) {
  GetAssayData(obj, assay = "RNA", slot = "counts")
})

rna_genes <- rownames(rna_counts)
map_dt <- unique(anno_dt[, .(gene_id, gene_id_dash, gene_name)])

# --------------------------------------------------
# 2. map RNA rownames to link-compatible names
# priority:
#   a) match RNA gene to gene_id_dash
#   b) match RNA gene to gene_id
#   c) match RNA gene to gene_name
# resulting display:
#   prefer gene_name if available, else gene_id
# --------------------------------------------------
rna_map <- data.table(rna_feature = rna_genes)

# match gene_id_dash
m1 <- merge(
  rna_map,
  map_dt,
  by.x = "rna_feature",
  by.y = "gene_id_dash",
  all.x = TRUE,
  sort = FALSE
)

# unmatched -> gene_id
idx1 <- which(is.na(m1$gene_id))
if (length(idx1) > 0) {
  sub1 <- rna_map[idx1]
  m2 <- merge(
    sub1,
    map_dt,
    by.x = "rna_feature",
    by.y = "gene_id",
    all.x = TRUE,
    sort = FALSE
  )
  m1$gene_id[idx1] <- m2$gene_id
  m1$gene_name[idx1] <- m2$gene_name
}

# unmatched -> gene_name
idx2 <- which(is.na(m1$gene_id) & is.na(m1$gene_name))
if (length(idx2) > 0) {
  sub2 <- rna_map[idx2]
  m3 <- merge(
    sub2,
    map_dt,
    by.x = "rna_feature",
    by.y = "gene_name",
    all.x = TRUE,
    sort = FALSE
  )
  m1$gene_id[idx2] <- m3$gene_id
  m1$gene_name[idx2] <- m3$gene_name
}

# target name for LinkPeaks
m1[, link_name := ifelse(!is.na(gene_name) & gene_name != "", gene_name, gene_id)]

# fallback: keep original if still missing
m1[is.na(link_name) | link_name == "", link_name := rna_feature]

# 去重：如果多个 RNA feature 映射到同一个 link_name，保留第一个
m1_unique <- m1[!duplicated(link_name)]

message("RNA features: ", nrow(m1))
message("Unique LINK names: ", nrow(m1_unique))
message("Mapped to non-original names: ", sum(m1_unique$link_name != m1_unique$rna_feature))

# 保存映射表
fwrite(m1_unique, file.path(outDir, "RNA_to_RNA_LINK_mapping.tsv"), sep = "\t")

# --------------------------------------------------
# 3. build RNA_LINK assay
# --------------------------------------------------
rna_counts_link <- rna_counts[m1_unique$rna_feature, , drop = FALSE]
rownames(rna_counts_link) <- m1_unique$link_name

obj[["RNA_LINK"]] <- CreateAssayObject(counts = rna_counts_link)

DefaultAssay(obj) <- "RNA_LINK"
obj <- NormalizeData(obj, verbose = FALSE)

# --------------------------------------------------
# 4. check overlap with ATAC annotation
# --------------------------------------------------
atac_anno <- Annotation(obj[["ATAC"]])

if (is.null(atac_anno)) {
  stop("ATAC annotation is missing.")
}

anno_names <- NULL
if ("gene_name" %in% colnames(mcols(atac_anno))) {
  anno_names <- unique(as.character(mcols(atac_anno)$gene_name))
} else if ("gene_id" %in% colnames(mcols(atac_anno))) {
  anno_names <- unique(as.character(mcols(atac_anno)$gene_id))
}

if (is.null(anno_names)) {
  stop("ATAC annotation has neither gene_name nor gene_id metadata.")
}

overlap_n <- sum(rownames(obj[["RNA_LINK"]]) %in% anno_names)
message("RNA_LINK features overlapping ATAC annotation gene names: ", overlap_n)

if (overlap_n == 0) {
  stop("No overlap between RNA_LINK rownames and ATAC annotation gene coordinates.")
}

# --------------------------------------------------
# 5. RegionStats
# --------------------------------------------------
DefaultAssay(obj) <- "ATAC"
message("Running RegionStats...")
obj <- RegionStats(
  object = obj,
  genome = BSgenome.FMdeer.Custom.v1
)

# --------------------------------------------------
# 6. LinkPeaks using RNA_LINK
# --------------------------------------------------
message("Running LinkPeaks with RNA_LINK...")
obj <- LinkPeaks(
  object = obj,
  peak.assay = "ATAC",
  expression.assay = "RNA_LINK",
  distance = 500000,
  min.cells = 10,
  method = "pearson"
)

saveRDS(obj, file = file.path(outDir, "FMdeer_multiome.wnn.final_annotated.linkpeaks.rds"))
message("Saved object with links.")
