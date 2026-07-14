suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# =========================================================
# 0. paths
# =========================================================
baseDir <- "06.multiome/04.seurat_signac"
objFile <- file.path(baseDir, "FMdeer_multiome.wnn.final_annotated.linkpeaks.rds")
candidateFile <- file.path(baseDir, "musk_function_genes.strict_confidence.re_scored.tsv")
finalMarkerDir <- file.path(baseDir, "final_celltype_markers")

figDir <- "06.multiome/05.figure5"
dir.create(figDir, recursive = TRUE, showWarnings = FALSE)

fgDir <- file.path(figDir, "Figure5FG_per_gene")
dir.create(fgDir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(objFile)) stop("Missing object: ", objFile)
if (!file.exists(candidateFile)) stop("Missing candidate file: ", candidateFile)

obj <- readRDS(objFile)
cand <- fread(candidateFile)

if (!all(c("RNA", "ATAC") %in% Assays(obj))) {
  stop("Object must contain RNA and ATAC assays.")
}
if (!"celltype" %in% colnames(obj@meta.data)) {
  stop("Object metadata must contain celltype.")
}
if (!"sample" %in% colnames(obj@meta.data)) {
  stop("Object metadata must contain sample.")
}

Idents(obj) <- "celltype"

# =========================================================
# helpers
# =========================================================
keep_existing_genes <- function(obj, genes, assay = "RNA") {
  if (!assay %in% Assays(obj)) return(character(0))
  genes[genes %in% rownames(obj[[assay]])]
}

get_fc_col <- function(dt) {
  for (cc in c("avg_log2FC", "avg_logFC", "log2FC", "avg_diff")) {
    if (cc %in% colnames(dt)) return(cc)
  }
  stop("No logFC column found.")
}

match_gene_to_annotation <- function(g, anno_df) {
  g2 <- gsub("-", "_", g)

  if ("gene_name" %in% colnames(anno_df)) {
    hit <- anno_df[gene_name == g]
    if (nrow(hit) > 0) return(hit[1])
  }

  if ("gene_id" %in% colnames(anno_df)) {
    hit <- anno_df[gene_id == g]
    if (nrow(hit) > 0) return(hit[1])

    hit <- anno_df[gene_id == g2]
    if (nrow(hit) > 0) return(hit[1])
  }

  if ("gene_name" %in% colnames(anno_df)) {
    hit <- anno_df[gene_name == g2]
    if (nrow(hit) > 0) return(hit[1])
  }

  NULL
}

file_manifest <- list()

add_manifest <- function(panel, file, source) {
  file_manifest[[length(file_manifest) + 1]] <<- data.table(
    panel = panel,
    file = file,
    source = source
  )
}

# =========================================================
# A. UMAPs
# =========================================================
pA1 <- DimPlot(
  obj,
  reduction = "wnn.umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + ggtitle("Cluster")

pA2 <- DimPlot(
  obj,
  reduction = "wnn.umap",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE
) + ggtitle("Cell type")

pA3 <- DimPlot(
  obj,
  reduction = "wnn.umap",
  group.by = "sample"
) + ggtitle("Sample")

pdf(file.path(figDir, "Figure5A_UMAP_cluster_celltype_sample.pdf"), width = 16, height = 12)
print((pA1 | pA2) / pA3)
dev.off()
add_manifest("A", file.path(figDir, "Figure5A_UMAP_cluster_celltype_sample.pdf"), "50_make_figure5_panels.R")

# =========================================================
# B. marker heatmap (RNA final celltype markers)
# =========================================================
rnaMarkerTopFile <- file.path(finalMarkerDir, "RNA_final_celltype_markers_top20.tsv")
if (file.exists(rnaMarkerTopFile)) {
  rna_markers <- fread(rnaMarkerTopFile)
  fc_col <- get_fc_col(rna_markers)

  top_heat <- rna_markers[order(cluster, -get(fc_col))]
  top_heat <- top_heat[, head(.SD, 3), by = cluster]

  heat_genes <- unique(top_heat$gene)
  heat_genes <- keep_existing_genes(obj, heat_genes, assay = "RNA")

  if (length(heat_genes) > 1) {
    DefaultAssay(obj) <- "RNA"
    pdf(file.path(figDir, "Figure5B_marker_heatmap.pdf"), width = 12, height = 10)
    print(
      DoHeatmap(
        obj,
        features = heat_genes,
        group.by = "celltype",
        assay = "RNA"
      ) + ggtitle("Representative marker genes by cell type")
    )
    dev.off()
    add_manifest("B", file.path(figDir, "Figure5B_marker_heatmap.pdf"), "50_make_figure5_panels.R")
  }
}

# =========================================================
# C. representative RNA/ACTIVITY feature plots
# pick biologically informative genes from current results
# =========================================================
rna_genes <- c("KCNC2", "PROX1", "COL1A1", "ACTG2", "IL7R", "VWF", "GABRR2")
rna_genes <- keep_existing_genes(obj, rna_genes, assay = "RNA")

act_genes <- c("KCNC2", "PROX1", "COL1A1", "ACTG2", "IL7R", "VWF", "GABRR2")
act_genes <- keep_existing_genes(obj, act_genes, assay = "ACTIVITY")

pdf(file.path(figDir, "Figure5C_featureplots_RNA_ACTIVITY.pdf"), width = 14, height = 10)

if (length(rna_genes) > 0) {
  DefaultAssay(obj) <- "RNA"
  for (g in rna_genes) {
    print(FeaturePlot(obj, features = g, reduction = "wnn.umap", order = TRUE) + ggtitle(paste("RNA", g)))
  }
}

if ("ACTIVITY" %in% Assays(obj) && length(act_genes) > 0) {
  DefaultAssay(obj) <- "ACTIVITY"
  for (g in act_genes) {
    print(FeaturePlot(obj, features = g, reduction = "wnn.umap", order = TRUE) + ggtitle(paste("ACTIVITY", g)))
  }
}

dev.off()
add_manifest("C", file.path(figDir, "Figure5C_featureplots_RNA_ACTIVITY.pdf"), "50_make_figure5_panels.R")

# =========================================================
# D. sample/celltype composition
# =========================================================
meta <- as.data.table(obj@meta.data)
comp <- meta[, .N, by = .(sample, celltype)]
setnames(comp, "N", "n_cells")
comp[, frac_in_sample := n_cells / sum(n_cells), by = sample]

pD1 <- ggplot(comp, aes(x = sample, y = frac_in_sample, fill = celltype)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Fraction", title = "Cell type composition by sample")

pD2 <- ggplot(comp, aes(x = sample, y = n_cells, fill = celltype)) +
  geom_col(width = 0.8) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample", y = "Cell number", title = "Cell type counts by sample")

pdf(file.path(figDir, "Figure5D_sample_celltype_composition.pdf"), width = 14, height = 10)
print(pD1 / pD2)
dev.off()
add_manifest("D", file.path(figDir, "Figure5D_sample_celltype_composition.pdf"), "50_make_figure5_panels.R")

# =========================================================
# E. chromatin accessibility / motif enrichment
# 当前不在此脚本重算 motif，直接记录现有 ArchR 输出文件
# =========================================================
archr_peak_heatmap <- "06.multiome/03.archr/ArchR_proj/03_markerPeaks_heatmap.pdf"
if (file.exists(archr_peak_heatmap)) {
  add_manifest("E", archr_peak_heatmap, "03_peak_calling_marker_analysis.R (existing result)")
}

# =========================================================
# F/G. peak-to-gene linkage / coverage plots
# choose a compact, representative set
# =========================================================
anno <- Annotation(obj[["ATAC"]])
anno_df <- as.data.table(as.data.frame(anno))

fg_genes <- c("KCNC2", "PROX1", "COL1A1", "ACTG2", "GABRR2")

status_dt <- data.table(
  input_gene = character(),
  matched_gene_name = character(),
  matched_gene_id = character(),
  success = logical(),
  mode = character(),
  output_file = character(),
  note = character()
)

pdf(file.path(figDir, "Figure5F_G_peak2gene_coverage_combined.pdf"), width = 12, height = 8)

for (g in fg_genes) {
  hit <- match_gene_to_annotation(g, anno_df)

  if (is.null(hit) || nrow(hit) == 0) {
    status_dt <- rbind(
      status_dt,
      data.table(
        input_gene = g,
        matched_gene_name = NA_character_,
        matched_gene_id = NA_character_,
        success = FALSE,
        mode = "no_match",
        output_file = NA_character_,
        note = "No annotation match"
      ),
      fill = TRUE
    )
    next
  }

  gene_name_use <- if ("gene_name" %in% colnames(hit)) as.character(hit$gene_name[[1]]) else NA_character_
  gene_id_use <- if ("gene_id" %in% colnames(hit)) as.character(hit$gene_id[[1]]) else NA_character_

  region_gr <- GRanges(
    seqnames = hit$seqnames[[1]],
    ranges = IRanges(start = hit$start[[1]], end = hit$end[[1]])
  )

  safe_name <- if (!is.na(gene_name_use) && gene_name_use != "") gene_name_use else g
  safe_name <- gsub("[^A-Za-z0-9_.-]", "_", safe_name)

  outfile <- file.path(fgDir, paste0("CoveragePlot_", safe_name, ".pdf"))

  ok <- FALSE
  mode <- NA_character_
  note <- ""

  tryCatch({
    p <- CoveragePlot(
      object = obj,
      region = region_gr,
      peaks = TRUE,
      links = TRUE,
      extend.upstream = 50000,
      extend.downstream = 50000,
      group.by = "celltype"
    )
    print(p + ggtitle(paste0("Coverage + links: ", g)))
    pdf(outfile, width = 12, height = 8)
    print(p + ggtitle(paste0("Coverage + links: ", g)))
    dev.off()
    ok <- TRUE
    mode <- "region+links"
    note <- "success"
  }, error = function(e) {
    note <<- e$message
  })

  if (!ok) {
    tryCatch({
      p <- CoveragePlot(
        object = obj,
        region = region_gr,
        peaks = TRUE,
        extend.upstream = 50000,
        extend.downstream = 50000,
        group.by = "celltype"
      )
      print(p + ggtitle(paste0("Coverage only: ", g)))
      pdf(outfile, width = 12, height = 8)
      print(p + ggtitle(paste0("Coverage only: ", g)))
      dev.off()
      ok <- TRUE
      mode <- "region_only"
      note <- "success_after_links_fallback"
    }, error = function(e) {
      note <<- paste(note, " | fallback:", e$message)
    })
  }

  status_dt <- rbind(
    status_dt,
    data.table(
      input_gene = g,
      matched_gene_name = gene_name_use,
      matched_gene_id = gene_id_use,
      success = ok,
      mode = mode,
      output_file = outfile,
      note = note
    ),
    fill = TRUE
  )
}

dev.off()

fwrite(status_dt, file.path(figDir, "Figure5F_G_status.tsv"), sep = "\t")
add_manifest("F/G", file.path(figDir, "Figure5F_G_peak2gene_coverage_combined.pdf"), "50_make_figure5_panels.R")

# =========================================================
# save manifest
# =========================================================
manifest_dt <- rbindlist(file_manifest, fill = TRUE)
fwrite(manifest_dt, file.path(figDir, "Figure5_file_manifest.tsv"), sep = "\t")

message("Done. Figure 5 outputs saved to: ", figDir)
