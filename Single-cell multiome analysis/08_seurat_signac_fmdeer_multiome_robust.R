suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(S4Vectors)
  library(rtracklayer)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(BSgenome)
  library(BSgenome.FMdeer.Custom.v1)
})

options(future.globals.maxSize = 16 * 1024^3)

# =========================================================
# 0. Paths and parameters
# =========================================================
outDir <- "06.multiome/04.seurat_signac"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

gtfFile <- "01.reference/FMdeer.primary.gtf"
bsgenomeObj <- BSgenome.FMdeer.Custom.v1

sample_dirs <- c(
  FMD3      = "06.multiome/cellranger_arc/FMD3/outs",
  FMD4      = "06.multiome/cellranger_arc/FMD4/outs",
  J21010063 = "06.multiome/cellranger_arc/J21010063/outs",
  J21090141 = "06.multiome/cellranger_arc/J21090141/outs"
)

valid_chr <- c(paste0("Chr", sprintf("%02d", 1:29)), "ChrX", "ChrY")

# QC thresholds
qc_rna_count_min    <- 500
qc_rna_feature_min  <- 200
qc_atac_count_min   <- 1000
qc_tss_min          <- 2
qc_nucleosome_max   <- 4

# clustering
rna_pcs     <- 1:30
atac_lsis   <- 2:30
resolution  <- 0.5

# =========================================================
# 1. Helpers
# =========================================================
check_sample_files <- function(sample_name, outdir) {
  h5file <- file.path(outdir, "filtered_feature_bc_matrix.h5")
  fragfile <- file.path(outdir, "atac_fragments.tsv.gz")
  fragidx <- file.path(outdir, "atac_fragments.tsv.gz.tbi")

  miss <- c(h5file, fragfile, fragidx)[!file.exists(c(h5file, fragfile, fragidx))]
  if (length(miss) > 0) {
    stop(
      "Missing files for sample ", sample_name, ":\n",
      paste0("  - ", miss, collapse = "\n")
    )
  }
}

guess_rna_assay_name <- function(x) {
  nms <- names(x)
  candidates <- c("Gene Expression", "GeneExpression", "RNA", "gex", "expression")
  hit <- candidates[candidates %in% nms]
  if (length(hit) > 0) return(hit[1])

  # fallback: choose first matrix not obviously peaks/atac/antibody
  for (nm in nms) {
    low <- tolower(nm)
    if (!grepl("peak|atac|antibody|adt|crispr", low)) {
      return(nm)
    }
  }
  stop("Cannot infer RNA assay from h5. Available assays: ", paste(nms, collapse = ", "))
}

guess_atac_assay_name <- function(x) {
  nms <- names(x)
  candidates <- c("Peaks", "ATAC", "Peak")
  hit <- candidates[candidates %in% nms]
  if (length(hit) > 0) return(hit[1])

  # fallback: choose assay with peak/atac keyword
  for (nm in nms) {
    low <- tolower(nm)
    if (grepl("peak|atac", low)) {
      return(nm)
    }
  }
  stop("Cannot infer ATAC assay from h5. Available assays: ", paste(nms, collapse = ", "))
}

safe_make_annotation <- function(gtfFile, valid_chr) {
  message("Importing GTF: ", gtfFile)
  gtf <- import(gtfFile)

  if (length(gtf) == 0) {
    stop("Imported GTF is empty.")
  }

  chr <- as.character(seqnames(gtf))
  keep <- !is.na(chr) & chr %in% valid_chr
  gtf <- gtf[keep]

  if (length(gtf) == 0) {
    stop("No GTF records remain after filtering to valid chromosomes.")
  }

  seqlevels(gtf) <- intersect(seqlevels(gtf), valid_chr)

  mc <- mcols(gtf)
  nms <- colnames(mc)

  # -----------------------------
  # ensure type
  # -----------------------------
  if (!"type" %in% nms) {
    if ("feature" %in% nms) {
      mc$type <- as.character(mc$feature)
    } else {
      mc$type <- "unknown"
    }
  }

  # -----------------------------
  # ensure gene_id
  # -----------------------------
  if (!"gene_id" %in% nms) {
    if ("gene" %in% nms) {
      mc$gene_id <- as.character(mc$gene)
    } else if ("ID" %in% nms) {
      mc$gene_id <- as.character(mc$ID)
    } else if ("Parent" %in% nms) {
      mc$gene_id <- as.character(mc$Parent)
    } else if ("Name" %in% nms) {
      mc$gene_id <- as.character(mc$Name)
    } else {
      mc$gene_id <- paste0("gene_", seq_along(gtf))
    }
  }

  # -----------------------------
  # ensure gene_name
  # -----------------------------
  if (!"gene_name" %in% nms) {
    if ("gene" %in% nms) {
      mc$gene_name <- as.character(mc$gene)
    } else if ("Name" %in% nms) {
      mc$gene_name <- as.character(mc$Name)
    } else if ("gene_id" %in% colnames(mc)) {
      mc$gene_name <- as.character(mc$gene_id)
    } else if ("ID" %in% nms) {
      mc$gene_name <- as.character(mc$ID)
    } else {
      mc$gene_name <- paste0("gene_", seq_along(gtf))
    }
  }

  # -----------------------------
  # ensure gene_biotype
  # -----------------------------
  if (!"gene_biotype" %in% nms) {
    if ("gene_type" %in% nms) {
      mc$gene_biotype <- as.character(mc$gene_type)
    } else if ("gene_biotype" %in% nms) {
      mc$gene_biotype <- as.character(mc$gene_biotype)
    } else if ("biotype" %in% nms) {
      mc$gene_biotype <- as.character(mc$biotype)
    } else {
      mc$gene_biotype <- "unknown"
    }
  }

  # force character
  mc$type <- as.character(mc$type)
  mc$gene_id <- as.character(mc$gene_id)
  mc$gene_name <- as.character(mc$gene_name)
  mc$gene_biotype <- as.character(mc$gene_biotype)

  mcols(gtf) <- mc

  message("Annotation seqnames:")
  print(sort(unique(as.character(seqnames(gtf)))))

  message("Annotation metadata columns:")
  print(colnames(mcols(gtf)))

  # 检查必需列
  required_cols <- c("gene_name", "gene_id", "gene_biotype", "type")
  miss <- required_cols[!required_cols %in% colnames(mcols(gtf))]
  if (length(miss) > 0) {
    stop("Annotation still missing required columns: ", paste(miss, collapse = ", "))
  }

  gtf
}


CreateArcSeuratObject <- function(sample_name, outdir, annotation_gr) {
  message("=================================================")
  message("Processing sample: ", sample_name)

  h5file <- file.path(outdir, "filtered_feature_bc_matrix.h5")
  fragpath <- file.path(outdir, "atac_fragments.tsv.gz")

  dat <- Read10X_h5(h5file)
  message("Assays found in h5 for ", sample_name, ": ", paste(names(dat), collapse = ", "))

  rna_name <- guess_rna_assay_name(dat)
  atac_name <- guess_atac_assay_name(dat)

  message("Using RNA assay:  ", rna_name)
  message("Using ATAC assay: ", atac_name)

  rna_counts <- dat[[rna_name]]
  atac_counts <- dat[[atac_name]]

  obj <- CreateSeuratObject(
    counts = rna_counts,
    assay = "RNA",
    project = sample_name
  )

  obj[["ATAC"]] <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    fragments = fragpath,
    annotation = annotation_gr,
    min.cells = 1,
    min.features = 1
  )

  obj$sample <- sample_name
  obj
}

write_cluster_summary <- function(obj, outFile) {
  dt <- as.data.table(table(obj$seurat_clusters, obj$sample))
  colnames(dt) <- c("cluster", "sample", "n_cells")
  fwrite(dt, outFile, sep = "\t")
}

write_meta_summary <- function(obj, outFile) {
  meta <- obj@meta.data
  meta$cell <- rownames(meta)
  fwrite(as.data.table(meta), outFile, sep = "\t")
}

# =========================================================
# 2. Check files
# =========================================================
invisible(mapply(check_sample_files, names(sample_dirs), sample_dirs))

# =========================================================
# 3. Build annotation
# =========================================================
annotation <- safe_make_annotation(gtfFile, valid_chr)

# =========================================================
# 4. Create sample objects
# =========================================================
objs <- mapply(
  FUN = CreateArcSeuratObject,
  sample_name = names(sample_dirs),
  outdir = sample_dirs,
  MoreArgs = list(annotation_gr = annotation),
  SIMPLIFY = FALSE
)

# =========================================================
# 5. Merge
# =========================================================
message("Merging all samples...")
obj <- merge(
  x = objs[[1]],
  y = objs[-1],
  add.cell.ids = names(sample_dirs),
  project = "FMdeer_multiome"
)

rm(objs)
gc()

saveRDS(obj, file = file.path(outDir, "FMdeer_multiome.merged.raw.rds"))

# =========================================================
# 6. Initial RNA preprocessing
# =========================================================
message("Initial RNA preprocessing...")
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 50)

# =========================================================
# 7. ATAC QC + reduction
# =========================================================
message("ATAC reduction...")
DefaultAssay(obj) <- "ATAC"
message("Skipping Signac NucleosomeSignal/TSSEnrichment on custom genome annotation.")
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = "q0")
obj <- RunSVD(obj)


pdf(file.path(outDir, "01_QC_violin_before_filter.pdf"), width = 16, height = 8)
print(VlnPlot(
  obj,
  features = c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  pt.size = 0,
  ncol = 5
))
dev.off()

pdf(file.path(outDir, "02_DepthCor_before_filter.pdf"), width = 7, height = 6)
print(DepthCor(obj))
dev.off()

# =========================================================
# 8. Filter cells
# =========================================================
message("Filtering cells by QC...")

obj <- subset(
  obj,
  subset =
    nCount_RNA > qc_rna_count_min &
    nFeature_RNA > qc_rna_feature_min &
    nCount_ATAC > qc_atac_count_min
)


message("Cells remaining after QC: ", ncol(obj))

pdf(file.path(outDir, "03_QC_violin_after_filter.pdf"), width = 16, height = 8)
print(VlnPlot(
  obj,
  features = c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  pt.size = 0,
  ncol = 5
))
dev.off()

# =========================================================
# 9. Recompute reductions after QC
# =========================================================
message("Recomputing reductions after filtering...")

DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 50)

DefaultAssay(obj) <- "ATAC"
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = "q0")
obj <- RunSVD(obj)

pdf(file.path(outDir, "04_DepthCor_after_filter.pdf"), width = 7, height = 6)
print(DepthCor(obj))
dev.off()

# =========================================================
# 10. WNN
# =========================================================
message("Running WNN multimodal integration...")

obj <- FindMultiModalNeighbors(
  object = obj,
  reduction.list = list("pca", "lsi"),
  dims.list = list(rna_pcs, atac_lsis),
  modality.weight.name = "RNA.weight"
)

obj <- RunUMAP(
  object = obj,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_"
)

obj <- FindClusters(
  object = obj,
  graph.name = "wsnn",
  algorithm = 3,
  resolution = resolution
)

# Optional single-modality UMAPs
obj <- RunUMAP(
  object = obj,
  reduction = "pca",
  dims = rna_pcs,
  reduction.name = "rna.umap",
  reduction.key = "rnaUMAP_"
)

obj <- RunUMAP(
  object = obj,
  reduction = "lsi",
  dims = atac_lsis,
  reduction.name = "atac.umap",
  reduction.key = "atacUMAP_"
)

pdf(file.path(outDir, "05_UMAPs.pdf"), width = 14, height = 10)
print(DimPlot(obj, reduction = "wnn.umap", group.by = "sample"))
print(DimPlot(obj, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE))
print(DimPlot(obj, reduction = "rna.umap", group.by = "seurat_clusters", label = TRUE))
print(DimPlot(obj, reduction = "atac.umap", group.by = "seurat_clusters", label = TRUE))
dev.off()

# =========================================================
# 11. Gene activity
# =========================================================
message("Computing gene activity...")
DefaultAssay(obj) <- "ATAC"
gene.activities <- GeneActivity(obj)

obj[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

DefaultAssay(obj) <- "ACTIVITY"
obj <- NormalizeData(obj)
obj <- ScaleData(obj)

# =========================================================
# 12. RegionStats
# =========================================================
message("Computing RegionStats...")
DefaultAssay(obj) <- "ATAC"
obj <- RegionStats(obj, genome = bsgenomeObj)

# =========================================================
# 13. Export summaries
# =========================================================
write_cluster_summary(obj, file.path(outDir, "cluster_by_sample_counts.tsv"))
write_meta_summary(obj, file.path(outDir, "cell_metadata.tsv"))

cluster_sizes <- as.data.table(table(obj$seurat_clusters))
colnames(cluster_sizes) <- c("cluster", "n_cells")
fwrite(cluster_sizes, file.path(outDir, "cluster_sizes.tsv"), sep = "\t")

# =========================================================
# 14. Save final object
# =========================================================
saveRDS(obj, file = file.path(outDir, "FMdeer_multiome.wnn.robust.rds"))

message("Done.")
message("Saved object: ", file.path(outDir, "FMdeer_multiome.wnn.robust.rds"))
