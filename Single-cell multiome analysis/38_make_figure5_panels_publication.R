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
# 0.1. hybrid vector/raster export settings
# =========================================================
# Large UMAP/FeaturePlot point layers make Illustrator/CorelDRAW/Photoshop very
# slow because every cell is stored as an editable vector point. Use ggrastr to
# rasterize only GeomPoint layers while keeping axes, labels, legends, titles,
# gene models, coverage tracks and other annotations as vector elements.
use_point_raster <- TRUE
point_raster_dpi <- 600
has_ggrastr <- requireNamespace("ggrastr", quietly = TRUE)
if (use_point_raster && !has_ggrastr) {
  message(
    "Package 'ggrastr' is not installed; UMAP/FeaturePlot points will remain vector. ",
    "Install with install.packages('ggrastr') if hybrid vector/raster PDFs are needed."
  )
}

# =========================================================
# 0. paths
# =========================================================
baseDir <- "06.multiome/04.seurat_signac"
objFile <- file.path(baseDir, "FMdeer_multiome.wnn.final_annotated.linkpeaks.rds")
finalMarkerDir <- file.path(baseDir, "final_celltype_markers")
geneAnnotFile <- "/data/xb/FMdeer/06.multiome/04.seurat_signac/FMdeer_unified_gene_annotation.unique.tsv"

pubDir <- "06.multiome/05.figure5_pub"
dir.create(pubDir, recursive = TRUE, showWarnings = FALSE)

fgDir <- file.path(pubDir, "Figure5FG_loci")
dir.create(fgDir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(objFile)) stop("Missing object: ", objFile)

obj <- readRDS(objFile)

if (!all(c("RNA", "ATAC") %in% Assays(obj))) {
  stop("Object must contain RNA and ATAC assays.")
}
if (!"celltype" %in% colnames(obj@meta.data)) stop("celltype missing")
if (!"sample" %in% colnames(obj@meta.data)) stop("sample missing")

obj_all <- obj
exclude_celltypes <- c("Unresolved_lineage")
keep_cells <- rownames(obj@meta.data)[
  !as.character(obj@meta.data$celltype) %in% exclude_celltypes
]
if (length(keep_cells) < ncol(obj)) {
  message(
    "Excluding celltype(s) from Figure 5 outputs: ",
    paste(exclude_celltypes, collapse = ", "),
    " (removed ", ncol(obj) - length(keep_cells), " cells)"
  )
  obj <- subset(obj, cells = keep_cells)
  obj$celltype <- droplevels(factor(as.character(obj$celltype)))
}

Idents(obj) <- "celltype"

# =========================================================
# 1. plotting theme
# =========================================================
theme_pub <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold")
    )
}

get_pastel_palette <- function(groups) {
  groups <- unique(as.character(groups))
  base_cols <- c(
    "#F4A6A6", "#F6C28B", "#F7E7A6", "#B8DEB2", "#9FD8CB",
    "#A7C7E7", "#B8B5E8", "#D8B4E2", "#F3B6D0", "#C9D6A3",
    "#BFD7EA", "#E6C9A8", "#C8E6C9", "#D7CCC8", "#F8BBD0",
    "#B2DFDB", "#D1C4E9", "#FFE0B2", "#C5CAE9", "#DCEDC8"
  )
  if (length(groups) > length(base_cols)) {
    base_cols <- grDevices::colorRampPalette(base_cols)(length(groups))
  } else {
    base_cols <- base_cols[seq_along(groups)]
  }
  setNames(base_cols, groups)
}

rasterize_point_layers <- function(p, dpi = point_raster_dpi) {
  if (!use_point_raster || !has_ggrastr) return(p)
  ggrastr::rasterise(p, layers = "Point", dpi = dpi)
}

save_pdf_png <- function(plot_obj, file_base, width = 8, height = 6, dpi = 600) {
  ggsave(
    filename = paste0(file_base, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height
  )
  ggsave(
    filename = paste0(file_base, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi
  )
}

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

load_gene_short_name_map <- function(annotation_file) {
  if (!file.exists(annotation_file)) {
    warning("Gene annotation file not found: ", annotation_file, ". Heatmap will use original gene IDs.")
    return(data.table(gene = character(), gene_label = character()))
  }

  ann <- fread(annotation_file, sep = "\t")
  id_col <- intersect(
    c("gene_id", "gene", "GeneID", "geneid", "id", "feature", "feature_id"),
    colnames(ann)
  )[1]
  label_col <- intersect(
    c("gene_name", "gene_short_name", "symbol", "GeneSymbol", "external_gene_name", "name"),
    colnames(ann)
  )[1]

  if (is.na(id_col) || is.na(label_col)) {
    warning(
      "Could not infer gene ID/name columns from annotation file. Columns are: ",
      paste(colnames(ann), collapse = ", "),
      ". Heatmap will use original gene IDs."
    )
    return(data.table(gene = character(), gene_label = character()))
  }

  out <- ann[
    !is.na(get(id_col)) & nzchar(as.character(get(id_col))),
    .(
      gene = as.character(get(id_col)),
      gene_label = as.character(get(label_col))
    )
  ]
  out <- out[is.na(gene_label) | !nzchar(gene_label), gene_label := gene]

  # Support both FMdeer_XXXXX and FMdeer-XXXXX naming conventions.
  out2 <- copy(out)
  out2[, gene := gsub("_", "-", gene)]
  unique(rbind(out, out2, fill = TRUE), by = "gene")
}

make_publication_marker_heatmap <- function(obj, genes, group_col = "celltype",
                                            assay = "RNA", gene_label_map = NULL) {
  DefaultAssay(obj) <- assay

  avg <- AverageExpression(
    obj,
    assays = assay,
    features = genes,
    group.by = group_col,
    slot = "data",
    verbose = FALSE
  )[[assay]]

  avg <- avg[intersect(genes, rownames(avg)), , drop = FALSE]
  if (nrow(avg) < 2 || ncol(avg) < 2) {
    stop("Too few genes or groups for Figure5B heatmap.")
  }

  z <- t(scale(t(as.matrix(avg))))
  z[is.na(z)] <- 0
  z <- pmax(pmin(z, 2), -2)

  heat_dt <- as.data.table(as.table(z))
  setnames(heat_dt, c("gene", "celltype", "z_score"))
  heat_dt[, gene := as.character(gene)]
  heat_dt[, celltype := as.character(celltype)]

  label_dt <- data.table(gene = rownames(z))
  if (!is.null(gene_label_map) && nrow(gene_label_map) > 0) {
    label_dt <- merge(label_dt, gene_label_map, by = "gene", all.x = TRUE)
  } else {
    label_dt[, gene_label := gene]
  }
  label_dt[is.na(gene_label) | !nzchar(gene_label), gene_label := gene]
  label_dt[, gene_label_unique := make.unique(gene_label)]

  heat_dt <- merge(heat_dt, label_dt[, .(gene, gene_label_unique)], by = "gene", all.x = TRUE)
  heat_dt[, gene_label_unique := factor(gene_label_unique, levels = rev(label_dt$gene_label_unique))]
  heat_dt[, celltype := factor(celltype, levels = colnames(z))]

  ggplot(heat_dt, aes(x = celltype, y = gene_label_unique, fill = z_score)) +
    geom_tile(color = "white", linewidth = 0.08) +
    scale_fill_gradient2(
      low = "#C9DFF2",
      mid = "#FAFAFA",
      high = "#D99A9A",
      midpoint = 0,
      limits = c(-2, 2),
      name = "Scaled\nexpression"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = "Top representative marker per cell type"
    ) +
    theme_pub() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
      axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 7, color = "black"),
      axis.text.y = element_text(size = 6.5, color = "black"),
      axis.ticks = element_blank(),
      legend.position = "right",
      legend.key.height = unit(0.45, "cm"),
      legend.key.width = unit(0.22, "cm"),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.margin = margin(8, 12, 18, 8)
    )
}

compress_peak_coordinate_labels <- function(labels, keep_each_side = 3L) {
  if (is.null(labels)) return(labels)
  labels <- as.character(labels)
  n <- length(labels)
  if (n <= keep_each_side * 2 + 1) return(labels)

  # Only compress coordinate-like labels; leave simple cluster labels alone.
  coord_like <- grepl("^(chr|Chr|CHR)[A-Za-z0-9._-]+[:-][0-9,]+[-:][0-9,]+", labels)
  if (mean(coord_like, na.rm = TRUE) < 0.25) return(labels)

  out <- rep("", n)
  left_idx <- seq_len(min(keep_each_side, n))
  right_idx <- seq(max(1, n - keep_each_side + 1), n)
  mid_idx <- ceiling(n / 2)
  out[left_idx] <- labels[left_idx]
  out[right_idx] <- labels[right_idx]
  out[mid_idx] <- "..."
  out
}

clean_archr_marker_heatmap <- function(ht) {
  clean_one_heatmap <- function(h) {
    if (is.null(h) || !inherits(h, "Heatmap")) return(h)

    # Compress long peak-coordinate labels across the top. Show only a few
    # coordinates at the left/right edges and one ellipsis in the middle.
    h <- tryCatch({
      if ("column_names_param" %in% slotNames(h) && !is.null(h@column_names_param)) {
        h@column_names_param$show <- TRUE
        if (!is.null(h@column_names_param$labels)) {
          h@column_names_param$labels <- compress_peak_coordinate_labels(
            h@column_names_param$labels,
            keep_each_side = 3L
          )
        }
        h@column_names_param$rot <- 60
        if (is.null(h@column_names_param$gp)) {
          h@column_names_param$gp <- grid::gpar(fontsize = 5)
        } else {
          h@column_names_param$gp$fontsize <- 5
        }
      }
      h
    }, error = function(e) {
      message("Skipping peak-coordinate label compression for one heatmap component: ", e$message)
      h
    })

    # ArchR cluster labels are often C1-C25. For Figure 5E, use simple 1-25
    # labels on the right side to save space and improve readability.
    h <- tryCatch({
      if ("row_names_param" %in% slotNames(h) &&
          !is.null(h@row_names_param) &&
          !is.null(h@row_names_param$labels)) {
        h@row_names_param$labels <- sub("^C(?=\\d+$)", "", h@row_names_param$labels, perl = TRUE)
      }
      h
    }, error = function(e) {
      message("Skipping row-label simplification for one heatmap component: ", e$message)
      h
    })
    h
  }

  if (is.null(ht)) return(ht)
  if (inherits(ht, "HeatmapList")) {
    ht <- tryCatch({
      for (i in seq_along(ht@ht_list)) {
        if (!is.null(ht@ht_list[[i]]) && inherits(ht@ht_list[[i]], "Heatmap")) {
          ht@ht_list[[i]] <- clean_one_heatmap(ht@ht_list[[i]])
        }
      }
      ht
    }, error = function(e) {
      message("Skipping HeatmapList-level cleanup: ", e$message)
      ht
    })
  } else if (inherits(ht, "Heatmap")) {
    ht <- clean_one_heatmap(ht)
  }
  ht
}

manifest_list <- list()
add_manifest <- function(panel, file, note) {
  manifest_list[[length(manifest_list) + 1]] <<- data.table(
    panel = panel,
    file = file,
    note = note
  )
}

# =========================================================
# 2. A panel: UMAP
# =========================================================
pA_cluster <- DimPlot(
  obj,
  reduction = "wnn.umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + ggtitle("Cluster identity") + theme_pub()
pA_cluster <- rasterize_point_layers(pA_cluster)

pA_celltype <- DimPlot(
  obj,
  reduction = "wnn.umap",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE
) + ggtitle("Annotated cell populations") + theme_pub()
pA_celltype <- rasterize_point_layers(pA_celltype)

pA_sample <- DimPlot(
  obj,
  reduction = "wnn.umap",
  group.by = "sample"
) + ggtitle("Sample origin") + theme_pub()
pA_sample <- rasterize_point_layers(pA_sample)

pdf(file.path(pubDir, "Figure5A_UMAP.pdf"), width = 16, height = 12)
print(
  (pA_cluster | pA_celltype) /
    (pA_sample | plot_spacer())
)
dev.off()
add_manifest("A", file.path(pubDir, "Figure5A_UMAP.pdf"), "Integrated UMAP by cluster/celltype/sample")

# =========================================================
# 3. B panel: marker heatmap
# =========================================================
rnaMarkerTopFile <- file.path(finalMarkerDir, "RNA_final_celltype_markers_top20.tsv")
if (file.exists(rnaMarkerTopFile)) {
  rna_markers <- fread(rnaMarkerTopFile)
  if ("cluster" %in% colnames(rna_markers)) {
    rna_markers <- rna_markers[
      !as.character(cluster) %in% exclude_celltypes
    ]
  }
  fc_col <- get_fc_col(rna_markers)

  top_heat <- rna_markers[order(cluster, -get(fc_col))]
  top_heat <- top_heat[, head(.SD, 1), by = cluster]
  heat_genes <- unique(top_heat$gene)
  heat_genes <- keep_existing_genes(obj, heat_genes, assay = "RNA")

  if (length(heat_genes) > 1) {
    gene_label_map <- load_gene_short_name_map(geneAnnotFile)
    pB <- make_publication_marker_heatmap(
      obj = obj,
      genes = heat_genes,
      group_col = "celltype",
      assay = "RNA",
      gene_label_map = gene_label_map
    )

    pdf(file.path(pubDir, "Figure5B_marker_heatmap.pdf"), width = 4.5, height = 6.5)
    print(pB)
    dev.off()
    add_manifest("B", file.path(pubDir, "Figure5B_marker_heatmap.pdf"), "Top one RNA marker per cell type heatmap")
  }
}

# =========================================================
# 4. C panel: representative feature plots
# =========================================================
rna_genes <- c("KCNC2", "PROX1", "COL1A1", "ACTG2", "IL7R", "VWF")
rna_genes <- keep_existing_genes(obj, rna_genes, assay = "RNA")

act_genes <- c("KCNC2", "PROX1", "COL1A1", "ACTG2", "IL7R", "VWF")
act_genes <- keep_existing_genes(obj, act_genes, assay = "ACTIVITY")

feature_plot_list <- list()

if (length(rna_genes) > 0) {
  DefaultAssay(obj) <- "RNA"
  for (g in rna_genes) {
    p_feature <- FeaturePlot(obj, features = g, reduction = "wnn.umap", order = TRUE) +
      ggtitle(paste("RNA", g)) +
      theme_pub()
    feature_plot_list[[length(feature_plot_list) + 1]] <- rasterize_point_layers(p_feature)
  }
}

if ("ACTIVITY" %in% Assays(obj) && length(act_genes) > 0) {
  DefaultAssay(obj) <- "ACTIVITY"
  for (g in act_genes) {
    p_feature <- FeaturePlot(obj, features = g, reduction = "wnn.umap", order = TRUE) +
      ggtitle(paste("Gene activity", g)) +
      theme_pub()
    feature_plot_list[[length(feature_plot_list) + 1]] <- rasterize_point_layers(p_feature)
  }
}

if (length(feature_plot_list) > 0) {
  pdf(file.path(pubDir, "Figure5C_featureplots.pdf"), width = 16, height = 12, onefile = TRUE)
  page_ids <- split(seq_along(feature_plot_list), ceiling(seq_along(feature_plot_list) / 4))
  for (idx in page_ids) {
    print(
      wrap_plots(feature_plot_list[idx], ncol = 2) +
        plot_annotation(title = "Representative RNA expression and gene activity")
    )
  }
  dev.off()
  add_manifest(
    "C",
    file.path(pubDir, "Figure5C_featureplots.pdf"),
    "Representative RNA and gene activity plots; four panels per page"
  )
}

# =========================================================
# 5. D panel: composition
# =========================================================
meta <- as.data.table(obj_all@meta.data)
meta <- meta[!is.na(celltype)]
comp <- meta[, .N, by = .(sample, celltype)]
setnames(comp, "N", "n_cells")
comp[, frac_in_sample := n_cells / sum(n_cells), by = sample]
comp[, celltype := factor(as.character(celltype), levels = unique(as.character(celltype)))]
composition_palette <- get_pastel_palette(levels(comp$celltype))
if ("Unresolved_lineage" %in% names(composition_palette)) {
  composition_palette["Unresolved_lineage"] <- "#D9D9D9"
}

pD1 <- ggplot(comp, aes(x = sample, y = frac_in_sample, fill = celltype)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = composition_palette, name = "Cell type") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Sample", y = "Fraction", title = "Cell population composition by sample") +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold")
  )

pD2 <- ggplot(comp, aes(x = sample, y = n_cells, fill = celltype)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = composition_palette, name = "Cell type") +
  labs(x = "Sample", y = "Cell number", title = "Cell counts by sample") +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold")
  )

pdf(file.path(pubDir, "Figure5D_composition.pdf"), width = 14, height = 10)
print(pD1 / pD2)
dev.off()
add_manifest("D", file.path(pubDir, "Figure5D_composition.pdf"), "Sample contributions across annotated populations")

# =========================================================
# 6. E panel: chromatin accessibility marker heatmap
# =========================================================
archr_marker_peaks_file <- "06.multiome/03.archr/ArchR_proj/markersPeaks.rds"
figure5e_file <- file.path(pubDir, "Figure5E_chromatin_accessibility_heatmap.pdf")
figure5e_peak_table <- file.path(pubDir, "Figure5E_chromatin_accessibility_heatmap_peak_labels.tsv")
if (file.exists(archr_marker_peaks_file)) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    warning("SummarizedExperiment is not installed; cannot regenerate Figure5E heatmap.")
  } else if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    warning("ComplexHeatmap is not installed; cannot regenerate Figure5E heatmap.")
  } else if (!requireNamespace("circlize", quietly = TRUE)) {
    warning("circlize is not installed; cannot regenerate Figure5E heatmap.")
  } else {
    markersPeaks <- readRDS(archr_marker_peaks_file)

    pick_marker_assay <- function(se, candidates) {
      available <- SummarizedExperiment::assayNames(se)
      hit <- candidates[candidates %in% available]
      if (length(hit) > 0) return(hit[[1]])
      available[[1]]
    }

    log2fc_assay <- pick_marker_assay(markersPeaks, c("Log2FC", "log2FC", "Log2FCMatrix"))
    fdr_assay <- pick_marker_assay(markersPeaks, c("FDR", "fdr", "padj", "qval"))
    log2fc <- as.matrix(SummarizedExperiment::assay(markersPeaks, log2fc_assay))
    fdr <- as.matrix(SummarizedExperiment::assay(markersPeaks, fdr_assay))

    cluster_names <- colnames(log2fc)
    if (is.null(cluster_names)) cluster_names <- paste0("C", seq_len(ncol(log2fc)))
    cluster_names_clean <- sub("^C(?=\\d+$)", "", cluster_names, perl = TRUE)

    selected_by_cluster <- lapply(seq_along(cluster_names), function(j) {
      ok <- which(is.finite(fdr[, j]) & is.finite(log2fc[, j]) & fdr[, j] <= 0.01 & log2fc[, j] >= 1)
      if (length(ok) == 0) return(integer(0))
      ok <- ok[order(log2fc[ok, j], decreasing = TRUE)]
      head(ok, 4L)
    })
    names(selected_by_cluster) <- cluster_names

    selected <- unique(unlist(selected_by_cluster, use.names = FALSE))
    if (length(selected) == 0) {
      warning("No marker peaks passed FDR <= 0.01 and Log2FC >= 1 for Figure5E.")
    } else {
      if (length(selected) > 140L) selected <- selected[seq_len(140L)]

      owner <- rep(NA_character_, length(selected))
      names(owner) <- selected
      for (cl in names(selected_by_cluster)) {
        idx <- intersect(selected_by_cluster[[cl]], selected)
        owner[as.character(idx)] <- cl
      }
      owner[is.na(owner)] <- "Marker peak"

      mat <- t(log2fc[selected, , drop = FALSE])
      rownames(mat) <- cluster_names_clean
      plot_peak_ids <- paste0("P", seq_along(selected))
      colnames(mat) <- plot_peak_ids

      z <- t(scale(t(mat)))
      z[!is.finite(z)] <- 0
      z <- pmax(pmin(z, 2), -2)

      rr <- tryCatch(SummarizedExperiment::rowRanges(markersPeaks), error = function(e) NULL)
      peak_label <- paste0("PeakID_", selected)
      full_coordinate <- rep(NA_character_, length(selected))
      label_source <- rep("peak_id_fallback", length(selected))
      if (!is.null(rr) && length(rr) == nrow(log2fc)) {
        chr_all <- as.character(GenomicRanges::seqnames(rr))
        start_all <- IRanges::start(rr)
        end_all <- IRanges::end(rr)
        midpoint_all <- (as.numeric(start_all) + as.numeric(end_all)) / 2
        valid_coord <- grepl("^(chr|Chr|CHR)[A-Za-z0-9._-]+$", chr_all) &
          is.finite(midpoint_all) &
          as.numeric(end_all) >= as.numeric(start_all)
        if (any(valid_coord)) {
          peak_label <- ifelse(
            midpoint_all[selected] >= 1e6,
            sprintf("%s:%.2f Mb", chr_all[selected], midpoint_all[selected] / 1e6),
            sprintf("%s:%.1f kb", chr_all[selected], midpoint_all[selected] / 1e3)
          )
          full_coordinate <- paste0(chr_all[selected], ":", start_all[selected], "-", end_all[selected])
          label_source <- rep("genomic_coordinate", length(selected))
        }
      }

      fwrite(
        data.table(
          column_index = seq_along(selected),
          plot_label = plot_peak_ids,
          marker_cluster = owner[as.character(selected)],
          peak_label = peak_label,
          full_coordinate = full_coordinate,
          label_source = label_source,
          peak_source_index = selected
        ),
        figure5e_peak_table,
        sep = "\t"
      )

      group_palette <- setNames(
        grDevices::colorRampPalette(c(
          "#F7B7A3", "#F3DFA2", "#B8E0D2", "#A8DADC", "#BFD7EA",
          "#CDB4DB", "#FFC8DD", "#D9D9D9"
        ))(length(unique(owner))),
        unique(owner)
      )

      top_anno <- ComplexHeatmap::HeatmapAnnotation(
        cluster = owner[as.character(selected)],
        col = list(cluster = group_palette),
        annotation_name_gp = grid::gpar(fontsize = 7),
        annotation_legend_param = list(
          title = "Marker cluster",
          title_gp = grid::gpar(fontsize = 7, fontface = "bold"),
          labels_gp = grid::gpar(fontsize = 6),
          nrow = 2
        ),
        simple_anno_size = grid::unit(2.2, "mm")
      )

      heatmapPeaks <- ComplexHeatmap::Heatmap(
        z,
        name = "Column Z-score",
        col = circlize::colorRamp2(c(-2, 0, 2), c("#C9DFF2", "#FAFAFA", "#D99A9A")),
        top_annotation = top_anno,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        row_names_side = "right",
        row_names_gp = grid::gpar(fontsize = 7),
        column_names_side = "top",
        column_labels = plot_peak_ids,
        column_names_rot = 45,
        column_names_gp = grid::gpar(fontsize = 5.2),
        column_title = "Marker peaks (P labels; see peak label table)",
        column_title_gp = grid::gpar(fontsize = 7, fontface = "bold"),
        border = TRUE,
        rect_gp = grid::gpar(col = NA),
        heatmap_legend_param = list(
          title_gp = grid::gpar(fontsize = 7, fontface = "bold"),
          labels_gp = grid::gpar(fontsize = 6),
          legend_height = grid::unit(22, "mm")
        )
      )

      pdf(figure5e_file, width = 7.2, height = 4.8, useDingbats = FALSE)
      ComplexHeatmap::draw(
        heatmapPeaks,
        heatmap_legend_side = "bottom",
        annotation_legend_side = "bottom",
        padding = grid::unit(c(2, 2, 2, 2), "mm")
      )
      dev.off()
      add_manifest(
        "E",
        figure5e_file,
        "Custom chromatin accessibility marker peak heatmap; generated without ArchR::plotMarkerHeatmap to avoid column annotation layout errors"
      )
    }
  }
} else {
  warning("ArchR marker peak object not found: ", archr_marker_peaks_file)
}

# =========================================================
# 7. F/G panel: representative loci
# =========================================================
anno <- Annotation(obj[["ATAC"]])
anno_df <- as.data.table(as.data.frame(anno))

# Representative loci for Figure 5F/G.
# KCNC2 is the main example used in the assembled Figure 5. SMPDL3A and NRCAM
# are added here as musk secretion-related loci for the full/supplementary
# peak-to-gene and coverage locus PDFs.
musky_secretion_loci <- c("SMPDL3A", "NRCAM")
panelF_genes <- c("KCNC2", "PROX1", "COL1A1", musky_secretion_loci)
panelG_genes <- c("KCNC2", "GABRR2", "ACTG2", musky_secretion_loci)

plot_locus_set <- function(genes, file_prefix, panel_label) {
  status <- data.table(
    input_gene = character(),
    matched_gene_name = character(),
    matched_gene_id = character(),
    success = logical(),
    mode = character(),
    output_file = character(),
    note = character()
  )

  combined_pdf <- file.path(pubDir, paste0(file_prefix, ".pdf"))
  pdf(combined_pdf, width = 12, height = 8)

  for (g in genes) {
    hit <- match_gene_to_annotation(g, anno_df)

    if (is.null(hit) || nrow(hit) == 0) {
      status <- rbind(
        status,
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
    outfile <- file.path(fgDir, paste0(file_prefix, "_", safe_name, ".pdf"))

    ok <- FALSE
    mode <- NA_character_
    note <- ""

    # 关键：CoveragePlot 前强制切回 ChromatinAssay
    DefaultAssay(obj) <- "ATAC"

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
      print(p + ggtitle(paste0(panel_label, ": ", g)))
      pdf(outfile, width = 12, height = 8)
      print(p + ggtitle(paste0(panel_label, ": ", g)))
      dev.off()
      ok <- TRUE
      mode <- "region+links"
      note <- "success"
    }, error = function(e) {
      note <<- e$message
    })

    if (!ok) {
      DefaultAssay(obj) <- "ATAC"
      tryCatch({
        p <- CoveragePlot(
          object = obj,
          region = region_gr,
          peaks = TRUE,
          extend.upstream = 50000,
          extend.downstream = 50000,
          group.by = "celltype"
        )
        print(p + ggtitle(paste0(panel_label, ": ", g, " (coverage only)")))
        pdf(outfile, width = 12, height = 8)
        print(p + ggtitle(paste0(panel_label, ": ", g, " (coverage only)")))
        dev.off()
        ok <- TRUE
        mode <- "region_only"
        note <- "success_after_links_fallback"
      }, error = function(e) {
        note <<- paste(note, " | fallback:", e$message)
      })
    }

    status <- rbind(
      status,
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
  fwrite(status, file.path(pubDir, paste0(file_prefix, "_status.tsv")), sep = "\t")

  add_manifest(panel_label, combined_pdf, paste0(panel_label, " combined coverage/linkage plots"))
}

plot_locus_set(panelF_genes, "Figure5F_peak2gene_loci", "F")
plot_locus_set(panelG_genes, "Figure5G_coverage_loci", "G")

# =========================================================
# 8. save manifest
# =========================================================
manifest_dt <- rbindlist(manifest_list, fill = TRUE)
fwrite(manifest_dt, file.path(pubDir, "Figure5_pub_manifest.tsv"), sep = "\t")

message("Done. Publication-oriented Figure 5 outputs saved to: ", pubDir)
