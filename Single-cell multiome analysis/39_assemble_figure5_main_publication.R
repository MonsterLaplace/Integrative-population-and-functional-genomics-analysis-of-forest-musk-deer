suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(grid)
})

# =========================================================
# Assemble Figure 5 main layout
# =========================================================
# Main figure design:
#   - A block: cluster/celltype/sample UMAPs
#   - C block: KCNC2 feature plot placed at A block bottom-right
#   - B block: marker heatmap placed to the right of A/C
#   - D block: composition below A/C
#   - E block: chromatin accessibility heatmap to the right of D
#   - F/G block: one KCNC2 locus plot only, because F and G KCNC2 are redundant
#
# Other generated panels remain available as supplementary files.

baseDir <- "06.multiome/04.seurat_signac"
pubDir <- "06.multiome/05.figure5_pub"
objFile <- file.path(baseDir, "FMdeer_multiome.wnn.final_annotated.linkpeaks.rds")

figureA_file <- file.path(pubDir, "Figure5A_UMAP.pdf")
figureB_file <- file.path(pubDir, "Figure5B_marker_heatmap.pdf")
figureD_file <- file.path(pubDir, "Figure5D_composition.pdf")
figureE_file <- file.path(pubDir, "Figure5E_chromatin_accessibility_heatmap.pdf")
figureF_kcnc2_file <- file.path(pubDir, "Figure5FG_loci", "Figure5F_peak2gene_loci_KCNC2.pdf")
figureG_kcnc2_file <- file.path(pubDir, "Figure5FG_loci", "Figure5G_coverage_loci_KCNC2.pdf")

out_pdf <- file.path(pubDir, "Figure5_main_assembled.pdf")
out_png <- file.path(pubDir, "Figure5_main_assembled.png")
out_manifest <- file.path(pubDir, "Figure5_main_assembled_manifest.tsv")

dir.create(pubDir, recursive = TRUE, showWarnings = FALSE)

required_pkgs <- c("magick", "cowplot")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing package(s): ", paste(missing_pkgs, collapse = ", "),
    ". Install on the server first, e.g. install.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
    "))."
  )
}

use_point_raster <- TRUE
point_raster_dpi <- 600
has_ggrastr <- requireNamespace("ggrastr", quietly = TRUE)
if (use_point_raster && !has_ggrastr) {
  message("Package 'ggrastr' not installed; UMAP/FeaturePlot points remain vector.")
}

pdf_render_density <- 400
main_width <- 16
main_height <- 17
main_png_dpi <- 500
trim_imported_pdf_panels <- TRUE
trim_fuzz_percent <- 6
trim_border_px <- 10

theme_pub <- function() {
  theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black", size = 7),
      axis.title = element_text(color = "black", size = 8),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
      legend.title = element_text(face = "bold", size = 7),
      legend.text = element_text(size = 6),
      plot.margin = margin(2, 2, 2, 2)
    )
}

rasterize_point_layers <- function(p, dpi = point_raster_dpi) {
  if (!use_point_raster || !has_ggrastr) return(p)
  ggrastr::rasterise(p, layers = "Point", dpi = dpi)
}

render_pdf_page_to_png <- function(file, page = 1, density = pdf_render_density) {
  if (!file.exists(file)) stop("Missing panel file: ", file)
  pdftoppm <- Sys.which("pdftoppm")
  if (!nzchar(pdftoppm)) {
    stop(
      "Cannot render PDF because system command 'pdftoppm' is not available. ",
      "Install poppler-utils on Ubuntu, e.g. sudo apt-get install -y poppler-utils."
    )
  }

  tmp_prefix <- tempfile(pattern = "figure5_pdf_page_")
  cmd_status <- system2(
    pdftoppm,
    args = c(
      "-f", as.character(page),
      "-singlefile",
      "-png",
      "-r", as.character(density),
      shQuote(normalizePath(file, mustWork = TRUE)),
      shQuote(tmp_prefix)
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  png_file <- paste0(tmp_prefix, ".png")
  if (!file.exists(png_file)) {
    stop(
      "pdftoppm failed to render PDF page from: ", file,
      "\nOutput:\n", paste(cmd_status, collapse = "\n")
    )
  }
  png_file
}

read_pdf_as_plot <- function(file, page = 1, density = pdf_render_density,
                             trim = trim_imported_pdf_panels,
                             fuzz = trim_fuzz_percent,
                             border_px = trim_border_px) {
  png_file <- render_pdf_page_to_png(file, page = page, density = density)
  img <- magick::image_read(png_file)
  if (trim) {
    img <- magick::image_trim(img, fuzz = fuzz)
    if (!is.null(border_px) && border_px > 0) {
      img <- magick::image_border(img, color = "white", geometry = paste0(border_px, "x", border_px))
    }
  }
  grob <- rasterGrob(as.raster(img), interpolate = TRUE)
  ggdraw() + draw_grob(grob, x = 0, y = 0, width = 1, height = 1)
}

add_panel_label <- function(p, label, x = 0.006, y = 0.994, size = 17) {
  ggdraw(p) +
    draw_label(
      label,
      x = x,
      y = y,
      hjust = 0,
      vjust = 1,
      fontface = "bold",
      size = size
    )
}

make_A_C_block <- function(obj) {
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

  kcnc2_assay <- if ("ACTIVITY" %in% Assays(obj) && "KCNC2" %in% rownames(obj[["ACTIVITY"]])) {
    "ACTIVITY"
  } else if ("RNA" %in% Assays(obj) && "KCNC2" %in% rownames(obj[["RNA"]])) {
    "RNA"
  } else {
    NA_character_
  }

  if (is.na(kcnc2_assay)) {
    pC <- ggdraw() +
      draw_label("KCNC2 not found in RNA or ACTIVITY assay", x = 0.5, y = 0.5, size = 10)
  } else {
    DefaultAssay(obj) <- kcnc2_assay
    pC <- FeaturePlot(
      obj,
      features = "KCNC2",
      reduction = "wnn.umap",
      order = TRUE
    ) +
      ggtitle(ifelse(kcnc2_assay == "ACTIVITY", "Gene activity KCNC2", "RNA KCNC2")) +
      theme_pub()
    pC <- rasterize_point_layers(pC)
  }
  pC <- add_panel_label(pC, "C", size = 15)

  pA <- (pA_cluster | pA_celltype) / (pA_sample | pC)
  add_panel_label(pA, "A", size = 17)
}

if (!file.exists(objFile)) stop("Missing object: ", objFile)
obj <- readRDS(objFile)
if (!"celltype" %in% colnames(obj@meta.data)) stop("celltype missing in object metadata.")
if (!"sample" %in% colnames(obj@meta.data)) stop("sample missing in object metadata.")

exclude_celltypes <- c("Unresolved_lineage")
keep_cells <- rownames(obj@meta.data)[
  !as.character(obj@meta.data$celltype) %in% exclude_celltypes
]
if (length(keep_cells) < ncol(obj)) {
  message(
    "Excluding celltype(s) from assembled Figure 5: ",
    paste(exclude_celltypes, collapse = ", "),
    " (removed ", ncol(obj) - length(keep_cells), " cells)"
  )
  obj <- subset(obj, cells = keep_cells)
  obj$celltype <- droplevels(factor(as.character(obj$celltype)))
}

message("Building A/C block from Seurat object...")
panel_AC <- make_A_C_block(obj)

message("Reading existing panel PDFs...")
panel_B <- add_panel_label(read_pdf_as_plot(figureB_file), "B", size = 17)
panel_D <- add_panel_label(read_pdf_as_plot(figureD_file), "D", size = 17)
panel_E <- add_panel_label(read_pdf_as_plot(figureE_file), "E", size = 17)

fg_file <- if (file.exists(figureF_kcnc2_file)) figureF_kcnc2_file else figureG_kcnc2_file
if (!file.exists(fg_file)) {
  stop("Missing KCNC2 locus file. Checked: ", figureF_kcnc2_file, " and ", figureG_kcnc2_file)
}
panel_FG <- add_panel_label(read_pdf_as_plot(fg_file), "F/G", size = 17)

# CNS-style fixed canvas layout.
# Absolute coordinates make panel sizing more predictable than nested plot_grid,
# especially when imported PDF panels contain different internal aspect ratios.
main_fig <- ggdraw() +
  theme(plot.background = element_rect(fill = "white", colour = NA)) +
  draw_plot(panel_AC, x = 0.010, y = 0.610, width = 0.470, height = 0.380) +
  draw_plot(panel_B,  x = 0.500, y = 0.600, width = 0.485, height = 0.390) +
  draw_plot(panel_D,  x = 0.010, y = 0.365, width = 0.470, height = 0.225) +
  draw_plot(panel_E,  x = 0.500, y = 0.360, width = 0.485, height = 0.235) +
  draw_plot(panel_FG, x = 0.010, y = 0.020, width = 0.975, height = 0.325)

message("Saving assembled Figure 5...")
ggsave(out_pdf, main_fig, width = main_width, height = main_height, units = "in", bg = "white")
ggsave(out_png, main_fig, width = main_width, height = main_height, units = "in", dpi = main_png_dpi, bg = "white")

all_panel_files <- list.files(pubDir, pattern = "^Figure5.*\\.(pdf|tsv)$", recursive = TRUE, full.names = TRUE)
main_files <- normalizePath(
  c(figureA_file, figureB_file, figureD_file, figureE_file, fg_file, out_pdf, out_png),
  winslash = "/",
  mustWork = FALSE
)
manifest <- data.table(
  file = normalizePath(all_panel_files, winslash = "/", mustWork = FALSE)
)
manifest[, role := fifelse(file %in% main_files, "main_figure_input_or_output", "supplementary")]
manifest[, note := fifelse(
  grepl("Figure5C_featureplots", file),
  "Only KCNC2 is used in the main figure; full feature plot PDF can be supplementary.",
  fifelse(
    grepl("Figure5F|Figure5G", file) & !file %in% main_files,
    "Non-KCNC2 locus panels can be supplementary.",
    ""
  )
)]
fwrite(manifest, out_manifest, sep = "\t")

message("Done. Main Figure 5 saved to: ", out_pdf)
