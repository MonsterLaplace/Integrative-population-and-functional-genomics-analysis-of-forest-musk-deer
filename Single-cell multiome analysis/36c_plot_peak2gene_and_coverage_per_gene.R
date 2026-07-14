suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
})

outDir <- "06.multiome/04.seurat_signac"
objFile <- file.path(outDir, "FMdeer_multiome.wnn.final_annotated.linkpeaks.rds")

# 使用 36b 成功状态表
statusFile <- file.path(
  outDir,
  "peak2gene_coverage_plots_fallback",
  "36b_peak2gene_coverage_candidate_genes.status.tsv"
)

plotDir <- file.path(outDir, "peak2gene_coverage_plots_per_gene")
dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(objFile)) stop("Object not found: ", objFile)
if (!file.exists(statusFile)) stop("Status file not found: ", statusFile)

obj <- readRDS(objFile)
status_dt <- fread(statusFile)

anno <- Annotation(obj[["ATAC"]])
if (is.null(anno) || length(anno) == 0) {
  stop("ATAC annotation missing.")
}

anno_df <- as.data.table(as.data.frame(anno))

# --------------------------------------------------
# helper: match gene to annotation
# --------------------------------------------------
match_gene_to_annotation <- function(g, anno_df) {
  g2 <- gsub("-", "_", g)

  # 1. exact gene_name
  if ("gene_name" %in% colnames(anno_df)) {
    hit <- anno_df[gene_name == g]
    if (nrow(hit) > 0) return(hit[1])
  }

  # 2. gene_id
  if ("gene_id" %in% colnames(anno_df)) {
    hit <- anno_df[gene_id == g]
    if (nrow(hit) > 0) return(hit[1])

    hit <- anno_df[gene_id == g2]
    if (nrow(hit) > 0) return(hit[1])
  }

  # 3. gene_name underscore fallback
  if ("gene_name" %in% colnames(anno_df)) {
    hit <- anno_df[gene_name == g2]
    if (nrow(hit) > 0) return(hit[1])
  }

  NULL
}

# --------------------------------------------------
# only keep successful genes from 36b
# --------------------------------------------------
ok_genes <- status_dt[success == TRUE, input_gene]
ok_genes <- unique(ok_genes)

message("Genes to export individually:")
print(ok_genes)

out_status <- data.table(
  input_gene = character(),
  matched_gene_name = character(),
  matched_gene_id = character(),
  output_file = character(),
  success = logical(),
  mode = character(),
  note = character()
)

for (g in ok_genes) {
  message("Exporting gene: ", g)

  hit <- match_gene_to_annotation(g, anno_df)

  if (is.null(hit) || nrow(hit) == 0) {
    out_status <- rbind(
      out_status,
      data.table(
        input_gene = g,
        matched_gene_name = NA_character_,
        matched_gene_id = NA_character_,
        output_file = NA_character_,
        success = FALSE,
        mode = "no_match",
        note = "No matching gene found in ATAC annotation"
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

  outfile <- file.path(plotDir, paste0("CoveragePlot_", safe_name, ".pdf"))

  ok <- FALSE
  mode <- NA_character_
  note <- ""

  # 1) preferred: region + links
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

    pdf(outfile, width = 12, height = 8)
    print(p + ggtitle(paste0("Coverage + links: ", g, " (", gene_name_use, ")")))
    dev.off()

    ok <- TRUE
    mode <- "region+links"
    note <- "success"
  }, error = function(e) {
    note <<- e$message
  })

  # 2) fallback: region only
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

      pdf(outfile, width = 12, height = 8)
      print(p + ggtitle(paste0("Coverage only: ", g, " (", gene_name_use, ")")))
      dev.off()

      ok <- TRUE
      mode <- "region_only"
      note <- "success_after_links_fallback"
    }, error = function(e) {
      note <<- paste(note, " | fallback:", e$message)
    })
  }

  out_status <- rbind(
    out_status,
    data.table(
      input_gene = g,
      matched_gene_name = gene_name_use,
      matched_gene_id = gene_id_use,
      output_file = outfile,
      success = ok,
      mode = mode,
      note = note
    ),
    fill = TRUE
  )
}

fwrite(
  out_status,
  file.path(plotDir, "36c_peak2gene_coverage_per_gene.status.tsv"),
  sep = "\t"
)

message("Done. Per-gene PDFs saved to: ", plotDir)
