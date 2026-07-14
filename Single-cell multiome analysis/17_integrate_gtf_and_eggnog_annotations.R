suppressPackageStartupMessages({
  library(rtracklayer)
  library(data.table)
})

# =========================================================
# files
# =========================================================
gtfFile <- "01.reference/FMdeer.primary.gtf"
eggnogFile <- "01.reference/FMdeer.emapper.annotations"

outDir <- "06.multiome/04.seurat_signac"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

markerDir <- file.path(outDir, "markers")

# =========================================================
# helper: ensure a column exists
# =========================================================
ensure_col <- function(dt, col, default = NA_character_) {
  if (!col %in% colnames(dt)) {
    dt[, (col) := default]
  }
  dt
}

# =========================================================
# helper: safely read eggNOG emapper annotations
# keeps #query header, removes ## comment lines
# =========================================================
read_eggnog_annotations <- function(eggnogFile) {
  x <- readLines(eggnogFile)

  # remove comment lines beginning with ##
  x <- x[!grepl("^##", x)]

  if (length(x) == 0) {
    stop("eggNOG annotation file is empty after removing comment lines.")
  }

  # convert #query -> query
  x[1] <- sub("^#query", "query", x[1])

  tf <- tempfile(fileext = ".tsv")
  writeLines(x, tf)

  dt <- fread(
    tf,
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE
  )

  if (!"query" %in% colnames(dt)) {
    stop("eggNOG file missing 'query' column after parsing.")
  }

  dt
}

# =========================================================
# helper: annotate one marker table using unique gene annotation
# matching priority:
#   1) gene_id_dash
#   2) gene_id
#   3) gene_name
# =========================================================
annotate_marker_file <- function(infile, outfile, anno_unique_dt) {
  if (!file.exists(infile)) {
    message("Skip missing marker file: ", infile)
    return(invisible(NULL))
  }

  dt <- fread(infile)
  if (!"gene" %in% colnames(dt)) {
    message("No 'gene' column in marker file: ", infile)
    fwrite(dt, outfile, sep = "\t")
    return(invisible(NULL))
  }

  dt[, row_id__tmp := .I]

  # 1) match by gene_id_dash
  m1 <- merge(
    dt,
    anno_unique_dt,
    by.x = "gene",
    by.y = "gene_id_dash",
    all.x = TRUE,
    sort = FALSE
  )

  # recover original order as much as possible
  if ("row_id__tmp" %in% colnames(m1)) {
    setorder(m1, row_id__tmp)
  }

  # 2) unmatched -> gene_id
  idx1 <- which(is.na(m1$display_name) | m1$display_name == "")
  if (length(idx1) > 0) {
    sub1 <- dt[idx1]

    m2 <- merge(
      sub1,
      anno_unique_dt,
      by.x = "gene",
      by.y = "gene_id",
      all.x = TRUE,
      sort = FALSE
    )

    if ("row_id__tmp" %in% colnames(m2)) {
      setorder(m2, row_id__tmp)
    }

    fill_cols <- intersect(colnames(m2), colnames(m1))
    fill_cols <- setdiff(fill_cols, c("gene", "row_id__tmp"))

    for (cc in fill_cols) {
      m1[[cc]][idx1] <- m2[[cc]]
    }
  }

  # 3) still unmatched -> gene_name
  idx2 <- which(is.na(m1$display_name) | m1$display_name == "")
  if (length(idx2) > 0) {
    sub2 <- dt[idx2]

    m3 <- merge(
      sub2,
      anno_unique_dt,
      by.x = "gene",
      by.y = "gene_name",
      all.x = TRUE,
      sort = FALSE
    )

    if ("row_id__tmp" %in% colnames(m3)) {
      setorder(m3, row_id__tmp)
    }

    fill_cols <- intersect(colnames(m3), colnames(m1))
    fill_cols <- setdiff(fill_cols, c("gene", "row_id__tmp"))

    for (cc in fill_cols) {
      m1[[cc]][idx2] <- m3[[cc]]
    }
  }

  # display label
  if (!"display_label" %in% colnames(m1)) {
    m1[, display_label := fifelse(
      !is.na(display_name) & display_name != "",
      display_name,
      gene
    )]
  }

  # remove temp
  if ("row_id__tmp" %in% colnames(m1)) {
    m1[, row_id__tmp := NULL]
  }

  fwrite(m1, outfile, sep = "\t")
  message("Saved annotated marker file: ", outfile)

  invisible(NULL)
}

# =========================================================
# 1. read GTF
# =========================================================
message("Reading GTF: ", gtfFile)
gtf <- import(gtfFile)

if (length(gtf) == 0) {
  stop("Imported GTF is empty.")
}

gtf_dt <- as.data.table(as.data.frame(gtf))

# unify feature/type column
if (!"type" %in% colnames(gtf_dt)) {
  if ("feature" %in% colnames(gtf_dt)) {
    setnames(gtf_dt, "feature", "type")
  } else {
    gtf_dt[, type := NA_character_]
  }
}

# ensure common columns
gtf_dt <- ensure_col(gtf_dt, "gene_id")
gtf_dt <- ensure_col(gtf_dt, "gene")
gtf_dt <- ensure_col(gtf_dt, "gene_name")
gtf_dt <- ensure_col(gtf_dt, "gene_biotype")
gtf_dt <- ensure_col(gtf_dt, "description")
gtf_dt <- ensure_col(gtf_dt, "product")
gtf_dt <- ensure_col(gtf_dt, "protein_id")
gtf_dt <- ensure_col(gtf_dt, "transcript_id")
gtf_dt <- ensure_col(gtf_dt, "locus_tag")
gtf_dt <- ensure_col(gtf_dt, "type")

# fill gene_name from gene if needed
gtf_dt[is.na(gene_name) | gene_name == "", gene_name := as.character(gene)]
gtf_dt[is.na(gene_name) | gene_name == "", gene_name := as.character(gene_id)]

# fill gene_biotype from possible alternatives
if (all(is.na(gtf_dt$gene_biotype) | gtf_dt$gene_biotype == "")) {
  if ("gene_type" %in% colnames(gtf_dt)) {
    gtf_dt[is.na(gene_biotype) | gene_biotype == "", gene_biotype := as.character(gene_type)]
  } else if ("biotype" %in% colnames(gtf_dt)) {
    gtf_dt[is.na(gene_biotype) | gene_biotype == "", gene_biotype := as.character(biotype)]
  }
}
gtf_dt[is.na(gene_biotype) | gene_biotype == "", gene_biotype := "unknown"]

# normalize ids
gtf_dt[, gene_id := as.character(gene_id)]
gtf_dt[, gene_id_dash := gsub("_", "-", gene_id)]
gtf_dt[, gene_name := as.character(gene_name)]
gtf_dt[, description := as.character(description)]
gtf_dt[, product := as.character(product)]
gtf_dt[, protein_id := as.character(protein_id)]
gtf_dt[, transcript_id := as.character(transcript_id)]
gtf_dt[, locus_tag := as.character(locus_tag)]
gtf_dt[, type := as.character(type)]

message("GTF columns:")
print(colnames(gtf_dt))

# =========================================================
# 2. gene-level table from GTF
# =========================================================
gene_dt0 <- gtf_dt[type == "gene"]

if (nrow(gene_dt0) == 0) {
  gene_dt <- unique(gtf_dt[, .(
    gene_id, gene_id_dash, gene_name, gene_biotype, description, product, locus_tag
  )])
} else {
  gene_dt <- unique(gene_dt0[, .(
    gene_id, gene_id_dash, gene_name, gene_biotype, description, product, locus_tag
  )])
}

# remove rows with missing gene_id
gene_dt <- gene_dt[!is.na(gene_id) & gene_id != ""]

# one row per gene
gene_dt <- gene_dt[order(gene_id), .SD[1], by = gene_id]

# =========================================================
# 3. protein-to-gene map from GTF
# =========================================================
protein_map_dt <- unique(
  gtf_dt[!is.na(protein_id) & protein_id != "", .(
    protein_id,
    gene_id,
    gene_id_dash,
    transcript_id,
    gene_name,
    product,
    description
  )]
)

if (nrow(protein_map_dt) == 0) {
  warning("No protein_id records found in GTF; eggNOG merge may be mostly empty.")
}

# one row per protein
protein_map_dt <- protein_map_dt[order(protein_id), .SD[1], by = protein_id]

# =========================================================
# 4. read eggNOG
# =========================================================
message("Reading eggNOG annotations: ", eggnogFile)
eggnog_dt <- read_eggnog_annotations(eggnogFile)

# rename to avoid ambiguity
if ("Description" %in% colnames(eggnog_dt)) {
  setnames(eggnog_dt, "Description", "eggNOG_Description")
}
if ("Preferred_name" %in% colnames(eggnog_dt)) {
  setnames(eggnog_dt, "Preferred_name", "eggNOG_Preferred_name")
}

message("eggNOG columns:")
print(colnames(eggnog_dt))

# one row per query protein
eggnog_dt <- eggnog_dt[order(query), .SD[1], by = query]

# =========================================================
# 5. merge GTF protein map with eggNOG by protein_id
# =========================================================
anno_dt <- merge(
  protein_map_dt,
  eggnog_dt,
  by.x = "protein_id",
  by.y = "query",
  all.x = TRUE,
  sort = FALSE
)

# merge gene-level metadata
anno_dt <- merge(
  anno_dt,
  gene_dt,
  by = c("gene_id", "gene_id_dash"),
  all.x = TRUE,
  suffixes = c("", ".gene"),
  sort = FALSE
)

# fill from gene-level columns if needed
if ("gene_name.gene" %in% colnames(anno_dt)) {
  anno_dt[is.na(gene_name) | gene_name == "", gene_name := gene_name.gene]
}
if ("description.gene" %in% colnames(anno_dt)) {
  anno_dt[is.na(description) | description == "", description := description.gene]
}
if ("product.gene" %in% colnames(anno_dt)) {
  anno_dt[is.na(product) | product == "", product := product.gene]
}
if ("gene_biotype" %in% colnames(anno_dt) == FALSE && "gene_biotype.gene" %in% colnames(anno_dt)) {
  anno_dt[, gene_biotype := gene_biotype.gene]
}
if ("locus_tag.gene" %in% colnames(anno_dt)) {
  anno_dt[is.na(locus_tag) | locus_tag == "", locus_tag := locus_tag.gene]
}

anno_dt <- ensure_col(anno_dt, "gene_biotype", "unknown")
anno_dt <- ensure_col(anno_dt, "locus_tag")

# display name
if (!"eggNOG_Preferred_name" %in% colnames(anno_dt)) {
  anno_dt[, eggNOG_Preferred_name := NA_character_]
}
if (!"eggNOG_Description" %in% colnames(anno_dt)) {
  anno_dt[, eggNOG_Description := NA_character_]
}

anno_dt[, display_name := fifelse(
  !is.na(gene_name) & gene_name != "",
  gene_name,
  fifelse(
    !is.na(eggNOG_Preferred_name) & eggNOG_Preferred_name != "",
    eggNOG_Preferred_name,
    gene_id
  )
)]

anno_dt[, display_description := fifelse(
  !is.na(description) & description != "",
  description,
  fifelse(
    !is.na(product) & product != "",
    product,
    eggNOG_Description
  )
)]

# =========================================================
# 6. build final unified annotation tables
# =========================================================

# full table: may contain multiple proteins/transcripts per gene
keep_cols_full <- intersect(c(
  "gene_id",
  "gene_id_dash",
  "gene_name",
  "gene_biotype",
  "locus_tag",
  "protein_id",
  "transcript_id",
  "product",
  "description",
  "eggNOG_Preferred_name",
  "eggNOG_Description",
  "seed_ortholog",
  "evalue",
  "score",
  "GOs",
  "EC",
  "KEGG_ko",
  "KEGG_Pathway",
  "KEGG_Module",
  "KEGG_Reaction",
  "BRITE",
  "KEGG_TC",
  "CAZy",
  "BiGG_Reaction",
  "PFAMs",
  "display_name",
  "display_description"
), colnames(anno_dt))

gene_anno_dt <- unique(anno_dt[, ..keep_cols_full])

# unique gene-level table: one row per gene
# this is the table used for marker annotation
gene_anno_unique_dt <- gene_anno_dt[
  order(gene_id, protein_id, transcript_id),
  .SD[1],
  by = .(gene_id, gene_id_dash)
]

# if any gene in gene_dt did not appear in protein/eggnog merge,
# append plain GTF-only rows so every gene can still be annotated
missing_gene_ids <- setdiff(gene_dt$gene_id, gene_anno_unique_dt$gene_id)

if (length(missing_gene_ids) > 0) {
  missing_dt <- gene_dt[gene_id %in% missing_gene_ids]

  missing_dt[, `:=`(
    protein_id = NA_character_,
    transcript_id = NA_character_,
    eggNOG_Preferred_name = NA_character_,
    eggNOG_Description = NA_character_,
    seed_ortholog = NA_character_,
    evalue = NA_real_,
    score = NA_real_,
    GOs = NA_character_,
    EC = NA_character_,
    KEGG_ko = NA_character_,
    KEGG_Pathway = NA_character_,
    KEGG_Module = NA_character_,
    KEGG_Reaction = NA_character_,
    BRITE = NA_character_,
    KEGG_TC = NA_character_,
    CAZy = NA_character_,
    BiGG_Reaction = NA_character_,
    PFAMs = NA_character_,
    display_name = fifelse(!is.na(gene_name) & gene_name != "", gene_name, gene_id),
    display_description = fifelse(
      !is.na(description) & description != "",
      description,
      product
    )
  )]

  need_cols <- union(colnames(gene_anno_unique_dt), colnames(missing_dt))
  for (cc in setdiff(need_cols, colnames(gene_anno_unique_dt))) {
    gene_anno_unique_dt[, (cc) := NA]
  }
  for (cc in setdiff(need_cols, colnames(missing_dt))) {
    missing_dt[, (cc) := NA]
  }

  gene_anno_unique_dt <- rbindlist(
    list(gene_anno_unique_dt[, ..need_cols], missing_dt[, ..need_cols]),
    use.names = TRUE,
    fill = TRUE
  )

  gene_anno_unique_dt <- gene_anno_unique_dt[
    order(gene_id, protein_id, transcript_id),
    .SD[1],
    by = .(gene_id, gene_id_dash)
  ]
}

# =========================================================
# 7. save annotation tables
# =========================================================
fwrite(gene_dt, file.path(outDir, "FMdeer_gtf_gene_table.tsv"), sep = "\t")
fwrite(protein_map_dt, file.path(outDir, "FMdeer_gtf_protein_gene_map.tsv"), sep = "\t")
fwrite(gene_anno_dt, file.path(outDir, "FMdeer_unified_gene_annotation.tsv"), sep = "\t")
fwrite(gene_anno_unique_dt, file.path(outDir, "FMdeer_unified_gene_annotation.unique.tsv"), sep = "\t")

message("Saved unified annotation tables:")
message("  ", file.path(outDir, "FMdeer_unified_gene_annotation.tsv"))
message("  ", file.path(outDir, "FMdeer_unified_gene_annotation.unique.tsv"))

# =========================================================
# 8. annotate marker tables with unique gene-level annotation
# =========================================================
annotate_marker_file(
  infile = file.path(markerDir, "RNA_markers_top20.tsv"),
  outfile = file.path(markerDir, "RNA_markers_top20.unified.tsv"),
  anno_unique_dt = gene_anno_unique_dt
)

annotate_marker_file(
  infile = file.path(markerDir, "RNA_markers_all.tsv"),
  outfile = file.path(markerDir, "RNA_markers_all.unified.tsv"),
  anno_unique_dt = gene_anno_unique_dt
)

annotate_marker_file(
  infile = file.path(markerDir, "ACTIVITY_markers_top20.tsv"),
  outfile = file.path(markerDir, "ACTIVITY_markers_top20.unified.tsv"),
  anno_unique_dt = gene_anno_unique_dt
)

annotate_marker_file(
  infile = file.path(markerDir, "ACTIVITY_markers_all.tsv"),
  outfile = file.path(markerDir, "ACTIVITY_markers_all.unified.tsv"),
  anno_unique_dt = gene_anno_unique_dt
)

message("Done.")
