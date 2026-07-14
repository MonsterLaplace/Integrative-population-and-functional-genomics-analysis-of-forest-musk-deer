#!/usr/bin/env Rscript
# Create sample and pseudobulk-group metadata for the three analysis levels.

suppressPackageStartupMessages({
  library(data.table)
})

outDir <- "06.multiome/04.seurat_signac"
pbDir <- file.path(outDir, "pseudobulk")
primaryFile <- file.path(pbDir, "RNA_pseudobulk_group_metadata.filtered.tsv")
if (!file.exists(primaryFile)) {
  stop("Missing primary pseudobulk metadata: ", primaryFile, ". Run script 24 first.")
}

make_sample_metadata <- function(samples) {
  x <- data.table(sample = sort(unique(samples)))
  x[, `:=`(
    age_label = NA_character_, age_month = NA_real_, age_group = NA_character_,
    tissue = "musk_gland", species = "Moschus berezovskii", notes = NA_character_
  )]
  x[sample == "FMD3", `:=`(age_label = "11m", age_month = 11, age_group = "young")]
  x[sample == "FMD4", `:=`(age_label = "3y", age_month = 36, age_group = "adult")]
  x[sample == "J21010063", `:=`(age_label = "6m", age_month = 6, age_group = "young")]
  x[sample == "J21090141", `:=`(age_label = "9y", age_month = 108, age_group = "adult")]
  x[, group := age_group] # compatible with downstream edgeR/DESeq2 scripts
  x
}

metadata_files <- data.table(
  file = c(
    "RNA_pseudobulk_group_metadata.filtered.tsv",
    "RNA_pseudobulk_epithelial_minor_group_metadata.filtered.tsv",
    "RNA_pseudobulk_cell_state_group_metadata.filtered.tsv"
  ),
  annotation_level = c("celltype_major_inferred", "celltype_minor", "cell_state"),
  output = c(
    "pseudobulk_metadata_major_celltype.tsv",
    "pseudobulk_metadata_epithelial_minor.tsv",
    "pseudobulk_metadata_cell_state.tsv"
  )
)

available <- metadata_files[file.exists(file.path(pbDir, file))]
if (nrow(available) == 0) stop("No filtered pseudobulk metadata files found.")
all_groups <- rbindlist(lapply(available$file, function(f) fread(file.path(pbDir, f), sep = "\t")), fill = TRUE)
if (!"sample" %in% names(all_groups)) stop("Pseudobulk metadata missing 'sample' column.")

sample_meta <- make_sample_metadata(all_groups$sample)
fwrite(sample_meta, file.path(pbDir, "sample_metadata_for_pseudobulk.tsv"), sep = "\t", quote = FALSE)

for (i in seq_len(nrow(available))) {
  spec <- available[i]
  pb_meta <- fread(file.path(pbDir, spec$file), sep = "\t")
  if (!all(c("pb_group", "sample", "annotation", "n_cells") %in% names(pb_meta))) {
    stop("Unexpected pseudobulk metadata columns in: ", spec$file)
  }
  pb_meta[, annotation_level := spec$annotation_level]
  pb_meta <- merge(pb_meta, sample_meta, by = "sample", all.x = TRUE, sort = FALSE)
  setcolorder(pb_meta, c("pb_group", "sample", "annotation", "annotation_level", "n_cells",
                          "age_label", "age_month", "age_group", "group", "tissue", "species", "notes"))
  fwrite(pb_meta, file.path(pbDir, spec$output), sep = "\t", quote = FALSE)
  message("Saved pseudobulk design metadata: ", file.path(pbDir, spec$output))
}

message("Saved shared sample metadata: ", file.path(pbDir, "sample_metadata_for_pseudobulk.tsv"))
print(sample_meta)
