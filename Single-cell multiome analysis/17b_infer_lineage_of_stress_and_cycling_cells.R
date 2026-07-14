#!/usr/bin/env Rscript
# 17_infer_lineage_of_stress_and_cycling_cells.R
# Infer the parent lineage of state-dominated clusters 4 (stress) and 31
# (cycling) from nearest stable-lineage neighbours in WNN-UMAP space.
#
# Input : 06.multiome/04.seurat_signac/FMdeer_multiome.wnn.final_annotated.rds
# Output: 06.multiome/04.seurat_signac/state_lineage_inference/
#
# This script does not overwrite the original celltype_major annotation.
# It writes celltype_major_inferred, lineage_inference_confidence, and
# celltype_with_state as additional metadata columns.

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
})

if (!requireNamespace("RANN", quietly = TRUE)) {
  stop("Package 'RANN' is required. Install it once with install.packages('RANN').")
}

# -------------------------------
# User-adjustable settings
# -------------------------------
outDir <- "06.multiome/04.seurat_signac"
inputFile <- file.path(outDir, "FMdeer_multiome.wnn.final_annotated.rds")
inferenceDir <- file.path(outDir, "state_lineage_inference")
state_clusters <- c("4", "31")
reduction_name <- "wnn.umap"
n_neighbors <- 50L
min_confidence <- 0.60

dir.create(inferenceDir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(inputFile)) stop("Annotated object not found: ", inputFile)

obj <- readRDS(inputFile)
required_columns <- c("seurat_clusters", "celltype_major", "cell_state")
if (!all(required_columns %in% colnames(obj@meta.data))) {
  stop("Object is missing: ", paste(setdiff(required_columns, colnames(obj@meta.data)), collapse = ", "),
       ". Run 19_apply_cluster_annotation_from_table.R first.")
}
if (!reduction_name %in% Reductions(obj)) stop("Reduction not found: ", reduction_name)

obj$cluster_id_for_inference <- as.character(obj$seurat_clusters)
state_cells <- colnames(obj)[obj$cluster_id_for_inference %in% state_clusters]
reference_cells <- colnames(obj)[
  !(obj$cluster_id_for_inference %in% state_clusters) &
    !is.na(obj$celltype_major) &
    as.character(obj$celltype_major) != "Unresolved_lineage"
]

if (length(state_cells) == 0) stop("No cells found in state clusters: ", paste(state_clusters, collapse = ", "))
if (length(reference_cells) < 10) stop("Too few stable-lineage reference cells.")

emb <- Embeddings(obj, reduction = reduction_name)
reference_cells <- intersect(reference_cells, rownames(emb))
state_cells <- intersect(state_cells, rownames(emb))
if (length(state_cells) == 0 || length(reference_cells) < 10) {
  stop("State or reference cells are absent from the selected reduction.")
}

k_use <- min(as.integer(n_neighbors), length(reference_cells))
nn <- RANN::nn2(
  data = emb[reference_cells, , drop = FALSE],
  query = emb[state_cells, , drop = FALSE],
  k = k_use
)

reference_labels <- as.character(obj$celltype_major[reference_cells])
cell_votes <- rbindlist(lapply(seq_along(state_cells), function(i) {
  neighbour_labels <- reference_labels[nn$nn.idx[i, ]]
  vote_table <- sort(table(neighbour_labels), decreasing = TRUE)
  data.table(
    cell = state_cells[i],
    original_cluster = as.character(obj$cluster_id_for_inference[state_cells[i]]),
    cell_state = as.character(obj$cell_state[state_cells[i]]),
    inferred_lineage = names(vote_table)[1],
    lineage_inference_confidence = as.numeric(vote_table[1]) / sum(vote_table),
    n_neighbors = k_use,
    mean_neighbor_distance = mean(nn$nn.dists[i, ])
  )
}))
cell_votes[, inference_status := fifelse(
  lineage_inference_confidence >= min_confidence,
  "accepted", "low_confidence"
)]

# Preserve original major labels for stable cells; infer only state-cluster cells.
obj$celltype_major_inferred <- as.character(obj$celltype_major)
obj$lineage_inference_confidence <- NA_real_
obj$lineage_inference_status <- "not_applicable"

obj$celltype_major_inferred[cell_votes$cell] <- fifelse(
  cell_votes$inference_status == "accepted",
  cell_votes$inferred_lineage,
  "Unresolved_lineage"
)
obj$lineage_inference_confidence[cell_votes$cell] <- cell_votes$lineage_inference_confidence
obj$lineage_inference_status[cell_votes$cell] <- cell_votes$inference_status

obj$celltype_with_state <- ifelse(
  as.character(obj$cell_state) == "None_detected",
  obj$celltype_major_inferred,
  paste(obj$celltype_major_inferred, as.character(obj$cell_state), sep = " | ")
)

# Cluster-level consensus is an audit summary only; cell-level confidence is
# retained because state clusters can contain more than one parent lineage.
cluster_summary <- cell_votes[, {
  lineage_votes <- sort(table(inferred_lineage), decreasing = TRUE)
  .(
    n_cells = .N,
    cluster_consensus_lineage = names(lineage_votes)[1],
    cluster_consensus_fraction = as.numeric(lineage_votes[1]) / sum(lineage_votes),
    median_cell_confidence = median(lineage_inference_confidence),
    mean_cell_confidence = mean(lineage_inference_confidence),
    n_accepted = sum(inference_status == "accepted"),
    fraction_accepted = mean(inference_status == "accepted")
  )
}, by = .(original_cluster, cell_state)]

fwrite(cell_votes, file.path(inferenceDir, "state_cell_lineage_neighbour_votes.tsv"), sep = "\t")
fwrite(cluster_summary, file.path(inferenceDir, "state_cluster_lineage_inference_summary.tsv"), sep = "\t")

metadata_export <- as.data.table(obj@meta.data, keep.rownames = "cell")
fwrite(metadata_export, file.path(inferenceDir, "cell_metadata_with_state_lineage_inference.tsv"), sep = "\t")

if (reduction_name %in% Reductions(obj)) {
  p <- DimPlot(obj, reduction = reduction_name, group.by = "celltype_with_state", label = TRUE, repel = TRUE) +
    ggtitle("WNN UMAP: inferred lineage of stress and cycling cells")
  ggsave(file.path(inferenceDir, "UMAP_state_lineage_inference.pdf"), p, width = 14, height = 9)
}

outputFile <- file.path(outDir, "FMdeer_multiome.wnn.final_annotated.state_lineage_inferred.rds")
saveRDS(obj, outputFile)
message("Done. Inferred object saved to: ", outputFile)
message("Accepted assignments require neighbour-vote confidence >= ", min_confidence)
