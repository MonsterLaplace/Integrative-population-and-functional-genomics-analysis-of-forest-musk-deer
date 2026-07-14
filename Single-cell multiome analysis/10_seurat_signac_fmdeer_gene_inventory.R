suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

outDir <- "06.multiome/04.seurat_signac"
obj <- readRDS(file.path(outDir, "FMdeer_multiome.wnn.robust.rds"))

rna_genes <- data.table(gene = rownames(obj[["RNA"]]), assay = "RNA")
activity_genes <- data.table(gene = rownames(obj[["ACTIVITY"]]), assay = "ACTIVITY")

fwrite(rna_genes, file.path(outDir, "RNA_gene_names.tsv"), sep = "\t")
fwrite(activity_genes, file.path(outDir, "ACTIVITY_gene_names.tsv"), sep = "\t")

common_genes <- intersect(rna_genes$gene, activity_genes$gene)
fwrite(data.table(gene = common_genes), file.path(outDir, "common_gene_names.tsv"), sep = "\t")

message("Saved gene inventories.")
