library(tidyverse)
library(stringr)

extract_attr <- function(x, key) {
  pattern <- paste0(key, ' "([^"]*)"')
  m <- stringr::str_match(x, pattern)
  m[,2]
}

# 1. read gtf
gtf <- read.delim(
  "01.reference/complete.genomic.gtf",
  comment.char = "#",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = ""
)

colnames(gtf) <- c(
  "seqname", "source", "feature", "start", "end",
  "score", "strand", "frame", "attribute"
)

gtf_parsed <- gtf %>%
  mutate(
    gene_id = extract_attr(attribute, "gene_id"),
    transcript_id = extract_attr(attribute, "transcript_id"),
    gene_name = extract_attr(attribute, "gene"),
    protein_id = extract_attr(attribute, "protein_id"),
    description = extract_attr(attribute, "description")
  )

gene_info <- gtf_parsed %>%
  filter(feature == "gene") %>%
  select(gene_id, gene_name, description) %>%
  distinct()

protein_map <- gtf_parsed %>%
  filter(!is.na(protein_id), protein_id != "", !is.na(gene_id), gene_id != "") %>%
  select(gene_id, transcript_id, protein_id) %>%
  distinct()

gtf_map <- protein_map %>%
  left_join(gene_info, by = "gene_id") %>%
  distinct()

# 2. read eggnog
egg <- read.delim(
  "01.reference/FMdeer.emapper.annotations",
  comment.char = "#",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = ""
)

colnames(egg) <- c(
  "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
  "max_annot_lvl", "COG_category", "Description", "Preferred_name",
  "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
  "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC",
  "CAZy", "BiGG_Reaction", "PFAMs"
)

# 3. merge
anno_all <- gtf_map %>%
  left_join(egg, by = c("protein_id" = "query"))

gene_annotation <- anno_all %>%
  mutate(
    final_description = coalesce(Description, description),
    final_symbol = coalesce(gene_name, Preferred_name)
  ) %>%
  select(
    gene_id, final_symbol, gene_name, Preferred_name,
    transcript_id, protein_id, final_description,
    GOs, KEGG_ko, KEGG_Pathway, PFAMs, COG_category
  ) %>%
  distinct()

write.table(
  gene_annotation,
  file = "FMdeer_gene_annotation.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 4. GO TERM2GENE
go_term2gene <- anno_all %>%
  select(gene_id, GOs) %>%
  filter(!is.na(GOs), GOs != "-", GOs != "") %>%
  separate_rows(GOs, sep = ",") %>%
  rename(GO_ID = GOs) %>%
  filter(GO_ID != "") %>%
  distinct()

go_term2gene_cp <- go_term2gene %>%
  select(GO_ID, gene_id)

write.table(
  go_term2gene_cp,
  file = "FMdeer_TERM2GENE_GO.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 5. KEGG KO TERM2GENE
keggko_term2gene_cp <- anno_all %>%
  select(gene_id, KEGG_ko) %>%
  filter(!is.na(KEGG_ko), KEGG_ko != "-", KEGG_ko != "") %>%
  separate_rows(KEGG_ko, sep = ",") %>%
  distinct() %>%
  select(KEGG_ko, gene_id)

write.table(
  keggko_term2gene_cp,
  file = "FMdeer_TERM2GENE_KEGGKO.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# 6. KEGG pathway TERM2GENE
pathway_term2gene_cp <- anno_all %>%
  select(gene_id, KEGG_Pathway) %>%
  filter(!is.na(KEGG_Pathway), KEGG_Pathway != "-", KEGG_Pathway != "") %>%
  separate_rows(KEGG_Pathway, sep = ",") %>%
  distinct() %>%
  select(KEGG_Pathway, gene_id)

write.table(
  pathway_term2gene_cp,
  file = "FMdeer_TERM2GENE_KEGGPathway.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)
