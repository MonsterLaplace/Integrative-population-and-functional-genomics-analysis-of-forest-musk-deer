library(tidyverse)

# ================================================================
# Extract genes overlapping high-confidence selection windows
# ================================================================
gene_file <- "gene_coordinates.tsv"
autosome_window_file <- "selection_autosome_candidates_highconf.tsv"
sexchr_window_file <- "selection_sexchr_candidates_highconf.tsv"

# Convert common chromosome conventions to the convention used in the window
# tables, e.g. 1/chr1/Chr01 -> Chr01 and X/chrX -> ChrX.
normalize_chr <- function(x) {
  x <- trimws(as.character(x))
  x <- stringr::str_replace(x, stringr::regex("^chromosome", ignore_case = TRUE), "")
  x <- stringr::str_replace(x, stringr::regex("^chr", ignore_case = TRUE), "")
  x <- trimws(x)
  dplyr::case_when(
    toupper(x) == "X" ~ "ChrX",
    toupper(x) == "Y" ~ "ChrY",
    grepl("^[0-9]+$", x) ~ sprintf("Chr%02d", as.integer(x)),
    TRUE ~ paste0("Chr", x)
  )
}

read_windows <- function(path, chromosome_group) {
  if (!file.exists(path)) stop("Cannot find window result: ", path)
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  needed <- c("CHROM", "start", "end")
  if (!all(needed %in% names(x))) {
    stop("Window file must contain these columns: ", paste(needed, collapse = ", "), "; file: ", path)
  }
  x %>%
    mutate(CHROM = normalize_chr(CHROM),
           start = as.numeric(start), end = as.numeric(end),
           chromosome_group = chromosome_group,
           window_id = paste(chromosome_group, CHROM, start, end, sep = "_"))
}

# The supplied coordinate file has: gene, CHROM, start, end.
genes <- read.delim(gene_file, check.names = FALSE, stringsAsFactors = FALSE)
needed_gene <- c("gene", "CHROM", "start", "end")
if (!all(needed_gene %in% names(genes))) {
  stop("Gene coordinate file must contain: ", paste(needed_gene, collapse = ", "))
}
genes <- genes %>%
  transmute(gene, CHROM = normalize_chr(CHROM),
            gene_start = as.numeric(start), gene_end = as.numeric(end)) %>%
  filter(!is.na(gene_start), !is.na(gene_end), gene_start <= gene_end)

autosome_windows <- read_windows(autosome_window_file, "autosome")
sexchr_windows <- read_windows(sexchr_window_file, "sexchr")

# One gene is retained for every high-confidence window that overlaps it.
# Overlap condition: gene_start <= window_end AND gene_end >= window_start.
extract_overlaps <- function(windows, genes) {
  hit_template <- tibble(
    chromosome_group = character(), window_id = character(), window_CHROM = character(),
    window_start = numeric(), window_end = numeric(), gene = character(),
    gene_CHROM = character(), gene_start = numeric(), gene_end = numeric(),
    overlap_start = numeric(), overlap_end = numeric(), overlap_bp = numeric()
  )
  if (nrow(windows) == 0) return(hit_template)
  hits <- purrr::map_dfr(seq_len(nrow(windows)), function(i) {
    w <- windows[i, , drop = FALSE]
    genes %>%
      filter(CHROM == w$CHROM,
             gene_start <= w$end,
             gene_end >= w$start) %>%
      transmute(
        chromosome_group = w$chromosome_group,
        window_id = w$window_id,
        window_CHROM = w$CHROM,
        window_start = w$start,
        window_end = w$end,
        gene, gene_CHROM = CHROM, gene_start, gene_end,
        overlap_start = pmax(gene_start, w$start),
        overlap_end = pmin(gene_end, w$end),
        overlap_bp = overlap_end - overlap_start + 1
      )
  })
  if (ncol(hits) == 0) hit_template else hits
}

autosome_hits <- extract_overlaps(autosome_windows, genes)
sexchr_hits <- extract_overlaps(sexchr_windows, genes)
all_hits <- bind_rows(autosome_hits, sexchr_hits)

# Window-to-gene mapping: a gene can occur more than once if it overlaps more
# than one selected window.  This preserves the exact supporting windows.
write.table(autosome_hits, "selection_autosome_highconf_window_genes.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sexchr_hits, "selection_sexchr_highconf_window_genes.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all_hits, "selection_highconf_window_gene_overlaps.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Non-redundant gene lists, convenient for enrichment and annotation.
autosome_genes <- autosome_hits %>% distinct(gene, gene_CHROM, gene_start, gene_end)
sexchr_genes <- sexchr_hits %>% distinct(gene, gene_CHROM, gene_start, gene_end)
all_genes <- all_hits %>% distinct(gene, gene_CHROM, gene_start, gene_end)
write.table(autosome_genes, "selection_autosome_highconf_genes_unique.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sexchr_genes, "selection_sexchr_highconf_genes_unique.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(all_genes, "selection_highconf_genes_unique.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

summary_df <- tibble(
  chromosome_group = c("autosome", "sexchr", "combined"),
  highconf_windows = c(nrow(autosome_windows), nrow(sexchr_windows), nrow(autosome_windows) + nrow(sexchr_windows)),
  gene_window_overlaps = c(nrow(autosome_hits), nrow(sexchr_hits), nrow(all_hits)),
  unique_genes = c(nrow(autosome_genes), nrow(sexchr_genes), nrow(all_genes))
)
write.table(summary_df, "selection_highconf_gene_summary.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
print(summary_df)
cat("Finished extracting genes from high-confidence windows.\n")
