suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(clusterProfiler)
  library(fgsea)
  library(GO.db)
  library(AnnotationDbi)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

# =========================
# 0. 参数设置
# =========================
gtf_file <- "01.reference/complete.genomic.gtf"
egg_file <- "01.reference/FMdeer.emapper.annotations"

deseq_files <- c(
  tissue = "05.bulkRNA/05.deseq2/res_tissue.tsv",
  interaction = "05.bulkRNA/05.deseq2/res_interaction.tsv",
  gland_age = "05.bulkRNA/05.deseq2/res_gland_age.tsv"
)

outdir <- "05.bulkRNA/06.enrichment"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

padj_cutoff <- 0.05
lfc_cutoff  <- 1
min_gset_size <- 10
max_gset_size <- 500

# =========================
# 1. 工具函数
# =========================
extract_attr <- function(x, key) {
  pattern <- paste0(key, ' "([^"]*)"')
  m <- stringr::str_match(x, pattern)
  m[, 2]
}

parse_ratio <- function(x) {
  sapply(strsplit(x, "/"), function(z) {
    if (length(z) == 2) {
      as.numeric(z[1]) / as.numeric(z[2])
    } else {
      NA_real_
    }
  })
}

# 把 list 列安全转成字符列再写出
safe_write_table <- function(df, file) {
  if (is.null(df)) {
    df <- data.frame()
  }

  df2 <- as.data.frame(df, stringsAsFactors = FALSE)

  for (nm in colnames(df2)) {
    if (is.list(df2[[nm]])) {
      df2[[nm]] <- vapply(
        df2[[nm]],
        function(x) {
          if (length(x) == 0 || all(is.na(x))) {
            ""
          } else {
            paste(as.character(x), collapse = ",")
          }
        },
        character(1)
      )
    }
  }

  write.table(
    df2,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

plot_empty <- function(title_text, outfile) {
  p <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 1, y = 1, label = "No significant enrichment", size = 6) +
    ggplot2::theme_void() +
    ggplot2::ggtitle(title_text)
  ggplot2::ggsave(outfile, p, width = 8, height = 5)
}

save_barplot <- function(enrich_df, title_text, outfile, top_n = 15, term_col = "TERM") {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) {
    plot_empty(title_text, outfile)
    return(invisible(NULL))
  }

  plot_df <- enrich_df %>%
    dplyr::filter(!is.na(p.adjust)) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = top_n)

  if (nrow(plot_df) == 0) {
    plot_empty(title_text, outfile)
    return(invisible(NULL))
  }

  if (!term_col %in% colnames(plot_df)) {
    term_col <- "ID"
  }

  plot_df[[term_col]] <- factor(plot_df[[term_col]], levels = rev(unique(plot_df[[term_col]])))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[term_col]], y = -log10(p.adjust))) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "", y = "-log10(adj.P)", title = title_text)

  ggplot2::ggsave(outfile, p, width = 9, height = 6)
}

save_dotplot <- function(enrich_df, title_text, outfile, top_n = 20, term_col = "TERM") {
  if (is.null(enrich_df) || nrow(enrich_df) == 0) {
    plot_empty(title_text, outfile)
    return(invisible(NULL))
  }

  plot_df <- enrich_df %>%
    dplyr::filter(!is.na(p.adjust)) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = top_n)

  if (nrow(plot_df) == 0) {
    plot_empty(title_text, outfile)
    return(invisible(NULL))
  }

  if (!term_col %in% colnames(plot_df)) {
    term_col <- "ID"
  }

  if ("GeneRatio" %in% colnames(plot_df)) {
    plot_df$GeneRatio_num <- parse_ratio(plot_df$GeneRatio)
  } else {
    plot_df$GeneRatio_num <- NA_real_
  }

  plot_df[[term_col]] <- factor(plot_df[[term_col]], levels = rev(unique(plot_df[[term_col]])))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = GeneRatio_num, y = .data[[term_col]])) +
    ggplot2::geom_point(ggplot2::aes(size = Count, color = p.adjust)) +
    ggplot2::scale_color_gradient(low = "red", high = "blue", trans = "reverse") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "GeneRatio", y = "", title = title_text)

  ggplot2::ggsave(outfile, p, width = 9, height = 7)
}

# =========================
# 2. 读取并解析 GTF
# =========================
message("Step 1/6: Reading GTF...")

gtf <- read.delim(
  gtf_file,
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
  dplyr::mutate(
    gene_id = extract_attr(attribute, "gene_id"),
    transcript_id = extract_attr(attribute, "transcript_id"),
    gene_name = extract_attr(attribute, "gene"),
    protein_id = extract_attr(attribute, "protein_id"),
    gtf_description = extract_attr(attribute, "description")
  )

gene_info <- gtf_parsed %>%
  dplyr::filter(feature == "gene") %>%
  dplyr::select(gene_id, gene_name, gtf_description) %>%
  dplyr::distinct()

protein_map <- gtf_parsed %>%
  dplyr::filter(!is.na(gene_id), gene_id != "",
                !is.na(protein_id), protein_id != "") %>%
  dplyr::select(gene_id, transcript_id, protein_id) %>%
  dplyr::distinct()

gtf_map <- protein_map %>%
  dplyr::left_join(gene_info, by = "gene_id") %>%
  dplyr::distinct()

safe_write_table(gtf_map, file.path(outdir, "FMdeer_gtf_gene_protein_map.tsv"))

# =========================
# 3. 读取 eggNOG 注释
# =========================
message("Step 2/6: Reading eggNOG annotations...")

egg <- read.delim(
  egg_file,
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

# =========================
# 4. 合并注释
# =========================
message("Step 3/6: Building merged annotation...")

anno_all <- gtf_map %>%
  dplyr::left_join(egg, by = c("protein_id" = "query"))

gene_annotation <- anno_all %>%
  dplyr::mutate(
    final_symbol = dplyr::coalesce(gene_name, Preferred_name),
    final_description = dplyr::coalesce(Description, gtf_description)
  ) %>%
  dplyr::select(
    gene_id,
    final_symbol,
    gene_name,
    Preferred_name,
    transcript_id,
    protein_id,
    final_description,
    GOs,
    KEGG_ko,
    KEGG_Pathway,
    PFAMs,
    COG_category
  ) %>%
  dplyr::distinct()

safe_write_table(gene_annotation, file.path(outdir, "FMdeer_gene_annotation.tsv"))

# =========================
# 5. 构建 TERM2GENE / TERM2NAME
# =========================
message("Step 4/6: Building TERM2GENE...")

term2gene_go <- anno_all %>%
  dplyr::select(gene_id, GOs) %>%
  dplyr::filter(!is.na(gene_id), gene_id != "",
                !is.na(GOs), GOs != "", GOs != "-") %>%
  tidyr::separate_rows(GOs, sep = ",") %>%
  dplyr::rename(GO_ID = GOs) %>%
  dplyr::filter(GO_ID != "") %>%
  dplyr::distinct() %>%
  dplyr::select(GO_ID, gene_id)

safe_write_table(term2gene_go, file.path(outdir, "FMdeer_TERM2GENE_GO.tsv"))

message("Fetching GO term names from GO.db ...")

go_name <- AnnotationDbi::select(
  GO.db,
  keys = unique(term2gene_go$GO_ID),
  columns = c("TERM", "ONTOLOGY"),
  keytype = "GOID"
) %>%
  dplyr::distinct()

safe_write_table(go_name, file.path(outdir, "FMdeer_TERM2NAME_GO.tsv"))

term2gene_keggko <- anno_all %>%
  dplyr::select(gene_id, KEGG_ko) %>%
  dplyr::filter(!is.na(gene_id), gene_id != "",
                !is.na(KEGG_ko), KEGG_ko != "", KEGG_ko != "-") %>%
  tidyr::separate_rows(KEGG_ko, sep = ",") %>%
  dplyr::distinct() %>%
  dplyr::select(KEGG_ko, gene_id)

safe_write_table(term2gene_keggko, file.path(outdir, "FMdeer_TERM2GENE_KEGGKO.tsv"))

term2gene_keggpath <- anno_all %>%
  dplyr::select(gene_id, KEGG_Pathway) %>%
  dplyr::filter(!is.na(gene_id), gene_id != "",
                !is.na(KEGG_Pathway), KEGG_Pathway != "", KEGG_Pathway != "-") %>%
  tidyr::separate_rows(KEGG_Pathway, sep = ",") %>%
  dplyr::distinct() %>%
  dplyr::select(KEGG_Pathway, gene_id)

safe_write_table(term2gene_keggpath, file.path(outdir, "FMdeer_TERM2GENE_KEGGPathway.tsv"))

pathways_go <- split(term2gene_go$gene_id, term2gene_go$GO_ID)

# =========================
# 6. DESeq2 富集
# =========================
message("Step 5/6: Running enrichment on DESeq2 results...")

run_ora <- function(genes, bg_genes, term2gene, term_name_df = NULL, join_left = "ID", join_right = NULL) {
  if (length(genes) == 0) {
    return(data.frame())
  }

  ego <- tryCatch(
    clusterProfiler::enricher(
      gene = genes,
      universe = bg_genes,
      TERM2GENE = term2gene,
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    ),
    error = function(e) NULL
  )

  if (is.null(ego)) {
    return(data.frame())
  }

  df <- as.data.frame(ego)
  if (nrow(df) == 0) {
    return(df)
  }

  if (!is.null(term_name_df)) {
    if (is.null(join_right)) {
      join_right <- join_left
    }
    df <- dplyr::left_join(df, term_name_df, by = setNames(join_right, join_left))
  }

  df
}

run_one_deseq <- function(res_file, prefix) {
  message("Processing: ", prefix)

  subdir <- file.path(outdir, prefix)
  dir.create(subdir, showWarnings = FALSE, recursive = TRUE)

  res <- read.delim(
    res_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = ""
  )

  required_cols <- c("gene", "log2FoldChange", "stat", "pvalue", "padj")
  miss_cols <- setdiff(required_cols, colnames(res))
  if (length(miss_cols) > 0) {
    stop("Missing required columns in ", res_file, ": ", paste(miss_cols, collapse = ", "))
  }

  res <- res %>%
    dplyr::distinct(gene, .keep_all = TRUE)

  res_anno <- res %>%
    dplyr::left_join(gene_annotation, by = c("gene" = "gene_id"))

  safe_write_table(res_anno, file.path(subdir, paste0(prefix, "_annotated.tsv")))

  bg_genes <- res %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::pull(gene) %>%
    unique()

  deg_all <- res %>%
    dplyr::filter(!is.na(padj), padj < padj_cutoff) %>%
    dplyr::pull(gene) %>%
    unique()

  deg_up <- res %>%
    dplyr::filter(!is.na(padj), padj < padj_cutoff, log2FoldChange > lfc_cutoff) %>%
    dplyr::pull(gene) %>%
    unique()

  deg_down <- res %>%
    dplyr::filter(!is.na(padj), padj < padj_cutoff, log2FoldChange < -lfc_cutoff) %>%
    dplyr::pull(gene) %>%
    unique()

  stat_summary <- data.frame(
    comparison = prefix,
    bg_n = length(bg_genes),
    deg_all_n = length(deg_all),
    deg_up_n = length(deg_up),
    deg_down_n = length(deg_down),
    stringsAsFactors = FALSE
  )
  safe_write_table(stat_summary, file.path(subdir, paste0(prefix, "_DEG_summary.tsv")))

  # GO ORA
  go_all_df <- run_ora(deg_all, bg_genes, term2gene_go, go_name, "ID", "GOID")
  go_up_df <- run_ora(deg_up, bg_genes, term2gene_go, go_name, "ID", "GOID")
  go_down_df <- run_ora(deg_down, bg_genes, term2gene_go, go_name, "ID", "GOID")

  safe_write_table(go_all_df, file.path(subdir, paste0(prefix, "_GO_ORA_all.tsv")))
  safe_write_table(go_up_df, file.path(subdir, paste0(prefix, "_GO_ORA_up.tsv")))
  safe_write_table(go_down_df, file.path(subdir, paste0(prefix, "_GO_ORA_down.tsv")))

  save_barplot(go_all_df, paste0(prefix, " GO ORA (all)"),
               file.path(subdir, paste0(prefix, "_GO_ORA_all_barplot.pdf")))
  save_dotplot(go_all_df, paste0(prefix, " GO ORA (all)"),
               file.path(subdir, paste0(prefix, "_GO_ORA_all_dotplot.pdf")))

  save_barplot(go_up_df, paste0(prefix, " GO ORA (up)"),
               file.path(subdir, paste0(prefix, "_GO_ORA_up_barplot.pdf")))
  save_dotplot(go_up_df, paste0(prefix, " GO ORA (up)"),
               file.path(subdir, paste0(prefix, "_GO_ORA_up_dotplot.pdf")))

  save_barplot(go_down_df, paste0(prefix, " GO ORA (down)"),
               file.path(subdir, paste0(prefix, "_GO_ORA_down_barplot.pdf")))
  save_dotplot(go_down_df, paste0(prefix, " GO ORA (down)"),
               file.path(subdir, paste0(prefix, "_GO_ORA_down_dotplot.pdf")))

  # KEGG KO ORA
  keggko_all_df <- run_ora(deg_all, bg_genes, term2gene_keggko, NULL, "ID", NULL)
  safe_write_table(keggko_all_df, file.path(subdir, paste0(prefix, "_KEGGKO_ORA_all.tsv")))

  save_barplot(keggko_all_df, paste0(prefix, " KEGG KO ORA"),
               file.path(subdir, paste0(prefix, "_KEGGKO_ORA_barplot.pdf")),
               term_col = "Description")
  save_dotplot(keggko_all_df, paste0(prefix, " KEGG KO ORA"),
               file.path(subdir, paste0(prefix, "_KEGGKO_ORA_dotplot.pdf")),
               term_col = "Description")

  # GO GSEA，使用 fgseaMultilevel（不写 nperm）
  tmp <- res %>%
    dplyr::filter(!is.na(stat)) %>%
    dplyr::distinct(gene, .keep_all = TRUE)

  geneList <- tmp$stat
  names(geneList) <- tmp$gene
  geneList <- sort(geneList, decreasing = TRUE)

  fgsea_go <- tryCatch(
    fgsea::fgsea(
      pathways = pathways_go,
      stats = geneList,
      minSize = min_gset_size,
      maxSize = max_gset_size
    ) %>% as.data.frame(),
    error = function(e) data.frame()
  )

  if (nrow(fgsea_go) > 0) {
    fgsea_go <- fgsea_go %>%
      dplyr::left_join(go_name, by = c("pathway" = "GOID")) %>%
      dplyr::arrange(padj)
  }

  safe_write_table(fgsea_go, file.path(subdir, paste0(prefix, "_GO_fgsea.tsv")))

  if (nrow(fgsea_go) > 0) {
    fgsea_top <- fgsea_go %>%
      dplyr::filter(!is.na(padj)) %>%
      dplyr::arrange(padj) %>%
      dplyr::slice_head(n = 20)

    if (nrow(fgsea_top) > 0) {
      fgsea_top$term_label <- ifelse(is.na(fgsea_top$TERM), fgsea_top$pathway, fgsea_top$TERM)
      fgsea_top$term_label <- factor(fgsea_top$term_label, levels = rev(unique(fgsea_top$term_label)))

      p_fg <- ggplot2::ggplot(fgsea_top, ggplot2::aes(x = NES, y = term_label)) +
        ggplot2::geom_point(ggplot2::aes(size = size, color = padj)) +
        ggplot2::scale_color_gradient(low = "red", high = "blue", trans = "reverse") +
        ggplot2::theme_bw() +
        ggplot2::labs(title = paste0(prefix, " GO GSEA"), x = "NES", y = "")

      ggplot2::ggsave(file.path(subdir, paste0(prefix, "_GO_fgsea_dotplot.pdf")), p_fg, width = 9, height = 7)
    } else {
      plot_empty(paste0(prefix, " GO GSEA"), file.path(subdir, paste0(prefix, "_GO_fgsea_dotplot.pdf")))
    }
  } else {
    plot_empty(paste0(prefix, " GO GSEA"), file.path(subdir, paste0(prefix, "_GO_fgsea_dotplot.pdf")))
  }

  invisible(list(
    res = res,
    res_anno = res_anno,
    go_all = go_all_df,
    go_up = go_up_df,
    go_down = go_down_df,
    keggko_all = keggko_all_df,
    fgsea_go = fgsea_go
  ))
}

all_results <- list()
for (nm in names(deseq_files)) {
  all_results[[nm]] <- run_one_deseq(deseq_files[[nm]], nm)
}

# =========================
# 7. 汇总
# =========================
run_merged_up_down_enrichment <- function() {
  message("Step 6/7: Running enrichment on merged UP/DOWN genes...")

  subdir <- file.path(outdir, "merged_up_down")
  dir.create(subdir, showWarnings = FALSE, recursive = TRUE)

  merged_res <- dplyr::bind_rows(lapply(names(deseq_files), function(nm) {
    res_file <- deseq_files[[nm]]
    res <- read.delim(
      res_file,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      quote = ""
    )

    required_cols <- c("gene", "log2FoldChange", "padj")
    miss_cols <- setdiff(required_cols, colnames(res))
    if (length(miss_cols) > 0) {
      stop("Missing required columns in ", res_file, ": ", paste(miss_cols, collapse = ", "))
    }

    res %>%
      dplyr::distinct(gene, .keep_all = TRUE) %>%
      dplyr::mutate(
        comparison = nm,
        direction = dplyr::case_when(
          !is.na(padj) & padj < padj_cutoff & log2FoldChange >= lfc_cutoff ~ "UP",
          !is.na(padj) & padj < padj_cutoff & log2FoldChange <= -lfc_cutoff ~ "DOWN",
          TRUE ~ "NS"
        )
      )
  }))

  bg_genes <- merged_res %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::pull(gene) %>%
    unique()

  merged_up_detail <- merged_res %>%
    dplyr::filter(direction == "UP") %>%
    dplyr::select(gene, comparison, direction, log2FoldChange, padj) %>%
    dplyr::arrange(gene, comparison)

  merged_down_detail <- merged_res %>%
    dplyr::filter(direction == "DOWN") %>%
    dplyr::select(gene, comparison, direction, log2FoldChange, padj) %>%
    dplyr::arrange(gene, comparison)

  merged_up_genes <- merged_up_detail %>%
    dplyr::distinct(gene) %>%
    dplyr::arrange(gene)

  merged_down_genes <- merged_down_detail %>%
    dplyr::distinct(gene) %>%
    dplyr::arrange(gene)

  merged_up_summary <- merged_up_detail %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      n_comparisons = dplyr::n_distinct(comparison),
      comparisons = paste(sort(unique(comparison)), collapse = ";"),
      max_log2FoldChange = max(log2FoldChange, na.rm = TRUE),
      min_padj = min(padj, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_comparisons), min_padj, gene)

  merged_down_summary <- merged_down_detail %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      n_comparisons = dplyr::n_distinct(comparison),
      comparisons = paste(sort(unique(comparison)), collapse = ";"),
      min_log2FoldChange = min(log2FoldChange, na.rm = TRUE),
      min_padj = min(padj, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_comparisons), min_padj, gene)

  safe_write_table(merged_up_genes, file.path(subdir, "merged_UP_gene_ids.tsv"))
  safe_write_table(merged_down_genes, file.path(subdir, "merged_DOWN_gene_ids.tsv"))
  safe_write_table(merged_up_summary, file.path(subdir, "merged_UP_gene_ids_with_sources.tsv"))
  safe_write_table(merged_down_summary, file.path(subdir, "merged_DOWN_gene_ids_with_sources.tsv"))
  safe_write_table(merged_up_detail, file.path(subdir, "merged_UP_gene_long_table.tsv"))
  safe_write_table(merged_down_detail, file.path(subdir, "merged_DOWN_gene_long_table.tsv"))

  safe_write_table(
    data.frame(
      gene_set = c("merged_UP", "merged_DOWN"),
      bg_n = length(bg_genes),
      gene_n = c(nrow(merged_up_genes), nrow(merged_down_genes)),
      stringsAsFactors = FALSE
    ),
    file.path(subdir, "merged_UP_DOWN_gene_summary.tsv")
  )

  go_up_df <- run_ora(merged_up_genes$gene, bg_genes, term2gene_go, go_name, "ID", "GOID")
  go_down_df <- run_ora(merged_down_genes$gene, bg_genes, term2gene_go, go_name, "ID", "GOID")
  keggko_up_df <- run_ora(merged_up_genes$gene, bg_genes, term2gene_keggko, NULL, "ID", NULL)
  keggko_down_df <- run_ora(merged_down_genes$gene, bg_genes, term2gene_keggko, NULL, "ID", NULL)
  keggpath_up_df <- run_ora(merged_up_genes$gene, bg_genes, term2gene_keggpath, NULL, "ID", NULL)
  keggpath_down_df <- run_ora(merged_down_genes$gene, bg_genes, term2gene_keggpath, NULL, "ID", NULL)

  safe_write_table(go_up_df, file.path(subdir, "merged_UP_GO_ORA.tsv"))
  safe_write_table(go_down_df, file.path(subdir, "merged_DOWN_GO_ORA.tsv"))
  safe_write_table(keggko_up_df, file.path(subdir, "merged_UP_KEGGKO_ORA.tsv"))
  safe_write_table(keggko_down_df, file.path(subdir, "merged_DOWN_KEGGKO_ORA.tsv"))
  safe_write_table(keggpath_up_df, file.path(subdir, "merged_UP_KEGGPathway_ORA.tsv"))
  safe_write_table(keggpath_down_df, file.path(subdir, "merged_DOWN_KEGGPathway_ORA.tsv"))

  save_barplot(go_up_df, "merged UP GO ORA", file.path(subdir, "merged_UP_GO_ORA_barplot.pdf"))
  save_dotplot(go_up_df, "merged UP GO ORA", file.path(subdir, "merged_UP_GO_ORA_dotplot.pdf"))
  save_barplot(go_down_df, "merged DOWN GO ORA", file.path(subdir, "merged_DOWN_GO_ORA_barplot.pdf"))
  save_dotplot(go_down_df, "merged DOWN GO ORA", file.path(subdir, "merged_DOWN_GO_ORA_dotplot.pdf"))

  save_barplot(keggko_up_df, "merged UP KEGG KO ORA",
               file.path(subdir, "merged_UP_KEGGKO_ORA_barplot.pdf"),
               term_col = "Description")
  save_dotplot(keggko_up_df, "merged UP KEGG KO ORA",
               file.path(subdir, "merged_UP_KEGGKO_ORA_dotplot.pdf"),
               term_col = "Description")
  save_barplot(keggko_down_df, "merged DOWN KEGG KO ORA",
               file.path(subdir, "merged_DOWN_KEGGKO_ORA_barplot.pdf"),
               term_col = "Description")
  save_dotplot(keggko_down_df, "merged DOWN KEGG KO ORA",
               file.path(subdir, "merged_DOWN_KEGGKO_ORA_dotplot.pdf"),
               term_col = "Description")

  save_barplot(keggpath_up_df, "merged UP KEGG Pathway ORA",
               file.path(subdir, "merged_UP_KEGGPathway_ORA_barplot.pdf"),
               term_col = "Description")
  save_dotplot(keggpath_up_df, "merged UP KEGG Pathway ORA",
               file.path(subdir, "merged_UP_KEGGPathway_ORA_dotplot.pdf"),
               term_col = "Description")
  save_barplot(keggpath_down_df, "merged DOWN KEGG Pathway ORA",
               file.path(subdir, "merged_DOWN_KEGGPathway_ORA_barplot.pdf"),
               term_col = "Description")
  save_dotplot(keggpath_down_df, "merged DOWN KEGG Pathway ORA",
               file.path(subdir, "merged_DOWN_KEGGPathway_ORA_dotplot.pdf"),
               term_col = "Description")

  safe_write_table(
    data.frame(
      gene_set = c("merged_UP", "merged_DOWN", "merged_UP", "merged_DOWN", "merged_UP", "merged_DOWN"),
      database = c("GO", "GO", "KEGGKO", "KEGGKO", "KEGGPathway", "KEGGPathway"),
      input_gene_n = c(
        nrow(merged_up_genes), nrow(merged_down_genes),
        nrow(merged_up_genes), nrow(merged_down_genes),
        nrow(merged_up_genes), nrow(merged_down_genes)
      ),
      enriched_term_n = c(
        nrow(go_up_df), nrow(go_down_df),
        nrow(keggko_up_df), nrow(keggko_down_df),
        nrow(keggpath_up_df), nrow(keggpath_down_df)
      ),
      stringsAsFactors = FALSE
    ),
    file.path(subdir, "merged_UP_DOWN_GO_KEGG_enrichment_summary.tsv")
  )

  invisible(list(
    merged_up = merged_up_genes,
    merged_down = merged_down_genes,
    go_up = go_up_df,
    go_down = go_down_df,
    keggko_up = keggko_up_df,
    keggko_down = keggko_down_df,
    keggpath_up = keggpath_up_df,
    keggpath_down = keggpath_down_df
  ))
}

merged_results <- run_merged_up_down_enrichment()

message("Step 7/7: Writing summary...")

summary_df <- dplyr::bind_rows(lapply(names(deseq_files), function(nm) {
  subdir <- file.path(outdir, nm)
  f <- file.path(subdir, paste0(nm, "_DEG_summary.tsv"))
  read.delim(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
}))

safe_write_table(summary_df, file.path(outdir, "all_comparisons_DEG_summary.tsv"))

message("All done.")
message("Results written to: ", outdir)
