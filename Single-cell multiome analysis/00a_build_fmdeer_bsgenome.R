suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
  library(GenomeInfoDb)
  library(IRanges)
  library(S4Vectors)
  library(BSgenome)
  library(BSgenome.FMdeer.Custom.v1)
})

message("===== Start building custom ArchR annotation for FMdeer =====")

projDir   <- "06.multiome/03.archr"
fastaFile <- "01.reference/FMdeer.clean.fa"
faiFile   <- "01.reference/FMdeer.clean.fa.fai"
gtfFile   <- "01.reference/complete.clean.gtf"

bsgenomePkg <- "BSgenome.FMdeer.Custom.v1"

dir.create(projDir, recursive = TRUE, showWarnings = FALSE)

geneAnnoRDS   <- file.path(projDir, "FMdeer_geneAnnotation.rds")
genomeAnnoRDS <- file.path(projDir, "FMdeer_genomeAnnotation.rds")

requiredFiles <- c(fastaFile, faiFile, gtfFile)
missingFiles <- requiredFiles[!file.exists(requiredFiles)]
if (length(missingFiles) > 0) {
  stop("Missing files:\n", paste0("  - ", missingFiles, collapse = "\n"))
}

bsobj <- getBSgenome(bsgenomePkg)
message("Loaded BSgenome: ", bsgenomePkg)

# =========================================================
# 1. 读取 FAI
# =========================================================
fai <- fread(faiFile, header = FALSE)
colnames(fai)[1:min(5, ncol(fai))] <- c("seqnames", "seqlengths", "offset", "linebases", "linewidth")[1:min(5, ncol(fai))]
fai <- as.data.frame(fai)
fai$seqnames <- as.character(fai$seqnames)
fai$seqlengths <- as.numeric(fai$seqlengths)

chromSizes <- fai[, c("seqnames", "seqlengths"), drop = FALSE]
chromSizes <- chromSizes[!is.na(chromSizes$seqnames) & !is.na(chromSizes$seqlengths), , drop = FALSE]

message("Chromosomes from FAI:")
print(chromSizes$seqnames)

make_seqinfo_from_chr <- function(seqs, chromSizes) {
  seqs <- unique(as.character(seqs))
  seqs <- seqs[!is.na(seqs)]
  seqs <- seqs[seqs %in% chromSizes$seqnames]
  Seqinfo(
    seqnames   = seqs,
    seqlengths = chromSizes$seqlengths[match(seqs, chromSizes$seqnames)]
  )
}

clean_granges <- function(gr, chromSizes, objectName = "GRanges") {
  if (length(gr) == 0) return(gr)

  chr <- as.character(seqnames(gr))
  keep <- !is.na(chr) & chr %in% chromSizes$seqnames

  if (sum(!keep) > 0) {
    message(objectName, ": removing ", sum(!keep), " records with invalid seqnames.")
  }

  gr <- gr[keep]

  if (length(gr) == 0) {
    return(gr)
  }

  chr2 <- unique(as.character(seqnames(gr)))
  chr2 <- chr2[!is.na(chr2)]
  chr2 <- chr2[chr2 %in% chromSizes$seqnames]

  gr <- keepSeqlevels(gr, value = chr2, pruning.mode = "coarse")
  seqlevels(gr) <- chr2
  seqinfo(gr) <- make_seqinfo_from_chr(chr2, chromSizes)

  gr
}

# =========================================================
# 2. 导入 GTF
# =========================================================
gtf <- import(gtfFile)

if (length(gtf) == 0) {
  stop("Imported GTF is empty.")
}

gtf <- gtf[as.character(seqnames(gtf)) %in% chromSizes$seqnames]
if (length(gtf) == 0) {
  stop("No GTF records remain after filtering by FAI seqnames.")
}

gtf <- clean_granges(gtf, chromSizes, "gtf")

if (!"type" %in% colnames(mcols(gtf))) {
  stop("No 'type' column found in imported GTF.")
}

# =========================================================
# 3. 提取 gene / exon
# =========================================================
genesGR <- gtf[mcols(gtf)$type == "gene"]
exonsGR <- gtf[mcols(gtf)$type == "exon"]

if (length(genesGR) == 0) {
  txGR <- gtf[mcols(gtf)$type %in% c("transcript", "mRNA")]
  if (length(txGR) > 0 && "gene_id" %in% colnames(mcols(txGR))) {
    splitTx <- split(txGR, mcols(txGR)$gene_id)
    genesGR <- unlist(range(splitTx))
    mcols(genesGR)$gene_id <- names(splitTx)
  } else if (length(exonsGR) > 0 && "gene_id" %in% colnames(mcols(exonsGR))) {
    splitExon <- split(exonsGR, mcols(exonsGR)$gene_id)
    genesGR <- unlist(range(splitExon))
    mcols(genesGR)$gene_id <- names(splitExon)
  } else {
    stop("Cannot infer genes from transcript/exon records.")
  }
}

# =========================================================
# 4. 标准化元数据
# =========================================================
normalize_gene_metadata <- function(gr) {
  mc <- mcols(gr)
  nms <- colnames(mc)

  if (!"gene_id" %in% nms) {
    if ("ID" %in% nms) mc$gene_id <- mc$ID
    else if ("gene" %in% nms) mc$gene_id <- mc$gene
    else if ("Name" %in% nms) mc$gene_id <- mc$Name
    else mc$gene_id <- paste0("gene_", seq_along(gr))
  }

  if (!"gene_name" %in% nms) {
    if ("gene" %in% nms) mc$gene_name <- mc$gene
    else if ("Name" %in% nms) mc$gene_name <- mc$Name
    else if ("ID" %in% nms) mc$gene_name <- mc$ID
    else mc$gene_name <- mc$gene_id
  }

  mc$gene_id   <- as.character(mc$gene_id)
  mc$gene_name <- as.character(mc$gene_name)
  mc$symbol    <- mc$gene_name
  mcols(gr) <- mc
  gr
}

normalize_exon_metadata <- function(gr) {
  mc <- mcols(gr)
  nms <- colnames(mc)

  if (!"gene_id" %in% nms) {
    if ("Parent" %in% nms) mc$gene_id <- mc$Parent
    else if ("gene" %in% nms) mc$gene_id <- mc$gene
    else mc$gene_id <- NA_character_
  }

  if (!"gene_name" %in% nms) {
    if ("gene" %in% nms) mc$gene_name <- mc$gene
    else if ("Name" %in% nms) mc$gene_name <- mc$Name
    else mc$gene_name <- mc$gene_id
  }

  mc$gene_id   <- as.character(mc$gene_id)
  mc$gene_name <- as.character(mc$gene_name)
  mc$symbol    <- mc$gene_name
  mcols(gr) <- mc
  gr
}

genesGR <- normalize_gene_metadata(genesGR)
exonsGR <- normalize_exon_metadata(exonsGR)

# 清理非法染色体
genesGR <- clean_granges(genesGR, chromSizes, "genesGR")
exonsGR <- clean_granges(exonsGR, chromSizes, "exonsGR")

if (length(exonsGR) > 0) {
  exonsGR <- exonsGR[!is.na(mcols(exonsGR)$gene_id)]
  exonsGR <- clean_granges(exonsGR, chromSizes, "exonsGR_after_geneid_filter")
}

if (length(genesGR) == 0) {
  stop("No genes remain after cleaning.")
}

# =========================================================
# 5. 构建 TSS
# =========================================================
tssGR <- resize(genesGR, width = 1, fix = "start")
minusIdx <- which(as.character(strand(genesGR)) == "-")
if (length(minusIdx) > 0) {
  tssGR[minusIdx] <- resize(genesGR[minusIdx], width = 1, fix = "end")
}

mcols(tssGR)$gene_id   <- mcols(genesGR)$gene_id
mcols(tssGR)$gene_name <- mcols(genesGR)$gene_name
mcols(tssGR)$symbol    <- mcols(genesGR)$symbol

tssGR <- clean_granges(tssGR, chromSizes, "tssGR")

if (length(tssGR) == 0) {
  stop("No TSS remain after cleaning.")
}

# =========================================================
# 6. 构建 genomeAnnotation
# =========================================================
chromSizesGR <- GRanges(
  seqnames = chromSizes$seqnames,
  ranges   = IRanges(start = 1, end = chromSizes$seqlengths)
)
seqlevels(chromSizesGR) <- chromSizes$seqnames
seqinfo(chromSizesGR) <- make_seqinfo_from_chr(chromSizes$seqnames, chromSizes)

blacklistGR <- GRanges()

genomeAnnotation <- SimpleList(
  genome = bsgenomePkg,
  chromSizes = chromSizesGR,
  blacklist = blacklistGR
)

# =========================================================
# 7. 构建 geneAnnotation
# =========================================================
message("Final seqnames in genesGR:")
print(sort(unique(as.character(seqnames(genesGR)))))

message("Final seqnames in exonsGR:")
print(sort(unique(as.character(seqnames(exonsGR)))))

message("Final seqnames in tssGR:")
print(sort(unique(as.character(seqnames(tssGR)))))

if (any(is.na(as.character(seqnames(genesGR))))) stop("genesGR still contains NA seqnames.")
if (any(is.na(as.character(seqnames(exonsGR))))) stop("exonsGR still contains NA seqnames.")
if (any(is.na(as.character(seqnames(tssGR))))) stop("tssGR still contains NA seqnames.")

geneAnnotation <- createGeneAnnotation(
  genes = genesGR,
  exons = exonsGR,
  TSS   = tssGR
)

saveRDS(geneAnnotation, geneAnnoRDS)
saveRDS(genomeAnnotation, genomeAnnoRDS)

message("Saved geneAnnotation to: ", geneAnnoRDS)
message("Saved genomeAnnotation to: ", genomeAnnoRDS)
message("===== Done =====")
