if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
  stop("GenomicRanges is required to build example data")
}
if (!requireNamespace("IRanges", quietly = TRUE)) {
  stop("IRanges is required to build example data")
}

query <- GenomicRanges::GRanges(
  seqnames = c("chr1", "chr1", "chr1", "chr2", "chr2", "chr3"),
  ranges = IRanges::IRanges(start = c(100, 180, 350, 40, 300, 10), width = c(60, 90, 70, 40, 80, 25)),
  strand = c("+", "+", "-", "*", "+", "-")
)
S4Vectors::mcols(query)$query_id <- paste0("q", seq_along(query))
S4Vectors::mcols(query)$score <- c(8L, 10L, 6L, 3L, 9L, 5L)

subject <- GenomicRanges::GRanges(
  seqnames = c("chr1", "chr1", "chr1", "chr2", "chr2", "chr3", "chr3"),
  ranges = IRanges::IRanges(start = c(90, 210, 320, 20, 330, 1, 80), width = c(80, 100, 50, 70, 90, 40, 35)),
  strand = c("+", "-", "-", "+", "*", "-", "+")
)
S4Vectors::mcols(subject)$gene_id <- c("gA", "gB", "gC", "gD", "gE", "gF", "gG")
S4Vectors::mcols(subject)$biotype <- c("promoter", "enhancer", "enhancer", "promoter", "gene_body", "promoter", "enhancer")

fast_ranges_example <- list(query = query, subject = subject)

if (!dir.exists("data")) dir.create("data", recursive = TRUE)
save(fast_ranges_example, file = "data/fast_ranges_example.rda", compress = "xz")

if (!dir.exists("inst/extdata")) dir.create("inst/extdata", recursive = TRUE)

query_bed <- data.frame(
  seqnames = as.character(GenomicRanges::seqnames(query)),
  start0 = IRanges::start(query) - 1L,
  end = IRanges::end(query),
  name = S4Vectors::mcols(query)$query_id,
  score = S4Vectors::mcols(query)$score,
  strand = as.character(GenomicRanges::strand(query)),
  stringsAsFactors = FALSE
)

subject_bed <- data.frame(
  seqnames = as.character(GenomicRanges::seqnames(subject)),
  start0 = IRanges::start(subject) - 1L,
  end = IRanges::end(subject),
  name = S4Vectors::mcols(subject)$gene_id,
  score = 0L,
  strand = as.character(GenomicRanges::strand(subject)),
  stringsAsFactors = FALSE
)

utils::write.table(query_bed, file = "inst/extdata/query_peaks.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
utils::write.table(subject_bed, file = "inst/extdata/subject_genes.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
