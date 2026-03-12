#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(fastRanges)
  library(GenomicRanges)
  library(IRanges)
})

make_random_granges <- function(n, seqlevels = paste0("chr", 1:6), seed = 1L) {
  set.seed(seed)
  starts <- sample.int(5e6, n, replace = TRUE)
  widths <- sample(50:600, n, replace = TRUE)
  GRanges(
    seqnames = sample(seqlevels, n, replace = TRUE),
    ranges = IRanges(start = starts, width = widths),
    strand = sample(c("+", "-", "*"), n, replace = TRUE)
  )
}

query_n <- as.integer(Sys.getenv("FASTRANGES_BENCH_QUERY", unset = "200000"))
subject_n <- as.integer(Sys.getenv("FASTRANGES_BENCH_SUBJECT", unset = "500000"))
threads <- as.integer(Sys.getenv("FASTRANGES_BENCH_THREADS", unset = "4"))

cat("Building synthetic benchmark data...\n")
query <- make_random_granges(query_n, seed = 1001L)
subject <- make_random_granges(subject_n, seed = 2001L)
index <- fast_build_index(subject)

cat(sprintf("query_n=%d subject_n=%d threads=%d\n", query_n, subject_n, threads))

cat("Running GenomicRanges::findOverlaps()...\n")
base_time <- system.time({
  h_base <- GenomicRanges::findOverlaps(query, subject, ignore.strand = FALSE)
})

cat("Running fastRanges::fast_find_overlaps()...\n")
fast_time <- system.time({
  h_fast <- fast_find_overlaps(query, index, ignore_strand = FALSE, threads = threads)
})

cat("\nResults\n")
print(rbind(
  genomic_ranges = base_time[1:3],
  fast_ranges = fast_time[1:3]
))

cat(sprintf("\nHits: base=%d fast=%d\n", length(h_base), length(h_fast)))
