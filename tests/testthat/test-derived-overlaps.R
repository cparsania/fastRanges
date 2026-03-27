library(GenomicRanges)

canon_hits <- function(h) {
  x <- cbind(S4Vectors::queryHits(h), S4Vectors::subjectHits(h))
  if (nrow(x) == 0L) {
    return(x)
  }
  x[order(x[, 1], x[, 2]), , drop = FALSE]
}

test_that("fast_self_overlaps matches Bioconductor semantics", {
  x <- IRanges::IRanges(start = c(1, 3, 10, 11), end = c(5, 8, 12, 14))

  ref_all <- IRanges::findOverlaps(x, x)
  qh <- S4Vectors::queryHits(ref_all)
  sh <- S4Vectors::subjectHits(ref_all)
  keep <- qh != sh & qh < sh
  ref <- S4Vectors::Hits(from = qh[keep], to = sh[keep], nLnode = length(x), nRnode = length(x))
  got <- fast_self_overlaps(x)

  expect_equal(canon_hits(got), canon_hits(ref))
})

test_that("fast_cluster_overlaps assigns connected components", {
  x <- IRanges::IRanges(start = c(1, 3, 10, 11, 20), end = c(5, 8, 12, 14, 22))
  cl <- fast_cluster_overlaps(x)

  expect_length(cl, length(x))
  expect_equal(cl[1], cl[2])
  expect_equal(cl[3], cl[4])
  expect_true(cl[1] != cl[3])
  expect_true(cl[3] != cl[5])

  df <- fast_cluster_overlaps(x, return = "data.frame")
  expect_true(all(c("range_id", "cluster_id", "cluster_size") %in% names(df)))
  expect_equal(nrow(df), length(x))
})

test_that("fast_window_count_overlaps returns expected window counts", {
  q <- IRanges::IRanges(start = c(1, 12), end = c(10, 20))
  s <- IRanges::IRanges(start = c(2, 5, 15), end = c(3, 6, 16))

  got <- fast_window_count_overlaps(q, s, window_width = 5L, step_width = 5L)

  expect_equal(got$start, c(1L, 6L, 11L, 16L))
  expect_equal(got$end, c(5L, 10L, 15L, 20L))
  expect_equal(got$overlap_count, c(2L, 1L, 1L, 1L))
})

test_that("overlap iterator collects the same hits as direct call", {
  data(fast_ranges_example, package = "fastRanges")
  q <- fast_ranges_example$query
  s <- fast_ranges_example$subject

  iter <- fast_find_overlaps_iter(q, s, chunk_size = 25L, threads = 2)
  expect_true(fast_iter_has_next(iter))

  h_iter <- fast_iter_collect(iter)
  expect_false(fast_iter_has_next(iter))

  fast_iter_reset(iter)
  expect_true(fast_iter_has_next(iter))

  h_direct <- fast_find_overlaps(q, s, threads = 2)
  expect_equal(canon_hits(h_iter), canon_hits(h_direct))
})

test_that("index save/load/stats round-trip works", {
  s <- IRanges::IRanges(start = c(3, 9, 18, 30), end = c(6, 14, 22, 29))
  idx <- fast_build_index(s)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  fast_save_index(idx, path)
  idx2 <- fast_load_index(path)

  expect_s3_class(idx2, "fast_ranges_index")
  expect_equal(idx2$subject_n, idx$subject_n)
  expect_identical(idx$has_empty_ranges, TRUE)
  expect_identical(idx2$has_empty_ranges, TRUE)
  expect_identical(idx$has_circular_sequences, FALSE)
  expect_identical(idx2$has_circular_sequences, FALSE)

  stats <- fast_index_stats(idx2)
  expect_true(all(c("subject_n", "partition_n", "seqlevel_n", "index_size_mb") %in% names(stats)))

  detailed <- fast_index_stats(idx2, detailed = TRUE)
  expect_true(all(c("summary", "partitions") %in% names(detailed)))

  legacy_idx <- idx
  legacy_idx$has_empty_ranges <- NULL
  legacy_idx$has_circular_sequences <- NULL
  saveRDS(legacy_idx, path)
  idx3 <- fast_load_index(path)
  expect_identical(idx3$has_empty_ranges, TRUE)
  expect_identical(idx3$has_circular_sequences, FALSE)
})

test_that("fast_gaps matches GenomicRanges when GRanges bounds are implicit", {
  x <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr2"),
    ranges = IRanges::IRanges(start = c(1L, 5L, 10L), end = c(3L, 8L, 12L)),
    strand = c("+", "+", "*")
  )
  seqlengths(x) <- c(chr1 = 20L, chr2 = 20L)

  expect_equal(
    fast_gaps(x, ignore_strand = TRUE),
    GenomicRanges::gaps(x, ignore.strand = TRUE)
  )
})

test_that("self overlaps and clustering handle empty and isolated inputs", {
  empty_ir <- IRanges::IRanges()
  empty_hits <- fast_self_overlaps(empty_ir)
  empty_clusters <- fast_cluster_overlaps(empty_ir)
  empty_df <- fast_cluster_overlaps(empty_ir, return = "data.frame")

  expect_length(empty_hits, 0L)
  expect_identical(empty_clusters, integer())
  expect_equal(nrow(empty_df), 0L)

  isolated <- IRanges::IRanges(start = c(1L, 20L, 40L), width = c(3L, 3L, 3L))
  cl <- fast_cluster_overlaps(isolated)
  expect_identical(cl, c(1L, 2L, 3L))
})

test_that("window count handles no-hit and strand-sensitive cases", {
  q <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(1L, 20L), end = c(10L, 30L)),
    strand = c("+", "-")
  )
  s <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(50L, 60L), end = c(55L, 65L)),
    strand = c("+", "-")
  )

  got_zero <- fast_window_count_overlaps(q, s, window_width = 10L, step_width = 10L, threads = 2)
  expect_true(all(got_zero$overlap_count == 0L))

  s2 <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(2L, 21L), end = c(4L, 22L)),
    strand = c("-", "+")
  )
  ref <- GenomicRanges::countOverlaps(
    GenomicRanges::GRanges("chr1", IRanges::IRanges(start = c(1L, 11L, 21L), width = 10L), strand = "*"),
    s2,
    ignore.strand = TRUE
  )
  got <- fast_window_count_overlaps(
    q,
    s2,
    window_width = 10L,
    step_width = 10L,
    ignore_strand = TRUE,
    threads = 2
  )

  expect_equal(got$overlap_count, as.integer(ref))
})
