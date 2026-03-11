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
  s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
  idx <- fast_build_index(s)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  fast_save_index(idx, path)
  idx2 <- fast_load_index(path)

  expect_s3_class(idx2, "fast_ranges_index")
  expect_equal(idx2$subject_n, idx$subject_n)

  stats <- fast_index_stats(idx2)
  expect_true(all(c("subject_n", "partition_n", "seqlevel_n", "index_size_mb") %in% names(stats)))

  detailed <- fast_index_stats(idx2, detailed = TRUE)
  expect_true(all(c("summary", "partitions") %in% names(detailed)))
})
