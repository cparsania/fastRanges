library(GenomicRanges)

canon_hits <- function(h) {
  x <- cbind(S4Vectors::queryHits(h), S4Vectors::subjectHits(h))
  if (nrow(x) == 0L) return(x)
  x[order(x[, 1], x[, 2]), , drop = FALSE]
}

test_that("fast_find_overlaps matches GenomicRanges for GRanges", {
  data(fast_ranges_example, package = "fastRanges")
  q <- fast_ranges_example$query
  s <- fast_ranges_example$subject

  ref <- GenomicRanges::findOverlaps(q, s, ignore.strand = FALSE)
  got <- fast_find_overlaps(q, s, threads = 2)

  expect_equal(canon_hits(got), canon_hits(ref))
})

test_that("ignore_strand behavior matches GenomicRanges", {
  data(fast_ranges_example, package = "fastRanges")
  q <- fast_ranges_example$query
  s <- fast_ranges_example$subject

  ref <- GenomicRanges::findOverlaps(q, s, ignore.strand = TRUE)
  got <- fast_find_overlaps(q, s, ignore_strand = TRUE, threads = 2)

  expect_equal(canon_hits(got), canon_hits(ref))
})

test_that("index-backed queries match direct subject queries", {
  data(fast_ranges_example, package = "fastRanges")
  q <- fast_ranges_example$query
  s <- fast_ranges_example$subject

  idx <- fast_build_index(s)

  from_subject <- fast_find_overlaps(q, s, threads = 2)
  from_index <- fast_find_overlaps(q, idx, threads = 2)

  expect_equal(canon_hits(from_subject), canon_hits(from_index))
})

test_that("thread counts produce deterministic output", {
  data(fast_ranges_example, package = "fastRanges")
  q <- fast_ranges_example$query
  s <- fast_ranges_example$subject

  h1 <- fast_find_overlaps(q, s, threads = 1, deterministic = TRUE)
  h4 <- fast_find_overlaps(q, s, threads = 4, deterministic = TRUE)

  expect_identical(S4Vectors::queryHits(h1), S4Vectors::queryHits(h4))
  expect_identical(S4Vectors::subjectHits(h1), S4Vectors::subjectHits(h4))
})

test_that("IRanges support matches IRanges::findOverlaps", {
  q <- IRanges::IRanges(start = c(1, 10, 20), end = c(5, 15, 25))
  s <- IRanges::IRanges(start = c(3, 8, 21), end = c(4, 12, 30))

  ref <- IRanges::findOverlaps(q, s)
  got <- fast_find_overlaps(q, s, threads = 2)

  expect_equal(canon_hits(got), canon_hits(ref))
})

test_that("indexed IRanges hits are thread-invariant when deterministic is FALSE", {
  set.seed(7)
  q <- IRanges::IRanges(
    start = sample.int(20000L, 600L, replace = TRUE),
    width = sample.int(150L, 600L, replace = TRUE)
  )
  s <- IRanges::IRanges(
    start = sample.int(20000L, 2000L, replace = TRUE),
    width = sample.int(150L, 2000L, replace = TRUE)
  )
  idx <- fast_build_index(s)

  h1 <- fast_find_overlaps(q, idx, threads = 1, deterministic = FALSE)
  h8 <- fast_find_overlaps(q, idx, threads = 8, deterministic = FALSE)

  expect_equal(canon_hits(h1), canon_hits(h8))
})

test_that("count and any wrappers agree with overlap hits", {
  data(fast_ranges_example, package = "fastRanges")
  q <- fast_ranges_example$query
  s <- fast_ranges_example$subject

  h <- fast_find_overlaps(q, s, threads = 2)
  counts <- fast_count_overlaps(q, s, threads = 2)
  any_hits <- fast_overlaps_any(q, s, threads = 2)

  expect_equal(counts, tabulate(S4Vectors::queryHits(h), nbins = length(q)))
  expect_equal(any_hits, counts > 0L)
})
