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

test_that("indexed count kernel matches hit tabulation across threads", {
  set.seed(11)
  q <- IRanges::IRanges(
    start = sample.int(30000L, 800L, replace = TRUE),
    width = sample.int(180L, 800L, replace = TRUE)
  )
  s <- IRanges::IRanges(
    start = sample.int(30000L, 3000L, replace = TRUE),
    width = sample.int(180L, 3000L, replace = TRUE)
  )
  idx <- fast_build_index(s)

  c1 <- fast_count_overlaps(q, idx, threads = 1, deterministic = FALSE)
  c8 <- fast_count_overlaps(q, idx, threads = 8, deterministic = FALSE)
  h <- fast_find_overlaps(q, idx, threads = 4, deterministic = FALSE)
  href <- tabulate(S4Vectors::queryHits(h), nbins = length(q))

  expect_equal(c1, href)
  expect_equal(c8, href)
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

test_that("within semantics with max_gap match GenomicRanges", {
  q <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = 10L, end = 20L))
  s <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(
      start = c(5L, 5L, 8L, 10L, 11L),
      end = c(20L, 25L, 23L, 25L, 20L)
    )
  )

  ref0 <- GenomicRanges::findOverlaps(
    q, s,
    type = "within",
    maxgap = 0L,
    minoverlap = 0L,
    ignore.strand = FALSE
  )
  got0 <- fast_find_overlaps(
    q, s,
    type = "within",
    max_gap = 0L,
    min_overlap = 0L,
    ignore_strand = FALSE,
    threads = 2
  )

  ref5 <- GenomicRanges::findOverlaps(
    q, s,
    type = "within",
    maxgap = 5L,
    minoverlap = 0L,
    ignore.strand = FALSE
  )
  got5 <- fast_find_overlaps(
    q, s,
    type = "within",
    max_gap = 5L,
    min_overlap = 0L,
    ignore_strand = FALSE,
    threads = 2
  )

  expect_equal(canon_hits(got0), canon_hits(ref0))
  expect_equal(canon_hits(got5), canon_hits(ref5))
})

test_that("equal semantics with max_gap match GenomicRanges", {
  q <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = 10L, end = 20L))
  s <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(
      start = c(5L, 5L, 10L, 11L, 10L),
      end = c(20L, 21L, 25L, 20L, 20L)
    )
  )

  ref0 <- GenomicRanges::findOverlaps(
    q, s,
    type = "equal",
    maxgap = 0L,
    minoverlap = 0L,
    ignore.strand = FALSE
  )
  got0 <- fast_find_overlaps(
    q, s,
    type = "equal",
    max_gap = 0L,
    min_overlap = 0L,
    ignore_strand = FALSE,
    threads = 2
  )

  ref5 <- GenomicRanges::findOverlaps(
    q, s,
    type = "equal",
    maxgap = 5L,
    minoverlap = 0L,
    ignore.strand = FALSE
  )
  got5 <- fast_find_overlaps(
    q, s,
    type = "equal",
    max_gap = 5L,
    min_overlap = 0L,
    ignore_strand = FALSE,
    threads = 2
  )

  expect_equal(canon_hits(got0), canon_hits(ref0))
  expect_equal(canon_hits(got5), canon_hits(ref5))
})

test_that("empty GRanges inputs return empty-compatible outputs", {
  q <- GenomicRanges::GRanges()
  s <- GenomicRanges::GRanges()

  h <- fast_find_overlaps(q, s, threads = 2)
  counts <- fast_count_overlaps(q, s, threads = 2)
  any_hits <- fast_overlaps_any(q, s, threads = 2)

  expect_s4_class(h, "Hits")
  expect_length(h, 0L)
  expect_identical(as.integer(counts), integer())
  expect_identical(as.logical(any_hits), logical())
})

test_that("single-range no-hit cases match GenomicRanges", {
  q <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = 10L, end = 20L), strand = "+")
  s <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = 30L, end = 40L), strand = "+")

  ref <- GenomicRanges::findOverlaps(q, s, ignore.strand = FALSE)
  got <- fast_find_overlaps(q, s, threads = 2)

  expect_equal(canon_hits(got), canon_hits(ref))
  expect_identical(as.integer(fast_count_overlaps(q, s, threads = 2)), 0L)
  expect_identical(as.logical(fast_overlaps_any(q, s, threads = 2)), FALSE)
})

test_that("no shared seqlevels produce zero hits", {
  q <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = c(1L, 10L), width = 5L), strand = c("+", "-"))
  s <- GenomicRanges::GRanges("chr9", IRanges::IRanges(start = c(1L, 10L), width = 5L), strand = c("+", "-"))

  got <- suppressWarnings(fast_find_overlaps(q, s, threads = 2))
  expect_length(got, 0L)
  expect_identical(as.integer(fast_count_overlaps(q, s, threads = 2)), c(0L, 0L))
  expect_identical(as.logical(fast_overlaps_any(q, s, threads = 2)), c(FALSE, FALSE))
})

test_that("all-star strand inputs behave like ignore_strand TRUE baseline", {
  q <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = c(5L, 15L), width = c(5L, 5L)),
    strand = "*"
  )
  s <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = c(1L, 15L), width = c(10L, 5L)),
    strand = "*"
  )

  ref <- GenomicRanges::findOverlaps(q, s, ignore.strand = TRUE)
  got_false <- fast_find_overlaps(q, s, ignore_strand = FALSE, threads = 2)
  got_true <- fast_find_overlaps(q, s, ignore_strand = TRUE, threads = 2)

  expect_equal(canon_hits(got_false), canon_hits(ref))
  expect_equal(canon_hits(got_true), canon_hits(ref))
})

test_that("deterministic indexed GRanges output is stable across threads", {
  set.seed(99)
  q <- GenomicRanges::GRanges(
    seqnames = sample(c("chr1", "chr2"), 300L, replace = TRUE),
    ranges = IRanges::IRanges(
      start = sample.int(5000L, 300L, replace = TRUE),
      width = sample.int(50L, 300L, replace = TRUE)
    ),
    strand = sample(c("+", "-", "*"), 300L, replace = TRUE)
  )
  s <- GenomicRanges::GRanges(
    seqnames = sample(c("chr1", "chr2"), 900L, replace = TRUE),
    ranges = IRanges::IRanges(
      start = sample.int(5000L, 900L, replace = TRUE),
      width = sample.int(50L, 900L, replace = TRUE)
    ),
    strand = sample(c("+", "-", "*"), 900L, replace = TRUE)
  )
  idx <- fast_build_index(s)

  h1 <- fast_find_overlaps(q, idx, threads = 1, deterministic = TRUE)
  h8 <- fast_find_overlaps(q, idx, threads = 8, deterministic = TRUE)

  expect_identical(S4Vectors::queryHits(h1), S4Vectors::queryHits(h8))
  expect_identical(S4Vectors::subjectHits(h1), S4Vectors::subjectHits(h8))
})
