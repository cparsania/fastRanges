library(GenomicRanges)

test_that("fast_nearest returns expected columns", {
  data(fast_ranges_example, package = "fastRanges")
  q <- fast_ranges_example$query
  s <- fast_ranges_example$subject

  x <- fast_nearest(q, s)

  expect_true(all(c("query_id", "subject_id", "distance") %in% names(x)))
  expect_true(nrow(as.data.frame(x)) > 0)
})

test_that("nearest helper variants return Bioconductor-compatible outputs", {
  data(fast_ranges_example, package = "fastRanges")
  q <- fast_ranges_example$query
  s <- fast_ranges_example$subject

  expect_equal(fast_distance_to_nearest(q, s), fast_nearest(q, s))
  expect_equal(
    fast_precede(q, s, ignore_strand = TRUE),
    GenomicRanges::precede(q, s, ignore.strand = TRUE)
  )
  expect_equal(
    fast_follow(q, s, ignore_strand = TRUE),
    GenomicRanges::follow(q, s, ignore.strand = TRUE)
  )
})

test_that("nearest family handles zero-hit and seqlevel-mismatch cases", {
  q <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = c(1L, 10L), width = 5L), strand = c("+", "-"))
  s <- GenomicRanges::GRanges("chr9", IRanges::IRanges(start = c(1L, 10L), width = 5L), strand = c("+", "-"))

  dn <- suppressWarnings(fast_distance_to_nearest(q, s, ignore_strand = FALSE, threads = 2))
  nr <- suppressWarnings(fast_nearest(q, s, ignore_strand = FALSE, threads = 2))

  expect_equal(nrow(as.data.frame(dn)), 0L)
  expect_equal(nrow(as.data.frame(nr)), 0L)
  expect_identical(
    suppressWarnings(fast_precede(q, s, ignore_strand = FALSE, threads = 2)),
    suppressWarnings(GenomicRanges::precede(q, s, ignore.strand = FALSE))
  )
  expect_identical(
    suppressWarnings(fast_follow(q, s, ignore_strand = FALSE, threads = 2)),
    suppressWarnings(GenomicRanges::follow(q, s, ignore.strand = FALSE))
  )
})

test_that("nearest family respects all-star strand edge case", {
  q <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = c(10L, 30L), width = 5L), strand = "*")
  s <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = c(1L, 20L), width = 5L), strand = "*")

  got <- fast_distance_to_nearest(q, s, ignore_strand = FALSE, threads = 2)
  ref <- GenomicRanges::distanceToNearest(q, s, ignore.strand = FALSE)

  expect_identical(as.integer(got$query_id), as.integer(S4Vectors::queryHits(ref)))
  expect_identical(as.integer(got$subject_id), as.integer(S4Vectors::subjectHits(ref)))
  expect_equal(as.numeric(got$distance), as.numeric(S4Vectors::mcols(ref)$distance))
  expect_equal(
    fast_precede(q, s, ignore_strand = FALSE, threads = 2),
    GenomicRanges::precede(q, s, ignore.strand = FALSE)
  )
  expect_equal(
    fast_follow(q, s, ignore_strand = FALSE, threads = 2),
    GenomicRanges::follow(q, s, ignore.strand = FALSE)
  )
})
