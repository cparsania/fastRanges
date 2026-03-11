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
