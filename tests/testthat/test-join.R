library(GenomicRanges)

test_that("inner join row count equals overlap count", {
  data(fast_ranges_example, package = "fastRanges")
  q <- fast_ranges_example$query
  s <- fast_ranges_example$subject

  h <- fast_find_overlaps(q, s, threads = 2)
  joined <- fast_overlap_join(q, s, join = "inner", threads = 2)

  expect_equal(nrow(joined), length(h))
  expect_true(all(c("query_id", "subject_id") %in% names(joined)))
})

test_that("left join keeps all query ranges", {
  data(fast_ranges_example, package = "fastRanges")
  q <- fast_ranges_example$query
  s <- fast_ranges_example$subject

  joined <- fast_overlap_join(q, s, join = "left", threads = 2)

  expect_true(all(seq_len(length(q)) %in% joined$query_id))
  expect_true(any(is.na(joined$subject_id)) || nrow(joined) >= length(q))
})

test_that("join prefixes are applied consistently", {
  data(fast_ranges_example, package = "fastRanges")
  q <- fast_ranges_example$query
  s <- fast_ranges_example$subject

  joined <- fast_overlap_join(
    q,
    s,
    join = "inner",
    query_prefix = "q_",
    subject_prefix = "s_",
    threads = 2
  )

  expect_true(any(startsWith(names(joined), "q_")))
  expect_true(any(startsWith(names(joined), "s_")))
})

test_that("join family handles zero-hit and empty-query edge cases", {
  q <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = c(1L, 10L), width = 5L), strand = c("+", "-"))
  s <- GenomicRanges::GRanges("chr2", IRanges::IRanges(start = c(1L, 10L), width = 5L), strand = c("+", "-"))

  inner <- suppressWarnings(fast_inner_overlap_join(q, s, threads = 2))
  left <- suppressWarnings(fast_left_overlap_join(q, s, threads = 2))
  semi <- suppressWarnings(fast_semi_overlap_join(q, s, threads = 2))
  anti <- suppressWarnings(fast_anti_overlap_join(q, s, threads = 2))

  expect_equal(nrow(inner), 0L)
  expect_equal(nrow(left), length(q))
  expect_equal(nrow(semi), 0L)
  expect_equal(sort(anti$query_id), seq_len(length(q)))

  q_empty <- GenomicRanges::GRanges()
  inner_empty <- fast_inner_overlap_join(q_empty, s, threads = 2)
  left_empty <- fast_left_overlap_join(q_empty, s, threads = 2)
  expect_equal(nrow(inner_empty), 0L)
  expect_equal(nrow(left_empty), 0L)
})

test_that("deterministic join output is stable across thread counts", {
  set.seed(101)
  q <- GenomicRanges::GRanges(
    seqnames = sample(c("chr1", "chr2"), 120L, replace = TRUE),
    ranges = IRanges::IRanges(
      start = sample.int(2000L, 120L, replace = TRUE),
      width = sample.int(40L, 120L, replace = TRUE)
    ),
    strand = sample(c("+", "-", "*"), 120L, replace = TRUE)
  )
  s <- GenomicRanges::GRanges(
    seqnames = sample(c("chr1", "chr2"), 400L, replace = TRUE),
    ranges = IRanges::IRanges(
      start = sample.int(2000L, 400L, replace = TRUE),
      width = sample.int(40L, 400L, replace = TRUE)
    ),
    strand = sample(c("+", "-", "*"), 400L, replace = TRUE)
  )

  j1 <- fast_overlap_join(q, s, join = "inner", threads = 1, deterministic = TRUE)
  j8 <- fast_overlap_join(q, s, join = "inner", threads = 8, deterministic = TRUE)

  expect_identical(j1$query_id, j8$query_id)
  expect_identical(j1$subject_id, j8$subject_id)
})
