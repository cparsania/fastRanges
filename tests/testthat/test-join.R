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
