library(GenomicRanges)

test_that("fast_count_overlaps_by_group matches reference tabulation", {
  q <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 10, 20), width = 5))
  s <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(2, 9, 18, 22), width = 5))
  S4Vectors::mcols(s)$grp <- c("A", "B", "A", NA_character_)

  got <- fast_count_overlaps_by_group(q, s, group_col = "grp", include_na_group = TRUE)

  h <- GenomicRanges::findOverlaps(q, s)
  qh <- S4Vectors::queryHits(h)
  sh <- S4Vectors::subjectHits(h)
  grp <- as.character(S4Vectors::mcols(s)$grp[sh])
  grp[is.na(grp)] <- "<NA>"
  lev <- sort(unique(grp))
  ref <- matrix(0L, nrow = length(q), ncol = length(lev), dimnames = list(NULL, lev))
  tab <- table(qh, factor(grp, levels = lev))
  ref[as.integer(rownames(tab)), ] <- tab

  expect_equal(got[, colnames(ref), drop = FALSE], ref)
})

test_that("fast_overlap_aggregate count and sum are correct", {
  q <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 10, 20), width = 5))
  s <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(2, 9, 18, 22), width = 5))
  S4Vectors::mcols(s)$score <- c(2, 5, 1, 4)

  count_got <- fast_overlap_aggregate(q, s, fun = "count")
  sum_got <- fast_overlap_aggregate(q, s, value_col = "score", fun = "sum")

  h <- GenomicRanges::findOverlaps(q, s)
  qh <- S4Vectors::queryHits(h)
  sh <- S4Vectors::subjectHits(h)
  scores <- S4Vectors::mcols(s)$score[sh]

  count_ref <- tabulate(qh, nbins = length(q))
  sum_ref <- rep(NA_real_, length(q))
  for (i in seq_len(length(q))) {
    vals <- scores[qh == i]
    if (length(vals) > 0L) {
      sum_ref[i] <- sum(vals)
    }
  }

  expect_equal(as.integer(count_got), count_ref)
  expect_equal(sum_got, sum_ref)
})

test_that("range operation wrappers match IRanges implementations", {
  x <- IRanges::IRanges(start = c(1, 4, 10), end = c(5, 8, 12))
  y <- IRanges::IRanges(start = c(3, 20), end = c(8, 21))

  expect_equal(fast_reduce(x), IRanges::reduce(x))
  expect_equal(fast_disjoin(x), IRanges::disjoin(x))
  expect_equal(fast_gaps(x, start = 1L, end = 15L), IRanges::gaps(x, start = 1L, end = 15L))
  expect_equal(fast_range_union(x, y), IRanges::union(x, y))
  expect_equal(fast_range_intersect(x, y), IRanges::intersect(x, y))
  expect_equal(fast_range_setdiff(x, y), IRanges::setdiff(x, y))
})

test_that("join variant wrappers return expected subsets", {
  q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
  s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))

  inner <- fast_inner_overlap_join(q, s)
  left <- fast_left_overlap_join(q, s)
  semi <- fast_semi_overlap_join(q, s)
  anti <- fast_anti_overlap_join(q, s)

  expect_equal(inner, fast_overlap_join(q, s, join = "inner"))
  expect_equal(left, fast_overlap_join(q, s, join = "left"))
  expect_true(all(semi$overlap_count > 0L))
  expect_true(all(anti$overlap_count == 0L))
  expect_setequal(c(semi$query_id, anti$query_id), seq_len(length(q)))
})

test_that("coverage wrappers return expected structure", {
  x <- IRanges::IRanges(start = c(1, 4, 10), end = c(5, 8, 12))

  got_cov <- fast_coverage(x)
  ref_cov <- IRanges::coverage(x)
  expect_equal(got_cov, ref_cov)

  tiles <- fast_tile_coverage(x, tile_width = 4L, step_width = 4L)
  expect_true(all(c("start", "end", "coverage_sum") %in% names(tiles)))
  expect_true(nrow(tiles) > 0L)
})
