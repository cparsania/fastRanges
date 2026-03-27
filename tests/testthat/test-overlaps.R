library(GenomicRanges)

canon_hits <- function(h) {
  x <- cbind(S4Vectors::queryHits(h), S4Vectors::subjectHits(h))
  if (nrow(x) == 0L) return(x)
  x[order(x[, 1], x[, 2]), , drop = FALSE]
}

validate_select_result <- function(selected, ref_all_hits, query_n, mode = c("first", "last", "arbitrary")) {
  mode <- match.arg(mode)
  expect_type(selected, "integer")
  expect_length(selected, query_n)

  qh <- S4Vectors::queryHits(ref_all_hits)
  sh <- S4Vectors::subjectHits(ref_all_hits)
  ref_choices <- split(sh, qh)
  has_hit <- rep.int(FALSE, query_n)
  if (length(ref_choices) > 0L) {
    has_hit[as.integer(names(ref_choices))] <- TRUE
  }

  expect_identical(is.na(selected), !has_hit)

  for (qid in which(has_hit)) {
    choices <- ref_choices[[as.character(qid)]]
    got <- selected[[qid]]
    expect_false(is.na(got))
    if (identical(mode, "first")) {
      expect_identical(as.integer(got), as.integer(min(choices)))
    } else if (identical(mode, "last")) {
      expect_identical(as.integer(got), as.integer(max(choices)))
    } else {
      expect_true(got %in% choices)
    }
  }
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

test_that("select semantics match GenomicRanges for GRanges", {
  data(fast_ranges_example, package = "fastRanges")
  q <- fast_ranges_example$query
  s <- fast_ranges_example$subject
  ref_all <- GenomicRanges::findOverlaps(q, s, ignore.strand = FALSE)

  for (sel in c("first", "last")) {
    ref <- GenomicRanges::findOverlaps(q, s, select = sel, ignore.strand = FALSE)
    got <- fast_find_overlaps(q, s, select = sel, threads = 2)
    expect_identical(got, ref)
  }

  got_arb <- fast_find_overlaps(q, s, select = "arbitrary", threads = 2)
  validate_select_result(got_arb, ref_all, length(q), mode = "arbitrary")
})

test_that("select semantics match IRanges for IRanges inputs", {
  q <- IRanges::IRanges(start = c(1, 10, 20, 40), width = c(5, 4, 3, 2))
  s <- IRanges::IRanges(start = c(3, 9, 18, 41), width = c(4, 6, 5, 2))
  ref_all <- IRanges::findOverlaps(q, s)

  for (sel in c("first", "last")) {
    ref <- IRanges::findOverlaps(q, s, select = sel)
    got <- fast_find_overlaps(q, s, select = sel, threads = 2)
    expect_identical(got, ref)
  }

  got_arb <- fast_find_overlaps(q, s, select = "arbitrary", threads = 2)
  validate_select_result(got_arb, ref_all, length(q), mode = "arbitrary")
})

test_that("select returns NA for queries with no hits", {
  q <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = c(1L, 100L), width = c(5L, 5L)))
  s <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = 3L, width = 2L))

  got <- fast_find_overlaps(q, s, select = "first", threads = 2)

  expect_identical(got, c(1L, NA_integer_))
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

test_that("width-0 GRanges semantics match GenomicRanges for direct and indexed paths", {
  q <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = c(5L, 8L, 12L, 20L), end = c(4L, 10L, 11L, 23L))
  )
  s <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = c(3L, 5L, 8L, 12L, 19L), end = c(6L, 4L, 10L, 11L, 21L))
  )
  idx <- fast_build_index(s)

  ref_all <- GenomicRanges::findOverlaps(q, s)
  got_direct <- fast_find_overlaps(q, s, threads = 2)
  got_index <- fast_find_overlaps(q, idx, threads = 2)
  expect_equal(canon_hits(got_direct), canon_hits(ref_all))
  expect_equal(canon_hits(got_index), canon_hits(ref_all))

  ref_first <- GenomicRanges::findOverlaps(q, s, select = "first")
  ref_last <- GenomicRanges::findOverlaps(q, s, select = "last")
  expect_identical(fast_find_overlaps(q, s, select = "first", threads = 2), ref_first)
  expect_identical(fast_find_overlaps(q, idx, select = "first", threads = 2), ref_first)
  expect_identical(fast_find_overlaps(q, s, select = "last", threads = 2), ref_last)
  expect_identical(fast_find_overlaps(q, idx, select = "last", threads = 2), ref_last)

  got_arb_direct <- fast_find_overlaps(q, s, select = "arbitrary", threads = 2)
  got_arb_index <- fast_find_overlaps(q, idx, select = "arbitrary", threads = 2)
  validate_select_result(got_arb_direct, ref_all, length(q), mode = "arbitrary")
  validate_select_result(got_arb_index, ref_all, length(q), mode = "arbitrary")

  ref_count <- GenomicRanges::countOverlaps(q, s)
  expect_identical(as.integer(fast_count_overlaps(q, s, threads = 2)), as.integer(ref_count))
  expect_identical(as.integer(fast_count_overlaps(q, idx, threads = 2)), as.integer(ref_count))
  expect_identical(as.logical(fast_overlaps_any(q, s, threads = 2)), ref_count > 0L)
  expect_identical(as.logical(fast_overlaps_any(q, idx, threads = 2)), ref_count > 0L)
})

test_that("width-0 IRanges semantics match IRanges for direct path", {
  q <- IRanges::IRanges(start = c(5L, 8L, 12L, 20L), end = c(4L, 10L, 11L, 23L))
  s <- IRanges::IRanges(start = c(3L, 5L, 8L, 12L, 19L), end = c(6L, 4L, 10L, 11L, 21L))
  idx <- fast_build_index(s)

  ref_all <- IRanges::findOverlaps(q, s)
  got_direct <- fast_find_overlaps(q, s, threads = 2)
  got_index <- fast_find_overlaps(q, idx, threads = 2)
  expect_equal(canon_hits(got_direct), canon_hits(ref_all))
  expect_equal(canon_hits(got_index), canon_hits(ref_all))

  ref_first <- IRanges::findOverlaps(q, s, select = "first")
  ref_last <- IRanges::findOverlaps(q, s, select = "last")
  expect_identical(fast_find_overlaps(q, s, select = "first", threads = 2), ref_first)
  expect_identical(fast_find_overlaps(q, idx, select = "first", threads = 2), ref_first)
  expect_identical(fast_find_overlaps(q, s, select = "last", threads = 2), ref_last)
  expect_identical(fast_find_overlaps(q, idx, select = "last", threads = 2), ref_last)

  got_arb_direct <- fast_find_overlaps(q, s, select = "arbitrary", threads = 2)
  got_arb_index <- fast_find_overlaps(q, idx, select = "arbitrary", threads = 2)
  validate_select_result(got_arb_direct, ref_all, length(q), mode = "arbitrary")
  validate_select_result(got_arb_index, ref_all, length(q), mode = "arbitrary")

  ref_count <- IRanges::countOverlaps(q, s)
  expect_identical(as.integer(fast_count_overlaps(q, s, threads = 2)), as.integer(ref_count))
  expect_identical(as.integer(fast_count_overlaps(q, idx, threads = 2)), as.integer(ref_count))
  expect_identical(as.logical(fast_overlaps_any(q, s, threads = 2)), ref_count > 0L)
  expect_identical(as.logical(fast_overlaps_any(q, idx, threads = 2)), ref_count > 0L)
})

test_that("circular GRanges are rejected explicitly", {
  si <- GenomeInfoDb::Seqinfo("chrC", seqlengths = 100L, isCircular = TRUE)
  q <- suppressWarnings(GenomicRanges::GRanges("chrC", IRanges::IRanges(start = 95L, width = 10L), seqinfo = si))
  s <- suppressWarnings(GenomicRanges::GRanges("chrC", IRanges::IRanges(start = 2L, width = 8L), seqinfo = si))

  expect_error(
    fast_find_overlaps(q, s, threads = 2),
    "circular sequences, which are not currently supported"
  )
  expect_error(
    fast_count_overlaps(q, s, threads = 2),
    "circular sequences, which are not currently supported"
  )
  expect_error(
    fast_build_index(s),
    "circular sequences, which are not currently supported"
  )
})

test_that("indexed subject with circular flag is rejected explicitly", {
  s <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = c(1L, 10L), width = 5L))
  q <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = c(2L, 12L), width = 3L))
  idx <- fast_build_index(s)
  idx$has_circular_sequences <- TRUE

  expect_error(
    fast_find_overlaps(q, idx, threads = 2),
    "circular sequences, which are not currently supported"
  )
  expect_error(
    fast_count_overlaps(q, idx, threads = 2),
    "circular sequences, which are not currently supported"
  )
})

test_that("GRangesList is rejected explicitly", {
  q <- GenomicRanges::GRangesList(
    a = GenomicRanges::GRanges("chr1", IRanges::IRanges(start = c(1L, 10L), width = c(3L, 4L))),
    b = GenomicRanges::GRanges("chr1", IRanges::IRanges(start = 20L, width = 5L))
  )
  s <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = c(2L, 12L, 19L), width = c(2L, 3L, 4L)))

  expect_error(
    fast_find_overlaps(q, s, threads = 2),
    "GRangesList is not currently supported"
  )
  expect_error(
    fast_count_overlaps(q, s, threads = 2),
    "GRangesList is not currently supported"
  )
  expect_error(
    fast_build_index(q),
    "GRangesList is not currently supported"
  )
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
