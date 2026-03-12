#' Reduce Overlapping Ranges
#'
#' Merge overlapping or adjacent ranges.
#'
#' @param x An `IRanges` or `GRanges` object.
#' @param ignore_strand Logical scalar. Ignored for non-genomic ranges.
#' @param min_gap_width Integer scalar controlling when nearby ranges should be
#'   merged.
#'
#'   `1` merges overlapping or directly adjacent ranges.
#'
#'   Larger values require larger gaps before two ranges are kept separate.
#' @param with_revmap Logical scalar. If `TRUE`, include reverse mapping from
#'   each output range back to contributing input ranges.
#'
#' @details
#' `fast_reduce()` simplifies a range set by merging runs of overlapping or
#' near-adjacent intervals.
#'
#' @return An object of the same range class as `x`.
#' @export
#'
#' @examples
#' x <- IRanges::IRanges(start = c(1, 4, 10), end = c(5, 8, 12))
#' fast_reduce(x)
fast_reduce <- function(
    x,
    ignore_strand = FALSE,
    min_gap_width = 1L,
    with_revmap = FALSE) {
  .assert_supported_ranges(x, "x")
  .assert_scalar_logical(ignore_strand, "ignore_strand")
  .assert_scalar_integerish(min_gap_width, "min_gap_width")
  .assert_scalar_logical(with_revmap, "with_revmap")

  if (.is_granges(x)) {
    return(GenomicRanges::reduce(
      x,
      ignore.strand = ignore_strand,
      min.gapwidth = as.integer(min_gap_width),
      with.revmap = with_revmap
    ))
  }

  IRanges::reduce(
    x,
    min.gapwidth = as.integer(min_gap_width),
    with.revmap = with_revmap
  )
}

#' Disjoin Ranges
#'
#' Return non-overlapping segments induced by input ranges.
#'
#' @inheritParams fast_reduce
#'
#' @details
#' `fast_disjoin()` cuts the covered span of `x` into the smallest
#' non-overlapping pieces.
#'
#' @return An object of the same range class as `x`.
#' @export
#'
#' @examples
#' x <- IRanges::IRanges(start = c(1, 4, 10), end = c(5, 8, 12))
#' fast_disjoin(x)
fast_disjoin <- function(
    x,
    ignore_strand = FALSE,
    with_revmap = FALSE) {
  .assert_supported_ranges(x, "x")
  .assert_scalar_logical(ignore_strand, "ignore_strand")
  .assert_scalar_logical(with_revmap, "with_revmap")

  if (.is_granges(x)) {
    return(GenomicRanges::disjoin(x, ignore.strand = ignore_strand, with.revmap = with_revmap))
  }

  IRanges::disjoin(x, with.revmap = with_revmap)
}

#' Gaps Between Ranges
#'
#' Compute uncovered regions between ranges.
#'
#' @param x An `IRanges` or `GRanges` object.
#' @param start,end Optional integer bounds for non-genomic ranges.
#'   These are most useful for `IRanges`. For `GRanges`, sequence lengths are
#'   usually taken from `seqinfo(x)`.
#' @param ignore_strand Logical scalar. Ignored for non-genomic ranges.
#'
#' @details
#' `fast_gaps()` returns the regions not covered by `x` inside the requested
#' bounds.
#'
#' @return An object of the same range class as `x`.
#' @export
#'
#' @examples
#' x <- IRanges::IRanges(start = c(1, 4, 10), end = c(5, 8, 12))
#' fast_gaps(x, start = 1L, end = 15L)
fast_gaps <- function(
    x,
    start = NULL,
    end = NULL,
    ignore_strand = FALSE) {
  .assert_supported_ranges(x, "x")
  .assert_scalar_logical(ignore_strand, "ignore_strand")

  if (.is_granges(x)) {
    args <- list(x, ignore.strand = ignore_strand)
    if (!is.null(start)) {
      args$start <- start
    }
    if (!is.null(end)) {
      args$end <- end
    }
    return(do.call(GenomicRanges::gaps, args))
  }

  IRanges::gaps(x, start = start, end = end)
}

#' Union of Two Range Sets
#'
#' Compute range-wise union.
#'
#' @param x,y `IRanges` or `GRanges` objects of compatible class.
#' @param ignore_strand Logical scalar. Ignored for non-genomic ranges.
#'
#' @details
#' `fast_range_union()` returns the combined interval coverage of `x` and `y`.
#'
#' @return An object of the same range class as `x` and `y`.
#' @export
#'
#' @examples
#' x <- IRanges::IRanges(start = c(1, 10), end = c(5, 12))
#' y <- IRanges::IRanges(start = c(3, 20), end = c(8, 21))
#' fast_range_union(x, y)
fast_range_union <- function(x, y, ignore_strand = FALSE) {
  .assert_supported_ranges(x, "x")
  .assert_supported_ranges(y, "y")
  .assert_scalar_logical(ignore_strand, "ignore_strand")

  if (.is_granges(x) || .is_granges(y)) {
    return(GenomicRanges::union(x, y, ignore.strand = ignore_strand))
  }
  IRanges::union(x, y)
}

#' Intersection of Two Range Sets
#'
#' Compute range-wise intersection.
#'
#' @inheritParams fast_range_union
#'
#' @details
#' `fast_range_intersect()` keeps only the coordinate span shared by `x` and
#' `y`.
#'
#' @return An object of the same range class as `x` and `y`.
#' @export
#'
#' @examples
#' x <- IRanges::IRanges(start = c(1, 10), end = c(5, 12))
#' y <- IRanges::IRanges(start = c(3, 20), end = c(8, 21))
#' fast_range_intersect(x, y)
fast_range_intersect <- function(x, y, ignore_strand = FALSE) {
  .assert_supported_ranges(x, "x")
  .assert_supported_ranges(y, "y")
  .assert_scalar_logical(ignore_strand, "ignore_strand")

  if (.is_granges(x) || .is_granges(y)) {
    return(GenomicRanges::intersect(x, y, ignore.strand = ignore_strand))
  }
  IRanges::intersect(x, y)
}

#' Set Difference of Two Range Sets
#'
#' Compute ranges in `x` that are not covered by `y`.
#'
#' @inheritParams fast_range_union
#'
#' @details
#' `fast_range_setdiff()` subtracts the covered span of `y` from `x`.
#'
#' @return An object of the same range class as `x`.
#' @export
#'
#' @examples
#' x <- IRanges::IRanges(start = c(1, 10), end = c(5, 12))
#' y <- IRanges::IRanges(start = c(3, 20), end = c(8, 21))
#' fast_range_setdiff(x, y)
fast_range_setdiff <- function(x, y, ignore_strand = FALSE) {
  .assert_supported_ranges(x, "x")
  .assert_supported_ranges(y, "y")
  .assert_scalar_logical(ignore_strand, "ignore_strand")

  if (.is_granges(x) || .is_granges(y)) {
    return(GenomicRanges::setdiff(x, y, ignore.strand = ignore_strand))
  }
  IRanges::setdiff(x, y)
}
