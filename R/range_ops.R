#' Reduce Overlapping Ranges
#'
#' Merge overlapping or adjacent ranges.
#'
#' @param x An `IRanges` or `GRanges` object.
#' @param ignore_strand Logical scalar. Ignored for non-genomic ranges.
#' @param min_gap_width Integer scalar minimum gap width for separating ranges.
#' @param with_revmap Logical scalar. If `TRUE`, include reverse mapping.
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
#' @param ignore_strand Logical scalar. Ignored for non-genomic ranges.
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
    return(GenomicRanges::gaps(x, start = start, end = end, ignore.strand = ignore_strand))
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
