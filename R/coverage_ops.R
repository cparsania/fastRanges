#' Coverage Across Ranges
#'
#' Compute per-position coverage for input ranges.
#'
#' @param x An `IRanges` or `GRanges` object.
#' @param shift Passed to `coverage()`.
#' @param width Passed to `coverage()`.
#' @param weight Passed to `coverage()`.
#' @param method Coverage method.
#' @param threads Integer scalar thread count. Reserved for API consistency.
#'
#' @return `Rle` (for `IRanges`) or `RleList` (for `GRanges`).
#' @export
#'
#' @examples
#' x <- IRanges::IRanges(start = c(1, 4, 10), end = c(5, 8, 12))
#' fast_coverage(x)
fast_coverage <- function(
    x,
    shift = 0L,
    width = NULL,
    weight = 1L,
    method = c("auto", "sort", "hash"),
    threads = fast_default_threads()) {
  .assert_supported_ranges(x, "x")
  .normalize_threads(threads)
  method <- match.arg(method)

  if (.is_granges(x)) {
    return(GenomicRanges::coverage(x, shift = shift, width = width, weight = weight, method = method))
  }

  IRanges::coverage(x, shift = shift, width = width, weight = weight, method = method)
}

#' Tile-Based Coverage Summary
#'
#' Aggregate coverage into fixed-width tiles.
#'
#' @param x An `IRanges` or `GRanges` object.
#' @param tile_width Integer scalar tile width.
#' @param step_width Integer scalar step width.
#' @param shift Passed to `coverage()`.
#' @param width Passed to `coverage()`.
#' @param weight Passed to `coverage()`.
#' @param method Coverage method.
#' @param threads Integer scalar thread count. Reserved for API consistency.
#'
#' @return A `data.frame` with tile coordinates and `coverage_sum`.
#' @export
#'
#' @examples
#' x <- IRanges::IRanges(start = c(1, 4, 10), end = c(5, 8, 12))
#' fast_tile_coverage(x, tile_width = 5L)
fast_tile_coverage <- function(
    x,
    tile_width,
    step_width = tile_width,
    shift = 0L,
    width = NULL,
    weight = 1L,
    method = c("auto", "sort", "hash"),
    threads = fast_default_threads()) {
  .assert_supported_ranges(x, "x")
  .assert_positive_integerish(tile_width, "tile_width")
  .assert_positive_integerish(step_width, "step_width")
  .normalize_threads(threads)
  method <- match.arg(method)

  cov <- fast_coverage(
    x = x,
    shift = shift,
    width = width,
    weight = weight,
    method = method,
    threads = threads
  )

  if (.is_granges(x)) {
    cov_names <- names(cov)
    if (is.null(cov_names)) {
      cov_names <- paste0("seq", seq_along(cov))
    }

    out <- vector("list", length(cov))
    for (i in seq_along(cov)) {
      n <- length(cov[[i]])
      starts <- seq.int(1L, n, by = as.integer(step_width))
      ends <- pmin(starts + as.integer(tile_width) - 1L, n)
      sums <- as.numeric(IRanges::viewSums(IRanges::Views(cov[[i]], start = starts, end = ends)))
      out[[i]] <- data.frame(
        seqnames = cov_names[i],
        start = starts,
        end = ends,
        coverage_sum = sums,
        stringsAsFactors = FALSE
      )
    }
    return(do.call(rbind, out))
  }

  n <- length(cov)
  starts <- seq.int(1L, n, by = as.integer(step_width))
  ends <- pmin(starts + as.integer(tile_width) - 1L, n)
  sums <- as.numeric(IRanges::viewSums(IRanges::Views(cov, start = starts, end = ends)))
  data.frame(
    start = starts,
    end = ends,
    coverage_sum = sums,
    stringsAsFactors = FALSE
  )
}
