#' Nearest Subject Range per Query
#'
#' Compute nearest-neighbor mapping from query ranges to subject ranges.
#'
#' @param query An `IRanges` or `GRanges` query object.
#' @param subject An `IRanges` or `GRanges` subject object.
#' @param ignore_strand Logical scalar. Ignored for non-genomic ranges.
#' @param threads Integer scalar thread count.
#'
#' @return A `S4Vectors::DataFrame` with `query_id`, `subject_id`, and
#'   `distance` columns.
#' @export
#'
#' @examples
#' q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' fast_nearest(q, s)
fast_nearest <- function(
    query,
    subject,
    ignore_strand = FALSE,
    threads = fast_default_threads()) {
  fast_distance_to_nearest(
    query = query,
    subject = subject,
    ignore_strand = ignore_strand,
    threads = threads
  )
}

#' Distance to Nearest Subject Range
#'
#' Compute nearest-neighbor mapping from query ranges to subject ranges.
#'
#' @inheritParams fast_nearest
#' @return A `S4Vectors::DataFrame` with `query_id`, `subject_id`, and
#'   `distance` columns.
#' @export
#'
#' @examples
#' q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' fast_distance_to_nearest(q, s)
fast_distance_to_nearest <- function(
    query,
    subject,
    ignore_strand = FALSE,
    threads = fast_default_threads()) {
  .assert_supported_ranges(query, "query")
  .assert_supported_ranges(subject, "subject")
  .assert_scalar_logical(ignore_strand, "ignore_strand")
  .normalize_threads(threads)

  nh <- if (.is_granges(query) || .is_granges(subject)) {
    GenomicRanges::distanceToNearest(query, subject, ignore.strand = ignore_strand)
  } else {
    IRanges::distanceToNearest(query, subject)
  }

  S4Vectors::DataFrame(
    query_id = S4Vectors::queryHits(nh),
    subject_id = S4Vectors::subjectHits(nh),
    distance = S4Vectors::mcols(nh)$distance
  )
}

#' Precede Query Ranges
#'
#' Return the index of the first subject range that is strictly after each
#' query range.
#'
#' @inheritParams fast_nearest
#' @return Integer vector of subject indices, with `NA` for unmatched queries.
#' @export
#'
#' @examples
#' q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' fast_precede(q, s)
fast_precede <- function(
    query,
    subject,
    ignore_strand = FALSE,
    threads = fast_default_threads()) {
  .assert_supported_ranges(query, "query")
  .assert_supported_ranges(subject, "subject")
  .assert_scalar_logical(ignore_strand, "ignore_strand")
  .normalize_threads(threads)

  if (.is_granges(query) || .is_granges(subject)) {
    return(GenomicRanges::precede(query, subject, ignore.strand = ignore_strand))
  }
  IRanges::precede(query, subject)
}

#' Follow Query Ranges
#'
#' Return the index of the first subject range that is strictly before each
#' query range.
#'
#' @inheritParams fast_nearest
#' @return Integer vector of subject indices, with `NA` for unmatched queries.
#' @export
#'
#' @examples
#' q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' fast_follow(q, s)
fast_follow <- function(
    query,
    subject,
    ignore_strand = FALSE,
    threads = fast_default_threads()) {
  .assert_supported_ranges(query, "query")
  .assert_supported_ranges(subject, "subject")
  .assert_scalar_logical(ignore_strand, "ignore_strand")
  .normalize_threads(threads)

  if (.is_granges(query) || .is_granges(subject)) {
    return(GenomicRanges::follow(query, subject, ignore.strand = ignore_strand))
  }
  IRanges::follow(query, subject)
}
