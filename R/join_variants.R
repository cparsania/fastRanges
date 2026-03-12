#' Inner Overlap Join
#'
#' Convenience wrapper for `fast_overlap_join(..., join = "inner")`.
#'
#' @inheritParams fast_overlap_join
#' @inheritSection fast_find_overlaps Overlap semantics
#' @return A `data.frame` overlap join result.
#' @export
#'
#' @examples
#' q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' fast_inner_overlap_join(q, s)
fast_inner_overlap_join <- function(
    query,
    subject,
    max_gap = -1L,
    min_overlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore_strand = FALSE,
    threads = fast_default_threads(),
    deterministic = TRUE,
    query_prefix = "query_",
    subject_prefix = "subject_") {
  fast_overlap_join(
    query = query,
    subject = subject,
    join = "inner",
    max_gap = max_gap,
    min_overlap = min_overlap,
    type = match.arg(type),
    ignore_strand = ignore_strand,
    threads = threads,
    deterministic = deterministic,
    query_prefix = query_prefix,
    subject_prefix = subject_prefix
  )
}

#' Left Overlap Join
#'
#' Convenience wrapper for `fast_overlap_join(..., join = "left")`.
#'
#' @inheritParams fast_inner_overlap_join
#' @inheritSection fast_find_overlaps Overlap semantics
#' @return A `data.frame` overlap join result.
#' @export
#'
#' @examples
#' q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' fast_left_overlap_join(q, s)
fast_left_overlap_join <- function(
    query,
    subject,
    max_gap = -1L,
    min_overlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore_strand = FALSE,
    threads = fast_default_threads(),
    deterministic = TRUE,
    query_prefix = "query_",
    subject_prefix = "subject_") {
  fast_overlap_join(
    query = query,
    subject = subject,
    join = "left",
    max_gap = max_gap,
    min_overlap = min_overlap,
    type = match.arg(type),
    ignore_strand = ignore_strand,
    threads = threads,
    deterministic = deterministic,
    query_prefix = query_prefix,
    subject_prefix = subject_prefix
  )
}

#' Semi Overlap Join
#'
#' Return query rows that overlap at least one subject range.
#'
#' @inheritParams fast_find_overlaps
#' @inheritSection fast_find_overlaps Overlap semantics
#' @param query_prefix Prefix applied to query columns.
#'
#' @details
#' This is similar to a SQL `SEMI JOIN`.
#'
#' It keeps only query rows that have at least one overlap hit.
#'
#' `overlap_count` tells you how many subject ranges matched each retained
#' query row.
#'
#' @return A `data.frame` containing matching query rows and overlap counts.
#' @export
#'
#' @examples
#' q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' fast_semi_overlap_join(q, s)
fast_semi_overlap_join <- function(
    query,
    subject,
    max_gap = -1L,
    min_overlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore_strand = FALSE,
    threads = fast_default_threads(),
    deterministic = TRUE,
    query_prefix = "query_") {
  .assert_supported_ranges(query, "query")
  counts <- fast_count_overlaps(
    query = query,
    subject = subject,
    max_gap = max_gap,
    min_overlap = min_overlap,
    type = match.arg(type),
    ignore_strand = ignore_strand,
    threads = threads,
    deterministic = deterministic
  )

  keep <- counts > 0L
  out <- as.data.frame(query)
  names(out) <- paste0(query_prefix, names(out))
  out <- out[keep, , drop = FALSE]
  out$query_id <- which(keep)
  out$overlap_count <- counts[keep]
  out[, c("query_id", "overlap_count", setdiff(names(out), c("query_id", "overlap_count"))), drop = FALSE]
}

#' Anti Overlap Join
#'
#' Return query rows that do not overlap any subject ranges.
#'
#' @inheritParams fast_semi_overlap_join
#' @inheritSection fast_find_overlaps Overlap semantics
#'
#' @details
#' This is similar to a SQL `ANTI JOIN`.
#'
#' It keeps only query rows with zero overlap hits.
#'
#' `overlap_count` is always `0` in the returned table and is included to keep
#' the result grammar parallel to `fast_semi_overlap_join()`.
#'
#' @return A `data.frame` containing query rows with zero overlaps.
#' @export
#'
#' @examples
#' q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' fast_anti_overlap_join(q, s)
fast_anti_overlap_join <- function(
    query,
    subject,
    max_gap = -1L,
    min_overlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore_strand = FALSE,
    threads = fast_default_threads(),
    deterministic = TRUE,
    query_prefix = "query_") {
  .assert_supported_ranges(query, "query")
  counts <- fast_count_overlaps(
    query = query,
    subject = subject,
    max_gap = max_gap,
    min_overlap = min_overlap,
    type = match.arg(type),
    ignore_strand = ignore_strand,
    threads = threads,
    deterministic = deterministic
  )

  keep <- counts == 0L
  out <- as.data.frame(query)
  names(out) <- paste0(query_prefix, names(out))
  out <- out[keep, , drop = FALSE]
  out$query_id <- which(keep)
  out$overlap_count <- counts[keep]
  out[, c("query_id", "overlap_count", setdiff(names(out), c("query_id", "overlap_count"))), drop = FALSE]
}
