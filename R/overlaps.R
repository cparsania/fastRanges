#' @useDynLib fastRanges, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

#' @keywords internal
.prepare_overlap_inputs <- function(
    query,
    subject,
    max_gap,
    min_overlap,
    type,
    threads) {
  .assert_supported_ranges(query, "query")
  type <- match.arg(type, c("any", "start", "end", "within", "equal"))
  .assert_scalar_integerish(max_gap, "max_gap")
  .assert_scalar_integerish(min_overlap, "min_overlap")
  threads <- .normalize_threads(threads)

  if (inherits(subject, "fast_ranges_index")) {
    idx <- subject
    subject_n <- idx$subject_n
    seq_map <- idx$seq_map
  } else {
    .assert_supported_ranges(subject, "subject")
    idx <- fast_build_index(subject)
    subject_n <- length(subject)
    seq_map <- idx$seq_map
  }

  q <- .encode_query(query, seq_map)

  list(
    idx = idx,
    q = q,
    subject_n = subject_n,
    type = type,
    max_gap = as.integer(max_gap),
    min_overlap = as.integer(min_overlap),
    threads = as.integer(threads)
  )
}

#' Find Overlaps with Deterministic Multithreading
#'
#' Compute overlap pairs between `query` and `subject` using a multithreaded
#' C++ backend. The result is a `Hits` object compatible with Bioconductor
#' workflows.
#'
#' @param query An `IRanges` or `GRanges` query object.
#' @param subject An `IRanges`/`GRanges` object or a `fast_ranges_index`.
#'   Use `fast_build_index(subject)` when the same subject is reused across
#'   many overlap queries.
#' @template overlap_shared_args
#' @template overlap_type_details
#'
#' @details
#' This is the core matching function in `fastRanges`.
#'
#' Think of it as answering the question: "for each query range, which subject
#' ranges satisfy my overlap rule?"
#'
#' The return value is a `Hits` object. The important pieces are:
#'
#' `queryHits(h)` gives the row numbers from `query`.
#'
#' `subjectHits(h)` gives the matching row numbers from `subject`.
#'
#' If you need only counts or a yes/no answer, prefer `fast_count_overlaps()`
#' or `fast_overlaps_any()`, because they return simpler summaries.
#'
#' @return A `S4Vectors::Hits` object.
#' @export
#'
#' @examples
#' q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' h <- fast_find_overlaps(q, s, threads = 1)
#' length(h)
fast_find_overlaps <- function(
    query,
    subject,
    max_gap = -1L,
    min_overlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore_strand = FALSE,
    threads = fast_default_threads(),
    deterministic = TRUE) {
  inputs <- .prepare_overlap_inputs(
    query = query,
    subject = subject,
    max_gap = max_gap,
    min_overlap = min_overlap,
    type = type,
    threads = threads
  )

  hit_idx <- cpp_find_overlaps_indexed(
    q_start = inputs$q$start,
    q_end = inputs$q$end,
    q_seq = inputs$q$seq_id,
    q_strand = inputs$q$strand,
    s_start = inputs$idx$subject_start,
    s_end = inputs$idx$subject_end,
    s_seq = inputs$idx$subject_seq,
    s_strand = inputs$idx$subject_strand,
    s_original = inputs$idx$subject_original_index,
    partition_keys = inputs$idx$partition_keys,
    partition_starts = inputs$idx$partition_starts,
    partition_ends = inputs$idx$partition_ends,
    max_gap = inputs$max_gap,
    min_overlap = inputs$min_overlap,
    type = inputs$type,
    ignore_strand = isTRUE(ignore_strand),
    threads = inputs$threads,
    deterministic = isTRUE(deterministic)
  )

  S4Vectors::Hits(
    from = hit_idx$query_hits,
    to = hit_idx$subject_hits,
    nLnode = length(query),
    nRnode = inputs$subject_n,
    sort.by.query = FALSE
  )
}

#' Count Overlaps
#'
#' Count subject overlaps per query range.
#'
#' `deterministic` does not change returned counts for this summary output.
#'
#' @inheritParams fast_find_overlaps
#' @inheritSection fast_find_overlaps Overlap semantics
#'
#' @details
#' Returns one integer per query range.
#'
#' A value of `0` means that query had no matching subjects.
#'
#' A value of `5` means that query matched five subject ranges under the chosen
#' overlap rule.
#'
#' @return Integer vector with one element per query range.
#' @export
#'
#' @examples
#' q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' fast_count_overlaps(q, s, threads = 1)
fast_count_overlaps <- function(
    query,
    subject,
    max_gap = -1L,
    min_overlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore_strand = FALSE,
    threads = fast_default_threads(),
    deterministic = TRUE) {
  .assert_scalar_logical(deterministic, "deterministic")
  inputs <- .prepare_overlap_inputs(
    query = query,
    subject = subject,
    max_gap = max_gap,
    min_overlap = min_overlap,
    type = type,
    threads = threads
  )

  cpp_count_overlaps_indexed(
    q_start = inputs$q$start,
    q_end = inputs$q$end,
    q_seq = inputs$q$seq_id,
    q_strand = inputs$q$strand,
    s_start = inputs$idx$subject_start,
    s_end = inputs$idx$subject_end,
    s_seq = inputs$idx$subject_seq,
    s_strand = inputs$idx$subject_strand,
    partition_keys = inputs$idx$partition_keys,
    partition_starts = inputs$idx$partition_starts,
    partition_ends = inputs$idx$partition_ends,
    max_gap = inputs$max_gap,
    min_overlap = inputs$min_overlap,
    type = inputs$type,
    ignore_strand = isTRUE(ignore_strand),
    threads = inputs$threads
  )
}

#' Overlap Existence per Query
#'
#' Return `TRUE` for queries that overlap at least one subject range.
#'
#' `deterministic` does not change returned logical values for this summary
#' output.
#'
#' @inheritParams fast_find_overlaps
#' @inheritSection fast_find_overlaps Overlap semantics
#'
#' @details
#' Returns one logical value per query range.
#'
#' `TRUE` means at least one subject range matched.
#'
#' `FALSE` means no subject range matched.
#'
#' @return Logical vector with one element per query range.
#' @export
#'
#' @examples
#' q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' fast_overlaps_any(q, s, threads = 1)
fast_overlaps_any <- function(
    query,
    subject,
    max_gap = -1L,
    min_overlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore_strand = FALSE,
    threads = fast_default_threads(),
    deterministic = TRUE) {
  fast_count_overlaps(
    query = query,
    subject = subject,
    max_gap = max_gap,
    min_overlap = min_overlap,
    type = type,
    ignore_strand = ignore_strand,
    threads = threads,
    deterministic = deterministic
  ) > 0L
}
