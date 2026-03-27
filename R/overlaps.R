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
    idx <- ._add_block_index_fields(subject)
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
#' @param select Character scalar controlling whether all hits are returned or
#'   a single subject match is selected per query.
#'
#'   Use `"all"` to return a `Hits` object.
#'
#'   Use `"first"`, `"last"`, or `"arbitrary"` to return an integer vector
#'   with one selected subject index per query and `NA` for queries with no hit.
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
#' Compatibility notes:
#'
#' `fastRanges` aims to stay close to Bioconductor overlap semantics for
#' supported inputs, and its outputs are routinely validated against
#' `IRanges::findOverlaps()` / `GenomicRanges::findOverlaps()`.
#'
#' Currently supported core input types are `IRanges` and `GRanges`.
#'
#' Empty-range semantics are delegated to Bioconductor-compatible reference
#' behavior.
#'
#' Circular genomic sequences are not currently supported and will raise an
#' explicit error.
#'
#' `GRangesList` inputs are not currently supported and will raise an explicit
#' error.
#'
#' Performance notes:
#'
#' For one-off overlap calls, use the raw `subject`.
#'
#' For repeated-query or throughput-oriented workloads, build a reusable index
#' once with `fast_build_index(subject)` and pass that index as `subject`.
#'
#' For maximum multithreaded throughput, consider `deterministic = FALSE` when
#' output order is not important.
#'
#' @return
#' If `select = "all"`, a `S4Vectors::Hits` object.
#'
#' Otherwise, an integer vector of length `length(query)` containing one
#' selected subject index per query and `NA` when no subject matched.
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
    select = c("all", "first", "last", "arbitrary"),
    max_gap = -1L,
    min_overlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore_strand = FALSE,
    threads = fast_default_threads(),
    deterministic = TRUE) {
  select <- match.arg(select)
  .assert_supported_ranges(query, "query")
  .assert_no_circular_ranges(query, "query")
  if (inherits(subject, "fast_ranges_index")) {
    if (isTRUE(subject$has_circular_sequences)) {
      stop("`subject` contains circular sequences, which are not currently supported", call. = FALSE)
    }
  } else {
    .assert_supported_ranges(subject, "subject")
    .assert_no_circular_ranges(subject, "subject")
  }
  query_has_empty <- .has_empty_ranges(query)
  subject_has_empty <- if (inherits(subject, "fast_ranges_index")) {
    isTRUE(subject$has_empty_ranges)
  } else {
    .has_empty_ranges(subject)
  }
  if (query_has_empty || subject_has_empty) {
    subject_ref <- if (inherits(subject, "fast_ranges_index")) .restore_subject_from_index(subject) else subject
    return(.find_overlaps_reference(
      query = query,
      subject = subject_ref,
      select = select,
      max_gap = max_gap,
      min_overlap = min_overlap,
      type = type,
      ignore_strand = ignore_strand
    ))
  }
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
    block_starts = inputs$idx$block_starts,
    block_ends = inputs$idx$block_ends,
    block_first_start = inputs$idx$block_first_start,
    block_max_end = inputs$idx$block_max_end,
    partition_block_starts = inputs$idx$partition_block_starts,
    partition_block_ends = inputs$idx$partition_block_ends,
    max_gap = inputs$max_gap,
    min_overlap = inputs$min_overlap,
    type = inputs$type,
    ignore_strand = isTRUE(ignore_strand),
    threads = inputs$threads,
    deterministic = isTRUE(deterministic)
  )

  hits <- S4Vectors::Hits(
    from = hit_idx$query_hits,
    to = hit_idx$subject_hits,
    nLnode = length(query),
    nRnode = inputs$subject_n,
    sort.by.query = FALSE
  )

  if (identical(select, "all")) {
    return(hits)
  }

  .select_from_hits(hits, select = select, query_n = length(query))
}

#' @keywords internal
.select_from_hits <- function(hits, select, query_n) {
  select <- match.arg(select, c("first", "last", "arbitrary"))
  out <- rep.int(NA_integer_, query_n)
  if (length(hits) == 0L || query_n == 0L) {
    return(out)
  }

  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)
  ord <- order(qh, sh)
  qh <- qh[ord]
  sh <- sh[ord]

  if (identical(select, "last")) {
    keep <- !duplicated(qh, fromLast = TRUE)
  } else {
    keep <- !duplicated(qh)
  }

  out[qh[keep]] <- sh[keep]
  out
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
  .assert_supported_ranges(query, "query")
  .assert_no_circular_ranges(query, "query")
  if (inherits(subject, "fast_ranges_index")) {
    if (isTRUE(subject$has_circular_sequences)) {
      stop("`subject` contains circular sequences, which are not currently supported", call. = FALSE)
    }
  } else {
    .assert_supported_ranges(subject, "subject")
    .assert_no_circular_ranges(subject, "subject")
  }
  query_has_empty <- .has_empty_ranges(query)
  subject_has_empty <- if (inherits(subject, "fast_ranges_index")) {
    isTRUE(subject$has_empty_ranges)
  } else {
    .has_empty_ranges(subject)
  }
  if (query_has_empty || subject_has_empty) {
    subject_ref <- if (inherits(subject, "fast_ranges_index")) .restore_subject_from_index(subject) else subject
    return(.count_overlaps_reference(
      query = query,
      subject = subject_ref,
      max_gap = max_gap,
      min_overlap = min_overlap,
      type = type,
      ignore_strand = ignore_strand
    ))
  }
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
    block_starts = inputs$idx$block_starts,
    block_ends = inputs$idx$block_ends,
    block_first_start = inputs$idx$block_first_start,
    block_max_end = inputs$idx$block_max_end,
    partition_block_starts = inputs$idx$partition_block_starts,
    partition_block_ends = inputs$idx$partition_block_ends,
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
