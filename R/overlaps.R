#' @useDynLib fastRanges, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

#' Find Overlaps with Deterministic Multithreading
#'
#' Compute overlap pairs between `query` and `subject` using a multithreaded
#' C++ backend. The result is a `Hits` object compatible with Bioconductor
#' workflows.
#'
#' @param query An `IRanges` or `GRanges` query object.
#' @param subject An `IRanges`/`GRanges` object or a `fast_ranges_index`.
#' @param max_gap Integer scalar maximum allowed gap. A value of `-1`
#'   requires direct overlap.
#' @param min_overlap Integer scalar minimum overlap width.
#' @param type Overlap mode: one of `"any"`, `"start"`, `"end"`,
#'   `"within"`, or `"equal"`.
#' @param ignore_strand Logical scalar. Ignored for non-genomic ranges.
#' @param threads Integer scalar thread count.
#' @param deterministic Logical scalar. If `TRUE`, output is sorted by
#'   `(query_id, subject_id)`.
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
  .assert_supported_ranges(query, "query")
  type <- match.arg(type)
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

  hit_idx <- cpp_find_overlaps_indexed(
    q_start = q$start,
    q_end = q$end,
    q_seq = q$seq_id,
    q_strand = q$strand,
    s_start = idx$subject_start,
    s_end = idx$subject_end,
    s_seq = idx$subject_seq,
    s_strand = idx$subject_strand,
    s_original = idx$subject_original_index,
    partition_keys = idx$partition_keys,
    partition_starts = idx$partition_starts,
    partition_ends = idx$partition_ends,
    max_gap = as.integer(max_gap),
    min_overlap = as.integer(min_overlap),
    type = type,
    ignore_strand = isTRUE(ignore_strand),
    threads = as.integer(threads),
    deterministic = isTRUE(deterministic)
  )

  S4Vectors::Hits(
    from = hit_idx$query_hits,
    to = hit_idx$subject_hits,
    nLnode = length(query),
    nRnode = subject_n,
    sort.by.query = FALSE
  )
}

#' Count Overlaps
#'
#' Count subject overlaps per query range.
#'
#' @inheritParams fast_find_overlaps
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
  hits <- fast_find_overlaps(
    query = query,
    subject = subject,
    max_gap = max_gap,
    min_overlap = min_overlap,
    type = match.arg(type),
    ignore_strand = ignore_strand,
    threads = threads,
    deterministic = deterministic
  )

  tabulate(S4Vectors::queryHits(hits), nbins = length(query))
}

#' Overlap Existence per Query
#'
#' Return `TRUE` for queries that overlap at least one subject range.
#'
#' @inheritParams fast_find_overlaps
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
    type = match.arg(type),
    ignore_strand = ignore_strand,
    threads = threads,
    deterministic = deterministic
  ) > 0L
}
