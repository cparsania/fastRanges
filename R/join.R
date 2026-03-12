#' Join Ranges by Overlap
#'
#' Build a metadata-preserving overlap join with a consistent, tidy tabular
#' output grammar.
#'
#' @inheritParams fast_find_overlaps
#' @inheritSection fast_find_overlaps Overlap semantics
#' @param join Join mode.
#'
#'   `"inner"` returns one row per overlap hit and drops queries with no hit.
#'
#'   `"left"` keeps every query at least once. Queries with no hit get
#'   `NA` values in subject columns.
#' @param query_prefix Prefix added to columns derived from `query`.
#'   This helps you see which output columns came from the query object.
#' @param subject_prefix Prefix added to columns derived from `subject`.
#'   This helps you see which output columns came from the subject object.
#'
#' @details
#' The join family turns overlap hits into beginner-friendly tabular output.
#'
#' Output always starts with `query_id` and `subject_id`.
#'
#' `query_id` is the 1-based row index from `query`.
#'
#' `subject_id` is the 1-based row index from `subject`, or `NA` for unmatched
#' queries in a left join.
#'
#' The remaining columns are prefixed copies of the original range columns and
#' metadata columns.
#'
#' @return A `data.frame` with overlap ids and prefixed query/subject columns.
#' @export
#'
#' @examples
#' q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' fast_overlap_join(q, s, join = "inner", threads = 1)
fast_overlap_join <- function(
    query,
    subject,
    join = c("inner", "left"),
    max_gap = -1L,
    min_overlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore_strand = FALSE,
    threads = fast_default_threads(),
    deterministic = TRUE,
    query_prefix = "query_",
    subject_prefix = "subject_") {
  if (inherits(subject, "fast_ranges_index")) {
    stop("`fast_overlap_join()` requires raw `subject` ranges for metadata output", call. = FALSE)
  }

  join <- match.arg(join)
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

  query_df <- as.data.frame(query)
  subject_df <- as.data.frame(subject)
  names(query_df) <- paste0(query_prefix, names(query_df))
  names(subject_df) <- paste0(subject_prefix, names(subject_df))

  q_idx <- S4Vectors::queryHits(hits)
  s_idx <- S4Vectors::subjectHits(hits)

  if (join == "inner") {
    if (length(q_idx) == 0L) {
      empty <- data.frame(query_id = integer(), subject_id = integer())
      return(cbind(empty, query_df[FALSE, , drop = FALSE], subject_df[FALSE, , drop = FALSE]))
    }

    out <- data.frame(query_id = q_idx, subject_id = s_idx)
    return(cbind(out, query_df[q_idx, , drop = FALSE], subject_df[s_idx, , drop = FALSE]))
  }

  hit_map <- data.frame(query_id = q_idx, subject_id = s_idx)
  all_query <- data.frame(query_id = seq_len(length(query)))
  if (nrow(hit_map) == 0L) {
    joined <- transform(all_query, subject_id = NA_integer_)
  } else {
    joined <- merge(all_query, hit_map, by = "query_id", all.x = TRUE, sort = FALSE)
  }

  cbind(
    joined,
    query_df[joined$query_id, , drop = FALSE],
    subject_df[joined$subject_id, , drop = FALSE]
  )
}
