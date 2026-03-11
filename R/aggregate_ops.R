#' Count Overlaps by Subject Group
#'
#' Count overlaps per query and per subject metadata group.
#'
#' @inheritParams fast_find_overlaps
#' @param group_col Subject metadata column name used for grouping.
#' @param include_na_group Logical scalar. If `TRUE`, missing group values are
#'   counted as `"<NA>"`.
#'
#' @return Integer matrix with one row per query and one column per group.
#' @export
#'
#' @examples
#' q <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 10), width = 5))
#' s <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(2, 9), width = 5))
#' S4Vectors::mcols(s)$grp <- c("A", "B")
#' fast_count_overlaps_by_group(q, s, group_col = "grp")
fast_count_overlaps_by_group <- function(
    query,
    subject,
    group_col,
    max_gap = -1L,
    min_overlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore_strand = FALSE,
    threads = fast_default_threads(),
    deterministic = TRUE,
    include_na_group = FALSE) {
  .assert_supported_ranges(query, "query")
  .assert_supported_ranges(subject, "subject")
  .assert_scalar_logical(include_na_group, "include_na_group")

  if (!is.character(group_col) || length(group_col) != 1L || is.na(group_col) || nchar(group_col) == 0L) {
    stop("`group_col` must be a non-empty character scalar", call. = FALSE)
  }

  subj_mcols <- S4Vectors::mcols(subject)
  if (!(group_col %in% names(subj_mcols))) {
    stop(sprintf("`group_col` '%s' is not present in `mcols(subject)`", group_col), call. = FALSE)
  }

  groups_raw <- subj_mcols[[group_col]]
  groups <- as.character(groups_raw)
  if (include_na_group) {
    groups[is.na(groups)] <- "<NA>"
  }

  valid <- !is.na(groups)
  group_levels <- sort(unique(groups[valid]))
  if (length(group_levels) == 0L) {
    return(matrix(0L, nrow = length(query), ncol = 0L))
  }

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

  out <- matrix(0L, nrow = length(query), ncol = length(group_levels))
  colnames(out) <- group_levels

  if (length(hits) == 0L) {
    return(out)
  }

  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)
  hit_groups <- groups[sh]
  keep <- !is.na(hit_groups)

  if (!any(keep)) {
    return(out)
  }

  qh <- qh[keep]
  hit_groups <- hit_groups[keep]
  tab <- table(qh, factor(hit_groups, levels = group_levels))
  out[as.integer(rownames(tab)), ] <- out[as.integer(rownames(tab)), , drop = FALSE] + tab
  out
}

#' Aggregate Subject Metadata Over Overlaps
#'
#' Aggregate a numeric subject metadata column across overlaps for each query.
#'
#' @inheritParams fast_find_overlaps
#' @param value_col Subject metadata column name containing numeric values.
#' @param fun Aggregation function: one of `"count"`, `"sum"`, `"mean"`,
#'   `"min"`, or `"max"`.
#' @param na_rm Logical scalar. If `TRUE`, remove missing values in
#'   aggregation.
#'
#' @return Numeric vector with one value per query.
#' @export
#'
#' @examples
#' q <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 10), width = 5))
#' s <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(2, 9), width = 5))
#' S4Vectors::mcols(s)$score <- c(2, 5)
#' fast_overlap_aggregate(q, s, value_col = "score", fun = "sum")
fast_overlap_aggregate <- function(
    query,
    subject,
    value_col = NULL,
    fun = c("count", "sum", "mean", "min", "max"),
    max_gap = -1L,
    min_overlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore_strand = FALSE,
    threads = fast_default_threads(),
    deterministic = TRUE,
    na_rm = TRUE) {
  .assert_supported_ranges(query, "query")
  .assert_supported_ranges(subject, "subject")
  .assert_scalar_logical(na_rm, "na_rm")
  fun <- match.arg(fun)

  if (fun != "count") {
    if (!is.character(value_col) || length(value_col) != 1L || is.na(value_col) || nchar(value_col) == 0L) {
      stop("`value_col` must be a non-empty character scalar when `fun` is not 'count'", call. = FALSE)
    }

    subj_mcols <- S4Vectors::mcols(subject)
    if (!(value_col %in% names(subj_mcols))) {
      stop(sprintf("`value_col` '%s' is not present in `mcols(subject)`", value_col), call. = FALSE)
    }
    values <- as.numeric(subj_mcols[[value_col]])
  } else {
    values <- NULL
  }

  if (fun == "count") {
    return(as.numeric(fast_count_overlaps(
      query = query,
      subject = subject,
      max_gap = max_gap,
      min_overlap = min_overlap,
      type = match.arg(type),
      ignore_strand = ignore_strand,
      threads = threads,
      deterministic = deterministic
    )))
  }

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

  out <- rep(NA_real_, length(query))
  if (length(hits) == 0L) {
    return(out)
  }

  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)
  x <- values[sh]

  if (na_rm) {
    keep <- !is.na(x)
    qh <- qh[keep]
    x <- x[keep]
  }

  if (length(x) == 0L) {
    return(out)
  }

  split_x <- split(x, qh)
  agg_fun <- switch(
    fun,
    sum = function(z) sum(z, na.rm = na_rm),
    mean = function(z) mean(z, na.rm = na_rm),
    min = function(z) min(z, na.rm = na_rm),
    max = function(z) max(z, na.rm = na_rm)
  )

  aggregated <- vapply(split_x, agg_fun, numeric(1))
  out[as.integer(names(aggregated))] <- aggregated
  out
}
