#' Self Overlaps
#'
#' Find overlaps within a single range object.
#'
#' @param x An `IRanges` or `GRanges` object.
#' @param drop_self Logical scalar. If `TRUE`, self-hits (`i` vs `i`) are
#'   removed.
#' @param drop_redundant Logical scalar. If `TRUE`, redundant pairs are removed
#'   for self-comparisons by keeping only `query_id < subject_id`.
#' @inheritParams fast_find_overlaps
#' @inheritSection fast_find_overlaps Overlap semantics
#'
#' @return A `S4Vectors::Hits` object.
#' @export
#'
#' @examples
#' x <- IRanges::IRanges(start = c(1, 3, 10), end = c(5, 8, 12))
#' fast_self_overlaps(x)
fast_self_overlaps <- function(
    x,
    max_gap = -1L,
    min_overlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore_strand = FALSE,
    threads = fast_default_threads(),
    deterministic = TRUE,
    drop_self = TRUE,
    drop_redundant = TRUE) {
  .assert_supported_ranges(x, "x")
  .assert_scalar_logical(drop_self, "drop_self")
  .assert_scalar_logical(drop_redundant, "drop_redundant")

  # Always compute unsorted hits, then apply self/redundant filters and
  # optional deterministic ordering once. This avoids unnecessary sorting.
  hits <- fast_find_overlaps(
    query = x,
    subject = x,
    max_gap = max_gap,
    min_overlap = min_overlap,
    type = match.arg(type),
    ignore_strand = ignore_strand,
    threads = threads,
    deterministic = FALSE
  )

  if (length(hits) == 0L) {
    return(hits)
  }

  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)

  if (drop_redundant) {
    keep <- qh < sh
    qh <- qh[keep]
    sh <- sh[keep]
  } else if (drop_self) {
    keep <- qh != sh
    qh <- qh[keep]
    sh <- sh[keep]
  }

  if (isTRUE(deterministic) && length(qh) > 1L) {
    ord <- order(qh, sh)
    qh <- qh[ord]
    sh <- sh[ord]
  }

  S4Vectors::Hits(
    from = qh,
    to = sh,
    nLnode = length(x),
    nRnode = length(x),
    sort.by.query = FALSE
  )
}

#' Cluster Overlapping Ranges
#'
#' Assign each range to an overlap-connected component.
#'
#' Two ranges are assigned to the same cluster when they are connected by a
#' chain of overlaps under the provided overlap settings.
#'
#' @inheritParams fast_self_overlaps
#' @inheritSection fast_find_overlaps Overlap semantics
#' @param return One of `"vector"` or `"data.frame"`.
#'
#' @return If `return = "vector"`, an integer vector of cluster IDs with one
#'   element per range in `x`. If `return = "data.frame"`, a `data.frame` with
#'   `range_id`, `cluster_id`, and `cluster_size`.
#' @export
#'
#' @examples
#' x <- IRanges::IRanges(start = c(1, 3, 10, 11, 20), end = c(5, 8, 12, 14, 22))
#' fast_cluster_overlaps(x)
fast_cluster_overlaps <- function(
    x,
    max_gap = -1L,
    min_overlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore_strand = FALSE,
    threads = fast_default_threads(),
    deterministic = TRUE,
    return = c("vector", "data.frame")) {
  .assert_supported_ranges(x, "x")
  return <- match.arg(return)

  n <- length(x)
  if (n == 0L) {
    if (return == "vector") {
      return(integer())
    }
    return(data.frame(
      range_id = integer(),
      cluster_id = integer(),
      cluster_size = integer(),
      stringsAsFactors = FALSE
    ))
  }

  hits <- fast_self_overlaps(
    x = x,
    max_gap = max_gap,
    min_overlap = min_overlap,
    type = match.arg(type),
    ignore_strand = ignore_strand,
    threads = threads,
    deterministic = deterministic,
    drop_self = TRUE,
    drop_redundant = TRUE
  )

  parent <- seq_len(n)
  rank <- integer(n)

  find_root <- function(i) {
    while (parent[i] != i) {
      parent[i] <<- parent[parent[i]]
      i <- parent[i]
    }
    i
  }

  union_sets <- function(a, b) {
    root_a <- find_root(a)
    root_b <- find_root(b)
    if (root_a == root_b) {
      return(invisible(NULL))
    }

    if (rank[root_a] < rank[root_b]) {
      parent[root_a] <<- root_b
    } else if (rank[root_a] > rank[root_b]) {
      parent[root_b] <<- root_a
    } else {
      parent[root_b] <<- root_a
      rank[root_a] <<- rank[root_a] + 1L
    }
    invisible(NULL)
  }

  if (length(hits) > 0L) {
    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)
    for (i in seq_along(qh)) {
      union_sets(qh[i], sh[i])
    }
  }

  roots <- vapply(seq_len(n), find_root, integer(1))
  unique_roots <- unique(roots)
  cluster_id <- match(roots, unique_roots)
  cluster_size <- tabulate(cluster_id, nbins = length(unique_roots))

  if (return == "vector") {
    return(as.integer(cluster_id))
  }

  data.frame(
    range_id = seq_len(n),
    cluster_id = as.integer(cluster_id),
    cluster_size = as.integer(cluster_size[cluster_id]),
    stringsAsFactors = FALSE
  )
}

#' Windowed Overlap Counts
#'
#' Count overlaps in sliding windows across the coordinate span of `query`.
#'
#' @param query An `IRanges` or `GRanges` query object. Windows are generated
#'   from `min(start(query))` to `max(end(query))`, separately by chromosome
#'   for `GRanges`.
#' @param subject An `IRanges` or `GRanges` subject object.
#' @param window_width Integer scalar window width.
#' @param step_width Integer scalar window step.
#' @inheritParams fast_find_overlaps
#' @inheritSection fast_find_overlaps Overlap semantics
#'
#' @return A `data.frame` containing window coordinates and overlap counts.
#' @export
#'
#' @examples
#' q <- IRanges::IRanges(start = c(1, 12), end = c(10, 20))
#' s <- IRanges::IRanges(start = c(2, 5, 15), end = c(3, 6, 16))
#' fast_window_count_overlaps(q, s, window_width = 5L, step_width = 5L)
fast_window_count_overlaps <- function(
    query,
    subject,
    window_width,
    step_width = window_width,
    max_gap = -1L,
    min_overlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore_strand = FALSE,
    threads = fast_default_threads(),
    deterministic = TRUE) {
  .assert_supported_ranges(query, "query")
  .assert_supported_ranges(subject, "subject")
  .assert_positive_integerish(window_width, "window_width")
  .assert_positive_integerish(step_width, "step_width")

  if (length(query) == 0L) {
    if (.is_granges(query)) {
      return(data.frame(
        window_id = integer(),
        seqnames = character(),
        start = integer(),
        end = integer(),
        overlap_count = integer(),
        stringsAsFactors = FALSE
      ))
    }
    return(data.frame(
      window_id = integer(),
      start = integer(),
      end = integer(),
      overlap_count = integer(),
      stringsAsFactors = FALSE
    ))
  }

  make_windows <- function(start, end) {
    starts <- seq.int(start, end, by = as.integer(step_width))
    IRanges::IRanges(start = starts, width = as.integer(window_width))
  }

  if (.is_granges(query)) {
    seq_chr <- as.character(GenomicRanges::seqnames(query))
    split_idx <- split(seq_len(length(query)), seq_chr)
    seq_levels <- names(split_idx)

    pieces <- vector("list", length(seq_levels))
    for (i in seq_along(seq_levels)) {
      idx <- split_idx[[i]]
      s <- min(IRanges::start(query[idx]))
      e <- max(IRanges::end(query[idx]))
      w <- make_windows(s, e)
      pieces[[i]] <- GenomicRanges::GRanges(
        seqnames = seq_levels[i],
        ranges = w,
        strand = "*"
      )
    }

    windows <- do.call(c, pieces)
    counts <- fast_count_overlaps(
      query = windows,
      subject = subject,
      max_gap = max_gap,
      min_overlap = min_overlap,
      type = match.arg(type),
      ignore_strand = ignore_strand,
      threads = threads,
      deterministic = deterministic
    )

    return(data.frame(
      window_id = seq_len(length(windows)),
      seqnames = as.character(GenomicRanges::seqnames(windows)),
      start = IRanges::start(windows),
      end = IRanges::end(windows),
      overlap_count = as.integer(counts),
      stringsAsFactors = FALSE
    ))
  }

  window_ir <- make_windows(min(IRanges::start(query)), max(IRanges::end(query)))
  counts <- fast_count_overlaps(
    query = window_ir,
    subject = subject,
    max_gap = max_gap,
    min_overlap = min_overlap,
    type = match.arg(type),
    ignore_strand = ignore_strand,
    threads = threads,
    deterministic = deterministic
  )

  data.frame(
    window_id = seq_len(length(window_ir)),
    start = IRanges::start(window_ir),
    end = IRanges::end(window_ir),
    overlap_count = as.integer(counts),
    stringsAsFactors = FALSE
  )
}
