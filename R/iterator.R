#' Create an Overlap Iterator
#'
#' Create an iterator that computes overlaps in query chunks.
#'
#' @param query An `IRanges` or `GRanges` query object.
#' @param subject An `IRanges`/`GRanges` object or a `fast_ranges_index`.
#' @param chunk_size Integer scalar number of query ranges per chunk.
#' @inheritParams fast_find_overlaps
#' @inheritSection fast_find_overlaps Overlap semantics
#'
#' @return A `fast_ranges_iter` iterator object.
#' @export
#'
#' @examples
#' q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' it <- fast_find_overlaps_iter(q, s, chunk_size = 2L)
#' fast_iter_has_next(it)
fast_find_overlaps_iter <- function(
    query,
    subject,
    chunk_size = 50000L,
    max_gap = -1L,
    min_overlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore_strand = FALSE,
    threads = fast_default_threads(),
    deterministic = TRUE) {
  .assert_supported_ranges(query, "query")
  .assert_positive_integerish(chunk_size, "chunk_size")
  .assert_scalar_integerish(max_gap, "max_gap")
  .assert_scalar_integerish(min_overlap, "min_overlap")
  .assert_scalar_logical(ignore_strand, "ignore_strand")
  .assert_scalar_logical(deterministic, "deterministic")
  threads <- .normalize_threads(threads)
  type <- match.arg(type)

  subject_n <- if (inherits(subject, "fast_ranges_index")) {
    .assert_fast_ranges_index(subject, "subject")
    as.integer(subject$subject_n)
  } else {
    .assert_supported_ranges(subject, "subject")
    as.integer(length(subject))
  }

  state <- new.env(parent = emptyenv())
  state$position <- 1L

  structure(
    list(
      query = query,
      subject = subject,
      query_n = as.integer(length(query)),
      subject_n = subject_n,
      chunk_size = as.integer(chunk_size),
      max_gap = as.integer(max_gap),
      min_overlap = as.integer(min_overlap),
      type = type,
      ignore_strand = isTRUE(ignore_strand),
      threads = as.integer(threads),
      deterministic = isTRUE(deterministic),
      state = state
    ),
    class = "fast_ranges_iter"
  )
}

#' Check if an Iterator Has Remaining Chunks
#'
#' @param iter A `fast_ranges_iter` object.
#'
#' @return Logical scalar indicating whether more query chunks remain.
#' @export
#'
#' @examples
#' q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' it <- fast_find_overlaps_iter(q, s, chunk_size = 2L)
#' fast_iter_has_next(it)
fast_iter_has_next <- function(iter) {
  .assert_fast_ranges_iter(iter)
  iter$state$position <= iter$query_n
}

#' Advance an Overlap Iterator
#'
#' Compute overlaps for the next query chunk.
#'
#' @param iter A `fast_ranges_iter` object.
#'
#' @return A `S4Vectors::Hits` object for the next chunk, with global query
#'   indices.
#' @export
#'
#' @examples
#' q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' it <- fast_find_overlaps_iter(q, s, chunk_size = 2L)
#' fast_iter_next(it)
fast_iter_next <- function(iter) {
  .assert_fast_ranges_iter(iter)

  if (!fast_iter_has_next(iter)) {
    return(S4Vectors::Hits(
      from = integer(),
      to = integer(),
      nLnode = iter$query_n,
      nRnode = iter$subject_n,
      sort.by.query = FALSE
    ))
  }

  from <- iter$state$position
  to <- min(iter$query_n, from + iter$chunk_size - 1L)
  iter$state$position <- to + 1L

  q_chunk <- iter$query[from:to]
  chunk_hits <- fast_find_overlaps(
    query = q_chunk,
    subject = iter$subject,
    max_gap = iter$max_gap,
    min_overlap = iter$min_overlap,
    type = iter$type,
    ignore_strand = iter$ignore_strand,
    threads = iter$threads,
    deterministic = iter$deterministic
  )

  qh <- S4Vectors::queryHits(chunk_hits)
  sh <- S4Vectors::subjectHits(chunk_hits)
  if (length(qh) > 0L) {
    qh <- qh + from - 1L
  }

  S4Vectors::Hits(
    from = qh,
    to = sh,
    nLnode = iter$query_n,
    nRnode = iter$subject_n,
    sort.by.query = FALSE
  )
}

#' Reset an Overlap Iterator
#'
#' Reset iterator position to the first query chunk.
#'
#' @param iter A `fast_ranges_iter` object.
#'
#' @return Invisibly returns `iter`.
#' @export
#'
#' @examples
#' q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' it <- fast_find_overlaps_iter(q, s, chunk_size = 2L)
#' fast_iter_reset(it)
fast_iter_reset <- function(iter) {
  .assert_fast_ranges_iter(iter)
  iter$state$position <- 1L
  invisible(iter)
}

#' Collect All Iterator Chunks
#'
#' Materialize all overlap chunks from an iterator into a single `Hits` object.
#'
#' @param iter A `fast_ranges_iter` object.
#'
#' @return A `S4Vectors::Hits` object.
#' @export
#'
#' @examples
#' q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' it <- fast_find_overlaps_iter(q, s, chunk_size = 2L)
#' fast_iter_collect(it)
fast_iter_collect <- function(iter) {
  .assert_fast_ranges_iter(iter)

  qh <- integer()
  sh <- integer()
  while (fast_iter_has_next(iter)) {
    chunk_hits <- fast_iter_next(iter)
    qh <- c(qh, S4Vectors::queryHits(chunk_hits))
    sh <- c(sh, S4Vectors::subjectHits(chunk_hits))
  }

  if (isTRUE(iter$deterministic) && length(qh) > 1L) {
    ord <- order(qh, sh)
    qh <- qh[ord]
    sh <- sh[ord]
  }

  S4Vectors::Hits(
    from = qh,
    to = sh,
    nLnode = iter$query_n,
    nRnode = iter$subject_n,
    sort.by.query = FALSE
  )
}

#' @export
print.fast_ranges_iter <- function(x, ...) {
  remaining <- max(0L, x$query_n - x$state$position + 1L)
  cat("<fast_ranges_iter>\n", sep = "")
  cat("  query ranges:", x$query_n, "\n", sep = " ")
  cat("  subject ranges:", x$subject_n, "\n", sep = " ")
  cat("  chunk size:", x$chunk_size, "\n", sep = " ")
  cat("  remaining query ranges:", remaining, "\n", sep = " ")
  invisible(x)
}

#' @keywords internal
.assert_fast_ranges_iter <- function(iter) {
  if (!inherits(iter, "fast_ranges_iter")) {
    stop("`iter` must be a `fast_ranges_iter` object", call. = FALSE)
  }
}
