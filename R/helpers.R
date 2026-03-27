#' @keywords internal
`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

#' @keywords internal
.is_granges <- function(x) {
  isTRUE(methods::is(x, "GenomicRanges"))
}

#' @keywords internal
.is_iranges <- function(x) {
  isTRUE(methods::is(x, "IntegerRanges"))
}

#' @keywords internal
.assert_supported_ranges <- function(x, arg_name) {
  if (isTRUE(methods::is(x, "GenomicRangesList"))) {
    stop(sprintf("`%s` must not be a GRangesList; GRangesList is not currently supported", arg_name), call. = FALSE)
  }
  if (!(.is_granges(x) || .is_iranges(x))) {
    stop(sprintf("`%s` must inherit from IntegerRanges or GenomicRanges", arg_name), call. = FALSE)
  }
}

#' @keywords internal
.assert_scalar_integerish <- function(x, arg_name) {
  if (length(x) != 1L || is.na(x) || !is.numeric(x)) {
    stop(sprintf("`%s` must be a single non-missing numeric value", arg_name), call. = FALSE)
  }
}

#' @keywords internal
.assert_scalar_logical <- function(x, arg_name) {
  if (length(x) != 1L || is.na(x) || !is.logical(x)) {
    stop(sprintf("`%s` must be a single non-missing logical value", arg_name), call. = FALSE)
  }
}

#' @keywords internal
.assert_scalar_character <- function(x, arg_name) {
  if (length(x) != 1L || is.na(x) || !is.character(x) || nchar(x) == 0L) {
    stop(sprintf("`%s` must be a single non-empty character value", arg_name), call. = FALSE)
  }
}

#' @keywords internal
.assert_positive_integerish <- function(x, arg_name) {
  .assert_scalar_integerish(x, arg_name)
  if (x < 1) {
    stop(sprintf("`%s` must be >= 1", arg_name), call. = FALSE)
  }
}

#' @keywords internal
.normalize_threads <- function(threads) {
  .assert_scalar_integerish(threads, "threads")
  as.integer(max(1L, threads))
}

#' @keywords internal
.assert_fast_ranges_index <- function(x, arg_name = "index") {
  if (!inherits(x, "fast_ranges_index")) {
    stop(sprintf("`%s` must be a `fast_ranges_index` object", arg_name), call. = FALSE)
  }
}

#' @keywords internal
.empty_like_query <- function(query) {
  if (.is_granges(query)) {
    return(query[FALSE])
  }
  if (.is_iranges(query)) {
    return(query[FALSE])
  }
  stop("Unsupported query class", call. = FALSE)
}

#' @keywords internal
.has_empty_ranges <- function(x) {
  isTRUE(any(IRanges::width(x) == 0L))
}

#' @keywords internal
.has_circular_sequences <- function(x) {
  if (!.is_granges(x) || length(x) == 0L) {
    return(FALSE)
  }
  si <- GenomeInfoDb::seqinfo(x)
  flags <- GenomeInfoDb::isCircular(si)
  used <- unique(as.character(GenomicRanges::seqnames(x)))
  if (length(used) == 0L) {
    return(FALSE)
  }
  any(flags[used] %in% TRUE)
}

#' @keywords internal
.assert_no_circular_ranges <- function(x, arg_name) {
  if (.has_circular_sequences(x)) {
    stop(sprintf("`%s` contains circular sequences, which are not currently supported", arg_name), call. = FALSE)
  }
}

#' @keywords internal
.restore_subject_from_index <- function(index) {
  .assert_fast_ranges_index(index, "index")

  n <- as.integer(index$subject_n)
  if (n == 0L) {
    if ("GRanges" %in% index$subject_class || "GenomicRanges" %in% index$subject_class) {
      return(GenomicRanges::GRanges())
    }
    return(IRanges::IRanges())
  }

  ord <- as.integer(index$subject_original_index)
  start <- integer(n)
  end <- integer(n)
  seq_id <- integer(n)
  strand_id <- integer(n)

  start[ord] <- as.integer(index$subject_start)
  end[ord] <- as.integer(index$subject_end)
  seq_id[ord] <- as.integer(index$subject_seq)
  strand_id[ord] <- as.integer(index$subject_strand)

  if ("GRanges" %in% index$subject_class || "GenomicRanges" %in% index$subject_class) {
    seq_levels <- index$seq_levels
    seqnames <- seq_levels[seq_id]
    strand <- rep.int("*", n)
    strand[strand_id == 1L] <- "+"
    strand[strand_id == 2L] <- "-"
    return(GenomicRanges::GRanges(
      seqnames = seqnames,
      ranges = IRanges::IRanges(start = start, end = end),
      strand = strand
    ))
  }

  IRanges::IRanges(start = start, end = end)
}

#' @keywords internal
.find_overlaps_reference <- function(
    query,
    subject,
    select = c("all", "first", "last", "arbitrary"),
    max_gap = -1L,
    min_overlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore_strand = FALSE) {
  select <- match.arg(select)
  type <- match.arg(type)

  if (.is_granges(query) || .is_granges(subject)) {
    return(GenomicRanges::findOverlaps(
      query, subject,
      maxgap = max_gap,
      minoverlap = min_overlap,
      type = type,
      select = select,
      ignore.strand = ignore_strand
    ))
  }

  IRanges::findOverlaps(
    query, subject,
    maxgap = max_gap,
    minoverlap = min_overlap,
    type = type,
    select = select
  )
}

#' @keywords internal
.count_overlaps_reference <- function(
    query,
    subject,
    max_gap = -1L,
    min_overlap = 0L,
    type = c("any", "start", "end", "within", "equal"),
    ignore_strand = FALSE) {
  type <- match.arg(type)

  if (.is_granges(query) || .is_granges(subject)) {
    return(GenomicRanges::countOverlaps(
      query, subject,
      maxgap = max_gap,
      minoverlap = min_overlap,
      type = type,
      ignore.strand = ignore_strand
    ))
  }

  IRanges::countOverlaps(
    query, subject,
    maxgap = max_gap,
    minoverlap = min_overlap,
    type = type
  )
}

#' Default Thread Count
#'
#' Returns the default thread count used by `fastRanges` overlap routines.
#'
#' The default is controlled by `getOption("fastRanges.threads")` and falls
#' back to `1L`.
#'
#' @details
#' This helper is mainly useful when you want package-wide thread control
#' without passing `threads =` to every call.
#'
#' Example:
#'
#' `options(fastRanges.threads = 8L)` sets the default thread count for later
#' calls that do not specify `threads` explicitly.
#'
#' @return Integer scalar thread count.
#' @export
#'
#' @examples
#' fast_default_threads()
#' old_threads <- getOption("fastRanges.threads")
#' options(fastRanges.threads = 3L)
#' fast_default_threads()
#' options(fastRanges.threads = old_threads)
fast_default_threads <- function() {
  .normalize_threads(getOption("fastRanges.threads", 1L))
}
