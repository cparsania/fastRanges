#' @keywords internal
`%||%` <- function(lhs, rhs) {
  if (is.null(lhs)) rhs else lhs
}

#' @keywords internal
.is_granges <- function(x) {
  methods::is(x, "GenomicRanges")
}

#' @keywords internal
.is_iranges <- function(x) {
  methods::is(x, "IntegerRanges")
}

#' @keywords internal
.assert_supported_ranges <- function(x, arg_name) {
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
