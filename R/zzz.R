#' @keywords internal
.onLoad <- function(libname, pkgname) {
  if (is.null(getOption("fastRanges.threads"))) {
    options(fastRanges.threads = 1L)
  }
}
