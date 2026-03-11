#' Build a Reusable Subject Index
#'
#' Build a sorted subject index that can be reused across repeated overlap
#' queries.
#'
#' @param subject An `IRanges` or `GRanges` object.
#'
#' @return A `fast_ranges_index` object.
#' @export
#'
#' @examples
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' idx <- fast_build_index(s)
#' idx
fast_build_index <- function(subject) {
  .assert_supported_ranges(subject, "subject")

  sorted <- .build_sorted_subject_vectors(subject)
  structure(
    c(sorted, list(
      subject_n = length(subject),
      subject_class = class(subject)
    )),
    class = "fast_ranges_index"
  )
}

#' @export
print.fast_ranges_index <- function(x, ...) {
  cat("<fast_ranges_index>\n", sep = "")
  cat("  subject class:", paste(x$subject_class, collapse = ", "), "\n", sep = " ")
  cat("  subject ranges:", x$subject_n, "\n", sep = " ")
  cat("  partitions:", length(x$partition_keys), "\n", sep = " ")
  invisible(x)
}
