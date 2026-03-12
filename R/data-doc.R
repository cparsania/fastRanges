#' Example Genomic Ranges for Documentation, Tests, and Tutorials
#'
#' `fast_ranges_example` is a small reproducible dataset shipped with the
#' package. It is intended for documentation examples, teaching, unit tests,
#' and quick smoke checks of overlap workflows.
#'
#' The object is a named list with two `GRanges` components:
#'
#' - `query`: six example genomic intervals that behave like user-supplied
#'   peaks, windows, or regions of interest.
#' - `subject`: seven example genomic intervals that behave like annotation
#'   features such as genes, promoters, enhancers, or other reference ranges.
#'
#' Metadata columns are included so examples can demonstrate joins,
#' grouped summaries, and overlap aggregation:
#'
#' - `query` contains `query_id` and `score`
#' - `subject` contains `gene_id` and `biotype`
#'
#' The same records are also distributed as plain BED files in
#' `inst/extdata/query_peaks.bed` and `inst/extdata/subject_genes.bed`. Use the
#' in-memory dataset when you want a ready-to-run example in R, and use the BED
#' files when you want to demonstrate file import or command-line workflows.
#'
#' This dataset is intentionally small and synthetic. It is designed for
#' examples and tests, not as a biological reference resource.
#'
#' @name fast_ranges_example
#' @docType data
#' @format A named list with two components:
#' \describe{
#'   \item{query}{A `GRanges` object with 6 genomic intervals and metadata
#'   columns `query_id` and `score`.}
#'   \item{subject}{A `GRanges` object with 7 genomic intervals and metadata
#'   columns `gene_id` and `biotype`.}
#' }
#' @source Generated from `data-raw/make_example_data.R`, which also writes
#' matching BED files to `inst/extdata/`.
#' @examples
#' data("fast_ranges_example", package = "fastRanges")
#' query <- fast_ranges_example$query
#' subject <- fast_ranges_example$subject
#'
#' query
#' subject
#'
#' fast_find_overlaps(query, subject, threads = 1)
#' @keywords datasets

NULL
