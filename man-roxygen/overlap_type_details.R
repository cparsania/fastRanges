#' @section Overlap semantics:
#' Core interval semantics (ASCII schematic):
#'
#' ```text
#' type = "any"
#' query  :    |------|
#' subject: |------|
#'
#' type = "within"
#' subject: |--------------|
#' query  :    |------|
#'
#' type = "start"
#' query  :    |------|
#' subject:    |------------|
#'
#' type = "end"
#' query  :       |------|
#' subject: |------------|
#'
#' type = "equal"
#' query  :    |------|
#' subject:    |------|
#'
#' gap / min_overlap controls
#' query  : |------|      |------| : subject
#'                < gap >
#' ```
#'
#' The middle distance is the gap. A hit is allowed when this distance is
#' `<= max_gap` (for `max_gap >= 0`), and overlap width is `>= min_overlap`.
#'
#' This argument grammar is intentionally aligned with Bioconductor overlap
#' APIs (`IRanges` / `GenomicRanges`).
