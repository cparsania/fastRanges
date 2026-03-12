#' @section Overlap semantics:
#' `query` is the range set you ask about. `subject` is the range set you
#' compare it against.
#'
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
#' Beginner-friendly interpretation:
#'
#' `type = "any"` asks "do these ranges touch or overlap closely enough to
#' count?"
#'
#' `type = "start"` and `type = "end"` are useful when interval boundaries are
#' biologically meaningful, for example transcription start or end sites.
#'
#' `type = "within"` asks whether each query lies inside a subject interval.
#'
#' `type = "equal"` asks whether query and subject describe the same interval,
#' optionally with endpoint tolerance when `max_gap >= 0`.
#'
#' This argument grammar is intentionally aligned with Bioconductor overlap
#' APIs (`IRanges` / `GenomicRanges`).
