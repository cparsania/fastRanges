#' @param max_gap Integer scalar controlling how far apart two ranges may be
#'   and still count as a hit.
#'
#'   Use `-1` to require a true overlap.
#'
#'   Use `0` to allow touching ranges for `"any"` and to keep Bioconductor's
#'   default tolerance behavior for the other overlap modes.
#'
#'   Use positive values when you want "nearby" ranges to count as matches even
#'   if they do not overlap directly.
#'
#'   Units are bases. The meaning is intentionally aligned with
#'   `IRanges::findOverlaps()` / `GenomicRanges::findOverlaps()`.
#' @param min_overlap Integer scalar minimum overlap width, in bases.
#'
#'   `0` is the least strict setting.
#'
#'   Larger values require wider shared overlap and therefore return fewer hits.
#'
#'   This argument matters only when the selected `type` allows an actual
#'   overlap width to be measured.
#' @param type Character scalar describing what "match" means.
#'
#'   `"any"` matches any overlap that satisfies `max_gap` / `min_overlap`.
#'
#'   `"start"` matches ranges with compatible start coordinates.
#'
#'   `"end"` matches ranges with compatible end coordinates.
#'
#'   `"within"` matches queries contained inside subjects.
#'
#'   `"equal"` matches queries and subjects with the same interval, or with
#'   start/end differences no larger than `max_gap` when tolerance is allowed.
#' @param ignore_strand Logical scalar controlling strand handling for genomic
#'   ranges.
#'
#'   For `GRanges`, `FALSE` means `"+"`, `"-"`, and `"*"` are interpreted
#'   using standard Bioconductor strand rules.
#'
#'   `TRUE` means strand is ignored and only genomic coordinates are compared.
#'
#'   For `IRanges`, this argument has no effect because there is no strand.
#' @param threads Integer scalar number of worker threads to use.
#'
#'   Use `1` for the most conservative behavior and easiest debugging.
#'
#'   Use larger values on multicore machines when throughput matters.
#'
#'   For repeated-query workloads, combine a prebuilt index from
#'   `fast_build_index(subject)` with a thread count chosen empirically on your
#'   hardware.
#' @param deterministic Logical scalar controlling output order.
#'
#'   `TRUE` returns a stable order, which is useful for testing, reproducible
#'   reports, and direct comparison across thread counts.
#'
#'   `FALSE` allows the implementation to return hits in an unspecified order,
#'   which can be slightly faster for large jobs.
