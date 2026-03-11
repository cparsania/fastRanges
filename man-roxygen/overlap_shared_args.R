#' @param max_gap Integer scalar maximum allowed gap between query and subject.
#'   A value of `-1` requires direct overlap according to `type`. Values
#'   `>= 0` allow near matches separated by up to `max_gap` bases.
#' @param min_overlap Integer scalar minimum overlap width (in bases) required
#'   for a hit. Larger values are stricter.
#' @param type Overlap mode:
#'   `"any"`: any qualifying overlap (or gap-constrained hit when
#'   `max_gap >= 0`);
#'   `"start"`: matching start coordinates;
#'   `"end"`: matching end coordinates;
#'   `"within"`: query fully contained in subject;
#'   `"equal"`: query and subject have identical start and end.
#' @param ignore_strand Logical scalar. For `GRanges`, if `FALSE`, comparisons
#'   follow strand-aware Bioconductor semantics. Ignored for non-genomic
#'   ranges.
#' @param threads Integer scalar thread count. For repeated-query workloads,
#'   prefer `fast_build_index(subject)` and tune `threads` empirically for your
#'   hardware (scaling can plateau at high core counts).
#' @param deterministic Logical scalar. If `TRUE`, output is sorted by
#'   `(query_id, subject_id)`, giving stable ordering across thread counts.
#'   If `FALSE`, output order is unspecified and throughput is typically higher.
