# Cluster Overlapping Ranges

Assign each range to an overlap-connected component.

## Usage

``` r
fast_cluster_overlaps(
  x,
  max_gap = -1L,
  min_overlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  ignore_strand = FALSE,
  threads = fast_default_threads(),
  deterministic = TRUE,
  return = c("vector", "data.frame")
)
```

## Arguments

- x:

  An `IRanges` or `GRanges` object.

- max_gap:

  Integer scalar controlling how far apart two ranges may be and still
  count as a hit.

  Use `-1` to require a true overlap.

  Use `0` to allow touching ranges for `"any"` and to keep
  Bioconductor's default tolerance behavior for the other overlap modes.

  Use positive values when you want "nearby" ranges to count as matches
  even if they do not overlap directly.

  Units are bases. The meaning is intentionally aligned with
  [`IRanges::findOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html)
  /
  [`GenomicRanges::findOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html).

- min_overlap:

  Integer scalar minimum overlap width, in bases.

  `0` is the least strict setting.

  Larger values require wider shared overlap and therefore return fewer
  hits.

  This argument matters only when the selected `type` allows an actual
  overlap width to be measured.

- type:

  Character scalar describing what "match" means.

  `"any"` matches any overlap that satisfies `max_gap` / `min_overlap`.

  `"start"` matches ranges with compatible start coordinates.

  `"end"` matches ranges with compatible end coordinates.

  `"within"` matches queries contained inside subjects.

  `"equal"` matches queries and subjects with the same interval, or with
  start/end differences no larger than `max_gap` when tolerance is
  allowed.

- ignore_strand:

  Logical scalar controlling strand handling for genomic ranges.

  For `GRanges`, `FALSE` means `"+"`, `"-"`, and `"*"` are interpreted
  using standard Bioconductor strand rules.

  `TRUE` means strand is ignored and only genomic coordinates are
  compared.

  For `IRanges`, this argument has no effect because there is no strand.

- threads:

  Integer scalar number of worker threads to use.

  Use `1` for the most conservative behavior and easiest debugging.

  Use larger values on multicore machines when throughput matters.

  For repeated-query workloads, combine a prebuilt index from
  `fast_build_index(subject)` with a thread count chosen empirically on
  your hardware.

- deterministic:

  Logical scalar controlling output order.

  `TRUE` returns a stable order, which is useful for testing,
  reproducible reports, and direct comparison across thread counts.

  `FALSE` allows the implementation to return hits in an unspecified
  order, which can be slightly faster for large jobs.

- return:

  One of `"vector"` or `"data.frame"`.

## Value

If `return = "vector"`, an integer vector of cluster IDs with one
element per range in `x`. If `return = "data.frame"`, a `data.frame`
with `range_id`, `cluster_id`, and `cluster_size`.

## Details

Two ranges are assigned to the same cluster when they are connected by a
chain of overlaps under the provided overlap settings.

Two ranges can end up in the same cluster even if they do not overlap
directly, as long as they are linked by intermediate overlaps.

Example: if range 1 overlaps range 2, and range 2 overlaps range 3, then
all three belong to one cluster.

## Overlap semantics

`query` is the range set you ask about. `subject` is the range set you
compare it against.

Core interval semantics (ASCII schematic):

The middle distance is the gap. A hit is allowed when this distance is
`<= max_gap` (for `max_gap >= 0`), and overlap width is
`>= min_overlap`.

Beginner-friendly interpretation:

`type = "any"` asks "do these ranges touch or overlap closely enough to
count?"

`type = "start"` and `type = "end"` are useful when interval boundaries
are biologically meaningful, for example transcription start or end
sites.

`type = "within"` asks whether each query lies inside a subject
interval.

`type = "equal"` asks whether query and subject describe the same
interval, optionally with endpoint tolerance when `max_gap >= 0`.

This argument grammar is intentionally aligned with Bioconductor overlap
APIs (`IRanges` / `GenomicRanges`).

## Examples

``` r
x <- IRanges::IRanges(start = c(1, 3, 10, 11, 20), end = c(5, 8, 12, 14, 22))
fast_cluster_overlaps(x)
#> [1] 1 1 2 2 3
```
