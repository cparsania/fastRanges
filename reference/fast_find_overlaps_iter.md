# Create an Overlap Iterator

Create an iterator that computes overlaps in query chunks.

## Usage

``` r
fast_find_overlaps_iter(
  query,
  subject,
  chunk_size = 50000L,
  max_gap = -1L,
  min_overlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  ignore_strand = FALSE,
  threads = fast_default_threads(),
  deterministic = TRUE
)
```

## Arguments

- query:

  An `IRanges` or `GRanges` query object.

- subject:

  An `IRanges`/`GRanges` object or a `fast_ranges_index`.

- chunk_size:

  Integer scalar number of query ranges to process in one iterator step.

  Smaller values use less memory and give more frequent progress points.

  Larger values usually improve throughput but each call to
  [`fast_iter_next()`](https://cparsania.github.io/fastRanges/reference/fast_iter_next.md)
  does more work.

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

  `fastRanges` is optimized for large and throughput-oriented workloads.
  For one-off or small jobs, Bioconductor's native overlap routines may
  be competitive.

- deterministic:

  Logical scalar controlling output order.

  `TRUE` returns a stable order, which is useful for testing,
  reproducible reports, and direct comparison across thread counts.

  `FALSE` allows the implementation to return hits in an unspecified
  order, which can be noticeably faster for large multithreaded jobs
  because it avoids extra global ordering work.

## Value

A `fast_ranges_iter` iterator object.

## Details

Use the iterator API when `query` is large and you do not want to
materialize all overlap hits in one call.

Typical workflow:

create the iterator with `fast_find_overlaps_iter()`

inspect progress with
[`fast_iter_has_next()`](https://cparsania.github.io/fastRanges/reference/fast_iter_has_next.md)

pull one chunk with
[`fast_iter_next()`](https://cparsania.github.io/fastRanges/reference/fast_iter_next.md)

or collect all remaining chunks with
[`fast_iter_collect()`](https://cparsania.github.io/fastRanges/reference/fast_iter_collect.md)

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
q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
it <- fast_find_overlaps_iter(q, s, chunk_size = 2L)
fast_iter_has_next(it)
#> [1] TRUE
```
