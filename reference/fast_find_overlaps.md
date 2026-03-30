# Find Overlaps with Deterministic Multithreading

Compute overlap pairs between `query` and `subject` using a
multithreaded C++ backend. The result is a `Hits` object compatible with
Bioconductor workflows.

## Usage

``` r
fast_find_overlaps(
  query,
  subject,
  select = c("all", "first", "last", "arbitrary"),
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

  An `IRanges`/`GRanges` object or a `fast_ranges_index`. Use
  `fast_build_index(subject)` when the same subject is reused across
  many overlap queries.

- select:

  Character scalar controlling whether all hits are returned or a single
  subject match is selected per query.

  Use `"all"` to return a `Hits` object.

  Use `"first"`, `"last"`, or `"arbitrary"` to return an integer vector
  with one selected subject index per query and `NA` for queries with no
  hit.

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

If `select = "all"`, a
[`S4Vectors::Hits`](https://rdrr.io/pkg/S4Vectors/man/Hits-class.html)
object.

Otherwise, an integer vector of length `length(query)` containing one
selected subject index per query and `NA` when no subject matched.

## Details

This is the core matching function in `fastRanges`.

Think of it as answering the question: "for each query range, which
subject ranges satisfy my overlap rule?"

The return value is a `Hits` object. The important pieces are:

`queryHits(h)` gives the row numbers from `query`.

`subjectHits(h)` gives the matching row numbers from `subject`.

If you need only counts or a yes/no answer, prefer
[`fast_count_overlaps()`](https://cparsania.github.io/fastRanges/reference/fast_count_overlaps.md)
or
[`fast_overlaps_any()`](https://cparsania.github.io/fastRanges/reference/fast_overlaps_any.md),
because they return simpler summaries.

Compatibility notes:

`fastRanges` aims to stay close to Bioconductor overlap semantics for
supported inputs, and its outputs are routinely validated against
[`IRanges::findOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html)
/
[`GenomicRanges::findOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html).

Currently supported core input types are `IRanges` and `GRanges`.

Empty-range semantics are delegated to Bioconductor-compatible reference
behavior.

Circular genomic sequences are not currently supported and will raise an
explicit error.

`GRangesList` inputs are not currently supported and will raise an
explicit error.

Performance notes:

For one-off overlap calls, use the raw `subject`.

For repeated-query or throughput-oriented workloads, build a reusable
index once with `fast_build_index(subject)` and pass that index as
`subject`.

For maximum multithreaded throughput, consider `deterministic = FALSE`
when output order is not important.

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
h <- fast_find_overlaps(q, s, threads = 1)
length(h)
#> [1] 3
```
