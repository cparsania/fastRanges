# Coverage Across Ranges

Compute per-position coverage for input ranges.

## Usage

``` r
fast_coverage(
  x,
  shift = 0L,
  width = NULL,
  weight = 1L,
  method = c("auto", "sort", "hash"),
  threads = fast_default_threads()
)
```

## Arguments

- x:

  An `IRanges` or `GRanges` object.

- shift:

  Passed to
  [`coverage()`](https://rdrr.io/pkg/IRanges/man/coverage-methods.html).
  Use this to shift intervals before computing coverage.

- width:

  Passed to
  [`coverage()`](https://rdrr.io/pkg/IRanges/man/coverage-methods.html).
  Use this to force the output length.

- weight:

  Passed to
  [`coverage()`](https://rdrr.io/pkg/IRanges/man/coverage-methods.html).
  This can be a scalar or vector and is useful when intervals should
  contribute weighted counts instead of `1`.

- method:

  Coverage method.

- threads:

  Integer scalar thread count. Reserved for API consistency.

## Value

`Rle` (for `IRanges`) or `RleList` (for `GRanges`).

## Details

Coverage answers the question: "how many ranges cover each position?"

For `IRanges`, the result is one `Rle`.

For `GRanges`, the result is an `RleList`, usually one coverage track
per sequence level.

## Examples

``` r
x <- IRanges::IRanges(start = c(1, 4, 10), end = c(5, 8, 12))
fast_coverage(x)
#> integer-Rle of length 12 with 5 runs
#>   Lengths: 3 2 3 1 3
#>   Values : 1 2 1 0 1
```
