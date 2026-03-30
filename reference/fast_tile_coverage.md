# Tile-Based Coverage Summary

Aggregate coverage into fixed-width tiles.

## Usage

``` r
fast_tile_coverage(
  x,
  tile_width,
  step_width = tile_width,
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

- tile_width:

  Integer scalar tile width.

- step_width:

  Integer scalar step width. Use a value smaller than `tile_width` for
  overlapping tiles and the same value for adjacent tiles.

- shift:

  Passed to
  [`coverage()`](https://rdrr.io/pkg/IRanges/man/coverage-methods.html).

- width:

  Passed to
  [`coverage()`](https://rdrr.io/pkg/IRanges/man/coverage-methods.html).

- weight:

  Passed to
  [`coverage()`](https://rdrr.io/pkg/IRanges/man/coverage-methods.html).

- method:

  Coverage method.

- threads:

  Integer scalar thread count. Reserved for API consistency.

## Value

A `data.frame` with tile coordinates and `coverage_sum`.

## Details

`fast_tile_coverage()` converts base-resolution coverage into a simpler
table of fixed-width summaries.

Each row in the result is one tile.

`coverage_sum` is the sum of coverage values across that tile.

For `GRanges`, the output also includes `seqnames`.

## Examples

``` r
x <- IRanges::IRanges(start = c(1, 4, 10), end = c(5, 8, 12))
fast_tile_coverage(x, tile_width = 5L)
#>   start end coverage_sum
#> 1     1   5            7
#> 2     6  10            4
#> 3    11  12            2
```
