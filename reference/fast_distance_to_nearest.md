# Distance to Nearest Subject Range

Compute nearest-neighbor mapping from query ranges to subject ranges.

## Usage

``` r
fast_distance_to_nearest(
  query,
  subject,
  ignore_strand = FALSE,
  threads = fast_default_threads()
)
```

## Arguments

- query:

  An `IRanges` or `GRanges` query object.

- subject:

  An `IRanges` or `GRanges` subject object.

- ignore_strand:

  Logical scalar. Ignored for non-genomic ranges.

- threads:

  Integer scalar thread count. Included for API consistency.

## Value

A
[`S4Vectors::DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
with `query_id`, `subject_id`, and `distance` columns.

## Details

This function currently returns the same object shape as
[`fast_nearest()`](https://cparsania.github.io/fastRanges/reference/fast_nearest.md).
It is provided so users can choose the more explicit name when they care
about the distance column.

## Examples

``` r
q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
fast_distance_to_nearest(q, s)
#> DataFrame with 3 rows and 3 columns
#>    query_id subject_id  distance
#>   <integer>  <integer> <integer>
#> 1         1          1         0
#> 2         2          2         0
#> 3         3          3         0
```
