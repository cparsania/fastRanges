# Follow Query Ranges

Return the index of the first subject range that is strictly before each
query range.

## Usage

``` r
fast_follow(
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

Integer vector of subject indices, with `NA` for unmatched queries.

## Details

For each query range, return the index of the first subject range that
comes strictly before the query.

## Examples

``` r
q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
fast_follow(q, s)
#> [1] NA  1  2
```
