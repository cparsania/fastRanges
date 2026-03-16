# Disjoin Ranges

Return non-overlapping segments induced by input ranges.

## Usage

``` r
fast_disjoin(x, ignore_strand = FALSE, with_revmap = FALSE)
```

## Arguments

- x:

  An `IRanges` or `GRanges` object.

- ignore_strand:

  Logical scalar. Ignored for non-genomic ranges.

- with_revmap:

  Logical scalar. If `TRUE`, include reverse mapping from each output
  range back to contributing input ranges.

## Value

An object of the same range class as `x`.

## Details

`fast_disjoin()` cuts the covered span of `x` into the smallest
non-overlapping pieces.

## Examples

``` r
x <- IRanges::IRanges(start = c(1, 4, 10), end = c(5, 8, 12))
fast_disjoin(x)
#> IRanges object with 4 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]         1         3         3
#>   [2]         4         5         2
#>   [3]         6         8         3
#>   [4]        10        12         3
```
