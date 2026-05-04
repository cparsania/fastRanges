# Gaps Between Ranges

Compute uncovered regions between ranges.

## Usage

``` r
fast_gaps(x, start = NULL, end = NULL, ignore_strand = FALSE)
```

## Arguments

- x:

  An `IRanges` or `GRanges` object.

- start, end:

  Optional integer bounds for non-genomic ranges. These are most useful
  for `IRanges`. For `GRanges`, sequence lengths are usually taken from
  `seqinfo(x)`.

- ignore_strand:

  Logical scalar. Ignored for non-genomic ranges.

## Value

An object of the same range class as `x`.

## Details

`fast_gaps()` returns the regions not covered by `x` inside the
requested bounds.

## Examples

``` r
x <- IRanges::IRanges(start = c(1, 4, 10), end = c(5, 8, 12))
fast_gaps(x, start = 1L, end = 15L)
#> IRanges object with 2 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]         9         9         1
#>   [2]        13        15         3
```
