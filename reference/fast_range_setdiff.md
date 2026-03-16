# Set Difference of Two Range Sets

Compute ranges in `x` that are not covered by `y`.

## Usage

``` r
fast_range_setdiff(x, y, ignore_strand = FALSE)
```

## Arguments

- x, y:

  `IRanges` or `GRanges` objects of compatible class.

- ignore_strand:

  Logical scalar. Ignored for non-genomic ranges.

## Value

An object of the same range class as `x`.

## Details

`fast_range_setdiff()` subtracts the covered span of `y` from `x`.

## Examples

``` r
x <- IRanges::IRanges(start = c(1, 10), end = c(5, 12))
y <- IRanges::IRanges(start = c(3, 20), end = c(8, 21))
fast_range_setdiff(x, y)
#> IRanges object with 2 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]         1         2         2
#>   [2]        10        12         3
```
