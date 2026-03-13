# Intersection of Two Range Sets

Compute range-wise intersection.

## Usage

``` r
fast_range_intersect(x, y, ignore_strand = FALSE)
```

## Arguments

- x, y:

  `IRanges` or `GRanges` objects of compatible class.

- ignore_strand:

  Logical scalar. Ignored for non-genomic ranges.

## Value

An object of the same range class as `x` and `y`.

## Details

`fast_range_intersect()` keeps only the coordinate span shared by `x`
and `y`.

## Examples

``` r
x <- IRanges::IRanges(start = c(1, 10), end = c(5, 12))
y <- IRanges::IRanges(start = c(3, 20), end = c(8, 21))
fast_range_intersect(x, y)
#> IRanges object with 1 range and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]         3         5         3
```
