# Union of Two Range Sets

Compute range-wise union.

## Usage

``` r
fast_range_union(x, y, ignore_strand = FALSE)
```

## Arguments

- x, y:

  `IRanges` or `GRanges` objects of compatible class.

- ignore_strand:

  Logical scalar. Ignored for non-genomic ranges.

## Value

An object of the same range class as `x` and `y`.

## Details

`fast_range_union()` returns the combined interval coverage of `x` and
`y`.

## Examples

``` r
x <- IRanges::IRanges(start = c(1, 10), end = c(5, 12))
y <- IRanges::IRanges(start = c(3, 20), end = c(8, 21))
fast_range_union(x, y)
#> IRanges object with 3 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]         1         8         8
#>   [2]        10        12         3
#>   [3]        20        21         2
```
