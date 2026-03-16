# Reduce Overlapping Ranges

Merge overlapping or adjacent ranges.

## Usage

``` r
fast_reduce(x, ignore_strand = FALSE, min_gap_width = 1L, with_revmap = FALSE)
```

## Arguments

- x:

  An `IRanges` or `GRanges` object.

- ignore_strand:

  Logical scalar. Ignored for non-genomic ranges.

- min_gap_width:

  Integer scalar controlling when nearby ranges should be merged.

  `1` merges overlapping or directly adjacent ranges.

  Larger values require larger gaps before two ranges are kept separate.

- with_revmap:

  Logical scalar. If `TRUE`, include reverse mapping from each output
  range back to contributing input ranges.

## Value

An object of the same range class as `x`.

## Details

`fast_reduce()` simplifies a range set by merging runs of overlapping or
near-adjacent intervals.

## Examples

``` r
x <- IRanges::IRanges(start = c(1, 4, 10), end = c(5, 8, 12))
fast_reduce(x)
#> IRanges object with 2 ranges and 0 metadata columns:
#>           start       end     width
#>       <integer> <integer> <integer>
#>   [1]         1         8         8
#>   [2]        10        12         3
```
