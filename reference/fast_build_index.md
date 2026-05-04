# Build a Reusable Subject Index

Build a sorted subject index that can be reused across repeated overlap
queries.

## Usage

``` r
fast_build_index(subject)
```

## Arguments

- subject:

  An `IRanges` or `GRanges` object.

## Value

A `fast_ranges_index` object.

## Details

The index stores a sorted representation of `subject` that is optimized
for repeated overlap queries.

Build the index once, then pass it as `subject` to
[`fast_find_overlaps()`](https://cparsania.github.io/fastRanges/reference/fast_find_overlaps.md),
[`fast_count_overlaps()`](https://cparsania.github.io/fastRanges/reference/fast_count_overlaps.md),
or other overlap-summary functions.

Use the raw `subject` object, not the index, when you need subject
metadata in the output table, for example with overlap joins.

## Examples

``` r
s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
idx <- fast_build_index(s)
idx
#> <fast_ranges_index>
#>   subject class: IRanges 
#>   subject ranges: 3 
#>   partitions: 1 
```
