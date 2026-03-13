# Check if an Iterator Has Remaining Chunks

Check if an Iterator Has Remaining Chunks

## Usage

``` r
fast_iter_has_next(iter)
```

## Arguments

- iter:

  A `fast_ranges_iter` object.

## Value

Logical scalar indicating whether more query chunks remain.

## Examples

``` r
q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
it <- fast_find_overlaps_iter(q, s, chunk_size = 2L)
fast_iter_has_next(it)
#> [1] TRUE
```
