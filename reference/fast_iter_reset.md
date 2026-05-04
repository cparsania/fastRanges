# Reset an Overlap Iterator

Reset iterator position to the first query chunk.

## Usage

``` r
fast_iter_reset(iter)
```

## Arguments

- iter:

  A `fast_ranges_iter` object.

## Value

Invisibly returns `iter`.

## Details

After reset, the next call to
[`fast_iter_next()`](https://cparsania.github.io/fastRanges/reference/fast_iter_next.md)
starts again from the first query chunk.

## Examples

``` r
q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
it <- fast_find_overlaps_iter(q, s, chunk_size = 2L)
fast_iter_reset(it)
```
