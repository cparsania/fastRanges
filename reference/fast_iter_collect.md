# Collect All Iterator Chunks

Materialize all overlap chunks from an iterator into a single `Hits`
object.

## Usage

``` r
fast_iter_collect(iter)
```

## Arguments

- iter:

  A `fast_ranges_iter` object.

## Value

A [`S4Vectors::Hits`](https://rdrr.io/pkg/S4Vectors/man/Hits-class.html)
object.

## Details

This collects only the chunks that have not yet been consumed.

If you want all hits from the beginning after some chunks were already
read, call `fast_iter_reset(iter)` first and then
`fast_iter_collect(iter)`.

## Examples

``` r
q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
it <- fast_find_overlaps_iter(q, s, chunk_size = 2L)
fast_iter_collect(it)
#> Hits object with 3 hits and 0 metadata columns:
#>            from        to
#>       <integer> <integer>
#>   [1]         1         1
#>   [2]         2         2
#>   [3]         3         3
#>   -------
#>   nLnode: 3 / nRnode: 3
```
