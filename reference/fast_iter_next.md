# Advance an Overlap Iterator

Compute overlaps for the next query chunk.

## Usage

``` r
fast_iter_next(iter)
```

## Arguments

- iter:

  A `fast_ranges_iter` object.

## Value

A [`S4Vectors::Hits`](https://rdrr.io/pkg/S4Vectors/man/Hits-class.html)
object for the next chunk, with global query indices.

## Details

Each call advances the iterator state.

The returned `Hits` object uses global query indices, not chunk-local
indices, so you can combine chunk outputs safely.

## Examples

``` r
q <- IRanges::IRanges(start = c(1, 10, 20), width = c(5, 4, 3))
s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
it <- fast_find_overlaps_iter(q, s, chunk_size = 2L)
fast_iter_next(it)
#> Hits object with 2 hits and 0 metadata columns:
#>            from        to
#>       <integer> <integer>
#>   [1]         1         1
#>   [2]         2         2
#>   -------
#>   nLnode: 3 / nRnode: 3
```
