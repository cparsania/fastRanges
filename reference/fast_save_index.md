# Save a Reusable Subject Index

Save a `fast_ranges_index` object to disk for reuse across sessions.

## Usage

``` r
fast_save_index(index, path, compress = TRUE)
```

## Arguments

- index:

  A `fast_ranges_index` object created by
  [`fast_build_index()`](https://cparsania.github.io/fastRanges/reference/fast_build_index.md).

- path:

  Output file path for the serialized index.

- compress:

  Logical scalar. If `TRUE`, uses xz compression.

## Value

Invisibly returns `path`.

## Details

Saving an index is useful when index construction is expensive and the
same subject set will be queried again in a later R session.

## Examples

``` r
s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
idx <- fast_build_index(s)
f <- tempfile(fileext = ".rds")
fast_save_index(idx, f)
unlink(f)
```
