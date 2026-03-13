# Load a Reusable Subject Index

Load a `fast_ranges_index` object saved with
[`fast_save_index()`](https://cparsania.github.io/fastRanges/reference/fast_save_index.md).

## Usage

``` r
fast_load_index(path)
```

## Arguments

- path:

  File path to a serialized index.

## Value

A `fast_ranges_index` object.

## Details

This function validates that the file contains the fields required by
`fastRanges` before returning the object.

## Examples

``` r
s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
idx <- fast_build_index(s)
f <- tempfile(fileext = ".rds")
fast_save_index(idx, f)
idx2 <- fast_load_index(f)
unlink(f)
print(idx2)
#> <fast_ranges_index>
#>   subject class: IRanges 
#>   subject ranges: 3 
#>   partitions: 1 
```
