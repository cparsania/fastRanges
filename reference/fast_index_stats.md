# Index Summary Statistics

Summarize index size and partition structure.

## Usage

``` r
fast_index_stats(index, detailed = FALSE)
```

## Arguments

- index:

  A `fast_ranges_index` object.

- detailed:

  Logical scalar. If `TRUE`, returns partition-level details.

## Value

By default, a one-row
[`S4Vectors::DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
with summary fields. If `detailed = TRUE`, returns a list with `summary`
and `partitions`.

## Details

Use this function to inspect how the subject was partitioned internally
and to get a rough sense of memory footprint.

`subject_n` is the number of indexed ranges.

`partition_n` is the number of internal partitions, usually driven by
sequence structure.

`index_size_mb` is the in-memory object size, not the serialized file
size.

## Examples

``` r
s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
idx <- fast_build_index(s)
fast_index_stats(idx)
#> DataFrame with 1 row and 6 columns
#>   subject_n partition_n   block_n seqlevel_n mean_partition_size index_size_mb
#>   <integer>   <integer> <integer>  <integer>           <numeric>     <numeric>
#> 1         3           1         1          1                   3    0.00354004
```
