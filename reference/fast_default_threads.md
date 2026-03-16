# Default Thread Count

Returns the default thread count used by `fastRanges` overlap routines.

## Usage

``` r
fast_default_threads()
```

## Value

Integer scalar thread count.

## Details

The default is controlled by `getOption("fastRanges.threads")` and falls
back to `1L`.

This helper is mainly useful when you want package-wide thread control
without passing `threads =` to every call.

Example:

`options(fastRanges.threads = 8L)` sets the default thread count for
later calls that do not specify `threads` explicitly.

## Examples

``` r
fast_default_threads()
#> [1] 1
old_threads <- getOption("fastRanges.threads")
options(fastRanges.threads = 3L)
fast_default_threads()
#> [1] 3
options(fastRanges.threads = old_threads)
```
