# Example Genomic Ranges for Documentation, Tests, and Tutorials

`fast_ranges_example` is a small reproducible dataset shipped with the
package. It is intended for documentation examples, teaching, unit
tests, and quick smoke checks of overlap workflows.

## Format

A named list with two components:

- query:

  A `GRanges` object with 6 genomic intervals and metadata columns
  `query_id` and `score`.

- subject:

  A `GRanges` object with 7 genomic intervals and metadata columns
  `gene_id` and `biotype`.

## Source

Generated from `data-raw/make_example_data.R`, which also writes
matching BED files to `inst/extdata/`.

## Details

The object is a named list with two `GRanges` components:

- `query`: six example genomic intervals that behave like user-supplied
  peaks, windows, or regions of interest.

- `subject`: seven example genomic intervals that behave like annotation
  features such as genes, promoters, enhancers, or other reference
  ranges.

Metadata columns are included so examples can demonstrate joins, grouped
summaries, and overlap aggregation:

- `query` contains `query_id` and `score`

- `subject` contains `gene_id` and `biotype`

The same records are also distributed as plain BED files in
`inst/extdata/query_peaks.bed` and `inst/extdata/subject_genes.bed`. Use
the in-memory dataset when you want a ready-to-run example in R, and use
the BED files when you want to demonstrate file import or command-line
workflows.

This dataset is intentionally small and synthetic. It is designed for
examples and tests, not as a biological reference resource.

## Examples

``` r
data("fast_ranges_example", package = "fastRanges")
query <- fast_ranges_example$query
subject <- fast_ranges_example$subject

query
#> GRanges object with 6 ranges and 2 metadata columns:
#>       seqnames    ranges strand |    query_id     score
#>          <Rle> <IRanges>  <Rle> | <character> <integer>
#>   [1]     chr1   100-159      + |          q1         8
#>   [2]     chr1   180-269      + |          q2        10
#>   [3]     chr1   350-419      - |          q3         6
#>   [4]     chr2     40-79      * |          q4         3
#>   [5]     chr2   300-379      + |          q5         9
#>   [6]     chr3     10-34      - |          q6         5
#>   -------
#>   seqinfo: 3 sequences from an unspecified genome; no seqlengths
subject
#> GRanges object with 7 ranges and 2 metadata columns:
#>       seqnames    ranges strand |     gene_id     biotype
#>          <Rle> <IRanges>  <Rle> | <character> <character>
#>   [1]     chr1    90-169      + |          gA    promoter
#>   [2]     chr1   210-309      - |          gB    enhancer
#>   [3]     chr1   320-369      - |          gC    enhancer
#>   [4]     chr2     20-89      + |          gD    promoter
#>   [5]     chr2   330-419      * |          gE   gene_body
#>   [6]     chr3      1-40      - |          gF    promoter
#>   [7]     chr3    80-114      + |          gG    enhancer
#>   -------
#>   seqinfo: 3 sequences from an unspecified genome; no seqlengths

fast_find_overlaps(query, subject, threads = 1)
#> Hits object with 5 hits and 0 metadata columns:
#>            from        to
#>       <integer> <integer>
#>   [1]         1         1
#>   [2]         3         3
#>   [3]         4         4
#>   [4]         5         5
#>   [5]         6         6
#>   -------
#>   nLnode: 6 / nRnode: 7
```
