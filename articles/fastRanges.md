# fastRanges: A Practical Introduction to Genomic Interval Analysis

## Overview

`fastRanges` provides high-throughput interval operations for
Bioconductor workflows built on `IRanges` and `GRanges`.

For genomic data, the central data type is `GRanges`. Each row is one
genomic interval with:

- a sequence name such as `"chr1"`
- a genomic range (`start`, `end`, or `width`)
- a strand (`"+"`, `"-"`, or `"*"`)

Typical analyses ask questions such as:

- Which peaks overlap genes or promoters?
- How many query ranges overlap each subject annotation?
- Which subject is nearest to each query?
- How can the same subject annotation be reused across many batches of
  query intervals?

`fastRanges` is designed to answer those questions while keeping its
grammar close to the Bioconductor overlap APIs in `IRanges` and
`GenomicRanges`.

``` r
library(fastRanges)
library(GenomicRanges)
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
library(S4Vectors)

data("fast_ranges_example", package = "fastRanges")
query <- fast_ranges_example$query
subject <- fast_ranges_example$subject
```

## A Mental Model

The package uses the same two-object language as Bioconductor:

- `query`: the intervals you are asking about
- `subject`: the intervals you compare the query against

If a query peak overlaps three subject genes, then:

- [`fast_find_overlaps()`](https://cparsania.github.io/fastRanges/reference/fast_find_overlaps.md)
  returns three hit pairs
- [`fast_count_overlaps()`](https://cparsania.github.io/fastRanges/reference/fast_count_overlaps.md)
  returns `3` for that query row
- [`fast_overlaps_any()`](https://cparsania.github.io/fastRanges/reference/fast_overlaps_any.md)
  returns `TRUE` for that query row

That distinction is the main decision point when choosing a function.

## Which Function Should You Use?

| Goal                                         | Function                                                                                                                                                                                                                                                 | Main output         |
|----------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------|
| Get every matching pair                      | [`fast_find_overlaps()`](https://cparsania.github.io/fastRanges/reference/fast_find_overlaps.md)                                                                                                                                                         | `Hits`              |
| Count matches per query                      | [`fast_count_overlaps()`](https://cparsania.github.io/fastRanges/reference/fast_count_overlaps.md)                                                                                                                                                       | integer vector      |
| Ask only yes/no                              | [`fast_overlaps_any()`](https://cparsania.github.io/fastRanges/reference/fast_overlaps_any.md)                                                                                                                                                           | logical vector      |
| Build a tabular overlap join                 | [`fast_overlap_join()`](https://cparsania.github.io/fastRanges/reference/fast_overlap_join.md)                                                                                                                                                           | `data.frame`        |
| Keep only matching query rows                | [`fast_semi_overlap_join()`](https://cparsania.github.io/fastRanges/reference/fast_semi_overlap_join.md)                                                                                                                                                 | `data.frame`        |
| Keep only non-matching query rows            | [`fast_anti_overlap_join()`](https://cparsania.github.io/fastRanges/reference/fast_anti_overlap_join.md)                                                                                                                                                 | `data.frame`        |
| Reuse the same subject many times            | [`fast_build_index()`](https://cparsania.github.io/fastRanges/reference/fast_build_index.md)                                                                                                                                                             | `fast_ranges_index` |
| Find nearest subject                         | [`fast_nearest()`](https://cparsania.github.io/fastRanges/reference/fast_nearest.md)                                                                                                                                                                     | `DataFrame`         |
| Summarize overlap counts by subject metadata | [`fast_count_overlaps_by_group()`](https://cparsania.github.io/fastRanges/reference/fast_count_overlaps_by_group.md)                                                                                                                                     | matrix              |
| Summarize subject scores over overlaps       | [`fast_overlap_aggregate()`](https://cparsania.github.io/fastRanges/reference/fast_overlap_aggregate.md)                                                                                                                                                 | numeric vector      |
| Merge or transform intervals themselves      | [`fast_reduce()`](https://cparsania.github.io/fastRanges/reference/fast_reduce.md), [`fast_disjoin()`](https://cparsania.github.io/fastRanges/reference/fast_disjoin.md), [`fast_gaps()`](https://cparsania.github.io/fastRanges/reference/fast_gaps.md) | ranges              |
| Summarize coverage in bins                   | [`fast_tile_coverage()`](https://cparsania.github.io/fastRanges/reference/fast_tile_coverage.md)                                                                                                                                                         | `data.frame`        |

## The Example Data

The package ships with a small reproducible example:

``` r
names(fast_ranges_example)
#> [1] "query"   "subject"
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
```

In a real analysis, `query` might be ChIP-seq peaks and `subject` might
be genes, exons, enhancers, promoters, or other genomic annotations.

This example is intentionally small so every function in the vignette
can run quickly. The structure is:

- `query`: 6 genomic intervals with metadata columns `query_id` and
  `score`
- `subject`: 7 genomic intervals with metadata columns `gene_id` and
  `biotype`

If you want file-based examples instead of in-memory objects, the same
records are also installed as BED files:

``` r
system.file("extdata", "query_peaks.bed", package = "fastRanges")
#> [1] "/home/runner/work/_temp/Library/fastRanges/extdata/query_peaks.bed"
system.file("extdata", "subject_genes.bed", package = "fastRanges")
#> [1] "/home/runner/work/_temp/Library/fastRanges/extdata/subject_genes.bed"
```

That split is deliberate:

- use `fast_ranges_example` for package examples, tutorials, and tests
- use the BED files when demonstrating import pipelines or command-line
  inputs

## Core Overlap Operations

### Return all hit pairs

``` r
hits <- fast_find_overlaps(query, subject, threads = 2)
hits
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

This result is a `Hits` object. The important parts are:

- `queryHits(hits)`: row indices from `query`
- `subjectHits(hits)`: matching row indices from `subject`

``` r
cbind(
  query_id = S4Vectors::queryHits(hits),
  subject_id = S4Vectors::subjectHits(hits)
)
#>      query_id subject_id
#> [1,]        1          1
#> [2,]        3          3
#> [3,]        4          4
#> [4,]        5          5
#> [5,]        6          6
```

### Count matches per query

``` r
counts <- fast_count_overlaps(query, subject, threads = 2)
counts
#> [1] 1 0 1 1 1 1
```

There is one value per query row:

- `0` means no subject matched
- `2` means that query matched two subject ranges

### Ask only whether a match exists

``` r
any_hits <- fast_overlaps_any(query, subject, threads = 2)
any_hits
#> [1]  TRUE FALSE  TRUE  TRUE  TRUE  TRUE
```

This is the lightest-weight summary when you need only a logical answer.

## Important Arguments

Most overlap-style functions use the same argument grammar:

- `query`
- `subject`
- `max_gap`
- `min_overlap`
- `type`
- `ignore_strand`
- `threads`
- `deterministic`

### `type`

`type` controls what “match” means.

``` text
type = "any"
query  :    |------|
subject: |------|
=> any qualifying overlap is a hit

type = "within"
subject: |--------------|
query  :    |------|
=> query must lie inside subject

type = "start"
query  :    |------|
subject:    |------------|
=> start coordinates must match or be within tolerance

type = "end"
query  :       |------|
subject: |------------|
=> end coordinates must match or be within tolerance

type = "equal"
query  :    |------|
subject:    |------|
=> start and end coordinates must match, or be within tolerance
```

### `max_gap`

`max_gap` controls how far apart two ranges may be and still count as a
hit.

- `-1` means true overlap is required
- `0` allows Bioconductor’s zero-gap tolerance behavior
- positive values allow nearby, non-overlapping ranges to count in the
  modes that support that interpretation

### `min_overlap`

`min_overlap` is the minimum overlap width, in bases.

- `0` is the least strict setting
- larger values require wider shared overlap

### `ignore_strand`

For `GRanges`:

- `FALSE` uses strand-aware Bioconductor semantics
- `TRUE` ignores strand and compares only genomic coordinates

For `IRanges`, strand does not exist, so this argument has no effect.

### `threads` and `deterministic`

- `threads` controls parallel execution
- `deterministic = TRUE` keeps stable ordering of hit pairs
- `deterministic = FALSE` can be slightly faster for large jobs when
  order is unimportant

## Bioconductor Compatibility

`fastRanges` is intended to stay close to Bioconductor overlap
semantics.

The easiest way to see that is to compare outputs directly:

``` r
ref <- GenomicRanges::findOverlaps(query, subject, ignore.strand = FALSE)

identical(
  cbind(S4Vectors::queryHits(hits), S4Vectors::subjectHits(hits)),
  cbind(S4Vectors::queryHits(ref), S4Vectors::subjectHits(ref))
)
#> [1] TRUE
```

That compatibility is important because `fastRanges` is designed to fit
into existing `IRanges` / `GRanges` workflows rather than replace those
classes.

## Return Types You Will See Most Often

``` r
ret_classes <- c(
  hits = class(fast_find_overlaps(query, subject, threads = 1))[1],
  count = class(fast_count_overlaps(query, subject, threads = 1))[1],
  any = class(fast_overlaps_any(query, subject, threads = 1))[1],
  nearest = class(fast_nearest(query, subject))[1],
  join = class(fast_overlap_join(query, subject, threads = 1))[1],
  coverage = class(fast_coverage(query))[1]
)

ret_classes
#>            hits           count             any         nearest            join 
#>          "Hits"       "integer"       "logical"        "DFrame"    "data.frame" 
#>        coverage 
#> "SimpleRleList"
```

Common return shapes:

- `Hits`: pairwise overlap relationships
- integer vector: one summary value per query
- logical vector: one yes/no value per query
- `DataFrame`: nearest-neighbor tables
- `data.frame`: join or coverage summary tables
- `GRanges` / `IRanges`: transformed intervals
- `RleList` or `Rle`: base-resolution coverage

## Reusing the Same Subject with an Index

If you will query the same `subject` repeatedly, build an index once and
reuse it.

``` r
subject_index <- fast_build_index(subject)
subject_index
#> <fast_ranges_index>
#>   subject class: GRanges 
#>   subject ranges: 7 
#>   partitions: 3
```

Now use that index in later overlap calls:

``` r
hits_from_index <- fast_find_overlaps(query, subject_index, threads = 2)

identical(
  cbind(S4Vectors::queryHits(hits_from_index), S4Vectors::subjectHits(hits_from_index)),
  cbind(S4Vectors::queryHits(hits), S4Vectors::subjectHits(hits))
)
#> [1] TRUE
```

When should you build an index?

- yes: many query batches against the same subject
- yes: benchmarking repeated overlap workloads
- no: one small one-off overlap call
- no: when you need subject metadata directly in a join output

## Overlap Joins

The join family turns overlap hits into beginner-friendly tables.

``` r
joined_inner <- fast_overlap_join(query, subject, join = "inner", threads = 2)
head(joined_inner)
#>   query_id subject_id query_seqnames query_start query_end query_width
#> 1        1          1           chr1         100       159          60
#> 3        3          3           chr1         350       419          70
#> 4        4          4           chr2          40        79          40
#> 5        5          5           chr2         300       379          80
#> 6        6          6           chr3          10        34          25
#>   query_strand query_query_id query_score subject_seqnames subject_start
#> 1            +             q1           8             chr1            90
#> 3            -             q3           6             chr1           320
#> 4            *             q4           3             chr2            20
#> 5            +             q5           9             chr2           330
#> 6            -             q6           5             chr3             1
#>   subject_end subject_width subject_strand subject_gene_id subject_biotype
#> 1         169            80              +              gA        promoter
#> 3         369            50              -              gC        enhancer
#> 4          89            70              +              gD        promoter
#> 5         419            90              *              gE       gene_body
#> 6          40            40              -              gF        promoter
```

Important columns:

- `query_id`: row index from `query`
- `subject_id`: matching row index from `subject`
- `query_*`: columns copied from `query`
- `subject_*`: columns copied from `subject`

Left join keeps all query rows:

``` r
joined_left <- fast_left_overlap_join(query, subject, threads = 2)
head(joined_left)
#>   query_id subject_id query_seqnames query_start query_end query_width
#> 1        1          1           chr1         100       159          60
#> 3        3          3           chr1         350       419          70
#> 4        4          4           chr2          40        79          40
#> 5        5          5           chr2         300       379          80
#> 6        6          6           chr3          10        34          25
#> 2        2         NA           chr1         180       269          90
#>   query_strand query_query_id query_score subject_seqnames subject_start
#> 1            +             q1           8             chr1            90
#> 3            -             q3           6             chr1           320
#> 4            *             q4           3             chr2            20
#> 5            +             q5           9             chr2           330
#> 6            -             q6           5             chr3             1
#> 2            +             q2          10             <NA>            NA
#>   subject_end subject_width subject_strand subject_gene_id subject_biotype
#> 1         169            80              +              gA        promoter
#> 3         369            50              -              gC        enhancer
#> 4          89            70              +              gD        promoter
#> 5         419            90              *              gE       gene_body
#> 6          40            40              -              gF        promoter
#> 2          NA            NA           <NA>            <NA>            <NA>
```

Semi and anti joins keep only query rows:

``` r
semi_tbl <- fast_semi_overlap_join(query, subject, threads = 2)
anti_tbl <- fast_anti_overlap_join(query, subject, threads = 2)

head(semi_tbl)
#>   query_id overlap_count query_seqnames query_start query_end query_width
#> 1        1             1           chr1         100       159          60
#> 3        3             1           chr1         350       419          70
#> 4        4             1           chr2          40        79          40
#> 5        5             1           chr2         300       379          80
#> 6        6             1           chr3          10        34          25
#>   query_strand query_query_id query_score
#> 1            +             q1           8
#> 3            -             q3           6
#> 4            *             q4           3
#> 5            +             q5           9
#> 6            -             q6           5
head(anti_tbl)
#>   query_id overlap_count query_seqnames query_start query_end query_width
#> 2        2             0           chr1         180       269          90
#>   query_strand query_query_id query_score
#> 2            +             q2          10
```

## Nearest and Directional Queries

These functions are useful when direct overlap is not the only
biological question.

``` r
nearest_tbl <- fast_nearest(query, subject)
head(as.data.frame(nearest_tbl))
#>   query_id subject_id distance
#> 1        1          1        0
#> 2        2          1       10
#> 3        3          3        0
#> 4        4          4        0
#> 5        5          5        0
#> 6        6          6        0
```

Interpretation:

- `subject_id` is the nearest subject for each matched query
- `distance = 0` means overlap
- positive distance means the nearest subject is separated by a gap

Directional helpers return subject row indices:

``` r
fast_precede(query, subject)
#> [1] NA NA  2  5 NA NA
fast_follow(query, subject)
#> [1] NA  1 NA NA  4 NA
```

## Grouped Overlap Summaries

Many analyses need more than raw hits. They need counts or scores
summarized by subject annotation.

``` r
set.seed(1)
S4Vectors::mcols(subject)$group <- sample(
  c("promoter", "enhancer"),
  length(subject),
  replace = TRUE
)
S4Vectors::mcols(subject)$score <- seq_len(length(subject))
```

### Count overlaps by group

``` r
by_group <- fast_count_overlaps_by_group(
  query,
  subject,
  group_col = "group",
  threads = 2
)

by_group
#>      enhancer promoter
#> [1,]        0        1
#> [2,]        0        0
#> [3,]        0        1
#> [4,]        0        1
#> [5,]        1        0
#> [6,]        0        1
```

Rows are query ranges. Columns are subject groups.

### Aggregate a subject score over overlaps

``` r
sum_score <- fast_overlap_aggregate(
  query,
  subject,
  value_col = "score",
  fun = "sum",
  threads = 2
)

sum_score
#> [1]  1 NA  3  4  5  6
```

This pattern is useful when subject intervals carry signal, annotation
scores, or other numeric metadata.

## Self-Overlaps, Clusters, and Sliding Windows

Self-overlaps compare a range set to itself:

``` r
self_hits <- fast_self_overlaps(query, threads = 2)
self_hits
#> Hits object with 0 hits and 0 metadata columns:
#>         from        to
#>    <integer> <integer>
#>   -------
#>   nLnode: 6 / nRnode: 6
```

Clusters identify connected components in the self-overlap graph:

``` r
clusters <- fast_cluster_overlaps(query, return = "data.frame")
head(clusters)
#>   range_id cluster_id cluster_size
#> 1        1          1            1
#> 2        2          2            1
#> 3        3          3            1
#> 4        4          4            1
#> 5        5          5            1
#> 6        6          6            1
```

Sliding-window summaries convert interval data into a regular table:

``` r
window_counts <- fast_window_count_overlaps(
  query,
  subject,
  window_width = 1000L,
  step_width = 500L,
  threads = 2
)
#> Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)

head(window_counts)
#>   window_id seqnames start  end overlap_count
#> 1         1     chr1   100 1099             3
#> 2         2     chr2    40 1039             2
#> 3         3     chr3    10 1009             2
```

## Range Algebra

These functions work on the ranges themselves rather than overlap hit
pairs.

``` r
query_reduced <- fast_reduce(query)
query_disjoint <- fast_disjoin(query)
query_gaps <- fast_gaps(query)

c(
  original = length(query),
  reduced = length(query_reduced),
  disjoint = length(query_disjoint),
  gaps = length(query_gaps)
)
#> original  reduced disjoint     gaps 
#>        6        6        6        6
```

You can also combine two range sets directly:

``` r
union_ranges <- fast_range_union(query, subject)
intersect_ranges <- fast_range_intersect(query, subject)
setdiff_ranges <- fast_range_setdiff(query, subject)

c(
  union = length(union_ranges),
  intersect = length(intersect_ranges),
  setdiff = length(setdiff_ranges)
)
#>     union intersect   setdiff 
#>        10         3         4
```

## Coverage and Binned Coverage

Coverage answers: “how many ranges cover each base position?”

``` r
query_cov <- fast_coverage(query)
query_cov
#> RleList of length 3
#> $chr1
#> integer-Rle of length 419 with 6 runs
#>   Lengths: 99 60 20 90 80 70
#>   Values :  0  1  0  1  0  1
#> 
#> $chr2
#> integer-Rle of length 379 with 4 runs
#>   Lengths:  39  40 220  80
#>   Values :   0   1   0   1
#> 
#> $chr3
#> integer-Rle of length 34 with 2 runs
#>   Lengths:  9 25
#>   Values :  0  1
```

For plotting and downstream summaries, tile-based coverage is often
easier to work with:

``` r
query_tiles <- fast_tile_coverage(
  query,
  tile_width = 1000L,
  step_width = 1000L
)

head(query_tiles)
#>   seqnames start end coverage_sum
#> 1     chr1     1 419          220
#> 2     chr2     1 379          120
#> 3     chr3     1  34           25
```

## Iterating over Large Query Sets

Use the iterator API when you do not want to materialize all hits at
once.

``` r
iter <- fast_find_overlaps_iter(query, subject_index, chunk_size = 2L, threads = 2)

fast_iter_has_next(iter)
#> [1] TRUE
chunk1 <- fast_iter_next(iter)
chunk1
#> Hits object with 1 hit and 0 metadata columns:
#>            from        to
#>       <integer> <integer>
#>   [1]         1         1
#>   -------
#>   nLnode: 6 / nRnode: 7
fast_iter_has_next(iter)
#> [1] TRUE
```

To start over:

``` r
fast_iter_reset(iter)
all_iter_hits <- fast_iter_collect(iter)
length(all_iter_hits)
#> [1] 5
```

[`fast_iter_collect()`](https://cparsania.github.io/fastRanges/reference/fast_iter_collect.md)
gathers the remaining chunks from the current iterator position. If you
already consumed some chunks and want everything again, reset the
iterator first.

## Saving and Loading an Index

``` r
tmp_index <- tempfile(fileext = ".rds")

fast_save_index(subject_index, tmp_index)
loaded_index <- fast_load_index(tmp_index)
fast_index_stats(loaded_index)
#> DataFrame with 1 row and 6 columns
#>   subject_n partition_n   block_n seqlevel_n mean_partition_size index_size_mb
#>   <integer>   <integer> <integer>  <integer>           <numeric>     <numeric>
#> 1         7           3         3          3             2.33333    0.00395966

unlink(tmp_index)
```

This is useful when index construction is expensive and the same subject
set will be queried again in a later session.

## A Small Benchmark Example

The chunk below is only a small demonstration. Publication-scale
benchmarking should use the Quarto workflow in `inst/benchmarks`.

``` r
base_time <- system.time({
  h_base <- GenomicRanges::findOverlaps(query, subject, ignore.strand = FALSE)
})

fast_time <- system.time({
  h_fast <- fast_find_overlaps(query, subject_index, ignore_strand = FALSE, threads = 2)
})

rbind(
  genomic_ranges = base_time[1:3],
  fast_ranges = fast_time[1:3]
)
#>                user.self sys.self elapsed
#> genomic_ranges     0.010        0   0.011
#> fast_ranges        0.003        0   0.002

c(base_hits = length(h_base), fast_hits = length(h_fast))
#> base_hits fast_hits 
#>         5         5
```

## Benchmark Resources

The vignette example above is intentionally small. For large benchmark
studies, repeated-query experiments, and figure-rich interpretation, use
the benchmark resources maintained in the repository:

- [Benchmark
  summary](https://github.com/cparsania/fastRanges/blob/main/inst/benchmarks/README.md)
- [Benchmark
  runner](https://github.com/cparsania/fastRanges/blob/main/inst/benchmarks/benchmark_bioc.qmd)
- [Benchmark interpretation
  report](https://github.com/cparsania/fastRanges/blob/main/inst/benchmarks/benchmark_result_interpretation.qmd)

This separation is deliberate:

- the vignette stays fast and reproducible on ordinary machines
- the benchmark runner can be executed on larger servers
- the interpretation report can be rendered later from saved result
  tables

## Package Example Files

The package also ships example BED files and benchmarking helpers:

``` r
c(
  list.files(system.file("extdata", package = "fastRanges")),
  list.files(system.file("benchmarks", package = "fastRanges"))
)
#> [1] "query_peaks.bed"    "subject_genes.bed"  "01"                
#> [4] "README.md"          "validate_outputs.R"
```

## Session Info

``` r
sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#> [1] GenomicRanges_1.62.1 Seqinfo_1.0.0        IRanges_2.44.0      
#> [4] S4Vectors_0.48.0     BiocGenerics_0.56.0  generics_0.1.4      
#> [7] fastRanges_0.99.0   
#> 
#> loaded via a namespace (and not attached):
#>  [1] cli_3.6.5         knitr_1.51        rlang_1.1.7       xfun_0.56        
#>  [5] textshaping_1.0.5 jsonlite_2.0.0    htmltools_0.5.9   ragg_1.5.1       
#>  [9] sass_0.4.10       rmarkdown_2.30    evaluate_1.0.5    jquerylib_0.1.4  
#> [13] fastmap_1.2.0     yaml_2.3.12       lifecycle_1.0.5   compiler_4.5.3   
#> [17] fs_1.6.7          Rcpp_1.1.1        XVector_0.50.0    systemfonts_1.3.2
#> [21] digest_0.6.39     R6_2.6.1          bslib_0.10.0      tools_4.5.3      
#> [25] pkgdown_2.2.0     cachem_1.1.0      desc_1.4.3
```
