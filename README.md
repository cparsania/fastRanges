# fastRanges

[![R-CMD-check](https://github.com/cparsania/fastRanges/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cparsania/fastRanges/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/cparsania/fastRanges/branch/main/graph/badge.svg)](https://app.codecov.io/gh/cparsania/fastRanges)
[![pkgdown site](https://img.shields.io/badge/pkgdown-site-blue)](https://cparsania.github.io/fastRanges/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-coming_soon-lightgrey)](https://bioconductor.org/)

`fastRanges` is a high-performance, Bioconductor-style package for deterministic multithreaded interval operations on `IRanges` and `GRanges`.

Website: <https://cparsania.github.io/fastRanges/>  
Source: <https://github.com/cparsania/fastRanges>

## Background 

In Bioconductor, genomic intervals are represented by `GRanges` (from
`GenomicRanges`), which encode chromosome (`seqnames`), coordinates (`ranges`),
and strand. These objects are central to analyses such as:

- ChIP-seq and ATAC-seq peak-to-gene annotation,
- overlap of variants with exons, promoters, or regulatory elements,
- transcriptome and epigenome region-level summarization.

Core overlap computations are performed with `findOverlaps()`, a Bioconductor
API provided through `IRanges` and extended for genomic coordinates in
`GenomicRanges`.

## Problem 

At modern dataset scale, interval overlap is frequently the dominant runtime
cost. This is especially pronounced when analyses require:

- very large overlap queries (millions of intervals),
- repeated query batches against the same subject set,
- downstream counting and aggregation that revisit overlap results.

These patterns are common in both exploratory analyses and production
bioinformatics pipelines.

## Limitations of Current Approaches at Scale

Bioconductor tools are mature and statistically reliable, but practical
bottlenecks can emerge for high-throughput overlap workloads:

- repeated `findOverlaps()` calls may re-incur non-trivial indexing/traversal
  costs when subjects are reused,
- multithreading behavior is not uniform across all overlap-related operations,
- overlap detection, grouping, and aggregation are often executed as separate
  stages, increasing end-to-end runtime.

## fastRanges Approach

`fastRanges` addresses these constraints while remaining Bioconductor-oriented:

- deterministic multithreaded overlap search,
- reusable subject indexing for repeated-query workflows,
- derived overlap summaries (counts, grouped counts, aggregation) in the same
  package,
- direct support for `IRanges`/`GRanges` classes and Bioconductor-style outputs.

## Bioconductor Compatibility

### Does fastRanges achieve the same things as current tools?

For core overlap tasks, yes. `fastRanges` is built to produce equivalent
overlap relationships to Bioconductor baselines for supported modes, while
adding reusable indexing and multithreaded execution.

Current scope differences:

- focus is on deterministic overlap computation and high-throughput workflows,
- not every optional argument from every Bioconductor overlap helper is exposed
  in every `fastRanges` wrapper.

### Argument conventions

`fastRanges` keeps argument names aligned with Bioconductor intent:

- `query`, `subject`, `max_gap`, `min_overlap`, `type`, `ignore_strand`
- join/nearest/range operations follow familiar Bioconductor semantics

Package-specific additions are explicit:

- `threads` for parallel execution
- `deterministic` for reproducible output ordering

### Return types

- `fast_find_overlaps()`, `fast_self_overlaps()`, iterator collectors:
  `S4Vectors::Hits`
- `fast_count_overlaps()`: integer vector
- `fast_overlaps_any()`: logical vector
- `fast_nearest()`, `fast_distance_to_nearest()`: `S4Vectors::DataFrame`
- `fast_precede()`, `fast_follow()`: integer vector (with `NA` when unmatched)
- `fast_overlap_join()` and join variants: `data.frame`
- `fast_count_overlaps_by_group()`: integer matrix
- `fast_overlap_aggregate()`: numeric vector
- range algebra helpers: `IRanges`/`GRanges` objects
- `fast_coverage()`: `Rle` or `RleList`
- `fast_tile_coverage()`: `data.frame`
- `fast_index_stats()`: `S4Vectors::DataFrame` (or list with detailed output)

## Installation

### Bioconductor (when available)

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("fastRanges")
```

### GitHub `main`

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("cparsania/fastRanges", ref = "main")
```

## Quick Start

```r
library(fastRanges)
library(GenomicRanges)

data("fast_ranges_example", package = "fastRanges")
query <- fast_ranges_example$query
subject <- fast_ranges_example$subject

hits <- fast_find_overlaps(query, subject, threads = 4)
counts <- fast_count_overlaps(query, subject, threads = 4)
joined <- fast_overlap_join(query, subject, threads = 4)

# Build once, reuse many times
subject_index <- fast_build_index(subject)
hits_indexed <- fast_find_overlaps(query, subject_index, threads = 4)

S4Vectors::mcols(subject)$type <- sample(c("A", "B"), length(subject), replace = TRUE)
S4Vectors::mcols(subject)$score <- seq_len(length(subject))
group_counts <- fast_count_overlaps_by_group(query, subject, group_col = "type", threads = 4)
sum_score <- fast_overlap_aggregate(query, subject, value_col = "score", fun = "sum", threads = 4)
```

## Function Guide

### Overlap Core

- `fast_find_overlaps()`: return overlap pairs as `Hits`
- `fast_count_overlaps()`: per-query overlap counts
- `fast_overlaps_any()`: per-query logical overlap flag
- `fast_build_index()`: build reusable subject index

### Join Operations

- `fast_overlap_join()`: inner/left overlap join
- `fast_inner_overlap_join()`, `fast_left_overlap_join()`
- `fast_semi_overlap_join()`, `fast_anti_overlap_join()`

### Nearest and Directional Queries

- `fast_nearest()`, `fast_distance_to_nearest()`
- `fast_precede()`, `fast_follow()`

### Overlap Analytics

- `fast_count_overlaps_by_group()`
- `fast_overlap_aggregate()`
- `fast_window_count_overlaps()`
- `fast_self_overlaps()`, `fast_cluster_overlaps()`

### Range Algebra

- `fast_reduce()`, `fast_disjoin()`, `fast_gaps()`
- `fast_range_union()`, `fast_range_intersect()`, `fast_range_setdiff()`

### Coverage and Binning

- `fast_coverage()`
- `fast_tile_coverage()`

### Index Persistence

- `fast_save_index()`, `fast_load_index()`, `fast_index_stats()`

### Chunked Iteration

- `fast_find_overlaps_iter()`
- `fast_iter_has_next()`, `fast_iter_next()`
- `fast_iter_reset()`, `fast_iter_collect()`

## Benchmarking

For benchmarking and result interpretation, use the following entry points:

- [Benchmark Summary](inst/benchmarks/README.md)
- [Benchmark Runner](inst/benchmarks/benchmark_bioc.qmd)
- [Benchmark Interpretation Report](inst/benchmarks/benchmark_result_interpretation.qmd)

Recommended workflow:

1. Run the large benchmark on the target machine with `benchmark_bioc.qmd`.
2. Review the GitHub-friendly benchmark summary in `inst/benchmarks/README.md`.
3. Render `benchmark_result_interpretation.qmd` for figure-rich post hoc
   analysis from saved result tables.

Render with:

```bash
quarto render inst/benchmarks/benchmark_bioc.qmd -P max_threads:96
```

Interpret saved results with:

```bash
quarto render inst/benchmarks/benchmark_result_interpretation.qmd
```
