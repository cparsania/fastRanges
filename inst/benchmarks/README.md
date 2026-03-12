# fastRanges Benchmark Guide

This directory contains reproducible benchmark assets for `fastRanges`.

## Files

- `benchmark_bioc.qmd`: publication-style benchmark comparing `fastRanges` and Bioconductor overlap engines.
- `benchmark_overlaps.R`: quick command-line benchmark script for a single overlap workload.

## Choosing the Right Execution Mode

Use these rules for production and benchmarking:

1. `fast_find_overlaps(query, subject, ...)` (direct subject)
   - Best for one-off calls.
   - Includes subject indexing cost on each call.
2. `idx <- fast_build_index(subject)` then `fast_find_overlaps(query, idx, ...)`
   - Best for repeated-query workloads.
   - Avoids rebuilding index for each query batch.
3. `deterministic = TRUE`
   - Use when stable output order is required across thread counts.
4. `deterministic = FALSE`
   - Use for throughput benchmarking when order does not matter.
5. `threads`
   - Set based on physical cores and observed saturation on your machine.
   - Do not assume `max cores` is always fastest.

## Fair Comparison Protocol

For fair engine-vs-engine comparisons:

1. Reuse the same `query`/`subject` objects across engines.
2. Verify hit counts are identical before comparing runtimes.
3. Use median runtime across multiple iterations.
4. For repeated-query scenarios, compare indexed mode for all engines where available.
5. Report hardware and thread grid explicitly.

## 96-Core Server Runs

The benchmark profile is tuned for high-core machines and uses thread points:

- `2, 8, 32, 64, 96` plus automatic baseline `1` for speedup normalization.

Run:

```bash
quarto render inst/benchmarks/benchmark_bioc.qmd -P max_threads:96
```

Optional overrides:

```bash
quarto render inst/benchmarks/benchmark_bioc.qmd \
  -P max_threads:96 \
  -P iterations:3 \
  -P repeated_queries:8
```

## Output Artifacts

- Rendered HTML report:
  - at your render location (Quarto default)
- Figure exports for manuscript panels:
  - `./benchmark_result/figures/*.png`
- Tabular exports for downstream/custom plotting:
  - `./benchmark_result/results/*.csv`
  - `./benchmark_result/results/benchmark_table_manifest.csv` (index of all exported tables)

## Notes

- Large defaults are intentional and may require substantial RAM.
- For CI or local smoke tests, reduce `iterations` and object sizes via Quarto parameters.
