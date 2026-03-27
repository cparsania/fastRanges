# Reviewer Response TODO

Ordered from highest to lowest priority, based on what most directly affects
Bioconductor acceptance and how much the current package contract differs from
`findOverlaps()`.

## Priority 1: Core Compatibility / Semantics

- [x] Add `select=` support to `fast_find_overlaps()`
  - Target parity with Bioconductor-style values:
    - `all`
    - `first`
    - `last`
    - `arbitrary`
  - Default should remain `select = "all"` so current `Hits` behavior is preserved.
  - Added tests against `GenomicRanges::findOverlaps()` / `IRanges::findOverlaps()`.

- [x] Fix overlap behavior involving empty ranges
  - Match Bioconductor semantics for width-0 / empty ranges.
  - Added explicit regression tests for `IRanges` and `GRanges`.
  - Added reviewer-focused validation examples in `validate_outputs.R`.

- [x] Handle circular sequences safely
  - Detect circular `Seqinfo` and fail explicitly with a clear error.
  - Added rejection for raw overlap/count calls and index construction.
  - Added index flag and save/load backward-compatible upgrade.
  - Do not silently return incorrect results.

- [x] Handle `GRangesList` safely
  - Detect it and stop with a precise unsupported-feature message.
  - Centralized the rejection in shared input validation.
  - Added explicit tests for overlap/count/index entry points.
  - Document current limitation clearly.

## Priority 2: Benchmark / Performance Clarity

- [ ] Run reviewer-response benchmark script and capture exact outputs
  - Use reviewer’s synthetic workload.
  - Include:
    - default `fast_find_overlaps()`
    - `threads = 1`
    - `threads > 1`
    - `deterministic = TRUE/FALSE`
    - indexed vs direct
    - `findOverlaps()`
    - `findOverlaps(GNCList(...))`

- [ ] Quantify preprocessing cost explicitly
  - `fast_build_index(subject)`
  - `GenomicRanges::GNCList(subject)`
  - Include runtime and object size.

- [ ] Measure memory overhead in a reproducible way
  - Prefer scripted measurement over ad hoc `top`.
  - If possible, use `peakRAM` or another stable R-side measurement.
  - Compare:
    - direct fastRanges
    - indexed fastRanges
    - `findOverlaps()`
    - `GNCList()`

- [ ] Confirm whether default single-thread performance is slower on some systems
  - Benchmark `threads = 1, deterministic = TRUE`.
  - Separate “default behavior” from “recommended throughput mode”.
  - Avoid mixing these in the response.

## Priority 3: Documentation / Positioning

- [ ] Add an explicit compatibility section to documentation
  - README
  - vignette
  - `?fast_find_overlaps`
  - State clearly:
    - what matches `findOverlaps()`
    - what does not
    - what is unsupported
    - what is planned

- [ ] Clarify intended usage modes
  - Direct mode:
    - one-off queries
  - Indexed mode:
    - repeated-query / throughput workflows
  - Make this distinction unavoidable in docs and examples.

- [ ] Document the `deterministic` tradeoff prominently
  - Explain that `deterministic = TRUE` may reduce multithread speedup.
  - Show both modes in examples or benchmark docs.
  - Tell users when `deterministic = FALSE` is appropriate.

- [ ] Reframe performance claims more narrowly
  - Avoid implying universal replacement of `findOverlaps()`.
  - Position `fastRanges` as strongest for:
    - multithreaded workloads
    - large subject/query sets
    - repeated-query pipelines
    - indexed reuse

- [ ] Switch vignette styling to `BiocStyle`
  - Add `BiocStyle` integration.
  - Rebuild vignette and check for rendering issues.

## Priority 4: Defensive API Behavior

- [ ] Add explicit unsupported-feature checks in input validation
  - circular sequences
  - `GRangesList`
  - possibly other currently unsupported GenomicRanges variants
  - Prefer clear errors over silent semantic drift.

- [ ] Add reviewer-facing validation script outputs to repository workflow
  - Generate CSVs that can be turned into a response QMD.
  - Keep this separate from ordinary package tests.

## Priority 5: Tests

- [ ] Add focused tests for reviewer concerns
  - empty ranges
  - circular sequences
  - `GRangesList`
  - `select=`
  - deterministic ordering differences

- [x] Extend randomized validation to cover `select=`
  - `first` / `last` checked for exact identity with Bioconductor.
  - `arbitrary` checked for semantic validity against the full reference hit set.
  - Added for both `GRanges` and `IRanges`.

- [ ] Add platform-sensitive performance sanity checks where feasible
  - Not strict runtime thresholds.
  - Just enough to catch accidental regressions in direct vs indexed mode.

## Priority 6: Performance Engineering

- [ ] Investigate deterministic sorting overhead
  - Profile cost of final hit ordering.
  - Explore cheaper merge/sort strategies.
  - Decide whether this is worth changing before submission.

- [ ] Investigate memory overhead of hit materialization
  - Profile temporary allocations.
  - Check whether per-thread hit buffers or merge steps dominate.
  - Optimize only if the gain is meaningful and low-risk.

- [ ] Reassess block-index optimization branch
  - Current evidence suggests it may worsen some workloads.
  - Do not keep it unless benchmarks justify it.
  - If continuing, do so on a separate optimization branch.

## Priority 7: Reviewer Response Materials

- [ ] Prepare a point-by-point response draft
  - Thank reviewer for technical feedback.
  - Acknowledge current non-drop-in behavior honestly.
  - Separate:
    - semantics / API scope
    - default behavior
    - throughput mode

- [ ] Create a reviewer-response QMD from scripted outputs
  - Semantic compatibility table
  - reviewer workload benchmark table
  - preprocessing cost table
  - memory table
  - concise narrative per concern

- [ ] Make sure every claim in the response is backed by:
  - a test
  - a benchmark CSV
  - or an explicit documented limitation

## Suggested Execution Order

1. run reviewer-response script
2. update docs + vignette + `BiocStyle`
3. prepare response QMD
4. revisit performance engineering only if needed
