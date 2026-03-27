#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(fastRanges)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})

get_arg <- function(args, key, default = NULL) {
  prefix <- paste0("--", key, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0L) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

get_flag <- function(args, key, default = FALSE) {
  val <- get_arg(args, key, default = NULL)
  if (is.null(val)) return(isTRUE(default))
  tolower(val) %in% c("1", "true", "t", "yes", "y")
}

safe_pkg_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    return(NA_character_)
  }
  as.character(utils::packageVersion(pkg))
}

canon_hits <- function(h) {
  x <- cbind(S4Vectors::queryHits(h), S4Vectors::subjectHits(h))
  if (nrow(x) == 0L) return(x)
  x[order(x[, 1], x[, 2]), , drop = FALSE]
}

equal_hits <- function(x, y) {
  identical(canon_hits(x), canon_hits(y))
}

capture_error <- function(expr) {
  tryCatch(
    list(ok = TRUE, value = force(expr), message = NA_character_),
    error = function(e) {
      list(ok = FALSE, value = NULL, message = conditionMessage(e))
    }
  )
}

make_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, quote = TRUE, na = "")
  invisible(path)
}

session_pkg_versions <- function() {
  data.frame(
    package = c(
      "R", "fastRanges", "GenomicRanges", "IRanges", "S4Vectors",
      "GenomeInfoDb", "BiocGenerics"
    ),
    version = c(
      paste(R.version$major, R.version$minor, sep = "."),
      safe_pkg_version("fastRanges"),
      safe_pkg_version("GenomicRanges"),
      safe_pkg_version("IRanges"),
      safe_pkg_version("S4Vectors"),
      safe_pkg_version("GenomeInfoDb"),
      safe_pkg_version("BiocGenerics")
    ),
    stringsAsFactors = FALSE
  )
}

record_check <- function(
    concern,
    scenario,
    status,
    ref_supported,
    fast_supported,
    ref_metric = NA_character_,
    fast_metric = NA_character_,
    message = NA_character_) {
  data.frame(
    concern = concern,
    scenario = scenario,
    status = status,
    reference_supported = ref_supported,
    fastRanges_supported = fast_supported,
    reference_metric = ref_metric,
    fastRanges_metric = fast_metric,
    message = message,
    stringsAsFactors = FALSE
  )
}

time_engine <- function(label, expr_fun, iters = 3L) {
  out <- vector("list", iters)
  for (i in seq_len(iters)) {
    gc()
    tm <- system.time({
      hits <- expr_fun()
    })
    out[[i]] <- data.frame(
      engine = label,
      iter = i,
      elapsed_sec = unname(tm["elapsed"]),
      hits = length(hits),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

time_object_build <- function(label, expr_fun, iters = 3L) {
  out <- vector("list", iters)
  for (i in seq_len(iters)) {
    gc()
    tm <- system.time({
      obj <- expr_fun()
    })
    out[[i]] <- data.frame(
      engine = label,
      iter = i,
      elapsed_sec = unname(tm["elapsed"]),
      object_size_mb = as.numeric(utils::object.size(obj)) / 1024^2,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

summarise_runtime <- function(x, baseline_engine) {
  agg <- aggregate(elapsed_sec ~ engine, data = x, FUN = median)
  q25 <- aggregate(elapsed_sec ~ engine, data = x, FUN = function(z) as.numeric(stats::quantile(z, 0.25)))
  q75 <- aggregate(elapsed_sec ~ engine, data = x, FUN = function(z) as.numeric(stats::quantile(z, 0.75)))
  names(q25)[2] <- "q25_sec"
  names(q75)[2] <- "q75_sec"
  out <- Reduce(function(a, b) merge(a, b, by = "engine", sort = FALSE), list(agg, q25, q75))
  names(out)[2] <- "median_sec"
  baseline <- out$median_sec[out$engine == baseline_engine]
  if (length(baseline) == 1L && !is.na(baseline)) {
    out$speedup_vs_baseline <- baseline / out$median_sec
  } else {
    out$speedup_vs_baseline <- NA_real_
  }
  out[order(out$median_sec), ]
}

make_reviewer_ranges <- function() {
  set.seed(333)
  gr1 <- GRanges(
    sample(LETTERS, 3e6, replace = TRUE),
    IRanges(
      sample(5000000L, 3e6, replace = TRUE),
      width = sample(10:500, 3e6, replace = TRUE)
    )
  )
  gr2 <- GRanges(
    sample(LETTERS, 1e6, replace = TRUE),
    IRanges(
      sample(5000000L, 1e6, replace = TRUE),
      width = sample(44000L, 1e6, replace = TRUE)
    )
  )
  list(query = gr1, subject = gr2)
}

run_semantic_checks <- function(threads) {
  out <- list()

  q_empty <- GRanges("chr1", IRanges(start = c(5L, 8L), width = c(0L, 3L)))
  s_empty <- GRanges("chr1", IRanges(start = c(5L, 8L), width = c(0L, 2L)))
  ref_empty <- capture_error(GenomicRanges::findOverlaps(q_empty, s_empty))
  got_empty <- capture_error(fastRanges::fast_find_overlaps(q_empty, s_empty, threads = threads))
  out[[length(out) + 1L]] <- record_check(
    concern = "empty_ranges",
    scenario = "GRanges width-0 overlaps",
    status = if (ref_empty$ok && got_empty$ok && equal_hits(ref_empty$value, got_empty$value)) "match" else "mismatch",
    ref_supported = ref_empty$ok,
    fast_supported = got_empty$ok,
    ref_metric = if (ref_empty$ok) as.character(length(ref_empty$value)) else NA_character_,
    fast_metric = if (got_empty$ok) as.character(length(got_empty$value)) else NA_character_,
    message = if (!got_empty$ok) got_empty$message else if (!ref_empty$ok) ref_empty$message else NA_character_
  )

  si <- GenomeInfoDb::Seqinfo("chrC", seqlengths = 100L, isCircular = TRUE)
  q_circ <- suppressWarnings(GRanges("chrC", IRanges(start = 95L, width = 10L), seqinfo = si))
  s_circ <- suppressWarnings(GRanges("chrC", IRanges(start = 2L, width = 8L), seqinfo = si))
  ref_circ <- capture_error(GenomicRanges::findOverlaps(q_circ, s_circ))
  got_circ <- capture_error(fastRanges::fast_find_overlaps(q_circ, s_circ, threads = threads))
  circ_status <- if (ref_circ$ok && got_circ$ok && equal_hits(ref_circ$value, got_circ$value)) "match" else if (ref_circ$ok && !got_circ$ok) "unsupported_or_error" else "mismatch"
  out[[length(out) + 1L]] <- record_check(
    concern = "circular_sequences",
    scenario = "circular seqinfo wrap-around overlap",
    status = circ_status,
    ref_supported = ref_circ$ok,
    fast_supported = got_circ$ok,
    ref_metric = if (ref_circ$ok) as.character(length(ref_circ$value)) else NA_character_,
    fast_metric = if (got_circ$ok) as.character(length(got_circ$value)) else NA_character_,
    message = if (!got_circ$ok) got_circ$message else if (!ref_circ$ok) ref_circ$message else NA_character_
  )

  q_grl <- GenomicRanges::GRangesList(
    a = GRanges("chr1", IRanges(start = c(1L, 10L), width = c(3L, 4L))),
    b = GRanges("chr1", IRanges(start = 20L, width = 5L))
  )
  s_gr <- GRanges("chr1", IRanges(start = c(2L, 12L, 19L), width = c(2L, 3L, 4L)))
  ref_grl <- capture_error(GenomicRanges::findOverlaps(q_grl, s_gr))
  got_grl <- capture_error(fastRanges::fast_find_overlaps(q_grl, s_gr, threads = threads))
  out[[length(out) + 1L]] <- record_check(
    concern = "GRangesList_support",
    scenario = "GRangesList query vs GRanges subject",
    status = if (ref_grl$ok && got_grl$ok) "supported" else "unsupported_or_error",
    ref_supported = ref_grl$ok,
    fast_supported = got_grl$ok,
    ref_metric = if (ref_grl$ok) as.character(length(ref_grl$value)) else NA_character_,
    fast_metric = if (got_grl$ok) as.character(length(got_grl$value)) else NA_character_,
    message = if (!got_grl$ok) got_grl$message else NA_character_
  )

  has_select <- "select" %in% names(formals(fastRanges::fast_find_overlaps))
  out[[length(out) + 1L]] <- record_check(
    concern = "select_argument",
    scenario = "presence of select= API argument",
    status = if (has_select) "supported" else "missing",
    ref_supported = TRUE,
    fast_supported = has_select,
    ref_metric = "select in findOverlaps",
    fast_metric = if (has_select) "select in fast_find_overlaps" else "select missing",
    message = if (has_select) NA_character_ else "fast_find_overlaps() has no select= argument"
  )

  do.call(rbind, out)
}

maybe_peakram <- function(expr_fun) {
  if (!requireNamespace("peakRAM", quietly = TRUE)) {
    return(list(
      available = FALSE,
      peak_ram_mb = NA_real_,
      elapsed_sec = NA_real_
    ))
  }
  res <- peakRAM::peakRAM({
    hits <- expr_fun()
  })
  list(
    available = TRUE,
    peak_ram_mb = as.numeric(res$Peak_RAM_Used_MiB[1]),
    elapsed_sec = as.numeric(res$Elapsed_Time_sec[1])
  )
}

run_peak_ram_checks <- function(gr1, gr2, idx, threads_max) {
  engines <- list(
    list(name = "fastRanges_default", fn = function() fastRanges::fast_find_overlaps(gr1, gr2)),
    list(name = "fastRanges_direct_t1_det_true", fn = function() fastRanges::fast_find_overlaps(gr1, gr2, threads = 1L, deterministic = TRUE)),
    list(name = "fastRanges_direct_tmax_det_false", fn = function() fastRanges::fast_find_overlaps(gr1, gr2, threads = threads_max, deterministic = FALSE)),
    list(name = "fastRanges_indexed_t1_det_true", fn = function() fastRanges::fast_find_overlaps(gr1, idx, threads = 1L, deterministic = TRUE)),
    list(name = "fastRanges_indexed_tmax_det_false", fn = function() fastRanges::fast_find_overlaps(gr1, idx, threads = threads_max, deterministic = FALSE)),
    list(name = "GenomicRanges_findOverlaps", fn = function() GenomicRanges::findOverlaps(gr1, gr2)),
    list(name = "GenomicRanges_GNCList", fn = function() GenomicRanges::findOverlaps(gr1, GenomicRanges::GNCList(gr2)))
  )

  out <- lapply(engines, function(engine) {
    gc()
    res <- maybe_peakram(engine$fn)
    data.frame(
      engine = engine$name,
      peakram_available = res$available,
      peak_ram_mb = res$peak_ram_mb,
      elapsed_sec = res$elapsed_sec,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

args <- commandArgs(trailingOnly = TRUE)
out_dir <- make_dir(get_arg(args, "out_dir", "inst/benchmarks/reviewer_response_results"))
iters <- as.integer(get_arg(args, "iters", "3"))
threads_arg <- get_arg(args, "threads", NULL)
run_peak_ram <- get_flag(args, "run_peak_ram", default = FALSE)

physical_cores <- parallel::detectCores(logical = FALSE)
if (is.na(physical_cores) || physical_cores < 1L) physical_cores <- 1L

threads_grid <- if (is.null(threads_arg)) {
  unique(as.integer(c(1L, 8L, physical_cores)))
} else {
  unique(as.integer(strsplit(threads_arg, ",", fixed = TRUE)[[1L]]))
}
threads_grid <- sort(unique(pmax(1L, pmin(threads_grid, physical_cores))))
if (!1L %in% threads_grid) threads_grid <- sort(unique(c(1L, threads_grid)))

cat("Reviewer response checks started\n")
cat("out_dir=", out_dir, "\n", sep = "")
cat("iters=", iters, " threads=", paste(threads_grid, collapse = ","), " run_peak_ram=", run_peak_ram, "\n", sep = "")

metadata <- data.frame(
  timestamp = as.character(Sys.time()),
  platform = paste(R.version$platform, R.version$os),
  physical_cores = physical_cores,
  fastRanges_default_threads = fastRanges::fast_default_threads(),
  iterations = iters,
  threads_grid = paste(threads_grid, collapse = ","),
  stringsAsFactors = FALSE
)
write_csv(metadata, file.path(out_dir, "run_metadata.csv"))
write_csv(session_pkg_versions(), file.path(out_dir, "package_versions.csv"))

semantic_checks <- run_semantic_checks(threads = max(threads_grid))
write_csv(semantic_checks, file.path(out_dir, "semantic_checks.csv"))

range_data <- make_reviewer_ranges()
gr1 <- range_data$query
gr2 <- range_data$subject

prep_bench <- rbind(
  time_object_build("fastRanges_build_index", function() fastRanges::fast_build_index(gr2), iters = iters),
  time_object_build("GenomicRanges_GNCList_build", function() GenomicRanges::GNCList(gr2), iters = iters)
)
write_csv(prep_bench, file.path(out_dir, "preprocessing_raw.csv"))
write_csv(summarise_runtime(prep_bench, baseline_engine = "GenomicRanges_GNCList_build"), file.path(out_dir, "preprocessing_summary.csv"))

idx <- fastRanges::fast_build_index(gr2)
gnc <- GenomicRanges::GNCList(gr2)

reviewer_workload <- data.frame(
  query_n = length(gr1),
  subject_n = length(gr2),
  query_mean_width = mean(IRanges::width(gr1)),
  subject_mean_width = mean(IRanges::width(gr2)),
  query_seqlevels = length(unique(as.character(GenomicRanges::seqnames(gr1)))),
  subject_seqlevels = length(unique(as.character(GenomicRanges::seqnames(gr2)))),
  stringsAsFactors = FALSE
)
write_csv(reviewer_workload, file.path(out_dir, "reviewer_workload.csv"))

perf_runs <- list()

perf_runs[[length(perf_runs) + 1L]] <- time_engine(
  "fastRanges_default",
  function() fastRanges::fast_find_overlaps(gr1, gr2),
  iters = iters
)

for (th in threads_grid) {
  for (det in c(TRUE, FALSE)) {
    det_label <- if (det) "det_true" else "det_false"
    perf_runs[[length(perf_runs) + 1L]] <- time_engine(
      paste0("fastRanges_direct_t", th, "_", det_label),
      function() fastRanges::fast_find_overlaps(gr1, gr2, threads = th, deterministic = det),
      iters = iters
    )
    perf_runs[[length(perf_runs) + 1L]] <- time_engine(
      paste0("fastRanges_indexed_t", th, "_", det_label),
      function() fastRanges::fast_find_overlaps(gr1, idx, threads = th, deterministic = det),
      iters = iters
    )
  }
}

perf_runs[[length(perf_runs) + 1L]] <- time_engine(
  "GenomicRanges_findOverlaps",
  function() GenomicRanges::findOverlaps(gr1, gr2),
  iters = iters
)

perf_runs[[length(perf_runs) + 1L]] <- time_engine(
  "GenomicRanges_GNCList",
  function() GenomicRanges::findOverlaps(gr1, gnc),
  iters = iters
)

performance_raw <- do.call(rbind, perf_runs)
write_csv(performance_raw, file.path(out_dir, "performance_raw.csv"))

if (length(unique(performance_raw$hits)) != 1L) {
  warning("Hit count mismatch detected across reviewer benchmark engines")
}

performance_summary <- summarise_runtime(
  performance_raw,
  baseline_engine = "GenomicRanges_findOverlaps"
)
write_csv(performance_summary, file.path(out_dir, "performance_summary.csv"))

determinism_focus <- subset(
  performance_summary,
  grepl("^fastRanges_(direct|indexed)_t", engine) | engine %in% c("GenomicRanges_findOverlaps", "GenomicRanges_GNCList", "fastRanges_default")
)
write_csv(determinism_focus, file.path(out_dir, "determinism_thread_summary.csv"))

if (run_peak_ram) {
  peak_ram <- run_peak_ram_checks(gr1, gr2, idx, threads_max = max(threads_grid))
  write_csv(peak_ram, file.path(out_dir, "peak_ram_summary.csv"))
} else {
  peak_ram <- data.frame(
    engine = character(),
    peakram_available = logical(),
    peak_ram_mb = numeric(),
    elapsed_sec = numeric(),
    stringsAsFactors = FALSE
  )
  write_csv(peak_ram, file.path(out_dir, "peak_ram_summary.csv"))
}

manifest <- data.frame(
  file = c(
    "run_metadata.csv",
    "package_versions.csv",
    "semantic_checks.csv",
    "preprocessing_raw.csv",
    "preprocessing_summary.csv",
    "reviewer_workload.csv",
    "performance_raw.csv",
    "performance_summary.csv",
    "determinism_thread_summary.csv",
    "peak_ram_summary.csv"
  ),
  description = c(
    "Benchmark run metadata",
    "Package versions used in the run",
    "Reviewer concern compatibility checks",
    "Index / preprocessing benchmark raw timings",
    "Median preprocessing summary and speedup vs GNCList build",
    "Reviewer workload dimensions",
    "Per-iteration benchmark timings on the reviewer workload",
    "Median runtime summary and speedup vs GenomicRanges::findOverlaps()",
    "Thread and deterministic comparison summary for fastRanges",
    "Optional peak RAM summary if peakRAM is installed"
  ),
  stringsAsFactors = FALSE
)
write_csv(manifest, file.path(out_dir, "manifest.csv"))

cat("Reviewer response checks completed\n")
cat("Results written to: ", out_dir, "\n", sep = "")
