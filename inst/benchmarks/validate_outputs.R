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

canon_hits <- function(h) {
  x <- cbind(S4Vectors::queryHits(h), S4Vectors::subjectHits(h))
  if (nrow(x) == 0L) return(x)
  x[order(x[, 1], x[, 2]), , drop = FALSE]
}

validate_select_result <- function(selected, ref_all_hits, query_n, mode = c("first", "last", "arbitrary")) {
  mode <- match.arg(mode)
  if (!is.integer(selected)) {
    return(FALSE)
  }
  if (length(selected) != query_n) {
    return(FALSE)
  }

  qh <- S4Vectors::queryHits(ref_all_hits)
  sh <- S4Vectors::subjectHits(ref_all_hits)
  ref_choices <- split(sh, qh)
  has_hit <- rep.int(FALSE, query_n)
  if (length(ref_choices) > 0L) {
    has_hit[as.integer(names(ref_choices))] <- TRUE
  }

  if (!identical(is.na(selected), !has_hit)) {
    return(FALSE)
  }
  if (!any(has_hit)) {
    return(TRUE)
  }

  for (qid in which(has_hit)) {
    choices <- ref_choices[[as.character(qid)]]
    if (is.null(choices) || length(choices) == 0L) {
      return(FALSE)
    }
    got <- selected[[qid]]
    if (is.na(got)) {
      return(FALSE)
    }
    if (identical(mode, "first") && !identical(as.integer(got), as.integer(min(choices)))) {
      return(FALSE)
    }
    if (identical(mode, "last") && !identical(as.integer(got), as.integer(max(choices)))) {
      return(FALSE)
    }
    if (identical(mode, "arbitrary") && !(got %in% choices)) {
      return(FALSE)
    }
  }

  TRUE
}

same_partition <- function(x, y) {
  if (length(x) != length(y)) {
    return(FALSE)
  }
  if (length(x) == 0L) {
    return(TRUE)
  }
  identical(outer(x, x, `==`), outer(y, y, `==`))
}

strip_rownames <- function(x) {
  if (is.data.frame(x)) {
    rownames(x) <- NULL
  }
  x
}

equal_df <- function(x, y, sort_by = NULL, tolerance = 1e-10) {
  x <- strip_rownames(as.data.frame(x, stringsAsFactors = FALSE))
  y <- strip_rownames(as.data.frame(y, stringsAsFactors = FALSE))
  if (!is.null(sort_by) && length(sort_by) > 0L) {
    sort_by <- intersect(sort_by, intersect(names(x), names(y)))
    if (length(sort_by) > 0L) {
      ox <- do.call(order, x[sort_by])
      oy <- do.call(order, y[sort_by])
      x <- x[ox, , drop = FALSE]
      y <- y[oy, , drop = FALSE]
    }
  }
  isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE))
}

manual_group_counts <- function(query, subject, group_col, include_na_group = FALSE, hits = NULL) {
  groups <- as.character(S4Vectors::mcols(subject)[[group_col]])
  if (include_na_group) {
    groups[is.na(groups)] <- "<NA>"
  }
  valid <- !is.na(groups)
  levels <- sort(unique(groups[valid]))
  out <- matrix(0L, nrow = length(query), ncol = length(levels))
  colnames(out) <- levels
  if (length(levels) == 0L) return(out)

  if (is.null(hits)) {
    hits <- GenomicRanges::findOverlaps(query, subject, ignore.strand = FALSE)
  }
  if (length(hits) == 0L) return(out)

  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)
  g <- groups[sh]
  keep <- !is.na(g)
  if (!any(keep)) return(out)

  qh <- qh[keep]
  g <- g[keep]
  tab <- table(qh, factor(g, levels = levels))
  out[as.integer(rownames(tab)), ] <- tab
  out
}

manual_overlap_aggregate <- function(
    query,
    subject,
    value_col,
    fun,
    na_rm = TRUE,
    hits = NULL,
    ignore_strand = FALSE) {
  if (fun == "count") {
    if (!is.null(hits)) {
      return(as.numeric(tabulate(S4Vectors::queryHits(hits), nbins = length(query))))
    }
    return(as.numeric(GenomicRanges::countOverlaps(query, subject, ignore.strand = ignore_strand)))
  }
  values <- as.numeric(S4Vectors::mcols(subject)[[value_col]])
  out <- rep(NA_real_, length(query))

  if (is.null(hits)) {
    hits <- GenomicRanges::findOverlaps(query, subject, ignore.strand = ignore_strand)
  }
  if (length(hits) == 0L) return(out)

  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)
  x <- values[sh]
  if (na_rm) {
    keep <- !is.na(x)
    qh <- qh[keep]
    x <- x[keep]
  }
  if (length(x) == 0L) return(out)

  split_x <- split(x, qh)
  agg_fun <- switch(
    fun,
    sum = function(z) sum(z, na.rm = na_rm),
    mean = function(z) mean(z, na.rm = na_rm),
    min = function(z) min(z, na.rm = na_rm),
    max = function(z) max(z, na.rm = na_rm)
  )
  agg <- vapply(split_x, agg_fun, numeric(1))
  out[as.integer(names(agg))] <- agg
  out
}

baseline_window_count <- function(query, subject, window_width, step_width, ignore_strand = FALSE) {
  seq_chr <- as.character(GenomicRanges::seqnames(query))
  split_idx <- split(seq_len(length(query)), seq_chr)

  pieces <- vector("list", length(split_idx))
  for (i in seq_along(split_idx)) {
    idx <- split_idx[[i]]
    s <- min(IRanges::start(query[idx]))
    e <- max(IRanges::end(query[idx]))
    starts <- seq.int(s, e, by = as.integer(step_width))
    pieces[[i]] <- GenomicRanges::GRanges(
      seqnames = names(split_idx)[i],
      ranges = IRanges::IRanges(start = starts, width = as.integer(window_width)),
      strand = "*"
    )
  }
  windows <- do.call(c, pieces)
  counts <- GenomicRanges::countOverlaps(windows, subject, ignore.strand = ignore_strand)
  data.frame(
    seqnames = as.character(GenomicRanges::seqnames(windows)),
    start = IRanges::start(windows),
    end = IRanges::end(windows),
    overlap_count = as.integer(counts),
    stringsAsFactors = FALSE
  )
}

baseline_tile_coverage <- function(
    x,
    tile_width,
    step_width = tile_width,
    shift = 0L,
    width = NULL,
    weight = 1L,
    method = c("auto", "sort", "hash")) {
  method <- match.arg(method)
  cov <- if (methods::is(x, "GenomicRanges")) {
    GenomicRanges::coverage(x, shift = shift, width = width, weight = weight, method = method)
  } else {
    IRanges::coverage(x, shift = shift, width = width, weight = weight, method = method)
  }

  if (methods::is(x, "GenomicRanges")) {
    cov_names <- names(cov)
    if (is.null(cov_names)) {
      cov_names <- paste0("seq", seq_along(cov))
    }

    out <- vector("list", length(cov))
    for (i in seq_along(cov)) {
      n <- length(cov[[i]])
      starts <- seq.int(1L, n, by = as.integer(step_width))
      ends <- pmin(starts + as.integer(tile_width) - 1L, n)
      sums <- as.numeric(IRanges::viewSums(IRanges::Views(cov[[i]], start = starts, end = ends)))
      out[[i]] <- data.frame(
        seqnames = cov_names[i],
        start = starts,
        end = ends,
        coverage_sum = sums,
        stringsAsFactors = FALSE
      )
    }
    return(do.call(rbind, out))
  }

  n <- length(cov)
  starts <- seq.int(1L, n, by = as.integer(step_width))
  ends <- pmin(starts + as.integer(tile_width) - 1L, n)
  sums <- as.numeric(IRanges::viewSums(IRanges::Views(cov, start = starts, end = ends)))
  data.frame(
    start = starts,
    end = ends,
    coverage_sum = sums,
    stringsAsFactors = FALSE
  )
}

cluster_from_hits <- function(n, hits) {
  if (n == 0L) {
    return(integer())
  }

  parent <- seq_len(n)
  rank <- integer(n)

  find_root <- function(i) {
    while (parent[i] != i) {
      parent[i] <<- parent[parent[i]]
      i <- parent[i]
    }
    i
  }

  union_sets <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (ra == rb) {
      return(invisible(NULL))
    }
    if (rank[ra] < rank[rb]) {
      parent[ra] <<- rb
    } else if (rank[ra] > rank[rb]) {
      parent[rb] <<- ra
    } else {
      parent[rb] <<- ra
      rank[ra] <<- rank[ra] + 1L
    }
    invisible(NULL)
  }

  if (length(hits) > 0L) {
    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)
    for (i in seq_along(qh)) {
      union_sets(qh[i], sh[i])
    }
  }

  roots <- vapply(seq_len(n), find_root, integer(1))
  match(roots, unique(roots))
}

baseline_overlap_join <- function(
    query,
    subject,
    join = c("inner", "left"),
    hits,
    query_prefix = "query_",
    subject_prefix = "subject_") {
  join <- match.arg(join)
  query_df <- as.data.frame(query)
  subject_df <- as.data.frame(subject)
  names(query_df) <- paste0(query_prefix, names(query_df))
  names(subject_df) <- paste0(subject_prefix, names(subject_df))

  q_idx <- S4Vectors::queryHits(hits)
  s_idx <- S4Vectors::subjectHits(hits)

  if (join == "inner") {
    if (length(q_idx) == 0L) {
      empty <- data.frame(query_id = integer(), subject_id = integer())
      return(cbind(empty, query_df[FALSE, , drop = FALSE], subject_df[FALSE, , drop = FALSE]))
    }
    out <- data.frame(query_id = q_idx, subject_id = s_idx)
    return(cbind(out, query_df[q_idx, , drop = FALSE], subject_df[s_idx, , drop = FALSE]))
  }

  parts <- vector("list", length(query))
  for (i in seq_len(length(query))) {
    sh_i <- s_idx[q_idx == i]
    if (length(sh_i) == 0L) {
      parts[[i]] <- cbind(
        data.frame(query_id = i, subject_id = NA_integer_),
        query_df[i, , drop = FALSE],
        subject_df[NA_integer_, , drop = FALSE]
      )
    } else {
      parts[[i]] <- cbind(
        data.frame(query_id = rep.int(i, length(sh_i)), subject_id = sh_i),
        query_df[rep.int(i, length(sh_i)), , drop = FALSE],
        subject_df[sh_i, , drop = FALSE]
      )
    }
  }
  do.call(rbind, parts)
}

baseline_semi_join <- function(query, counts, query_prefix = "query_") {
  keep <- counts > 0L
  out <- as.data.frame(query)
  names(out) <- paste0(query_prefix, names(out))
  out <- out[keep, , drop = FALSE]
  out$query_id <- which(keep)
  out$overlap_count <- counts[keep]
  out[, c("query_id", "overlap_count", setdiff(names(out), c("query_id", "overlap_count"))), drop = FALSE]
}

baseline_anti_join <- function(query, counts, query_prefix = "query_") {
  keep <- counts == 0L
  out <- as.data.frame(query)
  names(out) <- paste0(query_prefix, names(out))
  out <- out[keep, , drop = FALSE]
  out$query_id <- which(keep)
  out$overlap_count <- counts[keep]
  out[, c("query_id", "overlap_count", setdiff(names(out), c("query_id", "overlap_count"))), drop = FALSE]
}

make_granges <- function(n, seed, seqlevels = c("chr1", "chr2", "chr3")) {
  set.seed(seed)
  starts <- sample.int(50000L, n, replace = TRUE)
  widths <- sample.int(400L, n, replace = TRUE) + 30L
  GenomicRanges::GRanges(
    seqnames = sample(seqlevels, n, replace = TRUE),
    ranges = IRanges::IRanges(start = starts, width = widths),
    strand = sample(c("+", "-", "*"), n, replace = TRUE, prob = c(0.45, 0.45, 0.10))
  )
}

make_iranges <- function(n, seed) {
  set.seed(seed)
  starts <- sample.int(50000L, n, replace = TRUE)
  widths <- sample.int(400L, n, replace = TRUE) + 30L
  IRanges::IRanges(start = starts, width = widths)
}

args <- commandArgs(trailingOnly = TRUE)
n_iter <- as.integer(get_arg(args, "iters", Sys.getenv("FASTRANGES_VALIDATE_ITERS", "5")))
seed <- as.integer(get_arg(args, "seed", Sys.getenv("FASTRANGES_VALIDATE_SEED", "42")))
threads <- as.integer(get_arg(args, "threads", Sys.getenv("FASTRANGES_VALIDATE_THREADS", "8")))
q_gr_n <- as.integer(get_arg(args, "q_gr_n", Sys.getenv("FASTRANGES_VALIDATE_Q_GR_N", "250")))
s_gr_n <- as.integer(get_arg(args, "s_gr_n", Sys.getenv("FASTRANGES_VALIDATE_S_GR_N", "900")))
q_ir_n <- as.integer(get_arg(args, "q_ir_n", Sys.getenv("FASTRANGES_VALIDATE_Q_IR_N", "220")))
s_ir_n <- as.integer(get_arg(args, "s_ir_n", Sys.getenv("FASTRANGES_VALIDATE_S_IR_N", "880")))
threads <- max(1L, threads)
q_gr_n <- max(100L, q_gr_n)
s_gr_n <- max(100L, s_gr_n)
q_ir_n <- max(100L, q_ir_n)
s_ir_n <- max(100L, s_ir_n)

set.seed(seed)

message("fastRanges validation started")
message("iters=", n_iter, " seed=", seed, " threads=", threads)
message("GRanges sizes: query=", q_gr_n, " subject=", s_gr_n)
message("IRanges sizes: query=", q_ir_n, " subject=", s_ir_n)

types <- c("any", "start", "end", "within", "equal")
max_gaps <- c(-1L, 0L, 5L)
min_overlaps <- c(0L, 1L)

# Keep the validator aligned with Bioconductor's accepted parameter grammar.
# For type = "any", at least one of max_gap and min_overlap must be default.
is_valid_bioc_combo <- function(type, max_gap, min_overlap) {
  if (identical(type, "any") && (max_gap != -1L) && (min_overlap != 0L)) {
    return(FALSE)
  }
  TRUE
}

# Basic exported helper and index API checks.
old_threads <- getOption("fastRanges.threads")
on.exit(options(fastRanges.threads = old_threads), add = TRUE)
options(fastRanges.threads = 3L)
if (!identical(fastRanges::fast_default_threads(), 3L)) {
  stop("Mismatch: fast_default_threads option handling")
}
options(fastRanges.threads = old_threads)

idx_subject <- make_granges(120L, seed + 5000L)
idx <- fastRanges::fast_build_index(idx_subject)
if (!inherits(idx, "fast_ranges_index") || !identical(as.integer(idx$subject_n), as.integer(length(idx_subject)))) {
  stop("Mismatch: fast_build_index structure")
}

stats_basic <- fastRanges::fast_index_stats(idx)
if (!all(c("subject_n", "partition_n", "seqlevel_n", "mean_partition_size", "index_size_mb") %in% names(stats_basic))) {
  stop("Mismatch: fast_index_stats summary columns")
}
if (!identical(as.integer(stats_basic$subject_n), as.integer(length(idx_subject)))) {
  stop("Mismatch: fast_index_stats subject_n")
}

stats_detailed <- fastRanges::fast_index_stats(idx, detailed = TRUE)
if (!all(c("summary", "partitions") %in% names(stats_detailed))) {
  stop("Mismatch: fast_index_stats detailed structure")
}

idx_path <- tempfile(fileext = ".rds")
on.exit(unlink(idx_path), add = TRUE)
fastRanges::fast_save_index(idx, idx_path)
idx_loaded <- fastRanges::fast_load_index(idx_path)
if (!isTRUE(all.equal(idx, idx_loaded, check.attributes = FALSE))) {
  stop("Mismatch: fast_save_index/fast_load_index roundtrip")
}

# Range algebra and coverage wrappers against Bioconductor references.
xi <- IRanges::IRanges(start = c(1L, 4L, 10L, 20L), end = c(5L, 8L, 12L, 24L))
yi <- IRanges::IRanges(start = c(3L, 7L, 22L), end = c(6L, 14L, 26L))
if (!isTRUE(all.equal(fastRanges::fast_reduce(xi), IRanges::reduce(xi), check.attributes = FALSE))) {
  stop("Mismatch: fast_reduce IRanges")
}
if (!isTRUE(all.equal(fastRanges::fast_disjoin(xi), IRanges::disjoin(xi), check.attributes = FALSE))) {
  stop("Mismatch: fast_disjoin IRanges")
}
if (!isTRUE(all.equal(fastRanges::fast_gaps(xi, start = 1L, end = 30L), IRanges::gaps(xi, start = 1L, end = 30L), check.attributes = FALSE))) {
  stop("Mismatch: fast_gaps IRanges")
}
if (!isTRUE(all.equal(fastRanges::fast_range_union(xi, yi), IRanges::union(xi, yi), check.attributes = FALSE))) {
  stop("Mismatch: fast_range_union IRanges")
}
if (!isTRUE(all.equal(fastRanges::fast_range_intersect(xi, yi), IRanges::intersect(xi, yi), check.attributes = FALSE))) {
  stop("Mismatch: fast_range_intersect IRanges")
}
if (!isTRUE(all.equal(fastRanges::fast_range_setdiff(xi, yi), IRanges::setdiff(xi, yi), check.attributes = FALSE))) {
  stop("Mismatch: fast_range_setdiff IRanges")
}

xg <- GenomicRanges::GRanges(
  seqnames = c("chr1", "chr1", "chr2", "chr2"),
  ranges = IRanges::IRanges(start = c(1L, 4L, 10L, 20L), end = c(5L, 8L, 12L, 24L)),
  strand = c("+", "+", "-", "*")
)
yg <- GenomicRanges::GRanges(
  seqnames = c("chr1", "chr2", "chr2"),
  ranges = IRanges::IRanges(start = c(3L, 7L, 22L), end = c(6L, 14L, 26L)),
  strand = c("+", "-", "*")
)
seqlengths(xg) <- c(chr1 = 30L, chr2 = 30L)
seqlengths(yg) <- c(chr1 = 30L, chr2 = 30L)
if (!isTRUE(all.equal(
  fastRanges::fast_reduce(xg, ignore_strand = TRUE),
  GenomicRanges::reduce(xg, ignore.strand = TRUE),
  check.attributes = FALSE
))) {
  stop("Mismatch: fast_reduce GRanges")
}
if (!isTRUE(all.equal(
  fastRanges::fast_disjoin(xg, ignore_strand = TRUE),
  GenomicRanges::disjoin(xg, ignore.strand = TRUE),
  check.attributes = FALSE
))) {
  stop("Mismatch: fast_disjoin GRanges")
}
if (!isTRUE(all.equal(
  fastRanges::fast_gaps(xg, ignore_strand = TRUE),
  GenomicRanges::gaps(xg, ignore.strand = TRUE),
  check.attributes = FALSE
))) {
  stop("Mismatch: fast_gaps GRanges")
}
if (!isTRUE(all.equal(
  suppressWarnings(fastRanges::fast_range_union(xg, yg, ignore_strand = TRUE)),
  suppressWarnings(GenomicRanges::union(xg, yg, ignore.strand = TRUE)),
  check.attributes = FALSE
))) {
  stop("Mismatch: fast_range_union GRanges")
}
if (!isTRUE(all.equal(
  suppressWarnings(fastRanges::fast_range_intersect(xg, yg, ignore_strand = TRUE)),
  suppressWarnings(GenomicRanges::intersect(xg, yg, ignore.strand = TRUE)),
  check.attributes = FALSE
))) {
  stop("Mismatch: fast_range_intersect GRanges")
}
if (!isTRUE(all.equal(
  suppressWarnings(fastRanges::fast_range_setdiff(xg, yg, ignore_strand = TRUE)),
  suppressWarnings(GenomicRanges::setdiff(xg, yg, ignore.strand = TRUE)),
  check.attributes = FALSE
))) {
  stop("Mismatch: fast_range_setdiff GRanges")
}

cov_i_got <- fastRanges::fast_coverage(xi, method = "sort", threads = threads)
cov_i_ref <- IRanges::coverage(xi, method = "sort")
if (!isTRUE(all.equal(cov_i_got, cov_i_ref, check.attributes = FALSE))) {
  stop("Mismatch: fast_coverage IRanges")
}

tile_i_got <- fastRanges::fast_tile_coverage(xi, tile_width = 5L, step_width = 3L, method = "sort", threads = threads)
tile_i_ref <- baseline_tile_coverage(xi, tile_width = 5L, step_width = 3L, method = "sort")
if (!equal_df(tile_i_got, tile_i_ref, sort_by = c("start", "end"))) {
  stop("Mismatch: fast_tile_coverage IRanges")
}

dn_i_got <- fastRanges::fast_distance_to_nearest(xi, yi, threads = threads)
dn_i_ref <- IRanges::distanceToNearest(xi, yi)
if (!identical(as.integer(dn_i_got$query_id), as.integer(S4Vectors::queryHits(dn_i_ref))) ||
    !identical(as.integer(dn_i_got$subject_id), as.integer(S4Vectors::subjectHits(dn_i_ref))) ||
    !isTRUE(all.equal(as.numeric(dn_i_got$distance), as.numeric(S4Vectors::mcols(dn_i_ref)$distance), check.attributes = FALSE))) {
  stop("Mismatch: fast_distance_to_nearest IRanges")
}
if (!isTRUE(all.equal(fastRanges::fast_nearest(xi, yi, threads = threads), dn_i_got, check.attributes = FALSE))) {
  stop("Mismatch: fast_nearest IRanges")
}
if (!identical(fastRanges::fast_precede(xi, yi, threads = threads), IRanges::precede(xi, yi))) {
  stop("Mismatch: fast_precede IRanges")
}
if (!identical(fastRanges::fast_follow(xi, yi, threads = threads), IRanges::follow(xi, yi))) {
  stop("Mismatch: fast_follow IRanges")
}

q0g <- GenomicRanges::GRanges(
  "chr1",
  IRanges::IRanges(start = c(5L, 8L, 12L, 20L), end = c(4L, 10L, 11L, 23L))
)
s0g <- GenomicRanges::GRanges(
  "chr1",
  IRanges::IRanges(start = c(3L, 5L, 8L, 12L, 19L), end = c(6L, 4L, 10L, 11L, 21L))
)
idx0g <- fastRanges::fast_build_index(s0g)
ref0g <- GenomicRanges::findOverlaps(q0g, s0g)
got0g_direct <- fastRanges::fast_find_overlaps(q0g, s0g, threads = threads)
got0g_index <- fastRanges::fast_find_overlaps(q0g, idx0g, threads = threads)
if (!isTRUE(all.equal(canon_hits(got0g_direct), canon_hits(ref0g), check.attributes = FALSE))) {
  stop("Mismatch: fast_find_overlaps direct empty-width GRanges")
}
if (!isTRUE(all.equal(canon_hits(got0g_index), canon_hits(ref0g), check.attributes = FALSE))) {
  stop("Mismatch: fast_find_overlaps indexed empty-width GRanges")
}
for (sel in c("first", "last")) {
  ref0g_sel <- GenomicRanges::findOverlaps(q0g, s0g, select = sel)
  got0g_sel_direct <- fastRanges::fast_find_overlaps(q0g, s0g, select = sel, threads = threads)
  got0g_sel_index <- fastRanges::fast_find_overlaps(q0g, idx0g, select = sel, threads = threads)
  if (!identical(got0g_sel_direct, ref0g_sel)) {
    stop("Mismatch: fast_find_overlaps direct empty-width GRanges select=", sel)
  }
  if (!identical(got0g_sel_index, ref0g_sel)) {
    stop("Mismatch: fast_find_overlaps indexed empty-width GRanges select=", sel)
  }
}
if (!validate_select_result(fastRanges::fast_find_overlaps(q0g, s0g, select = "arbitrary", threads = threads), ref0g, length(q0g), mode = "arbitrary")) {
  stop("Mismatch: fast_find_overlaps direct empty-width GRanges select=arbitrary")
}
if (!validate_select_result(fastRanges::fast_find_overlaps(q0g, idx0g, select = "arbitrary", threads = threads), ref0g, length(q0g), mode = "arbitrary")) {
  stop("Mismatch: fast_find_overlaps indexed empty-width GRanges select=arbitrary")
}
ref0g_count <- GenomicRanges::countOverlaps(q0g, s0g)
if (!identical(as.integer(fastRanges::fast_count_overlaps(q0g, s0g, threads = threads)), as.integer(ref0g_count))) {
  stop("Mismatch: fast_count_overlaps direct empty-width GRanges")
}
if (!identical(as.integer(fastRanges::fast_count_overlaps(q0g, idx0g, threads = threads)), as.integer(ref0g_count))) {
  stop("Mismatch: fast_count_overlaps indexed empty-width GRanges")
}
if (!identical(as.logical(fastRanges::fast_overlaps_any(q0g, s0g, threads = threads)), ref0g_count > 0L)) {
  stop("Mismatch: fast_overlaps_any direct empty-width GRanges")
}
if (!identical(as.logical(fastRanges::fast_overlaps_any(q0g, idx0g, threads = threads)), ref0g_count > 0L)) {
  stop("Mismatch: fast_overlaps_any indexed empty-width GRanges")
}

q0i <- IRanges::IRanges(start = c(5L, 8L, 12L, 20L), end = c(4L, 10L, 11L, 23L))
s0i <- IRanges::IRanges(start = c(3L, 5L, 8L, 12L, 19L), end = c(6L, 4L, 10L, 11L, 21L))
idx0i <- fastRanges::fast_build_index(s0i)
ref0i <- IRanges::findOverlaps(q0i, s0i)
got0i_direct <- fastRanges::fast_find_overlaps(q0i, s0i, threads = threads)
got0i_index <- fastRanges::fast_find_overlaps(q0i, idx0i, threads = threads)
if (!isTRUE(all.equal(canon_hits(got0i_direct), canon_hits(ref0i), check.attributes = FALSE))) {
  stop("Mismatch: fast_find_overlaps direct empty-width IRanges")
}
if (!isTRUE(all.equal(canon_hits(got0i_index), canon_hits(ref0i), check.attributes = FALSE))) {
  stop("Mismatch: fast_find_overlaps indexed empty-width IRanges")
}
for (sel in c("first", "last")) {
  ref0i_sel <- IRanges::findOverlaps(q0i, s0i, select = sel)
  got0i_sel_direct <- fastRanges::fast_find_overlaps(q0i, s0i, select = sel, threads = threads)
  got0i_sel_index <- fastRanges::fast_find_overlaps(q0i, idx0i, select = sel, threads = threads)
  if (!identical(got0i_sel_direct, ref0i_sel)) {
    stop("Mismatch: fast_find_overlaps direct empty-width IRanges select=", sel)
  }
  if (!identical(got0i_sel_index, ref0i_sel)) {
    stop("Mismatch: fast_find_overlaps indexed empty-width IRanges select=", sel)
  }
}
if (!validate_select_result(fastRanges::fast_find_overlaps(q0i, s0i, select = "arbitrary", threads = threads), ref0i, length(q0i), mode = "arbitrary")) {
  stop("Mismatch: fast_find_overlaps direct empty-width IRanges select=arbitrary")
}
if (!validate_select_result(fastRanges::fast_find_overlaps(q0i, idx0i, select = "arbitrary", threads = threads), ref0i, length(q0i), mode = "arbitrary")) {
  stop("Mismatch: fast_find_overlaps indexed empty-width IRanges select=arbitrary")
}
ref0i_count <- IRanges::countOverlaps(q0i, s0i)
if (!identical(as.integer(fastRanges::fast_count_overlaps(q0i, s0i, threads = threads)), as.integer(ref0i_count))) {
  stop("Mismatch: fast_count_overlaps direct empty-width IRanges")
}
if (!identical(as.integer(fastRanges::fast_count_overlaps(q0i, idx0i, threads = threads)), as.integer(ref0i_count))) {
  stop("Mismatch: fast_count_overlaps indexed empty-width IRanges")
}
if (!identical(as.logical(fastRanges::fast_overlaps_any(q0i, s0i, threads = threads)), ref0i_count > 0L)) {
  stop("Mismatch: fast_overlaps_any direct empty-width IRanges")
}
if (!identical(as.logical(fastRanges::fast_overlaps_any(q0i, idx0i, threads = threads)), ref0i_count > 0L)) {
  stop("Mismatch: fast_overlaps_any indexed empty-width IRanges")
}

for (iter in seq_len(n_iter)) {
  qg <- make_granges(q_gr_n, seed + iter * 10L)
  sg <- make_granges(s_gr_n, seed + iter * 10L + 1L)
  S4Vectors::mcols(sg)$grp <- sample(c("A", "B", "C", NA_character_), length(sg), replace = TRUE)
  S4Vectors::mcols(sg)$score <- stats::rnorm(length(sg), mean = 2, sd = 1)

  idx <- fastRanges::fast_build_index(sg)

  for (ignore_strand in c(FALSE, TRUE)) {
    for (type in types) {
      for (max_gap in max_gaps) {
        for (min_overlap in min_overlaps) {
          if (!is_valid_bioc_combo(type, max_gap, min_overlap)) {
            next
          }

          ref <- GenomicRanges::findOverlaps(
            qg, sg,
            maxgap = max_gap,
            minoverlap = min_overlap,
            type = type,
            ignore.strand = ignore_strand
          )
          got_direct <- fastRanges::fast_find_overlaps(
            qg, sg,
            max_gap = max_gap,
            min_overlap = min_overlap,
            type = type,
            ignore_strand = ignore_strand,
            threads = threads
          )
          got_index <- fastRanges::fast_find_overlaps(
            qg, idx,
            max_gap = max_gap,
            min_overlap = min_overlap,
            type = type,
            ignore_strand = ignore_strand,
            threads = threads
          )
          if (!isTRUE(all.equal(canon_hits(got_direct), canon_hits(ref), check.attributes = FALSE))) {
            stop("Mismatch: fast_find_overlaps direct vs GenomicRanges at iter=", iter,
                 " ignore_strand=", ignore_strand,
                 " type=", type, " max_gap=", max_gap, " min_overlap=", min_overlap)
          }
          if (!isTRUE(all.equal(canon_hits(got_index), canon_hits(ref), check.attributes = FALSE))) {
            stop("Mismatch: fast_find_overlaps index vs GenomicRanges at iter=", iter,
                 " ignore_strand=", ignore_strand,
                 " type=", type, " max_gap=", max_gap, " min_overlap=", min_overlap)
          }

          for (sel in c("first", "last")) {
            ref_sel <- GenomicRanges::findOverlaps(
              qg, sg,
              maxgap = max_gap,
              minoverlap = min_overlap,
              type = type,
              ignore.strand = ignore_strand,
              select = sel
            )
            got_sel_direct <- fastRanges::fast_find_overlaps(
              qg, sg,
              select = sel,
              max_gap = max_gap,
              min_overlap = min_overlap,
              type = type,
              ignore_strand = ignore_strand,
              threads = threads
            )
            got_sel_index <- fastRanges::fast_find_overlaps(
              qg, idx,
              select = sel,
              max_gap = max_gap,
              min_overlap = min_overlap,
              type = type,
              ignore_strand = ignore_strand,
              threads = threads
            )
            if (!identical(got_sel_direct, ref_sel)) {
              stop("Mismatch: fast_find_overlaps select=", sel, " direct vs GenomicRanges at iter=", iter,
                   " ignore_strand=", ignore_strand,
                   " type=", type, " max_gap=", max_gap, " min_overlap=", min_overlap)
            }
            if (!identical(got_sel_index, ref_sel)) {
              stop("Mismatch: fast_find_overlaps select=", sel, " index vs GenomicRanges at iter=", iter,
                   " ignore_strand=", ignore_strand,
                   " type=", type, " max_gap=", max_gap, " min_overlap=", min_overlap)
            }
          }

          got_sel_arb_direct <- fastRanges::fast_find_overlaps(
            qg, sg,
            select = "arbitrary",
            max_gap = max_gap,
            min_overlap = min_overlap,
            type = type,
            ignore_strand = ignore_strand,
            threads = threads
          )
          got_sel_arb_index <- fastRanges::fast_find_overlaps(
            qg, idx,
            select = "arbitrary",
            max_gap = max_gap,
            min_overlap = min_overlap,
            type = type,
            ignore_strand = ignore_strand,
            threads = threads
          )
          if (!validate_select_result(got_sel_arb_direct, ref, length(qg), mode = "arbitrary")) {
            stop("Mismatch: fast_find_overlaps select=arbitrary direct semantic validity at iter=", iter,
                 " ignore_strand=", ignore_strand,
                 " type=", type, " max_gap=", max_gap, " min_overlap=", min_overlap)
          }
          if (!validate_select_result(got_sel_arb_index, ref, length(qg), mode = "arbitrary")) {
            stop("Mismatch: fast_find_overlaps select=arbitrary index semantic validity at iter=", iter,
                 " ignore_strand=", ignore_strand,
                 " type=", type, " max_gap=", max_gap, " min_overlap=", min_overlap)
          }

          cnt <- fastRanges::fast_count_overlaps(
            qg, idx,
            max_gap = max_gap,
            min_overlap = min_overlap,
            type = type,
            ignore_strand = ignore_strand,
            threads = threads
          )
          ref_cnt <- as.integer(GenomicRanges::countOverlaps(
            qg, sg,
            maxgap = max_gap,
            minoverlap = min_overlap,
            type = type,
            ignore.strand = ignore_strand
          ))
          if (!identical(as.integer(cnt), ref_cnt)) {
            stop("Mismatch: fast_count_overlaps at iter=", iter,
                 " ignore_strand=", ignore_strand,
                 " type=", type, " max_gap=", max_gap, " min_overlap=", min_overlap)
          }

          any_got <- fastRanges::fast_overlaps_any(
            qg, idx,
            max_gap = max_gap,
            min_overlap = min_overlap,
            type = type,
            ignore_strand = ignore_strand,
            threads = threads
          )
          if (!identical(as.logical(any_got), ref_cnt > 0L)) {
            stop("Mismatch: fast_overlaps_any at iter=", iter,
                 " ignore_strand=", ignore_strand,
                 " type=", type, " max_gap=", max_gap, " min_overlap=", min_overlap)
          }
        }
      }
    }
  }

  # Iterator API.
  iter_direct <- fastRanges::fast_find_overlaps_iter(
    qg, sg,
    chunk_size = 41L,
    threads = threads,
    deterministic = TRUE
  )
  iter_index <- fastRanges::fast_find_overlaps_iter(
    qg, idx,
    chunk_size = 41L,
    threads = threads,
    deterministic = TRUE
  )
  if (!fastRanges::fast_iter_has_next(iter_direct) || !fastRanges::fast_iter_has_next(iter_index)) {
    stop("Mismatch: fast_find_overlaps_iter initial state at iter=", iter)
  }
  iter_first <- fastRanges::fast_iter_next(iter_direct)
  iter_remaining <- fastRanges::fast_iter_collect(iter_direct)
  fastRanges::fast_iter_reset(iter_direct)
  iter_all <- fastRanges::fast_iter_collect(iter_direct)
  ref_all_default <- GenomicRanges::findOverlaps(qg, sg, ignore.strand = FALSE)
  if (!isTRUE(all.equal(canon_hits(iter_all), canon_hits(ref_all_default), check.attributes = FALSE))) {
    stop("Mismatch: fast_iter_collect direct at iter=", iter)
  }
  if (length(iter_first) == 0L && length(ref_all_default) > 0L) {
    stop("Mismatch: fast_iter_next empty first chunk at iter=", iter)
  }
  iter_combined <- S4Vectors::Hits(
    from = c(S4Vectors::queryHits(iter_first), S4Vectors::queryHits(iter_remaining)),
    to = c(S4Vectors::subjectHits(iter_first), S4Vectors::subjectHits(iter_remaining)),
    nLnode = length(qg),
    nRnode = length(sg),
    sort.by.query = FALSE
  )
  if (!isTRUE(all.equal(canon_hits(iter_combined), canon_hits(ref_all_default), check.attributes = FALSE))) {
    stop("Mismatch: fast_iter_collect after partial iteration at iter=", iter)
  }
  iter_all_index <- fastRanges::fast_iter_collect(iter_index)
  if (!isTRUE(all.equal(canon_hits(iter_all_index), canon_hits(ref_all_default), check.attributes = FALSE))) {
    stop("Mismatch: fast_iter_collect indexed at iter=", iter)
  }

  for (ignore_strand in c(FALSE, TRUE)) {
    # Join family.
    ref_join_hits <- GenomicRanges::findOverlaps(qg, sg, ignore.strand = ignore_strand)
    join_inner_ref <- baseline_overlap_join(qg, sg, join = "inner", hits = ref_join_hits)
    join_left_ref <- baseline_overlap_join(qg, sg, join = "left", hits = ref_join_hits)
    counts_default <- as.integer(GenomicRanges::countOverlaps(qg, sg, ignore.strand = ignore_strand))
    semi_ref <- baseline_semi_join(qg, counts_default)
    anti_ref <- baseline_anti_join(qg, counts_default)

    join_inner_got <- fastRanges::fast_overlap_join(
      qg, sg,
      join = "inner",
      ignore_strand = ignore_strand,
      threads = threads,
      deterministic = TRUE
    )
    join_left_got <- fastRanges::fast_overlap_join(
      qg, sg,
      join = "left",
      ignore_strand = ignore_strand,
      threads = threads,
      deterministic = TRUE
    )
    inner_wrap_got <- fastRanges::fast_inner_overlap_join(
      qg, sg,
      ignore_strand = ignore_strand,
      threads = threads,
      deterministic = TRUE
    )
    left_wrap_got <- fastRanges::fast_left_overlap_join(
      qg, sg,
      ignore_strand = ignore_strand,
      threads = threads,
      deterministic = TRUE
    )
    semi_got <- fastRanges::fast_semi_overlap_join(
      qg, sg,
      ignore_strand = ignore_strand,
      threads = threads,
      deterministic = TRUE
    )
    anti_got <- fastRanges::fast_anti_overlap_join(
      qg, sg,
      ignore_strand = ignore_strand,
      threads = threads,
      deterministic = TRUE
    )

    if (!equal_df(join_inner_got, join_inner_ref, sort_by = c("query_id", "subject_id"))) {
      stop("Mismatch: fast_overlap_join(inner) at iter=", iter, " ignore_strand=", ignore_strand)
    }
    if (!equal_df(join_left_got, join_left_ref, sort_by = c("query_id", "subject_id"))) {
      stop("Mismatch: fast_overlap_join(left) at iter=", iter, " ignore_strand=", ignore_strand)
    }
    if (!equal_df(inner_wrap_got, join_inner_ref, sort_by = c("query_id", "subject_id"))) {
      stop("Mismatch: fast_inner_overlap_join at iter=", iter, " ignore_strand=", ignore_strand)
    }
    if (!equal_df(left_wrap_got, join_left_ref, sort_by = c("query_id", "subject_id"))) {
      stop("Mismatch: fast_left_overlap_join at iter=", iter, " ignore_strand=", ignore_strand)
    }
    if (!equal_df(semi_got, semi_ref, sort_by = "query_id")) {
      stop("Mismatch: fast_semi_overlap_join at iter=", iter, " ignore_strand=", ignore_strand)
    }
    if (!equal_df(anti_got, anti_ref, sort_by = "query_id")) {
      stop("Mismatch: fast_anti_overlap_join at iter=", iter, " ignore_strand=", ignore_strand)
    }

    # Group counts and aggregates.
    got_grp <- fastRanges::fast_count_overlaps_by_group(
      qg, sg,
      group_col = "grp",
      include_na_group = TRUE,
      ignore_strand = ignore_strand,
      threads = threads,
      deterministic = FALSE
    )
    ref_grp <- manual_group_counts(qg, sg, "grp", include_na_group = TRUE, hits = ref_join_hits)
    if (!isTRUE(all.equal(unname(got_grp[, colnames(ref_grp), drop = FALSE]), unname(ref_grp)))) {
      stop("Mismatch: fast_count_overlaps_by_group at iter=", iter, " ignore_strand=", ignore_strand)
    }

    for (fun in c("count", "sum", "mean", "min", "max")) {
      got_agg <- fastRanges::fast_overlap_aggregate(
        qg, sg,
        value_col = if (fun == "count") NULL else "score",
        fun = fun,
        ignore_strand = ignore_strand,
        threads = threads,
        deterministic = FALSE
      )
      ref_agg <- manual_overlap_aggregate(
        qg, sg,
        value_col = "score",
        fun = fun,
        hits = ref_join_hits,
        ignore_strand = ignore_strand
      )
      if (!isTRUE(all.equal(got_agg, ref_agg, tolerance = 1e-10, check.attributes = FALSE))) {
        stop("Mismatch: fast_overlap_aggregate(", fun, ") at iter=", iter, " ignore_strand=", ignore_strand)
      }
    }

    # Self overlaps.
    ref_self_all <- GenomicRanges::findOverlaps(qg, qg, ignore.strand = ignore_strand)
    rq <- S4Vectors::queryHits(ref_self_all)
    rs <- S4Vectors::subjectHits(ref_self_all)
    keep <- rq != rs & rq < rs
    ref_self <- S4Vectors::Hits(from = rq[keep], to = rs[keep], nLnode = length(qg), nRnode = length(qg))
    got_self <- fastRanges::fast_self_overlaps(
      qg,
      ignore_strand = ignore_strand,
      threads = threads,
      deterministic = FALSE
    )
    if (!isTRUE(all.equal(canon_hits(got_self), canon_hits(ref_self), check.attributes = FALSE))) {
      stop("Mismatch: fast_self_overlaps at iter=", iter, " ignore_strand=", ignore_strand)
    }

    # Window counts.
    ww <- 1000L
    sw <- 500L
    got_w <- fastRanges::fast_window_count_overlaps(
      qg, sg,
      window_width = ww,
      step_width = sw,
      ignore_strand = ignore_strand,
      threads = threads,
      deterministic = FALSE
    )
    ref_w <- baseline_window_count(qg, sg, window_width = ww, step_width = sw, ignore_strand = ignore_strand)
    if (!identical(got_w$overlap_count, ref_w$overlap_count)) {
      stop("Mismatch: fast_window_count_overlaps at iter=", iter, " ignore_strand=", ignore_strand)
    }

    # Nearest/distance/precede/follow.
    got_dn <- fastRanges::fast_distance_to_nearest(qg, sg, threads = threads, ignore_strand = ignore_strand)
    ref_dn <- GenomicRanges::distanceToNearest(qg, sg, ignore.strand = ignore_strand)
    if (!identical(as.integer(got_dn$query_id), as.integer(S4Vectors::queryHits(ref_dn))) ||
        !identical(as.integer(got_dn$subject_id), as.integer(S4Vectors::subjectHits(ref_dn))) ||
        !isTRUE(all.equal(as.numeric(got_dn$distance), as.numeric(S4Vectors::mcols(ref_dn)$distance), check.attributes = FALSE))) {
      stop("Mismatch: fast_distance_to_nearest at iter=", iter, " ignore_strand=", ignore_strand)
    }
    got_nearest <- fastRanges::fast_nearest(qg, sg, threads = threads, ignore_strand = ignore_strand)
    if (!isTRUE(all.equal(got_nearest, got_dn, check.attributes = FALSE))) {
      stop("Mismatch: fast_nearest at iter=", iter, " ignore_strand=", ignore_strand)
    }

    if (!identical(
      fastRanges::fast_precede(qg, sg, ignore_strand = ignore_strand, threads = threads),
      GenomicRanges::precede(qg, sg, ignore.strand = ignore_strand)
    )) {
      stop("Mismatch: fast_precede at iter=", iter, " ignore_strand=", ignore_strand)
    }
    if (!identical(
      fastRanges::fast_follow(qg, sg, ignore_strand = ignore_strand, threads = threads),
      GenomicRanges::follow(qg, sg, ignore.strand = ignore_strand)
    )) {
      stop("Mismatch: fast_follow at iter=", iter, " ignore_strand=", ignore_strand)
    }

    # Clustering derived from self-overlap graph.
    cl_got <- fastRanges::fast_cluster_overlaps(
      qg,
      ignore_strand = ignore_strand,
      threads = threads,
      deterministic = TRUE
    )
    cl_ref <- cluster_from_hits(length(qg), ref_self)
    if (!same_partition(cl_got, cl_ref)) {
      stop("Mismatch: fast_cluster_overlaps vector at iter=", iter, " ignore_strand=", ignore_strand)
    }
    cl_df_got <- fastRanges::fast_cluster_overlaps(
      qg,
      ignore_strand = ignore_strand,
      threads = threads,
      deterministic = TRUE,
      return = "data.frame"
    )
    if (!identical(as.integer(cl_df_got$range_id), seq_len(length(qg))) ||
        !same_partition(cl_df_got$cluster_id, cl_ref) ||
        !identical(
          as.integer(cl_df_got$cluster_size),
          as.integer(tabulate(cl_ref, nbins = max(cl_ref))[cl_ref])
        )) {
      stop("Mismatch: fast_cluster_overlaps data.frame at iter=", iter, " ignore_strand=", ignore_strand)
    }
  }

  # Coverage and tile coverage.
  cov_got <- fastRanges::fast_coverage(qg, threads = threads)
  cov_ref <- GenomicRanges::coverage(qg)
  if (!isTRUE(all.equal(cov_got, cov_ref, check.attributes = FALSE))) {
    stop("Mismatch: fast_coverage at iter=", iter)
  }
  tc <- fastRanges::fast_tile_coverage(qg, tile_width = 500L, step_width = 500L, threads = threads)
  tc_ref <- baseline_tile_coverage(qg, tile_width = 500L, step_width = 500L)
  if (!equal_df(tc, tc_ref, sort_by = c("seqnames", "start", "end"))) {
    stop("Mismatch: fast_tile_coverage at iter=", iter)
  }

  # IRanges equivalence for core overlap.
  qi <- make_iranges(q_ir_n, seed + iter * 100L)
  si <- make_iranges(s_ir_n, seed + iter * 100L + 1L)
  ref_i <- IRanges::findOverlaps(qi, si)
  got_i <- fastRanges::fast_find_overlaps(qi, si, threads = threads)
  if (!isTRUE(all.equal(canon_hits(got_i), canon_hits(ref_i), check.attributes = FALSE))) {
    stop("Mismatch: fast_find_overlaps IRanges at iter=", iter)
  }
  for (sel in c("first", "last")) {
    ref_sel_i <- IRanges::findOverlaps(qi, si, select = sel)
    got_sel_i <- fastRanges::fast_find_overlaps(qi, si, select = sel, threads = threads)
    if (!identical(got_sel_i, ref_sel_i)) {
      stop("Mismatch: fast_find_overlaps IRanges select=", sel, " at iter=", iter)
    }
  }
  got_sel_i_arb <- fastRanges::fast_find_overlaps(qi, si, select = "arbitrary", threads = threads)
  if (!validate_select_result(got_sel_i_arb, ref_i, length(qi), mode = "arbitrary")) {
    stop("Mismatch: fast_find_overlaps IRanges select=arbitrary semantic validity at iter=", iter)
  }

  message("iter ", iter, "/", n_iter, " OK")
}

message("All fastRanges output checks passed.")
