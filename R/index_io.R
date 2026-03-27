._add_block_index_fields <- function(index) {
  if (is.null(index$has_empty_ranges)) {
    index$has_empty_ranges <- any((as.integer(index$subject_end) - as.integer(index$subject_start) + 1L) == 0L)
  }
  if (is.null(index$has_circular_sequences)) {
    index$has_circular_sequences <- FALSE
  }

  if (!is.null(index$block_starts) &&
      !is.null(index$block_ends) &&
      !is.null(index$block_first_start) &&
      !is.null(index$block_max_end) &&
      !is.null(index$partition_block_starts) &&
      !is.null(index$partition_block_ends)) {
    return(index)
  }

  block_size <- 2048L
  partition_starts <- as.integer(index$partition_starts)
  partition_ends <- as.integer(index$partition_ends)
  subject_start <- as.integer(index$subject_start)
  subject_end <- as.integer(index$subject_end)

  block_starts <- integer()
  block_ends <- integer()
  block_first_start <- integer()
  block_max_end <- integer()
  partition_block_starts <- integer(length(partition_starts))
  partition_block_ends <- integer(length(partition_starts))

  next_block <- 1L
  for (i in seq_along(partition_starts)) {
    p_start <- partition_starts[[i]]
    p_end <- partition_ends[[i]]
    starts <- seq.int(p_start, p_end, by = block_size)
    ends <- pmin.int(starts + block_size - 1L, p_end)
    n_blocks <- length(starts)
    if (n_blocks == 0L) {
      partition_block_starts[[i]] <- 1L
      partition_block_ends[[i]] <- 0L
      next
    }
    max_end <- vapply(
      seq_len(n_blocks),
      function(j) max(subject_end[starts[[j]]:ends[[j]]]),
      integer(1)
    )
    idx <- seq.int(next_block, length.out = n_blocks)
    partition_block_starts[[i]] <- idx[[1]]
    partition_block_ends[[i]] <- idx[[length(idx)]]
    next_block <- next_block + n_blocks
    block_starts <- c(block_starts, starts)
    block_ends <- c(block_ends, ends)
    block_first_start <- c(block_first_start, subject_start[starts])
    block_max_end <- c(block_max_end, max_end)
  }

  index$block_size <- as.integer(block_size)
  index$block_starts <- as.integer(block_starts)
  index$block_ends <- as.integer(block_ends)
  index$block_first_start <- as.integer(block_first_start)
  index$block_max_end <- as.integer(block_max_end)
  index$partition_block_starts <- as.integer(partition_block_starts)
  index$partition_block_ends <- as.integer(partition_block_ends)
  index
}

#' Save a Reusable Subject Index
#'
#' Save a `fast_ranges_index` object to disk for reuse across sessions.
#'
#' @param index A `fast_ranges_index` object created by `fast_build_index()`.
#' @param path Output file path for the serialized index.
#' @param compress Logical scalar. If `TRUE`, uses xz compression.
#'
#' @details
#' Saving an index is useful when index construction is expensive and the same
#' subject set will be queried again in a later R session.
#'
#' @return Invisibly returns `path`.
#' @export
#'
#' @examples
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' idx <- fast_build_index(s)
#' f <- tempfile(fileext = ".rds")
#' fast_save_index(idx, f)
#' unlink(f)
fast_save_index <- function(index, path, compress = TRUE) {
  .assert_fast_ranges_index(index, "index")
  .assert_scalar_character(path, "path")
  .assert_scalar_logical(compress, "compress")

  out_dir <- dirname(path)
  if (!dir.exists(out_dir)) {
    stop("Directory for `path` does not exist: ", out_dir, call. = FALSE)
  }

  saveRDS(
    object = index,
    file = path,
    compress = if (isTRUE(compress)) "xz" else FALSE
  )

  invisible(path)
}

#' Load a Reusable Subject Index
#'
#' Load a `fast_ranges_index` object saved with `fast_save_index()`.
#'
#' @param path File path to a serialized index.
#'
#' @details
#' This function validates that the file contains the fields required by
#' `fastRanges` before returning the object.
#'
#' @return A `fast_ranges_index` object.
#' @export
#'
#' @examples
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' idx <- fast_build_index(s)
#' f <- tempfile(fileext = ".rds")
#' fast_save_index(idx, f)
#' idx2 <- fast_load_index(f)
#' unlink(f)
#' print(idx2)
fast_load_index <- function(path) {
  .assert_scalar_character(path, "path")
  if (!file.exists(path)) {
    stop("`path` does not exist: ", path, call. = FALSE)
  }

  index <- readRDS(path)
  .assert_fast_ranges_index(index, "index")

  required_fields <- c(
    "subject_start", "subject_end", "subject_seq", "subject_strand",
    "subject_original_index", "partition_keys", "partition_starts",
    "partition_ends", "seq_map", "seq_levels", "subject_n", "subject_class"
  )
  missing_fields <- setdiff(required_fields, names(index))
  if (length(missing_fields) > 0L) {
    stop(
      "`path` does not contain a valid `fast_ranges_index` object. Missing fields: ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  index <- ._add_block_index_fields(index)
  index
}

#' Index Summary Statistics
#'
#' Summarize index size and partition structure.
#'
#' @param index A `fast_ranges_index` object.
#' @param detailed Logical scalar. If `TRUE`, returns partition-level details.
#'
#' @details
#' Use this function to inspect how the subject was partitioned internally and
#' to get a rough sense of memory footprint.
#'
#' `subject_n` is the number of indexed ranges.
#'
#' `partition_n` is the number of internal partitions, usually driven by
#' sequence structure.
#'
#' `index_size_mb` is the in-memory object size, not the serialized file size.
#'
#' @return By default, a one-row `S4Vectors::DataFrame` with summary fields.
#'   If `detailed = TRUE`, returns a list with `summary` and `partitions`.
#' @export
#'
#' @examples
#' s <- IRanges::IRanges(start = c(3, 9, 18), width = c(4, 6, 5))
#' idx <- fast_build_index(s)
#' fast_index_stats(idx)
fast_index_stats <- function(index, detailed = FALSE) {
  .assert_fast_ranges_index(index, "index")
  .assert_scalar_logical(detailed, "detailed")

  partition_n <- length(index$partition_keys)
  mean_partition_size <- if (partition_n == 0L) 0 else index$subject_n / partition_n

  summary_df <- S4Vectors::DataFrame(
    subject_n = as.integer(index$subject_n),
    partition_n = as.integer(partition_n),
    block_n = as.integer(length(index$block_starts %||% integer())),
    seqlevel_n = as.integer(length(index$seq_map)),
    mean_partition_size = as.numeric(mean_partition_size),
    index_size_mb = as.numeric(utils::object.size(index)) / (1024^2)
  )

  if (!isTRUE(detailed)) {
    return(summary_df)
  }

  if (partition_n == 0L) {
    partition_df <- data.frame(
      seq_id = integer(),
      seqname = character(),
      partition_start = integer(),
      partition_end = integer(),
      n_ranges = integer(),
      stringsAsFactors = FALSE
    )
  } else {
    id_to_name <- names(index$seq_map)
    partition_df <- data.frame(
      seq_id = as.integer(index$partition_keys),
      seqname = as.character(id_to_name[as.integer(index$partition_keys)]),
      partition_start = as.integer(index$partition_starts),
      partition_end = as.integer(index$partition_ends),
      n_ranges = as.integer(index$partition_ends - index$partition_starts + 1L),
      block_start = as.integer(index$partition_block_starts %||% integer(length(index$partition_keys))),
      block_end = as.integer(index$partition_block_ends %||% integer(length(index$partition_keys))),
      stringsAsFactors = FALSE
    )
    partition_df$n_blocks <- pmax.int(0L, partition_df$block_end - partition_df$block_start + 1L)
  }

  list(summary = summary_df, partitions = partition_df)
}
