#' @keywords internal
.encode_strand <- function(x) {
  if (!.is_granges(x)) {
    return(integer(length(x)))
  }

  strand_chr <- as.character(GenomicRanges::strand(x))
  out <- integer(length(strand_chr))
  out[strand_chr == "+"] <- 1L
  out[strand_chr == "-"] <- 2L
  out
}

#' @keywords internal
.encode_subject <- function(subject) {
  .assert_supported_ranges(subject, "subject")

  start <- as.integer(IRanges::start(subject))
  end <- as.integer(IRanges::end(subject))
  strand <- .encode_strand(subject)

  if (.is_granges(subject)) {
    seq_chr <- as.character(GenomicRanges::seqnames(subject))
    seq_levels <- unique(seq_chr)
    seq_map <- stats::setNames(seq_along(seq_levels), seq_levels)
    seq_id <- as.integer(unname(seq_map[seq_chr]))
  } else {
    seq_levels <- "range"
    seq_map <- stats::setNames(1L, seq_levels)
    seq_id <- rep.int(1L, length(subject))
  }

  list(
    start = start,
    end = end,
    seq_id = seq_id,
    strand = strand,
    seq_map = seq_map,
    seq_levels = seq_levels
  )
}

#' @keywords internal
.encode_query <- function(query, seq_map) {
  .assert_supported_ranges(query, "query")

  start <- as.integer(IRanges::start(query))
  end <- as.integer(IRanges::end(query))
  strand <- .encode_strand(query)

  if (.is_granges(query)) {
    seq_chr <- as.character(GenomicRanges::seqnames(query))
    seq_id <- as.integer(unname(seq_map[seq_chr]))
    seq_id[is.na(seq_id)] <- 0L
  } else {
    seq_id <- rep.int(1L, length(query))
  }

  list(
    start = start,
    end = end,
    seq_id = seq_id,
    strand = strand
  )
}

#' @keywords internal
.build_sorted_subject_vectors <- function(subject) {
  encoded <- .encode_subject(subject)

  n <- length(encoded$start)
  if (n == 0L) {
    return(list(
      subject_start = integer(),
      subject_end = integer(),
      subject_seq = integer(),
      subject_strand = integer(),
      subject_original_index = integer(),
      partition_keys = integer(),
      partition_starts = integer(),
      partition_ends = integer(),
      seq_map = encoded$seq_map,
      seq_levels = encoded$seq_levels
    ))
  }

  ord <- order(encoded$seq_id, encoded$start, encoded$end)
  subject_start <- encoded$start[ord]
  subject_end <- encoded$end[ord]
  subject_seq <- encoded$seq_id[ord]
  subject_strand <- encoded$strand[ord]

  seq_rle <- rle(subject_seq)
  partition_ends <- cumsum(seq_rle$lengths)
  partition_starts <- partition_ends - seq_rle$lengths + 1L

  list(
    subject_start = subject_start,
    subject_end = subject_end,
    subject_seq = subject_seq,
    subject_strand = subject_strand,
    subject_original_index = as.integer(ord),
    partition_keys = as.integer(seq_rle$values),
    partition_starts = as.integer(partition_starts),
    partition_ends = as.integer(partition_ends),
    seq_map = encoded$seq_map,
    seq_levels = encoded$seq_levels
  )
}
