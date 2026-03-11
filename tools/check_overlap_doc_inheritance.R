#!/usr/bin/env Rscript

read_lines <- function(path) {
  x <- readLines(path, warn = FALSE)
  if (!length(x)) character() else x
}

extract_doc_block <- function(lines, fn_line) {
  i <- fn_line - 1L
  block <- character()
  while (i >= 1L) {
    line <- lines[[i]]
    if (grepl("^#'", line)) {
      block <- c(line, block)
      i <- i - 1L
      next
    }
    if (grepl("^\\s*$", line)) {
      i <- i - 1L
      next
    }
    break
  }
  block
}

r_files <- list.files("R", pattern = "[.]R$", full.names = TRUE)
issues <- character()

for (f in r_files) {
  lines <- read_lines(f)
  if (!length(lines)) {
    next
  }

  fn_lines <- grep("<-\\s*function\\s*\\(", lines)
  for (ln in fn_lines) {
    header <- lines[[ln]]
    fn_name <- sub("^\\s*([.A-Za-z0-9_]+)\\s*<-\\s*function\\s*\\(.*$", "\\1", header)
    if (identical(fn_name, header)) {
      next
    }

    scan_to <- min(length(lines), ln + 40L)
    sig <- paste(lines[ln:scan_to], collapse = "\n")
    has_overlap_type <- grepl('type\\s*=\\s*c\\("any",\\s*"start",\\s*"end",\\s*"within",\\s*"equal"\\)', sig)
    if (!has_overlap_type) {
      next
    }

    if (fn_name == "fast_find_overlaps") {
      next
    }

    doc <- extract_doc_block(lines, ln)
    has_inherit <- any(grepl(
      "@inheritParams\\s+(fast_find_overlaps|fast_overlap_join|fast_inner_overlap_join|fast_semi_overlap_join|fast_self_overlaps)",
      doc
    ))
    has_template <- any(grepl("@template\\s+overlap_shared_args", doc))

    if (!(has_inherit || has_template)) {
      issues <- c(issues, sprintf("%s:%d (%s)", f, ln, fn_name))
    }
  }
}

if (length(issues)) {
  cat("Functions with overlap `type` argument must inherit or template shared docs.\n", sep = "")
  cat("Add `@inheritParams fast_find_overlaps` (or approved inherited chain) or `@template overlap_shared_args`.\n\n", sep = "")
  cat(paste0(" - ", issues, collapse = "\n"), "\n", sep = "")
  quit(status = 1L)
}

cat("Overlap argument documentation inheritance check: OK\n")
