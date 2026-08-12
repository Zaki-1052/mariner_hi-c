# tads/scripts/cre_blacklist_utils.R
# Shared functions for Math1-cre transgene artifact blacklist filtering

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(readr)
})

find_repo_root <- function(start = getwd()) {
  d <- normalizePath(start)
  for (i in 1:6) {
    if (file.exists(file.path(d, "CLAUDE.md")) &&
        dir.exists(file.path(d, "tads")) &&
        dir.exists(file.path(d, "loops"))) {
      return(d)
    }
    parent <- dirname(d)
    if (parent == d) break
    d <- parent
  }
  env_root <- Sys.getenv("MARINER_REPO_ROOT", "")
  if (nzchar(env_root) && dir.exists(env_root)) return(normalizePath(env_root))
  stop("Could not locate mariner_hi-c repo root from: ", start,
       "\n  Set MARINER_REPO_ROOT env var to override")
}

load_cre_blacklist <- function(repo_root) {
  bl_path <- file.path(repo_root, "tads", "cre_artifact_blacklist.bed")
  if (!file.exists(bl_path)) stop("Cre blacklist not found: ", bl_path)
  bl <- import.bed(bl_path)
  cat("  Loaded cre blacklist:", length(bl), "regions,",
      sum(width(bl)), "bp total\n")
  for (i in seq_along(bl)) {
    cat("    ", as.character(seqnames(bl)[i]), ":",
        start(bl)[i], "-", end(bl)[i],
        " (", mcols(bl)$name[i], ")\n", sep = "")
  }
  reduce(bl)
}

normalize_chr <- function(x) {
  ifelse(grepl("^chr", x), x, paste0("chr", x))
}

flag_single <- function(df, chr_col, start_col, end_col, blacklist,
                        resolution_bp = NA, needs_chr_prefix = FALSE) {
  chrs <- df[[chr_col]]
  if (needs_chr_prefix) chrs <- normalize_chr(chrs)
  starts <- as.numeric(df[[start_col]])
  if (!is.na(resolution_bp)) {
    ends <- starts + resolution_bp - 1
  } else {
    ends <- as.numeric(df[[end_col]])
  }
  valid <- !is.na(chrs) & !is.na(starts) & !is.na(ends)
  result <- rep(FALSE, nrow(df))
  if (sum(valid) == 0) return(result)
  gr <- GRanges(
    seqnames = chrs[valid],
    ranges = IRanges(start = starts[valid], end = ends[valid])
  )
  result[valid] <- countOverlaps(gr, blacklist) > 0
  result
}

flag_paired <- function(df, chr1_col, start1_col, end1_col,
                        chr2_col, start2_col, end2_col, blacklist,
                        needs_chr_prefix = FALSE) {
  c1 <- df[[chr1_col]]
  c2 <- df[[chr2_col]]
  if (needs_chr_prefix) {
    c1 <- normalize_chr(c1)
    c2 <- normalize_chr(c2)
  }
  s1 <- as.numeric(df[[start1_col]])
  e1 <- as.numeric(df[[end1_col]])
  s2 <- as.numeric(df[[start2_col]])
  e2 <- as.numeric(df[[end2_col]])

  valid <- !is.na(c1) & !is.na(s1) & !is.na(e1) &
           !is.na(c2) & !is.na(s2) & !is.na(e2)
  result <- rep(FALSE, nrow(df))
  if (sum(valid) == 0) return(result)

  gr1 <- GRanges(seqnames = c1[valid], ranges = IRanges(s1[valid], e1[valid]))
  gr2 <- GRanges(seqnames = c2[valid], ranges = IRanges(s2[valid], e2[valid]))
  result[valid] <- (countOverlaps(gr1, blacklist) > 0) |
                   (countOverlaps(gr2, blacklist) > 0)
  result
}

identify_region <- function(df, hit_mask, chr_col, start_col, end_col,
                            blacklist_raw, needs_chr_prefix = FALSE) {
  if (sum(hit_mask) == 0) return(character(0))
  sub <- df[hit_mask, , drop = FALSE]
  chrs <- sub[[chr_col]]
  if (needs_chr_prefix) chrs <- normalize_chr(chrs)
  starts <- as.numeric(sub[[start_col]])
  ends <- as.numeric(sub[[end_col]])
  valid <- !is.na(chrs) & !is.na(starts) & !is.na(ends)
  regions <- rep(NA_character_, nrow(sub))
  if (sum(valid) == 0) return(regions)
  gr <- GRanges(seqnames = chrs[valid], ranges = IRanges(starts[valid], ends[valid]))
  fo <- findOverlaps(gr, blacklist_raw)
  regions[valid][queryHits(fo)] <- mcols(blacklist_raw)$name[subjectHits(fo)]
  regions
}

write_cre_outputs <- function(df, in_blacklist, out_path, label,
                              gr_original = NULL) {
  stem <- tools::file_path_sans_ext(out_path)
  ext <- tools::file_ext(out_path)
  filtered_path <- paste0(stem, "_crefiltered.", ext)
  removed_path <- paste0(stem, "_cre_removed.", ext)

  if (!is.null(gr_original) && ext == "bed") {
    export.bed(gr_original[!in_blacklist], filtered_path)
    if (any(in_blacklist)) export.bed(gr_original[in_blacklist], removed_path)
  } else {
    write_tsv(df[!in_blacklist, ], filtered_path)
    if (any(in_blacklist)) write_tsv(df[in_blacklist, ], removed_path)
  }

  cat("    Wrote", sum(!in_blacklist), "rows ->", basename(filtered_path), "\n")
  if (any(in_blacklist)) {
    cat("    Wrote", sum(in_blacklist), "removed ->", basename(removed_path), "\n")
  }

  invisible(NULL)
}
