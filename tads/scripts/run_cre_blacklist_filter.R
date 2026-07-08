# tads/scripts/run_cre_blacklist_filter.R
# Manifest-driven engine for Math1-cre transgene artifact blacklist filtering
#
# Usage:
#   Rscript tads/scripts/run_cre_blacklist_filter.R <mode> <module>
#   mode:   "audit" (dry run, counts only) or "filter" (writes *_crefiltered files)
#   module: "all", "loops", "stripes", "tads", "cmpts", "peaks", "abc", "biomodal"

# ============================================================================
# Setup
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_cre_blacklist_filter.R <audit|filter> <module>\n",
       "  module: all, loops, stripes, tads, cmpts, peaks, abc, biomodal")
}

MODE <- match.arg(args[1], c("audit", "filter"))
MODULE <- args[2]
VALID_MODULES <- c("all", "loops", "stripes", "tads", "cmpts", "peaks", "abc", "biomodal")
if (!MODULE %in% VALID_MODULES) {
  stop("Invalid module: ", MODULE, "\n  Valid: ", paste(VALID_MODULES, collapse = ", "))
}

get_script_dir <- function() {
  # Works with Rscript (commandArgs), source(), and interactive
  for (frame in rev(sys.frames())) {
    if (!is.null(frame$ofile)) return(dirname(normalizePath(frame$ofile)))
  }
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", ca, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  getwd()
}
script_dir <- get_script_dir()
source(file.path(script_dir, "cre_blacklist_utils.R"))

REPO_ROOT <- find_repo_root(script_dir)
cat("=== Math1-cre Transgene Blacklist Filter ===\n")
cat("Mode:", MODE, "\n")
cat("Module:", MODULE, "\n")
cat("Repo root:", REPO_ROOT, "\n\n")

# ============================================================================
# Load blacklist
# ============================================================================

cat("Loading cre artifact blacklist...\n")
bl_path <- file.path(REPO_ROOT, "tads", "cre_artifact_blacklist.bed")
blacklist_raw <- import.bed(bl_path)
blacklist <- reduce(blacklist_raw)
cat("\n")

# ============================================================================
# Load and filter manifest
# ============================================================================

manifest_path <- file.path(script_dir, "cre_blacklist_manifest.tsv")
if (!file.exists(manifest_path)) stop("Manifest not found: ", manifest_path)

manifest <- read_tsv(manifest_path, show_col_types = FALSE,
                     col_types = cols(.default = "c"))

if (MODULE != "all") {
  manifest <- manifest %>% filter(module == MODULE)
}

if (nrow(manifest) == 0) {
  stop("No manifest entries for module: ", MODULE)
}

cat("Manifest entries to process:", nrow(manifest), "\n\n")

# ============================================================================
# Process each manifest entry
# ============================================================================

results <- list()

for (i in seq_len(nrow(manifest))) {
  entry <- manifest[i, ]
  file_pattern <- file.path(REPO_ROOT, entry$file_relative)

  # Expand globs
  if (grepl("\\*", file_pattern)) {
    files <- Sys.glob(file_pattern)
    if (length(files) == 0) {
      cat("[SKIP]", entry$module, "/", entry$label,
          "- no files matching:", entry$file_relative, "\n")
      next
    }
  } else {
    files <- file_pattern
    if (!file.exists(files)) {
      cat("[SKIP]", entry$module, "/", entry$label,
          "- file not found:", entry$file_relative, "\n")
      next
    }
  }

  for (fpath in files) {
    rel_path <- sub(paste0(REPO_ROOT, "/"), "", fpath)
    file_label <- if (length(files) > 1) {
      paste0(entry$label, ":", basename(fpath))
    } else {
      entry$label
    }

    cat("[", entry$module, "/", entry$timepoint, "] ", file_label, "\n", sep = "")

    has_header <- entry$has_header == "TRUE"
    needs_chr_prefix <- entry$needs_chr_prefix == "TRUE"
    resolution_bp <- if (is.na(entry$resolution_bp) || entry$resolution_bp == "NA") NA else as.numeric(entry$resolution_bp)

    # ---- Read file ----
    gr_original <- NULL
    if (entry$mode == "bed") {
      gr_original <- import.bed(fpath)
      df <- as.data.frame(gr_original)
      hit_mask <- countOverlaps(gr_original, blacklist) > 0

    } else {
      df <- tryCatch(
        read_tsv(fpath, show_col_types = FALSE, col_names = has_header),
        error = function(e) {
          cat("    ERROR reading file:", conditionMessage(e), "\n")
          return(NULL)
        }
      )
      if (is.null(df) || nrow(df) == 0) {
        cat("    Empty or unreadable, skipping\n")
        next
      }

      # ---- Flag overlaps ----
      if (entry$mode == "paired") {
        required_cols <- c(entry$chr1_col, entry$start1_col, entry$end1_col,
                           entry$chr2_col, entry$start2_col, entry$end2_col)
        missing <- setdiff(required_cols, colnames(df))
        if (length(missing) > 0) {
          cat("    ERROR: missing columns:", paste(missing, collapse = ", "), "\n")
          next
        }
        hit_mask <- flag_paired(df,
                                entry$chr1_col, entry$start1_col, entry$end1_col,
                                entry$chr2_col, entry$start2_col, entry$end2_col,
                                blacklist, needs_chr_prefix)
      } else {
        required_cols <- c(entry$chr1_col, entry$start1_col)
        if (!is.na(resolution_bp)) {
          # Point + resolution mode (tads)
        } else {
          required_cols <- c(required_cols, entry$end1_col)
        }
        missing <- setdiff(required_cols, colnames(df))
        if (length(missing) > 0) {
          cat("    ERROR: missing columns:", paste(missing, collapse = ", "), "\n")
          next
        }
        hit_mask <- flag_single(df,
                                entry$chr1_col, entry$start1_col,
                                if (is.na(resolution_bp)) entry$end1_col else entry$start1_col,
                                blacklist, resolution_bp, needs_chr_prefix)
      }
    }

    n_total <- nrow(df)
    n_hit <- sum(hit_mask)
    pct <- round(100 * n_hit / n_total, 2)

    cat("    Total:", n_total, " | Overlap:", n_hit,
        " (", pct, "%)\n", sep = "")

    # ---- Region breakdown ----
    if (n_hit > 0) {
      chr_col_for_id <- if (entry$mode == "bed") "seqnames" else entry$chr1_col
      start_col_for_id <- if (entry$mode == "bed") "start" else entry$start1_col
      end_col_for_id <- if (entry$mode == "bed") "end" else {
        if (is.na(resolution_bp)) entry$end1_col else entry$start1_col
      }
      regions <- identify_region(df, hit_mask,
                                 chr_col_for_id, start_col_for_id, end_col_for_id,
                                 blacklist_raw, needs_chr_prefix)
      if (entry$mode == "paired" && entry$mode != "bed") {
        regions2 <- identify_region(df, hit_mask,
                                    entry$chr2_col, entry$start2_col, entry$end2_col,
                                    blacklist_raw, needs_chr_prefix)
        regions <- ifelse(is.na(regions), regions2,
                          ifelse(is.na(regions2), regions,
                                 paste(regions, regions2, sep = "+")))
      }
      region_counts <- sort(table(regions, useNA = "ifany"), decreasing = TRUE)
      for (rn in names(region_counts)) {
        cat("      ", ifelse(is.na(rn), "<unresolved>", rn),
            ": ", region_counts[rn], "\n", sep = "")
      }
    }

    # ---- Write filtered output (filter mode only) ----
    if (MODE == "filter" && n_hit > 0) {
      write_cre_outputs(df, hit_mask, fpath, file_label, gr_original)
    } else if (MODE == "filter" && n_hit == 0) {
      cat("    No overlaps - no filtered file needed\n")
    }

    # Accumulate for report
    results[[length(results) + 1]] <- list(
      module = entry$module,
      timepoint = entry$timepoint,
      label = file_label,
      file = rel_path,
      n_total = n_total,
      n_hit = n_hit,
      pct = pct
    )
  }
}

# ============================================================================
# Write summary report
# ============================================================================

report_dir <- file.path(REPO_ROOT, "tads", "results")
if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)

report_tsv_path <- file.path(report_dir, "cre_blacklist_audit_report.tsv")
report_txt_path <- file.path(report_dir, "cre_blacklist_audit_report.txt")

# Build summary data frame from audit entries
audit_rows <- lapply(results, function(r) {
  if (!is.null(r$n_total)) {
    data.frame(
      module = r$module, timepoint = r$timepoint, label = r$label,
      file = r$file, n_total = r$n_total, n_hit = r$n_hit, pct = r$pct,
      stringsAsFactors = FALSE
    )
  }
})
audit_df <- do.call(rbind, Filter(Negate(is.null), audit_rows))

if (!is.null(audit_df) && nrow(audit_df) > 0) {
  write_tsv(audit_df, report_tsv_path)
  cat("\nAudit report (TSV):", report_tsv_path, "\n")

  # Human-readable summary
  sink(report_txt_path)
  cat("Math1-cre Transgene Artifact Blacklist - Audit Report\n")
  cat("=====================================================\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Mode:", MODE, "\n")
  cat("Module:", MODULE, "\n")
  cat("Blacklist: tads/cre_artifact_blacklist.bed (3 regions)\n")
  cat("  chr6:64,725,000-64,735,000  Atoh1 locus\n")
  cat("  chr8:49,660,000-49,835,000  Math1-cre insertion site\n")
  cat("  chr8:30,000,000-43,000,000  Haplotype block\n\n")

  for (mod in unique(audit_df$module)) {
    mod_df <- audit_df[audit_df$module == mod, ]
    total_rows <- sum(mod_df$n_total)
    total_hits <- sum(mod_df$n_hit)
    cat(sprintf("== %s ==  (%d files, %d total rows, %d overlapping [%.1f%%])\n",
                toupper(mod), nrow(mod_df), total_rows, total_hits,
                100 * total_hits / max(total_rows, 1)))
    for (j in seq_len(nrow(mod_df))) {
      r <- mod_df[j, ]
      marker <- if (r$n_hit > 0) "*" else " "
      cat(sprintf("  %s %-40s %6d rows, %4d hit (%5.1f%%)\n",
                  marker, r$label, r$n_total, r$n_hit, r$pct))
    }
    cat("\n")
  }

  total_all <- sum(audit_df$n_total)
  total_hit <- sum(audit_df$n_hit)
  cat(sprintf("TOTAL: %d rows across %d files, %d overlapping (%.2f%%)\n",
              total_all, nrow(audit_df), total_hit, 100 * total_hit / max(total_all, 1)))
  sink()

  cat("Audit report (TXT):", report_txt_path, "\n")
} else {
  cat("\nNo files were processed.\n")
}

# ============================================================================
# Final summary
# ============================================================================

cat("\n=== Done ===\n")
if (MODE == "audit") {
  cat("Audit complete. Review reports before running:\n")
  cat("  Rscript tads/scripts/run_cre_blacklist_filter.R filter", MODULE, "\n")
} else {
  n_filtered <- sum(sapply(results, function(r) !is.null(r$n_hit) && r$n_hit > 0))
  cat("Filter complete.", n_filtered, "files had rows removed.\n")
  cat("Look for *_crefiltered.* and *_cre_removed.* files alongside originals.\n")
}
