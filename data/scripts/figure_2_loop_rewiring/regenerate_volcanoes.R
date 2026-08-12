#!/usr/bin/env Rscript
# data/scripts/figure_2_loop_rewiring/regenerate_volcanoes.R
# Regenerate every volcano artifact using the fixed create_publication_volcano().
#
# Scope:
#   - Both timepoints (late, early)
#   - All five volcano variants per timepoint (5kb, 10kb, 25kb, merged multi-res, non-redundant)
#   - Overwrites outputs/{timepoint}/visualizations/volcano/* and promotes the merged
#     hero plot to data/plots/figure_2_loop_rewiring/2A_{timepoint}_volcano_merged.{pdf,svg,jpg}
#
# Why this exists:
#   The shared create_publication_volcano() in visualizations.R previously derived
#   its "significant" counts from an upstream `significant` column that flags FDR<0.05
#   only (ignoring the |logFC|>0.3 cutoff drawn on the plot) AND placed the up/down
#   annotations on the wrong sides. Both defects are now fixed at the source; this
#   driver regenerates every output so no stale file remains.
#
# Usage:
#   Rscript data/scripts/figure_2_loop_rewiring/regenerate_volcanoes.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(EnhancedVolcano)
  library(yaml)
})

# Reuse the existing multi-format saver (pdf + svg + jpg)
source("data/scripts/_shared/multi_format_output.R")

# Honor the repo's base_dir convention (matches visualizations.R setup block)
if (file.exists("config/paths_config.yaml")) {
  config <- yaml::read_yaml("config/paths_config.yaml")
  base_dir <- config$project$base_dir
  if (dir.exists(base_dir)) setwd(base_dir)
}

# Pull the (fixed) create_publication_volcano() from the bundled visualizations.R
# without triggering its top-level workflow. Brace-balanced extraction — robust
# to the function body length and internal nesting.
viz_file <- "data/scripts/figure_2_loop_rewiring/visualizations.R"
viz_lines <- readLines(viz_file)
fn_start_line <- grep("^create_publication_volcano <- function", viz_lines)
if (length(fn_start_line) != 1) {
  stop(sprintf("Expected exactly one create_publication_volcano definition in %s; found %d",
               viz_file, length(fn_start_line)))
}
depth <- 0L
fn_end_line <- NA_integer_
for (i in seq(fn_start_line, length(viz_lines))) {
  chars <- strsplit(viz_lines[i], "", fixed = TRUE)[[1]]
  depth <- depth + sum(chars == "{") - sum(chars == "}")
  if (depth == 0L && i > fn_start_line) {
    fn_end_line <- i
    break
  }
}
if (is.na(fn_end_line)) stop("Could not locate closing brace of create_publication_volcano")
eval(parse(text = paste(viz_lines[fn_start_line:fn_end_line], collapse = "\n")))
cat(sprintf("Loaded create_publication_volcano() from %s (lines %d-%d)\n\n",
            viz_file, fn_start_line, fn_end_line))

# -----------------------------------------------------------------------------
# Regenerate
# -----------------------------------------------------------------------------

timepoints <- list(
  late  = "outputs/250402-late_outputs",
  early = "outputs/250831-early_outputs"
)

curated_dir <- "data/plots/figure_2_loop_rewiring"
dir.create(curated_dir, recursive = TRUE, showWarnings = FALSE)

for (tp in names(timepoints)) {
  tp_root <- timepoints[[tp]]
  vol_dir <- file.path(tp_root, "visualizations/volcano")
  dir.create(vol_dir, recursive = TRUE, showWarnings = FALSE)

  cat("========================================\n")
  cat(sprintf("Timepoint: %s  (%s)\n", tp, tp_root))
  cat("========================================\n\n")

  # Per-resolution volcanoes
  for (res_kb in c(5, 10, 25)) {
    create_publication_volcano(
      sprintf("%s/edgeR_results_res_%dkb/primary_analysis/all_results_primary.tsv",
              tp_root, res_kb),
      file.path(vol_dir, sprintf("volcano_%dkb.pdf", res_kb)),
      sprintf("%dkb", res_kb)
    )
  }

  # Merged multi-resolution (the hero figure promoted to data/plots/)
  merged_all_results_file <- file.path(tp_root, "merged_loops/merged_all_results.tsv")
  if (file.exists(merged_all_results_file)) {
    create_publication_volcano(
      merged_all_results_file,
      file.path(vol_dir, "volcano_merged_multiresolution.pdf"),
      "Multi-Resolution"
    )
  } else {
    cat(sprintf("  ! Skipping merged multi-res (not found: %s)\n\n", merged_all_results_file))
  }

  # Non-redundant (stringent characterized set)
  characterized_file <- file.path(tp_root, "merged_loops/characterized_loops.tsv")
  if (file.exists(characterized_file)) {
    create_publication_volcano(
      characterized_file,
      file.path(vol_dir, "volcano_nonredundant.pdf"),
      "Non-Redundant"
    )
  } else {
    cat(sprintf("  ! Skipping non-redundant (not found: %s)\n\n", characterized_file))
  }

  # Promote merged hero figure to the curated location, overwriting stale files.
  # save_multiformat_ggplot() writes into a subfolder named after the figure.
  src_subdir <- file.path(vol_dir, "volcano_merged_multiresolution")
  if (dir.exists(src_subdir)) {
    for (ext in c("pdf", "svg", "jpg")) {
      src <- file.path(src_subdir, paste0("volcano_merged_multiresolution.", ext))
      dst <- file.path(curated_dir, sprintf("2A_%s_volcano_merged.%s", tp, ext))
      if (file.exists(src)) {
        ok <- file.copy(src, dst, overwrite = TRUE)
        cat(sprintf("  %s  %s -> %s\n", ifelse(ok, "\u2713", "\u2717"), src, dst))
      }
    }
    cat("\n")
  }
}

cat("Done. Regenerated volcano artifacts for all timepoints.\n")
