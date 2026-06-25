# scripts/viz_sections/_figure_config.R
# =============================================================================
# Shared utilities for consolidated section scripts (sections 79-88).
#
# Sourced by each consolidated section script AFTER _shared_config.R.
# By the time this file is sourced, the following are already in scope
# and MUST NOT be reloaded:
#   - BASE_DIR, REPO_ROOT, OUTPUT_DIR, TABLES_DIR  (path constants)
#   - COLORS, theme_biomodal(), CHROMATIN_STATE_ORDER/_COLORS
#   - save_multiformat_ggplot()/_base()/_pheatmap()  (multi_format_output.R)
#   - mc_dmr, hmc_dmr, region_dmrs, bioqc, upstream, SAMPLES, KEY_GENES
#
# Working directory MUST be downstream/ (BASE_DIR == getwd()).
#
# This file provides, in order:
#   1. theme_pub(base_size = 8)   - publication theme extending theme_biomodal()
#   2. FIGURE_DIR + FIGURE_TABLES_DIR (auto-created)
#   3. save_figure(plot, name, width_mm, height_mm) - mm -> in, 300 DPI wrapper
#   4. read_table_tsv(filename)   - reads a pre-computed TSV from TABLES_DIR
#   5. add_panel_labels(plot)     - patchwork tag_levels = "A"
#   6. stat_label(name, value, fmt) - formatted on-panel annotation helper
#   7. PUB_COLORS                 - canonical colour re-exports from COLORS
#
# Output formats: the save helpers emit PDF + SVG + JPG ONLY (never PNG).
# =============================================================================

# -----------------------------------------------------------------------------
# Fail loudly if the prerequisite shared config has not been sourced first.
# -----------------------------------------------------------------------------
if (!exists("BASE_DIR") || !exists("TABLES_DIR") || !exists("COLORS") ||
    !exists("theme_biomodal") || !exists("save_multiformat_ggplot")) {
  stop(
    "_figure_config.R requires scripts/viz_sections/_shared_config.R to be ",
    "sourced first (provides BASE_DIR, TABLES_DIR, COLORS, theme_biomodal(), ",
    "and save_multiformat_ggplot())."
  )
}

# patchwork supplies plot_annotation()/wrap_elements() used by the figure
# scripts; it is listed as a _shared_config.R dependency but require it here so
# add_panel_labels() fails loudly rather than silently producing a non-tagged
# plot if the package is unavailable.
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("_figure_config.R requires the 'patchwork' package for panel assembly.")
}

# =============================================================================
# 1. theme_pub(): publication theme extending theme_biomodal()
# =============================================================================

#' Publication-grade ggplot2 theme for consolidated figures.
#'
#' Extends theme_biomodal() with a smaller default base size, a boxed panel
#' border, thin axis ticks, compact legend keys, and a bold patchwork panel
#' tag. Helvetica/Arial is requested as the base family for journal output;
#' R silently falls back to the default sans serif if the font is unavailable,
#' which keeps the script portable across machines.
#'
#' @param base_size Base font size in points (default 8 for multi-panel print).
#' @return A ggplot2 theme object.
theme_pub <- function(base_size = 8) {
  theme_biomodal(base_size = base_size) +
    ggplot2::theme(
      text             = ggplot2::element_text(family = "Helvetica"),
      panel.border     = ggplot2::element_rect(colour = "black",
                                               fill = NA, linewidth = 0.5),
      axis.ticks       = ggplot2::element_line(linewidth = 0.3),
      legend.key.size  = ggplot2::unit(3, "mm"),
      plot.tag         = ggplot2::element_text(face = "bold", size = 12),
      strip.background = ggplot2::element_blank()
    )
}

# =============================================================================
# 2. FIGURE_DIR: output directory for consolidated figures (auto-created)
# =============================================================================

# Consolidated sections output to the same plots/visualizations/ tree as all
# other sections. TABLES_DIR is already created by _shared_config.R.
FIGURE_DIR        <- OUTPUT_DIR
FIGURE_TABLES_DIR <- TABLES_DIR

# =============================================================================
# 3. save_figure(): mm-dimensioned wrapper around save_multiformat_ggplot()
# =============================================================================

#' Save a consolidated figure in PDF + SVG + JPG under plots/figures/<name>/.
#'
#' Journals specify figure dimensions in millimetres; save_multiformat_ggplot()
#' (and the underlying ggsave) expect inches. This wrapper converts mm -> in,
#' fixes 300 DPI, and routes output into a per-figure subfolder under
#' FIGURE_DIR via the utility's default use_subfolders = TRUE behaviour. The
#' subfolder and three format files are named after `name`
#' (e.g. plots/figures/figure1_methylation_phenotype/figure1_methylation_phenotype.{pdf,svg,jpg}).
#'
#' @param plot      A ggplot/patchwork object.
#' @param name      Figure basename WITHOUT extension (becomes the subfolder).
#' @param width_mm  Figure width in millimetres.
#' @param height_mm Figure height in millimetres.
#' @return Invisibly returns the plot object (passthrough from the utility).
save_figure <- function(plot, name, width_mm, height_mm) {
  if (is.null(plot)) {
    stop("save_figure(): 'plot' is NULL for figure '", name, "'.")
  }
  mm_to_in <- 1 / 25.4
  base_path <- file.path(FIGURE_DIR, name)
  save_multiformat_ggplot(
    plot          = plot,
    base_path     = base_path,
    width         = width_mm  * mm_to_in,
    height        = height_mm * mm_to_in,
    dpi           = 300,
    verbose       = TRUE,
    use_subfolders = TRUE
  )
}

# =============================================================================
# 4. read_table_tsv(): load a pre-computed section TSV from TABLES_DIR
# =============================================================================

#' Read a pre-computed section table from plots/visualizations/tables/.
#'
#' The 255 section TSVs are the canonical data source for every figure panel.
#' This reader uses strict, lossless options: tab-separated, no quoting (some
#' tables contain unbalanced quote characters in gene/annotation fields),
#' character-preserving, and verbatim column names (check.names = FALSE keeps
#' names like '5mC_pct' intact). Fails loudly if the file is absent.
#'
#' @param filename Bare filename (e.g. "64_global_methylation_summary.tsv").
#' @return A data.frame.
read_table_tsv <- function(filename) {
  path <- file.path(TABLES_DIR, filename)
  if (!file.exists(path)) {
    stop("read_table_tsv(): table not found: ", path)
  }
  utils::read.table(
    path,
    header           = TRUE,
    sep              = "\t",
    stringsAsFactors = FALSE,
    quote            = "",
    check.names      = FALSE,
    comment.char     = ""
  )
}

# =============================================================================
# 5. add_panel_labels(): patchwork panel lettering (A, B, C, ...)
# =============================================================================

#' Add bold A/B/C... panel tags to a patchwork composite.
#'
#' The tag appearance (bold, size 12) is inherited from theme_pub()'s
#' plot.tag; this helper only assigns the tag *sequence*.
#'
#' @param plot A patchwork composite (or a ggplot wrapped via wrap_plots()).
#' @return The plot annotated with tag_levels = "A".
add_panel_labels <- function(plot) {
  plot + patchwork::plot_annotation(tag_levels = "A")
}

# =============================================================================
# 6. stat_label(): formatted statistic annotation helper
# =============================================================================

#' Format a single statistic for on-panel annotation.
#'
#' Produces strings like "OR = 4.71", "AUC = 0.793", or "slope = -0.96" with
#' consistent numeric formatting across all figures. P-values are best passed
#' with fmt = "%.2e" (e.g. stat_label("p", 9.43e-49, "%.2e")). Every value
#' supplied here MUST match the corresponding source TSV.
#'
#' @param name  Statistic label (e.g. "OR", "AUC", "p", "slope").
#' @param value Numeric value to format.
#' @param fmt   sprintf format string for the value (default "%.2f").
#' @return A single character string "name = value".
stat_label <- function(name, value, fmt = "%.2f") {
  sprintf("%s = %s", name, sprintf(fmt, value))
}

# =============================================================================
# 7. PUB_COLORS: canonical colour re-exports for cross-figure consistency
# =============================================================================

# Re-export the load-bearing colours from COLORS so every figure script uses
# identical hex values without re-typing them. Per PLAN.md verification:
# 5mC must always be #E41A1C and 5hmC #377EB8 (the COLORS$methylation key,
# NOT the COLORS$direction key). Hyper/Hypo come from COLORS$direction;
# Control/Mutant from COLORS$condition.
PUB_COLORS <- list(
  mC      = COLORS$methylation[["5mC"]],            # "#E41A1C"
  hmC     = COLORS$methylation[["5hmC"]],           # "#377EB8"
  Hyper   = COLORS$direction[["Hypermethylated"]],  # "#D7191C"
  Hypo    = COLORS$direction[["Hypomethylated"]],   # "#2C7BB6"
  Control = COLORS$condition[["Control"]],          # "#2166AC"
  Mutant  = COLORS$condition[["Mutant"]]            # "#B2182B"
)

# -----------------------------------------------------------------------------
# Sourced-confirmation message (mirrors multi_format_output.R convention).
# -----------------------------------------------------------------------------
cat("Loaded _figure_config.R - consolidated section utilities available:\n")
cat("  - theme_pub(base_size = 8)\n")
cat("  - FIGURE_DIR:", FIGURE_DIR, "\n")
cat("  - save_figure(plot, name, width_mm, height_mm)  [PDF + SVG + JPG]\n")
cat("  - read_table_tsv(filename)  [from", TABLES_DIR, "]\n")
cat("  - add_panel_labels(plot), stat_label(name, value, fmt), PUB_COLORS\n")
