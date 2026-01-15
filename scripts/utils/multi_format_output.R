# scripts/utils/multi_format_output.R
# Utility functions for multi-format plot output (PDF + SVG + JPEG)
# Author: Zakir Alibhai
# Date: 2026-01-14
#
# Purpose:
#   Provides helper functions to output plots in three formats simultaneously:
#   - PDF (publication standard, vector)
#   - SVG (Illustrator-friendly, editable vector)
#   - JPEG (Google Slides, presentations)
#
# Usage:
#   source("scripts/utils/multi_format_output.R")
#
#   # For ggplot2 objects:
#   save_multiformat_ggplot(my_plot, "path/to/output/filename", width = 10, height = 8)
#
#   # For base R graphics:
#   save_multiformat_base(quote({ plot(x, y); lines(x, y) }), "path/to/output/filename", width = 10, height = 8)
#
# Note: Requires svglite package for high-quality SVG output
#   install.packages("svglite")

# Load svglite for Illustrator-friendly SVG output
if (!requireNamespace("svglite", quietly = TRUE)) {
  message("Installing svglite package for SVG output...")
  install.packages("svglite", repos = "https://cloud.r-project.org")
}
library(svglite)

#' Save ggplot object in multiple formats (PDF, SVG, JPEG)
#'
#' @param plot ggplot object to save
#' @param base_path Output path WITHOUT file extension (e.g., "outputs/volcano" not "outputs/volcano.pdf")
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 8)
#' @param dpi DPI for JPEG output (default: 300 for print quality)
#' @param verbose Print confirmation messages (default: TRUE)
#'
#' @return Invisibly returns the plot object
#'
#' @examples
#' p <- ggplot(mtcars, aes(mpg, hp)) + geom_point()
#' save_multiformat_ggplot(p, "outputs/my_plot", width = 8, height = 6)
#' # Creates: my_plot.pdf, my_plot.svg, my_plot.jpg
save_multiformat_ggplot <- function(plot, base_path, width = 10, height = 8, dpi = 300, verbose = TRUE) {
  # Ensure output directory exists
  output_dir <- dirname(base_path)
  if (!dir.exists(output_dir) && output_dir != ".") {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # PDF (publication standard)
  pdf_path <- paste0(base_path, ".pdf")
  ggplot2::ggsave(pdf_path, plot, width = width, height = height)

 # SVG (Illustrator-friendly with editable text)
  svg_path <- paste0(base_path, ".svg")
  ggplot2::ggsave(svg_path, plot, width = width, height = height, device = svglite::svglite)

  # JPEG (presentations/slides)
  jpg_path <- paste0(base_path, ".jpg")
  ggplot2::ggsave(jpg_path, plot, width = width, height = height, dpi = dpi, device = "jpeg")

  if (verbose) {
    cat(sprintf("  Saved: %s.{pdf,svg,jpg}\n", basename(base_path)))
  }

  invisible(plot)
}


#' Save base R graphics in multiple formats (PDF, SVG, JPEG)
#'
#' @param plot_expr A quoted expression containing the plotting code.
#'   Use quote({...}) to wrap multiple plotting commands.
#' @param base_path Output path WITHOUT file extension
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 8)
#' @param dpi DPI for JPEG output (default: 300)
#' @param verbose Print confirmation messages (default: TRUE)
#'
#' @return NULL (invisibly)
#'
#' @examples
#' save_multiformat_base(
#'   quote({
#'     plot(1:10, main = "Test Plot")
#'     abline(h = 5, col = "red")
#'   }),
#'   "outputs/base_plot",
#'   width = 8, height = 6
#' )
save_multiformat_base <- function(plot_expr, base_path, width = 10, height = 8, dpi = 300, verbose = TRUE) {
  # Ensure output directory exists
  output_dir <- dirname(base_path)
  if (!dir.exists(output_dir) && output_dir != ".") {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # PDF
  pdf_path <- paste0(base_path, ".pdf")
  pdf(pdf_path, width = width, height = height)
  tryCatch({
    eval(plot_expr)
  }, finally = {
    dev.off()
  })

  # SVG (using svglite for better Illustrator compatibility)
  svg_path <- paste0(base_path, ".svg")
  svglite::svglite(svg_path, width = width, height = height)
  tryCatch({
    eval(plot_expr)
  }, finally = {
    dev.off()
  })

  # JPEG
  jpg_path <- paste0(base_path, ".jpg")
  jpeg(jpg_path, width = width * dpi, height = height * dpi, res = dpi, quality = 95)
  tryCatch({
    eval(plot_expr)
  }, finally = {
    dev.off()
  })

  if (verbose) {
    cat(sprintf("  Saved: %s.{pdf,svg,jpg}\n", basename(base_path)))
  }

  invisible(NULL)
}


#' Save pheatmap in multiple formats (PDF, SVG, JPEG)
#'
#' pheatmap is special because it draws directly and doesn't return a ggplot object.
#' This wrapper handles pheatmap's unique behavior.
#'
#' @param pheatmap_call A quoted pheatmap() call expression
#' @param base_path Output path WITHOUT file extension
#' @param width Plot width in inches (default: 8)
#' @param height Plot height in inches (default: 10)
#' @param dpi DPI for JPEG output (default: 300)
#' @param verbose Print confirmation messages (default: TRUE)
#'
#' @return NULL (invisibly)
#'
#' @examples
#' save_multiformat_pheatmap(
#'   quote(pheatmap(mat, main = "Heatmap", cluster_rows = TRUE)),
#'   "outputs/heatmap",
#'   width = 8, height = 10
#' )
save_multiformat_pheatmap <- function(pheatmap_call, base_path, width = 8, height = 10, dpi = 300, verbose = TRUE) {
  # Ensure output directory exists
  output_dir <- dirname(base_path)
  if (!dir.exists(output_dir) && output_dir != ".") {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # PDF
  pdf_path <- paste0(base_path, ".pdf")
  pdf(pdf_path, width = width, height = height)
  tryCatch({
    eval(pheatmap_call)
  }, finally = {
    dev.off()
  })

  # SVG
  svg_path <- paste0(base_path, ".svg")
  svglite::svglite(svg_path, width = width, height = height)
  tryCatch({
    eval(pheatmap_call)
  }, finally = {
    dev.off()
  })

  # JPEG
  jpg_path <- paste0(base_path, ".jpg")
  jpeg(jpg_path, width = width * dpi, height = height * dpi, res = dpi, quality = 95)
  tryCatch({
    eval(pheatmap_call)
  }, finally = {
    dev.off()
  })

  if (verbose) {
    cat(sprintf("  Saved: %s.{pdf,svg,jpg}\n", basename(base_path)))
  }

  invisible(NULL)
}


#' Get the path to the multi_format_output.R utility
#'
#' Helper to find this file from any script location
#' @return Character string with the path to this file
get_multiformat_util_path <- function() {
  # Try multiple locations
  candidates <- c(
    "scripts/utils/multi_format_output.R",
    "../scripts/utils/multi_format_output.R",
    "../../scripts/utils/multi_format_output.R",
    file.path(Sys.getenv("MARINER_BASE"), "scripts/utils/multi_format_output.R")
  )

  for (path in candidates) {
    if (file.exists(path)) {
      return(normalizePath(path))
    }
  }

  stop("Could not find multi_format_output.R utility. Set MARINER_BASE environment variable.")
}


# Print confirmation when sourced
cat("Loaded multi_format_output.R - Functions available:\n")
cat("  - save_multiformat_ggplot(plot, base_path, width, height, dpi)\n")
cat("  - save_multiformat_base(plot_expr, base_path, width, height, dpi)\n")
cat("  - save_multiformat_pheatmap(pheatmap_call, base_path, width, height, dpi)\n")
