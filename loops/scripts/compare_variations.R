# scripts/compare_variations.R
# Compare edgeR results across normalization and buffer size variations

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

option_list <- list(
  make_option("--output-dir", type = "character",
              default = "outputs/250402-late_outputs/variations/comparison",
              help = "Output directory for comparison results")
)
parsed <- parse_args(OptionParser(option_list = option_list))
output_dir <- parsed[["output-dir"]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

config <- yaml::read_yaml("config/paths_config.yaml")
base_dir <- config$project$base_dir
setwd(base_dir)

cat("\n========================================\n")
cat("Variation Comparison\n")
cat("========================================\n\n")

# Define all conditions to compare
conditions <- list(
  list(name = "baseline_TMM_5x5",
       resolutions = c(5, 10, 25),
       base = "outputs/250402-late_outputs"),
  list(name = "cisnorm_5x5",
       resolutions = c(5, 10, 25),
       base = "outputs/250402-late_outputs/variations/cisnorm"),
  list(name = "buffer_1x1_TMM",
       resolutions = c(25),
       base = "outputs/250402-late_outputs/variations/buffer_1x1"),
  list(name = "buffer_3x3_TMM",
       resolutions = c(25),
       base = "outputs/250402-late_outputs/variations/buffer_3x3")
)

results <- data.frame(
  condition = character(),
  resolution = character(),
  loops_tested = integer(),
  sig_fdr05 = integer(),
  sig_strict = integer(),
  up_fdr05 = integer(),
  down_fdr05 = integer(),
  median_abslogfc_sig = numeric(),
  stringsAsFactors = FALSE
)

for (cond in conditions) {
  for (res_kb in cond$resolutions) {
    results_file <- file.path(cond$base,
                              sprintf("edgeR_results_res_%dkb/primary_analysis/all_results_primary.tsv", res_kb))

    if (!file.exists(results_file)) {
      cat(sprintf("  SKIP: %s %dkb — file not found: %s\n", cond$name, res_kb, results_file))
      next
    }

    df <- read.table(results_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    sig05 <- df[!is.na(df$FDR) & df$FDR < 0.05, ]
    sig_strict <- df[!is.na(df$FDR) & df$FDR < 0.03 & abs(df$logFC) > 0.3, ]
    up <- sig05[sig05$logFC > 0, ]
    down <- sig05[sig05$logFC < 0, ]

    results <- rbind(results, data.frame(
      condition = cond$name,
      resolution = sprintf("%dkb", res_kb),
      loops_tested = nrow(df),
      sig_fdr05 = nrow(sig05),
      sig_strict = nrow(sig_strict),
      up_fdr05 = nrow(up),
      down_fdr05 = nrow(down),
      median_abslogfc_sig = if (nrow(sig05) > 0) median(abs(sig05$logFC)) else NA,
      stringsAsFactors = FALSE
    ))

    cat(sprintf("  %s %dkb: %d tested, %d sig (FDR<0.05), %d strict, %d up / %d down\n",
                cond$name, res_kb, nrow(df), nrow(sig05), nrow(sig_strict),
                nrow(up), nrow(down)))
  }
}

cat("\n========================================\n")
cat("Comparison Table\n")
cat("========================================\n\n")
print(results, row.names = FALSE)

out_file <- file.path(output_dir, "variation_comparison.tsv")
write.table(results, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nSaved to: %s\n", out_file))

# Overlap analysis between baseline and cis-norm (25kb)
cat("\n========================================\n")
cat("Overlap Analysis: Baseline vs Cis-norm (all resolutions)\n")
cat("========================================\n\n")

for (res_kb in c(5, 10, 25)) {
  baseline_file <- sprintf("outputs/250402-late_outputs/edgeR_results_res_%dkb/primary_analysis/all_results_primary.tsv", res_kb)
  cisnorm_file <- sprintf("outputs/250402-late_outputs/variations/cisnorm/edgeR_results_res_%dkb/primary_analysis/all_results_primary.tsv", res_kb)

  if (!file.exists(baseline_file) || !file.exists(cisnorm_file)) {
    cat(sprintf("  SKIP %dkb — missing file\n", res_kb))
    next
  }

  bl <- read.table(baseline_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  cn <- read.table(cisnorm_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  bl_sig <- bl$FDR < 0.05 & !is.na(bl$FDR)
  cn_sig <- cn$FDR < 0.05 & !is.na(cn$FDR)

  shared <- sum(bl_sig & cn_sig)
  bl_only <- sum(bl_sig & !cn_sig)
  cn_only <- sum(!bl_sig & cn_sig)

  cat(sprintf("  %dkb: shared=%d, baseline-only=%d, cisnorm-only=%d\n",
              res_kb, shared, bl_only, cn_only))

  # Direction concordance among shared
  if (shared > 0) {
    concordant <- sum(sign(bl$logFC[bl_sig & cn_sig]) == sign(cn$logFC[bl_sig & cn_sig]))
    cat(sprintf("    Direction concordance: %d/%d (%.1f%%)\n",
                concordant, shared, 100 * concordant / shared))
  }
}

cat("\nDone.\n")
