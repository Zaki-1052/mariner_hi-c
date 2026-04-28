#!/usr/bin/env Rscript
# /expanse/lustre/projects/csd940/zalibhai/mariner/scripts/diagnose_extraction.R
# Purpose: Diagnose pullHicMatrices function signature and parameters

cat("\n=== Diagnosing pullHicMatrices Function ===\n\n")

# Load libraries
library(mariner)
library(InteractionSet)

# Check mariner version
cat("Mariner version:", as.character(packageVersion("mariner")), "\n\n")

# Check function signature
cat("Function signature for pullHicMatrices:\n")
cat("----------------------------------------\n")
args_list <- args(pullHicMatrices)
print(args_list)

# Check formal arguments
cat("\nFormal arguments:\n")
formals_list <- formals(pullHicMatrices)
cat("Available parameters:", paste(names(formals_list), collapse=", "), "\n\n")

# Check if specific parameters exist
params_to_check <- c("outDir", "BPPARAM", "blockSize", "onDisk", "norm", "matrix")
cat("Parameter availability check:\n")
for (param in params_to_check) {
  exists_param <- param %in% names(formals_list)
  cat(sprintf("  %s: %s\n", param, ifelse(exists_param, "✓ EXISTS", "✗ NOT FOUND")))
}

# Check methods available for pullHicMatrices
cat("\n\nMethods for pullHicMatrices:\n")
methods_list <- methods("pullHicMatrices")
if (length(methods_list) > 0) {
  print(methods_list)
} else {
  cat("No specific methods found (using default)\n")
}

# Try to get help documentation programmatically
cat("\n\nSearching for documentation...\n")
help_obj <- help("pullHicMatrices", package = "mariner")
if (length(help_obj) > 0) {
  cat("Documentation found - run ?pullHicMatrices in R for details\n")
}

# Load test data to try a minimal extraction
cat("\n\nTesting minimal extraction with discovered parameters...\n")
buffered <- readRDS("outputs/test/04_buffered.rds")

# Get just first 5 loops for ultra-quick test
test_loops <- buffered[1:5]

# Define files
hicFiles <- c(
  ctrl = "/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/merged/hic/resorted_ctrl.hic",
  mut = "/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/merged/hic/resorted_mut.hic"
)

# Try extraction with minimal parameters
cat("Attempting extraction with basic parameters only...\n")
tryCatch({
  # Try most basic call
  result1 <- pullHicMatrices(
    x = test_loops,
    files = hicFiles,
    binSize = 5000
  )
  cat("✓ Basic extraction successful!\n")
  cat("  Result class:", class(result1), "\n")
  cat("  Dimensions:", paste(dim(counts(result1)), collapse=" x "), "\n")
}, error = function(e) {
  cat("✗ Basic extraction failed:", e$message, "\n")
})

# Try with normalization
cat("\nTrying with VC normalization...\n")
tryCatch({
  result2 <- pullHicMatrices(
    x = test_loops,
    files = hicFiles,
    binSize = 5000,
    norm = "VC"
  )
  cat("✓ Extraction with norm='VC' successful!\n")
}, error = function(e) {
  cat("✗ Extraction with norm failed:", e$message, "\n")
})

# Try with other potential parameters
cat("\nTrying with blockSize...\n")
tryCatch({
  if ("blockSize" %in% names(formals_list)) {
    result3 <- pullHicMatrices(
      x = test_loops,
      files = hicFiles,
      binSize = 5000,
      blockSize = 5
    )
    cat("✓ blockSize parameter works!\n")
  } else {
    cat("- blockSize not in function signature\n")
  }
}, error = function(e) {
  cat("✗ blockSize caused error:", e$message, "\n")
})

cat("\n=== Diagnosis Complete ===\n")
