# scripts/print_versions.R
# Print installed versions of all R/Bioconductor packages used in the
# loops and stripes pipelines.  Run from any directory:
#   Rscript scripts/print_versions.R 2>&1 | tee software_versions.txt

cat("=== R session info ===\n")
cat("R version:", paste(R.version$major, R.version$minor, sep = "."), "\n")
cat("Platform:", R.version$platform, "\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

pkgs <- c(
  # Bioconductor core (loops pipeline)
  "mariner", "edgeR", "limma", "InteractionSet", "GenomicRanges",
  "IRanges", "GenomicFeatures", "SummarizedExperiment",
  "HDF5Array", "DelayedArray", "strawr",
  # Annotation
  "ChIPseeker", "clusterProfiler", "enrichplot", "DOSE",
  "TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db",
  "AnnotationDbi", "rtracklayer",
  # Visualization
  "EnhancedVolcano", "ggplot2", "patchwork", "pheatmap",
  "viridis", "RColorBrewer", "scales", "svglite", "ggrepel",
  # Data
  "tidyverse", "dplyr", "tidyr", "readr", "readxl",
  "yaml", "jsonlite", "gtools",
  # Permutation / region analysis
  "regioneReloaded", "VennDiagram",
  # BiocManager itself
  "BiocManager"
)

cat("=== Package versions ===\n")
cat(sprintf("%-45s %s\n", "Package", "Version"))
cat(paste(rep("-", 60), collapse = ""), "\n")

for (pkg in sort(pkgs)) {
  v <- tryCatch(
    as.character(packageVersion(pkg)),
    error = function(e) "NOT INSTALLED"
  )
  cat(sprintf("%-45s %s\n", pkg, v))
}

cat("\n=== Bioconductor version ===\n")
tryCatch(
  cat("BiocManager::version():", as.character(BiocManager::version()), "\n"),
  error = function(e) cat("BiocManager not available\n")
)

cat("\n=== sessionInfo() ===\n")
print(sessionInfo())
