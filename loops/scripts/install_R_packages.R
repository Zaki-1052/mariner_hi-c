#!/usr/bin/env Rscript
# install_R_packages.R
# Comprehensive R package installation script for mariner pipeline
# Author: Zakir Alibhai
# Date: 2025-11-10

cat("\n========================================\n")
cat("Mariner Pipeline: R Package Installation\n")
cat("========================================\n\n")

# =============================================================================
# 1. BIOCONDUCTOR SETUP
# =============================================================================

cat("Step 1: Installing BiocManager...\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
  cat("  ✓ BiocManager installed\n\n")
} else {
  cat("  ✓ BiocManager already installed\n\n")
}

# Set Bioconductor version (3.18 for R 4.3.x, 3.19 for R 4.4.x)
BiocManager::install(version = "3.19", ask = FALSE, update = FALSE)

# =============================================================================
# 2. CRAN PACKAGES
# =============================================================================

cat("Step 2: Installing CRAN packages...\n\n")

cran_packages <- c(
  "yaml",              # Configuration file parsing (edgeR.R)
  "ggplot2",           # Core plotting (all visualization scripts)
  "dplyr",             # Data manipulation (edgeR.R, visualizations.R)
  "tibble",            # Modern data frames (edgeR.R)
  "tidyr",             # Data tidying (compare_resolutions.R, downstream_analysis.R)
  "pheatmap",          # Heatmap plotting (compare_resolutions.R, visualizations.R)
  "patchwork",         # Plot composition (compare_resolutions.R, downstream_analysis.R)
  "RColorBrewer",      # Color palettes (visualizations.R, apa_analysis.R)
  "viridis",           # Perceptually uniform color scales (apa_analysis.R)
  "VennDiagram",       # Venn diagram generation (compare_resolutions.R)
  "Matrix",            # Sparse matrix operations (aggregate.R)
  "scales",            # Scale functions for ggplot2 (visualizations.R)
  "igraph",            # Graph construction and layout (network_analysis.R)
  "ggraph",            # ggplot2-based network visualization (network_analysis.R)
  "tidygraph",         # Tidy graph manipulation (network_analysis.R)
  "ggforce",           # Extended ggplot2 geoms incl. mark_hull (network_analysis.R)
  "ggnewscale"         # Multiple fill/color scales in one plot (network_analysis.R)
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", dependencies = TRUE)
  } else {
    cat(sprintf("  ✓ %s already installed\n", pkg))
  }
}

cat("\n✓ CRAN packages complete\n\n")

# =============================================================================
# 3. BIOCONDUCTOR CORE INFRASTRUCTURE
# =============================================================================

cat("Step 3: Installing Bioconductor core infrastructure...\n\n")

bioc_core <- c(
  "BiocGenerics",      # S4 generic functions
  "S4Vectors",         # S4 vector infrastructure
  "IRanges",           # Integer range infrastructure
  "GenomeInfoDb",      # Genome metadata
  "GenomicRanges",     # Genomic ranges (all scripts)
  "InteractionSet",    # Hi-C interaction data (all scripts)
  "SummarizedExperiment",  # Experiment data containers (apa_analysis.R)
  "DelayedArray",      # Delayed array operations (extract_counts.R, aggregate.R)
  "HDF5Array"          # HDF5-backed arrays (extract_counts.R, aggregate.R, apa_analysis.R)
)

for (pkg in bioc_core) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    cat(sprintf("  ✓ %s already installed\n", pkg))
  }
}

cat("\n✓ Bioconductor core complete\n\n")

# =============================================================================
# 4. BIOCONDUCTOR HI-C ANALYSIS
# =============================================================================

cat("Step 4: Installing Hi-C analysis packages...\n\n")

bioc_hic <- c(
  "mariner",           # Hi-C loop analysis (all pipeline stages)
  "strawr"             # Straw API for .hic files (extract_counts.R, apa_analysis.R)
)

for (pkg in bioc_hic) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    cat(sprintf("  ✓ %s already installed\n", pkg))
  }
}

cat("\n✓ Hi-C packages complete\n\n")

# =============================================================================
# 5. BIOCONDUCTOR DIFFERENTIAL ANALYSIS
# =============================================================================

cat("Step 5: Installing differential analysis packages...\n\n")

bioc_diff <- c(
  "edgeR",             # Differential expression/loop analysis (aggregate.R, edgeR.R)
  "limma"              # Linear models (dependency for edgeR)
)

for (pkg in bioc_diff) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    cat(sprintf("  ✓ %s already installed\n", pkg))
  }
}

cat("\n✓ Differential analysis packages complete\n\n")

# =============================================================================
# 6. BIOCONDUCTOR GENOME ANNOTATION
# =============================================================================

cat("Step 6: Installing genome annotation packages...\n\n")

bioc_anno <- c(
  "rtracklayer",       # Import/export genomic data (downstream_analysis.R)
  "GenomicFeatures",   # Gene/transcript models (downstream_analysis.R)
  "TxDb.Mmusculus.UCSC.mm10.knownGene",  # Mouse mm10 gene annotations
  "org.Mm.eg.db",      # Mouse gene ID mappings
  "AnnotationDbi"      # Annotation database interface
)

for (pkg in bioc_anno) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    cat(sprintf("  ✓ %s already installed\n", pkg))
  }
}

cat("\n✓ Genome annotation packages complete\n\n")

# =============================================================================
# 7. BIOCONDUCTOR FUNCTIONAL ENRICHMENT
# =============================================================================

cat("Step 7: Installing functional enrichment packages...\n\n")

bioc_enrich <- c(
  "ChIPseeker",        # ChIP-seq peak annotation (visualizations.R)
  "clusterProfiler",   # GO/KEGG enrichment analysis (visualizations.R)
  "enrichplot",        # Enrichment visualization (visualizations.R)
  "DOSE",              # Disease ontology enrichment (visualizations.R)
  "GO.db",             # Gene Ontology database
  "KEGG.db"            # KEGG pathway database
)

for (pkg in bioc_enrich) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    cat(sprintf("  ✓ %s already installed\n", pkg))
  }
}

cat("\n✓ Functional enrichment packages complete\n\n")

# =============================================================================
# 8. ENHANCED VISUALIZATION
# =============================================================================

cat("Step 8: Installing enhanced visualization packages...\n\n")

# EnhancedVolcano may be on Bioconductor or GitHub
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  cat("  Installing EnhancedVolcano...\n")
  BiocManager::install("EnhancedVolcano", ask = FALSE, update = FALSE)
} else {
  cat("  ✓ EnhancedVolcano already installed\n")
}

cat("\n✓ Enhanced visualization complete\n\n")

# =============================================================================
# 9. VERIFICATION
# =============================================================================

cat("\n========================================\n")
cat("Installation Verification\n")
cat("========================================\n\n")

all_packages <- c(
  cran_packages,
  bioc_core,
  bioc_hic,
  bioc_diff,
  bioc_anno,
  bioc_enrich,
  "EnhancedVolcano"
)

missing_packages <- character()

for (pkg in all_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
    cat(sprintf("  ✗ %s - NOT INSTALLED\n", pkg))
  } else {
    cat(sprintf("  ✓ %s\n", pkg))
  }
}

cat("\n========================================\n")

if (length(missing_packages) == 0) {
  cat("✓ All packages installed successfully!\n")
  cat("========================================\n\n")
  cat("You can now run the mariner pipeline:\n")
  cat("  Rscript scripts/prep_loops.R 5000\n")
  cat("  Rscript scripts/extract_counts.R 5000\n")
  cat("  Rscript scripts/aggregate.R 5000\n")
  cat("  Rscript scripts/edgeR.R 5000\n")
  cat("  Rscript scripts/compare_resolutions.R\n")
  cat("  Rscript scripts/downstream_analysis.R\n")
  cat("  Rscript scripts/visualizations.R\n\n")
} else {
  cat("✗ Installation incomplete\n")
  cat("========================================\n\n")
  cat("Missing packages:\n")
  for (pkg in missing_packages) {
    cat(sprintf("  - %s\n", pkg))
  }
  cat("\nPlease install missing packages manually:\n")
  cat("  BiocManager::install(c(")
  cat(paste(sprintf('"%s"', missing_packages), collapse = ", "))
  cat("))\n\n")
}

# =============================================================================
# 10. SAVE SESSION INFO
# =============================================================================

session_file <- "R_package_installation_session_info.txt"
cat(sprintf("Saving session info to: %s\n", session_file))

sink(session_file)
cat("========================================\n")
cat("R Package Installation Session Info\n")
cat("========================================\n")
cat(sprintf("Date: %s\n\n", Sys.time()))
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("BiocManager version: %s\n\n", packageVersion("BiocManager")))

cat("Installed package versions:\n")
cat("---------------------------\n\n")

for (pkg in sort(all_packages)) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("%-40s %s\n", pkg, packageVersion(pkg)))
  }
}

cat("\n========================================\n")
cat("Full sessionInfo():\n")
cat("========================================\n\n")
print(sessionInfo())
sink()

cat(sprintf("\n✓ Session info saved: %s\n\n", session_file))

cat("========================================\n")
cat("Installation Complete\n")
cat("========================================\n\n")
