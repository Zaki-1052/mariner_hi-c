# ML/cmpts/scripts/utils/install_calder2_deps.R
# Installs all CALDER2 dependencies and the package itself from local source.
# Called by A0_setup_calder2_env.sh after conda env creation.
#
# Usage: Rscript install_calder2_deps.R <path_to_calder2_source>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript install_calder2_deps.R <calder2_source_dir>")

calder2_src <- args[1]
if (!file.exists(file.path(calder2_src, "DESCRIPTION"))) {
    stop(paste("CALDER2 DESCRIPTION not found at:", calder2_src))
}

cat("── Installing BiocManager + Bioconductor packages ──\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install(c("GenomicRanges", "rhdf5"), update = FALSE, ask = FALSE)

cat("\n── Installing CRAN dependencies ──\n")
cran_pkgs <- c(
    "strawr", "data.table", "ape", "dendextend", "fitdistrplus",
    "igraph", "Matrix", "rARPACK", "factoextra", "fields",
    "ggplot2", "optparse", "R.utils", "doParallel", "foreach"
)
install.packages(cran_pkgs, repos = "https://cloud.r-project.org")

cat("\n── Installing CALDER2 from local source ──\n")
install.packages(calder2_src, repos = NULL, type = "source")

cat("\n── Verification ──\n")
library(CALDER)
cat("CALDER version:", as.character(packageVersion("CALDER")), "\n")
ref <- system.file("extdata", "mm10_all_sub_compartments.bed", package = "CALDER")
cat("mm10 reference:", ref, "\n")
