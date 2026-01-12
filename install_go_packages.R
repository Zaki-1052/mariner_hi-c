# install_go_packages.R
# Install packages needed for GO enrichment analysis

cat("Installing required packages for GO enrichment...\n")

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Install clusterProfiler and org.Mm.eg.db
BiocManager::install(c("clusterProfiler", "org.Mm.eg.db"),
                     update = FALSE,
                     ask = FALSE)

cat("Installation complete.\n")
