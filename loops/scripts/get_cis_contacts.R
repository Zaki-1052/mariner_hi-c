# scripts/get_cis_contacts.R
# Compute total valid cis contacts per sample from .hic files
# Uses strawr::straw() at 2.5Mb resolution to sum all intrachromosomal contacts
# Output: config/cis_contacts.yaml

suppressPackageStartupMessages({
  library(strawr)
  library(yaml)
})

cat("\n========================================\n")
cat("Valid Cis Contact Extraction\n")
cat("========================================\n\n")

config <- yaml::read_yaml("config/paths_config.yaml")

hic_files <- c(
  ctrl_M1 = config$hic_files$ctrl_M1,
  ctrl_M2 = config$hic_files$ctrl_M2,
  ctrl_M3 = config$hic_files$ctrl_M3,
  mut_M1 = config$hic_files$mut_M1,
  mut_M2 = config$hic_files$mut_M2,
  mut_M3 = config$hic_files$mut_M3
)

cat("Verifying .hic files...\n")
for (name in names(hic_files)) {
  if (!file.exists(hic_files[name])) {
    stop(sprintf("Missing .hic file for %s: %s", name, hic_files[name]))
  }
}
cat(sprintf("  All %d files found\n\n", length(hic_files)))

coarse_res <- 2500000
cis_contacts <- numeric(length(hic_files))
names(cis_contacts) <- names(hic_files)

for (sample_name in names(hic_files)) {
  hic_path <- hic_files[sample_name]
  cat(sprintf("Processing %s...\n", sample_name))

  chroms <- readHicChroms(hic_path)
  standard <- chroms[grep("^chr[0-9XY]+$", chroms$name), ]
  cat(sprintf("  %d standard chromosomes\n", nrow(standard)))

  total <- 0
  for (chrom in standard$name) {
    contacts <- straw("NONE", hic_path, chrom, chrom, "BP", coarse_res)
    chrom_total <- sum(contacts$counts)
    total <- total + chrom_total
  }

  cis_contacts[sample_name] <- total
  cat(sprintf("  Valid cis contacts: %.0f\n", total))
}

cat("\n========================================\n")
cat("Summary\n")
cat("========================================\n")
for (name in names(cis_contacts)) {
  cat(sprintf("  %s: %15.0f\n", name, cis_contacts[name]))
}

ratio <- max(cis_contacts) / min(cis_contacts)
cat(sprintf("\n  Max/min ratio: %.2f\n", ratio))

if (any(cis_contacts == 0)) {
  stop("One or more samples has zero valid cis contacts")
}

output <- list(cis_contacts = as.list(cis_contacts))
yaml::write_yaml(output, "config/cis_contacts.yaml")
cat(sprintf("\nSaved to: config/cis_contacts.yaml\n"))
