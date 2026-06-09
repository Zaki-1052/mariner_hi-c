# biomodal/downstream/scripts/prepare_mecp2_regions.R
# Prepare MeCP2 peaks as a modality-compatible regions BED file
# for running non-CG (CHG/CHH) feature extraction at MeCP2 binding sites.
#
# Modality Regions_Directory format:
#   Chromosome\tStart\tEnd\tAnnotation\tName (no chr prefix, gzipped)
#
# Run from downstream/ directory:
#   Rscript scripts/prepare_mecp2_regions.R
#
# Then scp the output to HPC:
#   scp modality/mecp2_regions/mecp2_peaks.annotation.bed.gz \
#       expanse:/expanse/lustre/projects/csd940/zalibhai/biomodal/modality/mecp2_regions/

REPO_ROOT <- normalizePath(file.path(getwd(), "../.."))

mecp2_up_path <- file.path(REPO_ROOT, "peaks/mecp2/MeCP2_up.bed")
mecp2_down_path <- file.path(REPO_ROOT, "peaks/mecp2/MeCP2_down.bed")

stopifnot(file.exists(mecp2_up_path))
stopifnot(file.exists(mecp2_down_path))

cat("Reading MeCP2 peak files...\n")
up <- read.table(mecp2_up_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
down <- read.table(mecp2_down_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(up) <- c("chr", "start", "end")
colnames(down) <- c("chr", "start", "end")

up$Annotation <- "MeCP2_Up"
down$Annotation <- "MeCP2_Down"

combined <- rbind(up, down)

# Strip chr prefix for modality format
combined$Chromosome <- gsub("^chr", "", combined$chr)
combined$Name <- paste0(combined$chr, ":", combined$start, "-", combined$end)

out_df <- combined[, c("Chromosome", "start", "end", "Annotation", "Name")]
colnames(out_df) <- c("Chromosome", "Start", "End", "Annotation", "Name")

# Sort by chromosome and start
out_df <- out_df[order(out_df$Chromosome, out_df$Start), ]

cat(sprintf("  MeCP2 Up: %d peaks\n", nrow(up)))
cat(sprintf("  MeCP2 Down: %d peaks\n", nrow(down)))
cat(sprintf("  Total: %d peaks\n", nrow(out_df)))

# Write to regions directory
out_dir <- file.path(getwd(), "modality/mecp2_regions")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_path <- file.path(out_dir, "mecp2_peaks.annotation.bed.gz")
gz <- gzfile(out_path, "w")
write.table(out_df, gz, sep = "\t", row.names = FALSE, quote = FALSE)
close(gz)

cat(sprintf("\nWrote: %s\n", out_path))
cat("\nNext steps:\n")
cat("  1. scp modality/mecp2_regions/mecp2_peaks.annotation.bed.gz \\\n")
cat("       expanse:/expanse/lustre/projects/csd940/zalibhai/biomodal/modality/mecp2_regions/\n")
cat("  2. sbatch scripts/run_modality_mecp2.sb CHG\n")
cat("  3. sbatch scripts/run_modality_mecp2.sb CHH\n")
