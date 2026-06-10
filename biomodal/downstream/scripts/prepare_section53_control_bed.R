# biomodal/downstream/scripts/prepare_section53_control_bed.R
# Generate matched control regions for Section 53 genome-wide baseline comparison.
# Shuffles 8,886 regions (same count and size as MeCP2 peaks) to random
# genomic positions, excluding blacklist and MeCP2 peaks themselves.
#
# Output format matches modality Regions_Directory BED:
#   Chromosome\tStart\tEnd\tAnnotation\tName (Ensembl chromosomes, gzipped)
#
# Run from downstream/ directory:
#   Rscript scripts/prepare_section53_control_bed.R
#
# Then: git add/commit/push → HPC git pull → sbatch run_modality_control.sb CHG/CHH

suppressPackageStartupMessages({
  library(GenomicRanges)
})

set.seed(42)

REPO_ROOT <- normalizePath(file.path(getwd(), "../.."))

# ---------------------------------------------------------------------------
# Input paths
# ---------------------------------------------------------------------------
MECP2_BED <- file.path(getwd(), "modality/mecp2/mecp2_peaks.annotation.bed.gz")
CHROM_SIZES <- file.path(REPO_ROOT, "peaks/mm10.chrom.sizes")
BLACKLIST <- file.path(REPO_ROOT, "tads/mm10-blacklist.v2.bed")

stopifnot(file.exists(MECP2_BED))
stopifnot(file.exists(CHROM_SIZES))
stopifnot(file.exists(BLACKLIST))

# Output
OUT_DIR <- file.path(getwd(), "modality/mecp2")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PATH <- file.path(OUT_DIR, "control_peaks.annotation.bed.gz")

# ---------------------------------------------------------------------------
# Load MeCP2 peaks
# ---------------------------------------------------------------------------
cat("Loading MeCP2 peaks...\n")
mecp2 <- read.table(gzfile(MECP2_BED), header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)
cat(sprintf("  %d MeCP2 peaks loaded\n", nrow(mecp2)))

mecp2_gr <- GRanges(
  seqnames = mecp2$Chromosome,
  ranges = IRanges(start = mecp2$Start, end = mecp2$End)
)

# ---------------------------------------------------------------------------
# Load chromosome sizes (UCSC → Ensembl)
# ---------------------------------------------------------------------------
cat("Loading chromosome sizes...\n")
chroms <- read.table(CHROM_SIZES, header = FALSE, sep = "\t",
                     stringsAsFactors = FALSE, col.names = c("chr", "size"))
chroms$chr_ensembl <- gsub("^chr", "", chroms$chr)
chroms <- chroms[!grepl("_", chroms$chr), ]
chroms <- chroms[chroms$chr != "chrM", ]
cat(sprintf("  %d chromosomes\n", nrow(chroms)))

# ---------------------------------------------------------------------------
# Load blacklist (UCSC → Ensembl)
# ---------------------------------------------------------------------------
cat("Loading blacklist...\n")
bl <- read.table(BLACKLIST, header = FALSE, sep = "\t",
                 stringsAsFactors = FALSE)
bl_gr <- GRanges(
  seqnames = gsub("^chr", "", bl$V1),
  ranges = IRanges(start = bl$V2, end = bl$V3)
)
cat(sprintf("  %d blacklist regions\n", length(bl_gr)))

exclude_gr <- c(mecp2_gr, bl_gr)

# ---------------------------------------------------------------------------
# Generate shuffled control regions
# ---------------------------------------------------------------------------
cat("Generating shuffled control regions...\n")

peak_widths <- width(mecp2_gr)
peak_chroms <- as.character(seqnames(mecp2_gr))

control_list <- vector("list", length(mecp2_gr))
max_attempts <- 1000

for (i in seq_along(mecp2_gr)) {
  chr <- peak_chroms[i]
  w <- peak_widths[i]
  chr_size <- chroms$size[chroms$chr_ensembl == chr]
  stopifnot(length(chr_size) == 1)

  placed <- FALSE
  for (attempt in seq_len(max_attempts)) {
    s <- sample.int(chr_size - w, 1)
    candidate <- GRanges(seqnames = chr, ranges = IRanges(start = s, end = s + w - 1))
    if (sum(countOverlaps(candidate, exclude_gr)) == 0) {
      control_list[[i]] <- candidate
      exclude_gr <- c(exclude_gr, candidate)
      placed <- TRUE
      break
    }
  }
  stopifnot(placed)

  if (i %% 1000 == 0) cat(sprintf("  %d / %d placed\n", i, length(mecp2_gr)))
}

control_gr <- do.call(c, control_list)
cat(sprintf("  %d control regions generated\n", length(control_gr)))

# ---------------------------------------------------------------------------
# Format output BED
# ---------------------------------------------------------------------------
out_df <- data.frame(
  Chromosome = as.character(seqnames(control_gr)),
  Start = start(control_gr),
  End = end(control_gr),
  Annotation = "Control",
  Name = paste0("chr", as.character(seqnames(control_gr)), ":",
                start(control_gr), "-", end(control_gr)),
  stringsAsFactors = FALSE
)
out_df <- out_df[order(out_df$Chromosome, out_df$Start), ]

gz <- gzfile(OUT_PATH, "w")
write.table(out_df, gz, sep = "\t", row.names = FALSE, quote = FALSE)
close(gz)

cat(sprintf("\nWrote: %s (%d regions)\n", OUT_PATH, nrow(out_df)))
cat("\nNext steps:\n")
cat("  1. git add modality/mecp2/control_peaks.annotation.bed.gz\n")
cat("  2. git commit -m 'add shuffled control regions for section 53'\n")
cat("  3. git push → HPC git pull\n")
cat("  4. sbatch run_modality_control.sb CHG\n")
cat("  5. sbatch run_modality_control.sb CHH\n")
