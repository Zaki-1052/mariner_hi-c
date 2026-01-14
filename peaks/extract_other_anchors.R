# peaks/extract_other_anchors.R
# Extract unique "Other" category anchors from early timepoint to BED format
# For CTCF motif scanning with HOMER
# Created: 2026-01-13

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(tidyverse)
})

# ============================================================================
# Configuration
# ============================================================================

INPUT_FILE <- "loop_annotation_extended/early/extended_characterized_loops.tsv"
OUTPUT_DIR <- "loop_annotation_extended/early"
OUTPUT_BED <- file.path(OUTPUT_DIR, "other_anchors.bed")

# ============================================================================
# Load Data
# ============================================================================

cat("=== Extracting 'Other' Category Anchors for CTCF Motif Scanning ===\n\n")

cat("Loading early timepoint characterized loops...\n")
loops <- read.table(INPUT_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  Total loops: %d\n", nrow(loops)))

# ============================================================================
# Filter for Up Loops with "Other" Anchors
# ============================================================================

cat("\n=== Filtering for 'Other' Anchors in Up Loops ===\n")

up_loops <- loops %>% filter(direction == "up_in_mutant")
cat(sprintf("  Up loops: %d\n", nrow(up_loops)))

# Extract anchor1 coordinates where type is "Other"
other_anchor1 <- up_loops %>%
  filter(anchor1_type_extended == "Other") %>%
  select(chr = chr1, start = start1, end = end1, loop_id) %>%
  mutate(anchor = "anchor1")

cat(sprintf("  Anchor1 'Other': %d\n", nrow(other_anchor1)))

# Extract anchor2 coordinates where type is "Other"
other_anchor2 <- up_loops %>%
  filter(anchor2_type_extended == "Other") %>%
  select(chr = chr2, start = start2, end = end2, loop_id) %>%
  mutate(anchor = "anchor2")

cat(sprintf("  Anchor2 'Other': %d\n", nrow(other_anchor2)))

# Combine all "Other" anchors
all_other_anchors <- bind_rows(other_anchor1, other_anchor2)
cat(sprintf("  Total 'Other' anchor entries: %d\n", nrow(all_other_anchors)))

# ============================================================================
# Deduplicate Anchors
# ============================================================================

cat("\n=== Deduplicating Anchors ===\n")

# Create unique anchor ID for deduplication
all_other_anchors <- all_other_anchors %>%
  mutate(anchor_id = sprintf("%s:%d-%d", chr, start, end))

# Get unique anchors (keep first occurrence)
unique_anchors <- all_other_anchors %>%
  distinct(anchor_id, .keep_all = TRUE)

cat(sprintf("  Unique 'Other' anchors: %d\n", nrow(unique_anchors)))

# Show which loops share anchors (if any)
duplicated_anchors <- all_other_anchors %>%
  group_by(anchor_id) %>%
  filter(n() > 1) %>%
  summarize(n_loops = n(), loops = paste(unique(loop_id), collapse = ", "))

if (nrow(duplicated_anchors) > 0) {
  cat(sprintf("  Anchors shared by multiple loops: %d\n", nrow(duplicated_anchors)))
}

# ============================================================================
# Create BED File
# ============================================================================

cat("\n=== Creating BED File ===\n")

# Format for BED: chr, start, end, name, score, strand
# Name = anchor_id for tracking back to loops
bed_output <- unique_anchors %>%
  select(chr, start, end, name = anchor_id) %>%
  mutate(
    score = 0,
    strand = "."
  ) %>%
  arrange(chr, start)

# Write BED file (no header for standard BED format)
write.table(bed_output, OUTPUT_BED, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

cat(sprintf("  Saved: %s\n", OUTPUT_BED))
cat(sprintf("  Total regions: %d\n", nrow(bed_output)))

# ============================================================================
# Summary Statistics
# ============================================================================

cat("\n=== Summary Statistics ===\n")

# Region size distribution
region_sizes <- bed_output$end - bed_output$start
cat(sprintf("  Region size range: %d - %d bp\n", min(region_sizes), max(region_sizes)))
cat(sprintf("  Median region size: %d bp\n", median(region_sizes)))

# Chromosome distribution
chr_counts <- table(bed_output$chr)
cat("\n  Chromosome distribution:\n")
for (chr in names(sort(chr_counts, decreasing = TRUE))[1:5]) {
  cat(sprintf("    %s: %d anchors\n", chr, chr_counts[chr]))
}

# ============================================================================
# Create Mapping File (anchor_id to loop_id)
# ============================================================================

cat("\n=== Creating Anchor-to-Loop Mapping ===\n")

mapping_file <- file.path(OUTPUT_DIR, "other_anchors_mapping.tsv")
mapping <- all_other_anchors %>%
  select(anchor_id, loop_id, anchor)

write.table(mapping, mapping_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: %s\n", mapping_file))

cat("\n=== DONE ===\n")
cat("\nNext step: Run HOMER on HPC:\n")
cat("  sbatch scripts/scan_ctcf_motifs.sb\n")
