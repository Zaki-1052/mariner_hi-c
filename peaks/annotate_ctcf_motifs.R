# peaks/annotate_ctcf_motifs.R
# Integrate HOMER CTCF motif scanning results with existing loop annotations
# Generates concordance matrix comparing motif-based vs ChIP-seq-based CTCF detection
# Created: 2026-01-13

suppressPackageStartupMessages({
  library(tidyverse)
})

# ============================================================================
# Configuration
# ============================================================================

# Input files
HOMER_OUTPUT <- "loop_annotation_extended/early/homer_ctcf_annotate.txt"
ANCHOR_MAPPING <- "loop_annotation_extended/early/other_anchors_mapping.tsv"
EXISTING_LOOPS <- "loop_annotation_extended/early/other_category_up_loops.tsv"
FULL_LOOPS <- "loop_annotation_extended/early/other_up_loops_with_genes.tsv"

# Output files
OUTPUT_DIR <- "loop_annotation_extended/early"
VALIDATION_REPORT <- file.path(OUTPUT_DIR, "ctcf_motif_validation.txt")
UPDATED_LOOPS <- file.path(OUTPUT_DIR, "other_category_up_loops.tsv")
UPDATED_SUMMARY <- file.path(OUTPUT_DIR, "other_up_loops_summary.txt")

# ============================================================================
# Load HOMER Results
# ============================================================================

cat("=== Integrating CTCF Motif Results ===\n\n")

cat("Loading HOMER annotatePeaks output...\n")

# HOMER annotatePeaks output is tab-delimited with header
# Key columns: PeakID (col 1), and motif counts (varies by -m option)
# When using -nmotifs with -m ctcf, the CTCF motif count is in a column

if (!file.exists(HOMER_OUTPUT)) {
  stop(sprintf("HOMER output not found: %s\nRun scan_ctcf_motifs.sb first.",
               HOMER_OUTPUT))
}

homer_raw <- read.table(HOMER_OUTPUT, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE, comment.char = "",
                        quote = "", fill = TRUE)

cat(sprintf("  HOMER output rows: %d\n", nrow(homer_raw)))
cat(sprintf("  Columns: %s\n", paste(head(colnames(homer_raw), 10), collapse = ", ")))

# ============================================================================
# Parse CTCF Motif Counts
# ============================================================================

cat("\n=== Parsing CTCF Motif Counts ===\n")
# HOMER annotatePeaks column names can vary
# Look for column containing "ctcf" (case insensitive)
ctcf_cols <- grep("ctcf", colnames(homer_raw), ignore.case = TRUE, value = TRUE)
cat(sprintf("  CTCF-related columns found: %s\n", paste(ctcf_cols, collapse = ", ")))

# The first column is typically PeakID/Chr:Start-End format
# Extract anchor_id from first column
homer_data <- homer_raw %>%
  mutate(
    # First column is peak identifier - standardize naming
    anchor_id = homer_raw[[1]]
  )

# Find CTCF motif count column
# HOMER uses format like "CTCF(Zf)/..." for the motif name
if (length(ctcf_cols) > 0) {
  # Use first CTCF column found (should be motif count)
  ctcf_col <- ctcf_cols[1]
  homer_data$ctcf_motif_count <- homer_raw[[ctcf_col]]
  cat(sprintf("  Using column '%s' for CTCF motif counts\n", ctcf_col))
} else {
  # If no explicit CTCF column, check for general motif density columns
  cat("  WARNING: No explicit CTCF column found. Checking for motif density...\n")

  # Look for columns that might contain motif information
  motif_cols <- grep("motif|Motif", colnames(homer_raw), value = TRUE)
  if (length(motif_cols) > 0) {
    cat(sprintf("  Found motif columns: %s\n", paste(motif_cols, collapse = ", ")))
    # Use first motif column as fallback
    homer_data$ctcf_motif_count <- homer_raw[[motif_cols[1]]]
  } else {
    # Set to NA if no motif information found
    cat("  WARNING: No motif columns found. Setting counts to NA.\n")
    homer_data$ctcf_motif_count <- NA
  }
}

# Convert to numeric (HOMER may output as character)
homer_data$ctcf_motif_count <- as.numeric(homer_data$ctcf_motif_count)

# Determine has_CTCF_motif based on count > 0
homer_data$has_CTCF_motif <- !is.na(homer_data$ctcf_motif_count) &
                              homer_data$ctcf_motif_count > 0

# Summary of motif detection
n_with_motif <- sum(homer_data$has_CTCF_motif, na.rm = TRUE)
n_total <- nrow(homer_data)
cat(sprintf("\n  Anchors with CTCF motif: %d / %d (%.1f%%)\n",
            n_with_motif, n_total, 100 * n_with_motif / n_total))

# ============================================================================
# Load Anchor Mapping and Existing Loop Data
# ============================================================================

cat("\n=== Loading Mapping and Loop Data ===\n")

# Load anchor-to-loop mapping
mapping <- read.table(ANCHOR_MAPPING, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
cat(sprintf("  Anchor mapping entries: %d\n", nrow(mapping)))

# Load existing loop annotations (has CTCF ChIP-seq overlap)
existing_loops <- read.table(EXISTING_LOOPS, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE)
cat(sprintf("  Existing loop annotations: %d\n", nrow(existing_loops)))

# ============================================================================
# Join Motif Results with Mapping
# ============================================================================

cat("\n=== Joining Motif Results with Loop Data ===\n")

# Clean anchor_id in homer_data if needed (remove quotes, standardize format)
homer_data$anchor_id <- gsub("\"", "", homer_data$anchor_id)

# Join motif results with mapping
mapping_with_motif <- mapping %>%
  left_join(
    homer_data %>% select(anchor_id, ctcf_motif_count, has_CTCF_motif),
    by = "anchor_id"
  )

cat(sprintf("  Mapped entries: %d\n", nrow(mapping_with_motif)))

# Check for unmapped anchors
n_unmapped <- sum(is.na(mapping_with_motif$has_CTCF_motif))
if (n_unmapped > 0) {
  cat(sprintf("  WARNING: %d anchors could not be matched to HOMER output\n", n_unmapped))
}

# ============================================================================
# Aggregate Motif Info per Loop
# ============================================================================

cat("\n=== Aggregating Motif Info per Loop ===\n")

# For each loop, determine if anchor1 and/or anchor2 has CTCF motif
loop_motif_info <- mapping_with_motif %>%
  group_by(loop_id) %>%
  summarize(
    anchor1_has_CTCF_motif = any(anchor == "anchor1" & has_CTCF_motif, na.rm = TRUE),
    anchor2_has_CTCF_motif = any(anchor == "anchor2" & has_CTCF_motif, na.rm = TRUE),
    anchor1_CTCF_motif_count = max(
      ifelse(anchor == "anchor1", ctcf_motif_count, 0), na.rm = TRUE
    ),
    anchor2_CTCF_motif_count = max(
      ifelse(anchor == "anchor2", ctcf_motif_count, 0), na.rm = TRUE
    ),
    .groups = "drop"
  )

# Handle infinite values from max() on empty sets
loop_motif_info <- loop_motif_info %>%
  mutate(
    anchor1_CTCF_motif_count = ifelse(
      is.infinite(anchor1_CTCF_motif_count), 0, anchor1_CTCF_motif_count
    ),
    anchor2_CTCF_motif_count = ifelse(
      is.infinite(anchor2_CTCF_motif_count), 0, anchor2_CTCF_motif_count
    )
  )

cat(sprintf("  Loops with motif info: %d\n", nrow(loop_motif_info)))

# ============================================================================
# Merge with Existing Loop Data
# ============================================================================

cat("\n=== Merging with Existing Loop Annotations ===\n")

# Join motif info to existing loops
updated_loops <- existing_loops %>%
  left_join(loop_motif_info, by = "loop_id")

cat(sprintf("  Updated loops: %d\n", nrow(updated_loops)))

# ============================================================================
# Build Concordance Matrix
# ============================================================================

cat("\n=== Building Concordance Matrix ===\n")

# For "Other" anchors, compare:
# - has_CTCF_motif (in silico prediction)
# - overlaps_CTCF (late ChIP-seq proxy)

# Separate analysis for anchor1 and anchor2

# Anchor1 concordance (only for "Other" anchor1)
anchor1_concordance <- updated_loops %>%
  filter(anchor1_type == "Other") %>%
  mutate(
    motif_status = ifelse(anchor1_has_CTCF_motif, "Motif+", "Motif-"),
    chip_status = ifelse(anchor1_overlaps_CTCF, "ChIP+", "ChIP-")
  ) %>%
  count(motif_status, chip_status) %>%
  pivot_wider(names_from = chip_status, values_from = n, values_fill = 0)

cat("\nAnchor1 'Other' Concordance:\n")
print(anchor1_concordance)

# Anchor2 concordance (only for "Other" anchor2)
anchor2_concordance <- updated_loops %>%
  filter(anchor2_type == "Other") %>%
  mutate(
    motif_status = ifelse(anchor2_has_CTCF_motif, "Motif+", "Motif-"),
    chip_status = ifelse(anchor2_overlaps_CTCF, "ChIP+", "ChIP-")
  ) %>%
  count(motif_status, chip_status) %>%
  pivot_wider(names_from = chip_status, values_from = n, values_fill = 0)

cat("\nAnchor2 'Other' Concordance:\n")
print(anchor2_concordance)

# Combined concordance (all "Other" anchors)
all_other_anchor1 <- updated_loops %>%
  filter(anchor1_type == "Other") %>%
  transmute(
    has_motif = anchor1_has_CTCF_motif,
    overlaps_chip = anchor1_overlaps_CTCF
  )

all_other_anchor2 <- updated_loops %>%
  filter(anchor2_type == "Other") %>%
  transmute(
    has_motif = anchor2_has_CTCF_motif,
    overlaps_chip = anchor2_overlaps_CTCF
  )

all_other_anchors <- bind_rows(all_other_anchor1, all_other_anchor2)

combined_concordance <- all_other_anchors %>%
  mutate(
    motif_status = ifelse(has_motif, "Motif+", "Motif-"),
    chip_status = ifelse(overlaps_chip, "ChIP+", "ChIP-")
  ) %>%
  count(motif_status, chip_status) %>%
  pivot_wider(names_from = chip_status, values_from = n, values_fill = 0)

cat("\nCombined 'Other' Anchors Concordance:\n")
print(combined_concordance)

# ============================================================================
# Calculate Summary Statistics
# ============================================================================

cat("\n=== Summary Statistics ===\n")

n_total_other <- nrow(all_other_anchors)
n_motif_pos <- sum(all_other_anchors$has_motif, na.rm = TRUE)
n_chip_pos <- sum(all_other_anchors$overlaps_chip, na.rm = TRUE)
n_both_pos <- sum(all_other_anchors$has_motif & all_other_anchors$overlaps_chip,
                  na.rm = TRUE)
n_either_pos <- sum(all_other_anchors$has_motif | all_other_anchors$overlaps_chip,
                    na.rm = TRUE)
n_neither <- sum(!all_other_anchors$has_motif & !all_other_anchors$overlaps_chip,
                 na.rm = TRUE)

cat(sprintf("Total 'Other' anchors: %d\n", n_total_other))
cat(sprintf("  CTCF motif+ (in silico): %d (%.1f%%)\n",
            n_motif_pos, 100 * n_motif_pos / n_total_other))
cat(sprintf("  Late CTCF ChIP+: %d (%.1f%%)\n",
            n_chip_pos, 100 * n_chip_pos / n_total_other))
cat(sprintf("  Both motif+ and ChIP+: %d (%.1f%%)\n",
            n_both_pos, 100 * n_both_pos / n_total_other))
cat(sprintf("  Either motif+ or ChIP+: %d (%.1f%%)\n",
            n_either_pos, 100 * n_either_pos / n_total_other))
cat(sprintf("  Neither (true 'Other'): %d (%.1f%%)\n",
            n_neither, 100 * n_neither / n_total_other))

# ============================================================================
# Save Updated Loop File
# ============================================================================

cat("\n=== Saving Updated Files ===\n")

write.table(updated_loops, UPDATED_LOOPS, sep = "\t", quote = FALSE,
            row.names = FALSE)
cat(sprintf("  Updated: %s\n", UPDATED_LOOPS))

# ============================================================================
# Generate Validation Report
# ============================================================================

report_text <- sprintf("
================================================================================
CTCF MOTIF VALIDATION FOR EARLY 'OTHER' ANCHORS
Generated: %s
================================================================================

OVERVIEW
--------
This analysis validates early timepoint 'Other' category anchors by scanning
for CTCF DNA binding motifs using HOMER. Since the early timepoint lacks CTCF
ChIP-seq data, many potential CTCF sites are classified as 'Other'.

We compare two approaches:
1. CTCF motif scanning (in silico prediction)
2. Late timepoint CTCF ChIP-seq overlap (proxy)

RESULTS SUMMARY
---------------
Total 'Other' anchors analyzed: %d

                        Count       Percent
CTCF motif+ (in silico):   %d         %.1f%%
Late CTCF ChIP+:           %d         %.1f%%
Both motif+ AND ChIP+:     %d         %.1f%%
Either motif+ OR ChIP+:    %d         %.1f%%
Neither (true 'Other'):    %d         %.1f%%

CONCORDANCE MATRIX
------------------
                    Late CTCF ChIP+    Late CTCF ChIP-
CTCF Motif+         %d (%.1f%%)          %d (%.1f%%)
CTCF Motif-         %d (%.1f%%)          %d (%.1f%%)

INTERPRETATION
--------------
- %.1f%% of 'Other' anchors have CTCF evidence (motif or ChIP)
- %.1f%% are concordant (both motif+ and ChIP+ or both negative)
- %.1f%% remain truly 'Other' (no CTCF evidence)

These 'Other' anchors may represent:
1. Cohesin-only binding sites (CTCF-independent loop anchors)
2. Novel structural elements
3. Weak or indirect CTCF binding not detectable by motif scanning

TECHNICAL NOTES
---------------
- HOMER annotatePeaks.pl with -m ctcf -nmotifs
- Motif threshold: default HOMER settings
- ChIP-seq data: late timepoint (P21+) as proxy for early (P13)
- Some timepoint-specific CTCF binding may not be captured by either method

FILES GENERATED
---------------
- %s (this report)
- %s (updated with motif columns)

================================================================================
",
  Sys.time(),
  n_total_other,
  n_motif_pos, 100 * n_motif_pos / n_total_other,
  n_chip_pos, 100 * n_chip_pos / n_total_other,
  n_both_pos, 100 * n_both_pos / n_total_other,
  n_either_pos, 100 * n_either_pos / n_total_other,
  n_neither, 100 * n_neither / n_total_other,
  # Concordance matrix values - need to extract from combined_concordance
  ifelse("ChIP+" %in% colnames(combined_concordance),
         combined_concordance$`ChIP+`[combined_concordance$motif_status == "Motif+"], 0),
  ifelse("ChIP+" %in% colnames(combined_concordance),
         100 * combined_concordance$`ChIP+`[combined_concordance$motif_status == "Motif+"] / n_total_other, 0),
  ifelse("ChIP-" %in% colnames(combined_concordance),
         combined_concordance$`ChIP-`[combined_concordance$motif_status == "Motif+"], 0),
  ifelse("ChIP-" %in% colnames(combined_concordance),
         100 * combined_concordance$`ChIP-`[combined_concordance$motif_status == "Motif+"] / n_total_other, 0),
  ifelse("ChIP+" %in% colnames(combined_concordance),
         combined_concordance$`ChIP+`[combined_concordance$motif_status == "Motif-"], 0),
  ifelse("ChIP+" %in% colnames(combined_concordance),
         100 * combined_concordance$`ChIP+`[combined_concordance$motif_status == "Motif-"] / n_total_other, 0),
  ifelse("ChIP-" %in% colnames(combined_concordance),
         combined_concordance$`ChIP-`[combined_concordance$motif_status == "Motif-"], 0),
  ifelse("ChIP-" %in% colnames(combined_concordance),
         100 * combined_concordance$`ChIP-`[combined_concordance$motif_status == "Motif-"] / n_total_other, 0),
  # Interpretation percentages
  100 * n_either_pos / n_total_other,
  100 * (n_both_pos + n_neither) / n_total_other,
  100 * n_neither / n_total_other,
  # File names
  VALIDATION_REPORT,
  UPDATED_LOOPS
)

writeLines(report_text, VALIDATION_REPORT)
cat(sprintf("  Saved: %s\n", VALIDATION_REPORT))

# ============================================================================
# Append to Existing Summary
# ============================================================================

cat("\n=== Appending to Summary File ===\n")

# Read existing summary
if (file.exists(UPDATED_SUMMARY)) {
  existing_summary <- readLines(UPDATED_SUMMARY)

  # Add motif validation section
  motif_section <- sprintf("

CTCF MOTIF VALIDATION (Added %s)
--------------------------------
Total 'Other' anchors: %d
  - CTCF motif+ (in silico): %d (%.1f%%)
  - Late CTCF ChIP+ (proxy): %d (%.1f%%)
  - Either motif+ or ChIP+: %d (%.1f%%)
  - Neither (true 'Other'): %d (%.1f%%)

See ctcf_motif_validation.txt for full concordance matrix.
",
    format(Sys.time(), "%Y-%m-%d"),
    n_total_other,
    n_motif_pos, 100 * n_motif_pos / n_total_other,
    n_chip_pos, 100 * n_chip_pos / n_total_other,
    n_either_pos, 100 * n_either_pos / n_total_other,
    n_neither, 100 * n_neither / n_total_other
  )

  # Append to summary
  writeLines(c(existing_summary, motif_section), UPDATED_SUMMARY)
  cat(sprintf("  Appended motif section to: %s\n", UPDATED_SUMMARY))
} else {
  cat(sprintf("  Summary file not found: %s\n", UPDATED_SUMMARY))
}

cat("\n=== DONE ===\n")
