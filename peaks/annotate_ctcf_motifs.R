# peaks/annotate_ctcf_motifs.R
# Integrate CTCF motif counts (from bedtools) with existing loop annotations
# Generates concordance matrix comparing motif-based vs ChIP-seq-based CTCF detection
# Created: 2026-01-13
# Revised: 2026-01-13 (simplified to use bedtools output instead of HOMER)

suppressPackageStartupMessages({
  library(tidyverse)
})

# ============================================================================
# Configuration
# ============================================================================

# Input files
CTCF_COUNTS_FILE <- "loop_annotation_extended/early/anchors_ctcf_motif_count.bed"
ANCHOR_MAPPING <- "loop_annotation_extended/early/other_anchors_mapping.tsv"
EXISTING_LOOPS <- "loop_annotation_extended/early/other_category_up_loops.tsv"

# Output files
OUTPUT_DIR <- "loop_annotation_extended/early"
VALIDATION_REPORT <- file.path(OUTPUT_DIR, "ctcf_motif_validation.txt")
UPDATED_LOOPS <- file.path(OUTPUT_DIR, "other_category_up_loops.tsv")
UPDATED_SUMMARY <- file.path(OUTPUT_DIR, "other_up_loops_summary.txt")

# ============================================================================
# Load Bedtools Output
# ============================================================================

cat("=== Integrating CTCF Motif Results ===\n\n")

cat("Loading bedtools intersection output...\n")

if (!file.exists(CTCF_COUNTS_FILE)) {
  stop(sprintf("Bedtools output not found: %s\nRun intersect_ctcf_motifs.sh on HPC first.",
               CTCF_COUNTS_FILE))
}

# Bedtools output format: chr start end anchor_id score strand ctcf_count
ctcf_data <- read.table(CTCF_COUNTS_FILE, header = FALSE, sep = "\t",
                        stringsAsFactors = FALSE)
colnames(ctcf_data) <- c("chr", "start", "end", "anchor_id", "score", "strand",
                         "ctcf_motif_count")

cat(sprintf("  Loaded %d anchors with CTCF motif counts\n", nrow(ctcf_data)))

# Derive has_CTCF_motif from count
ctcf_data$has_CTCF_motif <- ctcf_data$ctcf_motif_count > 0

# Summary
n_with_motif <- sum(ctcf_data$has_CTCF_motif)
n_total <- nrow(ctcf_data)
cat(sprintf("  Anchors with CTCF motif: %d / %d (%.1f%%)\n",
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

# Join motif data with mapping
mapping_with_motif <- mapping %>%
  left_join(
    ctcf_data %>% select(anchor_id, ctcf_motif_count, has_CTCF_motif),
    by = "anchor_id"
  )

cat(sprintf("  Mapped entries: %d\n", nrow(mapping_with_motif)))

# Check for unmapped anchors
n_unmapped <- sum(is.na(mapping_with_motif$has_CTCF_motif))
if (n_unmapped > 0) {
  cat(sprintf("  WARNING: %d anchors could not be matched\n", n_unmapped))
}

# ============================================================================
# Aggregate Motif Info per Loop
# ============================================================================

cat("\n=== Aggregating Motif Info per Loop ===\n")

# For each loop, determine if anchor1 and/or anchor2 has CTCF motif
loop_motif_info <- mapping_with_motif %>%
  group_by(loop_id) %>%
  summarize(
    anchor1_has_CTCF_motif = any(anchor == "anchor1" & has_CTCF_motif,
                                 na.rm = TRUE),
    anchor2_has_CTCF_motif = any(anchor == "anchor2" & has_CTCF_motif,
                                 na.rm = TRUE),
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
# - has_CTCF_motif (in silico, from pre-computed motifs)
# - overlaps_CTCF (late ChIP-seq proxy)

# Collect all "Other" anchor data
all_other_anchor1 <- updated_loops %>%
  filter(anchor1_type == "Other") %>%
  transmute(
    has_motif = anchor1_has_CTCF_motif,
    overlaps_chip = anchor1_overlaps_CTCF,
    motif_count = anchor1_CTCF_motif_count
  )

all_other_anchor2 <- updated_loops %>%
  filter(anchor2_type == "Other") %>%
  transmute(
    has_motif = anchor2_has_CTCF_motif,
    overlaps_chip = anchor2_overlaps_CTCF,
    motif_count = anchor2_CTCF_motif_count
  )

all_other_anchors <- bind_rows(all_other_anchor1, all_other_anchor2)

# Build concordance table
concordance <- all_other_anchors %>%
  mutate(
    motif_status = ifelse(has_motif, "Motif+", "Motif-"),
    chip_status = ifelse(overlaps_chip, "ChIP+", "ChIP-")
  ) %>%
  count(motif_status, chip_status) %>%
  pivot_wider(names_from = chip_status, values_from = n, values_fill = 0)

cat("\nConcordance Matrix (all 'Other' anchors):\n")
print(concordance)

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

# Extract concordance values safely
get_concordance_val <- function(df, motif, chip) {
  row <- df[df$motif_status == motif, ]
  if (nrow(row) == 0) return(0)
  if (chip %in% colnames(row)) return(row[[chip]])
  return(0)
}

motif_pos_chip_pos <- get_concordance_val(concordance, "Motif+", "ChIP+")
motif_pos_chip_neg <- get_concordance_val(concordance, "Motif+", "ChIP-")
motif_neg_chip_pos <- get_concordance_val(concordance, "Motif-", "ChIP+")
motif_neg_chip_neg <- get_concordance_val(concordance, "Motif-", "ChIP-")

report_text <- sprintf("
================================================================================
CTCF MOTIF VALIDATION FOR EARLY 'OTHER' ANCHORS
Generated: %s
================================================================================

OVERVIEW
--------
This analysis validates early timepoint 'Other' category anchors by checking
for CTCF DNA binding motifs. Since the early timepoint lacks CTCF ChIP-seq
data, many potential CTCF sites are classified as 'Other'.

We compare two approaches:
1. CTCF motif presence (from HOMER pre-computed genome-wide motif scan)
2. Late timepoint CTCF ChIP-seq overlap (proxy)

METHOD
------
- CTCF motifs extracted from: homer.KnownMotifs.mm10.191020.bed.gz
- Used canonical CTCF(Zf) motif only (excluded Satellite variants)
- bedtools intersect to count motifs per anchor

RESULTS SUMMARY
---------------
Total 'Other' anchors analyzed: %d

                          Count       Percent
CTCF motif+ (in silico):     %d         %.1f%%
Late CTCF ChIP+:             %d         %.1f%%
Both motif+ AND ChIP+:       %d         %.1f%%
Either motif+ OR ChIP+:      %d         %.1f%%
Neither (true 'Other'):      %d         %.1f%%

CONCORDANCE MATRIX
------------------
                    Late CTCF ChIP+    Late CTCF ChIP-
CTCF Motif+              %d                  %d
CTCF Motif-              %d                  %d

INTERPRETATION
--------------
- %.1f%% of 'Other' anchors have CTCF evidence (motif or ChIP)
- %.1f%% are truly unmarked structural elements (neither motif nor ChIP)

Categories:
- Motif+ ChIP+ : Strong CTCF site (motif present, ChIP binding confirmed)
- Motif+ ChIP- : Motif present but no late ChIP binding
                 (may be timepoint-specific or weak binding)
- Motif- ChIP+ : ChIP binding without canonical motif
                 (non-canonical binding or indirect recruitment)
- Motif- ChIP- : True 'Other' - structural loop without CTCF

BIOLOGICAL NOTES
----------------
- High Motif+ ChIP- may indicate early-specific CTCF binding not present in late
- Motif- ChIP+ anchors may use non-canonical CTCF binding or cohesin recruitment
- True 'Other' anchors may represent cohesin-only loops or novel structural elements

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
  motif_pos_chip_pos, motif_pos_chip_neg,
  motif_neg_chip_pos, motif_neg_chip_neg,
  100 * n_either_pos / n_total_other,
  100 * n_neither / n_total_other,
  VALIDATION_REPORT,
  UPDATED_LOOPS
)

writeLines(report_text, VALIDATION_REPORT)
cat(sprintf("  Saved: %s\n", VALIDATION_REPORT))

# ============================================================================
# Append to Existing Summary
# ============================================================================

cat("\n=== Appending to Summary File ===\n")

if (file.exists(UPDATED_SUMMARY)) {
  existing_summary <- readLines(UPDATED_SUMMARY)

  motif_section <- sprintf("

CTCF MOTIF VALIDATION (Added %s)
--------------------------------
Method: HOMER pre-computed motifs (homer.KnownMotifs.mm10.191020.bed.gz)
        Canonical CTCF(Zf) motif only, bedtools intersection

Total 'Other' anchors: %d
  - CTCF motif+ (in silico): %d (%.1f%%)
  - Late CTCF ChIP+ (proxy): %d (%.1f%%)
  - Either motif+ or ChIP+: %d (%.1f%%)
  - Neither (true 'Other'): %d (%.1f%%)

See ctcf_motif_validation.txt for full concordance matrix.
",
    format(Sys.time(), "%%Y-%%m-%%d"),
    n_total_other,
    n_motif_pos, 100 * n_motif_pos / n_total_other,
    n_chip_pos, 100 * n_chip_pos / n_total_other,
    n_either_pos, 100 * n_either_pos / n_total_other,
    n_neither, 100 * n_neither / n_total_other
  )

  writeLines(c(existing_summary, motif_section), UPDATED_SUMMARY)
  cat(sprintf("  Appended motif section to: %s\n", UPDATED_SUMMARY))
} else {
  cat(sprintf("  Summary file not found: %s\n", UPDATED_SUMMARY))
}

cat("\n=== DONE ===\n")
