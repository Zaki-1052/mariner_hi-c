# peaks/extract_other_up_loops.R
# Extract significant "Other" category up loops from early timepoint for PI review
# Created: 2026-01-13

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(tidyverse)
})

# ============================================================================
# Configuration
# ============================================================================

INPUT_FILE <- "loop_annotation_extended/early/extended_characterized_loops.tsv"
CTCF_FILE <- "CTCF.bed"
OUTPUT_DIR <- "loop_annotation_extended/early"

# ============================================================================
# Helper Functions
# ============================================================================

get_nearest_gene <- function(chr, start, end, genes_gr, gene_symbols) {
  # Create GRanges for the anchor
  anchor_gr <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))

  # Find nearest gene
  nearest_idx <- nearest(anchor_gr, genes_gr)

  if (is.na(nearest_idx)) {
    return(list(gene = NA, distance = NA))
  }

  nearest_gene <- genes_gr[nearest_idx]
  distance <- distance(anchor_gr, nearest_gene)
  gene_id <- names(nearest_gene)

  # Convert to symbol
  symbol <- tryCatch({
    gene_symbols[gene_id]
  }, error = function(e) gene_id)

  if (is.na(symbol)) symbol <- gene_id

  return(list(gene = symbol, distance = distance))
}

# ============================================================================
# Load Data
# ============================================================================

cat("=== Extracting 'Other' Category Up Loops from Early Timepoint ===\n\n")

# Load early timepoint loops
cat("Loading early timepoint characterized loops...\n")
loops <- read.table(INPUT_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  Total loops: %d\n", nrow(loops)))

# Load CTCF peaks for overlap check (narrowPeak format)
cat("Loading CTCF peaks for overlap validation...\n")
ctcf_df <- read.table(CTCF_FILE, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
ctcf_peaks <- GRanges(
  seqnames = ctcf_df$V1,
  ranges = IRanges(start = ctcf_df$V2, end = ctcf_df$V3)
)
cat(sprintf("  CTCF peaks: %d\n", length(ctcf_peaks)))

# Load gene annotations
cat("Loading gene annotations...\n")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- genes(txdb)
gene_symbols <- mapIds(org.Mm.eg.db, keys = names(genes), column = "SYMBOL",
                       keytype = "ENTREZID", multiVals = "first")
cat(sprintf("  Gene annotations loaded: %d genes\n", length(genes)))

# ============================================================================
# Filter for "Other" Up Loops
# ============================================================================

cat("\n=== Filtering for 'Other' Category Up Loops ===\n")

# Filter for up_in_mutant
up_loops <- loops %>% filter(direction == "up_in_mutant")
cat(sprintf("  Up loops (up_in_mutant): %d\n", nrow(up_loops)))

# Filter for loops with at least one "Other" anchor
other_up_loops <- up_loops %>%
  filter(anchor1_type_extended == "Other" | anchor2_type_extended == "Other")
cat(sprintf("  Other up loops (at least one anchor): %d\n", nrow(other_up_loops)))

# Add breakdown category
other_up_loops <- other_up_loops %>%
  mutate(
    other_category = case_when(
      anchor1_type_extended == "Other" & anchor2_type_extended == "Other" ~ "both_other",
      TRUE ~ "single_other"
    )
  )

# Count breakdown
both_other <- sum(other_up_loops$other_category == "both_other")
single_other <- sum(other_up_loops$other_category == "single_other")
cat(sprintf("    - Both anchors Other: %d\n", both_other))
cat(sprintf("    - Single anchor Other: %d\n", single_other))

# ============================================================================
# CTCF Overlap Analysis
# ============================================================================

cat("\n=== CTCF Overlap Analysis ===\n")
cat("Checking if 'Other' anchors overlap with CTCF peaks...\n")

# Function to check CTCF overlap
check_ctcf_overlap <- function(chr, start, end, ctcf_gr) {
  anchor_gr <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))
  overlaps <- countOverlaps(anchor_gr, ctcf_gr)
  return(overlaps > 0)
}

# Check CTCF overlap for each anchor
other_up_loops <- other_up_loops %>%
  rowwise() %>%
  mutate(
    anchor1_overlaps_CTCF = check_ctcf_overlap(chr1, start1, end1, ctcf_peaks),
    anchor2_overlaps_CTCF = check_ctcf_overlap(chr2, start2, end2, ctcf_peaks)
  ) %>%
  ungroup()

# Summarize CTCF overlaps for "Other" anchors only
other_anchor1 <- other_up_loops %>% filter(anchor1_type_extended == "Other")
other_anchor2 <- other_up_loops %>% filter(anchor2_type_extended == "Other")

ctcf_in_other_anchor1 <- sum(other_anchor1$anchor1_overlaps_CTCF)
ctcf_in_other_anchor2 <- sum(other_anchor2$anchor2_overlaps_CTCF)

cat(sprintf("  Anchor1 'Other' with CTCF overlap: %d / %d (%.1f%%)\n",
            ctcf_in_other_anchor1, nrow(other_anchor1),
            100 * ctcf_in_other_anchor1 / max(nrow(other_anchor1), 1)))
cat(sprintf("  Anchor2 'Other' with CTCF overlap: %d / %d (%.1f%%)\n",
            ctcf_in_other_anchor2, nrow(other_anchor2),
            100 * ctcf_in_other_anchor2 / max(nrow(other_anchor2), 1)))

# ============================================================================
# Gene Annotation
# ============================================================================

cat("\n=== Adding Gene Annotations ===\n")

# Get nearest genes for each anchor
other_up_loops <- other_up_loops %>%
  rowwise() %>%
  mutate(
    anchor1_gene_info = list(get_nearest_gene(chr1, start1, end1, genes, gene_symbols)),
    anchor2_gene_info = list(get_nearest_gene(chr2, start2, end2, genes, gene_symbols)),
    anchor1_nearest_gene = anchor1_gene_info$gene,
    anchor1_gene_distance = anchor1_gene_info$distance,
    anchor2_nearest_gene = anchor2_gene_info$gene,
    anchor2_gene_distance = anchor2_gene_info$distance
  ) %>%
  select(-anchor1_gene_info, -anchor2_gene_info) %>%
  ungroup()

cat("  Gene annotations added.\n")

# ============================================================================
# Sort by Significance
# ============================================================================

cat("\n=== Sorting by Significance ===\n")

other_up_loops <- other_up_loops %>%
  arrange(FDR, desc(abs(logFC)))

cat("  Sorted by FDR (ascending), then |logFC| (descending).\n")

# ============================================================================
# Generate Output
# ============================================================================

cat("\n=== Generating Output Files ===\n")

# Select columns for PI table
pi_table <- other_up_loops %>%
  mutate(
    anchor1_coords = sprintf("%s:%d-%d", chr1, start1, end1),
    anchor2_coords = sprintf("%s:%d-%d", chr2, start2, end2)
  ) %>%
  select(
    loop_id,
    anchor1_coords, anchor2_coords,
    logFC, FDR, category,
    anchor1_type = anchor1_type_extended,
    anchor2_type = anchor2_type_extended,
    loop_type = loop_type_extended,
    other_category,
    anchor1_nearest_gene, anchor1_gene_distance,
    anchor2_nearest_gene, anchor2_gene_distance,
    anchor1_overlaps_CTCF, anchor2_overlaps_CTCF,
    loop_distance, resolution_kb
  )

# Save main output file
output_file <- file.path(OUTPUT_DIR, "other_category_up_loops.tsv")
write.table(pi_table, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: %s\n", output_file))

# Save full annotation file (with all original columns)
full_output_file <- file.path(OUTPUT_DIR, "other_up_loops_with_genes.tsv")
write.table(other_up_loops, full_output_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: %s\n", full_output_file))

# ============================================================================
# Generate Summary Report
# ============================================================================

summary_file <- file.path(OUTPUT_DIR, "other_up_loops_summary.txt")

summary_text <- sprintf("
================================================================================
SUMMARY: 'Other' Category Up Loops from Early Timepoint
Generated: %s
================================================================================

DATA OVERVIEW
-------------
Total early timepoint loops: %d
Up loops (up_in_mutant): %d
'Other' up loops (at least one anchor): %d
  - Both anchors 'Other': %d (%.1f%%)
  - Single anchor 'Other': %d (%.1f%%)

LOOP TYPE BREAKDOWN
-------------------
%s

CTCF OVERLAP ANALYSIS
---------------------
IMPORTANT: Early timepoint lacks CTCF ChIP-seq data, so 'Other' category
includes potential CTCF sites. We checked overlap with late CTCF peaks to
estimate how many 'Other' anchors may be CTCF-mediated.

Anchor1 'Other' overlapping CTCF: %d / %d (%.1f%%)
Anchor2 'Other' overlapping CTCF: %d / %d (%.1f%%)

KEY CAVEATS
-----------
1. Early 'Other' category is INFLATED (~37%%) due to missing CTCF data
   - Late timepoint with CTCF data shows only ~9%% 'Other'
   - An estimated 60-70%% of early 'Other' may be CTCF-mediated

2. The high 'Repressed_Promoter' in early (35%%) vs late (10%%) suggests
   early cerebellum has more Polycomb-repressed chromatin architecture

3. Multiple strong 'Other-Other' up loops cluster at chr8:71Mb region

TOP 10 MOST SIGNIFICANT 'OTHER' UP LOOPS
----------------------------------------
%s

================================================================================
",
  Sys.time(),
  nrow(loops),
  nrow(up_loops),
  nrow(other_up_loops),
  both_other, 100 * both_other / nrow(other_up_loops),
  single_other, 100 * single_other / nrow(other_up_loops),
  paste(capture.output(print(table(other_up_loops$loop_type_extended))), collapse = "\n"),
  ctcf_in_other_anchor1, nrow(other_anchor1), 100 * ctcf_in_other_anchor1 / max(nrow(other_anchor1), 1),
  ctcf_in_other_anchor2, nrow(other_anchor2), 100 * ctcf_in_other_anchor2 / max(nrow(other_anchor2), 1),
  paste(capture.output(print(head(pi_table[, c("loop_id", "anchor1_coords", "anchor2_coords",
                                                "logFC", "FDR", "loop_type",
                                                "anchor1_nearest_gene", "anchor2_nearest_gene")], 10),
                             width = 200)), collapse = "\n")
)

writeLines(summary_text, summary_file)
cat(sprintf("  Saved: %s\n", summary_file))

# ============================================================================
# Console Output for Quick Review
# ============================================================================

cat("\n")
cat("================================================================================\n")
cat("TOP 10 MOST SIGNIFICANT 'OTHER' UP LOOPS (for PI)\n")
cat("================================================================================\n\n")

top10 <- head(pi_table, 10)
for (i in 1:nrow(top10)) {
  row <- top10[i, ]
  cat(sprintf("%d. %s\n", i, row$loop_id))
  cat(sprintf("   %s <-> %s\n", row$anchor1_coords, row$anchor2_coords))
  cat(sprintf("   logFC: %.3f | FDR: %.2e | Category: %s\n", row$logFC, row$FDR, row$category))
  cat(sprintf("   Loop type: %s | Other type: %s\n", row$loop_type, row$other_category))
  cat(sprintf("   Anchor1 gene: %s (%s bp) | CTCF: %s\n",
              row$anchor1_nearest_gene, format(row$anchor1_gene_distance, big.mark = ","),
              ifelse(row$anchor1_overlaps_CTCF, "YES", "no")))
  cat(sprintf("   Anchor2 gene: %s (%s bp) | CTCF: %s\n",
              row$anchor2_nearest_gene, format(row$anchor2_gene_distance, big.mark = ","),
              ifelse(row$anchor2_overlaps_CTCF, "YES", "no")))
  cat("\n")
}

cat("================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n")
cat(sprintf("Total 'Other' up loops: %d\n", nrow(other_up_loops)))
cat(sprintf("  - Both anchors 'Other': %d\n", both_other))
cat(sprintf("  - Single anchor 'Other': %d\n", single_other))
cat(sprintf("\nCTCF overlap in 'Other' anchors: %.1f%% (Anchor1), %.1f%% (Anchor2)\n",
            100 * ctcf_in_other_anchor1 / max(nrow(other_anchor1), 1),
            100 * ctcf_in_other_anchor2 / max(nrow(other_anchor2), 1)))
cat("\nOutput files saved to:\n")
cat(sprintf("  - %s\n", output_file))
cat(sprintf("  - %s\n", full_output_file))
cat(sprintf("  - %s\n", summary_file))

cat("\n=== DONE ===\n")
