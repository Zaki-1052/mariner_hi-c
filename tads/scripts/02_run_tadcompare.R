# 02_run_tadcompare.R
# TADCompare differential boundary analysis: Control vs BAP1-KO mutant
# Uses merged contact matrices at 25kb resolution

suppressPackageStartupMessages({
    library(TADCompare)
    library(dplyr)
    library(readr)
})

# === Configuration ===
RESOLUTION <- 25000
BASE_DIR <- "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads"
INPUT_DIR <- file.path(BASE_DIR, "data/extracted/merged")
OUTPUT_DIR <- file.path(BASE_DIR, "results/tadcompare")

# Mouse autosomes
CHROMOSOMES <- paste0("chr", 1:19)

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=========================================\n")
cat("TADCompare Analysis: Control vs BAP1-KO\n")
cat("=========================================\n")
cat("Resolution:", RESOLUTION, "bp\n")
cat("Input directory:", INPUT_DIR, "\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("Chromosomes:", length(CHROMOSOMES), "\n")
cat("=========================================\n\n")

# === Function: Run TADCompare on one chromosome ===
run_tadcompare_chr <- function(chrom) {
    cat("Processing", chrom, "...\n")
    
    # Input file paths
    ctrl_file <- file.path(INPUT_DIR, paste0("ctrl_merged.", chrom, ".25kb.txt"))
    mut_file <- file.path(INPUT_DIR, paste0("mut_merged.", chrom, ".25kb.txt"))
    
    # Check files exist
    if (!file.exists(ctrl_file)) {
        cat("  ERROR: Control file not found:", ctrl_file, "\n")
        return(NULL)
    }
    if (!file.exists(mut_file)) {
        cat("  ERROR: Mutant file not found:", mut_file, "\n")
        return(NULL)
    }
    
    # Read sparse 3-column matrices
    # Format: region_i, region_j, contacts
    ctrl_mat <- tryCatch(
        read.table(ctrl_file, header = FALSE, col.names = c("i", "j", "counts")),
        error = function(e) NULL
    )
    mut_mat <- tryCatch(
        read.table(mut_file, header = FALSE, col.names = c("i", "j", "counts")),
        error = function(e) NULL
    )
    
    if (is.null(ctrl_mat) || is.null(mut_mat)) {
        cat("  ERROR: Failed to read matrix files\n")
        return(NULL)
    }
    
    # Check for sufficient data
    if (nrow(ctrl_mat) < 100 || nrow(mut_mat) < 100) {
        cat("  WARNING: Insufficient contacts (ctrl:", nrow(ctrl_mat), 
            ", mut:", nrow(mut_mat), ")\n")
        return(NULL)
    }
    
    cat("  Contacts - Control:", nrow(ctrl_mat), ", Mutant:", nrow(mut_mat), "\n")
    
    # Run TADCompare
    # Input: sparse 3-column matrices, TADCompare converts internally
    result <- tryCatch({
        TADCompare(
            cont_mat1 = ctrl_mat,
            cont_mat2 = mut_mat,
            resolution = RESOLUTION
        )
    }, error = function(e) {
        cat("  ERROR in TADCompare:", conditionMessage(e), "\n")
        return(NULL)
    })
    
    if (is.null(result)) {
        return(NULL)
    }
    
    # Extract TAD_Frame (boundaries with classifications)
    tad_frame <- result$TAD_Frame
    
    if (is.null(tad_frame) || nrow(tad_frame) == 0) {
        cat("  WARNING: No boundaries detected\n")
        return(NULL)
    }
    
    # Add chromosome column
    tad_frame$chr <- chrom
    
    # Reorder columns for clarity
    tad_frame <- tad_frame %>%
        select(chr, Boundary, everything())
    
    # Summary stats
    n_total <- nrow(tad_frame)
    n_diff <- sum(tad_frame$Differential == "Differential", na.rm = TRUE)
    
    cat("  Boundaries detected:", n_total, "\n")
    cat("  Differential boundaries:", n_diff, "\n")
    
    # Type breakdown (excluding NA)
    type_counts <- table(tad_frame$Type, useNA = "no")
    if (length(type_counts) > 0) {
        cat("  Types:", paste(names(type_counts), type_counts, sep = "=", collapse = ", "), "\n")
    }
    
    cat("  ✓ Complete\n\n")
    
    return(tad_frame)
}

# === Main analysis loop ===
cat("Starting TADCompare analysis across all chromosomes...\n\n")

results_list <- list()
failed_chroms <- c()

for (chrom in CHROMOSOMES) {
    result <- run_tadcompare_chr(chrom)
    
    if (!is.null(result)) {
        results_list[[chrom]] <- result
    } else {
        failed_chroms <- c(failed_chroms, chrom)
    }
}

# === Combine results ===
cat("=========================================\n")
cat("Combining results...\n")

if (length(results_list) == 0) {
    cat("ERROR: No results to combine. Check input files.\n")
    quit(status = 1)
}

# Combine all chromosomes
all_boundaries <- bind_rows(results_list)

cat("Total boundaries across all chromosomes:", nrow(all_boundaries), "\n")

# === Summary statistics ===
cat("\n=========================================\n")
cat("Summary Statistics\n")
cat("=========================================\n")

# Overall differential status
diff_summary <- all_boundaries %>%
    count(Differential) %>%
    mutate(pct = round(n / sum(n) * 100, 1))

cat("\nDifferential status:\n")
print(diff_summary)

# Type breakdown
type_summary <- all_boundaries %>%
    filter(!is.na(Type)) %>%
    count(Type) %>%
    arrange(desc(n)) %>%
    mutate(pct = round(n / sum(n) * 100, 1))

cat("\nBoundary type classification:\n")
print(type_summary)

# Enrichment direction
enrichment_summary <- all_boundaries %>%
    filter(Differential == "Differential") %>%
    count(Enriched_In) %>%
    mutate(pct = round(n / sum(n) * 100, 1))

cat("\nDifferential boundaries by enrichment:\n")
print(enrichment_summary)

# Per-chromosome summary
chr_summary <- all_boundaries %>%
    group_by(chr) %>%
    summarise(
        total = n(),
        differential = sum(Differential == "Differential", na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(pct_diff = round(differential / total * 100, 1))

cat("\nPer-chromosome summary:\n")
print(chr_summary, n = 19)

# === Save outputs ===
cat("\n=========================================\n")
cat("Saving outputs...\n")

# Full results table
output_file <- file.path(OUTPUT_DIR, "tadcompare_all_boundaries.tsv")
write_tsv(all_boundaries, output_file)
cat("Full results:", output_file, "\n")

# Differential boundaries only
diff_boundaries <- all_boundaries %>%
    filter(Differential == "Differential")

output_diff <- file.path(OUTPUT_DIR, "tadcompare_differential_only.tsv")
write_tsv(diff_boundaries, output_diff)
cat("Differential only:", output_diff, "\n")

# BED format for IGV (all boundaries)
# BED: chr, start, end, name, score, strand
bed_all <- all_boundaries %>%
    mutate(
        start = Boundary,
        end = Boundary + RESOLUTION,
        name = paste(Type, Enriched_In, sep = "_"),
        score = round(abs(Gap_Score) * 100),  # Scale for visualization
        strand = "."
    ) %>%
    select(chr, start, end, name, score, strand, 
           Gap_Score, TAD_Score1, TAD_Score2, Differential, Enriched_In, Type)

output_bed <- file.path(OUTPUT_DIR, "tadcompare_boundaries.bed")
write_tsv(bed_all, output_bed, col_names = FALSE)
cat("BED file (all):", output_bed, "\n")

# BED format for differential boundaries only
bed_diff <- bed_all %>%
    filter(Differential == "Differential")

output_bed_diff <- file.path(OUTPUT_DIR, "tadcompare_differential.bed")
write_tsv(bed_diff, output_bed_diff, col_names = FALSE)
cat("BED file (differential):", output_bed_diff, "\n")

# Summary statistics file
summary_output <- file.path(OUTPUT_DIR, "tadcompare_summary.txt")
sink(summary_output)
cat("TADCompare Analysis Summary\n")
cat("===========================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Resolution:", RESOLUTION, "bp\n")
cat("Chromosomes analyzed:", length(results_list), "/", length(CHROMOSOMES), "\n")
if (length(failed_chroms) > 0) {
    cat("Failed chromosomes:", paste(failed_chroms, collapse = ", "), "\n")
}
cat("\nDifferential Status:\n")
print(diff_summary)
cat("\nBoundary Types:\n")
print(type_summary)
cat("\nEnrichment Direction:\n")
print(enrichment_summary)
cat("\nPer-Chromosome:\n")
print(chr_summary, n = 19)
sink()
cat("Summary file:", summary_output, "\n")

# === Done ===
cat("\n=========================================\n")
cat("✅ TADCompare analysis complete\n")
cat("=========================================\n")
cat("Chromosomes processed:", length(results_list), "\n")
if (length(failed_chroms) > 0) {
    cat("Failed chromosomes:", paste(failed_chroms, collapse = ", "), "\n")
}
cat("\nNext step: Run 03_run_consensus.R for replicate robustness check\n")
