# 04_analyze_shifts.R
# Calculate shift distances and generate final annotated outputs
# Requires: bedtools closest output files from 04_postprocess.sb
# Usage: Rscript scripts/04_analyze_shifts.R [timepoint]
#   timepoint: "late" (default) or "early"

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

# === Parse Arguments ===
args <- commandArgs(trailingOnly = TRUE)
TIMEPOINT <- if (length(args) >= 1) args[1] else "late"

# Validate timepoint
if (!TIMEPOINT %in% c("late", "early")) {
    stop("ERROR: Unknown timepoint '", TIMEPOINT, "'. Use 'late' or 'early'.")
}

cat("Timepoint:", TIMEPOINT, "\n\n")

# === Configuration ===
RESOLUTION <- 25000
BASE_DIR <- "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads"
TMP_DIR <- file.path(BASE_DIR, "results", TIMEPOINT, "final/tmp")
FINAL_DIR <- file.path(BASE_DIR, "results", TIMEPOINT, "final")
CONSENSUS_DIR <- file.path(BASE_DIR, "results", TIMEPOINT, "consensus")
TADCOMPARE_DIR <- file.path(BASE_DIR, "results", TIMEPOINT, "tadcompare")

cat("=========================================\n")
cat("Shift Distance Analysis\n")
cat("=========================================\n")

# === Load bedtools closest output ===
cat("Loading distance files from bedtools...\n")

# bedtools closest output: chr, start, end, id, chr2, start2, end2, id2, distance
col_names <- c("chr", "start", "end", "id", "chr2", "start2", "end2", "id2", "distance")

ctrl_to_mut <- tryCatch(
    read.table(file.path(TMP_DIR, "ctrl_to_mut_distances.txt"), 
               header = FALSE, sep = "\t", col.names = col_names),
    error = function(e) {
        cat("  No ctrl→mut distances found\n")
        data.frame()
    }
)

mut_to_ctrl <- tryCatch(
    read.table(file.path(TMP_DIR, "mut_to_ctrl_distances.txt"),
               header = FALSE, sep = "\t", col.names = col_names),
    error = function(e) {
        cat("  No mut→ctrl distances found\n
")
        data.frame()
    }
)

cat("  Control→Mutant pairs:", nrow(ctrl_to_mut), "\n")
cat("  Mutant→Control pairs:", nrow(mut_to_ctrl), "\n")

# === Load TADCompare results ===
cat("\nLoading TADCompare results...\n")

if (file.exists(file.path(CONSENSUS_DIR, "tadcompare_with_robustness.tsv"))) {
    tadcompare <- read_tsv(file.path(CONSENSUS_DIR, "tadcompare_with_robustness.tsv"), 
                           show_col_types = FALSE)
    has_robustness <- TRUE
    cat("  Using robustness-annotated file\n")
} else {
    tadcompare <- read_tsv(file.path(TADCOMPARE_DIR, "tadcompare_all_boundaries.tsv"),
                           show_col_types = FALSE)
    has_robustness <- FALSE
    cat("  Using original TADCompare file\n")
}

cat("  Total boundaries:", nrow(tadcompare), "\n")

# === Merge distances back to main results ===
cat("\nMerging shift distances...\n")

# Create distance lookup tables
ctrl_distances <- ctrl_to_mut %>%
    select(chr, boundary = start, shift_distance = distance)

mut_distances <- mut_to_ctrl %>%
    select(chr, boundary = start, shift_distance = distance)

all_distances <- bind_rows(ctrl_distances, mut_distances)

# Join to main results
tadcompare_final <- tadcompare %>%
    left_join(all_distances, by = c("chr", "Boundary" = "boundary"))

# Only keep shift distance for "Shifted" boundaries
tadcompare_final <- tadcompare_final %>%
    mutate(
        shift_distance = ifelse(Type == "Shifted", shift_distance, NA_real_),
        shift_distance_kb = shift_distance / 1000
    )

# === Shift Distance Statistics ===
cat("\n=========================================\n")
cat("Shift Distance Summary\n")
cat("=========================================\n")

shifted_only <- tadcompare_final %>% 
    filter(Type == "Shifted", !is.na(shift_distance))

cat("\nShifted boundaries with distances:", nrow(shifted_only), "/", 
    sum(tadcompare_final$Type == "Shifted", na.rm = TRUE), "\n")

if (nrow(shifted_only) > 0) {
    cat("\nDistribution (kb):\n")
    cat("  Min:    ", round(min(shifted_only$shift_distance_kb), 1), "\n")
    cat("  Q1:     ", round(quantile(shifted_only$shift_distance_kb, 0.25), 1), "\n")
    cat("  Median: ", round(median(shifted_only$shift_distance_kb), 1), "\n")
    cat("  Mean:   ", round(mean(shifted_only$shift_distance_kb), 1), "\n")
    cat("  Q3:     ", round(quantile(shifted_only$shift_distance_kb, 0.75), 1), "\n")
    cat("  Max:    ", round(max(shifted_only$shift_distance_kb), 1), "\n")
    
    # Binned distribution
    cat("\nBinned distribution:\n")
    bins <- shifted_only %>%
        mutate(bin = cut(shift_distance_kb, 
                        breaks = c(0, 25, 50, 100, 200, 500, Inf),
                        labels = c("0-25kb", "25-50kb", "50-100kb", 
                                  "100-200kb", "200-500kb", ">500kb"))) %>%
        count(bin) %>%
        mutate(pct = round(n / sum(n) * 100, 1))
    print(bins)
    
    # By enrichment direction
    cat("\nBy enrichment:\n")
    by_enrichment <- shifted_only %>%
        group_by(Enriched_In) %>%
        summarise(
            n = n(),
            median_kb = round(median(shift_distance_kb), 1),
            mean_kb = round(mean(shift_distance_kb), 1),
            .groups = "drop"
        )
    print(by_enrichment)
}

# === Save outputs ===
cat("\n=========================================\n")
cat("Saving outputs...\n")

# Full annotated TSV
output_tsv <- file.path(FINAL_DIR, "tadcompare_final_annotated.tsv")
write_tsv(tadcompare_final, output_tsv)
cat("Full results:", output_tsv, "\n")

# BED format helper
make_bed <- function(df) {
    df %>%
        mutate(
            start = Boundary,
            end = Boundary + RESOLUTION,
            name = paste(Type, Enriched_In, sep = "_"),
            score = pmin(round(abs(Gap_Score) * 100), 1000),
            strand = "."
        ) %>%
        select(chr, start, end, name, score, strand,
               Gap_Score, TAD_Score1, TAD_Score2, Differential, 
               Enriched_In, Type, shift_distance_kb)
}

# All boundaries BED
bed_all <- make_bed(tadcompare_final)
write_tsv(bed_all, file.path(FINAL_DIR, "tadcompare_final.bed"), col_names = FALSE)
cat("All boundaries BED:", file.path(FINAL_DIR, "tadcompare_final.bed"), "\n")

# Differential only BED
bed_diff <- bed_all %>% filter(Differential == "Differential")
write_tsv(bed_diff, file.path(FINAL_DIR, "differential_boundaries_final.bed"), col_names = FALSE)
cat("Differential BED:", file.path(FINAL_DIR, "differential_boundaries_final.bed"), "\n")

# Shifted only BED
bed_shifted <- bed_all %>% filter(Type == "Shifted")
write_tsv(bed_shifted, file.path(FINAL_DIR, "shifted_boundaries.bed"), col_names = FALSE)
cat("Shifted BED:", file.path(FINAL_DIR, "shifted_boundaries.bed"), "\n")

# Summary text file
summary_file <- file.path(FINAL_DIR, "analysis_summary.txt")
sink(summary_file)
cat("TADCompare Final Analysis Summary\n")
cat("=================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Resolution:", RESOLUTION, "bp\n\n")
cat("Total boundaries:", nrow(tadcompare_final), "\n")
cat("Differential:", sum(tadcompare_final$Differential == "Differential", na.rm = TRUE), "\n")
cat("Shifted:", sum(tadcompare_final$Type == "Shifted", na.rm = TRUE), "\n")
if (nrow(shifted_only) > 0) {
    cat("\nShift distances (kb):\n")
    cat("  Median:", round(median(shifted_only$shift_distance_kb), 1), "\n")
    cat("  Mean:", round(mean(shifted_only$shift_distance_kb), 1), "\n")
}
sink()
cat("Summary:", summary_file, "\n")

cat("\n=========================================\n")
cat("✅ Shift distance analysis complete\n")
cat("=========================================\n")
