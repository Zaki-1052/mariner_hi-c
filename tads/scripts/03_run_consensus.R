# 03_run_consensus.R
# ConsensusTADs analysis for replicate robustness assessment
# Identifies consistently-detected boundaries across biological replicates

suppressPackageStartupMessages({
    library(TADCompare)
    library(dplyr)
    library(readr)
})

# === Configuration ===
RESOLUTION <- 25000
BASE_DIR <- "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads"
INPUT_DIR <- file.path(BASE_DIR, "data/extracted/replicates")
OUTPUT_DIR <- file.path(BASE_DIR, "results/consensus")
TADCOMPARE_DIR <- file.path(BASE_DIR, "results/tadcompare")

CHROMOSOMES <- paste0("chr", 1:19)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=========================================\n")
cat("ConsensusTADs: Replicate Robustness Analysis\n")
cat("=========================================\n")
cat("Resolution:", RESOLUTION, "bp\n")
cat("Input directory:", INPUT_DIR, "\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("=========================================\n\n")

# === Function: Load sparse matrix ===
load_sparse_matrix <- function(filepath) {
    if (!file.exists(filepath)) return(NULL)
    tryCatch(
        read.table(filepath, header = FALSE, col.names = c("i", "j", "counts")),
        error = function(e) NULL
    )
}

# === Function: Run ConsensusTADs for one condition/chromosome ===
run_consensus_chr <- function(chrom, condition, rep_names) {
    # Load all replicate matrices for this chromosome
    matrices <- lapply(rep_names, function(rep) {
        filepath <- file.path(INPUT_DIR, paste0(rep, ".", chrom, ".25kb.txt"))
        load_sparse_matrix(filepath)
    })
    names(matrices) <- rep_names
    
    # Check all matrices loaded
    loaded <- !sapply(matrices, is.null)
    if (sum(loaded) < 2) {
        cat("    WARNING: Only", sum(loaded), "replicates loaded, need at least 2\n")
        return(NULL)
    }
    
    # Filter to successfully loaded matrices
    matrices <- matrices[loaded]
    
    # Check for sufficient data in each
    sufficient <- sapply(matrices, function(m) nrow(m) > 100)
    if (sum(sufficient) < 2) {
        cat("    WARNING: Insufficient contacts in replicates\n")
        return(NULL)
    }
    matrices <- matrices[sufficient]
    
    cat("    Replicates used:", paste(names(matrices), collapse = ", "), "\n")
    
    # Run ConsensusTADs
    result <- tryCatch({
        ConsensusTADs(matrices, resolution = RESOLUTION)
    }, error = function(e) {
        cat("    ERROR:", conditionMessage(e), "\n")
        return(NULL)
    })
    
    if (is.null(result)) return(NULL)
    
    # Extract consensus boundaries (significant ones)
    consensus <- result$Consensus
    if (is.null(consensus) || nrow(consensus) == 0) {
        cat("    No consensus boundaries detected\n")
        return(NULL)
    }
    
    # Add metadata
    consensus$chr <- chrom
    consensus$condition <- condition
    consensus$n_replicates <- length(matrices)
    
    cat("    Consensus boundaries:", nrow(consensus), "\n")
    
    return(consensus)
}

# === Run ConsensusTADs for Control replicates ===
cat("Processing CONTROL replicates...\n")
cat("================================\n")

ctrl_reps <- c("ctrl_M1", "ctrl_M2", "ctrl_M3")
ctrl_results <- list()

for (chrom in CHROMOSOMES) {
    cat("  ", chrom, ":\n")
    result <- run_consensus_chr(chrom, "control", ctrl_reps)
    if (!is.null(result)) {
        ctrl_results[[chrom]] <- result
    }
}

cat("\n")

# === Run ConsensusTADs for Mutant replicates ===
cat("Processing MUTANT replicates...\n")
cat("===============================\n")

mut_reps <- c("mut_M1", "mut_M2", "mut_M3")
mut_results <- list()

for (chrom in CHROMOSOMES) {
    cat("  ", chrom, ":\n")
    result <- run_consensus_chr(chrom, "mutant", mut_reps)
    if (!is.null(result)) {
        mut_results[[chrom]] <- result
    }
}

cat("\n")

# === Combine results ===
cat("=========================================\n")
cat("Combining and analyzing results...\n")
cat("=========================================\n")

ctrl_consensus <- bind_rows(ctrl_results)
mut_consensus <- bind_rows(mut_results)

cat("Control consensus boundaries:", nrow(ctrl_consensus), "\n")
cat("Mutant consensus boundaries:", nrow(mut_consensus), "\n")

# === Load TADCompare results for integration ===
cat("\nIntegrating with TADCompare results...\n")

tadcompare_file <- file.path(TADCOMPARE_DIR, "tadcompare_all_boundaries.tsv")
if (file.exists(tadcompare_file)) {
    tadcompare_results <- read_tsv(tadcompare_file, show_col_types = FALSE)
    cat("TADCompare boundaries loaded:", nrow(tadcompare_results), "\n")
} else {
    cat("WARNING: TADCompare results not found at", tadcompare_file, "\n")
    tadcompare_results <- NULL
}

# === Create robustness scores ===
# A boundary is "robust" if it appears in consensus for that condition

if (!is.null(tadcompare_results)) {
    # Create lookup keys for consensus boundaries
    # Use a tolerance window (resolution/2) for matching
    tolerance <- RESOLUTION / 2
    
    ctrl_consensus_coords <- ctrl_consensus %>%
        mutate(key = paste(chr, Coordinate, sep = "_")) %>%
        pull(key)
    
    mut_consensus_coords <- mut_consensus %>%
        mutate(key = paste(chr, Coordinate, sep = "_")) %>%
        pull(key)
    
    # Function to check if boundary is near any consensus boundary
    is_robust <- function(chr_val, boundary_val, consensus_df, tol) {
        subset_df <- consensus_df %>% filter(chr == chr_val)
        if (nrow(subset_df) == 0) return(FALSE)
        any(abs(subset_df$Coordinate - boundary_val) <= tol)
    }
    
    # Add robustness flags to TADCompare results
    tadcompare_annotated <- tadcompare_results %>%
        rowwise() %>%
        mutate(
            robust_in_ctrl = is_robust(chr, Boundary, ctrl_consensus, tolerance),
            robust_in_mut = is_robust(chr, Boundary, mut_consensus, tolerance),
            robustness = case_when(
                robust_in_ctrl & robust_in_mut ~ "Both",
                robust_in_ctrl ~ "Control_only",
                robust_in_mut ~ "Mutant_only",
                TRUE ~ "Neither"
            )
        ) %>%
        ungroup()
    
    # Summary of robustness
    cat("\nRobustness assessment:\n")
    robustness_summary <- tadcompare_annotated %>%
        count(robustness) %>%
        mutate(pct = round(n / sum(n) * 100, 1))
    print(robustness_summary)
    
    # Robustness by differential status
    cat("\nRobustness by differential status:\n")
    robust_by_diff <- tadcompare_annotated %>%
        group_by(Differential) %>%
        summarise(
            total = n(),
            robust_both = sum(robustness == "Both"),
            robust_ctrl = sum(robustness == "Control_only"),
            robust_mut = sum(robustness == "Mutant_only"),
            not_robust = sum(robustness == "Neither"),
            pct_robust = round((robust_both + robust_ctrl + robust_mut) / total * 100, 1),
            .groups = "drop"
        )
    print(robust_by_diff)
    
    # High-confidence differential boundaries (differential AND robust)
    high_conf_diff <- tadcompare_annotated %>%
        filter(Differential == "Differential", robustness != "Neither")
    
    cat("\nHigh-confidence differential boundaries:", nrow(high_conf_diff), "\n")
    cat("(Differential AND detected in replicate consensus)\n")
}

# === Save outputs ===
cat("\n=========================================\n")
cat("Saving outputs...\n")

# Control consensus
ctrl_out <- file.path(OUTPUT_DIR, "consensus_control.tsv")
write_tsv(ctrl_consensus, ctrl_out)
cat("Control consensus:", ctrl_out, "\n")

# Mutant consensus
mut_out <- file.path(OUTPUT_DIR, "consensus_mutant.tsv")
write_tsv(mut_consensus, mut_out)
cat("Mutant consensus:", mut_out, "\n")

# Annotated TADCompare results (if available)
if (!is.null(tadcompare_results)) {
    annotated_out <- file.path(OUTPUT_DIR, "tadcompare_with_robustness.tsv")
    write_tsv(tadcompare_annotated, annotated_out)
    cat("Annotated TADCompare:", annotated_out, "\n")
    
    # High-confidence differential boundaries
    highconf_out <- file.path(OUTPUT_DIR, "high_confidence_differential.tsv")
    write_tsv(high_conf_diff, highconf_out)
    cat("High-confidence differential:", highconf_out, "\n")
    
    # BED file for high-confidence differential
    bed_highconf <- high_conf_diff %>%
        mutate(
            start = Boundary,
            end = Boundary + RESOLUTION,
            name = paste(Type, Enriched_In, robustness, sep = "_"),
            score = round(abs(Gap_Score) * 100),
            strand = "."
        ) %>%
        select(chr, start, end, name, score, strand,
               Gap_Score, TAD_Score1, TAD_Score2, Type, Enriched_In, robustness)
    
    bed_out <- file.path(OUTPUT_DIR, "high_confidence_differential.bed")
    write_tsv(bed_highconf, bed_out, col_names = FALSE)
    cat("High-confidence BED:", bed_out, "\n")
}

# === Summary file ===
summary_out <- file.path(OUTPUT_DIR, "consensus_summary.txt")
sink(summary_out)
cat("ConsensusTADs Analysis Summary\n")
cat("==============================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Resolution:", RESOLUTION, "bp\n\n")
cat("Consensus Boundaries:\n")
cat("  Control:", nrow(ctrl_consensus), "\n")
cat("  Mutant:", nrow(mut_consensus), "\n\n")
if (!is.null(tadcompare_results)) {
    cat("Robustness Assessment:\n")
    print(robustness_summary)
    cat("\nBy Differential Status:\n")
    print(robust_by_diff)
    cat("\nHigh-confidence differential:", nrow(high_conf_diff), "\n")
}
sink()
cat("Summary:", summary_out, "\n")

# === Done ===
cat("\n=========================================\n")
cat("✅ ConsensusTADs analysis complete\n")
cat("=========================================\n")
cat("\nKey outputs:\n")
cat("  - consensus_control.tsv: Robust boundaries in control\n")
cat("  - consensus_mutant.tsv: Robust boundaries in mutant\n")
cat("  - tadcompare_with_robustness.tsv: TADCompare + robustness flags\n")
cat("  - high_confidence_differential.tsv/bed: Filtered high-confidence set\n")
cat("\nNext step: Run 04_postprocess.sh for shift distance calculation\n")
