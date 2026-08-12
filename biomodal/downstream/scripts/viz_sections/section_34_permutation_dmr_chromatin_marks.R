# biomodal/downstream/scripts/viz_sections/section_34_permutation_dmr_chromatin_marks.R
# Section 34: DMR x Chromatin Mark Interval Permutation Testing
# Validates Fisher's exact tests from sections 12a, 14a, 14b, 19f using
# regioneReloaded crosswisePermTest to assess spatial co-occurrence of DMR
# intervals with chromatin mark peaks while accounting for genomic clustering.
#
# Figures:
#   34a: Crosswise z-score heatmap (2 x 8 matrix)
#   34b: Fisher-vs-permutation comparison forest plot
#   34c: Local z-score curve for strongest association
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_34_permutation_dmr_chromatin_marks.R

source("scripts/viz_sections/_shared_config.R")

# Skip if cached results exist (avoids re-run from run_all_sections.sh glob)
# FORCE_RERUN env var (set by run_permutation.sb) bypasses this check
.cache_check <- file.path(TABLES_DIR, "permutation_34_dmr_marks.rds")
if (file.exists(.cache_check) && Sys.getenv("FORCE_RERUN", unset = "") == "") {
  cat("Section 34: Cached RDS found at", .cache_check, "-- skipping.\n")
  cat("  To re-run: use run_permutation.sb or export FORCE_RERUN=1\n")
  quit(save = "no", status = 0)
}
.force_rerun <- Sys.getenv("FORCE_RERUN", unset = "") != ""

cat("\n")
cat("================================================================================\n")
cat("SECTION 34: DMR x CHROMATIN MARK INTERVAL PERMUTATION TESTING\n")
cat("================================================================================\n\n")

# =============================================================================
# BLOCK 1: SECTION-SPECIFIC PACKAGES
# =============================================================================

suppressPackageStartupMessages({
  library(regioneR)
  library(regioneReloaded)
  library(BSgenome.Mmusculus.UCSC.mm10)
})
patch_chooseHclustMet()

cat(sprintf("regioneR version: %s\n", packageVersion("regioneR")))
cat(sprintf("regioneReloaded version: %s\n", packageVersion("regioneReloaded")))


# =============================================================================
# BLOCK 2: PERMUTATION PARAMETERS
# =============================================================================

# Accept ntimes from command line (e.g., Rscript section_34_... 5000)
args <- commandArgs(trailingOnly = TRUE)
PERM_NTIMES <- if (length(args) >= 1 && !is.na(as.numeric(args[1]))) as.numeric(args[1]) else 5000

# Use SLURM_CPUS_PER_TASK if available, else default to 8
slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
PERM_CORES <- if (slurm_cpus != "") as.integer(slurm_cpus) else 8L

PERM_PER_CHR  <- TRUE
PERM_RANFUN   <- "randomizeRegions"
PERM_EVFUN    <- "numOverlaps"
PERM_SEED     <- 42

SECTION_OUTPUT <- file.path(OUTPUT_DIR, "34_permutation_dmr_chromatin_marks")
dir.create(SECTION_OUTPUT, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("\nPermutation parameters: ntimes=%d, cores=%d, ranFUN=%s, per.chr=%s\n",
            PERM_NTIMES, PERM_CORES, PERM_RANFUN, PERM_PER_CHR))

# =============================================================================
# BLOCK 3: GENOME OBJECT
# =============================================================================

cat("\nConstructing mm10 genome object (standard chromosomes)...\n")

genome_full <- getGenomeAndMask("mm10")$genome
standard_chrs <- paste0("chr", c(1:19, "X"))
genome <- genome_full[seqnames(genome_full) %in% standard_chrs]
stopifnot(length(genome) == 20)

cat(sprintf("  Genome: %d chromosomes, total length %.1f Gb\n",
            length(genome), sum(as.numeric(width(genome))) / 1e9))

# =============================================================================
# BLOCK 4: ALIST CONSTRUCTION (2 DMR direction sets)
# =============================================================================

cat("\nConstructing Alist (DMR direction sets)...\n")

hyper_gr <- dmr_to_granges(mc_dmr %>% dplyr::filter(significant, mod_difference > 0))
hypo_gr  <- dmr_to_granges(mc_dmr %>% dplyr::filter(significant, mod_difference < 0))

stopifnot(length(hyper_gr) > 0)
stopifnot(length(hypo_gr) > 0)

# Restrict to standard chromosomes and harmonize seqlevels with genome
restrict_to_standard <- function(gr) {
  gr <- gr[seqnames(gr) %in% standard_chrs]
  seqlevels(gr) <- seqlevels(genome)
  gr
}

hyper_gr <- restrict_to_standard(hyper_gr)
hypo_gr  <- restrict_to_standard(hypo_gr)

Alist <- list(
  "mC Hyper DMRs" = hyper_gr,
  "mC Hypo DMRs"  = hypo_gr
)

cat(sprintf("  mC Hyper DMRs: %d regions\n", length(hyper_gr)))
cat(sprintf("  mC Hypo DMRs:  %d regions\n", length(hypo_gr)))

# =============================================================================
# BLOCK 5: BLIST CONSTRUCTION (8 chromatin mark peak sets)
# =============================================================================

cat("\nConstructing Blist (chromatin mark peak sets)...\n")

# Derive differential peaks from condition-specific peak sets
# Gained = mut-only, Lost = ctrl-only (matching section_14 pattern)
derive_differential_peaks <- function(ctrl_gr, mut_gr) {
  mut_hits <- countOverlaps(mut_gr, ctrl_gr)
  ctrl_hits <- countOverlaps(ctrl_gr, mut_gr)
  list(
    gained = mut_gr[mut_hits == 0],
    lost   = ctrl_gr[ctrl_hits == 0],
    shared = mut_gr[mut_hits > 0]
  )
}

# ATAC differential peaks
atac_up   <- load_chip_peaks(ATAC_FILES$up, "ATAC Up")
atac_down <- load_chip_peaks(ATAC_FILES$down, "ATAC Down")
stopifnot(!is.null(atac_up), !is.null(atac_down))

# K119ub condition-specific peaks
k119ub_ctrl <- load_chip_peaks(K119UB_FILES$ctrl, "K119ub Ctrl")
k119ub_mut  <- load_chip_peaks(K119UB_FILES$mut, "K119ub Mut")
stopifnot(!is.null(k119ub_ctrl), !is.null(k119ub_mut))

# Derived gained/lost K119ub peaks
diff_k119ub   <- derive_differential_peaks(k119ub_ctrl, k119ub_mut)
k119ub_gained <- diff_k119ub$gained
k119ub_lost   <- diff_k119ub$lost

cat(sprintf("  K119ub Gained (mut-only): %d peaks\n", length(k119ub_gained)))
cat(sprintf("  K119ub Lost (ctrl-only):  %d peaks\n", length(k119ub_lost)))

# H3K27ac condition-specific peaks
h3k27ac_ctrl <- load_chip_peaks(H3K27AC_FILES$ctrl, "H3K27ac Ctrl")
h3k27ac_mut  <- load_chip_peaks(H3K27AC_FILES$mut, "H3K27ac Mut")
stopifnot(!is.null(h3k27ac_ctrl), !is.null(h3k27ac_mut))

# Build Blist with seqlevel harmonization
Blist <- list(
  "ATAC Up"       = restrict_to_standard(atac_up),
  "ATAC Down"     = restrict_to_standard(atac_down),
  "K119ub Ctrl"   = restrict_to_standard(k119ub_ctrl),
  "K119ub Mut"    = restrict_to_standard(k119ub_mut),
  "K119ub Gained" = restrict_to_standard(k119ub_gained),
  "K119ub Lost"   = restrict_to_standard(k119ub_lost),
  "H3K27ac Ctrl"  = restrict_to_standard(h3k27ac_ctrl),
  "H3K27ac Mut"   = restrict_to_standard(h3k27ac_mut)
)

# Validate all Blist entries have regions
for (nm in names(Blist)) {
  stopifnot(length(Blist[[nm]]) > 0)
  cat(sprintf("  %s: %d regions (after chr filter)\n", nm, length(Blist[[nm]])))
}

# =============================================================================
# BLOCK 6: RDS CACHE + crosswisePermTest
# =============================================================================

cache_path <- file.path(TABLES_DIR, "permutation_34_dmr_marks.rds")

if (file.exists(cache_path) && !.force_rerun) {
  cat("\nLoading cached permutation results from:", cache_path, "\n")
  cw_34 <- readRDS(cache_path)
} else {
  cat(sprintf("\nRunning crosswisePermTest (ntimes=%d, cores=%d)...\n",
              PERM_NTIMES, PERM_CORES))
  cat("  Estimated runtime: 20-30 minutes\n")

  set.seed(PERM_SEED)
  cw_34 <- crosswisePermTest(
    Alist          = Alist,
    Blist          = Blist,
    genome         = genome,
    ranFUN         = PERM_RANFUN,
    evFUN          = PERM_EVFUN,
    ntimes         = PERM_NTIMES,
    mc.cores       = PERM_CORES,
    per.chromosome = PERM_PER_CHR
  )

  cat("  crosswisePermTest complete. Making crosswise matrix...\n")

  # symm_matrix=FALSE: non-square Alist(2) x Blist(8) design
  # hc.method="average": single method (patched chooseHclustMet handles <=2 dim)
  cw_34 <- makeCrosswiseMatrix(cw_34, pvcut = 1,
                                symm_matrix = FALSE, hc.method = "average")

  saveRDS(cw_34, cache_path)
  cat("  Saved cache:", cache_path, "\n")
}

# =============================================================================
# BLOCK 7: EXTRACT PER-CELL PERMUTATION RESULTS
# =============================================================================

cat("\nExtracting per-cell permutation results...\n")

# Access multiOverlaps slot: list of data.frames, one per Alist entry
multi_overlaps <- cw_34@multiOverlaps

perm_results <- do.call(rbind, lapply(names(multi_overlaps), function(a_name) {
  df <- multi_overlaps[[a_name]]
  df$dmr_set <- a_name
  df
}))

# Create readable test ID for each cell
perm_results$test_id <- paste(perm_results$dmr_set, "x", perm_results$name)

cat(sprintf("  Extracted %d pairwise test results\n", nrow(perm_results)))
cat(sprintf("  Z-score range: [%.2f, %.2f]\n",
            min(perm_results$z_score), max(perm_results$z_score)))

# =============================================================================
# BLOCK 8: FIGURE 34a - Crosswise Z-Score Heatmap
# =============================================================================

cat("\n--- FIGURE 34a: Crosswise Z-Score Heatmap ---\n")

p_34a <- plotCrosswiseMatrix(cw_34, matrix_type = "association")

# plotCrosswiseMatrix returns a ggplot object
if (is.ggplot(p_34a)) {
  save_multiformat_ggplot(p_34a, file.path(SECTION_OUTPUT, "34a_crosswise_dmr_x_marks"),
                          width = 12, height = 6)
} else {
  save_multiformat_base(
    quote(plotCrosswiseMatrix(cw_34, matrix_type = "association")),
    file.path(SECTION_OUTPUT, "34a_crosswise_dmr_x_marks"),
    width = 12, height = 6
  )
}

cat("  Saved: 34a_crosswise_dmr_x_marks\n")

# =============================================================================
# BLOCK 9: RECOMPUTE ORIGINAL FISHER TESTS INLINE
# =============================================================================

cat("\n--- Recomputing original Fisher tests for comparison ---\n")

# Test 1: Section 12a - DMR direction x ATAC direction
fisher_12a_mat <- matrix(
  c(sum(countOverlaps(Alist[["mC Hyper DMRs"]], Blist[["ATAC Up"]]) > 0),
    sum(countOverlaps(Alist[["mC Hyper DMRs"]], Blist[["ATAC Down"]]) > 0),
    sum(countOverlaps(Alist[["mC Hypo DMRs"]], Blist[["ATAC Up"]]) > 0),
    sum(countOverlaps(Alist[["mC Hypo DMRs"]], Blist[["ATAC Down"]]) > 0)),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("Hyper", "Hypo"), c("ATAC Up", "ATAC Down"))
)
fisher_12a <- fisher.test(fisher_12a_mat)
cat(sprintf("  12a (ATAC):    OR = %.2f, p = %.2e\n",
            fisher_12a$estimate, fisher_12a$p.value))

# Test 2: Section 14a - DMR direction x K119ub condition
fisher_14a_mat <- matrix(
  c(sum(countOverlaps(Alist[["mC Hyper DMRs"]], Blist[["K119ub Ctrl"]]) > 0),
    sum(countOverlaps(Alist[["mC Hyper DMRs"]], Blist[["K119ub Mut"]]) > 0),
    sum(countOverlaps(Alist[["mC Hypo DMRs"]], Blist[["K119ub Ctrl"]]) > 0),
    sum(countOverlaps(Alist[["mC Hypo DMRs"]], Blist[["K119ub Mut"]]) > 0)),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("Hyper", "Hypo"), c("K119ub Ctrl", "K119ub Mut"))
)
fisher_14a <- fisher.test(fisher_14a_mat)
cat(sprintf("  14a (K119ub):  OR = %.2f, p = %.2e\n",
            fisher_14a$estimate, fisher_14a$p.value))

# Test 3: Section 14b - DMR direction x K119ub differential
fisher_14b_mat <- matrix(
  c(sum(countOverlaps(Alist[["mC Hyper DMRs"]], Blist[["K119ub Gained"]]) > 0),
    sum(countOverlaps(Alist[["mC Hyper DMRs"]], Blist[["K119ub Lost"]]) > 0),
    sum(countOverlaps(Alist[["mC Hypo DMRs"]], Blist[["K119ub Gained"]]) > 0),
    sum(countOverlaps(Alist[["mC Hypo DMRs"]], Blist[["K119ub Lost"]]) > 0)),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("Hyper", "Hypo"), c("K119ub Gained", "K119ub Lost"))
)
fisher_14b <- fisher.test(fisher_14b_mat)
cat(sprintf("  14b (K119ub diff): OR = %.2f, p = %.2e\n",
            fisher_14b$estimate, fisher_14b$p.value))

# Test 4: Section 19f - DMR direction x H3K27ac condition
fisher_19f_mat <- matrix(
  c(sum(countOverlaps(Alist[["mC Hyper DMRs"]], Blist[["H3K27ac Ctrl"]]) > 0),
    sum(countOverlaps(Alist[["mC Hyper DMRs"]], Blist[["H3K27ac Mut"]]) > 0),
    sum(countOverlaps(Alist[["mC Hypo DMRs"]], Blist[["H3K27ac Ctrl"]]) > 0),
    sum(countOverlaps(Alist[["mC Hypo DMRs"]], Blist[["H3K27ac Mut"]]) > 0)),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("Hyper", "Hypo"), c("H3K27ac Ctrl", "H3K27ac Mut"))
)
fisher_19f <- fisher.test(fisher_19f_mat)
cat(sprintf("  19f (H3K27ac): OR = %.2f, p = %.2e\n",
            fisher_19f$estimate, fisher_19f$p.value))

# =============================================================================
# BLOCK 10: BUILD 16-ROW COMPARISON TABLE
# =============================================================================

cat("\n--- Building Fisher vs Permutation comparison table ---\n")

# Map each crosswise cell to its parent Fisher test
cell_mapping <- data.frame(
  dmr_set = rep(c("mC Hyper DMRs", "mC Hypo DMRs"), each = 8),
  mark_set = rep(c("ATAC Up", "ATAC Down", "K119ub Ctrl", "K119ub Mut",
                    "K119ub Gained", "K119ub Lost", "H3K27ac Ctrl", "H3K27ac Mut"), 2),
  original_section = rep(c("12", "12", "14", "14", "14", "14", "19", "19"), 2),
  original_figure = rep(c("12a", "12a", "14a", "14a", "14b", "14b", "19f", "19f"), 2),
  stringsAsFactors = FALSE
)

# Attach Fisher results to each cell via parent test
fisher_lookup <- list(
  "12a" = fisher_12a,
  "14a" = fisher_14a,
  "14b" = fisher_14b,
  "19f" = fisher_19f
)

cell_mapping$fisher_or <- sapply(cell_mapping$original_figure, function(fig) {
  fisher_lookup[[fig]]$estimate
})
cell_mapping$fisher_p <- sapply(cell_mapping$original_figure, function(fig) {
  fisher_lookup[[fig]]$p.value
})

# Attach permutation results to each cell
cell_mapping$perm_z <- NA_real_
cell_mapping$perm_p <- NA_real_
cell_mapping$perm_norm_z <- NA_real_
cell_mapping$perm_n_hits <- NA_real_

for (i in seq_len(nrow(cell_mapping))) {
  match_idx <- which(perm_results$dmr_set == cell_mapping$dmr_set[i] &
                     perm_results$name == cell_mapping$mark_set[i])
  if (length(match_idx) == 1) {
    cell_mapping$perm_z[i]       <- perm_results$z_score[match_idx]
    cell_mapping$perm_p[i]       <- perm_results$p_value[match_idx]
    cell_mapping$perm_norm_z[i]  <- perm_results$norm_zscore[match_idx]
    cell_mapping$perm_n_hits[i]  <- perm_results$n_hits[match_idx]
  }
}

# Classify concordance
classify_concordance <- function(fisher_p, perm_p, threshold = 0.05) {
  dplyr::case_when(
    fisher_p < threshold & perm_p < threshold ~ "Confirmed",
    fisher_p < threshold & perm_p >= threshold ~ "Weakened",
    fisher_p >= threshold & perm_p < threshold ~ "Strengthened",
    TRUE ~ "Concordant NS"
  )
}

cell_mapping$concordance <- classify_concordance(cell_mapping$fisher_p, cell_mapping$perm_p)
cell_mapping$test_id <- paste(cell_mapping$dmr_set, "x", cell_mapping$mark_set)

# Save comparison table
write.table(cell_mapping, file.path(TABLES_DIR, "permutation_34_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: permutation_34_comparison.tsv\n")

# Print concordance summary
conc_table <- table(cell_mapping$concordance)
cat("\n  Concordance Summary:\n")
for (cat_name in names(conc_table)) {
  cat(sprintf("    %s: %d/%d\n", cat_name, conc_table[cat_name], nrow(cell_mapping)))
}

# =============================================================================
# BLOCK 11: FIGURE 34b - Fisher vs Permutation Forest Plot
# =============================================================================

cat("\n--- FIGURE 34b: Fisher vs Permutation Comparison ---\n")

# Convert Fisher p to z-score equivalent for visual comparison
cell_mapping$fisher_z_equiv <- qnorm(1 - cell_mapping$fisher_p / 2) *
  sign(log2(cell_mapping$fisher_or))

# Build long-form data for dual-point forest plot
comparison_long <- cell_mapping %>%
  dplyr::select(test_id, mark_set, dmr_set, concordance,
                fisher_z_equiv, perm_z, original_figure) %>%
  tidyr::pivot_longer(
    cols = c(fisher_z_equiv, perm_z),
    names_to = "method",
    values_to = "z_score"
  ) %>%
  dplyr::mutate(
    method = ifelse(method == "fisher_z_equiv", "Fisher", "Permutation"),
    mark_group = dplyr::case_when(
      grepl("ATAC", mark_set) ~ "ATAC-seq",
      grepl("K119ub", mark_set) ~ "H2AK119ub",
      grepl("H3K27ac", mark_set) ~ "H3K27ac"
    )
  )

concordance_colors <- c(
  "Confirmed"     = "#2CA02C",
  "Weakened"      = "#D62728",
  "Strengthened"  = "#FF7F0E",
  "Concordant NS" = "grey60"
)

p_34b <- ggplot(comparison_long,
                aes(x = z_score,
                    y = reorder(test_id, z_score),
                    shape = method,
                    color = concordance)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = "dotted", color = "grey70") +
  geom_point(size = 3, alpha = 0.9) +
  scale_shape_manual(values = c("Fisher" = 16, "Permutation" = 17), name = "Method") +
  scale_color_manual(values = concordance_colors, name = "Concordance") +
  facet_grid(mark_group ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Fisher vs Permutation Test Concordance",
    subtitle = sprintf("Section 34: DMR intervals x chromatin mark peaks (%s permutations)",
                       format(PERM_NTIMES, big.mark = ",")),
    x = "Z-score (Fisher: qnorm equivalent; Permutation: raw z-score)",
    y = ""
  ) +
  theme_biomodal() +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0, face = "bold"),
    axis.text.y = element_text(size = 8)
  )

save_multiformat_ggplot(p_34b, file.path(SECTION_OUTPUT, "34b_fisher_vs_permutation_comparison"),
                        width = 14, height = 10)
cat("  Saved: 34b_fisher_vs_permutation_comparison\n")

# =============================================================================
# BLOCK 12: FIGURE 34c - Local Z-Score Curve
# =============================================================================

cat("\n--- FIGURE 34c: Local Z-Score Curve ---\n")

# Find strongest association
strongest_idx <- which.max(abs(perm_results$z_score))
strongest_a <- perm_results$dmr_set[strongest_idx]
strongest_b <- perm_results$name[strongest_idx]

cat(sprintf("  Strongest association: %s x %s (z = %.2f, p = %.2e)\n",
            strongest_a, strongest_b,
            perm_results$z_score[strongest_idx],
            perm_results$p_value[strongest_idx]))

lz_cache_path <- file.path(TABLES_DIR, "permutation_34_local_zscore.rds")

if (file.exists(lz_cache_path) && !.force_rerun) {
  cat("  Loading cached local z-score results...\n")
  mlz_34 <- readRDS(lz_cache_path)
} else {
  cat("  Running multiLocalZscore (window=50000, step=1000)...\n")

  set.seed(PERM_SEED)
  mlz_34 <- multiLocalZscore(
    A       = Alist[[strongest_a]],
    Blist   = Blist,
    ranFUN  = PERM_RANFUN,
    evFUN   = PERM_EVFUN,
    window  = 50000,
    step    = 1000,
    genome  = genome
  )
  mlz_34 <- makeLZMatrix(mlz_34)

  saveRDS(mlz_34, lz_cache_path)
  cat("  Saved local z-score cache:", lz_cache_path, "\n")
}

p_34c <- plotSingleLZ(mLZ = mlz_34, RS = strongest_b, smoothing = TRUE)

if (is.ggplot(p_34c)) {
  save_multiformat_ggplot(p_34c, file.path(SECTION_OUTPUT, "34c_local_zscore_strongest"),
                          width = 10, height = 6)
} else {
  save_multiformat_base(
    quote(plotSingleLZ(mLZ = mlz_34, RS = strongest_b, smoothing = TRUE)),
    file.path(SECTION_OUTPUT, "34c_local_zscore_strongest"),
    width = 10, height = 6
  )
}

cat("  Saved: 34c_local_zscore_strongest\n")

# =============================================================================
# BLOCK 13: SECTION SUMMARY
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("SECTION 34 SUMMARY\n")
cat("================================================================================\n\n")

cat(sprintf("Permutation parameters:\n"))
cat(sprintf("  ntimes = %d, ranFUN = %s, per.chromosome = %s, seed = %d\n",
            PERM_NTIMES, PERM_RANFUN, PERM_PER_CHR, PERM_SEED))
cat(sprintf("  Alist: %d sets, Blist: %d sets -> %d pairwise tests\n\n",
            length(Alist), length(Blist), length(Alist) * length(Blist)))

cat("Original Fisher tests:\n")
cat(sprintf("  12a (ATAC):        OR = %.2f, p = %.2e\n", fisher_12a$estimate, fisher_12a$p.value))
cat(sprintf("  14a (K119ub):      OR = %.2f, p = %.2e\n", fisher_14a$estimate, fisher_14a$p.value))
cat(sprintf("  14b (K119ub diff): OR = %.2f, p = %.2e\n", fisher_14b$estimate, fisher_14b$p.value))
cat(sprintf("  19f (H3K27ac):     OR = %.2f, p = %.2e\n", fisher_19f$estimate, fisher_19f$p.value))

cat("\nPermutation concordance:\n")
for (cat_name in c("Confirmed", "Weakened", "Strengthened", "Concordant NS")) {
  n <- sum(cell_mapping$concordance == cat_name)
  cat(sprintf("  %-15s %d/%d (%.0f%%)\n", paste0(cat_name, ":"), n,
              nrow(cell_mapping), 100 * n / nrow(cell_mapping)))
}

cat(sprintf("\nStrongest permutation association:\n"))
cat(sprintf("  %s x %s\n", strongest_a, strongest_b))
cat(sprintf("  z = %.2f, p = %.2e\n",
            perm_results$z_score[strongest_idx],
            perm_results$p_value[strongest_idx]))

cat("\nOutputs:\n")
cat("  34a: Crosswise z-score heatmap (2 x 8)\n")
cat("  34b: Fisher vs permutation comparison (forest plot)\n")
cat("  34c: Local z-score curve (strongest association)\n")
cat("  Table: permutation_34_comparison.tsv\n")
cat("  Cache: permutation_34_dmr_marks.rds\n")

cat("\nSection 34 complete.\n")
