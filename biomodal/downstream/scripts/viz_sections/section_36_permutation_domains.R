# biomodal/downstream/scripts/viz_sections/section_36_permutation_domains.R
# Section 36: Domain-Level Permutation Testing (Compartments + Polycomb)
# Standalone script - sources shared config for all dependencies and data
#
# Validates existing Fisher's exact tests from sections 29 and 30 using
# regioneReloaded genomic permutation tests. Addresses reviewer concerns
# about spatial non-independence of megabase-scale compartment domains and
# Polycomb regions by randomizing interval positions while preserving size,
# count, and chromosome assignment.
#
# Two sub-analyses:
#   36A: DMR direction x A/B compartments + shifts (4x4 matrix, validates section 29)
#   36B: DMR direction x Polycomb domains (4x2 matrix, validates section 30)
#
# Figures produced:
#   36a - Crosswise heatmap: DMR direction x compartment (4x4)
#   36b - Crosswise heatmap: DMR direction x Polycomb domains (4x2)
#   36c - Local z-score curve: mC Hyper DMRs at A Compartment (+/- 50kb)
#   36d - Fisher vs permutation comparison dot-plot
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_36_permutation_domains.R

# =============================================================================
# SETUP
# =============================================================================

source("scripts/viz_sections/_shared_config.R")

# Skip if cached results exist (avoids re-run from run_all_sections.sh glob)
# FORCE_RERUN env var (set by run_permutation.sb) bypasses this check
.cache_check <- file.path(TABLES_DIR, "permutation_36_domains.rds")
if (file.exists(.cache_check) && Sys.getenv("FORCE_RERUN", unset = "") == "") {
  cat("Section 36: Cached RDS found at", .cache_check, "-- skipping.\n")
  cat("  To re-run: use run_permutation.sb or export FORCE_RERUN=1\n")
  quit(save = "no", status = 0)
}
.force_rerun <- Sys.getenv("FORCE_RERUN", unset = "") != ""

suppressPackageStartupMessages({
  library(regioneR)
  library(regioneReloaded)
  library(BSgenome.Mmusculus.UCSC.mm10)
})

cat("\n")
cat("===========================================================================\n")
cat("  SECTION 36: DOMAIN-LEVEL PERMUTATION (COMPARTMENTS + POLYCOMB)\n")
cat("===========================================================================\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

SECTION_DIR <- "36_permutation_domains"
SECTION_OUTPUT <- file.path(OUTPUT_DIR, SECTION_DIR)
dir.create(SECTION_OUTPUT, recursive = TRUE, showWarnings = FALSE)

# Permutation parameters (consistent across sections 34-37)
# Accept ntimes from command line (e.g., Rscript section_36_... 5000)
args <- commandArgs(trailingOnly = TRUE)
PERM_NTIMES <- if (length(args) >= 1 && !is.na(as.numeric(args[1]))) as.numeric(args[1]) else 5000

# Use SLURM_CPUS_PER_TASK if available, else default to 8
slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
PERM_CORES <- if (slurm_cpus != "") as.integer(slurm_cpus) else 8L
PERM_PER_CHR  <- TRUE    # Circular permutation within chromosomes
PERM_RANFUN   <- "randomizeRegions"  # Preserves interval sizes exactly
PERM_EVFUN    <- "numOverlaps"       # Default evaluation function
PERM_SEED     <- 42      # Reproducibility

# Cache file path
CACHE_PATH <- file.path(TABLES_DIR, "permutation_36_domains.rds")

# HOMER compartment file (LATE timepoint)
COMPARTMENT_FILE <- file.path(BASE_DIR, "../../tads/tad-pc-analysis/inputs/late/diffPC/diffcompartments.txt")

# Compartment shift thresholds (matching section 29)
SHIFT_FDR  <- 0.05
SHIFT_DIFF <- 0.30

# =============================================================================
# GENOME OBJECT
# =============================================================================

cat("Constructing mm10 genome object (standard chromosomes only)...\n")
genome_full <- getGenomeAndMask("mm10")$genome
standard_chrs <- paste0("chr", c(1:19, "X"))
genome <- genome_full[seqnames(genome_full) %in% standard_chrs]
cat(sprintf("  Genome: %d chromosomes, %s total bp\n",
            length(genome), format(sum(width(genome)), big.mark = ",")))

# =============================================================================
# LOAD REGION SETS
# =============================================================================

cat("\n--- Loading region sets ---\n\n")

# --- DMR GRanges (Alist, shared by 36A and 36B) ---

cat("Constructing DMR GRanges (4 direction sets)...\n")
stopifnot("mc_dmr not loaded from shared config" = exists("mc_dmr") && nrow(mc_dmr) > 0)
stopifnot("hmc_dmr not loaded from shared config" = exists("hmc_dmr") && nrow(hmc_dmr) > 0)

mc_hyper_gr  <- dmr_to_granges(mc_dmr  %>% dplyr::filter(significant, mod_difference > 0))
mc_hypo_gr   <- dmr_to_granges(mc_dmr  %>% dplyr::filter(significant, mod_difference < 0))
hmc_hypo_gr  <- dmr_to_granges(hmc_dmr %>% dplyr::filter(significant, mod_difference < 0))
hmc_hyper_gr <- dmr_to_granges(hmc_dmr %>% dplyr::filter(significant, mod_difference > 0))

cat(sprintf("  mC Hyper DMRs:  %d regions\n", length(mc_hyper_gr)))
cat(sprintf("  mC Hypo DMRs:   %d regions\n", length(mc_hypo_gr)))
cat(sprintf("  hmC Hypo DMRs:  %d regions\n", length(hmc_hypo_gr)))
cat(sprintf("  hmC Hyper DMRs: %d regions\n", length(hmc_hyper_gr)))

stopifnot("No mC Hyper DMRs" = length(mc_hyper_gr) > 0)
stopifnot("No mC Hypo DMRs" = length(mc_hypo_gr) > 0)
stopifnot("No hmC Hypo DMRs" = length(hmc_hypo_gr) > 0)
stopifnot("No hmC Hyper DMRs" = length(hmc_hyper_gr) > 0)

# --- Compartment GRanges (Blist for 36A) ---

cat("\nLoading HOMER compartment data...\n")
stopifnot("Compartment file not found" = file.exists(COMPARTMENT_FILE))

comp_raw <- read.table(COMPARTMENT_FILE, sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE, check.names = FALSE,
                        comment.char = "", quote = "", fill = TRUE)
cat(sprintf("  Loaded %d compartment bins\n", nrow(comp_raw)))

# Identify PC1 columns (same pattern as section 29)
pc1_cols <- grep("bedGraph avg over given bp", names(comp_raw), value = TRUE)
stopifnot("Expected 6 PC1 sample columns" = length(pc1_cols) == 6)
ctrl_pc1_cols <- pc1_cols[1:3]

# Identify difference and adj.pvalue columns
difference_col <- grep("ctrl vs\\. mut Difference",
                        names(comp_raw), value = TRUE, ignore.case = TRUE)[1]
adj_pvalue_col <- grep("ctrl vs\\. mut adj\\. p-value",
                        names(comp_raw), value = TRUE, ignore.case = TRUE)[1]
stopifnot("Difference column not found" = !is.na(difference_col))
stopifnot("Adj. p-value column not found" = !is.na(adj_pvalue_col))

cat(sprintf("  PC1 columns: %d ctrl\n", length(ctrl_pc1_cols)))
cat(sprintf("  Difference column: %s\n", difference_col))
cat(sprintf("  Adj. p-value column: %s\n", adj_pvalue_col))

# Compute mean control PC1
mean_ctrl_pc1 <- rowMeans(
  data.frame(
    as.numeric(comp_raw[[ctrl_pc1_cols[1]]]),
    as.numeric(comp_raw[[ctrl_pc1_cols[2]]]),
    as.numeric(comp_raw[[ctrl_pc1_cols[3]]])
  ),
  na.rm = TRUE
)

comp_difference <- as.numeric(comp_raw[[difference_col]])
comp_adj_pvalue <- as.numeric(comp_raw[[adj_pvalue_col]])

# Construct GRanges for all compartment bins
comp_gr <- GRanges(
  seqnames = comp_raw[["Chr"]],
  ranges = IRanges(start = comp_raw[["Start"]], end = comp_raw[["End"]]),
  mean_ctrl_pc1 = mean_ctrl_pc1,
  difference = comp_difference,
  adj_pvalue = comp_adj_pvalue
)

# Filter to standard chromosomes
comp_gr <- comp_gr[seqnames(comp_gr) %in% standard_chrs]
cat(sprintf("  After filtering to standard chromosomes: %d bins\n", length(comp_gr)))

# Split into A/B compartments and shifts
a_comp       <- comp_gr[comp_gr$mean_ctrl_pc1 > 0]
b_comp       <- comp_gr[comp_gr$mean_ctrl_pc1 <= 0]
b_to_a_shift <- comp_gr[!is.na(comp_gr$adj_pvalue) &
                           comp_gr$adj_pvalue < SHIFT_FDR &
                           comp_gr$difference > SHIFT_DIFF]
a_to_b_shift <- comp_gr[!is.na(comp_gr$adj_pvalue) &
                           comp_gr$adj_pvalue < SHIFT_FDR &
                           comp_gr$difference < -SHIFT_DIFF]

cat(sprintf("  A Compartment:  %d bins\n", length(a_comp)))
cat(sprintf("  B Compartment:  %d bins\n", length(b_comp)))
cat(sprintf("  B->A Shift:     %d bins\n", length(b_to_a_shift)))
cat(sprintf("  A->B Shift:     %d bins\n", length(a_to_b_shift)))

stopifnot("No A compartment bins" = length(a_comp) > 0)
stopifnot("No B compartment bins" = length(b_comp) > 0)

if (length(b_to_a_shift) < 50) {
  cat("  WARNING: B->A shift has < 50 bins — permutation may have low power\n")
}
if (length(a_to_b_shift) < 50) {
  cat("  WARNING: A->B shift has < 50 bins — permutation may have low power\n")
}

# --- Polycomb GRanges (Blist for 36B) ---

cat("\nLoading Polycomb peak sets...\n")
h3k27me3_gr <- load_chip_peaks(CHIP_PEAK_FILES$h3k27me3, "H3K27me3")
bivalent_gr <- load_chip_peaks(CHIP_PEAK_FILES$bivalent, "Bivalent")
stopifnot("H3K27me3 peaks not loaded" = !is.null(h3k27me3_gr))
stopifnot("Bivalent peaks not loaded" = !is.null(bivalent_gr))

# =============================================================================
# BUILD ALIST / BLIST FOR EACH SUB-ANALYSIS
# =============================================================================

cat("\n--- Assembling region set lists ---\n\n")

# Shared Alist (4 DMR direction sets)
Alist_36 <- list(
  "mC Hyper DMRs"  = mc_hyper_gr,
  "mC Hypo DMRs"   = mc_hypo_gr,
  "hmC Hypo DMRs"  = hmc_hypo_gr,
  "hmC Hyper DMRs" = hmc_hyper_gr
)

# 36A: DMR direction x compartments
Blist_36a <- list(
  "A Compartment" = a_comp,
  "B Compartment" = b_comp
)
# Add shift sets only if they have regions
if (length(b_to_a_shift) > 0) {
  Blist_36a[["B->A Shift"]] <- b_to_a_shift
}
if (length(a_to_b_shift) > 0) {
  Blist_36a[["A->B Shift"]] <- a_to_b_shift
}
cat(sprintf("36A: %d x %d = %d pairwise tests\n",
            length(Alist_36), length(Blist_36a),
            length(Alist_36) * length(Blist_36a)))

# 36B: DMR direction x Polycomb domains
Blist_36b <- list(
  "H3K27me3 Peaks" = h3k27me3_gr,
  "Bivalent Peaks"  = bivalent_gr
)
cat(sprintf("36B: %d x %d = %d pairwise tests\n",
            length(Alist_36), length(Blist_36b),
            length(Alist_36) * length(Blist_36b)))

cat(sprintf("\nTotal pairwise permutation tests: %d\n",
            length(Alist_36) * length(Blist_36a) +
            length(Alist_36) * length(Blist_36b)))


# =============================================================================
# PERMUTATION TESTS (with RDS caching)
# =============================================================================

if (file.exists(CACHE_PATH) && !.force_rerun) {
  cat("\n--- Loading cached permutation results ---\n")
  cat(sprintf("  Cache: %s\n", CACHE_PATH))
  perm_36 <- readRDS(CACHE_PATH)
  cw_36a <- perm_36$cw_36a
  cw_36b <- perm_36$cw_36b
  cat("  Loaded 2 crosswise results from cache.\n")
} else {
  cat("\n--- Running permutation tests ---\n")
  cat(sprintf("  Parameters: ntimes=%d, cores=%d, ranFUN=%s, per.chromosome=%s\n",
              PERM_NTIMES, PERM_CORES, PERM_RANFUN, PERM_PER_CHR))

  # 36A: DMR direction x compartments
  cat(sprintf("\n[36A] DMR direction x A/B compartments (4x%d)...\n", length(Blist_36a)))
  set.seed(PERM_SEED)
  cw_36a <- crosswisePermTest(
    Alist          = Alist_36,
    Blist          = Blist_36a,
    genome         = genome,
    ranFUN         = PERM_RANFUN,
    evFUN          = PERM_EVFUN,
    ntimes         = PERM_NTIMES,
    mc.cores       = PERM_CORES,
    per.chromosome = PERM_PER_CHR
  )
  # symm_matrix=FALSE + hc.method="average": non-square design, <=2 Blist
  # elements makes chooseHclustMet crash (cophenetic corr undefined for 1x1 dist)
  cw_36a <- makeCrosswiseMatrix(cw_36a, pvcut = 1,
                                 symm_matrix = FALSE, hc.method = "average")
  cat("  36A complete.\n")

  # 36B: DMR direction x Polycomb domains
  cat("\n[36B] DMR direction x Polycomb domains (4x2)...\n")
  set.seed(PERM_SEED)
  cw_36b <- crosswisePermTest(
    Alist          = Alist_36,
    Blist          = Blist_36b,
    genome         = genome,
    ranFUN         = PERM_RANFUN,
    evFUN          = PERM_EVFUN,
    ntimes         = PERM_NTIMES,
    mc.cores       = PERM_CORES,
    per.chromosome = PERM_PER_CHR
  )
  cw_36b <- makeCrosswiseMatrix(cw_36b, pvcut = 1,
                                 symm_matrix = FALSE, hc.method = "average")
  cat("  36B complete.\n")

  # Save cache
  cat("\nSaving permutation results to cache...\n")
  perm_36 <- list(cw_36a = cw_36a, cw_36b = cw_36b)
  saveRDS(perm_36, CACHE_PATH)
  cat(sprintf("  Saved: %s\n", CACHE_PATH))
}

# =============================================================================
# FIGURE 36a: CROSSWISE HEATMAP — DMR DIRECTION x COMPARTMENT
# =============================================================================

cat("\n--- Generating Figure 36a: DMR x compartment heatmap ---\n")

p_36a <- plotCrosswiseMatrix(cw_36a, matrix_type = "association") +
  ggtitle("Permutation Test: DMR Direction x A/B Compartment",
          subtitle = sprintf("randomizeRegions, %d permutations, per.chromosome=TRUE",
                             PERM_NTIMES)) +
  theme_biomodal()

save_multiformat_ggplot(
  p_36a,
  file.path(SECTION_OUTPUT, "36a_crosswise_dmr_x_compartment"),
  width = 10, height = 8
)

# =============================================================================
# FIGURE 36b: CROSSWISE HEATMAP — DMR DIRECTION x POLYCOMB DOMAINS
# =============================================================================

cat("--- Generating Figure 36b: DMR x Polycomb heatmap ---\n")

p_36b <- plotCrosswiseMatrix(cw_36b, matrix_type = "association") +
  ggtitle("Permutation Test: DMR Direction x Polycomb Domains",
          subtitle = sprintf("randomizeRegions, %d permutations, per.chromosome=TRUE",
                             PERM_NTIMES)) +
  theme_biomodal()

save_multiformat_ggplot(
  p_36b,
  file.path(SECTION_OUTPUT, "36b_crosswise_dmr_x_polycomb"),
  width = 8, height = 8
)

# =============================================================================
# FIGURE 36c: LOCAL Z-SCORE — mC HYPER x A COMPARTMENT
# =============================================================================

cat("--- Generating Figure 36c: Local z-score (mC Hyper x A Compartment) ---\n")
cat("  Running multiLocalZscore (window=50kb, step=1kb)...\n")

set.seed(PERM_SEED)
mlz_36 <- multiLocalZscore(
  A        = mc_hyper_gr,
  Blist    = Blist_36a,
  ranFUN   = PERM_RANFUN,
  evFUN    = PERM_EVFUN,
  genome   = genome,
  window   = 50000,
  step     = 1000,
  mc.cores = PERM_CORES,
  ntimes   = PERM_NTIMES
)
mlz_36 <- makeLZMatrix(mlz_36)

p_36c <- plotSingleLZ(mlz_36, RS = "A Compartment", smoothing = TRUE) +
  ggtitle("Local Z-Score: mC Hyper DMRs at A Compartment",
          subtitle = "+/- 50kb window, 1kb step (broad = domain-level, not focal)") +
  theme_biomodal()

save_multiformat_ggplot(
  p_36c,
  file.path(SECTION_OUTPUT, "36c_local_zscore_compartment"),
  width = 10, height = 6
)

# =============================================================================
# FIGURE 36d: FISHER vs PERMUTATION COMPARISON TABLE
# =============================================================================

cat("--- Generating Figure 36d: Fisher vs permutation comparison ---\n")

# Extract permutation results from crosswise objects
extract_perm_results <- function(cw_obj) {
  mo <- getMultiEvaluation(cw_obj)
  results <- do.call(rbind, lapply(names(mo), function(rs1_name) {
    df <- mo[[rs1_name]]
    data.frame(
      RS1           = rs1_name,
      RS2           = df$name,
      perm_zscore   = df$norm_zscore,
      perm_pvalue   = df$adj_p_value,
      stringsAsFactors = FALSE
    )
  }))
  return(results)
}

perm_results_36a <- extract_perm_results(cw_36a)
perm_results_36b <- extract_perm_results(cw_36b)

# --- Re-compute Fisher tests inline (sections are independent) ---

# 36A: For each DMR pair contrast x each compartment target
# mC contrast: mC Hyper vs mC Hypo
# hmC contrast: hmC Hypo vs hmC Hyper (hmC Hypo is the expected direction)
fisher_36a <- data.frame()
for (target_name in names(Blist_36a)) {
  target_gr <- Blist_36a[[target_name]]

  # mC Hyper vs mC Hypo
  ov_hyper <- sum(countOverlaps(mc_hyper_gr, target_gr) > 0)
  no_hyper <- length(mc_hyper_gr) - ov_hyper
  ov_hypo  <- sum(countOverlaps(mc_hypo_gr, target_gr) > 0)
  no_hypo  <- length(mc_hypo_gr) - ov_hypo

  ft_hyper <- fisher.test(matrix(c(ov_hyper, no_hyper, ov_hypo, no_hypo), nrow = 2))
  fisher_36a <- rbind(fisher_36a, data.frame(
    RS1 = "mC Hyper DMRs", RS2 = target_name,
    fisher_or = as.numeric(ft_hyper$estimate),
    fisher_p = ft_hyper$p.value, stringsAsFactors = FALSE
  ))

  ft_hypo <- fisher.test(matrix(c(ov_hypo, no_hypo, ov_hyper, no_hyper), nrow = 2))
  fisher_36a <- rbind(fisher_36a, data.frame(
    RS1 = "mC Hypo DMRs", RS2 = target_name,
    fisher_or = as.numeric(ft_hypo$estimate),
    fisher_p = ft_hypo$p.value, stringsAsFactors = FALSE
  ))

  # hmC Hypo vs hmC Hyper
  ov_hmc_hypo  <- sum(countOverlaps(hmc_hypo_gr, target_gr) > 0)
  no_hmc_hypo  <- length(hmc_hypo_gr) - ov_hmc_hypo
  ov_hmc_hyper <- sum(countOverlaps(hmc_hyper_gr, target_gr) > 0)
  no_hmc_hyper <- length(hmc_hyper_gr) - ov_hmc_hyper

  ft_hmc_hypo <- fisher.test(matrix(c(ov_hmc_hypo, no_hmc_hypo,
                                        ov_hmc_hyper, no_hmc_hyper), nrow = 2))
  fisher_36a <- rbind(fisher_36a, data.frame(
    RS1 = "hmC Hypo DMRs", RS2 = target_name,
    fisher_or = as.numeric(ft_hmc_hypo$estimate),
    fisher_p = ft_hmc_hypo$p.value, stringsAsFactors = FALSE
  ))

  ft_hmc_hyper <- fisher.test(matrix(c(ov_hmc_hyper, no_hmc_hyper,
                                         ov_hmc_hypo, no_hmc_hypo), nrow = 2))
  fisher_36a <- rbind(fisher_36a, data.frame(
    RS1 = "hmC Hyper DMRs", RS2 = target_name,
    fisher_or = as.numeric(ft_hmc_hyper$estimate),
    fisher_p = ft_hmc_hyper$p.value, stringsAsFactors = FALSE
  ))
}

# 36B: Same pattern for Polycomb targets
fisher_36b <- data.frame()
for (target_name in names(Blist_36b)) {
  target_gr <- Blist_36b[[target_name]]

  # mC Hyper vs mC Hypo
  ov_hyper <- sum(countOverlaps(mc_hyper_gr, target_gr) > 0)
  no_hyper <- length(mc_hyper_gr) - ov_hyper
  ov_hypo  <- sum(countOverlaps(mc_hypo_gr, target_gr) > 0)
  no_hypo  <- length(mc_hypo_gr) - ov_hypo

  ft_hyper <- fisher.test(matrix(c(ov_hyper, no_hyper, ov_hypo, no_hypo), nrow = 2))
  fisher_36b <- rbind(fisher_36b, data.frame(
    RS1 = "mC Hyper DMRs", RS2 = target_name,
    fisher_or = as.numeric(ft_hyper$estimate),
    fisher_p = ft_hyper$p.value, stringsAsFactors = FALSE
  ))

  ft_hypo <- fisher.test(matrix(c(ov_hypo, no_hypo, ov_hyper, no_hyper), nrow = 2))
  fisher_36b <- rbind(fisher_36b, data.frame(
    RS1 = "mC Hypo DMRs", RS2 = target_name,
    fisher_or = as.numeric(ft_hypo$estimate),
    fisher_p = ft_hypo$p.value, stringsAsFactors = FALSE
  ))

  # hmC Hypo vs hmC Hyper
  ov_hmc_hypo  <- sum(countOverlaps(hmc_hypo_gr, target_gr) > 0)
  no_hmc_hypo  <- length(hmc_hypo_gr) - ov_hmc_hypo
  ov_hmc_hyper <- sum(countOverlaps(hmc_hyper_gr, target_gr) > 0)
  no_hmc_hyper <- length(hmc_hyper_gr) - ov_hmc_hyper

  ft_hmc_hypo <- fisher.test(matrix(c(ov_hmc_hypo, no_hmc_hypo,
                                        ov_hmc_hyper, no_hmc_hyper), nrow = 2))
  fisher_36b <- rbind(fisher_36b, data.frame(
    RS1 = "hmC Hypo DMRs", RS2 = target_name,
    fisher_or = as.numeric(ft_hmc_hypo$estimate),
    fisher_p = ft_hmc_hypo$p.value, stringsAsFactors = FALSE
  ))

  ft_hmc_hyper <- fisher.test(matrix(c(ov_hmc_hyper, no_hmc_hyper,
                                         ov_hmc_hypo, no_hmc_hypo), nrow = 2))
  fisher_36b <- rbind(fisher_36b, data.frame(
    RS1 = "hmC Hyper DMRs", RS2 = target_name,
    fisher_or = as.numeric(ft_hmc_hyper$estimate),
    fisher_p = ft_hmc_hyper$p.value, stringsAsFactors = FALSE
  ))
}

# Merge Fisher + permutation results
merge_fisher_perm <- function(fisher_df, perm_df, sub_analysis) {
  merged <- dplyr::inner_join(fisher_df, perm_df, by = c("RS1", "RS2"))
  merged$sub_analysis <- sub_analysis
  merged
}

comparison_36a <- merge_fisher_perm(fisher_36a, perm_results_36a, "36A")
comparison_36b <- merge_fisher_perm(fisher_36b, perm_results_36b, "36B")

comparison_all <- rbind(comparison_36a, comparison_36b)

# Classify concordance
comparison_all$concordance <- dplyr::case_when(
  comparison_all$fisher_p < 0.05 & comparison_all$perm_pvalue < 0.05 ~ "Confirmed",
  comparison_all$fisher_p < 0.05 & comparison_all$perm_pvalue >= 0.05 ~ "Weakened",
  comparison_all$fisher_p >= 0.05 & comparison_all$perm_pvalue < 0.05 ~ "Strengthened",
  TRUE ~ "Concordant NS"
)

# Map to original sections
comparison_all$original_section <- dplyr::case_when(
  comparison_all$sub_analysis == "36A" ~ "29",
  comparison_all$sub_analysis == "36B" ~ "30",
  TRUE ~ "other"
)

# Create test_id column
comparison_all$test_id <- paste(comparison_all$RS1, "x", comparison_all$RS2)

cat(sprintf("  Comparison table: %d tests\n", nrow(comparison_all)))
cat("  Concordance breakdown:\n")
conc_counts <- table(comparison_all$concordance)
for (cat_name in names(conc_counts)) {
  cat(sprintf("    %s: %d\n", cat_name, conc_counts[cat_name]))
}

# Plot comparison dot-plot
concordance_colors <- c(
  "Confirmed"     = "#2ca02c",
  "Weakened"      = "#d62728",
  "Strengthened"  = "#ff7f0e",
  "Concordant NS" = "#7f7f7f"
)

p_36d <- ggplot(comparison_all,
                aes(x = perm_zscore, y = reorder(test_id, perm_zscore),
                    color = concordance)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = "dotted", color = "grey70") +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("OR=%.1f", fisher_or)),
            hjust = -0.2, size = 2.5, show.legend = FALSE) +
  scale_color_manual(values = concordance_colors, name = "Concordance") +
  facet_wrap(~sub_analysis, scales = "free_y", ncol = 1) +
  labs(
    title = "Fisher vs Permutation Test Comparison (Domains)",
    subtitle = "Points = permutation z-score; labels = Fisher OR",
    x = "Permutation Normalized Z-Score",
    y = NULL
  ) +
  theme_biomodal() +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 8)
  )

save_multiformat_ggplot(
  p_36d,
  file.path(SECTION_OUTPUT, "36d_fisher_vs_permutation_table"),
  width = 14, height = max(10, nrow(comparison_all) * 0.4)
)

# =============================================================================
# EXPORT TABLES
# =============================================================================

cat("\n--- Exporting tables ---\n")

comparison_export <- comparison_all[, c("test_id", "sub_analysis", "original_section",
                                         "RS1", "RS2", "fisher_or", "fisher_p",
                                         "perm_zscore", "perm_pvalue", "concordance")]
write.table(
  comparison_export,
  file.path(TABLES_DIR, "permutation_36_comparison.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat(sprintf("  Saved: %s\n", file.path(TABLES_DIR, "permutation_36_comparison.tsv")))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("===========================================================================\n")
cat("  SECTION 36 SUMMARY\n")
cat("===========================================================================\n\n")

cat(sprintf("Sub-analysis 36A (DMR x Compartment):  %d x %d = %d tests\n",
            length(Alist_36), length(Blist_36a),
            length(Alist_36) * length(Blist_36a)))
cat(sprintf("Sub-analysis 36B (DMR x Polycomb):     %d x %d = %d tests\n",
            length(Alist_36), length(Blist_36b),
            length(Alist_36) * length(Blist_36b)))
cat(sprintf("\nTotal permutation tests: %d\n", nrow(comparison_all)))
cat(sprintf("Permutations per test:   %d\n", PERM_NTIMES))
cat(sprintf("Randomization function:  %s\n", PERM_RANFUN))

cat("\nConcordance with Fisher's exact tests:\n")
for (cat_name in c("Confirmed", "Weakened", "Strengthened", "Concordant NS")) {
  n <- sum(comparison_all$concordance == cat_name)
  pct <- 100 * n / nrow(comparison_all)
  cat(sprintf("  %-16s %3d (%5.1f%%)\n", cat_name, n, pct))
}

cat("\n--- Output Locations ---\n")
cat(sprintf("  Figures: %s/\n", SECTION_OUTPUT))
cat(sprintf("  Tables:  %s/permutation_36_*.tsv\n", TABLES_DIR))
cat(sprintf("  Cache:   %s\n", CACHE_PATH))

cat("\n=== Section 36 Complete ===\n\n")
