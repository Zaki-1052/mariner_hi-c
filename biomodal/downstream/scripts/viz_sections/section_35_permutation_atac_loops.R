# biomodal/downstream/scripts/viz_sections/section_35_permutation_atac_loops.R
# Section 35: ATAC x Chromatin Features + Loop Anchor Permutation Testing
# Standalone script - sources shared config for all dependencies and data
#
# Validates existing Fisher's exact tests from sections 13b, 13c, 13d, 27a/b/e,
# and 31a using regioneReloaded genomic permutation tests. Addresses reviewer
# concerns about spatial non-independence of genomic features by randomizing
# interval positions while preserving size, count, and chromosome assignment.
#
# Three sub-analyses:
#   35A: ATAC peaks x 6 ChIP marks (2x6 matrix, validates section 13b)
#   35B: ATAC peaks x 7 chromatin states (2x7 matrix, validates section 13c)
#   35C: Loop anchors x chromatin features (2x6 matrix, validates 13d, 27a/b/e, 31a)
#
# Figures produced:
#   35a - Crosswise heatmap: ATAC direction x 6 ChIP marks
#   35b - Crosswise heatmap: ATAC direction x 7 chromatin states
#   35c - Crosswise heatmap: Loop anchor direction x chromatin features
#   35d - Fisher vs permutation comparison table
#   35e - Local z-score curve for strongest loop anchor association
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_35_permutation_atac_loops.R

# =============================================================================
# SETUP
# =============================================================================

source("scripts/viz_sections/_shared_config.R")

# Skip if cached results exist (avoids re-run from run_all_sections.sh glob)
# FORCE_RERUN env var (set by run_permutation.sb) bypasses this check
.cache_check <- file.path(TABLES_DIR, "permutation_35_atac_loops.rds")
if (file.exists(.cache_check) && Sys.getenv("FORCE_RERUN", unset = "") == "") {
  cat("Section 35: Cached RDS found at", .cache_check, "-- skipping.\n")
  cat("  To re-run: use run_permutation.sb or export FORCE_RERUN=1\n")
  quit(save = "no", status = 0)
}
.force_rerun <- Sys.getenv("FORCE_RERUN", unset = "") != ""

suppressPackageStartupMessages({
  library(regioneR)
  library(regioneReloaded)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(ChIPseeker)
})

cat("\n")
cat("===========================================================================\n")
cat("  SECTION 35: ATAC x CHROMATIN FEATURES + LOOP ANCHOR PERMUTATION\n")
cat("===========================================================================\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

SECTION_DIR <- "35_permutation_atac_loops"
SECTION_OUTPUT <- file.path(OUTPUT_DIR, SECTION_DIR)
dir.create(SECTION_OUTPUT, recursive = TRUE, showWarnings = FALSE)

# Permutation parameters (consistent across sections 34-37)
# Accept ntimes from command line (e.g., Rscript section_35_... 5000)
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
CACHE_PATH <- file.path(TABLES_DIR, "permutation_35_atac_loops.rds")

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

# ATAC differential peaks (Alist for 35A and 35B)
cat("Loading ATAC differential peaks...\n")
atac_up   <- load_chip_peaks(ATAC_FILES$up, "ATAC Up")
atac_down <- load_chip_peaks(ATAC_FILES$down, "ATAC Down")
stopifnot(!is.null(atac_up), !is.null(atac_down))

# 6 ChIP mark peaks (Blist for 35A)
cat("\nLoading ChIP-seq peak sets...\n")
chip_ctcf     <- load_chip_peaks(CHIP_PEAK_FILES$ctcf, "CTCF")
chip_h3k27ac  <- load_chip_peaks(CHIP_PEAK_FILES$h3k27ac, "H3K27ac")
chip_h3k27me3 <- load_chip_peaks(CHIP_PEAK_FILES$h3k27me3, "H3K27me3")
chip_h3k4me1  <- load_chip_peaks(CHIP_PEAK_FILES$h3k4me1, "H3K4me1")
chip_h3k4me3  <- load_chip_peaks(CHIP_PEAK_FILES$h3k4me3, "H3K4me3")
chip_bivalent <- load_chip_peaks(CHIP_PEAK_FILES$bivalent, "Bivalent")
stopifnot(!is.null(chip_ctcf), !is.null(chip_h3k27ac), !is.null(chip_h3k27me3),
          !is.null(chip_h3k4me1), !is.null(chip_h3k4me3), !is.null(chip_bivalent))

# MeCP2 differential peaks (for 35C Blist)
cat("\nLoading MeCP2 differential peaks...\n")
mecp2_up   <- load_chip_peaks(MECP2_FILES$up, "MeCP2 Up")
mecp2_down <- load_chip_peaks(MECP2_FILES$down, "MeCP2 Down")
stopifnot(!is.null(mecp2_up), !is.null(mecp2_down))

# DMR GRanges for 35C Blist
cat("\nConstructing DMR GRanges...\n")
hyper_gr <- dmr_to_granges(mc_dmr %>% dplyr::filter(significant, mod_difference > 0))
hypo_gr  <- dmr_to_granges(mc_dmr %>% dplyr::filter(significant, mod_difference < 0))
cat(sprintf("  mC Hyper DMRs: %d regions\n", length(hyper_gr)))
cat(sprintf("  mC Hypo DMRs:  %d regions\n", length(hypo_gr)))

# Loop anchors (for 35C Alist)
cat("\nLoading Hi-C loop data...\n")
stopifnot(file.exists(LOOP_FILES$late))
loops <- read.table(LOOP_FILES$late, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
required_cols <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "direction")
stopifnot(all(required_cols %in% colnames(loops)))
cat(sprintf("  Loaded %d loops\n", nrow(loops)))

anchor1_gr <- GRanges(seqnames = loops$chr1,
                      ranges = IRanges(start = loops$start1, end = loops$end1))
anchor2_gr <- GRanges(seqnames = loops$chr2,
                      ranges = IRanges(start = loops$start2, end = loops$end2))

gained_mask <- loops$direction == "up_in_mutant"
lost_mask   <- loops$direction == "down_in_mutant"

gained_anchors <- reduce(c(anchor1_gr[gained_mask], anchor2_gr[gained_mask]))
lost_anchors   <- reduce(c(anchor1_gr[lost_mask], anchor2_gr[lost_mask]))
cat(sprintf("  Gained loop anchors: %d unique regions (from %d loops)\n",
            length(gained_anchors), sum(gained_mask)))
cat(sprintf("  Lost loop anchors:   %d unique regions (from %d loops)\n",
            length(lost_anchors), sum(lost_mask)))

# Chromatin state regions (for 35B Blist)
cat("\nConstructing chromatin state regions...\n")
all_dmr_gr <- dmr_to_granges(mc_dmr)
chip_peaks_list <- list(
  ctcf     = chip_ctcf,
  h3k27ac  = chip_h3k27ac,
  h3k27me3 = chip_h3k27me3,
  h3k4me1  = chip_h3k4me1,
  h3k4me3  = chip_h3k4me3,
  bivalent = chip_bivalent
)
chip_overlaps <- compute_chip_overlaps(all_dmr_gr, chip_peaks_list)

# Get distance to TSS via ChIPseeker
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
peak_anno <- annotatePeak(all_dmr_gr, TxDb = txdb, level = "gene",
                          annoDb = "org.Mm.eg.db", verbose = FALSE)
distance_to_tss <- abs(peak_anno@anno$distanceToTSS)

chromatin_states <- classify_chromatin_state(chip_overlaps, distance_to_tss)

# Split DMR GRanges by chromatin state
state_gr_list <- split(all_dmr_gr, chromatin_states)

cat("  Chromatin state region counts:\n")
for (state in CHROMATIN_STATE_ORDER) {
  n <- if (state %in% names(state_gr_list)) length(state_gr_list[[state]]) else 0
  cat(sprintf("    %s: %d regions\n", state, n))
}

# Filter out empty states
state_gr_list <- state_gr_list[sapply(state_gr_list, length) > 0]

# =============================================================================
# BUILD ALIST / BLIST FOR EACH SUB-ANALYSIS
# =============================================================================

cat("\n--- Assembling region set lists ---\n\n")

# 35A: ATAC direction x 6 ChIP marks
Alist_35a <- list("ATAC Up" = atac_up, "ATAC Down" = atac_down)
Blist_35a <- list(
  "CTCF"      = chip_ctcf,
  "H3K27ac"   = chip_h3k27ac,
  "H3K27me3"  = chip_h3k27me3,
  "H3K4me1"   = chip_h3k4me1,
  "H3K4me3"   = chip_h3k4me3,
  "Bivalent"  = chip_bivalent
)
cat(sprintf("35A: %d x %d = %d pairwise tests\n",
            length(Alist_35a), length(Blist_35a),
            length(Alist_35a) * length(Blist_35a)))

# 35B: ATAC direction x 7 chromatin states
Alist_35b <- Alist_35a
Blist_35b <- state_gr_list
cat(sprintf("35B: %d x %d = %d pairwise tests\n",
            length(Alist_35b), length(Blist_35b),
            length(Alist_35b) * length(Blist_35b)))

# 35C: Loop anchor direction x chromatin features
Alist_35c <- list(
  "Gained Loop Anchors" = gained_anchors,
  "Lost Loop Anchors"   = lost_anchors
)
Blist_35c <- list(
  "ATAC Up"       = atac_up,
  "ATAC Down"     = atac_down,
  "MeCP2 Up"      = mecp2_up,
  "MeCP2 Down"    = mecp2_down,
  "mC Hyper DMRs" = hyper_gr,
  "mC Hypo DMRs"  = hypo_gr
)
cat(sprintf("35C: %d x %d = %d pairwise tests\n",
            length(Alist_35c), length(Blist_35c),
            length(Alist_35c) * length(Blist_35c)))

cat(sprintf("\nTotal pairwise permutation tests: %d\n",
            length(Alist_35a) * length(Blist_35a) +
            length(Alist_35b) * length(Blist_35b) +
            length(Alist_35c) * length(Blist_35c)))

# =============================================================================
# PERMUTATION TESTS (with RDS caching)
# =============================================================================

if (file.exists(CACHE_PATH) && !.force_rerun) {
  cat("\n--- Loading cached permutation results ---\n")
  cat(sprintf("  Cache: %s\n", CACHE_PATH))
  perm_35 <- readRDS(CACHE_PATH)
  cw_35a <- perm_35$cw_35a
  cw_35b <- perm_35$cw_35b
  cw_35c <- perm_35$cw_35c
  cat("  Loaded 3 crosswise results from cache.\n")
} else {
  cat("\n--- Running permutation tests ---\n")
  cat(sprintf("  Parameters: ntimes=%d, cores=%d, ranFUN=%s, per.chromosome=%s\n",
              PERM_NTIMES, PERM_CORES, PERM_RANFUN, PERM_PER_CHR))

  # 35A: ATAC x ChIP marks
  cat("\n[35A] ATAC peaks x 6 ChIP marks (2x6)...\n")
  set.seed(PERM_SEED)
  cw_35a <- crosswisePermTest(
    Alist          = Alist_35a,
    Blist          = Blist_35a,
    genome         = genome,
    ranFUN         = PERM_RANFUN,
    evFUN          = PERM_EVFUN,
    ntimes         = PERM_NTIMES,
    mc.cores       = PERM_CORES,
    per.chromosome = PERM_PER_CHR
  )
  # symm_matrix=FALSE + hc.method="average": non-square design, <=2 Alist
  # elements makes chooseHclustMet crash (cophenetic corr undefined for 1x1 dist)
  cw_35a <- makeCrosswiseMatrix(cw_35a, pvcut = 1,
                                 symm_matrix = FALSE, hc.method = "average")
  cat("  35A complete.\n")

  # 35B: ATAC x chromatin states
  cat("\n[35B] ATAC peaks x 7 chromatin states (2x7)...\n")
  set.seed(PERM_SEED)
  cw_35b <- crosswisePermTest(
    Alist          = Alist_35b,
    Blist          = Blist_35b,
    genome         = genome,
    ranFUN         = PERM_RANFUN,
    evFUN          = PERM_EVFUN,
    ntimes         = PERM_NTIMES,
    mc.cores       = PERM_CORES,
    per.chromosome = PERM_PER_CHR
  )
  cw_35b <- makeCrosswiseMatrix(cw_35b, pvcut = 1,
                                 symm_matrix = FALSE, hc.method = "average")
  cat("  35B complete.\n")

  # 35C: Loop anchors x chromatin features
  cat("\n[35C] Loop anchors x chromatin features (2x6)...\n")
  set.seed(PERM_SEED)
  cw_35c <- crosswisePermTest(
    Alist          = Alist_35c,
    Blist          = Blist_35c,
    genome         = genome,
    ranFUN         = PERM_RANFUN,
    evFUN          = PERM_EVFUN,
    ntimes         = PERM_NTIMES,
    mc.cores       = PERM_CORES,
    per.chromosome = PERM_PER_CHR
  )
  cw_35c <- makeCrosswiseMatrix(cw_35c, pvcut = 1,
                                 symm_matrix = FALSE, hc.method = "average")
  cat("  35C complete.\n")

  # Save cache
  cat("\nSaving permutation results to cache...\n")
  perm_35 <- list(cw_35a = cw_35a, cw_35b = cw_35b, cw_35c = cw_35c)
  saveRDS(perm_35, CACHE_PATH)
  cat(sprintf("  Saved: %s\n", CACHE_PATH))
}

# =============================================================================
# FIGURE 35a: CROSSWISE HEATMAP — ATAC x ChIP MARKS (2x6)
# =============================================================================

cat("\n--- Generating Figure 35a: ATAC x ChIP marks heatmap ---\n")

p_35a <- plotCrosswiseMatrix(cw_35a, matrix_type = "association") +
  ggtitle("Permutation Test: ATAC Peaks x ChIP Marks",
          subtitle = sprintf("randomizeRegions, %d permutations, per.chromosome=TRUE",
                             PERM_NTIMES)) +
  theme_biomodal()

save_multiformat_ggplot(
  p_35a,
  file.path(SECTION_OUTPUT, "35a_crosswise_atac_x_chip"),
  width = 10, height = 6
)

# =============================================================================
# FIGURE 35b: CROSSWISE HEATMAP — ATAC x CHROMATIN STATES (2x7)
# =============================================================================

cat("--- Generating Figure 35b: ATAC x chromatin states heatmap ---\n")

p_35b <- plotCrosswiseMatrix(cw_35b, matrix_type = "association") +
  ggtitle("Permutation Test: ATAC Peaks x Chromatin States",
          subtitle = sprintf("randomizeRegions, %d permutations, per.chromosome=TRUE",
                             PERM_NTIMES)) +
  theme_biomodal()

save_multiformat_ggplot(
  p_35b,
  file.path(SECTION_OUTPUT, "35b_crosswise_atac_x_chromatin_state"),
  width = 12, height = 6
)

# =============================================================================
# FIGURE 35c: CROSSWISE HEATMAP — LOOP ANCHORS x FEATURES (2x6)
# =============================================================================

cat("--- Generating Figure 35c: Loop anchors x chromatin features heatmap ---\n")

p_35c <- plotCrosswiseMatrix(cw_35c, matrix_type = "association") +
  ggtitle("Permutation Test: Loop Anchors x Chromatin Features",
          subtitle = sprintf("randomizeRegions, %d permutations, per.chromosome=TRUE",
                             PERM_NTIMES)) +
  theme_biomodal()

save_multiformat_ggplot(
  p_35c,
  file.path(SECTION_OUTPUT, "35c_crosswise_loops_x_features"),
  width = 10, height = 6
)

# =============================================================================
# FIGURE 35d: FISHER vs PERMUTATION COMPARISON TABLE
# =============================================================================

cat("--- Generating Figure 35d: Fisher vs permutation comparison ---\n")

# Extract permutation results from all 3 crosswise objects
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

perm_results_35a <- extract_perm_results(cw_35a)
perm_results_35b <- extract_perm_results(cw_35b)
perm_results_35c <- extract_perm_results(cw_35c)

# --- Re-compute Fisher tests inline (sections are independent) ---

# Helper: build Fisher 2x2 from overlap counts
run_fisher_overlap <- function(query_gr, target_gr, bg_gr) {
  a <- sum(countOverlaps(query_gr, target_gr) > 0)
  b <- length(query_gr) - a
  c <- sum(countOverlaps(bg_gr, target_gr) > 0)
  d <- length(bg_gr) - c
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
  data.frame(
    fisher_or = as.numeric(ft$estimate),
    fisher_p  = ft$p.value,
    stringsAsFactors = FALSE
  )
}

# 35A Fisher tests: ATAC Up/Down x each ChIP mark
# Pattern from section 13b: for each mark, test ATAC-down enrichment vs ATAC-up
fisher_35a <- data.frame()
for (mark_name in names(Blist_35a)) {
  mark_gr <- Blist_35a[[mark_name]]

  # ATAC Up x mark
  ov_up <- sum(countOverlaps(atac_up, mark_gr) > 0)
  no_up <- length(atac_up) - ov_up
  ov_dn <- sum(countOverlaps(atac_down, mark_gr) > 0)
  no_dn <- length(atac_down) - ov_dn

  # Fisher: ATAC Up row
  ft_up <- fisher.test(matrix(c(ov_up, no_up, ov_dn, no_dn), nrow = 2))
  fisher_35a <- rbind(fisher_35a, data.frame(
    RS1 = "ATAC Up", RS2 = mark_name,
    fisher_or = as.numeric(ft_up$estimate),
    fisher_p = ft_up$p.value, stringsAsFactors = FALSE
  ))

  # Fisher: ATAC Down row
  ft_dn <- fisher.test(matrix(c(ov_dn, no_dn, ov_up, no_up), nrow = 2))
  fisher_35a <- rbind(fisher_35a, data.frame(
    RS1 = "ATAC Down", RS2 = mark_name,
    fisher_or = as.numeric(ft_dn$estimate),
    fisher_p = ft_dn$p.value, stringsAsFactors = FALSE
  ))
}

# 35B Fisher tests: ATAC Up/Down x each chromatin state
fisher_35b <- data.frame()
for (state in names(Blist_35b)) {
  state_gr <- Blist_35b[[state]]

  ov_up <- sum(countOverlaps(atac_up, state_gr) > 0)
  no_up <- length(atac_up) - ov_up
  ov_dn <- sum(countOverlaps(atac_down, state_gr) > 0)
  no_dn <- length(atac_down) - ov_dn

  ft_up <- fisher.test(matrix(c(ov_up, no_up, ov_dn, no_dn), nrow = 2))
  fisher_35b <- rbind(fisher_35b, data.frame(
    RS1 = "ATAC Up", RS2 = state,
    fisher_or = as.numeric(ft_up$estimate),
    fisher_p = ft_up$p.value, stringsAsFactors = FALSE
  ))

  ft_dn <- fisher.test(matrix(c(ov_dn, no_dn, ov_up, no_up), nrow = 2))
  fisher_35b <- rbind(fisher_35b, data.frame(
    RS1 = "ATAC Down", RS2 = state,
    fisher_or = as.numeric(ft_dn$estimate),
    fisher_p = ft_dn$p.value, stringsAsFactors = FALSE
  ))
}

# 35C Fisher tests: Loop anchors x features
fisher_35c <- data.frame()
for (feat_name in names(Blist_35c)) {
  feat_gr <- Blist_35c[[feat_name]]

  ov_gained <- sum(countOverlaps(gained_anchors, feat_gr) > 0)
  no_gained <- length(gained_anchors) - ov_gained
  ov_lost   <- sum(countOverlaps(lost_anchors, feat_gr) > 0)
  no_lost   <- length(lost_anchors) - ov_lost

  ft_gained <- fisher.test(matrix(c(ov_gained, no_gained, ov_lost, no_lost), nrow = 2))
  fisher_35c <- rbind(fisher_35c, data.frame(
    RS1 = "Gained Loop Anchors", RS2 = feat_name,
    fisher_or = as.numeric(ft_gained$estimate),
    fisher_p = ft_gained$p.value, stringsAsFactors = FALSE
  ))

  ft_lost <- fisher.test(matrix(c(ov_lost, no_lost, ov_gained, no_gained), nrow = 2))
  fisher_35c <- rbind(fisher_35c, data.frame(
    RS1 = "Lost Loop Anchors", RS2 = feat_name,
    fisher_or = as.numeric(ft_lost$estimate),
    fisher_p = ft_lost$p.value, stringsAsFactors = FALSE
  ))
}

# Merge Fisher + permutation results
merge_fisher_perm <- function(fisher_df, perm_df, sub_analysis) {
  merged <- dplyr::inner_join(fisher_df, perm_df, by = c("RS1", "RS2"))
  merged$sub_analysis <- sub_analysis
  merged
}

comparison_35a <- merge_fisher_perm(fisher_35a, perm_results_35a, "35A")
comparison_35b <- merge_fisher_perm(fisher_35b, perm_results_35b, "35B")
comparison_35c <- merge_fisher_perm(fisher_35c, perm_results_35c, "35C")

comparison_all <- rbind(comparison_35a, comparison_35b, comparison_35c)

# Classify concordance
comparison_all$concordance <- dplyr::case_when(
  comparison_all$fisher_p < 0.05 & comparison_all$perm_pvalue < 0.05 ~ "Confirmed",
  comparison_all$fisher_p < 0.05 & comparison_all$perm_pvalue >= 0.05 ~ "Weakened",
  comparison_all$fisher_p >= 0.05 & comparison_all$perm_pvalue < 0.05 ~ "Strengthened",
  TRUE ~ "Concordant NS"
)

# Add original section mapping
comparison_all$original_section <- dplyr::case_when(
  comparison_all$sub_analysis == "35A" ~ "13b",
  comparison_all$sub_analysis == "35B" ~ "13c",
  comparison_all$sub_analysis == "35C" &
    comparison_all$RS2 %in% c("ATAC Up", "ATAC Down") ~ "13d",
  comparison_all$sub_analysis == "35C" &
    comparison_all$RS2 %in% c("MeCP2 Up", "MeCP2 Down") ~ "31a",
  comparison_all$sub_analysis == "35C" &
    comparison_all$RS2 %in% c("mC Hyper DMRs", "mC Hypo DMRs") ~ "27a/b",
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

# Plot comparison table as formatted figure
concordance_colors <- c(
  "Confirmed"     = "#2ca02c",
  "Weakened"      = "#d62728",
  "Strengthened"  = "#ff7f0e",
  "Concordant NS" = "#7f7f7f"
)

p_35d <- ggplot(comparison_all,
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
    title = "Fisher vs Permutation Test Comparison",
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
  p_35d,
  file.path(SECTION_OUTPUT, "35d_fisher_vs_permutation_table"),
  width = 14, height = max(10, nrow(comparison_all) * 0.35)
)

# =============================================================================
# FIGURE 35e: LOCAL Z-SCORE FOR STRONGEST LOOP ANCHOR ASSOCIATION
# =============================================================================

cat("--- Generating Figure 35e: Local z-score for strongest loop anchor ---\n")

# Identify strongest association from 35C
perm_35c_abs <- perm_results_35c
perm_35c_abs$abs_z <- abs(perm_35c_abs$perm_zscore)
strongest_idx <- which.max(perm_35c_abs$abs_z)
strongest_rs1 <- perm_35c_abs$RS1[strongest_idx]
strongest_rs2 <- perm_35c_abs$RS2[strongest_idx]
cat(sprintf("  Strongest 35C association: %s x %s (z=%.2f)\n",
            strongest_rs1, strongest_rs2,
            perm_35c_abs$perm_zscore[strongest_idx]))

# Select the corresponding Alist element
strongest_a <- Alist_35c[[strongest_rs1]]

# Run multiLocalZscore with +/- 50kb window
cat("  Running multiLocalZscore (window=50kb, step=1kb)...\n")
set.seed(PERM_SEED)
mlz_35 <- multiLocalZscore(
  A        = strongest_a,
  Blist    = Blist_35c,
  ranFUN   = PERM_RANFUN,
  evFUN    = PERM_EVFUN,
  genome   = genome,
  window   = 50000,
  step     = 1000,
  mc.cores = PERM_CORES,
  ntimes   = PERM_NTIMES
)
mlz_35 <- makeLZMatrix(mlz_35)

p_35e <- plotSingleLZ(mlz_35, RS = strongest_rs2, smoothing = TRUE) +
  ggtitle(sprintf("Local Z-Score: %s at %s",
                  strongest_rs1, strongest_rs2),
          subtitle = "+/- 50kb window, 1kb step") +
  theme_biomodal()

save_multiformat_ggplot(
  p_35e,
  file.path(SECTION_OUTPUT, "35e_local_zscore_loop"),
  width = 10, height = 6
)

# =============================================================================
# EXPORT TABLES
# =============================================================================

cat("\n--- Exporting tables ---\n")

# Comparison table
comparison_export <- comparison_all[, c("test_id", "sub_analysis", "original_section",
                                         "RS1", "RS2", "fisher_or", "fisher_p",
                                         "perm_zscore", "perm_pvalue", "concordance")]
write.table(
  comparison_export,
  file.path(TABLES_DIR, "permutation_35_comparison.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat(sprintf("  Saved: %s\n", file.path(TABLES_DIR, "permutation_35_comparison.tsv")))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("===========================================================================\n")
cat("  SECTION 35 SUMMARY\n")
cat("===========================================================================\n\n")

cat(sprintf("Sub-analysis 35A (ATAC x ChIP):       %d x %d = %d tests\n",
            length(Alist_35a), length(Blist_35a),
            length(Alist_35a) * length(Blist_35a)))
cat(sprintf("Sub-analysis 35B (ATAC x States):      %d x %d = %d tests\n",
            length(Alist_35b), length(Blist_35b),
            length(Alist_35b) * length(Blist_35b)))
cat(sprintf("Sub-analysis 35C (Loops x Features):   %d x %d = %d tests\n",
            length(Alist_35c), length(Blist_35c),
            length(Alist_35c) * length(Blist_35c)))
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
cat(sprintf("  Tables:  %s/permutation_35_*.tsv\n", TABLES_DIR))
cat(sprintf("  Cache:   %s\n", CACHE_PATH))

cat("\n=== Section 35 Complete ===\n\n")
