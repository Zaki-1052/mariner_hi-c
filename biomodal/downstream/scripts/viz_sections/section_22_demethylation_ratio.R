# biomodal/downstream/scripts/viz_sections/section_22_demethylation_ratio.R
# Section 22: 5hmC/(5mC+5hmC) Demethylation Efficiency Ratio
# Standalone script - sources shared config for all dependencies and data
#
# Computes the TET conversion efficiency ratio per gene: the fraction of
# modified cytosines that have been oxidized to 5hmC. A decrease in BAP1-KO
# indicates impaired active demethylation via the TET pathway.
#
# Data sources:
#   1. Group-level means from DMR BED files (mc_dmr, hmc_dmr) — ~20,979 genes
#   2. Per-sample regional fractions from feature extraction TSVs — ~21,673 genes
#
# Figures:
#   22a: WT ratio density distribution
#   22b: WT vs KO paired violin+box comparison
#   22c: Delta-ratio histogram (sign-colored)
#   22d: 2-panel scatter (delta-ratio vs mc_diff, vs hmc_diff)
#   22e: Delta-ratio by DMR status (4 groups)
#   22f: Delta-ratio by chromatin state (7 categories)
#   22g: Top 30 genes lollipop (most negative delta-ratio)
#   22h: Per-sample ratio distributions (4 samples)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_22_demethylation_ratio.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 22 CONFIGURATION
# =============================================================================

# Feature extraction file paths (per-sample regional fractions, from _shared_config.R)
MC_EXTRACT_FILE  <- EXTRACT_PATHS$mc_regional_frac
HMC_EXTRACT_FILE <- EXTRACT_PATHS$hmc_regional_frac

# Minimum total methylation threshold (below this, ratio is meaningless)
MIN_TOTAL_METHYLATION <- 0.01

# Colors for delta-ratio sign
DELTA_COLORS <- c("Decreased" = "#B2182B", "Increased" = "#2166AC")

# DMR status colors and order
DMR_STATUS_COLORS <- c(
  "Both mC+hmC DMR" = "#984EA3",
  "mC DMR only"     = "#E41A1C",
  "hmC DMR only"    = "#377EB8",
  "Not DMR"         = "grey70"
)
DMR_STATUS_ORDER <- c("Both mC+hmC DMR", "mC DMR only", "hmC DMR only", "Not DMR")

# Per-sample colors (condition-based, lighter/darker for sex)
SAMPLE_COLORS <- c(
  "ctrl-M" = "#4393C3", "ctrl-F" = "#92C5DE",
  "mut-M"  = "#D6604D", "mut-F"  = "#F4A582"
)

# Helper: format p-value
fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

# Helper: paired Cliff's delta (sign-based effect size)
paired_cliffs_delta <- function(x, y) {
  diffs <- x - y
  diffs <- diffs[is.finite(diffs)]
  n <- length(diffs)
  d <- (sum(diffs > 0) - sum(diffs < 0)) / n
  magnitude <- if (abs(d) < 0.147) "negligible"
               else if (abs(d) < 0.330) "small"
               else if (abs(d) < 0.474) "medium"
               else "large"
  list(delta = d, magnitude = magnitude)
}

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

cat("\n")
cat("========================================================================\n")
cat("SECTION 22: 5hmC/(5mC+5hmC) DEMETHYLATION EFFICIENCY RATIO\n")
cat("========================================================================\n\n")

cat("Validating inputs...\n")
stopifnot("mc_dmr not loaded from shared config" = exists("mc_dmr") && !is.null(mc_dmr))
stopifnot("hmc_dmr not loaded from shared config" = exists("hmc_dmr") && !is.null(hmc_dmr))
stopifnot("mC extract file not found" = file.exists(MC_EXTRACT_FILE))
stopifnot("hmC extract file not found" = file.exists(HMC_EXTRACT_FILE))
cat("  All inputs validated.\n\n")

# =============================================================================
# STEP 1: LOAD PER-SAMPLE FEATURE EXTRACTION DATA
# =============================================================================

cat("--- Loading per-sample feature extraction data ---\n")

mc_extract <- read.table(gzfile(MC_EXTRACT_FILE), header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, na.strings = c("", "NA"))
cat(sprintf("  mC extract: %d genes x %d columns\n", nrow(mc_extract), ncol(mc_extract)))

hmc_extract <- read.table(gzfile(HMC_EXTRACT_FILE), header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, na.strings = c("", "NA"))
cat(sprintf("  hmC extract: %d genes x %d columns\n", nrow(hmc_extract), ncol(hmc_extract)))

# Rename sample columns by position (R converts hyphens to dots in headers)
# Column order: Chromosome, Start, End, ctrl-M, mut-M, ctrl-F, mut-F, Annotation, Name
names(mc_extract)[4:7]  <- c("mc_ctrl_M", "mc_mut_M", "mc_ctrl_F", "mc_mut_F")
names(hmc_extract)[4:7] <- c("hmc_ctrl_M", "hmc_mut_M", "hmc_ctrl_F", "hmc_mut_F")

mc_renamed <- mc_extract %>%
  dplyr::select(gene = Name, mc_ctrl_M, mc_mut_M, mc_ctrl_F, mc_mut_F)

hmc_renamed <- hmc_extract %>%
  dplyr::select(gene = Name, hmc_ctrl_M, hmc_mut_M, hmc_ctrl_F, hmc_mut_F)

# Deduplicate genes (keep first occurrence for genes with multiple entries)
mc_renamed  <- mc_renamed  %>% dplyr::distinct(gene, .keep_all = TRUE)
hmc_renamed <- hmc_renamed %>% dplyr::distinct(gene, .keep_all = TRUE)

# Join mC and hmC per-sample data
per_sample_raw <- dplyr::inner_join(mc_renamed, hmc_renamed, by = "gene")
cat(sprintf("  Joined: %d genes with both mC and hmC data\n", nrow(per_sample_raw)))

# Filter complete cases
n_before <- nrow(per_sample_raw)
per_sample <- per_sample_raw %>% dplyr::filter(complete.cases(.))
n_excluded <- n_before - nrow(per_sample)
cat(sprintf("  Excluded %d genes with NA values (%d retained)\n", n_excluded, nrow(per_sample)))

# =============================================================================
# STEP 2: COMPUTE GROUP-LEVEL RATIOS (from DMR BED)
# =============================================================================

cat("\n--- Computing group-level demethylation ratios ---\n")

# Join ALL genes from mc_dmr and hmc_dmr (not just significant)
# Deduplicate by gene name (keep first occurrence for multi-locus genes)
mc_for_join <- mc_dmr %>%
  dplyr::select(gene, chr, start, end,
                mc_ctrl = mean_mod_group1, mc_mut = mean_mod_group2,
                mc_diff = mod_difference, mc_q = dmr_qvalue,
                mc_sig = significant) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

hmc_for_join <- hmc_dmr %>%
  dplyr::select(gene,
                hmc_ctrl = mean_mod_group1, hmc_mut = mean_mod_group2,
                hmc_diff = mod_difference, hmc_q = dmr_qvalue,
                hmc_sig = significant) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

all_genes <- dplyr::inner_join(mc_for_join, hmc_for_join, by = "gene")
cat(sprintf("  Inner-joined mc_dmr (%d) + hmc_dmr (%d) = %d genes\n",
            nrow(mc_dmr), nrow(hmc_dmr), nrow(all_genes)))

# Compute total methylation and ratios
all_genes <- all_genes %>%
  dplyr::mutate(
    total_ctrl = mc_ctrl + hmc_ctrl,
    total_mut  = mc_mut + hmc_mut,
    ratio_ctrl = hmc_ctrl / total_ctrl,
    ratio_mut  = hmc_mut / total_mut,
    delta_ratio = ratio_mut - ratio_ctrl
  )

# Filter meaningless ratios (both denominators near zero) and non-finite values
n_before_filter <- nrow(all_genes)
ratio_data <- all_genes %>%
  dplyr::filter(!(total_ctrl < MIN_TOTAL_METHYLATION & total_mut < MIN_TOTAL_METHYLATION)) %>%
  dplyr::filter(is.finite(ratio_ctrl) & is.finite(ratio_mut) & is.finite(delta_ratio))
n_filtered <- n_before_filter - nrow(ratio_data)
cat(sprintf("  Filtered %d genes (near-zero denominators or non-finite ratios)\n", n_filtered))
cat(sprintf("  Valid ratio data: %d genes\n", nrow(ratio_data)))

# Classify DMR status
ratio_data <- ratio_data %>%
  dplyr::mutate(
    dmr_status = factor(dplyr::case_when(
      mc_sig & hmc_sig  ~ "Both mC+hmC DMR",
      mc_sig & !hmc_sig ~ "mC DMR only",
      !mc_sig & hmc_sig ~ "hmC DMR only",
      TRUE              ~ "Not DMR"
    ), levels = DMR_STATUS_ORDER)
  )

cat("  DMR status breakdown:\n")
for (s in DMR_STATUS_ORDER) {
  cat(sprintf("    %s: %d genes\n", s, sum(ratio_data$dmr_status == s)))
}

# =============================================================================
# STEP 3: COMPUTE PER-SAMPLE RATIOS
# =============================================================================

cat("\n--- Computing per-sample demethylation ratios ---\n")

per_sample_ratios <- per_sample %>%
  dplyr::mutate(
    ratio_ctrl_M = hmc_ctrl_M / (mc_ctrl_M + hmc_ctrl_M),
    ratio_mut_M  = hmc_mut_M  / (mc_mut_M  + hmc_mut_M),
    ratio_ctrl_F = hmc_ctrl_F / (mc_ctrl_F + hmc_ctrl_F),
    ratio_mut_F  = hmc_mut_F  / (mc_mut_F  + hmc_mut_F),
    ratio_ctrl_avg = (ratio_ctrl_M + ratio_ctrl_F) / 2,
    ratio_mut_avg  = (ratio_mut_M + ratio_mut_F) / 2
  )

# Filter genes where any per-sample denominator is near zero
per_sample_ratios <- per_sample_ratios %>%
  dplyr::filter(
    (mc_ctrl_M + hmc_ctrl_M) >= MIN_TOTAL_METHYLATION,
    (mc_mut_M  + hmc_mut_M)  >= MIN_TOTAL_METHYLATION,
    (mc_ctrl_F + hmc_ctrl_F) >= MIN_TOTAL_METHYLATION,
    (mc_mut_F  + hmc_mut_F)  >= MIN_TOTAL_METHYLATION
  ) %>%
  dplyr::filter(complete.cases(ratio_ctrl_M, ratio_mut_M, ratio_ctrl_F, ratio_mut_F))

cat(sprintf("  Valid per-sample ratios: %d genes\n", nrow(per_sample_ratios)))

sample_medians <- data.frame(
  sample = c("ctrl-M", "ctrl-F", "mut-M", "mut-F"),
  condition = c("Control", "Control", "Mutant", "Mutant"),
  median_ratio = c(
    median(per_sample_ratios$ratio_ctrl_M, na.rm = TRUE),
    median(per_sample_ratios$ratio_ctrl_F, na.rm = TRUE),
    median(per_sample_ratios$ratio_mut_M, na.rm = TRUE),
    median(per_sample_ratios$ratio_mut_F, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

cat("  Per-sample median ratios:\n")
for (i in 1:nrow(sample_medians)) {
  cat(sprintf("    %s: %.4f\n", sample_medians$sample[i], sample_medians$median_ratio[i]))
}

# =============================================================================
# STEP 4: CHROMATIN STATE ANNOTATION
# =============================================================================

cat("\n--- Annotating chromatin states ---\n")

chip_peaks <- list(
  ctcf     = load_chip_peaks(CHIP_PEAK_FILES$ctcf, "CTCF"),
  h3k27ac  = load_chip_peaks(CHIP_PEAK_FILES$h3k27ac, "H3K27ac"),
  h3k27me3 = load_chip_peaks(CHIP_PEAK_FILES$h3k27me3, "H3K27me3"),
  h3k4me1  = load_chip_peaks(CHIP_PEAK_FILES$h3k4me1, "H3K4me1"),
  h3k4me3  = load_chip_peaks(CHIP_PEAK_FILES$h3k4me3, "H3K4me3"),
  bivalent = load_chip_peaks(CHIP_PEAK_FILES$bivalent, "Bivalent")
)
stopifnot("Missing ChIP peak files" = all(!sapply(chip_peaks, is.null)))

# Convert to GRanges for overlap analysis
ratio_gr <- GRanges(
  seqnames = ratio_data$chr,
  ranges = IRanges(start = ratio_data$start, end = ratio_data$end)
)

# Compute TSS distances
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
gene_gr <- genes(txdb)
tss_gr <- resize(gene_gr, width = 1, fix = "start")
nearest_hits <- distanceToNearest(ratio_gr, tss_gr)
distance_to_tss <- rep(NA_real_, length(ratio_gr))
distance_to_tss[queryHits(nearest_hits)] <- mcols(nearest_hits)$distance

cat(sprintf("  TSS distance: %d genes with distance=0 (contain TSS)\n",
            sum(distance_to_tss == 0, na.rm = TRUE)))

# Compute overlaps and classify
chip_overlaps <- compute_chip_overlaps(ratio_gr, chip_peaks)
ratio_data$chromatin_state <- classify_chromatin_state(chip_overlaps, distance_to_tss, TSS_THRESHOLD)

cat("  Chromatin state breakdown:\n")
state_counts <- table(ratio_data$chromatin_state)
for (s in CHROMATIN_STATE_ORDER) {
  if (s %in% names(state_counts)) {
    cat(sprintf("    %s: %d genes\n", s, state_counts[s]))
  }
}

# =============================================================================
# STEP 5: CORE STATISTICS
# =============================================================================

cat("\n--- Core statistical tests ---\n")

# Paired Wilcoxon signed-rank: ratio_ctrl vs ratio_mut
wilcox_paired <- wilcox.test(ratio_data$ratio_ctrl, ratio_data$ratio_mut,
                             paired = TRUE, alternative = "two.sided")
cat(sprintf("  Wilcoxon signed-rank (paired): %s\n", fmt_p(wilcox_paired$p.value)))

# One-sample Wilcoxon: delta_ratio vs 0
wilcox_delta <- wilcox.test(ratio_data$delta_ratio, mu = 0, alternative = "two.sided")
cat(sprintf("  Wilcoxon (delta vs 0): %s\n", fmt_p(wilcox_delta$p.value)))

# Cliff's delta effect size
cd <- paired_cliffs_delta(ratio_data$ratio_ctrl, ratio_data$ratio_mut)
cat(sprintf("  Cliff's delta: %.4f (%s)\n", cd$delta, cd$magnitude))

# Proportion with decreased ratio
pct_decreased <- mean(ratio_data$delta_ratio < 0, na.rm = TRUE) * 100
pct_increased <- mean(ratio_data$delta_ratio > 0, na.rm = TRUE) * 100
cat(sprintf("  Genes with decreased ratio (impaired TET): %.1f%%\n", pct_decreased))
cat(sprintf("  Genes with increased ratio: %.1f%%\n", pct_increased))

# Spearman correlations
rho_mc  <- cor.test(ratio_data$delta_ratio, ratio_data$mc_diff, method = "spearman")
rho_hmc <- cor.test(ratio_data$delta_ratio, ratio_data$hmc_diff, method = "spearman")
cat(sprintf("  Spearman rho (delta-ratio vs mc_diff): %.3f (%s)\n",
            rho_mc$estimate, fmt_p(rho_mc$p.value)))
cat(sprintf("  Spearman rho (delta-ratio vs hmc_diff): %.3f (%s)\n",
            rho_hmc$estimate, fmt_p(rho_hmc$p.value)))

# =============================================================================
# FIGURE 22a: WT RATIO DENSITY DISTRIBUTION
# =============================================================================

cat("\n--- Figure 22a: WT demethylation ratio density ---\n")

ratio_median_wt <- median(ratio_data$ratio_ctrl, na.rm = TRUE)
ratio_mean_wt   <- mean(ratio_data$ratio_ctrl, na.rm = TRUE)
ratio_iqr_wt    <- IQR(ratio_data$ratio_ctrl, na.rm = TRUE)
ratio_q25_wt    <- quantile(ratio_data$ratio_ctrl, 0.25, na.rm = TRUE)
ratio_q75_wt    <- quantile(ratio_data$ratio_ctrl, 0.75, na.rm = TRUE)

p_22a <- ggplot(ratio_data, aes(x = ratio_ctrl)) +
  geom_density(fill = COLORS$condition["Control"], alpha = 0.5,
               color = "black", linewidth = 0.5) +
  geom_vline(xintercept = ratio_median_wt, linetype = "dashed",
             color = "black", linewidth = 0.8) +
  geom_vline(xintercept = ratio_mean_wt, linetype = "dotted",
             color = "grey40", linewidth = 0.6) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3.5,
           label = sprintf("N = %s genes\nMedian = %.4f\nMean = %.4f\nIQR = %.4f [%.4f\u2013%.4f]",
                           format(nrow(ratio_data), big.mark = ","),
                           ratio_median_wt, ratio_mean_wt, ratio_iqr_wt,
                           ratio_q25_wt, ratio_q75_wt)) +
  labs(
    title = "Wildtype 5hmC/(5mC+5hmC) Ratio Distribution",
    subtitle = "TET conversion efficiency across all genes (gene body CG context)",
    x = "5hmC / (5mC + 5hmC) ratio",
    y = "Density"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_22a, file.path(OUTPUT_DIR, "22a_wt_ratio_density"), 10, 7)

# =============================================================================
# FIGURE 22b: WT VS KO PAIRED VIOLIN+BOX
# =============================================================================

cat("--- Figure 22b: WT vs KO ratio comparison ---\n")

ratio_long <- ratio_data %>%
  dplyr::select(gene, ratio_ctrl, ratio_mut) %>%
  tidyr::pivot_longer(cols = c(ratio_ctrl, ratio_mut),
                      names_to = "condition", values_to = "ratio") %>%
  dplyr::mutate(
    condition = factor(
      ifelse(condition == "ratio_ctrl", "Control", "Mutant"),
      levels = c("Control", "Mutant")
    )
  )

p_22b <- ggplot(ratio_long, aes(x = condition, y = ratio, fill = condition)) +
  geom_violin(alpha = 0.6, width = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.8) +
  scale_fill_manual(values = COLORS$condition) +
  annotate("text", x = 1.5, y = max(ratio_long$ratio, na.rm = TRUE) * 0.98,
           label = sprintf("Wilcoxon signed-rank\n%s\nCliff's delta = %.3f (%s)",
                           fmt_p(wilcox_paired$p.value), cd$delta, cd$magnitude),
           size = 3.5, hjust = 0.5) +
  labs(
    title = "Demethylation Efficiency: Wildtype vs BAP1-KO",
    subtitle = "5hmC/(5mC+5hmC) ratio per gene (paired comparison)",
    x = "Condition",
    y = "5hmC / (5mC + 5hmC) ratio"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_22b, file.path(OUTPUT_DIR, "22b_wt_vs_ko_ratio"), 8, 8)

# =============================================================================
# FIGURE 22c: DELTA-RATIO HISTOGRAM
# =============================================================================

cat("--- Figure 22c: Delta-ratio distribution ---\n")

ratio_data <- ratio_data %>%
  dplyr::mutate(delta_sign = ifelse(delta_ratio < 0, "Decreased", "Increased"))

p_22c <- ggplot(ratio_data, aes(x = delta_ratio, fill = delta_sign)) +
  geom_histogram(bins = 100, color = "black", linewidth = 0.2, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = median(ratio_data$delta_ratio, na.rm = TRUE),
             linetype = "dotted", color = "grey30", linewidth = 0.6) +
  scale_fill_manual(values = DELTA_COLORS, name = "Direction") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 3.5,
           label = sprintf("%.1f%% decreased (impaired TET)\n%s\nMedian = %.4f",
                           pct_decreased, fmt_p(wilcox_delta$p.value),
                           median(ratio_data$delta_ratio, na.rm = TRUE))) +
  labs(
    title = "Distribution of Delta Demethylation Ratio (KO - WT)",
    subtitle = "Negative values indicate impaired TET conversion in BAP1-KO",
    x = expression(Delta ~ "ratio (Mutant - Control)"),
    y = "Number of genes"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_22c, file.path(OUTPUT_DIR, "22c_delta_ratio_histogram"), 10, 7)

# =============================================================================
# FIGURE 22d: 2-PANEL SCATTER (delta-ratio vs mc_diff, hmc_diff)
# =============================================================================

cat("--- Figure 22d: Delta-ratio vs methylation changes ---\n")

ratio_data <- ratio_data %>%
  dplyr::mutate(is_key_gene = gene %in% KEY_GENES)

# Left panel: delta-ratio vs mc_diff
p_22d_left <- ggplot(ratio_data, aes(x = mc_diff, y = delta_ratio)) +
  geom_point(alpha = 0.1, size = 0.5, color = "grey50") +
  geom_smooth(method = "lm", color = COLORS$methylation["5mC"],
              linewidth = 1, se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_point(data = dplyr::filter(ratio_data, is_key_gene),
             color = COLORS$methylation["5mC"], size = 2) +
  geom_text_repel(data = dplyr::filter(ratio_data, is_key_gene),
                  aes(label = gene), size = 3, max.overlaps = 20,
                  color = COLORS$methylation["5mC"], fontface = "bold") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 3.5,
           label = sprintf("rho = %.3f\n%s",
                           rho_mc$estimate, fmt_p(rho_mc$p.value))) +
  labs(title = "vs 5mC Difference",
       x = "5mC difference (KO - WT)",
       y = expression(Delta ~ "ratio (KO - WT)")) +
  theme_biomodal()

# Right panel: delta-ratio vs hmc_diff
p_22d_right <- ggplot(ratio_data, aes(x = hmc_diff, y = delta_ratio)) +
  geom_point(alpha = 0.1, size = 0.5, color = "grey50") +
  geom_smooth(method = "lm", color = COLORS$methylation["5hmC"],
              linewidth = 1, se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_point(data = dplyr::filter(ratio_data, is_key_gene),
             color = COLORS$methylation["5hmC"], size = 2) +
  geom_text_repel(data = dplyr::filter(ratio_data, is_key_gene),
                  aes(label = gene), size = 3, max.overlaps = 20,
                  color = COLORS$methylation["5hmC"], fontface = "bold") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 3.5,
           label = sprintf("rho = %.3f\n%s",
                           rho_hmc$estimate, fmt_p(rho_hmc$p.value))) +
  labs(title = "vs 5hmC Difference",
       x = "5hmC difference (KO - WT)",
       y = expression(Delta ~ "ratio (KO - WT)")) +
  theme_biomodal()

p_22d <- (p_22d_left | p_22d_right) +
  plot_annotation(
    title = "Delta Demethylation Ratio vs Individual Methylation Changes",
    subtitle = "Spearman correlation with key genes highlighted",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  )

save_multiformat_ggplot(p_22d, file.path(OUTPUT_DIR, "22d_delta_ratio_vs_diff"), 14, 7)

# =============================================================================
# FIGURE 22e: DELTA-RATIO BY DMR STATUS
# =============================================================================

cat("--- Figure 22e: Delta-ratio by DMR status ---\n")

# Kruskal-Wallis across DMR status groups
kw_dmr <- kruskal.test(delta_ratio ~ dmr_status, data = ratio_data)
cat(sprintf("  Kruskal-Wallis: %s\n", fmt_p(kw_dmr$p.value)))

# Pairwise Wilcoxon with BH correction
pw_dmr <- pairwise.wilcox.test(ratio_data$delta_ratio, ratio_data$dmr_status,
                               p.adjust.method = "BH")

# Per-group medians
dmr_medians <- ratio_data %>%
  dplyr::group_by(dmr_status) %>%
  dplyr::summarise(
    n = dplyr::n(),
    median_delta = median(delta_ratio, na.rm = TRUE),
    .groups = "drop"
  )

p_22e <- ggplot(ratio_data, aes(x = dmr_status, y = delta_ratio, fill = dmr_status)) +
  geom_violin(alpha = 0.6, width = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  scale_fill_manual(values = DMR_STATUS_COLORS) +
  geom_text(data = dmr_medians,
            aes(x = dmr_status, y = -Inf,
                label = sprintf("n=%s\nmed=%.4f", format(n, big.mark = ","), median_delta)),
            vjust = -0.3, size = 3) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.3, size = 3.5,
           label = sprintf("Kruskal-Wallis %s", fmt_p(kw_dmr$p.value))) +
  labs(
    title = "Delta Demethylation Ratio by DMR Status",
    subtitle = "Genes classified by differential methylation significance (q < 0.05)",
    x = "DMR Status",
    y = expression(Delta ~ "ratio (Mutant - Control)")
  ) +
  theme_biomodal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

save_multiformat_ggplot(p_22e, file.path(OUTPUT_DIR, "22e_delta_ratio_by_dmr_status"), 10, 8)

# =============================================================================
# FIGURE 22f: DELTA-RATIO BY CHROMATIN STATE
# =============================================================================

cat("--- Figure 22f: Delta-ratio by chromatin state ---\n")

# Per-state Wilcoxon vs 0 with BH correction
state_tests <- ratio_data %>%
  dplyr::group_by(chromatin_state) %>%
  dplyr::summarise(
    n = dplyr::n(),
    median_delta = median(delta_ratio, na.rm = TRUE),
    iqr_delta = IQR(delta_ratio, na.rm = TRUE),
    p_value = tryCatch(
      wilcox.test(delta_ratio, mu = 0)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    sig_label = dplyr::case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

cat("  Per-state Wilcoxon vs 0 (BH-corrected):\n")
for (i in 1:nrow(state_tests)) {
  cat(sprintf("    %s: median=%.4f, n=%d, padj=%s\n",
              state_tests$chromatin_state[i], state_tests$median_delta[i],
              state_tests$n[i], fmt_p(state_tests$p_adj[i])))
}

# Compute y-position for significance labels
y_max_f <- max(ratio_data$delta_ratio, na.rm = TRUE)

p_22f <- ggplot(ratio_data, aes(x = chromatin_state, y = delta_ratio,
                                fill = chromatin_state)) +
  geom_violin(alpha = 0.6, width = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  scale_fill_manual(values = CHROMATIN_STATE_COLORS) +
  geom_text(data = state_tests,
            aes(x = chromatin_state, y = y_max_f * 1.05, label = sig_label),
            size = 5, vjust = 0) +
  geom_text(data = state_tests,
            aes(x = chromatin_state, y = -Inf,
                label = sprintf("n=%s", format(n, big.mark = ","))),
            vjust = -0.3, size = 2.8) +
  labs(
    title = "Delta Demethylation Ratio by Chromatin State",
    subtitle = "Per-state Wilcoxon vs 0 (BH-corrected): *** p<0.001, ** p<0.01, * p<0.05",
    x = "Chromatin State",
    y = expression(Delta ~ "ratio (Mutant - Control)")
  ) +
  theme_biomodal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1))

save_multiformat_ggplot(p_22f, file.path(OUTPUT_DIR, "22f_delta_ratio_by_chromatin_state"), 12, 9)

# =============================================================================
# FIGURE 22g: TOP 30 GENES LOLLIPOP
# =============================================================================

cat("--- Figure 22g: Top 30 genes by most negative delta-ratio ---\n")

top30 <- ratio_data %>%
  dplyr::arrange(delta_ratio) %>%
  dplyr::slice_head(n = 30) %>%
  dplyr::mutate(
    is_key = gene %in% KEY_GENES,
    gene_label = factor(gene, levels = rev(gene))
  )

# Build face/color vectors aligned to factor level order (bottom to top)
level_order <- levels(top30$gene_label)
is_key_ordered <- top30$is_key[match(level_order, top30$gene)]
y_face  <- ifelse(is_key_ordered, "bold", "plain")
y_color <- ifelse(is_key_ordered, "red", "black")

p_22g <- ggplot(top30, aes(x = delta_ratio, y = gene_label, color = chromatin_state)) +
  geom_segment(aes(x = 0, xend = delta_ratio, y = gene_label, yend = gene_label),
               linewidth = 0.8) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = CHROMATIN_STATE_COLORS, name = "Chromatin State") +
  labs(
    title = "Top 30 Genes with Most Decreased Demethylation Ratio",
    subtitle = "Largest decrease in 5hmC/(5mC+5hmC) in BAP1-KO vs wildtype",
    x = expression(Delta ~ "ratio (Mutant - Control)"),
    y = NULL
  ) +
  theme_biomodal() +
  theme(axis.text.y = element_text(face = y_face, color = y_color, size = 9))

save_multiformat_ggplot(p_22g, file.path(OUTPUT_DIR, "22g_top30_decreased_ratio"), 10, 10)

# =============================================================================
# FIGURE 22h: PER-SAMPLE RATIO DISTRIBUTIONS
# =============================================================================

cat("--- Figure 22h: Per-sample ratio distributions ---\n")

per_sample_long <- per_sample_ratios %>%
  dplyr::select(gene, ratio_ctrl_M, ratio_ctrl_F, ratio_mut_M, ratio_mut_F) %>%
  tidyr::pivot_longer(
    cols = starts_with("ratio_"),
    names_to = "sample",
    values_to = "ratio",
    names_prefix = "ratio_"
  ) %>%
  dplyr::mutate(
    sample = factor(sample,
                    levels = c("ctrl_M", "ctrl_F", "mut_M", "mut_F"),
                    labels = c("ctrl-M", "ctrl-F", "mut-M", "mut-F")),
    condition = ifelse(grepl("ctrl", sample), "Control", "Mutant")
  )

p_22h <- ggplot(per_sample_long, aes(x = sample, y = ratio, fill = sample)) +
  geom_violin(alpha = 0.6, width = 0.8) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.8) +
  scale_fill_manual(values = SAMPLE_COLORS) +
  geom_text(data = sample_medians,
            aes(x = sample, y = -Inf, label = sprintf("med=%.4f", median_ratio)),
            vjust = -0.3, size = 3, inherit.aes = FALSE) +
  labs(
    title = "Per-Sample Demethylation Ratio Distributions",
    subtitle = "5hmC/(5mC+5hmC) per gene for each biological replicate",
    x = "Sample",
    y = "5hmC / (5mC + 5hmC) ratio"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_22h, file.path(OUTPUT_DIR, "22h_per_sample_ratio"), 10, 8)

# =============================================================================
# EXPORT TABLES
# =============================================================================

cat("\n--- Exporting tables ---\n")

# Table 1: Full per-gene table (group-level)
export_all <- ratio_data %>%
  dplyr::select(gene, chr, start, end,
                mc_ctrl, mc_mut, mc_diff, mc_q, mc_sig,
                hmc_ctrl, hmc_mut, hmc_diff, hmc_q, hmc_sig,
                total_ctrl, total_mut, ratio_ctrl, ratio_mut, delta_ratio,
                dmr_status, chromatin_state) %>%
  dplyr::arrange(delta_ratio)

write.table(export_all,
            file.path(TABLES_DIR, "demethylation_ratio_all_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved demethylation_ratio_all_genes.tsv (%d rows)\n", nrow(export_all)))

# Table 2: Per-sample ratios
export_per_sample <- per_sample_ratios %>%
  dplyr::select(gene,
                mc_ctrl_M, mc_ctrl_F, mc_mut_M, mc_mut_F,
                hmc_ctrl_M, hmc_ctrl_F, hmc_mut_M, hmc_mut_F,
                ratio_ctrl_M, ratio_ctrl_F, ratio_mut_M, ratio_mut_F,
                ratio_ctrl_avg, ratio_mut_avg) %>%
  dplyr::mutate(delta_ratio_avg = ratio_mut_avg - ratio_ctrl_avg) %>%
  dplyr::arrange(delta_ratio_avg)

write.table(export_per_sample,
            file.path(TABLES_DIR, "demethylation_ratio_per_sample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved demethylation_ratio_per_sample.tsv (%d rows)\n", nrow(export_per_sample)))

# Table 3: Top 50 most decreased
export_top <- ratio_data %>%
  dplyr::arrange(delta_ratio) %>%
  dplyr::slice_head(n = 50) %>%
  dplyr::select(gene, chr, start, end,
                mc_ctrl, mc_mut, mc_diff, mc_q, mc_sig,
                hmc_ctrl, hmc_mut, hmc_diff, hmc_q, hmc_sig,
                ratio_ctrl, ratio_mut, delta_ratio,
                dmr_status, chromatin_state)

write.table(export_top,
            file.path(TABLES_DIR, "demethylation_ratio_top_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved demethylation_ratio_top_genes.tsv (%d rows)\n", nrow(export_top)))

# Table 4: Per-chromatin-state summary
write.table(state_tests,
            file.path(TABLES_DIR, "demethylation_ratio_by_chromatin_state.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved demethylation_ratio_by_chromatin_state.tsv (%d rows)\n", nrow(state_tests)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("========================================================================\n")
cat("SECTION 22 SUMMARY: DEMETHYLATION EFFICIENCY RATIO\n")
cat("========================================================================\n\n")

cat(sprintf("Total genes in mc_dmr: %d\n", nrow(mc_dmr)))
cat(sprintf("Total genes in hmc_dmr: %d\n", nrow(hmc_dmr)))
cat(sprintf("Valid ratio data: %d genes\n", nrow(ratio_data)))

cat(sprintf("\nGroup-level ratio statistics:\n"))
cat(sprintf("  WT median ratio:  %.4f\n", median(ratio_data$ratio_ctrl, na.rm = TRUE)))
cat(sprintf("  KO median ratio:  %.4f\n", median(ratio_data$ratio_mut, na.rm = TRUE)))
cat(sprintf("  Delta median:     %.4f\n", median(ratio_data$delta_ratio, na.rm = TRUE)))
cat(sprintf("  Paired Wilcoxon:  %s\n", fmt_p(wilcox_paired$p.value)))
cat(sprintf("  Cliff's delta:    %.4f (%s)\n", cd$delta, cd$magnitude))
cat(sprintf("  %% decreased:      %.1f%%\n", pct_decreased))

cat(sprintf("\nCorrelations:\n"))
cat(sprintf("  delta-ratio vs mc_diff:  rho=%.3f (%s)\n",
            rho_mc$estimate, fmt_p(rho_mc$p.value)))
cat(sprintf("  delta-ratio vs hmc_diff: rho=%.3f (%s)\n",
            rho_hmc$estimate, fmt_p(rho_hmc$p.value)))

cat(sprintf("\nTop 5 genes (most decreased ratio):\n"))
top5 <- ratio_data %>% dplyr::arrange(delta_ratio) %>% dplyr::slice_head(n = 5)
for (i in 1:5) {
  cat(sprintf("  %d. %s: delta=%.4f (WT=%.4f, KO=%.4f) [%s]\n",
              i, top5$gene[i], top5$delta_ratio[i],
              top5$ratio_ctrl[i], top5$ratio_mut[i],
              as.character(top5$chromatin_state[i])))
}

cat(sprintf("\nPer-chromatin-state median delta-ratios:\n"))
for (i in 1:nrow(state_tests)) {
  cat(sprintf("  %s: %.4f (n=%d, %s)\n",
              state_tests$chromatin_state[i], state_tests$median_delta[i],
              state_tests$n[i], state_tests$sig_label[i]))
}

cat(sprintf("\nLiterature comparison (TET triple-KO):\n"))
cat("  Published TET1/2/3 triple-KO: near-complete 5hmC loss (ratio -> 0)\n")
cat("  BAP1-KO effect is indirect (via PRC1/H2AK119ub -> TET recruitment)\n")
cat(sprintf("  Observed WT median ratio: %.4f\n",
            median(ratio_data$ratio_ctrl, na.rm = TRUE)))
cat(sprintf("  Observed KO median ratio: %.4f\n",
            median(ratio_data$ratio_mut, na.rm = TRUE)))
cat("  Partial reduction consistent with indirect TET impairment, not ablation\n")

cat(sprintf("\nPer-sample analysis (%d genes):\n", nrow(per_sample_ratios)))
for (i in 1:nrow(sample_medians)) {
  cat(sprintf("  %s: median ratio = %.4f\n",
              sample_medians$sample[i], sample_medians$median_ratio[i]))
}

cat(sprintf("\nDelta-ratio column exported in demethylation_ratio_all_genes.tsv\n"))
cat("for reuse as response variable by other sections (Sections 1, 3, 8, 11).\n")

cat("\n--- Output files ---\n")
cat(sprintf("  Figures: %s/22{a-h}_*/\n", OUTPUT_DIR))
cat(sprintf("  Tables:  %s/demethylation_ratio_*.tsv (4 files)\n", TABLES_DIR))

cat("\n=== Section 22 complete ===\n")
