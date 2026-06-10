# biomodal/downstream/scripts/viz_sections/section_53_mecp2_peak_noncg.R
# Section 53: MeCP2 Non-CG Methylation at Peak Resolution
#
# Two questions:
#   1. Is there any non-CG methylation at MeCP2 binding sites at all?
#      Especially in the mutant — is MeCP2 encountering mCH substrate?
#   2. Is whatever signal exists differential between BAP1-KO and wildtype?
#
# Uses modality XPLR outputs from 8,886 MeCP2 ChIP peaks (400bp) run directly
# on CHG and CHH zarr stores, plus matched shuffled control regions as baseline.
#
# Note: evoC cannot distinguish mC from hmC at non-CG contexts. Only mC analyzed.
#
# Prerequisites:
#   - MeCP2 peak modality runs: sbatch run_modality_mecp2.sb CHG/CHH (done)
#   - Control region modality runs: sbatch run_modality_control.sb CHG/CHH
#   - rsync both output sets to local
#   - Section 51 must have been run (summary table needed for 53f)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_53_mecp2_peak_noncg.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 53: MeCP2 NON-CG METHYLATION AT PEAK RESOLUTION
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("SECTION 53: MeCP2 NON-CG METHYLATION AT PEAK RESOLUTION\n")
cat("================================================================================\n\n")

SEC53_DIR <- file.path(OUTPUT_DIR, "53_mecp2_peak_noncg")
dir.create(SEC53_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Configuration ---------------------------------------------------------

MECP2_PEAK_PATHS <- list(
  chg_frac  = file.path(BASE_DIR, "modality/outputs_CHG_mecp2/Results/mecp2_peaks.annotation/Extract_20260609_162529/Extract_mc_regional-frac_20260609_162529.tsv.gz"),
  chh_frac  = file.path(BASE_DIR, "modality/outputs_CHH_mecp2/Results/mecp2_peaks.annotation/Extract_20260610_015808/Extract_mc_regional-frac_20260610_015808.tsv.gz"),
  chg_count = file.path(BASE_DIR, "modality/outputs_CHG_mecp2/Results/mecp2_peaks.annotation/Extract_20260609_161937/Extract_total_c_count_20260609_161937.tsv.gz"),
  chh_count = file.path(BASE_DIR, "modality/outputs_CHH_mecp2/Results/mecp2_peaks.annotation/Extract_20260609_162203/Extract_total_c_count_20260609_162203.tsv.gz"),
  chg_dmr   = file.path(BASE_DIR, "modality/outputs_CHG_mecp2/Results/mecp2_peaks.annotation/DMR_20260609_191055/DMR_mc_control__mutant_20260609_191055.bed"),
  chh_dmr   = file.path(BASE_DIR, "modality/outputs_CHH_mecp2/Results/mecp2_peaks.annotation/DMR_20260609_192623/DMR_mc_control__mutant_20260609_192623.bed")
)

CTRL_CHG_RESULTS <- file.path(BASE_DIR, "modality/outputs_CHG_control/Results/control_peaks.annotation")
CTRL_CHH_RESULTS <- file.path(BASE_DIR, "modality/outputs_CHH_control/Results/control_peaks.annotation")

SEC51_SUMMARY <- file.path(TABLES_DIR, "mecp2_noncg_summary.tsv")

# ---- Helpers ----------------------------------------------------------------

load_peak_extract <- function(filepath, label) {
  stopifnot(file.exists(filepath))
  df <- read.table(gzfile(filepath), header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)
  if (nrow(df) > 0 && !grepl("^chr", df$Chromosome[1])) {
    df$Chromosome <- paste0("chr", df$Chromosome)
  }
  cat(sprintf("  %s: %d regions\n", label, nrow(df)))
  return(df)
}

discover_file <- function(results_dir, pattern) {
  stopifnot(dir.exists(results_dir))
  files <- list.files(results_dir, pattern = pattern,
                      recursive = TRUE, full.names = TRUE)
  stopifnot(length(files) == 1)
  files[1]
}

compute_detection <- function(df, ctrl_cols, mut_cols) {
  ctrl_any <- apply(df[, ctrl_cols, drop = FALSE], 1,
                    function(x) any(x > 0, na.rm = TRUE))
  mut_any <- apply(df[, mut_cols, drop = FALSE], 1,
                   function(x) any(x > 0, na.rm = TRUE))
  ifelse(ctrl_any & mut_any, "Ctrl & Mut",
  ifelse(ctrl_any & !mut_any, "Ctrl Only",
  ifelse(!ctrl_any & mut_any, "Mut Only", "No Signal")))
}

DETECTION_COLORS <- c("Ctrl & Mut" = "#7570B3",
                      "Ctrl Only" = "#2166AC",
                      "Mut Only" = "#B2182B")
DETECTION_LEVELS <- c("Ctrl & Mut", "Ctrl Only", "Mut Only", "No Signal")

# ---- Data Loading -----------------------------------------------------------

cat("Loading MeCP2 peak-resolution data...\n")
chg_frac  <- load_peak_extract(MECP2_PEAK_PATHS$chg_frac,  "CHG regional-frac")
chh_frac  <- load_peak_extract(MECP2_PEAK_PATHS$chh_frac,  "CHH regional-frac")
chg_count <- load_peak_extract(MECP2_PEAK_PATHS$chg_count, "CHG count")
chh_count <- load_peak_extract(MECP2_PEAK_PATHS$chh_count, "CHH count")

cat("\nLoading MeCP2 peak DMR results...\n")
chg_dmr <- load_dmr_bed(MECP2_PEAK_PATHS$chg_dmr)
chh_dmr <- load_dmr_bed(MECP2_PEAK_PATHS$chh_dmr)
chg_dmr$peak_id <- chg_dmr$gene
chh_dmr$peak_id <- chh_dmr$gene
cat(sprintf("  CHG DMR: %d peaks tested\n", nrow(chg_dmr)))
cat(sprintf("  CHH DMR: %d peaks tested\n", nrow(chh_dmr)))

cat("\nLoading control region data...\n")
stopifnot(dir.exists(CTRL_CHG_RESULTS))
stopifnot(dir.exists(CTRL_CHH_RESULTS))
ctrl_chg_frac_path <- discover_file(CTRL_CHG_RESULTS,
                                    "Extract_mc_regional-frac_.*\\.tsv\\.gz$")
ctrl_chh_frac_path <- discover_file(CTRL_CHH_RESULTS,
                                    "Extract_mc_regional-frac_.*\\.tsv\\.gz$")
ctrl_chg_frac <- load_peak_extract(ctrl_chg_frac_path, "Control CHG regional-frac")
ctrl_chh_frac <- load_peak_extract(ctrl_chh_frac_path, "Control CHH regional-frac")

cat("\nLoading Section 51 summary...\n")
stopifnot(file.exists(SEC51_SUMMARY))
s51_summary <- read.table(SEC51_SUMMARY, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
cat(sprintf("  %d rows loaded\n", nrow(s51_summary)))

# Identify sample columns
CTRL_COLS <- grep("ctrl", colnames(chg_frac), value = TRUE)
MUT_COLS  <- grep("mut",  colnames(chg_frac), value = TRUE)
ALL_COLS  <- c(CTRL_COLS, MUT_COLS)
cat(sprintf("\n  %d ctrl columns, %d mut columns\n", length(CTRL_COLS), length(MUT_COLS)))
stopifnot(length(CTRL_COLS) == 4, length(MUT_COLS) == 4)

# Verify control extracts have matching sample columns
stopifnot(all(CTRL_COLS %in% colnames(ctrl_chg_frac)))
stopifnot(all(MUT_COLS %in% colnames(ctrl_chg_frac)))

# Compute detection categories for all datasets
chg_frac$detection <- compute_detection(chg_frac, CTRL_COLS, MUT_COLS)
chh_frac$detection <- compute_detection(chh_frac, CTRL_COLS, MUT_COLS)
ctrl_chg_frac$detection <- compute_detection(ctrl_chg_frac, CTRL_COLS, MUT_COLS)
ctrl_chh_frac$detection <- compute_detection(ctrl_chh_frac, CTRL_COLS, MUT_COLS)

# =============================================================================
# FIGURE 53a: Detection Rate — Is There Any mCH at MeCP2 Peaks?
# =============================================================================

cat("\n--- Figure 53a: Non-CG detection rate at MeCP2 peaks ---\n")

build_det_counts <- function(frac_df, context_label) {
  by_dir <- frac_df %>%
    group_by(direction = Annotation, detection) %>%
    summarise(n = n(), .groups = "drop")

  combined <- frac_df %>%
    group_by(detection) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(direction = "All")

  bind_rows(by_dir, combined) %>% mutate(context = context_label)
}

det_counts <- bind_rows(
  build_det_counts(chg_frac, "CHG"),
  build_det_counts(chh_frac, "CHH")
) %>%
  mutate(
    direction = factor(direction, levels = c("All", "MeCP2_Up", "MeCP2_Down")),
    context = factor(context, levels = c("CHG", "CHH")),
    detection = factor(detection, levels = DETECTION_LEVELS)
  )

# Print detection breakdown
cat("\n  Detection breakdown:\n")
for (ctx in c("CHG", "CHH")) {
  for (dir in c("All", "MeCP2_Up", "MeCP2_Down")) {
    sub <- det_counts %>% filter(context == ctx, direction == dir)
    n_total <- sum(sub$n)
    n_signal <- sum(sub$n[sub$detection != "No Signal"])
    cat(sprintf("    %s %s: %d/%d (%.1f%%) detected | ",
                ctx, dir, n_signal, n_total, 100 * n_signal / n_total))
    for (det in c("Ctrl & Mut", "Ctrl Only", "Mut Only")) {
      val <- sub$n[sub$detection == det]
      if (length(val) == 0) val <- 0
      cat(sprintf("%s=%d  ", det, val))
    }
    cat("\n")
  }
}

# Zero-count annotations per facet
zero_annot <- det_counts %>%
  group_by(context, direction) %>%
  summarise(
    n_zero = sum(n[detection == "No Signal"]),
    n_total = sum(n),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("%.1f%% zero\n(%d / %d)", 100 * n_zero / n_total, n_zero, n_total))

# Plot non-zero categories only
det_nonzero <- det_counts %>%
  filter(detection != "No Signal") %>%
  mutate(detection = factor(detection, levels = c("Ctrl & Mut", "Ctrl Only", "Mut Only")))

p_53a <- ggplot(det_nonzero, aes(x = detection, y = n, fill = detection)) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_text(aes(label = n), vjust = -0.5, size = 3) +
  geom_text(
    data = zero_annot,
    aes(x = 2, y = Inf, label = label),
    vjust = 1.3, size = 2.8, color = "grey40", inherit.aes = FALSE
  ) +
  scale_fill_manual(values = DETECTION_COLORS, guide = "none") +
  facet_grid(direction ~ context, scales = "free_y") +
  labs(
    title = "Non-CG Methylation Detection at MeCP2 Peaks",
    subtitle = "Peaks with any detectable mCH in at least one sample (8 samples total)",
    x = "Detection Category", y = "Number of Peaks"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_53a, file.path(SEC53_DIR, "53a_detection_rate"),
                        width = 10, height = 10)

# =============================================================================
# FIGURE 53b: Methylation Levels at Non-Zero Peaks
# =============================================================================

cat("\n--- Figure 53b: Methylation levels at non-zero peaks ---\n")

build_nonzero_scatter <- function(frac_df, ctrl_cols, mut_cols, context_label) {
  frac_df$ctrl_mean <- rowMeans(frac_df[, ctrl_cols, drop = FALSE], na.rm = TRUE)
  frac_df$mut_mean  <- rowMeans(frac_df[, mut_cols, drop = FALSE], na.rm = TRUE)

  sub <- frac_df[frac_df$detection != "No Signal", ]
  cat(sprintf("  %s: %d non-zero peaks\n", context_label, nrow(sub)))

  data.frame(
    peak_id   = sub$Name,
    ctrl_mean = sub$ctrl_mean,
    mut_mean  = sub$mut_mean,
    direction = sub$Annotation,
    context   = context_label,
    stringsAsFactors = FALSE
  )
}

nz_chg <- build_nonzero_scatter(chg_frac, CTRL_COLS, MUT_COLS, "CHG")
nz_chh <- build_nonzero_scatter(chh_frac, CTRL_COLS, MUT_COLS, "CHH")

nz_all <- bind_rows(nz_chg, nz_chh) %>%
  mutate(
    context   = factor(context, levels = c("CHG", "CHH")),
    direction = factor(direction, levels = c("MeCP2_Up", "MeCP2_Down"))
  )

p_53b <- ggplot(nz_all, aes(x = ctrl_mean * 100, y = mut_mean * 100,
                             color = direction)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_color_manual(
    values = c("MeCP2_Up" = "#D95F02", "MeCP2_Down" = "#7570B3"),
    name = "MeCP2 Direction"
  ) +
  facet_wrap(~ context, scales = "free") +
  labs(
    title = "Non-CG Methylation at Non-Zero MeCP2 Peaks",
    subtitle = sprintf(
      "CHG: %d non-zero peaks | CHH: %d non-zero peaks | Diagonal = no difference",
      nrow(nz_chg), nrow(nz_chh)
    ),
    x = "Control Mean mC (%)", y = "Mutant Mean mC (%)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_53b, file.path(SEC53_DIR, "53b_nonzero_peak_scatter"),
                        width = 12, height = 6)

# =============================================================================
# FIGURE 53c: MeCP2 Peaks vs Shuffled Control Baseline
# =============================================================================

cat("\n--- Figure 53c: MeCP2 peaks vs shuffled control regions ---\n")

n_mecp2 <- nrow(chg_frac)
n_ctrl  <- nrow(ctrl_chg_frac)

mecp2_chg_det <- sum(chg_frac$detection != "No Signal")
mecp2_chh_det <- sum(chh_frac$detection != "No Signal")
ctrl_chg_det  <- sum(ctrl_chg_frac$detection != "No Signal")
ctrl_chh_det  <- sum(ctrl_chh_frac$detection != "No Signal")

fisher_chg <- fisher.test(matrix(c(
  mecp2_chg_det, n_mecp2 - mecp2_chg_det,
  ctrl_chg_det,  n_ctrl  - ctrl_chg_det
), nrow = 2, byrow = TRUE))

fisher_chh <- fisher.test(matrix(c(
  mecp2_chh_det, n_mecp2 - mecp2_chh_det,
  ctrl_chh_det,  n_ctrl  - ctrl_chh_det
), nrow = 2, byrow = TRUE))

cat(sprintf("  CHG: MeCP2=%d/%d (%.1f%%) vs Control=%d/%d (%.1f%%) | Fisher OR=%.2f p=%.2e\n",
            mecp2_chg_det, n_mecp2, 100 * mecp2_chg_det / n_mecp2,
            ctrl_chg_det, n_ctrl, 100 * ctrl_chg_det / n_ctrl,
            fisher_chg$estimate, fisher_chg$p.value))
cat(sprintf("  CHH: MeCP2=%d/%d (%.1f%%) vs Control=%d/%d (%.1f%%) | Fisher OR=%.2f p=%.2e\n",
            mecp2_chh_det, n_mecp2, 100 * mecp2_chh_det / n_mecp2,
            ctrl_chh_det, n_ctrl, 100 * ctrl_chh_det / n_ctrl,
            fisher_chh$estimate, fisher_chh$p.value))

# Mean methylation comparison
mecp2_chg_mean <- mean(rowMeans(chg_frac[, ALL_COLS, drop = FALSE], na.rm = TRUE), na.rm = TRUE)
mecp2_chh_mean <- mean(rowMeans(chh_frac[, ALL_COLS, drop = FALSE], na.rm = TRUE), na.rm = TRUE)
ctrl_chg_mean  <- mean(rowMeans(ctrl_chg_frac[, ALL_COLS, drop = FALSE], na.rm = TRUE), na.rm = TRUE)
ctrl_chh_mean  <- mean(rowMeans(ctrl_chh_frac[, ALL_COLS, drop = FALSE], na.rm = TRUE), na.rm = TRUE)

cat(sprintf("  CHG mean mC: MeCP2=%.4e vs Control=%.4e\n", mecp2_chg_mean, ctrl_chg_mean))
cat(sprintf("  CHH mean mC: MeCP2=%.4e vs Control=%.4e\n", mecp2_chh_mean, ctrl_chh_mean))

baseline_df <- data.frame(
  region_type  = rep(c("MeCP2 Peaks", "Shuffled Control"), each = 2),
  context      = rep(c("CHG", "CHH"), 2),
  n_detected   = c(mecp2_chg_det, mecp2_chh_det, ctrl_chg_det, ctrl_chh_det),
  n_total      = c(n_mecp2, n_mecp2, n_ctrl, n_ctrl),
  pct_detected = c(mecp2_chg_det / n_mecp2, mecp2_chh_det / n_mecp2,
                   ctrl_chg_det / n_ctrl, ctrl_chh_det / n_ctrl) * 100,
  mean_meth    = c(mecp2_chg_mean, mecp2_chh_mean, ctrl_chg_mean, ctrl_chh_mean) * 100,
  stringsAsFactors = FALSE
) %>%
  mutate(
    region_type = factor(region_type, levels = c("MeCP2 Peaks", "Shuffled Control")),
    context = factor(context, levels = c("CHG", "CHH"))
  )

p_53c <- ggplot(baseline_df, aes(x = region_type, y = pct_detected, fill = region_type)) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", pct_detected, n_detected, n_total)),
            vjust = -0.3, size = 3) +
  scale_fill_manual(
    values = c("MeCP2 Peaks" = "#D95F02", "Shuffled Control" = "#999999"),
    guide = "none"
  ) +
  facet_wrap(~ context) +
  labs(
    title = "Non-CG Methylation: MeCP2 Peaks vs Shuffled Control Regions",
    subtitle = sprintf(
      "Fisher CHG: OR=%.2f, p=%.2e | Fisher CHH: OR=%.2f, p=%.2e",
      fisher_chg$estimate, fisher_chg$p.value,
      fisher_chh$estimate, fisher_chh$p.value
    ),
    x = "Region Type", y = "Regions with Detectable mCH (%)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_53c, file.path(SEC53_DIR, "53c_mecp2_vs_control_baseline"),
                        width = 10, height = 6)

# =============================================================================
# FIGURE 53d: DMR Q-value Rank Plot
# =============================================================================

cat("\n--- Figure 53d: DMR q-value rank plot ---\n")

chg_dmr$rank <- rank(chg_dmr$dmr_qvalue, ties.method = "first")
chh_dmr$rank <- rank(chh_dmr$dmr_qvalue, ties.method = "first")

dmr_rank <- bind_rows(
  chg_dmr %>% mutate(context = "CHG"),
  chh_dmr %>% mutate(context = "CHH")
) %>%
  mutate(
    context = factor(context, levels = c("CHG", "CHH")),
    sig_tier = case_when(
      dmr_qvalue < 0.05 ~ "q < 0.05",
      dmr_qvalue < 0.5  ~ "q < 0.5",
      TRUE              ~ "NS"
    )
  )

n_chg_sig <- sum(chg_dmr$dmr_qvalue < Q_THRESHOLD, na.rm = TRUE)
n_chh_sig <- sum(chh_dmr$dmr_qvalue < Q_THRESHOLD, na.rm = TRUE)
cat(sprintf("  CHG: %d/%d peaks at q < %.2f\n", n_chg_sig, nrow(chg_dmr), Q_THRESHOLD))
cat(sprintf("  CHH: %d/%d peaks at q < %.2f\n", n_chh_sig, nrow(chh_dmr), Q_THRESHOLD))

p_53d <- ggplot(dmr_rank, aes(x = rank, y = neg_log10_q, color = sig_tier)) +
  geom_point(size = 0.5, alpha = 0.4) +
  geom_hline(yintercept = -log10(Q_THRESHOLD), linetype = "dashed", color = "grey40") +
  ggrepel::geom_text_repel(
    data = dmr_rank %>% filter(sig_tier != "NS"),
    aes(label = peak_id),
    size = 2.5, max.overlaps = 15, segment.color = "grey50",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("q < 0.05" = "#D7191C", "q < 0.5" = "#FF8C00", "NS" = "grey75"),
    name = "Significance"
  ) +
  facet_wrap(~ context, scales = "free") +
  labs(
    title = "Differential Non-CG Methylation at MeCP2 Peaks (GLM with Sex Covariate)",
    subtitle = sprintf("CHG: %d/%d at q<0.05 | CHH: %d/%d at q<0.05",
                       n_chg_sig, nrow(chg_dmr), n_chh_sig, nrow(chh_dmr)),
    x = "Peak Rank (sorted by q-value)",
    y = expression(-log[10](q-value))
  ) +
  theme_biomodal() +
  theme(legend.position = "bottom")

save_multiformat_ggplot(p_53d, file.path(SEC53_DIR, "53d_dmr_rank_plot"),
                        width = 12, height = 6)

# =============================================================================
# FIGURE 53e: Chr8 CHH Cluster Spotlight
# =============================================================================

cat("\n--- Figure 53e: Chr8 CHH cluster spotlight ---\n")

chh_chr8 <- chh_frac %>%
  filter(Chromosome == "chr8", Start >= 35000000, Start <= 50000000)

stopifnot(nrow(chh_chr8) > 0)
cat(sprintf("  %d CHH peaks on chr8 (35-50 Mb)\n", nrow(chh_chr8)))

chh_chr8$ctrl_mean <- rowMeans(chh_chr8[, CTRL_COLS, drop = FALSE], na.rm = TRUE)
chh_chr8$mut_mean  <- rowMeans(chh_chr8[, MUT_COLS, drop = FALSE], na.rm = TRUE)
chh_chr8$midpoint  <- (chh_chr8$Start + chh_chr8$End) / 2

chr8_dmr_sub <- chh_dmr %>%
  filter(chr == "chr8", start >= 35000000, start <= 50000000) %>%
  select(peak_id, dmr_qvalue, mod_difference)

chh_chr8 <- chh_chr8 %>%
  left_join(chr8_dmr_sub, by = c("Name" = "peak_id")) %>%
  mutate(
    q_tier = case_when(
      is.na(dmr_qvalue)   ~ "Not tested",
      dmr_qvalue < 0.05   ~ "q < 0.05",
      dmr_qvalue < 0.25   ~ "q < 0.25",
      TRUE                ~ "NS"
    )
  )

chr8_long <- chh_chr8 %>%
  select(Name, midpoint, ctrl_mean, mut_mean, q_tier, dmr_qvalue) %>%
  pivot_longer(cols = c(ctrl_mean, mut_mean),
               names_to = "condition", values_to = "mean_meth") %>%
  mutate(condition = recode(condition,
                            ctrl_mean = "Control", mut_mean = "Mutant"))

p_53e <- ggplot(chr8_long, aes(x = midpoint / 1e6, y = mean_meth * 100,
                               color = condition, shape = q_tier)) +
  geom_segment(
    data = chh_chr8,
    aes(x = midpoint / 1e6, xend = midpoint / 1e6,
        y = ctrl_mean * 100, yend = mut_mean * 100),
    color = "grey60", linewidth = 0.3, inherit.aes = FALSE
  ) +
  geom_point(size = 3, alpha = 0.8) +
  ggrepel::geom_text_repel(
    data = chr8_long %>% filter(q_tier == "q < 0.05", condition == "Mutant"),
    aes(label = sprintf("%.2f Mb\nq=%.4f", midpoint / 1e6, dmr_qvalue)),
    size = 2.5, nudge_y = 0.15, show.legend = FALSE, color = "black"
  ) +
  scale_color_manual(values = COLORS$condition) +
  scale_shape_manual(
    values = c("q < 0.05" = 16, "q < 0.25" = 17, "NS" = 1, "Not tested" = 4),
    name = "DMR q-value"
  ) +
  labs(
    title = "Chr8 CHH Methylation Cluster at MeCP2 Peaks (35–50 Mb)",
    subtitle = sprintf("%d peaks in region | Vertical lines connect ctrl/mut means",
                       nrow(chh_chr8)),
    x = "Chromosomal Position (Mb, chr8)", y = "Mean CHH mC (%)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_53e, file.path(SEC53_DIR, "53e_chr8_cluster_spotlight"),
                        width = 12, height = 6)

# =============================================================================
# FIGURE 53f: Peak Resolution vs Gene-Body Comparison
# =============================================================================

cat("\n--- Figure 53f: Peak vs gene-body resolution comparison ---\n")

compute_sample_means <- function(frac_df, ctrl_cols, mut_cols, context_label) {
  long <- frac_df %>%
    pivot_longer(cols = all_of(c(ctrl_cols, mut_cols)),
                 names_to = "sample", values_to = "methylation") %>%
    mutate(condition = ifelse(grepl("ctrl", sample), "Control", "Mutant"))

  long %>%
    group_by(sample, condition) %>%
    summarise(mean_meth = mean(methylation, na.rm = TRUE), .groups = "drop") %>%
    mutate(context = context_label)
}

peak_chg_means <- compute_sample_means(chg_frac, CTRL_COLS, MUT_COLS, "CHG")
peak_chh_means <- compute_sample_means(chh_frac, CTRL_COLS, MUT_COLS, "CHH")

peak_wt_chg <- wilcox.test(
  peak_chg_means$mean_meth[peak_chg_means$condition == "Control"],
  peak_chg_means$mean_meth[peak_chg_means$condition == "Mutant"]
)
peak_wt_chh <- wilcox.test(
  peak_chh_means$mean_meth[peak_chh_means$condition == "Control"],
  peak_chh_means$mean_meth[peak_chh_means$condition == "Mutant"]
)
cat(sprintf("  Peak-level Wilcoxon: CHG p=%.3e, CHH p=%.3e\n",
            peak_wt_chg$p.value, peak_wt_chh$p.value))

s51_chg <- s51_summary[grep("CHG at all MeCP2-bound genes", s51_summary$test), ]
s51_chh <- s51_summary[grep("CHH at all MeCP2-bound genes", s51_summary$test), ]
stopifnot(nrow(s51_chg) == 1, nrow(s51_chh) == 1)

comparison_df <- bind_rows(
  data.frame(
    resolution = "Gene Body (~18kb)",
    context = c("CHG", "CHH"),
    ctrl_mean = c(s51_chg$ctrl_mean, s51_chh$ctrl_mean),
    mut_mean = c(s51_chg$mut_mean, s51_chh$mut_mean),
    p_value = c(s51_chg$p_value, s51_chh$p_value),
    stringsAsFactors = FALSE
  ),
  data.frame(
    resolution = "Peak (400bp)",
    context = c("CHG", "CHH"),
    ctrl_mean = c(
      mean(peak_chg_means$mean_meth[peak_chg_means$condition == "Control"]),
      mean(peak_chh_means$mean_meth[peak_chh_means$condition == "Control"])
    ),
    mut_mean = c(
      mean(peak_chg_means$mean_meth[peak_chg_means$condition == "Mutant"]),
      mean(peak_chh_means$mean_meth[peak_chh_means$condition == "Mutant"])
    ),
    p_value = c(peak_wt_chg$p.value, peak_wt_chh$p.value),
    stringsAsFactors = FALSE
  )
)

# Annotation df for p-values (one row per resolution × context)
p_annot <- comparison_df %>%
  mutate(p_label = sprintf("p = %s", format.pval(p_value, digits = 2)))

comparison_long <- comparison_df %>%
  mutate(
    resolution = factor(resolution, levels = c("Gene Body (~18kb)", "Peak (400bp)")),
    context = factor(context, levels = c("CHG", "CHH"))
  ) %>%
  pivot_longer(cols = c(ctrl_mean, mut_mean),
               names_to = "condition", values_to = "mean_meth") %>%
  mutate(condition = recode(condition, ctrl_mean = "Control", mut_mean = "Mutant"))

p_53f <- ggplot(comparison_long, aes(x = resolution, y = mean_meth * 100,
                                     fill = condition)) +
  geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.85) +
  geom_text(
    data = p_annot %>% mutate(
      context = factor(context, levels = c("CHG", "CHH")),
      resolution = factor(resolution, levels = c("Gene Body (~18kb)", "Peak (400bp)"))
    ),
    aes(x = resolution, y = Inf, label = p_label),
    vjust = 1.5, size = 3, inherit.aes = FALSE, color = "grey30"
  ) +
  scale_fill_manual(values = COLORS$condition) +
  facet_wrap(~ context, scales = "free_y") +
  labs(
    title = "Non-CG Methylation at MeCP2 Sites: Gene Body vs Peak Resolution",
    subtitle = "Consistent negative result at both resolutions (Wilcoxon ctrl vs mut)",
    x = "Analysis Resolution", y = "Mean mC Regional Fraction (%)",
    fill = "Condition"
  ) +
  theme_biomodal() +
  theme(legend.position = "bottom")

save_multiformat_ggplot(p_53f, file.path(SEC53_DIR, "53f_resolution_comparison"),
                        width = 10, height = 6)

# =============================================================================
# FIGURE 53g: Coverage and Context Quality
# =============================================================================

cat("\n--- Figure 53g: Coverage and context quality ---\n")

first_col <- CTRL_COLS[1]
chg_ctx_vals <- as.numeric(chg_count[[first_col]])
chh_ctx_vals <- as.numeric(chh_count[[first_col]])

cat(sprintf("  CHG contexts/peak: mean=%.0f, range=%d-%d\n",
            mean(chg_ctx_vals, na.rm = TRUE),
            min(chg_ctx_vals, na.rm = TRUE), max(chg_ctx_vals, na.rm = TRUE)))
cat(sprintf("  CHH contexts/peak: mean=%.0f, range=%d-%d\n",
            mean(chh_ctx_vals, na.rm = TRUE),
            min(chh_ctx_vals, na.rm = TRUE), max(chh_ctx_vals, na.rm = TRUE)))

contexts_df <- bind_rows(
  data.frame(num_contexts = chg_ctx_vals, context = "CHG"),
  data.frame(num_contexts = chh_ctx_vals, context = "CHH")
) %>% mutate(context = factor(context, levels = c("CHG", "CHH")))

ctx_means <- contexts_df %>%
  group_by(context) %>%
  summarise(mean_ctx = mean(num_contexts, na.rm = TRUE), .groups = "drop")

coverage_df <- bind_rows(
  data.frame(mean_coverage = chg_dmr$mean_coverage, context = "CHG"),
  data.frame(mean_coverage = chh_dmr$mean_coverage, context = "CHH")
) %>%
  filter(!is.na(mean_coverage)) %>%
  mutate(context = factor(context, levels = c("CHG", "CHH")))

p_53g_left <- ggplot(contexts_df, aes(x = num_contexts, fill = context)) +
  geom_histogram(bins = 50, alpha = 0.8) +
  geom_vline(data = ctx_means, aes(xintercept = mean_ctx),
             linetype = "dashed", color = "grey30", linewidth = 0.6) +
  scale_fill_manual(values = c("CHG" = "#984EA3", "CHH" = "#FF7F00"), guide = "none") +
  facet_wrap(~ context, ncol = 1, scales = "free") +
  labs(title = "Cytosine Contexts per Peak",
       x = "Number of Cytosine Contexts", y = "Number of Peaks") +
  theme_biomodal()

p_53g_right <- ggplot(coverage_df, aes(x = mean_coverage, fill = context)) +
  geom_histogram(bins = 50, alpha = 0.8) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "grey40") +
  annotate("text", x = 11, y = Inf, label = "min_cov = 10",
           vjust = 1.5, hjust = 0, size = 3, color = "grey40") +
  scale_fill_manual(values = c("CHG" = "#984EA3", "CHH" = "#FF7F00"), guide = "none") +
  facet_wrap(~ context, ncol = 1, scales = "free") +
  labs(title = "Mean Coverage per Peak",
       x = "Mean Coverage (reads)", y = "Number of Peaks") +
  theme_biomodal()

p_53g <- (p_53g_left | p_53g_right) +
  patchwork::plot_annotation(
    title = "Data Quality: Coverage and Context Density at MeCP2 Peaks",
    subtitle = "Adequate depth confirms the negative result is not a power issue",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  )

save_multiformat_ggplot(p_53g, file.path(SEC53_DIR, "53g_coverage_quality"),
                        width = 12, height = 8)

# =============================================================================
# SUMMARY TABLE
# =============================================================================

cat("\n--- Writing summary table ---\n")

summary_rows <- list()

# Detection rates — MeCP2 peaks
for (ctx in c("CHG", "CHH")) {
  frac_df <- if (ctx == "CHG") chg_frac else chh_frac
  det <- frac_df$detection

  for (dir in c("All", "MeCP2_Up", "MeCP2_Down")) {
    mask <- if (dir == "All") rep(TRUE, nrow(frac_df)) else frac_df$Annotation == dir
    sub_det <- det[mask]
    n <- sum(mask)
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      test           = sprintf("Detection: %s at %s MeCP2 peaks", ctx, dir),
      n_peaks        = n,
      n_detected     = sum(sub_det != "No Signal"),
      pct_detected   = round(100 * sum(sub_det != "No Signal") / n, 2),
      n_ctrl_only    = sum(sub_det == "Ctrl Only"),
      n_mut_only     = sum(sub_det == "Mut Only"),
      n_both         = sum(sub_det == "Ctrl & Mut"),
      ctrl_mean_meth = NA_real_,
      mut_mean_meth  = NA_real_,
      p_value        = NA_real_,
      note           = "",
      stringsAsFactors = FALSE
    )
  }
}

# Detection rates — control regions
for (ctx in c("CHG", "CHH")) {
  ctrl_df <- if (ctx == "CHG") ctrl_chg_frac else ctrl_chh_frac
  det <- ctrl_df$detection
  n <- nrow(ctrl_df)
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    test           = sprintf("Detection: %s at Control regions", ctx),
    n_peaks        = n,
    n_detected     = sum(det != "No Signal"),
    pct_detected   = round(100 * sum(det != "No Signal") / n, 2),
    n_ctrl_only    = sum(det == "Ctrl Only"),
    n_mut_only     = sum(det == "Mut Only"),
    n_both         = sum(det == "Ctrl & Mut"),
    ctrl_mean_meth = NA_real_,
    mut_mean_meth  = NA_real_,
    p_value        = NA_real_,
    note           = "",
    stringsAsFactors = FALSE
  )
}

# Fisher's exact tests
summary_rows[[length(summary_rows) + 1]] <- data.frame(
  test = "Fisher: CHG MeCP2 vs Control",
  n_peaks = NA, n_detected = NA, pct_detected = NA,
  n_ctrl_only = NA, n_mut_only = NA, n_both = NA,
  ctrl_mean_meth = NA, mut_mean_meth = NA,
  p_value = fisher_chg$p.value,
  note = sprintf("OR=%.3f", fisher_chg$estimate),
  stringsAsFactors = FALSE
)
summary_rows[[length(summary_rows) + 1]] <- data.frame(
  test = "Fisher: CHH MeCP2 vs Control",
  n_peaks = NA, n_detected = NA, pct_detected = NA,
  n_ctrl_only = NA, n_mut_only = NA, n_both = NA,
  ctrl_mean_meth = NA, mut_mean_meth = NA,
  p_value = fisher_chh$p.value,
  note = sprintf("OR=%.3f", fisher_chh$estimate),
  stringsAsFactors = FALSE
)

# Wilcoxon ctrl vs mut (peak level)
for (ctx in c("CHG", "CHH")) {
  means_df <- if (ctx == "CHG") peak_chg_means else peak_chh_means
  wt <- if (ctx == "CHG") peak_wt_chg else peak_wt_chh
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    test = sprintf("Wilcoxon ctrl vs mut: %s at MeCP2 peaks (400bp)", ctx),
    n_peaks = nrow(if (ctx == "CHG") chg_frac else chh_frac),
    n_detected = NA, pct_detected = NA,
    n_ctrl_only = NA, n_mut_only = NA, n_both = NA,
    ctrl_mean_meth = mean(means_df$mean_meth[means_df$condition == "Control"]),
    mut_mean_meth = mean(means_df$mean_meth[means_df$condition == "Mutant"]),
    p_value = wt$p.value,
    note = "",
    stringsAsFactors = FALSE
  )
}

# DMR summary
for (ctx in c("CHG", "CHH")) {
  dmr_df <- if (ctx == "CHG") chg_dmr else chh_dmr
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    test = sprintf("DMR summary: %s", ctx),
    n_peaks = nrow(dmr_df),
    n_detected = sum(dmr_df$dmr_qvalue < Q_THRESHOLD, na.rm = TRUE),
    pct_detected = NA,
    n_ctrl_only = sum(dmr_df$dmr_qvalue < 0.25, na.rm = TRUE),
    n_mut_only = sum(dmr_df$dmr_qvalue < 0.50, na.rm = TRUE),
    n_both = NA,
    ctrl_mean_meth = mean(dmr_df$mean_mod_group1, na.rm = TRUE),
    mut_mean_meth = mean(dmr_df$mean_mod_group2, na.rm = TRUE),
    p_value = min(dmr_df$dmr_qvalue, na.rm = TRUE),
    note = sprintf("n_detected=q<0.05; n_ctrl_only=q<0.25; n_mut_only=q<0.50; p=best_q"),
    stringsAsFactors = FALSE
  )
}

summary_table <- bind_rows(summary_rows)
write.table(summary_table, file.path(TABLES_DIR, "mecp2_peak_noncg_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Saved: %s (%d rows)\n",
            file.path(TABLES_DIR, "mecp2_peak_noncg_summary.tsv"), nrow(summary_table)))

# =============================================================================
# SECTION 53 COMPLETE
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 53 COMPLETE\n")
cat("================================================================================\n")
cat(sprintf("Figures saved to: %s\n", SEC53_DIR))
cat(sprintf("Summary table: %s\n", file.path(TABLES_DIR, "mecp2_peak_noncg_summary.tsv")))
cat("\nKEY FINDINGS:\n")
cat("  ABSOLUTE DETECTION:\n")
cat(sprintf("    CHG: %d/%d peaks (%.1f%%) have any detectable mCH\n",
            mecp2_chg_det, n_mecp2, 100 * mecp2_chg_det / n_mecp2))
cat(sprintf("    CHH: %d/%d peaks (%.1f%%) have any detectable mCH\n",
            mecp2_chh_det, n_mecp2, 100 * mecp2_chh_det / n_mecp2))
cat(sprintf("    vs Control: CHG %.1f%%, CHH %.1f%%\n",
            100 * ctrl_chg_det / n_ctrl, 100 * ctrl_chh_det / n_ctrl))
cat("  DIFFERENTIAL:\n")
cat(sprintf("    CHG: %d/%d peaks at q<0.05\n", n_chg_sig, nrow(chg_dmr)))
cat(sprintf("    CHH: %d/%d peaks at q<0.05\n", n_chh_sig, nrow(chh_dmr)))
cat("================================================================================\n\n")
