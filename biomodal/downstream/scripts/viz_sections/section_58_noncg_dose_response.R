# biomodal/downstream/scripts/viz_sections/section_58_noncg_dose_response.R
# Section 58: Non-CG Methylation Dose-Response on MeCP2 Residual
#
# Tests whether Ecker CH (non-CG) methylation level scales the magnitude of
# CG-unexplained MeCP2 enrichment at MeCP2-Up peaks. Bins peaks by Ecker CH
# quartile and tests for monotonic trend in MeCP2 residual.
#
# Also includes a corrected version of the section 57b correlation plot
# (MeCP2-Up only, fixing the Simpson's paradox from pooling Up+Down).
#
# Adds quantile regression (median regression) as a robust alternative to the
# OLS residuals from section 56.
#
# Input:  tables/57_ecker_peak_signal.tsv (from section 57)
#         tables/56_mecp2_peak_annotated.tsv (from section 56, for cg_mc_log2fc)
# Usage:  cd downstream && Rscript scripts/viz_sections/section_58_noncg_dose_response.R

source("scripts/viz_sections/_shared_config.R")

library(quantreg)
library(clinfun)

# =============================================================================
# SECTION 58 CONFIGURATION
# =============================================================================

SEC58_DIR <- file.path(OUTPUT_DIR, "58_noncg_dose_response")
dir.create(SEC58_DIR, recursive = TRUE, showWarnings = FALSE)

SIGNAL_FILE <- file.path(TABLES_DIR, "57_ecker_peak_signal.tsv")
PEAK_FILE   <- file.path(TABLES_DIR, "56_mecp2_peak_annotated.tsv")

# =============================================================================
# LOAD DATA
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 58: NON-CG METHYLATION DOSE-RESPONSE ON MeCP2 RESIDUAL\n")
cat("================================================================================\n\n")

if (!file.exists(SIGNAL_FILE)) {
  stop("Section 57 signal file not found: ", SIGNAL_FILE,
       "\nRun section_57_ecker_noncg_validation.R first.")
}
if (!file.exists(PEAK_FILE)) {
  stop("Section 56 peak file not found: ", PEAK_FILE,
       "\nRun section_56_mecp2_peak_cg_correction.R first.")
}

sig <- read.table(SIGNAL_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("Loaded %d peaks from section 57 signal file\n", nrow(sig)))

pk <- read.table(PEAK_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                 fill = TRUE, comment.char = "")
cat(sprintf("Loaded %d peaks from section 56 annotated file\n", nrow(pk)))

peaks <- merge(sig, pk[, c("Chr", "Start", "End", "cg_mc_log2fc")],
               by.x = c("chr", "start", "end"),
               by.y = c("Chr", "Start", "End"),
               all.x = TRUE)
cat(sprintf("Merged: %d peaks with cg_mc_log2fc\n", sum(!is.na(peaks$cg_mc_log2fc))))

cat(sprintf("  MeCP2 Up:           %d\n", sum(peaks$mecp2_class == "MeCP2 Up")))
cat(sprintf("  MeCP2 Down:         %d\n", sum(peaks$mecp2_class == "MeCP2 Down")))
cat(sprintf("  Not Significant:    %d\n\n", sum(peaks$mecp2_class == "Not Significant")))

# =============================================================================
# QUANTILE REGRESSION (ROBUST ALTERNATIVE TO OLS)
# =============================================================================

cat("--- Quantile regression: MeCP2 Fold ~ CG mC log2FC (tau=0.5) ---\n\n")

valid <- !is.na(peaks$Fold) & !is.na(peaks$cg_mc_log2fc) & is.finite(peaks$cg_mc_log2fc)

ols_model <- lm(Fold ~ cg_mc_log2fc, data = peaks, subset = valid)
qr_model  <- rq(Fold ~ cg_mc_log2fc, data = peaks, subset = valid, tau = 0.5)
qr_summary <- summary(qr_model, se = "boot", R = 1000)

cat(sprintf("  OLS:  intercept=%.4f, slope=%.4f, R²=%.4f\n",
            coef(ols_model)[1], coef(ols_model)[2], summary(ols_model)$r.squared))
cat(sprintf("  QR:   intercept=%.4f, slope=%.4f\n",
            coef(qr_model)[1], coef(qr_model)[2]))

peaks$ols_residual <- peaks$residual
peaks$qr_residual <- NA_real_
peaks$qr_residual[valid] <- residuals(qr_model)

ols_qr_cor <- cor(peaks$ols_residual[valid], peaks$qr_residual[valid],
                  method = "spearman", use = "complete.obs")
cat(sprintf("  OLS vs QR residual Spearman correlation: %.4f\n\n", ols_qr_cor))

reg_compare <- data.frame(
  model = c("OLS", "Quantile (tau=0.5)"),
  intercept = c(coef(ols_model)[1], coef(qr_model)[1]),
  slope = c(coef(ols_model)[2], coef(qr_model)[2]),
  r_squared = c(summary(ols_model)$r.squared, NA),
  n = sum(valid),
  residual_cor_spearman = ols_qr_cor,
  stringsAsFactors = FALSE, row.names = NULL
)

reg_path <- file.path(TABLES_DIR, "58_noncg_regression_comparison.tsv")
write.table(reg_compare, reg_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Regression comparison: %s\n\n", reg_path))

# =============================================================================
# HELPER: QUARTILE DOSE-RESPONSE ANALYSIS
# =============================================================================

run_quartile_analysis <- function(df, ch_col, resid_col, label) {
  vals <- df[[ch_col]]
  resids <- df[[resid_col]]
  valid_mask <- !is.na(vals) & !is.na(resids) & is.finite(vals) & is.finite(resids)
  vals <- vals[valid_mask]
  resids <- resids[valid_mask]

  if (length(vals) < 20) {
    cat(sprintf("  %s: too few valid peaks (%d), skipping\n", label, length(vals)))
    return(NULL)
  }

  breaks <- quantile(vals, probs = c(0, 0.25, 0.5, 0.75, 1))
  breaks[1] <- breaks[1] - 1e-10
  quartile <- as.integer(cut(vals, breaks = breaks, labels = 1:4))

  rows <- list()
  for (q in 1:4) {
    idx <- which(quartile == q)
    if (length(idx) < 5) next
    r <- resids[idx]
    ch <- vals[idx]

    set.seed(42 + q)
    boot_medians <- replicate(1000, median(sample(r, replace = TRUE)))

    wt <- wilcox.test(r, mu = 0)
    rows[[q]] <- data.frame(
      quartile = q,
      quartile_label = c("Q1 (low CH)", "Q2", "Q3", "Q4 (high CH)")[q],
      n = length(idx),
      ch_range = sprintf("%.6f - %.6f", min(ch), max(ch)),
      median_ch = median(ch),
      mean_residual = mean(r),
      median_residual = median(r),
      sd_residual = sd(r),
      ci_lo = quantile(boot_medians, 0.025, names = FALSE),
      ci_hi = quantile(boot_medians, 0.975, names = FALSE),
      wilcox_p_vs_zero = wt$p.value,
      stringsAsFactors = FALSE, row.names = NULL
    )
    cat(sprintf("  %s Q%d: n=%d, median_CH=%.6f, median_resid=%.4f, p=%.2e\n",
                label, q, length(idx), median(ch), median(r), wt$p.value))
  }

  result_df <- do.call(rbind, rows)

  jt <- jonckheere.test(resids, quartile, alternative = "increasing")
  kw <- kruskal.test(resids ~ quartile)

  cat(sprintf("  %s Jonckheere (increasing): stat=%.1f, p=%.2e\n",
              label, jt$statistic, jt$p.value))
  cat(sprintf("  %s Kruskal-Wallis: chi²=%.2f, p=%.2e\n\n",
              label, kw$statistic, kw$p.value))

  jt_dec <- jonckheere.test(resids, quartile, alternative = "decreasing")

  tests_df <- data.frame(
    test = c("Jonckheere-Terpstra (increasing)", "Jonckheere-Terpstra (decreasing)",
             "Kruskal-Wallis"),
    statistic = c(jt$statistic, jt_dec$statistic, kw$statistic),
    p_value = c(jt$p.value, jt_dec$p.value, kw$p.value),
    n = length(vals),
    stringsAsFactors = FALSE, row.names = NULL
  )

  list(quartiles = result_df, tests = tests_df, raw_quartile = quartile,
       raw_resid = resids, raw_ch = vals)
}

# =============================================================================
# CORRECTED 57b CORRELATION (MeCP2-Up ONLY)
# =============================================================================

cat("--- Corrected 57b: Ecker CH vs residual (MeCP2-Up only) ---\n\n")

up_peaks <- peaks[peaks$mecp2_class == "MeCP2 Up" &
                  !is.na(peaks$ols_residual) & !is.na(peaks$ecker_ch), ]
cat(sprintf("  MeCP2-Up peaks with valid residual + Ecker CH: %d\n", nrow(up_peaks)))

up_cor <- cor.test(up_peaks$ols_residual, up_peaks$ecker_ch, method = "spearman")
cat(sprintf("  Spearman rho = %.4f, p = %.2e (MeCP2-Up only)\n", up_cor$estimate, up_cor$p.value))

all_sig <- peaks[peaks$mecp2_class != "Not Significant" &
                 !is.na(peaks$ols_residual) & !is.na(peaks$ecker_ch), ]
pooled_cor <- cor.test(all_sig$ols_residual, all_sig$ecker_ch, method = "spearman")
cat(sprintf("  Spearman rho = %.4f, p = %.2e (pooled Up+Down, for reference)\n\n",
            pooled_cor$estimate, pooled_cor$p.value))

up_peaks$peak_group <- ifelse(up_peaks$cg_class == "CG Increase",
                               "CG-Concordant", "Non-CG Candidate")

# =============================================================================
# ECKER CH QUARTILE ANALYSIS — MeCP2-Up
# =============================================================================

cat("--- Ecker CH quartile analysis: MeCP2-Up peaks ---\n")

up_ols <- run_quartile_analysis(up_peaks, "ecker_ch", "ols_residual", "Up_OLS")
up_qr  <- run_quartile_analysis(up_peaks, "ecker_ch", "qr_residual", "Up_QR")

up_ols$quartiles$residual_type <- "OLS"
up_qr$quartiles$residual_type  <- "Quantile Regression"
up_combined <- rbind(up_ols$quartiles, up_qr$quartiles)

up_path <- file.path(TABLES_DIR, "58_noncg_dose_response_up.tsv")
write.table(up_combined, up_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  MeCP2-Up quartile table: %s\n\n", up_path))

# =============================================================================
# NEGATIVE CONTROL — MeCP2-Down
# =============================================================================

cat("--- Ecker CH quartile analysis: MeCP2-Down peaks (negative control) ---\n")

down_peaks <- peaks[peaks$mecp2_class == "MeCP2 Down" &
                    !is.na(peaks$ols_residual) & !is.na(peaks$ecker_ch), ]
cat(sprintf("  MeCP2-Down peaks: %d\n", nrow(down_peaks)))

down_ols <- run_quartile_analysis(down_peaks, "ecker_ch", "ols_residual", "Down_OLS")

if (!is.null(down_ols)) {
  down_ols$quartiles$residual_type <- "OLS"
  down_path <- file.path(TABLES_DIR, "58_noncg_dose_response_down.tsv")
  write.table(down_ols$quartiles, down_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  MeCP2-Down quartile table: %s\n\n", down_path))
}

# =============================================================================
# COMBINED TREND TESTS
# =============================================================================

cat("--- Trend test summary ---\n")

all_tests <- rbind(
  cbind(group = "MeCP2-Up (OLS)", up_ols$tests),
  cbind(group = "MeCP2-Up (QR)", up_qr$tests)
)
if (!is.null(down_ols)) {
  all_tests <- rbind(all_tests, cbind(group = "MeCP2-Down (OLS)", down_ols$tests))
}

tests_path <- file.path(TABLES_DIR, "58_noncg_trend_tests.tsv")
write.table(all_tests, tests_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Trend tests: %s\n\n", tests_path))

# =============================================================================
# JOINT CG × CH STRATIFICATION
# =============================================================================

cat("--- Joint CG × CH stratification (MeCP2-Up) ---\n")

joint <- up_peaks[!is.na(up_peaks$ecker_cg) & !is.na(up_peaks$ecker_ch) &
                  !is.na(up_peaks$ols_residual), ]
cg_med <- median(joint$ecker_cg)
ch_med <- median(joint$ecker_ch)
cat(sprintf("  Ecker CG median: %.6f, Ecker CH median: %.6f\n", cg_med, ch_med))

joint$cg_group <- ifelse(joint$ecker_cg >= cg_med, "CG High", "CG Low")
joint$ch_group <- ifelse(joint$ecker_ch >= ch_med, "CH High", "CH Low")
joint$joint_group <- paste(joint$cg_group, joint$ch_group, sep = " / ")

joint_rows <- list()
for (cg in c("CG Low", "CG High")) {
  for (ch in c("CH Low", "CH High")) {
    idx <- which(joint$cg_group == cg & joint$ch_group == ch)
    if (length(idx) < 5) next
    r <- joint$ols_residual[idx]
    wt <- wilcox.test(r, mu = 0)
    joint_rows[[paste(cg, ch)]] <- data.frame(
      cg_group = cg, ch_group = ch,
      n = length(idx),
      mean_residual = mean(r),
      median_residual = median(r),
      sd_residual = sd(r),
      wilcox_p_vs_zero = wt$p.value,
      mean_ecker_cg = mean(joint$ecker_cg[idx]),
      mean_ecker_ch = mean(joint$ecker_ch[idx]),
      stringsAsFactors = FALSE, row.names = NULL
    )
    cat(sprintf("  %s / %s: n=%d, median_resid=%.4f, p=%.2e\n",
                cg, ch, length(idx), median(r), wt$p.value))
  }
}

joint_df <- do.call(rbind, joint_rows)
joint_path <- file.path(TABLES_DIR, "58_noncg_joint_stratification.tsv")
write.table(joint_df, joint_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Joint stratification: %s\n\n", joint_path))

# =============================================================================
# FIGURE 58a: CORRECTED 57b — ECKER CH vs RESIDUAL (MeCP2-Up ONLY)
# =============================================================================

cat("--- 58a: Corrected 57b correlation (MeCP2-Up only) ---\n")

p58a <- ggplot(up_peaks, aes(x = ecker_ch, y = ols_residual)) +
  geom_point(aes(color = peak_group), alpha = 0.15, size = 0.4) +
  geom_smooth(method = "lm", color = "#D73027", se = TRUE, linewidth = 0.8) +
  scale_color_manual(values = c("Non-CG Candidate" = "#984EA3",
                                 "CG-Concordant" = "#E41A1C")) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("rho = %.4f\np = %.2e\nn = %d\n(MeCP2-Up only)",
                            up_cor$estimate, up_cor$p.value, nrow(up_peaks)),
           hjust = 1.1, vjust = 1.3, size = 3.5) +
  labs(title = "MeCP2 Residual vs Ecker Non-CG Methylation (MeCP2-Up Only)",
       subtitle = "Corrected from 57b: restricting to Up peaks removes Simpson's paradox from pooling Up+Down",
       x = "Ecker CH methylation (wildtype cerebellum)",
       y = "MeCP2 residual (CG-unexplained fold change)",
       color = "Peak Class") +
  theme_biomodal() +
  guides(color = guide_legend(override.aes = list(alpha = 0.8, size = 2)))

save_multiformat_ggplot(p58a, file.path(SEC58_DIR, "58a_corrected_residual_vs_ecker_ch"),
                        width = 10, height = 8)
cat("  58a saved.\n")

# =============================================================================
# FIGURE 58b: MEDIAN RESIDUAL BY ECKER CH QUARTILE — MeCP2-Up (KEY FIGURE)
# =============================================================================

cat("--- 58b: Dose-response — MeCP2-Up ---\n")

q_df <- up_ols$quartiles
q_df$quartile_label <- factor(q_df$quartile_label,
  levels = c("Q1 (low CH)", "Q2", "Q3", "Q4 (high CH)"))

jt_label <- sprintf("Jonckheere (increasing): p = %.2e\nKruskal-Wallis: p = %.2e",
                     up_ols$tests$p_value[up_ols$tests$test == "Jonckheere-Terpstra (increasing)"],
                     up_ols$tests$p_value[up_ols$tests$test == "Kruskal-Wallis"])

p58b <- ggplot(q_df, aes(x = quartile_label, y = median_residual, fill = quartile_label)) +
  geom_col(alpha = 0.85, width = 0.6) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Q1 (low CH)" = "#F7F7F7", "Q2" = "#CCCCFF",
                                "Q3" = "#9999CC", "Q4 (high CH)" = "#663399")) +
  geom_text(aes(label = sprintf("n=%d", n)), vjust = ifelse(q_df$median_residual >= 0, -1.5, 1.5),
            size = 3) +
  annotate("text", x = Inf, y = Inf, label = jt_label,
           hjust = 1.1, vjust = 1.3, size = 3.5) +
  labs(title = "MeCP2 Residual by Ecker CH Quartile (MeCP2-Up Peaks)",
       subtitle = "Does non-CG methylation level scale MeCP2 recruitment strength?",
       x = "Ecker CH methylation quartile", y = "Median MeCP2 residual (OLS)") +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p58b, file.path(SEC58_DIR, "58b_dose_response_up"),
                        width = 10, height = 7)
cat("  58b saved.\n")

# =============================================================================
# FIGURE 58c: NEGATIVE CONTROL — MeCP2-Down
# =============================================================================

cat("--- 58c: Dose-response — MeCP2-Down (negative control) ---\n")

if (!is.null(down_ols) && nrow(down_ols$quartiles) >= 4) {
  dq_df <- down_ols$quartiles
  dq_df$quartile_label <- factor(dq_df$quartile_label,
    levels = c("Q1 (low CH)", "Q2", "Q3", "Q4 (high CH)"))

  jt_down_label <- sprintf("Jonckheere (increasing): p = %.2e\nKruskal-Wallis: p = %.2e",
    down_ols$tests$p_value[down_ols$tests$test == "Jonckheere-Terpstra (increasing)"],
    down_ols$tests$p_value[down_ols$tests$test == "Kruskal-Wallis"])

  p58c <- ggplot(dq_df, aes(x = quartile_label, y = median_residual, fill = quartile_label)) +
    geom_col(alpha = 0.85, width = 0.6) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = c("Q1 (low CH)" = "#F7F7F7", "Q2" = "#CCCCCC",
                                  "Q3" = "#999999", "Q4 (high CH)" = "#555555")) +
    geom_text(aes(label = sprintf("n=%d", n)), vjust = ifelse(dq_df$median_residual >= 0, -1.5, 1.5),
              size = 3) +
    annotate("text", x = Inf, y = Inf, label = jt_down_label,
             hjust = 1.1, vjust = 1.3, size = 3.5) +
    labs(title = "MeCP2 Residual by Ecker CH Quartile (MeCP2-Down Peaks)",
         subtitle = "Negative control: non-CG dose-response not expected for lost peaks",
         x = "Ecker CH methylation quartile", y = "Median MeCP2 residual (OLS)") +
    theme_biomodal() +
    theme(legend.position = "none")

  save_multiformat_ggplot(p58c, file.path(SEC58_DIR, "58c_dose_response_down"),
                          width = 10, height = 7)
  cat("  58c saved.\n")
}

# =============================================================================
# FIGURE 58d: SCATTER — ECKER CH vs RESIDUAL COLORED BY QUARTILE
# =============================================================================

cat("--- 58d: Scatter colored by CH quartile ---\n")

up_scatter <- up_peaks[!is.na(up_peaks$ecker_ch) & !is.na(up_peaks$ols_residual), ]
ch_breaks <- quantile(up_scatter$ecker_ch, probs = c(0, 0.25, 0.5, 0.75, 1))
ch_breaks[1] <- ch_breaks[1] - 1e-10
up_scatter$ch_quartile <- factor(
  c("Q1 (low CH)", "Q2", "Q3", "Q4 (high CH)")[
    as.integer(cut(up_scatter$ecker_ch, breaks = ch_breaks, labels = 1:4))],
  levels = c("Q1 (low CH)", "Q2", "Q3", "Q4 (high CH)")
)

p58d <- ggplot(up_scatter, aes(x = ecker_ch, y = ols_residual)) +
  geom_point(aes(color = ch_quartile), alpha = 0.15, size = 0.4) +
  geom_smooth(method = "loess", color = "#D73027", se = TRUE, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = c("Q1 (low CH)" = "#BCBDDC", "Q2" = "#9E9AC8",
                                 "Q3" = "#807DBA", "Q4 (high CH)" = "#6A51A3")) +
  labs(title = "MeCP2 Residual vs Ecker CH (MeCP2-Up, by Quartile)",
       subtitle = "LOESS trend line shows local relationship",
       x = "Ecker CH methylation (wildtype cerebellum)",
       y = "MeCP2 residual (OLS)",
       color = "CH Quartile") +
  theme_biomodal() +
  guides(color = guide_legend(override.aes = list(alpha = 0.8, size = 2)))

save_multiformat_ggplot(p58d, file.path(SEC58_DIR, "58d_scatter_by_quartile"),
                        width = 10, height = 8)
cat("  58d saved.\n")

# =============================================================================
# FIGURE 58e: JOINT CG × CH STRATIFICATION HEATMAP
# =============================================================================

cat("--- 58e: Joint CG × CH heatmap ---\n")

joint_df$cg_group <- factor(joint_df$cg_group, levels = c("CG Low", "CG High"))
joint_df$ch_group <- factor(joint_df$ch_group, levels = c("CH Low", "CH High"))

p58e <- ggplot(joint_df, aes(x = ch_group, y = cg_group, fill = median_residual)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.3f\nn=%d\np=%.1e",
                                 median_residual, n, wilcox_p_vs_zero)),
            size = 4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "Median\nResidual") +
  labs(title = "Joint CG × CH Stratification (MeCP2-Up Peaks)",
       subtitle = "Does CH effect on MeCP2 residual persist after controlling for CG level?",
       x = "Ecker CH methylation (median split)",
       y = "Ecker CG methylation (median split)") +
  theme_biomodal() +
  theme(panel.grid = element_blank())

save_multiformat_ggplot(p58e, file.path(SEC58_DIR, "58e_joint_cg_ch_heatmap"),
                        width = 8, height = 7)
cat("  58e saved.\n")

# =============================================================================
# FIGURE 58f: OLS vs QUANTILE REGRESSION RESIDUAL COMPARISON
# =============================================================================

cat("--- 58f: OLS vs QR residual comparison ---\n")

compare_df <- peaks[valid & !is.na(peaks$ols_residual) & !is.na(peaks$qr_residual), ]

set.seed(42)
if (nrow(compare_df) > 20000) {
  compare_df <- compare_df[sample(nrow(compare_df), 20000), ]
}

p58f <- ggplot(compare_df, aes(x = ols_residual, y = qr_residual)) +
  geom_point(alpha = 0.05, size = 0.3, color = "grey40") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#D73027") +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("Spearman rho = %.4f\nn = %d (sampled)",
                            ols_qr_cor, nrow(compare_df)),
           hjust = 1.1, vjust = 1.3, size = 3.5) +
  labs(title = "OLS vs Quantile Regression Residuals",
       subtitle = "Quantile regression (tau=0.5) as robust alternative to OLS",
       x = "OLS residual", y = "Quantile regression residual") +
  theme_biomodal()

save_multiformat_ggplot(p58f, file.path(SEC58_DIR, "58f_ols_vs_qr_comparison"),
                        width = 8, height = 8)
cat("  58f saved.\n")

# =============================================================================
# FIGURE 58g: COMPOSITE PANEL
# =============================================================================

cat("--- 58g: Composite panel ---\n")

composite_parts <- list()
if (exists("p58a")) composite_parts$a <- p58a + labs(title = NULL, subtitle = "A. Corrected Residual vs CH (Up Only)")
if (exists("p58b")) composite_parts$b <- p58b + labs(title = NULL, subtitle = "B. Dose-Response (CH Quartiles)")
if (exists("p58e")) composite_parts$e <- p58e + labs(title = NULL, subtitle = "C. Joint CG × CH Heatmap")

if (length(composite_parts) >= 2) {
  p58g <- wrap_plots(composite_parts, nrow = 1) +
    plot_annotation(
      title = "Non-CG Methylation Dose-Response on MeCP2 Residual",
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    )
  save_multiformat_ggplot(p58g, file.path(SEC58_DIR, "58g_composite"),
                          width = 22, height = 7)
  cat("  58g saved.\n")
}

# =============================================================================
# COMPLETE
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 58 COMPLETE\n")
cat("================================================================================\n")
cat(sprintf("Output directory: %s\n", SEC58_DIR))
cat(sprintf("Figures: %d panels generated\n",
            length(list.files(SEC58_DIR, recursive = TRUE, pattern = "\\.png$"))))

cat(sprintf("\nTSV outputs in %s:\n", TABLES_DIR))
for (f in list.files(TABLES_DIR, pattern = "^58_")) {
  cat(sprintf("  %s\n", f))
}

cat(sprintf("\n=== KEY RESULTS ===\n"))
cat(sprintf("  Corrected correlation (MeCP2-Up only):\n"))
cat(sprintf("    Spearman rho = %.4f, p = %.2e\n", up_cor$estimate, up_cor$p.value))
cat(sprintf("  (Compare to pooled Up+Down: rho = %.4f, p = %.2e)\n\n",
            pooled_cor$estimate, pooled_cor$p.value))

cat("  MeCP2-Up dose-response (OLS residual by Ecker CH quartile):\n")
for (i in seq_len(nrow(up_ols$quartiles))) {
  cat(sprintf("    %s: median_resid = %.4f [%.4f, %.4f]\n",
              up_ols$quartiles$quartile_label[i],
              up_ols$quartiles$median_residual[i],
              up_ols$quartiles$ci_lo[i],
              up_ols$quartiles$ci_hi[i]))
}
cat(sprintf("    Jonckheere (increasing) p = %.2e\n",
            up_ols$tests$p_value[up_ols$tests$test == "Jonckheere-Terpstra (increasing)"]))

if (!is.null(down_ols)) {
  cat("\n  MeCP2-Down negative control:\n")
  for (i in seq_len(nrow(down_ols$quartiles))) {
    cat(sprintf("    %s: median_resid = %.4f\n",
                down_ols$quartiles$quartile_label[i],
                down_ols$quartiles$median_residual[i]))
  }
  cat(sprintf("    Jonckheere (increasing) p = %.2e\n",
              down_ols$tests$p_value[down_ols$tests$test == "Jonckheere-Terpstra (increasing)"]))
}

cat(sprintf("\n  Joint CG × CH stratification:\n"))
for (i in seq_len(nrow(joint_df))) {
  cat(sprintf("    %s / %s: median_resid = %.4f, n = %d\n",
              joint_df$cg_group[i], joint_df$ch_group[i],
              joint_df$median_residual[i], joint_df$n[i]))
}

cat(sprintf("\n  Quantile regression vs OLS: residual Spearman r = %.4f\n", ols_qr_cor))
cat("================================================================================\n")
