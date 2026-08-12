# biomodal/downstream/scripts/viz_sections/section_57_ecker_noncg_validation.R
# Section 57: Ecker WGBS Validation of Non-CG MeCP2 Candidate Peaks
#
# Uses Joe Ecker's adult mouse cerebellum WGBS data (wildtype) to test whether
# the 2,726 MeCP2-Up peaks that lack CG methylation increases (non-CG candidates
# from section 56) sit at regions with higher non-CG (CH) methylation than the
# CG-concordant MeCP2-Up peaks.
#
# Input:  tables/56_mecp2_peak_annotated.tsv (from section 56)
#         + Ecker BigWigs (ECKER_BIGWIGS from _shared_config.R)
# Usage:  cd downstream && Rscript scripts/viz_sections/section_57_ecker_noncg_validation.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 57 CONFIGURATION
# =============================================================================

SEC57_DIR <- file.path(OUTPUT_DIR, "57_ecker_noncg_validation")
dir.create(SEC57_DIR, recursive = TRUE, showWarnings = FALSE)

PEAK_FILE <- file.path(TABLES_DIR, "56_mecp2_peak_annotated.tsv")

# =============================================================================
# LOAD DATA
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 57: ECKER WGBS NON-CG VALIDATION\n")
cat("================================================================================\n\n")

if (!file.exists(PEAK_FILE)) {
  stop("Section 56 peak data not found: ", PEAK_FILE,
       "\nRun section_56_mecp2_peak_cg_correction.R first.")
}

peaks <- read.table(PEAK_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                    fill = TRUE, comment.char = "")
cat(sprintf("Loaded %d peaks from section 56\n", nrow(peaks)))
cat(sprintf("  MeCP2 Up:  %d\n", sum(peaks$mecp2_class == "MeCP2 Up")))
cat(sprintf("  MeCP2 Down: %d\n", sum(peaks$mecp2_class == "MeCP2 Down")))
cat(sprintf("  Not Significant: %d\n", sum(peaks$mecp2_class == "Not Significant")))

up_peaks <- peaks[peaks$mecp2_class == "MeCP2 Up", ]
cat(sprintf("  CG-Concordant:   %d\n", sum(up_peaks$cg_class == "CG Increase")))
cat(sprintf("  Non-CG Candidate: %d\n\n",
            sum(up_peaks$cg_class %in% c("CG Decrease", "CG Unchanged"))))

for (nm in names(ECKER_BIGWIGS)) {
  if (!file.exists(ECKER_BIGWIGS[[nm]])) {
    stop("Ecker BigWig not found: ", ECKER_BIGWIGS[[nm]])
  }
}
cat("Ecker BigWig files validated.\n\n")

peak_gr <- GRanges(
  seqnames = peaks$Chr,
  ranges = IRanges(start = peaks$Start, end = peaks$End)
)

# =============================================================================
# EXTRACT ECKER SIGNAL AT EACH PEAK
# =============================================================================

extract_signal <- function(bw_path, peak_gr, label) {
  cat(sprintf("  %s (%d peaks)...", label, length(peak_gr)))

  bw <- rtracklayer::import.bw(bw_path, which = peak_gr)
  cov <- coverage(bw, weight = "score")

  ranges_by_chr <- split(ranges(peak_gr), seqnames(peak_gr))
  shared_chrs <- intersect(names(cov), names(ranges_by_chr))
  v <- Views(cov[shared_chrs], ranges_by_chr[shared_chrs])
  means <- viewMeans(v)

  result <- rep(NA_real_, length(peak_gr))
  chr_names <- as.character(seqnames(peak_gr))
  for (chr in names(means)) {
    idx <- which(chr_names == chr)
    result[idx] <- as.numeric(means[[chr]])
  }

  nonzero <- sum(result > 0, na.rm = TRUE)
  cat(sprintf(" median=%.6f, non-zero=%d/%d\n",
              median(result, na.rm = TRUE), nonzero, length(result)))
  result
}

cat("--- Extracting Ecker signals ---\n")
peaks$ecker_ch <- extract_signal(ECKER_BIGWIGS$ch, peak_gr, "Ecker CH (non-CG)")
peaks$ecker_cg <- extract_signal(ECKER_BIGWIGS$cg, peak_gr, "Ecker CG")

cat("\n--- Extracting evoC non-CG signals ---\n")
peaks$evoc_chg_ctrl <- extract_signal(METHYLATION_BIGWIGS$chg_mc_ctrl, peak_gr, "evoC CHG ctrl")
peaks$evoc_chg_mut  <- extract_signal(METHYLATION_BIGWIGS$chg_mc_mut, peak_gr, "evoC CHG mut")
peaks$evoc_chh_ctrl <- extract_signal(METHYLATION_BIGWIGS$chh_mc_ctrl, peak_gr, "evoC CHH ctrl")
peaks$evoc_chh_mut  <- extract_signal(METHYLATION_BIGWIGS$chh_mc_mut, peak_gr, "evoC CHH mut")

peaks$evoc_chg_mean <- (peaks$evoc_chg_ctrl + peaks$evoc_chg_mut) / 2
peaks$evoc_chh_mean <- (peaks$evoc_chh_ctrl + peaks$evoc_chh_mut) / 2
cat("\n")

# =============================================================================
# GROUP COMPARISONS
# =============================================================================

cat("--- Group comparisons ---\n\n")

peaks$peak_group <- peaks$mecp2_class
up_mask <- peaks$mecp2_class == "MeCP2 Up"
peaks$peak_group[up_mask & peaks$cg_class == "CG Increase"] <- "CG-Concordant"
peaks$peak_group[up_mask & peaks$cg_class != "CG Increase"] <- "Non-CG Candidate"

summary_rows <- list()
groups <- c("CG-Concordant", "Non-CG Candidate", "MeCP2 Down", "Not Significant")

for (grp in groups) {
  idx <- peaks$peak_group == grp
  if (sum(idx) < 5) next
  ch_vals <- peaks$ecker_ch[idx]
  cg_vals <- peaks$ecker_cg[idx]

  summary_rows[[grp]] <- data.frame(
    group = grp,
    n = sum(idx),
    median_ch = median(ch_vals, na.rm = TRUE),
    mean_ch = mean(ch_vals, na.rm = TRUE),
    pct_ch_detected = 100 * sum(ch_vals > 0, na.rm = TRUE) / sum(idx),
    median_cg = median(cg_vals, na.rm = TRUE),
    mean_cg = mean(cg_vals, na.rm = TRUE),
    stringsAsFactors = FALSE, row.names = NULL
  )
  cat(sprintf("  %s (n=%d): CH median=%.6f, CH detected=%.1f%%, CG median=%.4f\n",
              grp, sum(idx),
              median(ch_vals, na.rm = TRUE),
              100 * sum(ch_vals > 0, na.rm = TRUE) / sum(idx),
              median(cg_vals, na.rm = TRUE)))
}

summary_df <- do.call(rbind, summary_rows)

# Wilcoxon tests
noncg_ch <- peaks$ecker_ch[peaks$peak_group == "Non-CG Candidate"]
conc_ch  <- peaks$ecker_ch[peaks$peak_group == "CG-Concordant"]
down_ch  <- peaks$ecker_ch[peaks$peak_group == "MeCP2 Down"]

wt_main <- wilcox.test(noncg_ch, conc_ch)
cat(sprintf("\n  KEY TEST: Non-CG Candidate vs CG-Concordant (Ecker CH):\n"))
cat(sprintf("    median %.6f vs %.6f, W=%.0f, p=%.2e\n",
            median(noncg_ch, na.rm = TRUE), median(conc_ch, na.rm = TRUE),
            wt_main$statistic, wt_main$p.value))

noncg_cg <- peaks$ecker_cg[peaks$peak_group == "Non-CG Candidate"]
conc_cg  <- peaks$ecker_cg[peaks$peak_group == "CG-Concordant"]
wt_control <- wilcox.test(noncg_cg, conc_cg)
cat(sprintf("\n  CONTROL: Non-CG Candidate vs CG-Concordant (Ecker CG):\n"))
cat(sprintf("    median %.4f vs %.4f, W=%.0f, p=%.2e\n",
            median(noncg_cg, na.rm = TRUE), median(conc_cg, na.rm = TRUE),
            wt_control$statistic, wt_control$p.value))

summary_df$wilcox_ch_vs_concordant_p <- NA_real_
summary_df$wilcox_cg_vs_concordant_p <- NA_real_
for (grp in c("Non-CG Candidate", "MeCP2 Down")) {
  ch_v <- peaks$ecker_ch[peaks$peak_group == grp]
  cg_v <- peaks$ecker_cg[peaks$peak_group == grp]
  if (length(ch_v) >= 5) {
    summary_df$wilcox_ch_vs_concordant_p[summary_df$group == grp] <-
      wilcox.test(ch_v, conc_ch)$p.value
    summary_df$wilcox_cg_vs_concordant_p[summary_df$group == grp] <-
      wilcox.test(cg_v, conc_cg)$p.value
  }
}

# evoC non-CG comparisons
cat("\n--- evoC non-CG comparisons ---\n")
noncg_chg <- peaks$evoc_chg_mean[peaks$peak_group == "Non-CG Candidate"]
conc_chg  <- peaks$evoc_chg_mean[peaks$peak_group == "CG-Concordant"]
noncg_chh <- peaks$evoc_chh_mean[peaks$peak_group == "Non-CG Candidate"]
conc_chh  <- peaks$evoc_chh_mean[peaks$peak_group == "CG-Concordant"]

wt_evoc_chg <- wilcox.test(noncg_chg, conc_chg)
wt_evoc_chh <- wilcox.test(noncg_chh, conc_chh)

cat(sprintf("  evoC CHG: Non-CG Candidate median=%.6f vs CG-Concordant median=%.6f, p=%.2e\n",
            median(noncg_chg, na.rm = TRUE), median(conc_chg, na.rm = TRUE), wt_evoc_chg$p.value))
cat(sprintf("    Non-zero: %d/%d (%.1f%%) vs %d/%d (%.1f%%)\n",
            sum(noncg_chg > 0, na.rm = TRUE), sum(!is.na(noncg_chg)),
            100 * sum(noncg_chg > 0, na.rm = TRUE) / sum(!is.na(noncg_chg)),
            sum(conc_chg > 0, na.rm = TRUE), sum(!is.na(conc_chg)),
            100 * sum(conc_chg > 0, na.rm = TRUE) / sum(!is.na(conc_chg))))

cat(sprintf("  evoC CHH: Non-CG Candidate median=%.6f vs CG-Concordant median=%.6f, p=%.2e\n",
            median(noncg_chh, na.rm = TRUE), median(conc_chh, na.rm = TRUE), wt_evoc_chh$p.value))
cat(sprintf("    Non-zero: %d/%d (%.1f%%) vs %d/%d (%.1f%%)\n",
            sum(noncg_chh > 0, na.rm = TRUE), sum(!is.na(noncg_chh)),
            100 * sum(noncg_chh > 0, na.rm = TRUE) / sum(!is.na(noncg_chh)),
            sum(conc_chh > 0, na.rm = TRUE), sum(!is.na(conc_chh)),
            100 * sum(conc_chh > 0, na.rm = TRUE) / sum(!is.na(conc_chh))))

for (grp in groups) {
  idx <- peaks$peak_group == grp
  if (sum(idx) < 5) next
  summary_df$median_evoc_chg[summary_df$group == grp] <- median(peaks$evoc_chg_mean[idx], na.rm = TRUE)
  summary_df$median_evoc_chh[summary_df$group == grp] <- median(peaks$evoc_chh_mean[idx], na.rm = TRUE)
  summary_df$pct_evoc_chg_detected[summary_df$group == grp] <-
    100 * sum(peaks$evoc_chg_mean[idx] > 0, na.rm = TRUE) / sum(idx)
  summary_df$pct_evoc_chh_detected[summary_df$group == grp] <-
    100 * sum(peaks$evoc_chh_mean[idx] > 0, na.rm = TRUE) / sum(idx)
}

summary_path <- file.path(TABLES_DIR, "57_ecker_noncg_validation_summary.tsv")
write.table(summary_df, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  Summary table: %s\n\n", summary_path))

# Detection rate Fisher's test
noncg_detected <- sum(noncg_ch > 0, na.rm = TRUE)
noncg_total    <- sum(!is.na(noncg_ch))
conc_detected  <- sum(conc_ch > 0, na.rm = TRUE)
conc_total     <- sum(!is.na(conc_ch))

ft_mat <- matrix(c(noncg_detected, noncg_total - noncg_detected,
                    conc_detected, conc_total - conc_detected),
                 nrow = 2, dimnames = list(c("Detected", "Not"),
                                           c("Non-CG Cand", "CG-Conc")))
ft <- fisher.test(ft_mat)
cat(sprintf("  Fisher's (CH detection): OR=%.2f, p=%.2e\n\n", ft$estimate, ft$p.value))

# =============================================================================
# CORRELATION: RESIDUAL vs ECKER CH
# =============================================================================

cat("--- Correlation: MeCP2 residual vs Ecker CH ---\n")

sig_peaks <- peaks[peaks$mecp2_class != "Not Significant" &
                   !is.na(peaks$residual) & !is.na(peaks$ecker_ch), ]

if (nrow(sig_peaks) > 20) {
  cor_test <- cor.test(sig_peaks$residual, sig_peaks$ecker_ch, method = "spearman")
  cat(sprintf("  Spearman rho = %.4f, p = %.2e, n = %d\n",
              cor_test$estimate, cor_test$p.value, nrow(sig_peaks)))

  cor_df <- data.frame(
    metric = "residual_vs_ecker_ch",
    rho = as.numeric(cor_test$estimate),
    p_value = cor_test$p.value,
    n = nrow(sig_peaks),
    stringsAsFactors = FALSE
  )
  cor_path <- file.path(TABLES_DIR, "57_ecker_correlation.tsv")
  write.table(cor_df, cor_path, sep = "\t", quote = FALSE, row.names = FALSE)
}
cat("\n")

# =============================================================================
# SAVE PER-PEAK SIGNAL TSV
# =============================================================================

cat("--- Saving per-peak signal data ---\n")

out_df <- data.frame(
  chr = peaks$Chr, start = peaks$Start, end = peaks$End,
  mecp2_class = peaks$mecp2_class,
  peak_group = peaks$peak_group,
  cg_class = peaks$cg_class,
  Fold = peaks$Fold,
  FDR = peaks$FDR,
  residual = peaks$residual,
  ecker_ch = peaks$ecker_ch,
  ecker_cg = peaks$ecker_cg,
  evoc_chg_mean = peaks$evoc_chg_mean,
  evoc_chh_mean = peaks$evoc_chh_mean,
  stringsAsFactors = FALSE
)
if ("SYMBOL" %in% names(peaks)) out_df$SYMBOL <- peaks$SYMBOL

signal_path <- file.path(TABLES_DIR, "57_ecker_peak_signal.tsv")
write.table(out_df, signal_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Per-peak signal: %s (%d rows)\n\n", signal_path, nrow(out_df)))

# =============================================================================
# FIGURE 57a: KEY FIGURE — ECKER CH AT NON-CG vs CG-CONCORDANT
# =============================================================================

cat("--- 57a: Ecker CH at Non-CG Candidate vs CG-Concordant ---\n")

up_plot <- peaks[peaks$mecp2_class == "MeCP2 Up" & !is.na(peaks$ecker_ch), ]
up_plot$peak_group <- ifelse(up_plot$cg_class == "CG Increase",
                              "CG-Concordant", "Non-CG Candidate")
up_plot$peak_group <- factor(up_plot$peak_group,
                              levels = c("Non-CG Candidate", "CG-Concordant"))

stat_label <- sprintf("Wilcoxon p = %.2e", wt_main$p.value)

p57a <- ggplot(up_plot, aes(x = peak_group, y = ecker_ch, fill = peak_group)) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.2, alpha = 0.9) +
  scale_fill_manual(values = c("Non-CG Candidate" = "#984EA3",
                                "CG-Concordant" = "#E41A1C")) +
  annotate("text", x = 1.5, y = Inf, label = stat_label, vjust = 1.5, size = 3.5) +
  labs(title = "Ecker Non-CG (CH) Methylation at MeCP2-Up Peaks",
       subtitle = sprintf("Non-CG Candidates (n=%d) vs CG-Concordant (n=%d)\nEcker = wildtype adult cerebellum WGBS",
                           sum(up_plot$peak_group == "Non-CG Candidate"),
                           sum(up_plot$peak_group == "CG-Concordant")),
       x = NULL, y = "Ecker CH methylation (fractional)") +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p57a, file.path(SEC57_DIR, "57a_ecker_ch_noncg_vs_concordant"),
                        width = 8, height = 7)
cat("  57a saved.\n")

# =============================================================================
# FIGURE 57b: RESIDUAL vs ECKER CH CORRELATION
# =============================================================================

cat("--- 57b: Residual vs Ecker CH ---\n")

if (nrow(sig_peaks) > 20) {
  p57b <- ggplot(sig_peaks, aes(x = ecker_ch, y = residual)) +
    geom_point(aes(color = mecp2_class), alpha = 0.15, size = 0.4) +
    geom_smooth(method = "lm", color = "#D73027", se = TRUE, linewidth = 0.8) +
    scale_color_manual(values = c("MeCP2 Up" = "#D95F02", "MeCP2 Down" = "#7570B3")) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("rho = %.4f\np = %.2e\nn = %d",
                              cor_test$estimate, cor_test$p.value, nrow(sig_peaks)),
             hjust = 1.1, vjust = 1.3, size = 3.5) +
    labs(title = "MeCP2 Residual vs Ecker Non-CG Methylation",
         subtitle = "Does CG-unexplained MeCP2 enrichment correlate with non-CG methylation?",
         x = "Ecker CH methylation (wildtype cerebellum)",
         y = "MeCP2 residual (CG-unexplained fold change)",
         color = "MeCP2 Direction") +
    theme_biomodal() +
    guides(color = guide_legend(override.aes = list(alpha = 0.8, size = 2)))

  save_multiformat_ggplot(p57b, file.path(SEC57_DIR, "57b_residual_vs_ecker_ch"),
                          width = 10, height = 8)
  cat("  57b saved.\n")
}

# =============================================================================
# FIGURE 57c: DETECTION RATE BAR CHART
# =============================================================================

cat("--- 57c: Detection rate comparison ---\n")

det_groups <- c("Non-CG Candidate", "CG-Concordant", "MeCP2 Down")
det_data <- data.frame(
  group = character(0), detected = integer(0), total = integer(0),
  pct = numeric(0), stringsAsFactors = FALSE
)

for (grp in det_groups) {
  ch_v <- peaks$ecker_ch[peaks$peak_group == grp]
  ch_v <- ch_v[!is.na(ch_v)]
  det_data <- rbind(det_data, data.frame(
    group = grp,
    detected = sum(ch_v > 0),
    total = length(ch_v),
    pct = 100 * sum(ch_v > 0) / length(ch_v),
    stringsAsFactors = FALSE
  ))
}
det_data$group <- factor(det_data$group, levels = det_groups)

p57c <- ggplot(det_data, aes(x = group, y = pct, fill = group)) +
  geom_col(alpha = 0.85, width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", pct, detected, total)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Non-CG Candidate" = "#984EA3",
                                "CG-Concordant" = "#E41A1C",
                                "MeCP2 Down" = "#7570B3")) +
  labs(title = "Ecker CH Detection Rate at MeCP2 Peaks",
       subtitle = sprintf("Fisher's (Non-CG vs Concordant): OR=%.2f, p=%.2e",
                           ft$estimate, ft$p.value),
       x = NULL, y = "% peaks with Ecker CH > 0") +
  theme_biomodal() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, max(det_data$pct) * 1.2))

save_multiformat_ggplot(p57c, file.path(SEC57_DIR, "57c_detection_rate"),
                        width = 8, height = 7)
cat("  57c saved.\n")

# =============================================================================
# FIGURE 57d: CONTROL — ECKER CG AT NON-CG vs CG-CONCORDANT
# =============================================================================

cat("--- 57d: Control — Ecker CG comparison ---\n")

stat_label_cg <- sprintf("Wilcoxon p = %.2e", wt_control$p.value)

p57d <- ggplot(up_plot, aes(x = peak_group, y = ecker_cg, fill = peak_group)) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.2, alpha = 0.9) +
  scale_fill_manual(values = c("Non-CG Candidate" = "#984EA3",
                                "CG-Concordant" = "#E41A1C")) +
  annotate("text", x = 1.5, y = Inf, label = stat_label_cg, vjust = 1.5, size = 3.5) +
  labs(title = "Ecker CG Methylation at MeCP2-Up Peaks (Control)",
       subtitle = "Same comparison as 57a but for CG methylation\nIf NO difference here → CH difference is specific",
       x = NULL, y = "Ecker CG methylation (fractional)") +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p57d, file.path(SEC57_DIR, "57d_ecker_cg_control"),
                        width = 8, height = 7)
cat("  57d saved.\n")

# =============================================================================
# FIGURE 57e: THREE-WAY — ECKER CH AT UP vs DOWN vs NS
# =============================================================================

cat("--- 57e: Three-way Ecker CH comparison ---\n")

three_df <- peaks[!is.na(peaks$ecker_ch), ]
set.seed(42)
ns_idx <- which(three_df$mecp2_class == "Not Significant")
if (length(ns_idx) > 10000) {
  keep <- c(which(three_df$mecp2_class != "Not Significant"),
            sample(ns_idx, 10000))
  three_df <- three_df[keep, ]
}
three_df$mecp2_class <- factor(three_df$mecp2_class,
                                levels = c("MeCP2 Up", "MeCP2 Down", "Not Significant"))

wt_up_down <- wilcox.test(
  peaks$ecker_ch[peaks$mecp2_class == "MeCP2 Up"],
  peaks$ecker_ch[peaks$mecp2_class == "MeCP2 Down"])

p57e <- ggplot(three_df, aes(x = mecp2_class, y = ecker_ch, fill = mecp2_class)) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.2, alpha = 0.9) +
  scale_fill_manual(values = COLORS$mecp2) +
  annotate("text", x = 1.5, y = Inf,
           label = sprintf("Up vs Down: p = %.2e", wt_up_down$p.value),
           vjust = 1.5, size = 3.5) +
  labs(title = "Ecker Non-CG Methylation by MeCP2 Status",
       subtitle = "Do MeCP2-Up peaks sit at regions with more non-CG methylation?",
       x = NULL, y = "Ecker CH methylation (fractional)") +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p57e, file.path(SEC57_DIR, "57e_three_way_ecker_ch"),
                        width = 8, height = 7)
cat("  57e saved.\n")

# =============================================================================
# FIGURE 57f: SIDE-BY-SIDE — ECKER CH vs evoC CHG vs evoC CHH
# =============================================================================

cat("--- 57f: Side-by-side Ecker vs evoC non-CG ---\n")

sbs_data <- rbind(
  data.frame(
    dataset = "Ecker CH\n(WGBS)",
    peak_group = up_plot$peak_group,
    signal = up_plot$ecker_ch,
    stringsAsFactors = FALSE
  ),
  data.frame(
    dataset = "evoC CHG\n(DUET)",
    peak_group = up_plot$peak_group,
    signal = up_plot$evoc_chg_mean,
    stringsAsFactors = FALSE
  ),
  data.frame(
    dataset = "evoC CHH\n(DUET)",
    peak_group = up_plot$peak_group,
    signal = up_plot$evoc_chh_mean,
    stringsAsFactors = FALSE
  )
)
sbs_data$dataset <- factor(sbs_data$dataset,
  levels = c("Ecker CH\n(WGBS)", "evoC CHG\n(DUET)", "evoC CHH\n(DUET)"))

sbs_stats <- data.frame(
  dataset = c("Ecker CH\n(WGBS)", "evoC CHG\n(DUET)", "evoC CHH\n(DUET)"),
  p_label = sprintf("p = %.2e", c(wt_main$p.value, wt_evoc_chg$p.value, wt_evoc_chh$p.value)),
  stringsAsFactors = FALSE
)
sbs_stats$dataset <- factor(sbs_stats$dataset, levels = levels(sbs_data$dataset))

p57f <- ggplot(sbs_data, aes(x = peak_group, y = signal, fill = peak_group)) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.2, alpha = 0.9) +
  facet_wrap(~ dataset, scales = "free_y") +
  scale_fill_manual(values = c("Non-CG Candidate" = "#984EA3",
                                "CG-Concordant" = "#E41A1C")) +
  geom_text(data = sbs_stats, aes(x = 1.5, y = Inf, label = p_label),
            vjust = 1.5, size = 3, inherit.aes = FALSE) +
  labs(title = "Non-CG Methylation: Ecker WGBS vs evoC DUET",
       subtitle = "Same MeCP2-Up peaks compared across three non-CG datasets",
       x = NULL, y = "Non-CG methylation (fractional)") +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p57f, file.path(SEC57_DIR, "57f_ecker_vs_evoc_sidebyside"),
                        width = 14, height = 7)
cat("  57f saved.\n")

# =============================================================================
# FIGURE 57g: COMPOSITE
# =============================================================================

cat("--- 57g: Composite panel ---\n")

composite_parts <- list()
if (exists("p57a")) composite_parts$a <- p57a + labs(title = NULL, subtitle = "A. Ecker CH: Non-CG vs Concordant")
if (exists("p57b")) composite_parts$b <- p57b + labs(title = NULL, subtitle = "B. Residual vs Ecker CH")
if (exists("p57d")) composite_parts$d <- p57d + labs(title = NULL, subtitle = "C. Control: Ecker CG")

if (length(composite_parts) >= 2) {
  p57g <- wrap_plots(composite_parts, nrow = 1) +
    plot_annotation(
      title = "Ecker WGBS Validation: Non-CG Methylation at MeCP2 Candidate Peaks",
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    )
  save_multiformat_ggplot(p57g, file.path(SEC57_DIR, "57g_composite"),
                          width = 22, height = 7)
  cat("  57g saved.\n")
}

# =============================================================================
# COMPLETE
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 57 COMPLETE\n")
cat("================================================================================\n")
cat(sprintf("Output directory: %s\n", SEC57_DIR))

cat(sprintf("\nTSV outputs in %s:\n", TABLES_DIR))
for (f in list.files(TABLES_DIR, pattern = "^57_")) {
  cat(sprintf("  %s\n", f))
}

cat(sprintf("\n=== KEY RESULTS ===\n"))
cat(sprintf("  Non-CG Candidate vs CG-Concordant (Ecker CH):\n"))
cat(sprintf("    Median: %.6f vs %.6f\n",
            median(noncg_ch, na.rm = TRUE), median(conc_ch, na.rm = TRUE)))
cat(sprintf("    Wilcoxon p = %.2e\n", wt_main$p.value))
cat(sprintf("  Control (Ecker CG):\n"))
cat(sprintf("    Median: %.4f vs %.4f\n",
            median(noncg_cg, na.rm = TRUE), median(conc_cg, na.rm = TRUE)))
cat(sprintf("    Wilcoxon p = %.2e\n", wt_control$p.value))
if (exists("cor_test")) {
  cat(sprintf("  Residual-CH correlation: rho=%.4f, p=%.2e\n",
              cor_test$estimate, cor_test$p.value))
}
cat(sprintf("  CH detection: Non-CG=%.1f%%, Concordant=%.1f%% (Fisher p=%.2e)\n",
            100 * noncg_detected / noncg_total,
            100 * conc_detected / conc_total,
            ft$p.value))
cat("================================================================================\n")
