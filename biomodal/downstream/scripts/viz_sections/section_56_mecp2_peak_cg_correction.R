# biomodal/downstream/scripts/viz_sections/section_56_mecp2_peak_cg_correction.R
# Section 56: Peak-Level MeCP2 vs CG Methylation — Non-CG Binding Test
#
# For each MeCP2 DiffBind peak, extracts CG 5mC signal and tests whether
# MeCP2 enrichment exceeds what CG methylation change alone predicts.
# Peaks where MeCP2 goes UP but CG does NOT increase are candidate
# non-CG methylation binding sites.
#
# Input:  peaks/mecp2/MeCP2_annotated.txt (202,650 DiffBind peaks)
#         + merged CG mC BigWigs (METHYLATION_BIGWIGS from _shared_config.R)
# Usage:  cd downstream && Rscript scripts/viz_sections/section_56_mecp2_peak_cg_correction.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 56 CONFIGURATION
# =============================================================================

SEC56_DIR <- file.path(OUTPUT_DIR, "56_mecp2_peak_cg_correction")
dir.create(SEC56_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# LOAD MeCP2 DiffBind DATA
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 56: PEAK-LEVEL MeCP2 vs CG METHYLATION\n")
cat("================================================================================\n\n")

cat("--- Loading MeCP2 DiffBind results ---\n")
peaks <- load_diffbind_flex(MECP2_FILES$annotated, "MeCP2")

peak_gr <- GRanges(
  seqnames = peaks$Chr,
  ranges = IRanges(start = peaks$Start, end = peaks$End)
)
cat(sprintf("  Created GRanges for %d peaks\n\n", length(peak_gr)))

# =============================================================================
# EXTRACT CG 5mC SIGNAL AT EACH PEAK
# =============================================================================

cat("--- Extracting CG 5mC signal at MeCP2 peak coordinates ---\n")

extract_signal_at_peaks <- function(bw_path, peak_gr, label) {
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

cat("  CG mC ctrl:\n")
peaks$cg_mc_ctrl <- extract_signal_at_peaks(
  METHYLATION_BIGWIGS$cg_mc_ctrl, peak_gr, "CG_mC_ctrl")

cat("  CG mC mut:\n")
peaks$cg_mc_mut <- extract_signal_at_peaks(
  METHYLATION_BIGWIGS$cg_mc_mut, peak_gr, "CG_mC_mut")

all_nonzero <- c(peaks$cg_mc_ctrl[peaks$cg_mc_ctrl > 0],
                 peaks$cg_mc_mut[peaks$cg_mc_mut > 0])
cg_pseudo <- if (length(all_nonzero) > 0) {
  max(quantile(all_nonzero, 0.05, na.rm = TRUE), 1e-6)
} else { 1e-6 }

peaks$cg_mc_log2fc <- log2((peaks$cg_mc_mut + cg_pseudo) /
                            (peaks$cg_mc_ctrl + cg_pseudo))

cat(sprintf("\n  CG mC pseudocount: %.2e\n", cg_pseudo))
cat(sprintf("  CG mC log2FC: median=%.4f, mean=%.4f\n",
            median(peaks$cg_mc_log2fc, na.rm = TRUE),
            mean(peaks$cg_mc_log2fc, na.rm = TRUE)))
cat(sprintf("  Peaks with CG mC increase: %d (%.1f%%)\n",
            sum(peaks$cg_mc_log2fc > 0, na.rm = TRUE),
            100 * sum(peaks$cg_mc_log2fc > 0, na.rm = TRUE) / nrow(peaks)))

# =============================================================================
# PEAK-LEVEL REGRESSION
# =============================================================================

cat("\n--- Peak-level regression: MeCP2 Fold ~ CG mC log2FC ---\n")

valid <- !is.na(peaks$Fold) & !is.na(peaks$cg_mc_log2fc) &
         is.finite(peaks$cg_mc_log2fc)
cat(sprintf("  Valid peaks for regression: %d / %d\n", sum(valid), nrow(peaks)))

model <- lm(Fold ~ cg_mc_log2fc, data = peaks, subset = valid)
model_summary <- summary(model)

cat(sprintf("  R² = %.4f, adj-R² = %.4f\n",
            model_summary$r.squared, model_summary$adj.r.squared))
cat(sprintf("  Slope = %.4f, Intercept = %.4f\n",
            coef(model)[2], coef(model)[1]))

fstat <- model_summary$fstatistic
f_pval <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
cat(sprintf("  F = %.2f, p = %.2e\n", fstat[1], f_pval))

peaks$predicted_fold <- NA_real_
peaks$predicted_fold[valid] <- predict(model)
peaks$residual <- NA_real_
peaks$residual[valid] <- residuals(model)

# =============================================================================
# SAVE REGRESSION COEFFICIENTS
# =============================================================================

coef_df <- data.frame(
  term = names(coef(model)),
  estimate = as.numeric(coef(model)),
  std_error = model_summary$coefficients[, "Std. Error"],
  t_value = model_summary$coefficients[, "t value"],
  p_value = model_summary$coefficients[, "Pr(>|t|)"],
  r_squared = model_summary$r.squared,
  adj_r_squared = model_summary$adj.r.squared,
  f_statistic = fstat[1],
  f_pvalue = f_pval,
  n = sum(valid),
  stringsAsFactors = FALSE, row.names = NULL
)

reg_path <- file.path(TABLES_DIR, "56_mecp2_peak_regression.tsv")
write.table(coef_df, reg_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Regression table: %s\n\n", reg_path))

# =============================================================================
# PEAK CLASSIFICATION
# =============================================================================

cat("--- Peak classification ---\n")

peaks$mecp2_class <- "Not Significant"
sig <- peaks$FDR < 0.05
peaks$mecp2_class[sig & peaks$Fold > 0] <- "MeCP2 Up"
peaks$mecp2_class[sig & peaks$Fold < 0] <- "MeCP2 Down"

peaks$cg_class <- "CG Unchanged"
peaks$cg_class[peaks$cg_mc_log2fc > 0] <- "CG Increase"
peaks$cg_class[peaks$cg_mc_log2fc < 0] <- "CG Decrease"

up_peaks <- peaks[peaks$mecp2_class == "MeCP2 Up", ]
n_up <- nrow(up_peaks)

n_cg_concordant <- sum(up_peaks$cg_mc_log2fc > 0, na.rm = TRUE)
n_noncg_candidate <- sum(up_peaks$cg_mc_log2fc <= 0, na.rm = TRUE)

cat(sprintf("  MeCP2-Up peaks (FDR<0.05, Fold>0): %d\n", n_up))
cat(sprintf("    CG-Concordant (CG also increases):     %d (%.1f%%)\n",
            n_cg_concordant, 100 * n_cg_concordant / n_up))
cat(sprintf("    Non-CG Candidate (CG unchanged/down):  %d (%.1f%%)\n",
            n_noncg_candidate, 100 * n_noncg_candidate / n_up))

down_peaks <- peaks[peaks$mecp2_class == "MeCP2 Down", ]
n_down <- nrow(down_peaks)
cat(sprintf("\n  MeCP2-Down peaks: %d\n", n_down))
cat(sprintf("    CG decreases:  %d (%.1f%%)\n",
            sum(down_peaks$cg_mc_log2fc < 0, na.rm = TRUE),
            100 * sum(down_peaks$cg_mc_log2fc < 0, na.rm = TRUE) / n_down))
cat(sprintf("    CG increases:  %d (%.1f%%)\n",
            sum(down_peaks$cg_mc_log2fc > 0, na.rm = TRUE),
            100 * sum(down_peaks$cg_mc_log2fc > 0, na.rm = TRUE) / n_down))

# Residual statistics for MeCP2-up peaks
up_resid <- up_peaks$residual[!is.na(up_peaks$residual)]
if (length(up_resid) > 5) {
  wt <- wilcox.test(up_resid, mu = 0)
  cat(sprintf("\n  MeCP2-Up residual: mean=%.3f, median=%.3f, Wilcoxon p=%.2e\n",
              mean(up_resid), median(up_resid), wt$p.value))
}

# Save classification summary
class_summary <- data.frame(
  category = c("CG-Concordant", "Non-CG Candidate", "Total MeCP2-Up"),
  n_peaks = c(n_cg_concordant, n_noncg_candidate, n_up),
  pct = c(100 * n_cg_concordant / n_up, 100 * n_noncg_candidate / n_up, 100),
  mean_mecp2_fold = c(
    mean(up_peaks$Fold[up_peaks$cg_mc_log2fc > 0], na.rm = TRUE),
    mean(up_peaks$Fold[up_peaks$cg_mc_log2fc <= 0], na.rm = TRUE),
    mean(up_peaks$Fold, na.rm = TRUE)
  ),
  mean_cg_log2fc = c(
    mean(up_peaks$cg_mc_log2fc[up_peaks$cg_mc_log2fc > 0], na.rm = TRUE),
    mean(up_peaks$cg_mc_log2fc[up_peaks$cg_mc_log2fc <= 0], na.rm = TRUE),
    mean(up_peaks$cg_mc_log2fc, na.rm = TRUE)
  ),
  mean_residual = c(
    mean(up_peaks$residual[up_peaks$cg_mc_log2fc > 0], na.rm = TRUE),
    mean(up_peaks$residual[up_peaks$cg_mc_log2fc <= 0], na.rm = TRUE),
    mean(up_peaks$residual, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

class_path <- file.path(TABLES_DIR, "56_mecp2_peak_classification_summary.tsv")
write.table(class_summary, class_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  Classification summary: %s\n\n", class_path))

# =============================================================================
# CG DMR OVERLAP
# =============================================================================

cat("--- CG DMR overlap cross-tabulation ---\n")

peaks$dmr_overlap <- "No CG DMR"

if (!is.null(mc_dmr) && nrow(mc_dmr) > 0) {
  sig_hyper <- mc_dmr[mc_dmr$dmr_qvalue < 0.05 & mc_dmr$mod_fold_change > 0, ]
  sig_hypo  <- mc_dmr[mc_dmr$dmr_qvalue < 0.05 & mc_dmr$mod_fold_change < 0, ]

  if (nrow(sig_hyper) > 0) {
    hyper_gr <- GRanges(seqnames = sig_hyper$chr,
                        ranges = IRanges(start = sig_hyper$start, end = sig_hyper$end))
    peaks$dmr_overlap[countOverlaps(peak_gr, hyper_gr) > 0] <- "Hyper-DMR"
  }
  if (nrow(sig_hypo) > 0) {
    hypo_gr <- GRanges(seqnames = sig_hypo$chr,
                       ranges = IRanges(start = sig_hypo$start, end = sig_hypo$end))
    hypo_mask <- countOverlaps(peak_gr, hypo_gr) > 0
    peaks$dmr_overlap[hypo_mask & peaks$dmr_overlap == "No CG DMR"] <- "Hypo-DMR"
  }
}

crosstab <- table(peaks$mecp2_class, peaks$dmr_overlap)
cat("  Cross-tabulation:\n")
print(crosstab)

crosstab_df <- as.data.frame(crosstab)
colnames(crosstab_df) <- c("mecp2_class", "dmr_overlap", "count")
crosstab_path <- file.path(TABLES_DIR, "56_mecp2_peak_dmr_crosstab.tsv")
write.table(crosstab_df, crosstab_path, sep = "\t", quote = FALSE, row.names = FALSE)

if ("MeCP2 Up" %in% rownames(crosstab) && ncol(crosstab) >= 2) {
  up_dmr <- crosstab["MeCP2 Up", ]
  cat(sprintf("\n  MeCP2-Up peak DMR overlap:\n"))
  for (nm in names(up_dmr)) {
    cat(sprintf("    %s: %d (%.1f%%)\n", nm, up_dmr[nm],
                100 * up_dmr[nm] / sum(up_dmr)))
  }
}
cat("\n")

# =============================================================================
# SAVE PER-PEAK ANNOTATED TSV
# =============================================================================

cat("--- Saving per-peak annotated data ---\n")

out_cols <- c("Chr", "Start", "End", "Fold", "FDR", "Conc_mut", "Conc_ctrl")
if ("SYMBOL" %in% names(peaks)) out_cols <- c(out_cols, "SYMBOL")
if ("annotation" %in% names(peaks)) out_cols <- c(out_cols, "annotation")

peak_out <- peaks[, out_cols, drop = FALSE]
peak_out$cg_mc_ctrl <- peaks$cg_mc_ctrl
peak_out$cg_mc_mut <- peaks$cg_mc_mut
peak_out$cg_mc_log2fc <- peaks$cg_mc_log2fc
peak_out$predicted_fold <- peaks$predicted_fold
peak_out$residual <- peaks$residual
peak_out$mecp2_class <- peaks$mecp2_class
peak_out$cg_class <- peaks$cg_class
peak_out$dmr_overlap <- peaks$dmr_overlap

peak_path <- file.path(TABLES_DIR, "56_mecp2_peak_annotated.tsv")
write.table(peak_out, peak_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Per-peak data: %s (%d rows x %d cols)\n\n",
            peak_path, nrow(peak_out), ncol(peak_out)))

# Save non-CG candidates separately
noncg <- peak_out[peak_out$mecp2_class == "MeCP2 Up" &
                   peak_out$cg_mc_log2fc <= 0, ]
noncg <- noncg[order(noncg$FDR), ]
noncg_path <- file.path(TABLES_DIR, "56_mecp2_peak_noncg_candidates.tsv")
write.table(noncg, noncg_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Non-CG candidates: %s (%d peaks)\n\n", noncg_path, nrow(noncg)))

# =============================================================================
# FIGURE 56a: MeCP2 FOLD vs CG mC CHANGE SCATTER
# =============================================================================

cat("--- 56a: MeCP2 Fold vs CG mC change ---\n")

plot_df <- peaks[valid, ]
plot_df$significance <- factor(
  ifelse(plot_df$FDR < 0.05,
         ifelse(plot_df$Fold > 0, "MeCP2 Up", "MeCP2 Down"),
         "Not Significant"),
  levels = c("Not Significant", "MeCP2 Up", "MeCP2 Down")
)

r2_label <- sprintf("R² = %.4f\nSlope = %.2f\np = %.2e\nn = %d",
                     model_summary$r.squared, coef(model)[2], f_pval, sum(valid))

p56a <- ggplot(plot_df, aes(x = cg_mc_log2fc, y = Fold)) +
  geom_point(aes(color = significance), alpha = 0.05, size = 0.3) +
  geom_smooth(method = "lm", color = "#D73027", se = TRUE, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = c("Not Significant" = "grey80",
                                 "MeCP2 Up" = "#D95F02",
                                 "MeCP2 Down" = "#7570B3")) +
  annotate("text", x = Inf, y = Inf, label = r2_label,
           hjust = 1.1, vjust = 1.3, size = 3.5) +
  labs(title = "MeCP2 Fold Change vs CG 5mC Change at Peak Level",
       subtitle = "Each point = one MeCP2 DiffBind peak (n=202,650)",
       x = "CG 5mC log2FC (mutant/control)",
       y = "MeCP2 DiffBind Fold (positive = mutant enriched)",
       color = "MeCP2 Status") +
  theme_biomodal() +
  guides(color = guide_legend(override.aes = list(alpha = 0.8, size = 2)))

save_multiformat_ggplot(p56a, file.path(SEC56_DIR, "56a_mecp2_vs_cg_scatter"),
                        width = 10, height = 8)
cat("  56a saved.\n")

# =============================================================================
# FIGURE 56b: CLASSIFICATION BAR CHART (KEY FIGURE)
# =============================================================================

cat("--- 56b: MeCP2-Up peak classification ---\n")

bar_df <- data.frame(
  class = c("CG-Concordant\n(CG also increases)",
            "Non-CG Candidate\n(CG unchanged or decreases)"),
  count = c(n_cg_concordant, n_noncg_candidate),
  pct = c(100 * n_cg_concordant / n_up, 100 * n_noncg_candidate / n_up),
  stringsAsFactors = FALSE
)
bar_df$class <- factor(bar_df$class, levels = bar_df$class)

p56b <- ggplot(bar_df, aes(x = class, y = count, fill = class)) +
  geom_col(alpha = 0.85, width = 0.6) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", count, pct)),
            vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("#E41A1C", "#984EA3")) +
  labs(title = sprintf("Classification of %d MeCP2-Up Peaks (FDR<0.05)", n_up),
       subtitle = "Is the MeCP2 gain at each peak explained by CG methylation increase?",
       x = NULL, y = "Number of peaks") +
  theme_biomodal() +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, max(bar_df$count) * 1.15))

save_multiformat_ggplot(p56b, file.path(SEC56_DIR, "56b_peak_classification"),
                        width = 8, height = 7)
cat("  56b saved.\n")

# =============================================================================
# FIGURE 56c: CG mC LEVEL AT MeCP2 PEAK CATEGORIES
# =============================================================================

cat("--- 56c: CG mC levels by MeCP2 status ---\n")

level_df <- peaks[peaks$mecp2_class %in% c("MeCP2 Up", "MeCP2 Down"), ]

level_long <- rbind(
  data.frame(mecp2_class = level_df$mecp2_class,
             condition = "Control",
             cg_mc = level_df$cg_mc_ctrl,
             stringsAsFactors = FALSE),
  data.frame(mecp2_class = level_df$mecp2_class,
             condition = "Mutant",
             cg_mc = level_df$cg_mc_mut,
             stringsAsFactors = FALSE)
)
level_long$mecp2_class <- factor(level_long$mecp2_class,
                                  levels = c("MeCP2 Up", "MeCP2 Down"))
level_long$condition <- factor(level_long$condition, levels = c("Control", "Mutant"))

p_labels <- list()
for (cond in c("Control", "Mutant")) {
  up_cg <- level_long$cg_mc[level_long$mecp2_class == "MeCP2 Up" &
                             level_long$condition == cond]
  down_cg <- level_long$cg_mc[level_long$mecp2_class == "MeCP2 Down" &
                               level_long$condition == cond]
  if (length(up_cg) > 5 && length(down_cg) > 5) {
    wt <- wilcox.test(up_cg, down_cg)
    p_labels[[cond]] <- data.frame(
      condition = factor(cond, levels = c("Control", "Mutant")),
      label = sprintf("p = %.1e\nmedian: %.4f vs %.4f",
                       wt$p.value,
                       median(up_cg, na.rm = TRUE),
                       median(down_cg, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
    cat(sprintf("  %s — MeCP2-Up vs Down: median %.4f vs %.4f, p=%.2e\n",
                cond, median(up_cg, na.rm = TRUE),
                median(down_cg, na.rm = TRUE), wt$p.value))
  }
}
p_label_df <- do.call(rbind, p_labels)

p56c <- ggplot(level_long, aes(x = mecp2_class, y = cg_mc, fill = mecp2_class)) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.2, alpha = 0.9) +
  facet_wrap(~ condition) +
  scale_fill_manual(values = c("MeCP2 Up" = "#D95F02", "MeCP2 Down" = "#7570B3")) +
  geom_text(data = p_label_df,
            aes(x = 1.5, y = Inf, label = label),
            vjust = 1.3, size = 3.2, inherit.aes = FALSE) +
  labs(title = "CG 5mC at Significant MeCP2 Peaks: Control vs Mutant",
       subtitle = "MeCP2-Up peaks sit at lower-CG regions in both conditions",
       x = NULL, y = "CG 5mC mean signal") +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p56c, file.path(SEC56_DIR, "56c_cg_levels_by_mecp2"),
                        width = 10, height = 7)
cat("  56c saved.\n")

# =============================================================================
# FIGURE 56d: RESIDUAL DISTRIBUTION FOR SIGNIFICANT PEAKS
# =============================================================================

cat("--- 56d: Residual distributions ---\n")

sig_peaks <- peaks[peaks$mecp2_class != "Not Significant" & !is.na(peaks$residual), ]
sig_peaks$mecp2_class <- factor(sig_peaks$mecp2_class,
                                 levels = c("MeCP2 Up", "MeCP2 Down"))

up_r <- sig_peaks$residual[sig_peaks$mecp2_class == "MeCP2 Up"]
down_r <- sig_peaks$residual[sig_peaks$mecp2_class == "MeCP2 Down"]

stat_labels <- character(0)
if (length(up_r) > 5) {
  wt_up <- wilcox.test(up_r, mu = 0)
  stat_labels <- c(stat_labels,
    sprintf("Up vs 0: median=%.2f, p=%.2e", median(up_r), wt_up$p.value))
}
if (length(down_r) > 5) {
  wt_down <- wilcox.test(down_r, mu = 0)
  stat_labels <- c(stat_labels,
    sprintf("Down vs 0: median=%.2f, p=%.2e", median(down_r), wt_down$p.value))
}

p56d <- ggplot(sig_peaks, aes(x = residual, fill = mecp2_class)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("MeCP2 Up" = "#D95F02", "MeCP2 Down" = "#7570B3")) +
  annotate("text", x = Inf, y = Inf,
           label = paste(stat_labels, collapse = "\n"),
           hjust = 1.1, vjust = 1.3, size = 3.2) +
  labs(title = "MeCP2 Residual at Significant Peaks",
       subtitle = "Residual = MeCP2 fold change unexplained by CG methylation change",
       x = "Residual (actual Fold - CG-predicted Fold)", y = "Density",
       fill = "MeCP2 Direction") +
  theme_biomodal()

save_multiformat_ggplot(p56d, file.path(SEC56_DIR, "56d_residual_distributions"),
                        width = 10, height = 7)
cat("  56d saved.\n")

# =============================================================================
# FIGURE 56e: CG DMR CROSS-TABULATION
# =============================================================================

cat("--- 56e: CG DMR overlap ---\n")

ct_plot <- as.data.frame(crosstab)
colnames(ct_plot) <- c("mecp2_class", "dmr_overlap", "count")
ct_plot <- ct_plot[ct_plot$mecp2_class != "Not Significant", ]
ct_plot$mecp2_class <- factor(ct_plot$mecp2_class,
                               levels = c("MeCP2 Up", "MeCP2 Down"))
ct_plot$dmr_overlap <- factor(ct_plot$dmr_overlap,
                               levels = c("No CG DMR", "Hypo-DMR", "Hyper-DMR"))

if (nrow(ct_plot) > 0) {
  ct_totals <- aggregate(count ~ mecp2_class, data = ct_plot, FUN = sum)
  ct_plot <- merge(ct_plot, ct_totals, by = "mecp2_class", suffixes = c("", "_total"))
  ct_plot$pct <- 100 * ct_plot$count / ct_plot$count_total

  p56e <- ggplot(ct_plot, aes(x = mecp2_class, y = count, fill = dmr_overlap)) +
    geom_col(position = "dodge", alpha = 0.85) +
    geom_text(aes(label = sprintf("%d\n(%.0f%%)", count, pct)),
              position = position_dodge(0.9), vjust = -0.3, size = 3) +
    scale_fill_manual(values = c("No CG DMR" = "#4575B4",
                                  "Hypo-DMR" = "#74ADD1",
                                  "Hyper-DMR" = "#D73027")) +
    labs(title = "MeCP2 Peak Overlap with CG mC DMRs",
         subtitle = "Do MeCP2-Up peaks sit at genes with CG hypermethylation?",
         x = NULL, y = "Number of peaks", fill = "CG DMR Status") +
    theme_biomodal() +
    coord_cartesian(ylim = c(0, max(ct_plot$count) * 1.2))

  # Fisher's test: MeCP2 direction vs DMR overlap
  if ("MeCP2 Up" %in% rownames(crosstab) && "MeCP2 Down" %in% rownames(crosstab)) {
    sig_ct <- crosstab[c("MeCP2 Up", "MeCP2 Down"), ]
    if (ncol(sig_ct) >= 2) {
      ft <- fisher.test(sig_ct)
      cat(sprintf("  Fisher's test (MeCP2 direction x DMR overlap): p=%.2e\n", ft$p.value))
    }
  }

  save_multiformat_ggplot(p56e, file.path(SEC56_DIR, "56e_dmr_crosstab"),
                          width = 10, height = 7)
  cat("  56e saved.\n")
}

# =============================================================================
# FIGURE 56f: NON-CG CANDIDATE PEAK CHARACTERIZATION
# =============================================================================

cat("--- 56f: Non-CG candidate chromatin context ---\n")

up_for_chromatin <- peaks[peaks$mecp2_class == "MeCP2 Up" & !is.na(peaks$cg_mc_log2fc), ]
up_for_chromatin$peak_class <- ifelse(up_for_chromatin$cg_mc_log2fc > 0,
                                       "CG-Concordant", "Non-CG Candidate")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
tss <- promoters(txdb, upstream = 0, downstream = 1)

up_gr <- GRanges(seqnames = up_for_chromatin$Chr,
                 ranges = IRanges(start = up_for_chromatin$Start,
                                  end = up_for_chromatin$End))

dist_to_tss <- rep(NA_real_, length(up_gr))
nearest_idx <- nearest(up_gr, tss, ignore.strand = TRUE)
valid_nearest <- !is.na(nearest_idx)
if (any(valid_nearest)) {
  dist_to_tss[valid_nearest] <- distance(up_gr[valid_nearest],
                                          tss[nearest_idx[valid_nearest]])
}

load_bed_as_gr <- function(path) {
  bed <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  GRanges(seqnames = bed$V1, ranges = IRanges(start = bed$V2, end = bed$V3))
}

overlaps_df <- data.frame(
  H3K27ac_overlap  = countOverlaps(up_gr, load_bed_as_gr(CHIP_PEAK_FILES$h3k27ac)) > 0,
  H3K27me3_overlap = countOverlaps(up_gr, load_bed_as_gr(CHIP_PEAK_FILES$h3k27me3)) > 0,
  H3K4me1_overlap  = countOverlaps(up_gr, load_bed_as_gr(CHIP_PEAK_FILES$h3k4me1)) > 0,
  H3K4me3_overlap  = countOverlaps(up_gr, load_bed_as_gr(CHIP_PEAK_FILES$h3k4me3)) > 0,
  Bivalent_overlap = countOverlaps(up_gr, load_bed_as_gr(CHIP_PEAK_FILES$bivalent)) > 0
)

up_for_chromatin$chromatin_state <- classify_chromatin_state(
  overlaps_df, dist_to_tss, TSS_THRESHOLD)

chrom_summary <- up_for_chromatin %>%
  group_by(peak_class, chromatin_state) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(peak_class) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

chrom_summary$chromatin_state <- factor(chrom_summary$chromatin_state,
                                         levels = CHROMATIN_STATE_ORDER)

p56f <- ggplot(chrom_summary, aes(x = peak_class, y = pct, fill = chromatin_state)) +
  geom_col(position = "stack", alpha = 0.85) +
  scale_fill_manual(values = CHROMATIN_STATE_COLORS) +
  labs(title = "Chromatin Context of MeCP2-Up Peaks",
       subtitle = "Non-CG Candidate vs CG-Concordant peaks: different genomic contexts?",
       x = NULL, y = "Percentage of peaks", fill = "Chromatin State") +
  theme_biomodal()

save_multiformat_ggplot(p56f, file.path(SEC56_DIR, "56f_chromatin_context"),
                        width = 10, height = 7)
cat("  56f saved.\n")

chrom_path <- file.path(TABLES_DIR, "56_mecp2_peak_chromatin_summary.tsv")
write.table(chrom_summary, chrom_path, sep = "\t", quote = FALSE, row.names = FALSE)

# =============================================================================
# FIGURE 56g: COMPOSITE
# =============================================================================

cat("--- 56g: Composite panel ---\n")

composite_parts <- list()
if (exists("p56a")) composite_parts$a <- p56a + labs(title = NULL, subtitle = "A. MeCP2 vs CG Change")
if (exists("p56b")) composite_parts$b <- p56b + labs(title = NULL, subtitle = "B. Peak Classification")
if (exists("p56d")) composite_parts$d <- p56d + labs(title = NULL, subtitle = "C. Residual Distribution")

if (length(composite_parts) >= 2) {
  p56g <- wrap_plots(composite_parts, nrow = 1) +
    plot_annotation(
      title = "Peak-Level MeCP2 vs CG Methylation: Non-CG Binding Test",
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    )
  save_multiformat_ggplot(p56g, file.path(SEC56_DIR, "56g_composite"),
                          width = 22, height = 7)
  cat("  56g saved.\n")
}

# =============================================================================
# COMPLETE
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 56 COMPLETE\n")
cat("================================================================================\n")
cat(sprintf("Output directory: %s\n", SEC56_DIR))
cat(sprintf("Figures: %d panels generated\n",
            length(list.files(SEC56_DIR, recursive = TRUE, pattern = "\\.png$"))))

cat(sprintf("\nTSV outputs in %s:\n", TABLES_DIR))
for (f in list.files(TABLES_DIR, pattern = "^56_")) {
  cat(sprintf("  %s\n", f))
}

cat(sprintf("\n=== KEY RESULTS ===\n"))
cat(sprintf("  Regression R²: %.4f (slope=%.2f)\n",
            model_summary$r.squared, coef(model)[2]))
cat(sprintf("  MeCP2-Up peaks: %d total\n", n_up))
cat(sprintf("    CG-Concordant:    %d (%.1f%%) — CG increase may explain MeCP2 gain\n",
            n_cg_concordant, 100 * n_cg_concordant / n_up))
cat(sprintf("    Non-CG Candidate: %d (%.1f%%) — MeCP2 gain WITHOUT CG increase\n",
            n_noncg_candidate, 100 * n_noncg_candidate / n_up))
cat("================================================================================\n")
