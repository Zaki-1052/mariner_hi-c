# biomodal/downstream/scripts/viz_sections/section_55_mecp2_cg_corrected_proxy.R
# Section 55: MeCP2 as CG-Corrected Proxy for Non-CG Methylation in TADs
#
# Uses OLS regression to decompose MeCP2 ChIP enrichment into CG-explained and
# residual (candidate non-CG) components. Tests whether the residual is
# organized within TAD domains (variance ratio), enriched at boundaries, and
# whether it persists in TADs without CG methylation changes.
#
# Model A: MeCP2 ~ CG 5mC         (primary)
# Model B: MeCP2 ~ CG 5mC + CG 5hmC  (sensitivity)
#
# Requires: add_mecp2_tad_signal.sb to have been run on HPC first
# Input:  data/tad_methylation_signal_late.tsv (augmented with MeCP2 columns)
#         data/tad_mecp2_binlevel_variance.tsv (pre-computed on HPC)
# Usage:  cd downstream && Rscript scripts/viz_sections/section_55_mecp2_cg_corrected_proxy.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 55 CONFIGURATION
# =============================================================================

SEC55_DIR <- file.path(OUTPUT_DIR, "55_mecp2_cg_corrected_proxy")
dir.create(SEC55_DIR, recursive = TRUE, showWarnings = FALSE)

TAD_SIGNAL_FILE <- file.path(BASE_DIR, "data/tad_methylation_signal_late.tsv")

# =============================================================================
# LOAD DATA
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 55: MeCP2 CG-CORRECTED PROXY FOR NON-CG METHYLATION\n")
cat("================================================================================\n\n")

if (!file.exists(TAD_SIGNAL_FILE)) {
  stop("TAD signal file not found: ", TAD_SIGNAL_FILE,
       "\nRun compute_tad_methylation.R then add_mecp2_tad_signal.R first.")
}

tad <- read.table(TAD_SIGNAL_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("Loaded %d TAD domains from %s\n", nrow(tad), TAD_SIGNAL_FILE))

required_cols <- c("mecp2_ctrl_mean", "mecp2_mut_mean",
                   "cg_mc_ctrl_mean", "cg_mc_mut_mean",
                   "cg_hmc_ctrl_mean", "cg_hmc_mut_mean",
                   "cg_mc_log2fc")
missing <- setdiff(required_cols, names(tad))
if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "),
       "\nRun add_mecp2_tad_signal.R to add MeCP2 columns.")
}

# Load pre-computed bin-level variance data (from HPC add_mecp2_tad_signal.R)
BINVAR_FILE <- file.path(BASE_DIR, "data/tad_mecp2_binlevel_variance.tsv")
if (!file.exists(BINVAR_FILE)) {
  stop("Bin-level variance file not found: ", BINVAR_FILE,
       "\nRun add_mecp2_tad_signal.sb on Expanse first, then rsync both TSVs.")
}
binvar <- read.table(BINVAR_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("Loaded bin-level variance data: %d rows x %d cols\n\n",
            nrow(binvar), ncol(binvar)))

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

bootstrap_variance_ratio <- function(domain_means, within_vars, B = 1000, seed = 42) {
  valid <- !is.na(domain_means) & !is.na(within_vars) & within_vars > 0
  means_v <- domain_means[valid]
  vars_v  <- within_vars[valid]
  n <- length(means_v)
  if (n < 10) return(c(estimate = NA, ci_lo = NA, ci_hi = NA))

  set.seed(seed)
  boot_ratios <- replicate(B, {
    idx <- sample(n, replace = TRUE)
    var(means_v[idx]) / mean(vars_v[idx])
  })

  c(estimate = var(means_v) / mean(vars_v),
    ci_lo = quantile(boot_ratios, 0.025),
    ci_hi = quantile(boot_ratios, 0.975))
}

# =============================================================================
# FIT REGRESSION MODELS (TAD-LEVEL)
# =============================================================================

cat("--- Fitting regression models ---\n\n")

valid_rows <- complete.cases(tad[, c("mecp2_ctrl_mean", "cg_mc_ctrl_mean",
                                      "cg_hmc_ctrl_mean")])
cat(sprintf("  Complete cases for modelling: %d / %d\n", sum(valid_rows), nrow(tad)))

# Model A: MeCP2 ~ CG 5mC
model_a <- lm(mecp2_ctrl_mean ~ cg_mc_ctrl_mean, data = tad, subset = valid_rows)
model_a_summary <- summary(model_a)
cat(sprintf("\n  Model A: MeCP2 ~ CG_mC\n"))
cat(sprintf("    R² = %.4f, adj-R² = %.4f\n",
            model_a_summary$r.squared, model_a_summary$adj.r.squared))
cat(sprintf("    Slope = %.4f, Intercept = %.4f\n",
            coef(model_a)[2], coef(model_a)[1]))
cat(sprintf("    F = %.2f, p = %.2e\n",
            model_a_summary$fstatistic[1],
            pf(model_a_summary$fstatistic[1],
               model_a_summary$fstatistic[2],
               model_a_summary$fstatistic[3],
               lower.tail = FALSE)))

# Model B: MeCP2 ~ CG 5mC + CG 5hmC
model_b <- lm(mecp2_ctrl_mean ~ cg_mc_ctrl_mean + cg_hmc_ctrl_mean,
              data = tad, subset = valid_rows)
model_b_summary <- summary(model_b)
cat(sprintf("\n  Model B: MeCP2 ~ CG_mC + CG_hmC\n"))
cat(sprintf("    R² = %.4f, adj-R² = %.4f\n",
            model_b_summary$r.squared, model_b_summary$adj.r.squared))
cat(sprintf("    Slope(mC) = %.4f, Slope(hmC) = %.4f, Intercept = %.4f\n",
            coef(model_b)[2], coef(model_b)[3], coef(model_b)[1]))
cat(sprintf("    F = %.2f, p = %.2e\n\n",
            model_b_summary$fstatistic[1],
            pf(model_b_summary$fstatistic[1],
               model_b_summary$fstatistic[2],
               model_b_summary$fstatistic[3],
               lower.tail = FALSE)))

# Compute per-TAD residuals for ctrl
tad$resid_a_ctrl <- NA_real_
tad$resid_a_ctrl[valid_rows] <- residuals(model_a)

tad$resid_b_ctrl <- NA_real_
tad$resid_b_ctrl[valid_rows] <- residuals(model_b)

# Apply ctrl model to mut data
tad$predicted_a_mut <- predict(model_a,
  newdata = data.frame(cg_mc_ctrl_mean = tad$cg_mc_mut_mean))
tad$resid_a_mut <- tad$mecp2_mut_mean - tad$predicted_a_mut

tad$predicted_b_mut <- predict(model_b,
  newdata = data.frame(cg_mc_ctrl_mean = tad$cg_mc_mut_mean,
                       cg_hmc_ctrl_mean = tad$cg_hmc_mut_mean))
tad$resid_b_mut <- tad$mecp2_mut_mean - tad$predicted_b_mut

tad$predicted_a_ctrl <- predict(model_a,
  newdata = data.frame(cg_mc_ctrl_mean = tad$cg_mc_ctrl_mean))
tad$predicted_b_ctrl <- predict(model_b,
  newdata = data.frame(cg_mc_ctrl_mean = tad$cg_mc_ctrl_mean,
                       cg_hmc_ctrl_mean = tad$cg_hmc_ctrl_mean))

# =============================================================================
# SAVE REGRESSION COEFFICIENTS TSV
# =============================================================================

cat("--- Saving regression coefficient tables ---\n")

coef_rows <- list()

extract_model_stats <- function(mod, mod_summary, model_name) {
  fstat <- mod_summary$fstatistic
  data.frame(
    model = model_name,
    formula = paste(deparse(formula(mod)), collapse = " "),
    n = nrow(mod$model),
    r_squared = mod_summary$r.squared,
    adj_r_squared = mod_summary$adj.r.squared,
    residual_se = mod_summary$sigma,
    f_statistic = fstat[1],
    f_df1 = fstat[2],
    f_df2 = fstat[3],
    f_pvalue = pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
}

coef_table <- rbind(
  extract_model_stats(model_a, model_a_summary, "Model_A"),
  extract_model_stats(model_b, model_b_summary, "Model_B")
)

coef_detail <- rbind(
  data.frame(
    model = "Model_A",
    term = names(coef(model_a)),
    estimate = as.numeric(coef(model_a)),
    std_error = model_a_summary$coefficients[, "Std. Error"],
    t_value = model_a_summary$coefficients[, "t value"],
    p_value = model_a_summary$coefficients[, "Pr(>|t|)"],
    stringsAsFactors = FALSE, row.names = NULL
  ),
  data.frame(
    model = "Model_B",
    term = names(coef(model_b)),
    estimate = as.numeric(coef(model_b)),
    std_error = model_b_summary$coefficients[, "Std. Error"],
    t_value = model_b_summary$coefficients[, "t value"],
    p_value = model_b_summary$coefficients[, "Pr(>|t|)"],
    stringsAsFactors = FALSE, row.names = NULL
  )
)

coef_path <- file.path(TABLES_DIR, "55_mecp2_regression_coefficients.tsv")
write.table(coef_detail, coef_path, sep = "\t", quote = FALSE, row.names = FALSE)

model_summary_path <- file.path(TABLES_DIR, "55_mecp2_regression_model_summary.tsv")
write.table(coef_table, model_summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Coefficients: %s\n", coef_path))
cat(sprintf("  Model summary: %s\n\n", model_summary_path))

# =============================================================================
# LOAD PRE-COMPUTED BIN-LEVEL VARIANCE DATA
# =============================================================================

cat("--- Loading pre-computed bin-level variance data ---\n")

within_var_raw_ctrl     <- binvar$within_var_raw_ctrl
within_var_raw_mut      <- binvar$within_var_raw_mut
within_var_resid_a_ctrl <- binvar$within_var_resid_a_ctrl
within_var_resid_a_mut  <- binvar$within_var_resid_a_mut
within_var_resid_b_ctrl <- binvar$within_var_resid_b_ctrl
within_var_resid_b_mut  <- binvar$within_var_resid_b_mut

cat(sprintf("  HPC model A R²: %.4f\n", binvar$model_a_r_squared[1]))
cat(sprintf("  HPC model B R²: %.4f\n", binvar$model_b_r_squared[1]))
cat(sprintf("  Non-NA raw ctrl variances: %d / %d\n",
            sum(!is.na(within_var_raw_ctrl)), length(within_var_raw_ctrl)))
cat(sprintf("  Non-NA resid A ctrl variances: %d / %d\n\n",
            sum(!is.na(within_var_resid_a_ctrl)), length(within_var_resid_a_ctrl)))

# =============================================================================
# VARIANCE RATIO COMPUTATION
# =============================================================================

cat("--- Variance ratios (bootstrap, B=1000) ---\n")

vr_rows <- list()

compute_and_report_vr <- function(label, domain_means, within_vars) {
  boot <- bootstrap_variance_ratio(domain_means, within_vars)
  cat(sprintf("  %s: %.4f [%.4f, %.4f]\n", label,
              boot["estimate"], boot["ci_lo"], boot["ci_hi"]))
  data.frame(
    metric = label,
    variance_ratio = boot["estimate"],
    ci_lo = boot["ci_lo"],
    ci_hi = boot["ci_hi"],
    n = sum(!is.na(domain_means) & !is.na(within_vars) & within_vars > 0),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

vr_rows[[1]] <- compute_and_report_vr("CG_5mC_ctrl",
  tad$cg_mc_ctrl_mean, tad$cg_mc_ctrl_within_var)
vr_rows[[2]] <- compute_and_report_vr("CG_5mC_mut",
  tad$cg_mc_mut_mean, tad$cg_mc_mut_within_var)
vr_rows[[3]] <- compute_and_report_vr("MeCP2_raw_ctrl",
  tad$mecp2_ctrl_mean, within_var_raw_ctrl)
vr_rows[[4]] <- compute_and_report_vr("MeCP2_raw_mut",
  tad$mecp2_mut_mean, within_var_raw_mut)
vr_rows[[5]] <- compute_and_report_vr("MeCP2_resid_A_ctrl",
  tad$resid_a_ctrl, within_var_resid_a_ctrl)
vr_rows[[6]] <- compute_and_report_vr("MeCP2_resid_A_mut",
  tad$resid_a_mut, within_var_resid_a_mut)
vr_rows[[7]] <- compute_and_report_vr("MeCP2_resid_B_ctrl",
  tad$resid_b_ctrl, within_var_resid_b_ctrl)
vr_rows[[8]] <- compute_and_report_vr("MeCP2_resid_B_mut",
  tad$resid_b_mut, within_var_resid_b_mut)

vr_df <- do.call(rbind, vr_rows)
vr_path <- file.path(TABLES_DIR, "55_mecp2_variance_ratio_summary.tsv")
write.table(vr_df, vr_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  Variance ratio table: %s\n\n", vr_path))

# =============================================================================
# CG-STRATIFIED RESIDUAL ANALYSIS
# =============================================================================

cat("--- CG-stratified residual analysis ---\n")

tad$cg_mc_quartile <- NA_integer_
valid_fc <- !is.na(tad$cg_mc_log2fc)
if (sum(valid_fc) > 0) {
  breaks <- quantile(tad$cg_mc_log2fc[valid_fc], probs = c(0, 0.25, 0.5, 0.75, 1))
  breaks[1] <- breaks[1] - 1
  tad$cg_mc_quartile[valid_fc] <- as.integer(
    cut(tad$cg_mc_log2fc[valid_fc], breaks = breaks, labels = 1:4))
}

strat_rows <- list()
for (q in 1:4) {
  idx <- which(tad$cg_mc_quartile == q)
  if (length(idx) < 5) next
  resid_vals <- tad$resid_a_ctrl[idx]
  resid_vals <- resid_vals[!is.na(resid_vals)]
  if (length(resid_vals) < 5) next

  wt <- wilcox.test(resid_vals, mu = 0, alternative = "two.sided")
  pk <- if ("mecp2_peak_count" %in% names(tad)) mean(tad$mecp2_peak_count[idx]) else NA

  strat_rows[[q]] <- data.frame(
    quartile = q,
    quartile_label = c("Q1 (hypo)", "Q2", "Q3", "Q4 (hyper)")[q],
    n_tads = length(idx),
    mean_cg_log2fc = mean(tad$cg_mc_log2fc[idx], na.rm = TRUE),
    mean_resid_a = mean(resid_vals),
    median_resid_a = median(resid_vals),
    sd_resid_a = sd(resid_vals),
    wilcox_p_vs_zero = wt$p.value,
    wilcox_statistic = as.numeric(wt$statistic),
    mean_peak_count = pk,
    stringsAsFactors = FALSE, row.names = NULL
  )
  cat(sprintf("  Q%d: n=%d, median_resid=%.4f, p=%.2e\n",
              q, length(idx), median(resid_vals), wt$p.value))
}

strat_df <- do.call(rbind, strat_rows)
strat_path <- file.path(TABLES_DIR, "55_mecp2_cg_stratified_summary.tsv")
write.table(strat_df, strat_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  CG-stratified table: %s\n\n", strat_path))

# =============================================================================
# DMR-STRATIFIED RESIDUAL ANALYSIS
# =============================================================================

cat("--- DMR-stratified residual analysis ---\n")

domain_gr <- GRanges(
  seqnames = tad$chr,
  ranges = IRanges(start = tad$start, end = tad$end)
)

tad$dmr_class <- "No DMR"

if (file.exists(DATA_PATHS$mc_dmr)) {
  mc_dmr <- load_dmr_bed(DATA_PATHS$mc_dmr)
  sig_hyper <- mc_dmr[mc_dmr$dmr_qvalue < 0.05 & mc_dmr$mod_fold_change > 0, ]
  sig_hypo  <- mc_dmr[mc_dmr$dmr_qvalue < 0.05 & mc_dmr$mod_fold_change < 0, ]

  if (nrow(sig_hyper) > 0) {
    hyper_gr <- GRanges(seqnames = sig_hyper$chr,
                        ranges = IRanges(start = sig_hyper$start, end = sig_hyper$end))
    hyper_olap <- countOverlaps(domain_gr, hyper_gr) > 0
    tad$dmr_class[hyper_olap] <- "Hyper-DMR"
  }
  if (nrow(sig_hypo) > 0) {
    hypo_gr <- GRanges(seqnames = sig_hypo$chr,
                       ranges = IRanges(start = sig_hypo$start, end = sig_hypo$end))
    hypo_olap <- countOverlaps(domain_gr, hypo_gr) > 0
    tad$dmr_class[hypo_olap & tad$dmr_class == "No DMR"] <- "Hypo-DMR"
    tad$dmr_class[hyper_olap & hypo_olap] <- "Both"
  }

  cat(sprintf("  DMR classes: %s\n",
              paste(sprintf("%s=%d", names(table(tad$dmr_class)),
                            as.integer(table(tad$dmr_class))), collapse = ", ")))
} else {
  cat("  WARNING: mc_dmr file not found, skipping DMR stratification\n")
}

dmr_rows <- list()
for (cls in unique(tad$dmr_class)) {
  idx <- which(tad$dmr_class == cls)
  if (length(idx) < 5) next
  resid_vals <- tad$resid_a_ctrl[idx]
  resid_vals <- resid_vals[!is.na(resid_vals)]
  if (length(resid_vals) < 5) next

  dmr_rows[[cls]] <- data.frame(
    dmr_class = cls,
    n_tads = length(idx),
    mean_resid_a = mean(resid_vals),
    median_resid_a = median(resid_vals),
    sd_resid_a = sd(resid_vals),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

dmr_strat_df <- do.call(rbind, dmr_rows)

if (nrow(dmr_strat_df) >= 2) {
  classes <- unique(tad$dmr_class)
  pw_tests <- list()
  for (i in seq_along(classes)) {
    for (j in seq(i + 1, length(classes))) {
      if (j > length(classes)) break
      v1 <- tad$resid_a_ctrl[tad$dmr_class == classes[i] & !is.na(tad$resid_a_ctrl)]
      v2 <- tad$resid_a_ctrl[tad$dmr_class == classes[j] & !is.na(tad$resid_a_ctrl)]
      if (length(v1) >= 5 && length(v2) >= 5) {
        wt <- wilcox.test(v1, v2)
        pw_tests[[length(pw_tests) + 1]] <- data.frame(
          group1 = classes[i], group2 = classes[j],
          wilcox_p = wt$p.value, stringsAsFactors = FALSE)
        cat(sprintf("  %s vs %s: W=%.0f, p=%.2e\n",
                    classes[i], classes[j], wt$statistic, wt$p.value))
      }
    }
  }
  if (length(pw_tests) > 0) {
    pw_df <- do.call(rbind, pw_tests)
    dmr_strat_df <- merge(dmr_strat_df,
      data.frame(dmr_class = character(0), stringsAsFactors = FALSE),
      all.x = TRUE)
  }
}

dmr_path <- file.path(TABLES_DIR, "55_mecp2_dmr_stratified_summary.tsv")
write.table(dmr_strat_df, dmr_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  DMR-stratified table: %s\n\n", dmr_path))

# =============================================================================
# SAVE PER-TAD RESIDUALS TSV
# =============================================================================

cat("--- Saving per-TAD residual data ---\n")

resid_out <- data.frame(
  chr = tad$chr, start = tad$start, end = tad$end,
  width_kb = tad$width_kb,
  is_differential = tad$is_differential,
  mecp2_ctrl_mean = tad$mecp2_ctrl_mean,
  mecp2_mut_mean = tad$mecp2_mut_mean,
  cg_mc_ctrl_mean = tad$cg_mc_ctrl_mean,
  cg_mc_mut_mean = tad$cg_mc_mut_mean,
  cg_hmc_ctrl_mean = tad$cg_hmc_ctrl_mean,
  cg_hmc_mut_mean = tad$cg_hmc_mut_mean,
  predicted_a_ctrl = tad$predicted_a_ctrl,
  predicted_a_mut = tad$predicted_a_mut,
  resid_a_ctrl = tad$resid_a_ctrl,
  resid_a_mut = tad$resid_a_mut,
  predicted_b_ctrl = tad$predicted_b_ctrl,
  predicted_b_mut = tad$predicted_b_mut,
  resid_b_ctrl = tad$resid_b_ctrl,
  resid_b_mut = tad$resid_b_mut,
  cg_mc_log2fc = tad$cg_mc_log2fc,
  cg_mc_quartile = tad$cg_mc_quartile,
  dmr_class = tad$dmr_class,
  mecp2_log2fc = tad$mecp2_log2fc,
  mean_tad_score = tad$mean_tad_score,
  stringsAsFactors = FALSE
)
if ("mecp2_peak_count" %in% names(tad)) {
  resid_out$mecp2_peak_count <- tad$mecp2_peak_count
  resid_out$mecp2_up_count   <- tad$mecp2_up_count
  resid_out$mecp2_down_count <- tad$mecp2_down_count
}

resid_path <- file.path(TABLES_DIR, "55_mecp2_tad_residuals.tsv")
write.table(resid_out, resid_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Per-TAD residuals: %s (%d rows x %d cols)\n\n",
            resid_path, nrow(resid_out), ncol(resid_out)))

# =============================================================================
# FIGURE 55a: MeCP2 vs CG 5mC SCATTER + REGRESSION
# =============================================================================

cat("--- 55a: MeCP2 vs CG 5mC scatter ---\n")

plot_df_a <- tad[valid_rows, c("cg_mc_ctrl_mean", "mecp2_ctrl_mean", "is_differential")]
plot_df_a$is_differential <- factor(plot_df_a$is_differential,
                                     levels = c(FALSE, TRUE),
                                     labels = c("Non-Differential", "Differential"))

r2_label <- sprintf("R² = %.4f\nSlope = %.2f\np = %.2e\nn = %d",
                     model_a_summary$r.squared, coef(model_a)[2],
                     pf(model_a_summary$fstatistic[1],
                        model_a_summary$fstatistic[2],
                        model_a_summary$fstatistic[3],
                        lower.tail = FALSE),
                     nrow(model_a$model))

p55a <- ggplot(plot_df_a, aes(x = cg_mc_ctrl_mean, y = mecp2_ctrl_mean)) +
  geom_point(aes(color = is_differential), alpha = 0.15, size = 0.4) +
  geom_smooth(method = "lm", color = "#D73027", se = TRUE, linewidth = 0.8) +
  scale_color_manual(values = c("Non-Differential" = "grey60",
                                 "Differential" = "#D73027")) +
  annotate("text", x = Inf, y = Inf, label = r2_label,
           hjust = 1.1, vjust = 1.3, size = 3.5) +
  labs(title = "MeCP2 ChIP Enrichment vs CG 5mC per TAD",
       subtitle = "Model A: how much MeCP2 is explained by CG methylation?",
       x = "CG 5mC mean signal (ctrl)", y = "MeCP2 mean signal (ctrl)",
       color = "TAD type") +
  theme_biomodal()

save_multiformat_ggplot(p55a, file.path(SEC55_DIR, "55a_mecp2_vs_cg_scatter"),
                        width = 10, height = 8)
cat("  55a saved.\n")

# =============================================================================
# FIGURE 55b: MODEL DIAGNOSTICS
# =============================================================================

cat("--- 55b: Model diagnostics ---\n")

diag_df <- data.frame(
  fitted = fitted(model_a),
  residuals = residuals(model_a),
  std_residuals = rstandard(model_a)
)

p55b_left <- ggplot(diag_df, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.1, size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#D73027") +
  geom_smooth(method = "loess", color = "#377EB8", se = FALSE, linewidth = 0.8) +
  labs(title = "Residuals vs Fitted", x = "Fitted values", y = "Residuals") +
  theme_biomodal()

p55b_right <- ggplot(diag_df, aes(sample = std_residuals)) +
  stat_qq(alpha = 0.1, size = 0.3) +
  stat_qq_line(color = "#D73027", linewidth = 0.8) +
  labs(title = "Normal Q-Q", x = "Theoretical quantiles",
       y = "Standardized residuals") +
  theme_biomodal()

p55b <- p55b_left + p55b_right +
  plot_annotation(title = "Model A Diagnostics (MeCP2 ~ CG 5mC)")

save_multiformat_ggplot(p55b, file.path(SEC55_DIR, "55b_model_diagnostics"),
                        width = 14, height = 6)
cat("  55b saved.\n")

# =============================================================================
# FIGURE 55c: RESIDUAL DISTRIBUTIONS (CTRL vs MUT)
# =============================================================================

cat("--- 55c: Residual distributions ---\n")

resid_plot <- rbind(
  data.frame(condition = "Control", residual = tad$resid_a_ctrl[!is.na(tad$resid_a_ctrl)]),
  data.frame(condition = "Mutant",  residual = tad$resid_a_mut[!is.na(tad$resid_a_mut)])
)

ctrl_resid <- tad$resid_a_ctrl[!is.na(tad$resid_a_ctrl)]
mut_resid  <- tad$resid_a_mut[!is.na(tad$resid_a_mut)]

shift_test <- wilcox.test(mut_resid, ctrl_resid)
mut_vs_zero <- wilcox.test(mut_resid, mu = 0)

stat_label <- sprintf("Mut vs Ctrl: p = %.2e\nMut vs 0: p = %.2e",
                       shift_test$p.value, mut_vs_zero$p.value)

p55c <- ggplot(resid_plot, aes(x = residual, fill = condition)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Control" = "#2166AC", "Mutant" = "#B2182B")) +
  annotate("text", x = Inf, y = Inf, label = stat_label,
           hjust = 1.1, vjust = 1.3, size = 3.5) +
  labs(title = "MeCP2 Residual Distribution (Model A)",
       subtitle = "Positive residual = MeCP2 enrichment beyond CG methylation prediction",
       x = "MeCP2 residual (actual - CG-predicted)", y = "Density",
       fill = "Condition") +
  theme_biomodal()

save_multiformat_ggplot(p55c, file.path(SEC55_DIR, "55c_residual_distributions"),
                        width = 10, height = 7)
cat("  55c saved.\n")

# =============================================================================
# FIGURE 55d: VARIANCE RATIO COMPARISON (KEY FIGURE)
# =============================================================================

cat("--- 55d: Variance ratio comparison ---\n")

vr_plot <- vr_df[grep("ctrl", vr_df$metric), ]
vr_plot$track <- factor(
  c("CG 5mC", "MeCP2 (raw)", "MeCP2 residual\n(Model A)", "MeCP2 residual\n(Model B)"),
  levels = c("CG 5mC", "MeCP2 (raw)", "MeCP2 residual\n(Model A)", "MeCP2 residual\n(Model B)")
)

p55d <- ggplot(vr_plot, aes(x = track, y = variance_ratio, fill = track)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("CG 5mC" = "#E41A1C",
                                "MeCP2 (raw)" = "#984EA3",
                                "MeCP2 residual\n(Model A)" = "#6A3D9A",
                                "MeCP2 residual\n(Model B)" = "#B15928")) +
  labs(title = "Inter-TAD / Intra-TAD Variance Ratio",
       subtitle = "Does the non-CG component of MeCP2 show TAD organization?\n(dashed line: ratio=1, no TAD organization; ctrl condition)",
       x = NULL, y = "Variance ratio (inter / intra)") +
  theme_biomodal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10))

save_multiformat_ggplot(p55d, file.path(SEC55_DIR, "55d_variance_ratio_comparison"),
                        width = 10, height = 7)
cat("  55d saved.\n")

# =============================================================================
# FIGURE 55e: BOUNDARY CORRELATION
# =============================================================================

cat("--- 55e: Boundary strength correlation ---\n")

boundary_df <- tad[!is.na(tad$resid_a_ctrl) & !is.na(tad$mean_tad_score), ]

if (nrow(boundary_df) > 20) {
  rho <- cor(boundary_df$mean_tad_score, boundary_df$resid_a_ctrl,
             method = "spearman", use = "complete.obs")
  rho_test <- cor.test(boundary_df$mean_tad_score, boundary_df$resid_a_ctrl,
                       method = "spearman")

  p55e <- ggplot(boundary_df, aes(x = mean_tad_score, y = resid_a_ctrl)) +
    geom_point(alpha = 0.15, size = 0.5) +
    geom_smooth(method = "lm", color = "#D73027", se = TRUE, linewidth = 0.8) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("rho = %.3f\np = %.2e\nn = %d",
                             rho, rho_test$p.value, nrow(boundary_df)),
             hjust = 1.1, vjust = 1.3, size = 3.5) +
    labs(title = "MeCP2 Residual vs TAD Boundary Strength",
         subtitle = "Does non-CG MeCP2 binding correlate with TAD boundary strength?",
         x = "Mean TAD boundary score (TADCompare)",
         y = "MeCP2 residual (Model A, ctrl)") +
    theme_biomodal()

  save_multiformat_ggplot(p55e, file.path(SEC55_DIR, "55e_boundary_correlation"),
                          width = 10, height = 8)
  cat(sprintf("  rho=%.3f, p=%.2e\n", rho, rho_test$p.value))
  cat("  55e saved.\n")
}

# =============================================================================
# FIGURE 55f: CG-STRATIFIED RESIDUAL BOXPLOT
# =============================================================================

cat("--- 55f: CG-stratified residuals ---\n")

strat_plot <- tad[!is.na(tad$cg_mc_quartile) & !is.na(tad$resid_a_ctrl), ]
strat_plot$quartile_label <- factor(
  c("Q1\n(hypo)", "Q2", "Q3", "Q4\n(hyper)")[strat_plot$cg_mc_quartile],
  levels = c("Q1\n(hypo)", "Q2", "Q3", "Q4\n(hyper)")
)

quartile_stats <- strat_df[, c("quartile_label", "wilcox_p_vs_zero")]
quartile_stats$label <- sprintf("p = %.2e", quartile_stats$wilcox_p_vs_zero)
quartile_stats$quartile_label <- factor(
  c("Q1\n(hypo)", "Q2", "Q3", "Q4\n(hyper)"),
  levels = c("Q1\n(hypo)", "Q2", "Q3", "Q4\n(hyper)")
)

p55f <- ggplot(strat_plot, aes(x = quartile_label, y = resid_a_ctrl,
                                fill = quartile_label)) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text(data = quartile_stats,
            aes(x = quartile_label, y = Inf, label = label),
            vjust = 1.5, size = 3, inherit.aes = FALSE) +
  scale_fill_brewer(palette = "RdYlBu", direction = -1) +
  labs(title = "MeCP2 Residual by CG Methylation Change Quartile",
       subtitle = "Key test: are residuals non-zero even at Q2-Q3 (no CG change)?",
       x = "CG 5mC log2FC quartile", y = "MeCP2 residual (Model A, ctrl)") +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p55f, file.path(SEC55_DIR, "55f_cg_stratified_residuals"),
                        width = 10, height = 7)
cat("  55f saved.\n")

# =============================================================================
# FIGURE 55g: DMR-STRATIFIED RESIDUAL
# =============================================================================

cat("--- 55g: DMR-stratified residuals ---\n")

dmr_plot <- tad[!is.na(tad$resid_a_ctrl), ]
dmr_plot$dmr_class <- factor(dmr_plot$dmr_class,
  levels = c("No DMR", "Hypo-DMR", "Hyper-DMR", "Both"))
dmr_plot <- dmr_plot[!is.na(dmr_plot$dmr_class), ]

if (length(unique(dmr_plot$dmr_class)) >= 2) {
  class_counts <- table(dmr_plot$dmr_class)
  class_labels <- sprintf("%s\n(n=%d)", names(class_counts), as.integer(class_counts))
  names(class_labels) <- names(class_counts)

  p55g <- ggplot(dmr_plot, aes(x = dmr_class, y = resid_a_ctrl, fill = dmr_class)) +
    geom_violin(alpha = 0.6, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = c("No DMR" = "#4575B4", "Hypo-DMR" = "#74ADD1",
                                  "Hyper-DMR" = "#D73027", "Both" = "#F46D43")) +
    scale_x_discrete(labels = class_labels) +
    labs(title = "MeCP2 Residual by CG DMR Overlap",
         subtitle = "Does MeCP2 excess differ at TADs with hypermethylated genes?",
         x = "TAD contains CG mC DMR gene (q<0.05)",
         y = "MeCP2 residual (Model A, ctrl)") +
    theme_biomodal() +
    theme(legend.position = "none")

  save_multiformat_ggplot(p55g, file.path(SEC55_DIR, "55g_dmr_stratified_residuals"),
                          width = 10, height = 7)
  cat("  55g saved.\n")
}

# =============================================================================
# FIGURE 55h: MODEL A vs MODEL B VARIANCE RATIO COMPARISON
# =============================================================================

cat("--- 55h: Model A vs B comparison ---\n")

vr_compare <- vr_df
vr_compare$condition <- ifelse(grepl("ctrl", vr_compare$metric), "Control", "Mutant")
vr_compare$track <- gsub("_(ctrl|mut)", "", vr_compare$metric)
vr_compare$track <- factor(vr_compare$track,
  levels = c("CG_5mC", "MeCP2_raw", "MeCP2_resid_A", "MeCP2_resid_B"),
  labels = c("CG 5mC", "MeCP2\n(raw)", "Residual\n(Model A)", "Residual\n(Model B)"))

p55h <- ggplot(vr_compare, aes(x = track, y = variance_ratio, fill = condition)) +
  geom_col(alpha = 0.8, position = position_dodge(0.7), width = 0.6) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.2, position = position_dodge(0.7)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Control" = "#2166AC", "Mutant" = "#B2182B")) +
  labs(title = "TAD Variance Ratios: Model A vs Model B",
       subtitle = "Does including 5hmC in the correction change the result?",
       x = NULL, y = "Variance ratio (inter / intra)",
       fill = "Condition") +
  theme_biomodal() +
  theme(axis.text.x = element_text(size = 10))

save_multiformat_ggplot(p55h, file.path(SEC55_DIR, "55h_model_comparison"),
                        width = 12, height = 7)
cat("  55h saved.\n")

# =============================================================================
# FIGURE 55i: COMPOSITE PANEL
# =============================================================================

cat("--- 55i: Composite panel ---\n")

composite_parts <- list()
if (exists("p55a")) composite_parts$a <- p55a + labs(title = NULL, subtitle = "A. MeCP2 vs CG 5mC")
if (exists("p55d")) composite_parts$d <- p55d + labs(title = NULL, subtitle = "B. Variance Ratio")
if (exists("p55f")) composite_parts$f <- p55f + labs(title = NULL, subtitle = "C. CG-Stratified Residuals")

if (length(composite_parts) >= 2) {
  p55i <- wrap_plots(composite_parts, nrow = 1) +
    plot_annotation(
      title = "MeCP2 as CG-Corrected Non-CG Methylation Proxy",
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    )
  save_multiformat_ggplot(p55i, file.path(SEC55_DIR, "55i_composite"),
                          width = 20, height = 7)
  cat("  55i saved.\n")
}

# =============================================================================
# COMPLETE
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 55 COMPLETE\n")
cat("================================================================================\n")
cat(sprintf("Output directory: %s\n", SEC55_DIR))
cat(sprintf("Figures: %d panels generated\n",
            length(list.files(SEC55_DIR, recursive = TRUE, pattern = "\\.png$"))))
cat(sprintf("\nTSV outputs in %s:\n", TABLES_DIR))
for (f in list.files(TABLES_DIR, pattern = "^55_")) {
  cat(sprintf("  %s\n", f))
}
cat("\nKey results:\n")
cat(sprintf("  Model A R²: %.4f\n", model_a_summary$r.squared))
cat(sprintf("  Model B R²: %.4f\n", model_b_summary$r.squared))
for (i in seq_len(nrow(vr_df))) {
  cat(sprintf("  Variance ratio %s: %.4f [%.4f, %.4f]\n",
              vr_df$metric[i], vr_df$variance_ratio[i],
              vr_df$ci_lo[i], vr_df$ci_hi[i]))
}
cat("================================================================================\n")
