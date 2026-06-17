# biomodal/downstream/scripts/viz_sections/section_61_stoichiometry_mechanism.R
# Section 61: Stoichiometry Test & Mechanism Differentiation
#
# Tests whether 5mC↑/5hmC↓ at BAP1-KO genes is driven by:
#   (A) DNMT3A dehydroxymethylase — direct conversion of 5hmC→5mC (stoichiometric)
#   (B) TET inhibition — blocked TET + ongoing DNMT3A de novo methylation (independent)
#
# Also tests whether chromatin accessibility (ATAC) mediates K27me3→delta_ratio,
# defines the neuronal gene subset from GO enrichment, and compares BAP1-KO vs TET-KO.
#
# Figures:
#   61a: Total methylation (mC+hmC) ctrl vs mut violin by gene group
#   61b: Δ5mC vs Δ5hmC stoichiometry scatter with slope test
#   61c: R² bar — chromatin accessibility models vs section 25 baseline models
#   61d: Mediation scatter panels (K27me3→ATAC, ATAC→delta_ratio, bootstrap CI)
#   61e: 3-way Venn (Neuronal ∩ Coordinated ∩ MeCP2-Up)
#   61f: BAP1-KO vs TET-KO delta_ratio scatter at neuronal genes
#   61g: Composite
#
# Input:
#   diffbind_gene_level_all_marks.tsv, 59_quadrant_master.tsv,
#   mecp2_coordinated_genes.tsv, enrichment_go_bp.tsv, tet_ko_gene_signal.tsv
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_61_stoichiometry_mechanism.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(ggVennDiagram)
  library(patchwork)
})

# =============================================================================
# SECTION 61 CONFIGURATION
# =============================================================================

cat("================================================================================\n")
cat("SECTION 61: STOICHIOMETRY TEST & MECHANISM DIFFERENTIATION\n")
cat("  DNMT3A dehydroxymethylase vs TET inhibition\n")
cat("================================================================================\n\n")

SEC61_DIR <- file.path(OUTPUT_DIR, "61_stoichiometry_mechanism")
dir.create(SEC61_DIR, recursive = TRUE, showWarnings = FALSE)

DIFFBIND_PATH   <- file.path(TABLES_DIR, "diffbind_gene_level_all_marks.tsv")
QUADRANT_PATH   <- file.path(TABLES_DIR, "59_quadrant_master.tsv")
COORDINATED_PATH <- file.path(TABLES_DIR, "mecp2_coordinated_genes.tsv")
ENRICHMENT_PATH <- file.path(TABLES_DIR, "enrichment_go_bp.tsv")
TETKO_PATH      <- file.path("data", "tet_ko_gene_signal.tsv")

NEURONAL_GO_PATTERN <- "synap|neuron|axon"

fmt_p <- function(p) {
  if (is.na(p) || !is.finite(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

sig_stars <- function(q) {
  if (is.na(q) || !is.finite(q)) return("")
  if (q < 0.001) return("***")
  if (q < 0.01)  return("**")
  if (q < 0.05)  return("*")
  return("ns")
}

save_plot <- function(p, name, w = 10, h = 7) {
  save_multiformat_ggplot(p,
    base_path = file.path(SEC61_DIR, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

# =============================================================================
# DATA LOADING
# =============================================================================

cat("--- Loading data ---\n\n")

input_files <- c(DIFFBIND_PATH, QUADRANT_PATH, COORDINATED_PATH, ENRICHMENT_PATH, TETKO_PATH)
for (f in input_files) {
  if (!file.exists(f)) stop(paste0("Required file not found: ", f))
  cat(sprintf("  [OK] %s\n", f))
}

db <- read.table(DIFFBIND_PATH, header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, quote = "")
cat(sprintf("  DiffBind all-marks: %d genes\n", nrow(db)))

qm <- read.table(QUADRANT_PATH, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, quote = "")
cat(sprintf("  Quadrant master:    %d genes\n", nrow(qm)))

coord <- read.table(COORDINATED_PATH, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, quote = "")
cat(sprintf("  Coordinated genes:  %d genes\n", nrow(coord)))

go_bp <- read.table(ENRICHMENT_PATH, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, quote = "")
cat(sprintf("  GO BP terms:        %d terms\n", nrow(go_bp)))

tetko <- read.table(TETKO_PATH, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, quote = "")
cat(sprintf("  TET-KO signal:      %d genes\n", nrow(tetko)))

# --- Build master table ---

mecp2_cols <- c("gene", "mecp2_mean_fold", "mecp2_nearest_fdr", "mecp2_n_peaks",
                "gb_log2fc", "gb_signal_class")
mecp2_cols <- intersect(mecp2_cols, colnames(qm))
qm_slim <- qm[, mecp2_cols, drop = FALSE]

df <- merge(db, qm_slim, by = "gene", all.x = TRUE)

df$mecp2_status <- dplyr::case_when(
  !is.na(df$mecp2_nearest_fdr) & df$mecp2_nearest_fdr < Q_THRESHOLD & df$mecp2_mean_fold > 0 ~ "MeCP2 Up",
  !is.na(df$mecp2_nearest_fdr) & df$mecp2_nearest_fdr < Q_THRESHOLD & df$mecp2_mean_fold < 0 ~ "MeCP2 Down",
  TRUE ~ "NS"
)
df$mecp2_status <- factor(df$mecp2_status, levels = c("MeCP2 Up", "MeCP2 Down", "NS"))

df$is_coordinated <- df$gene %in% coord$gene

cat(sprintf("\n  Merged table: %d genes\n", nrow(df)))
cat(sprintf("  MeCP2 Up: %d | Down: %d | NS: %d\n",
            sum(df$mecp2_status == "MeCP2 Up"),
            sum(df$mecp2_status == "MeCP2 Down"),
            sum(df$mecp2_status == "NS")))
cat(sprintf("  Coordinated: %d genes\n", sum(df$is_coordinated)))

# --- Extract neuronal gene set from GO ---

cat("\n--- Extracting neuronal gene set from GO BP enrichment ---\n")

neuronal_terms <- go_bp[grepl(NEURONAL_GO_PATTERN, go_bp$Description, ignore.case = TRUE), ]
cat(sprintf("  Found %d neuronal GO terms (matching '%s')\n",
            nrow(neuronal_terms), NEURONAL_GO_PATTERN))

if (nrow(neuronal_terms) > 0) {
  neuronal_genes <- unique(unlist(strsplit(neuronal_terms$geneID, "/")))
  cat(sprintf("  Union: %d unique neuronal genes\n", length(neuronal_genes)))
  cat(sprintf("  Top terms: %s\n",
              paste(head(neuronal_terms$Description, 5), collapse = "; ")))
} else {
  neuronal_genes <- character(0)
  cat("  WARNING: No neuronal GO terms found\n")
}

df$is_neuronal <- df$gene %in% neuronal_genes

cat(sprintf("  Neuronal genes in master table: %d / %d\n",
            sum(df$is_neuronal), nrow(df)))

# =============================================================================
# 61a: TOTAL METHYLATION CONSERVATION TEST
# =============================================================================

cat("\n--- 61a: Total methylation conservation test ---\n")

df$delta_total <- df$total_mut - df$total_ctrl

gene_groups <- list(
  "All genes"      = rep(TRUE, nrow(df)),
  "MeCP2 Up"       = df$mecp2_status == "MeCP2 Up",
  "MeCP2 Down"     = df$mecp2_status == "MeCP2 Down",
  "Coordinated"    = df$is_coordinated,
  "Neuronal"       = df$is_neuronal,
  "Neuronal + Coord" = df$is_neuronal & df$is_coordinated
)

stoich_results <- data.frame()
for (grp_name in names(gene_groups)) {
  mask <- gene_groups[[grp_name]] & !is.na(df$delta_total)
  vals <- df$delta_total[mask]
  n <- length(vals)
  if (n < 3) next

  wt <- wilcox.test(vals, mu = 0, conf.int = TRUE)
  stoich_results <- rbind(stoich_results, data.frame(
    group = grp_name,
    n = n,
    mean_total_ctrl = mean(df$total_ctrl[mask], na.rm = TRUE),
    mean_total_mut  = mean(df$total_mut[mask], na.rm = TRUE),
    mean_delta_total = mean(vals),
    median_delta_total = median(vals),
    wilcox_p = wt$p.value,
    wilcox_estimate = as.numeric(wt$estimate),
    ci_lo = wt$conf.int[1],
    ci_hi = wt$conf.int[2],
    stringsAsFactors = FALSE
  ))
}

stoich_results$wilcox_q <- p.adjust(stoich_results$wilcox_p, method = "BH")
stoich_results$sig <- sapply(stoich_results$wilcox_q, sig_stars)

cat("\n  Total methylation conservation (Wilcoxon vs 0, BH-corrected):\n")
for (i in seq_len(nrow(stoich_results))) {
  r <- stoich_results[i, ]
  direction <- ifelse(r$mean_delta_total > 0, "INCREASES",
                ifelse(r$mean_delta_total < 0, "DECREASES", "FLAT"))
  cat(sprintf("    %-20s n=%5d  mean_Δ=%+.5f  %s  q=%.2e %s\n",
              r$group, r$n, r$mean_delta_total, direction,
              r$wilcox_q, r$sig))
}

# Plot 61a: violin of delta_total by gene group
plot_groups <- c("All genes", "MeCP2 Up", "Coordinated", "Neuronal")
violin_data <- do.call(rbind, lapply(plot_groups, function(grp) {
  mask <- gene_groups[[grp]] & !is.na(df$delta_total)
  data.frame(group = grp, delta_total = df$delta_total[mask],
             stringsAsFactors = FALSE)
}))
violin_data$group <- factor(violin_data$group, levels = plot_groups)

p_61a <- ggplot(violin_data, aes(x = group, y = delta_total)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_violin(aes(fill = group), alpha = 0.6, show.legend = FALSE) +
  geom_boxplot(width = 0.15, outlier.size = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Total methylation (5mC + 5hmC) change in BAP1-KO",
    subtitle = "Δ(total) = total_mut - total_ctrl | dashed line = conservation (dehydroxymethylase)",
    x = NULL, y = "Δ Total methylation (mut − ctrl)"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

save_plot(p_61a, "61a_total_methylation_conservation", w = 10, h = 7)

# =============================================================================
# 61b: STOICHIOMETRY SCATTER (Δ5mC vs Δ5hmC)
# =============================================================================

cat("\n--- 61b: Stoichiometry scatter ---\n")

valid_stoich <- !is.na(df$mc_diff) & !is.na(df$hmc_diff)

slope_model <- lm(mc_diff ~ hmc_diff, data = df, subset = valid_stoich)
slope_summary <- summary(slope_model)
slope_ci <- confint(slope_model, "hmc_diff", level = 0.95)

cat(sprintf("  Global slope: %.4f [%.4f, %.4f]\n",
            coef(slope_model)["hmc_diff"], slope_ci[1], slope_ci[2]))
cat(sprintf("  R² = %.4f, %s\n", slope_summary$r.squared, fmt_p(slope_summary$coefficients["hmc_diff", 4])))
cat(sprintf("  Reference slope for dehydroxymethylase: -1.0\n"))
cat(sprintf("  Observed slope is %s -1: %s\n",
            ifelse(slope_ci[1] > -1 | slope_ci[2] < -1, "significantly different from",
                   "NOT significantly different from"),
            ifelse(slope_ci[1] > -1 | slope_ci[2] < -1, "→ NOT stoichiometric",
                   "→ consistent with stoichiometric conversion")))

stoich_slope_row <- data.frame(
  test = "Global slope (mc_diff ~ hmc_diff)",
  slope = coef(slope_model)["hmc_diff"],
  slope_ci_lo = slope_ci[1],
  slope_ci_hi = slope_ci[2],
  r_squared = slope_summary$r.squared,
  p_value = slope_summary$coefficients["hmc_diff", 4],
  reference_slope = -1,
  differs_from_ref = slope_ci[1] > -1 | slope_ci[2] < -1,
  stringsAsFactors = FALSE
)

# Per-group slopes
for (grp in c("MeCP2 Up", "Coordinated", "Neuronal")) {
  mask <- gene_groups[[grp]] & valid_stoich
  if (sum(mask) < 10) next
  m <- lm(mc_diff ~ hmc_diff, data = df[mask, ])
  ci <- confint(m, "hmc_diff", level = 0.95)
  cat(sprintf("  %-18s slope=%.4f [%.4f, %.4f]  R²=%.4f  n=%d\n",
              paste0(grp, ":"), coef(m)["hmc_diff"], ci[1], ci[2],
              summary(m)$r.squared, sum(mask)))
  stoich_slope_row <- rbind(stoich_slope_row, data.frame(
    test = paste0(grp, " slope"),
    slope = coef(m)["hmc_diff"], slope_ci_lo = ci[1], slope_ci_hi = ci[2],
    r_squared = summary(m)$r.squared,
    p_value = summary(m)$coefficients["hmc_diff", 4],
    reference_slope = -1,
    differs_from_ref = ci[1] > -1 | ci[2] < -1,
    stringsAsFactors = FALSE
  ))
}

# Scatter plot
scatter_df <- df[valid_stoich, c("gene", "mc_diff", "hmc_diff", "mecp2_status")]
status_colors <- c("MeCP2 Up" = "#D73027", "MeCP2 Down" = "#4575B4", "NS" = "grey80")

p_61b <- ggplot(scatter_df, aes(x = hmc_diff, y = mc_diff)) +
  geom_point(aes(color = mecp2_status), alpha = 0.4, size = 0.8) +
  geom_smooth(method = "lm", color = "black", linewidth = 1, se = TRUE, alpha = 0.2) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed",
              color = "#E69F00", linewidth = 1) +
  scale_color_manual(values = status_colors, name = "MeCP2 Status") +
  labs(
    title = "Stoichiometry test: Δ5mC vs Δ5hmC",
    subtitle = sprintf(
      "Black = OLS (slope=%.3f, R²=%.4f) | Orange dashed = reference slope −1 (dehydroxymethylase)",
      coef(slope_model)["hmc_diff"], slope_summary$r.squared
    ),
    x = "Δ5hmC (mod_difference)", y = "Δ5mC (mod_difference)"
  ) +
  coord_fixed(ratio = 1) +
  theme_biomodal()

save_plot(p_61b, "61b_stoichiometry_scatter", w = 9, h = 8)

# =============================================================================
# 61c: DELTA-RATIO CHROMATIN ACCESSIBILITY MODELS
# =============================================================================

cat("\n--- 61c: Delta-ratio chromatin accessibility models ---\n")

model_df <- df[complete.cases(df[, c("delta_ratio", "atac_fold", "k27me3_fold", "k119ub_fold")]), ]
cat(sprintf("  Complete cases for chromatin models: %d genes\n", nrow(model_df)))

lm_atac     <- lm(delta_ratio ~ atac_fold, data = model_df)
lm_atac_k27 <- lm(delta_ratio ~ atac_fold + k27me3_fold, data = model_df)
lm_full_chr <- lm(delta_ratio ~ atac_fold + k27me3_fold + k119ub_fold, data = model_df)
lm_full_all <- lm(delta_ratio ~ atac_fold + k27me3_fold + k119ub_fold +
                     k27ac_fold, data = model_df[complete.cases(model_df[, "k27ac_fold"]), ])

chromatin_models <- list(
  "ATAC only"                = lm_atac,
  "ATAC + K27me3"            = lm_atac_k27,
  "ATAC + K27me3 + K119ub"   = lm_full_chr,
  "ATAC + K27me3 + K119ub + K27ac" = lm_full_all
)

model_comparison <- data.frame()
for (nm in names(chromatin_models)) {
  s <- summary(chromatin_models[[nm]])
  f_p <- pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)
  model_comparison <- rbind(model_comparison, data.frame(
    model = nm,
    n = nrow(chromatin_models[[nm]]$model),
    r_squared = s$r.squared,
    adj_r_squared = s$adj.r.squared,
    f_stat = s$fstatistic[1],
    f_p = f_p,
    aic = AIC(chromatin_models[[nm]]),
    bic = BIC(chromatin_models[[nm]]),
    stringsAsFactors = FALSE
  ))
}

cat("\n  Chromatin accessibility models for delta_ratio:\n")
for (i in seq_len(nrow(model_comparison))) {
  r <- model_comparison[i, ]
  cat(sprintf("    %-35s R²=%.4f  adj_R²=%.4f  AIC=%.0f  %s\n",
              r$model, r$r_squared, r$adj_r_squared, r$aic, fmt_p(r$f_p)))
}

# Print per-predictor coefficients for the full chromatin model
cat("\n  Coefficients (ATAC + K27me3 + K119ub):\n")
full_coefs <- summary(lm_full_chr)$coefficients
for (j in seq_len(nrow(full_coefs))) {
  cat(sprintf("    %-18s beta=%+.5f  SE=%.5f  t=%.2f  %s\n",
              rownames(full_coefs)[j], full_coefs[j, 1], full_coefs[j, 2],
              full_coefs[j, 3], fmt_p(full_coefs[j, 4])))
}

# Plot 61c: R² comparison bar
model_comparison$model <- factor(model_comparison$model, levels = model_comparison$model)

p_61c <- ggplot(model_comparison, aes(x = model, y = r_squared)) +
  geom_col(fill = "#4393C3", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.4f", r_squared)), vjust = -0.5, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Chromatin accessibility models predicting TET efficiency change",
    subtitle = "Response: delta_ratio = hmC/(mC+hmC) change | Higher R² = more variance explained",
    x = NULL, y = expression(R^2)
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

save_plot(p_61c, "61c_chromatin_model_r2", w = 10, h = 7)

# =============================================================================
# 61d: MEDIATION ANALYSIS (K27me3 → ATAC → delta_ratio)
# =============================================================================

cat("\n--- 61d: Mediation analysis ---\n")

med_df <- df[complete.cases(df[, c("delta_ratio", "atac_fold", "k27me3_fold")]), ]
cat(sprintf("  Mediation sample: %d genes\n", nrow(med_df)))

path_a  <- lm(atac_fold ~ k27me3_fold, data = med_df)
path_b  <- lm(delta_ratio ~ atac_fold, data = med_df)
path_c  <- lm(delta_ratio ~ k27me3_fold, data = med_df)
path_cp <- lm(delta_ratio ~ k27me3_fold + atac_fold, data = med_df)

coef_a  <- coef(path_a)["k27me3_fold"]
coef_b  <- coef(path_cp)["atac_fold"]
coef_c  <- coef(path_c)["k27me3_fold"]
coef_cp <- coef(path_cp)["k27me3_fold"]
indirect_obs <- coef_a * coef_b
proportion_mediated <- indirect_obs / coef_c

cat(sprintf("  Path a (K27me3 → ATAC):        beta = %+.5f, %s\n",
            coef_a, fmt_p(summary(path_a)$coefficients["k27me3_fold", 4])))
cat(sprintf("  Path b (ATAC → delta_ratio|K27me3): beta = %+.5f, %s\n",
            coef_b, fmt_p(summary(path_cp)$coefficients["atac_fold", 4])))
cat(sprintf("  Path c  (total: K27me3 → delta_ratio):   beta = %+.5f\n", coef_c))
cat(sprintf("  Path c' (direct: K27me3 → delta_ratio|ATAC): beta = %+.5f\n", coef_cp))
cat(sprintf("  Indirect effect (a × b): %.6f\n", indirect_obs))
cat(sprintf("  Proportion mediated:     %.1f%%\n", 100 * proportion_mediated))

# Bootstrap indirect effect
set.seed(42)
B <- 5000
n_med <- nrow(med_df)
boot_indirect <- numeric(B)
for (b in seq_len(B)) {
  idx <- sample(n_med, replace = TRUE)
  boot_data <- med_df[idx, ]
  a_b <- coef(lm(atac_fold ~ k27me3_fold, data = boot_data))["k27me3_fold"]
  cp_b <- lm(delta_ratio ~ k27me3_fold + atac_fold, data = boot_data)
  b_b <- coef(cp_b)["atac_fold"]
  boot_indirect[b] <- a_b * b_b
}

boot_ci <- quantile(boot_indirect, c(0.025, 0.975))
boot_p <- 2 * min(mean(boot_indirect > 0), mean(boot_indirect < 0))

cat(sprintf("  Bootstrap 95%% CI: [%.6f, %.6f]\n", boot_ci[1], boot_ci[2]))
cat(sprintf("  Bootstrap p ≈ %.4f\n", boot_p))
cat(sprintf("  Mediation %s (CI %s zero)\n",
            ifelse(boot_p < 0.05, "SIGNIFICANT", "not significant"),
            ifelse(boot_ci[1] > 0 | boot_ci[2] < 0, "excludes", "includes")))

mediation_results <- data.frame(
  path = c("a (K27me3→ATAC)", "b (ATAC→DR|K27me3)", "c (total)", "c' (direct)", "indirect (a×b)"),
  coefficient = c(coef_a, coef_b, coef_c, coef_cp, indirect_obs),
  ci_lo = c(confint(path_a, "k27me3_fold")[1], confint(path_cp, "atac_fold")[1],
            confint(path_c, "k27me3_fold")[1], confint(path_cp, "k27me3_fold")[1],
            boot_ci[1]),
  ci_hi = c(confint(path_a, "k27me3_fold")[2], confint(path_cp, "atac_fold")[2],
            confint(path_c, "k27me3_fold")[2], confint(path_cp, "k27me3_fold")[2],
            boot_ci[2]),
  p_value = c(summary(path_a)$coefficients["k27me3_fold", 4],
              summary(path_cp)$coefficients["atac_fold", 4],
              summary(path_c)$coefficients["k27me3_fold", 4],
              summary(path_cp)$coefficients["k27me3_fold", 4],
              boot_p),
  proportion_mediated = c(NA, NA, NA, NA, proportion_mediated),
  stringsAsFactors = FALSE
)

# Plot 61d: 3-panel mediation
p_med1 <- ggplot(med_df, aes(x = k27me3_fold, y = atac_fold)) +
  geom_point(alpha = 0.2, size = 0.5, color = "grey50") +
  geom_smooth(method = "lm", color = "#D73027", linewidth = 1, se = TRUE) +
  labs(title = "Path a: K27me3 → ATAC",
       subtitle = sprintf("β = %.4f, %s", coef_a,
                          fmt_p(summary(path_a)$coefficients["k27me3_fold", 4])),
       x = "H3K27me3 log2FC", y = "ATAC-seq log2FC") +
  theme_biomodal(base_size = 10)

p_med2 <- ggplot(med_df, aes(x = atac_fold, y = delta_ratio)) +
  geom_point(alpha = 0.2, size = 0.5, color = "grey50") +
  geom_smooth(method = "lm", color = "#4575B4", linewidth = 1, se = TRUE) +
  labs(title = "Path b: ATAC → TET efficiency",
       subtitle = sprintf("β = %.4f, %s", coef_b,
                          fmt_p(summary(path_cp)$coefficients["atac_fold", 4])),
       x = "ATAC-seq log2FC", y = "Delta ratio (hmC/(mC+hmC) change)") +
  theme_biomodal(base_size = 10)

boot_df <- data.frame(indirect = boot_indirect)
p_med3 <- ggplot(boot_df, aes(x = indirect)) +
  geom_histogram(bins = 80, fill = "#4393C3", alpha = 0.7, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = boot_ci, linetype = "dotted", color = "#D73027", linewidth = 0.8) +
  geom_vline(xintercept = indirect_obs, color = "#D73027", linewidth = 1) +
  labs(title = "Indirect effect bootstrap (B=5000)",
       subtitle = sprintf("a×b = %.5f [%.5f, %.5f]  %s",
                          indirect_obs, boot_ci[1], boot_ci[2],
                          ifelse(boot_p < 0.05, "SIGNIFICANT", "ns")),
       x = "Indirect effect (a × b)", y = "Count") +
  theme_biomodal(base_size = 10)

p_61d <- p_med1 + p_med2 + p_med3 +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "Mediation: K27me3 → ATAC accessibility → TET efficiency",
    subtitle = sprintf("Proportion mediated: %.1f%% | If significant → TET inhibition via chromatin closure",
                       100 * proportion_mediated),
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

save_plot(p_61d, "61d_mediation_panels", w = 16, h = 6)

# =============================================================================
# 61e: NEURONAL GENE SUBSET DEFINITION + VENN
# =============================================================================

cat("\n--- 61e: Neuronal gene subset definition ---\n")

mecp2_up_genes <- df$gene[df$mecp2_status == "MeCP2 Up"]
coordinated_genes <- coord$gene
neuronal_gene_set <- neuronal_genes

cat(sprintf("  MeCP2-Up:     %d genes\n", length(mecp2_up_genes)))
cat(sprintf("  Coordinated:  %d genes\n", length(coordinated_genes)))
cat(sprintf("  Neuronal GO:  %d genes\n", length(neuronal_gene_set)))

# Compute overlaps
n_neur_coord <- length(intersect(neuronal_gene_set, coordinated_genes))
n_neur_up    <- length(intersect(neuronal_gene_set, mecp2_up_genes))
n_coord_up   <- length(intersect(coordinated_genes, mecp2_up_genes))
n_all_three  <- length(Reduce(intersect, list(neuronal_gene_set, coordinated_genes, mecp2_up_genes)))

cat(sprintf("  Neuronal ∩ Coordinated: %d\n", n_neur_coord))
cat(sprintf("  Neuronal ∩ MeCP2-Up:    %d\n", n_neur_up))
cat(sprintf("  Coordinated ∩ MeCP2-Up: %d\n", n_coord_up))
cat(sprintf("  All three:              %d\n", n_all_three))

# Venn diagram
venn_list <- list(
  "Neuronal GO" = neuronal_gene_set,
  "Coordinated\n(mC↑/hmC↓)" = coordinated_genes,
  "MeCP2 Up" = mecp2_up_genes
)

p_61e <- ggVennDiagram(venn_list, label_alpha = 0, edge_size = 0.5) +
  scale_fill_gradient(low = "white", high = "#4393C3") +
  labs(title = "Overlap: Neuronal × Coordinated × MeCP2-Up gene sets") +
  theme(legend.position = "none",
        plot.title = element_text(size = 13, face = "bold"))

save_plot(p_61e, "61e_neuronal_venn", w = 9, h = 8)

# Save neuronal gene set with flags
neuronal_df <- data.frame(
  gene = neuronal_gene_set,
  is_coordinated = neuronal_gene_set %in% coordinated_genes,
  is_mecp2_up = neuronal_gene_set %in% mecp2_up_genes,
  stringsAsFactors = FALSE
)
neuronal_df <- merge(neuronal_df, df[, c("gene", "mc_diff", "hmc_diff", "delta_ratio",
                                          "mecp2_mean_fold", "atac_fold",
                                          "k27me3_fold", "k119ub_fold")],
                     by = "gene", all.x = TRUE)

write.table(neuronal_df, file.path(TABLES_DIR, "61_neuronal_gene_set.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Saved: 61_neuronal_gene_set.tsv (%d genes)\n", nrow(neuronal_df)))

# --- Repeat stoichiometry on neuronal subset ---

cat("\n  Stoichiometry on neuronal subset:\n")
neur_mask <- df$is_neuronal & !is.na(df$mc_diff) & !is.na(df$hmc_diff)
if (sum(neur_mask) >= 10) {
  neur_slope <- lm(mc_diff ~ hmc_diff, data = df[neur_mask, ])
  neur_ci <- confint(neur_slope, "hmc_diff", level = 0.95)
  cat(sprintf("    Neuronal slope: %.4f [%.4f, %.4f]  R²=%.4f  n=%d\n",
              coef(neur_slope)["hmc_diff"], neur_ci[1], neur_ci[2],
              summary(neur_slope)$r.squared, sum(neur_mask)))
}

# --- Repeat mediation on neuronal subset ---

cat("  Mediation on neuronal subset:\n")
neur_med <- df[df$is_neuronal & complete.cases(df[, c("delta_ratio", "atac_fold", "k27me3_fold")]), ]
if (nrow(neur_med) >= 20) {
  a_n <- coef(lm(atac_fold ~ k27me3_fold, data = neur_med))["k27me3_fold"]
  cp_n <- lm(delta_ratio ~ k27me3_fold + atac_fold, data = neur_med)
  b_n <- coef(cp_n)["atac_fold"]
  indirect_n <- a_n * b_n

  set.seed(42)
  boot_n <- replicate(5000, {
    idx <- sample(nrow(neur_med), replace = TRUE)
    bd <- neur_med[idx, ]
    a_b <- coef(lm(atac_fold ~ k27me3_fold, data = bd))["k27me3_fold"]
    b_b <- coef(lm(delta_ratio ~ k27me3_fold + atac_fold, data = bd))["atac_fold"]
    a_b * b_b
  })
  boot_ci_n <- quantile(boot_n, c(0.025, 0.975))
  boot_p_n <- 2 * min(mean(boot_n > 0), mean(boot_n < 0))

  cat(sprintf("    indirect = %.6f [%.6f, %.6f]  p≈%.4f  n=%d\n",
              indirect_n, boot_ci_n[1], boot_ci_n[2], boot_p_n, nrow(neur_med)))
} else {
  cat(sprintf("    Too few neuronal genes with complete data: n=%d\n", nrow(neur_med)))
}

# =============================================================================
# 61f: TET-KO COMPARISON AT NEURONAL GENES
# =============================================================================

cat("\n--- 61f: TET-KO comparison ---\n")

tetko_merged <- merge(df[, c("gene", "delta_ratio", "mecp2_status", "is_neuronal",
                              "is_coordinated", "mc_diff", "hmc_diff")],
                      tetko[, c("gene", "delta_ratio_tet")],
                      by = "gene")
cat(sprintf("  Genes with both BAP1-KO and TET-KO data: %d\n", nrow(tetko_merged)))
cat(sprintf("  Of which neuronal: %d\n", sum(tetko_merged$is_neuronal)))

tetko_groups <- list(
  "All genes" = rep(TRUE, nrow(tetko_merged)),
  "Neuronal"  = tetko_merged$is_neuronal,
  "Coordinated" = tetko_merged$is_coordinated
)

cat("\n  BAP1-KO vs TET-KO delta_ratio Spearman correlations:\n")
tetko_cor_results <- data.frame()
for (grp in names(tetko_groups)) {
  mask <- tetko_groups[[grp]] &
    !is.na(tetko_merged$delta_ratio) & !is.na(tetko_merged$delta_ratio_tet)
  n <- sum(mask)
  if (n < 5) next

  ct <- cor.test(tetko_merged$delta_ratio[mask], tetko_merged$delta_ratio_tet[mask],
                 method = "spearman")
  cat(sprintf("    %-15s rho=%.4f  %s  n=%d\n", grp, ct$estimate, fmt_p(ct$p.value), n))

  tetko_cor_results <- rbind(tetko_cor_results, data.frame(
    group = grp, n = n, rho = as.numeric(ct$estimate),
    p_value = ct$p.value, stringsAsFactors = FALSE
  ))
}

# Plot 61f: BAP1-KO vs TET-KO scatter
plot_tetko <- tetko_merged[!is.na(tetko_merged$delta_ratio) &
                            !is.na(tetko_merged$delta_ratio_tet), ]
plot_tetko$highlight <- ifelse(plot_tetko$is_neuronal, "Neuronal", "Other")

global_rho <- cor.test(plot_tetko$delta_ratio, plot_tetko$delta_ratio_tet,
                       method = "spearman")

p_61f <- ggplot(plot_tetko, aes(x = delta_ratio_tet, y = delta_ratio)) +
  geom_point(aes(color = highlight), alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8, se = TRUE) +
  scale_color_manual(values = c("Neuronal" = "#D73027", "Other" = "grey70"),
                     name = "Gene class") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  labs(
    title = "BAP1-KO vs TET-KO: TET efficiency change",
    subtitle = sprintf("Spearman ρ = %.3f, %s, n = %d | If correlated → shared TET pathway",
                       global_rho$estimate, fmt_p(global_rho$p.value), nrow(plot_tetko)),
    x = "TET-KO Δ ratio (hmC/(mC+hmC))",
    y = "BAP1-KO Δ ratio (hmC/(mC+hmC))"
  ) +
  theme_biomodal()

save_plot(p_61f, "61f_bap1_vs_tetko_scatter", w = 9, h = 8)

# =============================================================================
# SAVE TABLES
# =============================================================================

cat("\n--- Saving tables ---\n")

write.table(stoich_results, file.path(TABLES_DIR, "61_stoichiometry_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 61_stoichiometry_summary.tsv\n")

write.table(stoich_slope_row, file.path(TABLES_DIR, "61_stoichiometry_slopes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 61_stoichiometry_slopes.tsv\n")

write.table(model_comparison, file.path(TABLES_DIR, "61_delta_ratio_chromatin_models.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 61_delta_ratio_chromatin_models.tsv\n")

write.table(mediation_results, file.path(TABLES_DIR, "61_mediation_results.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 61_mediation_results.tsv\n")

write.table(tetko_cor_results, file.path(TABLES_DIR, "61_tetko_comparison.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 61_tetko_comparison.tsv\n")

# =============================================================================
# COMPOSITE
# =============================================================================

cat("\n--- 61g: Composite ---\n")

p_61g <- (p_61a | p_61b) / (p_61c | p_61e) +
  plot_annotation(
    title = "Section 61: Stoichiometry & Mechanism Differentiation",
    subtitle = "DNMT3A dehydroxymethylase vs TET inhibition",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

save_plot(p_61g, "61g_composite", w = 18, h = 14)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 61 SUMMARY\n")
cat("================================================================================\n\n")

cat("  KEY RESULTS:\n")
cat(sprintf("    Total methylation (mC+hmC): %s in BAP1-KO (%s)\n",
            ifelse(stoich_results$mean_delta_total[1] > 0, "INCREASES", "CONSERVED"),
            ifelse(stoich_results$mean_delta_total[1] > 0,
                   "→ favors TET inhibition + de novo DNMT3A",
                   "→ favors dehydroxymethylase conversion")))
cat(sprintf("    Δ5mC~Δ5hmC slope: %.3f (reference: -1.0)\n",
            stoich_slope_row$slope[1]))
cat(sprintf("    Slope %s -1 → %s\n",
            ifelse(stoich_slope_row$differs_from_ref[1], "≠", "≈"),
            ifelse(stoich_slope_row$differs_from_ref[1],
                   "NOT stoichiometric → independent processes",
                   "consistent with stoichiometric conversion")))
cat(sprintf("    Mediation (K27me3→ATAC→delta_ratio): %s (proportion=%.1f%%)\n",
            ifelse(boot_p < 0.05, "SIGNIFICANT", "not significant"),
            100 * proportion_mediated))
cat(sprintf("    TET-KO phenocopy: rho=%.3f (%s)\n",
            global_rho$estimate,
            ifelse(global_rho$p.value < 0.05, "SIGNIFICANT", "ns")))

cat(sprintf("\n  Neuronal gene set: %d genes (saved to 61_neuronal_gene_set.tsv)\n",
            nrow(neuronal_df)))
cat(sprintf("    Overlap with coordinated: %d\n", n_neur_coord))
cat(sprintf("    Overlap with MeCP2-Up:    %d\n", n_neur_up))
cat(sprintf("    Triple overlap:           %d\n", n_all_three))

cat("\n  Plots saved to:", SEC61_DIR, "\n")
cat("    61a: Total methylation conservation violin\n")
cat("    61b: Stoichiometry scatter (Δ5mC vs Δ5hmC)\n")
cat("    61c: Chromatin accessibility model R² comparison\n")
cat("    61d: Mediation panels (K27me3→ATAC→delta_ratio)\n")
cat("    61e: Neuronal gene Venn diagram\n")
cat("    61f: BAP1-KO vs TET-KO delta_ratio scatter\n")
cat("    61g: Composite\n")

cat("\n  Tables saved to:", TABLES_DIR, "\n")
cat("    61_stoichiometry_summary.tsv\n")
cat("    61_stoichiometry_slopes.tsv\n")
cat("    61_delta_ratio_chromatin_models.tsv\n")
cat("    61_mediation_results.tsv\n")
cat("    61_neuronal_gene_set.tsv\n")
cat("    61_tetko_comparison.tsv\n")

cat("\nSection 61 complete.\n")
