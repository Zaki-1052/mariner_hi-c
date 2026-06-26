# biomodal/downstream/scripts/viz_sections/section_90_5hmc_subcompartment_gradient_mechanism.R
# Section 90: 5hmC Subcompartment Analysis + 5mC/5hmC Gradient Comparison
#              + Continuous PC1 Regression + TET vs DNMT3A Mechanism Differentiation
#
# Extends Section 66 (which only plots 5mC DMR fractions by subcompartment)
# to analyze 5hmC DMRs by subcompartment, compare gradient steepness between
# the two modalities, overlay continuous Homer PC1 eigenvalues, and use
# mediation analysis to differentiate TET impediment from DNMT3A recruitment.
#
#   Part 1 (90a-b): 5hmC DMR fractions and directions by subcompartment
#   Part 2 (90c-d): 5mC vs 5hmC gradient comparison
#   Part 3 (90e-f): Continuous PC1 regression and slope comparison
#   Part 4 (90g-h): Mediation analysis and variance partition
#   Part 5 (90i):   B-compartment behavior
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_90_5hmc_subcompartment_gradient_mechanism.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 90: 5hmC SUBCOMPARTMENT + MECHANISM DIFFERENTIATION
# =============================================================================

cat("================================================================================\n")
cat("SECTION 90: 5hmC SUBCOMPARTMENT + MECHANISM DIFFERENTIATION\n")
cat("================================================================================\n\n")

# ---- Constants --------------------------------------------------------------

SUBCMPT_ORDER <- c("A.1", "A.2", "B.1", "B.2")
SUBCMPT_LABELS <- c("A.1" = "A.1\n(Strong Active)",
                     "A.2" = "A.2\n(Weak Active)",
                     "B.1" = "B.1\n(Fac. Het)",
                     "B.2" = "B.2\n(Constit. Het)")
SUBCMPT_COLORS <- c("A.1" = "#E41A1C", "A.2" = "#FF7F00",
                     "B.1" = "#4DAF4A", "B.2" = "#377EB8")

fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

# ---- Load CALDER2 subcompartment labels -------------------------------------

SUBCMPT_PATH <- file.path(REPO_ROOT, "ML/cmpts/outputs/calder2/late/250402_subcompartment_labels_100kb.tsv")

if (!file.exists(SUBCMPT_PATH)) {
  stop("CALDER2 subcompartment file not found: ", SUBCMPT_PATH)
}

cat("Loading CALDER2 subcompartment labels...\n")
subcmpt <- read.table(SUBCMPT_PATH, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
subcmpt <- subcmpt %>% dplyr::filter(!is.na(ctrl_label) & ctrl_label != "NA")
cat(sprintf("  %d non-NA bins loaded\n", nrow(subcmpt)))

# ---- Assign subcompartment to gene bodies (Section 66 pattern) --------------

assign_subcompartment <- function(dmr_df, subcmpt_df) {
  dmr_df$midpoint <- (dmr_df$start + dmr_df$end) %/% 2
  dmr_df$bin_start <- (dmr_df$midpoint %/% 100000) * 100000 + 1
  dmr_df$bin_end <- dmr_df$bin_start + 99999
  merged <- dmr_df %>%
    left_join(
      subcmpt_df %>% dplyr::select(chr, bin_start, bin_end, ctrl_label, mut_label, label_changed),
      by = c("chr", "bin_start", "bin_end")
    )
  merged
}

cat("\nAssigning subcompartments to gene bodies...\n")
mc_cmpt <- assign_subcompartment(mc_dmr, subcmpt) %>% dplyr::filter(!is.na(ctrl_label))
hmc_cmpt <- assign_subcompartment(hmc_dmr, subcmpt) %>% dplyr::filter(!is.na(ctrl_label))
cat(sprintf("  5mC: %d genes assigned to subcompartments\n", nrow(mc_cmpt)))
cat(sprintf("  5hmC: %d genes assigned to subcompartments\n", nrow(hmc_cmpt)))

# =============================================================================
# PART 1: 5hmC SUBCOMPARTMENT DMR ANALYSIS (mirrors Section 66a/66b)
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("PART 1: 5hmC Subcompartment DMR Analysis\n")
cat("================================================================================\n\n")

# ---- 5hmC subcompartment summary -------------------------------------------

hmc_summary <- hmc_cmpt %>%
  group_by(ctrl_label) %>%
  summarise(
    total = n(),
    sig = sum(significant),
    pct_sig = 100 * sum(significant) / n(),
    hyper = sum(significant & mod_difference > 0),
    hypo = sum(significant & mod_difference < 0),
    .groups = "drop"
  ) %>%
  dplyr::mutate(ctrl_label = factor(ctrl_label, levels = SUBCMPT_ORDER))

cat("  5hmC per-subcompartment DMR summary:\n")
for (i in seq_len(nrow(hmc_summary))) {
  row <- hmc_summary[i, ]
  cat(sprintf("    %s: %d total, %d sig (%.1f%%), %d hyper, %d hypo\n",
              row$ctrl_label, row$total, row$sig, row$pct_sig, row$hyper, row$hypo))
}

# Chi-squared test
hmc_chisq_mat <- matrix(
  c(hmc_summary$sig, hmc_summary$total - hmc_summary$sig),
  ncol = 2,
  dimnames = list(hmc_summary$ctrl_label, c("Sig", "NotSig"))
)
hmc_chisq <- chisq.test(hmc_chisq_mat)
cat(sprintf("\n  Chi-squared (5hmC, subcompartment x significance): X²=%.1f, df=%d, %s\n",
            hmc_chisq$statistic, hmc_chisq$parameter, fmt_p(hmc_chisq$p.value)))

# ---- Panel 90a: 5hmC DMR fractions by subcompartment -----------------------

cat("\nCreating Figure 90a: 5hmC significant DMRs by subcompartment...\n")

hmc_bar_df <- hmc_summary %>%
  tidyr::pivot_longer(cols = c(hyper, hypo), names_to = "direction", values_to = "count") %>%
  dplyr::mutate(
    pct = 100 * count / total,
    direction = factor(ifelse(direction == "hyper", "Hypermethylated", "Hypomethylated"),
                       levels = c("Hypermethylated", "Hypomethylated"))
  )

p_90a <- ggplot(hmc_bar_df, aes(x = ctrl_label, y = pct, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", count, pct)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
  scale_fill_manual(values = COLORS$direction, name = "Direction") +
  scale_x_discrete(labels = SUBCMPT_LABELS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Significant 5hmC DMRs by Subcompartment",
    subtitle = sprintf("Chi² %s | Compare with 5mC gradient in Section 66a",
                        fmt_p(hmc_chisq$p.value)),
    x = "CALDER2 Subcompartment (Control)", y = "% of Gene Bodies with Sig DMR (q<0.05)"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_90a,
                        file.path(OUTPUT_DIR, "90a_5hmc_subcompartment_dmr_fraction"),
                        width = 10, height = 7)

# ---- Panel 90b: 5hmC direction split ---------------------------------------

cat("Creating Figure 90b: 5hmC direction split by subcompartment...\n")

hmc_dir_df <- hmc_cmpt %>%
  dplyr::filter(significant) %>%
  group_by(ctrl_label, direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(ctrl_label) %>%
  dplyr::mutate(
    total_sig = sum(n),
    pct = 100 * n / total_sig,
    ctrl_label = factor(ctrl_label, levels = SUBCMPT_ORDER),
    direction = factor(direction, levels = c("Hypermethylated", "Hypomethylated"))
  )

p_90b <- ggplot(hmc_dir_df, aes(x = ctrl_label, y = pct, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%\n(%d)", pct, n)),
            position = position_stack(vjust = 0.5), size = 3.5,
            color = "white", fontface = "bold") +
  scale_fill_manual(values = COLORS$direction, name = "Direction") +
  scale_x_discrete(labels = SUBCMPT_LABELS) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title = "5hmC Change Direction by Subcompartment",
    subtitle = "Compare directionality with 5mC (Section 66b)",
    x = "CALDER2 Subcompartment (Control)", y = "% of Significant DMRs"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_90b,
                        file.path(OUTPUT_DIR, "90b_5hmc_direction_by_subcompartment"),
                        width = 10, height = 7)

# =============================================================================
# PART 2: 5mC vs 5hmC GRADIENT COMPARISON
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("PART 2: 5mC vs 5hmC Gradient Comparison\n")
cat("================================================================================\n\n")

# ---- Compute 5mC summary (same as Section 66) for comparison ---------------

mc_summary <- mc_cmpt %>%
  group_by(ctrl_label) %>%
  summarise(
    total = n(),
    sig = sum(significant),
    pct_sig = 100 * sum(significant) / n(),
    hyper = sum(significant & mod_difference > 0),
    hypo = sum(significant & mod_difference < 0),
    .groups = "drop"
  ) %>%
  dplyr::mutate(ctrl_label = factor(ctrl_label, levels = SUBCMPT_ORDER))

# ---- Panel 90c: Net direction score across compartments --------------------

cat("Creating Figure 90c: gradient comparison (net direction)...\n")

gradient_df <- rbind(
  mc_summary %>%
    dplyr::mutate(
      pct_hyper = 100 * hyper / sig,
      pct_hypo = 100 * hypo / sig,
      net_direction = pct_hyper - pct_hypo,
      modality = "5mC"
    ),
  hmc_summary %>%
    dplyr::mutate(
      pct_hyper = 100 * hyper / sig,
      pct_hypo = 100 * hypo / sig,
      net_direction = pct_hyper - pct_hypo,
      modality = "5hmC"
    )
) %>%
  dplyr::mutate(
    modality = factor(modality, levels = c("5mC", "5hmC")),
    subcmpt_rank = as.numeric(factor(ctrl_label, levels = SUBCMPT_ORDER))
  )

# Spearman of net direction vs compartment rank
mc_grad_cor <- cor.test(
  gradient_df$subcmpt_rank[gradient_df$modality == "5mC"],
  gradient_df$net_direction[gradient_df$modality == "5mC"],
  method = "spearman", exact = FALSE
)
hmc_grad_cor <- cor.test(
  gradient_df$subcmpt_rank[gradient_df$modality == "5hmC"],
  gradient_df$net_direction[gradient_df$modality == "5hmC"],
  method = "spearman", exact = FALSE
)

cat(sprintf("  5mC net direction gradient: rho=%.3f\n", mc_grad_cor$estimate))
cat(sprintf("  5hmC net direction gradient: rho=%.3f\n", hmc_grad_cor$estimate))

p_90c <- ggplot(gradient_df, aes(x = ctrl_label, y = net_direction,
                                  color = modality, group = modality)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  geom_text(aes(label = sprintf("%.0f%%", net_direction)),
            vjust = -1.2, size = 3.5, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("5mC" = COLORS$methylation["5mC"],
                                 "5hmC" = COLORS$methylation["5hmC"]),
                     name = "Modality") +
  scale_x_discrete(labels = SUBCMPT_LABELS) +
  labs(
    title = "Methylation Gradient: Net Direction (% Hyper - % Hypo) by Subcompartment",
    subtitle = "Which modality shows a steeper gradient from active to heterochromatic?",
    x = "CALDER2 Subcompartment (Control)",
    y = "Net Direction Score (% Hyper - % Hypo)"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_90c,
                        file.path(OUTPUT_DIR, "90c_gradient_comparison_net_direction"),
                        width = 12, height = 7)

# ---- Panel 90d: Median effect size by subcompartment -----------------------

cat("Creating Figure 90d: median effect size by subcompartment...\n")

effect_df <- rbind(
  mc_cmpt %>%
    group_by(ctrl_label) %>%
    summarise(median_diff = median(mod_difference),
              mean_diff = mean(mod_difference),
              .groups = "drop") %>%
    dplyr::mutate(modality = "5mC"),
  hmc_cmpt %>%
    group_by(ctrl_label) %>%
    summarise(median_diff = median(mod_difference),
              mean_diff = mean(mod_difference),
              .groups = "drop") %>%
    dplyr::mutate(modality = "5hmC")
) %>%
  dplyr::mutate(
    ctrl_label = factor(ctrl_label, levels = SUBCMPT_ORDER),
    modality = factor(modality, levels = c("5mC", "5hmC"))
  )

p_90d <- ggplot(effect_df, aes(x = ctrl_label, y = median_diff * 100,
                                fill = modality)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f%%", median_diff * 100)),
            position = position_dodge(width = 0.7),
            vjust = ifelse(effect_df$median_diff > 0, -0.5, 1.5),
            size = 3) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = c("5mC" = COLORS$methylation["5mC"],
                                "5hmC" = COLORS$methylation["5hmC"]),
                    name = "Modality") +
  scale_x_discrete(labels = SUBCMPT_LABELS) +
  labs(
    title = "Median Methylation Change by Subcompartment",
    subtitle = "Absolute effect sizes across the compartment spectrum",
    x = "CALDER2 Subcompartment (Control)",
    y = "Median mod_difference (%)"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_90d,
                        file.path(OUTPUT_DIR, "90d_median_effect_by_subcompartment"),
                        width = 12, height = 7)

# =============================================================================
# PART 3: CONTINUOUS PC1 REGRESSION
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("PART 3: Continuous PC1 Regression\n")
cat("================================================================================\n\n")

# ---- Load Homer PC1 and join to genes via Gene_Name -------------------------

PC1_PATH <- file.path(REPO_ROOT, "tads/tad-pc-analysis/output/compartment_analysis/compartment_all_annotated.tsv")

if (!file.exists(PC1_PATH)) {
  stop("Homer PC1 annotated file not found: ", PC1_PATH)
}

cat("Loading Homer PC1 eigenvalues for gene-level join...\n")
pc1_raw <- read.table(PC1_PATH, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                      comment.char = "", fill = TRUE, quote = "")

pc1_per_gene <- pc1_raw %>%
  dplyr::filter(!is.na(Gene_Name) & Gene_Name != "") %>%
  group_by(Gene_Name) %>%
  summarise(
    mean_ctrl_pc1 = mean(ctrl_avg_PC1, na.rm = TRUE),
    mean_mut_pc1 = mean(mut_avg_PC1, na.rm = TRUE),
    n_bins = n(),
    .groups = "drop"
  )
cat(sprintf("  %d unique genes with PC1 values\n", nrow(pc1_per_gene)))

# Join with DMR data
mc_dedup <- mc_dmr %>%
  dplyr::select(gene, mc_diff = mod_difference, mc_q = dmr_qvalue,
                mc_sig = significant, mc_direction = direction) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

hmc_dedup <- hmc_dmr %>%
  dplyr::select(gene, hmc_diff = mod_difference, hmc_q = dmr_qvalue,
                hmc_sig = significant, hmc_direction = direction) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

gene_df <- pc1_per_gene %>%
  inner_join(mc_dedup, by = c("Gene_Name" = "gene")) %>%
  inner_join(hmc_dedup, by = c("Gene_Name" = "gene"))

cat(sprintf("  Matched: %d genes (%.1f%% of DMR genes)\n",
            nrow(gene_df), 100 * nrow(gene_df) / nrow(mc_dedup)))

stopifnot("Gene match rate too low" = nrow(gene_df) / nrow(mc_dedup) > 0.50)

# Also add subcompartment labels for B-compartment analysis later
mc_subcmpt_map <- mc_cmpt %>%
  dplyr::select(gene, ctrl_label) %>%
  dplyr::distinct(gene, .keep_all = TRUE)
gene_df <- gene_df %>%
  left_join(mc_subcmpt_map, by = c("Gene_Name" = "gene"))

# ---- Panel 90e: PC1 vs methylation change scatters -------------------------

cat("Creating Figure 90e: PC1 regression scatters...\n")

# 5mC regression
mod_mc <- lm(mc_diff ~ mean_ctrl_pc1, data = gene_df)
mc_r2 <- summary(mod_mc)$r.squared
mc_beta <- coef(mod_mc)["mean_ctrl_pc1"]
mc_pval <- summary(mod_mc)$coefficients["mean_ctrl_pc1", 4]
mc_spearman <- cor.test(gene_df$mean_ctrl_pc1, gene_df$mc_diff,
                         method = "spearman", exact = FALSE)

cat(sprintf("  5mC ~ PC1: beta=%.5f, R²=%.4f, %s, Spearman rho=%.4f\n",
            mc_beta, mc_r2, fmt_p(mc_pval), mc_spearman$estimate))

# 5hmC regression
mod_hmc <- lm(hmc_diff ~ mean_ctrl_pc1, data = gene_df)
hmc_r2 <- summary(mod_hmc)$r.squared
hmc_beta <- coef(mod_hmc)["mean_ctrl_pc1"]
hmc_pval <- summary(mod_hmc)$coefficients["mean_ctrl_pc1", 4]
hmc_spearman <- cor.test(gene_df$mean_ctrl_pc1, gene_df$hmc_diff,
                          method = "spearman", exact = FALSE)

cat(sprintf("  5hmC ~ PC1: beta=%.5f, R²=%.4f, %s, Spearman rho=%.4f\n",
            hmc_beta, hmc_r2, fmt_p(hmc_pval), hmc_spearman$estimate))

p_90e_mc <- ggplot(gene_df, aes(x = mean_ctrl_pc1, y = mc_diff)) +
  geom_point(alpha = 0.1, size = 0.3, color = COLORS$methylation["5mC"]) +
  geom_smooth(method = "loess", se = TRUE, color = "black", linewidth = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "grey30",
              linewidth = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  annotate("text", x = max(gene_df$mean_ctrl_pc1) * 0.6,
           y = max(gene_df$mc_diff) * 0.85,
           label = sprintf("beta = %.4f\nR² = %.4f\nrho = %.3f",
                           mc_beta, mc_r2, mc_spearman$estimate),
           hjust = 0, size = 3.5, fontface = "italic") +
  labs(title = "5mC Change vs Control PC1",
       x = "Control PC1 Eigenvalue", y = "Δ5mC (mut − ctrl)") +
  theme_biomodal()

p_90e_hmc <- ggplot(gene_df, aes(x = mean_ctrl_pc1, y = hmc_diff)) +
  geom_point(alpha = 0.1, size = 0.3, color = COLORS$methylation["5hmC"]) +
  geom_smooth(method = "loess", se = TRUE, color = "black", linewidth = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "grey30",
              linewidth = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  annotate("text", x = max(gene_df$mean_ctrl_pc1) * 0.6,
           y = min(gene_df$hmc_diff) * 0.85,
           label = sprintf("beta = %.4f\nR² = %.4f\nrho = %.3f",
                           hmc_beta, hmc_r2, hmc_spearman$estimate),
           hjust = 0, size = 3.5, fontface = "italic") +
  labs(title = "5hmC Change vs Control PC1",
       x = "Control PC1 Eigenvalue", y = "Δ5hmC (mut − ctrl)") +
  theme_biomodal()

p_90e <- p_90e_mc | p_90e_hmc
p_90e <- p_90e + plot_annotation(
  title = "Continuous PC1 vs Methylation Change: Gene-Level Regression",
  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
)

save_multiformat_ggplot(p_90e,
                        file.path(OUTPUT_DIR, "90e_pc1_regression_scatter"),
                        width = 16, height = 7)

# ---- Panel 90f: Standardized slope comparison with bootstrap ----------------

cat("Creating Figure 90f: gradient steepness comparison...\n")

gene_df$z_pc1 <- as.numeric(scale(gene_df$mean_ctrl_pc1))
gene_df$z_mc  <- as.numeric(scale(gene_df$mc_diff))
gene_df$z_hmc <- as.numeric(scale(gene_df$hmc_diff))

std_beta_mc  <- coef(lm(z_mc ~ z_pc1, data = gene_df))["z_pc1"]
std_beta_hmc <- coef(lm(z_hmc ~ z_pc1, data = gene_df))["z_pc1"]

cat(sprintf("  Standardized betas:\n"))
cat(sprintf("    5mC:  std_beta = %+.4f\n", std_beta_mc))
cat(sprintf("    5hmC: std_beta = %+.4f\n", std_beta_hmc))
cat(sprintf("    |5hmC| - |5mC| = %+.4f\n", abs(std_beta_hmc) - abs(std_beta_mc)))

set.seed(42)
B <- 5000
n_genes <- nrow(gene_df)
boot_diff <- numeric(B)
boot_mc <- numeric(B)
boot_hmc <- numeric(B)

for (b in seq_len(B)) {
  idx <- sample(n_genes, replace = TRUE)
  bd <- gene_df[idx, ]
  b_mc  <- coef(lm(z_mc ~ z_pc1, data = bd))["z_pc1"]
  b_hmc <- coef(lm(z_hmc ~ z_pc1, data = bd))["z_pc1"]
  boot_mc[b] <- b_mc
  boot_hmc[b] <- b_hmc
  boot_diff[b] <- abs(b_hmc) - abs(b_mc)
}

boot_ci <- quantile(boot_diff, c(0.025, 0.975))
boot_p <- 2 * min(mean(boot_diff > 0), mean(boot_diff < 0))

cat(sprintf("  Bootstrap comparison (|std_beta_hmC| - |std_beta_mC|):\n"))
cat(sprintf("    Observed: %+.4f\n", abs(std_beta_hmc) - abs(std_beta_mc)))
cat(sprintf("    95%% CI: [%+.4f, %+.4f]\n", boot_ci[1], boot_ci[2]))
cat(sprintf("    p = %.4f\n", boot_p))

if (abs(std_beta_hmc) > abs(std_beta_mc)) {
  cat("  --> 5hmC has STEEPER compartment gradient → consistent with TET impediment\n")
} else {
  cat("  --> 5mC has STEEPER compartment gradient → consistent with DNMT3A recruitment\n")
}

# Interaction test
stacked <- rbind(
  data.frame(gene = gene_df$Gene_Name, pc1 = gene_df$mean_ctrl_pc1,
             delta = gene_df$mc_diff, modality = "5mC", stringsAsFactors = FALSE),
  data.frame(gene = gene_df$Gene_Name, pc1 = gene_df$mean_ctrl_pc1,
             delta = gene_df$hmc_diff, modality = "5hmC", stringsAsFactors = FALSE)
)
stacked$modality <- factor(stacked$modality, levels = c("5mC", "5hmC"))
m_interaction <- lm(delta ~ pc1 * modality, data = stacked)
int_coef <- summary(m_interaction)$coefficients
cat(sprintf("\n  Interaction test (delta ~ pc1 * modality):\n"))
cat(sprintf("    pc1:modality5hmC coefficient: %.5f, %s\n",
            int_coef["pc1:modality5hmC", 1],
            fmt_p(int_coef["pc1:modality5hmC", 4])))

slope_df <- data.frame(
  modality = c("5mC", "5hmC"),
  std_beta = c(std_beta_mc, std_beta_hmc),
  abs_std_beta = c(abs(std_beta_mc), abs(std_beta_hmc)),
  ci_lo = c(quantile(boot_mc, 0.025), quantile(boot_hmc, 0.025)),
  ci_hi = c(quantile(boot_mc, 0.975), quantile(boot_hmc, 0.975)),
  stringsAsFactors = FALSE
)

p_90f <- ggplot(slope_df, aes(x = modality, y = std_beta,
                               fill = modality)) +
  geom_col(width = 0.5, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15, linewidth = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = c("5mC" = COLORS$methylation["5mC"],
                                "5hmC" = COLORS$methylation["5hmC"]),
                    guide = "none") +
  annotate("text", x = 1.5, y = max(slope_df$ci_hi) * 1.3,
           label = sprintf("|hmC| - |mC| = %+.4f\n95%% CI [%+.4f, %+.4f]\np = %.4f",
                           abs(std_beta_hmc) - abs(std_beta_mc),
                           boot_ci[1], boot_ci[2], boot_p),
           size = 3.5, fontface = "italic") +
  labs(
    title = "Compartment Gradient Steepness: 5mC vs 5hmC",
    subtitle = "Standardized beta from PC1 regression (bootstrap 95% CI)",
    x = "Modality", y = "Standardized Beta (Δmethylation ~ PC1)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_90f,
                        file.path(OUTPUT_DIR, "90f_gradient_steepness_comparison"),
                        width = 8, height = 7)

# =============================================================================
# PART 4: MECHANISM DIFFERENTIATION
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("PART 4: Mechanism Differentiation (Mediation + Variance Partition)\n")
cat("================================================================================\n\n")

# ---- Panel 90g: Mediation analysis ------------------------------------------

cat("Running mediation analysis (following Section 61 bootstrap pattern)...\n")

run_mediation <- function(gene_df, exposure, mediator, outcome, label, B = 5000) {
  path_a  <- lm(as.formula(paste(mediator, "~", exposure)), data = gene_df)
  path_cp <- lm(as.formula(paste(outcome, "~", exposure, "+", mediator)), data = gene_df)
  path_c  <- lm(as.formula(paste(outcome, "~", exposure)), data = gene_df)

  coef_a  <- coef(path_a)[exposure]
  coef_b  <- coef(path_cp)[mediator]
  coef_c  <- coef(path_c)[exposure]
  coef_cp_val <- coef(path_cp)[exposure]
  indirect_obs <- coef_a * coef_b
  proportion_mediated <- indirect_obs / coef_c

  cat(sprintf("\n  %s:\n", label))
  cat(sprintf("    Path a (%s → %s):       beta = %+.6f, %s\n",
              exposure, mediator, coef_a,
              fmt_p(summary(path_a)$coefficients[exposure, 4])))
  cat(sprintf("    Path b (%s → %s|%s): beta = %+.6f, %s\n",
              mediator, outcome, exposure, coef_b,
              fmt_p(summary(path_cp)$coefficients[mediator, 4])))
  cat(sprintf("    Total (c):              beta = %+.6f\n", coef_c))
  cat(sprintf("    Direct (c'):            beta = %+.6f\n", coef_cp_val))
  cat(sprintf("    Indirect (a × b):       %.6f\n", indirect_obs))
  cat(sprintf("    Proportion mediated:    %.1f%%\n", 100 * proportion_mediated))

  # VIF check
  vif_val <- 1 / (1 - summary(lm(as.formula(paste(mediator, "~", exposure)),
                                   data = gene_df))$r.squared)
  cat(sprintf("    VIF (mediator~exposure): %.2f%s\n",
              vif_val, ifelse(vif_val > 5, " (WARNING: high multicollinearity)", "")))

  # Bootstrap
  set.seed(42)
  n <- nrow(gene_df)
  boot_indirect <- numeric(B)
  for (b_i in seq_len(B)) {
    idx <- sample(n, replace = TRUE)
    bd <- gene_df[idx, ]
    a_b <- coef(lm(as.formula(paste(mediator, "~", exposure)), data = bd))[exposure]
    cp_b <- lm(as.formula(paste(outcome, "~", exposure, "+", mediator)), data = bd)
    b_b <- coef(cp_b)[mediator]
    boot_indirect[b_i] <- a_b * b_b
  }

  boot_ci <- quantile(boot_indirect, c(0.025, 0.975))
  boot_p <- 2 * min(mean(boot_indirect > 0), mean(boot_indirect < 0))

  cat(sprintf("    Bootstrap 95%% CI: [%.6f, %.6f]\n", boot_ci[1], boot_ci[2]))
  cat(sprintf("    Bootstrap p = %.4f\n", boot_p))

  data.frame(
    model = label,
    path = c("a", "b", "c (total)", "c' (direct)", "indirect (a*b)"),
    coefficient = c(coef_a, coef_b, coef_c, coef_cp_val, indirect_obs),
    ci_lo = c(confint(path_a, exposure)[1], confint(path_cp, mediator)[1],
              confint(path_c, exposure)[1], confint(path_cp, exposure)[1],
              boot_ci[1]),
    ci_hi = c(confint(path_a, exposure)[2], confint(path_cp, mediator)[2],
              confint(path_c, exposure)[2], confint(path_cp, exposure)[2],
              boot_ci[2]),
    p_value = c(summary(path_a)$coefficients[exposure, 4],
                summary(path_cp)$coefficients[mediator, 4],
                summary(path_c)$coefficients[exposure, 4],
                summary(path_cp)$coefficients[exposure, 4],
                boot_p),
    proportion_mediated = c(NA, NA, NA, NA, proportion_mediated),
    vif = c(NA, NA, NA, NA, vif_val),
    stringsAsFactors = FALSE
  )
}

# Model 1 (TET-first): PC1 → Δ5hmC → Δ5mC
med_tet <- run_mediation(gene_df,
  exposure = "mean_ctrl_pc1", mediator = "hmc_diff", outcome = "mc_diff",
  label = "Model 1: TET-first (PC1 -> hmC -> mC)")

# Model 2 (DNMT3A-first): PC1 → Δ5mC → Δ5hmC
med_dnmt <- run_mediation(gene_df,
  exposure = "mean_ctrl_pc1", mediator = "mc_diff", outcome = "hmc_diff",
  label = "Model 2: DNMT3A-first (PC1 -> mC -> hmC)")

# Compare
prop_tet <- med_tet$proportion_mediated[med_tet$path == "indirect (a*b)"]
prop_dnmt <- med_dnmt$proportion_mediated[med_dnmt$path == "indirect (a*b)"]

cat(sprintf("\n  MEDIATION COMPARISON:\n"))
cat(sprintf("    TET-first proportion mediated:   %.1f%%\n", 100 * prop_tet))
cat(sprintf("    DNMT3A-first proportion mediated: %.1f%%\n", 100 * prop_dnmt))

if (abs(prop_tet) > abs(prop_dnmt)) {
  cat("  --> TET-first path shows STRONGER mediation\n")
} else {
  cat("  --> DNMT3A-first path shows STRONGER mediation\n")
}

med_compare <- data.frame(
  model = c("TET-first (PC1->hmC->mC)", "DNMT3A-first (PC1->mC->hmC)"),
  indirect_effect = c(med_tet$coefficient[5], med_dnmt$coefficient[5]),
  indirect_ci_lo = c(med_tet$ci_lo[5], med_dnmt$ci_lo[5]),
  indirect_ci_hi = c(med_tet$ci_hi[5], med_dnmt$ci_hi[5]),
  indirect_p = c(med_tet$p_value[5], med_dnmt$p_value[5]),
  proportion_mediated = c(prop_tet, prop_dnmt),
  vif = c(med_tet$vif[5], med_dnmt$vif[5]),
  stringsAsFactors = FALSE
)

p_90g <- ggplot(med_compare, aes(x = model, y = 100 * proportion_mediated,
                                  fill = model)) +
  geom_col(width = 0.5, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = 100 * indirect_ci_lo / indirect_effect * proportion_mediated,
                    ymax = 100 * indirect_ci_hi / indirect_effect * proportion_mediated),
                width = 0.15, linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%\np=%.4f", 100 * proportion_mediated, indirect_p)),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("TET-first (PC1->hmC->mC)" = COLORS$methylation["5hmC"],
                                "DNMT3A-first (PC1->mC->hmC)" = COLORS$methylation["5mC"]),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(
    title = "Mediation Analysis: Which Mechanism Mediates the PC1 → Methylation Link?",
    subtitle = "Higher proportion mediated = stronger evidence for that causal path",
    x = NULL, y = "Proportion of Total Effect Mediated (%)"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(size = 9))

save_multiformat_ggplot(p_90g,
                        file.path(OUTPUT_DIR, "90g_mediation_comparison"),
                        width = 10, height = 7)

# ---- Panel 90h: Variance partition -----------------------------------------

cat("\nComputing variance partition...\n")

partition_variance <- function(r2_mod1, r2_mod2, r2_both, label) {
  mod1_unique   <- r2_both - r2_mod2
  mod2_unique   <- r2_both - r2_mod1
  shared        <- r2_mod1 + r2_mod2 - r2_both
  unexplained   <- 1 - r2_both
  data.frame(
    outcome = label,
    component = c("PC1 unique", "Other modality unique", "Shared", "Unexplained"),
    fraction = c(mod1_unique, mod2_unique, shared, unexplained),
    stringsAsFactors = FALSE
  )
}

# For Δ5mC: how much do PC1 vs Δ5hmC explain?
r2_mc_pc1   <- summary(lm(mc_diff ~ mean_ctrl_pc1, data = gene_df))$r.squared
r2_mc_hmc   <- summary(lm(mc_diff ~ hmc_diff, data = gene_df))$r.squared
r2_mc_both  <- summary(lm(mc_diff ~ mean_ctrl_pc1 + hmc_diff, data = gene_df))$r.squared

# For Δ5hmC: how much do PC1 vs Δ5mC explain?
r2_hmc_pc1  <- summary(lm(hmc_diff ~ mean_ctrl_pc1, data = gene_df))$r.squared
r2_hmc_mc   <- summary(lm(hmc_diff ~ mc_diff, data = gene_df))$r.squared
r2_hmc_both <- summary(lm(hmc_diff ~ mean_ctrl_pc1 + mc_diff, data = gene_df))$r.squared

vp <- rbind(
  partition_variance(r2_mc_pc1, r2_mc_hmc, r2_mc_both, "Predicting Δ5mC"),
  partition_variance(r2_hmc_pc1, r2_hmc_mc, r2_hmc_both, "Predicting Δ5hmC")
)

cat("  Variance partition:\n")
for (i in seq_len(nrow(vp))) {
  cat(sprintf("    [%s] %-22s %.4f (%.1f%%)\n",
              vp$outcome[i], vp$component[i], vp$fraction[i], 100 * vp$fraction[i]))
}

vp$component <- factor(vp$component,
  levels = c("Unexplained", "Other modality unique", "Shared", "PC1 unique"))
vp_colors <- c("PC1 unique" = "#4393C3",
               "Shared" = "#FDB863",
               "Other modality unique" = "#D6604D",
               "Unexplained" = "grey85")

p_90h <- ggplot(vp, aes(x = outcome, y = fraction, fill = component)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(aes(label = ifelse(fraction > 0.02,
                sprintf("%.1f%%", 100 * fraction), "")),
            position = position_stack(vjust = 0.5), size = 3.5) +
  scale_fill_manual(values = vp_colors, name = "Component") +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Variance Partition: PC1 vs Other Modality",
    subtitle = "How much variance does compartment identity explain uniquely?",
    x = NULL, y = "Fraction of Total Variance"
  ) +
  theme_biomodal() +
  theme(legend.position = "right")

save_multiformat_ggplot(p_90h,
                        file.path(OUTPUT_DIR, "90h_variance_partition"),
                        width = 10, height = 7)

# =============================================================================
# PART 5: B-COMPARTMENT BEHAVIOR
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("PART 5: B-Compartment Behavior\n")
cat("================================================================================\n\n")

# ---- Panel 90i: B-compartment methylation changes ---------------------------

cat("Creating Figure 90i: B-compartment behavior...\n")

b_genes <- gene_df %>%
  dplyr::filter(ctrl_label %in% c("B.1", "B.2"))

cat(sprintf("  B-compartment genes: %d\n", nrow(b_genes)))
cat(sprintf("    mC sig: %d (%.1f%%)\n",
            sum(b_genes$mc_sig), 100 * mean(b_genes$mc_sig)))
cat(sprintf("    hmC sig: %d (%.1f%%)\n",
            sum(b_genes$hmc_sig), 100 * mean(b_genes$hmc_sig)))

# 2x2 contingency
b_contingency <- table(mC_sig = b_genes$mc_sig, hmC_sig = b_genes$hmc_sig)
cat("  B-compartment 2x2 (mC sig x hmC sig):\n")
print(b_contingency)
b_fisher <- fisher.test(b_contingency)
cat(sprintf("  Fisher's OR = %.2f, %s\n",
            b_fisher$estimate, fmt_p(b_fisher$p.value)))

b_violin_df <- rbind(
  data.frame(gene = b_genes$Gene_Name, subcmpt = b_genes$ctrl_label,
             delta = b_genes$mc_diff, modality = "Δ5mC", stringsAsFactors = FALSE),
  data.frame(gene = b_genes$Gene_Name, subcmpt = b_genes$ctrl_label,
             delta = b_genes$hmc_diff, modality = "Δ5hmC", stringsAsFactors = FALSE)
) %>%
  dplyr::mutate(modality = factor(modality, levels = c("Δ5mC", "Δ5hmC")))

# Wilcoxon within B compartment for each modality
mc_b_wilcox <- wilcox.test(b_genes$mc_diff, mu = 0)
hmc_b_wilcox <- wilcox.test(b_genes$hmc_diff, mu = 0)
cat(sprintf("  B-compartment Wilcoxon (vs 0):\n"))
cat(sprintf("    Δ5mC: median=%.5f, %s\n",
            median(b_genes$mc_diff), fmt_p(mc_b_wilcox$p.value)))
cat(sprintf("    Δ5hmC: median=%.5f, %s\n",
            median(b_genes$hmc_diff), fmt_p(hmc_b_wilcox$p.value)))

p_90i <- ggplot(b_violin_df, aes(x = modality, y = delta * 100, fill = modality)) +
  geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ subcmpt, labeller = labeller(
    subcmpt = c("B.1" = "B.1 (Fac. Het)", "B.2" = "B.2 (Constit. Het)"))) +
  scale_fill_manual(values = c("Δ5mC" = COLORS$methylation["5mC"],
                                "Δ5hmC" = COLORS$methylation["5hmC"]),
                    guide = "none") +
  labs(
    title = "Methylation Changes in B-Compartment Genes",
    subtitle = sprintf("B-cmpt Fisher OR (mC sig x hmC sig) = %.2f, %s",
                        b_fisher$estimate, fmt_p(b_fisher$p.value)),
    x = NULL, y = "Methylation Difference (%)"
  ) +
  theme_biomodal() +
  theme(strip.text = element_text(size = 10, face = "bold"))

save_multiformat_ggplot(p_90i,
                        file.path(OUTPUT_DIR, "90i_b_compartment_behavior"),
                        width = 10, height = 7)

# ---- Combined figures -------------------------------------------------------

cat("\nCreating combined figures...\n")

p_90_descriptive <- (p_90a | p_90b) /
  (p_90c | p_90d) +
  plot_annotation(
    title = "5hmC Subcompartment Analysis + Gradient Comparison",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )
save_multiformat_ggplot(p_90_descriptive,
                        file.path(OUTPUT_DIR, "90j_descriptive_composite"),
                        width = 16, height = 14)

p_90_mechanism <- (p_90e) /
  (p_90f | p_90g | p_90h) /
  p_90i +
  plot_layout(heights = c(1, 1, 0.8)) +
  plot_annotation(
    title = "Mechanism Differentiation: TET Impediment vs DNMT3A Recruitment",
    subtitle = "Compartment-dependent gradients, mediation, and variance partition",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  )
save_multiformat_ggplot(p_90_mechanism,
                        file.path(OUTPUT_DIR, "90j_mechanism_composite"),
                        width = 18, height = 18)

# ---- Save tables ------------------------------------------------------------

cat("\nSaving tables...\n")

# 5hmC subcompartment summary
write.table(hmc_summary,
            file.path(TABLES_DIR, "90_5hmc_subcompartment_dmr_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 90_5hmc_subcompartment_dmr_summary.tsv\n")

# Gradient comparison
write.table(gradient_df %>% dplyr::select(ctrl_label, modality, pct_sig,
                                            pct_hyper, pct_hypo, net_direction),
            file.path(TABLES_DIR, "90_gradient_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 90_gradient_comparison.tsv\n")

# PC1 regression coefficients
reg_coefs <- data.frame(
  modality = c("5mC", "5hmC"),
  beta = c(mc_beta, hmc_beta),
  std_beta = c(std_beta_mc, std_beta_hmc),
  r_squared = c(mc_r2, hmc_r2),
  spearman_rho = c(mc_spearman$estimate, hmc_spearman$estimate),
  p_value = c(mc_pval, hmc_pval),
  stringsAsFactors = FALSE
)
write.table(reg_coefs,
            file.path(TABLES_DIR, "90_pc1_regression_coefficients.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 90_pc1_regression_coefficients.tsv\n")

# Gradient steepness bootstrap
boot_stats <- data.frame(
  comparison = "|std_beta_hmC| - |std_beta_mC|",
  observed = abs(std_beta_hmc) - abs(std_beta_mc),
  boot_ci_lo = boot_ci[1],
  boot_ci_hi = boot_ci[2],
  boot_p = boot_p,
  interpretation = ifelse(abs(std_beta_hmc) > abs(std_beta_mc),
                          "5hmC steeper -> TET impediment", "5mC steeper -> DNMT3A recruitment"),
  stringsAsFactors = FALSE
)
write.table(boot_stats,
            file.path(TABLES_DIR, "90_gradient_steepness_bootstrap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 90_gradient_steepness_bootstrap.tsv\n")

# Interaction test
int_results <- data.frame(
  term = rownames(int_coef),
  estimate = int_coef[, 1],
  std_error = int_coef[, 2],
  t_value = int_coef[, 3],
  p_value = int_coef[, 4],
  stringsAsFactors = FALSE
)
write.table(int_results,
            file.path(TABLES_DIR, "90_interaction_test.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 90_interaction_test.tsv\n")

# Mediation results
write.table(med_tet,
            file.path(TABLES_DIR, "90_mediation_tet_first.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 90_mediation_tet_first.tsv\n")

write.table(med_dnmt,
            file.path(TABLES_DIR, "90_mediation_dnmt3a_first.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 90_mediation_dnmt3a_first.tsv\n")

write.table(med_compare,
            file.path(TABLES_DIR, "90_mediation_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 90_mediation_comparison.tsv\n")

# Variance partition
write.table(vp,
            file.path(TABLES_DIR, "90_variance_partition.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 90_variance_partition.tsv\n")

# B-compartment contingency
b_cont_df <- data.frame(
  mc_sig = c(TRUE, TRUE, FALSE, FALSE),
  hmc_sig = c(TRUE, FALSE, TRUE, FALSE),
  count = as.vector(b_contingency),
  fisher_or = b_fisher$estimate,
  fisher_p = b_fisher$p.value,
  stringsAsFactors = FALSE
)
write.table(b_cont_df,
            file.path(TABLES_DIR, "90_b_compartment_contingency.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 90_b_compartment_contingency.tsv\n")

# Per-gene PC1 assignment
gene_export <- gene_df %>%
  dplyr::select(Gene_Name, mean_ctrl_pc1, mean_mut_pc1, n_bins,
                ctrl_label, mc_diff, mc_sig, mc_direction,
                hmc_diff, hmc_sig, hmc_direction)
write.table(gene_export,
            file.path(TABLES_DIR, "90_per_gene_pc1_assignment.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: 90_per_gene_pc1_assignment.tsv (%d genes)\n", nrow(gene_export)))

# ---- Mechanism summary table ------------------------------------------------

cat("\n")
cat("================================================================================\n")
cat("MECHANISM DIFFERENTIATION SUMMARY\n")
cat("================================================================================\n\n")

mechanism_summary <- data.frame(
  test = c("Gradient steepness (|std_beta|)",
           "Interaction (PC1 x modality)",
           "Mediation proportion (TET-first)",
           "Mediation proportion (DNMT3A-first)",
           "Variance: PC1 unique for Δ5mC",
           "Variance: PC1 unique for Δ5hmC",
           "B-cmpt median Δ5mC",
           "B-cmpt median Δ5hmC"),
  mc_value = c(abs(std_beta_mc),
               NA,
               NA,
               100 * prop_dnmt,
               100 * (r2_mc_both - r2_mc_hmc),
               NA,
               100 * median(b_genes$mc_diff),
               NA),
  hmc_value = c(abs(std_beta_hmc),
                NA,
                100 * prop_tet,
                NA,
                NA,
                100 * (r2_hmc_both - r2_hmc_mc),
                NA,
                100 * median(b_genes$hmc_diff)),
  favors = c(ifelse(abs(std_beta_hmc) > abs(std_beta_mc), "TET impediment", "DNMT3A recruitment"),
             ifelse(int_coef["pc1:modality5hmC", 4] < 0.05, "Different slopes (sig interaction)", "Similar slopes"),
             ifelse(abs(prop_tet) > abs(prop_dnmt), "TET impediment", "DNMT3A recruitment"),
             ifelse(abs(prop_dnmt) > abs(prop_tet), "DNMT3A recruitment", "TET impediment"),
             "See values",
             "See values",
             "See values",
             "See values"),
  stringsAsFactors = FALSE
)

write.table(mechanism_summary,
            file.path(TABLES_DIR, "90_mechanism_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 90_mechanism_summary.tsv\n")

for (i in seq_len(nrow(mechanism_summary))) {
  row <- mechanism_summary[i, ]
  cat(sprintf("  %-45s mC=%.4f  hmC=%.4f  → %s\n",
              row$test,
              ifelse(is.na(row$mc_value), NA, row$mc_value),
              ifelse(is.na(row$hmc_value), NA, row$hmc_value),
              row$favors))
}

cat("\n")
cat("================================================================================\n")
cat("SECTION 90 COMPLETE\n")
cat("================================================================================\n")
