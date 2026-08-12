# biomodal/downstream/scripts/viz_sections/section_29_ab_compartment_methylation_mapping.R
# Section 29: A/B Compartment Methylation Mapping (TODO Task 2)
# Tests whether BAP1-KO methylation changes follow the Lopez-Moyado TET-KO
# signature: hypermethylation in A compartment (euchromatin) + hypomethylation
# in B compartment (heterochromatin), indicating DNMT3A redistribution.
# A match supports "convergent mechanisms" with TET-KO; a mismatch suggests
# a distinct BAP1-mediated pathway.
#
# Sub-analyses:
#   29a-b: mC/hmC mod_difference violins by A vs B compartment
#   29c-d: mC/hmC mod_difference violins by compartment shift (B->A, A->B, Stable)
#   29e:   Stacked bar — DMR direction proportions by compartment and shift
#   29f:   Scatter — PC1 vs mC mod_difference with loess + Spearman rho
#   29g:   Composite summary panel (29a + 29b + 29e)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_29_ab_compartment_methylation_mapping.R

source("scripts/viz_sections/_shared_config.R")

cat("\n")
cat("================================================================================\n")
cat("SECTION 29: A/B COMPARTMENT METHYLATION MAPPING\n")
cat("================================================================================\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================

# HOMER compartment file (LATE timepoint, 104,071 bins)
COMPARTMENT_FILE <- file.path(BASE_DIR, "../../tads/tad-pc-analysis/inputs/late/diffPC/diffcompartments.txt")

# Compartment shift thresholds (matching existing pipeline standard)
SHIFT_FDR <- 0.05
SHIFT_DIFF <- 0.30

# Section output directory
SECTION_DIR <- file.path(OUTPUT_DIR, "29_ab_compartment_methylation")
dir.create(SECTION_DIR, recursive = TRUE, showWarnings = FALSE)

# Colors
COMPARTMENT_COLORS <- c(A = "#E41A1C", B = "#2166AC")
SHIFT_COLORS <- c("B to A" = "#D7191C", "A to B" = "#2C7BB6", Stable = "grey70")

# Helper: format p-value
fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

cat("Validating inputs...\n")
stopifnot("Compartment file not found" = file.exists(COMPARTMENT_FILE))
stopifnot("mc_dmr not loaded from shared config" = exists("mc_dmr") && nrow(mc_dmr) > 0)
stopifnot("hmc_dmr not loaded from shared config" = exists("hmc_dmr") && nrow(hmc_dmr) > 0)
cat("  All inputs validated.\n\n")

# =============================================================================
# STEP 1 (Task 2a): LOAD HOMER DATA, ASSIGN GENES TO COMPARTMENTS
# =============================================================================

cat("================================================================================\n")
cat("STEP 1: Load HOMER compartment data and assign genes\n")
cat("================================================================================\n\n")

# Read HOMER diffcompartments.txt
comp_raw <- read.table(COMPARTMENT_FILE, sep = "\t", header = TRUE,
                       stringsAsFactors = FALSE, check.names = FALSE,
                       comment.char = "", quote = "")
cat(sprintf("  Loaded %d genomic bins\n", nrow(comp_raw)))

# Identify columns by grep pattern (following compartment_volcano_plot.R pattern)
pc1_cols <- grep("bedGraph avg over given bp", names(comp_raw), value = TRUE)
stopifnot("Expected 6 PC1 sample columns" = length(pc1_cols) == 6)
ctrl_pc1_cols <- pc1_cols[1:3]  # first 3 = ctrl
mut_pc1_cols <- pc1_cols[4:6]   # last 3 = mut

difference_col <- grep("ctrl vs\\. mut Difference",
                        names(comp_raw), value = TRUE, ignore.case = TRUE)[1]
adj_pvalue_col <- grep("ctrl vs\\. mut adj\\. p-value",
                        names(comp_raw), value = TRUE, ignore.case = TRUE)[1]
stopifnot("Difference column not found" = !is.na(difference_col))
stopifnot("Adj. p-value column not found" = !is.na(adj_pvalue_col))

cat(sprintf("  PC1 columns: %d ctrl, %d mut\n", length(ctrl_pc1_cols), length(mut_pc1_cols)))
cat(sprintf("  Difference column: %s\n", difference_col))
cat(sprintf("  Adj. p-value column: %s\n", adj_pvalue_col))

# Build clean working dataframe
comp_df <- data.frame(
  chr = comp_raw[["Chr"]],
  start = comp_raw[["Start"]],
  end = comp_raw[["End"]],
  gene_name = comp_raw[["Gene Name"]],
  dist_to_tss = as.numeric(comp_raw[["Distance to TSS"]]),
  ctrl_pc1_1 = as.numeric(comp_raw[[ctrl_pc1_cols[1]]]),
  ctrl_pc1_2 = as.numeric(comp_raw[[ctrl_pc1_cols[2]]]),
  ctrl_pc1_3 = as.numeric(comp_raw[[ctrl_pc1_cols[3]]]),
  mut_pc1_1 = as.numeric(comp_raw[[mut_pc1_cols[1]]]),
  mut_pc1_2 = as.numeric(comp_raw[[mut_pc1_cols[2]]]),
  mut_pc1_3 = as.numeric(comp_raw[[mut_pc1_cols[3]]]),
  difference = as.numeric(comp_raw[[difference_col]]),
  adj_pvalue = as.numeric(comp_raw[[adj_pvalue_col]]),
  stringsAsFactors = FALSE
)

# Compute control mean PC1 for A/B classification
comp_df$mean_ctrl_pc1 <- rowMeans(comp_df[, c("ctrl_pc1_1", "ctrl_pc1_2", "ctrl_pc1_3")],
                                   na.rm = TRUE)

# Classify compartment based on control mean PC1
comp_df$compartment <- ifelse(comp_df$mean_ctrl_pc1 > 0, "A", "B")

# Classify compartment shift
comp_df$shift <- case_when(
  comp_df$adj_pvalue < SHIFT_FDR & comp_df$difference > SHIFT_DIFF  ~ "B to A",
  comp_df$adj_pvalue < SHIFT_FDR & comp_df$difference < -SHIFT_DIFF ~ "A to B",
  TRUE ~ "Stable"
)

cat(sprintf("\n  Compartment assignment (control mean PC1):\n"))
cat(sprintf("    A compartment: %d bins (%.1f%%)\n",
            sum(comp_df$compartment == "A"),
            100 * mean(comp_df$compartment == "A")))
cat(sprintf("    B compartment: %d bins (%.1f%%)\n",
            sum(comp_df$compartment == "B"),
            100 * mean(comp_df$compartment == "B")))
cat(sprintf("\n  Compartment shifts (FDR < %.2f, |Diff| > %.2f):\n", SHIFT_FDR, SHIFT_DIFF))
cat(sprintf("    B to A: %d bins\n", sum(comp_df$shift == "B to A")))
cat(sprintf("    A to B: %d bins\n", sum(comp_df$shift == "A to B")))
cat(sprintf("    Stable: %d bins\n", sum(comp_df$shift == "Stable")))

# Deduplicate: one gene per bin — keep bin closest to TSS
comp_genes <- comp_df %>%
  dplyr::filter(!is.na(gene_name) & gene_name != "") %>%
  dplyr::group_by(gene_name) %>%
  dplyr::slice_min(abs(dist_to_tss), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

cat(sprintf("\n  After deduplication (closest bin per gene): %d genes\n", nrow(comp_genes)))

# Prepare DMR data (deduplicate by gene, keep first)
mc_dedup <- mc_dmr %>%
  dplyr::select(gene, mc_diff = mod_difference, mc_q = dmr_qvalue,
                mc_sig = significant, mc_direction = direction) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

hmc_dedup <- hmc_dmr %>%
  dplyr::select(gene, hmc_diff = mod_difference, hmc_q = dmr_qvalue,
                hmc_sig = significant, hmc_direction = direction) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

# Inner join compartment genes with mC and hmC DMRs
merged <- comp_genes %>%
  dplyr::inner_join(mc_dedup, by = c("gene_name" = "gene")) %>%
  dplyr::inner_join(hmc_dedup, by = c("gene_name" = "gene"))

# Report match rates
n_comp_genes <- length(unique(comp_genes$gene_name))
n_mc_genes <- length(unique(mc_dedup$gene))
n_overlap <- nrow(merged)
match_rate_mc <- n_overlap / n_mc_genes
match_rate_comp <- n_overlap / n_comp_genes

cat(sprintf("\n  Gene matching:\n"))
cat(sprintf("    Compartment genes: %d unique\n", n_comp_genes))
cat(sprintf("    mC DMR genes: %d unique\n", n_mc_genes))
cat(sprintf("    hmC DMR genes: %d unique\n", length(unique(hmc_dedup$gene))))
cat(sprintf("    Matched (inner join): %d genes\n", n_overlap))
cat(sprintf("    Match rate (of DMR genes): %.1f%%\n", 100 * match_rate_mc))
cat(sprintf("    Match rate (of compartment genes): %.1f%%\n", 100 * match_rate_comp))

stopifnot("Match rate too low — check gene name format" = match_rate_mc > 0.50)

# Set factor levels for plotting
merged$compartment <- factor(merged$compartment, levels = c("A", "B"))
merged$shift <- factor(merged$shift, levels = c("B to A", "A to B", "Stable"))

cat(sprintf("\n  Final merged dataset:\n"))
cat(sprintf("    A compartment: %d genes (%.1f%%)\n",
            sum(merged$compartment == "A"),
            100 * mean(merged$compartment == "A")))
cat(sprintf("    B compartment: %d genes (%.1f%%)\n",
            sum(merged$compartment == "B"),
            100 * mean(merged$compartment == "B")))

# =============================================================================
# STEP 2 (Task 2b): FISHER'S EXACT TESTS FOR COMPARTMENT ENRICHMENT
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("STEP 2: Fisher's exact tests — DMR direction x compartment enrichment\n")
cat("================================================================================\n\n")

run_compartment_fisher <- function(df, test_name, positive_set, compartment_expected) {
  # positive_set: logical vector (TRUE = in positive set)
  # compartment_expected: "A" or "B" — the compartment we expect enrichment in

  in_expected <- df$compartment == compartment_expected
  in_other <- !in_expected

  a <- sum(positive_set & in_expected)
  b <- sum(positive_set & in_other)
  c <- sum(!positive_set & in_expected)
  d <- sum(!positive_set & in_other)

  fisher_mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                       dimnames = list(c("Positive", "Negative"),
                                       c(compartment_expected,
                                         ifelse(compartment_expected == "A", "B", "A"))))
  fisher_res <- fisher.test(fisher_mat)

  pct_pos_expected <- ifelse((a + c) > 0, 100 * a / (a + c), 0)
  pct_pos_other <- ifelse((b + d) > 0, 100 * b / (b + d), 0)

  cat(sprintf("  %s:\n", test_name))
  cat(sprintf("    OR = %.2f (95%% CI: %.2f-%.2f), %s\n",
              fisher_res$estimate, fisher_res$conf.int[1], fisher_res$conf.int[2],
              fmt_p(fisher_res$p.value)))
  cat(sprintf("    %% positive in %s: %.1f%%, in %s: %.1f%%\n",
              compartment_expected, pct_pos_expected,
              ifelse(compartment_expected == "A", "B", "A"), pct_pos_other))
  cat(sprintf("    Contingency: a=%d, b=%d, c=%d, d=%d\n\n", a, b, c, d))

  data.frame(
    test = test_name,
    positive_set_n = sum(positive_set),
    expected_compartment = compartment_expected,
    a = a, b = b, c = c, d = d,
    odds_ratio = fisher_res$estimate,
    ci_lower = fisher_res$conf.int[1],
    ci_upper = fisher_res$conf.int[2],
    p_value = fisher_res$p.value,
    pct_in_expected = pct_pos_expected,
    pct_in_other = pct_pos_other,
    stringsAsFactors = FALSE
  )
}

# Four directional tests (TET-KO signature predictions)
fisher_results <- dplyr::bind_rows(
  # mC hypermethylation enriched in A? (TET-KO: DNMT3A shifts to euchromatin)
  run_compartment_fisher(merged, "mC hyper -> A enriched",
                         positive_set = merged$mc_sig & merged$mc_diff > 0,
                         compartment_expected = "A"),
  # mC hypomethylation enriched in B? (TET-KO: DNMT3A depleted from heterochromatin)
  run_compartment_fisher(merged, "mC hypo -> B enriched",
                         positive_set = merged$mc_sig & merged$mc_diff < 0,
                         compartment_expected = "B"),
  # hmC loss enriched in A? (TET impairment -> hmC loss where TET is active)
  run_compartment_fisher(merged, "hmC hypo -> A enriched",
                         positive_set = merged$hmc_sig & merged$hmc_diff < 0,
                         compartment_expected = "A"),
  # hmC gain enriched in B? (less likely, but tests reciprocal)
  run_compartment_fisher(merged, "hmC hyper -> B enriched",
                         positive_set = merged$hmc_sig & merged$hmc_diff > 0,
                         compartment_expected = "B")
)

# =============================================================================
# STEP 3 (Task 2c): DIFFERENTIAL COMPARTMENT SHIFT x DMR OVERLAP
# =============================================================================

cat("================================================================================\n")
cat("STEP 3: Compartment shift x DMR direction overlap\n")
cat("================================================================================\n\n")

# Only genes with shifts
shifted_genes <- merged %>% dplyr::filter(shift != "Stable")

if (nrow(shifted_genes) > 0) {
  cat(sprintf("  Genes at shifted compartments: %d\n", nrow(shifted_genes)))
  cat(sprintf("    B to A: %d, A to B: %d\n",
              sum(shifted_genes$shift == "B to A"),
              sum(shifted_genes$shift == "A to B")))

  # Fisher's: B->A genes enriched for mC hypermethylation?
  shift_fisher_results <- dplyr::bind_rows(
    run_compartment_fisher(
      df = merged %>% dplyr::mutate(compartment = factor(ifelse(shift == "B to A", "A", "B"),
                                                          levels = c("A", "B"))),
      test_name = "B->A shift -> mC hyper enriched",
      positive_set = merged$mc_sig & merged$mc_diff > 0,
      compartment_expected = "A"  # B->A bins
    ),
    run_compartment_fisher(
      df = merged %>% dplyr::mutate(compartment = factor(ifelse(shift == "A to B", "A", "B"),
                                                          levels = c("A", "B"))),
      test_name = "A->B shift -> mC hypo enriched",
      positive_set = merged$mc_sig & merged$mc_diff < 0,
      compartment_expected = "A"  # A->B bins
    )
  )

  # O/E ratios for shift x DMR direction crosstab
  cat("\n  Observed/Expected ratios (shift x DMR direction):\n")
  shift_crosstab <- merged %>%
    dplyr::filter(mc_sig) %>%
    dplyr::group_by(shift, mc_direction) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

  total_sig <- sum(shift_crosstab$n)
  shift_totals <- merged %>%
    dplyr::filter(mc_sig) %>%
    dplyr::count(shift, name = "shift_total")
  dir_totals <- merged %>%
    dplyr::filter(mc_sig) %>%
    dplyr::count(mc_direction, name = "dir_total")

  shift_crosstab <- shift_crosstab %>%
    dplyr::left_join(shift_totals, by = "shift") %>%
    dplyr::left_join(dir_totals, by = "mc_direction") %>%
    dplyr::mutate(
      expected = (shift_total * dir_total) / total_sig,
      oe_ratio = n / expected
    )

  for (i in seq_len(nrow(shift_crosstab))) {
    row <- shift_crosstab[i, ]
    cat(sprintf("    %s x %s: O=%d, E=%.1f, O/E=%.2f\n",
                row$shift, row$mc_direction, row$n, row$expected, row$oe_ratio))
  }

  fisher_results <- dplyr::bind_rows(fisher_results, shift_fisher_results)
} else {
  cat("  No genes at shifted compartments — skipping shift analysis.\n")
  shift_crosstab <- data.frame()
}

# =============================================================================
# STEP 4 (Task 2d): VISUALIZATIONS
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("STEP 4: Generating visualizations\n")
cat("================================================================================\n\n")

# --- 29a: mC mod_difference violin by A vs B ---
cat("  29a: mC violin by compartment...\n")

wilcox_mc <- wilcox.test(mc_diff ~ compartment, data = merged)

p_29a <- ggplot(merged, aes(x = compartment, y = mc_diff, fill = compartment)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = COMPARTMENT_COLORS) +
  labs(
    title = "5mC Change by A/B Compartment",
    subtitle = sprintf("Wilcoxon %s | A: n=%d, B: n=%d",
                       fmt_p(wilcox_mc$p.value),
                       sum(merged$compartment == "A"),
                       sum(merged$compartment == "B")),
    x = "Compartment (control)", y = "mC mod_difference (mutant - control)",
    fill = "Compartment"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_29a, file.path(SECTION_DIR, "29a_mc_violin_by_compartment"),
                        width = 8, height = 7)

# --- 29b: hmC mod_difference violin by A vs B ---
cat("  29b: hmC violin by compartment...\n")

wilcox_hmc <- wilcox.test(hmc_diff ~ compartment, data = merged)

p_29b <- ggplot(merged, aes(x = compartment, y = hmc_diff, fill = compartment)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = COMPARTMENT_COLORS) +
  labs(
    title = "5hmC Change by A/B Compartment",
    subtitle = sprintf("Wilcoxon %s | A: n=%d, B: n=%d",
                       fmt_p(wilcox_hmc$p.value),
                       sum(merged$compartment == "A"),
                       sum(merged$compartment == "B")),
    x = "Compartment (control)", y = "hmC mod_difference (mutant - control)",
    fill = "Compartment"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_29b, file.path(SECTION_DIR, "29b_hmc_violin_by_compartment"),
                        width = 8, height = 7)

# --- 29c: mC violin by compartment shift ---
cat("  29c: mC violin by compartment shift...\n")

# Pairwise Wilcoxon tests for shift groups
shift_pw_mc <- pairwise.wilcox.test(merged$mc_diff, merged$shift, p.adjust.method = "BH")

p_29c <- ggplot(merged, aes(x = shift, y = mc_diff, fill = shift)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = SHIFT_COLORS) +
  labs(
    title = "5mC Change by Compartment Shift",
    subtitle = sprintf("B->A vs Stable: %s | A->B vs Stable: %s",
                       fmt_p(shift_pw_mc$p.value["Stable", "B to A"]),
                       fmt_p(shift_pw_mc$p.value["Stable", "A to B"])),
    x = "Compartment shift", y = "mC mod_difference (mutant - control)",
    fill = "Shift"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_29c, file.path(SECTION_DIR, "29c_mc_violin_by_shift"),
                        width = 10, height = 7)

# --- 29d: hmC violin by compartment shift ---
cat("  29d: hmC violin by compartment shift...\n")

shift_pw_hmc <- pairwise.wilcox.test(merged$hmc_diff, merged$shift, p.adjust.method = "BH")

p_29d <- ggplot(merged, aes(x = shift, y = hmc_diff, fill = shift)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = SHIFT_COLORS) +
  labs(
    title = "5hmC Change by Compartment Shift",
    subtitle = sprintf("B->A vs Stable: %s | A->B vs Stable: %s",
                       fmt_p(shift_pw_hmc$p.value["Stable", "B to A"]),
                       fmt_p(shift_pw_hmc$p.value["Stable", "A to B"])),
    x = "Compartment shift", y = "hmC mod_difference (mutant - control)",
    fill = "Shift"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_29d, file.path(SECTION_DIR, "29d_hmc_violin_by_shift"),
                        width = 10, height = 7)

# --- 29e: Stacked bar — DMR direction proportions ---
cat("  29e: Stacked bar — DMR direction proportions...\n")

# By compartment (A/B)
bar_comp <- merged %>%
  dplyr::filter(mc_sig) %>%
  dplyr::count(compartment, mc_direction) %>%
  dplyr::group_by(compartment) %>%
  dplyr::mutate(pct = 100 * n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(panel = "By Compartment",
                x_label = as.character(compartment))

# By shift (B->A, A->B, Stable)
bar_shift <- merged %>%
  dplyr::filter(mc_sig) %>%
  dplyr::count(shift, mc_direction) %>%
  dplyr::group_by(shift) %>%
  dplyr::mutate(pct = 100 * n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(panel = "By Shift",
                x_label = as.character(shift))

bar_data <- dplyr::bind_rows(bar_comp, bar_shift)
bar_data$x_label <- factor(bar_data$x_label,
                            levels = c("A", "B", "B to A", "A to B", "Stable"))
bar_data$panel <- factor(bar_data$panel, levels = c("By Compartment", "By Shift"))

p_29e <- ggplot(bar_data, aes(x = x_label, y = pct, fill = mc_direction)) +
  geom_col(position = "stack", width = 0.7) +
  geom_text(aes(label = sprintf("n=%d", n)),
            position = position_stack(vjust = 0.5), size = 3) +
  facet_wrap(~panel, scales = "free_x") +
  scale_fill_manual(values = COLORS$direction) +
  labs(
    title = "mC DMR Direction Proportions by Compartment and Shift",
    subtitle = "Significant DMRs only (q < 0.05)",
    x = NULL, y = "Percentage of significant DMRs",
    fill = "DMR direction"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_29e, file.path(SECTION_DIR, "29e_dmr_direction_stacked_bar"),
                        width = 10, height = 7)

# --- 29f: Scatter — PC1 vs mC mod_difference ---
cat("  29f: Scatter — PC1 vs mC mod_difference...\n")

# Spearman correlation
spearman_res <- cor.test(merged$mean_ctrl_pc1, merged$mc_diff, method = "spearman",
                         exact = FALSE)

merged$mc_sig_label <- ifelse(merged$mc_sig, "Significant", "Not significant")

p_29f <- ggplot(merged, aes(x = mean_ctrl_pc1, y = mc_diff)) +
  geom_point(aes(color = mc_sig_label), alpha = 0.15, size = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "#333333", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("Significant" = "#E41A1C", "Not significant" = "grey70")) +
  annotate("text", x = max(merged$mean_ctrl_pc1) * 0.7,
           y = max(merged$mc_diff, na.rm = TRUE) * 0.9,
           label = sprintf("rho = %.3f\n%s",
                           spearman_res$estimate, fmt_p(spearman_res$p.value)),
           hjust = 0, size = 4, fontface = "italic") +
  annotate("text", x = max(merged$mean_ctrl_pc1) * 0.8, y = min(merged$mc_diff) * 0.5,
           label = "A compartment", color = COMPARTMENT_COLORS["A"],
           fontface = "bold", size = 4) +
  annotate("text", x = min(merged$mean_ctrl_pc1) * 0.8, y = min(merged$mc_diff) * 0.5,
           label = "B compartment", color = COMPARTMENT_COLORS["B"],
           fontface = "bold", size = 4) +
  labs(
    title = "Control PC1 vs 5mC Change (BAP1-KO - WT)",
    subtitle = sprintf("Spearman rho = %.3f, %s | n = %d genes",
                       spearman_res$estimate, fmt_p(spearman_res$p.value),
                       nrow(merged)),
    x = "Mean control PC1 (A/B compartment score)",
    y = "mC mod_difference (mutant - control)",
    color = NULL
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_29f, file.path(SECTION_DIR, "29f_pc1_vs_mc_scatter"),
                        width = 10, height = 9)

# --- 29g: Composite summary panel ---
cat("  29g: Composite summary panel...\n")

p_29g <- (p_29a + p_29b + p_29e) +
  plot_layout(widths = c(1, 1, 1.5)) +
  plot_annotation(
    title = "A/B Compartment Methylation Mapping: BAP1-KO vs TET-KO Prediction",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

save_multiformat_ggplot(p_29g, file.path(SECTION_DIR, "29g_composite_compartment_summary"),
                        width = 16, height = 6)

# =============================================================================
# STEP 5: EXPORT TABLES
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("STEP 5: Exporting tables\n")
cat("================================================================================\n\n")

# 1. Full gene assignment table
gene_export <- merged %>%
  dplyr::select(gene_name, chr, start, end, compartment, shift,
                mean_ctrl_pc1, difference, adj_pvalue, dist_to_tss,
                mc_diff, mc_q, mc_sig, mc_direction,
                hmc_diff, hmc_q, hmc_sig, hmc_direction)

write.table(gene_export,
            file.path(TABLES_DIR, "compartment_gene_assignment.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: compartment_gene_assignment.tsv (%d genes)\n", nrow(gene_export)))

# 2. Fisher's test results
write.table(fisher_results,
            file.path(TABLES_DIR, "compartment_fisher_tests.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: compartment_fisher_tests.tsv (%d tests)\n", nrow(fisher_results)))

# 3. Summary statistics per compartment/shift
summary_stats <- merged %>%
  dplyr::group_by(compartment) %>%
  dplyr::summarise(
    n_genes = dplyr::n(),
    mean_mc_diff = mean(mc_diff, na.rm = TRUE),
    median_mc_diff = median(mc_diff, na.rm = TRUE),
    mean_hmc_diff = mean(hmc_diff, na.rm = TRUE),
    median_hmc_diff = median(hmc_diff, na.rm = TRUE),
    n_mc_sig = sum(mc_sig),
    n_mc_hyper = sum(mc_sig & mc_diff > 0),
    n_mc_hypo = sum(mc_sig & mc_diff < 0),
    pct_mc_hyper = 100 * sum(mc_sig & mc_diff > 0) / sum(mc_sig),
    n_hmc_sig = sum(hmc_sig),
    n_hmc_hypo = sum(hmc_sig & hmc_diff < 0),
    .groups = "drop"
  ) %>%
  dplyr::mutate(grouping = "compartment", label = as.character(compartment))

shift_stats <- merged %>%
  dplyr::group_by(shift) %>%
  dplyr::summarise(
    n_genes = dplyr::n(),
    mean_mc_diff = mean(mc_diff, na.rm = TRUE),
    median_mc_diff = median(mc_diff, na.rm = TRUE),
    mean_hmc_diff = mean(hmc_diff, na.rm = TRUE),
    median_hmc_diff = median(hmc_diff, na.rm = TRUE),
    n_mc_sig = sum(mc_sig),
    n_mc_hyper = sum(mc_sig & mc_diff > 0),
    n_mc_hypo = sum(mc_sig & mc_diff < 0),
    pct_mc_hyper = ifelse(sum(mc_sig) > 0, 100 * sum(mc_sig & mc_diff > 0) / sum(mc_sig), NA),
    n_hmc_sig = sum(hmc_sig),
    n_hmc_hypo = sum(hmc_sig & hmc_diff < 0),
    .groups = "drop"
  ) %>%
  dplyr::mutate(grouping = "shift", label = as.character(shift))

all_summary <- dplyr::bind_rows(summary_stats, shift_stats)
write.table(all_summary,
            file.path(TABLES_DIR, "compartment_methylation_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: compartment_methylation_summary.tsv (%d rows)\n", nrow(all_summary)))

# =============================================================================
# STEP 6: SUMMARY PRINTOUT
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("SECTION 29 SUMMARY: A/B COMPARTMENT METHYLATION MAPPING\n")
cat("================================================================================\n\n")

cat(sprintf("Input: %d HOMER 25kb bins -> %d deduplicated genes -> %d matched to DMRs\n\n",
            nrow(comp_raw), n_comp_genes, n_overlap))

cat("--- Compartment Distribution ---\n")
cat(sprintf("  A compartment: %d genes (%.1f%%)\n",
            sum(merged$compartment == "A"), 100 * mean(merged$compartment == "A")))
cat(sprintf("  B compartment: %d genes (%.1f%%)\n",
            sum(merged$compartment == "B"), 100 * mean(merged$compartment == "B")))

cat("\n--- Key Fisher's Test Results ---\n")
for (i in seq_len(min(4, nrow(fisher_results)))) {
  row <- fisher_results[i, ]
  cat(sprintf("  %s: OR = %.2f (%.2f-%.2f), %s\n",
              row$test, row$odds_ratio, row$ci_lower, row$ci_upper,
              fmt_p(row$p_value)))
}

cat("\n--- Wilcoxon Tests (mC/hmC difference by compartment) ---\n")
cat(sprintf("  mC A vs B: %s (median A: %.4f, B: %.4f)\n",
            fmt_p(wilcox_mc$p.value),
            median(merged$mc_diff[merged$compartment == "A"]),
            median(merged$mc_diff[merged$compartment == "B"])))
cat(sprintf("  hmC A vs B: %s (median A: %.4f, B: %.4f)\n",
            fmt_p(wilcox_hmc$p.value),
            median(merged$hmc_diff[merged$compartment == "A"]),
            median(merged$hmc_diff[merged$compartment == "B"])))

cat(sprintf("\n--- PC1 x mC Correlation ---\n"))
cat(sprintf("  Spearman rho = %.3f, %s\n",
            spearman_res$estimate, fmt_p(spearman_res$p.value)))

# Biological interpretation
cat("\n--- Biological Interpretation ---\n")
mc_hyper_a_fisher <- fisher_results[fisher_results$test == "mC hyper -> A enriched", ]
mc_hypo_b_fisher <- fisher_results[fisher_results$test == "mC hypo -> B enriched", ]

tet_ko_match <- (mc_hyper_a_fisher$odds_ratio > 1 & mc_hyper_a_fisher$p_value < 0.05) ||
                (mc_hypo_b_fisher$odds_ratio > 1 & mc_hypo_b_fisher$p_value < 0.05)

if (tet_ko_match) {
  cat("  RESULT: BAP1-KO methylation changes show A/B compartment asymmetry\n")
  cat("  consistent with the Lopez-Moyado TET-KO signature.\n")
  cat("  This supports CONVERGENT MECHANISMS: BAP1 loss -> TET impairment ->\n")
  cat("  DNMT3A redistribution from B (heterochromatin) to A (euchromatin).\n")
} else {
  cat("  RESULT: BAP1-KO methylation changes do NOT show the expected TET-KO\n")
  cat("  compartment asymmetry pattern.\n")
  cat("  This suggests a DISTINCT BAP1-mediated pathway that does not involve\n")
  cat("  the same DNMT3A redistribution mechanism as TET-KO.\n")
}

cat(sprintf("\n  Plots saved to: %s\n", SECTION_DIR))
cat(sprintf("  Tables saved to: %s\n", TABLES_DIR))
cat("\nSection 29 complete.\n")
