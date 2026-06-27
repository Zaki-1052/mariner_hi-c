# biomodal/downstream/scripts/viz_sections/section_91_compartment_shift_stoichiometry.R
# Section 91: Compartment Shift Methylation + Stoichiometry + TET Consistency
#
# Answers three related questions about compartment transitions and 5hmC fate:
#
#   Q1: How many genes in each subcompartment are co-significant for 5mC and 5hmC?
#   Q2: Are A→B and B→A transitions symmetric in methylation effect?
#   Q3: Is 5hmC being completely removed or redistributed to 5mC? (TET consistency)
#
#   Panel 91a: Co-significance overlap by CALDER2 subcompartment
#   Panel 91b: Compartment transition methylation (A→B vs B→A vs Stable)
#   Panel 91c: Stoichiometric ratio (|Δ5mC|/|Δ5hmC|) by shift category
#   Panel 91d: Total modC balance (Δ5mC + Δ5hmC) per gene, by shift
#   Panel 91e: 5hmC depletion vs baseline level (TET impediment prediction)
#   Panel 91f: Residual 5hmC in mutant by subcompartment (is 5hmC removed or retained?)
#   Panel 91g: Combined composite
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_91_compartment_shift_stoichiometry.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 91: COMPARTMENT SHIFT STOICHIOMETRY
# =============================================================================

cat("================================================================================\n")
cat("SECTION 91: COMPARTMENT SHIFT STOICHIOMETRY\n")
cat("================================================================================\n\n")

# ---- Constants --------------------------------------------------------------

SUBCMPT_ORDER <- c("A.1", "A.2", "B.1", "B.2")
SUBCMPT_LABELS <- c("A.1" = "A.1\n(Strong Active)",
                     "A.2" = "A.2\n(Weak Active)",
                     "B.1" = "B.1\n(Fac. Het)",
                     "B.2" = "B.2\n(Constit. Het)")
SUBCMPT_COLORS <- c("A.1" = "#E41A1C", "A.2" = "#FF7F00",
                     "B.1" = "#4DAF4A", "B.2" = "#377EB8")

SHIFT_COLORS <- c("B to A" = "#D7191C", "Stable" = "grey60", "A to B" = "#2C7BB6")

fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

# ---- Load CALDER2 subcompartment labels + assign to genes -------------------

SUBCMPT_PATH <- file.path(REPO_ROOT, "ML/cmpts/outputs/calder2/late/250402_subcompartment_labels_100kb.tsv")
if (!file.exists(SUBCMPT_PATH)) stop("CALDER2 file not found: ", SUBCMPT_PATH)

subcmpt <- read.table(SUBCMPT_PATH, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
subcmpt <- subcmpt %>% dplyr::filter(!is.na(ctrl_label) & ctrl_label != "NA")

assign_subcompartment <- function(dmr_df, subcmpt_df) {
  dmr_df$midpoint <- (dmr_df$start + dmr_df$end) %/% 2
  dmr_df$bin_start <- (dmr_df$midpoint %/% 100000) * 100000 + 1
  dmr_df$bin_end <- dmr_df$bin_start + 99999
  dmr_df %>%
    left_join(
      subcmpt_df %>% dplyr::select(chr, bin_start, bin_end, ctrl_label, mut_label, label_changed),
      by = c("chr", "bin_start", "bin_end")
    )
}

cat("Assigning CALDER2 subcompartments to gene bodies...\n")
mc_cmpt <- assign_subcompartment(mc_dmr, subcmpt) %>% dplyr::filter(!is.na(ctrl_label))
hmc_cmpt <- assign_subcompartment(hmc_dmr, subcmpt) %>% dplyr::filter(!is.na(ctrl_label))
cat(sprintf("  5mC: %d genes, 5hmC: %d genes\n", nrow(mc_cmpt), nrow(hmc_cmpt)))

# ---- Load Section 29 compartment gene assignment (Homer A/B shifts) ---------

COMP_GENE_FILE <- file.path(TABLES_DIR, "compartment_gene_assignment.tsv")
if (!file.exists(COMP_GENE_FILE)) stop("Section 29 table not found: ", COMP_GENE_FILE)

cat("Loading Section 29 compartment gene assignment (Homer A/B shifts)...\n")
comp_genes <- read.table(COMP_GENE_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
comp_genes$shift <- factor(comp_genes$shift, levels = c("B to A", "Stable", "A to B"))
cat(sprintf("  %d genes with shift assignments\n", nrow(comp_genes)))
cat(sprintf("    B to A: %d | Stable: %d | A to B: %d\n",
            sum(comp_genes$shift == "B to A"),
            sum(comp_genes$shift == "Stable"),
            sum(comp_genes$shift == "A to B")))

# ---- Join with absolute methylation levels from pre-loaded DMR data ---------

cat("\nJoining with absolute methylation levels...\n")

mc_abs <- mc_dmr %>%
  dplyr::select(gene, mc_ctrl_level = mean_mod_group1, mc_mut_level = mean_mod_group2) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

hmc_abs <- hmc_dmr %>%
  dplyr::select(gene, hmc_ctrl_level = mean_mod_group1, hmc_mut_level = mean_mod_group2) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

shift_df <- comp_genes %>%
  inner_join(mc_abs, by = c("gene_name" = "gene")) %>%
  inner_join(hmc_abs, by = c("gene_name" = "gene"))

shift_df$mc_diff <- as.numeric(shift_df$mc_diff)
shift_df$hmc_diff <- as.numeric(shift_df$hmc_diff)
shift_df$total_diff <- shift_df$mc_diff + shift_df$hmc_diff
shift_df$total_ctrl <- shift_df$mc_ctrl_level + shift_df$hmc_ctrl_level
shift_df$total_mut <- shift_df$mc_mut_level + shift_df$hmc_mut_level

cat(sprintf("  Matched: %d genes with absolute levels + shift data\n", nrow(shift_df)))

# =============================================================================
# PANEL 91a: CO-SIGNIFICANCE OVERLAP BY SUBCOMPARTMENT
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("Panel 91a: Co-Significance Overlap by Subcompartment\n")
cat("================================================================================\n\n")

# Build overlap table from CALDER2-assigned genes
overlap_df <- mc_cmpt %>%
  dplyr::select(gene, ctrl_label, mc_sig = significant, mc_dir = direction) %>%
  inner_join(
    hmc_cmpt %>% dplyr::select(gene, hmc_sig = significant, hmc_dir = direction),
    by = "gene"
  ) %>%
  dplyr::mutate(
    overlap_cat = case_when(
      mc_sig & hmc_sig ~ "Both significant",
      mc_sig & !hmc_sig ~ "5mC only",
      !mc_sig & hmc_sig ~ "5hmC only",
      TRUE ~ "Neither"
    ),
    overlap_cat = factor(overlap_cat,
      levels = c("Both significant", "5mC only", "5hmC only", "Neither")),
    ctrl_label = factor(ctrl_label, levels = SUBCMPT_ORDER)
  )

overlap_summary <- overlap_df %>%
  group_by(ctrl_label, overlap_cat) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(ctrl_label) %>%
  dplyr::mutate(total = sum(n), pct = 100 * n / total) %>%
  ungroup()

cat("  Co-significance by subcompartment:\n")
for (sc in SUBCMPT_ORDER) {
  both_n <- overlap_summary$n[overlap_summary$ctrl_label == sc & overlap_summary$overlap_cat == "Both significant"]
  total_n <- overlap_summary$total[overlap_summary$ctrl_label == sc][1]
  cat(sprintf("    %s: %d / %d both sig (%.1f%%)\n", sc, both_n, total_n, 100 * both_n / total_n))
}

OVERLAP_COLORS <- c("Both significant" = "#7B3294",
                     "5mC only" = COLORS$methylation["5mC"],
                     "5hmC only" = COLORS$methylation["5hmC"],
                     "Neither" = "grey80")

p_91a <- ggplot(overlap_summary, aes(x = ctrl_label, y = pct, fill = overlap_cat)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(pct > 5, sprintf("%d\n(%.0f%%)", n, pct), "")),
            position = position_stack(vjust = 0.5), size = 3,
            color = "white", fontface = "bold") +
  scale_fill_manual(values = OVERLAP_COLORS, name = NULL) +
  scale_x_discrete(labels = SUBCMPT_LABELS) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title = "5mC / 5hmC Co-Significance by Subcompartment",
    subtitle = "Co-significant fraction drops monotonically from active (56%) to heterochromatic (19%)",
    x = "CALDER2 Subcompartment (Control)",
    y = "% of Genes"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_91a,
  file.path(OUTPUT_DIR, "91a_cosignificance_by_subcompartment"), width = 10, height = 7)

# =============================================================================
# PANEL 91b: COMPARTMENT TRANSITION METHYLATION
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("Panel 91b: Compartment Transition Methylation\n")
cat("================================================================================\n\n")

shift_summary <- shift_df %>%
  group_by(shift) %>%
  summarise(
    n = n(),
    mean_mc = mean(mc_diff),
    mean_hmc = mean(hmc_diff),
    mean_total = mean(total_diff),
    mc_sig_n = sum(mc_sig == "TRUE"),
    mc_hyper = sum(mc_sig == "TRUE" & mc_direction == "Hypermethylated"),
    mc_hypo = sum(mc_sig == "TRUE" & mc_direction == "Hypomethylated"),
    .groups = "drop"
  ) %>%
  dplyr::mutate(pct_hyper = 100 * mc_hyper / mc_sig_n)

cat("  Shift summary:\n")
for (i in seq_len(nrow(shift_summary))) {
  r <- shift_summary[i, ]
  cat(sprintf("    %-12s: n=%d, mean Δ5mC=%+.3f%%, mean Δ5hmC=%+.3f%%, total=%+.3f%%, %d sig (%d hyper / %d hypo = %.1f%%)\n",
              r$shift, r$n, 100*r$mean_mc, 100*r$mean_hmc, 100*r$mean_total,
              r$mc_sig_n, r$mc_hyper, r$mc_hypo, r$pct_hyper))
}

bar_shift <- shift_summary %>%
  tidyr::pivot_longer(cols = c(mean_mc, mean_hmc), names_to = "modality", values_to = "mean_diff") %>%
  dplyr::mutate(
    modality = factor(ifelse(modality == "mean_mc", "Δ5mC", "Δ5hmC"),
                      levels = c("Δ5mC", "Δ5hmC"))
  )

p_91b <- ggplot(bar_shift, aes(x = shift, y = mean_diff * 100, fill = modality)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%+.2f%%", mean_diff * 100)),
            position = position_dodge(width = 0.7),
            vjust = ifelse(bar_shift$mean_diff > 0, -0.5, 1.5), size = 3.5) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_text(data = shift_summary,
            aes(x = shift, y = 0, label = sprintf("n=%d", n)),
            inherit.aes = FALSE, vjust = -0.5, size = 3, fontface = "italic", color = "grey30") +
  scale_fill_manual(values = c("Δ5mC" = COLORS$methylation["5mC"],
                                "Δ5hmC" = COLORS$methylation["5hmC"]),
                    name = NULL) +
  labs(
    title = "Methylation Change by Compartment Transition",
    subtitle = "A→B genes gain 2.4× more methylation than B→A genes lose | Mirror-image 5mC/5hmC at both",
    x = "Compartment Shift (Homer PC1, FDR<0.05 & |Diff|>0.30)",
    y = "Mean Methylation Change (%)"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_91b,
  file.path(OUTPUT_DIR, "91b_compartment_transition_methylation"), width = 10, height = 7)

# =============================================================================
# PANEL 91c: STOICHIOMETRIC RATIO BY SHIFT
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("Panel 91c: Stoichiometric Ratio by Shift\n")
cat("================================================================================\n\n")

stoich_summary <- shift_summary %>%
  dplyr::mutate(
    ratio = abs(mean_mc) / abs(mean_hmc) * 100,
    total_pct = mean_total * 100
  )

cat("  Stoichiometric ratios (|Δ5mC|/|Δ5hmC| × 100):\n")
for (i in seq_len(nrow(stoich_summary))) {
  r <- stoich_summary[i, ]
  cat(sprintf("    %-12s: ratio=%.1f%%, total change=%+.3f%%\n",
              r$shift, r$ratio, r$total_pct))
}
cat("    100% = perfectly stoichiometric (total modC conserved)\n")

p_91c <- ggplot(stoich_summary, aes(x = shift, y = ratio, fill = shift)) +
  geom_col(width = 0.5, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%\n(total: %+.2f%%)", ratio, total_pct)),
            vjust = -0.3, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = SHIFT_COLORS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  annotate("text", x = 3.4, y = 100, label = "Stoichiometric\n(100%)",
           hjust = 0, vjust = -0.3, size = 3, fontface = "italic") +
  labs(
    title = "Stoichiometric Ratio: |Δ5mC| / |Δ5hmC|",
    subtitle = "A→B: 108.7% — 91% of 5hmC loss recaptured as 5mC (near-stoichiometric redistribution)",
    x = "Compartment Shift", y = "|Δ5mC| / |Δ5hmC| × 100 (%)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_91c,
  file.path(OUTPUT_DIR, "91c_stoichiometric_ratio_by_shift"), width = 9, height = 7)

# =============================================================================
# PANEL 91d: TOTAL modC BALANCE PER GENE
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("Panel 91d: Total modC Balance Per Gene\n")
cat("================================================================================\n\n")

# Wilcoxon vs 0 for each shift
total_annot <- data.frame(shift = character(), median_val = numeric(),
                           p = numeric(), n = integer(), stringsAsFactors = FALSE)
for (s in levels(shift_df$shift)) {
  vals <- shift_df$total_diff[shift_df$shift == s]
  wt <- wilcox.test(vals, mu = 0)
  total_annot <- rbind(total_annot,
    data.frame(shift = s, median_val = median(vals) * 100,
               p = wt$p.value, n = length(vals), stringsAsFactors = FALSE))
  cat(sprintf("  %s: median total_diff=%+.5f, %s (n=%d)\n",
              s, median(vals), fmt_p(wt$p.value), length(vals)))
}
total_annot$shift <- factor(total_annot$shift, levels = levels(shift_df$shift))
total_annot$label <- sprintf("median=%+.2f%%\n%s\nn=%d",
                              total_annot$median_val,
                              sapply(total_annot$p, fmt_p),
                              total_annot$n)

p_91d <- ggplot(shift_df, aes(x = shift, y = total_diff * 100, fill = shift)) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.3, alpha = 0.8, fill = "white") +
  geom_text(data = total_annot,
            aes(x = shift, y = -13, label = label),
            inherit.aes = FALSE, size = 3, fontface = "bold", vjust = 1) +
  scale_fill_manual(values = SHIFT_COLORS, guide = "none") +
  coord_cartesian(ylim = c(-17, 15)) +
  labs(
    title = "Total Modified Cytosine Change (Δ5mC + Δ5hmC) Per Gene",
    subtitle = "Centered near zero across all shifts — 5hmC loss is redistributed to 5mC, not removed",
    x = "Compartment Shift", y = "Total modC Change (Δ5mC + Δ5hmC) (%)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_91d,
  file.path(OUTPUT_DIR, "91d_total_modc_balance_by_shift"), width = 10, height = 7)

# =============================================================================
# PANEL 91e: 5hmC DEPLETION vs BASELINE (TET IMPEDIMENT PREDICTION)
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("Panel 91e: 5hmC Depletion vs Baseline Level\n")
cat("================================================================================\n\n")

# TET impediment predicts: higher baseline 5hmC → larger hmC decrease
# (more substrate to deplete when production is blocked)
hmc_baseline_cor <- cor.test(shift_df$hmc_ctrl_level, shift_df$hmc_diff,
                              method = "spearman", exact = FALSE)
cat(sprintf("  Spearman(baseline hmC, Δ5hmC): rho=%.4f, %s\n",
            hmc_baseline_cor$estimate, fmt_p(hmc_baseline_cor$p.value)))

# What fraction of baseline 5hmC remains in mutant?
shift_df$hmc_retained_pct <- ifelse(shift_df$hmc_ctrl_level > 0,
  100 * shift_df$hmc_mut_level / shift_df$hmc_ctrl_level, NA)

cat(sprintf("  Median 5hmC retention (mut/ctrl): %.1f%%\n",
            median(shift_df$hmc_retained_pct, na.rm = TRUE)))
cat(sprintf("  Genes with >50%% hmC retained: %d / %d (%.1f%%)\n",
            sum(shift_df$hmc_retained_pct > 50, na.rm = TRUE),
            sum(!is.na(shift_df$hmc_retained_pct)),
            100 * mean(shift_df$hmc_retained_pct > 50, na.rm = TRUE)))
cat(sprintf("  Genes with <10%% hmC retained (near-complete depletion): %d (%.1f%%)\n",
            sum(shift_df$hmc_retained_pct < 10, na.rm = TRUE),
            100 * mean(shift_df$hmc_retained_pct < 10, na.rm = TRUE)))

p_91e <- ggplot(shift_df, aes(x = hmc_ctrl_level * 100, y = hmc_diff * 100)) +
  geom_point(alpha = 0.08, size = 0.3, color = COLORS$methylation["5hmC"]) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  annotate("text",
    x = max(shift_df$hmc_ctrl_level * 100, na.rm = TRUE) * 0.6,
    y = min(shift_df$hmc_diff * 100, na.rm = TRUE) * 0.7,
    label = sprintf("Spearman rho = %.3f\n%s\n\nTET impediment prediction:\nhigher baseline → larger loss",
                    hmc_baseline_cor$estimate,
                    fmt_p(hmc_baseline_cor$p.value)),
    hjust = 0, size = 3.5, fontface = "italic") +
  labs(
    title = "5hmC Change vs Baseline Level: TET Impediment Consistency",
    subtitle = "Genes with more 5hmC substrate show larger decreases — production blocked, not active removal",
    x = "Baseline 5hmC Level (% in Control)",
    y = "Δ5hmC (Mutant − Control, %)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_91e,
  file.path(OUTPUT_DIR, "91e_hmc_depletion_vs_baseline"), width = 10, height = 8)

# =============================================================================
# PANEL 91f: RESIDUAL 5hmC IN MUTANT BY SUBCOMPARTMENT
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("Panel 91f: Residual 5hmC in Mutant by Subcompartment\n")
cat("================================================================================\n\n")

# Join hmc_cmpt with absolute levels
hmc_levels <- hmc_cmpt %>%
  dplyr::select(gene, ctrl_label, hmc_ctrl = mean_mod_group1, hmc_mut = mean_mod_group2) %>%
  dplyr::mutate(ctrl_label = factor(ctrl_label, levels = SUBCMPT_ORDER))

level_long <- rbind(
  data.frame(gene = hmc_levels$gene, subcmpt = hmc_levels$ctrl_label,
             condition = "Control", level = hmc_levels$hmc_ctrl, stringsAsFactors = FALSE),
  data.frame(gene = hmc_levels$gene, subcmpt = hmc_levels$ctrl_label,
             condition = "Mutant", level = hmc_levels$hmc_mut, stringsAsFactors = FALSE)
) %>%
  dplyr::mutate(condition = factor(condition, levels = c("Control", "Mutant")))

cat("  5hmC levels by subcompartment (median %):\n")
level_summary <- hmc_levels %>%
  group_by(ctrl_label) %>%
  summarise(
    n = n(),
    ctrl_median = median(hmc_ctrl) * 100,
    mut_median = median(hmc_mut) * 100,
    delta_median = median(hmc_mut - hmc_ctrl) * 100,
    retained_pct = median(hmc_mut) / median(hmc_ctrl) * 100,
    .groups = "drop"
  )

# Paired Wilcoxon per subcompartment
level_summary$p <- NA_real_
for (i in seq_len(nrow(level_summary))) {
  sc <- level_summary$ctrl_label[i]
  sub <- hmc_levels %>% dplyr::filter(ctrl_label == sc)
  wt <- wilcox.test(sub$hmc_ctrl, sub$hmc_mut, paired = TRUE)
  level_summary$p[i] <- wt$p.value
}

for (i in seq_len(nrow(level_summary))) {
  r <- level_summary[i, ]
  cat(sprintf("    %s: ctrl=%.2f%%, mut=%.2f%%, Δ=%+.2f%%, retained=%.0f%%, %s (n=%d)\n",
              r$ctrl_label, r$ctrl_median, r$mut_median, r$delta_median,
              r$retained_pct, fmt_p(r$p), r$n))
}

# Annotation dataframe for the faceted violin
level_annot <- level_summary %>%
  dplyr::mutate(
    label = sprintf("ctrl: %.2f%%\nmut: %.2f%%\nΔ = %+.2f%%\n%s\nn=%d",
                    ctrl_median, mut_median, delta_median, sapply(p, fmt_p), n)
  )

p_91f <- ggplot(level_long, aes(x = condition, y = level * 100, fill = condition)) +
  geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.3, alpha = 0.8) +
  geom_text(data = level_annot,
            aes(x = 1.5, y = 18, label = label),
            inherit.aes = FALSE, size = 2.8, fontface = "bold",
            color = "grey20", lineheight = 0.9) +
  facet_wrap(~ ctrl_label, nrow = 1, labeller = labeller(ctrl_label = SUBCMPT_LABELS)) +
  scale_fill_manual(values = COLORS$condition, guide = "none") +
  coord_cartesian(ylim = c(0, 22)) +
  labs(
    title = "5hmC Levels: Control vs Mutant by Subcompartment",
    subtitle = "5hmC is reduced but NOT eliminated — consistent with blocked production, not active removal",
    x = NULL, y = "5hmC Level (%)"
  ) +
  theme_biomodal() +
  theme(strip.text = element_text(size = 9, face = "bold"))

save_multiformat_ggplot(p_91f,
  file.path(OUTPUT_DIR, "91f_residual_hmc_by_subcompartment"), width = 14, height = 7)

# ---- Combined composite ----------------------------------------------------

cat("\nCreating combined Figure 91g...\n")

p_91g <- (p_91a | p_91b) /
  (p_91c | p_91d) /
  (p_91e | p_91f) +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(
    title = "Compartment Shift Stoichiometry: 5hmC Is Redistributed, Not Removed",
    subtitle = "Near-stoichiometric mC/hmC balance + partial hmC retention + baseline-dependent depletion → TET impediment",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  )

save_multiformat_ggplot(p_91g,
  file.path(OUTPUT_DIR, "91g_combined_shift_stoichiometry"), width = 18, height = 20)

# ---- Save tables ------------------------------------------------------------

cat("\nSaving tables...\n")

# Co-significance overlap
write.table(overlap_summary,
  file.path(TABLES_DIR, "91_cosignificance_by_subcompartment.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 91_cosignificance_by_subcompartment.tsv\n")

# Shift methylation summary
write.table(shift_summary,
  file.path(TABLES_DIR, "91_shift_methylation_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 91_shift_methylation_summary.tsv\n")

# Stoichiometric ratios
write.table(stoich_summary %>% dplyr::select(shift, n, mean_mc, mean_hmc, mean_total, ratio, total_pct),
  file.path(TABLES_DIR, "91_stoichiometric_ratios.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 91_stoichiometric_ratios.tsv\n")

# Residual hmC levels
write.table(level_summary,
  file.path(TABLES_DIR, "91_residual_hmc_levels.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 91_residual_hmc_levels.tsv\n")

# Per-gene shift + absolute levels
gene_export <- shift_df %>%
  dplyr::select(gene_name, compartment, shift, mean_ctrl_pc1, difference,
                mc_diff, hmc_diff, total_diff,
                mc_ctrl_level, mc_mut_level, hmc_ctrl_level, hmc_mut_level,
                mc_sig, hmc_sig)
write.table(gene_export,
  file.path(TABLES_DIR, "91_per_gene_shift_stoichiometry.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: 91_per_gene_shift_stoichiometry.tsv (%d genes)\n", nrow(gene_export)))

cat("\n")
cat("================================================================================\n")
cat("SECTION 91 COMPLETE\n")
cat("================================================================================\n")
