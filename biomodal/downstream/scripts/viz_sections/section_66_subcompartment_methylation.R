# biomodal/downstream/scripts/viz_sections/section_66_subcompartment_methylation.R
# Section 66: Methylation Stratified by Chromatin Compartment
# Standalone script - sources shared config for all dependencies and data
#
# Uses CALDER2 subcompartment labels (A.1/A.2/B.1/B.2) as proxy for chromatin
# state: B.2 = constitutive het, B.1 = facultative het, A.2 = weak active,
# A.1 = strong active. Also overlays H3K27me3 and H3K27ac peaks directly.
#
#   Panel 66a: % significant 5mC DMRs by subcompartment (hyper vs hypo)
#   Panel 66b: Direction split (hyper vs hypo fraction) for sig DMRs
#   Panel 66c: 5mC and 5hmC levels by subcompartment, ctrl vs mut (violin)
#   Panel 66d: Same analysis stratified by H3K27me3 / H3K27ac histone marks
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_66_subcompartment_methylation.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 66: METHYLATION STRATIFIED BY CHROMATIN COMPARTMENT
# =============================================================================

cat("================================================================================\n")
cat("SECTION 66: METHYLATION STRATIFIED BY CHROMATIN COMPARTMENT\n")
cat("================================================================================\n\n")

# ---- Load CALDER2 subcompartment labels --------------------------------------

SUBCMPT_PATH <- file.path(REPO_ROOT, "ML/cmpts/outputs/calder2/late/250402_subcompartment_labels_100kb.tsv")

if (!file.exists(SUBCMPT_PATH)) {
  stop("CALDER2 subcompartment file not found: ", SUBCMPT_PATH)
}

cat("Loading CALDER2 subcompartment labels (late/adult timepoint)...\n")
subcmpt <- read.table(SUBCMPT_PATH, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
subcmpt <- subcmpt %>% dplyr::filter(!is.na(ctrl_label) & ctrl_label != "NA")
cat(sprintf("  %d bins with labels (A.1=%d, A.2=%d, B.1=%d, B.2=%d)\n",
            nrow(subcmpt),
            sum(subcmpt$ctrl_label == "A.1"),
            sum(subcmpt$ctrl_label == "A.2"),
            sum(subcmpt$ctrl_label == "B.1"),
            sum(subcmpt$ctrl_label == "B.2")))

SUBCMPT_ORDER <- c("A.1", "A.2", "B.1", "B.2")
SUBCMPT_LABELS <- c("A.1" = "A.1\n(Strong Active)",
                     "A.2" = "A.2\n(Weak Active)",
                     "B.1" = "B.1\n(Fac. Het)",
                     "B.2" = "B.2\n(Constit. Het)")
SUBCMPT_COLORS <- c("A.1" = "#E41A1C", "A.2" = "#FF7F00",
                     "B.1" = "#4DAF4A", "B.2" = "#377EB8")

# ---- Assign subcompartment to each gene body ---------------------------------

cat("\nAssigning subcompartment to gene bodies by midpoint overlap...\n")

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

mc_dmr_cmpt <- assign_subcompartment(mc_dmr, subcmpt)
hmc_dmr_cmpt <- assign_subcompartment(hmc_dmr, subcmpt)

n_assigned <- sum(!is.na(mc_dmr_cmpt$ctrl_label))
n_total <- nrow(mc_dmr_cmpt)
cat(sprintf("  5mC: %d / %d genes assigned (%.1f%%)\n",
            n_assigned, n_total, 100 * n_assigned / n_total))

mc_dmr_cmpt <- mc_dmr_cmpt %>% dplyr::filter(!is.na(ctrl_label))
hmc_dmr_cmpt <- hmc_dmr_cmpt %>% dplyr::filter(!is.na(ctrl_label))

# Summary
cat("\n  Per-subcompartment gene counts (5mC):\n")
cmpt_summary <- mc_dmr_cmpt %>%
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

for (i in seq_len(nrow(cmpt_summary))) {
  row <- cmpt_summary[i, ]
  cat(sprintf("    %s: %d total, %d sig (%.1f%%), %d hyper, %d hypo\n",
              row$ctrl_label, row$total, row$sig, row$pct_sig, row$hyper, row$hypo))
}

# ---- Panel 66a: % significant DMRs by subcompartment (grouped bar) ----------

cat("\nCreating Figure 66a: % significant DMRs by subcompartment...\n")

bar_df <- cmpt_summary %>%
  tidyr::pivot_longer(cols = c(hyper, hypo), names_to = "direction", values_to = "count") %>%
  dplyr::mutate(
    pct = 100 * count / total,
    direction = factor(ifelse(direction == "hyper", "Hypermethylated", "Hypomethylated"),
                       levels = c("Hypermethylated", "Hypomethylated")),
    subcmpt_label = SUBCMPT_LABELS[as.character(ctrl_label)]
  )

# Chi-squared test: subcompartment x significance
chisq_mat <- matrix(
  c(cmpt_summary$sig, cmpt_summary$total - cmpt_summary$sig),
  ncol = 2,
  dimnames = list(cmpt_summary$ctrl_label, c("Sig", "NotSig"))
)
chisq_result <- chisq.test(chisq_mat)
cat(sprintf("  Chi-squared test (subcompartment x significance): X²=%.1f, df=%d, p=%.2e\n",
            chisq_result$statistic, chisq_result$parameter, chisq_result$p.value))

p_66a <- ggplot(bar_df, aes(x = ctrl_label, y = pct, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", count, pct)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
  scale_fill_manual(values = COLORS$direction, name = "Direction") +
  scale_x_discrete(labels = SUBCMPT_LABELS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Significant 5mC DMRs by Subcompartment",
    subtitle = sprintf("Constitutive het (B.2) is stable; active regions (A.1) most affected | Chi² p = %.2e",
                        chisq_result$p.value),
    x = "CALDER2 Subcompartment (Control)", y = "% of Gene Bodies with Sig DMR (q<0.05)"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_66a,
                        file.path(OUTPUT_DIR, "66a_subcompartment_dmr_fraction"),
                        width = 10, height = 7)

# ---- Panel 66b: Direction split (hyper vs hypo) for sig DMRs -----------------

cat("Creating Figure 66b: direction split by subcompartment...\n")

dir_df <- mc_dmr_cmpt %>%
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

p_66b <- ggplot(dir_df, aes(x = ctrl_label, y = pct, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%\n(%d)", pct, n)),
            position = position_stack(vjust = 0.5), size = 3.5,
            color = "white", fontface = "bold") +
  scale_fill_manual(values = COLORS$direction, name = "Direction") +
  scale_x_discrete(labels = SUBCMPT_LABELS) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title = "Methylation Change Direction by Subcompartment",
    subtitle = "A.1 (active): dominated by hypermethylation | B.2 (constit het): dominated by hypomethylation",
    x = "CALDER2 Subcompartment (Control)", y = "% of Significant DMRs"
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_multiformat_ggplot(p_66b,
                        file.path(OUTPUT_DIR, "66b_direction_by_subcompartment"),
                        width = 10, height = 7)

# ---- Panel 66c: 5mC and 5hmC levels by subcompartment (violin) --------------

cat("Creating Figure 66c: methylation levels by subcompartment...\n")

# 5mC violin
mc_violin <- rbind(
  data.frame(gene = mc_dmr_cmpt$gene, subcmpt = mc_dmr_cmpt$ctrl_label,
             condition = "Control", value = mc_dmr_cmpt$mean_mod_group1,
             stringsAsFactors = FALSE),
  data.frame(gene = mc_dmr_cmpt$gene, subcmpt = mc_dmr_cmpt$ctrl_label,
             condition = "Mutant", value = mc_dmr_cmpt$mean_mod_group2,
             stringsAsFactors = FALSE)
) %>%
  dplyr::mutate(
    subcmpt = factor(subcmpt, levels = SUBCMPT_ORDER),
    condition = factor(condition, levels = c("Control", "Mutant"))
  )

# Wilcoxon per subcompartment
cat("  5mC Wilcoxon tests per subcompartment:\n")
mc_wilcox_results <- data.frame(subcmpt = character(), p = numeric(), stringsAsFactors = FALSE)
for (sc in SUBCMPT_ORDER) {
  ctrl_vals <- mc_dmr_cmpt$mean_mod_group1[mc_dmr_cmpt$ctrl_label == sc]
  mut_vals <- mc_dmr_cmpt$mean_mod_group2[mc_dmr_cmpt$ctrl_label == sc]
  wt <- wilcox.test(ctrl_vals, mut_vals, paired = TRUE)
  mc_wilcox_results <- rbind(mc_wilcox_results,
                             data.frame(subcmpt = sc, p = wt$p.value, stringsAsFactors = FALSE))
  cat(sprintf("    %s: ctrl median=%.4f, mut median=%.4f, p=%.2e (n=%d)\n",
              sc, median(ctrl_vals), median(mut_vals), wt$p.value, length(ctrl_vals)))
}

mc_wilcox_annot <- mc_wilcox_results %>%
  dplyr::mutate(
    subcmpt = factor(subcmpt, levels = SUBCMPT_ORDER),
    label = ifelse(p < 0.001, sprintf("p = %.1e", p), sprintf("p = %.3f", p)),
    y_pos = Inf
  )

p_66c_mc <- ggplot(mc_violin, aes(x = condition, y = value * 100, fill = condition)) +
  geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.8) +
  geom_text(data = mc_wilcox_annot,
            aes(x = 1.5, y = y_pos, label = label),
            inherit.aes = FALSE, vjust = 1.5, size = 3, fontface = "italic") +
  facet_wrap(~ subcmpt, labeller = labeller(subcmpt = SUBCMPT_LABELS)) +
  scale_fill_manual(values = COLORS$condition, guide = "none") +
  coord_cartesian(ylim = c(0, 100)) +
  labs(title = "5mC Levels by Subcompartment: Control vs Mutant",
       x = "", y = "5mC Level (%)") +
  theme_biomodal() +
  theme(strip.text = element_text(size = 10, face = "bold"))

# 5hmC violin
hmc_violin <- rbind(
  data.frame(gene = hmc_dmr_cmpt$gene, subcmpt = hmc_dmr_cmpt$ctrl_label,
             condition = "Control", value = hmc_dmr_cmpt$mean_mod_group1,
             stringsAsFactors = FALSE),
  data.frame(gene = hmc_dmr_cmpt$gene, subcmpt = hmc_dmr_cmpt$ctrl_label,
             condition = "Mutant", value = hmc_dmr_cmpt$mean_mod_group2,
             stringsAsFactors = FALSE)
) %>%
  dplyr::mutate(
    subcmpt = factor(subcmpt, levels = SUBCMPT_ORDER),
    condition = factor(condition, levels = c("Control", "Mutant"))
  )

cat("  5hmC Wilcoxon tests per subcompartment:\n")
hmc_wilcox_results <- data.frame(subcmpt = character(), p = numeric(), stringsAsFactors = FALSE)
for (sc in SUBCMPT_ORDER) {
  ctrl_vals <- hmc_dmr_cmpt$mean_mod_group1[hmc_dmr_cmpt$ctrl_label == sc]
  mut_vals <- hmc_dmr_cmpt$mean_mod_group2[hmc_dmr_cmpt$ctrl_label == sc]
  wt <- wilcox.test(ctrl_vals, mut_vals, paired = TRUE)
  hmc_wilcox_results <- rbind(hmc_wilcox_results,
                              data.frame(subcmpt = sc, p = wt$p.value, stringsAsFactors = FALSE))
  cat(sprintf("    %s: ctrl median=%.4f, mut median=%.4f, p=%.2e (n=%d)\n",
              sc, median(ctrl_vals), median(mut_vals), wt$p.value, length(ctrl_vals)))
}

hmc_wilcox_annot <- hmc_wilcox_results %>%
  dplyr::mutate(
    subcmpt = factor(subcmpt, levels = SUBCMPT_ORDER),
    label = ifelse(p < 0.001, sprintf("p = %.1e", p), sprintf("p = %.3f", p)),
    y_pos = Inf
  )

p_66c_hmc <- ggplot(hmc_violin, aes(x = condition, y = value * 100, fill = condition)) +
  geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.8) +
  geom_text(data = hmc_wilcox_annot,
            aes(x = 1.5, y = y_pos, label = label),
            inherit.aes = FALSE, vjust = 1.5, size = 3, fontface = "italic") +
  facet_wrap(~ subcmpt, labeller = labeller(subcmpt = SUBCMPT_LABELS)) +
  scale_fill_manual(values = COLORS$condition, guide = "none") +
  coord_cartesian(ylim = c(0, 25)) +
  labs(title = "5hmC Levels by Subcompartment: Control vs Mutant",
       x = "", y = "5hmC Level (%)") +
  theme_biomodal() +
  theme(strip.text = element_text(size = 10, face = "bold"))

p_66c <- p_66c_mc / p_66c_hmc +
  plot_annotation(
    title = "Methylation Levels by Subcompartment",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

save_multiformat_ggplot(p_66c,
                        file.path(OUTPUT_DIR, "66c_methylation_levels_by_subcompartment"),
                        width = 14, height = 12)

# ---- Panel 66d: Histone mark overlay ----------------------------------------

cat("\nCreating Figure 66d: histone mark overlay...\n")

# Load H3K27me3 and H3K27ac peaks
k27me3_gr <- load_chip_peaks(CHIP_PEAK_FILES$h3k27me3, "H3K27me3")
k27ac_gr  <- load_chip_peaks(CHIP_PEAK_FILES$h3k27ac, "H3K27ac")

# Compute overlaps with gene body DMR coordinates
dmr_gr <- dmr_to_granges(mc_dmr_cmpt)

mc_dmr_cmpt$k27me3_overlap <- countOverlaps(dmr_gr, k27me3_gr) > 0
mc_dmr_cmpt$k27ac_overlap  <- countOverlaps(dmr_gr, k27ac_gr) > 0

mc_dmr_cmpt$mark_category <- case_when(
  mc_dmr_cmpt$k27me3_overlap & mc_dmr_cmpt$k27ac_overlap ~ "K27me3 + K27ac\n(Bivalent)",
  mc_dmr_cmpt$k27me3_overlap ~ "K27me3 only\n(Fac. Het)",
  mc_dmr_cmpt$k27ac_overlap  ~ "K27ac only\n(Active)",
  TRUE                        ~ "Neither"
)

cat("  Gene body mark classification:\n")
mark_counts <- table(mc_dmr_cmpt$mark_category)
for (nm in names(mark_counts)) {
  cat(sprintf("    %s: %d\n", gsub("\n", " ", nm), mark_counts[nm]))
}

# Summary by histone mark category
mark_summary <- mc_dmr_cmpt %>%
  group_by(mark_category) %>%
  summarise(
    total = n(),
    sig = sum(significant),
    pct_sig = 100 * sum(significant) / n(),
    hyper = sum(significant & mod_difference > 0),
    hypo = sum(significant & mod_difference < 0),
    .groups = "drop"
  )

cat("\n  DMR significance by histone mark:\n")
for (i in seq_len(nrow(mark_summary))) {
  row <- mark_summary[i, ]
  cat(sprintf("    %s: %d total, %d sig (%.1f%%), %d hyper, %d hypo\n",
              gsub("\n", " ", row$mark_category),
              row$total, row$sig, row$pct_sig, row$hyper, row$hypo))
}

# Grouped bar: % sig by mark, split hyper/hypo
mark_bar_df <- mark_summary %>%
  tidyr::pivot_longer(cols = c(hyper, hypo), names_to = "direction", values_to = "count") %>%
  dplyr::mutate(
    pct = 100 * count / total,
    direction = factor(ifelse(direction == "hyper", "Hypermethylated", "Hypomethylated"),
                       levels = c("Hypermethylated", "Hypomethylated")),
    mark_category = factor(mark_category,
                           levels = c("K27ac only\n(Active)",
                                      "K27me3 + K27ac\n(Bivalent)",
                                      "K27me3 only\n(Fac. Het)",
                                      "Neither"))
  )

MARK_COLORS <- c("K27ac only\n(Active)" = "#E41A1C",
                 "K27me3 + K27ac\n(Bivalent)" = "#984EA3",
                 "K27me3 only\n(Fac. Het)" = "#4DAF4A",
                 "Neither" = "grey60")

p_66d <- ggplot(mark_bar_df, aes(x = mark_category, y = pct, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", count, pct)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
  scale_fill_manual(values = COLORS$direction, name = "Direction") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Significant 5mC DMRs by Histone Mark Context",
    subtitle = "H3K27me3 (fac. het) genes mirror B.1 pattern; H3K27ac (active) genes mirror A.1",
    x = "Histone Mark Category", y = "% of Gene Bodies with Sig DMR (q<0.05)"
  ) +
  theme_biomodal() +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 9))

save_multiformat_ggplot(p_66d,
                        file.path(OUTPUT_DIR, "66d_histone_mark_dmr_fraction"),
                        width = 10, height = 7)

# ---- Combined figure ---------------------------------------------------------

cat("Creating combined Figure 66...\n")

p_66_combined <- ((p_66a | p_66b) / p_66c_mc / p_66d) +
  plot_annotation(
    title = "Methylation Changes Are Compartment-Specific",
    subtitle = "Constitutive het (B.2) is stable; active regions (A.1) show TET-blockade hypermethylation; B compartments show hypomethylation",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  )

save_multiformat_ggplot(p_66_combined,
                        file.path(OUTPUT_DIR, "66_subcompartment_methylation"),
                        width = 16, height = 18)

# ---- Save tables -------------------------------------------------------------

cat("\nSaving tables...\n")

# Subcompartment summary
write.table(cmpt_summary,
            file.path(TABLES_DIR, "66_subcompartment_dmr_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 66_subcompartment_dmr_summary.tsv\n")

# Histone mark summary
write.table(mark_summary,
            file.path(TABLES_DIR, "66_histone_mark_dmr_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 66_histone_mark_dmr_summary.tsv\n")

# Per-gene compartment assignment
gene_assignment <- mc_dmr_cmpt %>%
  dplyr::select(gene, chr, start, end, ctrl_label, mut_label, label_changed,
                k27me3_overlap, k27ac_overlap, mark_category,
                significant, direction, mod_difference, dmr_qvalue)
write.table(gene_assignment,
            file.path(TABLES_DIR, "66_per_gene_compartment_assignment.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 66_per_gene_compartment_assignment.tsv\n")

# Wilcoxon results
wilcox_combined <- rbind(
  mc_wilcox_results %>% dplyr::mutate(modality = "5mC"),
  hmc_wilcox_results %>% dplyr::mutate(modality = "5hmC")
)
write.table(wilcox_combined,
            file.path(TABLES_DIR, "66_subcompartment_wilcoxon.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 66_subcompartment_wilcoxon.tsv\n")

cat("\n")
cat("================================================================================\n")
cat("SECTION 66 COMPLETE\n")
cat("================================================================================\n")
