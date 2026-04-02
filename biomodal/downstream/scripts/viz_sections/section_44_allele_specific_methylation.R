# biomodal/downstream/scripts/viz_sections/section_44_allele_specific_methylation.R
# Section 44: Allele-Specific Methylation (ASM) Exploratory Analysis
# First downstream analysis of ASM data from DUET evoC pipeline
# Explores mC and hmC allelic imbalance at heterozygous SNV sites
# Standalone script - sources shared config for all dependencies and data

# Run from downstream/ directory
source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 44: ALLELE-SPECIFIC METHYLATION EXPLORATORY ANALYSIS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 44: ALLELE-SPECIFIC METHYLATION EXPLORATORY ANALYSIS\n")
cat("================================================================================\n\n")

# -----------------------------------------------------------------------------
# SETUP: ASM file paths and loader
# -----------------------------------------------------------------------------

ASM_DIR <- file.path(BASE_DIR, "../upstream/duet-1.5.0_evoC_Bap1_run_6bp/allele_specific_methylation")
stopifnot("ASM directory not found" = dir.exists(ASM_DIR))

# Build paths from SAMPLES (individual files have headers, unlike all_samples files)
ASM_FILES <- list(
  mC = setNames(
    file.path(ASM_DIR, paste0(SAMPLES$sample_id, ".mC.asm.csv")),
    SAMPLES$sample_id
  ),
  hmC = setNames(
    file.path(ASM_DIR, paste0(SAMPLES$sample_id, ".hmC.asm.csv")),
    SAMPLES$sample_id
  )
)

# Validate all 16 files exist
missing_mc  <- names(ASM_FILES$mC)[!file.exists(ASM_FILES$mC)]
missing_hmc <- names(ASM_FILES$hmC)[!file.exists(ASM_FILES$hmC)]
stopifnot("Missing mC ASM files" = length(missing_mc) == 0)
stopifnot("Missing hmC ASM files" = length(missing_hmc) == 0)

cat(sprintf("ASM directory: %s\n", ASM_DIR))
cat(sprintf("Files found: %d mC + %d hmC\n", length(ASM_FILES$mC), length(ASM_FILES$hmC)))

# Valid chromosomes (exclude scaffold contigs)
VALID_CHR <- paste0("chr", c(1:19, "X", "Y"))

# Section output directory
SECTION_DIR <- "44_allele_specific_methylation"
dir.create(file.path(OUTPUT_DIR, SECTION_DIR), recursive = TRUE, showWarnings = FALSE)

# ASM loader: reads individual sample CSV (with header), standardizes columns
load_asm_csv <- function(filepath, mod_type = "mC") {
  stopifnot(file.exists(filepath))
  df <- read.csv(filepath, stringsAsFactors = FALSE)

  # Standardize modification columns to generic names
  if (mod_type == "hmC") {
    colnames(df)[colnames(df) == "hmC_al1"] <- "mod_al1"
    colnames(df)[colnames(df) == "hmC_al2"] <- "mod_al2"
  } else {
    colnames(df)[colnames(df) == "mC_al1"] <- "mod_al1"
    colnames(df)[colnames(df) == "mC_al2"] <- "mod_al2"
  }

  # Add chr prefix to contig
  df$chr <- paste0("chr", df$contig)

  # Parse numeric fields (empty strings from non-PASS rows become NA)
  df$p_value <- suppressWarnings(as.numeric(df$p_value))
  df$corr_p_value <- suppressWarnings(as.numeric(df$corr_p_value))
  df$test_statistic <- suppressWarnings(as.numeric(df$test_statistic))
  df$methylation_diff <- suppressWarnings(as.numeric(df$methylation_diff))

  # Cap infinite test statistics
  df$test_statistic[is.infinite(df$test_statistic)] <- NA

  # Filter to standard chromosomes
  df <- df[df$chr %in% VALID_CHR, , drop = FALSE]

  df$mod_type <- mod_type
  return(df)
}

# -----------------------------------------------------------------------------
# DATA LOADING
# -----------------------------------------------------------------------------

cat("\nLoading ASM data from individual sample files...\n")

asm_mc_list <- lapply(names(ASM_FILES$mC), function(sid) {
  cat(sprintf("  Loading mC: %s\n", sid))
  load_asm_csv(ASM_FILES$mC[[sid]], mod_type = "mC")
})
asm_mc_all <- dplyr::bind_rows(asm_mc_list)
rm(asm_mc_list)

asm_hmc_list <- lapply(names(ASM_FILES$hmC), function(sid) {
  cat(sprintf("  Loading hmC: %s\n", sid))
  load_asm_csv(ASM_FILES$hmC[[sid]], mod_type = "hmC")
})
asm_hmc_all <- dplyr::bind_rows(asm_hmc_list)
rm(asm_hmc_list)

# Join sample metadata
asm_mc_all <- asm_mc_all %>%
  dplyr::left_join(SAMPLES, by = "sample_id")
asm_hmc_all <- asm_hmc_all %>%
  dplyr::left_join(SAMPLES, by = "sample_id")

# Subsets
asm_mc_pass <- asm_mc_all %>% dplyr::filter(filter == "PASS")
asm_hmc_pass <- asm_hmc_all %>% dplyr::filter(filter == "PASS")
asm_mc_sig <- asm_mc_pass %>% dplyr::filter(!is.na(corr_p_value) & corr_p_value < Q_THRESHOLD)
asm_hmc_sig <- asm_hmc_pass %>% dplyr::filter(!is.na(corr_p_value) & corr_p_value < Q_THRESHOLD)

cat(sprintf("\nmC:  %d total | %d PASS (%.1f%%) | %d significant (%.1f%% of PASS)\n",
            nrow(asm_mc_all), nrow(asm_mc_pass),
            100 * nrow(asm_mc_pass) / nrow(asm_mc_all),
            nrow(asm_mc_sig),
            100 * nrow(asm_mc_sig) / max(nrow(asm_mc_pass), 1)))
cat(sprintf("hmC: %d total | %d PASS (%.1f%%) | %d significant (%.1f%% of PASS)\n",
            nrow(asm_hmc_all), nrow(asm_hmc_pass),
            100 * nrow(asm_hmc_pass) / nrow(asm_hmc_all),
            nrow(asm_hmc_sig),
            100 * nrow(asm_hmc_sig) / max(nrow(asm_hmc_pass), 1)))

# =============================================================================
# FIGURE 44a: FILTER CATEGORY DISTRIBUTION PER SAMPLE
# =============================================================================

cat("\n--- Figure 44a: Filter Category Distribution Per Sample ---\n")

# Combine both mod types for filter overview
filter_mc <- asm_mc_all %>%
  dplyr::count(short_name, condition, filter, name = "count") %>%
  dplyr::mutate(mod_type = "5mC")
filter_hmc <- asm_hmc_all %>%
  dplyr::count(short_name, condition, filter, name = "count") %>%
  dplyr::mutate(mod_type = "5hmC")
filter_combined <- dplyr::bind_rows(filter_mc, filter_hmc)

# Order filter categories by total count
filter_order <- filter_combined %>%
  dplyr::group_by(filter) %>%
  dplyr::summarize(total = sum(count), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(total)) %>%
  dplyr::pull(filter)
filter_combined$filter <- factor(filter_combined$filter, levels = rev(filter_order))

# Order samples: controls then mutants
sample_order <- c("ctrl-F", "ctrl-M", "ctrl-F-B2", "ctrl-M-B2",
                  "mut-F", "mut-M", "mut-F-B2", "mut-M-B2")
filter_combined$short_name <- factor(filter_combined$short_name, levels = sample_order)

filter_palette <- c(
  "PASS" = "#4DAF4A",
  "LOW_METH_DIFF" = "#FF7F00",
  "NO_CPG_READS" = "#984EA3",
  "LOW_CPG_COVERAGE_AL1_AL2" = "#E41A1C",
  "LOW_CPG_COVERAGE_AL2" = "#A65628",
  "LOW_CPG_COVERAGE_AL1" = "#F781BF",
  "NO_CPG_COVERAGE_AL2" = "#999999",
  "NO_CPG_COVERAGE_AL1" = "#666666",
  "NO_READS" = "#377EB8"
)

p_filter <- ggplot(filter_combined, aes(x = short_name, y = count, fill = filter)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~mod_type, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = filter_palette, name = "Filter Status") +
  labs(
    title = "ASM Filter Category Distribution Per Sample",
    subtitle = "Heterozygous SNV sites | mC and hmC modification types",
    x = "Sample", y = "Number of Sites"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

save_multiformat_ggplot(p_filter, file.path(OUTPUT_DIR, SECTION_DIR, "44a_filter_distribution_per_sample"),
                        width = 14, height = 10)
cat("  Saved 44a_filter_distribution_per_sample\n")

# Release full data frames — subsequent figures use PASS/sig subsets
rm(asm_mc_all, asm_hmc_all, filter_mc, filter_hmc, filter_combined)
gc()

# =============================================================================
# FIGURE 44b: PASS SITE COUNTS PER SAMPLE BY CONDITION
# =============================================================================

cat("\n--- Figure 44b: PASS Site Counts Per Sample ---\n")

pass_counts_mc <- asm_mc_pass %>%
  dplyr::count(short_name, condition, name = "count") %>%
  dplyr::mutate(mod_type = "5mC")
pass_counts_hmc <- asm_hmc_pass %>%
  dplyr::count(short_name, condition, name = "count") %>%
  dplyr::mutate(mod_type = "5hmC")
pass_counts <- dplyr::bind_rows(pass_counts_mc, pass_counts_hmc)
pass_counts$short_name <- factor(pass_counts$short_name, levels = sample_order)

p_pass <- ggplot(pass_counts, aes(x = short_name, y = count, fill = condition)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = scales::comma(count)), vjust = -0.3, size = 3.2) +
  facet_wrap(~mod_type, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = COLORS$condition, name = "Condition") +
  labs(
    title = "PASS ASM Sites Per Sample",
    subtitle = "Sites with sufficient coverage for allele-specific evaluation",
    x = "Sample", y = "Number of PASS Sites"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_multiformat_ggplot(p_pass, file.path(OUTPUT_DIR, SECTION_DIR, "44b_pass_site_counts_per_sample"),
                        width = 10, height = 8)
cat("  Saved 44b_pass_site_counts_per_sample\n")

# =============================================================================
# FIGURE 44c: CHROMOSOME DISTRIBUTION OF SIGNIFICANT mC ASM SITES
# =============================================================================

cat("\n--- Figure 44c: Chromosome Distribution (mC ASM) ---\n")

# Count unique significant loci per chromosome (deduplicate across samples)
mc_sig_loci <- asm_mc_sig %>%
  dplyr::distinct(chr, pos_first_variant) %>%
  dplyr::count(chr, name = "sig_count")

mc_pass_loci <- asm_mc_pass %>%
  dplyr::distinct(chr, pos_first_variant) %>%
  dplyr::count(chr, name = "total_count")

chr_dist <- mc_pass_loci %>%
  dplyr::left_join(mc_sig_loci, by = "chr") %>%
  dplyr::mutate(
    sig_count = replace(sig_count, is.na(sig_count), 0),
    expected = total_count / sum(total_count) * sum(sig_count),
    obs_exp_ratio = sig_count / expected,
    chr_label = gsub("chr", "", chr),
    chr_num = suppressWarnings(as.numeric(chr_label)),
    chr_order = ifelse(is.na(chr_num), 100 + match(chr_label, c("X", "Y")), chr_num)
  ) %>%
  dplyr::arrange(chr_order) %>%
  dplyr::mutate(chr = factor(chr, levels = chr))

# Fisher's exact test per chromosome
n_sig_total <- sum(chr_dist$sig_count)
chr_dist$fisher_p <- NA_real_
chr_dist$fisher_or <- NA_real_
for (i in seq_len(nrow(chr_dist))) {
  sig_on   <- chr_dist$sig_count[i]
  sig_off  <- n_sig_total - sig_on
  total_on <- chr_dist$total_count[i]
  total_off <- sum(chr_dist$total_count) - total_on
  nonsig_on  <- total_on - sig_on
  nonsig_off <- total_off - sig_off

  mat <- matrix(c(sig_on, nonsig_on, sig_off, nonsig_off), nrow = 2)
  ft <- fisher.test(mat)
  chr_dist$fisher_p[i] <- ft$p.value
  chr_dist$fisher_or[i] <- ft$estimate
}

chr_dist$enrichment <- dplyr::case_when(
  chr_dist$fisher_p < 0.05 & chr_dist$obs_exp_ratio > 1 ~ "Enriched",
  chr_dist$fisher_p < 0.05 & chr_dist$obs_exp_ratio < 1 ~ "Depleted",
  TRUE ~ "Not significant"
)
chr_dist$fisher_label <- ifelse(
  chr_dist$fisher_p < 0.001,
  sprintf("p=%.1e", chr_dist$fisher_p),
  ifelse(chr_dist$fisher_p < 0.05,
         sprintf("p=%.3f", chr_dist$fisher_p), "")
)

# Log enrichment results
chr_enriched <- chr_dist %>% dplyr::filter(enrichment == "Enriched")
chr_depleted <- chr_dist %>% dplyr::filter(enrichment == "Depleted")
cat(sprintf("  Enriched chromosomes: %s\n",
            ifelse(nrow(chr_enriched) > 0,
                   paste(chr_enriched$chr, collapse = ", "), "none")))
cat(sprintf("  Depleted chromosomes: %s\n",
            ifelse(nrow(chr_depleted) > 0,
                   paste(chr_depleted$chr, collapse = ", "), "none")))

p_chr <- ggplot(chr_dist, aes(x = chr, y = sig_count)) +
  geom_bar(aes(fill = enrichment), stat = "identity", width = 0.7) +
  geom_point(aes(y = expected), color = "black", shape = 4, size = 3) +
  geom_text(aes(label = fisher_label), vjust = -0.5, size = 3, color = "#E41A1C",
            fontface = "bold") +
  scale_fill_manual(
    values = c("Enriched" = "#D7191C", "Depleted" = "#2C7BB6", "Not significant" = "grey70"),
    name = "Enrichment"
  ) +
  labs(
    title = "Chromosome Distribution of Significant mC ASM Sites",
    subtitle = sprintf("%d unique significant loci (corr_p < %.2f) | X = expected count",
                       sum(chr_dist$sig_count), Q_THRESHOLD),
    x = "Chromosome", y = "Significant ASM Loci"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_multiformat_ggplot(p_chr, file.path(OUTPUT_DIR, SECTION_DIR, "44c_asm_mc_chromosome_distribution"),
                        width = 12, height = 7)
cat("  Saved 44c_asm_mc_chromosome_distribution\n")

# =============================================================================
# FIGURE 44d: METHYLATION DIFFERENCE DISTRIBUTIONS
# =============================================================================

cat("\n--- Figure 44d: Methylation Difference Distributions ---\n")

# Panel A: Density of methylation_diff for significant mC ASM, by condition
p_density <- ggplot(asm_mc_sig, aes(x = methylation_diff, fill = condition)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0.3, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = COLORS$condition, name = "Condition") +
  labs(
    title = "mC Allelic Methylation Difference Distribution",
    subtitle = sprintf("Significant ASM sites (corr_p < %.2f) | Dashed = 0.3 PASS threshold",
                       Q_THRESHOLD),
    x = "Methylation Difference (|allele 1 - allele 2|)",
    y = "Density"
  ) +
  theme_biomodal() +
  theme(legend.position = "bottom")

# Panel B: Boxplot of abs(methylation_diff) per sample, colored by condition
p_box <- ggplot(asm_mc_sig, aes(x = factor(short_name, levels = sample_order),
                                y = methylation_diff, fill = condition)) +
  geom_boxplot(outlier.alpha = 0.3, outlier.size = 0.5) +
  scale_fill_manual(values = COLORS$condition, name = "Condition") +
  labs(
    title = "Per-Sample Allelic Methylation Bias Magnitude",
    subtitle = "Significant mC ASM sites",
    x = "Sample", y = "Methylation Difference"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

p_methdiff <- p_density / p_box +
  plot_annotation(
    title = "Allelic Methylation Difference: mC ASM",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )

save_multiformat_ggplot(p_methdiff, file.path(OUTPUT_DIR, SECTION_DIR, "44d_methylation_diff_distributions"),
                        width = 12, height = 10)
cat("  Saved 44d_methylation_diff_distributions\n")

# =============================================================================
# FIGURE 44e: MUTANT vs CONTROL ASM SITE COUNTS
# =============================================================================

cat("\n--- Figure 44e: Mutant vs Control ASM Comparison ---\n")

sig_per_sample <- asm_mc_sig %>%
  dplyr::count(sample_id, condition, short_name, name = "sig_count")

# Wilcoxon test
ctrl_counts <- sig_per_sample %>% dplyr::filter(condition == "Control") %>% dplyr::pull(sig_count)
mut_counts <- sig_per_sample %>% dplyr::filter(condition == "Mutant") %>% dplyr::pull(sig_count)
wt <- wilcox.test(ctrl_counts, mut_counts)

condition_means <- sig_per_sample %>%
  dplyr::group_by(condition) %>%
  dplyr::summarize(
    mean_count = mean(sig_count),
    sd_count = sd(sig_count),
    .groups = "drop"
  )

p_compare <- ggplot(sig_per_sample, aes(x = condition, y = sig_count, fill = condition)) +
  geom_bar(data = condition_means, aes(x = condition, y = mean_count),
           stat = "identity", width = 0.6, alpha = 0.7) +
  geom_errorbar(data = condition_means,
                aes(x = condition, y = mean_count,
                    ymin = mean_count - sd_count, ymax = mean_count + sd_count),
                width = 0.2, inherit.aes = FALSE) +
  geom_jitter(width = 0.1, size = 3, shape = 21, color = "black") +
  scale_fill_manual(values = COLORS$condition, name = "Condition") +
  labs(
    title = "Significant mC ASM Sites: Mutant vs Control",
    subtitle = sprintf("Wilcoxon p = %.4f | Mean mut = %.0f, ctrl = %.0f (%.1fx)",
                       wt$p.value,
                       condition_means$mean_count[condition_means$condition == "Mutant"],
                       condition_means$mean_count[condition_means$condition == "Control"],
                       condition_means$mean_count[condition_means$condition == "Mutant"] /
                         condition_means$mean_count[condition_means$condition == "Control"]),
    x = "", y = "Significant ASM Sites Per Sample"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

save_multiformat_ggplot(p_compare, file.path(OUTPUT_DIR, SECTION_DIR, "44e_mutant_vs_control_asm_comparison"),
                        width = 8, height = 7)
cat(sprintf("  Saved 44e_mutant_vs_control_asm_comparison (Wilcoxon p=%.4f)\n", wt$p.value))

# =============================================================================
# FIGURE 44f: P-VALUE DISTRIBUTIONS
# =============================================================================

cat("\n--- Figure 44f: P-value Distributions ---\n")

# Left panel: raw p-value histogram (uniform under null)
p_pval_hist <- ggplot(asm_mc_pass %>% dplyr::filter(!is.na(p_value)),
                      aes(x = p_value, fill = condition)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  facet_wrap(~condition, ncol = 2) +
  scale_fill_manual(values = COLORS$condition, name = "Condition") +
  labs(
    title = "Raw P-value Distribution (mC ASM)",
    subtitle = "PASS sites | Uniform distribution expected under null",
    x = "Raw P-value", y = "Count"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

# Right panel: -log10(corrected p-value) distribution
asm_mc_pass_plot <- asm_mc_pass %>%
  dplyr::filter(!is.na(corr_p_value)) %>%
  dplyr::mutate(neg_log10_corrp = pmin(-log10(corr_p_value), 50))

p_corrp <- ggplot(asm_mc_pass_plot, aes(x = neg_log10_corrp, fill = condition)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 10, linetype = "dotted", color = "#E41A1C") +
  scale_fill_manual(values = COLORS$condition, name = "Condition") +
  annotate("text", x = -log10(0.05) + 1, y = Inf, vjust = 2,
           label = "q = 0.05", size = 3, color = "grey40") +
  annotate("text", x = 11, y = Inf, vjust = 2,
           label = "Strong ASM\n(imprinting)", size = 3, color = "#E41A1C") +
  labs(
    title = expression(-log[10](corrected~p-value)~Distribution),
    subtitle = "PASS mC ASM sites | Capped at 50",
    x = expression(-log[10](corrected~p-value)), y = "Count"
  ) +
  theme_biomodal() +
  theme(legend.position = "bottom")

p_pval_combined <- p_pval_hist / p_corrp +
  plot_annotation(
    title = "Statistical Significance of mC Allele-Specific Methylation",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )

save_multiformat_ggplot(p_pval_combined, file.path(OUTPUT_DIR, SECTION_DIR, "44f_pvalue_distributions"),
                        width = 12, height = 10)
cat("  Saved 44f_pvalue_distributions\n")
rm(asm_mc_pass_plot)

# =============================================================================
# FIGURE 44g: ASM VOLCANO PLOT
# =============================================================================

cat("\n--- Figure 44g: ASM Volcano Plot ---\n")

volcano_data <- asm_mc_pass %>%
  dplyr::filter(!is.na(corr_p_value) & !is.na(methylation_diff)) %>%
  dplyr::mutate(
    neg_log10_corrp = pmin(-log10(corr_p_value), 300),
    sig_status = ifelse(corr_p_value < Q_THRESHOLD, "Significant", "Not Significant")
  )

# Count per quadrant for annotation
n_sig_ctrl <- sum(volcano_data$sig_status == "Significant" & volcano_data$condition == "Control")
n_sig_mut <- sum(volcano_data$sig_status == "Significant" & volcano_data$condition == "Mutant")

p_volcano <- ggplot(volcano_data, aes(x = methylation_diff, y = neg_log10_corrp)) +
  geom_point(aes(color = sig_status), alpha = 0.4, size = 0.8) +
  geom_hline(yintercept = -log10(Q_THRESHOLD), linetype = "dashed", color = "grey40") +
  facet_wrap(~condition, ncol = 2) +
  scale_color_manual(
    values = c("Significant" = "#E41A1C", "Not Significant" = "grey70"),
    name = ""
  ) +
  labs(
    title = "mC ASM Volcano: Allelic Methylation Difference vs Significance",
    subtitle = sprintf("PASS sites | Control: %s sig, Mutant: %s sig",
                       scales::comma(n_sig_ctrl), scales::comma(n_sig_mut)),
    x = "Methylation Difference Between Alleles",
    y = expression(-log[10](corrected~p-value))
  ) +
  coord_cartesian(ylim = c(0, min(50, max(volcano_data$neg_log10_corrp, na.rm = TRUE)))) +
  theme_biomodal() +
  theme(legend.position = "bottom")

save_multiformat_ggplot(p_volcano, file.path(OUTPUT_DIR, SECTION_DIR, "44g_asm_volcano_plot"),
                        width = 12, height = 7)
cat("  Saved 44g_asm_volcano_plot\n")
rm(volcano_data)

# =============================================================================
# FIGURE 44h: SHARED ASM SITES BETWEEN CONDITIONS (VENN)
# =============================================================================

cat("\n--- Figure 44h: Shared ASM Sites Venn ---\n")

# Unique locus key: chr:pos_first_variant
# A locus is "in condition" if significant in ANY sample of that condition
ctrl_loci <- asm_mc_sig %>%
  dplyr::filter(condition == "Control") %>%
  dplyr::mutate(locus = paste(chr, pos_first_variant, sep = ":")) %>%
  dplyr::pull(locus) %>%
  unique()

mut_loci <- asm_mc_sig %>%
  dplyr::filter(condition == "Mutant") %>%
  dplyr::mutate(locus = paste(chr, pos_first_variant, sep = ":")) %>%
  dplyr::pull(locus) %>%
  unique()

venn_sets <- list(
  "Control" = ctrl_loci,
  "Mutant" = mut_loci
)

n_ctrl_only <- length(setdiff(ctrl_loci, mut_loci))
n_mut_only <- length(setdiff(mut_loci, ctrl_loci))
n_shared <- length(intersect(ctrl_loci, mut_loci))

cat(sprintf("  Control-only: %d | Mutant-only: %d | Shared: %d\n",
            n_ctrl_only, n_mut_only, n_shared))

p_venn <- ggVennDiagram(venn_sets, label_alpha = 0, set_size = 5) +
  scale_fill_gradient(low = "white", high = "#E41A1C") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  labs(
    title = "Significant mC ASM Loci: Control vs Mutant",
    subtitle = sprintf("%s shared | %s control-only | %s mutant-only",
                       scales::comma(n_shared),
                       scales::comma(n_ctrl_only),
                       scales::comma(n_mut_only))
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

save_multiformat_ggplot(p_venn, file.path(OUTPUT_DIR, SECTION_DIR, "44h_shared_asm_sites_venn"),
                        width = 8, height = 7)
cat("  Saved 44h_shared_asm_sites_venn\n")

# =============================================================================
# FIGURE 44i: ASM-DMR OVERLAP
# =============================================================================

cat("\n--- Figure 44i: ASM-DMR Overlap ---\n")

stopifnot("mc_dmr not loaded" = exists("mc_dmr") && !is.null(mc_dmr))
stopifnot("hmc_dmr not loaded" = exists("hmc_dmr") && !is.null(hmc_dmr))

# Create GRanges for ASM sites (point ranges from pos_first_variant)
asm_sig_unique <- asm_mc_sig %>%
  dplyr::distinct(chr, pos_first_variant, .keep_all = TRUE)

asm_gr <- GenomicRanges::GRanges(
  seqnames = asm_sig_unique$chr,
  ranges = IRanges::IRanges(start = asm_sig_unique$pos_first_variant,
                            end = asm_sig_unique$pos_first_variant)
)

# DMR GRanges (significant only)
mc_sig_dmr <- mc_dmr %>% dplyr::filter(significant)
hmc_sig_dmr <- hmc_dmr %>% dplyr::filter(significant)

mc_hyper_gr <- dmr_to_granges(mc_sig_dmr %>% dplyr::filter(mod_difference > 0))
mc_hypo_gr <- dmr_to_granges(mc_sig_dmr %>% dplyr::filter(mod_difference <= 0))
hmc_hyper_gr <- dmr_to_granges(hmc_sig_dmr %>% dplyr::filter(mod_difference > 0))
hmc_hypo_gr <- dmr_to_granges(hmc_sig_dmr %>% dplyr::filter(mod_difference <= 0))

# Compute overlaps
n_asm <- length(asm_gr)
n_mc_hyper_ol <- length(unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(asm_gr, mc_hyper_gr))))
n_mc_hypo_ol <- length(unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(asm_gr, mc_hypo_gr))))
n_hmc_hyper_ol <- length(unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(asm_gr, hmc_hyper_gr))))
n_hmc_hypo_ol <- length(unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(asm_gr, hmc_hypo_gr))))

# Any DMR overlap
any_dmr_gr <- c(mc_hyper_gr, mc_hypo_gr, hmc_hyper_gr, hmc_hypo_gr)
n_any_dmr_ol <- length(unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(asm_gr, any_dmr_gr))))
n_no_dmr <- n_asm - n_any_dmr_ol

cat(sprintf("  ASM sites in mC hyper-DMRs: %d\n", n_mc_hyper_ol))
cat(sprintf("  ASM sites in mC hypo-DMRs:  %d\n", n_mc_hypo_ol))
cat(sprintf("  ASM sites in hmC hyper-DMRs: %d\n", n_hmc_hyper_ol))
cat(sprintf("  ASM sites in hmC hypo-DMRs: %d\n", n_hmc_hypo_ol))
cat(sprintf("  No DMR overlap: %d\n", n_no_dmr))

overlap_df <- data.frame(
  category = c("mC Hyper-DMR", "mC Hypo-DMR", "hmC Hyper-DMR", "hmC Hypo-DMR", "No DMR Overlap"),
  count = c(n_mc_hyper_ol, n_mc_hypo_ol, n_hmc_hyper_ol, n_hmc_hypo_ol, n_no_dmr),
  stringsAsFactors = FALSE
)
overlap_df$category <- factor(overlap_df$category, levels = overlap_df$category)
overlap_df$pct <- 100 * overlap_df$count / n_asm

p_overlap <- ggplot(overlap_df, aes(x = category, y = count, fill = category)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", count, pct)), vjust = -0.3, size = 3.5) +
  scale_fill_manual(
    values = c("mC Hyper-DMR" = "#D7191C", "mC Hypo-DMR" = "#2C7BB6",
               "hmC Hyper-DMR" = "#FDB863", "hmC Hypo-DMR" = "#5E4FA2",
               "No DMR Overlap" = "grey70"),
    name = ""
  ) +
  labs(
    title = "Significant ASM Sites Overlapping Gene Body DMRs",
    subtitle = sprintf("%s unique ASM loci | DMR q < %.2f",
                       scales::comma(n_asm), Q_THRESHOLD),
    x = "", y = "ASM Loci Count"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")

# Panel B: For genes with both DMR and ASM, scatter DMR effect vs ASM bias
# Annotate ASM sites with their overlapping DMR gene
mc_all_dmr_gr <- c(mc_hyper_gr, mc_hypo_gr)
ol_mc <- GenomicRanges::findOverlaps(asm_gr, mc_all_dmr_gr)

if (length(ol_mc) > 0) {
  asm_dmr_pairs <- data.frame(
    asm_idx = S4Vectors::queryHits(ol_mc),
    dmr_idx = S4Vectors::subjectHits(ol_mc),
    stringsAsFactors = FALSE
  )
  # Get DMR info
  mc_all_sig <- dplyr::bind_rows(
    mc_sig_dmr %>% dplyr::filter(mod_difference > 0),
    mc_sig_dmr %>% dplyr::filter(mod_difference <= 0)
  )
  asm_dmr_pairs$dmr_gene <- mc_all_sig$gene[asm_dmr_pairs$dmr_idx]
  asm_dmr_pairs$dmr_diff <- mc_all_sig$mod_difference[asm_dmr_pairs$dmr_idx]
  asm_dmr_pairs$asm_meth_diff <- asm_sig_unique$methylation_diff[asm_dmr_pairs$asm_idx]

  # Aggregate per gene: mean ASM methylation_diff
  gene_summary <- asm_dmr_pairs %>%
    dplyr::group_by(dmr_gene, dmr_diff) %>%
    dplyr::summarize(
      mean_asm_diff = mean(asm_meth_diff, na.rm = TRUE),
      n_asm_sites = dplyr::n(),
      .groups = "drop"
    )

  p_scatter <- ggplot(gene_summary, aes(x = dmr_diff * 100, y = mean_asm_diff)) +
    geom_point(aes(size = n_asm_sites), alpha = 0.5, color = COLORS$methylation[["5mC"]]) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey40") +
    scale_size_continuous(range = c(1, 6), name = "ASM Sites") +
    labs(
      title = "DMR Effect Size vs ASM Allelic Bias",
      subtitle = sprintf("%d genes with both DMR and ASM", nrow(gene_summary)),
      x = "mC DMR Difference (Mutant - Control, %)",
      y = "Mean ASM Methylation Difference"
    ) +
    theme_biomodal()

  p_dmr_overlap <- p_overlap + p_scatter +
    plot_annotation(
      title = "ASM and Gene Body DMR Integration",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    )
} else {
  cat("  No ASM-DMR overlaps found for scatter plot\n")
  p_dmr_overlap <- p_overlap +
    plot_annotation(
      title = "ASM and Gene Body DMR Integration",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    )
}

save_multiformat_ggplot(p_dmr_overlap, file.path(OUTPUT_DIR, SECTION_DIR, "44i_asm_dmr_overlap"),
                        width = 14, height = 7)
cat("  Saved 44i_asm_dmr_overlap\n")

# =============================================================================
# FIGURE 44j: mC vs hmC ASM COMPARISON
# =============================================================================

cat("\n--- Figure 44j: mC vs hmC ASM Comparison ---\n")

# Panel A: Side-by-side significant counts per sample
mc_sig_counts <- asm_mc_sig %>%
  dplyr::count(short_name, condition, name = "sig_count") %>%
  dplyr::mutate(mod_type = "5mC")
hmc_sig_counts <- asm_hmc_sig %>%
  dplyr::count(short_name, condition, name = "sig_count") %>%
  dplyr::mutate(mod_type = "5hmC")

# Ensure all samples present for hmC (some may have 0 significant sites)
all_samples_hmc <- SAMPLES %>%
  dplyr::select(short_name, condition) %>%
  dplyr::mutate(mod_type = "5hmC")
hmc_sig_counts <- all_samples_hmc %>%
  dplyr::left_join(hmc_sig_counts %>% dplyr::select(short_name, sig_count),
                   by = "short_name") %>%
  dplyr::mutate(sig_count = replace(sig_count, is.na(sig_count), 0))

sig_counts_combined <- dplyr::bind_rows(mc_sig_counts, hmc_sig_counts)
sig_counts_combined$short_name <- factor(sig_counts_combined$short_name, levels = sample_order)

p_mc_hmc_bar <- ggplot(sig_counts_combined, aes(x = short_name, y = sig_count, fill = condition)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = scales::comma(sig_count)), vjust = -0.3, size = 3) +
  facet_wrap(~mod_type, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = COLORS$condition, name = "Condition") +
  labs(
    title = "Significant ASM Sites: mC vs hmC",
    subtitle = sprintf("mC: %s total sig | hmC: %s total sig (exploratory, limited power)",
                       scales::comma(nrow(asm_mc_sig)),
                       scales::comma(nrow(asm_hmc_sig))),
    x = "Sample", y = "Significant ASM Sites"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Panel B: Venn of shared significant loci between mC and hmC
mc_sig_loci_set <- asm_mc_sig %>%
  dplyr::mutate(locus = paste(chr, pos_first_variant, sep = ":")) %>%
  dplyr::pull(locus) %>%
  unique()
hmc_sig_loci_set <- asm_hmc_sig %>%
  dplyr::mutate(locus = paste(chr, pos_first_variant, sep = ":")) %>%
  dplyr::pull(locus) %>%
  unique()

n_mc_only_loci <- length(setdiff(mc_sig_loci_set, hmc_sig_loci_set))
n_hmc_only_loci <- length(setdiff(hmc_sig_loci_set, mc_sig_loci_set))
n_both_loci <- length(intersect(mc_sig_loci_set, hmc_sig_loci_set))

if (length(hmc_sig_loci_set) > 0) {
  venn_mod <- list(
    "mC ASM" = mc_sig_loci_set,
    "hmC ASM" = hmc_sig_loci_set
  )
  p_venn_mod <- ggVennDiagram(venn_mod, label_alpha = 0, set_size = 5) +
    scale_fill_gradient(low = "white", high = "#377EB8") +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    labs(
      title = "mC vs hmC Significant ASM Loci Overlap",
      subtitle = sprintf("%d shared | %s mC-only | %d hmC-only",
                         n_both_loci,
                         scales::comma(n_mc_only_loci),
                         n_hmc_only_loci)
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))

  p_mc_hmc <- p_mc_hmc_bar + p_venn_mod +
    plot_layout(widths = c(2, 1)) +
    plot_annotation(
      title = "mC vs hmC Allele-Specific Methylation Comparison",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    )
} else {
  cat("  No significant hmC ASM loci — Venn not generated\n")
  p_mc_hmc <- p_mc_hmc_bar +
    plot_annotation(
      title = "mC vs hmC Allele-Specific Methylation Comparison",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    )
}

save_multiformat_ggplot(p_mc_hmc, file.path(OUTPUT_DIR, SECTION_DIR, "44j_mc_vs_hmc_asm_comparison"),
                        width = 14, height = 8)
cat(sprintf("  Saved 44j_mc_vs_hmc_asm_comparison (mC-only: %s, hmC-only: %d, both: %d)\n",
            scales::comma(n_mc_only_loci), n_hmc_only_loci, n_both_loci))

# =============================================================================
# TABLE EXPORTS
# =============================================================================

cat("\n--- Table Exports ---\n")

# Table 1: Per-sample significant mC ASM counts
table_per_sample <- sig_per_sample %>%
  dplyr::left_join(SAMPLES %>% dplyr::select(sample_id, sex, batch), by = "sample_id") %>%
  dplyr::arrange(condition, sample_id)

write.table(table_per_sample,
            file.path(TABLES_DIR, "asm_mc_significant_per_sample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved asm_mc_significant_per_sample.tsv (%d samples)\n", nrow(table_per_sample)))

# Table 2: Unique significant loci with per-condition counts
loci_summary <- asm_mc_sig %>%
  dplyr::mutate(locus = paste(chr, pos_first_variant, sep = ":")) %>%
  dplyr::group_by(locus, chr, pos_first_variant, ref_allele, alt_allele) %>%
  dplyr::summarize(
    n_samples = dplyr::n(),
    n_ctrl = sum(condition == "Control"),
    n_mut = sum(condition == "Mutant"),
    mean_meth_diff = mean(methylation_diff, na.rm = TRUE),
    min_corr_p = min(corr_p_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(min_corr_p)

write.table(loci_summary,
            file.path(TABLES_DIR, "asm_mc_significant_loci.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved asm_mc_significant_loci.tsv (%d unique loci)\n", nrow(loci_summary)))

# Table 3: ASM-DMR overlap gene summary
if (exists("gene_summary") && nrow(gene_summary) > 0) {
  write.table(gene_summary,
              file.path(TABLES_DIR, "asm_dmr_overlap_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Saved asm_dmr_overlap_summary.tsv (%d genes)\n", nrow(gene_summary)))
} else {
  cat("  No ASM-DMR overlap genes to export\n")
}

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("ALLELE-SPECIFIC METHYLATION EXPLORATORY ANALYSIS SUMMARY\n")
cat("================================================================================\n")
cat(sprintf("mC PASS ASM sites:            %s (across 8 samples)\n",
            scales::comma(nrow(asm_mc_pass))))
cat(sprintf("mC significant (q < %.2f):     %s (%.1f%% of PASS)\n",
            Q_THRESHOLD, scales::comma(nrow(asm_mc_sig)),
            100 * nrow(asm_mc_sig) / max(nrow(asm_mc_pass), 1)))
cat(sprintf("hmC PASS ASM sites:           %s\n",
            scales::comma(nrow(asm_hmc_pass))))
cat(sprintf("hmC significant (q < %.2f):    %s (%.1f%% of PASS)\n",
            Q_THRESHOLD, scales::comma(nrow(asm_hmc_sig)),
            100 * nrow(asm_hmc_sig) / max(nrow(asm_hmc_pass), 1)))
cat(sprintf("Mutant mean sig sites (mC):   %.0f\n",
            condition_means$mean_count[condition_means$condition == "Mutant"]))
cat(sprintf("Control mean sig sites (mC):  %.0f\n",
            condition_means$mean_count[condition_means$condition == "Control"]))
cat(sprintf("Ratio (mut/ctrl):             %.2fx\n",
            condition_means$mean_count[condition_means$condition == "Mutant"] /
              condition_means$mean_count[condition_means$condition == "Control"]))
cat(sprintf("Shared ctrl/mut loci:         %s\n", scales::comma(n_shared)))
cat(sprintf("Control-only loci:            %s\n", scales::comma(n_ctrl_only)))
cat(sprintf("Mutant-only loci:             %s\n", scales::comma(n_mut_only)))
cat(sprintf("ASM-DMR overlap loci:         %d / %s (%.1f%%)\n",
            n_any_dmr_ol, scales::comma(n_asm),
            100 * n_any_dmr_ol / max(n_asm, 1)))
cat("================================================================================\n")
cat("\nSection 44 complete.\n\n")
