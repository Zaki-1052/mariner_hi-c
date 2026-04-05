# biomodal/downstream/scripts/viz_sections/compare_shallow_vs_deep.R
# Standalone comparison: deep-8 (run-4, no sex covariate) vs sex-covariate (run-5, sex covariate)
# NOT a section_*.R — will not be picked up by run_all_sections.sh
# Run manually: cd downstream/ && Rscript scripts/viz_sections/compare_shallow_vs_deep.R

# Source shared config for packages, helpers, colors, and deep-8 DATA_PATHS
source("scripts/viz_sections/_shared_config.R")

cat("================================================================================\n")
cat("DEEP-8 (RUN-4, no sex covariate) vs SEX-COVARIATE (RUN-5) COMPARISON\n")
cat("================================================================================\n\n")

# =============================================================================
# DEEP-8 PATHS (run-4, 8 samples, no sex covariate, hardcoded timestamps)
# =============================================================================

SHALLOW_PATHS <- list(
  mc_dmr = file.path(BASE_DIR, "modality/outputs/run-4/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260329_203415/DMR_mc_control__mutant_20260329_203415.bed"),
  hmc_dmr = file.path(BASE_DIR, "modality/outputs/run-4/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260329_203415/DMR_hmc_control__mutant_20260329_203415.bed"),
  cpg_islands_mc = file.path(BASE_DIR, "modality/outputs/run-4/outputs_CG/Results/gencode.vM25.mouse.cpg_islands.annotation/DMR_20260329_202138/DMR_mc_control__mutant_20260329_202138.bed"),
  cpg_shores_mc = file.path(BASE_DIR, "modality/outputs/run-4/outputs_CG/Results/gencode.vM25.mouse.cpg_shores.annotation/DMR_20260329_202951/DMR_mc_control__mutant_20260329_202951.bed"),
  cpg_shelves_mc = file.path(BASE_DIR, "modality/outputs/run-4/outputs_CG/Results/gencode.vM25.mouse.cpg_shelves.annotation/DMR_20260329_202514/DMR_mc_control__mutant_20260329_202514.bed"),
  promoters_mc = file.path(BASE_DIR, "modality/outputs/run-4/outputs_CG/Results/gencode.vM25.mouse.promoters.annotation/DMR_20260329_203816/DMR_mc_control__mutant_20260329_203816.bed"),
  tss_mc = file.path(BASE_DIR, "modality/outputs/run-4/outputs_CG/Results/gencode.vM25.mouse.tss_region.annotation/DMR_20260329_204205/DMR_mc_control__mutant_20260329_204205.bed"),
  cpg_islands_hmc = file.path(BASE_DIR, "modality/outputs/run-4/outputs_CG/Results/gencode.vM25.mouse.cpg_islands.annotation/DMR_20260329_202138/DMR_hmc_control__mutant_20260329_202138.bed"),
  cpg_shores_hmc = file.path(BASE_DIR, "modality/outputs/run-4/outputs_CG/Results/gencode.vM25.mouse.cpg_shores.annotation/DMR_20260329_202951/DMR_hmc_control__mutant_20260329_202951.bed"),
  cpg_shelves_hmc = file.path(BASE_DIR, "modality/outputs/run-4/outputs_CG/Results/gencode.vM25.mouse.cpg_shelves.annotation/DMR_20260329_202514/DMR_hmc_control__mutant_20260329_202514.bed"),
  promoters_hmc = file.path(BASE_DIR, "modality/outputs/run-4/outputs_CG/Results/gencode.vM25.mouse.promoters.annotation/DMR_20260329_203816/DMR_hmc_control__mutant_20260329_203816.bed"),
  tss_hmc = file.path(BASE_DIR, "modality/outputs/run-4/outputs_CG/Results/gencode.vM25.mouse.tss_region.annotation/DMR_20260329_204205/DMR_hmc_control__mutant_20260329_204205.bed")
)

# Sex-covariate paths come from DATA_PATHS (loaded via _shared_config.R, now run-5)
DEEP_PATHS <- DATA_PATHS

# =============================================================================
# VALIDATE ALL PATHS (no fallbacks per project rules)
# =============================================================================

cat("Validating deep-8 (run-4, no sex covariate) paths...\n")
for (name in names(SHALLOW_PATHS)) {
  if (!file.exists(SHALLOW_PATHS[[name]])) {
    stop(sprintf("Deep-8 (run-4) file missing: %s\n  Path: %s", name, SHALLOW_PATHS[[name]]))
  }
}
cat("  All deep-8 (run-4) paths validated.\n")

cat("Validating sex-covariate (run-5) paths...\n")
deep_dmr_keys <- c("mc_dmr", "hmc_dmr",
                    "cpg_islands_mc", "cpg_shores_mc", "cpg_shelves_mc",
                    "promoters_mc", "tss_mc",
                    "cpg_islands_hmc", "cpg_shores_hmc", "cpg_shelves_hmc",
                    "promoters_hmc", "tss_hmc")
for (name in deep_dmr_keys) {
  if (!file.exists(DEEP_PATHS[[name]])) {
    stop(sprintf("Sex-covariate (run-5) file missing: %s\n  Path: %s", name, DEEP_PATHS[[name]]))
  }
}
cat("  All sex-covariate (run-5) paths validated.\n\n")

# Color palette for run comparison
RUN_COLORS <- c("No Sex Cov (run-4)" = "#FF8C00", "Sex Cov (run-5)" = "#1E90FF")

# Output subdirectory for comparison plots
COMP_DIR <- file.path(OUTPUT_DIR, "comparison_shallow_vs_deep")
dir.create(COMP_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# LOAD ALL DATA
# =============================================================================

cat("Loading deep-8 (run-4) and sex-covariate (run-5) DMR data...\n")

# Deduplicate DMR BED: some genes have multiple regions annotated to the same
# gene name. Keep the entry with the lowest q-value per gene.
dedup_by_gene <- function(df) {
  df %>%
    dplyr::group_by(gene) %>%
    dplyr::slice_min(dmr_qvalue, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
}

# Gene body DMRs
shallow_mc <- dedup_by_gene(load_dmr_bed(SHALLOW_PATHS$mc_dmr))
shallow_hmc <- dedup_by_gene(load_dmr_bed(SHALLOW_PATHS$hmc_dmr))
deep_mc <- dedup_by_gene(mc_dmr)  # mc_dmr loaded by _shared_config.R
deep_hmc <- dedup_by_gene(hmc_dmr)

stopifnot(!is.null(shallow_mc), !is.null(shallow_hmc),
          !is.null(deep_mc), !is.null(deep_hmc))

cat(sprintf("  No sex cov (run-4) gene body: %d mC genes, %d hmC genes (deduplicated)\n",
            nrow(shallow_mc), nrow(shallow_hmc)))
cat(sprintf("  Sex cov (run-5) gene body: %d mC genes, %d hmC genes (deduplicated)\n",
            nrow(deep_mc), nrow(deep_hmc)))

# Regional DMRs — organized as list keyed by region name
REGION_KEYS <- list(
  "Gene bodies"  = list(mc = "mc_dmr",         hmc = "hmc_dmr"),
  "CpG islands"  = list(mc = "cpg_islands_mc",  hmc = "cpg_islands_hmc"),
  "CpG shores"   = list(mc = "cpg_shores_mc",   hmc = "cpg_shores_hmc"),
  "CpG shelves"  = list(mc = "cpg_shelves_mc",  hmc = "cpg_shelves_hmc"),
  "Promoters"    = list(mc = "promoters_mc",     hmc = "promoters_hmc"),
  "TSS regions"  = list(mc = "tss_mc",           hmc = "tss_hmc")
)

# Load all regional data for both runs
shallow_regional <- list()
deep_regional <- list()
for (region in names(REGION_KEYS)) {
  mc_key <- REGION_KEYS[[region]]$mc
  hmc_key <- REGION_KEYS[[region]]$hmc
  shallow_regional[[region]] <- list(
    mc = dedup_by_gene(load_dmr_bed(SHALLOW_PATHS[[mc_key]])),
    hmc = dedup_by_gene(load_dmr_bed(SHALLOW_PATHS[[hmc_key]]))
  )
  deep_regional[[region]] <- list(
    mc = dedup_by_gene(load_dmr_bed(DEEP_PATHS[[mc_key]])),
    hmc = dedup_by_gene(load_dmr_bed(DEEP_PATHS[[hmc_key]]))
  )
}

cat("  All regional data loaded.\n\n")

# =============================================================================
# ANALYSIS A: DMR Count Comparison by Region
# =============================================================================

cat("--- Analysis A: DMR count comparison by region ---\n")

count_rows <- list()
for (region in names(REGION_KEYS)) {
  for (mod in c("mC", "hmC")) {
    mod_key <- ifelse(mod == "mC", "mc", "hmc")
    s_df <- shallow_regional[[region]][[mod_key]]
    d_df <- deep_regional[[region]][[mod_key]]
    count_rows[[length(count_rows) + 1]] <- data.frame(
      Region = region,
      Modification = paste0("5", mod),
      Run = "No Sex Cov (run-4)",
      Significant = sum(s_df$significant),
      Total = nrow(s_df),
      stringsAsFactors = FALSE
    )
    count_rows[[length(count_rows) + 1]] <- data.frame(
      Region = region,
      Modification = paste0("5", mod),
      Run = "Sex Cov (run-5)",
      Significant = sum(d_df$significant),
      Total = nrow(d_df),
      stringsAsFactors = FALSE
    )
  }
}
dmr_counts <- do.call(rbind, count_rows)
dmr_counts$Percentage <- 100 * dmr_counts$Significant / dmr_counts$Total
dmr_counts$Run <- factor(dmr_counts$Run, levels = c("No Sex Cov (run-4)", "Sex Cov (run-5)"))

# Save counts table
write.table(dmr_counts, file.path(TABLES_DIR, "dmr_counts_shallow_vs_deep.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

p_counts <- ggplot(dmr_counts, aes(x = Region, y = Significant, fill = Run)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = Significant),
            position = position_dodge(width = 0.8), vjust = -0.3, size = 3) +
  facet_wrap(~Modification, scales = "free_y") +
  scale_fill_manual(values = RUN_COLORS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Significant DMRs: No Sex Cov vs Sex Cov",
    subtitle = "DMR count by genomic region and modification type (q < 0.05)",
    x = "", y = "Number of Significant DMRs"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_multiformat_ggplot(p_counts, file.path(COMP_DIR, "20a_dmr_counts_by_region"),
                        width = 14, height = 8)
cat("  Saved 20a_dmr_counts_by_region\n")

# =============================================================================
# ANALYSIS B: Gene Status Classification
# =============================================================================

cat("\n--- Analysis B: Gene status classification ---\n")

classify_gene_status <- function(shallow_df, deep_df, mod_label) {
  # Full outer join on gene name
  merged <- full_join(
    shallow_df %>% dplyr::select(gene,
      shallow_q = dmr_qvalue, shallow_effect = mod_difference,
      shallow_coverage = mean_coverage),
    deep_df %>% dplyr::select(gene,
      deep_q = dmr_qvalue, deep_effect = mod_difference,
      deep_coverage = mean_coverage),
    by = "gene"
  )

  merged$status <- case_when(
    !is.na(merged$shallow_q) & !is.na(merged$deep_q) &
      merged$shallow_q < Q_THRESHOLD & merged$deep_q < Q_THRESHOLD ~ "Retained",
    !is.na(merged$shallow_q) & !is.na(merged$deep_q) &
      merged$shallow_q >= Q_THRESHOLD & merged$deep_q < Q_THRESHOLD ~ "Newly significant",
    !is.na(merged$shallow_q) & !is.na(merged$deep_q) &
      merged$shallow_q < Q_THRESHOLD & merged$deep_q >= Q_THRESHOLD ~ "Lost significance",
    !is.na(merged$shallow_q) & !is.na(merged$deep_q) &
      merged$shallow_q >= Q_THRESHOLD & merged$deep_q >= Q_THRESHOLD ~ "Never significant",
    is.na(merged$shallow_q) & !is.na(merged$deep_q) ~ "Run-5 only",
    !is.na(merged$shallow_q) & is.na(merged$deep_q) ~ "Run-4 only",
    TRUE ~ "Unknown"
  )
  merged$modification <- mod_label
  return(merged)
}

mc_status <- classify_gene_status(shallow_mc, deep_mc, "5mC")
hmc_status <- classify_gene_status(shallow_hmc, deep_hmc, "5hmC")

# Print summary
for (mod_df in list(mc_status, hmc_status)) {
  mod_label <- mod_df$modification[1]
  cat(sprintf("\n  %s gene status:\n", mod_label))
  status_tbl <- table(mod_df$status)
  for (s in names(status_tbl)) {
    cat(sprintf("    %-20s %d\n", s, status_tbl[s]))
  }
}

# Stacked bar chart
status_summary <- bind_rows(mc_status, hmc_status) %>%
  dplyr::filter(status != "Unknown") %>%
  group_by(modification, status) %>%
  summarise(n = n(), .groups = "drop")

status_order <- c("Retained", "Newly significant", "Lost significance",
                  "Never significant", "Run-5 only", "Run-4 only")
status_colors <- c(
  "Retained" = "#2CA02C",
  "Newly significant" = "#1E90FF",
  "Lost significance" = "#FF6347",
  "Never significant" = "grey70",
  "Run-5 only" = "#9467BD",
  "Run-4 only" = "#FF8C00"
)
status_summary$status <- factor(status_summary$status, levels = status_order)

p_status <- ggplot(status_summary, aes(x = modification, y = n, fill = status)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3.5) +
  scale_fill_manual(values = status_colors, name = "Gene Status") +
  labs(
    title = "Gene Status: No Sex Cov vs Sex Cov",
    subtitle = "How gene significance changed with sex covariate",
    x = "Modification", y = "Number of Genes"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_status, file.path(COMP_DIR, "20b_gene_status_summary"),
                        width = 10, height = 8)
cat("\n  Saved 20b_gene_status_summary\n")

# =============================================================================
# ANALYSIS C: Newly Significant Genes Table
# =============================================================================

cat("\n--- Analysis C: Newly significant genes ---\n")

newly_sig_mc <- mc_status %>%
  dplyr::filter(status == "Newly significant") %>%
  dplyr::arrange(deep_q) %>%
  dplyr::mutate(modification = "5mC")

newly_sig_hmc <- hmc_status %>%
  dplyr::filter(status == "Newly significant") %>%
  dplyr::arrange(deep_q) %>%
  dplyr::mutate(modification = "5hmC")

newly_sig <- bind_rows(newly_sig_mc, newly_sig_hmc) %>%
  dplyr::select(gene, modification, shallow_q, deep_q, shallow_effect, deep_effect,
                shallow_coverage, deep_coverage, status)

write.table(newly_sig, file.path(TABLES_DIR, "newly_significant_deep_seq.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  %d newly significant mC genes\n", nrow(newly_sig_mc)))
cat(sprintf("  %d newly significant hmC genes\n", nrow(newly_sig_hmc)))
cat(sprintf("  Saved newly_significant_deep_seq.tsv (%d rows)\n", nrow(newly_sig)))

# =============================================================================
# ANALYSIS D: Effect Size Scatter
# =============================================================================

cat("\n--- Analysis D: Effect size scatter ---\n")

# Genes present in both runs (shared)
mc_shared <- mc_status %>%
  dplyr::filter(!is.na(shallow_effect) & !is.na(deep_effect))
hmc_shared <- hmc_status %>%
  dplyr::filter(!is.na(shallow_effect) & !is.na(deep_effect))

scatter_data <- bind_rows(
  mc_shared %>% dplyr::select(gene, shallow_effect, deep_effect, status) %>%
    dplyr::mutate(modification = "5mC"),
  hmc_shared %>% dplyr::select(gene, shallow_effect, deep_effect, status) %>%
    dplyr::mutate(modification = "5hmC")
)

# Correlation per modification
mc_cor <- cor(mc_shared$shallow_effect, mc_shared$deep_effect,
              use = "complete.obs", method = "pearson")
hmc_cor <- cor(hmc_shared$shallow_effect, hmc_shared$deep_effect,
               use = "complete.obs", method = "pearson")

p_effect <- ggplot(scatter_data,
                   aes(x = shallow_effect * 100, y = deep_effect * 100)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = status), alpha = 0.4, size = 1.2) +
  scale_color_manual(values = status_colors, name = "Gene Status") +
  facet_wrap(~modification, scales = "free") +
  labs(
    title = "Effect Size Stability: No Sex Cov vs Sex Cov",
    subtitle = sprintf("5mC r=%.3f | 5hmC r=%.3f | Dashed line = identity",
                       mc_cor, hmc_cor),
    x = "No Sex Cov Effect Size (%)",
    y = "Sex Cov Effect Size (%)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_effect, file.path(COMP_DIR, "20d_effect_size_scatter"),
                        width = 14, height = 7)
cat(sprintf("  5mC correlation: %.4f (%d genes)\n", mc_cor, nrow(mc_shared)))
cat(sprintf("  5hmC correlation: %.4f (%d genes)\n", hmc_cor, nrow(hmc_shared)))
cat("  Saved 20d_effect_size_scatter\n")

# =============================================================================
# ANALYSIS E: Q-value Improvement
# =============================================================================

cat("\n--- Analysis E: Q-value improvement ---\n")

# For genes present in both runs, compute delta q-value
qval_data <- bind_rows(
  mc_shared %>% dplyr::select(gene, shallow_q, deep_q) %>%
    dplyr::mutate(modification = "5mC"),
  hmc_shared %>% dplyr::select(gene, shallow_q, deep_q) %>%
    dplyr::mutate(modification = "5hmC")
) %>%
  dplyr::mutate(
    delta_q = deep_q - shallow_q,
    q_improved = deep_q < shallow_q,
    shallow_nlog10 = -log10(pmax(shallow_q, 1e-300)),
    deep_nlog10 = -log10(pmax(deep_q, 1e-300))
  )

pct_improved <- 100 * mean(qval_data$q_improved, na.rm = TRUE)

# Histogram of delta q-values
p_delta_hist <- ggplot(qval_data, aes(x = delta_q, fill = modification)) +
  geom_histogram(bins = 80, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(values = COLORS$methylation, name = "Modification") +
  labs(
    title = "Change in Q-values: Sex Cov - No Sex Cov",
    subtitle = sprintf("%.1f%% of genes improved (negative = better in sex cov)", pct_improved),
    x = "Delta Q-value (sex cov - no sex cov)",
    y = "Number of Genes"
  ) +
  theme_biomodal()

# -log10(q) scatter
p_nlog10 <- ggplot(qval_data, aes(x = shallow_nlog10, y = deep_nlog10)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = modification), alpha = 0.3, size = 1) +
  scale_color_manual(values = COLORS$methylation, name = "Modification") +
  labs(
    title = "-log10(Q) Comparison",
    subtitle = "Points above diagonal = more significant in sex cov",
    x = "-log10(q) No Sex Cov",
    y = "-log10(q) Sex Cov"
  ) +
  theme_biomodal()

p_qval_combined <- p_delta_hist / p_nlog10 +
  plot_annotation(
    title = "Statistical Power: Effect of Sex Covariate",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )

save_multiformat_ggplot(p_qval_combined, file.path(COMP_DIR, "20e_qvalue_improvement"),
                        width = 12, height = 12)
cat(sprintf("  %.1f%% of genes have improved q-values in sex cov (run-5) vs no sex cov (run-4)\n", pct_improved))
cat("  Saved 20e_qvalue_improvement\n")

# =============================================================================
# ANALYSIS F: Coverage Comparison
# =============================================================================

cat("\n--- Analysis F: Coverage comparison ---\n")

cov_data <- mc_status %>%
  dplyr::filter(!is.na(shallow_coverage) & !is.na(deep_coverage)) %>%
  dplyr::mutate(
    sig_gain = status == "Newly significant",
    sig_label = ifelse(sig_gain, "Gained significance", "Other")
  )

cov_cor <- cor(cov_data$shallow_coverage, cov_data$deep_coverage,
               use = "complete.obs")
median_shallow <- median(cov_data$shallow_coverage, na.rm = TRUE)
median_deep <- median(cov_data$deep_coverage, na.rm = TRUE)

p_coverage <- ggplot(cov_data, aes(x = shallow_coverage, y = deep_coverage)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = sig_label), alpha = 0.3, size = 1.2) +
  scale_color_manual(
    values = c("Gained significance" = "#1E90FF", "Other" = "grey60"),
    name = "Status"
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Coverage Comparison: No Sex Cov vs Sex Cov (5mC)",
    subtitle = sprintf("Median coverage: %.1f (no sex cov) vs %.1f (sex cov) | r=%.3f",
                       median_shallow, median_deep, cov_cor),
    x = "Mean Coverage (No Sex Cov, log10)",
    y = "Mean Coverage (Sex Cov, log10)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_coverage, file.path(COMP_DIR, "20f_coverage_comparison"),
                        width = 10, height = 9)
cat(sprintf("  Median coverage: %.1f (no sex cov) vs %.1f (sex cov)\n",
            median_shallow, median_deep))
cat("  Saved 20f_coverage_comparison\n")

# =============================================================================
# ANALYSIS G: Coordinated Pattern Stability
# =============================================================================

cat("\n--- Analysis G: Coordinated pattern stability ---\n")

# Identify coordinated genes in each run
get_coordinated_genes <- function(mc_df, hmc_df) {
  mc_sig <- mc_df %>%
    dplyr::filter(significant) %>%
    dplyr::select(gene, mc_diff = mod_difference)
  hmc_sig <- hmc_df %>%
    dplyr::filter(significant) %>%
    dplyr::select(gene, hmc_diff = mod_difference)
  coord <- inner_join(mc_sig, hmc_sig, by = "gene") %>%
    dplyr::mutate(coordinated = mc_diff > 0 & hmc_diff < 0)
  return(coord)
}

shallow_coord <- get_coordinated_genes(shallow_mc, shallow_hmc)
deep_coord <- get_coordinated_genes(deep_mc, deep_hmc)

shallow_coord_genes <- shallow_coord %>% dplyr::filter(coordinated) %>% pull(gene)
deep_coord_genes <- deep_coord %>% dplyr::filter(coordinated) %>% pull(gene)

# Compute overlap
shared_coord <- intersect(shallow_coord_genes, deep_coord_genes)
shallow_only_coord <- setdiff(shallow_coord_genes, deep_coord_genes)
deep_only_coord <- setdiff(deep_coord_genes, shallow_coord_genes)

cat(sprintf("  No sex cov (run-4) coordinated genes: %d (%.1f%% of co-significant)\n",
            length(shallow_coord_genes),
            100 * mean(shallow_coord$coordinated)))
cat(sprintf("  Sex cov (run-5) coordinated genes: %d (%.1f%% of co-significant)\n",
            length(deep_coord_genes),
            100 * mean(deep_coord$coordinated)))
cat(sprintf("  Shared:       %d\n", length(shared_coord)))
cat(sprintf("  Run-4 only:   %d\n", length(shallow_only_coord)))
cat(sprintf("  Run-5 only:   %d\n", length(deep_only_coord)))

# Bar chart: coordinated pattern percentage comparison
pattern_pct <- data.frame(
  Run = factor(c("No Sex Cov (run-4)", "Sex Cov (run-5)"),
               levels = c("No Sex Cov (run-4)", "Sex Cov (run-5)")),
  Co_significant = c(nrow(shallow_coord), nrow(deep_coord)),
  Coordinated = c(length(shallow_coord_genes), length(deep_coord_genes)),
  Pct_coordinated = c(100 * mean(shallow_coord$coordinated),
                      100 * mean(deep_coord$coordinated))
)

p_pattern_bar <- ggplot(pattern_pct, aes(x = Run, y = Pct_coordinated, fill = Run)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", Pct_coordinated,
                                Coordinated, Co_significant)),
            vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = RUN_COLORS) +
  scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
  labs(
    title = "Coordinated Pattern Stability",
    subtitle = "Percentage of co-significant genes with mC(+)/hmC(-) pattern",
    x = "", y = "Percentage (%)"
  ) +
  theme_biomodal() +
  theme(legend.position = "none")

# Venn diagram of coordinated gene overlap
venn_list <- list(
  "No Sex Cov (run-4)" = shallow_coord_genes,
  "Sex Cov (run-5)" = deep_coord_genes
)

p_venn <- ggVennDiagram(venn_list, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#1E90FF") +
  labs(
    title = "Overlap of Coordinated Genes",
    subtitle = sprintf("%d shared, %d run-4-only, %d run-5-only",
                       length(shared_coord), length(shallow_only_coord),
                       length(deep_only_coord))
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

p_coord_combined <- p_pattern_bar | p_venn

save_multiformat_ggplot(p_coord_combined,
                        file.path(COMP_DIR, "20g_coordinated_pattern_stability"),
                        width = 14, height = 7)
cat("  Saved 20g_coordinated_pattern_stability\n")

# =============================================================================
# MASTER SUMMARY TABLE
# =============================================================================

cat("\n--- Creating master summary table ---\n")

# Build master table from mc_status and hmc_status
master <- full_join(
  mc_status %>% dplyr::select(
    gene,
    shallow_mc_q = shallow_q, deep_mc_q = deep_q,
    shallow_mc_effect = shallow_effect, deep_mc_effect = deep_effect,
    shallow_mc_coverage = shallow_coverage, deep_mc_coverage = deep_coverage,
    mc_status = status
  ),
  hmc_status %>% dplyr::select(
    gene,
    shallow_hmc_q = shallow_q, deep_hmc_q = deep_q,
    shallow_hmc_effect = shallow_effect, deep_hmc_effect = deep_effect,
    hmc_status = status
  ),
  by = "gene"
) %>%
  dplyr::mutate(
    mc_effect_delta = deep_mc_effect - shallow_mc_effect,
    hmc_effect_delta = deep_hmc_effect - shallow_hmc_effect
  ) %>%
  dplyr::arrange(deep_mc_q)

write.table(master,
            file.path(TABLES_DIR, "shallow_vs_deep_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("  Master table: %d genes\n", nrow(master)))
cat(sprintf("  Columns: %s\n", paste(colnames(master), collapse = ", ")))
cat("  Saved shallow_vs_deep_comparison.tsv\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("COMPARISON COMPLETE\n")
cat("================================================================================\n\n")
cat("Output plots:  ", COMP_DIR, "\n")
cat("Output tables: ", TABLES_DIR, "\n\n")
cat("Key findings:\n")
cat(sprintf("  mC significant: %d (no sex cov) -> %d (sex cov)\n",
            sum(shallow_mc$significant), sum(deep_mc$significant)))
cat(sprintf("  hmC significant: %d (no sex cov) -> %d (sex cov)\n",
            sum(shallow_hmc$significant), sum(deep_hmc$significant)))
cat(sprintf("  Newly significant mC genes: %d\n", nrow(newly_sig_mc)))
cat(sprintf("  Newly significant hmC genes: %d\n", nrow(newly_sig_hmc)))
cat(sprintf("  Effect size correlation: mC r=%.3f, hmC r=%.3f\n", mc_cor, hmc_cor))
cat(sprintf("  Coordinated pattern: %.1f%% (no sex cov) -> %.1f%% (sex cov)\n",
            100 * mean(shallow_coord$coordinated),
            100 * mean(deep_coord$coordinated)))
cat(sprintf("  Coverage: %.1fx (no sex cov) -> %.1fx (sex cov) median\n",
            median_shallow, median_deep))
cat("\n")
