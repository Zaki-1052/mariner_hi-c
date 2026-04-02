# biomodal/downstream/scripts/viz_sections/section_42_max_significance_gene_list.R
# Section 42: Extract genes at maximum significance (q-value floor)
# Identifies genes whose q-values hit the numerical floor (~0), appearing at
# the -log10(q) = 300 ceiling in volcano plots. Exports merged mC/hmC table.

# Run from downstream/ directory
source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 42: MAXIMUM SIGNIFICANCE GENE LIST
# =============================================================================

cat("================================================================================\n")
cat("SECTION 42: MAXIMUM SIGNIFICANCE GENE LIST\n")
cat("================================================================================\n\n")

# Threshold: q-values that reach the numerical floor used by load_dmr_bed()
# load_dmr_bed() caps at pmax(q, 1e-300) => neg_log10_q = 300 for q < 1e-300
Q_FLOOR <- 1e-300

# ---- Extract floor-significance genes from each modality --------------------

mc_floor <- mc_dmr %>%
  dplyr::filter(dmr_qvalue < Q_FLOOR) %>%
  dplyr::select(chr, start, end, gene, num_contexts, mean_coverage,
                mean_mod_group1, mean_mod_group2,
                mod_difference, dmr_qvalue, direction) %>%
  dplyr::rename(
    mc_difference  = mod_difference,
    mc_qvalue      = dmr_qvalue,
    mc_direction   = direction,
    mc_ctrl_mean   = mean_mod_group1,
    mc_mut_mean    = mean_mod_group2
  )

hmc_floor <- hmc_dmr %>%
  dplyr::filter(dmr_qvalue < Q_FLOOR) %>%
  dplyr::select(chr, start, end, gene, mod_difference, dmr_qvalue, direction,
                mean_mod_group1, mean_mod_group2) %>%
  dplyr::rename(
    hmc_difference = mod_difference,
    hmc_qvalue     = dmr_qvalue,
    hmc_direction  = direction,
    hmc_ctrl_mean  = mean_mod_group1,
    hmc_mut_mean   = mean_mod_group2
  )

cat(sprintf("  5mC genes at q-value floor: %d\n", nrow(mc_floor)))
cat(sprintf("  5hmC genes at q-value floor: %d\n", nrow(hmc_floor)))

# ---- Merge on gene name ----------------------------------------------------

# Join on gene + coordinates to handle genes with multiple isoforms (e.g. Arhgap26)
merged <- dplyr::full_join(mc_floor, hmc_floor, by = c("gene", "chr", "start", "end"))

# Flag which modalities hit the floor
merged$mc_at_floor  <- !is.na(merged$mc_difference)
merged$hmc_at_floor <- !is.na(merged$hmc_difference)

# For genes at floor in both: classify coordinated pattern
merged$coordinated <- dplyr::case_when(
  merged$mc_at_floor & merged$hmc_at_floor &
    merged$mc_direction == "Hypermethylated" &
    merged$hmc_direction == "Hypomethylated" ~ "mC_up_hmC_down",
  merged$mc_at_floor & merged$hmc_at_floor &
    merged$mc_direction == "Hypomethylated" &
    merged$hmc_direction == "Hypermethylated" ~ "mC_down_hmC_up",
  merged$mc_at_floor & merged$hmc_at_floor ~ "same_direction",
  TRUE ~ "single_modality"
)

# Convert differences to percentage for readability
merged$mc_diff_pct  <- round(merged$mc_difference * 100, 2)
merged$hmc_diff_pct <- round(merged$hmc_difference * 100, 2)

# Sort by absolute mC effect size (largest first), then hmC
merged <- merged %>%
  dplyr::arrange(desc(abs(mc_diff_pct)), desc(abs(hmc_diff_pct)))

# ---- Summary statistics -----------------------------------------------------

n_both  <- sum(merged$mc_at_floor & merged$hmc_at_floor)
n_mc_only  <- sum(merged$mc_at_floor & !merged$hmc_at_floor)
n_hmc_only <- sum(!merged$mc_at_floor & merged$hmc_at_floor)
n_coordinated <- sum(merged$coordinated == "mC_up_hmC_down", na.rm = TRUE)

cat(sprintf("\n  At floor in both mC and hmC: %d\n", n_both))
cat(sprintf("  At floor in mC only: %d\n", n_mc_only))
cat(sprintf("  At floor in hmC only: %d\n", n_hmc_only))
cat(sprintf("  Coordinated mC-up/hmC-down (TET block): %d / %d (%.1f%%)\n",
            n_coordinated, n_both,
            ifelse(n_both > 0, n_coordinated / n_both * 100, 0)))

# ---- Export tables -----------------------------------------------------------

# Full merged table
output_merged <- merged %>%
  dplyr::select(gene, chr, start, end,
                mc_diff_pct, mc_direction, mc_ctrl_mean, mc_mut_mean, mc_at_floor,
                hmc_diff_pct, hmc_direction, hmc_ctrl_mean, hmc_mut_mean, hmc_at_floor,
                coordinated, num_contexts, mean_coverage)

write.table(output_merged,
            file.path(TABLES_DIR, "max_significance_genes_merged.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  Saved: %s (%d genes)\n",
            file.path(TABLES_DIR, "max_significance_genes_merged.tsv"),
            nrow(output_merged)))

# mC-only list (sorted by effect size)
mc_only_table <- mc_floor %>%
  dplyr::mutate(mc_diff_pct = round(mc_difference * 100, 2)) %>%
  dplyr::arrange(desc(mc_diff_pct)) %>%
  dplyr::select(gene, chr, start, end, mc_diff_pct, mc_direction,
                mc_ctrl_mean, mc_mut_mean, num_contexts, mean_coverage)

write.table(mc_only_table,
            file.path(TABLES_DIR, "max_significance_genes_mc.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: %s (%d genes)\n",
            file.path(TABLES_DIR, "max_significance_genes_mc.tsv"),
            nrow(mc_only_table)))

# hmC-only list (sorted by effect size)
hmc_only_table <- hmc_floor %>%
  dplyr::mutate(hmc_diff_pct = round(hmc_difference * 100, 2)) %>%
  dplyr::arrange(hmc_diff_pct) %>%
  dplyr::select(gene, hmc_diff_pct, hmc_direction,
                hmc_ctrl_mean, hmc_mut_mean)

write.table(hmc_only_table,
            file.path(TABLES_DIR, "max_significance_genes_hmc.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: %s (%d genes)\n",
            file.path(TABLES_DIR, "max_significance_genes_hmc.tsv"),
            nrow(hmc_only_table)))

# Gene-only list (just names, one per line) for pathway tools / GREAT / etc.
gene_list <- sort(unique(merged$gene))
writeLines(gene_list,
           file.path(TABLES_DIR, "max_significance_gene_names.txt"))
cat(sprintf("  Saved: %s (%d unique gene names)\n",
            file.path(TABLES_DIR, "max_significance_gene_names.txt"),
            length(gene_list)))

cat("\n✅ Section 42 complete.\n")
