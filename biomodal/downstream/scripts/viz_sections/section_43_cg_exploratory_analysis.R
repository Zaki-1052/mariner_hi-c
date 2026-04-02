# biomodal/downstream/scripts/viz_sections/section_43_cg_exploratory_analysis.R
# Section 43: CG Context Exploratory Analysis
# Chromosome distribution, direction breakdown, mC/hmC overlap, and regional analyses
# for the primary CG methylation context (run-4, deep-seq, 8 samples)
# Standalone script - sources shared config for all dependencies and data

# Run from downstream/ directory
source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 43: CG CONTEXT EXPLORATORY ANALYSIS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 43: CG CONTEXT EXPLORATORY ANALYSIS\n")
cat("================================================================================\n\n")

# -----------------------------------------------------------------------------
# SETUP: Extract significant subsets from shared-config-loaded data
# -----------------------------------------------------------------------------

# mc_dmr and hmc_dmr are loaded by _shared_config.R (gene body CG DMRs)
stopifnot("mc_dmr not loaded" = exists("mc_dmr") && !is.null(mc_dmr))
stopifnot("hmc_dmr not loaded" = exists("hmc_dmr") && !is.null(hmc_dmr))

mc_sig <- mc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

hmc_sig <- hmc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

n_mc_hyper <- sum(mc_sig$mod_difference > 0)
n_mc_hypo  <- sum(mc_sig$mod_difference <= 0)
n_hmc_hyper <- sum(hmc_sig$mod_difference > 0)
n_hmc_hypo  <- sum(hmc_sig$mod_difference <= 0)

cat(sprintf("CG mC significant: %d genes (%d hyper, %d hypo)\n",
            nrow(mc_sig), n_mc_hyper, n_mc_hypo))
cat(sprintf("CG hmC significant: %d genes (%d hyper, %d hypo)\n",
            nrow(hmc_sig), n_hmc_hyper, n_hmc_hypo))

# Create section output subdirectory
SECTION_DIR <- "43_cg_exploratory"
dir.create(file.path(OUTPUT_DIR, SECTION_DIR), recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# FIGURE 43a: DIRECTION BREAKDOWN (mC + hmC)
# =============================================================================

cat("\n--- Figure 43a: Direction Breakdown ---\n")

direction_df <- data.frame(
  context = c("CG mC", "CG mC", "CG hmC", "CG hmC"),
  direction = c("Hypermethylated", "Hypomethylated", "Hypermethylated", "Hypomethylated"),
  count = c(n_mc_hyper, n_mc_hypo, n_hmc_hyper, n_hmc_hypo),
  stringsAsFactors = FALSE
)

direction_df <- direction_df %>%
  dplyr::group_by(context) %>%
  dplyr::mutate(
    total = sum(count),
    pct = count / total * 100,
    label = sprintf("%d (%.1f%%)", count, pct)
  ) %>%
  dplyr::ungroup()

p_direction <- ggplot(direction_df, aes(x = context, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = label), position = position_dodge(width = 0.7),
            vjust = -0.3, size = 4) +
  scale_fill_manual(values = COLORS$direction, name = "Direction") +
  labs(
    title = "CG Significant DMR Direction Breakdown",
    subtitle = "Gene body DMRs | BAP1-KO vs Control",
    x = "", y = "Number of Genes"
  ) +
  theme_biomodal() +
  theme(legend.position = "bottom")

save_multiformat_ggplot(p_direction, file.path(OUTPUT_DIR, SECTION_DIR, "43a_cg_direction_breakdown"),
                        width = 8, height = 6)
cat("  Saved 43a_cg_direction_breakdown\n")

# =============================================================================
# FIGURE 43b: mC vs hmC VENN DIAGRAM
# =============================================================================

cat("\n--- Figure 43b: mC vs hmC Venn Diagram ---\n")

mc_genes <- mc_sig$gene
hmc_genes <- hmc_sig$gene

venn_sets <- list(
  "5mC" = mc_genes,
  "5hmC" = hmc_genes
)

n_mc_only <- length(setdiff(mc_genes, hmc_genes))
n_hmc_only <- length(setdiff(hmc_genes, mc_genes))
n_both <- length(intersect(mc_genes, hmc_genes))

p_venn <- ggVennDiagram(venn_sets, label_alpha = 0, set_size = 5) +
  scale_fill_gradient(low = "white", high = "#E41A1C") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  labs(
    title = "CG mC vs hmC: Significant Gene Overlap",
    subtitle = sprintf("%d overlap | %d mC-only | %d hmC-only",
                       n_both, n_mc_only, n_hmc_only)
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

save_multiformat_ggplot(p_venn, file.path(OUTPUT_DIR, SECTION_DIR, "43b_cg_mc_hmc_venn"),
                        width = 8, height = 7)
cat(sprintf("  Saved 43b_cg_mc_hmc_venn (%d overlap, %d mC-only, %d hmC-only)\n",
            n_both, n_mc_only, n_hmc_only))

# =============================================================================
# FIGURE 43c: CHROMOSOME DISTRIBUTION - mC DMRs
# =============================================================================

cat("\n--- Figure 43c: Chromosome Distribution (mC) ---\n")

# Helper: compute chromosome distribution with Fisher's test
compute_chr_distribution <- function(sig_df, all_df, label = "mC") {
  chr_sig <- sig_df %>%
    dplyr::count(chr, name = "sig_count")

  chr_total <- all_df %>%
    dplyr::count(chr, name = "total_genes")

  chr_dist <- chr_total %>%
    dplyr::left_join(chr_sig, by = "chr") %>%
    dplyr::mutate(
      sig_count = replace(sig_count, is.na(sig_count), 0),
      sig_rate = sig_count / total_genes,
      expected = total_genes / sum(total_genes) * sum(sig_count),
      obs_exp_ratio = sig_count / expected,
      chr_label = gsub("chr", "", chr),
      chr_num = suppressWarnings(as.numeric(chr_label)),
      chr_order = ifelse(is.na(chr_num), 100 + match(chr_label, c("X", "Y", "M")), chr_num)
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
    total_on <- chr_dist$total_genes[i]
    total_off <- sum(chr_dist$total_genes) - total_on
    nonsig_on  <- total_on - sig_on
    nonsig_off <- total_off - sig_off

    mat <- matrix(c(sig_on, nonsig_on, sig_off, nonsig_off), nrow = 2)
    ft <- fisher.test(mat)
    chr_dist$fisher_p[i] <- ft$p.value
    chr_dist$fisher_or[i] <- ft$estimate
  }

  # Mark significantly enriched/depleted chromosomes
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

  return(chr_dist)
}

# Compute mC chromosome distribution
mc_chr_dist <- compute_chr_distribution(mc_sig, mc_dmr, "mC")

# Summary of enriched chromosomes
mc_enriched <- mc_chr_dist %>% dplyr::filter(enrichment == "Enriched")
mc_depleted <- mc_chr_dist %>% dplyr::filter(enrichment == "Depleted")
cat(sprintf("  mC enriched chromosomes: %s\n",
            ifelse(nrow(mc_enriched) > 0,
                   paste(mc_enriched$chr, collapse = ", "), "none")))
cat(sprintf("  mC depleted chromosomes: %s\n",
            ifelse(nrow(mc_depleted) > 0,
                   paste(mc_depleted$chr, collapse = ", "), "none")))

p_chr_mc <- ggplot(mc_chr_dist, aes(x = chr, y = sig_count)) +
  geom_bar(aes(fill = enrichment), stat = "identity", width = 0.7) +
  geom_point(aes(y = expected), color = "black", shape = 4, size = 3) +
  geom_text(aes(label = fisher_label), vjust = -0.5, size = 3, color = "#E41A1C",
            fontface = "bold") +
  scale_fill_manual(
    values = c("Enriched" = "#D7191C", "Depleted" = "#2C7BB6", "Not significant" = "grey70"),
    name = "Enrichment"
  ) +
  labs(
    title = "Chromosome Distribution of Significant CG mC DMRs",
    subtitle = sprintf("%d significant gene body DMRs | X = expected count (genome-wide rate)",
                       nrow(mc_sig)),
    x = "Chromosome", y = "Significant DMR Count"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_multiformat_ggplot(p_chr_mc, file.path(OUTPUT_DIR, SECTION_DIR, "43c_cg_mc_chromosome_distribution"),
                        width = 12, height = 7)
cat("  Saved 43c_cg_mc_chromosome_distribution\n")

# =============================================================================
# FIGURE 43d: CHROMOSOME DISTRIBUTION - hmC DMRs
# =============================================================================

cat("\n--- Figure 43d: Chromosome Distribution (hmC) ---\n")

hmc_chr_dist <- compute_chr_distribution(hmc_sig, hmc_dmr, "hmC")

hmc_enriched <- hmc_chr_dist %>% dplyr::filter(enrichment == "Enriched")
hmc_depleted <- hmc_chr_dist %>% dplyr::filter(enrichment == "Depleted")
cat(sprintf("  hmC enriched chromosomes: %s\n",
            ifelse(nrow(hmc_enriched) > 0,
                   paste(hmc_enriched$chr, collapse = ", "), "none")))
cat(sprintf("  hmC depleted chromosomes: %s\n",
            ifelse(nrow(hmc_depleted) > 0,
                   paste(hmc_depleted$chr, collapse = ", "), "none")))

p_chr_hmc <- ggplot(hmc_chr_dist, aes(x = chr, y = sig_count)) +
  geom_bar(aes(fill = enrichment), stat = "identity", width = 0.7) +
  geom_point(aes(y = expected), color = "black", shape = 4, size = 3) +
  geom_text(aes(label = fisher_label), vjust = -0.5, size = 3, color = "#E41A1C",
            fontface = "bold") +
  scale_fill_manual(
    values = c("Enriched" = "#D7191C", "Depleted" = "#2C7BB6", "Not significant" = "grey70"),
    name = "Enrichment"
  ) +
  labs(
    title = "Chromosome Distribution of Significant CG hmC DMRs",
    subtitle = sprintf("%d significant gene body DMRs | X = expected count (genome-wide rate)",
                       nrow(hmc_sig)),
    x = "Chromosome", y = "Significant DMR Count"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_multiformat_ggplot(p_chr_hmc, file.path(OUTPUT_DIR, SECTION_DIR, "43d_cg_hmc_chromosome_distribution"),
                        width = 12, height = 7)
cat("  Saved 43d_cg_hmc_chromosome_distribution\n")

# =============================================================================
# FIGURE 43e: mC vs hmC CHROMOSOME COMPARISON
# =============================================================================

cat("\n--- Figure 43e: mC vs hmC Chromosome Comparison ---\n")

chr_compare <- mc_chr_dist %>%
  dplyr::select(chr, mc_sig_count = sig_count, mc_total = total_genes) %>%
  dplyr::left_join(
    hmc_chr_dist %>%
      dplyr::select(chr, hmc_sig_count = sig_count, hmc_total = total_genes),
    by = "chr"
  ) %>%
  dplyr::mutate(
    mc_pct  = mc_sig_count / sum(mc_sig_count) * 100,
    hmc_pct = hmc_sig_count / sum(hmc_sig_count) * 100
  ) %>%
  tidyr::pivot_longer(cols = c(mc_pct, hmc_pct),
                      names_to = "context", values_to = "pct") %>%
  dplyr::mutate(context = ifelse(context == "mc_pct", "CG mC", "CG hmC"))

p_chr_compare <- ggplot(chr_compare, aes(x = chr, y = pct, fill = context)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("CG mC" = "#E41A1C", "CG hmC" = "#377EB8"), name = "Context") +
  labs(
    title = "Chromosome Distribution: CG mC vs hmC Significant DMRs",
    subtitle = "Normalized to percentage of each context's significant DMRs",
    x = "Chromosome", y = "% of Significant DMRs"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

save_multiformat_ggplot(p_chr_compare, file.path(OUTPUT_DIR, SECTION_DIR, "43e_cg_mc_vs_hmc_chr_comparison"),
                        width = 12, height = 7)
cat("  Saved 43e_cg_mc_vs_hmc_chr_comparison\n")

# =============================================================================
# FIGURE 43f: OBSERVED/EXPECTED RATIO HEATMAP
# =============================================================================

cat("\n--- Figure 43f: Observed/Expected Ratio Heatmap ---\n")

oe_matrix <- mc_chr_dist %>%
  dplyr::select(chr, mC = obs_exp_ratio) %>%
  dplyr::left_join(
    hmc_chr_dist %>% dplyr::select(chr, hmC = obs_exp_ratio),
    by = "chr"
  ) %>%
  tibble::column_to_rownames("chr") %>%
  as.matrix()

# Cap extreme values for visualization
oe_capped <- pmin(pmax(oe_matrix, 0.25), 2.5)

# Annotation for significance
sig_annot <- data.frame(
  mC_sig = mc_chr_dist$enrichment,
  hmC_sig = hmc_chr_dist$enrichment,
  row.names = as.character(mc_chr_dist$chr)
)

p_oe_heatmap <- pheatmap(
  oe_capped,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(0.25, 2.5, length.out = 101),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = round(oe_matrix, 2),
  number_format = "%.2f",
  fontsize_number = 9,
  main = "Chromosome Enrichment: Observed/Expected Ratio\nCG mC and hmC Significant DMRs",
  angle_col = 0,
  silent = TRUE
)

pdf(file.path(OUTPUT_DIR, SECTION_DIR, "43f_cg_chr_obs_exp_heatmap.pdf"), width = 6, height = 10)
grid::grid.newpage()
grid::grid.draw(p_oe_heatmap$gtable)
dev.off()

png(file.path(OUTPUT_DIR, SECTION_DIR, "43f_cg_chr_obs_exp_heatmap.png"),
    width = 6, height = 10, units = "in", res = 300)
grid::grid.newpage()
grid::grid.draw(p_oe_heatmap$gtable)
dev.off()

cat("  Saved 43f_cg_chr_obs_exp_heatmap\n")

# =============================================================================
# FIGURE 43g: CHROMOSOME DIRECTION SPLIT (mC)
# =============================================================================

cat("\n--- Figure 43g: Chromosome Direction Split (mC) ---\n")

mc_chr_direction <- mc_sig %>%
  dplyr::count(chr, direction, name = "count") %>%
  dplyr::mutate(
    chr_label = gsub("chr", "", chr),
    chr_num = suppressWarnings(as.numeric(chr_label)),
    chr_order = ifelse(is.na(chr_num), 100 + match(chr_label, c("X", "Y", "M")), chr_num)
  ) %>%
  dplyr::arrange(chr_order) %>%
  dplyr::mutate(chr = factor(chr, levels = unique(chr)))

# Compute hyper fraction per chromosome for annotation
hyper_frac <- mc_sig %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(
    n = dplyr::n(),
    n_hyper = sum(mod_difference > 0),
    hyper_pct = n_hyper / n * 100,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    chr_label = gsub("chr", "", chr),
    chr_num = suppressWarnings(as.numeric(chr_label)),
    chr_order = ifelse(is.na(chr_num), 100 + match(chr_label, c("X", "Y", "M")), chr_num)
  ) %>%
  dplyr::arrange(chr_order) %>%
  dplyr::mutate(chr = factor(chr, levels = unique(chr)))

global_hyper_pct <- 100 * n_mc_hyper / nrow(mc_sig)

p_chr_direction <- ggplot(mc_chr_direction, aes(x = chr, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = COLORS$direction, name = "Direction") +
  labs(
    title = "CG mC DMRs: Direction by Chromosome",
    subtitle = sprintf("Global: %.1f%% hypermethylated | Stacked by direction per chromosome",
                       global_hyper_pct),
    x = "Chromosome", y = "Significant DMR Count"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

save_multiformat_ggplot(p_chr_direction, file.path(OUTPUT_DIR, SECTION_DIR, "43g_cg_mc_chr_direction"),
                        width = 12, height = 7)
cat("  Saved 43g_cg_mc_chr_direction\n")

# =============================================================================
# FIGURE 43h: CHROMOSOME DIRECTION SPLIT (hmC)
# =============================================================================

cat("\n--- Figure 43h: Chromosome Direction Split (hmC) ---\n")

hmc_chr_direction <- hmc_sig %>%
  dplyr::count(chr, direction, name = "count") %>%
  dplyr::mutate(
    chr_label = gsub("chr", "", chr),
    chr_num = suppressWarnings(as.numeric(chr_label)),
    chr_order = ifelse(is.na(chr_num), 100 + match(chr_label, c("X", "Y", "M")), chr_num)
  ) %>%
  dplyr::arrange(chr_order) %>%
  dplyr::mutate(chr = factor(chr, levels = unique(chr)))

global_hmc_hyper_pct <- 100 * n_hmc_hyper / nrow(hmc_sig)

p_hmc_chr_direction <- ggplot(hmc_chr_direction, aes(x = chr, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = COLORS$direction, name = "Direction") +
  labs(
    title = "CG hmC DMRs: Direction by Chromosome",
    subtitle = sprintf("Global: %.1f%% hypermethylated | Stacked by direction per chromosome",
                       global_hmc_hyper_pct),
    x = "Chromosome", y = "Significant DMR Count"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

save_multiformat_ggplot(p_hmc_chr_direction, file.path(OUTPUT_DIR, SECTION_DIR, "43h_cg_hmc_chr_direction"),
                        width = 12, height = 7)
cat("  Saved 43h_cg_hmc_chr_direction\n")

# =============================================================================
# FIGURE 43i: REGIONAL CHROMOSOME ENRICHMENT
# =============================================================================

cat("\n--- Figure 43i: Regional Chromosome Enrichment ---\n")

# Compute chromosome distribution for each genomic region's mC DMRs
region_chr_results <- list()

for (region_name in names(region_dmrs)) {
  rdf <- region_dmrs[[region_name]]
  if (is.null(rdf) || sum(rdf$significant) < 10) next

  sig_df <- rdf %>%
    dplyr::filter(significant) %>%
    dplyr::distinct(gene, .keep_all = TRUE)

  chr_dist_r <- compute_chr_distribution(sig_df, rdf, region_name)

  enriched <- chr_dist_r %>%
    dplyr::filter(enrichment == "Enriched") %>%
    dplyr::select(chr, obs_exp_ratio, fisher_p)

  if (nrow(enriched) > 0) {
    enriched$region <- region_name
    region_chr_results[[region_name]] <- enriched
  }
}

if (length(region_chr_results) > 0) {
  region_enrichment <- dplyr::bind_rows(region_chr_results) %>%
    dplyr::mutate(
      neg_log10_p = -log10(fisher_p),
      label = sprintf("O/E=%.2f", obs_exp_ratio)
    )

  p_regional_chr <- ggplot(region_enrichment,
                           aes(x = chr, y = region, fill = obs_exp_ratio, size = neg_log10_p)) +
    geom_point(shape = 21, color = "black") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 1, name = "O/E Ratio") +
    scale_size_continuous(range = c(3, 10), name = expression(-log[10](p))) +
    labs(
      title = "Chromosome Enrichment Across Genomic Regions",
      subtitle = "Only showing significantly enriched chromosomes (Fisher p < 0.05)",
      x = "Chromosome", y = "Genomic Region"
    ) +
    theme_biomodal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_multiformat_ggplot(p_regional_chr,
                          file.path(OUTPUT_DIR, SECTION_DIR, "43i_cg_regional_chr_enrichment"),
                          width = 12, height = 7)
  cat("  Saved 43i_cg_regional_chr_enrichment\n")
} else {
  cat("  No significantly enriched chromosomes across regions (skipping)\n")
}

# =============================================================================
# FIGURE 43j: TOP DMR GENES BAR CHART (Top 50 by q-value)
# =============================================================================

cat("\n--- Figure 43j: Top 50 mC DMR Genes ---\n")

top_n <- 50
mc_top <- mc_sig %>%
  dplyr::arrange(dmr_qvalue) %>%
  head(top_n) %>%
  dplyr::mutate(
    in_hmc = gene %in% hmc_genes,
    hmc_status = ifelse(in_hmc, "Also hmC sig", "mC only"),
    gene = factor(gene, levels = rev(gene))
  )

p_top_genes <- ggplot(mc_top, aes(x = gene, y = mod_difference * 100, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, color = "black") +
  geom_point(aes(shape = hmc_status),
             size = 2, color = "black", show.legend = TRUE,
             position = position_nudge(y = ifelse(mc_top$mod_difference > 0, 0.2, -0.2))) +
  scale_fill_manual(values = COLORS$direction, name = "Direction") +
  scale_shape_manual(values = c("Also hmC sig" = 16, "mC only" = 4),
                     name = "hmC Status") +
  coord_flip() +
  labs(
    title = sprintf("Top %d Significant CG mC DMR Genes (by q-value)", top_n),
    subtitle = "BAP1-KO vs Control | Gene body methylation difference",
    x = "", y = "mC Difference (%)"
  ) +
  theme_biomodal() +
  theme(axis.text.y = element_text(size = 8))

save_multiformat_ggplot(p_top_genes, file.path(OUTPUT_DIR, SECTION_DIR, "43j_cg_top_dmr_genes"),
                        width = 10, height = 14)
cat("  Saved 43j_cg_top_dmr_genes\n")

# =============================================================================
# FIGURE 43k: SIGNIFICANCE vs GENE DENSITY (mC)
# =============================================================================

cat("\n--- Figure 43k: Significance Rate vs Gene Density ---\n")

# Compute per-chromosome significance rate and gene density
chr_rates <- mc_chr_dist %>%
  dplyr::select(chr, total_genes, sig_count, sig_rate, obs_exp_ratio, fisher_p, enrichment)

# Correlation between gene count and significance rate
cor_test <- cor.test(chr_rates$total_genes, chr_rates$sig_rate, method = "spearman")
cat(sprintf("  Gene count vs sig rate: rho=%.3f, p=%.3f\n",
            cor_test$estimate, cor_test$p.value))

p_density <- ggplot(chr_rates, aes(x = total_genes, y = sig_rate * 100)) +
  geom_point(aes(color = enrichment), size = 4) +
  geom_text_repel(aes(label = chr), size = 3.5, max.overlaps = 25) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed") +
  scale_color_manual(
    values = c("Enriched" = "#D7191C", "Depleted" = "#2C7BB6", "Not significant" = "grey50"),
    name = "Enrichment"
  ) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("Spearman rho = %.3f\np = %.3f",
                           cor_test$estimate, cor_test$p.value),
           hjust = 1.1, vjust = 1.5, size = 4.5, fontface = "italic") +
  labs(
    title = "CG mC: Significance Rate vs Total Genes per Chromosome",
    subtitle = "Are enriched chromosomes simply gene-dense?",
    x = "Total Genes Tested", y = "Significance Rate (%)"
  ) +
  theme_biomodal()

save_multiformat_ggplot(p_density, file.path(OUTPUT_DIR, SECTION_DIR, "43k_cg_mc_sig_rate_vs_gene_density"),
                        width = 10, height = 8)
cat("  Saved 43k_cg_mc_sig_rate_vs_gene_density\n")

# =============================================================================
# FIGURE 43l: EFFECT SIZE BY CHROMOSOME (mC)
# =============================================================================

cat("\n--- Figure 43l: Effect Size by Chromosome (mC) ---\n")

mc_sig_chr <- mc_sig %>%
  dplyr::mutate(
    chr_label = gsub("chr", "", chr),
    chr_num = suppressWarnings(as.numeric(chr_label)),
    chr_order = ifelse(is.na(chr_num), 100 + match(chr_label, c("X", "Y", "M")), chr_num)
  ) %>%
  dplyr::arrange(chr_order) %>%
  dplyr::mutate(chr = factor(chr, levels = unique(chr)))

# Kruskal-Wallis test for effect size differences across chromosomes
kw_test <- kruskal.test(mod_difference ~ chr, data = mc_sig_chr)
cat(sprintf("  Kruskal-Wallis test (effect size ~ chromosome): chi2=%.2f, p=%.4f\n",
            kw_test$statistic, kw_test$p.value))

p_effect_chr <- ggplot(mc_sig_chr, aes(x = chr, y = mod_difference * 100)) +
  geom_boxplot(aes(fill = chr), outlier.shape = NA, alpha = 0.7, show.legend = FALSE) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.8, color = "grey30") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("Kruskal-Wallis p = %.4f", kw_test$p.value),
           hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
  scale_fill_manual(values = rep(c("grey80", "grey90"), length.out = nlevels(mc_sig_chr$chr))) +
  labs(
    title = "CG mC Effect Size Distribution by Chromosome",
    subtitle = sprintf("%d significant DMRs | Do certain chromosomes show stronger effects?",
                       nrow(mc_sig_chr)),
    x = "Chromosome", y = "mC Difference (%)"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_multiformat_ggplot(p_effect_chr, file.path(OUTPUT_DIR, SECTION_DIR, "43l_cg_mc_effect_size_by_chr"),
                        width = 14, height = 7)
cat("  Saved 43l_cg_mc_effect_size_by_chr\n")

# =============================================================================
# TABLE EXPORTS
# =============================================================================

cat("\n--- Table Exports ---\n")

# Chromosome distribution summary table (mC)
mc_chr_export <- mc_chr_dist %>%
  dplyr::select(chr, total_genes, sig_count, sig_rate, expected, obs_exp_ratio,
                fisher_p, enrichment) %>%
  dplyr::mutate(
    sig_rate = round(sig_rate, 4),
    expected = round(expected, 1),
    obs_exp_ratio = round(obs_exp_ratio, 3),
    fisher_p = signif(fisher_p, 3)
  )

write.table(mc_chr_export, file.path(TABLES_DIR, "cg_mc_chromosome_distribution.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved cg_mc_chromosome_distribution.tsv (%d chromosomes)\n", nrow(mc_chr_export)))

# Chromosome distribution summary table (hmC)
hmc_chr_export <- hmc_chr_dist %>%
  dplyr::select(chr, total_genes, sig_count, sig_rate, expected, obs_exp_ratio,
                fisher_p, enrichment) %>%
  dplyr::mutate(
    sig_rate = round(sig_rate, 4),
    expected = round(expected, 1),
    obs_exp_ratio = round(obs_exp_ratio, 3),
    fisher_p = signif(fisher_p, 3)
  )

write.table(hmc_chr_export, file.path(TABLES_DIR, "cg_hmc_chromosome_distribution.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved cg_hmc_chromosome_distribution.tsv (%d chromosomes)\n", nrow(hmc_chr_export)))

# Top 50 mC genes with hmC cross-reference
hmc_for_join <- hmc_dmr %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(gene,
                hmc_difference = mod_difference,
                hmc_qvalue = dmr_qvalue,
                hmc_significant = significant,
                hmc_direction = direction)

top_genes_export <- mc_sig %>%
  dplyr::arrange(dmr_qvalue) %>%
  head(top_n) %>%
  dplyr::select(gene, chr, start, end, num_contexts, mean_coverage,
                mc_group1 = mean_mod_group1, mc_group2 = mean_mod_group2,
                mc_difference = mod_difference, mc_qvalue = dmr_qvalue,
                mc_direction = direction) %>%
  dplyr::left_join(hmc_for_join, by = "gene")

write.table(top_genes_export, file.path(TABLES_DIR, "cg_top50_mc_dmr_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved cg_top50_mc_dmr_genes.tsv (%d genes)\n", nrow(top_genes_export)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("CG EXPLORATORY ANALYSIS SUMMARY\n")
cat("================================================================================\n")
cat(sprintf("CG mC significant DMRs:   %d (q < %.2f)\n", nrow(mc_sig), Q_THRESHOLD))
cat(sprintf("  Hypermethylated:         %d (%.1f%%)\n", n_mc_hyper, 100 * n_mc_hyper / nrow(mc_sig)))
cat(sprintf("  Hypomethylated:          %d (%.1f%%)\n", n_mc_hypo, 100 * n_mc_hypo / nrow(mc_sig)))
cat(sprintf("CG hmC significant DMRs:   %d\n", nrow(hmc_sig)))
cat(sprintf("  Hypermethylated:         %d (%.1f%%)\n", n_hmc_hyper, 100 * n_hmc_hyper / nrow(hmc_sig)))
cat(sprintf("  Hypomethylated:          %d (%.1f%%)\n", n_hmc_hypo, 100 * n_hmc_hypo / nrow(hmc_sig)))
cat(sprintf("mC/hmC overlap:            %d genes\n", n_both))
cat(sprintf("mC-only:                   %d genes\n", n_mc_only))
cat(sprintf("hmC-only:                  %d genes\n", n_hmc_only))
cat("\nChromosome enrichment (mC, Fisher p < 0.05):\n")
if (nrow(mc_enriched) > 0) {
  for (i in seq_len(nrow(mc_enriched))) {
    cat(sprintf("  %s: O/E=%.2f, p=%.2e\n",
                mc_enriched$chr[i], mc_enriched$obs_exp_ratio[i], mc_enriched$fisher_p[i]))
  }
} else {
  cat("  No significantly enriched chromosomes\n")
}
if (nrow(mc_depleted) > 0) {
  cat("Chromosome depletion (mC, Fisher p < 0.05):\n")
  for (i in seq_len(nrow(mc_depleted))) {
    cat(sprintf("  %s: O/E=%.2f, p=%.2e\n",
                mc_depleted$chr[i], mc_depleted$obs_exp_ratio[i], mc_depleted$fisher_p[i]))
  }
}
cat(sprintf("\nGene density vs sig rate:  rho=%.3f, p=%.3f\n",
            cor_test$estimate, cor_test$p.value))
cat(sprintf("Effect size ~ chromosome:  Kruskal-Wallis p=%.4f\n", kw_test$p.value))
cat("================================================================================\n")

cat("Section 43 complete.\n\n")
