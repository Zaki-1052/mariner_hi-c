# biomodal/downstream/scripts/viz_sections/section_65_mecp2_methylation_scale.R
# Section 65: MeCP2 Binding Context — Methylation Scale vs MeCP2 Occupancy
# Standalone script - sources shared config for all dependencies and data
#
# Shows that methylation != MeCP2 binding: far more genes have significant 5mC
# changes than significant MeCP2 binding changes. MeCP2 is context-dependent.
#
#   Panel 65a: Cascade hierarchy (methylation scale vs MeCP2 scale)
#   Panel 65b: Gene-level Venn (mC hyper, MeCP2 sig up, coordinated)
#   Panel 65c: Quadrant scatter (all tested genes: mC change vs MeCP2 fold)
#   Panel 65d: Proportional summary bar (methylated vs MeCP2-bound)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_65_mecp2_methylation_scale.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# LOCAL HELPER FUNCTIONS (copied from section_11 for independence)
# =============================================================================

load_mecp2_annotated <- function(filepath) {
  if (!file.exists(filepath)) {
    warning("MeCP2 annotated file not found: ", filepath)
    return(NULL)
  }
  df <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                   fill = TRUE, quote = "")
  numeric_cols <- c("start", "end", "width", "Conc", "Conc_mut", "Conc_ctrl",
                    "Fold", "p.value", "FDR", "geneStart", "geneEnd",
                    "geneLength", "distanceToTSS")
  for (col in numeric_cols) {
    if (col %in% colnames(df)) df[[col]] <- as.numeric(df[[col]])
  }
  cat(sprintf("  Loaded %d MeCP2 peaks\n", nrow(df)))
  cat(sprintf("  Significant (FDR < 0.05): %d up, %d down\n",
              sum(df$FDR < 0.05 & df$Fold > 0, na.rm = TRUE),
              sum(df$FDR < 0.05 & df$Fold < 0, na.rm = TRUE)))
  return(df)
}

aggregate_mecp2_by_gene <- function(mecp2_df) {
  mecp2_df %>%
    dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    group_by(SYMBOL) %>%
    summarise(
      mean_fold = mean(Fold, na.rm = TRUE),
      max_fold = max(Fold, na.rm = TRUE),
      min_fold = min(Fold, na.rm = TRUE),
      nearest_fold = Fold[which.min(abs(distanceToTSS))],
      nearest_fdr = FDR[which.min(abs(distanceToTSS))],
      min_distance_tss = min(abs(distanceToTSS), na.rm = TRUE),
      n_peaks = n(),
      any_sig = any(FDR < 0.05),
      any_sig_up = any(FDR < 0.05 & Fold > 0),
      any_sig_down = any(FDR < 0.05 & Fold < 0),
      .groups = "drop"
    )
}

# =============================================================================
# SECTION 65: MeCP2 BINDING CONTEXT
# =============================================================================

cat("================================================================================\n")
cat("SECTION 65: MeCP2 BINDING CONTEXT — METHYLATION SCALE vs MeCP2 OCCUPANCY\n")
cat("================================================================================\n\n")

# ---- Load MeCP2 data --------------------------------------------------------

cat("Loading MeCP2 differential binding data...\n")
mecp2_annotated <- load_mecp2_annotated(MECP2_FILES$annotated)

if (is.null(mecp2_annotated)) {
  stop("MeCP2 annotated file not found. Cannot proceed.")
}

cat("\nAggregating MeCP2 peaks to gene level...\n")
mecp2_gene <- aggregate_mecp2_by_gene(mecp2_annotated)
cat(sprintf("  %d unique genes with MeCP2 peaks\n", nrow(mecp2_gene)))

# ---- Compute gene-level sets ------------------------------------------------

cat("\nComputing gene-level sets...\n")

all_tested_genes <- unique(mc_dmr$gene)
mc_sig_genes     <- mc_dmr %>% dplyr::filter(significant) %>% dplyr::pull(gene)
mc_hyper_genes   <- mc_dmr %>% dplyr::filter(significant & mod_difference > 0) %>% dplyr::pull(gene)
mc_hypo_genes    <- mc_dmr %>% dplyr::filter(significant & mod_difference < 0) %>% dplyr::pull(gene)

mecp2_all_genes    <- mecp2_gene$SYMBOL
mecp2_sig_genes    <- mecp2_gene %>% dplyr::filter(any_sig) %>% dplyr::pull(SYMBOL)
mecp2_sig_up_genes <- mecp2_gene %>% dplyr::filter(any_sig_up) %>% dplyr::pull(SYMBOL)

# Coordinated genes (mC up + hmC down, both sig)
mc_sig_coord <- mc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::select(gene, mc_diff = mod_difference)
hmc_sig_coord <- hmc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::select(gene, hmc_diff = mod_difference)
coord_genes <- inner_join(mc_sig_coord, hmc_sig_coord, by = "gene") %>%
  dplyr::filter(mc_diff > 0 & hmc_diff < 0) %>%
  dplyr::pull(gene)

# Intersections
both_mc_hyper_and_mecp2_up <- intersect(mc_hyper_genes, mecp2_sig_up_genes)
mc_hyper_with_any_mecp2    <- intersect(mc_hyper_genes, mecp2_all_genes)

cat(sprintf("  Gene bodies tested:           %d\n", length(all_tested_genes)))
cat(sprintf("  Significant 5mC change:       %d (%.1f%%)\n",
            length(mc_sig_genes), 100 * length(mc_sig_genes) / length(all_tested_genes)))
cat(sprintf("  5mC Hypermethylated:          %d\n", length(mc_hyper_genes)))
cat(sprintf("  5mC Hypomethylated:           %d\n", length(mc_hypo_genes)))
cat(sprintf("  Coordinated (mC up/hmC down): %d\n", length(coord_genes)))
cat(sprintf("  MeCP2 peaks (any gene):       %d\n", length(mecp2_all_genes)))
cat(sprintf("  MeCP2 sig change (any gene):  %d\n", length(mecp2_sig_genes)))
cat(sprintf("  MeCP2 sig UP:                 %d\n", length(mecp2_sig_up_genes)))
cat(sprintf("  mC hyper + MeCP2 sig up:      %d (%.1f%% of mC hyper)\n",
            length(both_mc_hyper_and_mecp2_up),
            100 * length(both_mc_hyper_and_mecp2_up) / length(mc_hyper_genes)))

# ---- Panel 65a: Cascade hierarchy -------------------------------------------

cat("\nCreating Figure 65a: cascade hierarchy...\n")

# Left panel: methylation hierarchy
meth_cascade <- data.frame(
  category = c("Gene Bodies\nTested", "Sig 5mC\n(q<0.05)", "mC Hyper\n(mut > ctrl)"),
  n = c(length(all_tested_genes), length(mc_sig_genes), length(mc_hyper_genes)),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    category = factor(category, levels = rev(category)),
    fill_group = "Methylation"
  )

# Right panel: MeCP2 hierarchy
mecp2_cascade <- data.frame(
  category = c("Genes with\nMeCP2 Peaks", "Sig MeCP2\n(FDR<0.05)", "MeCP2 Sig Up\n(mut > ctrl)"),
  n = c(length(mecp2_all_genes), length(mecp2_sig_genes), length(mecp2_sig_up_genes)),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    category = factor(category, levels = rev(category)),
    fill_group = "MeCP2"
  )

meth_fill <- c("#FCBBA1", "#FB6A4A", "#CB181D")
mecp2_fill <- c("#BCBDDC", "#807DBA", "#54278F")

p_65a_left <- ggplot(meth_cascade, aes(x = category, y = n)) +
  geom_bar(stat = "identity", width = 0.7, fill = rev(meth_fill)) +
  geom_text(aes(label = sprintf("%s\n(%.1f%%)", format(n, big.mark = ","),
                                100 * n / length(all_tested_genes))),
            hjust = -0.05, size = 3.5) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(title = "Methylation (5mC) Hierarchy",
       x = "", y = "Number of Genes") +
  theme_biomodal() +
  theme(plot.title = element_text(color = "#CB181D"))

p_65a_right <- ggplot(mecp2_cascade, aes(x = category, y = n)) +
  geom_bar(stat = "identity", width = 0.7, fill = rev(mecp2_fill)) +
  geom_text(aes(label = sprintf("%s\n(%.1f%%)", format(n, big.mark = ","),
                                100 * n / length(mecp2_all_genes))),
            hjust = -0.05, size = 3.5) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(title = "MeCP2 Binding Hierarchy",
       x = "", y = "Number of Genes") +
  theme_biomodal() +
  theme(plot.title = element_text(color = "#54278F"))

# Intersection annotation
intersect_label <- sprintf("Overlap: %d genes have both\nsig mC hyper AND sig MeCP2 up\n(%.1f%% of mC hyper genes)",
                           length(both_mc_hyper_and_mecp2_up),
                           100 * length(both_mc_hyper_and_mecp2_up) / length(mc_hyper_genes))

p_65a <- (p_65a_left | p_65a_right) +
  plot_annotation(
    title = "Scale of Methylation Changes vs MeCP2 Binding Changes",
    subtitle = intersect_label,
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30")
    )
  )

save_multiformat_ggplot(p_65a,
                        file.path(OUTPUT_DIR, "65a_cascade_hierarchy"),
                        width = 14, height = 7)

# ---- Panel 65b: Venn diagram ------------------------------------------------

cat("Creating Figure 65b: gene-level Venn diagram...\n")

venn_list <- list(
  "mC Hyper (q<0.05)" = mc_hyper_genes,
  "MeCP2 Sig Up (FDR<0.05)" = mecp2_sig_up_genes,
  "Coordinated (mC↑/hmC↓)" = coord_genes
)

p_65b <- ggVennDiagram(venn_list, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#E41A1C", name = "Count") +
  scale_color_manual(values = c("#CB181D", "#54278F", "#377EB8")) +
  labs(
    title = "Gene-Level Overlap: Methylation vs MeCP2 Binding",
    subtitle = "Most mC hyper genes do NOT show significant MeCP2 binding increase"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

save_multiformat_ggplot(p_65b,
                        file.path(OUTPUT_DIR, "65b_gene_level_venn"),
                        width = 10, height = 9)

# ---- Panel 65c: Quadrant scatter (all tested genes) -------------------------

cat("Creating Figure 65c: quadrant scatter (all tested genes)...\n")

# Join all tested genes with MeCP2 gene-level data
all_gene_join <- mc_dmr %>%
  left_join(
    mecp2_gene %>% dplyr::select(SYMBOL, nearest_fold, nearest_fdr, any_sig, any_sig_up),
    by = c("gene" = "SYMBOL")
  ) %>%
  dplyr::filter(!is.na(nearest_fold))

cat(sprintf("  %d genes with both mC and MeCP2 data\n", nrow(all_gene_join)))

# Classify into quadrants based on significance
all_gene_join <- all_gene_join %>%
  dplyr::mutate(
    mc_sig_up = significant & mod_difference > 0,
    mecp2_sig_up_flag = any_sig_up,
    quadrant = case_when(
      mc_sig_up & mecp2_sig_up_flag  ~ "mC↑ + MeCP2↑\n(Context match)",
      mc_sig_up & !mecp2_sig_up_flag ~ "mC↑, MeCP2 unchanged\n(Methylation w/o MeCP2)",
      !mc_sig_up & mecp2_sig_up_flag ~ "MeCP2↑, mC unchanged\n(MeCP2 w/o methylation)",
      TRUE                           ~ "Neither significant"
    ),
    quadrant = factor(quadrant, levels = c(
      "mC↑ + MeCP2↑\n(Context match)",
      "mC↑, MeCP2 unchanged\n(Methylation w/o MeCP2)",
      "MeCP2↑, mC unchanged\n(MeCP2 w/o methylation)",
      "Neither significant"
    ))
  )

# Count per quadrant
q_counts <- all_gene_join %>%
  count(quadrant) %>%
  dplyr::mutate(pct = sprintf("%.1f%%", 100 * n / nrow(all_gene_join)))

cat("  Quadrant counts:\n")
for (i in seq_len(nrow(q_counts))) {
  cat(sprintf("    %s: %d (%s)\n",
              gsub("\n", " ", q_counts$quadrant[i]),
              q_counts$n[i], q_counts$pct[i]))
}

quad_colors <- c(
  "mC↑ + MeCP2↑\n(Context match)"              = "#FF7F00",
  "mC↑, MeCP2 unchanged\n(Methylation w/o MeCP2)" = "#E41A1C",
  "MeCP2↑, mC unchanged\n(MeCP2 w/o methylation)" = "#7570B3",
  "Neither significant"                          = "grey75"
)

# Compute quadrant annotation positions
x_range <- range(all_gene_join$mod_difference * 100, na.rm = TRUE)
y_range <- range(all_gene_join$nearest_fold, na.rm = TRUE)

q_annot <- data.frame(
  quadrant = levels(all_gene_join$quadrant),
  x = c(x_range[2] * 0.7, x_range[2] * 0.7, x_range[1] * 0.7, x_range[1] * 0.7),
  y = c(y_range[2] * 0.85, y_range[1] * 0.85, y_range[2] * 0.85, y_range[1] * 0.85),
  stringsAsFactors = FALSE
) %>%
  left_join(q_counts, by = "quadrant")

p_65c <- ggplot(all_gene_join, aes(x = mod_difference * 100, y = nearest_fold)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_point(aes(color = quadrant), alpha = 0.4, size = 1.2) +
  geom_text(data = q_annot,
            aes(x = x, y = y, label = sprintf("n=%s\n(%s)", format(n, big.mark = ","), pct)),
            size = 3.5, fontface = "bold", color = "black") +
  scale_color_manual(values = quad_colors, name = "Category") +
  labs(
    title = "5mC Change vs MeCP2 Binding Change (All Tested Genes)",
    subtitle = "Most genes with increased methylation show NO significant MeCP2 binding change",
    x = "5mC Change (Mut - Ctrl, %)", y = "MeCP2 Fold Change (log2)"
  ) +
  theme_biomodal() +
  theme(legend.position = "right")

save_multiformat_ggplot(p_65c,
                        file.path(OUTPUT_DIR, "65c_allgene_quadrant_scatter"),
                        width = 12, height = 10)

# ---- Panel 65d: Proportional summary bar ------------------------------------

cat("Creating Figure 65d: proportional summary bar...\n")

prop_data <- data.frame(
  category = c(
    "Gene Bodies\nTested",
    "Sig 5mC Change\n(q<0.05)",
    "5mC Hyper\n(mut > ctrl)",
    "Hyper + Any\nMeCP2 Peak",
    "Hyper + Sig\nMeCP2 Up"
  ),
  n = c(
    length(all_tested_genes),
    length(mc_sig_genes),
    length(mc_hyper_genes),
    length(mc_hyper_with_any_mecp2),
    length(both_mc_hyper_and_mecp2_up)
  ),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    category = factor(category, levels = category),
    pct_of_tested = 100 * n / length(all_tested_genes),
    fill_val = seq_len(n())
  )

prop_fills <- c("#FEE5D9", "#FCAE91", "#FB6A4A", "#DE2D26", "#A50F15")

p_65d <- ggplot(prop_data, aes(x = category, y = n)) +
  geom_bar(stat = "identity", width = 0.7, fill = prop_fills) +
  geom_text(aes(label = sprintf("%s\n(%.1f%%)", format(n, big.mark = ","), pct_of_tested)),
            vjust = -0.3, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Methylation-to-MeCP2 Binding Funnel",
    subtitle = "Most methylated genes are NOT bound by MeCP2 — context-dependent recruitment",
    x = "", y = "Number of Genes"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(size = 9))

save_multiformat_ggplot(p_65d,
                        file.path(OUTPUT_DIR, "65d_proportional_funnel"),
                        width = 10, height = 8)

# ---- Combined figure ---------------------------------------------------------

cat("Creating combined Figure 65...\n")

p_65_combined <- ((p_65a_left | p_65a_right) / (p_65b | p_65d)) +
  plot_annotation(
    title = "MeCP2 is a Context-Dependent Methylation Reader",
    subtitle = "Far more genes show methylation changes than MeCP2 binding changes — methylation is necessary but not sufficient",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11)
    )
  )

save_multiformat_ggplot(p_65_combined,
                        file.path(OUTPUT_DIR, "65_mecp2_methylation_scale"),
                        width = 16, height = 14)

# ---- Statistics & tables -----------------------------------------------------

cat("\nRunning Fisher's exact test for context-dependence...\n")

# 2x2: sig mC (yes/no) x sig MeCP2 up (yes/no) at gene level
genes_in_both <- intersect(all_tested_genes, mecp2_all_genes)
fisher_df <- data.frame(
  gene = genes_in_both,
  mc_sig = genes_in_both %in% mc_sig_genes,
  mecp2_sig_up = genes_in_both %in% mecp2_sig_up_genes,
  stringsAsFactors = FALSE
)

fisher_mat <- table(
  `5mC Sig` = fisher_df$mc_sig,
  `MeCP2 Sig Up` = fisher_df$mecp2_sig_up
)
fisher_result <- fisher.test(fisher_mat)

cat(sprintf("  Fisher's exact test (sig mC x sig MeCP2 up):\n"))
cat(sprintf("    OR = %.3f, p = %.2e\n", fisher_result$estimate, fisher_result$p.value))
cat(sprintf("    Contingency table:\n"))
print(fisher_mat)

# Save cascade table
cat("\nSaving tables...\n")

cascade_table <- data.frame(
  category = c("Gene bodies tested",
               "Significant 5mC change (q<0.05)",
               "5mC Hypermethylated (mut > ctrl)",
               "5mC Hypomethylated (mut < ctrl)",
               "Coordinated (mC up / hmC down)",
               "Genes with any MeCP2 peak",
               "Genes with sig MeCP2 change (FDR<0.05)",
               "Genes with sig MeCP2 Up",
               "mC Hyper + any MeCP2 peak",
               "mC Hyper + sig MeCP2 Up"),
  n = c(length(all_tested_genes),
        length(mc_sig_genes),
        length(mc_hyper_genes),
        length(mc_hypo_genes),
        length(coord_genes),
        length(mecp2_all_genes),
        length(mecp2_sig_genes),
        length(mecp2_sig_up_genes),
        length(mc_hyper_with_any_mecp2),
        length(both_mc_hyper_and_mecp2_up)),
  pct_of_tested = c(100,
                    100 * length(mc_sig_genes) / length(all_tested_genes),
                    100 * length(mc_hyper_genes) / length(all_tested_genes),
                    100 * length(mc_hypo_genes) / length(all_tested_genes),
                    100 * length(coord_genes) / length(all_tested_genes),
                    100 * length(mecp2_all_genes) / length(all_tested_genes),
                    100 * length(mecp2_sig_genes) / length(all_tested_genes),
                    100 * length(mecp2_sig_up_genes) / length(all_tested_genes),
                    100 * length(mc_hyper_with_any_mecp2) / length(all_tested_genes),
                    100 * length(both_mc_hyper_and_mecp2_up) / length(all_tested_genes)),
  stringsAsFactors = FALSE
)

write.table(cascade_table,
            file.path(TABLES_DIR, "65_methylation_mecp2_cascade.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 65_methylation_mecp2_cascade.tsv\n")

# Per-gene overlap table
gene_overlap <- data.frame(
  gene = all_tested_genes,
  mc_significant = all_tested_genes %in% mc_sig_genes,
  mc_hyper = all_tested_genes %in% mc_hyper_genes,
  coordinated = all_tested_genes %in% coord_genes,
  mecp2_bound = all_tested_genes %in% mecp2_all_genes,
  mecp2_sig_change = all_tested_genes %in% mecp2_sig_genes,
  mecp2_sig_up = all_tested_genes %in% mecp2_sig_up_genes,
  stringsAsFactors = FALSE
)

write.table(gene_overlap,
            file.path(TABLES_DIR, "65_methylation_mecp2_gene_overlap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 65_methylation_mecp2_gene_overlap.tsv\n")

# Fisher test results
fisher_table <- data.frame(
  test = "sig_mC_x_sig_MeCP2_up",
  n_genes_tested = length(genes_in_both),
  both_sig = sum(fisher_df$mc_sig & fisher_df$mecp2_sig_up),
  mc_only = sum(fisher_df$mc_sig & !fisher_df$mecp2_sig_up),
  mecp2_only = sum(!fisher_df$mc_sig & fisher_df$mecp2_sig_up),
  neither = sum(!fisher_df$mc_sig & !fisher_df$mecp2_sig_up),
  fisher_OR = fisher_result$estimate,
  fisher_pvalue = fisher_result$p.value,
  stringsAsFactors = FALSE
)

write.table(fisher_table,
            file.path(TABLES_DIR, "65_fisher_context_dependence.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 65_fisher_context_dependence.tsv\n")

cat("\n")
cat("================================================================================\n")
cat("SECTION 65 COMPLETE\n")
cat("================================================================================\n")
