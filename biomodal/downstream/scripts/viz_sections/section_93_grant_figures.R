# biomodal/downstream/scripts/viz_sections/section_93_grant_figures.R
# Section 93: Grant Figure Composites — 4 figures for PI grant writing
#
# Generates 4 composite patchwork figures (3 panels each) + 12 individual panel SVGs.
# All data read from pre-loaded objects (mc_dmr, hmc_dmr) and pre-computed TSVs.
#
# Grant Figure 1: Bidirectional mC-up/hmC-down phenotype
# Grant Figure 2: Euchromatin silencing (K119ub → active chromatin hypermethylation)
# Grant Figure 3: Heterochromatin protection (Polycomb excluded from hypermethylation)
# Grant Figure 4: MeCP2 reads chromatin state, not methylation
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_93_grant_figures.R

source("scripts/viz_sections/_shared_config.R")

cat("================================================================================\n")
cat("SECTION 93: GRANT FIGURE COMPOSITES\n")
cat("================================================================================\n\n")

GRANT_DIR <- file.path(OUTPUT_DIR, "93_grant_figures")
dir.create(GRANT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# GRANT FIGURE 1: Bidirectional mC-up/hmC-down Phenotype
# =============================================================================

cat("── GRANT FIGURE 1: Bidirectional mC-up/hmC-down ──────────────────────────────────\n\n")

# ---- Panel A: Dual Volcano --------------------------------------------------

cat("  Panel A: Dual volcano plots...\n")

mc_volc <- mc_dmr %>%
  mutate(
    neg_log_q = -log10(pmax(dmr_qvalue, 1e-300)),
    color = ifelse(significant, as.character(direction), "Not Significant"),
    label = ifelse(gene %in% KEY_GENES & significant, gene, NA_character_),
    mod_type = "5mC"
  )

hmc_volc <- hmc_dmr %>%
  mutate(
    neg_log_q = -log10(pmax(dmr_qvalue, 1e-300)),
    color = case_when(
      !significant ~ "Not Significant",
      mod_difference > 0 ~ "Increased",
      TRUE ~ "Decreased"
    ),
    label = ifelse(gene %in% KEY_GENES & significant, gene, NA_character_),
    mod_type = "5hmC"
  )

n_mc_sig <- sum(mc_dmr$significant)
n_mc_hyper <- sum(mc_dmr$significant & mc_dmr$direction == "Hypermethylated")
n_hmc_sig <- sum(hmc_dmr$significant)
n_hmc_down <- sum(hmc_dmr$significant & hmc_dmr$mod_difference < 0)

volc_colors_mc <- c("Hypermethylated" = "#D7191C", "Hypomethylated" = "#2C7BB6",
                    "Not Significant" = "grey80")
volc_colors_hmc <- c("Increased" = "#D7191C", "Decreased" = "#2C7BB6",
                     "Not Significant" = "grey80")

p1a_mc <- ggplot(mc_volc, aes(x = mod_difference, y = neg_log_q, color = color)) +
  geom_point(size = 0.3, alpha = 0.4) +
  geom_hline(yintercept = -log10(Q_THRESHOLD), linetype = "dashed", color = "grey40") +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20, color = "black",
                  segment.color = "grey50", segment.size = 0.3, na.rm = TRUE) +
  scale_color_manual(values = volc_colors_mc, guide = "none") +
  labs(x = "5mC difference (mutant − control)",
       y = expression(-log[10](q)),
       subtitle = sprintf("5mC: %s sig (%d%% hyper)",
                          format(n_mc_sig, big.mark = ","),
                          round(100 * n_mc_hyper / n_mc_sig))) +
  theme_biomodal() +
  theme(plot.subtitle = element_text(size = 9))

p1a_hmc <- ggplot(hmc_volc, aes(x = mod_difference, y = neg_log_q, color = color)) +
  geom_point(size = 0.3, alpha = 0.4) +
  geom_hline(yintercept = -log10(Q_THRESHOLD), linetype = "dashed", color = "grey40") +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20, color = "black",
                  segment.color = "grey50", segment.size = 0.3, na.rm = TRUE) +
  scale_color_manual(values = volc_colors_hmc, guide = "none") +
  labs(x = "5hmC difference (mutant − control)",
       y = expression(-log[10](q)),
       subtitle = sprintf("5hmC: %s sig (%d%% decreased)",
                          format(n_hmc_sig, big.mark = ","),
                          round(100 * n_hmc_down / n_hmc_sig))) +
  theme_biomodal() +
  theme(plot.subtitle = element_text(size = 9))

p1a <- p1a_mc + p1a_hmc + plot_layout(nrow = 1)

save_multiformat_ggplot(p1a,
                        file.path(GRANT_DIR, "93a_dual_volcano"),
                        width = 12, height = 5)

# ---- Panel B: Quadrant Scatter -----------------------------------------------

cat("  Panel B: Quadrant scatter (coordinated mC-up/hmC-down)...\n")

coord <- read_tsv(file.path(TABLES_DIR, "coordinated_changes.tsv"),
                  show_col_types = FALSE)

coord <- coord %>%
  mutate(quadrant = case_when(
    mc_diff > 0 & hmc_diff < 0 ~ "Q4: mC up / hmC down",
    mc_diff < 0 & hmc_diff > 0 ~ "Q2: mC down / hmC up",
    mc_diff > 0 & hmc_diff > 0 ~ "Q3: same direction",
    mc_diff < 0 & hmc_diff < 0 ~ "Q1: same direction"
  ))

n_q4 <- sum(coord$quadrant == "Q4: mC up / hmC down", na.rm = TRUE)
n_total <- nrow(coord)
pct_coord <- round(100 * n_q4 / n_total, 1)

quad_colors <- c("Q4: mC up / hmC down" = "#D7191C",
                 "Q2: mC down / hmC up" = "#2C7BB6",
                 "Q1: same direction" = "grey60",
                 "Q3: same direction" = "grey60")

p1b <- ggplot(coord, aes(x = mc_diff, y = hmc_diff, color = quadrant)) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.3) +
  geom_point(size = 0.3, alpha = 0.4) +
  scale_color_manual(values = quad_colors) +
  annotate("text", x = max(coord$mc_diff) * 0.6, y = min(coord$hmc_diff) * 0.6,
           label = sprintf("Q4: %s genes\n(%s%%)", format(n_q4, big.mark = ","), pct_coord),
           size = 3.5, fontface = "bold", color = "#D7191C") +
  labs(x = "5mC difference", y = "5hmC difference",
       subtitle = "Co-significant genes: coordinated quadrant") +
  theme_biomodal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.subtitle = element_text(size = 9))

save_multiformat_ggplot(p1b,
                        file.path(GRANT_DIR, "93b_quadrant_scatter"),
                        width = 6, height = 6)

# ---- Panel C: Effect-Size Density -------------------------------------------

cat("  Panel C: Effect-size mirror density...\n")

mc_sig <- mc_dmr %>% filter(significant) %>%
  mutate(mod_type = "5mC", value = mod_difference)
hmc_sig <- hmc_dmr %>% filter(significant) %>%
  mutate(mod_type = "5hmC", value = mod_difference)
eff_data <- bind_rows(
  mc_sig %>% dplyr::select(mod_type, value),
  hmc_sig %>% dplyr::select(mod_type, value)
)

med_mc <- median(mc_sig$value)
med_hmc <- median(hmc_sig$value)

p1c <- ggplot(eff_data, aes(x = value, fill = mod_type)) +
  geom_density(alpha = 0.5, linewidth = 0.4) +
  geom_vline(xintercept = med_mc, linetype = "dashed", color = COLORS$methylation["5mC"],
             linewidth = 0.5) +
  geom_vline(xintercept = med_hmc, linetype = "dashed", color = COLORS$methylation["5hmC"],
             linewidth = 0.5) +
  scale_fill_manual(values = COLORS$methylation) +
  annotate("text", x = med_mc, y = Inf, vjust = 2,
           label = sprintf("median = %+.2f%%", med_mc), size = 2.8,
           color = COLORS$methylation["5mC"]) +
  annotate("text", x = med_hmc, y = Inf, vjust = 3.5,
           label = sprintf("median = %+.2f%%", med_hmc), size = 2.8,
           color = COLORS$methylation["5hmC"]) +
  labs(x = "Methylation difference (mutant − control, %)",
       y = "Density",
       subtitle = "Effect-size distributions (significant DMRs)") +
  theme_biomodal() +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_blank(),
        plot.subtitle = element_text(size = 9))

save_multiformat_ggplot(p1c,
                        file.path(GRANT_DIR, "93c_effect_size_density"),
                        width = 6, height = 5)

# ---- Figure 1 Composite -----------------------------------------------------

cat("  Compositing Grant Figure 1...\n")

fig1 <- (p1a_mc | p1a_hmc | p1b) / (p1c + plot_spacer() + plot_spacer()) +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 0.8))

save_multiformat_ggplot(fig1,
                        file.path(GRANT_DIR, "grant_fig1_bidirectional"),
                        width = 18, height = 10)

cat("  Done.\n\n")

# =============================================================================
# GRANT FIGURE 2: Euchromatin Silencing
# =============================================================================

cat("── GRANT FIGURE 2: Euchromatin Silencing ───────────────────────────────────\n\n")

# ---- Panel A: Chromatin State Direction Bar ----------------------------------

cat("  Panel A: Chromatin-state direction bar...\n")

chrom_summary <- read_tsv(file.path(TABLES_DIR, "chromatin_state_summary.tsv"),
                          show_col_types = FALSE)

chrom_long <- chrom_summary %>%
  dplyr::select(chromatin_state, Hypermethylated = count_Hypermethylated,
         Hypomethylated = count_Hypomethylated) %>%
  pivot_longer(c(Hypermethylated, Hypomethylated),
               names_to = "direction", values_to = "count") %>%
  group_by(chromatin_state) %>%
  mutate(frac = count / sum(count)) %>%
  ungroup() %>%
  mutate(chromatin_state = factor(chromatin_state, levels = rev(CHROMATIN_STATE_ORDER)))

p2a <- ggplot(chrom_long, aes(x = chromatin_state, y = frac, fill = direction)) +
  geom_col(position = "fill", width = 0.7) +
  geom_text(aes(label = ifelse(frac > 0.15, sprintf("%d%%", round(100 * frac)), "")),
            position = position_fill(vjust = 0.5), size = 2.8, color = "white") +
  coord_flip() +
  scale_fill_manual(values = c("Hypermethylated" = "#D7191C", "Hypomethylated" = "#2C7BB6")) +
  scale_y_continuous(labels = percent_format()) +
  labs(x = NULL, y = "Fraction of significant DMRs",
       subtitle = "Active_Promoter: 93% hyper; Repressed: 94% hypo") +
  theme_biomodal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.subtitle = element_text(size = 9))

save_multiformat_ggplot(p2a,
                        file.path(GRANT_DIR, "93d_chromatin_direction_bar"),
                        width = 7, height = 5)

# ---- Panel B: K119ub Logistic Forest Plot ------------------------------------

cat("  Panel B: K119ub logistic forest plot...\n")

coefs <- read_tsv(file.path(TABLES_DIR, "diffbind_logistic_model_coefficients.tsv"),
                  show_col_types = FALSE)

coefs <- coefs %>%
  mutate(display_name = factor(display_name, levels = rev(display_name)))

p2b <- ggplot(coefs, aes(x = or, y = display_name)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_pointrange(aes(xmin = or_lower, xmax = or_upper),
                  size = 0.5, linewidth = 0.8, color = "black") +
  geom_text(aes(label = sprintf("OR = %.2f %s", or, sig_label)),
            hjust = -0.15, size = 2.8) +
  scale_x_log10() +
  labs(x = "Odds Ratio (hypermethylation)", y = NULL,
       subtitle = "4-mark logistic model: K119ub dominant") +
  theme_biomodal() +
  theme(plot.subtitle = element_text(size = 9))

save_multiformat_ggplot(p2b,
                        file.path(GRANT_DIR, "93e_k119ub_forest"),
                        width = 7, height = 4)

# ---- Panel C: A-Compartment Enrichment Bars ----------------------------------

cat("  Panel C: A-compartment enrichment bars...\n")

comp_tests <- read_tsv(file.path(TABLES_DIR, "compartment_fisher_tests.tsv"),
                       show_col_types = FALSE)

comp_a <- comp_tests %>%
  filter(grepl("A enriched|A compartment", test, ignore.case = TRUE)) %>%
  mutate(test_short = gsub(" -> A enriched", "", test),
         test_short = factor(test_short, levels = rev(test_short)))

p2c <- ggplot(comp_a, aes(x = odds_ratio, y = test_short)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper),
                  size = 0.5, linewidth = 0.8, color = "#D7191C") +
  geom_text(aes(label = sprintf("OR = %.1f", odds_ratio)),
            hjust = -0.15, size = 3) +
  scale_x_log10() +
  labs(x = "Odds Ratio (enrichment in A compartment)", y = NULL,
       subtitle = "Hypermethylation targets euchromatin (A compartment)") +
  theme_biomodal() +
  theme(plot.subtitle = element_text(size = 9))

save_multiformat_ggplot(p2c,
                        file.path(GRANT_DIR, "93f_a_compartment_enrichment"),
                        width = 7, height = 4)

# ---- Figure 2 Composite -----------------------------------------------------

cat("  Compositing Grant Figure 2...\n")

fig2 <- p2a | p2b | p2c
fig2 <- fig2 + plot_annotation(tag_levels = "A")

save_multiformat_ggplot(fig2,
                        file.path(GRANT_DIR, "grant_fig2_euchromatin"),
                        width = 20, height = 5)

cat("  Done.\n\n")

# =============================================================================
# GRANT FIGURE 3: Heterochromatin Protection
# =============================================================================

cat("── GRANT FIGURE 3: Heterochromatin Protection ──────────────────────────────\n\n")

# ---- Panel A: Per-State Hyper Rate Lollipop ----------------------------------

cat("  Panel A: Per-state hypermethylation rate lollipop...\n")

state_enrich <- read_tsv(file.path(TABLES_DIR,
                         "polycomb_per_chromatin_state_enrichment.tsv"),
                         show_col_types = FALSE)

state_enrich <- state_enrich %>%
  filter(chromatin_state %in% CHROMATIN_STATE_ORDER) %>%
  mutate(chromatin_state = factor(chromatin_state, levels = rev(CHROMATIN_STATE_ORDER)))

p3a <- ggplot(state_enrich, aes(x = pct_hyper, y = chromatin_state)) +
  geom_segment(aes(x = 0, xend = pct_hyper, yend = chromatin_state,
                   color = chromatin_state),
               linewidth = 1) +
  geom_point(aes(color = chromatin_state), size = 3) +
  geom_text(aes(label = sprintf("%.1f%%", pct_hyper)),
            hjust = -0.3, size = 2.8) +
  scale_color_manual(values = CHROMATIN_STATE_COLORS, guide = "none") +
  scale_x_continuous(limits = c(0, max(state_enrich$pct_hyper) * 1.2)) +
  labs(x = "% of genes hypermethylated", y = NULL,
       subtitle = "Active Promoter 72% vs Polycomb 1.8%") +
  theme_biomodal() +
  theme(plot.subtitle = element_text(size = 9))

save_multiformat_ggplot(p3a,
                        file.path(GRANT_DIR, "93g_per_state_hyper_rate"),
                        width = 7, height = 5)

# ---- Panel B: Polycomb Exclusion Forest --------------------------------------

cat("  Panel B: Polycomb exclusion forest plot...\n")

poly_tests <- read_tsv(file.path(TABLES_DIR, "polycomb_fisher_tests.tsv"),
                       show_col_types = FALSE)

poly_plot <- poly_tests %>%
  mutate(test_short = gsub("Chromatin State × ", "", test),
         test_short = factor(test_short, levels = rev(test_short)))

p3b <- ggplot(poly_plot, aes(x = odds_ratio, y = test_short)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper),
                  size = 0.5, linewidth = 0.8, color = "#4daf4a") +
  geom_text(aes(label = sprintf("OR = %.3f", odds_ratio)),
            hjust = -0.15, size = 2.8) +
  scale_x_log10() +
  labs(x = "Odds Ratio (Polycomb target enrichment)", y = NULL,
       subtitle = "Polycomb EXCLUDED from hypermethylation (OR = 0.063)") +
  theme_biomodal() +
  theme(plot.subtitle = element_text(size = 9))

save_multiformat_ggplot(p3b,
                        file.path(GRANT_DIR, "93h_polycomb_exclusion_forest"),
                        width = 7, height = 4)

# ---- Panel C: B-Compartment Enrichment Bars ----------------------------------

cat("  Panel C: B-compartment enrichment bars...\n")

comp_b <- comp_tests %>%
  filter(grepl("B enriched|B compartment", test, ignore.case = TRUE)) %>%
  mutate(test_short = gsub(" -> B enriched", "", test),
         test_short = factor(test_short, levels = rev(test_short)))

p3c <- ggplot(comp_b, aes(x = odds_ratio, y = test_short)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper),
                  size = 0.5, linewidth = 0.8, color = "#2C7BB6") +
  geom_text(aes(label = sprintf("OR = %.2f", odds_ratio)),
            hjust = -0.15, size = 3) +
  scale_x_log10() +
  labs(x = "Odds Ratio (enrichment in B compartment)", y = NULL,
       subtitle = "Hypomethylation enriched in heterochromatin (B compartment)") +
  theme_biomodal() +
  theme(plot.subtitle = element_text(size = 9))

save_multiformat_ggplot(p3c,
                        file.path(GRANT_DIR, "93i_b_compartment_enrichment"),
                        width = 7, height = 4)

# ---- Figure 3 Composite -----------------------------------------------------

cat("  Compositing Grant Figure 3...\n")

fig3 <- p3a | p3b | p3c
fig3 <- fig3 + plot_annotation(tag_levels = "A")

save_multiformat_ggplot(fig3,
                        file.path(GRANT_DIR, "grant_fig3_heterochromatin"),
                        width = 20, height = 5)

cat("  Done.\n\n")

# =============================================================================
# GRANT FIGURE 4: MeCP2 Reads Chromatin State, Not Methylation
# =============================================================================

cat("── GRANT FIGURE 4: MeCP2 Reads Chromatin ──────────────────────────────────\n\n")

# ---- Panel A: R² Comparison Bars --------------------------------------------

cat("  Panel A: R² comparison (CG vs Chromatin)...\n")

model_r2 <- read_tsv(file.path(TABLES_DIR, "62_model_comparison_summary.tsv"),
                     show_col_types = FALSE)

model_r2_plot <- model_r2 %>%
  filter(model %in% c("CG only", "Chromatin only", "Full")) %>%
  mutate(model = factor(model, levels = c("CG only", "Chromatin only", "Full")))

model_colors <- c("CG only" = "#E41A1C", "Chromatin only" = "#756BB1", "Full" = "#333333")

p4a <- ggplot(model_r2_plot, aes(x = model, y = r_squared, fill = model)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("R² = %.3f", r_squared)),
            vjust = -0.5, size = 3) +
  facet_wrap(~type, nrow = 1) +
  scale_fill_manual(values = model_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = expression(R^2),
       subtitle = "CG methylation explains <2% of MeCP2; chromatin explains ~25%") +
  theme_biomodal() +
  theme(plot.subtitle = element_text(size = 9))

save_multiformat_ggplot(p4a,
                        file.path(GRANT_DIR, "93j_mecp2_r2_comparison"),
                        width = 9, height = 5)

# ---- Panel B: K119ub at MeCP2-Up-No-Methylation Genes -----------------------

cat("  Panel B: K119ub gained at MeCP2-no-meth genes (linchpin)...\n")

stats_67 <- read_tsv(file.path(TABLES_DIR, "67_statistics.tsv"),
                     show_col_types = FALSE)

fisher_row <- stats_67 %>% filter(grepl("Fisher", test))
our_pct <- fisher_row$our_pct_gained
rest_pct <- fisher_row$rest_pct_gained
fisher_or <- as.numeric(gsub("OR = ", "", fisher_row$statistic))
fisher_p <- fisher_row$p_value

bar_data <- data.frame(
  group = c("MeCP2-up\nno methylation Δ\n(359 genes)", "Genome\nbackground"),
  pct = c(our_pct, rest_pct),
  stringsAsFactors = FALSE
)

p4b <- ggplot(bar_data, aes(x = group, y = pct, fill = group)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("MeCP2-up\nno methylation Δ\n(359 genes)" = "#D95F02",
                               "Genome\nbackground" = "grey60"),
                    guide = "none") +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.1))) +
  annotate("text", x = 1.5, y = 90,
           label = sprintf("Fisher OR = %.2f\np = %.1e", fisher_or, fisher_p),
           size = 3, fontface = "italic") +
  labs(x = NULL, y = "% gaining H2AK119ub",
       subtitle = "MeCP2 follows K119ub even without methylation change") +
  theme_biomodal() +
  theme(plot.subtitle = element_text(size = 9))

save_multiformat_ggplot(p4b,
                        file.path(GRANT_DIR, "93k_k119ub_mecp2_no_meth"),
                        width = 5, height = 5)

# ---- Panel C: Mirror-Image Mark Profiles -------------------------------------

cat("  Panel C: Mirror-image mark profiles by MeCP2 status...\n")

mecp2_stats <- read_tsv(file.path(TABLES_DIR, "60_mecp2_status_stats.tsv"),
                        show_col_types = FALSE)

mark_labels <- c("k119ub_bw" = "H2AK119ub", "mc_bw" = "5mC", "hmc_bw" = "5hmC",
                 "k27ac_bw" = "H3K27ac", "k27me3_bw" = "H3K27me3",
                 "k119ub_db" = "K119ub\n(DiffBind)", "k27ac_db" = "K27ac\n(DiffBind)",
                 "k27me3_db" = "K27me3\n(DiffBind)")

mirror_data <- mecp2_stats %>%
  filter(mark %in% c("k119ub_bw", "mc_bw", "hmc_bw", "k27ac_bw", "k27me3_bw")) %>%
  mutate(
    mark_label = recode(mark, !!!mark_labels),
    mark_label = factor(mark_label, levels = c("H2AK119ub", "5mC", "5hmC",
                                                "H3K27ac", "H3K27me3")),
    group = factor(group, levels = c("MeCP2 Up", "MeCP2 Down", "Not Significant"))
  )

group_colors <- c("MeCP2 Up" = "#D95F02", "MeCP2 Down" = "#7570B3",
                  "Not Significant" = "grey60")

p4c <- ggplot(mirror_data %>% filter(group != "Not Significant"),
              aes(x = mark_label, y = mean, color = group)) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.3) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = 0, ymax = mean),
                 position = position_dodge(width = 0.5), linewidth = 0.8) +
  scale_color_manual(values = group_colors) +
  labs(x = NULL, y = "Mean signal change (mutant − control)",
       subtitle = "MeCP2-Up mirrors all marks; MeCP2-Down: K119ub FLAT") +
  theme_biomodal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.subtitle = element_text(size = 9))

save_multiformat_ggplot(p4c,
                        file.path(GRANT_DIR, "93l_mirror_profiles"),
                        width = 7, height = 5)

# ---- Figure 4 Composite -----------------------------------------------------

cat("  Compositing Grant Figure 4...\n")

fig4 <- p4a | p4b | p4c
fig4 <- fig4 + plot_annotation(tag_levels = "A")

save_multiformat_ggplot(fig4,
                        file.path(GRANT_DIR, "grant_fig4_mecp2_chromatin"),
                        width = 20, height = 5)

cat("  Done.\n\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("================================================================================\n")
cat("SECTION 93 COMPLETE\n")
cat("================================================================================\n")
cat(sprintf("Output directory: %s\n", GRANT_DIR))
cat("Generated:\n")
cat("  4 composite figures (grant_fig1–4)\n")
cat("  12 individual panels (93a–93l)\n")
cat("All in PDF + SVG + JPG via multi_format_output.R\n")
