# scripts/viz_sections/section_88_supplemental_evidence.R
# =============================================================================
# Figure S1: Supplemental Evidence (BAP1-KO cerebellum methylation paper)
#
# Supporting panels that corroborate the main story but are too detailed for
# the main figures. Six panels, 3x2 grid (design = "AABB\nCCDD\nEEFF"),
# 180 x 240 mm, saved PDF + SVG + JPG under plots/figures/figureS1_supplemental/.
#
#   Panel A (Sec 37) — Gene-level permutation validation: 15 label-shuffle
#                      permutation tests, all Confirmed (perm-z forest).
#                      Source: permutation_37_summary.tsv
#   Panel B (Sec 45) — Field 2019 chr8 cross-species concordance at gene bodies
#                      (Fisher p = 0.0089, hypergeometric p = 4.6e-6).
#                      Source: field_chr8_comparison_full.tsv,
#                              field_chr8_statistical_tests.tsv
#   Panel C (Sec 47) — CTCF loop-anchor mC-hypermethylation: lost anchors are
#                      enriched for hyper at flanking dynamic regions (OR = 3.28).
#                      Source: 47a_fisher_results.tsv
#                              (47a_dynamic_ctcf_anchor_methylation.tsv read for n)
#   Panel D (Sec 48) — CpG-island H2AK119ub DEPLETION at hypermethylated islands
#                      (mC-Hyper 17.2% vs Non-sig 30.2% K119ub-mutant overlap).
#                      Source: 48_cpg_island_ubiquitination_summary.tsv
#   Panel E (Sec 44) — Allele-specific methylation doubling (mutant 1.95x more
#                      significant mC ASM sites per sample).
#                      Source: asm_mc_significant_per_sample.tsv
#                              (asm_dmr_overlap_summary.tsv read for DMR linkage n)
#   Panel F (Sec 78) — Unbiased broad-neuronal stoichiometry self-correction:
#                      total methylation DECREASES at neuronal genes (-0.0022),
#                      rises only at coordinated / MeCP2-Up sets.
#                      Source: 78_total_methylation_summary.tsv
#
# Run from the downstream/ directory (getwd() == downstream).
# =============================================================================

# -----------------------------------------------------------------------------
# Shared infrastructure: _shared_config.R first (data + helpers + COLORS),
# then _figure_config.R (theme_pub, save_figure, read_table_tsv, ...).
# -----------------------------------------------------------------------------
source("scripts/viz_sections/_shared_config.R")
source("scripts/viz_sections/_figure_config.R")

library(patchwork)

SEC88_DIR <- file.path(OUTPUT_DIR, "88_supplemental_evidence")
dir.create(SEC88_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, w = 10, h = 7) {
  save_multiformat_ggplot(p,
    base_path = file.path(SEC88_DIR, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

cat("\n========== Section 88: Supplemental Evidence ==========\n")

# =============================================================================
# Panel A — Section 37: gene-level permutation validation (15/15 Confirmed)
# =============================================================================
# Forest of permutation z-scores for the 15 label-shuffle tests; sign of perm_z
# matches the direction of the Fisher OR (positive = enrichment, negative =
# depletion). Every test is "Confirmed" in the source TSV.
cat("\n--- Panel A: permutation validation (Sec 37) ---\n")

perm <- read_table_tsv("permutation_37_summary.tsv")

n_perm_total     <- nrow(perm)                                  # 15
n_perm_confirmed <- sum(perm$concordance == "Confirmed")        # 15
stopifnot(n_perm_total == 15L, n_perm_confirmed == 15L)

# Order tests by perm_z so the forest reads as a clean gradient; label each
# row with its test_id + short description.
perm$row_label <- sprintf("%s  %s", perm$test_id, perm$description)
perm <- perm[order(perm$perm_z), ]
perm$row_label <- factor(perm$row_label, levels = perm$row_label)

panel_a <- ggplot(perm, aes(x = perm_z, y = row_label, colour = source_section)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3,
             colour = "grey50") +
  geom_segment(aes(x = 0, xend = perm_z, yend = row_label), linewidth = 0.4) +
  geom_point(size = 1.6) +
  scale_colour_brewer(palette = "Paired", name = "Source\nsection") +
  labs(
    title = "Gene-level permutation validation",
    subtitle = sprintf("%d / %d enrichments Confirmed (chr-stratified label shuffle, 1e5 perms)",
                       n_perm_confirmed, n_perm_total),
    x = "Permutation z-score",
    y = NULL
  ) +
  theme_pub(base_size = 7) +
  theme(legend.position = "right",
        legend.key.size = unit(2.5, "mm"),
        axis.text.y = element_text(size = 5.5))

# =============================================================================
# Panel B — Section 45: Field 2019 chr8 cross-species concordance (gene body)
# =============================================================================
# Gene-body mC direction concordance of Field human BAP1 chr8 hotspot orthologs
# in our mouse cerebellum. 21 Concordant / 38 Discordant / 18 NS / 4 Not-in-data.
cat("\n--- Panel B: Field chr8 cross-species concordance (Sec 45) ---\n")

field      <- read_table_tsv("field_chr8_comparison_full.tsv")
field_stat <- read_table_tsv("field_chr8_statistical_tests.tsv")

# Tabulate the gene-body concordance categories straight from the master table.
conc_levels <- c("Concordant", "Discordant", "Non-significant", "Not in data")
field_counts <- as.data.frame(table(factor(field$gb_concordance, levels = conc_levels)),
                              stringsAsFactors = FALSE)
names(field_counts) <- c("category", "n")
stopifnot(field_counts$n[field_counts$category == "Concordant"] == 21,
          field_counts$n[field_counts$category == "Discordant"] == 38)

# Pull the two load-bearing p-values verbatim from the stats table.
fisher_gb_p <- field_stat$p_value[field_stat$test == "Fisher's exact (gene-body mC concordance)"]
hyper_gb_p  <- field_stat$p_value[field_stat$test == "Hypergeometric enrichment (gene-body mC)"]
stopifnot(length(fisher_gb_p) == 1, length(hyper_gb_p) == 1)

field_counts$category <- factor(field_counts$category, levels = conc_levels)
field_fill <- c(
  "Concordant"      = unname(COLORS$direction[["Hypermethylated"]]),
  "Discordant"      = unname(COLORS$direction[["Hypomethylated"]]),
  "Non-significant" = "grey70",
  "Not in data"     = "grey85"
)

panel_b <- ggplot(field_counts, aes(x = category, y = n, fill = category)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_text(aes(label = n), vjust = -0.3, size = 2.4) +
  scale_fill_manual(values = field_fill, guide = "none") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4, size = 2.2,
           label = sprintf("Fisher %s\nHypergeometric %s",
                           stat_label("p", fisher_gb_p, "%.4f"),
                           stat_label("p", hyper_gb_p, "%.1e"))) +
  labs(
    title = "Field chr8 hotspot concordance",
    subtitle = "Gene-body 5mC direction, Field orthologs (n = 81)",
    x = NULL, y = "Field ortholog genes"
  ) +
  theme_pub(base_size = 7) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# =============================================================================
# Panel C — Section 47: CTCF loop-anchor mC-hypermethylation (lost vs gained)
# =============================================================================
# Lost CTCF anchors are enriched for mC-hypermethylation at flanking dynamic
# CpG regions (shores+shelves) vs gained anchors. Percentages and OR are read
# directly from the canonical Fisher result row (Lost_vs_Gained, mC_hyper).
cat("\n--- Panel C: CTCF anchor methylation (Sec 47) ---\n")

ctcf_fisher <- read_table_tsv("47a_fisher_results.tsv")
ctcf_row <- ctcf_fisher[ctcf_fisher$test_label == "Lost_vs_Gained" &
                          ctcf_fisher$ref_label == "mC_hyper", ]
stopifnot(nrow(ctcf_row) == 1)

ctcf_or    <- ctcf_row$odds_ratio          # 3.277...
ctcf_p     <- ctcf_row$p_value             # 5.36e-24
ctcf_lost  <- ctcf_row$test_pct            # 22.32 %
ctcf_gain  <- ctcf_row$ref_pct             # 8.06 %
ctcf_lost_n <- ctcf_row$test_n             # 1084 dynamic regions
ctcf_gain_n <- ctcf_row$ref_n              # 1427 dynamic regions
stopifnot(round(ctcf_or, 2) == 3.28)

ctcf_df <- data.frame(
  anchor = factor(c("Lost CTCF\nanchor", "Gained CTCF\nanchor"),
                  levels = c("Lost CTCF\nanchor", "Gained CTCF\nanchor")),
  pct    = c(ctcf_lost, ctcf_gain),
  n      = c(ctcf_lost_n, ctcf_gain_n),
  stringsAsFactors = FALSE
)

panel_c <- ggplot(ctcf_df, aes(x = anchor, y = pct, fill = anchor)) +
  geom_col(width = 0.65, colour = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.3, size = 2.6) +
  scale_fill_manual(
    values = c("Lost CTCF\nanchor" = unname(COLORS$direction[["Hypermethylated"]]),
               "Gained CTCF\nanchor" = unname(COLORS$direction[["Hypomethylated"]])),
    guide = "none"
  ) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4, size = 2.4,
           label = sprintf("%s\n%s",
                           stat_label("OR", ctcf_or, "%.2f"),
                           stat_label("p", ctcf_p, "%.1e"))) +
  labs(
    title = "CTCF anchor 5mC hypermethylation",
    subtitle = "% mC-hyper at flanking dynamic CpG regions",
    x = NULL, y = "% regions hypermethylated"
  ) +
  theme_pub(base_size = 7)

# =============================================================================
# Panel D — Section 48: CpG-island H2AK119ub depletion at hyper islands
# =============================================================================
# Counter to the naive "K119ub gain -> island methylation" hypothesis: K119ub
# peak overlap is DEPLETED at hypermethylated CpG islands relative to
# non-significant islands.
cat("\n--- Panel D: CpG island K119ub depletion (Sec 48) ---\n")

cpgi <- read_table_tsv("48_cpg_island_ubiquitination_summary.tsv")

# Compute K119ub-mutant peak-overlap fraction per DMR category (0/1 flag column).
cpgi$dmr_category <- factor(cpgi$dmr_category,
                            levels = c("mC Hyper", "mC Hypo", "Non-sig"))
cpgi_overlap <- do.call(rbind, lapply(levels(cpgi$dmr_category), function(g) {
  sub <- cpgi[cpgi$dmr_category == g, ]
  data.frame(
    dmr_category = g,
    n            = nrow(sub),
    pct_k119ub   = 100 * mean(sub$k119ub_mut == 1),
    pct_gained   = 100 * mean(sub$k119ub_gained == 1),
    stringsAsFactors = FALSE
  )
}))
cpgi_overlap$dmr_category <- factor(cpgi_overlap$dmr_category,
                                    levels = c("mC Hyper", "mC Hypo", "Non-sig"))

# Confirm the depletion direction vs the documented numbers (17.2% vs 30.2%).
hyper_k119  <- cpgi_overlap$pct_k119ub[cpgi_overlap$dmr_category == "mC Hyper"]
nonsig_k119 <- cpgi_overlap$pct_k119ub[cpgi_overlap$dmr_category == "Non-sig"]
stopifnot(round(hyper_k119, 1) == 17.2, round(nonsig_k119, 1) == 30.2)

cpgi_fill <- c("mC Hyper" = unname(COLORS$direction[["Hypermethylated"]]),
               "mC Hypo"  = unname(COLORS$direction[["Hypomethylated"]]),
               "Non-sig"  = "grey70")

panel_d <- ggplot(cpgi_overlap, aes(x = dmr_category, y = pct_k119ub,
                                    fill = dmr_category)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f%%", pct_k119ub)),
            vjust = -0.3, size = 2.6) +
  scale_fill_manual(values = cpgi_fill, guide = "none") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4, size = 2.4,
           label = "K119ub depleted\nat hyper islands") +
  labs(
    title = "CpG-island H2AK119ub depletion",
    subtitle = "% islands with mutant K119ub peak overlap",
    x = NULL, y = "% islands with K119ub overlap"
  ) +
  theme_pub(base_size = 7)

# =============================================================================
# Panel E — Section 44: allele-specific methylation doubling (mutant 1.95x)
# =============================================================================
# Per-sample counts of significant mC ASM sites; mutant has ~2x more. Points =
# samples (shaped by sex), bars = condition means. Ratio annotated.
cat("\n--- Panel E: ASM doubling (Sec 44) ---\n")

asm <- read_table_tsv("asm_mc_significant_per_sample.tsv")
stopifnot(all(c("condition", "sig_count", "sex") %in% names(asm)))

asm$condition <- factor(asm$condition, levels = c("Control", "Mutant"))
asm_means <- aggregate(sig_count ~ condition, data = asm, FUN = mean)
ctrl_mean <- asm_means$sig_count[asm_means$condition == "Control"]
mut_mean  <- asm_means$sig_count[asm_means$condition == "Mutant"]
asm_ratio <- mut_mean / ctrl_mean
stopifnot(round(asm_ratio, 2) == 1.95)

# Linkage context: number of gene-body DMRs harbouring >=1 significant ASM site.
asm_overlap <- read_table_tsv("asm_dmr_overlap_summary.tsv")
n_asm_dmr_genes <- nrow(asm_overlap)   # 1408

panel_e <- ggplot(asm, aes(x = condition, y = sig_count)) +
  geom_col(data = asm_means, aes(x = condition, y = sig_count, fill = condition),
           width = 0.6, alpha = 0.55, colour = "black", linewidth = 0.2) +
  geom_point(aes(shape = sex), size = 1.8, colour = "black",
             position = position_jitter(width = 0.12, height = 0, seed = 1)) +
  scale_fill_manual(values = c("Control" = unname(COLORS$condition[["Control"]]),
                               "Mutant"  = unname(COLORS$condition[["Mutant"]])),
                    guide = "none") +
  scale_shape_manual(values = c("Female" = 16, "Male" = 17), name = "Sex") +
  scale_y_continuous(labels = scales::comma) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4, size = 2.4,
           label = sprintf("Mutant %.2fx\n(%s DMR genes\nhave ASM sites)",
                           asm_ratio, format(n_asm_dmr_genes, big.mark = ","))) +
  labs(
    title = "Allele-specific methylation doubling",
    subtitle = "Significant mC ASM sites per sample (n = 4 / condition)",
    x = NULL, y = "Significant mC ASM sites"
  ) +
  theme_pub(base_size = 7) +
  theme(legend.position = "right", legend.key.size = unit(2.5, "mm"))

# =============================================================================
# Panel F — Section 78: unbiased broad-neuronal stoichiometry self-correction
# =============================================================================
# Mean change in total methylation (mut - ctrl) by gene group with 95% CIs.
# Total methylation DECREASES at (unbiased GO-derived) neuronal genes and only
# rises where it should (coordinated, MeCP2-Up). Dashed line at zero.
cat("\n--- Panel F: broad-neuronal total-methylation forest (Sec 78) ---\n")

tot <- read_table_tsv("78_total_methylation_summary.tsv")
stopifnot(all(c("group", "mean_delta", "ci_lo", "ci_hi", "sig") %in% names(tot)))

# Confirm the headline self-correction value: broad neuronal mean = -0.0022.
neur_broad_delta <- tot$mean_delta[tot$group == "Neuronal (broad)"]
stopifnot(round(neur_broad_delta, 4) == -0.0022)

# Order groups for the forest (declining biological breadth, with the two
# methylation-defined "up" sets at the bottom for contrast).
group_order <- c("All genes", "Non-neuronal", "Neuronal (broad)", "Synapse/axon",
                 "Neur + Coord", "Coordinated", "MeCP2-Up")
tot <- tot[tot$group %in% group_order, ]
tot$group <- factor(tot$group, levels = rev(group_order))
tot$dir   <- ifelse(tot$mean_delta >= 0, "Increase", "Decrease")

panel_f <- ggplot(tot, aes(x = mean_delta, y = group, colour = dir)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3,
             colour = "grey50") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), width = 0.25,
                linewidth = 0.4, orientation = "y") +
  geom_point(size = 1.8) +
  geom_text(aes(label = sig), vjust = -0.9, size = 2.2, colour = "black") +
  scale_colour_manual(
    values = c("Decrease" = unname(COLORS$direction[["Hypomethylated"]]),
               "Increase" = unname(COLORS$direction[["Hypermethylated"]])),
    name = "Total mC"
  ) +
  labs(
    title = "Total methylation by gene set (unbiased)",
    subtitle = sprintf("Neuronal (broad) Delta = %s (decreases); rises at methylation-defined sets (coordinated / MeCP2-Up)",
                       sprintf("%+.4f", neur_broad_delta)),
    x = "Mean total-methylation change (mut - ctrl)",
    y = NULL
  ) +
  theme_pub(base_size = 7) +
  theme(legend.position = "right", legend.key.size = unit(2.5, "mm"))

# =============================================================================
# Composite assembly: 3x2 grid (design = "AABB\nCCDD\nEEFF"), 180 x 240 mm
# =============================================================================
cat("\n--- Assembling Figure S1 composite ---\n")

design_s1 <- "AABB\nCCDD\nEEFF"

sec88_composite <- panel_a + panel_b + panel_c + panel_d + panel_e + panel_f +
  plot_layout(design = design_s1, heights = c(1, 1, 1)) +
  plot_annotation(
    title = "Section 88: Supplemental evidence",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      plot.tag   = element_text(size = 12, face = "bold")
    )
  )

save_plot(panel_a, "88a_permutation",   w = 10, h = 7)
save_plot(panel_b, "88b_field_chr8",     w = 8, h = 6)
save_plot(panel_c, "88c_ctcf_anchor",    w = 8, h = 6)
save_plot(panel_d, "88d_cpgi_k119ub",    w = 8, h = 6)
save_plot(panel_e, "88e_asm_doubling",   w = 8, h = 6)
save_plot(panel_f, "88f_total_meth",     w = 10, h = 6)
save_plot(sec88_composite, "88_composite", w = 14, h = 16)

cat("\n========== Section 88 complete ==========\n")
cat(sprintf("  Panel A: %d/%d permutation tests Confirmed\n",
            n_perm_confirmed, n_perm_total))
cat(sprintf("  Panel B: Field chr8 gene-body Fisher p = %.4f, hypergeometric p = %.1e\n",
            fisher_gb_p, hyper_gb_p))
cat(sprintf("  Panel C: CTCF lost-vs-gained mC-hyper OR = %.2f (p = %.1e)\n",
            ctcf_or, ctcf_p))
cat(sprintf("  Panel D: CpG-island K119ub overlap %.1f%% (hyper) vs %.1f%% (non-sig)\n",
            hyper_k119, nonsig_k119))
cat(sprintf("  Panel E: ASM mutant/control ratio = %.2fx (%s DMR genes carry ASM)\n",
            asm_ratio, format(n_asm_dmr_genes, big.mark = ",")))
cat(sprintf("  Panel F: broad-neuronal total-methylation Delta = %+.4f\n",
            neur_broad_delta))
