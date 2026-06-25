# scripts/viz_sections/section_87_mechanism_summary.R
# =============================================================================
# FIGURE 7: Mechanism Summary Data Panel (cascade)
#
# Paper mapping: All R-sections (cascade summary). Data companion to the
# Illustrator model schematic. Each panel quantifies ONE link in the chain:
#   BAP1 loss -> H2AK119ub accumulation -> TET-block -> coordinated mC up /
#   hmC down -> MeCP2 redistribution -> convergence on neuronal genes.
#
# Panels (per approved plan "## Figure 7: Mechanism Summary Data Panel"):
#   A  BAP1 -> K119ub      Effect-size bar   k119ub_bigwig_signal_summary.tsv
#                          Global K119ub increase: All DMR genes median log2FC
#                          = +0.007 (p vs 0 = 1.81e-20); mC Up +0.059, mC Down -0.080.
#   B  K119ub -> mC        OR cascade waterfall  diffbind_logistic_model_coefficients.tsv
#                          4 marks ordered by OR; H2AK119ub OR = 4.71 (dominant).
#   C  TET-block evidence  AUC comparison    dnmt3a_model_comparison.tsv +
#                          baseline_hmc_model_comparison.tsv
#                          TET impediment AUC 0.793 >> DNMT3A recruitment 0.696
#                          (DeLong p = 9.43e-49, from dnmt3a_exclusive_model_comparison.tsv).
#   D  mC -> MeCP2         R-squared cascade  62_model_comparison_summary.tsv
#                          Chromatin-only R2 0.246 >> CG-only R2 0.017 (binding level).
#   E  Neuronal convergence  OR summary      74_pairwise_fisher.tsv +
#                          72_fisher_results.tsv
#                          Coordinated x MeCP2-Up OR = 5.16; K119ub top-decile
#                          neuronal-enrichment OR = 1.70.
#
# Layout (plan): design = "AABB\nCCDD\n#EE#", heights = c(1, 1, 0.6)
# Dimensions (plan, full-width short option): 180 x 140 mm
#
# Output: PDF + SVG + JPG into plots/figures/figure7_mechanism_summary/
#
# NOTE: MeCP2 is CUT&RUN (not an antibody-pulldown assay) — all MeCP2 labels
#       use "signal" / "binding", never the pulldown-assay term.
# =============================================================================

# -----------------------------------------------------------------------------
# Shared infrastructure (must be sourced first, from downstream/ as cwd).
# _shared_config.R provides COLORS, theme_biomodal(), TABLES_DIR, the save
# helpers, and the pre-loaded DMR objects. _figure_config.R adds theme_pub(),
# save_figure(), read_table_tsv(), add_panel_labels(), stat_label(), PUB_COLORS.
# -----------------------------------------------------------------------------
source("scripts/viz_sections/_shared_config.R")
source("scripts/viz_sections/_figure_config.R")

SEC87_DIR <- file.path(OUTPUT_DIR, "87_mechanism_summary")
dir.create(SEC87_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, w = 10, h = 7) {
  save_multiformat_ggplot(p,
    base_path = file.path(SEC87_DIR, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

cat("\n=== FIGURE 7: Mechanism Summary Data Panel (cascade) ===\n")

# =============================================================================
# PANEL A — BAP1 -> K119ub : global K119ub effect-size bar (Section 18)
# Source: k119ub_bigwig_signal_summary.tsv
# Shows per-group median K119ub log2FC (mut/ctrl). The "All DMR Genes"
# background bar (+0.007, p vs 0 = 1.81e-20) is the evidence for a global
# K119ub increase; mC-Up genes shift positive, mC-Down negative.
# =============================================================================

k119_sig <- read_table_tsv("k119ub_bigwig_signal_summary.tsv")

# Confirm the columns the panel relies on are present (fail loudly otherwise).
req_a <- c("met_group", "median_fc", "p_vs_zero")
if (!all(req_a %in% names(k119_sig))) {
  stop("Panel A: k119ub_bigwig_signal_summary.tsv missing columns: ",
       paste(setdiff(req_a, names(k119_sig)), collapse = ", "))
}

# Order groups as a cascade-readable sequence: background first, then the two
# directional DMR classes that carry the signal, then the hmC mirror groups.
panelA_order <- c("All DMR Genes", "mC Up", "mC Down", "hmC Down", "hmC Up")
k119_sig <- k119_sig[k119_sig$met_group %in% panelA_order, ]
k119_sig$met_group <- factor(k119_sig$met_group, levels = rev(panelA_order))

# Colour: the genome-wide background bar uses the canonical K119ub purple to
# anchor the "global increase" message; directional bars use the methylation
# hyper/hypo palette so the reader maps them onto the mC story.
k119_sig$bar_fill <- ifelse(
  k119_sig$met_group == "All DMR Genes", "Background (all DMR genes)",
  ifelse(grepl("Up$", k119_sig$met_group) & grepl("^mC", k119_sig$met_group),
         "mC hyper genes",
         ifelse(grepl("Down$", k119_sig$met_group) & grepl("^mC", k119_sig$met_group),
                "mC hypo genes", "hmC genes")))

panelA_cols <- c(
  "Background (all DMR genes)" = COLORS$k119ub[["K119ub Gained"]],  # "#756BB1"
  "mC hyper genes"            = COLORS$direction[["Hypermethylated"]],                   # "#D7191C"
  "mC hypo genes"             = COLORS$direction[["Hypomethylated"]],                    # "#2C7BB6"
  "hmC genes"                 = "grey60"
)

# Background median log2FC for the annotation (TSV-confirmed = 0.006956).
bg_median <- k119_sig$median_fc[k119_sig$met_group == "All DMR Genes"]
bg_p      <- k119_sig$p_vs_zero[k119_sig$met_group == "All DMR Genes"]

p_A <- ggplot(k119_sig,
              aes(x = median_fc, y = met_group, fill = bar_fill)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "black") +
  scale_fill_manual(values = panelA_cols, name = NULL) +
  labs(
    title = "BAP1 loss -> H2AK119ub",
    subtitle = sprintf("Global K119ub rise: background %s (p = %s)",
                       stat_label("median log2FC", bg_median, "%+.3f"),
                       formatC(bg_p, format = "e", digits = 2)),
    x = "K119ub signal log2FC (mut / ctrl)",
    y = NULL
  ) +
  theme_pub() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 6))

# =============================================================================
# PANEL B — K119ub -> mC : multi-mark OR waterfall (Section 33)
# Source: diffbind_logistic_model_coefficients.tsv
# Logistic model of hypermethylation on 4 differential DiffBind marks; bars
# ordered by OR magnitude. H2AK119ub OR = 4.71 is the dominant positive driver.
# =============================================================================

logit <- read_table_tsv("diffbind_logistic_model_coefficients.tsv")

req_b <- c("display_name", "or", "or_lower", "or_upper", "p_value")
if (!all(req_b %in% names(logit))) {
  stop("Panel B: diffbind_logistic_model_coefficients.tsv missing columns: ",
       paste(setdiff(req_b, names(logit)), collapse = ", "))
}

# Order marks by OR (largest at top) for the waterfall reading.
logit <- logit[order(logit$or), ]
logit$display_name <- factor(logit$display_name, levels = logit$display_name)

# Colour positive vs negative association so the K119ub driver (OR > 1) reads
# distinctly from the marks that suppress hypermethylation (OR < 1).
logit$assoc <- ifelse(logit$or > 1, "Promotes hyper (OR > 1)",
                      "Suppresses hyper (OR < 1)")
panelB_cols <- c("Promotes hyper (OR > 1)" = COLORS$k119ub[["K119ub Gained"]],
                 "Suppresses hyper (OR < 1)" = "grey65")

k119_or <- logit$or[logit$display_name == "H2AK119ub"]

p_B <- ggplot(logit,
              aes(x = or, y = display_name, fill = assoc)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(xmin = or_lower, xmax = or_upper),
                width = 0.25, linewidth = 0.4, colour = "black", orientation = "y") +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4,
             colour = "black") +
  geom_text(aes(label = sprintf("%.2f", or)),
            hjust = -0.15, size = 2.4) +
  scale_fill_manual(values = panelB_cols, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  labs(
    title = "H2AK119ub -> hypermethylation",
    subtitle = sprintf("Differential-mark logistic model: %s dominant",
                       stat_label("H2AK119ub OR", k119_or, "%.2f")),
    x = "Odds ratio for hypermethylation",
    y = NULL
  ) +
  theme_pub() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 6))

# =============================================================================
# PANEL C — TET-block evidence : model AUC comparison (Section 24 + 23)
# Source: dnmt3a_model_comparison.tsv (TET impediment vs DNMT3A recruitment)
#         baseline_hmc_model_comparison.tsv (baseline-5hmC predictor of DMR)
#         dnmt3a_exclusive_model_comparison.tsv (DeLong p for the contrast)
# TET impediment AUC 0.793 >> DNMT3A recruitment 0.696; DeLong p = 9.43e-49.
# =============================================================================

dnmt3a <- read_table_tsv("dnmt3a_model_comparison.tsv")
basehmc <- read_table_tsv("baseline_hmc_model_comparison.tsv")
delong_tab <- read_table_tsv("dnmt3a_exclusive_model_comparison.tsv")

req_c1 <- c("model", "auc", "auc_ci_lower", "auc_ci_upper")
if (!all(req_c1 %in% names(dnmt3a))) {
  stop("Panel C: dnmt3a_model_comparison.tsv missing columns: ",
       paste(setdiff(req_c1, names(dnmt3a)), collapse = ", "))
}
if (!all(req_c1 %in% names(basehmc))) {
  stop("Panel C: baseline_hmc_model_comparison.tsv missing columns: ",
       paste(setdiff(req_c1, names(basehmc)), collapse = ", "))
}
if (!"delong_p" %in% names(delong_tab)) {
  stop("Panel C: dnmt3a_exclusive_model_comparison.tsv missing 'delong_p'.")
}

# Pull the two adjudicating models from Section 24 and the baseline-5hmC
# substrate model (Model A) from Section 23 as the supporting TET-substrate
# evidence. Each row is selected by its exact model label from the TSV.
tet_row    <- dnmt3a[dnmt3a$model == "TET impediment", ]
dnmt3a_row <- dnmt3a[dnmt3a$model == "DNMT3A recruitment", ]
base_row   <- basehmc[basehmc$model == "A: Baseline 5hmC", ]
if (nrow(tet_row) != 1 || nrow(dnmt3a_row) != 1 || nrow(base_row) != 1) {
  stop("Panel C: expected exactly one row each for 'TET impediment', ",
       "'DNMT3A recruitment', and 'A: Baseline 5hmC'.")
}

panelC <- data.frame(
  model_label = c("Baseline 5hmC\n(TET substrate)",
                  "TET impediment\n(5hmC + ATAC)",
                  "DNMT3A recruitment\n(K119ub + ATAC + CpG)"),
  auc      = c(base_row$auc,          tet_row$auc,          dnmt3a_row$auc),
  ci_lower = c(base_row$auc_ci_lower, tet_row$auc_ci_lower, dnmt3a_row$auc_ci_lower),
  ci_upper = c(base_row$auc_ci_upper, tet_row$auc_ci_upper, dnmt3a_row$auc_ci_upper),
  family   = c("TET pathway", "TET pathway", "DNMT3A pathway"),
  stringsAsFactors = FALSE
)
panelC$model_label <- factor(panelC$model_label, levels = rev(panelC$model_label))

panelC_cols <- c("TET pathway"    = COLORS$methylation[["5hmC"]],   # 5hmC blue = TET pathway
                 "DNMT3A pathway" = COLORS$methylation[["5mC"]])    # 5mC red  = de-novo path

# DeLong p for the shared-model TET-vs-DNMT3A contrast (TSV-confirmed 9.43e-49).
delong_p <- unique(delong_tab$delong_p[delong_tab$model %in%
                     c("DNMT3A (shared)", "TET (shared)")])
delong_p <- delong_p[1]

p_C <- ggplot(panelC,
              aes(x = auc, y = model_label, fill = family)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper),
                width = 0.25, linewidth = 0.4, colour = "black", orientation = "y") +
  geom_vline(xintercept = 0.5, linetype = "dotted", linewidth = 0.4,
             colour = "grey40") +
  geom_text(aes(label = sprintf("%.3f", auc)),
            hjust = -0.15, size = 2.4) +
  scale_fill_manual(values = panelC_cols, name = NULL) +
  coord_cartesian(xlim = c(0.5, 0.92)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "TET impediment, not DNMT3A recruitment",
    subtitle = sprintf("%s vs %s  (DeLong p = %s)",
                       stat_label("TET AUC", tet_row$auc, "%.3f"),
                       stat_label("DNMT3A AUC", dnmt3a_row$auc, "%.3f"),
                       formatC(delong_p, format = "e", digits = 2)),
    x = "Hypermethylation prediction AUC",
    y = NULL
  ) +
  theme_pub() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 6))

# =============================================================================
# PANEL D — mC -> MeCP2 : R-squared cascade (Section 62)
# Source: 62_model_comparison_summary.tsv
# Chromatin marks explain MeCP2 binding (R2 0.246) an order of magnitude better
# than CG methylation alone (R2 0.017) — MeCP2 reads chromatin, not methylation.
#
# IMPORTANT: section_62 writes this TSV with quote = FALSE while its `label`
# column embeds a literal newline (label = paste0(type, "\n", model)). Each of
# the 6 real records therefore spans 2 physical lines: a full record line plus a
# ragged 1-field spillover line (the second half of `label`). The generic
# read_table_tsv() uses read.table() defaults that error on ragged rows, so this
# one malformed file is read directly with fill = TRUE (same lossless options
# otherwise). We then select the 6 real records by their exact `type` x `model`
# identity; the spillover rows carry NA r_squared and a `type` that is not a
# valid model level, so they are excluded. No data fabricated.
# =============================================================================

mecp2_r2_path <- file.path(TABLES_DIR, "62_model_comparison_summary.tsv")
if (!file.exists(mecp2_r2_path)) {
  stop("Panel D: table not found: ", mecp2_r2_path)
}
mecp2_r2_raw <- utils::read.table(
  mecp2_r2_path,
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE,
  quote            = "",
  check.names      = FALSE,
  comment.char     = "",
  fill             = TRUE
)

req_d <- c("type", "model", "r_squared")
if (!all(req_d %in% names(mecp2_r2_raw))) {
  stop("Panel D: 62_model_comparison_summary.tsv missing columns: ",
       paste(setdiff(req_d, names(mecp2_r2_raw)), collapse = ", "))
}

valid_types  <- c("Binding level", "Differential")
valid_models <- c("CG only", "Chromatin only", "Full (CG + Chromatin)")

mecp2_r2 <- mecp2_r2_raw[
  mecp2_r2_raw$type  %in% valid_types &
  mecp2_r2_raw$model %in% valid_models &
  !is.na(suppressWarnings(as.numeric(mecp2_r2_raw$r_squared))), ]
mecp2_r2$r_squared <- as.numeric(mecp2_r2$r_squared)

if (nrow(mecp2_r2) != 6) {
  stop("Panel D: expected 6 valid model rows after de-spilling, got ",
       nrow(mecp2_r2), ".")
}

mecp2_r2$type  <- factor(mecp2_r2$type, levels = valid_types)
mecp2_r2$model <- factor(mecp2_r2$model, levels = valid_models)

panelD_cols <- c(
  "CG only"               = COLORS$methylation[["5mC"]],                      # methylation red
  "Chromatin only"        = COLORS$k119ub[["K119ub Gained"]],  # chromatin purple
  "Full (CG + Chromatin)" = "grey45"
)

# Headline binding-level CG vs Chromatin R2 (TSV-confirmed 0.0168 vs 0.2456).
cg_bind   <- mecp2_r2$r_squared[mecp2_r2$type == "Binding level" &
                                 mecp2_r2$model == "CG only"]
chrom_bind <- mecp2_r2$r_squared[mecp2_r2$type == "Binding level" &
                                  mecp2_r2$model == "Chromatin only"]

p_D <- ggplot(mecp2_r2,
              aes(x = type, y = r_squared, fill = model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.3f", r_squared)),
            position = position_dodge(width = 0.8),
            vjust = -0.4, size = 2.2) +
  scale_fill_manual(values = panelD_cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "DNA methylation -> MeCP2 binding",
    subtitle = sprintf("MeCP2 reads chromatin: %s vs %s (binding level)",
                       stat_label("chromatin R²", chrom_bind, "%.3f"),
                       stat_label("CG R²", cg_bind, "%.3f")),
    x = NULL,
    y = expression(paste("Variance explained (", R^2, ")"))
  ) +
  theme_pub() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 6))

# =============================================================================
# PANEL E — Neuronal convergence : OR summary (Sections 74 + 72)
# Source: 74_pairwise_fisher.tsv (gene-set overlaps)
#         72_fisher_results.tsv (K119ub top-decile neuronal enrichment)
# Coordinated x MeCP2-Up OR = 5.16 (methylation tracks redistribution) >>
# Neuronal x MeCP2-Up OR = 1.73 (lineage); constitutive K119ub top-decile
# neuronal-enrichment OR = 1.70.
# =============================================================================

pw_fisher <- read_table_tsv("74_pairwise_fisher.tsv")
neu_fisher <- read_table_tsv("72_fisher_results.tsv")

req_e1 <- c("comparison", "odds_ratio", "ci_lo", "ci_hi", "p_adj")
if (!all(req_e1 %in% names(pw_fisher))) {
  stop("Panel E: 74_pairwise_fisher.tsv missing columns: ",
       paste(setdiff(req_e1, names(pw_fisher)), collapse = ", "))
}
req_e2 <- c("test", "OR", "ci_lo", "ci_hi", "q_value")
if (!all(req_e2 %in% names(neu_fisher))) {
  stop("Panel E: 72_fisher_results.tsv missing columns: ",
       paste(setdiff(req_e2, names(neu_fisher)), collapse = ", "))
}

# Select the three convergence ORs by their exact source-row identity.
coord_mecp2 <- pw_fisher[pw_fisher$comparison == "Coordinated × MeCP2-Up", ]
if (nrow(coord_mecp2) == 0) {  # tolerate ASCII 'x' if the TSV used it
  coord_mecp2 <- pw_fisher[grepl("Coordinated.*MeCP2-Up", pw_fisher$comparison), ]
}
neu_mecp2   <- pw_fisher[grepl("Neuronal.*MeCP2-Up", pw_fisher$comparison), ]
k119_decile <- neu_fisher[neu_fisher$test == "Ctrl top decile (D9)", ]

if (nrow(coord_mecp2) != 1 || nrow(neu_mecp2) != 1 || nrow(k119_decile) != 1) {
  stop("Panel E: could not uniquely select the three convergence rows ",
       "(Coordinated x MeCP2-Up, Neuronal x MeCP2-Up, Ctrl top decile D9).")
}

panelE <- data.frame(
  label = c("Coordinated mC-up/hmC-down\nx MeCP2-Up",
            "Broad neuronal\n× MeCP2-Up",
            "Constitutive K119ub\ntop decile × neuronal"),
  or    = c(coord_mecp2$odds_ratio, neu_mecp2$odds_ratio, k119_decile$OR),
  ci_lo = c(coord_mecp2$ci_lo,      neu_mecp2$ci_lo,      k119_decile$ci_lo),
  ci_hi = c(coord_mecp2$ci_hi,      neu_mecp2$ci_hi,      k119_decile$ci_hi),
  layer = c("Methylation-driven", "Lineage", "Lineage"),
  padj  = c(coord_mecp2$p_adj,      neu_mecp2$p_adj,      k119_decile$q_value),
  stringsAsFactors = FALSE
)
panelE$label <- factor(panelE$label, levels = rev(panelE$label))

panelE_cols <- c("Methylation-driven" = COLORS$mecp2[["MeCP2 Up"]],  # "#D95F02"
                 "Lineage"            = "grey60")

coord_or <- coord_mecp2$odds_ratio
neu_or   <- neu_mecp2$odds_ratio

p_E <- ggplot(panelE,
              aes(x = or, y = label, colour = layer)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4,
             colour = "black") +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
                width = 0.2, linewidth = 0.5, orientation = "y") +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("OR = %.2f", or)),
            hjust = -0.25, vjust = -0.7, size = 2.3, show.legend = FALSE) +
  scale_colour_manual(values = panelE_cols, name = NULL) +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.25))) +
  labs(
    title = "Convergence on neuronal genes",
    subtitle = sprintf("MeCP2 tracks methylation (%s) > lineage (%s)",
                       stat_label("OR", coord_or, "%.2f"),
                       stat_label("OR", neu_or, "%.2f")),
    x = "Odds ratio (log scale)",
    y = NULL
  ) +
  theme_pub() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 6))

# =============================================================================
# COMPOSE — patchwork assembly per plan layout
#   design  = "AABB\nCCDD\n#EE#"
#   heights = c(1, 1, 0.6)
# Panel tags A..E come from tag_levels = "A" (same call add_panel_labels()
# wraps); the title/subtitle are attached on the same plot_annotation() so the
# cascade reads as one argument. Tag appearance is inherited from theme_pub().
# =============================================================================

fig7_design <- "AABB\nCCDD\n#EE#"

fig7 <- (p_A + p_B + p_C + p_D + p_E) +
  plot_layout(design = fig7_design, heights = c(1, 1, 0.6))

fig7 <- fig7 +
  patchwork::plot_annotation(
    tag_levels = "A",
    title    = "Section 87: Quantitative cascade: BAP1 loss -> H2AK119ub -> TET-block -> methylation -> MeCP2 -> neuronal genes",
    subtitle = "Each panel quantifies one mechanistic link (deep-seq, 8 samples, sex covariate)",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 11, hjust = 0),
      plot.subtitle = element_text(size = 8, colour = "grey40", hjust = 0)
    )
  )

# =============================================================================
# SAVE — full-width short composite (plan: 180 x 140 mm) + the load-bearing
# sub-panels (B = OR waterfall, D = R2 cascade) for flexible figure assembly.
# save_figure() handles mm -> inches and emits PDF + SVG + JPG into a subfolder.
# =============================================================================

save_plot(p_A, "87a_k119ub_effect_size", w = 8, h = 6)
save_plot(p_B, "87b_k119ub_or_waterfall", w = 8, h = 6)
save_plot(p_C, "87c_tet_dnmt3a_auc", w = 8, h = 6)
save_plot(p_D, "87d_mecp2_r2_cascade", w = 8, h = 6)
save_plot(p_E, "87e_neuronal_convergence", w = 8, h = 6)
save_plot(fig7, "87_composite", w = 14, h = 10)

cat("\nSection 87 complete. Outputs in ", SEC87_DIR, "\n")
cat("Key stats (TSV-confirmed):\n")
cat(sprintf("  A  background K119ub median log2FC = %+.3f (p = %s)\n",
            bg_median, formatC(bg_p, format = "e", digits = 2)))
cat(sprintf("  B  H2AK119ub OR for hypermethylation = %.2f\n", k119_or))
cat(sprintf("  C  TET AUC = %.3f vs DNMT3A AUC = %.3f (DeLong p = %s)\n",
            tet_row$auc, dnmt3a_row$auc,
            formatC(delong_p, format = "e", digits = 2)))
cat(sprintf("  D  MeCP2 binding R2: chromatin = %.3f vs CG = %.3f\n",
            chrom_bind, cg_bind))
cat(sprintf("  E  Coordinated x MeCP2-Up OR = %.2f; K119ub top-decile OR = %.2f\n",
            coord_or, k119_decile$OR))
