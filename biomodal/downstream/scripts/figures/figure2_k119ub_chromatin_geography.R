# scripts/figures/figure2_k119ub_chromatin_geography.R
# =============================================================================
# FIGURE 2: H2AK119ub Drives Hypermethylation at Active Chromatin (paper R2)
# =============================================================================
# Question: Why does methylation change in BAP1-KO? H2AK119ub is the upstream
# driver, and its effect is restricted to active/euchromatic chromatin while
# Polycomb heterochromatin is protected.
#
# This script LEADS with the data the reader can see directly (the chromatin-
# state direction split, the compartment geography, the Polycomb exclusion) and
# uses the multi-mark logistic model as quantitative confirmation, exactly per
# the master plan and the TODOS framing decision.
#
# Panels (master plan "Figure 2", layout AABB / CCDD / #EE#):
#   A  Sec 10  Stacked fraction bar: per-state hyper/hypo composition
#              (chromatin_state_summary.tsv). Active_Promoter 93.0% hyper.
#   B  Sec 33  Forest plot: 4-mark logistic ORs predicting mC hypermethylation
#              (diffbind_logistic_model_coefficients.tsv). H2AK119ub OR=4.71.
#   C  Sec 29  Compartment enrichment OR bars, log x
#              (compartment_fisher_tests.tsv). mC-hyper -> A OR=13.64.
#   D  Sec 30  Polycomb exclusion / enrichment bars
#              (polycomb_fisher_tests.tsv). Polycomb x hyper OR=0.063.
#   E  Sec 17  K119ub "honest" breakdown as data: ~58% of DMR genes carry no
#              K119ub peak; among peak-bearing genes mC-Up gains +14.2 pp over
#              background (k119ub_honest_breakdown.tsv).
#
# Output: plots/figures/figure2_k119ub_chromatin_geography/ (PDF + SVG + JPG)
# Dimensions: 180 x 200 mm (per master plan).
#
# All H2AK119ub references denote the BAP1 PR-DUB substrate (histone mark).
# =============================================================================

# -----------------------------------------------------------------------------
# Shared infrastructure (must be sourced first; working dir == downstream/).
# -----------------------------------------------------------------------------
source("scripts/viz_sections/_shared_config.R")
source("scripts/figures/_figure_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(patchwork)
})

FIG2_NAME <- "figure2_k119ub_chromatin_geography"

# =============================================================================
# PANEL A: Chromatin-state direction split (Section 10)
# =============================================================================
# chromatin_state_summary.tsv columns:
#   chromatin_state, count_Hypermethylated, count_Hypomethylated, mean_/median_*
# We show, per state, the fraction of significant DMRs that are hyper vs hypo.
# Active_Promoter is 93.0% hyper (4562/4906); Repressed_Promoter 94.4% hypo.
build_panel_A <- function() {
  cs <- read_table_tsv("chromatin_state_summary.tsv")

  required <- c("chromatin_state", "count_Hypermethylated", "count_Hypomethylated")
  missing <- setdiff(required, colnames(cs))
  if (length(missing) > 0) {
    stop("Panel A: chromatin_state_summary.tsv missing columns: ",
         paste(missing, collapse = ", "))
  }

  cs_long <- cs %>%
    transmute(
      chromatin_state,
      Hypermethylated = count_Hypermethylated,
      Hypomethylated  = count_Hypomethylated,
      total = count_Hypermethylated + count_Hypomethylated
    ) %>%
    pivot_longer(
      cols = c("Hypermethylated", "Hypomethylated"),
      names_to = "direction", values_to = "count"
    ) %>%
    mutate(
      direction = factor(direction,
                         levels = c("Hypermethylated", "Hypomethylated")),
      chromatin_state = factor(chromatin_state,
                               levels = rev(CHROMATIN_STATE_ORDER))
    )

  # Per-state totals for the n= annotation at the right edge of each bar.
  state_totals <- cs %>%
    transmute(
      chromatin_state = factor(chromatin_state,
                               levels = rev(CHROMATIN_STATE_ORDER)),
      total = count_Hypermethylated + count_Hypomethylated,
      label = paste0("n=", scales::comma(total))
    )

  # Verified key stat: Active_Promoter hyper fraction = 4562/4906 = 93.0%.
  ap_hyper_pct <- with(cs[cs$chromatin_state == "Active_Promoter", ],
                       100 * count_Hypermethylated /
                         (count_Hypermethylated + count_Hypomethylated))

  ggplot(cs_long, aes(x = count, y = chromatin_state, fill = direction)) +
    geom_col(position = "fill", width = 0.78) +
    geom_text(
      data = state_totals,
      aes(x = 1.02, y = chromatin_state, label = label),
      inherit.aes = FALSE, hjust = 0, size = 2.1, colour = "grey25"
    ) +
    scale_fill_manual(values = COLORS$direction, name = NULL) +
    scale_x_continuous(
      labels = scales::percent, expand = expansion(mult = c(0, 0.22)),
      breaks = c(0, 0.5, 1)
    ) +
    labs(
      title = "DMR direction by chromatin state",
      subtitle = sprintf("Active_Promoter %.0f%% hypermethylated", ap_hyper_pct),
      x = "Fraction of significant mC DMRs", y = NULL
    ) +
    theme_pub() +
    theme(legend.position = "top")
}

# =============================================================================
# PANEL B: 4-mark logistic forest plot (Section 33)
# =============================================================================
# diffbind_logistic_model_coefficients.tsv columns:
#   term, estimate, or, or_lower, or_upper, p_value, display_name, sig_label
# H2AK119ub OR=4.707; H3K27me3 1.443; H3K27ac 0.480; ATAC-seq 0.221.
build_panel_B <- function() {
  fm <- read_table_tsv("diffbind_logistic_model_coefficients.tsv")

  required <- c("display_name", "or", "or_lower", "or_upper", "p_value", "sig_label")
  missing <- setdiff(required, colnames(fm))
  if (length(missing) > 0) {
    stop("Panel B: diffbind_logistic_model_coefficients.tsv missing columns: ",
         paste(missing, collapse = ", "))
  }

  fm <- fm %>%
    mutate(
      mark_dir = ifelse(or >= 1, "Promotes hypermethylation", "Protective"),
      display_name = forcats::fct_reorder(display_name, or)
    )

  # Verified key stat: H2AK119ub OR = 4.71 (4.70747).
  k119_or <- fm$or[fm$display_name == "H2AK119ub"]
  k119_lab <- sprintf("H2AK119ub %s", stat_label("OR", k119_or, "%.2f"))

  ggplot(fm, aes(x = or, y = display_name, colour = mark_dir)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = or_lower, xmax = or_upper),
                   height = 0.22, linewidth = 0.5) +
    geom_point(size = 2.2) +
    geom_text(aes(label = sprintf("%.2f", or)),
              vjust = -1.0, size = 2.1, show.legend = FALSE) +
    scale_x_log10(
      breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10),
      labels = c("0.1", "0.25", "0.5", "1", "2", "5", "10")
    ) +
    scale_colour_manual(
      values = c("Promotes hypermethylation" = COLORS$direction[["Hypermethylated"]],
                 "Protective" = COLORS$direction[["Hypomethylated"]]),
      name = NULL
    ) +
    annotate("text", x = min(fm$or_lower), y = 0.6, hjust = 0,
             label = k119_lab, size = 2.2, fontface = "bold",
             colour = COLORS$direction[["Hypermethylated"]]) +
    labs(
      title = "Mark log2FC predicts mC hypermethylation",
      subtitle = "4-mark logistic model (DiffBind), AUC = 0.818",
      x = "Odds ratio per unit log2FC (log scale)", y = NULL
    ) +
    theme_pub() +
    theme(legend.position = "top")
}

# =============================================================================
# PANEL C: A/B compartment enrichment OR bars (Section 29)
# =============================================================================
# compartment_fisher_tests.tsv columns:
#   test, positive_set_n, expected_compartment, a,b,c,d,
#   odds_ratio, ci_lower, ci_upper, p_value, pct_in_expected, pct_in_other
# We show the four static directional tests (not the dynamic shift tests, which
# belong conceptually with panel D's exclusion narrative). mC-hyper -> A OR=13.64.
build_panel_C <- function() {
  cf <- read_table_tsv("compartment_fisher_tests.tsv")

  required <- c("test", "odds_ratio", "ci_lower", "ci_upper", "p_value")
  missing <- setdiff(required, colnames(cf))
  if (length(missing) > 0) {
    stop("Panel C: compartment_fisher_tests.tsv missing columns: ",
         paste(missing, collapse = ", "))
  }

  static_tests <- c(
    "mC hyper -> A enriched",
    "hmC hypo -> A enriched",
    "mC hypo -> B enriched",
    "hmC hyper -> B enriched"
  )
  cf_use <- cf %>% filter(test %in% static_tests)
  if (nrow(cf_use) != length(static_tests)) {
    stop("Panel C: expected ", length(static_tests),
         " static compartment tests, found ", nrow(cf_use),
         " in compartment_fisher_tests.tsv")
  }

  pretty_labels <- c(
    "mC hyper -> A enriched"  = "mC hyper → A",
    "hmC hypo -> A enriched"  = "hmC loss → A",
    "mC hypo -> B enriched"   = "mC hypo → B",
    "hmC hyper -> B enriched" = "hmC gain → B"
  )

  cf_use <- cf_use %>%
    mutate(
      label = pretty_labels[test],
      target = ifelse(grepl("A enriched$", test), "A (euchromatin)",
                      "B (heterochromatin)"),
      label = forcats::fct_reorder(label, odds_ratio)
    )

  # Verified key stat: mC hyper -> A OR = 13.64 (13.6423).
  a_or <- cf_use$odds_ratio[cf_use$test == "mC hyper -> A enriched"]

  ggplot(cf_use, aes(x = odds_ratio, y = label, fill = target)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_col(width = 0.66) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                   height = 0.2, linewidth = 0.4, colour = "grey25") +
    geom_text(aes(label = sprintf("OR=%.1f", odds_ratio)),
              hjust = -0.15, size = 2.1) +
    scale_x_log10(
      breaks = c(1, 2, 5, 10, 20),
      labels = c("1", "2", "5", "10", "20"),
      expand = expansion(mult = c(0, 0.18))
    ) +
    scale_fill_manual(
      values = c("A (euchromatin)" = COLORS$direction[["Hypermethylated"]],
                 "B (heterochromatin)" = COLORS$direction[["Hypomethylated"]]),
      name = NULL
    ) +
    labs(
      title = "Methylation geography in A/B compartments",
      subtitle = sprintf("Hypermethylation enriched in A (OR = %.1f)", a_or),
      x = "Fisher odds ratio (log scale)", y = NULL
    ) +
    theme_pub() +
    theme(legend.position = "top")
}

# =============================================================================
# PANEL D: Polycomb exclusion (Section 30)
# =============================================================================
# polycomb_fisher_tests.tsv columns:
#   test, n_universe, n_polycomb, n_dmr, a,b,c,d, odds_ratio, ci_lower,
#   ci_upper, p_value, pct_polycomb_dmr, pct_non_polycomb_dmr, q_value
# We show the chromatin-state Polycomb tests for the four methylation
# directions: hyper excluded (OR=0.063), hypo enriched (OR=9.80), plus the
# mirrored hmC pair. This is the decisive falsification of the naive
# "Polycomb gets hypermethylated" prediction.
build_panel_D <- function() {
  pf <- read_table_tsv("polycomb_fisher_tests.tsv")

  required <- c("test", "odds_ratio", "ci_lower", "ci_upper", "p_value")
  missing <- setdiff(required, colnames(pf))
  if (length(missing) > 0) {
    stop("Panel D: polycomb_fisher_tests.tsv missing columns: ",
         paste(missing, collapse = ", "))
  }

  state_tests <- c(
    "Chromatin State × mC hyper",
    "Chromatin State × mC hypo",
    "Chromatin State × hmC hyper",
    "Chromatin State × hmC hypo"
  )
  pf_use <- pf %>% filter(test %in% state_tests)
  if (nrow(pf_use) != length(state_tests)) {
    stop("Panel D: expected ", length(state_tests),
         " chromatin-state Polycomb tests, found ", nrow(pf_use),
         " in polycomb_fisher_tests.tsv")
  }

  pretty_labels <- c(
    "Chromatin State × mC hyper"  = "mC hyper",
    "Chromatin State × mC hypo"   = "mC hypo",
    "Chromatin State × hmC hyper" = "hmC gain",
    "Chromatin State × hmC hypo"  = "hmC loss"
  )

  pf_use <- pf_use %>%
    mutate(
      label = pretty_labels[test],
      effect = ifelse(odds_ratio >= 1, "Enriched at Polycomb",
                      "Excluded from Polycomb"),
      label = forcats::fct_reorder(label, odds_ratio)
    )

  # Verified key stat: Polycomb x mC hyper OR = 0.063 (0.0633).
  pc_hyper_or <- pf_use$odds_ratio[pf_use$test == "Chromatin State × mC hyper"]

  ggplot(pf_use, aes(x = odds_ratio, y = label, fill = effect)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_col(width = 0.66) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                   height = 0.2, linewidth = 0.4, colour = "grey25") +
    geom_text(aes(label = ifelse(odds_ratio < 1,
                                 sprintf("OR=%.3f", odds_ratio),
                                 sprintf("OR=%.1f", odds_ratio)),
                  hjust = ifelse(odds_ratio < 1, 1.1, -0.15)),
              size = 2.1) +
    scale_x_log10(
      breaks = c(0.05, 0.1, 0.5, 1, 5, 10),
      labels = c("0.05", "0.1", "0.5", "1", "5", "10"),
      expand = expansion(mult = c(0.20, 0.18))
    ) +
    scale_fill_manual(
      values = c("Enriched at Polycomb" = COLORS$direction[["Hypomethylated"]],
                 "Excluded from Polycomb" = COLORS$direction[["Hypermethylated"]]),
      name = NULL
    ) +
    labs(
      title = "Polycomb chromatin is protected",
      subtitle = sprintf("Hypermethylation excluded from Polycomb (OR = %.3f)",
                         pc_hyper_or),
      x = "Fisher odds ratio (log scale)", y = NULL
    ) +
    theme_pub() +
    theme(legend.position = "top")
}

# =============================================================================
# PANEL E: K119ub "honest" breakdown as data (Section 17)
# =============================================================================
# k119ub_honest_breakdown.tsv columns:
#   met_group, total_genes, genes_with_peaks, genes_no_peaks, pct_no_peaks,
#   genes_gained, genes_equal, genes_lost, pct_gained_of_total, pct_gained_of_peaks
# Two facts shown as data, not as an apology:
#   (1) ~58% of DMR genes carry no K119ub peak at all (pct_no_peaks).
#   (2) Among peak-bearing genes, mC-Up gains K119ub at 47.8% vs 33.6%
#       background = +14.2 pp (pct_gained_of_peaks).
build_panel_E <- function() {
  kb <- read_table_tsv("k119ub_honest_breakdown.tsv")

  required <- c("met_group", "pct_no_peaks", "pct_gained_of_peaks")
  missing <- setdiff(required, colnames(kb))
  if (length(missing) > 0) {
    stop("Panel E: k119ub_honest_breakdown.tsv missing columns: ",
         paste(missing, collapse = ", "))
  }

  group_order <- c("mC Up", "mC Down", "hmC Down", "hmC Up", "All DMR Genes")
  kb_use <- kb %>% filter(met_group %in% group_order)
  if (!all(group_order %in% kb_use$met_group)) {
    stop("Panel E: k119ub_honest_breakdown.tsv missing expected met_group rows.")
  }

  # Verified key stats: All DMR Genes pct_no_peaks = 58.16 (=58.2%);
  # background pct_gained_of_peaks = 33.626; mC Up = 47.834 -> +14.2 pp.
  bg_gain  <- kb_use$pct_gained_of_peaks[kb_use$met_group == "All DMR Genes"]
  mcup_gain <- kb_use$pct_gained_of_peaks[kb_use$met_group == "mC Up"]
  no_peak_all <- kb_use$pct_no_peaks[kb_use$met_group == "All DMR Genes"]
  delta_pp <- mcup_gain - bg_gain

  kb_long <- kb_use %>%
    transmute(
      met_group = factor(met_group, levels = group_order),
      `% with no K119ub peak`            = pct_no_peaks,
      `% K119ub gained (peak-bearing)`   = pct_gained_of_peaks
    ) %>%
    pivot_longer(cols = -met_group, names_to = "metric", values_to = "pct") %>%
    mutate(metric = factor(metric,
                           levels = c("% with no K119ub peak",
                                      "% K119ub gained (peak-bearing)")))

  ggplot(kb_long, aes(x = met_group, y = pct, fill = metric)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.68) +
    geom_hline(yintercept = bg_gain, linetype = "dashed",
               colour = COLORS$k119ub[["K119ub Gained"]], linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.0f", pct)),
              position = position_dodge(width = 0.75),
              vjust = -0.4, size = 1.9) +
    annotate("text",
             x = 1, y = mcup_gain + 14,
             label = sprintf("mC Up gain +%.1f pp\nover background (%.0f%%)",
                             delta_pp, bg_gain),
             size = 2.0, hjust = 0.5, colour = "grey20") +
    scale_fill_manual(
      values = c("% with no K119ub peak" = "grey60",
                 "% K119ub gained (peak-bearing)" = COLORS$k119ub[["K119ub Gained"]]),
      name = NULL
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.20)),
                       limits = c(0, 100)) +
    labs(
      title = "K119ub coupling: real but partial",
      subtitle = sprintf("%.0f%% of DMR genes carry no K119ub peak",
                         no_peak_all),
      x = NULL, y = "Percent of genes"
    ) +
    theme_pub() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 25, hjust = 1))
}

# =============================================================================
# COMPOSE + SAVE
# =============================================================================
message("Building Figure 2 panels ...")
pA <- build_panel_A()
pB <- build_panel_B()
pC <- build_panel_C()
pD <- build_panel_D()
pE <- build_panel_E()

# Master-plan layout: AABB / CCDD / #EE# with heights 1.2, 1, 0.7.
fig2_design <- "AABB\nCCDD\n#EE#"

fig2 <- (pA + pB + pC + pD + pE) +
  plot_layout(
    design  = fig2_design,
    heights = c(1.2, 1, 0.7)
  ) +
  plot_annotation(
    tag_levels = "A",
    title = "Figure 2  |  H2AK119ub drives hypermethylation at active chromatin",
    theme = theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0)
    )
  )

message("Saving Figure 2 (PDF + SVG + JPG) to ", file.path(FIGURE_DIR, FIG2_NAME))
save_figure(fig2, FIG2_NAME, width_mm = 180, height_mm = 200)

message("Figure 2 complete.")
