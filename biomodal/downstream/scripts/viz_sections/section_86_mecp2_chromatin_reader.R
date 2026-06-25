# scripts/viz_sections/section_86_mecp2_chromatin_reader.R
# =============================================================================
# Figure 5: MeCP2 Reads Chromatin State, Not Methylation  (paper section R4)
#
# Central claim: in BAP1-KO cerebellum, MeCP2 binding change tracks chromatin
# state -- above all the Polycomb mark H2AK119ub -- and NOT CG DNA methylation.
# CG methylation is necessary but not sufficient: most hypermethylated genes do
# not recruit more MeCP2, and the genes that gain MeCP2 without any methylation
# change do so because K119ub accumulates there. MeCP2 is a chromatin reader.
#
# All MeCP2 data is CUT&RUN. Labels use "binding"/"signal"/"mark" throughout.
#
# Panels (per approved master plan, Figure 5 table + layout):
#   A  Sec 62  R^2 comparison bar       62_model_comparison_summary.tsv
#             CG-only vs Chromatin-only vs Full, for Binding-level and
#             Differential models. CG R2=0.017 << Chromatin R2=0.246 (binding).
#   B  Sec 62  Variance partition       62_variance_partition.tsv
#             Stacked CG-unique / Shared / Chromatin-unique / Unexplained, for
#             Binding-level and Differential. Chromatin-unique 24.3% >> CG 1.5%.
#   C  Sec 67  K119ub at MeCP2 loci     59_quadrant_master.tsv (+ 67_statistics.tsv)
#             Per-gene gene-body K119ub log2FC by MeCP2 status (Up/Down/NS).
#             MeCP2-Up loci gain K119ub; Fisher OR=3.15 (MeCP2-up-no-meth genes).
#   D  Sec 60  Mirror-image profiles    60_mecp2_status_stats.tsv (wide -> long)
#             Per-mark mean change for MeCP2-Up vs Down vs NS. The FLAT K119ub
#             at MeCP2-Down (mean -0.04, ns) is the visual punchline.
#   E  Sec 71  K119ub vs ratio variance 71_variance_partition.tsv
#             Stacked K119ub-unique / Shared / Ratio-unique / Unexplained.
#             K119ub 7.3% vs methylation-ratio 0.0% of MeCP2-fold variance.
#
# Layout: design = "AABB\nCCDD\nEEEE"; heights c(1, 1, 0.8); 180 x 200 mm.
#
# Data source: pre-computed section TSVs in plots/visualizations/tables/.
# No new computation beyond in-script reshaping (60 wide->long) and the
# documented MeCP2-status partition of 59_quadrant_master (Up/Down/NS).
# =============================================================================

# -----------------------------------------------------------------------------
# Shared infrastructure. Working directory MUST be downstream/ (BASE_DIR=getwd()).
# _shared_config.R provides COLORS, theme_biomodal(), mc_dmr/hmc_dmr, helpers,
# and (via multi_format_output.R) the save_multiformat_* functions.
# _figure_config.R provides theme_pub(), save_figure(), read_table_tsv(),
# add_panel_labels(), stat_label(), FIGURE_DIR, FIGURE_TABLES_DIR, PUB_COLORS.
# -----------------------------------------------------------------------------
source("scripts/viz_sections/_shared_config.R")
source("scripts/viz_sections/_figure_config.R")

SEC86_DIR <- file.path(OUTPUT_DIR, "86_mecp2_chromatin_reader")
dir.create(SEC86_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, w = 10, h = 7) {
  save_multiformat_ggplot(p,
    base_path = file.path(SEC86_DIR, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(patchwork)
})

cat("\n=== Figure 5: MeCP2 Reads Chromatin State, Not Methylation ===\n")

# -----------------------------------------------------------------------------
# Canonical colours used across this figure (verbatim from the COLORS contract).
# K119ub gets its own load-bearing accent so it reads as the dominant signal in
# panels D/E; CG methylation marks use the methylation palette (5mC red /
# 5hmC blue). MeCP2 status uses the mecp2 palette (Up orange / Down purple).
# -----------------------------------------------------------------------------
COL_K119UB <- COLORS$k119ub[["K119ub Gained"]]   # "#756BB1"
COL_5MC    <- COLORS$methylation[["5mC"]]         # "#E41A1C"
COL_5HMC   <- COLORS$methylation[["5hmC"]]        # "#377EB8"
COL_K27AC  <- COLORS$h3k27ac[["H3K27ac Gained"]]  # "#FF7F00"
COL_K27ME3 <- "#1B7837"                           # repressive green (distinct)

MECP2_STATUS_COLORS <- c(
  "MeCP2 Up"        = COLORS$mecp2[["MeCP2 Up"]],        # "#D95F02"
  "MeCP2 Down"      = COLORS$mecp2[["MeCP2 Down"]],      # "#7570B3"
  "Not Significant" = COLORS$mecp2[["Not Significant"]]  # "grey70"
)

# Stacked-variance component palette (shared shape across panels B and E).
VARIANCE_COLORS <- c(
  "Chromatin unique" = COL_K119UB,
  "K119ub unique"    = COL_K119UB,
  "CG unique"        = COL_5MC,
  "Ratio unique"     = COL_5MC,
  "Shared"           = "#999999",
  "Unexplained"      = "grey85"
)

# =============================================================================
# PANEL A (Sec 62): model R^2 comparison -- CG vs Chromatin vs Full.
# Reads 62_model_comparison_summary.tsv. Key stat: binding-level CG R2=0.017 vs
# Chromatin R2=0.246 -- chromatin explains MeCP2 binding ~15x better than CG.
# =============================================================================
build_panel_a <- function() {
  # NOTE: 62_model_comparison_summary.tsv carries a trailing `label` column that
  # section_62 builds as paste0(type, "\n", model) -- i.e. it embeds a literal
  # newline. Written with quote = FALSE, every logical record is split across two
  # physical lines (10 fields then 1 field), so the shared read_table_tsv()
  # (read.table with fill = FALSE) aborts with "line 2 did not have 10 elements".
  # Read with fill = TRUE so the ragged continuation lines are absorbed; the
  # filter(model %in% names(model_map)) below then drops them cleanly (they carry
  # an empty `model`). Only type/model/r_squared are used downstream.
  tbl_path <- file.path(TABLES_DIR, "62_model_comparison_summary.tsv")
  if (!file.exists(tbl_path)) {
    stop("Panel A: table not found: ", tbl_path)
  }
  tbl <- utils::read.table(
    tbl_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
    quote = "", check.names = FALSE, comment.char = "", fill = TRUE
  )

  required <- c("type", "model", "r_squared")
  missing <- setdiff(required, names(tbl))
  if (length(missing) > 0) {
    stop("Panel A: 62_model_comparison_summary.tsv missing columns: ",
         paste(missing, collapse = ", "))
  }

  # Compact, ordered factor labels for model family and analysis level.
  model_map <- c(
    "CG only"               = "CG only",
    "Chromatin only"        = "Chromatin only",
    "Full (CG + Chromatin)" = "Full"
  )
  df <- tbl %>%
    filter(.data$model %in% names(model_map)) %>%
    mutate(
      model_lab = factor(model_map[.data$model],
                         levels = c("CG only", "Chromatin only", "Full")),
      type      = factor(.data$type,
                         levels = c("Binding level", "Differential")),
      r_squared = as.numeric(.data$r_squared)
    )
  if (nrow(df) != 6L) {
    stop("Panel A: expected 6 model rows (2 levels x 3 models), got ", nrow(df))
  }

  # Confirm the load-bearing key stats appear (binding-level CG vs Chromatin).
  bind_cg   <- df$r_squared[df$type == "Binding level" & df$model_lab == "CG only"]
  bind_chr  <- df$r_squared[df$type == "Binding level" & df$model_lab == "Chromatin only"]
  ann_a <- sprintf("Binding level:\nCG %s vs Chromatin %s",
                   stat_label("R²", bind_cg, "%.3f"),
                   stat_label("R²", bind_chr, "%.3f"))

  model_fill <- c("CG only" = COL_5MC, "Chromatin only" = COL_K119UB, "Full" = "#542788")

  ggplot(df, aes(x = .data$type, y = .data$r_squared, fill = .data$model_lab)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75,
             colour = "black", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.3f", .data$r_squared)),
              position = position_dodge(width = 0.8),
              vjust = -0.4, size = 2.1) +
    scale_fill_manual(values = model_fill, name = "Predictors") +
    scale_y_continuous(limits = c(0, 0.30), expand = expansion(mult = c(0, 0.08))) +
    annotate("text", x = 1.5, y = 0.285, label = ann_a,
             size = 2.2, fontface = "italic", hjust = 0.5, lineheight = 0.95) +
    labs(
      title = "Chromatin marks explain MeCP2 binding, CG methylation does not",
      x = NULL,
      y = expression("Model "*R^2)
    ) +
    theme_pub() +
    theme(legend.position = "top",
          plot.title = element_text(size = 8))
}

# =============================================================================
# PANEL B (Sec 62): variance partition (stacked) for binding + differential.
# Reads 62_variance_partition.tsv. Key stat: Chromatin-unique 24.3% >> CG 1.5%.
# =============================================================================
build_panel_b <- function() {
  tbl <- read_table_tsv("62_variance_partition.tsv")

  required <- c("type", "component", "fraction")
  missing <- setdiff(required, names(tbl))
  if (length(missing) > 0) {
    stop("Panel B: 62_variance_partition.tsv missing columns: ",
         paste(missing, collapse = ", "))
  }

  comp_levels <- c("CG unique", "Shared", "Chromatin unique", "Unexplained")
  df <- tbl %>%
    mutate(
      fraction  = as.numeric(.data$fraction),
      component = factor(.data$component, levels = comp_levels),
      type      = factor(.data$type, levels = c("Binding level", "Differential"))
    )
  if (any(is.na(df$component))) {
    stop("Panel B: unexpected variance component label(s): ",
         paste(unique(tbl$component), collapse = ", "))
  }

  chr_uniq <- df$fraction[df$type == "Binding level" & df$component == "Chromatin unique"]
  cg_uniq  <- df$fraction[df$type == "Binding level" & df$component == "CG unique"]
  ann_b <- sprintf("Binding level:\nChromatin %s vs CG %s",
                   stat_label("unique", chr_uniq * 100, "%.1f%%"),
                   stat_label("unique", cg_uniq * 100, "%.1f%%"))

  # Label only the explained (non-Unexplained) components to avoid clutter.
  df <- df %>%
    mutate(lab = ifelse(.data$component == "Unexplained", "",
                        sprintf("%.1f%%", .data$fraction * 100)))

  ggplot(df, aes(x = .data$type, y = .data$fraction, fill = .data$component)) +
    geom_col(width = 0.65, colour = "black", linewidth = 0.2) +
    geom_text(aes(label = .data$lab),
              position = position_stack(vjust = 0.5), size = 2.0) +
    scale_fill_manual(
      values = c("CG unique" = COL_5MC, "Shared" = "#999999",
                 "Chromatin unique" = COL_K119UB, "Unexplained" = "grey85"),
      name = "Variance"
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.02))) +
    annotate("text", x = 1.5, y = 1.02, label = ann_b,
             size = 2.2, fontface = "italic", hjust = 0.5, vjust = 0,
             lineheight = 0.95) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    labs(
      title = "Variance in MeCP2 binding partitioned by predictor class",
      x = NULL,
      y = "Fraction of MeCP2 variance"
    ) +
    theme_pub() +
    theme(legend.position = "top",
          plot.title = element_text(size = 8),
          plot.margin = margin(14, 5, 5, 5))
}

# =============================================================================
# PANEL C (Sec 67): gene-body K119ub log2FC partitioned by MeCP2 status.
# Reads 59_quadrant_master.tsv for per-gene K119ub (gb_log2fc) + MeCP2 status,
# and 67_statistics.tsv for the OR=3.15 annotation (MeCP2-up-no-methylation set).
# MeCP2-Up loci gain K119ub; MeCP2-Down loci are flat (matches panel D + Sec 60).
# =============================================================================
build_panel_c <- function() {
  master <- read_table_tsv("59_quadrant_master.tsv")

  required <- c("gb_log2fc", "mecp2_mean_fold", "mecp2_nearest_fdr")
  missing <- setdiff(required, names(master))
  if (length(missing) > 0) {
    stop("Panel C: 59_quadrant_master.tsv missing columns: ",
         paste(missing, collapse = ", "))
  }

  # Documented MeCP2-status rule (Sec 59): significant nearest-peak MeCP2 change
  # with positive concentration-weighted fold = Up, negative = Down, else NS.
  df <- master %>%
    mutate(
      gb_log2fc       = suppressWarnings(as.numeric(.data$gb_log2fc)),
      mecp2_mean_fold = suppressWarnings(as.numeric(.data$mecp2_mean_fold)),
      mecp2_nearest_fdr = suppressWarnings(as.numeric(.data$mecp2_nearest_fdr))
    ) %>%
    filter(is.finite(.data$gb_log2fc)) %>%
    mutate(
      mecp2_status = case_when(
        is.finite(.data$mecp2_nearest_fdr) & .data$mecp2_nearest_fdr < Q_THRESHOLD &
          is.finite(.data$mecp2_mean_fold) & .data$mecp2_mean_fold > 0 ~ "MeCP2 Up",
        is.finite(.data$mecp2_nearest_fdr) & .data$mecp2_nearest_fdr < Q_THRESHOLD &
          is.finite(.data$mecp2_mean_fold) & .data$mecp2_mean_fold < 0 ~ "MeCP2 Down",
        TRUE ~ "Not Significant"
      )
    )

  status_levels <- c("MeCP2 Up", "MeCP2 Down", "Not Significant")
  df$mecp2_status <- factor(df$mecp2_status, levels = status_levels)

  n_by_status <- table(df$mecp2_status)
  if (n_by_status[["MeCP2 Up"]] == 0L || n_by_status[["MeCP2 Down"]] == 0L) {
    stop("Panel C: MeCP2 Up/Down partition is empty -- check 59 master schema.")
  }
  cat("  Panel C MeCP2-status partition (K119ub-quantifiable genes): Up=",
      n_by_status[["MeCP2 Up"]], " Down=", n_by_status[["MeCP2 Down"]],
      " NS=", n_by_status[["Not Significant"]], "\n", sep = "")

  # OR for the MeCP2-up-no-methylation cohort (Sec 67 key result), annotated.
  stat67 <- read_table_tsv("67_statistics.tsv")
  fisher_row <- stat67[grepl("^Fisher", stat67$test), , drop = FALSE]
  if (nrow(fisher_row) != 1L) {
    stop("Panel C: could not locate the single Fisher row in 67_statistics.tsv")
  }
  fisher_or <- as.numeric(sub("OR = ", "", fisher_row$statistic[1]))
  fisher_p  <- as.numeric(fisher_row$p_value[1])
  if (!is.finite(fisher_or)) {
    stop("Panel C: failed to parse Fisher OR from 67_statistics.tsv")
  }

  # Per-group medians for the on-panel summary (genome-wide K119ub gain at Up).
  med_up   <- median(df$gb_log2fc[df$mecp2_status == "MeCP2 Up"])
  med_down <- median(df$gb_log2fc[df$mecp2_status == "MeCP2 Down"])

  n_labels <- sprintf("n=%d", as.integer(n_by_status[status_levels]))
  names(n_labels) <- status_levels

  ann_c <- sprintf(
    "MeCP2-up-without-methylation genes:\n%s gain K119ub (p = %.1e)",
    stat_label("OR", fisher_or, "%.2f"), fisher_p
  )

  y_top <- quantile(df$gb_log2fc, 0.995, na.rm = TRUE)
  y_bot <- quantile(df$gb_log2fc, 0.005, na.rm = TRUE)

  ggplot(df, aes(x = .data$mecp2_status, y = .data$gb_log2fc,
                 fill = .data$mecp2_status)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50",
               linewidth = 0.3) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, alpha = 0.85) +
    geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.3,
                 fill = "white", alpha = 0.9) +
    scale_fill_manual(values = MECP2_STATUS_COLORS, guide = "none") +
    scale_x_discrete(labels = function(x) paste0(x, "\n", n_labels[x])) +
    annotate("text", x = 2, y = y_top * 1.18, label = ann_c,
             size = 2.2, fontface = "italic", hjust = 0.5, lineheight = 0.95) +
    coord_cartesian(ylim = c(y_bot, y_top * 1.3)) +
    labs(
      title = "Gene-body H2AK119ub change by MeCP2 binding status",
      x = "MeCP2 (CUT&RUN) status",
      y = expression("Gene-body H2AK119ub "*log[2]*"FC")
    ) +
    theme_pub() +
    theme(plot.title = element_text(size = 8))
}

# =============================================================================
# PANEL D (Sec 60): mirror-image epigenetic profiles by MeCP2 status.
# Reshapes 60_mecp2_status_stats.tsv (wide -> long) to plot the mean change of
# 5 marks (K119ub BigWig, 5mC, 5hmC, K27ac, K27me3) for MeCP2 Up / Down / NS.
# Punchline: K119ub is GAINED at MeCP2-Up but FLAT (ns) at MeCP2-Down.
# 95% CI bars are derived from each row's mean and one-sample t-statistic
# (SE = mean / t), the only dispersion recoverable from this summary table.
# =============================================================================
build_panel_d <- function() {
  tbl <- read_table_tsv("60_mecp2_status_stats.tsv")

  required <- c("group", "mark", "mean", "ttest_t", "wilcox_q")
  missing <- setdiff(required, names(tbl))
  if (length(missing) > 0) {
    stop("Panel D: 60_mecp2_status_stats.tsv missing columns: ",
         paste(missing, collapse = ", "))
  }

  # Use the BigWig K119ub row (k119ub_bw) -- the genome-wide signal measure that
  # is FLAT at MeCP2-Down -- plus the CG and active/repressive histone marks.
  mark_map <- c(
    "k119ub_bw" = "H2AK119ub",
    "mc"        = "5mC",
    "hmc"       = "5hmC",
    "k27ac"     = "H3K27ac",
    "k27me3"    = "H3K27me3"
  )
  mark_levels <- c("H2AK119ub", "5mC", "5hmC", "H3K27ac", "H3K27me3")

  df <- tbl %>%
    filter(.data$mark %in% names(mark_map)) %>%
    mutate(
      mark_lab = factor(mark_map[.data$mark], levels = mark_levels),
      group    = factor(.data$group,
                        levels = c("MeCP2 Up", "MeCP2 Down", "Not Significant")),
      mean     = as.numeric(.data$mean),
      ttest_t  = as.numeric(.data$ttest_t),
      wilcox_q = as.numeric(.data$wilcox_q),
      se       = ifelse(.data$ttest_t == 0, NA_real_,
                        abs(.data$mean / .data$ttest_t)),
      ci_lo    = .data$mean - 1.96 * .data$se,
      ci_hi    = .data$mean + 1.96 * .data$se,
      # Significance star from the BH-adjusted one-sample Wilcoxon q-value.
      sig      = ifelse(.data$wilcox_q < 0.05, "*", "ns")
    )
  if (nrow(df) != 15L) {
    stop("Panel D: expected 15 rows (5 marks x 3 MeCP2 groups), got ", nrow(df))
  }

  # Confirm the punchline value: K119ub flat (ns) at MeCP2-Down.
  k119_down <- df %>% filter(.data$mark_lab == "H2AK119ub", .data$group == "MeCP2 Down")
  cat("  Panel D K119ub at MeCP2-Down: mean=", round(k119_down$mean, 3),
      " (wilcox_q=", signif(k119_down$wilcox_q, 3), ", ", k119_down$sig, ")\n", sep = "")

  ann_d <- sprintf("H2AK119ub is gained at MeCP2-Up\nbut flat at MeCP2-Down (%s)",
                   stat_label("mean", k119_down$mean, "%.2f"))

  dodge <- position_dodge(width = 0.7)

  ggplot(df, aes(x = .data$mark_lab, y = .data$mean,
                 colour = .data$group, group = .data$group)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50",
               linewidth = 0.3) +
    geom_linerange(aes(ymin = .data$ci_lo, ymax = .data$ci_hi),
                   position = dodge, linewidth = 0.5) +
    geom_point(position = dodge, size = 2.2) +
    geom_text(aes(label = .data$sig,
                  y = .data$ci_hi + 0.06 * sign(pmax(.data$mean, 0) + 0.01)),
              position = dodge, size = 2.4, show.legend = FALSE,
              vjust = 0) +
    scale_colour_manual(values = MECP2_STATUS_COLORS, name = "MeCP2 status") +
    annotate("text", x = 4, y = max(df$ci_hi, na.rm = TRUE) * 0.95,
             label = ann_d, size = 2.2, fontface = "italic",
             hjust = 0.5, lineheight = 0.95) +
    labs(
      title = "Mirror-image epigenetic profiles define MeCP2 status",
      x = "Epigenetic mark (mutant - control)",
      y = "Mean change at gene bodies"
    ) +
    theme_pub() +
    theme(legend.position = "top",
          plot.title = element_text(size = 8),
          axis.text.x = element_text(angle = 30, hjust = 1))
}

# =============================================================================
# PANEL E (Sec 71): variance of MeCP2 fold partitioned into K119ub vs ratio.
# Reads 71_variance_partition.tsv. Key stat: K119ub unique 7.3% vs methylation-
# ratio unique 0.0% -- MeCP2 responds to the Polycomb mark, not the mC/hmC skew.
# =============================================================================
build_panel_e <- function() {
  tbl <- read_table_tsv("71_variance_partition.tsv")

  required <- c("component", "fraction")
  missing <- setdiff(required, names(tbl))
  if (length(missing) > 0) {
    stop("Panel E: 71_variance_partition.tsv missing columns: ",
         paste(missing, collapse = ", "))
  }

  comp_levels <- c("K119ub unique", "Shared", "Ratio unique", "Unexplained")
  df <- tbl %>%
    mutate(
      fraction  = as.numeric(.data$fraction),
      component = factor(.data$component, levels = comp_levels)
    )
  if (any(is.na(df$component))) {
    stop("Panel E: unexpected variance component label(s): ",
         paste(unique(tbl$component), collapse = ", "))
  }

  k119_frac  <- df$fraction[df$component == "K119ub unique"]
  ratio_frac <- df$fraction[df$component == "Ratio unique"]
  ann_e <- sprintf("%s vs methylation-ratio %s\nof MeCP2-fold variance",
                   stat_label("K119ub", k119_frac * 100, "%.1f%%"),
                   stat_label("", ratio_frac * 100, "%.1f%%"))

  # Single stacked bar; label every explained component (the ratio slice is tiny
  # but must be visible/annotated to make the K119ub-dominance argument).
  df <- df %>%
    mutate(
      x   = "MeCP2 fold",
      lab = ifelse(.data$component == "Unexplained", "",
                   sprintf("%s %.1f%%", .data$component, .data$fraction * 100))
    )

  ggplot(df, aes(x = .data$x, y = .data$fraction, fill = .data$component)) +
    geom_col(width = 0.4, colour = "black", linewidth = 0.2) +
    scale_fill_manual(
      values = c("K119ub unique" = COL_K119UB, "Shared" = "#999999",
                 "Ratio unique" = COL_5MC, "Unexplained" = "grey85"),
      name = "Variance"
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.02))) +
    annotate("text", x = 1.45, y = 0.55, label = ann_e,
             size = 2.4, fontface = "italic", hjust = 0, lineheight = 0.95) +
    coord_flip(clip = "off") +
    labs(
      title = "MeCP2 binding tracks H2AK119ub, not the 5mC/5hmC methylation skew",
      x = NULL,
      y = "Fraction of MeCP2-fold variance"
    ) +
    theme_pub() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 8),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}

# =============================================================================
# Build panels.
# =============================================================================
cat("Building panel A (Sec 62 model R^2 comparison)...\n")
panel_a <- build_panel_a()

cat("Building panel B (Sec 62 variance partition)...\n")
panel_b <- build_panel_b()

cat("Building panel C (Sec 67 K119ub by MeCP2 status)...\n")
panel_c <- build_panel_c()

cat("Building panel D (Sec 60 mirror-image profiles)...\n")
panel_d <- build_panel_d()

cat("Building panel E (Sec 71 K119ub vs ratio variance)...\n")
panel_e <- build_panel_e()

# =============================================================================
# Compose with patchwork using the plan's design + panel tags.
# design = "AABB\nCCDD\nEEEE"; heights = c(1, 1, 0.8).
# =============================================================================
cat("Composing Section 86 (design AABB / CCDD / EEEE)...\n")

design <- "AABB\nCCDD\nEEEE"

sec86_composite <- panel_a + panel_b + panel_c + panel_d + panel_e +
  plot_layout(design = design, heights = c(1, 1, 0.8)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Section 86: MeCP2 reads chromatin state, not DNA methylation",
    theme = theme(plot.title = element_text(face = "bold", size = 11, hjust = 0))
  )

cat("Saving Section 86 panels + composite...\n")
save_plot(panel_a, "86a_model_r2",        w = 8, h = 6)
save_plot(panel_b, "86b_variance",        w = 8, h = 6)
save_plot(panel_c, "86c_k119ub_by_mecp2", w = 8, h = 6)
save_plot(panel_d, "86d_mirror_profiles", w = 8, h = 6)
save_plot(panel_e, "86e_ratio_variance",  w = 14, h = 5)
save_plot(sec86_composite, "86_composite", w = 14, h = 12)

cat("\n=== Section 86 complete ===\n")
cat("Output: ", SEC86_DIR, "\n", sep = "")
