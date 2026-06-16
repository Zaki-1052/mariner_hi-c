# biomodal/downstream/scripts/viz_sections/section_54_tad_methylation_organization.R
# Section 54: Non-CG vs CG Methylation TAD Organization
#
# Tests whether non-CG methylation (CHG/CHH) is more organized within TADs
# than CpG methylation, using inter/intra-TAD variance ratios, intra-TAD CV,
# boundary insulation scores, and meta-TAD profiles.
#
# Input: data/tad_methylation_signal_late.tsv (from compute_tad_methylation.R)
# Usage: cd downstream && Rscript scripts/viz_sections/section_54_tad_methylation_organization.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 54 CONFIGURATION
# =============================================================================

SEC54_DIR <- file.path(OUTPUT_DIR, "54_tad_methylation_organization")
dir.create(SEC54_DIR, recursive = TRUE, showWarnings = FALSE)

TAD_SIGNAL_FILE <- file.path(BASE_DIR, "data/tad_methylation_signal_late.tsv")

CONTEXTS <- c("cg_mc", "cg_hmc", "chg_mc", "chh_mc", "mecp2")
CONTEXT_LABELS <- c(
  cg_mc  = "CG 5mC",
  cg_hmc = "CG 5hmC",
  chg_mc = "CHG 5mC",
  chh_mc = "CHH 5mC",
  mecp2  = "MeCP2"
)

CONTEXT_COLORS <- c(
  "CG 5mC"  = "#E41A1C",
  "CG 5hmC" = "#377EB8",
  "CHG 5mC" = "#FF7F00",
  "CHH 5mC" = "#4DAF4A",
  "MeCP2"   = "#984EA3"
)

# =============================================================================
# LOAD DATA
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 54: TAD METHYLATION ORGANIZATION\n")
cat("================================================================================\n\n")

if (!file.exists(TAD_SIGNAL_FILE)) {
  stop("TAD signal file not found: ", TAD_SIGNAL_FILE,
       "\nRun compute_tad_methylation.R on Expanse first, then rsync the TSV.")
}

tad <- read.table(TAD_SIGNAL_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("Loaded %d TAD domains from %s\n", nrow(tad), TAD_SIGNAL_FILE))
cat(sprintf("  Differential: %d (%.1f%%)\n",
            sum(tad$is_differential), 100 * sum(tad$is_differential) / nrow(tad)))
cat(sprintf("  Width range: %.0f - %.0f kb (median: %.0f kb)\n\n",
            min(tad$width_kb), max(tad$width_kb), median(tad$width_kb)))

# =============================================================================
# FIGURE 54a: META-TAD METHYLATION PROFILE (Z-SCORED)
# =============================================================================

cat("--- 54a: Meta-TAD profile ---\n")

build_meta_profile <- function(tad_df, context, condition = "ctrl") {
  mean_col <- paste0(context, "_", condition, "_mean")
  if (!(mean_col %in% names(tad_df))) return(NULL)
  vals <- tad_df[[mean_col]]
  valid <- !is.na(vals)
  data.frame(
    context = CONTEXT_LABELS[context],
    condition = condition,
    value = vals[valid],
    stringsAsFactors = FALSE
  )
}

profile_data <- do.call(rbind, lapply(CONTEXTS, function(ctx) {
  rbind(
    build_meta_profile(tad, ctx, "ctrl"),
    build_meta_profile(tad, ctx, "mut")
  )
}))

if (!is.null(profile_data) && nrow(profile_data) > 0) {
  profile_summary <- profile_data %>%
    group_by(context, condition) %>%
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      median_val = median(value, na.rm = TRUE),
      sd_val = sd(value, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )

  p54a <- ggplot(profile_data, aes(x = context, y = value, fill = condition)) +
    geom_violin(alpha = 0.7, scale = "width", position = position_dodge(0.8)) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, position = position_dodge(0.8),
                 alpha = 0.9) +
    scale_fill_manual(values = c("ctrl" = "#2166AC", "mut" = "#B2182B"),
                      labels = c("ctrl" = "Control", "mut" = "Mutant")) +
    facet_wrap(~ context, scales = "free_y", nrow = 1) +
    labs(title = "Mean Methylation per TAD Domain",
         subtitle = "Per-context distribution across TAD domains",
         x = NULL, y = "Mean methylation signal", fill = "Condition") +
    theme_biomodal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  save_multiformat_ggplot(p54a, file.path(SEC54_DIR, "54a_metatad_signal"),
                          width = 12, height = 7)
  cat("  54a saved.\n")
}

# =============================================================================
# FIGURE 54b: BOUNDARY INSULATION SCORES BY CONTEXT
# =============================================================================

cat("--- 54b: Boundary insulation ---\n")

ins_data <- do.call(rbind, lapply(CONTEXTS, function(ctx) {
  ctrl_col <- paste0(ctx, "_ctrl_insulation")
  mut_col  <- paste0(ctx, "_mut_insulation")
  rows <- list()
  if (ctrl_col %in% names(tad)) {
    valid <- !is.na(tad[[ctrl_col]])
    rows[[1]] <- data.frame(context = CONTEXT_LABELS[ctx],
                            condition = "Control",
                            insulation = tad[[ctrl_col]][valid])
  }
  if (mut_col %in% names(tad)) {
    valid <- !is.na(tad[[mut_col]])
    rows[[2]] <- data.frame(context = CONTEXT_LABELS[ctx],
                            condition = "Mutant",
                            insulation = tad[[mut_col]][valid])
  }
  do.call(rbind, rows)
}))

if (!is.null(ins_data) && nrow(ins_data) > 0) {
  p54b <- ggplot(ins_data, aes(x = context, y = insulation, fill = context)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.9) +
    facet_wrap(~ condition) +
    scale_fill_manual(values = CONTEXT_COLORS) +
    labs(title = "TAD Boundary Insulation Score by Methylation Context",
         subtitle = "|inside - outside| signal at domain boundaries (50kb flanks)",
         x = NULL, y = "Boundary insulation score") +
    theme_biomodal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))

  # Pairwise Wilcoxon: CG vs CHG, CG vs CHH (ctrl only)
  ctrl_ins <- ins_data[ins_data$condition == "Control", ]
  if (nrow(ctrl_ins) > 0) {
    for (alt_ctx in c("CHG 5mC", "CHH 5mC")) {
      cg_vals <- ctrl_ins$insulation[ctrl_ins$context == "CG 5mC"]
      alt_vals <- ctrl_ins$insulation[ctrl_ins$context == alt_ctx]
      n_compare <- min(length(cg_vals), length(alt_vals))
      if (n_compare > 10) {
        wt <- wilcox.test(cg_vals[1:n_compare], alt_vals[1:n_compare])
        cat(sprintf("  CG 5mC vs %s insulation: W=%.0f, p=%.2e\n",
                    alt_ctx, wt$statistic, wt$p.value))
      }
    }
  }

  save_multiformat_ggplot(p54b, file.path(SEC54_DIR, "54b_boundary_insulation"),
                          width = 10, height = 7)
  cat("  54b saved.\n")
}

# =============================================================================
# FIGURE 54c: INTRA-TAD COEFFICIENT OF VARIATION
# =============================================================================

cat("--- 54c: Intra-TAD CV ---\n")

cv_data <- do.call(rbind, lapply(CONTEXTS, function(ctx) {
  ctrl_col <- paste0(ctx, "_ctrl_cv")
  mut_col  <- paste0(ctx, "_mut_cv")
  rows <- list()
  if (ctrl_col %in% names(tad)) {
    valid <- !is.na(tad[[ctrl_col]])
    rows[[1]] <- data.frame(context = CONTEXT_LABELS[ctx],
                            condition = "Control",
                            cv = tad[[ctrl_col]][valid])
  }
  if (mut_col %in% names(tad)) {
    valid <- !is.na(tad[[mut_col]])
    rows[[2]] <- data.frame(context = CONTEXT_LABELS[ctx],
                            condition = "Mutant",
                            cv = tad[[mut_col]][valid])
  }
  do.call(rbind, rows)
}))

if (!is.null(cv_data) && nrow(cv_data) > 0) {
  p54c <- ggplot(cv_data, aes(x = context, y = cv, fill = context)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.9) +
    facet_wrap(~ condition) +
    scale_fill_manual(values = CONTEXT_COLORS) +
    labs(title = "Intra-TAD Coefficient of Variation",
         subtitle = "Lower CV = more homogeneous signal within TADs",
         x = NULL, y = "Coefficient of variation (SD/mean)") +
    theme_biomodal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))

  ctrl_cv <- cv_data[cv_data$condition == "Control", ]
  for (alt_ctx in c("CHG 5mC", "CHH 5mC")) {
    cg_vals <- ctrl_cv$cv[ctrl_cv$context == "CG 5mC"]
    alt_vals <- ctrl_cv$cv[ctrl_cv$context == alt_ctx]
    n_compare <- min(length(cg_vals), length(alt_vals))
    if (n_compare > 10) {
      wt <- wilcox.test(cg_vals[1:n_compare], alt_vals[1:n_compare])
      cat(sprintf("  CG 5mC vs %s CV: W=%.0f, p=%.2e\n",
                  alt_ctx, wt$statistic, wt$p.value))
    }
  }

  save_multiformat_ggplot(p54c, file.path(SEC54_DIR, "54c_intratad_cv"),
                          width = 10, height = 7)
  cat("  54c saved.\n")
}

# =============================================================================
# FIGURE 54d: INTER/INTRA TAD VARIANCE RATIO (KEY FIGURE)
# =============================================================================

cat("--- 54d: Inter/Intra variance ratio ---\n")

compute_variance_ratio <- function(domain_means, within_vars) {
  valid <- !is.na(domain_means) & !is.na(within_vars) & within_vars > 0
  if (sum(valid) < 10) return(NA)
  var(domain_means[valid]) / mean(within_vars[valid])
}

bootstrap_variance_ratio <- function(domain_means, within_vars, B = 1000, seed = 42) {
  valid <- !is.na(domain_means) & !is.na(within_vars) & within_vars > 0
  means_v <- domain_means[valid]
  vars_v  <- within_vars[valid]
  n <- length(means_v)
  if (n < 10) return(c(estimate = NA, ci_lo = NA, ci_hi = NA))

  set.seed(seed)
  boot_ratios <- replicate(B, {
    idx <- sample(n, replace = TRUE)
    var(means_v[idx]) / mean(vars_v[idx])
  })

  c(estimate = var(means_v) / mean(vars_v),
    ci_lo = unname(quantile(boot_ratios, 0.025)),
    ci_hi = unname(quantile(boot_ratios, 0.975)))
}

vr_results <- do.call(rbind, lapply(CONTEXTS, function(ctx) {
  rows <- list()
  for (cond in c("ctrl", "mut")) {
    mean_col <- paste0(ctx, "_", cond, "_mean")
    var_col  <- paste0(ctx, "_", cond, "_within_var")
    if (mean_col %in% names(tad) && var_col %in% names(tad)) {
      boot <- bootstrap_variance_ratio(tad[[mean_col]], tad[[var_col]])
      rows[[cond]] <- data.frame(
        context = CONTEXT_LABELS[ctx],
        condition = ifelse(cond == "ctrl", "Control", "Mutant"),
        ratio = boot["estimate"],
        ci_lo = boot["ci_lo"],
        ci_hi = boot["ci_hi"],
        stringsAsFactors = FALSE, row.names = NULL
      )
    }
  }
  do.call(rbind, rows)
}))

if (!is.null(vr_results) && nrow(vr_results) > 0 && any(!is.na(vr_results$ratio))) {
  vr_results$context <- factor(vr_results$context, levels = unname(CONTEXT_LABELS))

  p54d <- ggplot(vr_results, aes(x = context, y = ratio, fill = context)) +
    geom_col(alpha = 0.8, position = position_dodge(0.8)) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                  width = 0.2, position = position_dodge(0.8)) +
    facet_wrap(~ condition) +
    scale_fill_manual(values = CONTEXT_COLORS) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    labs(title = "Inter-TAD / Intra-TAD Variance Ratio",
         subtitle = "Higher ratio = methylation more organized by TAD structure\n(dashed line: ratio=1, no TAD organization)",
         x = NULL, y = "Variance ratio (inter / intra)") +
    theme_biomodal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))

  save_multiformat_ggplot(p54d, file.path(SEC54_DIR, "54d_variance_ratio"),
                          width = 10, height = 7)

  cat("  Variance ratios:\n")
  for (i in seq_len(nrow(vr_results))) {
    cat(sprintf("    %s (%s): %.4f [%.4f, %.4f]\n",
                vr_results$context[i], vr_results$condition[i],
                vr_results$ratio[i], vr_results$ci_lo[i], vr_results$ci_hi[i]))
  }
  cat("  54d saved.\n")
}

# =============================================================================
# FIGURE 54e: BOUNDARY SCORE VS METHYLATION INSULATION CORRELATION
# =============================================================================

cat("--- 54e: Boundary score correlation ---\n")

cor_data <- do.call(rbind, lapply(CONTEXTS, function(ctx) {
  ins_col <- paste0(ctx, "_ctrl_insulation")
  if (!(ins_col %in% names(tad))) return(NULL)
  valid <- !is.na(tad$mean_tad_score) & !is.na(tad[[ins_col]])
  data.frame(
    context = CONTEXT_LABELS[ctx],
    tad_score = tad$mean_tad_score[valid],
    insulation = tad[[ins_col]][valid],
    stringsAsFactors = FALSE
  )
}))

if (!is.null(cor_data) && nrow(cor_data) > 0) {
  cor_stats <- cor_data %>%
    group_by(context) %>%
    summarise(
      rho = cor(tad_score, insulation, method = "spearman", use = "complete.obs"),
      p_val = cor.test(tad_score, insulation, method = "spearman")$p.value,
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(label = sprintf("rho = %.3f\np = %.1e\nn = %d", rho, p_val, n))

  p54e <- ggplot(cor_data, aes(x = tad_score, y = insulation)) +
    geom_point(alpha = 0.15, size = 0.5) +
    geom_smooth(method = "lm", color = "#D73027", se = TRUE, linewidth = 0.8) +
    geom_text(data = cor_stats,
              aes(x = Inf, y = Inf, label = label),
              hjust = 1.1, vjust = 1.3, size = 3, inherit.aes = FALSE) +
    facet_wrap(~ context, scales = "free") +
    labs(title = "TAD Boundary Strength vs Methylation Insulation",
         subtitle = "Spearman correlation between TAD score and boundary methylation transition",
         x = "Mean TAD boundary score (TADCompare)", y = "Methylation insulation score") +
    theme_biomodal()

  save_multiformat_ggplot(p54e, file.path(SEC54_DIR, "54e_boundary_correlation"),
                          width = 12, height = 10)
  cat("  54e saved.\n")
}

# =============================================================================
# FIGURE 54f: DIFFERENTIAL TAD × METHYLATION OVERLAY
# =============================================================================

cat("--- 54f: Differential TAD overlay ---\n")

diff_data <- do.call(rbind, lapply(CONTEXTS, function(ctx) {
  fc_col <- paste0(ctx, "_log2fc")
  if (!(fc_col %in% names(tad))) return(NULL)
  valid <- !is.na(tad[[fc_col]])
  data.frame(
    context = CONTEXT_LABELS[ctx],
    is_differential = ifelse(tad$is_differential[valid], "Differential", "Non-Differential"),
    log2fc = tad[[fc_col]][valid],
    stringsAsFactors = FALSE
  )
}))

if (!is.null(diff_data) && nrow(diff_data) > 0) {
  p54f <- ggplot(diff_data, aes(x = is_differential, y = log2fc, fill = is_differential)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~ context, scales = "free_y") +
    scale_fill_manual(values = c("Differential" = "#D73027", "Non-Differential" = "#4575B4")) +
    labs(title = "Methylation Change at Differential vs Non-Differential TADs",
         subtitle = "log2(mutant/control) per TAD domain",
         x = NULL, y = "log2FC methylation") +
    theme_biomodal() +
    theme(legend.position = "none")

  for (ctx in CONTEXTS) {
    fc_col <- paste0(ctx, "_log2fc")
    if (fc_col %in% names(tad)) {
      diff_vals <- tad[[fc_col]][tad$is_differential & !is.na(tad[[fc_col]])]
      nondiff_vals <- tad[[fc_col]][!tad$is_differential & !is.na(tad[[fc_col]])]
      if (length(diff_vals) > 5 && length(nondiff_vals) > 5) {
        wt <- wilcox.test(diff_vals, nondiff_vals)
        cat(sprintf("  %s diff vs non-diff log2FC: W=%.0f, p=%.2e\n",
                    CONTEXT_LABELS[ctx], wt$statistic, wt$p.value))
      }
    }
  }

  save_multiformat_ggplot(p54f, file.path(SEC54_DIR, "54f_differential_overlay"),
                          width = 12, height = 8)
  cat("  54f saved.\n")
}

# =============================================================================
# FIGURE 54g: BOUNDARY TYPE STRATIFICATION
# =============================================================================

cat("--- 54g: Boundary type stratification ---\n")

BOUNDARY_TYPES <- c("Shifted", "Strength Change", "Split", "Merge",
                     "Complex", "Non-Differential")

btype_data <- do.call(rbind, lapply(CONTEXTS, function(ctx) {
  ins_col <- paste0(ctx, "_ctrl_insulation")
  if (!(ins_col %in% names(tad))) return(NULL)
  valid <- !is.na(tad[[ins_col]])
  left_type <- tad$boundary_left_type[valid]
  left_type[!(left_type %in% BOUNDARY_TYPES)] <- "Other"
  data.frame(
    context = CONTEXT_LABELS[ctx],
    boundary_type = left_type,
    insulation = tad[[ins_col]][valid],
    stringsAsFactors = FALSE
  )
}))

if (!is.null(btype_data) && nrow(btype_data) > 0) {
  btype_data$boundary_type <- factor(btype_data$boundary_type,
                                      levels = c(BOUNDARY_TYPES, "Other"))

  p54g <- ggplot(btype_data, aes(x = boundary_type, y = insulation,
                                  fill = boundary_type)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
    facet_wrap(~ context, scales = "free_y") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Methylation Insulation by TAD Boundary Type",
         subtitle = "How methylation transition varies with type of TAD boundary change",
         x = "Boundary type (TADCompare)", y = "Insulation score") +
    theme_biomodal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  save_multiformat_ggplot(p54g, file.path(SEC54_DIR, "54g_boundary_type"),
                          width = 14, height = 8)
  cat("  54g saved.\n")
}

# =============================================================================
# FIGURE 54h: CROSS-CONTEXT LOG2FC CORRELATION
# =============================================================================

cat("--- 54h: Cross-context FC correlation ---\n")

cross_plots <- list()
for (alt_ctx in c("chg_mc", "chh_mc", "mecp2")) {
  cg_col  <- "cg_mc_log2fc"
  alt_col <- paste0(alt_ctx, "_log2fc")
  if (!(cg_col %in% names(tad) && alt_col %in% names(tad))) next

  valid <- !is.na(tad[[cg_col]]) & !is.na(tad[[alt_col]])
  if (sum(valid) < 20) next

  plot_df <- data.frame(
    cg_fc  = tad[[cg_col]][valid],
    alt_fc = tad[[alt_col]][valid],
    is_differential = tad$is_differential[valid]
  )

  rho <- cor(plot_df$cg_fc, plot_df$alt_fc, method = "spearman")
  pval <- cor.test(plot_df$cg_fc, plot_df$alt_fc, method = "spearman")$p.value

  cross_plots[[alt_ctx]] <- ggplot(plot_df, aes(x = cg_fc, y = alt_fc,
                                                 color = is_differential)) +
    geom_point(alpha = 0.2, size = 0.5) +
    geom_smooth(method = "lm", color = "#D73027", se = TRUE, linewidth = 0.8,
                inherit.aes = FALSE, aes(x = cg_fc, y = alt_fc)) +
    scale_color_manual(values = c("TRUE" = "#D73027", "FALSE" = "grey60"),
                       labels = c("TRUE" = "Differential", "FALSE" = "Non-Diff")) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("rho = %.3f\np = %.1e", rho, pval),
             hjust = 1.1, vjust = 1.3, size = 3.5) +
    labs(x = "CG 5mC log2FC", y = paste0(CONTEXT_LABELS[alt_ctx], " log2FC"),
         color = "TAD type") +
    theme_biomodal()

  cat(sprintf("  CG vs %s log2FC: rho=%.3f, p=%.2e\n",
              CONTEXT_LABELS[alt_ctx], rho, pval))
}

if (length(cross_plots) > 0) {
  p54h <- wrap_plots(cross_plots, nrow = 1) +
    plot_annotation(
      title = "Cross-Context Methylation Change Correlation per TAD",
      subtitle = "Do TADs with CG methylation changes also show non-CG changes?"
    )
  save_multiformat_ggplot(p54h, file.path(SEC54_DIR, "54h_cross_context_fc"),
                          width = 12, height = 6)
  cat("  54h saved.\n")
}

# =============================================================================
# FIGURE 54i: SUMMARY STATISTICS TABLE
# =============================================================================

cat("--- 54i: Summary table ---\n")

summary_rows <- list()
for (ctx in CONTEXTS) {
  for (cond in c("ctrl", "mut")) {
    mean_col <- paste0(ctx, "_", cond, "_mean")
    cv_col   <- paste0(ctx, "_", cond, "_cv")
    ins_col  <- paste0(ctx, "_", cond, "_insulation")
    var_col  <- paste0(ctx, "_", cond, "_within_var")

    row <- data.frame(
      context = CONTEXT_LABELS[ctx],
      condition = ifelse(cond == "ctrl", "Control", "Mutant"),
      stringsAsFactors = FALSE
    )

    if (mean_col %in% names(tad))
      row$median_signal <- median(tad[[mean_col]], na.rm = TRUE)
    if (cv_col %in% names(tad))
      row$median_cv <- median(tad[[cv_col]], na.rm = TRUE)
    if (ins_col %in% names(tad))
      row$median_insulation <- median(tad[[ins_col]], na.rm = TRUE)
    if (mean_col %in% names(tad) && var_col %in% names(tad)) {
      vr <- compute_variance_ratio(tad[[mean_col]], tad[[var_col]])
      row$variance_ratio <- vr
    }

    summary_rows[[length(summary_rows) + 1]] <- row
  }
}

summary_df <- do.call(rbind, summary_rows)

summary_path <- file.path(TABLES_DIR, "54_tad_methylation_organization_summary.tsv")
write.table(summary_df, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Summary table written: %s\n", summary_path))
print(summary_df)

# =============================================================================
# FIGURE 54j: COMPOSITE PANEL
# =============================================================================

cat("--- 54j: Composite panel ---\n")

composite_parts <- list()
if (exists("p54a")) composite_parts$a <- p54a + labs(title = NULL, subtitle = "A. Signal per TAD")
if (exists("p54c")) composite_parts$c <- p54c + labs(title = NULL, subtitle = "B. Intra-TAD CV")
if (exists("p54d")) composite_parts$d <- p54d + labs(title = NULL, subtitle = "C. Variance Ratio")

if (length(composite_parts) >= 2) {
  p54j <- wrap_plots(composite_parts, nrow = 1) +
    plot_annotation(
      title = "Non-CG vs CG Methylation: TAD Organization",
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    )
  save_multiformat_ggplot(p54j, file.path(SEC54_DIR, "54j_composite"),
                          width = 18, height = 7)
  cat("  54j saved.\n")
}

# =============================================================================
# COMPLETE
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 54 COMPLETE\n")
cat(sprintf("Output: %s\n", SEC54_DIR))
cat(sprintf("Figures: %d panels generated\n",
            length(list.files(SEC54_DIR, recursive = TRUE, pattern = "\\.png$"))))
cat("================================================================================\n")
