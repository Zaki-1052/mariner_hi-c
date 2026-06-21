# biomodal/downstream/scripts/viz_sections/section_72g_neuronal_chromatin_remodeling.R
# Section 72g: Neuronal K119ub-High Genes Show Disproportionate Chromatin Remodeling
#
# Tests whether neuronal genes with constitutively high K119ub experience
# stronger chromatin state changes upon BAP1-KO: ATAC down, K27ac down, K27me3 up.
# No MeCP2 filtering. Independent confirmation of the theory:
#   BAP1 loss -> K119ub accumulation -> chromatin remodeling strongest at
#   neuronal genes -> MeCP2 redistribution follows chromatin, not methylation
#
# Panels:
#   72g-a: Multi-mark chromatin change (ATAC, K27ac, K27me3) neuronal vs non-neuronal
#   72g-b: 4-group comparison: K119ub level x neuronal status per mark
#   72g-c: Chromatin state distribution of K119ub-high neuronal vs non-neuronal
#   72g-d: Dose-response: K119ub decile vs mark change, neuronal vs non-neuronal
#   72g-e: Interaction models: mark_fold ~ k119ub_signal * is_neuronal
#   72g:   Composite
#
# Input:
#   data/k119ub_gene_signal.tsv, tables/diffbind_gene_level_all_marks.tsv,
#   tables/72_neuronal_gene_set_go_derived.tsv (from section 72)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_72g_neuronal_chromatin_remodeling.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(patchwork)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

cat("================================================================================\n")
cat("SECTION 72g: NEURONAL K119ub-HIGH GENES — CHROMATIN REMODELING\n")
cat("  Does BAP1-KO drive heterochromatin shift at neuronal K119ub loci?\n")
cat("================================================================================\n\n")

SEC72G_DIR <- file.path(OUTPUT_DIR, "72g_neuronal_chromatin_remodeling")
dir.create(SEC72G_DIR, recursive = TRUE, showWarnings = FALSE)

K119UB_SIGNAL_PATH <- file.path(BASE_DIR, "data/k119ub_gene_signal.tsv")
DIFFBIND_PATH      <- file.path(TABLES_DIR, "diffbind_gene_level_all_marks.tsv")
NEURONAL_PATH      <- file.path(TABLES_DIR, "72_neuronal_gene_set_go_derived.tsv")

stopifnot(
  "k119ub_gene_signal.tsv not found" = file.exists(K119UB_SIGNAL_PATH),
  "diffbind_gene_level_all_marks.tsv not found" = file.exists(DIFFBIND_PATH),
  "72_neuronal_gene_set_go_derived.tsv not found (run section 72 first)" = file.exists(NEURONAL_PATH)
)

fmt_p <- function(p) {
  if (is.na(p) || !is.finite(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

save_plot <- function(p, name, w = 10, h = 7) {
  save_multiformat_ggplot(p,
    base_path = file.path(SEC72G_DIR, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

MARK_COLORS <- c("ATAC" = "#E6AB02", "K27ac" = "#FF7F00", "K27me3" = "#756BB1")

# =============================================================================
# DATA LOADING AND MERGE
# =============================================================================

cat("--- Loading data ---\n\n")

k119ub <- read.table(K119UB_SIGNAL_PATH, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, quote = "")
k119ub_q <- k119ub[k119ub$gb_signal_class == "quantifiable", ]
cat(sprintf("  K119ub quantifiable genes: %d\n", nrow(k119ub_q)))

db <- read.table(DIFFBIND_PATH, header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, quote = "")
cat(sprintf("  DiffBind all-marks genes:  %d\n", nrow(db)))

neuronal_genes <- read.table(NEURONAL_PATH, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, quote = "")$gene
cat(sprintf("  Neuronal gene set:         %d genes\n", length(neuronal_genes)))

df <- merge(
  db[, c("gene", "chr", "start", "end", "chromatin_state",
         "atac_fold", "k27ac_fold", "k27me3_fold", "k119ub_fold",
         "mc_diff", "hmc_diff")],
  k119ub_q[, c("symbol", "gb_ctrl_signal", "gb_mut_signal", "gb_log2fc")],
  by.x = "gene", by.y = "symbol"
)

df$is_neuronal <- df$gene %in% neuronal_genes
df$group <- ifelse(df$is_neuronal, "Neuronal", "Non-neuronal")

ctrl_q3 <- quantile(df$gb_ctrl_signal, 0.75)
df$k119ub_level <- ifelse(df$gb_ctrl_signal >= ctrl_q3, "K119ub High", "K119ub Low")
df$k119ub_level <- factor(df$k119ub_level, levels = c("K119ub Low", "K119ub High"))

df$ctrl_decile <- as.integer(cut(df$gb_ctrl_signal,
                                  breaks = quantile(df$gb_ctrl_signal, probs = seq(0, 1, 0.1)),
                                  include.lowest = TRUE, labels = FALSE))

cat(sprintf("\n  Merged universe: %d genes\n", nrow(df)))
cat(sprintf("  Neuronal: %d (%.1f%%)\n", sum(df$is_neuronal),
            100 * sum(df$is_neuronal) / nrow(df)))
cat(sprintf("  K119ub High (ctrl Q3 >= %.3f): %d genes\n", ctrl_q3, sum(df$k119ub_level == "K119ub High")))
cat(sprintf("  Non-NA counts: ATAC=%d, K27ac=%d, K27me3=%d\n",
            sum(!is.na(df$atac_fold)), sum(!is.na(df$k27ac_fold)), sum(!is.na(df$k27me3_fold))))

# 4-group counts
for (lev in c("K119ub High", "K119ub Low")) {
  for (grp in c("Neuronal", "Non-neuronal")) {
    n <- sum(df$k119ub_level == lev & df$group == grp)
    cat(sprintf("  %s + %s: %d genes\n", lev, grp, n))
  }
}

# =============================================================================
# 72g-a: MULTI-MARK CHROMATIN CHANGE — NEURONAL vs NON-NEURONAL
# =============================================================================

cat("\n--- 72g-a: Multi-mark chromatin change ---\n")

mark_long <- rbind(
  data.frame(mark = "ATAC", fold = df$atac_fold, group = df$group, stringsAsFactors = FALSE),
  data.frame(mark = "K27ac", fold = df$k27ac_fold, group = df$group, stringsAsFactors = FALSE),
  data.frame(mark = "K27me3", fold = df$k27me3_fold, group = df$group, stringsAsFactors = FALSE)
)
mark_long <- mark_long[!is.na(mark_long$fold), ]
mark_long$mark <- factor(mark_long$mark, levels = c("ATAC", "K27ac", "K27me3"))
mark_long$group <- factor(mark_long$group, levels = c("Neuronal", "Non-neuronal"))

mark_stats <- data.frame()
for (m in c("ATAC", "K27ac", "K27me3")) {
  col <- switch(m, ATAC = "atac_fold", K27ac = "k27ac_fold", K27me3 = "k27me3_fold")
  neur_vals <- df[[col]][df$is_neuronal & !is.na(df[[col]])]
  other_vals <- df[[col]][!df$is_neuronal & !is.na(df[[col]])]

  wt <- wilcox.test(neur_vals, other_vals)
  mark_stats <- rbind(mark_stats, data.frame(
    mark = m,
    n_neuronal = length(neur_vals),
    n_other = length(other_vals),
    median_neuronal = median(neur_vals),
    median_other = median(other_vals),
    delta_median = median(neur_vals) - median(other_vals),
    wilcox_p = wt$p.value,
    stringsAsFactors = FALSE
  ))
  cat(sprintf("  %s: neuronal median=%+.4f (n=%d), other=%+.4f (n=%d), %s\n",
              m, median(neur_vals), length(neur_vals),
              median(other_vals), length(other_vals), fmt_p(wt$p.value)))
}

p_72ga <- ggplot(mark_long, aes(x = group, y = fold, fill = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_violin(alpha = 0.6, show.legend = FALSE, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.3, show.legend = FALSE) +
  facet_wrap(~ mark, scales = "free_y") +
  scale_fill_manual(values = c("Neuronal" = "#756BB1", "Non-neuronal" = "grey70")) +
  labs(
    title = "Chromatin mark changes in BAP1-KO: neuronal vs non-neuronal genes",
    subtitle = "DiffBind log2FC per mark | No MeCP2 filter",
    x = NULL, y = "DiffBind log2 fold change (mut/ctrl)"
  ) +
  theme_biomodal() +
  theme(strip.text = element_text(size = 12, face = "bold"))

save_plot(p_72ga, "72ga_multimark_neuronal_vs_other", w = 12, h = 7)

# =============================================================================
# 72g-b: 4-GROUP COMPARISON — K119ub LEVEL x NEURONAL STATUS
# =============================================================================

cat("\n--- 72g-b: 4-group comparison ---\n")

df$quad_group <- paste0(df$k119ub_level, "\n", df$group)
df$quad_group <- factor(df$quad_group, levels = c(
  "K119ub Low\nNon-neuronal", "K119ub Low\nNeuronal",
  "K119ub High\nNon-neuronal", "K119ub High\nNeuronal"
))

quad_long <- rbind(
  data.frame(mark = "ATAC", fold = df$atac_fold, quad = df$quad_group,
             k119ub_level = df$k119ub_level, group = df$group, stringsAsFactors = FALSE),
  data.frame(mark = "K27ac", fold = df$k27ac_fold, quad = df$quad_group,
             k119ub_level = df$k119ub_level, group = df$group, stringsAsFactors = FALSE),
  data.frame(mark = "K27me3", fold = df$k27me3_fold, quad = df$quad_group,
             k119ub_level = df$k119ub_level, group = df$group, stringsAsFactors = FALSE)
)
quad_long <- quad_long[!is.na(quad_long$fold), ]
quad_long$mark <- factor(quad_long$mark, levels = c("ATAC", "K27ac", "K27me3"))

quad_stats <- data.frame()
for (m in c("ATAC", "K27ac", "K27me3")) {
  col <- switch(m, ATAC = "atac_fold", K27ac = "k27ac_fold", K27me3 = "k27me3_fold")

  for (lev in c("K119ub High", "K119ub Low")) {
    neur <- df[[col]][df$k119ub_level == lev & df$is_neuronal & !is.na(df[[col]])]
    other <- df[[col]][df$k119ub_level == lev & !df$is_neuronal & !is.na(df[[col]])]

    if (length(neur) < 5 || length(other) < 5) next
    wt <- wilcox.test(neur, other)

    quad_stats <- rbind(quad_stats, data.frame(
      mark = m, k119ub_level = lev,
      n_neuronal = length(neur), n_other = length(other),
      median_neuronal = median(neur), median_other = median(other),
      wilcox_p = wt$p.value,
      stringsAsFactors = FALSE
    ))
    cat(sprintf("  %s at %s: neuronal=%+.4f (n=%d) vs other=%+.4f (n=%d) %s\n",
                m, lev, median(neur), length(neur),
                median(other), length(other), fmt_p(wt$p.value)))
  }
}

quad_summary <- aggregate(fold ~ mark + quad + k119ub_level + group,
                          data = quad_long, FUN = median)

p_72gb <- ggplot(quad_long, aes(x = quad, y = fold, fill = interaction(k119ub_level, group))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_violin(alpha = 0.5, show.legend = FALSE, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.2, show.legend = FALSE) +
  facet_wrap(~ mark, scales = "free_y") +
  scale_fill_manual(values = c(
    "K119ub Low.Non-neuronal" = "grey80",
    "K119ub Low.Neuronal" = "#BCBDDC",
    "K119ub High.Non-neuronal" = "#FEC44F",
    "K119ub High.Neuronal" = "#D73027"
  )) +
  labs(
    title = "Chromatin remodeling: K119ub level x neuronal status",
    subtitle = sprintf("K119ub High = top quartile ctrl signal (>= %.3f) | No MeCP2 filter", ctrl_q3),
    x = NULL, y = "DiffBind log2FC"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(size = 8, angle = 20, hjust = 1),
        strip.text = element_text(size = 12, face = "bold"))

save_plot(p_72gb, "72gb_4group_k119ub_x_neuronal", w = 14, h = 7)

# =============================================================================
# 72g-c: CHROMATIN STATE DISTRIBUTION — K119ub-HIGH BY NEURONAL STATUS
# =============================================================================

cat("\n--- 72g-c: Chromatin state distribution ---\n")

k119_high <- df[df$k119ub_level == "K119ub High", ]

state_counts <- as.data.frame(table(
  group = k119_high$group,
  state = k119_high$chromatin_state
))

state_totals <- aggregate(Freq ~ group, data = state_counts, FUN = sum)
state_counts <- merge(state_counts, state_totals, by = "group", suffixes = c("", "_total"))
state_counts$fraction <- state_counts$Freq / state_counts$Freq_total

state_order <- c("Active_Promoter", "Active_Enhancer", "Bivalent_Promoter",
                 "Poised_Enhancer", "Repressed_Promoter", "Polycomb", "Other")
state_counts$state <- factor(state_counts$state,
                              levels = intersect(state_order, unique(state_counts$state)))

cat("  Chromatin state distribution of K119ub-high genes:\n")
for (grp in c("Neuronal", "Non-neuronal")) {
  sub <- state_counts[state_counts$group == grp, ]
  sub <- sub[order(-sub$fraction), ]
  cat(sprintf("    %s (n=%d):\n", grp, sum(sub$Freq)))
  for (i in seq_len(nrow(sub))) {
    if (sub$Freq[i] > 0) {
      cat(sprintf("      %-22s %4d (%.1f%%)\n",
                  sub$state[i], sub$Freq[i], 100 * sub$fraction[i]))
    }
  }
}

# Fisher test for each non-trivial state
cat("\n  Fisher tests (neuronal vs non-neuronal among K119ub-high):\n")
chromatin_fisher <- data.frame()
for (s in levels(state_counts$state)) {
  n_neur_in <- sum(k119_high$is_neuronal & k119_high$chromatin_state == s)
  n_neur_out <- sum(k119_high$is_neuronal & k119_high$chromatin_state != s)
  n_other_in <- sum(!k119_high$is_neuronal & k119_high$chromatin_state == s)
  n_other_out <- sum(!k119_high$is_neuronal & k119_high$chromatin_state != s)

  if (n_neur_in + n_other_in < 5) next

  ft <- fisher.test(matrix(c(n_neur_in, n_other_in, n_neur_out, n_other_out), nrow = 2))
  cat(sprintf("    %-22s neuronal=%d/%d (%.1f%%) vs other=%d/%d (%.1f%%)  OR=%.2f  %s\n",
              s, n_neur_in, sum(k119_high$is_neuronal),
              100 * n_neur_in / sum(k119_high$is_neuronal),
              n_other_in, sum(!k119_high$is_neuronal),
              100 * n_other_in / sum(!k119_high$is_neuronal),
              ft$estimate, fmt_p(ft$p.value)))

  chromatin_fisher <- rbind(chromatin_fisher, data.frame(
    state = s, OR = as.numeric(ft$estimate), p_value = ft$p.value,
    n_neuronal = n_neur_in, n_other = n_other_in,
    stringsAsFactors = FALSE
  ))
}

state_colors_extended <- c(CHROMATIN_STATE_COLORS, "Other" = "grey70")

p_72gc <- ggplot(state_counts[state_counts$Freq > 0, ],
                 aes(x = group, y = fraction, fill = state)) +
  geom_col(alpha = 0.85, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = state_colors_extended, name = "Chromatin State") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Chromatin state composition: K119ub-high genes",
    subtitle = sprintf("K119ub top quartile | Neuronal: %d genes, Non-neuronal: %d genes",
                       sum(k119_high$is_neuronal), sum(!k119_high$is_neuronal)),
    x = NULL, y = "Fraction of genes"
  ) +
  theme_biomodal() +
  theme(legend.position = "right")

save_plot(p_72gc, "72gc_chromatin_state_k119ub_high", w = 10, h = 7)

# =============================================================================
# 72g-d: DOSE-RESPONSE — K119ub DECILE vs MARK CHANGE BY NEURONAL STATUS
# =============================================================================

cat("\n--- 72g-d: Dose-response by K119ub decile ---\n")

decile_mark_summary <- data.frame()
for (m in c("ATAC", "K27ac", "K27me3")) {
  col <- switch(m, ATAC = "atac_fold", K27ac = "k27ac_fold", K27me3 = "k27me3_fold")

  for (d in 1:10) {
    for (grp in c(TRUE, FALSE)) {
      mask <- df$ctrl_decile == d & df$is_neuronal == grp & !is.na(df[[col]])
      n <- sum(mask)
      if (n < 3) next

      vals <- df[[col]][mask]
      decile_mark_summary <- rbind(decile_mark_summary, data.frame(
        mark = m, decile = d,
        group = ifelse(grp, "Neuronal", "Non-neuronal"),
        n = n, median_fold = median(vals),
        mean_fold = mean(vals),
        se = sd(vals) / sqrt(n),
        stringsAsFactors = FALSE
      ))
    }
  }
}

decile_mark_summary$mark <- factor(decile_mark_summary$mark, levels = c("ATAC", "K27ac", "K27me3"))

p_72gd <- ggplot(decile_mark_summary,
                 aes(x = decile, y = median_fold, color = group, shape = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_ribbon(aes(ymin = median_fold - se, ymax = median_fold + se, fill = group),
              alpha = 0.15, color = NA) +
  facet_wrap(~ mark, scales = "free_y") +
  scale_x_continuous(breaks = 1:10, labels = paste0("D", 1:10)) +
  scale_color_manual(values = c("Neuronal" = "#756BB1", "Non-neuronal" = "grey50"),
                     name = "Gene class") +
  scale_fill_manual(values = c("Neuronal" = "#756BB1", "Non-neuronal" = "grey50"),
                    guide = "none") +
  scale_shape_manual(values = c("Neuronal" = 16, "Non-neuronal" = 1), guide = "none") +
  labs(
    title = "Chromatin mark change across K119ub signal deciles",
    subtitle = "Median DiffBind log2FC per decile | D1 = lowest K119ub, D10 = highest",
    x = "K119ub ctrl signal decile", y = "Median log2FC (mut/ctrl)"
  ) +
  theme_biomodal() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 7))

save_plot(p_72gd, "72gd_dose_response_decile_marks", w = 14, h = 6)

# Print D10 vs D1 for neuronal genes
cat("\n  Neuronal gene mark changes: D10 (highest K119ub) vs D1 (lowest):\n")
for (m in c("ATAC", "K27ac", "K27me3")) {
  d1 <- decile_mark_summary[decile_mark_summary$mark == m &
                             decile_mark_summary$group == "Neuronal" &
                             decile_mark_summary$decile == 1, ]
  d10 <- decile_mark_summary[decile_mark_summary$mark == m &
                              decile_mark_summary$group == "Neuronal" &
                              decile_mark_summary$decile == 10, ]
  if (nrow(d1) > 0 && nrow(d10) > 0) {
    cat(sprintf("    %s: D1=%+.4f (n=%d) -> D10=%+.4f (n=%d)\n",
                m, d1$median_fold, d1$n, d10$median_fold, d10$n))
  }
}

# =============================================================================
# 72g-e: INTERACTION MODELS — mark_fold ~ k119ub_signal * is_neuronal
# =============================================================================

cat("\n--- 72g-e: Interaction models ---\n")

interaction_results <- data.frame()
for (m in c("ATAC", "K27ac", "K27me3")) {
  col <- switch(m, ATAC = "atac_fold", K27ac = "k27ac_fold", K27me3 = "k27me3_fold")

  model_df <- df[!is.na(df[[col]]), ]
  if (nrow(model_df) < 50) next

  # Additive model
  lm_add <- lm(as.formula(paste0(col, " ~ gb_ctrl_signal + is_neuronal")), data = model_df)
  # Interaction model
  lm_int <- lm(as.formula(paste0(col, " ~ gb_ctrl_signal * is_neuronal")), data = model_df)

  s_add <- summary(lm_add)
  s_int <- summary(lm_int)
  anova_p <- anova(lm_add, lm_int)$`Pr(>F)`[2]

  int_coef <- s_int$coefficients
  int_row <- int_coef[grep(":", rownames(int_coef)), , drop = FALSE]

  cat(sprintf("\n  %s ~ k119ub_ctrl_signal * is_neuronal (n=%d):\n", m, nrow(model_df)))
  cat(sprintf("    Additive R2 = %.4f\n", s_add$r.squared))
  cat(sprintf("    Interaction R2 = %.4f (ANOVA %s)\n", s_int$r.squared, fmt_p(anova_p)))

  for (j in seq_len(nrow(int_coef))) {
    cat(sprintf("    %-35s beta=%+.5f  SE=%.5f  %s\n",
                rownames(int_coef)[j], int_coef[j, 1], int_coef[j, 2],
                fmt_p(int_coef[j, 4])))
  }

  interaction_results <- rbind(interaction_results, data.frame(
    mark = m, n = nrow(model_df),
    r2_additive = s_add$r.squared,
    r2_interaction = s_int$r.squared,
    anova_p = anova_p,
    k119ub_beta = s_int$coefficients["gb_ctrl_signal", 1],
    k119ub_p = s_int$coefficients["gb_ctrl_signal", 4],
    neuronal_beta = s_int$coefficients["is_neuronalTRUE", 1],
    neuronal_p = s_int$coefficients["is_neuronalTRUE", 4],
    interaction_beta = if (nrow(int_row) > 0) int_row[1, 1] else NA,
    interaction_p = if (nrow(int_row) > 0) int_row[1, 4] else NA,
    stringsAsFactors = FALSE
  ))
}

write.table(interaction_results, file.path(TABLES_DIR, "72g_interaction_models.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n  Saved: 72g_interaction_models.tsv\n")

# Interaction coefficient forest plot
if (nrow(interaction_results) > 0) {
  int_plot_df <- interaction_results[!is.na(interaction_results$interaction_beta), ]
  int_plot_df$mark <- factor(int_plot_df$mark, levels = c("ATAC", "K27ac", "K27me3"))
  int_plot_df$sig <- ifelse(int_plot_df$interaction_p < 0.05, "Significant", "NS")
  int_plot_df$label <- sapply(int_plot_df$interaction_p, fmt_p)

  p_72ge <- ggplot(int_plot_df, aes(x = interaction_beta, y = mark)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = sig), size = 5, show.legend = FALSE) +
    scale_color_manual(values = c("Significant" = "#D73027", "NS" = "grey50")) +
    geom_text(aes(label = sprintf("beta=%+.4f\n%s", interaction_beta, label)),
              hjust = -0.15, size = 3.5) +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.5))) +
    labs(
      title = "Interaction: K119ub signal x neuronal status -> chromatin change",
      subtitle = "Positive beta = neuronal genes respond MORE to K119ub (stronger remodeling)",
      x = "Interaction coefficient (K119ub x Neuronal)", y = NULL
    ) +
    theme_biomodal()

  save_plot(p_72ge, "72ge_interaction_coefficients", w = 10, h = 5)
}

# =============================================================================
# SAVE TABLES
# =============================================================================

cat("\n--- Saving tables ---\n")

write.table(mark_stats, file.path(TABLES_DIR, "72g_mark_stats_neuronal_vs_other.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 72g_mark_stats_neuronal_vs_other.tsv\n")

write.table(quad_stats, file.path(TABLES_DIR, "72g_4group_stats.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 72g_4group_stats.tsv\n")

write.table(decile_mark_summary, file.path(TABLES_DIR, "72g_decile_mark_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 72g_decile_mark_summary.tsv\n")

if (nrow(chromatin_fisher) > 0) {
  write.table(chromatin_fisher, file.path(TABLES_DIR, "72g_chromatin_state_fisher.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  Saved: 72g_chromatin_state_fisher.tsv\n")
}

# =============================================================================
# COMPOSITE
# =============================================================================

cat("\n--- Composite ---\n")

p_composite <- (p_72ga | p_72gc) / (p_72gb) / (p_72gd) +
  plot_layout(heights = c(1, 1, 0.8)) +
  plot_annotation(
    title = "Section 72g: Chromatin Remodeling at Neuronal K119ub Loci",
    subtitle = "BAP1-KO heterochromatin shift is strongest at neuronal genes with high constitutive K119ub",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

save_plot(p_composite, "72g_composite", w = 18, h = 16)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n================================================================================\n")
cat("SECTION 72g COMPLETE: Neuronal Chromatin Remodeling\n")
cat("================================================================================\n\n")

cat(sprintf("  Universe: %d genes (K119ub quantifiable + DiffBind overlap)\n", nrow(df)))
cat(sprintf("  Neuronal: %d (%.1f%%)\n", sum(df$is_neuronal),
            100 * sum(df$is_neuronal) / nrow(df)))

cat("\n  72g-a: Mark changes, neuronal vs non-neuronal (all genes):\n")
for (i in seq_len(nrow(mark_stats))) {
  r <- mark_stats[i, ]
  direction <- ifelse(r$delta_median > 0, "neuronal MORE positive", "neuronal MORE negative")
  cat(sprintf("    %s: delta_median=%+.4f (%s)  %s\n",
              r$mark, r$delta_median, direction, fmt_p(r$wilcox_p)))
}

cat("\n  72g-b: K119ub-High neuronal genes (key test):\n")
for (i in seq_len(nrow(quad_stats))) {
  r <- quad_stats[i, ]
  if (r$k119ub_level == "K119ub High") {
    cat(sprintf("    %s at K119ub-High: neuronal=%+.4f vs other=%+.4f  %s\n",
                r$mark, r$median_neuronal, r$median_other, fmt_p(r$wilcox_p)))
  }
}

cat("\n  72g-e: Interaction tests (K119ub x neuronal -> mark change):\n")
for (i in seq_len(nrow(interaction_results))) {
  r <- interaction_results[i, ]
  sig <- ifelse(!is.na(r$interaction_p) && r$interaction_p < 0.05, "SIGNIFICANT", "ns")
  cat(sprintf("    %s: interaction beta=%+.5f  %s  %s\n",
              r$mark,
              ifelse(is.na(r$interaction_beta), 0, r$interaction_beta),
              fmt_p(ifelse(is.na(r$interaction_p), 1, r$interaction_p)),
              sig))
}

cat("\n  THEORY TEST:\n")
cat("    If neuronal K119ub-high genes show K27ac down + K27me3 up + ATAC down\n")
cat("    more than other genes, it confirms: BAP1 loss -> K119ub at neuronal loci\n")
cat("    -> heterochromatin shift -> MeCP2 responds to chromatin, not methylation.\n")

cat("\n  Plots saved to:", SEC72G_DIR, "\n")
cat("  Tables saved to:", TABLES_DIR, "\n")
cat("\nSection 72g complete.\n")
