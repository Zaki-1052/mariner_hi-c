# biomodal/downstream/scripts/viz_sections/section_64_global_methylation_levels.R
# Section 64: Global CpG Methylation Levels — BAP1-KO Mutant vs Wildtype Control
# Standalone script - sources shared config for all dependencies and data
#
# Plots genome-wide autosomal CpG methylation from the DUET upstream summary:
#   Panel 64a: Per-sample dot+bar by modality (5mC, 5hmC, modC)
#   Panel 64b: Matched-replicate delta (mut - ctrl)
#   Panel 64c: Composition stacked bar (5mC + 5hmC + unmethylated = 100%)
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_64_global_methylation_levels.R

source("scripts/viz_sections/_shared_config.R")

# =============================================================================
# SECTION 64: GLOBAL CpG METHYLATION LEVELS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 64: GLOBAL CpG METHYLATION LEVELS\n")
cat("================================================================================\n\n")

if (is.null(upstream)) {
  stop("Upstream summary CSV not loaded. Cannot proceed.")
}

# ---- Data preparation --------------------------------------------------------

cat("Preparing global methylation data...\n")

meth_data <- upstream %>%
  dplyr::select(
    sample_id,
    mc  = modality_summary_cg_autosomes_mc,
    hmc = modality_summary_cg_autosomes_hmc,
    modc = modality_summary_cg_autosomes_modc,
    unmethylated = modality_summary_cg_autosomes_c
  ) %>%
  dplyr::mutate(
    short_name = gsub("evoC-Bap1-", "", sample_id),
    condition = ifelse(grepl("ctrl", sample_id), "Control", "Mutant")
  ) %>%
  left_join(
    SAMPLES %>% dplyr::select(sample_id, sex, batch),
    by = "sample_id"
  )

cat(sprintf("  %d samples: %d Control, %d Mutant\n",
            nrow(meth_data),
            sum(meth_data$condition == "Control"),
            sum(meth_data$condition == "Mutant")))

# Pivot to long format for faceted plotting
meth_long <- meth_data %>%
  tidyr::pivot_longer(
    cols = c(mc, hmc, modc),
    names_to = "modality_key",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    modality = factor(
      dplyr::recode(modality_key, mc = "5mC", hmc = "5hmC", modc = "modC"),
      levels = c("5mC", "5hmC", "modC")
    )
  )

# Per-condition summary stats
condition_summary <- meth_long %>%
  group_by(modality, condition) %>%
  summarise(
    mean_val = mean(value),
    sd_val = sd(value),
    se_val = sd(value) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

cat("\nPer-condition means:\n")
for (mod in levels(meth_long$modality)) {
  ctrl_mean <- condition_summary$mean_val[condition_summary$modality == mod & condition_summary$condition == "Control"]
  mut_mean <- condition_summary$mean_val[condition_summary$modality == mod & condition_summary$condition == "Mutant"]
  cat(sprintf("  %s: Control = %.3f%%  Mutant = %.3f%%  Delta = %+.3f%%\n",
              mod, ctrl_mean, mut_mean, mut_mean - ctrl_mean))
}

# ---- Statistics --------------------------------------------------------------

cat("\nRunning statistical tests (paired Wilcoxon, matched by sex+batch)...\n")

# Match ctrl-mut pairs by sex and batch
ctrl_samples <- meth_data %>% dplyr::filter(condition == "Control") %>% arrange(sex, batch)
mut_samples  <- meth_data %>% dplyr::filter(condition == "Mutant") %>% arrange(sex, batch)

stats_results <- data.frame(
  modality = character(),
  ctrl_mean = numeric(), ctrl_sd = numeric(),
  mut_mean = numeric(), mut_sd = numeric(),
  delta_mean = numeric(),
  wilcox_p = numeric(), wilcox_V = numeric(),
  ttest_p = numeric(),
  cohens_d = numeric(),
  stringsAsFactors = FALSE
)

modality_cols <- c(mc = "mc", hmc = "hmc", modc = "modc")
modality_labels <- c(mc = "5mC", hmc = "5hmC", modc = "modC")

for (key in names(modality_cols)) {
  ctrl_vals <- ctrl_samples[[modality_cols[key]]]
  mut_vals  <- mut_samples[[modality_cols[key]]]

  wt <- wilcox.test(ctrl_vals, mut_vals, paired = TRUE)
  tt <- t.test(ctrl_vals, mut_vals, paired = TRUE)

  pooled_sd <- sqrt((sd(ctrl_vals)^2 + sd(mut_vals)^2) / 2)
  d <- (mean(mut_vals) - mean(ctrl_vals)) / pooled_sd

  stats_results <- rbind(stats_results, data.frame(
    modality = modality_labels[key],
    ctrl_mean = mean(ctrl_vals), ctrl_sd = sd(ctrl_vals),
    mut_mean = mean(mut_vals), mut_sd = sd(mut_vals),
    delta_mean = mean(mut_vals) - mean(ctrl_vals),
    wilcox_p = wt$p.value, wilcox_V = as.numeric(wt$statistic),
    ttest_p = tt$p.value,
    cohens_d = d,
    stringsAsFactors = FALSE
  ))

  cat(sprintf("  %s: delta = %+.3f%%, Wilcoxon p = %.4f, t-test p = %.4f, Cohen's d = %.2f\n",
              modality_labels[key],
              mean(mut_vals) - mean(ctrl_vals),
              wt$p.value, tt$p.value, d))
}

# ---- Panel 64a: Per-sample dot+bar, faceted by modality ---------------------

cat("\nCreating Figure 64a: per-sample dot+bar plot...\n")

# Y-axis zoom ranges per modality
zoom_ranges <- list(
  "5mC"  = c(71.5, 73.0),
  "5hmC" = c(9.4, 10.5),
  "modC" = c(82.0, 83.0)
)

# Build per-facet p-value annotation df
p_annot <- stats_results %>%
  dplyr::mutate(
    modality = factor(modality, levels = c("5mC", "5hmC", "modC")),
    label = ifelse(wilcox_p < 0.001,
                   sprintf("p = %.2e", wilcox_p),
                   sprintf("p = %.3f", wilcox_p)),
    y_pos = sapply(as.character(modality), function(m) zoom_ranges[[m]][2] * 0.995)
  )

p_64a <- ggplot(meth_long, aes(x = condition, y = value)) +
  geom_bar(data = condition_summary,
           aes(x = condition, y = mean_val, fill = condition),
           stat = "identity", width = 0.6, alpha = 0.7) +
  geom_errorbar(data = condition_summary,
                aes(x = condition, y = mean_val,
                    ymin = mean_val - se_val, ymax = mean_val + se_val),
                width = 0.2, linewidth = 0.6) +
  geom_jitter(aes(shape = sex), width = 0.12, size = 2.5, alpha = 0.9) +
  scale_fill_manual(values = COLORS$condition, name = "Condition") +
  scale_shape_manual(values = c("Female" = 16, "Male" = 17), name = "Sex") +
  facet_wrap(~ modality, scales = "free_y") +
  labs(
    title = "Genome-Wide CpG Methylation: BAP1-KO vs Wildtype",
    subtitle = "DUET evoC autosomal CpG (8 samples, paired Wilcoxon test)",
    x = "", y = "Methylation (%)"
  ) +
  theme_biomodal() +
  theme(legend.position = "right",
        strip.text = element_text(size = 12, face = "bold"))

# Apply zoomed y-axes per facet via ggh4x if available, otherwise use coord_cartesian fallback
# coord_cartesian doesn't support per-facet limits, so we use scale_y_continuous with oob
p_64a_list <- list()
for (mod in c("5mC", "5hmC", "modC")) {
  mod_data <- meth_long %>% dplyr::filter(modality == mod)
  mod_summary <- condition_summary %>% dplyr::filter(modality == mod)
  mod_stats <- stats_results %>% dplyr::filter(modality == mod)

  p_label <- ifelse(mod_stats$wilcox_p < 0.001,
                    sprintf("p = %.2e", mod_stats$wilcox_p),
                    sprintf("p = %.3f", mod_stats$wilcox_p))

  p_facet <- ggplot(mod_data, aes(x = condition, y = value)) +
    geom_bar(data = mod_summary,
             aes(x = condition, y = mean_val, fill = condition),
             stat = "identity", width = 0.6, alpha = 0.7) +
    geom_errorbar(data = mod_summary,
                  aes(x = condition, y = mean_val,
                      ymin = mean_val - se_val, ymax = mean_val + se_val),
                  width = 0.2, linewidth = 0.6) +
    geom_jitter(aes(shape = sex), width = 0.12, size = 3, alpha = 0.9) +
    annotate("text", x = 1.5, y = zoom_ranges[[mod]][2] * 0.998,
             label = p_label, size = 3.5, fontface = "italic") +
    annotate("segment", x = 1, xend = 2,
             y = zoom_ranges[[mod]][2] * 0.99,
             yend = zoom_ranges[[mod]][2] * 0.99,
             linewidth = 0.4) +
    scale_fill_manual(values = COLORS$condition, guide = "none") +
    scale_shape_manual(values = c("Female" = 16, "Male" = 17), guide = "none") +
    coord_cartesian(ylim = zoom_ranges[[mod]]) +
    labs(title = mod, x = "", y = "Methylation (%)") +
    theme_biomodal() +
    theme(plot.title = element_text(size = 14, face = "bold"))

  p_64a_list[[mod]] <- p_facet
}

# Add legend manually from one panel
legend_plot <- ggplot(meth_long, aes(x = condition, y = value, fill = condition, shape = sex)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.6, alpha = 0.7) +
  geom_point(size = 3) +
  scale_fill_manual(values = COLORS$condition, name = "Condition") +
  scale_shape_manual(values = c("Female" = 16, "Male" = 17), name = "Sex") +
  theme_biomodal()

p_64a_combined <- (p_64a_list[["5mC"]] | p_64a_list[["5hmC"]] | p_64a_list[["modC"]]) +
  plot_annotation(
    title = "Genome-Wide CpG Methylation: BAP1-KO vs Wildtype",
    subtitle = "DUET evoC autosomal CpG | Bars = condition mean +/- SE | Dots = individual samples",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  )

save_multiformat_ggplot(p_64a_combined,
                        file.path(OUTPUT_DIR, "64a_global_methylation_dotbar"),
                        width = 14, height = 6)

# ---- Panel 64b: Matched-replicate delta plot --------------------------------

cat("Creating Figure 64b: matched-replicate delta plot...\n")

delta_data <- data.frame(
  pair = c("F-B1", "M-B1", "F-B2", "M-B2"),
  stringsAsFactors = FALSE
)

for (key in names(modality_cols)) {
  col <- modality_cols[key]
  delta_data[[modality_labels[key]]] <- mut_samples[[col]] - ctrl_samples[[col]]
}

delta_long <- delta_data %>%
  tidyr::pivot_longer(
    cols = c("5mC", "5hmC", "modC"),
    names_to = "modality",
    values_to = "delta"
  ) %>%
  dplyr::mutate(
    modality = factor(modality, levels = c("5mC", "5hmC", "modC")),
    direction = ifelse(delta > 0, "Increased", "Decreased")
  )

p_64b <- ggplot(delta_long, aes(x = pair, y = delta, color = direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_segment(aes(xend = pair, y = 0, yend = delta), linewidth = 1.5) +
  geom_point(size = 4) +
  geom_text(aes(label = sprintf("%+.2f", delta)),
            vjust = ifelse(delta_long$delta > 0, -1, 1.5), size = 3) +
  scale_color_manual(values = c("Increased" = "#B2182B", "Decreased" = "#2166AC"),
                     name = "Direction\n(Mut - Ctrl)") +
  facet_wrap(~ modality, scales = "free_y") +
  labs(
    title = "Per-Replicate Methylation Change (Mutant - Control)",
    subtitle = "Matched by sex and batch | 5mC consistently up, 5hmC consistently down",
    x = "Replicate Pair", y = "Delta (%)"
  ) +
  theme_biomodal() +
  theme(strip.text = element_text(size = 12, face = "bold"))

save_multiformat_ggplot(p_64b,
                        file.path(OUTPUT_DIR, "64b_replicate_delta"),
                        width = 12, height = 6)

# ---- Panel 64c: Composition stacked bar -------------------------------------

cat("Creating Figure 64c: composition stacked bar...\n")

comp_data <- meth_data %>%
  group_by(condition) %>%
  summarise(
    `5mC` = mean(mc),
    `5hmC` = mean(hmc),
    Unmethylated = mean(unmethylated),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = c("5mC", "5hmC", "Unmethylated"),
    names_to = "component",
    values_to = "pct"
  ) %>%
  dplyr::mutate(
    component = factor(component, levels = c("Unmethylated", "5hmC", "5mC"))
  )

comp_colors <- c("5mC" = COLORS$methylation[["5mC"]],
                 "5hmC" = COLORS$methylation[["5hmC"]],
                 "Unmethylated" = "grey70")

p_64c <- ggplot(comp_data, aes(x = condition, y = pct, fill = component)) +
  geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_stack(vjust = 0.5), size = 3.5,
            color = ifelse(comp_data$component == "5mC", "white", "black")) +
  scale_fill_manual(values = comp_colors, name = "Component") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title = "CpG Methylation Composition",
    subtitle = "Autosomal mean | TET blockade: 5mC up, 5hmC down, total stable",
    x = "", y = "Fraction of CpG Sites (%)"
  ) +
  theme_biomodal() +
  theme(legend.position = "right")

save_multiformat_ggplot(p_64c,
                        file.path(OUTPUT_DIR, "64c_composition_stacked"),
                        width = 7, height = 7)

# ---- Combined figure ---------------------------------------------------------

cat("Creating combined Figure 64...\n")

p_64_combined <- p_64a_combined /
  (p_64b | p_64c) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Global CpG Methylation Levels: BAP1-KO vs Wildtype",
    subtitle = "Consistent with TET-pathway blockade: 5mC↑, 5hmC↓, total modC unchanged",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11)
    )
  )

save_multiformat_ggplot(p_64_combined,
                        file.path(OUTPUT_DIR, "64_global_methylation_levels"),
                        width = 16, height = 12)

# ---- Save tables -------------------------------------------------------------

cat("\nSaving summary tables...\n")

# Per-sample summary
sample_summary <- meth_data %>%
  dplyr::select(sample_id, condition, sex, batch, mc, hmc, modc, unmethylated) %>%
  dplyr::rename(`5mC_pct` = mc, `5hmC_pct` = hmc, `modC_pct` = modc,
                unmethylated_pct = unmethylated)

write.table(sample_summary,
            file.path(TABLES_DIR, "64_global_methylation_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 64_global_methylation_summary.tsv\n")

write.table(stats_results,
            file.path(TABLES_DIR, "64_statistics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Saved: 64_statistics.tsv\n")

cat("\n")
cat("================================================================================\n")
cat("SECTION 64 COMPLETE\n")
cat("================================================================================\n")
