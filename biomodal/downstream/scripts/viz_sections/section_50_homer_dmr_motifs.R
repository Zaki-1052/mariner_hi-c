# biomodal/downstream/scripts/viz_sections/section_50_homer_dmr_motifs.R
# Section 50: HOMER Motif Enrichment — DMR comparisons (A5 only)
# A1-A4 (genome background) produced no/artefactual signal.
# A5 (coordinated mC-up/hmC-down vs discordant) has genuine enrichment.

source("scripts/viz_sections/_shared_config.R")
source("../../scripts/utils/multi_format_output.R")

cat("================================================================================\n")
cat("SECTION 50: HOMER MOTIF ENRICHMENT (DMRs)\n")
cat("================================================================================\n\n")

HOMER_DIR <- file.path(BASE_DIR, "results/homer_motif_enrichment")
SECTION_DIR <- file.path(OUTPUT_DIR, "50_homer_dmr_motifs")
dir.create(SECTION_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Load A5 results
# =============================================================================

cat("Loading A5 (coordinated vs discordant) knownResults...\n")

fpath <- file.path(HOMER_DIR, "A5_coordinated_vs_discordant", "knownResults.txt")

a5 <- read.delim(fpath, sep = "\t", header = TRUE, check.names = FALSE,
                 stringsAsFactors = FALSE)
colnames(a5) <- c("motif_name", "consensus", "pvalue", "log_pvalue",
                   "qvalue", "target_count", "target_pct",
                   "bg_count", "bg_pct")

a5$target_pct <- as.numeric(gsub("%", "", a5$target_pct))
a5$bg_pct     <- as.numeric(gsub("%", "", a5$bg_pct))
a5$pvalue     <- as.numeric(a5$pvalue)
a5$qvalue     <- as.numeric(a5$qvalue)

a5$tf_name <- trimws(sub("\\(.*", "", a5$motif_name))
a5$family  <- sub(".*?\\(([^)]+)\\).*", "\\1", a5$motif_name)
a5$family[!grepl("\\(", a5$motif_name)] <- "Other"

a5$fold_enrichment <- ifelse(a5$bg_pct > 0.05, a5$target_pct / a5$bg_pct, NA_real_)
a5$neg_log10_p <- -log10(pmax(a5$pvalue, 1e-300))

n_sig <- sum(a5$qvalue < 0.05)
cat(sprintf("  Total motifs: %d, Significant (q < 0.05): %d\n", nrow(a5), n_sig))

# =============================================================================
# Plot 1: A5 dot plot — top 25 motifs
# =============================================================================

cat("Plot 1: A5 dot plot (top 25 significant motifs)...\n")

dot_data <- a5 %>%
  filter(qvalue < 0.05) %>%
  group_by(tf_name) %>%
  slice_min(order_by = pvalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  slice_min(order_by = pvalue, n = 25, with_ties = FALSE) %>%
  mutate(tf_name = factor(tf_name, levels = rev(tf_name[order(pvalue)])))

p_a5 <- ggplot(dot_data, aes(x = neg_log10_p, y = tf_name)) +
  geom_point(aes(size = target_pct, fill = family),
             shape = 21, color = "black", stroke = 0.3) +
  scale_size_continuous(range = c(3, 8), name = "% Target\nwith Motif") +
  scale_fill_brewer(palette = "Set2", name = "TF Family") +
  labs(
    x = expression(-log[10](p-value)),
    y = NULL,
    title = "DMR Motif Enrichment: Coordinated (mC up/hmC down) vs Discordant",
    subtitle = "A5 comparison, q < 0.05. A1-A4 (genome background) produced no significant motifs."
  ) +
  theme_biomodal() +
  theme(axis.text.y = element_text(size = 10))

save_multiformat_ggplot(p_a5, file.path(SECTION_DIR, "homer_a5_coordinated_dotplot"),
                        width = 12, height = 10)
cat("  Saved: homer_a5_coordinated_dotplot\n\n")

# =============================================================================
# Plot 2: A5 family breakdown
# =============================================================================

cat("Plot 2: A5 TF family counts...\n")

family_counts <- a5 %>%
  filter(qvalue < 0.05) %>%
  count(family, name = "n_motifs") %>%
  arrange(desc(n_motifs)) %>%
  head(10) %>%
  mutate(family = factor(family, levels = rev(family)))

p_fam <- ggplot(family_counts, aes(x = n_motifs, y = family)) +
  geom_col(fill = "#66C2A5", width = 0.7) +
  geom_text(aes(label = n_motifs), hjust = -0.3, size = 4) +
  labs(
    x = "Significant Motifs (q < 0.05)",
    y = NULL,
    title = "TF Families Enriched in Coordinated DMRs (A5)",
    subtitle = "Zinc-finger (Zf) family dominates coordinated mC-up/hmC-down regions"
  ) +
  theme_biomodal() +
  expand_limits(x = max(family_counts$n_motifs) * 1.15)

save_multiformat_ggplot(p_fam, file.path(SECTION_DIR, "homer_a5_family_counts"),
                        width = 10, height = 6)
cat("  Saved: homer_a5_family_counts\n\n")

# =============================================================================
# Table exports
# =============================================================================

cat("Exporting tables...\n")

sig_table <- a5 %>%
  filter(qvalue < 0.05) %>%
  dplyr::select(tf_name, family, consensus, pvalue, qvalue,
                target_pct, bg_pct, fold_enrichment) %>%
  arrange(pvalue)

write_tsv(sig_table, file.path(SECTION_DIR, "homer_a5_significant_motifs.tsv"))
cat(sprintf("  Saved: homer_a5_significant_motifs.tsv (%d rows)\n", nrow(sig_table)))

cat("\n================================================================================\n")
cat("Section 50 complete.\n")
cat(sprintf("Output directory: %s\n", SECTION_DIR))
cat("================================================================================\n")
