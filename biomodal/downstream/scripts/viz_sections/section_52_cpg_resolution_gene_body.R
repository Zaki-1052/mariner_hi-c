# biomodal/downstream/scripts/viz_sections/section_52_cpg_resolution_gene_body.R
# Section 52: Single-CpG Resolution Gene Body Analysis
# Analyzes sub-gene-feature methylation patterns from modality regional-frac
# extraction on custom BED files (exons, introns, UTRs, splice sites).
#
# Prerequisites:
#   - Run prepare_section52_beds.R → commit/push → HPC pull → sbatch
#   - rsync outputs back: modality/outputs/outputs_CG_section52/
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_52_cpg_resolution_gene_body.R

source("scripts/viz_sections/_shared_config.R")

library(dunn.test)

# =============================================================================
# SECTION 52 CONFIGURATION
# =============================================================================

cat("\n")
cat("========================================================================\n")
cat("SECTION 52: SINGLE-CpG RESOLUTION GENE BODY ANALYSIS\n")
cat("========================================================================\n\n")

SEC52_DIR <- file.path(OUTPUT_DIR, "52_cpg_resolution_gene_body")
dir.create(SEC52_DIR, recursive = TRUE, showWarnings = FALSE)

SECTION52_BASE <- file.path(BASE_DIR, "modality/outputs/outputs_CG_section52/Results")

# Auto-discover extract TSV files (timestamp-agnostic)
discover_extract_tsv <- function(region_dir, field) {
  pattern <- paste0("Extract_", field, "_regional-frac_.*\\.tsv\\.gz$")
  search_dir <- file.path(SECTION52_BASE, region_dir)
  stopifnot(dir.exists(search_dir))
  files <- list.files(search_dir, pattern = pattern,
                      recursive = TRUE, full.names = TRUE)
  stopifnot(length(files) == 1)
  files[1]
}

FEATURE_MC_PATH  <- discover_extract_tsv("sub_gene_features.annotation", "mc")
FEATURE_HMC_PATH <- discover_extract_tsv("sub_gene_features.annotation", "hmc")
METAGENE_MC_PATH  <- discover_extract_tsv("metagene_bins.annotation", "mc")
METAGENE_HMC_PATH <- discover_extract_tsv("metagene_bins.annotation", "hmc")

cat(sprintf("  Feature mC:  %s\n", basename(FEATURE_MC_PATH)))
cat(sprintf("  Feature hmC: %s\n", basename(FEATURE_HMC_PATH)))
cat(sprintf("  Metagene mC:  %s\n", basename(METAGENE_MC_PATH)))
cat(sprintf("  Metagene hmC: %s\n", basename(METAGENE_HMC_PATH)))

# Strand lookup for metagene orientation
STRAND_PATH <- file.path(BASE_DIR, "modality/section52_beds/gene_strand_lookup.tsv")
stopifnot(file.exists(STRAND_PATH))

# Feature type order for plots (5' to 3' gene structure)
FEATURE_ORDER <- c("5UTR", "Exon", "SpliceSite_Donor", "SpliceSite_Acceptor",
                   "Intron", "3UTR")

FEATURE_COLORS <- c(
  "5UTR" = "#66C2A5",
  "Exon" = "#FC8D62",
  "SpliceSite_Donor" = "#8DA0CB",
  "SpliceSite_Acceptor" = "#A6D854",
  "Intron" = "#E78AC3",
  "3UTR" = "#FFD92F"
)

# =============================================================================
# LOAD DATA
# =============================================================================

cat("\nLoading extract TSV files...\n")

load_extract <- function(path) {
  df <- read.table(gzfile(path), header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, na.strings = c("", "NA"))
  cat(sprintf("  %s: %d rows x %d cols\n", basename(path), nrow(df), ncol(df)))
  df
}

feat_mc  <- load_extract(FEATURE_MC_PATH)
feat_hmc <- load_extract(FEATURE_HMC_PATH)
meta_mc  <- load_extract(METAGENE_MC_PATH)
meta_hmc <- load_extract(METAGENE_HMC_PATH)

strand_lookup <- read.table(STRAND_PATH, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE)

# =============================================================================
# IDENTIFY SAMPLE COLUMNS
# =============================================================================

# Sample columns are between Start/End and Annotation/Name
# Pattern: columns containing "ctrl" or "mut"
identify_sample_cols <- function(df) {
  cols <- colnames(df)
  sample_idx <- grep("ctrl|mut", cols, ignore.case = TRUE)
  cols[sample_idx]
}

sample_cols <- identify_sample_cols(feat_mc)
ctrl_cols <- sample_cols[grep("ctrl", sample_cols, ignore.case = TRUE)]
mut_cols <- sample_cols[grep("mut", sample_cols, ignore.case = TRUE)]

cat(sprintf("  %d ctrl columns, %d mut columns\n",
            length(ctrl_cols), length(mut_cols)))
stopifnot(length(ctrl_cols) >= 2)
stopifnot(length(mut_cols) >= 2)

# =============================================================================
# COMPUTE DELTAS FOR SUB-GENE FEATURES
# =============================================================================

cat("\nComputing per-feature methylation deltas...\n")

compute_deltas <- function(mc_df, hmc_df, ctrl_cols, mut_cols) {
  mc_df$ctrl_mc <- rowMeans(mc_df[, ctrl_cols, drop = FALSE], na.rm = FALSE)
  mc_df$mut_mc <- rowMeans(mc_df[, mut_cols, drop = FALSE], na.rm = FALSE)
  mc_df$delta_mc <- mc_df$mut_mc - mc_df$ctrl_mc

  hmc_df$ctrl_hmc <- rowMeans(hmc_df[, ctrl_cols, drop = FALSE], na.rm = FALSE)
  hmc_df$mut_hmc <- rowMeans(hmc_df[, mut_cols, drop = FALSE], na.rm = FALSE)
  hmc_df$delta_hmc <- hmc_df$mut_hmc - hmc_df$ctrl_hmc

  # Note: hmc columns have "hmc" in the name; ctrl/mut cols differ
  # We join on Name which is the unique feature identifier
  mc_out <- mc_df[, c("Chromosome", "Start", "End", "Annotation", "Name",
                       "ctrl_mc", "mut_mc", "delta_mc")]
  hmc_out <- hmc_df[, c("Name", "ctrl_hmc", "mut_hmc", "delta_hmc")]

  merged <- merge(mc_out, hmc_out, by = "Name")
  merged
}

# The hmC extract has different column names (hmc instead of mc in suffix)
hmc_sample_cols <- identify_sample_cols(feat_hmc)
hmc_ctrl_cols <- hmc_sample_cols[grep("ctrl", hmc_sample_cols, ignore.case = TRUE)]
hmc_mut_cols <- hmc_sample_cols[grep("mut", hmc_sample_cols, ignore.case = TRUE)]

feat_mc$ctrl_mc <- rowMeans(feat_mc[, ctrl_cols, drop = FALSE], na.rm = FALSE)
feat_mc$mut_mc <- rowMeans(feat_mc[, mut_cols, drop = FALSE], na.rm = FALSE)
feat_mc$delta_mc <- feat_mc$mut_mc - feat_mc$ctrl_mc

feat_hmc$ctrl_hmc <- rowMeans(feat_hmc[, hmc_ctrl_cols, drop = FALSE], na.rm = FALSE)
feat_hmc$mut_hmc <- rowMeans(feat_hmc[, hmc_mut_cols, drop = FALSE], na.rm = FALSE)
feat_hmc$delta_hmc <- feat_hmc$mut_hmc - feat_hmc$ctrl_hmc

mc_slim <- feat_mc[, c("Chromosome", "Start", "End", "Annotation", "Name",
                        "ctrl_mc", "mut_mc", "delta_mc")]
hmc_slim <- feat_hmc[, c("Name", "ctrl_hmc", "mut_hmc", "delta_hmc")]

feat_data <- merge(mc_slim, hmc_slim, by = "Name")

# Parse gene name and feature type from Name column
# Name format: {Gene}_{FeatureType}_{rank}
# FeatureType can be multi-part: SpliceSite_Donor, SpliceSite_Acceptor
# Strategy: match the known feature types and extract gene as prefix
parse_name <- function(name_col) {
  gene <- character(length(name_col))
  feature <- character(length(name_col))

  for (ft in c("SpliceSite_Donor", "SpliceSite_Acceptor",
               "5UTR", "3UTR", "Exon", "Intron")) {
    pattern <- paste0("^(.+)_", ft, "_[0-9]+$")
    matches <- grepl(pattern, name_col)
    gene[matches] <- sub(pattern, "\\1", name_col[matches])
    feature[matches] <- ft
  }

  data.frame(gene = gene, feature_type = feature, stringsAsFactors = FALSE)
}

parsed <- parse_name(feat_data$Name)
feat_data$gene <- parsed$gene
feat_data$feature_type <- parsed$feature_type

# Filter to known feature types and set factor order
feat_data <- feat_data[feat_data$feature_type %in% FEATURE_ORDER, ]
feat_data$feature_type <- factor(feat_data$feature_type, levels = FEATURE_ORDER)

# Remove rows with NA deltas
feat_data <- feat_data[complete.cases(feat_data[, c("delta_mc", "delta_hmc")]), ]

cat(sprintf("  %d feature intervals with complete data across %d genes.\n",
            nrow(feat_data), length(unique(feat_data$gene))))

# Feature-type breakdown
cat("  Feature-type counts:\n")
print(table(feat_data$feature_type))

# =============================================================================
# FIGURE 52a: Feature Distribution Overview
# =============================================================================

cat("\n--- Figure 52a: Feature distribution overview ---\n")

fig_52a_dir <- file.path(SEC52_DIR, "52a_feature_distribution")
dir.create(fig_52a_dir, recursive = TRUE, showWarnings = FALSE)

feat_summary <- feat_data %>%
  dplyr::group_by(feature_type) %>%
  dplyr::summarise(n_regions = dplyr::n(),
                   n_genes = dplyr::n_distinct(gene),
                   .groups = "drop")

p_52a <- ggplot(feat_summary, aes(x = feature_type, y = n_regions,
                                   fill = feature_type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_regions), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = FEATURE_COLORS, guide = "none") +
  labs(title = "Sub-Gene Feature Distribution",
       subtitle = sprintf("%d features across %d max-significance genes",
                          nrow(feat_data), length(unique(feat_data$gene))),
       x = "Gene Feature", y = "Number of Regions") +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_multiformat_ggplot(
  p_52a, file.path(fig_52a_dir, "52a_feature_distribution"),
  width = 8, height = 6
)

# =============================================================================
# FIGURE 52b: Delta mC by Feature Type
# =============================================================================

cat("\n--- Figure 52b: Delta mC by feature type ---\n")

fig_52b_dir <- file.path(SEC52_DIR, "52b_delta_mc_by_feature")
dir.create(fig_52b_dir, recursive = TRUE, showWarnings = FALSE)

p_52b <- ggplot(feat_data, aes(x = feature_type, y = delta_mc * 100,
                                fill = feature_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  scale_fill_manual(values = FEATURE_COLORS, guide = "none") +
  labs(title = expression(Delta * "5mC by Gene Feature"),
       subtitle = "Mutant - Control (positive = hypermethylated in KO)",
       x = "Gene Feature", y = expression(Delta * "5mC (%)")) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_multiformat_ggplot(
  p_52b, file.path(fig_52b_dir, "52b_delta_mc_by_feature"),
  width = 9, height = 6
)

# Kruskal-Wallis test for feature type effect
kw_mc <- kruskal.test(delta_mc ~ feature_type, data = feat_data)
cat(sprintf("  Kruskal-Wallis (delta_mc ~ feature): chi-sq=%.2f, p=%s\n",
            kw_mc$statistic, format.pval(kw_mc$p.value, digits = 3)))

cat("\n  Dunn's post-hoc pairwise tests (delta_mC, BH-adjusted):\n")
mc_dunn <- dunn.test(feat_data$delta_mc, feat_data$feature_type,
                     method = "bh", kw = FALSE, table = FALSE, list = TRUE)
mc_posthoc <- data.frame(
  comparison = mc_dunn$comparisons,
  Z = round(mc_dunn$Z, 3),
  p_raw = signif(mc_dunn$P, 3),
  p_adj = signif(mc_dunn$P.adjusted, 3),
  stringsAsFactors = FALSE
)
mc_posthoc <- mc_posthoc[order(mc_posthoc$p_adj), ]
mc_sig <- mc_posthoc[mc_posthoc$p_adj < 0.05, ]
cat(sprintf("  %d / %d pairs significant at q < 0.05:\n",
            nrow(mc_sig), nrow(mc_posthoc)))
if (nrow(mc_sig) > 0) {
  for (i in seq_len(nrow(mc_sig))) {
    cat(sprintf("    %-40s Z=%7.3f  q=%.2e\n",
                mc_sig$comparison[i], mc_sig$Z[i], mc_sig$p_adj[i]))
  }
} else {
  cat("    (none)\n")
}
cat("  Non-significant pairs:\n")
mc_ns <- mc_posthoc[mc_posthoc$p_adj >= 0.05, ]
for (i in seq_len(nrow(mc_ns))) {
  cat(sprintf("    %-40s Z=%7.3f  q=%.3f\n",
              mc_ns$comparison[i], mc_ns$Z[i], mc_ns$p_adj[i]))
}

# =============================================================================
# FIGURE 52c: Delta hmC by Feature Type
# =============================================================================

cat("\n--- Figure 52c: Delta hmC by feature type ---\n")

fig_52c_dir <- file.path(SEC52_DIR, "52c_delta_hmc_by_feature")
dir.create(fig_52c_dir, recursive = TRUE, showWarnings = FALSE)

p_52c <- ggplot(feat_data, aes(x = feature_type, y = delta_hmc * 100,
                                fill = feature_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  scale_fill_manual(values = FEATURE_COLORS, guide = "none") +
  labs(title = expression(Delta * "5hmC by Gene Feature"),
       subtitle = "Mutant - Control (negative = hypo-hydroxymethylated in KO)",
       x = "Gene Feature", y = expression(Delta * "5hmC (%)")) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_multiformat_ggplot(
  p_52c, file.path(fig_52c_dir, "52c_delta_hmc_by_feature"),
  width = 9, height = 6
)

kw_hmc <- kruskal.test(delta_hmc ~ feature_type, data = feat_data)
cat(sprintf("  Kruskal-Wallis (delta_hmc ~ feature): chi-sq=%.2f, p=%s\n",
            kw_hmc$statistic, format.pval(kw_hmc$p.value, digits = 3)))

cat("\n  Dunn's post-hoc pairwise tests (delta_hmC, BH-adjusted):\n")
hmc_dunn <- dunn.test(feat_data$delta_hmc, feat_data$feature_type,
                      method = "bh", kw = FALSE, table = FALSE, list = TRUE)
hmc_posthoc <- data.frame(
  comparison = hmc_dunn$comparisons,
  Z = round(hmc_dunn$Z, 3),
  p_raw = signif(hmc_dunn$P, 3),
  p_adj = signif(hmc_dunn$P.adjusted, 3),
  stringsAsFactors = FALSE
)
hmc_posthoc <- hmc_posthoc[order(hmc_posthoc$p_adj), ]
hmc_sig <- hmc_posthoc[hmc_posthoc$p_adj < 0.05, ]
cat(sprintf("  %d / %d pairs significant at q < 0.05:\n",
            nrow(hmc_sig), nrow(hmc_posthoc)))
if (nrow(hmc_sig) > 0) {
  for (i in seq_len(nrow(hmc_sig))) {
    cat(sprintf("    %-40s Z=%7.3f  q=%.2e\n",
                hmc_sig$comparison[i], hmc_sig$Z[i], hmc_sig$p_adj[i]))
  }
} else {
  cat("    (none)\n")
}
cat("  Non-significant pairs:\n")
hmc_ns <- hmc_posthoc[hmc_posthoc$p_adj >= 0.05, ]
for (i in seq_len(nrow(hmc_ns))) {
  cat(sprintf("    %-40s Z=%7.3f  q=%.3f\n",
              hmc_ns$comparison[i], hmc_ns$Z[i], hmc_ns$p_adj[i]))
}

# =============================================================================
# FIGURE 52d: KEY_GENES Individual Locus Panels
# =============================================================================

cat("\n--- Figure 52d: KEY_GENES locus panels ---\n")

fig_52d_dir <- file.path(SEC52_DIR, "52d_key_gene_loci")
dir.create(fig_52d_dir, recursive = TRUE, showWarnings = FALSE)

key_data <- feat_data[feat_data$gene %in% KEY_GENES, ]
key_data <- key_data[order(key_data$gene, key_data$Start), ]

cat(sprintf("  %d features across %d KEY_GENES found.\n",
            nrow(key_data), length(unique(key_data$gene))))

if (nrow(key_data) > 0) {
  gene_panels <- list()
  genes_present <- intersect(KEY_GENES, unique(key_data$gene))

  for (g in genes_present) {
    gd <- key_data[key_data$gene == g, ]
    gd$midpoint <- (gd$Start + gd$End) / 2

    p_mc <- ggplot(gd, aes(x = midpoint / 1e6, y = delta_mc * 100,
                            color = feature_type)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
                 linewidth = 0.3) +
      geom_point(size = 2.5) +
      geom_segment(aes(xend = midpoint / 1e6, yend = 0), linewidth = 0.5) +
      scale_color_manual(values = FEATURE_COLORS) +
      labs(title = g, x = "Position (Mb)",
           y = expression(Delta * "mC (%)"),
           color = "Feature") +
      theme_biomodal() +
      theme(legend.position = "none",
            plot.title = element_text(size = 11, face = "bold"))

    gene_panels[[g]] <- p_mc
  }

  if (length(gene_panels) >= 2) {
    # Combine with shared legend
    combined <- patchwork::wrap_plots(gene_panels, ncol = 2) +
      patchwork::plot_layout(guides = "collect") +
      patchwork::plot_annotation(
        title = expression("Sub-Gene " * Delta * "5mC by Genomic Position"),
        subtitle = "KEY_GENES: lollipop at each sub-gene feature",
        theme = theme(plot.title = element_text(size = 14, face = "bold"))
      )

    # Add a standalone legend
    legend_data <- data.frame(
      x = 1:length(FEATURE_ORDER),
      y = 1:length(FEATURE_ORDER),
      feature_type = factor(FEATURE_ORDER, levels = FEATURE_ORDER)
    )
    p_legend <- ggplot(legend_data, aes(x, y, color = feature_type)) +
      geom_point(size = 3) +
      scale_color_manual(values = FEATURE_COLORS, name = "Feature") +
      theme_void() +
      theme(legend.position = "right")
    legend_grob <- cowplot::get_legend(p_legend)

    final_52d <- combined + patchwork::inset_element(
      legend_grob, left = 0.85, right = 1, top = 1, bottom = 0.7
    )

    save_multiformat_ggplot(
      combined, file.path(fig_52d_dir, "52d_key_gene_loci"),
      width = 14, height = ceiling(length(gene_panels) / 2) * 3.5
    )
  }
}

# =============================================================================
# FIGURE 52e: Metagene Profile
# =============================================================================

cat("\n--- Figure 52e: Metagene profile ---\n")

fig_52e_dir <- file.path(SEC52_DIR, "52e_metagene_profile")
dir.create(fig_52e_dir, recursive = TRUE, showWarnings = FALSE)

# Compute metagene deltas
meta_mc_sample_cols <- identify_sample_cols(meta_mc)
meta_mc_ctrl <- meta_mc_sample_cols[grep("ctrl", meta_mc_sample_cols, ignore.case = TRUE)]
meta_mc_mut <- meta_mc_sample_cols[grep("mut", meta_mc_sample_cols, ignore.case = TRUE)]

meta_hmc_sample_cols <- identify_sample_cols(meta_hmc)
meta_hmc_ctrl <- meta_hmc_sample_cols[grep("ctrl", meta_hmc_sample_cols, ignore.case = TRUE)]
meta_hmc_mut <- meta_hmc_sample_cols[grep("mut", meta_hmc_sample_cols, ignore.case = TRUE)]

meta_mc$ctrl_mc <- rowMeans(meta_mc[, meta_mc_ctrl, drop = FALSE], na.rm = FALSE)
meta_mc$mut_mc <- rowMeans(meta_mc[, meta_mc_mut, drop = FALSE], na.rm = FALSE)
meta_mc$delta_mc <- meta_mc$mut_mc - meta_mc$ctrl_mc

meta_hmc$ctrl_hmc <- rowMeans(meta_hmc[, meta_hmc_ctrl, drop = FALSE], na.rm = FALSE)
meta_hmc$mut_hmc <- rowMeans(meta_hmc[, meta_hmc_mut, drop = FALSE], na.rm = FALSE)
meta_hmc$delta_hmc <- meta_hmc$mut_hmc - meta_hmc$ctrl_hmc

# Parse gene and bin number from Name
meta_mc$gene <- sub("_bin_[0-9]+$", "", meta_mc$Name)
meta_mc$bin <- as.integer(sub(".*_bin_", "", meta_mc$Name))

meta_hmc$gene <- sub("_bin_[0-9]+$", "", meta_hmc$Name)
meta_hmc$bin <- as.integer(sub(".*_bin_", "", meta_hmc$Name))

# Reverse bin order for minus-strand genes
minus_genes <- strand_lookup$gene[strand_lookup$strand == "-"]

meta_mc$bin_oriented <- meta_mc$bin
meta_mc$bin_oriented[meta_mc$gene %in% minus_genes] <-
  101 - meta_mc$bin[meta_mc$gene %in% minus_genes]

meta_hmc$bin_oriented <- meta_hmc$bin
meta_hmc$bin_oriented[meta_hmc$gene %in% minus_genes] <-
  101 - meta_hmc$bin[meta_hmc$gene %in% minus_genes]

# Aggregate across genes: mean and SE per bin
mc_profile <- meta_mc %>%
  dplyr::filter(!is.na(delta_mc)) %>%
  dplyr::group_by(bin_oriented) %>%
  dplyr::summarise(
    mean_delta = mean(delta_mc, na.rm = TRUE),
    se_delta = sd(delta_mc, na.rm = TRUE) / sqrt(dplyr::n()),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(modification = "5mC")

hmc_profile <- meta_hmc %>%
  dplyr::filter(!is.na(delta_hmc)) %>%
  dplyr::group_by(bin_oriented) %>%
  dplyr::summarise(
    mean_delta = mean(delta_hmc, na.rm = TRUE),
    se_delta = sd(delta_hmc, na.rm = TRUE) / sqrt(dplyr::n()),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(modification = "5hmC")

profile_data <- rbind(mc_profile, hmc_profile)

p_52e <- ggplot(profile_data, aes(x = bin_oriented, y = mean_delta * 100,
                                   color = modification, fill = modification)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = (mean_delta - se_delta) * 100,
                  ymax = (mean_delta + se_delta) * 100),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("5mC" = "#E41A1C", "5hmC" = "#377EB8")) +
  scale_fill_manual(values = c("5mC" = "#E41A1C", "5hmC" = "#377EB8")) +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("5'", "25%", "50%", "75%", "3'")) +
  labs(title = "Metagene Methylation Change Profile",
       subtitle = sprintf("Mean +/- SE across %d max-significance gene bodies",
                          length(unique(meta_mc$gene))),
       x = "Relative Gene Body Position (5' to 3')",
       y = expression(Delta * " Methylation (%)"),
       color = "Modification", fill = "Modification") +
  theme_biomodal() +
  theme(legend.position = c(0.85, 0.85))

save_multiformat_ggplot(
  p_52e, file.path(fig_52e_dir, "52e_metagene_profile"),
  width = 10, height = 6
)

cat(sprintf("  mC profile: %d bins with data (mean n=%.0f genes/bin)\n",
            nrow(mc_profile), mean(mc_profile$n)))
cat(sprintf("  hmC profile: %d bins with data (mean n=%.0f genes/bin)\n",
            nrow(hmc_profile), mean(hmc_profile$n)))

# =============================================================================
# SUMMARY TABLE
# =============================================================================

cat("\n--- Writing summary table ---\n")

summary_table <- feat_data %>%
  dplyr::select(gene, feature_type, Name, Chromosome, Start, End,
                ctrl_mc, mut_mc, delta_mc,
                ctrl_hmc, mut_hmc, delta_hmc) %>%
  dplyr::mutate(
    delta_mc_pct = round(delta_mc * 100, 2),
    delta_hmc_pct = round(delta_hmc * 100, 2)
  ) %>%
  dplyr::arrange(gene, Start)

table_path <- file.path(TABLES_DIR, "section52_feature_methylation.tsv")
write.table(summary_table, table_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Wrote: %s (%d rows)\n", table_path, nrow(summary_table)))

# Dunn's post-hoc results table
mc_posthoc$modification <- "mC"
hmc_posthoc$modification <- "hmC"
posthoc_combined <- rbind(mc_posthoc, hmc_posthoc)
posthoc_path <- file.path(TABLES_DIR, "section52_dunn_posthoc.tsv")
write.table(posthoc_combined, posthoc_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Wrote: %s (%d rows)\n", posthoc_path, nrow(posthoc_combined)))

# =============================================================================
# DONE
# =============================================================================

cat("\n========================================================================\n")
cat("SECTION 52 COMPLETE\n")
cat(sprintf("  Figures: %s\n", SEC52_DIR))
cat(sprintf("  Table:   %s\n", table_path))
cat("========================================================================\n\n")
