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
# FIGURE 52f: Dunn's Post-Hoc Significance Heatmaps
# =============================================================================

cat("\n--- Figure 52f: Dunn's post-hoc significance heatmaps ---\n")

fig_52f_dir <- file.path(SEC52_DIR, "52f_dunn_posthoc")
dir.create(fig_52f_dir, recursive = TRUE, showWarnings = FALSE)

build_pairwise_matrix <- function(posthoc_df, feature_levels) {
  n <- length(feature_levels)
  mat <- matrix(NA, nrow = n, ncol = n,
                dimnames = list(feature_levels, feature_levels))
  for (i in seq_len(nrow(posthoc_df))) {
    parts <- trimws(strsplit(posthoc_df$comparison[i], " - ")[[1]])
    if (length(parts) == 2 && all(parts %in% feature_levels)) {
      mat[parts[1], parts[2]] <- posthoc_df$p_adj[i]
      mat[parts[2], parts[1]] <- posthoc_df$p_adj[i]
    }
  }
  diag(mat) <- 1
  mat
}

mc_mat <- build_pairwise_matrix(mc_posthoc, levels(feat_data$feature_type))
hmc_mat <- build_pairwise_matrix(hmc_posthoc, levels(feat_data$feature_type))

tile_from_matrix <- function(mat, title_text, subtitle_text) {
  df <- expand.grid(
    Feature1 = rownames(mat), Feature2 = colnames(mat),
    stringsAsFactors = FALSE
  )
  df$q_value <- mapply(function(r, c) mat[r, c], df$Feature1, df$Feature2)
  # Lower triangle only (including diagonal)
  idx <- match(df$Feature1, rownames(mat))
  idy <- match(df$Feature2, colnames(mat))
  df <- df[idx >= idy, ]
  df$neg_log10_q <- -log10(pmax(df$q_value, 1e-10))
  df$label <- ifelse(
    is.na(df$q_value), "",
    ifelse(df$Feature1 == df$Feature2, "",
           ifelse(df$q_value < 0.001, sprintf("%.1e", df$q_value),
                  sprintf("%.3f", df$q_value))))
  df$sig <- ifelse(is.na(df$q_value) | df$Feature1 == df$Feature2, "",
                   ifelse(df$q_value < 0.001, "***",
                          ifelse(df$q_value < 0.01, "**",
                                 ifelse(df$q_value < 0.05, "*", ""))))
  df$display <- paste0(df$label, ifelse(df$sig == "", "", paste0("\n", df$sig)))
  df$Feature1 <- factor(df$Feature1, levels = rev(rownames(mat)))
  df$Feature2 <- factor(df$Feature2, levels = colnames(mat))

  ggplot(df[df$Feature1 != df$Feature2, ],
         aes(x = Feature2, y = Feature1, fill = neg_log10_q)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = display), size = 3, lineheight = 0.85) +
    scale_fill_gradient2(
      low = "grey90", mid = "#FDB863", high = "#B2182B",
      midpoint = -log10(0.05), limits = c(0, NA),
      name = expression(-log[10] * "(q)")) +
    labs(title = title_text, subtitle = subtitle_text,
         x = NULL, y = NULL) +
    theme_biomodal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid = element_blank())
}

p_52f_mc <- tile_from_matrix(
  mc_mat,
  expression("Dunn's Post-Hoc: " * Delta * "5mC"),
  "BH-adjusted q-values; * < 0.05, ** < 0.01, *** < 0.001"
)

p_52f_hmc <- tile_from_matrix(
  hmc_mat,
  expression("Dunn's Post-Hoc: " * Delta * "5hmC"),
  "BH-adjusted q-values; * < 0.05, ** < 0.01, *** < 0.001"
)

p_52f_combined <- p_52f_mc + p_52f_hmc +
  patchwork::plot_layout(ncol = 2, guides = "collect")

save_multiformat_ggplot(
  p_52f_combined, file.path(fig_52f_dir, "52f_dunn_posthoc"),
  width = 16, height = 7
)

# Also save individual panels
save_multiformat_ggplot(
  p_52f_mc, file.path(fig_52f_dir, "52f_dunn_mc"),
  width = 8, height = 7
)
save_multiformat_ggplot(
  p_52f_hmc, file.path(fig_52f_dir, "52f_dunn_hmc"),
  width = 8, height = 7
)

# =============================================================================
# FIGURE 52g: Chromatin Mark Overlay on Sub-Gene Features
# =============================================================================

cat("\n--- Figure 52g: Chromatin mark overlay ---\n")

fig_52g_dir <- file.path(SEC52_DIR, "52g_chromatin_overlay")
dir.create(fig_52g_dir, recursive = TRUE, showWarnings = FALSE)

load_peaks_gr <- function(bed_path) {
  stopifnot(file.exists(bed_path))
  df <- read.table(bed_path, header = FALSE, sep = "\t",
                   stringsAsFactors = FALSE)[, 1:3]
  colnames(df) <- c("chr", "start", "end")
  GRanges(seqnames = df$chr, ranges = IRanges(df$start, df$end))
}

k27ac_peaks  <- load_peaks_gr(CHIP_PEAK_FILES$h3k27ac)
k4me1_peaks  <- load_peaks_gr(CHIP_PEAK_FILES$h3k4me1)
k27me3_peaks <- load_peaks_gr(CHIP_PEAK_FILES$h3k27me3)
cat(sprintf("  Loaded: H3K27ac=%d, H3K4me1=%d, H3K27me3=%d peaks\n",
            length(k27ac_peaks), length(k4me1_peaks), length(k27me3_peaks)))

feat_gr <- GRanges(
  seqnames = paste0("chr", feat_data$Chromosome),
  ranges = IRanges(feat_data$Start, feat_data$End)
)

feat_data$h3k27ac  <- overlapsAny(feat_gr, k27ac_peaks)
feat_data$h3k4me1  <- overlapsAny(feat_gr, k4me1_peaks)
feat_data$h3k27me3 <- overlapsAny(feat_gr, k27me3_peaks)

cat("\n  ChIP-seq peak overlap by feature type (%):\n")
overlap_summary <- feat_data %>%
  dplyr::group_by(feature_type) %>%
  dplyr::summarise(
    n = dplyr::n(),
    H3K27ac = round(100 * mean(h3k27ac), 1),
    H3K4me1 = round(100 * mean(h3k4me1), 1),
    H3K27me3 = round(100 * mean(h3k27me3), 1),
    .groups = "drop"
  )
print(as.data.frame(overlap_summary))

# --- 52g panel A: Overlap fraction per feature type (stacked bar) ---
overlap_long <- feat_data %>%
  dplyr::select(feature_type, h3k27ac, h3k4me1, h3k27me3) %>%
  tidyr::pivot_longer(cols = c(h3k27ac, h3k4me1, h3k27me3),
                      names_to = "mark", values_to = "overlaps") %>%
  dplyr::group_by(feature_type, mark) %>%
  dplyr::summarise(pct = 100 * mean(overlaps), .groups = "drop") %>%
  dplyr::mutate(mark = dplyr::recode(mark,
    h3k27ac = "H3K27ac", h3k4me1 = "H3K4me1", h3k27me3 = "H3K27me3"))

MARK_COLORS <- c("H3K27ac" = "#FF7F00", "H3K4me1" = "#33A02C", "H3K27me3" = "#6A3D9A")

p_52g_bar <- ggplot(overlap_long,
                    aes(x = feature_type, y = pct, fill = mark)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = MARK_COLORS, name = "ChIP Mark") +
  labs(title = "Chromatin Mark Overlap by Gene Feature",
       subtitle = "% of sub-gene intervals overlapping ChIP-seq peaks",
       x = "Gene Feature", y = "% Overlapping") +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_multiformat_ggplot(
  p_52g_bar, file.path(fig_52g_dir, "52g_overlap_by_feature"),
  width = 10, height = 6
)

# --- 52g panel B: Introns split by H3K27ac — delta mC and hmC ---
intron_data <- feat_data[feat_data$feature_type == "Intron", ]
intron_data$enhancer <- ifelse(intron_data$h3k27ac, "H3K27ac+", "Unmarked")

n_enh <- sum(intron_data$enhancer == "H3K27ac+")
n_unm <- sum(intron_data$enhancer == "Unmarked")
cat(sprintf("\n  Introns: %d H3K27ac+ (enhancer), %d unmarked\n", n_enh, n_unm))

intron_long <- rbind(
  data.frame(enhancer = intron_data$enhancer,
             delta = intron_data$delta_mc * 100,
             modification = "5mC", stringsAsFactors = FALSE),
  data.frame(enhancer = intron_data$enhancer,
             delta = intron_data$delta_hmc * 100,
             modification = "5hmC", stringsAsFactors = FALSE)
)
intron_long$modification <- factor(intron_long$modification,
                                   levels = c("5mC", "5hmC"))

ENH_COLORS <- c("H3K27ac+" = "#FF7F00", "Unmarked" = "grey70")

p_52g_intron <- ggplot(intron_long,
                       aes(x = enhancer, y = delta, fill = enhancer)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_boxplot(outlier.size = 0.6, alpha = 0.85) +
  facet_wrap(~ modification, scales = "free_y") +
  scale_fill_manual(values = ENH_COLORS, guide = "none") +
  labs(title = "Intronic Methylation Change by Enhancer Status",
       subtitle = sprintf("H3K27ac+ (n=%d) vs unmarked (n=%d) introns",
                          n_enh, n_unm),
       x = NULL, y = expression(Delta * " Methylation (%)")) +
  theme_biomodal()

save_multiformat_ggplot(
  p_52g_intron, file.path(fig_52g_dir, "52g_intron_enhancer_split"),
  width = 9, height = 6
)

# Wilcoxon rank-sum tests: H3K27ac+ vs unmarked introns
for (mod_name in c("delta_mc", "delta_hmc")) {
  marked <- intron_data[[mod_name]][intron_data$h3k27ac]
  unmarked <- intron_data[[mod_name]][!intron_data$h3k27ac]
  wt <- wilcox.test(marked, unmarked)
  cat(sprintf("  Wilcoxon (intron %s): H3K27ac+ median=%.2f%% vs unmarked median=%.2f%%, p=%s\n",
              mod_name,
              median(marked, na.rm = TRUE) * 100,
              median(unmarked, na.rm = TRUE) * 100,
              format.pval(wt$p.value, digits = 3)))
}

# --- 52g panel C: All features split by H3K27ac ---
feat_data$enhancer <- ifelse(feat_data$h3k27ac, "H3K27ac+", "Unmarked")

all_long <- rbind(
  data.frame(feature_type = feat_data$feature_type,
             enhancer = feat_data$enhancer,
             delta = feat_data$delta_hmc * 100,
             modification = "5hmC", stringsAsFactors = FALSE),
  data.frame(feature_type = feat_data$feature_type,
             enhancer = feat_data$enhancer,
             delta = feat_data$delta_mc * 100,
             modification = "5mC", stringsAsFactors = FALSE)
)

p_52g_all <- ggplot(all_long[all_long$modification == "5hmC", ],
                    aes(x = feature_type, y = delta, fill = enhancer)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_boxplot(outlier.size = 0.5, alpha = 0.85,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = ENH_COLORS, name = "H3K27ac") +
  labs(title = expression("All Features: " * Delta * "5hmC by Enhancer Status"),
       subtitle = "Does H3K27ac overlap explain the feature-type differences?",
       x = "Gene Feature", y = expression(Delta * "5hmC (%)")) +
  theme_biomodal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_multiformat_ggplot(
  p_52g_all, file.path(fig_52g_dir, "52g_all_features_enhancer_split"),
  width = 11, height = 6
)

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
