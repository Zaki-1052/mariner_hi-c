# abc/scripts/step14_loop_type_at_enhancer_subsets.R
# Loop type decomposition at enhancer phenotypic subsets.
#
# Breaks down 7-category loop types (Active_Promoter, Repressed_Promoter,
# Bivalent_Promoter, Polycomb, Active_Enhancer, Poised_Enhancer, Other)
# present at each enhancer class (Activity_Lost, K119ub_Only, Activity_Gain,
# Stable). Focuses on the "partner anchor" — what the enhancer connects to.
#
# Key question: if enhancer-TSS (promoter) changes are minimal at K119ub_Only
# sites, what other loop types are present and contributing to the logFC?
#
# Inputs (all relative to abc/ working directory):
#   enhancer_subsets/enhancer_classes_*.tsv    -- 4 enhancer class files
#   characterized_loops.tsv                    -- differential loops with
#                                                 7-cat type annotations
#
# Outputs:
#   results/loop_type_at_subsets/              -- TSVs + plots (PDF + SVG + JPG)
#
# Usage:
#   cd abc && Rscript scripts/step14_loop_type_at_enhancer_subsets.R

# =============================================================================
# CONFIGURATION
# =============================================================================

ENHANCER_FILES <- c(
  Activity_Lost = "enhancer_subsets/enhancer_classes_activity_lost.tsv",
  K119ub_Only   = "enhancer_subsets/enhancer_classes_k119ub_only.tsv",
  Activity_Gain = "enhancer_subsets/enhancer_classes_activity_gain.tsv",
  Stable        = "enhancer_subsets/enhancer_classes_stable.tsv"
)

CHAR_LOOPS_FILE <- "characterized_loops.tsv"

OUTPUT_DIR <- "results/loop_type_at_subsets"

CLASS_COLORS <- c(
  Activity_Lost = "#2166AC",
  K119ub_Only   = "#B2182B",
  Activity_Gain = "#F4A582",
  Stable        = "#D1E5F0"
)
CLASS_ORDER <- c("Activity_Lost", "K119ub_Only", "Activity_Gain", "Stable")

PARTNER_COLORS <- c(
  Active_Promoter    = "#E41A1C",
  Repressed_Promoter = "#377EB8",
  Bivalent_Promoter  = "#FF7F00",
  Polycomb           = "#984EA3",
  Active_Enhancer    = "#4DAF4A",
  Poised_Enhancer    = "#A6D854",
  Other              = "#999999"
)
PARTNER_ORDER <- c("Active_Promoter", "Repressed_Promoter", "Bivalent_Promoter",
                   "Polycomb", "Active_Enhancer", "Poised_Enhancer", "Other")

FUNCTIONAL_GROUPS <- c(
  Active_Promoter    = "Promoter",
  Repressed_Promoter = "Promoter",
  Bivalent_Promoter  = "Promoter",
  Polycomb           = "Polycomb",
  Active_Enhancer    = "Enhancer",
  Poised_Enhancer    = "Enhancer",
  Other              = "Structural"
)
FUNC_COLORS <- c(
  Promoter   = "#E41A1C",
  Polycomb   = "#984EA3",
  Enhancer   = "#4DAF4A",
  Structural = "#999999"
)
FUNC_ORDER <- c("Promoter", "Polycomb", "Enhancer", "Structural")

cat("================================================================================\n")
cat("STEP 14: LOOP TYPE DECOMPOSITION AT ENHANCER SUBSETS\n")
cat("================================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# =============================================================================
# LOAD PACKAGES
# =============================================================================

cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(ggplot2)
  library(scales)
})
cat("Packages loaded.\n\n")

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

cat("Validating input files...\n")
for (nm in names(ENHANCER_FILES)) {
  stopifnot(file.exists(ENHANCER_FILES[[nm]]))
  cat(sprintf("  [OK] %s: %s\n", nm, ENHANCER_FILES[[nm]]))
}
stopifnot(file.exists(CHAR_LOOPS_FILE))
cat(sprintf("  [OK] Characterized loops: %s\n", CHAR_LOOPS_FILE))

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("\nOutput directory: %s\n\n", OUTPUT_DIR))

# =============================================================================
# HELPERS
# =============================================================================

theme_pub <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

UTIL_PATH <- "../scripts/utils/multi_format_output.R"
stopifnot(file.exists(UTIL_PATH))
source(UTIL_PATH)

save_plot <- function(p, name, w = 7, h = 6) {
  save_multiformat_ggplot(p, base_path = file.path(OUTPUT_DIR, name),
                          width = w, height = h, dpi = 300,
                          verbose = TRUE, use_subfolders = TRUE)
}

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")

# --- Enhancer classes ---
enh_list <- lapply(names(ENHANCER_FILES), function(cls) {
  df <- read.table(ENHANCER_FILES[[cls]], sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE)
  df$enhancer_class <- cls
  df
})
enh_all <- do.call(rbind, enh_list)
enh_all$enhancer_class <- factor(enh_all$enhancer_class, levels = CLASS_ORDER)
cat(sprintf("  Enhancers: %d total across %d classes\n", nrow(enh_all),
            length(CLASS_ORDER)))
for (cls in CLASS_ORDER) {
  cat(sprintf("    %s: %d\n", cls, sum(enh_all$enhancer_class == cls)))
}

# --- Characterized loops (differential, with 7-category type annotations) ---
char_loops <- read.table(CHAR_LOOPS_FILE, sep = "\t", header = TRUE,
                         stringsAsFactors = FALSE)
cat(sprintf("  Characterized loops: %d\n", nrow(char_loops)))

# Validate required columns
required_cols <- c("chr1", "start1", "end1", "chr2", "start2", "end2",
                   "logFC", "FDR", "direction", "anchor1_type", "anchor2_type",
                   "loop_type")
missing <- setdiff(required_cols, colnames(char_loops))
if (length(missing) > 0) {
  stop(sprintf("Missing columns in characterized_loops.tsv: %s",
               paste(missing, collapse = ", ")))
}

# Validate anchor types are from the 7-category system
valid_types <- c("Active_Promoter", "Repressed_Promoter", "Bivalent_Promoter",
                 "Polycomb", "Active_Enhancer", "Poised_Enhancer", "Other")
actual_types <- unique(c(char_loops$anchor1_type, char_loops$anchor2_type))
unexpected <- setdiff(actual_types, valid_types)
if (length(unexpected) > 0) {
  stop(sprintf("Unexpected anchor types: %s", paste(unexpected, collapse = ", ")))
}

cat(sprintf("  Anchor types present: %s\n", paste(sort(actual_types), collapse = ", ")))
cat(sprintf("  Unique loop types: %d\n", length(unique(char_loops$loop_type))))
cat("\n")


# #############################################################################
# PART A: ENHANCER-LOOP OVERLAP WITH PARTNER ANCHOR CLASSIFICATION
# #############################################################################

cat("=== PART A: Enhancer-Loop Overlap & Partner Classification ===\n\n")

# Build GRanges
enh_gr <- GRanges(
  seqnames = enh_all$chr,
  ranges = IRanges(start = enh_all$start, end = enh_all$end),
  enhancer_class = enh_all$enhancer_class,
  enh_idx = seq_len(nrow(enh_all))
)

loop_anchor1_gr <- GRanges(
  seqnames = char_loops$chr1,
  ranges = IRanges(start = char_loops$start1, end = char_loops$end1),
  loop_idx = seq_len(nrow(char_loops))
)

loop_anchor2_gr <- GRanges(
  seqnames = char_loops$chr2,
  ranges = IRanges(start = char_loops$start2, end = char_loops$end2),
  loop_idx = seq_len(nrow(char_loops))
)

cat("GRanges built.\n")

# Find overlaps: enhancer <-> loop anchors
hits_a1 <- findOverlaps(enh_gr, loop_anchor1_gr, type = "any")
hits_a2 <- findOverlaps(enh_gr, loop_anchor2_gr, type = "any")

cat(sprintf("  Enhancer <-> anchor1 overlaps: %d\n", length(hits_a1)))
cat(sprintf("  Enhancer <-> anchor2 overlaps: %d\n", length(hits_a2)))

# Build enhancer-loop pair table
# When enhancer overlaps anchor1 -> partner is anchor2 (and vice versa)
pairs_a1 <- data.frame(
  enh_idx       = queryHits(hits_a1),
  loop_idx      = subjectHits(hits_a1),
  enh_on_anchor = 1L,
  partner_type  = char_loops$anchor2_type[subjectHits(hits_a1)],
  own_type      = char_loops$anchor1_type[subjectHits(hits_a1)],
  stringsAsFactors = FALSE
)

pairs_a2 <- data.frame(
  enh_idx       = queryHits(hits_a2),
  loop_idx      = subjectHits(hits_a2),
  enh_on_anchor = 2L,
  partner_type  = char_loops$anchor1_type[subjectHits(hits_a2)],
  own_type      = char_loops$anchor2_type[subjectHits(hits_a2)],
  stringsAsFactors = FALSE
)

enh_loop_pairs <- rbind(pairs_a1, pairs_a2)

# Deduplicate: if enhancer overlaps BOTH anchors of same loop, keep first
enh_loop_pairs <- enh_loop_pairs[
  !duplicated(enh_loop_pairs[, c("enh_idx", "loop_idx")]), ]

cat(sprintf("  Unique enhancer-loop pairs (after dedup): %d\n",
            nrow(enh_loop_pairs)))

# Attach metadata
enh_loop_pairs$enhancer_class <- enh_all$enhancer_class[enh_loop_pairs$enh_idx]
enh_loop_pairs$loop_type      <- char_loops$loop_type[enh_loop_pairs$loop_idx]
enh_loop_pairs$logFC           <- char_loops$logFC[enh_loop_pairs$loop_idx]
enh_loop_pairs$FDR             <- char_loops$FDR[enh_loop_pairs$loop_idx]
enh_loop_pairs$direction       <- char_loops$direction[enh_loop_pairs$loop_idx]

# Factor levels
enh_loop_pairs$enhancer_class <- factor(enh_loop_pairs$enhancer_class,
                                         levels = CLASS_ORDER)
enh_loop_pairs$partner_type <- factor(enh_loop_pairs$partner_type,
                                       levels = PARTNER_ORDER)

# Functional grouping
enh_loop_pairs$partner_group <- factor(
  FUNCTIONAL_GROUPS[as.character(enh_loop_pairs$partner_type)],
  levels = FUNC_ORDER
)

cat(sprintf("  Enhancers with >=1 loop: %d / %d (%.1f%%)\n",
            length(unique(enh_loop_pairs$enh_idx)), nrow(enh_all),
            100 * length(unique(enh_loop_pairs$enh_idx)) / nrow(enh_all)))

# Per-class summary
cat("\n  Per-class enhancer-loop pair counts:\n")
for (cls in CLASS_ORDER) {
  sub <- enh_loop_pairs[enh_loop_pairs$enhancer_class == cls, ]
  n_enh <- length(unique(sub$enh_idx))
  n_total_enh <- sum(enh_all$enhancer_class == cls)
  cat(sprintf("    %s: %d pairs, %d unique enhancers (%.1f%% of %d)\n",
              cls, nrow(sub), n_enh, 100 * n_enh / n_total_enh, n_total_enh))
}

# Sanity check: overall logFC by class (should match step11 Part C)
cat("\n  Sanity check - overall median logFC per class:\n")
for (cls in CLASS_ORDER) {
  sub <- enh_loop_pairs[enh_loop_pairs$enhancer_class == cls, ]
  cat(sprintf("    %s: median logFC = %.4f (n=%d)\n",
              cls, median(sub$logFC), nrow(sub)))
}

# Save detailed pairs table
detailed_pairs <- enh_loop_pairs
detailed_pairs$enh_chr   <- enh_all$chr[detailed_pairs$enh_idx]
detailed_pairs$enh_start <- enh_all$start[detailed_pairs$enh_idx]
detailed_pairs$enh_end   <- enh_all$end[detailed_pairs$enh_idx]
detailed_pairs$loop_chr1   <- char_loops$chr1[detailed_pairs$loop_idx]
detailed_pairs$loop_start1 <- char_loops$start1[detailed_pairs$loop_idx]
detailed_pairs$loop_end1   <- char_loops$end1[detailed_pairs$loop_idx]
detailed_pairs$loop_chr2   <- char_loops$chr2[detailed_pairs$loop_idx]
detailed_pairs$loop_start2 <- char_loops$start2[detailed_pairs$loop_idx]
detailed_pairs$loop_end2   <- char_loops$end2[detailed_pairs$loop_idx]

write.table(detailed_pairs,
            file.path(OUTPUT_DIR, "detailed_enhancer_loop_pairs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  Saved: %s\n",
            file.path(OUTPUT_DIR, "detailed_enhancer_loop_pairs.tsv")))


# #############################################################################
# PART B: LOOP TYPE COMPOSITION BY ENHANCER CLASS
# #############################################################################

cat("\n=== PART B: Loop Type Composition ===\n\n")

# --- Composition table (7-cat partner types) ---
comp_table <- as.data.frame(table(
  enhancer_class = enh_loop_pairs$enhancer_class,
  partner_type = enh_loop_pairs$partner_type
))
comp_table$proportion <- ave(comp_table$Freq, comp_table$enhancer_class,
                              FUN = function(x) x / sum(x))

# Print summary
cat("  Partner type composition by enhancer class:\n")
for (cls in CLASS_ORDER) {
  sub <- comp_table[comp_table$enhancer_class == cls, ]
  sub <- sub[order(-sub$proportion), ]
  n_total <- sum(sub$Freq)
  cat(sprintf("\n  %s (n=%d loops):\n", cls, n_total))
  for (i in seq_len(nrow(sub))) {
    if (sub$Freq[i] > 0) {
      cat(sprintf("    %-22s %4d (%5.1f%%)\n",
                  as.character(sub$partner_type[i]),
                  sub$Freq[i], 100 * sub$proportion[i]))
    }
  }
}

write.table(comp_table, file.path(OUTPUT_DIR, "partner_type_composition.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Chi-squared test ---
contingency <- xtabs(Freq ~ enhancer_class + partner_type, data = comp_table)
chi_test <- chisq.test(contingency, simulate.p.value = TRUE, B = 10000)
cat(sprintf("\n  Chi-squared (partner type ~ enhancer class): X2=%.1f, %s\n",
            chi_test$statistic, fmt_p(chi_test$p.value)))

# --- Plot 01: Stacked bar — partner type by enhancer class (proportions) ---
cat("\nGenerating Part B plots...\n")

p01 <- ggplot(comp_table[comp_table$Freq > 0, ],
              aes(x = enhancer_class, y = Freq, fill = partner_type)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = PARTNER_COLORS) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    x = "Enhancer class", y = "Proportion of loops",
    fill = "Partner anchor type",
    title = "Loop partner composition by enhancer class",
    subtitle = sprintf("7-category extended classification | Chi-sq %s",
                       fmt_p(chi_test$p.value))
  ) +
  theme_pub +
  guides(fill = guide_legend(nrow = 2))
save_plot(p01, "01_partner_type_stacked_bar", w = 8, h = 6)

# --- Plot 02: Simplified functional grouping stacked bar ---
func_comp <- as.data.frame(table(
  enhancer_class = enh_loop_pairs$enhancer_class,
  partner_group = enh_loop_pairs$partner_group
))
func_comp$proportion <- ave(func_comp$Freq, func_comp$enhancer_class,
                             FUN = function(x) x / sum(x))

n_labels <- aggregate(Freq ~ enhancer_class, data = func_comp, FUN = sum)
colnames(n_labels) <- c("enhancer_class", "n_total")

p02 <- ggplot(func_comp[func_comp$Freq > 0, ],
              aes(x = enhancer_class, y = Freq, fill = partner_group)) +
  geom_col(position = "fill") +
  geom_text(data = n_labels, aes(x = enhancer_class, y = 1.03,
            label = paste0("n=", format(n_total, big.mark = ","))),
            inherit.aes = FALSE, size = 3.2) +
  scale_fill_manual(values = FUNC_COLORS) +
  scale_y_continuous(labels = percent_format(),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = "Enhancer class", y = "Proportion of loops",
    fill = "Partner category",
    title = "Functional loop partner categories by enhancer class",
    subtitle = paste("Promoter (Active/Repressed/Bivalent) | Polycomb |",
                     "Enhancer (Active/Poised) | Structural (Other)")
  ) +
  theme_pub
save_plot(p02, "02_functional_group_stacked_bar", w = 8, h = 6)

# --- Plot 03: Heatmap — proportions (partner_type x enhancer_class) ---
heat_data <- comp_table[, c("enhancer_class", "partner_type",
                             "proportion", "Freq")]
heat_data$label <- ifelse(
  heat_data$Freq > 0,
  sprintf("%d\n(%.0f%%)", heat_data$Freq, 100 * heat_data$proportion),
  ""
)

p03 <- ggplot(heat_data, aes(x = enhancer_class, y = partner_type,
                              fill = proportion)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, NA),
                      labels = percent_format()) +
  scale_y_discrete(limits = rev(PARTNER_ORDER)) +
  labs(
    x = "Enhancer class", y = "Partner anchor type",
    fill = "Proportion",
    title = "Partner type proportions by enhancer class"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p03, "03_partner_type_heatmap", w = 8, h = 7)


# #############################################################################
# PART C: logFC BY LOOP TYPE x ENHANCER CLASS
# #############################################################################

cat("\n=== PART C: logFC by Loop Type x Enhancer Class ===\n\n")

# --- Build summary statistics table ---
logfc_summary <- do.call(rbind, lapply(CLASS_ORDER, function(cls) {
  do.call(rbind, lapply(PARTNER_ORDER, function(pt) {
    sub <- enh_loop_pairs[enh_loop_pairs$enhancer_class == cls &
                            enh_loop_pairs$partner_type == pt, ]
    n <- nrow(sub)
    if (n >= 3) {
      wt <- tryCatch(wilcox.test(sub$logFC, mu = 0),
                     error = function(e) list(p.value = NA))
      data.frame(
        enhancer_class = cls,
        partner_type = pt,
        n = n,
        median_logFC = median(sub$logFC),
        mean_logFC = mean(sub$logFC),
        sd_logFC = sd(sub$logFC),
        wilcox_p = wt$p.value,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        enhancer_class = cls,
        partner_type = pt,
        n = n,
        median_logFC = ifelse(n > 0, median(sub$logFC), NA_real_),
        mean_logFC = ifelse(n > 0, mean(sub$logFC), NA_real_),
        sd_logFC = NA_real_,
        wilcox_p = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }))
}))

# Print formatted table
cat("  logFC by partner type x enhancer class (n >= 5):\n\n")
cat(sprintf("  %-15s %-22s %5s %10s %12s\n",
            "Class", "Partner Type", "n", "med logFC", "Wilcoxon p"))
cat(paste(rep("-", 72), collapse = ""), "\n")
for (i in seq_len(nrow(logfc_summary))) {
  row <- logfc_summary[i, ]
  if (row$n >= 5) {
    cat(sprintf("  %-15s %-22s %5d %10.4f %12s\n",
                row$enhancer_class, row$partner_type, row$n,
                row$median_logFC, fmt_p(row$wilcox_p)))
  }
}

write.table(logfc_summary,
            file.path(OUTPUT_DIR, "loop_logfc_by_type_and_class.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Plot 04: Faceted boxplot — logFC by partner type, faceted by class ---
cat("\nGenerating Part C plots...\n")

# Filter to partner types with sufficient data
pt_counts <- table(enh_loop_pairs$partner_type)
viable_pts <- names(pt_counts[pt_counts >= 5])
plot_data <- enh_loop_pairs[enh_loop_pairs$partner_type %in% viable_pts, ]
plot_data$partner_type <- factor(plot_data$partner_type, levels = PARTNER_ORDER)

p04 <- ggplot(plot_data, aes(x = partner_type, y = logFC, fill = partner_type)) +
  geom_boxplot(outlier.size = 0.3, notch = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ enhancer_class, scales = "free_x") +
  scale_fill_manual(values = PARTNER_COLORS) +
  labs(
    x = "Partner anchor type", y = "Loop logFC (KO/WT)",
    title = "Loop strength change by partner type and enhancer class",
    subtitle = "Dashed line = no change"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )
save_plot(p04, "04_logfc_by_partner_faceted", w = 12, h = 8)

# --- Plot 05: Heatmap — median logFC (partner_type x enhancer_class) ---
heat_logfc <- logfc_summary[logfc_summary$n >= 5, ]
heat_logfc$enhancer_class <- factor(heat_logfc$enhancer_class, levels = CLASS_ORDER)
heat_logfc$partner_type <- factor(heat_logfc$partner_type, levels = PARTNER_ORDER)

# Significance markers
heat_logfc$sig <- ""
heat_logfc$sig[!is.na(heat_logfc$wilcox_p) &
                 heat_logfc$wilcox_p < 0.05] <- "*"
heat_logfc$sig[!is.na(heat_logfc$wilcox_p) &
                 heat_logfc$wilcox_p < 0.001] <- "**"
heat_logfc$sig[!is.na(heat_logfc$wilcox_p) &
                 heat_logfc$wilcox_p < 1e-10] <- "***"

heat_logfc$label <- sprintf("%.3f\n(n=%d)%s",
                             heat_logfc$median_logFC,
                             heat_logfc$n,
                             heat_logfc$sig)

p05 <- ggplot(heat_logfc, aes(x = enhancer_class, y = partner_type,
                               fill = median_logFC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 2.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-0.2, 0.2), oob = squish,
                       name = "Median\nlogFC") +
  scale_y_discrete(limits = rev(PARTNER_ORDER)) +
  labs(
    x = "Enhancer class", y = "Partner anchor type",
    title = "Median loop logFC by partner type x enhancer class",
    subtitle = "Wilcoxon vs 0: * p<0.05, ** p<0.001, *** p<1e-10"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p05, "05_logfc_heatmap", w = 8, h = 7)

# --- Plot 06: Contribution analysis ---
# For each class: contribution_i = (n_i / N_total) * median_logFC_i
contrib <- logfc_summary[logfc_summary$n >= 3, ]
contrib <- do.call(rbind, lapply(CLASS_ORDER, function(cls) {
  sub <- contrib[contrib$enhancer_class == cls, ]
  n_total <- sum(sub$n)
  sub$weight <- sub$n / n_total
  sub$contribution <- sub$weight * sub$median_logFC
  sub
}))

cat("\n  Contribution analysis (weighted logFC by partner type):\n")
for (cls in CLASS_ORDER) {
  sub <- contrib[contrib$enhancer_class == cls, ]
  sub <- sub[order(-abs(sub$contribution)), ]
  total_contrib <- sum(sub$contribution, na.rm = TRUE)
  cat(sprintf("\n  %s (total weighted logFC = %.4f):\n", cls, total_contrib))
  for (i in seq_len(min(nrow(sub), 7))) {
    row <- sub[i, ]
    if (!is.na(row$contribution)) {
      cat(sprintf("    %-22s weight=%.3f x median=%.4f = %.5f\n",
                  row$partner_type, row$weight, row$median_logFC,
                  row$contribution))
    }
  }
}

write.table(contrib, file.path(OUTPUT_DIR, "contribution_analysis.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

contrib_plot <- contrib[!is.na(contrib$contribution), ]
contrib_plot$enhancer_class <- factor(contrib_plot$enhancer_class,
                                       levels = CLASS_ORDER)
contrib_plot$partner_type <- factor(contrib_plot$partner_type,
                                     levels = PARTNER_ORDER)

p06 <- ggplot(contrib_plot, aes(x = partner_type, y = contribution,
                                 fill = partner_type)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ enhancer_class) +
  scale_fill_manual(values = PARTNER_COLORS) +
  labs(
    x = "Partner anchor type",
    y = "Weighted logFC contribution\n(weight x median logFC)",
    title = "Loop type contribution to overall logFC by enhancer class",
    subtitle = "Weight = fraction of loops with that partner type"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )
save_plot(p06, "06_contribution_analysis", w = 12, h = 8)


# #############################################################################
# PART D: KEY COMPARISONS AND STATISTICS
# #############################################################################

cat("\n=== PART D: Key Comparisons ===\n\n")

# --- K119ub_Only vs Stable: Wilcoxon by partner type ---
cat("  K119ub_Only vs Stable logFC (by partner type):\n")
ks_comparisons <- list()
for (pt in PARTNER_ORDER) {
  k_vals <- enh_loop_pairs$logFC[enh_loop_pairs$enhancer_class == "K119ub_Only" &
                                   enh_loop_pairs$partner_type == pt]
  s_vals <- enh_loop_pairs$logFC[enh_loop_pairs$enhancer_class == "Stable" &
                                   enh_loop_pairs$partner_type == pt]
  if (length(k_vals) >= 5 && length(s_vals) >= 5) {
    wt <- wilcox.test(k_vals, s_vals)
    cat(sprintf("    %-22s K119ub(n=%d, med=%.4f) vs Stable(n=%d, med=%.4f): %s\n",
                pt, length(k_vals), median(k_vals),
                length(s_vals), median(s_vals), fmt_p(wt$p.value)))
    ks_comparisons[[pt]] <- list(
      partner_type = pt,
      k_n = length(k_vals), k_median = median(k_vals),
      s_n = length(s_vals), s_median = median(s_vals),
      p = wt$p.value
    )
  } else {
    cat(sprintf("    %-22s insufficient data (K=%d, S=%d)\n",
                pt, length(k_vals), length(s_vals)))
  }
}

# --- K119ub_Only vs Activity_Lost: Fisher's exact for partner type enrichment ---
cat("\n  K119ub_Only vs Activity_Lost: partner type enrichment (Fisher's exact):\n")
k_all <- enh_loop_pairs[enh_loop_pairs$enhancer_class == "K119ub_Only", ]
al_all <- enh_loop_pairs[enh_loop_pairs$enhancer_class == "Activity_Lost", ]

fisher_results <- list()
for (pt in PARTNER_ORDER) {
  k_yes <- sum(k_all$partner_type == pt)
  k_no  <- nrow(k_all) - k_yes
  al_yes <- sum(al_all$partner_type == pt)
  al_no  <- nrow(al_all) - al_yes

  mat <- matrix(c(k_yes, k_no, al_yes, al_no), nrow = 2, byrow = TRUE)
  ft <- fisher.test(mat)
  cat(sprintf("    %-22s K119ub=%d/%d (%.1f%%) vs ActLost=%d/%d (%.1f%%) OR=%.2f %s\n",
              pt, k_yes, nrow(k_all), 100 * k_yes / nrow(k_all),
              al_yes, nrow(al_all), 100 * al_yes / nrow(al_all),
              ft$estimate, fmt_p(ft$p.value)))
  fisher_results[[pt]] <- list(
    partner_type = pt,
    k_count = k_yes, k_frac = k_yes / nrow(k_all),
    al_count = al_yes, al_frac = al_yes / nrow(al_all),
    OR = ft$estimate, p = ft$p.value
  )
}

# --- Functional group comparison ---
cat("\n  Functional group proportions by enhancer class:\n")
func_table <- table(enh_loop_pairs$enhancer_class, enh_loop_pairs$partner_group)
for (cls in CLASS_ORDER) {
  total <- sum(func_table[cls, ])
  props <- 100 * func_table[cls, ] / total
  cat(sprintf("    %s: Promoter=%.1f%%, Polycomb=%.1f%%, Enhancer=%.1f%%, Structural=%.1f%%\n",
              cls, props["Promoter"], props["Polycomb"],
              props["Enhancer"], props["Structural"]))
}

chi_func <- chisq.test(func_table, simulate.p.value = TRUE, B = 10000)
cat(sprintf("\n  Chi-squared (functional group ~ enhancer class): X2=%.1f, %s\n",
            chi_func$statistic, fmt_p(chi_func$p.value)))


# #############################################################################
# PART E: SUMMARY
# #############################################################################

cat("\n=== PART E: Summary ===\n\n")

summary_file <- file.path(OUTPUT_DIR, "summary.txt")
sink(summary_file)

cat("================================================================================\n")
cat("STEP 14: LOOP TYPE DECOMPOSITION AT ENHANCER SUBSETS\n")
cat("================================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("INPUTS:\n")
cat(sprintf("  Characterized loops: %d\n", nrow(char_loops)))
cat(sprintf("  Enhancers: %d total\n", nrow(enh_all)))
for (cls in CLASS_ORDER) {
  cat(sprintf("    %s: %d\n", cls, sum(enh_all$enhancer_class == cls)))
}

cat(sprintf("\nTOTAL ENHANCER-LOOP PAIRS: %d\n", nrow(enh_loop_pairs)))

cat("\nPARTNER TYPE COMPOSITION:\n")
for (cls in CLASS_ORDER) {
  sub <- comp_table[comp_table$enhancer_class == cls, ]
  sub <- sub[order(-sub$proportion), ]
  n_total <- sum(sub$Freq)
  cat(sprintf("\n  %s (n=%d loops):\n", cls, n_total))
  for (i in seq_len(nrow(sub))) {
    if (sub$Freq[i] > 0) {
      cat(sprintf("    %-22s %4d (%5.1f%%)\n",
                  as.character(sub$partner_type[i]),
                  sub$Freq[i], 100 * sub$proportion[i]))
    }
  }
}

cat("\nFUNCTIONAL GROUP PROPORTIONS:\n")
for (cls in CLASS_ORDER) {
  total <- sum(func_table[cls, ])
  props <- 100 * func_table[cls, ] / total
  cat(sprintf("  %s: Promoter=%.1f%%, Polycomb=%.1f%%, Enhancer=%.1f%%, Structural=%.1f%%\n",
              cls, props["Promoter"], props["Polycomb"],
              props["Enhancer"], props["Structural"]))
}

cat(sprintf("\nChi-squared (partner type ~ class): %s\n", fmt_p(chi_test$p.value)))
cat(sprintf("Chi-squared (functional group ~ class): %s\n", fmt_p(chi_func$p.value)))

cat("\nLOGFC BY PARTNER TYPE (K119ub_Only, n >= 5):\n")
k_summary <- logfc_summary[logfc_summary$enhancer_class == "K119ub_Only" &
                              logfc_summary$n >= 5, ]
k_summary <- k_summary[order(-abs(k_summary$median_logFC)), ]
for (i in seq_len(nrow(k_summary))) {
  row <- k_summary[i, ]
  cat(sprintf("  %-22s n=%d, median logFC=%.4f, %s\n",
              row$partner_type, row$n, row$median_logFC, fmt_p(row$wilcox_p)))
}

cat("\nCONTRIBUTION ANALYSIS (K119ub_Only):\n")
k_contrib <- contrib[contrib$enhancer_class == "K119ub_Only", ]
k_contrib <- k_contrib[order(-abs(k_contrib$contribution)), ]
total_k <- sum(k_contrib$contribution, na.rm = TRUE)
cat(sprintf("  Total weighted logFC: %.4f\n", total_k))
for (i in seq_len(nrow(k_contrib))) {
  row <- k_contrib[i, ]
  if (!is.na(row$contribution)) {
    cat(sprintf("  %-22s weight=%.3f x med=%.4f -> contrib=%.5f (%.0f%%)\n",
                row$partner_type, row$weight, row$median_logFC,
                row$contribution,
                ifelse(total_k != 0, 100 * row$contribution / total_k, 0)))
  }
}

cat("\nK119ub_Only vs Activity_Lost PARTNER ENRICHMENT:\n")
for (pt in names(fisher_results)) {
  r <- fisher_results[[pt]]
  cat(sprintf("  %-22s K119ub=%.1f%% vs ActLost=%.1f%%, OR=%.2f, %s\n",
              r$partner_type, 100 * r$k_frac, 100 * r$al_frac,
              r$OR, fmt_p(r$p)))
}

sink()
cat(sprintf("  Summary written to: %s\n", summary_file))

cat("\n================================================================================\n")
cat("STEP 14 COMPLETE\n")
cat("================================================================================\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat("Plots: 6 panels (PDF + SVG + JPG)\n")
cat("Tables: 4 TSVs + 1 summary\n")
cat("================================================================================\n")
