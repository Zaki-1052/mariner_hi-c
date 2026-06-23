# biomodal/downstream/scripts/viz_sections/section_75_mecp2_signal_reconciliation.R
# Section 75: MeCP2 Signal Direction Reconciliation
#
# Reconciles the apparent paradox: gene-body MeCP2 BigWig signal drops genome-wide
# in the mutant (section 65), but DiffBind identifies 7,686 UP peaks vs only 1,200
# DOWN peaks. Also compares K119ub vs MeCP2 GSEA results — K119ub has 115 sig
# neuronal terms while MeCP2 has only "synapse assembly."
#
# Panels:
#   75a: Peak-to-gene aggregation — UP/DOWN/NS peak counts per gene
#   75b: Gene-body MeCP2 fold by peak status (UP-peak genes vs DOWN vs no sig)
#   75c: GSEA GO term comparison — K119ub vs MeCP2 neuronal term counts
#   75d: MeCP2 UP vs DOWN peak annotation distribution
#   75e: Composite
#
# Input:
#   peaks/mecp2/MeCP2_annotated.txt, tables/59_quadrant_master.tsv,
#   tables/61k_gsea_mecp2_go_bp.tsv, tables/61k_gsea_k119ub_go_bp.tsv,
#   tables/72_neuronal_gene_set_go_derived.tsv
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_75_mecp2_signal_reconciliation.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(patchwork)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

cat("================================================================================\n")
cat("SECTION 75: MeCP2 SIGNAL DIRECTION RECONCILIATION\n")
cat("  Peak-level UP vs gene-body signal drop\n")
cat("================================================================================\n\n")

SEC75_DIR <- file.path(OUTPUT_DIR, "75_mecp2_signal_reconciliation")
dir.create(SEC75_DIR, recursive = TRUE, showWarnings = FALSE)

QUADRANT_PATH  <- file.path(TABLES_DIR, "59_quadrant_master.tsv")
NEURONAL_PATH  <- file.path(TABLES_DIR, "72_neuronal_gene_set_go_derived.tsv")
GSEA_MECP2     <- file.path(TABLES_DIR, "61k_gsea_mecp2_go_bp.tsv")
GSEA_K119UB    <- file.path(TABLES_DIR, "61k_gsea_k119ub_go_bp.tsv")

stopifnot(
  "MeCP2_annotated.txt not found"  = file.exists(MECP2_FILES$annotated),
  "59_quadrant_master.tsv not found" = file.exists(QUADRANT_PATH),
  "61k_gsea_mecp2_go_bp.tsv not found" = file.exists(GSEA_MECP2),
  "61k_gsea_k119ub_go_bp.tsv not found" = file.exists(GSEA_K119UB)
)

NEURONAL_PATTERN <- "synap|neuron|axon|dendrit|nervous"

fmt_p <- function(p) {
  if (is.na(p) || !is.finite(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

save_plot <- function(p, name, w = 10, h = 7) {
  save_multiformat_ggplot(p,
    base_path = file.path(SEC75_DIR, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

# =============================================================================
# DATA LOADING
# =============================================================================

cat("--- Loading data ---\n\n")

# Peak-level MeCP2 DiffBind
mecp2_peaks <- read.table(MECP2_FILES$annotated, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, quote = "", fill = TRUE)
cat(sprintf("  MeCP2 peaks:        %d\n", nrow(mecp2_peaks)))

# Normalize p-value column
if ("p.value" %in% colnames(mecp2_peaks)) {
  # already correct
} else if ("p-value" %in% colnames(mecp2_peaks)) {
  mecp2_peaks$p.value <- mecp2_peaks[["p-value"]]
}

mecp2_peaks$peak_status <- dplyr::case_when(
  mecp2_peaks$FDR < Q_THRESHOLD & mecp2_peaks$Fold > 0 ~ "UP",
  mecp2_peaks$FDR < Q_THRESHOLD & mecp2_peaks$Fold < 0 ~ "DOWN",
  TRUE ~ "NS"
)
mecp2_peaks$peak_status <- factor(mecp2_peaks$peak_status,
                                   levels = c("UP", "DOWN", "NS"))

cat(sprintf("  Peak-level: UP=%d, DOWN=%d, NS=%d\n",
            sum(mecp2_peaks$peak_status == "UP"),
            sum(mecp2_peaks$peak_status == "DOWN"),
            sum(mecp2_peaks$peak_status == "NS")))

# Quadrant master for gene-body info
qm <- read.table(QUADRANT_PATH, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, quote = "")
cat(sprintf("  Quadrant master:    %d genes\n", nrow(qm)))

# Neuronal genes
neuronal_genes <- read.table(NEURONAL_PATH, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, quote = "")$gene

# GSEA results
gsea_mecp2  <- read.table(GSEA_MECP2, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE, quote = "")
gsea_k119ub <- read.table(GSEA_K119UB, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE, quote = "")
cat(sprintf("  GSEA terms: MeCP2=%d, K119ub=%d\n", nrow(gsea_mecp2), nrow(gsea_k119ub)))

# =============================================================================
# AGGREGATE PEAKS TO GENE LEVEL
# =============================================================================

cat("\n--- Aggregating peaks to gene level ---\n\n")

peak_by_gene <- mecp2_peaks %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(
    n_peaks     = dplyr::n(),
    n_up        = sum(peak_status == "UP"),
    n_down      = sum(peak_status == "DOWN"),
    n_ns        = sum(peak_status == "NS"),
    mean_fold   = mean(Fold, na.rm = TRUE),
    min_fdr     = min(FDR, na.rm = TRUE),
    .groups = "drop"
  )

peak_by_gene$gene_peak_class <- dplyr::case_when(
  peak_by_gene$n_up > 0  ~ "Has UP peak(s)",
  peak_by_gene$n_down > 0 ~ "Has DOWN peak(s) only",
  TRUE ~ "No significant peaks"
)

n_genes_up   <- sum(peak_by_gene$gene_peak_class == "Has UP peak(s)")
n_genes_down <- sum(peak_by_gene$gene_peak_class == "Has DOWN peak(s) only")
n_genes_ns   <- sum(peak_by_gene$gene_peak_class == "No significant peaks")

cat(sprintf("  Genes with at least 1 MeCP2 peak: %d\n", nrow(peak_by_gene)))
cat(sprintf("  Genes with UP peak(s):     %d (carrying %d UP peaks)\n",
            n_genes_up, sum(peak_by_gene$n_up)))
cat(sprintf("  Genes with DOWN only:      %d (carrying %d DOWN peaks)\n",
            n_genes_down, sum(peak_by_gene$n_down[peak_by_gene$gene_peak_class == "Has DOWN peak(s) only"])))
cat(sprintf("  Genes with no sig peaks:   %d\n", n_genes_ns))

# Merge with quadrant master
gene_df <- merge(
  qm[, c("gene", "mecp2_mean_fold", "mecp2_nearest_fdr", "mecp2_n_peaks")],
  peak_by_gene,
  by.x = "gene", by.y = "SYMBOL", all.x = TRUE
)
gene_df$gene_peak_class[is.na(gene_df$gene_peak_class)] <- "No MeCP2 peaks"
gene_df$gene_peak_class <- factor(gene_df$gene_peak_class,
  levels = c("Has UP peak(s)", "Has DOWN peak(s) only",
             "No significant peaks", "No MeCP2 peaks"))

# Save
write.table(peak_by_gene,
            file.path(TABLES_DIR, "75_mecp2_peak_gene_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 75_mecp2_peak_gene_summary.tsv\n")

# =============================================================================
# 75a: PEAK DISTRIBUTION ACROSS GENES
# =============================================================================

cat("\n--- 75a: Peak distribution across genes ---\n")

class_counts <- as.data.frame(table(gene_df$gene_peak_class))
colnames(class_counts) <- c("class", "n_genes")

p_75a <- ggplot(class_counts, aes(x = class, y = n_genes, fill = class)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = format(n_genes, big.mark = ",")),
            vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = c(
    "Has UP peak(s)" = "#D95F02",
    "Has DOWN peak(s) only" = "#7570B3",
    "No significant peaks" = "grey70",
    "No MeCP2 peaks" = "grey90"
  )) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  labs(
    title = "MeCP2 peak status per gene",
    subtitle = sprintf("Peak-level: %s UP, %s DOWN across %s unique genes with peaks",
                        format(sum(mecp2_peaks$peak_status == "UP"), big.mark = ","),
                        format(sum(mecp2_peaks$peak_status == "DOWN"), big.mark = ","),
                        format(nrow(peak_by_gene), big.mark = ",")),
    x = NULL, y = "Number of genes"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(size = 10))

save_plot(p_75a, "75a_peak_distribution_per_gene", w = 9, h = 7)

# =============================================================================
# 75b: GENE-BODY MeCP2 FOLD BY PEAK STATUS
# =============================================================================

cat("\n--- 75b: Gene-body MeCP2 fold by peak status ---\n")

plot_df <- gene_df[!is.na(gene_df$mecp2_mean_fold) &
                    gene_df$gene_peak_class != "No MeCP2 peaks", ]

for (cls in levels(plot_df$gene_peak_class)) {
  vals <- plot_df$mecp2_mean_fold[plot_df$gene_peak_class == cls]
  if (length(vals) > 0) {
    cat(sprintf("  %s: n=%d, median fold=%.4f\n", cls, length(vals), median(vals)))
  }
}

# Pairwise Wilcoxon between UP-peak genes and the other two groups
up_vals   <- plot_df$mecp2_mean_fold[plot_df$gene_peak_class == "Has UP peak(s)"]
down_vals <- plot_df$mecp2_mean_fold[plot_df$gene_peak_class == "Has DOWN peak(s) only"]
ns_vals   <- plot_df$mecp2_mean_fold[plot_df$gene_peak_class == "No significant peaks"]

if (length(up_vals) > 0 && length(down_vals) > 0) {
  wt_up_down <- wilcox.test(up_vals, down_vals)
  cat(sprintf("  UP vs DOWN: %s\n", fmt_p(wt_up_down$p.value)))
}
if (length(up_vals) > 0 && length(ns_vals) > 0) {
  wt_up_ns <- wilcox.test(up_vals, ns_vals)
  cat(sprintf("  UP vs NS: %s\n", fmt_p(wt_up_ns$p.value)))
}

p_75b <- ggplot(plot_df, aes(x = gene_peak_class, y = mecp2_mean_fold,
                              fill = gene_peak_class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_violin(alpha = 0.6, show.legend = FALSE, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = c(
    "Has UP peak(s)" = "#D95F02",
    "Has DOWN peak(s) only" = "#7570B3",
    "No significant peaks" = "grey70"
  )) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  labs(
    title = "Gene-body MeCP2 signal by peak status",
    subtitle = "mecp2_mean_fold from quadrant master (concentration-weighted gene-level aggregate)",
    x = NULL, y = "MeCP2 mean fold (log2 mut/ctrl)"
  ) +
  theme_biomodal() +
  theme(axis.text.x = element_text(size = 10))

save_plot(p_75b, "75b_genebody_signal_by_peak_status", w = 9, h = 7)

# =============================================================================
# 75c: GSEA GO TERM COMPARISON — K119ub vs MeCP2
# =============================================================================

cat("\n--- 75c: GSEA term comparison ---\n")

gsea_mecp2$sig  <- gsea_mecp2$p.adjust < Q_THRESHOLD
gsea_k119ub$sig <- gsea_k119ub$p.adjust < Q_THRESHOLD

gsea_mecp2$is_neuronal  <- grepl(NEURONAL_PATTERN, gsea_mecp2$Description, ignore.case = TRUE)
gsea_k119ub$is_neuronal <- grepl(NEURONAL_PATTERN, gsea_k119ub$Description, ignore.case = TRUE)

# Counts
counts_df <- data.frame(
  dataset = rep(c("K119ub GSEA", "MeCP2 GSEA"), each = 2),
  category = rep(c("Neuronal terms", "Other terms"), 2),
  count = c(
    sum(gsea_k119ub$sig & gsea_k119ub$is_neuronal),
    sum(gsea_k119ub$sig & !gsea_k119ub$is_neuronal),
    sum(gsea_mecp2$sig & gsea_mecp2$is_neuronal),
    sum(gsea_mecp2$sig & !gsea_mecp2$is_neuronal)
  ),
  stringsAsFactors = FALSE
)

cat(sprintf("  K119ub sig terms:  %d total (%d neuronal)\n",
            sum(gsea_k119ub$sig), sum(gsea_k119ub$sig & gsea_k119ub$is_neuronal)))
cat(sprintf("  MeCP2 sig terms:   %d total (%d neuronal)\n",
            sum(gsea_mecp2$sig), sum(gsea_mecp2$sig & gsea_mecp2$is_neuronal)))

# List the specific MeCP2 neuronal terms
mecp2_neur_terms <- gsea_mecp2[gsea_mecp2$sig & gsea_mecp2$is_neuronal, ]
if (nrow(mecp2_neur_terms) > 0) {
  cat("  MeCP2 neuronal terms:\n")
  for (i in seq_len(nrow(mecp2_neur_terms))) {
    cat(sprintf("    - %s (NES=%.2f, q=%.3e)\n",
                mecp2_neur_terms$Description[i],
                mecp2_neur_terms$NES[i],
                mecp2_neur_terms$p.adjust[i]))
  }
}

# Fisher test: neuronal fraction among sig terms
#   K119ub: neuronal_sig / total_sig  vs  MeCP2: neuronal_sig / total_sig
k_neur <- sum(gsea_k119ub$sig & gsea_k119ub$is_neuronal)
k_other <- sum(gsea_k119ub$sig & !gsea_k119ub$is_neuronal)
m_neur <- sum(gsea_mecp2$sig & gsea_mecp2$is_neuronal)
m_other <- sum(gsea_mecp2$sig & !gsea_mecp2$is_neuronal)

if ((k_neur + k_other) > 0 && (m_neur + m_other) > 0) {
  ft_terms <- fisher.test(matrix(c(k_neur, m_neur, k_other, m_other), nrow = 2))
  cat(sprintf("  Fisher (neuronal fraction K119ub vs MeCP2): OR=%.2f, %s\n",
              as.numeric(ft_terms$estimate), fmt_p(ft_terms$p.value)))
}

# Save
write.table(counts_df,
            file.path(TABLES_DIR, "75_go_term_comparison.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

counts_df$dataset <- factor(counts_df$dataset, levels = c("K119ub GSEA", "MeCP2 GSEA"))
counts_df$category <- factor(counts_df$category, levels = c("Neuronal terms", "Other terms"))

p_75c <- ggplot(counts_df, aes(x = dataset, y = count, fill = category)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5),
            size = 5, fontface = "bold", color = "white") +
  scale_fill_manual(values = c("Neuronal terms" = "#756BB1", "Other terms" = "grey70")) +
  labs(
    title = "Significant GSEA GO BP terms: K119ub vs MeCP2",
    subtitle = sprintf("K119ub: %d sig (%d neuronal) | MeCP2: %d sig (%d neuronal)",
                        sum(gsea_k119ub$sig), k_neur, sum(gsea_mecp2$sig), m_neur),
    x = NULL, y = "Number of significant GO terms (q < 0.05)",
    fill = NULL
  ) +
  theme_biomodal() +
  theme(legend.position = "top")

save_plot(p_75c, "75c_gsea_term_comparison", w = 8, h = 7)

# =============================================================================
# 75d: PEAK ANNOTATION DISTRIBUTION — UP vs DOWN
# =============================================================================

cat("\n--- 75d: Peak annotation distribution ---\n")

sig_peaks <- mecp2_peaks[mecp2_peaks$peak_status %in% c("UP", "DOWN"), ]

sig_peaks$annotation_broad <- dplyr::case_when(
  grepl("^Promoter", sig_peaks$annotation)  ~ "Promoter",
  grepl("Intron", sig_peaks$annotation)      ~ "Intron",
  grepl("Exon|UTR", sig_peaks$annotation)    ~ "Exon/UTR",
  grepl("Intergenic", sig_peaks$annotation)  ~ "Distal Intergenic",
  grepl("Downstream", sig_peaks$annotation)  ~ "Downstream",
  TRUE ~ "Other"
)

ann_order <- c("Promoter", "Exon/UTR", "Intron", "Distal Intergenic", "Downstream", "Other")
sig_peaks$annotation_broad <- factor(sig_peaks$annotation_broad, levels = ann_order)

ann_counts <- sig_peaks %>%
  dplyr::count(peak_status, annotation_broad) %>%
  dplyr::group_by(peak_status) %>%
  dplyr::mutate(pct = 100 * n / sum(n)) %>%
  dplyr::ungroup()

cat("  UP peak annotations:\n")
for (i in seq_len(nrow(ann_counts[ann_counts$peak_status == "UP", ]))) {
  r <- ann_counts[ann_counts$peak_status == "UP", ][i, ]
  cat(sprintf("    %-20s %5d (%.1f%%)\n", as.character(r$annotation_broad), r$n, r$pct))
}
cat("  DOWN peak annotations:\n")
for (i in seq_len(nrow(ann_counts[ann_counts$peak_status == "DOWN", ]))) {
  r <- ann_counts[ann_counts$peak_status == "DOWN", ][i, ]
  cat(sprintf("    %-20s %5d (%.1f%%)\n", as.character(r$annotation_broad), r$n, r$pct))
}

# Chi-squared test on annotation x direction (drop empty margins)
ann_table <- table(droplevels(sig_peaks$annotation_broad),
                   droplevels(sig_peaks$peak_status))
ann_table <- ann_table[rowSums(ann_table) > 0, colSums(ann_table) > 0, drop = FALSE]
chi_test <- chisq.test(ann_table, simulate.p.value = TRUE, B = 10000)
cat(sprintf("  Chi-squared (annotation x direction): X2=%.1f, %s (simulated)\n",
            chi_test$statistic, fmt_p(chi_test$p.value)))

write.table(ann_counts,
            file.path(TABLES_DIR, "75_peak_annotation_distribution.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

p_75d <- ggplot(ann_counts, aes(x = peak_status, y = pct, fill = annotation_broad)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "MeCP2 UP vs DOWN peak annotation",
    subtitle = sprintf("Chi² %s | UP=%s peaks, DOWN=%s peaks",
                        fmt_p(chi_test$p.value),
                        format(sum(sig_peaks$peak_status == "UP"), big.mark = ","),
                        format(sum(sig_peaks$peak_status == "DOWN"), big.mark = ",")),
    x = "Peak direction (FDR < 0.05)", y = "Percentage of peaks",
    fill = "Annotation"
  ) +
  theme_biomodal() +
  theme(legend.position = "right")

save_plot(p_75d, "75d_peak_annotation_distribution", w = 9, h = 7)

# =============================================================================
# 75e: COMPOSITE
# =============================================================================

cat("\n--- 75e: Composite ---\n")

p_75e <- (p_75a | p_75b) / (p_75c | p_75d) +
  plot_annotation(
    title = "Section 75: MeCP2 Signal Direction Reconciliation",
    subtitle = "Peak-level DiffBind shows 7,686 UP peaks, but gene-body signal drops genome-wide",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

save_plot(p_75e, "75_composite", w = 18, h = 14)

cat("\n================================================================================\n")
cat("SECTION 75 COMPLETE\n")
cat("================================================================================\n")
