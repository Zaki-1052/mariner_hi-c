# biomodal/downstream/scripts/viz_sections/section_76_triple_overlap_synapse_chromatin.R
# Section 76: Triple-Overlap & Synapse-Specific Chromatin Analysis
#
# Addresses three questions from Jai:
#   (1) Add explicit p-values to the 73a neuronal vs non-neuronal chromatin violins
#   (2) Repeat ATAC/K27ac/K27me3 analysis for the triple-overlap set
#       (MeCP2-Up AND Coordinated AND Neuronal)
#   (3) Test whether synapse/axon-specific genes (GO: synap|axon only) show
#       stronger chromatin effects than broader neuronal genes
#
# Panels:
#   76a: 73a-style violins (ATAC/K27ac/K27me3) with p-values and n
#   76b: 4-group comparison (triple-overlap / neuronal-only / coordinated-only / rest)
#   76c: Synapse/axon vs broader neuronal vs non-neuronal chromatin
#   76d: K119ub decile enrichment for triple-overlap and synapse/axon genes
#   76e: Composite
#
# Input:
#   tables/diffbind_gene_level_all_marks.tsv, tables/59_quadrant_master.tsv,
#   tables/72_neuronal_gene_set_go_derived.tsv, tables/mecp2_coordinated_genes.tsv,
#   data/k119ub_gene_signal.tsv, org.Mm.eg.db + GO.db
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_76_triple_overlap_synapse_chromatin.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(GO.db)
  library(patchwork)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

cat("================================================================================\n")
cat("SECTION 76: TRIPLE-OVERLAP & SYNAPSE-SPECIFIC CHROMATIN\n")
cat("  P-values, overlap characterization, synapse/axon specificity\n")
cat("================================================================================\n\n")

SEC76_DIR <- file.path(OUTPUT_DIR, "76_triple_overlap_synapse_chromatin")
dir.create(SEC76_DIR, recursive = TRUE, showWarnings = FALSE)

DIFFBIND_PATH    <- file.path(TABLES_DIR, "diffbind_gene_level_all_marks.tsv")
QUADRANT_PATH    <- file.path(TABLES_DIR, "59_quadrant_master.tsv")
NEURONAL_PATH    <- file.path(TABLES_DIR, "72_neuronal_gene_set_go_derived.tsv")
COORDINATED_PATH <- file.path(TABLES_DIR, "mecp2_coordinated_genes.tsv")
K119UB_PATH      <- file.path(BASE_DIR, "data/k119ub_gene_signal.tsv")

stopifnot(
  "diffbind_gene_level_all_marks.tsv not found" = file.exists(DIFFBIND_PATH),
  "59_quadrant_master.tsv not found" = file.exists(QUADRANT_PATH),
  "72_neuronal_gene_set_go_derived.tsv not found" = file.exists(NEURONAL_PATH),
  "mecp2_coordinated_genes.tsv not found" = file.exists(COORDINATED_PATH),
  "k119ub_gene_signal.tsv not found" = file.exists(K119UB_PATH)
)

SYNAPSE_AXON_PATTERN <- "synap|axon"

MARK_COLORS <- c("ATAC" = "#E6AB02", "K27ac" = "#FF7F00", "K27me3" = "#756BB1")

fmt_p <- function(p) {
  if (is.na(p) || !is.finite(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

save_plot <- function(p, name, w = 10, h = 7) {
  save_multiformat_ggplot(p,
    base_path = file.path(SEC76_DIR, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

# =============================================================================
# DATA LOADING AND GENE SET DERIVATION
# =============================================================================

cat("--- Loading data ---\n\n")

db <- read.table(DIFFBIND_PATH, header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, quote = "")
cat(sprintf("  DiffBind all-marks: %d genes\n", nrow(db)))

qm <- read.table(QUADRANT_PATH, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, quote = "")
cat(sprintf("  Quadrant master:    %d genes\n", nrow(qm)))

neuronal_genes <- read.table(NEURONAL_PATH, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, quote = "")$gene
cat(sprintf("  Neuronal genes:     %d\n", length(neuronal_genes)))

coord_genes <- read.table(COORDINATED_PATH, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, quote = "")$gene
cat(sprintf("  Coordinated genes:  %d\n", length(coord_genes)))

k119ub <- read.table(K119UB_PATH, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, quote = "")
k119ub_q <- k119ub[k119ub$gb_signal_class == "quantifiable", ]
cat(sprintf("  K119ub quantifiable:%d genes\n", nrow(k119ub_q)))

# --- Derive synapse/axon gene set (stricter than neuronal) ---

cat("\n--- Deriving synapse/axon gene set from org.Mm.eg.db ---\n")
cat(sprintf("  Pattern: %s\n", SYNAPSE_AXON_PATTERN))

all_go_bp <- AnnotationDbi::select(GO.db,
  keys = AnnotationDbi::keys(GO.db, keytype = "GOID"),
  keytype = "GOID",
  columns = c("GOID", "TERM", "ONTOLOGY")
)
bp_terms <- all_go_bp[all_go_bp$ONTOLOGY == "BP" & !is.na(all_go_bp$TERM), ]

synapse_terms <- bp_terms[grepl(SYNAPSE_AXON_PATTERN, bp_terms$TERM, ignore.case = TRUE), ]
cat(sprintf("  Synapse/axon GO BP terms: %d\n", nrow(synapse_terms)))

syn_entrez <- AnnotationDbi::select(org.Mm.eg.db,
  keys = synapse_terms$GOID,
  keytype = "GOALL",
  columns = c("ENTREZID", "SYMBOL")
)
syn_entrez <- syn_entrez[!is.na(syn_entrez$SYMBOL) & syn_entrez$SYMBOL != "", ]
synapse_genes <- unique(syn_entrez$SYMBOL)
cat(sprintf("  Synapse/axon genes:  %d\n", length(synapse_genes)))

write.table(
  data.frame(gene = synapse_genes, stringsAsFactors = FALSE),
  file.path(TABLES_DIR, "76_synapse_axon_gene_set.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
cat("  Saved: 76_synapse_axon_gene_set.tsv\n")

# --- Build merged analysis table ---

cat("\n--- Building merged table ---\n")

mecp2_up_genes <- qm$gene[!is.na(qm$mecp2_nearest_fdr) &
                           qm$mecp2_nearest_fdr < Q_THRESHOLD &
                           !is.na(qm$mecp2_mean_fold) &
                           qm$mecp2_mean_fold > 0]

df <- merge(
  db[, c("gene", "chr", "start", "end", "chromatin_state",
         "atac_fold", "k27ac_fold", "k27me3_fold", "k119ub_fold")],
  k119ub_q[, c("symbol", "gb_ctrl_signal", "gb_mut_signal")],
  by.x = "gene", by.y = "symbol", all.x = TRUE
)

df$is_neuronal    <- df$gene %in% neuronal_genes
df$is_synapse     <- df$gene %in% synapse_genes
df$is_coordinated <- df$gene %in% coord_genes
df$is_mecp2_up    <- df$gene %in% mecp2_up_genes
df$is_triple      <- df$is_neuronal & df$is_coordinated & df$is_mecp2_up

# K119ub deciles (for genes with signal)
df_with_signal <- df[!is.na(df$gb_ctrl_signal), ]
df$ctrl_decile <- NA_integer_
df$ctrl_decile[!is.na(df$gb_ctrl_signal)] <- as.integer(
  cut(df$gb_ctrl_signal[!is.na(df$gb_ctrl_signal)],
      breaks = quantile(df$gb_ctrl_signal[!is.na(df$gb_ctrl_signal)],
                        probs = seq(0, 1, 0.1)),
      include.lowest = TRUE, labels = FALSE)
)

cat(sprintf("\n  Merged table:       %d genes\n", nrow(df)))
cat(sprintf("  Neuronal:           %d (%.1f%%)\n", sum(df$is_neuronal),
            100 * mean(df$is_neuronal)))
cat(sprintf("  Synapse/axon:       %d (%.1f%%)\n", sum(df$is_synapse),
            100 * mean(df$is_synapse)))
cat(sprintf("  Coordinated:        %d (%.1f%%)\n", sum(df$is_coordinated),
            100 * mean(df$is_coordinated)))
cat(sprintf("  MeCP2 Up:           %d\n", sum(df$is_mecp2_up)))
cat(sprintf("  Triple overlap:     %d\n", sum(df$is_triple)))

triple_genes <- df$gene[df$is_triple]
cat(sprintf("  Triple-overlap genes: %s\n",
            if (length(triple_genes) <= 20) paste(triple_genes, collapse = ", ")
            else paste(c(head(triple_genes, 15), sprintf("... +%d more", length(triple_genes) - 15)), collapse = ", ")))

# Save triple-overlap table
write.table(
  df[, c("gene", "chr", "start", "end", "chromatin_state",
         "atac_fold", "k27ac_fold", "k27me3_fold", "k119ub_fold",
         "gb_ctrl_signal", "gb_mut_signal",
         "is_neuronal", "is_synapse", "is_coordinated", "is_mecp2_up", "is_triple",
         "ctrl_decile")],
  file.path(TABLES_DIR, "76_triple_overlap_genes.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
cat("  Saved: 76_triple_overlap_genes.tsv\n")

# =============================================================================
# 76a: NEURONAL vs NON-NEURONAL VIOLINS WITH P-VALUES
# =============================================================================

cat("\n--- 76a: Neuronal vs non-neuronal chromatin with p-values ---\n")

df$group <- ifelse(df$is_neuronal, "Neuronal", "Non-neuronal")
df$group <- factor(df$group, levels = c("Neuronal", "Non-neuronal"))

mark_long <- rbind(
  data.frame(mark = "ATAC", fold = df$atac_fold, group = df$group, stringsAsFactors = FALSE),
  data.frame(mark = "K27ac", fold = df$k27ac_fold, group = df$group, stringsAsFactors = FALSE),
  data.frame(mark = "K27me3", fold = df$k27me3_fold, group = df$group, stringsAsFactors = FALSE)
)
mark_long <- mark_long[!is.na(mark_long$fold), ]
mark_long$mark <- factor(mark_long$mark, levels = c("ATAC", "K27ac", "K27me3"))

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
    wilcox_p = wt$p.value,
    stringsAsFactors = FALSE
  ))
  cat(sprintf("  %s: neuronal median=%+.4f (n=%d), other=%+.4f (n=%d), %s\n",
              m, median(neur_vals), length(neur_vals),
              median(other_vals), length(other_vals), fmt_p(wt$p.value)))
}

# Build p-value annotation data frame for geom_text
p_labels <- data.frame(
  mark = factor(mark_stats$mark, levels = c("ATAC", "K27ac", "K27me3")),
  label = sprintf("p = %s\nn(neur) = %s",
                   sapply(mark_stats$wilcox_p, function(p) {
                     if (p < 2.2e-16) "< 2.2e-16" else sprintf("%.2e", p)
                   }),
                   format(mark_stats$n_neuronal, big.mark = ",")),
  stringsAsFactors = FALSE
)

# Y position for annotation: top of each facet
y_positions <- mark_long %>%
  dplyr::group_by(mark) %>%
  dplyr::summarise(ymax = quantile(fold, 0.99, na.rm = TRUE), .groups = "drop")
p_labels <- merge(p_labels, y_positions, by = "mark")

p_76a <- ggplot(mark_long, aes(x = group, y = fold, fill = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_violin(alpha = 0.6, show.legend = FALSE, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.3, show.legend = FALSE) +
  geom_text(data = p_labels, aes(x = 1.5, y = ymax, label = label),
            inherit.aes = FALSE, size = 3.2, color = "black", vjust = -0.3) +
  facet_wrap(~ mark, scales = "free_y") +
  scale_fill_manual(values = c("Neuronal" = "#756BB1", "Non-neuronal" = "grey70")) +
  labs(
    title = "Chromatin mark changes in BAP1-KO: neuronal vs non-neuronal",
    subtitle = "DiffBind log2FC per mark | Wilcoxon rank-sum p-values shown",
    x = NULL, y = "DiffBind log2 fold change (mut/ctrl)"
  ) +
  theme_biomodal() +
  theme(strip.text = element_text(size = 12, face = "bold"))

save_plot(p_76a, "76a_neuronal_chromatin_with_pvalues", w = 13, h = 7)

# =============================================================================
# 76b: 4-GROUP COMPARISON (TRIPLE / NEURONAL-ONLY / COORD-ONLY / REST)
# =============================================================================

cat("\n--- 76b: 4-group chromatin comparison ---\n")

df$overlap_group <- dplyr::case_when(
  df$is_triple ~ "Triple overlap",
  df$is_neuronal & !df$is_triple ~ "Neuronal only",
  df$is_coordinated & !df$is_triple ~ "Coordinated only",
  TRUE ~ "Rest of genome"
)
df$overlap_group <- factor(df$overlap_group,
  levels = c("Triple overlap", "Neuronal only", "Coordinated only", "Rest of genome"))

for (grp in levels(df$overlap_group)) {
  cat(sprintf("  %s: n=%d\n", grp, sum(df$overlap_group == grp)))
}

quad_long <- rbind(
  data.frame(mark = "ATAC", fold = df$atac_fold, group = df$overlap_group, stringsAsFactors = FALSE),
  data.frame(mark = "K27ac", fold = df$k27ac_fold, group = df$overlap_group, stringsAsFactors = FALSE),
  data.frame(mark = "K27me3", fold = df$k27me3_fold, group = df$overlap_group, stringsAsFactors = FALSE)
)
quad_long <- quad_long[!is.na(quad_long$fold), ]
quad_long$mark <- factor(quad_long$mark, levels = c("ATAC", "K27ac", "K27me3"))

# Kruskal-Wallis + pairwise Wilcoxon per mark
quad_stats <- data.frame()
for (m in c("ATAC", "K27ac", "K27me3")) {
  sub <- quad_long[quad_long$mark == m, ]
  kt <- kruskal.test(fold ~ group, data = sub)
  cat(sprintf("  %s Kruskal-Wallis: %s\n", m, fmt_p(kt$p.value)))

  pairs <- combn(levels(sub$group), 2, simplify = FALSE)
  for (pr in pairs) {
    v1 <- sub$fold[sub$group == pr[1]]
    v2 <- sub$fold[sub$group == pr[2]]
    if (length(v1) >= 3 && length(v2) >= 3) {
      wt <- wilcox.test(v1, v2)
      quad_stats <- rbind(quad_stats, data.frame(
        mark = m,
        group_a = pr[1], group_b = pr[2],
        n_a = length(v1), n_b = length(v2),
        median_a = median(v1), median_b = median(v2),
        wilcox_p = wt$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
}
quad_stats$p_adj <- p.adjust(quad_stats$wilcox_p, method = "BH")

write.table(quad_stats,
            file.path(TABLES_DIR, "76_chromatin_stats.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 76_chromatin_stats.tsv\n")

p_76b <- ggplot(quad_long, aes(x = group, y = fold, fill = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_violin(alpha = 0.6, show.legend = FALSE, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.3, show.legend = FALSE) +
  facet_wrap(~ mark, scales = "free_y") +
  scale_fill_manual(values = c(
    "Triple overlap" = "#D95F02",
    "Neuronal only" = "#756BB1",
    "Coordinated only" = "#1B9E77",
    "Rest of genome" = "grey80"
  )) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  labs(
    title = "Chromatin remodeling by gene set membership",
    subtitle = sprintf("Triple = Neuronal ∩ Coordinated ∩ MeCP2-Up (%d genes)",
                        sum(df$is_triple)),
    x = NULL, y = "DiffBind log2 fold change (mut/ctrl)"
  ) +
  theme_biomodal() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 9))

save_plot(p_76b, "76b_four_group_chromatin", w = 14, h = 7)

# =============================================================================
# 76c: SYNAPSE/AXON vs BROADER NEURONAL vs NON-NEURONAL
# =============================================================================

cat("\n--- 76c: Synapse/axon vs broader neuronal ---\n")

df$neuro_class <- dplyr::case_when(
  df$is_synapse ~ "Synapse/axon",
  df$is_neuronal & !df$is_synapse ~ "Broader neuronal",
  TRUE ~ "Non-neuronal"
)
df$neuro_class <- factor(df$neuro_class,
  levels = c("Synapse/axon", "Broader neuronal", "Non-neuronal"))

for (cls in levels(df$neuro_class)) {
  cat(sprintf("  %s: n=%d\n", cls, sum(df$neuro_class == cls)))
}

syn_long <- rbind(
  data.frame(mark = "ATAC", fold = df$atac_fold, class = df$neuro_class, stringsAsFactors = FALSE),
  data.frame(mark = "K27ac", fold = df$k27ac_fold, class = df$neuro_class, stringsAsFactors = FALSE),
  data.frame(mark = "K27me3", fold = df$k27me3_fold, class = df$neuro_class, stringsAsFactors = FALSE)
)
syn_long <- syn_long[!is.na(syn_long$fold), ]
syn_long$mark <- factor(syn_long$mark, levels = c("ATAC", "K27ac", "K27me3"))

# Pairwise Wilcoxon for synapse vs broader-neuronal
syn_stats <- data.frame()
for (m in c("ATAC", "K27ac", "K27me3")) {
  sub <- syn_long[syn_long$mark == m, ]
  pairs <- combn(levels(sub$class), 2, simplify = FALSE)
  for (pr in pairs) {
    v1 <- sub$fold[sub$class == pr[1]]
    v2 <- sub$fold[sub$class == pr[2]]
    if (length(v1) >= 3 && length(v2) >= 3) {
      wt <- wilcox.test(v1, v2)
      syn_stats <- rbind(syn_stats, data.frame(
        mark = m,
        group_a = pr[1], group_b = pr[2],
        n_a = length(v1), n_b = length(v2),
        median_a = median(v1), median_b = median(v2),
        wilcox_p = wt$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
}
syn_stats$p_adj <- p.adjust(syn_stats$wilcox_p, method = "BH")

for (i in seq_len(nrow(syn_stats))) {
  r <- syn_stats[i, ]
  cat(sprintf("  %s: %s vs %s — Δ median=%+.4f, %s (adj %s)\n",
              r$mark, r$group_a, r$group_b,
              r$median_a - r$median_b,
              fmt_p(r$wilcox_p), fmt_p(r$p_adj)))
}

write.table(syn_stats,
            file.path(TABLES_DIR, "76_synapse_vs_neuronal_stats.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 76_synapse_vs_neuronal_stats.tsv\n")

p_76c <- ggplot(syn_long, aes(x = class, y = fold, fill = class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_violin(alpha = 0.6, show.legend = FALSE, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.3, show.legend = FALSE) +
  facet_wrap(~ mark, scales = "free_y") +
  scale_fill_manual(values = c(
    "Synapse/axon" = "#D95F02",
    "Broader neuronal" = "#756BB1",
    "Non-neuronal" = "grey70"
  )) +
  labs(
    title = "Chromatin remodeling: synapse/axon vs broader neuronal",
    subtitle = sprintf("Synapse/axon: %d genes (GO: synap|axon) | Broader neuronal: %d",
                        sum(df$is_synapse), sum(df$is_neuronal & !df$is_synapse)),
    x = NULL, y = "DiffBind log2 fold change (mut/ctrl)"
  ) +
  theme_biomodal() +
  theme(strip.text = element_text(size = 12, face = "bold"))

save_plot(p_76c, "76c_synapse_vs_neuronal_chromatin", w = 13, h = 7)

# =============================================================================
# 76d: K119ub DECILE ENRICHMENT — TRIPLE-OVERLAP AND SYNAPSE/AXON
# =============================================================================

cat("\n--- 76d: K119ub decile enrichment ---\n")

df_dec <- df[!is.na(df$ctrl_decile), ]

# Where do triple-overlap and synapse genes sit in K119ub deciles?
dec_summary <- df_dec %>%
  dplyr::group_by(ctrl_decile) %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_triple = sum(is_triple),
    n_synapse = sum(is_synapse),
    n_neuronal = sum(is_neuronal),
    pct_triple = 100 * sum(is_triple) / dplyr::n(),
    pct_synapse = 100 * sum(is_synapse) / dplyr::n(),
    pct_neuronal = 100 * sum(is_neuronal) / dplyr::n(),
    .groups = "drop"
  )

cat("  K119ub decile membership:\n")
for (i in seq_len(nrow(dec_summary))) {
  r <- dec_summary[i, ]
  cat(sprintf("    D%d: n=%d | triple=%d (%.1f%%) | synapse=%d (%.1f%%) | neuronal=%d (%.1f%%)\n",
              r$ctrl_decile, r$n_total, r$n_triple, r$pct_triple,
              r$n_synapse, r$pct_synapse, r$n_neuronal, r$pct_neuronal))
}

# Fisher test: top decile (D10) enrichment for each gene set
fisher_decile <- data.frame()
for (set_name in c("Triple overlap", "Synapse/axon", "Neuronal")) {
  flag <- switch(set_name,
    "Triple overlap" = df_dec$is_triple,
    "Synapse/axon" = df_dec$is_synapse,
    "Neuronal" = df_dec$is_neuronal
  )
  is_top <- df_dec$ctrl_decile == 10

  a <- sum(flag & is_top)
  b <- sum(!flag & is_top)
  c <- sum(flag & !is_top)
  d <- sum(!flag & !is_top)

  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
  fisher_decile <- rbind(fisher_decile, data.frame(
    gene_set = set_name,
    in_top_decile = a,
    total_in_set = sum(flag),
    odds_ratio = as.numeric(ft$estimate),
    ci_lo = ft$conf.int[1],
    ci_hi = ft$conf.int[2],
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  ))
  cat(sprintf("  Top decile %s: %d/%d, OR=%.2f [%.2f, %.2f], %s\n",
              set_name, a, sum(flag), as.numeric(ft$estimate),
              ft$conf.int[1], ft$conf.int[2], fmt_p(ft$p.value)))
}

write.table(fisher_decile,
            file.path(TABLES_DIR, "76_top_decile_fisher.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved: 76_top_decile_fisher.tsv\n")

# Forest plot of ORs
fisher_decile$gene_set <- factor(fisher_decile$gene_set,
  levels = c("Triple overlap", "Synapse/axon", "Neuronal"))

p_76d <- ggplot(fisher_decile, aes(x = gene_set, y = log2(odds_ratio))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 4, color = "#D95F02") +
  geom_errorbar(aes(ymin = log2(ci_lo), ymax = log2(ci_hi)),
                width = 0.2, linewidth = 0.8, color = "#D95F02") +
  geom_text(aes(label = sprintf("OR=%.2f\n%s", odds_ratio,
                                 sapply(p_value, function(p) {
                                   if (p < 0.001) "***"
                                   else if (p < 0.01) "**"
                                   else if (p < 0.05) "*"
                                   else "ns"
                                 }))),
            vjust = -1.2, size = 3.5) +
  coord_flip() +
  labs(
    title = "K119ub top-decile (D10) enrichment",
    subtitle = "Are triple-overlap / synapse genes disproportionately K119ub-high?",
    x = NULL, y = "log2(Odds Ratio)"
  ) +
  theme_biomodal()

save_plot(p_76d, "76d_k119ub_decile_enrichment", w = 9, h = 5)

# =============================================================================
# 76e: COMPOSITE
# =============================================================================

cat("\n--- 76e: Composite ---\n")

p_76e <- (p_76a) / (p_76b | p_76c) / (p_76d) +
  plot_layout(heights = c(3, 3, 2)) +
  plot_annotation(
    title = "Section 76: Triple-Overlap & Synapse-Specific Chromatin Analysis",
    subtitle = sprintf("Neuronal: %s | Synapse/axon: %s | Coordinated: %s | MeCP2-Up: %d | Triple: %d",
                        format(sum(df$is_neuronal), big.mark = ","),
                        format(sum(df$is_synapse), big.mark = ","),
                        format(sum(df$is_coordinated), big.mark = ","),
                        sum(df$is_mecp2_up),
                        sum(df$is_triple)),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5)
    )
  )

save_plot(p_76e, "76_composite", w = 20, h = 20)

cat("\n================================================================================\n")
cat("SECTION 76 COMPLETE\n")
cat("================================================================================\n")
