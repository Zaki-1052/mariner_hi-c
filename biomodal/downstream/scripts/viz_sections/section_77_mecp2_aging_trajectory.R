# biomodal/downstream/scripts/viz_sections/section_77_mecp2_aging_trajectory.R
# Section 77: MeCP2 Developmental Trajectory (Young vs Adult)
#
# Analyzes MeCP2 binding changes from young to adult in both control and mutant,
# comparing the aging trajectories to identify loci that uniquely super-increase
# in the mutant. From mecp2.png: ctrl young-vs-adult has ~13,752 FDR<0.05 loci,
# mut young-vs-adult has ~33,763 FDR<0.05 loci (3-fold more).
#
# Panels:
#   77a: Overview — stacked bar of ctrl vs mut aging (UP/DOWN/NS)
#   77b: Overlap Venn — ctrl aging-UP peaks vs mut aging-UP peaks
#   77c: Mut-specific aging genes — GO enrichment dotplot
#   77d: Fold comparison — scatter of ctrl vs mut aging fold at shared peaks
#   77e: Composite
#
# Input (from Jai, placed in downstream/peaks/):
#   MeCP2_ctrl_adultvsyoung_diffbind_results.txt
#   MeCP2_mut_adultvsyoung_diffbind_results.txt
#
# These are 11-column DiffBind results (no gene annotation):
#   seqnames, start, end, width, strand, Conc, Conc_young, Conc_adult,
#   Fold (log2 young/adult — INVERTED: negative = adult-higher = age-increased)
# The script negates Fold after loading so positive = age-increased.
# Gene annotation is added via TxDb.Mmusculus.UCSC.mm10.knownGene.
#
# If files are not found, this script skips gracefully.
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_77_mecp2_aging_trajectory.R

source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(ggVennDiagram)
  library(patchwork)
  library(GenomicRanges)
  library(ChIPseeker)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

cat("================================================================================\n")
cat("SECTION 77: MeCP2 DEVELOPMENTAL TRAJECTORY (YOUNG vs ADULT)\n")
cat("  Ctrl vs mut aging — shared vs mut-unique loci\n")
cat("================================================================================\n\n")

SEC77_DIR <- file.path(OUTPUT_DIR, "77_mecp2_aging_trajectory")
dir.create(SEC77_DIR, recursive = TRUE, showWarnings = FALSE)

CTRL_AGING_PATH <- file.path(BASE_DIR, "peaks/MeCP2_ctrl_adultvsyoung_diffbind_results.txt")
MUT_AGING_PATH  <- file.path(BASE_DIR, "peaks/MeCP2_mut_adultvsyoung_diffbind_results.txt")
NEURONAL_PATH   <- file.path(TABLES_DIR, "72_neuronal_gene_set_go_derived.tsv")

# --- Check for required data files ---

data_available <- file.exists(CTRL_AGING_PATH) && file.exists(MUT_AGING_PATH)

if (!data_available) {
  cat("  *** SKIPPING SECTION 77 ***\n")
  cat("  Required young-vs-adult MeCP2 DiffBind files not found:\n")
  if (!file.exists(CTRL_AGING_PATH)) cat(sprintf("    MISSING: %s\n", CTRL_AGING_PATH))
  if (!file.exists(MUT_AGING_PATH))  cat(sprintf("    MISSING: %s\n", MUT_AGING_PATH))
  cat("\n  These files should be provided by Jai.\n")
  cat("  Expected at:\n")
  cat(sprintf("    %s\n", CTRL_AGING_PATH))
  cat(sprintf("    %s\n", MUT_AGING_PATH))
  cat("\n  Re-run this script once the files are available.\n")
  cat("================================================================================\n")
  cat("SECTION 77 SKIPPED (data not available)\n")
  cat("================================================================================\n")
  quit(save = "no", status = 0)
}

TXDB <- TxDb.Mmusculus.UCSC.mm10.knownGene

NEURONAL_PATTERN <- "synap|neuron|axon|dendrit|nervous"

fmt_p <- function(p) {
  if (is.na(p) || !is.finite(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

save_plot <- function(p, name, w = 10, h = 7) {
  save_multiformat_ggplot(p,
    base_path = file.path(SEC77_DIR, name),
    width = w, height = h, dpi = 300,
    verbose = TRUE, use_subfolders = TRUE
  )
}

# =============================================================================
# DATA LOADING
# =============================================================================

cat("--- Loading data ---\n\n")

load_aging_diffbind <- function(filepath, label) {
  df <- read.table(filepath, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, quote = "", fill = TRUE)

  # Normalize p-value column
  if (!"p.value" %in% colnames(df) && "p-value" %in% colnames(df)) {
    df$p.value <- df[["p-value"]]
  }

  # Ensure required columns exist
  required <- c("Fold", "FDR")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf("Missing columns in %s: %s", filepath, paste(missing, collapse = ", ")))
  }

  # Negate Fold: raw file has log2(young/adult), we want log2(adult/young)
  # so positive = age-increased MeCP2 binding
  df$Fold <- -df$Fold

  df$aging_status <- dplyr::case_when(
    df$FDR < Q_THRESHOLD & df$Fold > 0 ~ "UP",
    df$FDR < Q_THRESHOLD & df$Fold < 0 ~ "DOWN",
    TRUE ~ "NS"
  )

  n_up   <- sum(df$aging_status == "UP")
  n_down <- sum(df$aging_status == "DOWN")
  n_ns   <- sum(df$aging_status == "NS")
  cat(sprintf("  %s: %d peaks (UP=%d, DOWN=%d, NS=%d)\n",
              label, nrow(df), n_up, n_down, n_ns))

  # If no SYMBOL column, annotate peaks with nearest gene
  if (!"SYMBOL" %in% colnames(df)) {
    cat(sprintf("  Annotating %s peaks with nearest gene (TxDb mm10)...\n", label))
    chr_col <- if ("seqnames" %in% colnames(df)) "seqnames" else "chr"
    start_col <- if ("start" %in% colnames(df)) "start" else "Start"
    end_col   <- if ("end" %in% colnames(df)) "end" else "End"

    gr <- GRanges(
      seqnames = df[[chr_col]],
      ranges = IRanges(start = df[[start_col]], end = df[[end_col]])
    )

    anno <- annotatePeak(gr, TxDb = TXDB, annoDb = "org.Mm.eg.db",
                         level = "gene", verbose = FALSE)
    anno_df <- as.data.frame(anno)

    # annotatePeak can drop peaks that fall outside known chromosomes;
    # match back by position to handle row count mismatch
    if (nrow(anno_df) == nrow(df)) {
      df$SYMBOL <- anno_df$SYMBOL
      df$distanceToTSS <- anno_df$distanceToTSS
      df$annotation <- anno_df$annotation
    } else {
      cat(sprintf("  Warning: annotatePeak returned %d rows for %d peaks — matching by coordinates\n",
                  nrow(anno_df), nrow(df)))
      df$SYMBOL <- NA_character_
      df$distanceToTSS <- NA_integer_
      df$annotation <- NA_character_

      anno_key <- paste(anno_df$seqnames, anno_df$start, anno_df$end, sep = ":")
      df_key   <- paste(df[[chr_col]], df[[start_col]], df[[end_col]], sep = ":")
      match_idx <- match(df_key, anno_key)
      matched <- !is.na(match_idx)
      df$SYMBOL[matched]       <- anno_df$SYMBOL[match_idx[matched]]
      df$distanceToTSS[matched] <- anno_df$distanceToTSS[match_idx[matched]]
      df$annotation[matched]   <- anno_df$annotation[match_idx[matched]]
    }

    n_annotated <- sum(!is.na(df$SYMBOL) & df$SYMBOL != "")
    cat(sprintf("  Annotated: %d / %d peaks mapped to a gene symbol\n",
                n_annotated, nrow(df)))
  }

  df
}

ctrl_aging <- load_aging_diffbind(CTRL_AGING_PATH, "Ctrl aging")
mut_aging  <- load_aging_diffbind(MUT_AGING_PATH, "Mut aging")

neuronal_genes <- read.table(NEURONAL_PATH, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, quote = "")$gene

# =============================================================================
# 77a: AGING OVERVIEW — STACKED BAR
# =============================================================================

cat("\n--- 77a: Aging peak overview ---\n")

aging_summary <- data.frame(
  genotype = rep(c("Control", "Mutant"), each = 3),
  status = rep(c("UP", "DOWN", "NS"), 2),
  count = c(
    sum(ctrl_aging$aging_status == "UP"),
    sum(ctrl_aging$aging_status == "DOWN"),
    sum(ctrl_aging$aging_status == "NS"),
    sum(mut_aging$aging_status == "UP"),
    sum(mut_aging$aging_status == "DOWN"),
    sum(mut_aging$aging_status == "NS")
  ),
  stringsAsFactors = FALSE
)
aging_summary$genotype <- factor(aging_summary$genotype, levels = c("Control", "Mutant"))
aging_summary$status <- factor(aging_summary$status, levels = c("UP", "DOWN", "NS"))

write.table(aging_summary,
            file.path(TABLES_DIR, "77_aging_peak_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

p_77a <- ggplot(aging_summary, aes(x = genotype, y = count, fill = status)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = format(count, big.mark = ",")),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("UP" = "#D95F02", "DOWN" = "#7570B3", "NS" = "grey70")) +
  labs(
    title = "MeCP2 binding changes: young → adult",
    subtitle = "DiffBind young vs adult within each genotype",
    x = NULL, y = "Number of peaks", fill = "Direction"
  ) +
  theme_biomodal()

save_plot(p_77a, "77a_aging_overview", w = 9, h = 7)

# =============================================================================
# 77b: OVERLAP — CTRL vs MUT AGING-UP PEAKS
# =============================================================================

cat("\n--- 77b: Ctrl vs mut aging-UP overlap ---\n")

get_coord_col <- function(df, candidates) {
  for (c in candidates) if (c %in% colnames(df)) return(c)
  stop(sprintf("None of %s found in columns", paste(candidates, collapse = "/")))
}
chr_col   <- get_coord_col(ctrl_aging, c("seqnames", "Chr", "chr"))
start_col <- get_coord_col(ctrl_aging, c("start", "Start"))
end_col   <- get_coord_col(ctrl_aging, c("end", "End"))

ctrl_up <- ctrl_aging[ctrl_aging$aging_status == "UP", ]
mut_up  <- mut_aging[mut_aging$aging_status == "UP", ]

cat(sprintf("  Ctrl aging-UP peaks: %d\n", nrow(ctrl_up)))
cat(sprintf("  Mut aging-UP peaks:  %d\n", nrow(mut_up)))

# Build GRanges for overlap
ctrl_gr <- GRanges(
  seqnames = ctrl_up[[chr_col]],
  ranges = IRanges(start = ctrl_up[[start_col]], end = ctrl_up[[end_col]])
)
mut_gr <- GRanges(
  seqnames = mut_up[[chr_col]],
  ranges = IRanges(start = mut_up[[start_col]], end = mut_up[[end_col]])
)

# Find overlaps
ctrl_hits <- queryHits(findOverlaps(ctrl_gr, mut_gr))
mut_hits  <- subjectHits(findOverlaps(ctrl_gr, mut_gr))

n_ctrl_shared <- length(unique(ctrl_hits))
n_mut_shared  <- length(unique(mut_hits))
n_ctrl_unique <- nrow(ctrl_up) - n_ctrl_shared
n_mut_unique  <- nrow(mut_up) - n_mut_shared

cat(sprintf("  Shared (ctrl peaks with mut overlap): %d / %d (%.1f%%)\n",
            n_ctrl_shared, nrow(ctrl_up), 100 * n_ctrl_shared / nrow(ctrl_up)))
cat(sprintf("  Ctrl-unique: %d\n", n_ctrl_unique))
cat(sprintf("  Mut-unique:  %d\n", n_mut_unique))

# Gene-level overlap for Venn
ctrl_up_genes <- unique(ctrl_up$SYMBOL[!is.na(ctrl_up$SYMBOL) & ctrl_up$SYMBOL != ""])
mut_up_genes  <- unique(mut_up$SYMBOL[!is.na(mut_up$SYMBOL) & mut_up$SYMBOL != ""])

cat(sprintf("  Ctrl aging-UP genes: %d\n", length(ctrl_up_genes)))
cat(sprintf("  Mut aging-UP genes:  %d\n", length(mut_up_genes)))
cat(sprintf("  Gene overlap:        %d\n", length(intersect(ctrl_up_genes, mut_up_genes))))

overlap_info <- data.frame(
  category = c("Ctrl aging-UP peaks", "Mut aging-UP peaks",
                "Ctrl peaks with mut overlap", "Mut peaks with ctrl overlap",
                "Ctrl-unique peaks", "Mut-unique peaks",
                "Ctrl aging-UP genes", "Mut aging-UP genes", "Gene overlap"),
  count = c(nrow(ctrl_up), nrow(mut_up),
            n_ctrl_shared, n_mut_shared,
            n_ctrl_unique, n_mut_unique,
            length(ctrl_up_genes), length(mut_up_genes),
            length(intersect(ctrl_up_genes, mut_up_genes))),
  stringsAsFactors = FALSE
)
write.table(overlap_info,
            file.path(TABLES_DIR, "77_aging_overlap.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Venn diagram at gene level
venn_list <- list(
  "Ctrl aging-UP" = ctrl_up_genes,
  "Mut aging-UP" = mut_up_genes
)

p_77b <- ggVennDiagram(venn_list, label = "count", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#D95F02", guide = "none") +
  labs(
    title = "Aging-increased MeCP2 binding: gene-level overlap",
    subtitle = sprintf("Ctrl: %s genes | Mut: %s genes | Shared: %s",
                        format(length(ctrl_up_genes), big.mark = ","),
                        format(length(mut_up_genes), big.mark = ","),
                        format(length(intersect(ctrl_up_genes, mut_up_genes)), big.mark = ","))
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )

save_plot(p_77b, "77b_aging_overlap_venn", w = 8, h = 7)

# =============================================================================
# 77c: MUT-SPECIFIC AGING GENES — GO ENRICHMENT
# =============================================================================

cat("\n--- 77c: Mut-specific aging GO enrichment ---\n")

mut_unique_genes <- setdiff(mut_up_genes, ctrl_up_genes)
cat(sprintf("  Mut-unique aging genes: %d\n", length(mut_unique_genes)))

# Are mut-unique genes enriched for neuronal genes?
n_mut_unique_neuronal <- sum(mut_unique_genes %in% neuronal_genes)
cat(sprintf("  Neuronal among mut-unique: %d / %d (%.1f%%)\n",
            n_mut_unique_neuronal, length(mut_unique_genes),
            100 * n_mut_unique_neuronal / max(1, length(mut_unique_genes))))

# GO enrichment
mut_unique_entrez <- tryCatch({
  AnnotationDbi::select(org.Mm.eg.db,
    keys = mut_unique_genes,
    keytype = "SYMBOL",
    columns = "ENTREZID"
  )$ENTREZID
}, error = function(e) character(0))
mut_unique_entrez <- unique(mut_unique_entrez[!is.na(mut_unique_entrez)])

if (length(mut_unique_entrez) >= 10) {
  ego <- enrichGO(
    gene = mut_unique_entrez,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  ego_df <- as.data.frame(ego)
  cat(sprintf("  Significant GO BP terms: %d\n", nrow(ego_df)))

  if (nrow(ego_df) > 0) {
    ego_df$is_neuronal <- grepl(NEURONAL_PATTERN, ego_df$Description, ignore.case = TRUE)
    n_neur_terms <- sum(ego_df$is_neuronal)
    cat(sprintf("  Neuronal GO terms among sig: %d / %d\n", n_neur_terms, nrow(ego_df)))

    write.table(ego_df,
                file.path(TABLES_DIR, "77_mut_specific_aging_go.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)

    p_77c <- dotplot(ego, showCategory = min(20, nrow(ego_df))) +
      labs(title = "GO enrichment: mut-specific aging genes",
           subtitle = sprintf("%d genes gaining MeCP2 only in mutant aging",
                               length(mut_unique_genes))) +
      theme_biomodal()

    save_plot(p_77c, "77c_mut_specific_go_enrichment", w = 10, h = 9)
  } else {
    cat("  No significant GO terms — creating placeholder plot\n")
    p_77c <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = sprintf("No significant GO BP terms\n(mut-unique: %d genes)",
                                length(mut_unique_genes)),
               size = 6, hjust = 0.5) +
      theme_void() +
      labs(title = "GO enrichment: mut-specific aging genes")
    save_plot(p_77c, "77c_mut_specific_go_enrichment", w = 10, h = 9)
  }
} else {
  cat(sprintf("  Only %d Entrez IDs — too few for GO enrichment\n", length(mut_unique_entrez)))
  p_77c <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = sprintf("Too few mut-unique genes for enrichment\n(%d Entrez IDs)",
                              length(mut_unique_entrez)),
             size = 6, hjust = 0.5) +
    theme_void() +
    labs(title = "GO enrichment: mut-specific aging genes")
  save_plot(p_77c, "77c_mut_specific_go_enrichment", w = 10, h = 9)
}

# =============================================================================
# 77d: FOLD COMPARISON AT SHARED PEAKS
# =============================================================================

cat("\n--- 77d: Fold comparison at shared loci ---\n")

# Match shared peaks: for each ctrl aging-UP peak that overlaps a mut aging-UP peak,
# pair them and compare fold changes
shared_ctrl_idx <- unique(ctrl_hits)
shared_mut_idx  <- unique(mut_hits)

if (length(shared_ctrl_idx) >= 10) {
  # For each shared ctrl peak, take the first overlapping mut peak
  ol <- findOverlaps(ctrl_gr, mut_gr)
  ol_df <- data.frame(ctrl_idx = queryHits(ol), mut_idx = subjectHits(ol))
  # Deduplicate: one ctrl peak → one mut peak (first match)
  ol_dedup <- ol_df[!duplicated(ol_df$ctrl_idx), ]

  ctrl_sym <- if ("SYMBOL" %in% colnames(ctrl_up)) ctrl_up$SYMBOL else rep(NA_character_, nrow(ctrl_up))
  shared_df <- data.frame(
    ctrl_fold = ctrl_up$Fold[ol_dedup$ctrl_idx],
    mut_fold  = mut_up$Fold[ol_dedup$mut_idx],
    gene = ctrl_sym[ol_dedup$ctrl_idx],
    stringsAsFactors = FALSE
  )
  shared_df <- shared_df[!is.na(shared_df$ctrl_fold) & !is.na(shared_df$mut_fold), ]

  cat(sprintf("  Paired shared peaks: %d\n", nrow(shared_df)))
  cat(sprintf("  Median ctrl aging fold:  %.4f\n", median(shared_df$ctrl_fold)))
  cat(sprintf("  Median mut aging fold:   %.4f\n", median(shared_df$mut_fold)))

  wt_shared <- wilcox.test(shared_df$mut_fold, shared_df$ctrl_fold, paired = TRUE)
  cat(sprintf("  Paired Wilcoxon (mut > ctrl): %s\n", fmt_p(wt_shared$p.value)))

  shared_df$is_neuronal <- shared_df$gene %in% neuronal_genes

  write.table(shared_df,
              file.path(TABLES_DIR, "77_shared_peak_fold_comparison.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # Scatter plot
  fold_range <- range(c(shared_df$ctrl_fold, shared_df$mut_fold), na.rm = TRUE)

  p_77d <- ggplot(shared_df, aes(x = ctrl_fold, y = mut_fold)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = is_neuronal), alpha = 0.3, size = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "#D95F02", linewidth = 0.8) +
    scale_color_manual(values = c("TRUE" = "#756BB1", "FALSE" = "grey50"),
                        labels = c("TRUE" = "Neuronal", "FALSE" = "Other"),
                        name = NULL) +
    coord_fixed(xlim = fold_range, ylim = fold_range) +
    labs(
      title = "MeCP2 aging fold: ctrl vs mut at shared peaks",
      subtitle = sprintf("n=%s shared peaks | Paired Wilcoxon %s\nAbove diagonal = mut ages more than ctrl",
                          format(nrow(shared_df), big.mark = ","),
                          fmt_p(wt_shared$p.value)),
      x = "Ctrl aging fold (log2 adult/young)",
      y = "Mut aging fold (log2 adult/young)"
    ) +
    theme_biomodal() +
    theme(legend.position = c(0.85, 0.15))

  save_plot(p_77d, "77d_shared_peak_fold_comparison", w = 9, h = 9)

} else {
  cat("  Too few shared peaks for fold comparison\n")
  p_77d <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "Too few shared peaks for comparison", size = 6) +
    theme_void() +
    labs(title = "MeCP2 aging fold: ctrl vs mut at shared peaks")
  save_plot(p_77d, "77d_shared_peak_fold_comparison", w = 9, h = 9)
}

# =============================================================================
# 77e: COMPOSITE
# =============================================================================

cat("\n--- 77e: Composite ---\n")

p_77e <- (p_77a | p_77b) / (p_77c | p_77d) +
  plot_annotation(
    title = "Section 77: MeCP2 Developmental Trajectory (Young vs Adult)",
    subtitle = "Mutants show ~2.5x more significant age-related MeCP2 changes than controls",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

save_plot(p_77e, "77_composite", w = 18, h = 16)

cat("\n================================================================================\n")
cat("SECTION 77 COMPLETE\n")
cat("================================================================================\n")
