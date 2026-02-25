# biomodal/downstream/scripts/viz_sections/section_21_discordant_mc_hmc_analysis.R
# Section 21: Discordant mC dn/hmC up Gene Characterization
# Standalone script - sources shared config for all dependencies and data
#
# Exploratory characterization of the minority discordant (mC dn/hmC up) genes.
# Main analysis: Discordant (Q2) vs Coordinated (Q4).
# Supplementary: All 4 quadrants compared.
#
# Analysis dimensions (modeled after abc/scripts/step13_discordant_gene_analysis.R):
#   1. Methylation effect size
#   2. RNA-seq expression
#   3. ATAC-seq accessibility
#   4. H2AK119ub signal
#   5. H3K27ac status
#   6. Chromatin state (7-category)
#   7. MeCP2 binding
#   8. Hi-C loop involvement
#   9. GO enrichment
#
# Figures:
#   21a: Composite 3x3 grid (Discordant vs Coordinated)
#   21b: mc_diff vs hmc_diff scatter colored by concordance
#   21c: mc_diff vs log2FC per group
#   21d: GO enrichment dot plot
#   21e: 4-quadrant comprehensive panel
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_21_discordant_mc_hmc_analysis.R

# Run from downstream/ directory
source("scripts/viz_sections/_shared_config.R")

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(readxl)
})

# =============================================================================
# SECTION 21 CONFIGURATION
# =============================================================================

RNA_SEQ_FILE <- "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx"
K119UB_GENE_SIGNAL_FILE <- file.path(BASE_DIR, "data/k119ub_gene_signal.tsv")

CONCORDANCE_COLORS <- c("Discordant (Q2)" = "#E41A1C", "Coordinated (Q4)" = "#4DAF4A")
QUADRANT_COLORS <- c("Q1: mC dn/hmC dn" = "#377EB8",
                      "Q2: mC dn/hmC up" = "#E41A1C",
                      "Q3: mC up/hmC up" = "#FF7F00",
                      "Q4: mC up/hmC dn" = "#4DAF4A")

# Helper for formatting p-values
fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.2e", p)
}

# =============================================================================
# VALIDATE INPUTS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 21: DISCORDANT mC dn/hmC up GENE CHARACTERIZATION\n")
cat("================================================================================\n\n")

stopifnot(file.exists(RNA_SEQ_FILE))
stopifnot(file.exists(K119UB_GENE_SIGNAL_FILE))
cat("  Input files validated.\n")

# =============================================================================
# RECOMPUTE COORDINATED CHANGES + DEFINE QUADRANTS
# =============================================================================

cat("\nRecomputing coordinated mC/hmC changes and defining quadrants...\n")

mc_sig <- mc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::select(gene, mc_diff = mod_difference, mc_q = dmr_qvalue,
                chr, start, end)

hmc_sig <- hmc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::select(gene, hmc_diff = mod_difference, hmc_q = dmr_qvalue)

coordinated_all <- inner_join(mc_sig, hmc_sig, by = "gene") %>%
  mutate(
    quadrant = case_when(
      mc_diff < 0 & hmc_diff < 0 ~ "Q1: mC dn/hmC dn",
      mc_diff < 0 & hmc_diff > 0 ~ "Q2: mC dn/hmC up",
      mc_diff > 0 & hmc_diff > 0 ~ "Q3: mC up/hmC up",
      mc_diff > 0 & hmc_diff < 0 ~ "Q4: mC up/hmC dn"
    ),
    combined_effect = abs(mc_diff) + abs(hmc_diff)
  )

# Define comparison groups
disc_genes <- coordinated_all %>% dplyr::filter(quadrant == "Q2: mC dn/hmC up")
coord_genes <- coordinated_all %>% dplyr::filter(quadrant == "Q4: mC up/hmC dn")

cat(sprintf("  Q1 (mC dn/hmC dn): %d genes\n", sum(coordinated_all$quadrant == "Q1: mC dn/hmC dn")))
cat(sprintf("  Q2 (mC dn/hmC up, DISCORDANT): %d genes\n", nrow(disc_genes)))
cat(sprintf("  Q3 (mC up/hmC up): %d genes\n", sum(coordinated_all$quadrant == "Q3: mC up/hmC up")))
cat(sprintf("  Q4 (mC up/hmC dn, COORDINATED): %d genes\n", nrow(coord_genes)))

# =============================================================================
# LOAD ALL DATA SOURCES
# =============================================================================

cat("\nLoading data sources...\n")

# RNA-seq
cat("  Loading RNA-seq...\n")
rna_all <- read_excel(RNA_SEQ_FILE) %>%
  dplyr::select(gene = ensembl_gene_id, log2FC = log2FoldChange, padj, baseMean) %>%
  dplyr::filter(!is.na(gene) & gene != "")
cat(sprintf("    %d genes\n", nrow(rna_all)))

# H2AK119ub gene signal
cat("  Loading K119ub gene signal...\n")
k119ub_signal <- read.table(K119UB_GENE_SIGNAL_FILE, sep = "\t", header = TRUE,
                            stringsAsFactors = FALSE)
cat(sprintf("    %d genes\n", nrow(k119ub_signal)))

# ATAC-seq differential peaks
cat("  Loading ATAC-seq differential peaks...\n")
atac_up_gr <- load_chip_peaks(ATAC_FILES$up, "ATAC Up")
atac_down_gr <- load_chip_peaks(ATAC_FILES$down, "ATAC Down")
stopifnot(!is.null(atac_up_gr) && !is.null(atac_down_gr))

# H3K27ac condition-specific peaks
cat("  Loading H3K27ac condition-specific peaks...\n")
k27ac_ctrl_gr <- load_chip_peaks(H3K27AC_FILES$ctrl, "H3K27ac Ctrl")
k27ac_mut_gr <- load_chip_peaks(H3K27AC_FILES$mut, "H3K27ac Mut")
stopifnot(!is.null(k27ac_ctrl_gr) && !is.null(k27ac_mut_gr))

# K119ub condition-specific peaks
cat("  Loading K119ub condition-specific peaks...\n")
k119ub_ctrl_gr <- load_chip_peaks(K119UB_FILES$ctrl, "K119ub Ctrl")
k119ub_mut_gr <- load_chip_peaks(K119UB_FILES$mut, "K119ub Mut")
stopifnot(!is.null(k119ub_ctrl_gr) && !is.null(k119ub_mut_gr))

# MeCP2 differential peaks
cat("  Loading MeCP2 differential peaks...\n")
mecp2_up_gr <- load_chip_peaks(MECP2_FILES$up, "MeCP2 Up")
mecp2_down_gr <- load_chip_peaks(MECP2_FILES$down, "MeCP2 Down")
stopifnot(!is.null(mecp2_up_gr) && !is.null(mecp2_down_gr))

# ChIP-seq peaks for chromatin state
cat("  Loading ChIP-seq peaks for chromatin state...\n")
chip_peaks <- list()
for (mark in names(CHIP_PEAK_FILES)) {
  chip_peaks[[mark]] <- load_chip_peaks(CHIP_PEAK_FILES[[mark]], mark)
}
stopifnot(all(!sapply(chip_peaks, is.null)))

# Hi-C loops
cat("  Loading Hi-C loops...\n")
loops <- read.table(LOOP_FILES$late, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("    %d loops\n", nrow(loops)))

# TxDb for gene annotation
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

cat("  All data sources loaded.\n")

# =============================================================================
# ANNOTATE PEAKS TO GENES (ATAC, MeCP2)
# =============================================================================

cat("\nAnnotating peaks to genes via ChIPseeker...\n")

# ATAC peaks to genes
cat("  ATAC Up peaks...\n")
atac_up_anno <- suppressMessages(annotatePeak(atac_up_gr, TxDb = txdb,
                                               annoDb = "org.Mm.eg.db", level = "gene"))
atac_up_gene <- as.data.frame(atac_up_anno) %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  group_by(SYMBOL) %>%
  summarise(n_atac_up = n(), .groups = "drop")

cat("  ATAC Down peaks...\n")
atac_down_anno <- suppressMessages(annotatePeak(atac_down_gr, TxDb = txdb,
                                                 annoDb = "org.Mm.eg.db", level = "gene"))
atac_down_gene <- as.data.frame(atac_down_anno) %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  group_by(SYMBOL) %>%
  summarise(n_atac_down = n(), .groups = "drop")

atac_gene <- full_join(atac_up_gene, atac_down_gene, by = "SYMBOL") %>%
  mutate(
    n_atac_up = replace_na(n_atac_up, 0),
    n_atac_down = replace_na(n_atac_down, 0),
    net_atac = n_atac_up - n_atac_down
  )
cat(sprintf("    %d genes with ATAC data\n", nrow(atac_gene)))

# MeCP2 peaks to genes
cat("  MeCP2 Up peaks...\n")
mecp2_up_anno <- suppressMessages(annotatePeak(mecp2_up_gr, TxDb = txdb,
                                                annoDb = "org.Mm.eg.db", level = "gene"))
mecp2_up_gene <- as.data.frame(mecp2_up_anno) %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  group_by(SYMBOL) %>%
  summarise(n_mecp2_up = n(), .groups = "drop")

cat("  MeCP2 Down peaks...\n")
mecp2_down_anno <- suppressMessages(annotatePeak(mecp2_down_gr, TxDb = txdb,
                                                  annoDb = "org.Mm.eg.db", level = "gene"))
mecp2_down_gene <- as.data.frame(mecp2_down_anno) %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  group_by(SYMBOL) %>%
  summarise(n_mecp2_down = n(), .groups = "drop")

mecp2_gene <- full_join(mecp2_up_gene, mecp2_down_gene, by = "SYMBOL") %>%
  mutate(
    n_mecp2_up = replace_na(n_mecp2_up, 0),
    n_mecp2_down = replace_na(n_mecp2_down, 0),
    net_mecp2 = n_mecp2_up - n_mecp2_down,
    mecp2_status = case_when(
      n_mecp2_up > 0 & n_mecp2_down == 0 ~ "MeCP2 Up Only",
      n_mecp2_down > 0 & n_mecp2_up == 0 ~ "MeCP2 Down Only",
      n_mecp2_up > 0 & n_mecp2_down > 0 ~ "MeCP2 Both",
      TRUE ~ "No MeCP2"
    )
  )
cat(sprintf("    %d genes with MeCP2 data\n", nrow(mecp2_gene)))

# =============================================================================
# DERIVE H3K27ac GAINED/LOST STATUS PER GENE
# =============================================================================

cat("\nDeriving H3K27ac gained/lost status...\n")

# Gained = mut-only, Lost = ctrl-only
k27ac_mut_hits <- countOverlaps(k27ac_mut_gr, k27ac_ctrl_gr)
k27ac_ctrl_hits <- countOverlaps(k27ac_ctrl_gr, k27ac_mut_gr)
k27ac_gained_gr <- k27ac_mut_gr[k27ac_mut_hits == 0]
k27ac_lost_gr <- k27ac_ctrl_gr[k27ac_ctrl_hits == 0]

cat(sprintf("  H3K27ac gained: %d peaks\n", length(k27ac_gained_gr)))
cat(sprintf("  H3K27ac lost: %d peaks\n", length(k27ac_lost_gr)))

# Annotate to genes
k27ac_gained_anno <- suppressMessages(annotatePeak(k27ac_gained_gr, TxDb = txdb,
                                                    annoDb = "org.Mm.eg.db", level = "gene"))
k27ac_gained_gene <- as.data.frame(k27ac_gained_anno) %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  group_by(SYMBOL) %>%
  summarise(n_k27ac_gained = n(), .groups = "drop")

k27ac_lost_anno <- suppressMessages(annotatePeak(k27ac_lost_gr, TxDb = txdb,
                                                  annoDb = "org.Mm.eg.db", level = "gene"))
k27ac_lost_gene <- as.data.frame(k27ac_lost_anno) %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  group_by(SYMBOL) %>%
  summarise(n_k27ac_lost = n(), .groups = "drop")

k27ac_gene <- full_join(k27ac_gained_gene, k27ac_lost_gene, by = "SYMBOL") %>%
  mutate(
    n_k27ac_gained = replace_na(n_k27ac_gained, 0),
    n_k27ac_lost = replace_na(n_k27ac_lost, 0),
    net_k27ac = n_k27ac_gained - n_k27ac_lost,
    k27ac_status = case_when(
      n_k27ac_gained > 0 & n_k27ac_lost == 0 ~ "Gained Only",
      n_k27ac_lost > 0 & n_k27ac_gained == 0 ~ "Lost Only",
      n_k27ac_gained > 0 & n_k27ac_lost > 0 ~ "Both",
      TRUE ~ "No Change"
    )
  )
cat(sprintf("  %d genes with H3K27ac differential data\n", nrow(k27ac_gene)))

# =============================================================================
# CHROMATIN STATE CLASSIFICATION
# =============================================================================

cat("\nClassifying chromatin state at gene bodies...\n")

# Create GRanges from coordinated gene coordinates
coord_all_gr <- GRanges(
  seqnames = coordinated_all$chr,
  ranges = IRanges(start = coordinated_all$start, end = coordinated_all$end)
)

# Compute ChIP overlaps
chip_overlaps <- compute_chip_overlaps(coord_all_gr, chip_peaks)

# Gene bodies contain their own TSS, so distance_to_tss = 0
# This means classifications will be promoter-centric (Active/Repressed/Bivalent_Promoter + Other)
distance_to_tss <- rep(0, nrow(coordinated_all))

coordinated_all$chromatin_state <- classify_chromatin_state(
  chip_overlaps, distance_to_tss, TSS_THRESHOLD
)

cat("  Chromatin state distribution:\n")
print(table(coordinated_all$chromatin_state))

# =============================================================================
# HI-C LOOP INVOLVEMENT
# =============================================================================

cat("\nDetermining Hi-C loop involvement...\n")

# Create GRanges for loop anchors (both anchors)
anchor1_gr <- GRanges(seqnames = loops$chr1,
                      ranges = IRanges(start = loops$start1, end = loops$end1))
anchor2_gr <- GRanges(seqnames = loops$chr2,
                      ranges = IRanges(start = loops$start2, end = loops$end2))

# Extend anchors by 5kb to capture nearby genes
anchor1_ext <- resize(anchor1_gr, width = width(anchor1_gr) + 10000, fix = "center")
anchor2_ext <- resize(anchor2_gr, width = width(anchor2_gr) + 10000, fix = "center")

# Check overlap of gene TSS with loop anchors
# Use gene start as approximate TSS
gene_tss_gr <- GRanges(
  seqnames = coordinated_all$chr,
  ranges = IRanges(start = coordinated_all$start, width = 1)
)

at_anchor1 <- countOverlaps(gene_tss_gr, anchor1_ext) > 0
at_anchor2 <- countOverlaps(gene_tss_gr, anchor2_ext) > 0
coordinated_all$loop_involved <- at_anchor1 | at_anchor2

# Also get loop direction for involved genes
loop_info_per_gene <- data.frame(
  gene = coordinated_all$gene,
  at_anchor1 = at_anchor1,
  at_anchor2 = at_anchor2,
  loop_involved = at_anchor1 | at_anchor2,
  stringsAsFactors = FALSE
)

cat(sprintf("  Genes at loop anchors: %d / %d (%.1f%%)\n",
            sum(coordinated_all$loop_involved), nrow(coordinated_all),
            100 * mean(coordinated_all$loop_involved)))

# =============================================================================
# BUILD MASTER GENE TABLE
# =============================================================================

cat("\nBuilding master gene characteristics table...\n")

master <- coordinated_all %>%
  # RNA-seq
  left_join(rna_all, by = "gene") %>%
  # ATAC
  left_join(atac_gene %>% dplyr::rename(gene = SYMBOL), by = "gene") %>%
  mutate(
    n_atac_up = replace_na(n_atac_up, 0),
    n_atac_down = replace_na(n_atac_down, 0),
    net_atac = replace_na(net_atac, 0)
  ) %>%
  # K119ub gene signal
  left_join(k119ub_signal %>% dplyr::select(gene = symbol, k119ub_gb_log2fc = gb_log2fc,
                                             k119ub_pr_log2fc = pr_log2fc,
                                             k119ub_gb_class = gb_signal_class),
            by = "gene") %>%
  # H3K27ac
  left_join(k27ac_gene %>% dplyr::rename(gene = SYMBOL), by = "gene") %>%
  mutate(
    n_k27ac_gained = replace_na(n_k27ac_gained, 0),
    n_k27ac_lost = replace_na(n_k27ac_lost, 0),
    net_k27ac = replace_na(net_k27ac, 0),
    k27ac_status = replace_na(k27ac_status, "No Change")
  ) %>%
  # MeCP2
  left_join(mecp2_gene %>% dplyr::rename(gene = SYMBOL), by = "gene") %>%
  mutate(
    n_mecp2_up = replace_na(n_mecp2_up, 0),
    n_mecp2_down = replace_na(n_mecp2_down, 0),
    net_mecp2 = replace_na(net_mecp2, 0),
    mecp2_status = replace_na(mecp2_status, "No MeCP2")
  )

# Split into disc vs coord for analysis
disc <- master %>% dplyr::filter(quadrant == "Q2: mC dn/hmC up")
coord <- master %>% dplyr::filter(quadrant == "Q4: mC up/hmC dn")

cat(sprintf("  Discordant (Q2): %d genes\n", nrow(disc)))
cat(sprintf("  Coordinated (Q4): %d genes\n", nrow(coord)))

# =============================================================================
# ANALYSIS 1: METHYLATION EFFECT SIZE
# =============================================================================

cat("\n--- Analysis 1: Methylation effect size ---\n")

wt_mc <- wilcox.test(abs(disc$mc_diff), abs(coord$mc_diff))
wt_hmc <- wilcox.test(abs(disc$hmc_diff), abs(coord$hmc_diff))

cat(sprintf("  |mc_diff|  - Disc median: %.4f, Coord median: %.4f, %s\n",
            median(abs(disc$mc_diff)), median(abs(coord$mc_diff)), fmt_p(wt_mc$p.value)))
cat(sprintf("  |hmc_diff| - Disc median: %.4f, Coord median: %.4f, %s\n",
            median(abs(disc$hmc_diff)), median(abs(coord$hmc_diff)), fmt_p(wt_hmc$p.value)))

p1_mc <- ggplot(bind_rows(
    disc %>% mutate(group = "Discordant (Q2)"),
    coord %>% mutate(group = "Coordinated (Q4)")
  ), aes(x = group, y = abs(mc_diff) * 100, fill = group)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.08, size = 0.4) +
  scale_fill_manual(values = CONCORDANCE_COLORS, guide = "none") +
  labs(title = "|mC diff|", subtitle = fmt_p(wt_mc$p.value),
       x = NULL, y = "|mC difference| (%)") +
  theme_biomodal(base_size = 10)

p1_hmc <- ggplot(bind_rows(
    disc %>% mutate(group = "Discordant (Q2)"),
    coord %>% mutate(group = "Coordinated (Q4)")
  ), aes(x = group, y = abs(hmc_diff) * 100, fill = group)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.08, size = 0.4) +
  scale_fill_manual(values = CONCORDANCE_COLORS, guide = "none") +
  labs(title = "|hmC diff|", subtitle = fmt_p(wt_hmc$p.value),
       x = NULL, y = "|hmC difference| (%)") +
  theme_biomodal(base_size = 10)

# =============================================================================
# ANALYSIS 2: RNA-seq EXPRESSION
# =============================================================================

cat("\n--- Analysis 2: RNA-seq expression ---\n")

disc_lfc <- disc$log2FC[!is.na(disc$log2FC)]
coord_lfc <- coord$log2FC[!is.na(coord$log2FC)]

wt_lfc <- wilcox.test(disc_lfc, coord_lfc)
cat(sprintf("  log2FC - Disc median: %.4f, Coord median: %.4f, %s\n",
            median(disc_lfc), median(coord_lfc), fmt_p(wt_lfc$p.value)))

# Fisher's: fraction up vs down regulated
disc_up <- sum(disc$padj < 0.05 & disc$log2FC > 0, na.rm = TRUE)
disc_dn <- sum(disc$padj < 0.05 & disc$log2FC < 0, na.rm = TRUE)
coord_up <- sum(coord$padj < 0.05 & coord$log2FC > 0, na.rm = TRUE)
coord_dn <- sum(coord$padj < 0.05 & coord$log2FC < 0, na.rm = TRUE)

fisher_expr <- tryCatch(
  fisher.test(matrix(c(disc_up, disc_dn, coord_up, coord_dn), nrow = 2, byrow = TRUE)),
  error = function(e) list(estimate = NA, p.value = NA)
)
cat(sprintf("  Fisher's (Up/Down): OR = %.2f, %s\n",
            ifelse(is.na(fisher_expr$estimate), NA, fisher_expr$estimate),
            fmt_p(fisher_expr$p.value)))

p2_lfc <- ggplot(bind_rows(
    disc %>% dplyr::filter(!is.na(log2FC)) %>% mutate(group = "Discordant (Q2)"),
    coord %>% dplyr::filter(!is.na(log2FC)) %>% mutate(group = "Coordinated (Q4)")
  ), aes(x = group, y = log2FC, fill = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.08, size = 0.4) +
  scale_fill_manual(values = CONCORDANCE_COLORS, guide = "none") +
  labs(title = "RNA-seq log2FC", subtitle = fmt_p(wt_lfc$p.value),
       x = NULL, y = "log2FC") +
  theme_biomodal(base_size = 10)

# =============================================================================
# ANALYSIS 3: ATAC-seq ACCESSIBILITY
# =============================================================================

cat("\n--- Analysis 3: ATAC-seq accessibility ---\n")

wt_atac <- wilcox.test(disc$net_atac, coord$net_atac)
cat(sprintf("  net_atac - Disc median: %.1f, Coord median: %.1f, %s\n",
            median(disc$net_atac), median(coord$net_atac), fmt_p(wt_atac$p.value)))

# Fisher's: any ATAC up vs any ATAC down
disc_atac_up <- sum(disc$n_atac_up > 0)
disc_atac_dn <- sum(disc$n_atac_down > 0)
coord_atac_up <- sum(coord$n_atac_up > 0)
coord_atac_dn <- sum(coord$n_atac_down > 0)

fisher_atac <- tryCatch(
  fisher.test(matrix(c(disc_atac_up, disc_atac_dn, coord_atac_up, coord_atac_dn),
                     nrow = 2, byrow = TRUE)),
  error = function(e) list(estimate = NA, p.value = NA)
)
cat(sprintf("  Fisher's (ATAC Up/Down): OR = %.2f, %s\n",
            ifelse(is.na(fisher_atac$estimate), NA, fisher_atac$estimate),
            fmt_p(fisher_atac$p.value)))

p3_atac <- ggplot(bind_rows(
    disc %>% mutate(group = "Discordant (Q2)"),
    coord %>% mutate(group = "Coordinated (Q4)")
  ), aes(x = group, y = net_atac, fill = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.08, size = 0.4) +
  scale_fill_manual(values = CONCORDANCE_COLORS, guide = "none") +
  labs(title = "Net ATAC Change", subtitle = fmt_p(wt_atac$p.value),
       x = NULL, y = "net ATAC (up - down)") +
  theme_biomodal(base_size = 10)

# =============================================================================
# ANALYSIS 4: H2AK119ub SIGNAL
# =============================================================================

cat("\n--- Analysis 4: H2AK119ub signal ---\n")

disc_k119 <- disc$k119ub_gb_log2fc[!is.na(disc$k119ub_gb_log2fc)]
coord_k119 <- coord$k119ub_gb_log2fc[!is.na(coord$k119ub_gb_log2fc)]

wt_k119 <- wilcox.test(disc_k119, coord_k119)
cat(sprintf("  K119ub gb_log2fc - Disc median: %.4f (n=%d), Coord median: %.4f (n=%d), %s\n",
            median(disc_k119), length(disc_k119),
            median(coord_k119), length(coord_k119), fmt_p(wt_k119$p.value)))

p4_k119 <- ggplot(bind_rows(
    disc %>% dplyr::filter(!is.na(k119ub_gb_log2fc)) %>% mutate(group = "Discordant (Q2)"),
    coord %>% dplyr::filter(!is.na(k119ub_gb_log2fc)) %>% mutate(group = "Coordinated (Q4)")
  ), aes(x = group, y = k119ub_gb_log2fc, fill = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.08, size = 0.4) +
  scale_fill_manual(values = CONCORDANCE_COLORS, guide = "none") +
  labs(title = "K119ub Gene Body Signal", subtitle = fmt_p(wt_k119$p.value),
       x = NULL, y = "K119ub log2FC (mut/ctrl)") +
  theme_biomodal(base_size = 10)

# =============================================================================
# ANALYSIS 5: H3K27ac STATUS
# =============================================================================

cat("\n--- Analysis 5: H3K27ac status ---\n")

# Construct table: group x gained/lost
disc_k27ac_gained <- sum(disc$n_k27ac_gained > 0)
disc_k27ac_lost <- sum(disc$n_k27ac_lost > 0)
coord_k27ac_gained <- sum(coord$n_k27ac_gained > 0)
coord_k27ac_lost <- sum(coord$n_k27ac_lost > 0)

fisher_k27ac <- tryCatch(
  fisher.test(matrix(c(disc_k27ac_gained, disc_k27ac_lost,
                       coord_k27ac_gained, coord_k27ac_lost),
                     nrow = 2, byrow = TRUE)),
  error = function(e) list(estimate = NA, p.value = NA)
)
cat(sprintf("  Disc: %d gained, %d lost | Coord: %d gained, %d lost\n",
            disc_k27ac_gained, disc_k27ac_lost, coord_k27ac_gained, coord_k27ac_lost))
cat(sprintf("  Fisher's (Gained/Lost): %s\n", fmt_p(fisher_k27ac$p.value)))

k27ac_bar_data <- data.frame(
  group = rep(c("Discordant (Q2)", "Coordinated (Q4)"), each = 2),
  status = rep(c("Gained", "Lost"), 2),
  count = c(disc_k27ac_gained, disc_k27ac_lost, coord_k27ac_gained, coord_k27ac_lost),
  total = rep(c(nrow(disc), nrow(coord)), each = 2)
) %>%
  mutate(pct = 100 * count / total)

p5_k27ac <- ggplot(k27ac_bar_data, aes(x = group, y = pct, fill = status)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, color = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%d", count)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 2.8) +
  scale_fill_manual(values = c("Gained" = "#FF7F00", "Lost" = "#1F78B4"),
                    name = "H3K27ac") +
  labs(title = "H3K27ac Gained/Lost", subtitle = fmt_p(fisher_k27ac$p.value),
       x = NULL, y = "% of genes") +
  theme_biomodal(base_size = 10) +
  theme(legend.position = "bottom", legend.key.size = unit(0.4, "cm"))

# =============================================================================
# ANALYSIS 6: CHROMATIN STATE
# =============================================================================

cat("\n--- Analysis 6: Chromatin state ---\n")

chrom_disc <- table(master$chromatin_state[master$quadrant == "Q2: mC dn/hmC up"])
chrom_coord <- table(master$chromatin_state[master$quadrant == "Q4: mC up/hmC dn"])

# Build combined table for Fisher's test
chrom_tab <- rbind(
  Discordant = chrom_disc[CHROMATIN_STATE_ORDER],
  Coordinated = chrom_coord[CHROMATIN_STATE_ORDER]
)
chrom_tab[is.na(chrom_tab)] <- 0

# Remove columns with all zeros
chrom_tab <- chrom_tab[, colSums(chrom_tab) > 0, drop = FALSE]

fisher_chrom <- tryCatch(
  fisher.test(chrom_tab, simulate.p.value = TRUE, B = 10000),
  error = function(e) list(p.value = NA)
)
cat(sprintf("  Fisher's (chromatin state overall): %s\n", fmt_p(fisher_chrom$p.value)))

chrom_bar_data <- bind_rows(
  master %>% dplyr::filter(quadrant == "Q2: mC dn/hmC up") %>%
    count(chromatin_state) %>% mutate(group = "Discordant (Q2)",
                                      total = nrow(disc)),
  master %>% dplyr::filter(quadrant == "Q4: mC up/hmC dn") %>%
    count(chromatin_state) %>% mutate(group = "Coordinated (Q4)",
                                      total = nrow(coord))
) %>%
  mutate(pct = 100 * n / total)

p6_chrom <- ggplot(chrom_bar_data, aes(x = group, y = pct, fill = chromatin_state)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = CHROMATIN_STATE_COLORS, name = "Chromatin State") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(title = "Chromatin State", subtitle = fmt_p(fisher_chrom$p.value),
       x = NULL, y = "% of genes") +
  theme_biomodal(base_size = 10) +
  theme(legend.position = "bottom", legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 7))

# =============================================================================
# ANALYSIS 7: MeCP2 BINDING
# =============================================================================

cat("\n--- Analysis 7: MeCP2 binding ---\n")

disc_mecp2_up <- sum(disc$n_mecp2_up > 0)
disc_mecp2_dn <- sum(disc$n_mecp2_down > 0)
coord_mecp2_up <- sum(coord$n_mecp2_up > 0)
coord_mecp2_dn <- sum(coord$n_mecp2_down > 0)

fisher_mecp2 <- tryCatch(
  fisher.test(matrix(c(disc_mecp2_up, disc_mecp2_dn,
                       coord_mecp2_up, coord_mecp2_dn),
                     nrow = 2, byrow = TRUE)),
  error = function(e) list(estimate = NA, p.value = NA)
)
cat(sprintf("  Disc: %d MeCP2 Up, %d MeCP2 Down | Coord: %d MeCP2 Up, %d MeCP2 Down\n",
            disc_mecp2_up, disc_mecp2_dn, coord_mecp2_up, coord_mecp2_dn))
cat(sprintf("  Fisher's (MeCP2 Up/Down): %s\n", fmt_p(fisher_mecp2$p.value)))

mecp2_bar_data <- data.frame(
  group = rep(c("Discordant (Q2)", "Coordinated (Q4)"), each = 2),
  status = rep(c("MeCP2 Up", "MeCP2 Down"), 2),
  count = c(disc_mecp2_up, disc_mecp2_dn, coord_mecp2_up, coord_mecp2_dn),
  total = rep(c(nrow(disc), nrow(coord)), each = 2)
) %>%
  mutate(pct = 100 * count / total)

p7_mecp2 <- ggplot(mecp2_bar_data, aes(x = group, y = pct, fill = status)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6, color = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%d", count)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 2.8) +
  scale_fill_manual(values = c("MeCP2 Up" = "#D95F02", "MeCP2 Down" = "#7570B3"),
                    name = "MeCP2") +
  labs(title = "MeCP2 Binding", subtitle = fmt_p(fisher_mecp2$p.value),
       x = NULL, y = "% of genes") +
  theme_biomodal(base_size = 10) +
  theme(legend.position = "bottom", legend.key.size = unit(0.4, "cm"))

# =============================================================================
# ANALYSIS 8: HI-C LOOP INVOLVEMENT
# =============================================================================

cat("\n--- Analysis 8: Hi-C loop involvement ---\n")

disc_loop <- sum(disc$loop_involved)
disc_no_loop <- nrow(disc) - disc_loop
coord_loop <- sum(coord$loop_involved)
coord_no_loop <- nrow(coord) - coord_loop

fisher_loop <- tryCatch(
  fisher.test(matrix(c(disc_loop, disc_no_loop, coord_loop, coord_no_loop),
                     nrow = 2, byrow = TRUE)),
  error = function(e) list(estimate = NA, p.value = NA)
)
cat(sprintf("  Disc: %d (%.1f%%) at loop anchors | Coord: %d (%.1f%%) at loop anchors\n",
            disc_loop, 100 * disc_loop / nrow(disc),
            coord_loop, 100 * coord_loop / nrow(coord)))
cat(sprintf("  Fisher's (Loop/No Loop): %s\n", fmt_p(fisher_loop$p.value)))

loop_bar_data <- data.frame(
  group = rep(c("Discordant (Q2)", "Coordinated (Q4)"), each = 2),
  status = rep(c("At Loop Anchor", "Not at Loop"), 2),
  count = c(disc_loop, disc_no_loop, coord_loop, coord_no_loop),
  total = rep(c(nrow(disc), nrow(coord)), each = 2)
) %>%
  mutate(pct = 100 * count / total)

p8_loop <- ggplot(loop_bar_data, aes(x = group, y = pct, fill = status)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7,
           color = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.0f%%\n(n=%d)", pct, count)),
            position = position_stack(vjust = 0.5), size = 2.8) +
  scale_fill_manual(values = c("At Loop Anchor" = "#984EA3", "Not at Loop" = "grey80"),
                    name = "") +
  labs(title = "Loop Involvement", subtitle = fmt_p(fisher_loop$p.value),
       x = NULL, y = "% of genes") +
  theme_biomodal(base_size = 10) +
  theme(legend.position = "bottom", legend.key.size = unit(0.4, "cm"))

# =============================================================================
# ANALYSIS 9: GO ENRICHMENT
# =============================================================================

cat("\n--- Analysis 9: GO enrichment ---\n")

# Convert discordant gene symbols to Entrez IDs
disc_entrez <- tryCatch(
  bitr(disc$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db),
  error = function(e) data.frame(SYMBOL = character(), ENTREZID = character())
)

coord_entrez <- tryCatch(
  bitr(coord$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db),
  error = function(e) data.frame(SYMBOL = character(), ENTREZID = character())
)

cat(sprintf("  Discordant genes mapped to Entrez: %d / %d\n",
            nrow(disc_entrez), nrow(disc)))
cat(sprintf("  Coordinated genes mapped to Entrez: %d / %d\n",
            nrow(coord_entrez), nrow(coord)))

ego <- NULL
if (nrow(disc_entrez) >= 5 && nrow(coord_entrez) >= 5) {
  ego <- tryCatch(
    enrichGO(
      gene = disc_entrez$ENTREZID,
      universe = coord_entrez$ENTREZID,
      OrgDb = org.Mm.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.2,
      readable = TRUE
    ),
    error = function(e) {
      cat(sprintf("  GO enrichment failed: %s\n", e$message))
      NULL
    }
  )

  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    cat(sprintf("  GO terms enriched: %d\n", nrow(as.data.frame(ego))))
  } else {
    cat("  No significant GO enrichment found.\n")
    ego <- NULL
  }
} else {
  cat("  Too few genes for GO enrichment analysis.\n")
}

# =============================================================================
# FIGURE 21a: COMPOSITE 3x3 GRID
# =============================================================================

cat("\n--- Assembling Figure 21a: Composite panel ---\n")

# Create a placeholder for position 9 if no GO enrichment
p9_placeholder <- ggplot() +
  annotate("text", x = 0.5, y = 0.5,
           label = ifelse(is.null(ego), "GO: No enrichment\nfound", "See 21d"),
           size = 4, color = "grey50") +
  theme_void() +
  labs(title = "GO Enrichment")

p_21a <- (p1_mc | p1_hmc | p2_lfc) /
         (p3_atac | p4_k119 | p5_k27ac) /
         (p6_chrom | p7_mecp2 | p8_loop) +
  plot_annotation(
    title = "Discordant Gene Characterization: Q2 (mC dn/hmC up) vs Q4 (mC up/hmC dn)",
    subtitle = sprintf("Discordant: %d genes | Coordinated: %d genes",
                       nrow(disc), nrow(coord)),
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 11, color = "grey40"))
  )

save_multiformat_ggplot(p_21a, file.path(OUTPUT_DIR, "21a_discordant_composite"),
                        width = 16, height = 16)

# Save individual panels
save_multiformat_ggplot(p1_mc, file.path(OUTPUT_DIR, "21a_panel_mc_diff"), width = 5, height = 5)
save_multiformat_ggplot(p1_hmc, file.path(OUTPUT_DIR, "21a_panel_hmc_diff"), width = 5, height = 5)
save_multiformat_ggplot(p2_lfc, file.path(OUTPUT_DIR, "21a_panel_log2fc"), width = 5, height = 5)
save_multiformat_ggplot(p3_atac, file.path(OUTPUT_DIR, "21a_panel_atac"), width = 5, height = 5)
save_multiformat_ggplot(p4_k119, file.path(OUTPUT_DIR, "21a_panel_k119ub"), width = 5, height = 5)
save_multiformat_ggplot(p5_k27ac, file.path(OUTPUT_DIR, "21a_panel_k27ac"), width = 6, height = 5)
save_multiformat_ggplot(p6_chrom, file.path(OUTPUT_DIR, "21a_panel_chromatin"), width = 6, height = 5)
save_multiformat_ggplot(p7_mecp2, file.path(OUTPUT_DIR, "21a_panel_mecp2"), width = 6, height = 5)
save_multiformat_ggplot(p8_loop, file.path(OUTPUT_DIR, "21a_panel_loop"), width = 6, height = 5)

# =============================================================================
# FIGURE 21b: mc_diff vs hmc_diff SCATTER COLORED BY CONCORDANCE
# =============================================================================

cat("\n--- Creating Figure 21b: mC vs hmC scatter ---\n")

# All significant genes as background
all_bg <- inner_join(
  mc_dmr %>% dplyr::select(gene, mc_diff = mod_difference),
  hmc_dmr %>% dplyr::select(gene, hmc_diff = mod_difference),
  by = "gene"
)

scatter_21b_data <- bind_rows(
  all_bg %>%
    dplyr::filter(!(gene %in% coordinated_all$gene)) %>%
    mutate(group = "Background"),
  coordinated_all %>%
    dplyr::filter(quadrant == "Q4: mC up/hmC dn") %>%
    dplyr::select(gene, mc_diff, hmc_diff) %>%
    mutate(group = "Coordinated (Q4)"),
  coordinated_all %>%
    dplyr::filter(quadrant == "Q2: mC dn/hmC up") %>%
    dplyr::select(gene, mc_diff, hmc_diff) %>%
    mutate(group = "Discordant (Q2)")
) %>%
  mutate(group = factor(group, levels = c("Background", "Coordinated (Q4)", "Discordant (Q2)")))

p_21b <- ggplot(scatter_21b_data, aes(x = mc_diff * 100, y = hmc_diff * 100)) +
  geom_point(data = scatter_21b_data %>% dplyr::filter(group == "Background"),
             color = "grey80", alpha = 0.3, size = 1) +
  geom_point(data = scatter_21b_data %>% dplyr::filter(group == "Coordinated (Q4)"),
             aes(color = group), alpha = 0.5, size = 1.5) +
  geom_point(data = scatter_21b_data %>% dplyr::filter(group == "Discordant (Q2)"),
             aes(color = group), alpha = 0.7, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = c("Coordinated (Q4)" = "#4DAF4A",
                                 "Discordant (Q2)" = "#E41A1C"),
                     name = "Group") +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           alpha = 0.05, fill = "#E41A1C") +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
           alpha = 0.05, fill = "#4DAF4A") +
  annotate("text", x = -8, y = 8, label = "DISCORDANT\nmC dn / hmC up",
           color = "#E41A1C", fontface = "bold", size = 3.5) +
  annotate("text", x = 12, y = -10, label = "COORDINATED\nmC up / hmC dn",
           color = "#4DAF4A", fontface = "bold", size = 3.5) +
  labs(
    title = "Coordinated vs Discordant mC/hmC Changes",
    subtitle = sprintf("All genes significant in both mC and hmC (n = %d) | Q4: %d, Q2: %d",
                       nrow(coordinated_all), nrow(coord), nrow(disc)),
    x = "5mC Change (%)", y = "5hmC Change (%)"
  ) +
  theme_biomodal() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))

save_multiformat_ggplot(p_21b, file.path(OUTPUT_DIR, "21b_mc_hmc_concordance_scatter"),
                        width = 10, height = 9)

# =============================================================================
# FIGURE 21c: mc_diff vs log2FC per group
# =============================================================================

cat("\n--- Creating Figure 21c: mC vs expression per group ---\n")

scatter_21c_data <- master %>%
  dplyr::filter(quadrant %in% c("Q2: mC dn/hmC up", "Q4: mC up/hmC dn") &
                !is.na(log2FC)) %>%
  mutate(group = ifelse(quadrant == "Q2: mC dn/hmC up",
                        "Discordant (Q2)", "Coordinated (Q4)"))

# Per-group correlations
cor_disc <- cor.test(
  scatter_21c_data$mc_diff[scatter_21c_data$group == "Discordant (Q2)"] * 100,
  scatter_21c_data$log2FC[scatter_21c_data$group == "Discordant (Q2)"],
  method = "spearman"
)
cor_coord <- cor.test(
  scatter_21c_data$mc_diff[scatter_21c_data$group == "Coordinated (Q4)"] * 100,
  scatter_21c_data$log2FC[scatter_21c_data$group == "Coordinated (Q4)"],
  method = "spearman"
)

cat(sprintf("  Spearman (Discordant): rho = %.3f, p = %.2e\n",
            cor_disc$estimate, cor_disc$p.value))
cat(sprintf("  Spearman (Coordinated): rho = %.3f, p = %.2e\n",
            cor_coord$estimate, cor_coord$p.value))

p_21c <- ggplot(scatter_21c_data, aes(x = mc_diff * 100, y = log2FC, color = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 0.8) +
  facet_wrap(~group, scales = "free_x") +
  scale_color_manual(values = CONCORDANCE_COLORS, guide = "none") +
  labs(
    title = "5mC Change vs Gene Expression by Concordance Group",
    subtitle = sprintf("Discordant: rho=%.3f (%s) | Coordinated: rho=%.3f (%s)",
                       cor_disc$estimate, fmt_p(cor_disc$p.value),
                       cor_coord$estimate, fmt_p(cor_coord$p.value)),
    x = "5mC Change (%)", y = "RNA-seq log2FC"
  ) +
  theme_biomodal() +
  theme(strip.text = element_text(size = 12, face = "bold"))

save_multiformat_ggplot(p_21c, file.path(OUTPUT_DIR, "21c_mc_vs_expression_per_group"),
                        width = 12, height = 7)

# =============================================================================
# FIGURE 21d: GO ENRICHMENT DOT PLOT
# =============================================================================

cat("\n--- Creating Figure 21d: GO enrichment ---\n")

if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  p_21d <- dotplot(ego, showCategory = 15) +
    labs(title = "GO Enrichment: Discordant (Q2) Genes",
         subtitle = sprintf("Against coordinated (Q4) background | %d terms",
                            nrow(as.data.frame(ego)))) +
    theme_biomodal()

  save_multiformat_ggplot(p_21d, file.path(OUTPUT_DIR, "21d_go_enrichment_discordant"),
                          width = 10, height = 9)
  cat("  GO dot plot saved.\n")
} else {
  cat("  No GO enrichment to plot.\n")
}

# =============================================================================
# FIGURE 21e: ALL 4 QUADRANTS COMPREHENSIVE PANEL
# =============================================================================

cat("\n--- Creating Figure 21e: All 4 quadrants ---\n")

# For the 4-quadrant figure, compute key dimensions for all quadrants
quad_data <- master %>%
  dplyr::filter(!is.na(log2FC)) %>%
  mutate(quadrant = factor(quadrant, levels = names(QUADRANT_COLORS)))

# Panel e1: log2FC by quadrant
pe1 <- ggplot(quad_data, aes(x = quadrant, y = log2FC, fill = quadrant)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  scale_fill_manual(values = QUADRANT_COLORS, guide = "none") +
  labs(title = "RNA-seq log2FC", x = NULL, y = "log2FC") +
  theme_biomodal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

# Panel e2: net ATAC by quadrant
pe2 <- ggplot(master %>% mutate(quadrant = factor(quadrant, levels = names(QUADRANT_COLORS))),
              aes(x = quadrant, y = net_atac, fill = quadrant)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  scale_fill_manual(values = QUADRANT_COLORS, guide = "none") +
  labs(title = "Net ATAC Change", x = NULL, y = "net ATAC") +
  theme_biomodal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

# Panel e3: K119ub by quadrant
pe3 <- ggplot(master %>%
                dplyr::filter(!is.na(k119ub_gb_log2fc)) %>%
                mutate(quadrant = factor(quadrant, levels = names(QUADRANT_COLORS))),
              aes(x = quadrant, y = k119ub_gb_log2fc, fill = quadrant)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  scale_fill_manual(values = QUADRANT_COLORS, guide = "none") +
  labs(title = "K119ub Gene Body Signal", x = NULL, y = "K119ub log2FC") +
  theme_biomodal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

# Panel e4: Chromatin state by quadrant
chrom_quad_data <- master %>%
  mutate(quadrant = factor(quadrant, levels = names(QUADRANT_COLORS))) %>%
  group_by(quadrant) %>%
  mutate(total = n()) %>%
  ungroup() %>%
  count(quadrant, chromatin_state, total) %>%
  mutate(pct = 100 * n / total)

pe4 <- ggplot(chrom_quad_data, aes(x = quadrant, y = pct, fill = chromatin_state)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = CHROMATIN_STATE_COLORS, name = "Chromatin") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(title = "Chromatin State", x = NULL, y = "% of genes") +
  theme_biomodal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        legend.position = "bottom", legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 6))

p_21e <- (pe1 | pe2) / (pe3 | pe4) +
  plot_annotation(
    title = "All 4 Quadrants: Key Epigenomic Dimensions",
    subtitle = sprintf("Q1: %d | Q2 (disc): %d | Q3: %d | Q4 (coord): %d",
                       sum(coordinated_all$quadrant == "Q1: mC dn/hmC dn"),
                       nrow(disc),
                       sum(coordinated_all$quadrant == "Q3: mC up/hmC up"),
                       nrow(coord)),
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 11, color = "grey40"))
  )

save_multiformat_ggplot(p_21e, file.path(OUTPUT_DIR, "21e_all_quadrants_comprehensive"),
                        width = 14, height = 12)

# =============================================================================
# OUTPUT TABLE
# =============================================================================

cat("\nSaving output tables...\n")

export_df <- master %>%
  dplyr::select(
    gene, mc_diff, hmc_diff, combined_effect, quadrant,
    log2FC, padj, baseMean,
    net_atac, n_atac_up, n_atac_down,
    k119ub_gb_log2fc, k119ub_pr_log2fc, k119ub_gb_class,
    n_k27ac_gained, n_k27ac_lost, k27ac_status,
    chromatin_state,
    n_mecp2_up, n_mecp2_down, mecp2_status,
    loop_involved
  ) %>%
  arrange(quadrant, desc(combined_effect))

write.table(export_df, file.path(TABLES_DIR, "discordant_gene_characteristics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: discordant_gene_characteristics.tsv (%d rows)\n", nrow(export_df)))

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("SECTION 21 COMPLETE -- DISCORDANT GENE CHARACTERIZATION SUMMARY\n")
cat("================================================================================\n\n")

cat(sprintf("Genes significant in both mC and hmC: %d\n", nrow(coordinated_all)))
cat(sprintf("  Q1 (mC dn/hmC dn): %d\n", sum(coordinated_all$quadrant == "Q1: mC dn/hmC dn")))
cat(sprintf("  Q2 (mC dn/hmC up, DISCORDANT): %d\n", nrow(disc)))
cat(sprintf("  Q3 (mC up/hmC up): %d\n", sum(coordinated_all$quadrant == "Q3: mC up/hmC up")))
cat(sprintf("  Q4 (mC up/hmC dn, COORDINATED): %d\n", nrow(coord)))

cat("\nKey comparisons (Discordant Q2 vs Coordinated Q4):\n")
cat(sprintf("  |mC diff|:       Wilcoxon %s\n", fmt_p(wt_mc$p.value)))
cat(sprintf("  |hmC diff|:      Wilcoxon %s\n", fmt_p(wt_hmc$p.value)))
cat(sprintf("  log2FC:          Wilcoxon %s\n", fmt_p(wt_lfc$p.value)))
cat(sprintf("  Expression Up/Dn: Fisher's %s\n", fmt_p(fisher_expr$p.value)))
cat(sprintf("  net ATAC:        Wilcoxon %s\n", fmt_p(wt_atac$p.value)))
cat(sprintf("  ATAC Up/Down:    Fisher's %s\n", fmt_p(fisher_atac$p.value)))
cat(sprintf("  K119ub gb sig:   Wilcoxon %s\n", fmt_p(wt_k119$p.value)))
cat(sprintf("  H3K27ac G/L:     Fisher's %s\n", fmt_p(fisher_k27ac$p.value)))
cat(sprintf("  Chromatin state: Fisher's %s\n", fmt_p(fisher_chrom$p.value)))
cat(sprintf("  MeCP2 Up/Down:   Fisher's %s\n", fmt_p(fisher_mecp2$p.value)))
cat(sprintf("  Loop involved:   Fisher's %s\n", fmt_p(fisher_loop$p.value)))

if (!is.null(ego)) {
  cat(sprintf("  GO terms enriched: %d (padj < 0.1)\n", nrow(as.data.frame(ego))))
}

cat(sprintf("\nOutputs:\n"))
cat(sprintf("  Figures: %s/21*\n", OUTPUT_DIR))
cat(sprintf("  Table: %s/discordant_gene_characteristics.tsv\n", TABLES_DIR))
cat("================================================================================\n")
