# biomodal/downstream/scripts/viz_sections/section_37_permutation_gene_level.R
# Section 37: Gene-Level Label-Shuffle Permutation Tests
# Standalone script - sources shared config for all dependencies and data
#
# Validates gene-level Fisher's exact test / O/E enrichment results from
# sections 12e, 15a-c, 19g, 20d, 27c, 29 Step 3, 31b, and 33c using
# chromosome-stratified label-shuffle permutation tests.
#
# Why not regioneReloaded? These analyses cross-tabulate gene-level binary
# labels (e.g., "mC Up" vs "mC Down" x "ATAC Up" vs "ATAC Down"). The unit
# is a gene symbol with two categorical assignments, not a genomic interval.
# regioneReloaded tests spatial co-occurrence of region sets, which is a
# different question answered in sections 34-36.
#
# Methodology: For each 2x2 gene-level contingency table, shuffle the column
# labels within chromosomes (preserving per-chromosome marginal counts and
# genomic neighborhood structure). Compute Fisher OR per permutation to build
# an empirical null distribution. Compare observed OR to null for an honest
# p-value that accounts for genomic non-independence.
#
# Figures:
#   37a: Forest plot of all permutation z-scores (color by source section)
#   37b: Null distribution histograms for 4 strongest effects (2x2 panel)
#   37c: Fisher-vs-permutation comparison table
#   37d: Scatter of observed log2(OR) vs permutation z-score
#
# Run from downstream/ directory:
#   Rscript scripts/viz_sections/section_37_permutation_gene_level.R

# =============================================================================
# BLOCK 0: SETUP
# =============================================================================

source("scripts/viz_sections/_shared_config.R")

# Skip if cached results exist (avoids re-run from run_all_sections.sh glob)
# FORCE_RERUN env var (set by run_permutation.sb) bypasses this check
.cache_check <- file.path(TABLES_DIR, "permutation_37_gene_level.rds")
if (file.exists(.cache_check) && Sys.getenv("FORCE_RERUN", unset = "") == "") {
  cat("Section 37: Cached RDS found at", .cache_check, "-- skipping.\n")
  cat("  To re-run: use run_permutation.sb or export FORCE_RERUN=1\n")
  quit(save = "no", status = 0)
}
.force_rerun <- Sys.getenv("FORCE_RERUN", unset = "") != ""

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(readxl)
  library(gridExtra)
})

cat("\n")
cat("===========================================================================\n")
cat("  SECTION 37: GENE-LEVEL LABEL-SHUFFLE PERMUTATION TESTS\n")
cat("===========================================================================\n\n")

# Configuration: accept ntimes from command line (e.g., Rscript section_37_... 10000)
args <- commandArgs(trailingOnly = TRUE)
PERM_NTIMES <- if (length(args) >= 1 && !is.na(as.numeric(args[1]))) as.numeric(args[1]) else 10000
PERM_SEED   <- 42

SECTION_DIR <- file.path(OUTPUT_DIR, "37_permutation_gene_level")
dir.create(SECTION_DIR, recursive = TRUE, showWarnings = FALSE)

FDR_THRESHOLD <- 0.05
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# RNA-seq and compartment file paths
RNA_SEQ_FILE <- "/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx"
COMPARTMENT_FILE <- file.path(BASE_DIR, "../../tads/tad-pc-analysis/inputs/late/diffPC/diffcompartments.txt")

# Compartment shift thresholds (matching section 29)
SHIFT_FDR  <- 0.05
SHIFT_DIFF <- 0.30

# GREAT-style regulatory domain parameters (matching section 27)
GREAT_UPSTREAM      <- 5000
GREAT_DOWNSTREAM    <- 1000
GREAT_MAX_EXTENSION <- 100000

# Section 27 anchor gene export
ANCHOR_GENE_TABLE <- file.path(TABLES_DIR, "anchor_gene_associations_great.tsv")

# Loop file from section 27
LOOPS_FILE <- file.path(BASE_DIR, "../../outputs/250402-late_outputs/merged_loops/characterized_loops.tsv")

# Input validation
cat("--- Input Validation ---\n\n")

stopifnot("ATAC Up BED not found" = file.exists(ATAC_FILES$up))
stopifnot("ATAC Down BED not found" = file.exists(ATAC_FILES$down))
stopifnot("K119ub ctrl BED not found" = file.exists(K119UB_FILES$ctrl))
stopifnot("K119ub mut BED not found" = file.exists(K119UB_FILES$mut))
stopifnot("H3K27ac ctrl BED not found" = file.exists(H3K27AC_FILES$ctrl))
stopifnot("H3K27ac mut BED not found" = file.exists(H3K27AC_FILES$mut))
stopifnot("MeCP2 annotated not found" = file.exists(MECP2_FILES$annotated))
stopifnot("DiffBind ATAC not found" = file.exists(DIFFBIND_FILES$atac))
stopifnot("DiffBind K27ac not found" = file.exists(DIFFBIND_FILES$k27ac))
stopifnot("DiffBind K27me3 not found" = file.exists(DIFFBIND_FILES$k27me3))
stopifnot("DiffBind K119ub not found" = file.exists(DIFFBIND_FILES$k119ub))
stopifnot("RNA-seq Excel not found" = file.exists(RNA_SEQ_FILE))
stopifnot("HOMER compartment file not found" = file.exists(COMPARTMENT_FILE))
stopifnot("Loop file not found" = file.exists(LOOPS_FILE))
stopifnot("mc_dmr not loaded" = exists("mc_dmr") && nrow(mc_dmr) > 0)
stopifnot("hmc_dmr not loaded" = exists("hmc_dmr") && nrow(hmc_dmr) > 0)

cat("  All core input files validated.\n")

# Check anchor gene table (from section 27) -- will re-derive if missing
has_anchor_table <- file.exists(ANCHOR_GENE_TABLE)
if (has_anchor_table) {
  cat(sprintf("  [OK] Anchor gene table: %s\n", basename(ANCHOR_GENE_TABLE)))
} else {
  cat("  [--] Anchor gene table not found -- will re-derive inline for tests 37-09/37-11\n")
}

cat("\n")

# =============================================================================
# BLOCK 1: CORE FUNCTION -- permute_gene_2x2()
# =============================================================================

#' Chromosome-stratified permutation test for 2x2 gene-level contingency tables
#'
#' Shuffles col_label assignments within each chromosome, preserving:
#'   - Per-chromosome marginal counts
#'   - Total count of each label
#'   - Genomic neighborhood structure (genes on same chr stay together)
#'
#' @param gene_df data.frame with columns: gene, chr, row_label, col_label
#' @param ntimes Number of permutations (default 10000)
#' @return list with observed_or, null_distribution, empirical_p, z_score, ci_95
permute_gene_2x2 <- function(gene_df, ntimes = PERM_NTIMES) {
  stopifnot(all(c("gene", "chr", "row_label", "col_label") %in% colnames(gene_df)))
  stopifnot(nrow(gene_df) > 0)

  # Observed Fisher OR
  obs_table <- table(gene_df$row_label, gene_df$col_label)
  obs_fisher <- fisher.test(obs_table)
  obs_log2or <- log2(pmax(obs_fisher$estimate, 1e-10))

  # Pre-split by chromosome for faster shuffling
  chr_list <- split(seq_len(nrow(gene_df)), gene_df$chr)

  # Permutation null: shuffle col_label within chromosomes
  null_log2or <- numeric(ntimes)
  col_labels <- gene_df$col_label
  row_labels <- gene_df$row_label

  for (i in seq_len(ntimes)) {
    shuffled_col <- col_labels
    for (idx in chr_list) {
      shuffled_col[idx] <- sample(col_labels[idx])
    }
    perm_table <- table(row_labels, shuffled_col)
    perm_fisher <- fisher.test(perm_table)
    null_log2or[i] <- log2(pmax(perm_fisher$estimate, 1e-10))
  }

  # Replace any Inf/-Inf in null (from zero cells)
  null_log2or[is.infinite(null_log2or)] <- NA_real_
  null_log2or <- null_log2or[!is.na(null_log2or)]

  if (length(null_log2or) < 100) {
    warning("Fewer than 100 valid permutations -- results unreliable")
  }

  # Empirical p-value (two-sided, Phipson & Smyth correction)
  empirical_p <- (sum(abs(null_log2or) >= abs(obs_log2or)) + 1) / (length(null_log2or) + 1)

  # Z-score
  null_sd <- sd(null_log2or)
  if (is.na(null_sd) || null_sd == 0) {
    z_score <- NA_real_
    warning("Zero-variance null distribution -- z-score set to NA")
  } else {
    z_score <- (obs_log2or - mean(null_log2or)) / null_sd
  }

  # 95% CI from null distribution
  ci_95 <- quantile(null_log2or, c(0.025, 0.975))

  list(
    observed_or     = obs_fisher$estimate,
    observed_log2or = obs_log2or,
    fisher_p        = obs_fisher$p.value,
    null_distribution = null_log2or,
    empirical_p     = empirical_p,
    z_score         = z_score,
    ci_95           = ci_95
  )
}

# =============================================================================
# BLOCK 2: REUSABLE DATA-PREP HELPERS
# =============================================================================

#' Prepare DMR row labels (gene, chr, row_label)
#' @param dmr_df DMR data frame (mc_dmr or hmc_dmr)
#' @param mod_type "mc" or "hmc"
#' @return data.frame with gene, chr, row_label
prep_dmr_row <- function(dmr_df, mod_type) {
  prefix <- ifelse(mod_type == "mc", "mC", "hmC")
  dmr_df %>%
    dplyr::filter(significant) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    dplyr::transmute(
      gene,
      chr,
      row_label = ifelse(mod_difference > 0,
                         paste(prefix, "Up"),
                         paste(prefix, "Down"))
    )
}

#' Prepare column labels from binary up/down BED files (ATAC pattern)
#' @param up_path Path to "up" BED file
#' @param down_path Path to "down" BED file
#' @param up_label Label for genes with more up peaks
#' @param down_label Label for genes with more down peaks
#' @param txdb_obj TxDb for ChIPseeker annotation
#' @return data.frame with gene, col_label
prep_binary_bed_col <- function(up_path, down_path, up_label, down_label, txdb_obj) {
  # Load and annotate up peaks
  up_gr <- load_chip_peaks(up_path, up_label)
  up_anno <- suppressMessages(annotatePeak(up_gr, TxDb = txdb_obj,
                                            annoDb = "org.Mm.eg.db", level = "gene"))
  up_df <- as.data.frame(up_anno) %>%
    dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(n_up = dplyr::n(), .groups = "drop")

  # Load and annotate down peaks
  down_gr <- load_chip_peaks(down_path, down_label)
  down_anno <- suppressMessages(annotatePeak(down_gr, TxDb = txdb_obj,
                                              annoDb = "org.Mm.eg.db", level = "gene"))
  down_df <- as.data.frame(down_anno) %>%
    dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(n_down = dplyr::n(), .groups = "drop")

  # Merge and compute net direction
  merged <- dplyr::full_join(up_df, down_df, by = "SYMBOL") %>%
    dplyr::mutate(
      n_up = tidyr::replace_na(n_up, 0),
      n_down = tidyr::replace_na(n_down, 0),
      net = n_up - n_down
    ) %>%
    dplyr::filter(net != 0) %>%
    dplyr::transmute(
      gene = SYMBOL,
      col_label = ifelse(net > 0, up_label, down_label)
    )

  cat(sprintf("    %s: %d up-genes, %d down-genes\n",
              up_label, sum(merged$col_label == up_label),
              sum(merged$col_label == down_label)))
  merged
}

#' Prepare column labels from condition-specific BED files (K119ub/H3K27ac pattern)
#' @param ctrl_path Path to control BED file
#' @param mut_path Path to mutant BED file
#' @param gained_label Label for genes with net gain
#' @param lost_label Label for genes with net loss
#' @param txdb_obj TxDb for ChIPseeker annotation
#' @return data.frame with gene, col_label
prep_condition_col <- function(ctrl_path, mut_path, gained_label, lost_label, txdb_obj) {
  # Load and annotate ctrl peaks
  ctrl_gr <- load_chip_peaks(ctrl_path, "ctrl")
  ctrl_anno <- suppressMessages(annotatePeak(ctrl_gr, TxDb = txdb_obj,
                                              annoDb = "org.Mm.eg.db", level = "gene"))
  ctrl_df <- as.data.frame(ctrl_anno) %>%
    dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(n_ctrl = dplyr::n(), .groups = "drop")

  # Load and annotate mut peaks
  mut_gr <- load_chip_peaks(mut_path, "mut")
  mut_anno <- suppressMessages(annotatePeak(mut_gr, TxDb = txdb_obj,
                                             annoDb = "org.Mm.eg.db", level = "gene"))
  mut_df <- as.data.frame(mut_anno) %>%
    dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(n_mut = dplyr::n(), .groups = "drop")

  # Merge and compute net direction
  merged <- dplyr::full_join(ctrl_df, mut_df, by = "SYMBOL") %>%
    dplyr::mutate(
      n_ctrl = tidyr::replace_na(n_ctrl, 0),
      n_mut = tidyr::replace_na(n_mut, 0),
      net = n_mut - n_ctrl
    ) %>%
    dplyr::filter(net != 0) %>%
    dplyr::transmute(
      gene = SYMBOL,
      col_label = ifelse(net > 0, gained_label, lost_label)
    )

  cat(sprintf("    %s: %d gained-genes, %d lost-genes\n",
              gained_label, sum(merged$col_label == gained_label),
              sum(merged$col_label == lost_label)))
  merged
}

#' Prepare MeCP2 column labels (nearest-to-TSS peak per gene)
#' @return data.frame with gene, col_label
prep_mecp2_col <- function() {
  mecp2_raw <- read.table(MECP2_FILES$annotated, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE, fill = TRUE,
                           quote = "", comment.char = "")

  # Ensure numeric
  mecp2_raw$Fold <- as.numeric(mecp2_raw$Fold)
  mecp2_raw$FDR <- as.numeric(mecp2_raw$FDR)
  mecp2_raw$distanceToTSS <- as.numeric(mecp2_raw$distanceToTSS)

  mecp2_gene <- mecp2_raw %>%
    dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(
      nearest_fold = Fold[which.min(abs(distanceToTSS))],
      nearest_fdr  = FDR[which.min(abs(distanceToTSS))],
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(nearest_fold) & nearest_fold != 0) %>%
    dplyr::transmute(
      gene = SYMBOL,
      col_label = ifelse(nearest_fold > 0, "MeCP2 Up", "MeCP2 Down")
    )

  cat(sprintf("    MeCP2: %d up-genes, %d down-genes\n",
              sum(mecp2_gene$col_label == "MeCP2 Up"),
              sum(mecp2_gene$col_label == "MeCP2 Down")))
  mecp2_gene
}

#' Prepare DiffBind column labels (nearest-to-TSS peak per gene, FDR-filtered)
#' @param diffbind_path Path to DiffBind TSV
#' @param up_label Label for upregulated
#' @param down_label Label for downregulated
#' @param txdb_obj TxDb for annotation
#' @return data.frame with gene, col_label
prep_diffbind_col <- function(diffbind_path, up_label, down_label, txdb_obj) {
  df <- read.table(diffbind_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  stopifnot(all(c("Summit_Chr", "Summit_Start", "Summit_End", "Fold", "FDR") %in% colnames(df)))

  valid <- !is.na(df$Summit_Chr) & !is.na(df$Summit_Start) & !is.na(df$Summit_End)
  df <- df[valid, ]

  gr <- GRanges(
    seqnames = df$Summit_Chr,
    ranges = IRanges(start = df$Summit_Start, end = df$Summit_End),
    Fold = df$Fold,
    FDR = df$FDR
  )

  anno <- suppressMessages(annotatePeak(gr, TxDb = txdb_obj,
                                         annoDb = "org.Mm.eg.db", level = "gene"))
  anno_df <- as.data.frame(anno)

  gene_summary <- anno_df %>%
    dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(
      nearest_fold = Fold[which.min(abs(distanceToTSS))],
      nearest_fdr  = FDR[which.min(abs(distanceToTSS))],
      .groups = "drop"
    ) %>%
    dplyr::filter(nearest_fdr < FDR_THRESHOLD & nearest_fold != 0) %>%
    dplyr::transmute(
      gene = SYMBOL,
      col_label = ifelse(nearest_fold > 0, up_label, down_label)
    )

  cat(sprintf("    %s: %d up-genes, %d down-genes\n",
              up_label, sum(gene_summary$col_label == up_label),
              sum(gene_summary$col_label == down_label)))
  gene_summary
}

#' Prepare expression column labels from RNA-seq
#' @param rna_path Path to RNA-seq Excel file
#' @return data.frame with gene, col_label
prep_expression_col <- function(rna_path) {
  rna_all <- read_excel(rna_path) %>%
    dplyr::select(gene = ensembl_gene_id, log2FC = log2FoldChange, padj) %>%
    dplyr::filter(!is.na(gene) & gene != "" & !is.na(padj) & !is.na(log2FC)) %>%
    dplyr::filter(padj < 0.05 & abs(log2FC) > 0.3) %>%
    dplyr::transmute(
      gene,
      col_label = ifelse(log2FC > 0, "Expr Up", "Expr Down")
    )

  cat(sprintf("    Expression: %d up-genes, %d down-genes\n",
              sum(rna_all$col_label == "Expr Up"),
              sum(rna_all$col_label == "Expr Down")))
  rna_all
}

#' Assemble gene_df from row and column data frames
#' Row df provides gene, chr, row_label; col df provides gene, col_label
#' @return data.frame with gene, chr, row_label, col_label (complete cases only)
assemble_gene_df <- function(row_df, col_df) {
  df <- dplyr::inner_join(row_df, col_df, by = "gene")
  df <- df[complete.cases(df[, c("gene", "chr", "row_label", "col_label")]), ]
  df
}

# =============================================================================
# BLOCK 3: PRE-COMPUTE COLUMN/ROW DATA & ASSEMBLE TEST gene_dfs
# =============================================================================

cat("================================================================================\n")
cat("PRE-COMPUTING GENE-LEVEL DIRECTION ASSIGNMENTS\n")
cat("================================================================================\n\n")

# --- Row data (DMR direction) ---
cat("--- DMR row labels ---\n")
mc_row  <- prep_dmr_row(mc_dmr, "mc")
cat(sprintf("  mC row: %d genes (%d Up, %d Down)\n",
            nrow(mc_row), sum(mc_row$row_label == "mC Up"), sum(mc_row$row_label == "mC Down")))
hmc_row <- prep_dmr_row(hmc_dmr, "hmc")
cat(sprintf("  hmC row: %d genes (%d Up, %d Down)\n",
            nrow(hmc_row), sum(hmc_row$row_label == "hmC Up"), sum(hmc_row$row_label == "hmC Down")))

# --- Column data (mark direction) ---
cat("\n--- Mark column labels ---\n")

cat("  ATAC binary BED:\n")
atac_col <- prep_binary_bed_col(ATAC_FILES$up, ATAC_FILES$down,
                                 "ATAC Up", "ATAC Down", txdb)

cat("  MeCP2:\n")
mecp2_col <- prep_mecp2_col()

cat("  K119ub condition-specific:\n")
k119ub_col <- prep_condition_col(K119UB_FILES$ctrl, K119UB_FILES$mut,
                                  "K119ub Gained", "K119ub Lost", txdb)

cat("  H3K27ac condition-specific:\n")
h3k27ac_col <- prep_condition_col(H3K27AC_FILES$ctrl, H3K27AC_FILES$mut,
                                   "H3K27ac Gained", "H3K27ac Lost", txdb)

cat("  DiffBind ATAC:\n")
db_atac_col <- prep_diffbind_col(DIFFBIND_FILES$atac, "ATAC DiffBind Up", "ATAC DiffBind Down", txdb)

cat("  DiffBind K27ac:\n")
db_k27ac_col <- prep_diffbind_col(DIFFBIND_FILES$k27ac, "K27ac DiffBind Up", "K27ac DiffBind Down", txdb)

cat("  DiffBind K27me3:\n")
db_k27me3_col <- prep_diffbind_col(DIFFBIND_FILES$k27me3, "K27me3 DiffBind Up", "K27me3 DiffBind Down", txdb)

cat("  DiffBind K119ub:\n")
db_k119ub_col <- prep_diffbind_col(DIFFBIND_FILES$k119ub, "K119ub DiffBind Up", "K119ub DiffBind Down", txdb)

cat("  RNA-seq expression:\n")
expr_col <- prep_expression_col(RNA_SEQ_FILE)

# --- Assemble standard tests (DMR row x mark col) ---
cat("\n--- Assembling test gene_dfs ---\n")

test_defs <- list(
  list(id = "37-01", section = "12e", desc = "mC direction x ATAC direction",
       gene_df = assemble_gene_df(mc_row, atac_col)),
  list(id = "37-02", section = "15a", desc = "hmC direction x MeCP2 direction",
       gene_df = assemble_gene_df(hmc_row, mecp2_col)),
  list(id = "37-03", section = "15b", desc = "hmC direction x ATAC direction",
       gene_df = assemble_gene_df(hmc_row, atac_col)),
  list(id = "37-04", section = "15c", desc = "hmC direction x K119ub direction",
       gene_df = assemble_gene_df(hmc_row, k119ub_col)),
  list(id = "37-05", section = "19g", desc = "mC direction x H3K27ac direction",
       gene_df = assemble_gene_df(mc_row, h3k27ac_col)),
  list(id = "37-06", section = "19g", desc = "hmC direction x H3K27ac direction",
       gene_df = assemble_gene_df(hmc_row, h3k27ac_col)),
  list(id = "37-07a", section = "33c", desc = "mC direction x ATAC DiffBind",
       gene_df = assemble_gene_df(mc_row, db_atac_col)),
  list(id = "37-07b", section = "33c", desc = "mC direction x K27ac DiffBind",
       gene_df = assemble_gene_df(mc_row, db_k27ac_col)),
  list(id = "37-07c", section = "33c", desc = "mC direction x K27me3 DiffBind",
       gene_df = assemble_gene_df(mc_row, db_k27me3_col)),
  list(id = "37-07d", section = "33c", desc = "mC direction x K119ub DiffBind",
       gene_df = assemble_gene_df(mc_row, db_k119ub_col)),
  list(id = "37-08", section = "20d", desc = "mC direction x expression direction",
       gene_df = assemble_gene_df(mc_row, expr_col))
)

# --- Special test 37-09: K119ub direction x Hyper/Not Hyper at loop anchors ---
cat("\n--- Preparing special tests (37-09, 37-10, 37-11) ---\n")

# Build anchor-gene mapping (reuse section 27 pattern)
cat("  Building GREAT-style anchor-gene associations...\n")
loops <- read.table(LOOPS_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("    Loaded %d loops\n", nrow(loops)))

# Build regulatory domains
build_gene_regulatory_domains <- function() {
  genes_gr <- genes(txdb)
  tss_pos <- ifelse(as.character(strand(genes_gr)) == "-",
                    end(genes_gr), start(genes_gr))
  gene_info <- data.frame(
    entrez_id = names(genes_gr),
    chr = as.character(seqnames(genes_gr)),
    tss = tss_pos,
    strand = as.character(strand(genes_gr)),
    stringsAsFactors = FALSE
  )
  gene_info <- gene_info %>%
    dplyr::arrange(chr, tss) %>%
    dplyr::mutate(
      basal_start = ifelse(strand == "+", tss - GREAT_UPSTREAM, tss - GREAT_DOWNSTREAM),
      basal_end = ifelse(strand == "+", tss + GREAT_DOWNSTREAM, tss + GREAT_UPSTREAM),
      max_start = tss - GREAT_MAX_EXTENSION,
      max_end = tss + GREAT_MAX_EXTENSION
    ) %>%
    dplyr::group_by(chr) %>%
    dplyr::mutate(
      prev_basal_end = dplyr::lag(basal_end, default = -Inf),
      next_basal_start = dplyr::lead(basal_start, default = Inf),
      reg_start = pmax(max_start, prev_basal_end, 1),
      reg_end = pmin(max_end, next_basal_start)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      reg_start = ifelse(reg_end < reg_start, basal_start, reg_start),
      reg_end = ifelse(reg_end < reg_start, basal_end, reg_end),
      reg_end = pmax(reg_end, reg_start)
    ) %>%
    dplyr::select(entrez_id, chr, tss, strand, reg_start, reg_end)
  gene_info
}

gene_domains <- build_gene_regulatory_domains()
gene_domains_gr <- GRanges(
  seqnames = gene_domains$chr,
  ranges = IRanges(start = gene_domains$reg_start, end = gene_domains$reg_end),
  entrez_id = gene_domains$entrez_id
)

# Entrez -> Symbol mapping
entrez_ids <- unique(gene_domains$entrez_id)
entrez_to_symbol <- suppressMessages(
  AnnotationDbi::select(org.Mm.eg.db, keys = entrez_ids,
                        columns = "SYMBOL", keytype = "ENTREZID")
)
entrez_map <- setNames(entrez_to_symbol$SYMBOL, entrez_to_symbol$ENTREZID)

# Build anchor GRanges from differential loops (all loops in this file are
# significant -- they come from edgeR output with filters already applied)
anchor1_gr <- GRanges(seqnames = loops$anchor1_chr,
                      ranges = IRanges(start = loops$anchor1_start, end = loops$anchor1_end),
                      loop_direction = loops$direction,
                      loop_logFC = loops$logFC)
anchor2_gr <- GRanges(seqnames = loops$anchor2_chr,
                      ranges = IRanges(start = loops$anchor2_start, end = loops$anchor2_end),
                      loop_direction = loops$direction,
                      loop_logFC = loops$logFC)
all_anchors <- c(anchor1_gr, anchor2_gr)
mcols(all_anchors)$anchor_id <- paste0(rep(loops$loop_id, 2), c(rep("_A1", nrow(loops)), rep("_A2", nrow(loops))))

# Associate anchors with genes
anchor_gene_overlaps <- findOverlaps(all_anchors, gene_domains_gr, ignore.strand = TRUE)
anchor_gene_assoc <- data.frame(
  anchor_idx = queryHits(anchor_gene_overlaps),
  entrez_id = gene_domains_gr$entrez_id[subjectHits(anchor_gene_overlaps)],
  anchor_chr = as.character(seqnames(all_anchors)[queryHits(anchor_gene_overlaps)]),
  loop_direction = all_anchors$loop_direction[queryHits(anchor_gene_overlaps)],
  stringsAsFactors = FALSE
)
anchor_gene_assoc$symbol <- entrez_map[anchor_gene_assoc$entrez_id]
anchor_gene_assoc <- anchor_gene_assoc %>%
  dplyr::filter(!is.na(symbol) & symbol != "")

cat(sprintf("    Anchor-gene associations: %d\n", nrow(anchor_gene_assoc)))

# Test 37-09: K119ub direction x Hyper/Not Hyper at loop anchors
# Get K119ub overlap per anchor gene
k119ub_ctrl_gr <- load_chip_peaks(K119UB_FILES$ctrl, "K119ub ctrl (37-09)")
k119ub_mut_gr  <- load_chip_peaks(K119UB_FILES$mut, "K119ub mut (37-09)")

anchor_genes_unique <- anchor_gene_assoc %>%
  dplyr::distinct(symbol, .keep_all = TRUE)

# Count K119ub overlaps per anchor gene's regulatory domain
anchor_gene_gr <- gene_domains_gr[match(anchor_genes_unique$entrez_id, gene_domains_gr$entrez_id)]
n_k119_ctrl <- countOverlaps(anchor_gene_gr, k119ub_ctrl_gr)
n_k119_mut  <- countOverlaps(anchor_gene_gr, k119ub_mut_gr)
anchor_genes_unique$k119_net <- n_k119_mut - n_k119_ctrl

# Row: K119ub Gained/Lost; Col: Hyper/Not Hyper
mc_sig_genes <- mc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  dplyr::transmute(gene, is_hyper = mod_difference > 0)

gene_df_09 <- anchor_genes_unique %>%
  dplyr::filter(k119_net != 0) %>%
  dplyr::inner_join(mc_sig_genes, by = c("symbol" = "gene")) %>%
  dplyr::transmute(
    gene = symbol,
    chr = anchor_chr,
    row_label = ifelse(k119_net > 0, "K119ub Gained", "K119ub Lost"),
    col_label = ifelse(is_hyper, "Hyper", "Not Hyper")
  )

cat(sprintf("  Test 37-09: %d genes\n", nrow(gene_df_09)))

test_defs <- c(test_defs, list(
  list(id = "37-09", section = "27c", desc = "K119ub direction x Hyper at loop anchors",
       gene_df = gene_df_09)
))

# --- Special test 37-10: Compartment shift x mC direction ---
cat("  Loading HOMER compartment data for test 37-10...\n")

comp_raw <- read.table(COMPARTMENT_FILE, sep = "\t", header = TRUE,
                       stringsAsFactors = FALSE, check.names = FALSE,
                       comment.char = "", quote = "")

pc1_cols <- grep("bedGraph avg over given bp", names(comp_raw), value = TRUE)
stopifnot("Expected 6 PC1 sample columns" = length(pc1_cols) == 6)
ctrl_pc1_cols <- pc1_cols[1:3]

difference_col <- grep("ctrl vs\\. mut Difference",
                        names(comp_raw), value = TRUE, ignore.case = TRUE)[1]
adj_pvalue_col <- grep("ctrl vs\\. mut adj\\. p-value",
                        names(comp_raw), value = TRUE, ignore.case = TRUE)[1]
stopifnot("Difference column not found" = !is.na(difference_col))
stopifnot("Adj. p-value column not found" = !is.na(adj_pvalue_col))

comp_df <- data.frame(
  chr = comp_raw[["Chr"]],
  start = comp_raw[["Start"]],
  end = comp_raw[["End"]],
  gene_name = comp_raw[["Gene Name"]],
  dist_to_tss = as.numeric(comp_raw[["Distance to TSS"]]),
  ctrl_pc1_1 = as.numeric(comp_raw[[ctrl_pc1_cols[1]]]),
  ctrl_pc1_2 = as.numeric(comp_raw[[ctrl_pc1_cols[2]]]),
  ctrl_pc1_3 = as.numeric(comp_raw[[ctrl_pc1_cols[3]]]),
  difference = as.numeric(comp_raw[[difference_col]]),
  adj_pvalue = as.numeric(comp_raw[[adj_pvalue_col]]),
  stringsAsFactors = FALSE
)

comp_df$mean_ctrl_pc1 <- rowMeans(comp_df[, c("ctrl_pc1_1", "ctrl_pc1_2", "ctrl_pc1_3")],
                                   na.rm = TRUE)

comp_df$shift <- dplyr::case_when(
  comp_df$adj_pvalue < SHIFT_FDR & comp_df$difference > SHIFT_DIFF  ~ "B to A",
  comp_df$adj_pvalue < SHIFT_FDR & comp_df$difference < -SHIFT_DIFF ~ "A to B",
  TRUE ~ "Stable"
)

# Deduplicate: one gene per bin, keep closest to TSS
comp_genes <- comp_df %>%
  dplyr::filter(!is.na(gene_name) & gene_name != "") %>%
  dplyr::group_by(gene_name) %>%
  dplyr::slice_min(abs(dist_to_tss), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

cat(sprintf("    Compartment genes: %d (B->A: %d, A->B: %d, Stable: %d)\n",
            nrow(comp_genes),
            sum(comp_genes$shift == "B to A"),
            sum(comp_genes$shift == "A to B"),
            sum(comp_genes$shift == "Stable")))

# Join with mC DMR data
mc_dedup <- mc_dmr %>%
  dplyr::filter(significant) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  dplyr::transmute(gene, mc_direction = ifelse(mod_difference > 0, "mC Hyper", "mC Hypo"))

comp_mc <- comp_genes %>%
  dplyr::inner_join(mc_dedup, by = c("gene_name" = "gene"))

# 37-10a: B->A vs Stable x mC Hyper/Hypo
gene_df_10a <- comp_mc %>%
  dplyr::filter(shift %in% c("B to A", "Stable")) %>%
  dplyr::transmute(
    gene = gene_name,
    chr,
    row_label = shift,
    col_label = mc_direction
  )

# 37-10b: A->B vs Stable x mC Hyper/Hypo
gene_df_10b <- comp_mc %>%
  dplyr::filter(shift %in% c("A to B", "Stable")) %>%
  dplyr::transmute(
    gene = gene_name,
    chr,
    row_label = shift,
    col_label = mc_direction
  )

cat(sprintf("  Test 37-10a (B->A vs Stable): %d genes\n", nrow(gene_df_10a)))
cat(sprintf("  Test 37-10b (A->B vs Stable): %d genes\n", nrow(gene_df_10b)))

test_defs <- c(test_defs, list(
  list(id = "37-10a", section = "29", desc = "B->A shift vs Stable x mC direction",
       gene_df = gene_df_10a),
  list(id = "37-10b", section = "29", desc = "A->B shift vs Stable x mC direction",
       gene_df = gene_df_10b)
))

# --- Special test 37-11: Loop Gained/Lost x MeCP2 Sig Up/Other ---
mecp2_raw <- read.table(MECP2_FILES$annotated, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, fill = TRUE,
                         quote = "", comment.char = "")
mecp2_raw$Fold <- as.numeric(mecp2_raw$Fold)
mecp2_raw$FDR <- as.numeric(mecp2_raw$FDR)
mecp2_raw$distanceToTSS <- as.numeric(mecp2_raw$distanceToTSS)

mecp2_gene_sig <- mecp2_raw %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(
    nearest_fold = Fold[which.min(abs(distanceToTSS))],
    nearest_fdr  = FDR[which.min(abs(distanceToTSS))],
    .groups = "drop"
  ) %>%
  dplyr::transmute(
    gene = SYMBOL,
    mecp2_status = ifelse(nearest_fdr < 0.05 & nearest_fold > 0, "MeCP2 Sig Up", "Other")
  )

# Anchor genes with loop direction
anchor_loop_dir <- anchor_gene_assoc %>%
  dplyr::distinct(symbol, .keep_all = TRUE) %>%
  dplyr::filter(loop_direction %in% c("up_in_mutant", "down_in_mutant")) %>%
  dplyr::transmute(
    gene = symbol,
    chr = anchor_chr,
    row_label = ifelse(loop_direction == "up_in_mutant", "Loop Gained", "Loop Lost")
  )

gene_df_11 <- anchor_loop_dir %>%
  dplyr::inner_join(mecp2_gene_sig, by = "gene")

# Rename mecp2_status to col_label
gene_df_11 <- gene_df_11 %>%
  dplyr::rename(col_label = mecp2_status)

cat(sprintf("  Test 37-11: %d genes\n", nrow(gene_df_11)))

test_defs <- c(test_defs, list(
  list(id = "37-11", section = "31b", desc = "Loop direction x MeCP2 Sig Up at anchors",
       gene_df = gene_df_11)
))

# Log all test sizes
cat("\n--- Test summary ---\n")
for (td in test_defs) {
  n <- nrow(td$gene_df)
  cat(sprintf("  %s (%s): %d genes", td$id, td$section, n))
  if (n > 0) {
    tab <- table(td$gene_df$row_label, td$gene_df$col_label)
    cat(sprintf(" [%s]\n", paste(dim(tab), collapse = "x")))
  } else {
    cat(" [SKIPPED - no genes]\n")
  }
}

# =============================================================================
# BLOCK 4: RDS CACHE + PERMUTATION RUNNER
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("RUNNING PERMUTATION TESTS\n")
cat("================================================================================\n\n")

cache_path <- file.path(TABLES_DIR, "permutation_37_gene_level.rds")

if (file.exists(cache_path) && !.force_rerun) {
  cat("Loading cached permutation results...\n")
  all_results <- readRDS(cache_path)
  cat(sprintf("  Loaded %d cached test results\n", length(all_results)))
} else {
  all_results <- list()

  for (td in test_defs) {
    gene_df <- td$gene_df

    cat(sprintf("\n--- Test %s: %s (validates section %s) ---\n",
                td$id, td$desc, td$section))

    if (nrow(gene_df) < 20) {
      cat(sprintf("  WARNING: Only %d genes -- skipping (need >= 20)\n", nrow(gene_df)))
      all_results[[td$id]] <- list(
        test_id = td$id, source_section = td$section, description = td$desc,
        n_genes = nrow(gene_df), skipped = TRUE, skip_reason = "Too few genes"
      )
      next
    }

    cat(sprintf("  Total genes: %d\n", nrow(gene_df)))
    cat(sprintf("  Chromosomes: %d\n", length(unique(gene_df$chr))))

    obs_tab <- table(gene_df$row_label, gene_df$col_label)
    cat("  Contingency table:\n")
    print(obs_tab)

    set.seed(PERM_SEED)
    result <- permute_gene_2x2(gene_df, ntimes = PERM_NTIMES)

    result$test_id          <- td$id
    result$source_section   <- td$section
    result$description      <- td$desc
    result$n_genes          <- nrow(gene_df)
    result$n_chromosomes    <- length(unique(gene_df$chr))
    result$contingency_table <- obs_tab
    result$skipped          <- FALSE

    all_results[[td$id]] <- result

    cat(sprintf("  Fisher OR = %.2f, p = %.2e\n", result$observed_or, result$fisher_p))
    cat(sprintf("  Permutation z = %.2f, empirical p = %.4f\n", result$z_score, result$empirical_p))
  }

  saveRDS(all_results, cache_path)
  cat(sprintf("\nSaved permutation cache: %s\n", cache_path))
}

# =============================================================================
# BLOCK 5: POST-PROCESSING -- CONCORDANCE CLASSIFICATION
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("POST-PROCESSING: CONCORDANCE CLASSIFICATION\n")
cat("================================================================================\n\n")

# Build summary data.frame from results
summary_rows <- lapply(names(all_results), function(id) {
  r <- all_results[[id]]
  if (isTRUE(r$skipped)) {
    return(data.frame(
      test_id = r$test_id, source_section = r$source_section,
      description = r$description, n_genes = r$n_genes,
      fisher_or = NA_real_, fisher_log2or = NA_real_, fisher_p = NA_real_,
      perm_z = NA_real_, perm_empirical_p = NA_real_,
      null_ci_lo = NA_real_, null_ci_hi = NA_real_,
      concordance = "Skipped",
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    test_id = r$test_id,
    source_section = r$source_section,
    description = r$description,
    n_genes = r$n_genes,
    fisher_or = r$observed_or,
    fisher_log2or = r$observed_log2or,
    fisher_p = r$fisher_p,
    perm_z = r$z_score,
    perm_empirical_p = r$empirical_p,
    null_ci_lo = r$ci_95[[1]],
    null_ci_hi = r$ci_95[[2]],
    concordance = NA_character_,
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, summary_rows)

# Classify concordance
summary_df$concordance <- dplyr::case_when(
  summary_df$concordance == "Skipped" ~ "Skipped",
  is.na(summary_df$perm_z) ~ "Indeterminate",
  summary_df$fisher_p < 0.05 & summary_df$perm_empirical_p < 0.05 &
    sign(summary_df$fisher_log2or) == sign(summary_df$perm_z) ~ "Confirmed",
  summary_df$fisher_p < 0.05 & summary_df$perm_empirical_p >= 0.05 ~ "Weakened",
  summary_df$fisher_p >= 0.05 & summary_df$perm_empirical_p < 0.05 ~ "Strengthened",
  TRUE ~ "Concordant NS"
)

# Print summary
cat("Concordance classification:\n")
print(table(summary_df$concordance))
cat("\n")

# Full table
cat("Full results:\n")
print(summary_df %>%
        dplyr::select(test_id, source_section, n_genes, fisher_or, fisher_p,
                      perm_z, perm_empirical_p, concordance) %>%
        dplyr::mutate(
          fisher_or = round(fisher_or, 2),
          fisher_p = signif(fisher_p, 3),
          perm_z = round(perm_z, 2),
          perm_empirical_p = round(perm_empirical_p, 4)
        ),
      row.names = FALSE)

# =============================================================================
# BLOCK 6: FIGURES
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("GENERATING FIGURES\n")
cat("================================================================================\n\n")

# Filter to non-skipped tests for plotting
plot_df <- summary_df %>%
  dplyr::filter(concordance != "Skipped" & concordance != "Indeterminate")

# Concordance colors
CONCORDANCE_COLORS <- c(
  "Confirmed" = "#2CA02C",
  "Weakened" = "#D62728",
  "Strengthened" = "#FF7F0E",
  "Concordant NS" = "#7F7F7F"
)

# --- Figure 37a: Forest plot of permutation z-scores ---
cat("Creating Figure 37a: Forest plot of permutation z-scores...\n")

forest_df <- plot_df %>%
  dplyr::mutate(
    test_label = paste0(test_id, " (s", source_section, ")"),
    sig_star = ifelse(perm_empirical_p < 0.05, "*", ""),
    sig_star = ifelse(perm_empirical_p < 0.01, "**", sig_star),
    sig_star = ifelse(perm_empirical_p < 0.001, "***", sig_star)
  ) %>%
  dplyr::arrange(abs(perm_z))

forest_df$test_label <- factor(forest_df$test_label, levels = forest_df$test_label)

# Compute CI in z-score space for display
forest_df$z_ci_lo <- (forest_df$null_ci_lo - 0) / ifelse(abs(forest_df$perm_z) > 0, abs(forest_df$fisher_log2or / forest_df$perm_z), 1)
forest_df$z_ci_hi <- (forest_df$null_ci_hi - 0) / ifelse(abs(forest_df$perm_z) > 0, abs(forest_df$fisher_log2or / forest_df$perm_z), 1)

p_37a <- ggplot(forest_df, aes(x = perm_z, y = test_label, color = source_section)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = "dotted", color = "grey70") +
  geom_point(size = 3) +
  geom_text(aes(label = sig_star), nudge_y = 0.25, size = 5, show.legend = FALSE) +
  geom_text(aes(label = sprintf("n=%d", n_genes)), x = min(forest_df$perm_z, na.rm = TRUE) - 0.5,
            hjust = 1, size = 2.5, show.legend = FALSE) +
  scale_color_brewer(palette = "Set2", name = "Source\nSection") +
  labs(
    title = "Gene-Level Permutation Z-Scores",
    subtitle = sprintf("%d tests, %s permutations per test | * p<0.05  ** p<0.01  *** p<0.001",
                       nrow(forest_df), format(PERM_NTIMES, big.mark = ",")),
    x = "Permutation Z-Score",
    y = NULL
  ) +
  theme_biomodal() +
  theme(axis.text.y = element_text(size = 9))

save_multiformat_ggplot(p_37a, file.path(SECTION_DIR, "37a_zscore_forest_plot"),
                        width = 12, height = max(8, nrow(forest_df) * 0.5 + 2))

# --- Figure 37b: Null distribution histograms for top 4 ---
cat("Creating Figure 37b: Null distribution histograms...\n")

top4_ids <- plot_df %>%
  dplyr::arrange(desc(abs(perm_z))) %>%
  dplyr::slice_head(n = 4) %>%
  dplyr::pull(test_id)

hist_plots <- lapply(top4_ids, function(id) {
  r <- all_results[[id]]
  null_df <- data.frame(log2or = r$null_distribution)
  obs_val <- r$observed_log2or

  # Determine tail region for shading
  tail_threshold <- abs(obs_val)

  ggplot(null_df, aes(x = log2or)) +
    geom_histogram(bins = 60, fill = "grey70", color = "grey50", linewidth = 0.2) +
    geom_vline(xintercept = obs_val, color = "#D62728", linewidth = 1) +
    geom_vline(xintercept = -obs_val, color = "#D62728", linewidth = 0.5, linetype = "dashed") +
    annotate("rect",
             xmin = tail_threshold, xmax = max(null_df$log2or, obs_val) + 0.5,
             ymin = -Inf, ymax = Inf, fill = "#D62728", alpha = 0.15) +
    annotate("rect",
             xmin = min(null_df$log2or, -obs_val) - 0.5, xmax = -tail_threshold,
             ymin = -Inf, ymax = Inf, fill = "#D62728", alpha = 0.15) +
    labs(
      title = sprintf("%s (z=%.1f, p=%.3f)", r$test_id, r$z_score, r$empirical_p),
      subtitle = r$description,
      x = "Null log2(OR)",
      y = "Count"
    ) +
    theme_biomodal(base_size = 10) +
    theme(plot.title = element_text(size = 11))
})

p_37b <- wrap_plots(hist_plots, ncol = 2, nrow = 2) +
  plot_annotation(
    title = "Null Distributions for Top 4 Strongest Effects",
    subtitle = "Grey = null, red line = observed, shaded = empirical p-value tail area",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                  plot.subtitle = element_text(size = 11, hjust = 0.5))
  )

save_multiformat_ggplot(p_37b, file.path(SECTION_DIR, "37b_null_distribution_top4"),
                        width = 14, height = 10)

# --- Figure 37c: Comparison table ---
cat("Creating Figure 37c: Fisher vs permutation comparison table...\n")

table_df <- summary_df %>%
  dplyr::transmute(
    `Test ID` = test_id,
    Section = source_section,
    N = n_genes,
    `Fisher OR` = ifelse(is.na(fisher_or), "—", sprintf("%.2f", fisher_or)),
    `Fisher p` = ifelse(is.na(fisher_p), "—", ifelse(fisher_p < 1e-10, sprintf("%.1e", fisher_p), sprintf("%.3e", fisher_p))),
    `Perm z` = ifelse(is.na(perm_z), "—", sprintf("%.2f", perm_z)),
    `Perm p` = ifelse(is.na(perm_empirical_p), "—", sprintf("%.4f", perm_empirical_p)),
    Result = concordance
  )

# Create a grob table
tt <- gridExtra::ttheme_default(
  core = list(fg_params = list(fontsize = 9)),
  colhead = list(fg_params = list(fontsize = 10, fontface = "bold"))
)
table_grob <- gridExtra::tableGrob(table_df, rows = NULL, theme = tt)

# Wrap in ggplot for save_multiformat_ggplot
p_37c <- ggplot() +
  annotation_custom(table_grob) +
  labs(title = "Fisher vs Permutation Comparison",
       subtitle = sprintf("Gene-level 2x2 tests | %s permutations | chr-stratified shuffle",
                          format(PERM_NTIMES, big.mark = ","))) +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        plot.margin = margin(10, 10, 10, 10))

save_multiformat_ggplot(p_37c, file.path(SECTION_DIR, "37c_fisher_vs_permutation_table"),
                        width = 14, height = max(6, nrow(table_df) * 0.4 + 3))

# --- Figure 37d: Scatter log2(OR) vs z-score ---
cat("Creating Figure 37d: log2(OR) vs z-score scatter...\n")

p_37d <- ggplot(plot_df, aes(x = fisher_log2or, y = perm_z)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "lm", se = TRUE, color = "grey40", linewidth = 0.5, alpha = 0.2) +
  geom_point(aes(color = concordance), size = 3) +
  geom_text_repel(aes(label = test_id, color = concordance),
                  size = 3, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = CONCORDANCE_COLORS, name = "Concordance") +
  labs(
    title = "Fisher log2(OR) vs Permutation Z-Score",
    subtitle = "Points along diagonal = Fisher and permutation agree in magnitude and direction",
    x = "Observed log2(OR) from Fisher's Exact Test",
    y = "Permutation Z-Score"
  ) +
  theme_biomodal()

# Add correlation annotation
if (nrow(plot_df) >= 3) {
  cor_test <- cor.test(plot_df$fisher_log2or, plot_df$perm_z, method = "spearman")
  p_37d <- p_37d +
    annotate("text", x = max(plot_df$fisher_log2or, na.rm = TRUE),
             y = min(plot_df$perm_z, na.rm = TRUE),
             label = sprintf("Spearman rho = %.2f, p = %.2e", cor_test$estimate, cor_test$p.value),
             hjust = 1, vjust = 0, size = 3.5, color = "grey30")
}

save_multiformat_ggplot(p_37d, file.path(SECTION_DIR, "37d_log2or_vs_zscore_scatter"),
                        width = 10, height = 9)

# =============================================================================
# BLOCK 7: EXPORTS + SUMMARY
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("EXPORTING RESULTS\n")
cat("================================================================================\n\n")

# Save summary TSV
write.table(summary_df, file.path(TABLES_DIR, "permutation_37_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Saved: %s\n", file.path(TABLES_DIR, "permutation_37_summary.tsv")))

# Final summary
cat("\n")
cat("===========================================================================\n")
cat("  SECTION 37 COMPLETE\n")
cat("===========================================================================\n\n")

cat(sprintf("Total tests: %d\n", nrow(summary_df)))
cat(sprintf("  Confirmed:     %d\n", sum(summary_df$concordance == "Confirmed")))
cat(sprintf("  Weakened:      %d\n", sum(summary_df$concordance == "Weakened")))
cat(sprintf("  Strengthened:  %d\n", sum(summary_df$concordance == "Strengthened")))
cat(sprintf("  Concordant NS: %d\n", sum(summary_df$concordance == "Concordant NS")))
cat(sprintf("  Skipped:       %d\n", sum(summary_df$concordance == "Skipped")))
cat(sprintf("  Indeterminate: %d\n", sum(summary_df$concordance == "Indeterminate")))

cat(sprintf("\nOutputs:\n"))
cat(sprintf("  Cache: %s\n", cache_path))
cat(sprintf("  Table: %s\n", file.path(TABLES_DIR, "permutation_37_summary.tsv")))
cat(sprintf("  Plots: %s/\n", SECTION_DIR))
cat("\nDone.\n")
