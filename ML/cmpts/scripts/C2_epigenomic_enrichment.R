# ML/cmpts/scripts/C2_epigenomic_enrichment.R
# Stage C2: Publication-quality epigenomic enrichment heatmaps (ComplexHeatmap).
#
# Reads pre-computed A3 enrichment/differential matrices and bin-level signals.
# Produces ComplexHeatmap-based figures with row/column annotations, category
# splitting, and Wilcoxon significance overlays. No BigWig extraction needed.
#
# Usage:
#   Rscript C2_epigenomic_enrichment.R <data_root> <code_root>
#     <data_root> : HPC data directory (kept for CLI consistency)
#     <code_root> : repo directory (e.g. /expanse/.../mariner_hi-c/ML/cmpts)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript C2_epigenomic_enrichment.R <data_root> <code_root>")
}

DATA_ROOT <- args[1]
CODE_ROOT <- args[2]

# -- Constants ----------------------------------------------------------------

TPS <- c("250402", "250831")
TP_LABELS <- c("250402" = "Late (adult)", "250831" = "Early (P12)")
TP_DIRS   <- c("250402" = "late", "250831" = "early")

LABEL_ORDER <- c("A.1", "A.2", "B.1", "B.2")
LABEL_COLORS <- c(
  "A.1" = "#e41a1c",
  "A.2" = "#ff7f00",
  "B.1" = "#4daf4a",
  "B.2" = "#377eb8"
)

MARK_ORDER <- c("H3K27ac", "H3K4me3", "ATAC", "H3K36me3", "RNA",
                "DNAmethylation", "H3K27me1", "H2AK119ub", "H3K27me3")

MARK_CATEGORIES <- c(
  "H3K27ac"        = "Active",
  "H3K4me3"        = "Active",
  "ATAC"           = "Active",
  "H3K36me3"       = "Gene body",
  "RNA"            = "Gene body",
  "DNAmethylation" = "Methylation",
  "H3K27me1"       = "Methylation",
  "H2AK119ub"      = "Repressive",
  "H3K27me3"       = "Repressive"
)

CATEGORY_ORDER  <- c("Active", "Gene body", "Methylation", "Repressive")
CATEGORY_COLORS <- c(
  "Active"      = "#e41a1c",
  "Gene body"   = "#ff7f00",
  "Methylation" = "#4daf4a",
  "Repressive"  = "#984ea3"
)

COMPARTMENT_SPLIT <- c("A.1" = "A", "A.2" = "A", "B.1" = "B", "B.2" = "B")

OUT_DIR   <- file.path(CODE_ROOT, "outputs", "calder2", "enrichment")
UTIL_PATH <- file.path(CODE_ROOT, "scripts", "utils", "multi_format_output.R")

# -- Libraries ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggplot2)
})
source(UTIL_PATH)

# -- Header -------------------------------------------------------------------

start_time <- proc.time()

cat("===========================================\n")
cat("C2: Epigenomic Enrichment (Publication)\n")
cat("===========================================\n")
cat(sprintf("CODE_ROOT:  %s\n", CODE_ROOT))
cat(sprintf("Output dir: %s\n", OUT_DIR))
cat(sprintf("Timepoints: %s\n", paste(TPS, collapse = ", ")))
cat(sprintf("Marks:      %d\n", length(MARK_ORDER)))
cat(sprintf("Start:      %s\n", date()))
cat("===========================================\n\n")

# -- Pre-flight validation ----------------------------------------------------

cat("=== Pre-flight validation ===\n")

enrich_paths <- list()
diff_paths   <- list()
signal_paths <- list()

for (tp in TPS) {
  tp_dir <- file.path(CODE_ROOT, "outputs", "calder2", TP_DIRS[tp])

  ep <- file.path(tp_dir, sprintf("%s_enrichment_matrix.tsv", tp))
  dp <- file.path(tp_dir, sprintf("%s_differential_matrix.tsv", tp))
  sp <- file.path(tp_dir, sprintf("%s_bin_signals.tsv", tp))

  if (!file.exists(ep)) stop(sprintf("Missing: %s", ep))
  if (!file.exists(dp)) stop(sprintf("Missing: %s", dp))
  if (!file.exists(sp)) stop(sprintf("Missing: %s", sp))

  enrich_paths[[tp]] <- ep
  diff_paths[[tp]]   <- dp
  signal_paths[[tp]] <- sp

  cat(sprintf("  OK: %s enrichment (%s bytes), differential (%s bytes), signals (%s bytes)\n",
              tp,
              format(file.info(ep)$size, big.mark = ","),
              format(file.info(dp)$size, big.mark = ","),
              format(file.info(sp)$size, big.mark = ",")))
}

if (!file.exists(UTIL_PATH)) stop(sprintf("Missing utility: %s", UTIL_PATH))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("  Pre-flight OK.\n\n")

# -- Function definitions -----------------------------------------------------

reshape_to_matrix <- function(dt, value_col = "log2_fold") {
  wide <- dcast(dt, mark ~ subcompartment, value.var = value_col)
  mat <- as.matrix(wide[, LABEL_ORDER, with = FALSE])
  rownames(mat) <- wide$mark
  mat <- mat[MARK_ORDER, , drop = FALSE]
  mat
}

make_row_annotation <- function(show_legend = TRUE) {
  cats <- factor(MARK_CATEGORIES[MARK_ORDER], levels = CATEGORY_ORDER)
  rowAnnotation(
    Category = cats,
    col = list(Category = CATEGORY_COLORS),
    show_annotation_name = TRUE,
    show_legend = show_legend,
    annotation_name_gp = gpar(fontsize = 9),
    annotation_legend_param = list(title = "Mark category")
  )
}

make_col_annotation <- function(n_bins_vec, show_legend = TRUE) {
  subcomps <- factor(LABEL_ORDER, levels = LABEL_ORDER)
  HeatmapAnnotation(
    Subcompartment = subcomps,
    `n bins` = anno_barplot(
      n_bins_vec,
      gp = gpar(fill = LABEL_COLORS[LABEL_ORDER]),
      height = unit(1.5, "cm"),
      axis_param = list(gp = gpar(fontsize = 7))
    ),
    col = list(Subcompartment = LABEL_COLORS),
    show_annotation_name = TRUE,
    show_legend = show_legend,
    annotation_name_gp = gpar(fontsize = 9),
    annotation_legend_param = list(title = "Subcompartment")
  )
}

build_enrichment_heatmap <- function(mat, fold_mat, title, n_bins_vec,
                                     show_row_anno = TRUE,
                                     show_col_anno = TRUE,
                                     show_legend = TRUE,
                                     col_cap = NULL) {
  finite_vals <- mat[is.finite(mat)]
  if (is.null(col_cap)) {
    col_cap <- min(5, ceiling(max(abs(finite_vals), na.rm = TRUE)))
  }

  col_fun <- colorRamp2(c(-col_cap, 0, col_cap), c("#2166AC", "white", "#B2182B"))

  row_split <- factor(MARK_CATEGORIES[rownames(mat)], levels = CATEGORY_ORDER)
  col_split <- factor(COMPARTMENT_SPLIT[colnames(mat)], levels = c("A", "B"))

  row_anno <- if (show_row_anno) make_row_annotation(show_legend) else NULL
  col_anno <- if (show_col_anno) make_col_annotation(n_bins_vec, show_legend) else NULL

  fold_for_label <- fold_mat

  Heatmap(
    mat,
    name             = "log2_fold",
    col              = col_fun,
    na_col           = "grey80",
    column_title     = title,
    column_title_gp  = gpar(fontsize = 12, fontface = "bold"),
    cluster_rows     = FALSE,
    cluster_columns  = FALSE,
    row_split        = row_split,
    column_split     = col_split,
    row_gap          = unit(2, "mm"),
    column_gap       = unit(2, "mm"),
    row_title_gp     = gpar(fontsize = 9),
    column_names_gp  = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp     = gpar(fontsize = 9),
    left_annotation  = row_anno,
    top_annotation   = col_anno,
    show_heatmap_legend = show_legend,
    heatmap_legend_param = list(
      title = expression(log[2]~"fold enrichment"),
      title_gp = gpar(fontsize = 9),
      labels_gp = gpar(fontsize = 8)
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      val <- fold_for_label[i, j]
      if (is.finite(val)) {
        grid.text(sprintf("%.2f", val), x, y,
                  gp = gpar(fontsize = 7, col = ifelse(abs(val) > col_cap * 0.7, "white", "black")))
      }
    },
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 1)
  )
}

build_differential_heatmap <- function(mat, pval_mat, title,
                                       show_row_anno = TRUE,
                                       show_legend = TRUE,
                                       col_cap = NULL) {
  finite_vals <- mat[is.finite(mat)]
  if (is.null(col_cap)) {
    col_cap <- min(2, ceiling(max(abs(finite_vals), na.rm = TRUE) * 10) / 10)
  }

  col_fun <- colorRamp2(c(-col_cap, 0, col_cap), c("#2166AC", "white", "#B2182B"))

  row_split <- factor(MARK_CATEGORIES[rownames(mat)], levels = CATEGORY_ORDER)
  col_split <- factor(COMPARTMENT_SPLIT[colnames(mat)], levels = c("A", "B"))

  row_anno <- if (show_row_anno) make_row_annotation(show_legend) else NULL

  Heatmap(
    mat,
    name             = "log2_fc",
    col              = col_fun,
    na_col           = "grey80",
    column_title     = title,
    column_title_gp  = gpar(fontsize = 12, fontface = "bold"),
    cluster_rows     = FALSE,
    cluster_columns  = FALSE,
    row_split        = row_split,
    column_split     = col_split,
    row_gap          = unit(2, "mm"),
    column_gap       = unit(2, "mm"),
    row_title_gp     = gpar(fontsize = 9),
    column_names_gp  = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp     = gpar(fontsize = 9),
    left_annotation  = row_anno,
    show_heatmap_legend = show_legend,
    heatmap_legend_param = list(
      title = expression(log[2]~"(mut / ctrl)"),
      title_gp = gpar(fontsize = 9),
      labels_gp = gpar(fontsize = 8)
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      val <- mat[i, j]
      pv  <- pval_mat[i, j]
      if (is.finite(val)) {
        label <- sprintf("%.2f", val)
        if (!is.na(pv) && pv < 0.05) label <- paste0(label, "*")
        grid.text(label, x, y,
                  gp = gpar(fontsize = 7, col = ifelse(abs(val) > col_cap * 0.7, "white", "black")))
      }
    },
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 1)
  )
}

compute_wilcoxon <- function(bin_signals_dt, tp) {
  results <- list()
  for (mark in MARK_ORDER) {
    ctrl_col <- paste0(mark, "_ctrl")
    mut_col  <- paste0(mark, "_mut")
    if (!ctrl_col %in% names(bin_signals_dt)) next

    for (subcomp in LABEL_ORDER) {
      bins <- bin_signals_dt[ctrl_label == subcomp]
      ctrl_vals <- bins[[ctrl_col]]
      mut_vals  <- bins[[mut_col]]

      valid <- is.finite(ctrl_vals) & is.finite(mut_vals)
      ctrl_vals <- ctrl_vals[valid]
      mut_vals  <- mut_vals[valid]

      if (length(ctrl_vals) < 10) {
        results[[length(results) + 1L]] <- data.table(
          mark = mark, subcompartment = subcomp,
          n_bins = length(ctrl_vals), p_value = NA_real_, tp = tp
        )
        next
      }

      wt <- suppressWarnings(wilcox.test(ctrl_vals, mut_vals, paired = TRUE))
      results[[length(results) + 1L]] <- data.table(
        mark = mark, subcompartment = subcomp,
        n_bins = length(ctrl_vals), p_value = wt$p.value, tp = tp
      )
    }
  }
  out <- rbindlist(results)
  out[, padj := p.adjust(p_value, method = "BH"), by = tp]
  out
}

pval_to_matrix <- function(stat_dt) {
  wide <- dcast(stat_dt, mark ~ subcompartment, value.var = "padj")
  mat <- as.matrix(wide[, LABEL_ORDER, with = FALSE])
  rownames(mat) <- wide$mark
  mat[MARK_ORDER, , drop = FALSE]
}

# -- Phase 1: Load and reshape data -------------------------------------------

cat("=== Phase 1: Loading A3 enrichment data ===\n")

enrich_data <- list()
diff_data   <- list()
signal_data <- list()
n_bins_data <- list()

for (tp in TPS) {
  enrich_data[[tp]] <- fread(enrich_paths[[tp]])
  diff_data[[tp]]   <- fread(diff_paths[[tp]])
  signal_data[[tp]] <- fread(signal_paths[[tp]])

  cat(sprintf("  %s (%s): enrichment=%d rows, differential=%d rows, signals=%d bins\n",
              tp, TP_LABELS[tp], nrow(enrich_data[[tp]]),
              nrow(diff_data[[tp]]), nrow(signal_data[[tp]])))
}

cat("\n  Reshaping to matrices...\n")

enrich_matrices <- list()
fold_matrices   <- list()
diff_matrices   <- list()

for (tp in TPS) {
  for (cond in c("ctrl", "mut")) {
    key <- paste(tp, cond, sep = "_")
    subset_dt <- enrich_data[[tp]][condition == cond]
    enrich_matrices[[key]] <- reshape_to_matrix(subset_dt, "log2_fold")
    fold_matrices[[key]]   <- reshape_to_matrix(subset_dt, "fold_enrichment")

    nbins <- subset_dt[mark == MARK_ORDER[1], .(subcompartment, n_bins)]
    nbins <- nbins[match(LABEL_ORDER, subcompartment)]
    n_bins_data[[key]] <- nbins$n_bins
  }

  diff_matrices[[tp]] <- reshape_to_matrix(diff_data[[tp]], "log2_fc")
}

cat(sprintf("  Built %d enrichment matrices + %d differential matrices.\n",
            length(enrich_matrices), length(diff_matrices)))

# -- Phase 2: Wilcoxon statistics ---------------------------------------------

cat("\n=== Phase 2: Wilcoxon rank-sum tests (ctrl vs mut) ===\n")

stat_list <- list()
pval_matrices <- list()

for (tp in TPS) {
  stat_dt <- compute_wilcoxon(signal_data[[tp]], tp)
  stat_list[[tp]] <- stat_dt
  pval_matrices[[tp]] <- pval_to_matrix(stat_dt)

  n_sig <- sum(stat_dt$padj < 0.05, na.rm = TRUE)
  cat(sprintf("  %s: %d/%d mark-subcompartment pairs significant (BH padj < 0.05)\n",
              tp, n_sig, nrow(stat_dt)))
}

all_stats <- rbindlist(stat_list)

# -- Phase 3: Build and render figures ----------------------------------------

cat("\n=== Phase 3: Generating figures ===\n")

# ---- Figure 1 & 2: Per-TP ctrl enrichment heatmaps ----

for (tp in TPS) {
  key <- paste(tp, "ctrl", sep = "_")
  fig_name <- sprintf("c2_enrichment_ctrl_%s", TP_DIRS[tp])

  cat(sprintf("  %s: %s\n", fig_name, TP_LABELS[tp]))

  ht <- build_enrichment_heatmap(
    mat        = enrich_matrices[[key]],
    fold_mat   = fold_matrices[[key]],
    title      = sprintf("Ctrl Subcompartment Enrichment: %s", TP_LABELS[tp]),
    n_bins_vec = n_bins_data[[key]]
  )

  save_multiformat_base(
    bquote(draw(.(ht),
                heatmap_legend_side = "right",
                annotation_legend_side = "right",
                padding = unit(c(2, 2, 2, 10), "mm"))),
    file.path(OUT_DIR, fig_name),
    width = 7, height = 7
  )
}

# ---- Figure 3: Main publication panel (2 TP x 2 conditions) ----

cat("  c2_enrichment_panel: combined 4-panel figure\n")

panel_cols <- character(0)
panel_mat  <- NULL
panel_fold <- NULL
panel_split <- character(0)
panel_labels <- character(0)
panel_subcomp_colors <- character(0)

for (tp in TPS) {
  for (cond in c("ctrl", "mut")) {
    key <- paste(tp, cond, sep = "_")
    m   <- enrich_matrices[[key]]
    f   <- fold_matrices[[key]]

    group_label <- sprintf("%s %s",
                           TP_LABELS[tp],
                           ifelse(cond == "ctrl", "Ctrl", "Mut"))
    new_names <- paste(group_label, colnames(m), sep = "\n")

    colnames(m) <- new_names
    colnames(f) <- new_names

    panel_mat   <- if (is.null(panel_mat)) m else cbind(panel_mat, m)
    panel_fold  <- if (is.null(panel_fold)) f else cbind(panel_fold, f)
    panel_split <- c(panel_split, rep(group_label, ncol(m)))
    panel_labels <- c(panel_labels, LABEL_ORDER)
    panel_subcomp_colors <- c(panel_subcomp_colors, LABEL_COLORS[LABEL_ORDER])
  }
}

panel_split_factor <- factor(panel_split, levels = unique(panel_split))

finite_vals <- panel_mat[is.finite(panel_mat)]
panel_cap <- min(5, ceiling(max(abs(finite_vals), na.rm = TRUE)))
panel_col_fun <- colorRamp2(c(-panel_cap, 0, panel_cap), c("#2166AC", "white", "#B2182B"))

row_split_panel <- factor(MARK_CATEGORIES[rownames(panel_mat)], levels = CATEGORY_ORDER)

panel_col_anno <- HeatmapAnnotation(
  Subcompartment = panel_labels,
  col = list(Subcompartment = LABEL_COLORS),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  annotation_legend_param = list(title = "Subcompartment")
)

panel_row_anno <- make_row_annotation(show_legend = TRUE)

fold_for_cell <- panel_fold

ht_panel <- Heatmap(
  panel_mat,
  name             = "log2_fold",
  col              = panel_col_fun,
  na_col           = "grey80",
  column_title     = "Epigenomic Mark Enrichment in CALDER2 Subcompartments",
  column_title_gp  = gpar(fontsize = 14, fontface = "bold"),
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  row_split        = row_split_panel,
  column_split     = panel_split_factor,
  row_gap          = unit(2, "mm"),
  column_gap       = unit(3, "mm"),
  row_title_gp     = gpar(fontsize = 9),
  column_labels    = panel_labels,
  column_names_gp  = gpar(fontsize = 9, fontface = "bold"),
  row_names_gp     = gpar(fontsize = 9),
  left_annotation  = panel_row_anno,
  top_annotation   = panel_col_anno,
  heatmap_legend_param = list(
    title = expression(log[2]~"fold enrichment"),
    title_gp = gpar(fontsize = 9),
    labels_gp = gpar(fontsize = 8)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- fold_for_cell[i, j]
    if (is.finite(val)) {
      grid.text(sprintf("%.1f", val), x, y,
                gp = gpar(fontsize = 6,
                          col = ifelse(abs(panel_mat[i, j]) > panel_cap * 0.7,
                                       "white", "black")))
    }
  },
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.5)
)

save_multiformat_base(
  bquote(draw(.(ht_panel),
              heatmap_legend_side = "bottom",
              annotation_legend_side = "bottom",
              padding = unit(c(2, 10, 2, 2), "mm"))),
  file.path(OUT_DIR, "c2_enrichment_panel"),
  width = 13, height = 8
)

# ---- Figure 4: Differential panel (2 TPs side by side) ----

cat("  c2_differential_panel: differential with significance\n")

diff_combined_mat  <- NULL
diff_combined_pval <- NULL
diff_split <- character(0)
diff_labels <- character(0)

for (tp in TPS) {
  m <- diff_matrices[[tp]]
  p <- pval_matrices[[tp]]

  group_label <- TP_LABELS[tp]
  new_names <- paste(group_label, colnames(m), sep = "\n")

  colnames(m) <- new_names
  colnames(p) <- new_names

  diff_combined_mat  <- if (is.null(diff_combined_mat)) m else cbind(diff_combined_mat, m)
  diff_combined_pval <- if (is.null(diff_combined_pval)) p else cbind(diff_combined_pval, p)
  diff_split  <- c(diff_split, rep(group_label, ncol(m)))
  diff_labels <- c(diff_labels, LABEL_ORDER)
}

diff_split_factor <- factor(diff_split, levels = unique(diff_split))

diff_finite <- diff_combined_mat[is.finite(diff_combined_mat)]
diff_cap <- min(2, ceiling(max(abs(diff_finite), na.rm = TRUE) * 10) / 10)
diff_col_fun <- colorRamp2(c(-diff_cap, 0, diff_cap), c("#2166AC", "white", "#B2182B"))

diff_row_split <- factor(MARK_CATEGORIES[rownames(diff_combined_mat)], levels = CATEGORY_ORDER)

diff_col_anno <- HeatmapAnnotation(
  Subcompartment = diff_labels,
  col = list(Subcompartment = LABEL_COLORS),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  annotation_legend_param = list(title = "Subcompartment")
)

diff_row_anno <- make_row_annotation(show_legend = TRUE)

diff_mat_for_cell  <- diff_combined_mat
diff_pval_for_cell <- diff_combined_pval

ht_diff <- Heatmap(
  diff_combined_mat,
  name             = "log2_fc",
  col              = diff_col_fun,
  na_col           = "grey80",
  column_title     = "Epigenomic Signal Change (Mut vs Ctrl) in Ctrl-Defined Subcompartments",
  column_title_gp  = gpar(fontsize = 12, fontface = "bold"),
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  row_split        = diff_row_split,
  column_split     = diff_split_factor,
  row_gap          = unit(2, "mm"),
  column_gap       = unit(3, "mm"),
  row_title_gp     = gpar(fontsize = 9),
  column_labels    = diff_labels,
  column_names_gp  = gpar(fontsize = 10, fontface = "bold"),
  row_names_gp     = gpar(fontsize = 9),
  left_annotation  = diff_row_anno,
  top_annotation   = diff_col_anno,
  heatmap_legend_param = list(
    title = expression(log[2]~"(mut / ctrl)"),
    title_gp = gpar(fontsize = 9),
    labels_gp = gpar(fontsize = 8)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- diff_mat_for_cell[i, j]
    pv  <- diff_pval_for_cell[i, j]
    if (is.finite(val)) {
      label <- sprintf("%.2f", val)
      if (!is.na(pv) && pv < 0.05) label <- paste0(label, "*")
      grid.text(label, x, y,
                gp = gpar(fontsize = 7,
                          col = ifelse(abs(val) > diff_cap * 0.7, "white", "black")))
    }
  },
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1)
)

save_multiformat_base(
  bquote(draw(.(ht_diff),
              heatmap_legend_side = "bottom",
              annotation_legend_side = "bottom",
              padding = unit(c(2, 10, 2, 2), "mm"))),
  file.path(OUT_DIR, "c2_differential_panel"),
  width = 10, height = 7
)

# ---- Figure 5: Enrichment dot plot (ggplot2) ----

cat("  c2_enrichment_dotplot: dot plot\n")

all_enrich <- rbindlist(enrich_data)
all_enrich[, mark := factor(mark, levels = rev(MARK_ORDER))]
all_enrich[, subcompartment := factor(subcompartment, levels = LABEL_ORDER)]
all_enrich[, condition := factor(ifelse(condition == "ctrl", "Ctrl", "Mut"),
                                 levels = c("Ctrl", "Mut"))]
all_enrich[, tp_label := factor(TP_LABELS[as.character(tp)], levels = TP_LABELS)]
all_enrich[, category := MARK_CATEGORIES[as.character(mark)]]

fold_cap <- 15
all_enrich[, fold_capped := pmin(fold_enrichment, fold_cap)]

log2_abs_max <- max(abs(all_enrich$log2_fold[is.finite(all_enrich$log2_fold)]),
                    na.rm = TRUE)
log2_cap <- min(5, ceiling(log2_abs_max))

p_dot <- ggplot(all_enrich, aes(x = subcompartment, y = mark)) +
  geom_point(aes(size = fold_capped, fill = log2_fold),
             shape = 21, color = "grey30", stroke = 0.3) +
  scale_size_continuous(
    name = "Fold enrichment",
    range = c(1, 10),
    breaks = c(1, 3, 5, 10, 15),
    limits = c(0, fold_cap)
  ) +
  scale_fill_gradient2(
    name = expression(log[2]~"fold"),
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-log2_cap, log2_cap),
    oob = scales::squish, na.value = "grey80"
  ) +
  facet_grid(rows = vars(), cols = vars(tp_label, condition)) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 9),
    strip.text  = element_text(face = "bold", size = 10),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.box       = "horizontal"
  ) +
  labs(x = NULL, y = NULL,
       title = "Epigenomic Mark Enrichment in CALDER2 Subcompartments")

save_multiformat_ggplot(p_dot, file.path(OUT_DIR, "c2_enrichment_dotplot"),
                        width = 12, height = 8)

# -- Phase 4: Write TSV outputs -----------------------------------------------

cat("\n=== Phase 4: Writing TSV outputs ===\n")

# 4a. Wide-format enrichment
wide_rows <- list()
for (tp in TPS) {
  for (cond in c("ctrl", "mut")) {
    key <- paste(tp, cond, sep = "_")
    mat <- fold_matrices[[key]]
    dt <- data.table(
      mark = rownames(mat),
      tp   = tp,
      condition = cond
    )
    for (sc in LABEL_ORDER) dt[, (sc) := mat[, sc]]
    wide_rows[[length(wide_rows) + 1L]] <- dt
  }
}
enrichment_wide <- rbindlist(wide_rows)
out_ew <- file.path(OUT_DIR, "enrichment_wide.tsv")
fwrite(enrichment_wide, out_ew, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("  Written: %s (%d rows)\n", basename(out_ew), nrow(enrichment_wide)))

# 4b. Wide-format differential
diff_rows <- list()
for (tp in TPS) {
  mat <- diff_matrices[[tp]]
  dt <- data.table(mark = rownames(mat), tp = tp)
  for (sc in LABEL_ORDER) dt[, (sc) := mat[, sc]]
  diff_rows[[length(diff_rows) + 1L]] <- dt
}
differential_wide <- rbindlist(diff_rows)
out_dw <- file.path(OUT_DIR, "differential_wide.tsv")
fwrite(differential_wide, out_dw, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("  Written: %s (%d rows)\n", basename(out_dw), nrow(differential_wide)))

# 4c. Enrichment statistics (Wilcoxon)
out_st <- file.path(OUT_DIR, "enrichment_statistics.tsv")
fwrite(all_stats, out_st, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("  Written: %s (%d rows)\n", basename(out_st), nrow(all_stats)))

# -- Verification -------------------------------------------------------------

cat("\n=== Verification ===\n")

checks_passed <- 0L
checks_total  <- 0L

# Check 1: H3K27ac ctrl enrichment monotonically decreasing (late)
checks_total <- checks_total + 1L
k27ac_late <- fold_matrices[["250402_ctrl"]]["H3K27ac", ]
if (k27ac_late["A.1"] > k27ac_late["A.2"] &&
    k27ac_late["A.2"] > k27ac_late["B.1"] &&
    k27ac_late["B.1"] > k27ac_late["B.2"]) {
  cat("  [PASS] H3K27ac ctrl gradient: A.1 > A.2 > B.1 > B.2 (late)\n")
  checks_passed <- checks_passed + 1L
} else {
  cat("  [WARN] H3K27ac ctrl gradient not monotonic (late)\n")
}

# Check 2: H3K27me3 B.1 enrichment >= 1.5 (late)
checks_total <- checks_total + 1L
k27me3_b1 <- fold_matrices[["250402_ctrl"]]["H3K27me3", "B.1"]
if (!is.na(k27me3_b1) && k27me3_b1 >= 1.5) {
  cat(sprintf("  [PASS] H3K27me3 B.1 fold = %.2f (>= 1.5, late)\n", k27me3_b1))
  checks_passed <- checks_passed + 1L
} else {
  cat(sprintf("  [WARN] H3K27me3 B.1 fold = %.2f (< 1.5, late)\n",
              ifelse(is.na(k27me3_b1), NA, k27me3_b1)))
}

# Check 3: Matrix dimensions
checks_total <- checks_total + 1L
dims_ok <- all(sapply(enrich_matrices, function(m) identical(dim(m), c(9L, 4L)))) &&
           all(sapply(diff_matrices, function(m) identical(dim(m), c(9L, 4L))))
if (dims_ok) {
  cat("  [PASS] All matrices are 9 x 4\n")
  checks_passed <- checks_passed + 1L
} else {
  cat("  [FAIL] Matrix dimension mismatch\n")
}

# Check 4: No NaN/Inf in enrichment
checks_total <- checks_total + 1L
has_bad <- any(sapply(enrich_matrices, function(m) any(is.nan(m) | is.infinite(m))))
if (!has_bad) {
  cat("  [PASS] No NaN/Inf in enrichment matrices\n")
  checks_passed <- checks_passed + 1L
} else {
  cat("  [WARN] NaN or Inf found in enrichment matrices\n")
}

# Check 5: Output file completeness
checks_total <- checks_total + 1L
expected_fig_dirs <- c(
  "c2_enrichment_ctrl_late", "c2_enrichment_ctrl_early",
  "c2_enrichment_panel", "c2_differential_panel",
  "c2_enrichment_dotplot"
)
fig_file_count <- 0L
for (fd in expected_fig_dirs) {
  fig_path <- file.path(OUT_DIR, fd)
  if (dir.exists(fig_path)) {
    fig_file_count <- fig_file_count + length(list.files(fig_path))
  }
}
expected_tsv_count <- 3L
actual_tsv_count <- sum(file.exists(c(out_ew, out_dw, out_st)))
if (fig_file_count == 20L && actual_tsv_count == expected_tsv_count) {
  cat(sprintf("  [PASS] %d figure files + %d TSVs\n", fig_file_count, actual_tsv_count))
  checks_passed <- checks_passed + 1L
} else {
  cat(sprintf("  [WARN] Expected 20 fig files + 3 TSVs, got %d + %d\n",
              fig_file_count, actual_tsv_count))
}

# Check 6: Wilcoxon p-values sensible
checks_total <- checks_total + 1L
pvals <- all_stats$p_value[!is.na(all_stats$p_value)]
if (length(pvals) > 0 && any(pvals < 0.05) && !all(pvals == 0)) {
  cat(sprintf("  [PASS] Wilcoxon: %d tests, %d significant (padj < 0.05)\n",
              nrow(all_stats), sum(all_stats$padj < 0.05, na.rm = TRUE)))
  checks_passed <- checks_passed + 1L
} else {
  cat("  [WARN] Wilcoxon p-values may be degenerate\n")
}

# -- Summary ------------------------------------------------------------------

elapsed <- (proc.time() - start_time)["elapsed"]
cat(sprintf("\n===========================================\n"))
cat(sprintf("C2 COMPLETE: %d/%d checks passed\n", checks_passed, checks_total))
cat(sprintf("Runtime: %.1f seconds\n", elapsed))
cat(sprintf("Outputs: 3 TSVs + 5 figure sets (20 files) in %s\n", OUT_DIR))
cat(sprintf("Finished: %s\n", date()))
cat("===========================================\n")
