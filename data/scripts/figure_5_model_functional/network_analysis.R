# scripts/network_analysis.R
#
# Gene-Level Multi-Layer Structural Disruption Network (Figure 5C)
#
# Integrates TAD boundary, chromatin loop, and ABC enhancer-gene layers
# to build a gene-level network of structurally disrupted genes in BAP1-KO.
#
# Usage:
#   Rscript scripts/network_analysis.R                       # Late timepoint (default)
#   Rscript scripts/network_analysis.R --timepoint late      # Explicit late
#   Rscript scripts/network_analysis.R --min-layers 3        # Require all 3 layers
#   Rscript scripts/network_analysis.R --max-nodes 200       # Cap nodes for readability
#
# Output:
#   output/network_analysis/{timepoint}/
#     tables/  - Gene profile TSV, edge list TSV, centrality, GO, summary
#     plots/   - Network plot, companion diagnostic plots (PDF/SVG/JPG)

# ==============================================================================
# 1. PACKAGE LOADING
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggforce)
  library(ggnewscale)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(readxl)
})

source("data/scripts/_shared/multi_format_output.R") # Original: source("scripts/utils/multi_format_output.R")

# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

BASE_DIR <- getwd()

TIMEPOINT_CONFIG <- list(
  late = list(
    gene_summary      = file.path(BASE_DIR, "data/tsvs/figure_5_model_functional/5B_gene_level_summary.tsv"), # Original: abc/results/gene_level_summary.tsv
    loops_file        = file.path(BASE_DIR, "data/tsvs/figure_5_model_functional/5C_loops_with_gene_assignments.tsv"), # Original: abc/results/loops_with_gene_assignments.tsv
    loops_gene_cols   = c(anchor1 = "anchor1_gene", anchor2 = "anchor2_gene"),
    boundary_genes    = file.path(BASE_DIR, "data/tsvs/figure_5_model_functional/5A_boundary_genes.tsv"), # Original: tads/results/visualizations/late/enrichment/boundary_genes.tsv
    delta_abc_pairs   = file.path(BASE_DIR, "data/tsvs/figure_4_abc_analysis/4A_delta_abc_with_rnaseq.tsv"), # Original: abc/results/delta_abc_with_rnaseq.tsv
    rnaseq_file       = NULL,
    output_dir        = file.path(BASE_DIR, "data/plots/figure_5_model_functional"), # Original: output/network_analysis/late
    tsv_dir           = file.path(BASE_DIR, "data/tsvs/figure_5_model_functional"), # Original: output/network_analysis/late/tables
    label             = "Late Timepoint (Adult)",
    has_abc           = TRUE
  ),
  early = list(
    gene_summary      = NULL,
    loops_file        = file.path(BASE_DIR, "data/upstream/loop_calls/early_characterized_loops.tsv"), # Original: outputs/250831-early_outputs/merged_loops/characterized_loops.tsv
    loops_gene_cols   = c(anchor1 = "anchor1_nearest_gene", anchor2 = "anchor2_nearest_gene"),
    boundary_genes    = file.path(BASE_DIR, "tads/results/visualizations/early/enrichment/boundary_genes.tsv"), # Note: repo-relative path, not bundled in data/
    delta_abc_pairs   = NULL,
    rnaseq_file       = file.path(BASE_DIR, "data/upstream/rna_seq/young_rnaseq_results.xlsx"), # Original: tads/young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx
    output_dir        = file.path(BASE_DIR, "data/plots/figure_5_model_functional"), # Original: output/network_analysis/early
    tsv_dir           = file.path(BASE_DIR, "data/tsvs/figure_5_model_functional"), # Original: output/network_analysis/early/tables
    label             = "Early Timepoint (P13)",
    has_abc           = FALSE
  )
)

THRESHOLDS <- list(
  min_layers           = 2L,
  max_nodes            = 150L,
  abc_delta_min        = 0.05,
  baseMean_min         = 10,
  shared_enh_delta_min = 0.005,
  # NOTE: GO enrichment thresholds unused since switching to pre-defined

  # GO groups from Go_term_selction.xlsx. Retained for reference/rollback.
  go_pvalue_cutoff     = 0.05,
  go_qvalue_cutoff     = 0.1,
  top_go_terms         = 10L,
  min_go_genes         = 3L,
  label_top_n          = 10L
)

# Pre-defined GO grouping from curated gene sets
GO_GROUPS_CONFIG <- list(
  excel_path = file.path(BASE_DIR, "data/upstream/Go_term_selction.xlsx"), # Original: peaks/Go_term_selction.xlsx
  sheet      = "Sheet1",
  skip_rows  = 2L,
  categories = list(
    "Cell Fate Commitment" = list(col_idx = 5L),
    "Synapse Assembly"     = list(col_idx = 6L),
    "Axon Guidance"        = list(col_idx = 7L)
  )
)

COLORS <- list(
  expression_low  = "#4575b4",
  expression_mid  = "white",
  expression_high = "#d73027",
  edge_color      = "grey60",
  cluster_alpha   = 0.10
)

THEME_PUB <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10)
  )

# ==============================================================================
# 3. ARGUMENT PARSING
# ==============================================================================

parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  params <- list(timepoint = "late", min_layers = NULL, max_nodes = NULL)

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--timepoint" && i < length(args)) {
      params$timepoint <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--min-layers" && i < length(args)) {
      params$min_layers <- as.integer(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--max-nodes" && i < length(args)) {
      params$max_nodes <- as.integer(args[i + 1])
      i <- i + 2
    } else if (args[i] %in% c("--help", "-h")) {
      cat("Usage: Rscript scripts/network_analysis.R [OPTIONS]\n")
      cat("  --timepoint  late|early  (default: late)\n")
      cat("  --min-layers N           (default: 2)\n")
      cat("  --max-nodes  N           (default: 150)\n")
      quit(save = "no", status = 0)
    } else {
      i <- i + 1
    }
  }

  if (!params$timepoint %in% names(TIMEPOINT_CONFIG)) {
    stop(sprintf("Invalid timepoint '%s'. Must be: %s",
                 params$timepoint, paste(names(TIMEPOINT_CONFIG), collapse = ", ")))
  }

  if (!is.null(params$min_layers)) THRESHOLDS$min_layers <<- params$min_layers
  if (!is.null(params$max_nodes))  THRESHOLDS$max_nodes  <<- params$max_nodes

  return(params)
}

# ==============================================================================
# 4. DATA LOADING FUNCTIONS
# ==============================================================================

load_gene_summary <- function(path) {
  stopifnot(file.exists(path))
  df <- read_tsv(path, show_col_types = FALSE)
  required <- c("TargetGene", "max_delta_unnorm", "sum_delta_unnorm",
                 "log2FC", "padj", "baseMean", "dysregulated")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) stop(sprintf("gene_level_summary missing columns: %s", paste(missing, collapse = ", ")))
  cat(sprintf("  Gene summary: %d genes loaded\n", nrow(df)))
  df
}

load_loop_gene_pairs <- function(path, gene_cols) {
  stopifnot(file.exists(path))
  df <- read_tsv(path, show_col_types = FALSE)
  required <- c("loop_id", gene_cols[["anchor1"]], gene_cols[["anchor2"]], "logFC", "FDR", "direction")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) stop(sprintf("Loop file missing columns: %s", paste(missing, collapse = ", ")))
  # Standardize gene column names
  df <- df %>% dplyr::rename(anchor1_gene = !!gene_cols[["anchor1"]],
                              anchor2_gene = !!gene_cols[["anchor2"]])
  # Convert Entrez IDs to symbols if gene columns are numeric
  if (is.numeric(df$anchor1_gene)) {
    cat("  Converting Entrez IDs to gene symbols...\n")
    all_ids <- unique(c(as.character(df$anchor1_gene), as.character(df$anchor2_gene)))
    id_map <- AnnotationDbi::select(org.Mm.eg.db, keys = all_ids,
                                     columns = "SYMBOL", keytype = "ENTREZID")
    id_map <- id_map %>% filter(!is.na(SYMBOL)) %>% distinct(ENTREZID, .keep_all = TRUE)
    df <- df %>%
      mutate(anchor1_gene = as.character(anchor1_gene),
             anchor2_gene = as.character(anchor2_gene)) %>%
      left_join(id_map, by = c("anchor1_gene" = "ENTREZID")) %>%
      dplyr::rename(anchor1_symbol = SYMBOL) %>%
      left_join(id_map, by = c("anchor2_gene" = "ENTREZID")) %>%
      dplyr::rename(anchor2_symbol = SYMBOL) %>%
      mutate(anchor1_gene = coalesce(anchor1_symbol, anchor1_gene),
             anchor2_gene = coalesce(anchor2_symbol, anchor2_gene)) %>%
      dplyr::select(-anchor1_symbol, -anchor2_symbol)
    cat(sprintf("    Mapped %d/%d Entrez IDs to symbols\n", nrow(id_map), length(all_ids)))
  }
  cat(sprintf("  Loop-gene pairs: %d rows (%d unique loops)\n", nrow(df), n_distinct(df$loop_id)))
  df
}

load_rnaseq_excel <- function(path) {
  stopifnot(file.exists(path))
  df <- read_excel(path, sheet = "Output")
  df <- df %>%
    transmute(
      TargetGene = ensembl_gene_id,
      log2FC     = log2FoldChange,
      padj       = padj,
      baseMean   = baseMean
    ) %>%
    filter(!is.na(TargetGene))
  cat(sprintf("  RNA-seq (Excel): %d genes loaded\n", nrow(df)))
  df
}

load_boundary_genes <- function(path) {
  stopifnot(file.exists(path))
  df <- read_tsv(path, show_col_types = FALSE)
  required <- c("symbol", "enriched_in")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) stop(sprintf("boundary_genes missing columns: %s", paste(missing, collapse = ", ")))
  cat(sprintf("  Boundary genes: %d rows (%d unique symbols)\n", nrow(df), n_distinct(df$symbol)))
  df
}

load_delta_abc_pairs <- function(path) {
  stopifnot(file.exists(path))
  df <- read_tsv(path, show_col_types = FALSE)
  required <- c("chr", "start", "end", "TargetGene", "delta_unnorm")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) stop(sprintf("delta_abc_pairs missing columns: %s", paste(missing, collapse = ", ")))
  cat(sprintf("  Delta ABC pairs: %d enhancer-gene pairs\n", nrow(df)))
  df
}

# ==============================================================================
# 5. BUILD PER-GENE STRUCTURAL PROFILE
# ==============================================================================

build_boundary_layer <- function(boundary_genes) {
  boundary_genes %>%
    distinct(symbol, .keep_all = TRUE) %>%
    transmute(
      gene = symbol,
      has_boundary = TRUE,
      boundary_direction = enriched_in
    )
}

build_loop_layer <- function(loops) {
  # Deduplicate loops first (same loop_id appears multiple times for multiple gene assignments)
  unique_loops <- loops %>%
    distinct(loop_id, .keep_all = TRUE)

  bind_rows(
    unique_loops %>% transmute(gene = anchor1_gene, loop_logFC = logFC, loop_direction = direction),
    unique_loops %>% transmute(gene = anchor2_gene, loop_logFC = logFC, loop_direction = direction)
  ) %>%
    group_by(gene) %>%
    summarise(
      n_diff_loops       = n(),
      n_lost_loops       = sum(loop_direction == "down_in_mutant"),
      n_gained_loops     = sum(loop_direction == "up_in_mutant"),
      max_loop_logFC_abs = max(abs(loop_logFC)),
      mean_loop_logFC    = mean(loop_logFC),
      .groups = "drop"
    )
}

build_abc_layer <- function(gene_summary, abc_delta_min) {
  gene_summary %>%
    filter(abs(max_delta_unnorm) >= abc_delta_min) %>%
    transmute(
      gene = TargetGene,
      has_abc = TRUE
    )
}

safe_zscore <- function(x) {
  x_clean <- replace_na(x, 0)
  s <- sd(x_clean, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x_clean - mean(x_clean, na.rm = TRUE)) / s
}

build_gene_structural_profile <- function(gene_summary, loops, boundary_genes, abc_delta_min) {
  cat("  Building boundary layer...\n")
  boundary_layer <- build_boundary_layer(boundary_genes)

  cat("  Building loop layer...\n")
  loop_layer <- build_loop_layer(loops)

  cat("  Building ABC layer...\n")
  abc_layer <- build_abc_layer(gene_summary, abc_delta_min)

  cat("  Joining layers...\n")
  profile <- gene_summary %>%
    transmute(
      gene            = TargetGene,
      log2FC          = log2FC,
      padj            = padj,
      baseMean        = baseMean,
      dysregulated    = dysregulated,
      max_delta_unnorm = max_delta_unnorm,
      sum_delta_unnorm = sum_delta_unnorm,
      n_enhancers     = n_enhancers,
      n_gained_enh    = n_gained,
      n_lost_enh      = n_lost
    ) %>%
    left_join(boundary_layer, by = "gene") %>%
    left_join(loop_layer, by = "gene") %>%
    left_join(abc_layer, by = "gene") %>%
    mutate(
      has_boundary = replace_na(has_boundary, FALSE),
      has_abc      = replace_na(has_abc, FALSE),
      has_loops    = !is.na(n_diff_loops) & n_diff_loops > 0,
      n_layers     = as.integer(has_boundary) + as.integer(has_loops) + as.integer(has_abc)
    )

  cat("  Computing combined structural score...\n")
  profile <- profile %>%
    mutate(
      z_boundary      = ifelse(has_boundary, 1, 0),
      z_loop_count    = safe_zscore(replace_na(n_diff_loops, 0)),
      z_loop_strength = safe_zscore(replace_na(max_loop_logFC_abs, 0)),
      z_abc           = safe_zscore(replace_na(abs(max_delta_unnorm), 0)),
      combined_score  = z_boundary + z_loop_count + z_loop_strength + z_abc
    )

  cat(sprintf("  Layer distribution:\n"))
  layer_counts <- table(profile$n_layers)
  for (l in names(layer_counts)) {
    cat(sprintf("    %s layers: %d genes\n", l, layer_counts[l]))
  }

  profile
}

# ==============================================================================
# 6. FILTER GENES FOR NETWORK
# ==============================================================================

filter_network_genes <- function(profile, thresholds) {
  filtered <- profile %>%
    filter(n_layers >= thresholds$min_layers) %>%
    filter(baseMean >= thresholds$baseMean_min)

  cat(sprintf("  After min_layers>=%d and baseMean>=%d: %d genes\n",
              thresholds$min_layers, thresholds$baseMean_min, nrow(filtered)))

  if (nrow(filtered) > thresholds$max_nodes) {
    cat(sprintf("  Capping to top %d by combined_score\n", thresholds$max_nodes))
    filtered <- filtered %>%
      arrange(desc(combined_score)) %>%
      slice_head(n = thresholds$max_nodes)
  }

  cat(sprintf("  Final network nodes: %d genes\n", nrow(filtered)))
  filtered
}

# ==============================================================================
# 7. BUILD EDGE LIST
# ==============================================================================

build_loop_edges <- function(filtered_gene_set, loops) {
  # Deduplicate to unique loops
  unique_loops <- loops %>% distinct(loop_id, .keep_all = TRUE)

  edges <- unique_loops %>%
    filter(
      anchor1_gene %in% filtered_gene_set,
      anchor2_gene %in% filtered_gene_set,
      anchor1_gene != anchor2_gene
    ) %>%
    transmute(
      from = pmin(anchor1_gene, anchor2_gene),
      to   = pmax(anchor1_gene, anchor2_gene),
      weight = abs(logFC),
      edge_type = "loop",
      loop_direction = direction
    ) %>%
    group_by(from, to) %>%
    summarise(
      weight = max(weight),
      n_connecting_loops = n(),
      edge_type = "loop",
      .groups = "drop"
    )

  cat(sprintf("  Loop edges: %d (from %d loops connecting filtered genes)\n",
              nrow(edges), sum(edges$n_connecting_loops)))
  edges
}

build_shared_enhancer_edges <- function(filtered_gene_set, delta_abc_pairs, delta_min) {
  # Filter to enhancers linked to filtered genes with sufficient delta
  enh_filtered <- delta_abc_pairs %>%
    filter(
      TargetGene %in% filtered_gene_set,
      abs(delta_unnorm) >= delta_min
    ) %>%
    mutate(enhancer_id = paste0(chr, ":", start, "-", end))

  # Find enhancers linked to 2+ filtered genes
  multi_gene_enh <- enh_filtered %>%
    group_by(enhancer_id) %>%
    filter(n_distinct(TargetGene) >= 2) %>%
    ungroup()

  if (nrow(multi_gene_enh) == 0) {
    cat("  Shared enhancer edges: 0\n")
    return(tibble(from = character(), to = character(), weight = numeric(),
                  n_shared_enhancers = integer(), edge_type = character()))
  }

  # Build pairwise edges
  gene_per_enh <- multi_gene_enh %>%
    distinct(enhancer_id, TargetGene, .keep_all = TRUE) %>%
    dplyr::select(enhancer_id, TargetGene, delta_unnorm)

  edges <- gene_per_enh %>%
    inner_join(gene_per_enh, by = "enhancer_id", suffix = c("_from", "_to"),
               relationship = "many-to-many") %>%
    filter(TargetGene_from < TargetGene_to) %>%
    transmute(
      from = TargetGene_from,
      to = TargetGene_to,
      enh_weight = (abs(delta_unnorm_from) + abs(delta_unnorm_to)) / 2,
      edge_type = "shared_enhancer"
    ) %>%
    group_by(from, to) %>%
    summarise(
      weight = max(enh_weight),
      n_shared_enhancers = n(),
      edge_type = "shared_enhancer",
      .groups = "drop"
    )

  cat(sprintf("  Shared enhancer edges: %d (from %d multi-gene enhancers)\n",
              nrow(edges), n_distinct(multi_gene_enh$enhancer_id)))
  edges
}

# ==============================================================================
# 8. PRE-DEFINED GO-BASED NODE GROUPING
# ==============================================================================

assign_go_groups_from_excel <- function(gene_symbols, config = GO_GROUPS_CONFIG) {
  stopifnot(file.exists(config$excel_path))

  cat(sprintf("  Loading curated GO gene sets from: %s\n", basename(config$excel_path)))

  # Read sheet, skipping merged header rows to reach data
  raw <- read_excel(config$excel_path, sheet = config$sheet,
                    skip = config$skip_rows, col_names = FALSE)

  # Extract gene lists from configured columns
  gene_sets <- list()
  for (category_name in names(config$categories)) {
    col_idx <- config$categories[[category_name]]$col_idx
    genes <- raw[[col_idx]]
    genes <- genes[!is.na(genes)]
    genes <- trimws(as.character(genes))
    genes <- genes[!grepl("^Top sets", genes, ignore.case = TRUE)]
    gene_sets[[category_name]] <- genes
    cat(sprintf("    %s: %d genes in curated set\n", category_name, length(genes)))
  }

  # Check for overlaps (defensive -- currently none exist)
  all_genes_flat <- unlist(gene_sets)
  if (any(duplicated(all_genes_flat))) {
    dup_genes <- unique(all_genes_flat[duplicated(all_genes_flat)])
    cat(sprintf("  WARNING: %d genes in multiple categories: %s\n",
                length(dup_genes), paste(dup_genes, collapse = ", ")))
    cat("  Priority assignment: first listed category wins\n")
  }

  # Build gene-to-category mapping (first match wins for overlaps)
  gene_to_group <- tibble(gene = character(), go_group = character())
  assigned <- character()
  for (category_name in names(gene_sets)) {
    new_genes <- setdiff(gene_sets[[category_name]], assigned)
    gene_to_group <- bind_rows(
      gene_to_group,
      tibble(gene = new_genes, go_group = category_name)
    )
    assigned <- c(assigned, new_genes)
  }

  # Join with network genes
  result <- tibble(gene = gene_symbols) %>%
    left_join(gene_to_group, by = "gene") %>%
    mutate(go_group = replace_na(go_group, "Other"))

  # Report overlap with network
  for (category_name in names(gene_sets)) {
    n_in_network <- sum(result$go_group == category_name)
    n_in_set <- length(gene_sets[[category_name]])
    cat(sprintf("    %s: %d/%d genes found in network\n",
                category_name, n_in_network, n_in_set))
  }
  cat(sprintf("  Assigned %d genes to GO groups (%d to 'Other')\n",
              sum(result$go_group != "Other"), sum(result$go_group == "Other")))

  list(groups = result, ego = NULL, gene_sets = gene_sets)
}

# ==============================================================================
# 9. NETWORK CONSTRUCTION
# ==============================================================================

build_network <- function(nodes, edges) {
  genes_with_edges <- unique(c(edges$from, edges$to))
  n_isolated <- sum(!nodes$gene %in% genes_with_edges)
  if (n_isolated > 0) {
    cat(sprintf("  %d genes have no edges (isolated nodes)\n", n_isolated))
  }

  # Ensure edge endpoints exist in node list
  edges_clean <- edges %>%
    filter(from %in% nodes$gene, to %in% nodes$gene)

  # Build via igraph (handles character from/to matching to vertex names)
  node_df <- nodes %>% dplyr::rename(name = gene)
  ig <- graph_from_data_frame(
    d = edges_clean,
    vertices = node_df,
    directed = FALSE
  )

  g <- as_tbl_graph(ig) %>%
    mutate(
      degree      = centrality_degree(),
      betweenness = centrality_betweenness(normalized = TRUE),
      closeness   = centrality_closeness(normalized = TRUE)
    )

  n_comp <- igraph::components(ig)$no
  cat(sprintf("  Network: %d nodes, %d edges, %d components\n",
              igraph::vcount(ig), igraph::ecount(ig), n_comp))

  g
}

# ==============================================================================
# 10A. CLUSTERED LAYOUT HELPERS
# ==============================================================================

wrap_go_label <- function(label, width = 22) {
  str_wrap(label, width = width)
}

compute_clustered_layout <- function(g, go_groups_df,
                                     base_radius = 2.5,
                                     center_spacing = 10.0,
                                     other_ring_factor = 2.0) {
  set.seed(42)
  node_df <- as_tibble(g)

  # Identify GO clusters (exclude "Other")
  cluster_names <- go_groups_df %>%
    filter(go_group != "Other") %>%
    count(go_group, name = "n_genes") %>%
    filter(n_genes >= 1) %>%
    arrange(desc(n_genes))

  n_clusters <- nrow(cluster_names)
  if (n_clusters == 0) {
    stop("No GO clusters found — cannot compute clustered layout")
  }

  # Compute per-cluster radius (area proportional to gene count)
  cluster_names <- cluster_names %>%
    mutate(
      radius = base_radius * sqrt(n_genes / pi),
      radius = pmax(radius, 1.5)
    )

  # Position cluster centers on a circle
  center_ring_radius <- center_spacing * n_clusters / (2 * pi)
  center_angles <- seq(0, 2 * pi, length.out = n_clusters + 1)[1:n_clusters]

  cluster_centers <- cluster_names %>%
    mutate(
      cx = center_ring_radius * cos(center_angles),
      cy = center_ring_radius * sin(center_angles)
    )

  # Ensure no overlap: scale outward if any cluster circles intersect
  min_gap <- 2.0
  needs_scaling <- TRUE
  max_iter <- 10
  iter <- 0
  while (needs_scaling && iter < max_iter) {
    needs_scaling <- FALSE
    iter <- iter + 1
    for (i in seq_len(nrow(cluster_centers) - 1)) {
      for (j in (i + 1):nrow(cluster_centers)) {
        dist_ij <- sqrt(
          (cluster_centers$cx[i] - cluster_centers$cx[j])^2 +
          (cluster_centers$cy[i] - cluster_centers$cy[j])^2
        )
        min_dist <- cluster_centers$radius[i] + cluster_centers$radius[j] + min_gap
        if (dist_ij < min_dist) {
          scale <- (min_dist / dist_ij) * 1.05
          cluster_centers <- cluster_centers %>%
            mutate(cx = cx * scale, cy = cy * scale)
          needs_scaling <- TRUE
          break
        }
      }
      if (needs_scaling) break
    }
  }

  # Position genes within each cluster in a circle
  positions <- list()
  for (i in seq_len(nrow(cluster_centers))) {
    grp <- cluster_centers$go_group[i]
    genes_in_group <- go_groups_df$gene[go_groups_df$go_group == grp]
    n <- length(genes_in_group)
    r <- cluster_centers$radius[i]
    cx <- cluster_centers$cx[i]
    cy <- cluster_centers$cy[i]

    # Order genes by combined_score for consistent placement
    gene_order <- node_df %>%
      filter(name %in% genes_in_group) %>%
      arrange(desc(combined_score)) %>%
      pull(name)

    # Handle any genes not in node_df (shouldn't happen, but defensive)
    gene_order <- c(gene_order, setdiff(genes_in_group, gene_order))

    if (n == 1) {
      positions[[grp]] <- tibble(name = gene_order, x = cx, y = cy)
    } else {
      angles <- seq(0, 2 * pi, length.out = n + 1)[1:n]
      positions[[grp]] <- tibble(
        name = gene_order,
        x = cx + r * 0.80 * cos(angles),
        y = cy + r * 0.80 * sin(angles)
      )
    }
  }

  # Position "Other" genes on an outer ring
  other_genes <- go_groups_df$gene[go_groups_df$go_group == "Other"]
  if (length(other_genes) > 0) {
    outer_radius <- max(sqrt(cluster_centers$cx^2 + cluster_centers$cy^2)) +
                    max(cluster_centers$radius) * other_ring_factor + 3
    # Sort "Other" genes by combined_score
    other_order <- node_df %>%
      filter(name %in% other_genes) %>%
      arrange(desc(combined_score)) %>%
      pull(name)
    other_order <- c(other_order, setdiff(other_genes, other_order))

    other_angles <- seq(0, 2 * pi, length.out = length(other_order) + 1)[1:length(other_order)]
    positions[["Other"]] <- tibble(
      name = other_order,
      x = outer_radius * cos(other_angles),
      y = outer_radius * sin(other_angles)
    )
  }

  all_positions <- bind_rows(positions)

  list(
    node_positions = all_positions,
    cluster_centers = cluster_centers
  )
}

# ==============================================================================
# 10B. NETWORK VISUALIZATION (FIGURE 5C)
# ==============================================================================

plot_network <- function(g, thresholds, colors, cfg_label) {
  set.seed(42)

  node_df <- as_tibble(g)
  edge_df <- as_tibble(g, "edges")

  # Edge from/to are integer indices in tidygraph — map to node names
  clustered_genes <- node_df$name[node_df$go_group != "Other"]
  # Only keep "Other" genes that connect to at least one clustered gene
  edges_touching_cluster <- edge_df %>%
    filter(node_df$name[from] %in% clustered_genes |
           node_df$name[to] %in% clustered_genes)
  other_genes_near_cluster <- setdiff(
    unique(c(node_df$name[edges_touching_cluster$from],
             node_df$name[edges_touching_cluster$to])),
    clustered_genes
  )
  show_node <- (node_df$go_group != "Other") |
    (node_df$name %in% other_genes_near_cluster)

  g_show <- g %>%
    mutate(show = show_node) %>%
    filter(show)

  n_total <- igraph::vcount(as.igraph(g))
  n_shown <- igraph::vcount(as.igraph(g_show))
  n_edges <- igraph::ecount(as.igraph(g_show))
  n_clustered <- sum(as_tibble(g_show)$go_group != "Other")
  n_peripheral <- n_shown - n_clustered
  cat(sprintf("  Plotting %d nodes (%d in clusters, %d peripheral, %d removed)\n",
              n_shown, n_clustered, n_peripheral, n_total - n_shown))

  # Compute clustered layout
  go_groups_df <- as_tibble(g_show) %>% dplyr::select(gene = name, go_group)
  layout_result <- compute_clustered_layout(g_show, go_groups_df)

  # Join positions onto graph node attributes
  g_show <- g_show %>%
    left_join(layout_result$node_positions, by = "name") %>%
    mutate(
      node_size = pmax(abs(replace_na(max_delta_unnorm, 0)), 0.005),
      is_other = (go_group == "Other"),
      display_size = ifelse(is_other, node_size * 0.5, node_size),
      display_alpha = ifelse(is_other, 0.4, 1.0),
      # Label top genes in clusters + a few "Other"
      cluster_rank = ifelse(!is_other,
                            rank(-combined_score, ties.method = "first"),
                            Inf),
      other_rank = ifelse(is_other,
                          rank(-combined_score, ties.method = "first"),
                          Inf),
      show_label = (cluster_rank <= thresholds$label_top_n) | (other_rank <= 3),
      label_text = ifelse(show_label, name, NA_character_)
    )

  # Add edge-level attribute: de-emphasize edges involving "Other" genes
  show_node_df <- as_tibble(g_show)
  g_show <- g_show %>%
    activate(edges) %>%
    mutate(
      from_clustered = show_node_df$go_group[from] != "Other",
      to_clustered   = show_node_df$go_group[to] != "Other",
      edge_alpha = ifelse(from_clustered & to_clustered, 0.40, 0.15),
      edge_width_mult = ifelse(from_clustered & to_clustered, 1.0, 0.5)
    ) %>%
    activate(nodes)

  # Create manual layout using precomputed x, y
  layout <- create_layout(g_show, layout = "manual", x = x, y = y)

  # Cluster background data
  cluster_centers <- layout_result$cluster_centers
  circle_data <- cluster_centers %>%
    mutate(
      label = sapply(go_group, wrap_go_label, width = 22),
      label_y = cy + radius + 1.2
    )

  # Get Set2 colors for clusters
  n_clusters <- nrow(circle_data)
  cluster_colors <- scales::brewer_pal(palette = "Set2")(max(3, n_clusters))[1:n_clusters]
  names(cluster_colors) <- circle_data$go_group

  # --- Build plot layer by layer ---
  p <- ggraph(layout) +
    scale_x_continuous(expand = expansion(mult = 0.12)) +
    scale_y_continuous(expand = expansion(mult = 0.12))

  # Layer 1: Cluster background circles
  if (nrow(circle_data) > 0) {
    p <- p +
      ggforce::geom_circle(
        data = circle_data,
        aes(x0 = cx, y0 = cy, r = radius, fill = go_group),
        alpha = colors$cluster_alpha,
        color = "grey70",
        linewidth = 0.4,
        inherit.aes = FALSE
      ) +
      scale_fill_manual(
        values = cluster_colors,
        name = "GO Group",
        guide = guide_legend(order = 1)
      )
  }

  # Layer 2: Cluster labels above circles
  if (nrow(circle_data) > 0) {
    p <- p +
      geom_text(
        data = circle_data,
        aes(x = cx, y = label_y, label = label),
        size = 2.8,
        fontface = "bold",
        color = "grey25",
        lineheight = 0.9,
        inherit.aes = FALSE
      )
  }

  # Layer 3: Edges (uniform color, de-emphasize edges to "Other" genes)
  p <- p +
    geom_edge_link(
      aes(edge_width = weight, edge_alpha = edge_alpha),
      edge_colour = colors$edge_color,
      show.legend = FALSE
    ) +
    scale_edge_width(range = c(0.2, 1.5), guide = "none") +
    scale_edge_alpha_identity()

  # Layer 4: Nodes as ggforce circles (fill=log2FC, linetype=n_layers)
  # Compute radii in data coordinates
  node_df <- as.data.frame(layout)
  x_extent <- diff(range(node_df$x, na.rm = TRUE))
  y_extent <- diff(range(node_df$y, na.rm = TRUE))
  eff_w_mm <- 14 * 25.4 * (1 - 2 * 0.12)
  eff_h_mm <- 12 * 25.4 * (1 - 2 * 0.12)
  mm_per_data <- min(eff_w_mm / x_extent, eff_h_mm / y_extent)

  node_circles <- node_df %>%
    transmute(
      x0 = x, y0 = y,
      size_mm = scales::rescale(display_size, to = c(2, 12)),
      r = (size_mm / 2) / mm_per_data * 1.15,
      log2FC = log2FC,
      border_lty = factor(
        case_when(
          n_layers >= 3 ~ "3 layers",
          n_layers == 2 ~ "2 layers",
          TRUE ~ "1 layer"
        ),
        levels = c("1 layer", "2 layers", "3 layers")
      ),
      node_alpha = display_alpha,
      display_size = display_size
    )

  p <- p +
    ggnewscale::new_scale_fill() +
    ggforce::geom_circle(
      data = node_circles,
      aes(x0 = x0, y0 = y0, r = r, fill = log2FC,
          linetype = border_lty, alpha = node_alpha),
      colour = "black", linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    scale_fill_gradient2(
      low = colors$expression_low,
      mid = colors$expression_mid,
      high = colors$expression_high,
      midpoint = 0,
      name = "RNA-seq\nlog2FC",
      limits = c(-3, 3),
      oob = squish
    ) +
    scale_linetype_manual(
      values = c("1 layer" = "dotted", "2 layers" = "dashed", "3 layers" = "solid"),
      name = "Structural\nLayers",
      guide = guide_legend(order = 4)
    ) +
    scale_alpha_identity() +
    # Hidden point layer for the size legend
    geom_point(
      data = node_circles,
      aes(x = x0, y = y0, size = display_size),
      alpha = 0, inherit.aes = FALSE
    ) +
    scale_size_continuous(
      range = c(2, 12), name = "|max delta AxC|",
      guide = guide_legend(
        order = 5,
        override.aes = list(alpha = 1, shape = 21, fill = "white",
                            colour = "black", stroke = 0.5)
      )
    )

  # Layer 5: Gene labels
  p <- p +
    geom_node_text(
      aes(label = label_text),
      repel = TRUE,
      size = 2.5,
      max.overlaps = 30,
      segment.size = 0.2,
      segment.alpha = 0.4,
      na.rm = TRUE
    )

  # Layer 6: Theme and titles
  p <- p +
    theme_graph(base_family = "") +
    labs(
      title = sprintf("Multi-Layer Structural Disruption Network (BAP1-KO, %s)", cfg_label),
      subtitle = sprintf("%d genes in %d GO clusters + %d peripheral | %d edges",
                         n_clustered, n_clusters, n_peripheral, n_edges)
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10),
      legend.box = "vertical"
    )

  p
}

# ==============================================================================
# 11. COMPANION PLOTS
# ==============================================================================

plot_layer_distribution <- function(profile) {
  layer_df <- profile %>%
    count(n_layers) %>%
    mutate(
      pct = n / sum(n) * 100,
      label = sprintf("%d\n(%.1f%%)", n, pct)
    )

  ggplot(layer_df, aes(x = factor(n_layers), y = n, fill = factor(n_layers))) +
    geom_col(width = 0.7) +
    geom_text(aes(label = label), vjust = -0.3, size = 3.5) +
    scale_fill_manual(
      values = c("0" = "#d9d9d9", "1" = "#bdbdbd", "2" = "#fdae6b", "3" = "#e6550d"),
      guide = "none"
    ) +
    labs(
      title = "Gene Disruption Layer Distribution",
      subtitle = sprintf("Total: %d genes", nrow(profile)),
      x = "Number of Structural Layers Disrupted",
      y = "Number of Genes"
    ) +
    THEME_PUB +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
}

plot_combined_score_distribution <- function(filtered) {
  ggplot(filtered, aes(x = combined_score)) +
    geom_histogram(bins = 30, fill = "#636363", color = "white", linewidth = 0.3) +
    geom_vline(xintercept = min(filtered$combined_score), linetype = "dashed", color = "#e6550d") +
    labs(
      title = "Combined Structural Score Distribution (Filtered Genes)",
      subtitle = sprintf("n = %d | range: [%.2f, %.2f]",
                         nrow(filtered), min(filtered$combined_score), max(filtered$combined_score)),
      x = "Combined Structural Score",
      y = "Count"
    ) +
    THEME_PUB
}

plot_edge_type_breakdown <- function(edges) {
  if (nrow(edges) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No edges") + theme_void())
  }

  edge_summary <- edges %>%
    count(edge_type) %>%
    mutate(pct = n / sum(n) * 100)

  ggplot(edge_summary, aes(x = edge_type, y = n, fill = edge_type)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%d (%.1f%%)", n, pct)), vjust = -0.3, size = 4) +
    scale_fill_manual(values = c("loop" = "#333333", "shared_enhancer" = "#FF7F00"), guide = "none") +
    labs(
      title = "Edge Type Breakdown",
      x = "Edge Type",
      y = "Number of Edges"
    ) +
    THEME_PUB +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
}

plot_degree_distribution <- function(g) {
  deg_df <- as_tibble(g) %>% dplyr::select(name, degree)

  ggplot(deg_df, aes(x = degree)) +
    geom_histogram(binwidth = 1, fill = "#636363", color = "white", linewidth = 0.3) +
    labs(
      title = "Node Degree Distribution",
      subtitle = sprintf("Mean degree: %.1f | Max degree: %d",
                         mean(deg_df$degree), max(deg_df$degree)),
      x = "Degree (Number of Connections)",
      y = "Count"
    ) +
    THEME_PUB
}

plot_go_group_summary <- function(go_groups) {
  if (is.null(go_groups) || nrow(go_groups) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
           label = "No GO group data") + theme_void())
  }

  group_counts <- go_groups %>%
    count(go_group) %>%
    mutate(is_other = (go_group == "Other"))

  # Order by count but keep "Other" last
  category_order <- group_counts %>%
    filter(!is_other) %>%
    arrange(desc(n)) %>%
    pull(go_group)
  category_order <- c(category_order, "Other")
  group_counts <- group_counts %>%
    mutate(go_group = factor(go_group, levels = category_order))

  n_categories <- sum(!group_counts$is_other)
  fill_colors <- c(
    scales::brewer_pal(palette = "Set2")(max(3, n_categories))[seq_len(n_categories)],
    "grey80"
  )
  names(fill_colors) <- levels(group_counts$go_group)

  n_assigned <- sum(group_counts$n[!group_counts$is_other])
  n_other <- group_counts$n[group_counts$is_other]

  ggplot(group_counts, aes(x = go_group, y = n, fill = go_group)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n), vjust = -0.3, size = 4) +
    scale_fill_manual(values = fill_colors, guide = "none") +
    labs(
      title = "GO Category Membership (Pre-defined Gene Sets)",
      subtitle = sprintf("%d genes assigned to categories | %d in 'Other'",
                         n_assigned, n_other),
      x = NULL,
      y = "Number of Network Genes"
    ) +
    THEME_PUB +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme(axis.text.x = element_text(angle = 15, hjust = 1))
}

# ==============================================================================
# 12. SUMMARY STATISTICS
# ==============================================================================

write_summary_statistics <- function(profile, filtered, edges, g, output_dir) {
  ig <- as.igraph(g)
  node_df <- as_tibble(g)
  comp <- igraph::components(ig)

  hub_genes <- node_df %>%
    arrange(desc(degree)) %>%
    slice_head(n = 10) %>%
    dplyr::select(name, degree, betweenness, log2FC, combined_score, n_layers)

  loop_edges_n <- sum(edges$edge_type == "loop")
  enh_edges_n  <- sum(edges$edge_type == "shared_enhancer")

  summary_text <- paste0(
    "===== Network Analysis Summary =====\n",
    sprintf("Date: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),

    "--- Gene Universe ---\n",
    sprintf("Total genes in ABC summary: %d\n", nrow(profile)),
    sprintf("Genes near diff TAD boundary: %d\n", sum(profile$has_boundary)),
    sprintf("Genes at loop anchors: %d\n", sum(profile$has_loops)),
    sprintf("Genes with ABC changes (|delta_unnorm|>=%.3f): %d\n", THRESHOLDS$abc_delta_min, sum(profile$has_abc)),
    sprintf("Genes in 1 layer: %d\n", sum(profile$n_layers == 1)),
    sprintf("Genes in 2 layers: %d\n", sum(profile$n_layers == 2)),
    sprintf("Genes in 3 layers: %d\n\n", sum(profile$n_layers == 3)),

    "--- Network Nodes ---\n",
    sprintf("Filtered nodes: %d\n", nrow(filtered)),
    sprintf("Min combined score: %.3f\n", min(filtered$combined_score)),
    sprintf("Max combined score: %.3f\n", max(filtered$combined_score)),
    sprintf("Upregulated (log2FC > 0): %d\n", sum(filtered$log2FC > 0, na.rm = TRUE)),
    sprintf("Downregulated (log2FC < 0): %d\n", sum(filtered$log2FC < 0, na.rm = TRUE)),
    sprintf("DEGs (padj < 0.05): %d\n\n", sum(filtered$padj < 0.05, na.rm = TRUE)),

    "--- GO Clusters (Pre-defined from Go_term_selction.xlsx) ---\n",
    sprintf("Clustered genes: %d\n", sum(node_df$go_group != "Other")),
    sprintf("Peripheral genes (Other): %d\n", sum(node_df$go_group == "Other")),
    sprintf("Number of GO clusters: %d\n",
            n_distinct(node_df$go_group[node_df$go_group != "Other"])),
    paste0(
      paste(
        sapply(
          sort(unique(node_df$go_group[node_df$go_group != "Other"])),
          function(grp) sprintf("  - %s: %d genes", grp,
                                sum(node_df$go_group == grp))
        ),
        collapse = "\n"
      ),
      "\n\n"
    ),

    "--- Network Edges ---\n",
    sprintf("Loop edges: %d\n", loop_edges_n),
    sprintf("Shared enhancer edges: %d\n", enh_edges_n),
    sprintf("Total edges: %d\n", nrow(edges)),
    sprintf("Graph density: %.4f\n", edge_density(ig)),
    sprintf("Connected components: %d\n", comp$no),
    sprintf("Largest component: %d nodes\n", max(comp$csize)),
    sprintf("Mean degree: %.2f\n", mean(degree(ig))),
    sprintf("Max degree: %d\n\n", max(degree(ig))),

    "--- Top Hub Genes (by degree) ---\n"
  )

  for (i in seq_len(nrow(hub_genes))) {
    row <- hub_genes[i, ]
    summary_text <- paste0(summary_text, sprintf(
      "  %2d. %-15s degree=%d  betweenness=%.3f  log2FC=%+.3f  score=%.2f  layers=%d\n",
      i, row$name, row$degree, row$betweenness,
      ifelse(is.na(row$log2FC), 0, row$log2FC),
      row$combined_score, row$n_layers
    ))
  }

  summary_text <- paste0(summary_text, sprintf(
    "\n--- Parameters ---\nmin_layers: %d\nmax_nodes: %d\nabc_delta_min: %.3f\nbaseMean_min: %d\nshared_enh_delta_min: %.3f\n",
    THRESHOLDS$min_layers, THRESHOLDS$max_nodes, THRESHOLDS$abc_delta_min,
    THRESHOLDS$baseMean_min, THRESHOLDS$shared_enh_delta_min
  ))

  write(summary_text, file = file.path(output_dir, "5C_summary_report.txt")) # Original: file.path(output_dir, "tables", "summary_report.txt")
  cat(summary_text)
}

# ==============================================================================
# 13. MAIN PIPELINE
# ==============================================================================

main <- function() {
  cat("\n================================================\n")
  cat("Gene-Level Network Analysis (Figure 5C)\n")
  cat("================================================\n\n")
  cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

  # --- Parse arguments ---
  params <- parse_arguments()
  cfg <- TIMEPOINT_CONFIG[[params$timepoint]]
  cat(sprintf("Timepoint: %s\n", cfg$label))
  cat(sprintf("Min layers: %d | Max nodes: %d\n\n", THRESHOLDS$min_layers, THRESHOLDS$max_nodes))

  # --- Create output directories ---
  tables_dir <- cfg$tsv_dir # Original: file.path(cfg$output_dir, "tables")
  plots_dir  <- cfg$output_dir # Original: file.path(cfg$output_dir, "plots")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir,  recursive = TRUE, showWarnings = FALSE)

  # --- Validate inputs ---
  cat("Validating input files...\n")
  required_files <- list(
    loops_file     = cfg$loops_file,
    boundary_genes = cfg$boundary_genes
  )
  if (cfg$has_abc) {
    required_files$gene_summary    <- cfg$gene_summary
    required_files$delta_abc_pairs <- cfg$delta_abc_pairs
  }
  if (!is.null(cfg$rnaseq_file)) {
    required_files$rnaseq_file <- cfg$rnaseq_file
  }
  required_files$go_groups_excel <- GO_GROUPS_CONFIG$excel_path
  for (nm in names(required_files)) {
    stopifnot(file.exists(required_files[[nm]]))
    cat(sprintf("  [OK] %s\n", nm))
  }

  # ==========================================================================
  # SECTION 1: Load Data
  # ==========================================================================
  cat("\n=== Section 1: Loading input data ===\n")
  loops          <- load_loop_gene_pairs(cfg$loops_file, cfg$loops_gene_cols)
  boundary_genes <- load_boundary_genes(cfg$boundary_genes)

  if (cfg$has_abc) {
    gene_summary <- load_gene_summary(cfg$gene_summary)
    delta_abc    <- load_delta_abc_pairs(cfg$delta_abc_pairs)
  } else {
    # Build gene universe from RNA-seq Excel + loop/boundary genes
    rnaseq <- load_rnaseq_excel(cfg$rnaseq_file)
    # Create a minimal gene_summary-like table (no ABC columns)
    gene_summary <- rnaseq %>%
      mutate(
        max_delta_unnorm = NA_real_,
        sum_delta_unnorm = NA_real_,
        n_enhancers = NA_integer_,
        n_gained = NA_integer_,
        n_lost = NA_integer_,
        dysregulated = FALSE
      )
    delta_abc <- NULL
    cat("  ABC data: not available for this timepoint\n")
  }

  # ==========================================================================
  # SECTION 2: Build Gene Profile
  # ==========================================================================
  cat("\n=== Section 2: Building per-gene structural profile ===\n")
  abc_min <- if (cfg$has_abc) THRESHOLDS$abc_delta_min else Inf
  profile <- build_gene_structural_profile(gene_summary, loops, boundary_genes, abc_min)
  write_tsv(profile, file.path(tables_dir, "5C_gene_structural_profile_all.tsv")) # Original: "gene_structural_profile_all.tsv"
  cat(sprintf("  Saved: gene_structural_profile_all.tsv (%d genes)\n", nrow(profile)))

  # ==========================================================================
  # SECTION 3: Filter Network Genes
  # ==========================================================================
  cat("\n=== Section 3: Filtering genes for network ===\n")
  # For early (2 layers only), allow min_layers=1 if needed
  early_thresholds <- THRESHOLDS
  if (!cfg$has_abc && THRESHOLDS$min_layers > 2) {
    cat("  Note: max 2 layers available (no ABC). Setting min_layers=2\n")
    early_thresholds$min_layers <- 2L
  }
  filtered <- filter_network_genes(profile, early_thresholds)

  # ==========================================================================
  # SECTION 4: Expand network with all GO top-set genes
  # ==========================================================================
  cat("\n=== Section 4: Expanding network with GO top-set genes ===\n")
  go_result <- assign_go_groups_from_excel(filtered$gene)

  # Add all top-set genes from the profile that aren't already in filtered
  all_topset_genes <- unique(unlist(go_result$gene_sets))
  missing_from_network <- setdiff(all_topset_genes, filtered$gene)
  missing_in_profile <- profile %>% filter(gene %in% missing_from_network)
  cat(sprintf("  Top-set genes not in structural top-%d: %d\n",
              THRESHOLDS$max_nodes, length(missing_from_network)))
  cat(sprintf("  Of those, found in gene profile: %d\n", nrow(missing_in_profile)))

  if (nrow(missing_in_profile) > 0) {
    filtered <- bind_rows(filtered, missing_in_profile)
  }

  # Re-run GO assignment on the expanded gene list
  go_result <- assign_go_groups_from_excel(filtered$gene)
  filtered  <- filtered %>% left_join(go_result$groups, by = "gene")

  cat(sprintf("  Expanded network: %d total genes\n", nrow(filtered)))
  write_tsv(filtered, file.path(tables_dir, "5C_gene_structural_profile_filtered.tsv")) # Original: "gene_structural_profile_filtered.tsv"

  # Write GO group assignments
  go_assignment_df <- go_result$groups %>%
    filter(go_group != "Other") %>%
    arrange(go_group, gene)
  write_tsv(go_assignment_df, file.path(tables_dir, "5C_go_group_assignments.tsv")) # Original: "go_group_assignments.tsv"
  cat(sprintf("  Saved: go_group_assignments.tsv (%d assigned genes)\n", nrow(go_assignment_df)))

  # ==========================================================================
  # SECTION 5: Build Edges
  # ==========================================================================
  cat("\n=== Section 5: Building edge list ===\n")
  loop_edges <- build_loop_edges(filtered$gene, loops)

  if (!is.null(delta_abc)) {
    enh_edges <- build_shared_enhancer_edges(filtered$gene, delta_abc, THRESHOLDS$shared_enh_delta_min)
    all_edges <- bind_rows(
      loop_edges %>% dplyr::select(from, to, weight, edge_type),
      enh_edges  %>% dplyr::select(from, to, weight, edge_type)
    )
  } else {
    cat("  Shared enhancer edges: skipped (no ABC data)\n")
    all_edges <- loop_edges %>% dplyr::select(from, to, weight, edge_type)
  }
  write_tsv(all_edges, file.path(tables_dir, "5C_edge_list.tsv")) # Original: "edge_list.tsv"
  cat(sprintf("  Total edges: %d\n", nrow(all_edges)))

  # ==========================================================================
  # SECTION 6: Build Network
  # ==========================================================================
  cat("\n=== Section 6: Constructing network ===\n")
  g <- build_network(filtered, all_edges)

  centrality <- as_tibble(g) %>% dplyr::select(name, degree, betweenness, closeness)
  write_tsv(centrality, file.path(tables_dir, "5C_network_centrality_metrics.tsv")) # Original: "network_centrality_metrics.tsv"

  # ==========================================================================
  # SECTION 7: Visualize
  # ==========================================================================
  cat("\n=== Section 7: Generating plots ===\n")

  cat("  Rendering network figure...\n")
  p_network <- plot_network(g, THRESHOLDS, COLORS, cfg$label)
  save_multiformat_ggplot(p_network, file.path(plots_dir, "5C_network_figure"), # Original: "network_figure5c"
                          width = 14, height = 12)

  cat("  Rendering layer distribution...\n")
  p_layers <- plot_layer_distribution(profile)
  save_multiformat_ggplot(p_layers, file.path(plots_dir, "5C_layer_distribution"), # Original: "layer_distribution"
                          width = 8, height = 6)

  cat("  Rendering combined score distribution...\n")
  p_scores <- plot_combined_score_distribution(filtered)
  save_multiformat_ggplot(p_scores, file.path(plots_dir, "5C_combined_score_distribution"), # Original: "combined_score_distribution"
                          width = 8, height = 6)

  cat("  Rendering edge type breakdown...\n")
  p_edges <- plot_edge_type_breakdown(all_edges)
  save_multiformat_ggplot(p_edges, file.path(plots_dir, "5C_edge_type_breakdown"), # Original: "edge_type_breakdown"
                          width = 8, height = 6)

  cat("  Rendering degree distribution...\n")
  p_degree <- plot_degree_distribution(g)
  save_multiformat_ggplot(p_degree, file.path(plots_dir, "5C_degree_distribution"), # Original: "degree_distribution"
                          width = 8, height = 6)

  cat("  Rendering GO group summary...\n")
  p_go <- plot_go_group_summary(go_result$groups)
  save_multiformat_ggplot(p_go, file.path(plots_dir, "5C_go_group_summary"), # Original: "go_group_summary"
                          width = 10, height = 6)

  # ==========================================================================
  # SECTION 8: Summary
  # ==========================================================================
  cat("\n=== Section 8: Writing summary ===\n")
  write_summary_statistics(profile, filtered, all_edges, g, tables_dir) # Original: cfg$output_dir

  cat(sprintf("\nDone! %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("  Tables: %s\n", tables_dir))
  cat(sprintf("  Plots:  %s\n", plots_dir))
}

main()
