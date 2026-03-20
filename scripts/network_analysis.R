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
})

source("scripts/utils/multi_format_output.R")

# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

BASE_DIR <- getwd()

TIMEPOINT_CONFIG <- list(
  late = list(
    gene_summary      = file.path(BASE_DIR, "abc/results/gene_level_summary.tsv"),
    loops_genes       = file.path(BASE_DIR, "abc/results/loops_with_gene_assignments.tsv"),
    boundary_genes    = file.path(BASE_DIR, "tads/results/visualizations/late/enrichment/boundary_genes.tsv"),
    delta_abc_pairs   = file.path(BASE_DIR, "abc/results/delta_abc_with_rnaseq.tsv"),
    output_dir        = file.path(BASE_DIR, "output/network_analysis/late"),
    label             = "Late Timepoint (Adult)"
  ),
  early = list(
    gene_summary      = file.path(BASE_DIR, "abc/results/gene_level_summary.tsv"),
    loops_genes       = NA,
    boundary_genes    = file.path(BASE_DIR, "tads/results/visualizations/early/enrichment/boundary_genes.tsv"),
    delta_abc_pairs   = NA,
    output_dir        = file.path(BASE_DIR, "output/network_analysis/early"),
    label             = "Early Timepoint (P13)"
  )
)

THRESHOLDS <- list(
  min_layers           = 2L,
  max_nodes            = 150L,
  abc_delta_min        = 0.05,
  baseMean_min         = 10,
  shared_enh_delta_min = 0.005,
  go_pvalue_cutoff     = 0.05,
  go_qvalue_cutoff     = 0.1,
  top_go_terms         = 8L,
  label_top_n          = 25L
)

COLORS <- list(
  expression_low  = "#4575b4",
  expression_mid  = "white",
  expression_high = "#d73027",
  edge_loop       = "#333333",
  edge_enhancer   = "#FF7F00",
  hull_alpha      = 0.12
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

  if (params$timepoint == "early") {
    stop("Early timepoint data not yet available for network analysis.")
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

load_loop_gene_pairs <- function(path) {
  stopifnot(file.exists(path))
  df <- read_tsv(path, show_col_types = FALSE)
  required <- c("loop_id", "anchor1_gene", "anchor2_gene", "logFC", "FDR", "direction")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) stop(sprintf("loops_with_gene_assignments missing columns: %s", paste(missing, collapse = ", ")))
  cat(sprintf("  Loop-gene pairs: %d rows (%d unique loops)\n", nrow(df), n_distinct(df$loop_id)))
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
# 8. GO-BASED NODE GROUPING
# ==============================================================================

assign_go_groups <- function(gene_symbols, top_n) {
  mapping <- tryCatch({
    AnnotationDbi::select(
      org.Mm.eg.db,
      keys = unique(gene_symbols),
      columns = "ENTREZID",
      keytype = "SYMBOL"
    ) %>%
      filter(!is.na(ENTREZID)) %>%
      distinct(SYMBOL, .keep_all = TRUE)
  }, error = function(e) {
    cat(sprintf("  WARNING: Symbol-to-Entrez mapping failed: %s\n", e$message))
    return(tibble(SYMBOL = character(), ENTREZID = character()))
  })

  if (nrow(mapping) == 0) {
    cat("  WARNING: No gene mappings found. Assigning all to 'Other'\n")
    return(tibble(gene = gene_symbols, go_group = "Other"))
  }

  cat(sprintf("  Mapped %d/%d gene symbols to Entrez IDs\n", nrow(mapping), length(gene_symbols)))

  ego <- tryCatch({
    enrichGO(
      gene          = mapping$ENTREZID,
      OrgDb         = org.Mm.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = THRESHOLDS$go_pvalue_cutoff,
      qvalueCutoff  = THRESHOLDS$go_qvalue_cutoff,
      readable      = TRUE
    )
  }, error = function(e) {
    cat(sprintf("  WARNING: enrichGO failed: %s\n", e$message))
    return(NULL)
  })

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    cat("  WARNING: No significant GO terms found. Assigning all to 'Other'\n")
    return(list(
      groups = tibble(gene = gene_symbols, go_group = "Other"),
      ego = NULL
    ))
  }

  ego_df <- as.data.frame(ego)
  cat(sprintf("  Found %d significant GO BP terms\n", nrow(ego_df)))

  top_terms <- ego_df %>% slice_head(n = top_n)

  gene_to_term <- top_terms %>%
    mutate(genes = str_split(geneID, "/")) %>%
    unnest(genes) %>%
    group_by(genes) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(gene = genes, go_group = Description)

  result <- tibble(gene = gene_symbols) %>%
    left_join(gene_to_term, by = "gene") %>%
    mutate(go_group = replace_na(go_group, "Other"))

  cat(sprintf("  Assigned %d genes to GO groups (%d to 'Other')\n",
              sum(result$go_group != "Other"), sum(result$go_group == "Other")))

  list(groups = result, ego = ego)
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
# 10. NETWORK VISUALIZATION (FIGURE 5C)
# ==============================================================================

plot_network <- function(g, thresholds, colors, cfg_label) {
  set.seed(42)

  # Remove isolated nodes for cleaner visualization
  g_connected <- g %>%
    filter(centrality_degree() > 0)

  n_total <- igraph::vcount(as.igraph(g))
  n_shown <- igraph::vcount(as.igraph(g_connected))
  n_edges <- igraph::ecount(as.igraph(g_connected))
  cat(sprintf("  Plotting %d connected genes (%d isolated removed)\n", n_shown, n_total - n_shown))

  # Add computed node aesthetics
  g_connected <- g_connected %>%
    mutate(
      node_size = pmax(abs(replace_na(max_delta_unnorm, 0)), 0.005),
      border_width = case_when(
        n_layers == 3 ~ 1.5,
        n_layers == 2 ~ 0.7,
        TRUE          ~ 0.3
      ),
      show_label = rank(-combined_score, ties.method = "first") <= thresholds$label_top_n,
      label_text = ifelse(show_label, name, NA_character_)
    )

  # Compute layout: Fruchterman-Reingold for better clustering of small components
  layout <- create_layout(g_connected, layout = "fr")

  # Identify GO groups with 3+ genes for hull drawing
  hull_groups <- layout %>%
    filter(go_group != "Other") %>%
    count(go_group) %>%
    filter(n >= 3) %>%
    pull(go_group)

  hull_data <- layout %>% filter(go_group %in% hull_groups)

  # Build plot using the pre-computed layout
  p <- ggraph(layout) +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    scale_y_continuous(expand = expansion(mult = 0.15))

  # Hull layer (behind everything)
  if (nrow(hull_data) > 0) {
    p <- p +
      geom_mark_hull(
        data = hull_data,
        aes(x = x, y = y, fill = go_group),
        alpha = colors$hull_alpha,
        expand = unit(4, "mm"),
        radius = unit(4, "mm"),
        show.legend = TRUE
      ) +
      scale_fill_brewer(palette = "Set2", name = "GO Group")
  }

  p <- p +
    geom_edge_link(
      aes(edge_width = weight, edge_colour = edge_type),
      alpha = 0.35,
      show.legend = TRUE
    ) +
    scale_edge_colour_manual(
      values = c("loop" = colors$edge_loop, "shared_enhancer" = colors$edge_enhancer),
      name = "Edge Type"
    ) +
    scale_edge_width(range = c(0.2, 2.0), guide = "none") +
    new_scale_fill() +
    geom_node_point(
      aes(size = node_size, fill = log2FC, stroke = border_width),
      shape = 21,
      colour = "black"
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
    scale_size_continuous(range = c(2, 12), name = "|max delta AxC|") +
    geom_node_text(
      aes(label = label_text),
      repel = TRUE,
      size = 2.5,
      max.overlaps = 25,
      segment.size = 0.2,
      segment.alpha = 0.4,
      na.rm = TRUE
    ) +
    theme_graph(base_family = "") +
    labs(
      title = sprintf("Multi-Layer Structural Disruption Network (BAP1-KO, %s)", cfg_label),
      subtitle = sprintf("%d connected genes (%d total in >=%d layers) | %d edges",
                         n_shown, n_total, thresholds$min_layers, n_edges)
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

plot_go_enrichment <- function(ego) {
  if (is.null(ego)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No significant GO terms") + theme_void())
  }

  dotplot(ego, showCategory = 20) +
    labs(title = "GO Biological Process Enrichment (Network Genes)") +
    THEME_PUB +
    theme(axis.text.y = element_text(size = 9))
}

# ==============================================================================
# 12. SUMMARY STATISTICS
# ==============================================================================

write_summary_statistics <- function(profile, filtered, edges, g, ego, output_dir) {
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

  if (!is.null(ego)) {
    ego_df <- as.data.frame(ego)
    summary_text <- paste0(summary_text, sprintf(
      "\n--- GO Enrichment ---\nTotal significant terms: %d\nTop 5 terms:\n", nrow(ego_df)
    ))
    for (i in seq_len(min(5, nrow(ego_df)))) {
      summary_text <- paste0(summary_text, sprintf(
        "  %d. %s (p.adj=%.2e, %d genes)\n",
        i, ego_df$Description[i], ego_df$p.adjust[i], ego_df$Count[i]
      ))
    }
  }

  summary_text <- paste0(summary_text, sprintf(
    "\n--- Parameters ---\nmin_layers: %d\nmax_nodes: %d\nabc_delta_min: %.3f\nbaseMean_min: %d\nshared_enh_delta_min: %.3f\n",
    THRESHOLDS$min_layers, THRESHOLDS$max_nodes, THRESHOLDS$abc_delta_min,
    THRESHOLDS$baseMean_min, THRESHOLDS$shared_enh_delta_min
  ))

  write(summary_text, file = file.path(output_dir, "tables", "summary_report.txt"))
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
  tables_dir <- file.path(cfg$output_dir, "tables")
  plots_dir  <- file.path(cfg$output_dir, "plots")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_dir,  recursive = TRUE, showWarnings = FALSE)

  # --- Validate inputs ---
  cat("Validating input files...\n")
  input_files <- list(
    gene_summary    = cfg$gene_summary,
    loops_genes     = cfg$loops_genes,
    boundary_genes  = cfg$boundary_genes,
    delta_abc_pairs = cfg$delta_abc_pairs
  )
  for (nm in names(input_files)) {
    stopifnot(file.exists(input_files[[nm]]))
    cat(sprintf("  [OK] %s\n", nm))
  }

  # ==========================================================================
  # SECTION 1: Load Data
  # ==========================================================================
  cat("\n=== Section 1: Loading input data ===\n")
  gene_summary   <- load_gene_summary(cfg$gene_summary)
  loops          <- load_loop_gene_pairs(cfg$loops_genes)
  boundary_genes <- load_boundary_genes(cfg$boundary_genes)
  delta_abc      <- load_delta_abc_pairs(cfg$delta_abc_pairs)

  # ==========================================================================
  # SECTION 2: Build Gene Profile
  # ==========================================================================
  cat("\n=== Section 2: Building per-gene structural profile ===\n")
  profile <- build_gene_structural_profile(gene_summary, loops, boundary_genes, THRESHOLDS$abc_delta_min)
  write_tsv(profile, file.path(tables_dir, "gene_structural_profile_all.tsv"))
  cat(sprintf("  Saved: gene_structural_profile_all.tsv (%d genes)\n", nrow(profile)))

  # ==========================================================================
  # SECTION 3: Filter Network Genes
  # ==========================================================================
  cat("\n=== Section 3: Filtering genes for network ===\n")
  filtered <- filter_network_genes(profile, THRESHOLDS)
  write_tsv(filtered, file.path(tables_dir, "gene_structural_profile_filtered.tsv"))

  # ==========================================================================
  # SECTION 4: Build Edges
  # ==========================================================================
  cat("\n=== Section 4: Building edge list ===\n")
  loop_edges <- build_loop_edges(filtered$gene, loops)
  enh_edges  <- build_shared_enhancer_edges(filtered$gene, delta_abc, THRESHOLDS$shared_enh_delta_min)
  all_edges  <- bind_rows(
    loop_edges %>% dplyr::select(from, to, weight, edge_type),
    enh_edges  %>% dplyr::select(from, to, weight, edge_type)
  )
  write_tsv(all_edges, file.path(tables_dir, "edge_list.tsv"))
  cat(sprintf("  Total edges: %d\n", nrow(all_edges)))

  # ==========================================================================
  # SECTION 5: GO Grouping
  # ==========================================================================
  cat("\n=== Section 5: GO enrichment & gene grouping ===\n")
  go_result <- assign_go_groups(filtered$gene, THRESHOLDS$top_go_terms)
  filtered  <- filtered %>% left_join(go_result$groups, by = "gene")

  if (!is.null(go_result$ego)) {
    ego_df <- as.data.frame(go_result$ego)
    write_tsv(ego_df, file.path(tables_dir, "go_enrichment_results.tsv"))
  }

  # ==========================================================================
  # SECTION 6: Build Network
  # ==========================================================================
  cat("\n=== Section 6: Constructing network ===\n")
  g <- build_network(filtered, all_edges)

  centrality <- as_tibble(g) %>% dplyr::select(name, degree, betweenness, closeness)
  write_tsv(centrality, file.path(tables_dir, "network_centrality_metrics.tsv"))

  # ==========================================================================
  # SECTION 7: Visualize
  # ==========================================================================
  cat("\n=== Section 7: Generating plots ===\n")

  cat("  Rendering network figure...\n")
  p_network <- plot_network(g, THRESHOLDS, COLORS, cfg$label)
  save_multiformat_ggplot(p_network, file.path(plots_dir, "network_figure5c"),
                          width = 14, height = 12)

  cat("  Rendering layer distribution...\n")
  p_layers <- plot_layer_distribution(profile)
  save_multiformat_ggplot(p_layers, file.path(plots_dir, "layer_distribution"),
                          width = 8, height = 6)

  cat("  Rendering combined score distribution...\n")
  p_scores <- plot_combined_score_distribution(filtered)
  save_multiformat_ggplot(p_scores, file.path(plots_dir, "combined_score_distribution"),
                          width = 8, height = 6)

  cat("  Rendering edge type breakdown...\n")
  p_edges <- plot_edge_type_breakdown(all_edges)
  save_multiformat_ggplot(p_edges, file.path(plots_dir, "edge_type_breakdown"),
                          width = 8, height = 6)

  cat("  Rendering degree distribution...\n")
  p_degree <- plot_degree_distribution(g)
  save_multiformat_ggplot(p_degree, file.path(plots_dir, "degree_distribution"),
                          width = 8, height = 6)

  cat("  Rendering GO enrichment dotplot...\n")
  p_go <- plot_go_enrichment(go_result$ego)
  save_multiformat_ggplot(p_go, file.path(plots_dir, "go_enrichment_dotplot"),
                          width = 10, height = 8)

  # ==========================================================================
  # SECTION 8: Summary
  # ==========================================================================
  cat("\n=== Section 8: Writing summary ===\n")
  write_summary_statistics(profile, filtered, all_edges, g, go_result$ego, cfg$output_dir)

  cat(sprintf("\nDone! %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("  Tables: %s\n", tables_dir))
  cat(sprintf("  Plots:  %s\n", plots_dir))
}

main()
