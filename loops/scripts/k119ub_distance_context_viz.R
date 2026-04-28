# scripts/k119ub_distance_context_viz.R
# Quick exploratory viz: K119ub effect on loop strength by distance and chromatin context

library(tidyverse)

# --- Load and merge data ---
chars <- read_tsv("outputs/250402-late_outputs/merged_loops/characterized_loops.tsv",
                  show_col_types = FALSE)
k119  <- read_tsv("output/h2ak119ub_loop_integration/late/tables/loops_with_k119ub_annotation.tsv",
                  show_col_types = FALSE)

# Join on genomic coordinates (both files share these)
shared_cols <- intersect(names(chars), names(k119))
cat("Shared columns:", paste(shared_cols, collapse = ", "), "\n")

# Use the K119ub file as base (has all 39K loops) and pull anchor types from chars
df <- k119 %>%
  left_join(
    chars %>% select(loop_id, anchor1_type, anchor2_type, loop_type),
    by = "loop_id"
  )

cat("Rows:", nrow(df), "| With anchor types:", sum(!is.na(df$anchor1_type)), "\n")

# --- Classify chromatin context of each loop ---
polycomb_types <- c("Polycomb", "Repressed_Promoter", "Bivalent_Promoter")
active_types   <- c("Active_Enhancer", "Active_Promoter")

df <- df %>%
  mutate(
    any_polycomb = anchor1_type %in% polycomb_types | anchor2_type %in% polycomb_types,
    any_active   = anchor1_type %in% active_types   | anchor2_type %in% active_types,
    chromatin_context = case_when(
      any_polycomb & !any_active ~ "Polycomb",
      any_active & !any_polycomb ~ "Active",
      any_polycomb & any_active  ~ "Mixed",
      TRUE                       ~ "Other/Structural"
    ),
    # Distance bins
    dist_kb = loop_distance / 1000,
    dist_bin = cut(dist_kb,
                   breaks = c(0, 100, 500, 1000, Inf),
                   labels = c("<100kb", "100-500kb", "500kb-1Mb", ">1Mb"),
                   right = FALSE),
    # K119ub status
    k119ub_status = case_when(
      K119ub_up_either   ~ "K119ub gained",
      K119ub_down_either ~ "K119ub lost",
      TRUE               ~ "No K119ub change"
    )
  )

cat("\nChromatin context breakdown:\n")
print(table(df$chromatin_context, useNA = "always"))
cat("\nK119ub status:\n")
print(table(df$k119ub_status))

# --- PLOT 1: logFC vs distance, smoothed by chromatin context ---
p1 <- df %>%
  filter(chromatin_context %in% c("Polycomb", "Active")) %>%
  ggplot(aes(x = log10(dist_kb), y = logFC, color = chromatin_context)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_point(alpha = 0.03, size = 0.5) +
  geom_smooth(method = "loess", span = 0.4, linewidth = 1.5, se = TRUE) +
  scale_color_manual(values = c("Active" = "#E64B35", "Polycomb" = "#4DBBD5")) +
  scale_x_continuous(
    breaks = log10(c(50, 100, 500, 1000, 5000)),
    labels = c("50kb", "100kb", "500kb", "1Mb", "5Mb")
  ) +
  labs(
    title = "Loop strength change by distance and chromatin context",
    subtitle = "LOESS smooth | All loops (differential + unchanged)",
    x = "Loop distance (log scale)", y = "logFC (BAP1-KO vs WT)",
    color = "Anchor context"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

# --- PLOT 2: Same but faceted by K119ub status ---
p2 <- df %>%
  filter(chromatin_context %in% c("Polycomb", "Active"),
         k119ub_status != "K119ub lost") %>%
  ggplot(aes(x = log10(dist_kb), y = logFC, color = chromatin_context)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_point(alpha = 0.03, size = 0.5) +
  geom_smooth(method = "loess", span = 0.5, linewidth = 1.5, se = TRUE) +
  scale_color_manual(values = c("Active" = "#E64B35", "Polycomb" = "#4DBBD5")) +
  scale_x_continuous(
    breaks = log10(c(50, 100, 500, 1000, 5000)),
    labels = c("50kb", "100kb", "500kb", "1Mb", "5Mb")
  ) +
  facet_wrap(~k119ub_status) +
  labs(
    title = "K119ub gained vs unchanged: effect on loop strength by context",
    subtitle = "Does K119ub accumulation have opposite distance effects?",
    x = "Loop distance (log scale)", y = "logFC (BAP1-KO vs WT)",
    color = "Anchor context"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

# --- PLOT 3: Boxplot grid -- distance bin x chromatin context x K119ub ---
p3 <- df %>%
  filter(chromatin_context %in% c("Polycomb", "Active"),
         k119ub_status %in% c("K119ub gained", "No K119ub change"),
         !is.na(dist_bin)) %>%
  ggplot(aes(x = dist_bin, y = logFC, fill = k119ub_status)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_boxplot(outlier.alpha = 0.1, outlier.size = 0.3, width = 0.7) +
  facet_wrap(~chromatin_context) +
  scale_fill_manual(values = c("K119ub gained" = "#7E6148", "No K119ub change" = "#B09C85")) +
  labs(
    title = "The dual effect: K119ub impact depends on context AND distance",
    subtitle = "Polycomb: short-range UP, long-range DOWN | Active: DOWN everywhere",
    x = "Loop distance", y = "logFC (BAP1-KO vs WT)",
    fill = "K119ub status"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 30, hjust = 1))

# --- PLOT 4: Median logFC heatmap (the clearest summary) ---
summary_df <- df %>%
  filter(chromatin_context %in% c("Polycomb", "Active"),
         k119ub_status %in% c("K119ub gained", "No K119ub change"),
         !is.na(dist_bin)) %>%
  group_by(chromatin_context, dist_bin, k119ub_status) %>%
  summarise(
    median_logFC = median(logFC, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

cat("\n=== THE KEY TABLE ===\n")
summary_df %>%
  arrange(chromatin_context, k119ub_status, dist_bin) %>%
  print(n = 30)

p4 <- summary_df %>%
  ggplot(aes(x = dist_bin, y = k119ub_status, fill = median_logFC)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.3f\n(n=%d)", median_logFC, n)),
            size = 3.2, color = "black") +
  facet_wrap(~chromatin_context) +
  scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35",
                       midpoint = 0, limits = c(-0.15, 0.15)) +
  labs(
    title = "Median logFC: K119ub effect is context + distance dependent",
    subtitle = "Blue = loop weakening | Red = loop strengthening",
    x = "Loop distance", y = "", fill = "Median\nlogFC"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# --- Save ---
outdir <- "output/k119ub_context_distance_viz"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(outdir, "01_loess_by_context.pdf"), p1, width = 8, height = 6)
ggsave(file.path(outdir, "02_loess_by_k119ub_facet.pdf"), p2, width = 12, height = 6)
ggsave(file.path(outdir, "03_boxplot_grid.pdf"), p3, width = 10, height = 6)
ggsave(file.path(outdir, "04_median_heatmap.pdf"), p4, width = 10, height = 5)

cat("\nSaved 4 plots to:", outdir, "\n")
cat("Done!\n")
