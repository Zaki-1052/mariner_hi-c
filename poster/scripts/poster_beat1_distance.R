# poster/scripts/poster_beat1_distance.R
# Beat 1 poster figure: density (all loops) + CDF (K27me3-anchored loops) side-by-side
# Run from repo root: Rscript poster/scripts/poster_beat1_distance.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(scales)
  library(GenomicRanges)
})

BASE_DIR <- getwd()
source(file.path(BASE_DIR, "data/scripts/_shared/multi_format_output.R"))

# ---- Configuration ----------------------------------------------------------

LOOPS_TSV  <- file.path(BASE_DIR, "loops/outputs/250402-late_outputs/merged_loops/characterized_loops.tsv")
K27ME3_BED <- file.path(BASE_DIR, "peaks/beds/H3K27me3CerebellumLate1.bed")
OUT_BASE   <- file.path(BASE_DIR, "poster/figures/beat1_distance_composite")

stopifnot(file.exists(LOOPS_TSV), file.exists(K27ME3_BED))

COL_LOST   <- "#d73027"
COL_GAINED <- "#4575b4"

# ---- Load data --------------------------------------------------------------

cat("Loading loops...\n")
loops <- read_tsv(LOOPS_TSV, show_col_types = FALSE)

loops_dir <- loops %>%
  filter(direction %in% c("up_in_mutant", "down_in_mutant")) %>%
  mutate(
    direction_label = factor(
      ifelse(direction == "down_in_mutant", "Lost in BAP1-KO", "Gained in BAP1-KO"),
      levels = c("Lost in BAP1-KO", "Gained in BAP1-KO")
    ),
    loop_distance_kb = loop_distance / 1000
  )

n_lost  <- sum(loops_dir$direction == "down_in_mutant")
n_gained <- sum(loops_dir$direction == "up_in_mutant")
cat(sprintf("  %d directional loops (%d lost, %d gained)\n", nrow(loops_dir), n_lost, n_gained))

# ---- K27me3 subset ----------------------------------------------------------

cat("Computing K27me3 overlaps...\n")
k27me3_bed <- read.table(K27ME3_BED, sep = "\t", stringsAsFactors = FALSE)
k27me3_gr <- GRanges(seqnames = k27me3_bed$V1,
                     ranges = IRanges(start = k27me3_bed$V2, end = k27me3_bed$V3))

anchor1_gr <- GRanges(
  seqnames = loops_dir$anchor1_chr,
  ranges = IRanges(start = loops_dir$anchor1_start, end = loops_dir$anchor1_end)
)
anchor2_gr <- GRanges(
  seqnames = loops_dir$anchor2_chr,
  ranges = IRanges(start = loops_dir$anchor2_start, end = loops_dir$anchor2_end)
)

has_k27me3 <- (countOverlaps(anchor1_gr, k27me3_gr) > 0) |
              (countOverlaps(anchor2_gr, k27me3_gr) > 0)

loops_k27 <- loops_dir[has_k27me3, ]
n_k27_lost   <- sum(loops_k27$direction == "down_in_mutant")
n_k27_gained <- sum(loops_k27$direction == "up_in_mutant")
cat(sprintf("  K27me3-anchored: %d loops (%d lost, %d gained)\n",
            nrow(loops_k27), n_k27_lost, n_k27_gained))

# ---- Stats ------------------------------------------------------------------

med_all_lost   <- median(loops_dir$loop_distance_kb[loops_dir$direction == "down_in_mutant"])
med_all_gained <- median(loops_dir$loop_distance_kb[loops_dir$direction == "up_in_mutant"])

med_k27_lost   <- median(loops_k27$loop_distance_kb[loops_k27$direction == "down_in_mutant"])
med_k27_gained <- median(loops_k27$loop_distance_kb[loops_k27$direction == "up_in_mutant"])

cat(sprintf("  All loops median: lost = %.0f kb, gained = %.0f kb\n", med_all_lost, med_all_gained))
cat(sprintf("  K27me3 median:    lost = %.0f kb, gained = %.0f kb\n", med_k27_lost, med_k27_gained))

# ---- Poster theme -----------------------------------------------------------

theme_poster <- function() {
  theme_minimal(base_size = 20) +
    theme(
      plot.title       = element_text(face = "bold", size = 22, hjust = 0.5),
      plot.subtitle    = element_text(size = 16, hjust = 0.5, color = "gray40"),
      axis.title       = element_text(face = "bold", size = 18),
      axis.text        = element_text(size = 16),
      legend.text      = element_text(size = 16),
      legend.title     = element_text(size = 13, color = "gray40", hjust = 0.5),
      legend.title.position = "bottom",
      legend.position  = "top",
      panel.grid.minor = element_blank()
    )
}

# ---- Left panel: density (all loops) ----------------------------------------

cat("Building density panel...\n")

p_density <- ggplot(loops_dir, aes(x = loop_distance_kb, fill = direction_label)) +
  geom_density(alpha = 0.5, color = "black", linewidth = 0.4) +
  geom_vline(xintercept = med_all_lost, color = COL_LOST,
             linetype = "dashed", linewidth = 1.0) +
  geom_vline(xintercept = med_all_gained, color = COL_GAINED,
             linetype = "dashed", linewidth = 1.0) +
  annotate("text", x = med_all_gained * 0.38, y = 0.92,
           label = sprintf("%.0f kb", med_all_gained),
           color = COL_GAINED, fontface = "bold", size = 6, hjust = 1) +
  annotate("text", x = med_all_lost * 2.2, y = 0.92,
           label = sprintf("%.0f kb", med_all_lost),
           color = COL_LOST, fontface = "bold", size = 6, hjust = 0) +
  annotate("text", x = 15, y = 0.55,
           label = sprintf("n = %s", format(nrow(loops_dir), big.mark = ",")),
           size = 5, color = "gray40", hjust = 0, fontface = "italic") +
  scale_x_log10(
    labels = label_comma(),
    breaks = c(10, 100, 1000, 10000)
  ) +
  scale_fill_manual(
    values = c("Lost in BAP1-KO" = COL_LOST,
               "Gained in BAP1-KO" = COL_GAINED)
  ) +
  guides(fill = guide_legend(
    title = sprintf("n = %s lost  |  n = %s gained",
                    format(n_lost, big.mark = ","),
                    format(n_gained, big.mark = ","))
  )) +
  labs(
    title = "All differential loops",
    x = "Loop distance (kb, log scale)",
    y = "Density"
  ) +
  theme_poster()

# ---- Right panel: CDF (K27me3-anchored) ------------------------------------

cat("Building CDF panel...\n")

p_cdf <- ggplot(loops_k27, aes(x = loop_distance_kb, color = direction_label)) +
  stat_ecdf(geom = "step", linewidth = 1.8) +
  geom_vline(xintercept = med_k27_lost, color = COL_LOST,
             linetype = "dashed", linewidth = 1.0) +
  geom_vline(xintercept = med_k27_gained, color = COL_GAINED,
             linetype = "dashed", linewidth = 1.0) +
  annotate("text", x = med_k27_gained, y = 0.55,
           label = sprintf("%.0f kb", med_k27_gained),
           color = COL_GAINED, fontface = "bold", size = 6,
           hjust = 1.15, vjust = 0.5) +
  annotate("text", x = med_k27_lost, y = 0.55,
           label = sprintf("%.0f kb", med_k27_lost),
           color = COL_LOST, fontface = "bold", size = 6,
           hjust = -0.15, vjust = 0.5) +
  scale_x_log10(
    labels = label_comma(),
    breaks = c(10, 100, 1000, 10000)
  ) +
  scale_color_manual(
    values = c("Lost in BAP1-KO" = COL_LOST,
               "Gained in BAP1-KO" = COL_GAINED)
  ) +
  guides(color = guide_legend(
    title = sprintf("n = %s lost  |  n = %s gained",
                    format(n_k27_lost, big.mark = ","),
                    format(n_k27_gained, big.mark = ","))
  )) +
  labs(
    title = "Polycomb-marked loops",
    subtitle = sprintf("n = %s (K27me3-anchored)", format(nrow(loops_k27), big.mark = ",")),
    x = "Loop distance (kb, log scale)",
    y = "Cumulative proportion"
  ) +
  theme_poster()

# ---- Composite --------------------------------------------------------------

cat("Assembling composite...\n")
combined <- free(p_density) | free(p_cdf)
combined <- combined +
  plot_annotation(
    caption = "Filtering to Polycomb-marked loops amplifies the distance shift",
    theme = theme(
      plot.caption = element_text(hjust = 0.5, size = 18, face = "italic", color = "grey30")
    )
  )

save_multiformat_ggplot(combined, OUT_BASE, width = 16, height = 7, dpi = 300)

cat("\nDone.\n")
