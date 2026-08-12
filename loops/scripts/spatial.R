# spatial_analysis.R
library(mariner)
# Load the original merged loops for coordinates
merged <- readRDS("outputs/full/02_merged.rds")

# Extract coordinates
loop_data <- data.frame(
  chr = as.character(seqnames1(anchors(merged, type="first"))),
  start1 = start(anchors(merged, type="first")),
  end1 = end(anchors(merged, type="first")),
  start2 = start(anchors(merged, type="second")),
  end2 = end(anchors(merged, type="second")),
  ctrl_count = counts[,1],
  mut_count = counts[,2],
  log2FC = M
)

# Add distance for intrachromosomal loops
loop_data$distance <- ifelse(
  loop_data$start2 > loop_data$start1,
  loop_data$start2 - loop_data$start1,
  NA
)

# Plot distance vs signal strength
ggplot(loop_data[!is.na(loop_data$distance),], 
       aes(x = distance/1e6, y = ctrl_count + mut_count)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Loop Signal vs Genomic Distance",
       x = "Distance (Mb, log10)",
       y = "Total Count (ctrl + mut, log10)")

# Chromosome distribution
chr_summary <- loop_data %>%
  group_by(chr) %>%
  summarise(
    n_loops = n(),
    mean_ctrl = mean(ctrl_count),
    mean_mut = mean(mut_count),
    mean_fc = mean(log2FC, na.rm = TRUE)
  )

ggplot(chr_summary, aes(x = chr, y = n_loops)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Loop Distribution by Chromosome",
       x = "Chromosome", y = "Number of Loops")

