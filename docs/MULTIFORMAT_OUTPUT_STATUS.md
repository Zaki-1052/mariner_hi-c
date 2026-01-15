# Multi-Format Figure Output Status

**Goal:** Output all figures as PDF + SVG (Illustrator) + JPEG (Google Slides)

## Completed Scripts (Native Multi-Format Output)

These scripts now automatically output `.pdf`, `.svg`, and `.jpg` for every figure:

| Script | Plots | Notes |
|--------|-------|-------|
| `scripts/visualizations.R` | 10 | Volcano, enrichment, loop classification |
| `scripts/tad_volcano_plot.R` | 2 | Relaxed + standard thresholds |
| `scripts/compartment_volcano_plot.R` | 2 | Relaxed + standard thresholds |
| `scripts/apa_analysis.R` | 5 | APA heatmaps, enrichment |
| `scripts/downstream_analysis.R` | 4 | Distance, chromosome, gene proximity |
| `scripts/loop_distance_analysis.R` | 9 | Loop rewriting visualizations |
| `scripts/annotate_loops_extended.R` | 3 | Anchor type distributions |

**Total: ~35 plots with native multi-format output**

## Scripts Still Needing Modification (PDF Only)

### High Priority (Main Pipeline QC)

| Script | Plots | Pattern | Difficulty |
|--------|-------|---------|------------|
| `scripts/qc-val.R` | 9 | `pdf()` + `dev.off()` | Medium (pheatmap + base R) |
| `scripts/edgeR.R` | 7 | `pdf()` + `dev.off()` | Medium (edgeR plots) |
| `scripts/compare_resolutions.R` | 5 | `pdf()` + `dev.off()` | Medium (base R + Venn) |

### Medium Priority (Secondary Analyses)

| Script | Plots | Pattern | Difficulty |
|--------|-------|---------|------------|
| `stripes/scripts/stripe_visualizations.R` | 9 | `ggsave()` | Easy |
| `stripes/scripts/phase3_edgeR.R` | 5 | `pdf()` + `dev.off()` | Medium |
| `tads/scripts/tad_visualizations.R` | 21 | `ggsave()` | Easy |
| `peaks/scripts/annotate_loops_extended.R` | 3 | `ggsave()` | Easy |

**Total remaining: ~59 plots**

## How to Modify Remaining Scripts

### For `ggsave()` scripts (Easy)

Replace:
```r
ggsave(file.path(output_dir, "plot.pdf"), p, width = 10, height = 8)
```

With:
```r
source("scripts/utils/multi_format_output.R")
save_multiformat_ggplot(p, file.path(output_dir, "plot"), width = 10, height = 8)
```

### For `pdf()` + `dev.off()` scripts (Medium)

Add this helper at the top of the script:
```r
library(svglite)

save_multiformat <- function(plot_code, base_path, width, height, dpi = 300) {
  pdf(paste0(base_path, ".pdf"), width = width, height = height)
  tryCatch(eval(plot_code), finally = dev.off())
  svglite(paste0(base_path, ".svg"), width = width, height = height)
  tryCatch(eval(plot_code), finally = dev.off())
  jpeg(paste0(base_path, ".jpg"), width = width * dpi, height = height * dpi, res = dpi, quality = 95)
  tryCatch(eval(plot_code), finally = dev.off())
}
```

Then replace:
```r
pdf(file.path(dir, "plot.pdf"), width = 8, height = 6)
plot(x, y)
abline(h = 0)
dev.off()
```

With:
```r
save_multiformat(
  quote({
    plot(x, y)
    abline(h = 0)
  }),
  file.path(dir, "plot"), width = 8, height = 6)
```

## Workaround: ImageMagick Conversion

For PDF-only scripts, you can convert after generation:

```bash
# Convert all PDFs to JPEG (300 DPI)
find . -name "*.pdf" -exec sh -c 'convert -density 300 "$1" -quality 95 "${1%.pdf}.jpg"' _ {} \;
```

**Note:** PDF-to-SVG conversion produces poor quality. Native R SVG output is strongly recommended for Illustrator work.

## Required R Package

```r
install.packages("svglite")
```

## Running Regeneration

```bash
# On HPC with SLURM
sbatch scripts/regenerate_all_figures.sb

# Locally
bash scripts/regenerate_all_figures.sb
```
