# Multi-Format Figure Output Status

**Goal:** Output all figures as PDF + SVG (Illustrator) + JPEG (Google Slides)

**Status:** ALL SCRIPTS CONVERTED - 15 scripts, ~90+ plots now output multi-format

---

## Completed Scripts (Native Multi-Format Output)

All visualization scripts now automatically output `.pdf`, `.svg`, and `.jpg` for every figure:

### Main Pipeline Scripts

| Script | Plots | Notes |
|--------|-------|-------|
| `scripts/visualizations.R` | 10 | Volcano, enrichment, loop classification |
| `scripts/tad_volcano_plot.R` | 2 | Relaxed + standard thresholds |
| `scripts/compartment_volcano_plot.R` | 2 | Relaxed + standard thresholds |
| `scripts/apa_analysis.R` | 5 | APA heatmaps, enrichment |
| `scripts/downstream_analysis.R` | 4 | Distance, chromosome, gene proximity |
| `scripts/loop_distance_analysis.R` | 9 | Loop rewriting visualizations |
| `scripts/annotate_loops_extended.R` | 3 | Anchor type distributions |
| `scripts/qc-val.R` | 9 | QC correlation, heatmaps, diagnostics |
| `scripts/edgeR.R` | 7 | MDS, BCV, volcano, MA plots |
| `scripts/compare_resolutions.R` | 5 | Resolution comparison, Venn diagrams |

### Secondary Analysis Scripts

| Script | Plots | Notes |
|--------|-------|-------|
| `stripes/scripts/stripe_visualizations.R` | 9 | Stripe analysis visualizations |
| `stripes/scripts/phase3_edgeR.R` | 5 | Stripe differential analysis |
| `tads/scripts/tad_visualizations.R` | 18+ | TAD boundary visualizations |
| `peaks/scripts/annotate_loops_extended.R` | 3 | Peak annotation plots |

### Utility Script

| Script | Notes |
|--------|-------|
| `scripts/utils/multi_format_output.R` | Shared helper functions |

---

## Total: ~90+ plots with native multi-format output

---

## How It Works

### For `ggsave()` scripts

Uses the shared utility:
```r
source("scripts/utils/multi_format_output.R")
save_multiformat_ggplot(p, file.path(output_dir, "plot"), width = 10, height = 8)
```

### For `pdf()` + `dev.off()` scripts

Each script includes an inline helper:
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

---

## Required R Package

```r
install.packages("svglite")
```

---

## Running Regeneration

The script automatically creates symlinks to the correct output directories for each timepoint.

```bash
# Run for late timepoint only (saves to 25042-late_outputs/)
bash scripts/regenerate_all_figures.sb late

# Run for early timepoint only (saves to 250831-early_outputs/)
bash scripts/regenerate_all_figures.sb early

# Run for both timepoints sequentially (default)
bash scripts/regenerate_all_figures.sb
bash scripts/regenerate_all_figures.sb both

# On HPC with SLURM
sbatch scripts/regenerate_all_figures.sb
```

**How it works:** The script creates a symlink (`outputs -> {timepoint}_outputs/`) before running each batch of visualization scripts, ensuring SVG/JPEG files are saved alongside existing PDFs in the correct directories.

---

## Output Quality Notes

- **SVG:** Native svglite output provides editable text and clean structure for Adobe Illustrator
- **JPEG:** 300 DPI at 95% quality for Google Slides presentations
- **PDF:** Vector format for publications

---

## Conversion Completed

- Jan 14, 2026: All 15 visualization scripts converted to multi-format output
