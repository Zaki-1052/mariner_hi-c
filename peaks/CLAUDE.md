# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Module Overview

**Loop Anchor Annotation Module** - Annotates differential chromatin loop anchors with 8 chromatin state categories based on histone modification ChIP-seq data (H3K27ac, H3K27me3, H3K4me1, H3K4me3, CTCF). Part of the larger Mariner Hi-C Differential Chromatin Loop Analysis Pipeline.

**Input:** Differential loops from parent pipeline (`{timepoint}_outputs/merged_loops/non_redundant_loops.tsv`)
**Output:** Annotated loops with chromatin state classifications, loop types, and visualizations

## Running the Scripts

```bash
cd peaks/

# Main annotation (both early & late timepoints)
Rscript annotate_loops_extended.R

# Single timepoint
Rscript annotate_loops_extended.R --timepoint early
Rscript annotate_loops_extended.R --timepoint late
Rscript annotate_loops_extended.R --timepoint late_consensus
Rscript annotate_loops_extended.R --timepoint early_p12ctrl

# All 4 timepoints
Rscript annotate_loops_extended.R --all

# Extract "Other" category up loops (early only)
Rscript extract_other_up_loops.R
```

## 8-Category Chromatin State System (Priority Order)

Classification follows strict priority - earlier categories take precedence:

| Priority | Category | Definition |
|----------|----------|------------|
| 1 | Active_Promoter | H3K4me3+ AND NOT H3K27me3 AND ≤2kb from TSS |
| 2 | Repressed_Promoter | H3K27me3+ AND NOT H3K27ac AND ≤2kb from TSS |
| 3 | Bivalent_Promoter | Overlaps pre-computed K4me3+K27me3 intersection |
| 4 | Polycomb | H3K27me3+ AND >2kb from TSS |
| 5 | Active_Enhancer | H3K27ac+ AND >2kb from TSS |
| 6 | Poised_Enhancer | H3K4me1+ AND NOT H3K27ac AND NOT H3K27me3 AND >2kb |
| 7 | CTCF_Site | CTCF+ AND not classified above |
| 8 | Other | Default (no marks) |

Loop types are 36 combinations of anchor types (e.g., `Active_Promoter-Active_Enhancer`).

## Key Data Files

**ChIP-seq Peaks (in `beds/`):**
- `H3K27ac{Tissue}{Timepoint}{Rep}.bed` - Active enhancer/promoter mark
- `H3K27me3{Tissue}{Timepoint}{Rep}.bed` - Polycomb repression mark
- `H3K4me1{Tissue}{Timepoint}{Rep}.bed` - Enhancer mark
- `H3K4me3{Tissue}{Timepoint}{Rep}.bed` - Active promoter mark
- `Bivalent_{Tissue}_{Timepoint}.bed` - Pre-computed K4me3+K27me3 overlap
- `../CTCF.bed` - Structural/insulator sites (32,487 peaks)

**Generate bivalent files:**
```bash
bedtools intersect -a H3K4me3.bed -b H3K27me3.bed -u > Bivalent.bed
```

## Core Script Functions (`annotate_loops_extended.R`)

| Function | Purpose |
|----------|---------|
| `load_chip_peaks()` | Load BED files into GRanges |
| `annotate_chip_overlaps_extended()` | Compute anchor-peak overlaps |
| `classify_anchor_type_extended()` | Apply priority-based classification |
| `classify_loop_type_extended()` | Combine anchor types into loop types |

## Output Structure

```
loop_annotation_extended/{timepoint}/
├── extended_characterized_loops.tsv   # Full annotation table
├── anchor_type_summary.tsv            # Anchor type counts
├── loop_type_summary.tsv              # Loop type counts
├── summary_statistics.txt             # Human-readable summary
└── plots/
    ├── loop_type_piechart_comparison.pdf
    ├── anchor_type_distribution.pdf
    └── loop_type_by_direction.pdf
```

## Adding New Timepoints/Peak Files

1. Add BED files to `beds/` with naming convention: `{Mark}{Tissue}{Timepoint}{Rep}.bed`
2. Generate bivalent file via bedtools intersection
3. Update `PEAK_FILES` list in `annotate_loops_extended.R`:
```r
PEAK_FILES <- list(
  new_timepoint = list(
    h3k27ac = "beds/H3K27acNew.bed",
    h3k27me3 = "beds/H3K27me3New.bed",
    h3k4me1 = "beds/H3K4me1New.bed",
    h3k4me3 = "beds/H3K4me3New.bed",
    bivalent = "beds/Bivalent_New.bed",
    ctcf = "CTCF.bed"
  )
)
```
4. Add entries to `DEFAULT_INPUT_FILES` and `DEFAULT_OUTPUT_DIRS`

## Dependencies

**R/Bioconductor:**
- GenomicRanges, InteractionSet, rtracklayer
- TxDb.Mmusculus.UCSC.mm10.knownGene, org.Mm.eg.db
- tidyverse, ggplot2, patchwork, RColorBrewer

**Runtime:** 2-5 min/timepoint | **Memory:** 4-8 GB RAM

## QC Checkpoints

- **Bivalent count** should be subset of H3K4me3 peaks
- **Active_Promoter < H3K4me3** (not all H3K4me3 near TSS)
- **>50% "Other"** may indicate missing peak files
- **Anchor1 ≈ Anchor2** distributions expected

## Known Caveats

- **Early timepoint lacks CTCF data** - "Other" category is ~36.5% (vs 9.3% in late)
- ~27% of early "Other" anchors overlap late CTCF peaks (validation)
- 82% of early up loops have "Other" anchors (structural/CTCF-mediated)
