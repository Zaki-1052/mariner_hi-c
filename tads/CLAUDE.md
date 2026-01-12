# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**TADCompare Differential TAD Boundary Analysis Pipeline** - A bioinformatics pipeline for comparing TAD (Topologically Associating Domain) boundaries between control (wildtype) and BAP1-KO mutant mouse cerebellum Hi-C samples using the TADCompare R/Bioconductor package.

**Key Technologies:**
- R/Bioconductor: TADCompare, GenomicRanges, rtracklayer, clusterProfiler
- Juicer `straw` tool for .hic matrix extraction
- Bash/SLURM for HPC job orchestration
- bedtools for genomic distance calculations

**Experimental Context:**
- Organism: Mouse (mm10 genome)
- Conditions: Control vs BAP1-KO mutant
- Replicates: 3 biological replicates per condition
- Resolution: 25kb
- Timepoint: Late (250402)

## Running the Pipeline

### Full Pipeline (HPC with SLURM)

```bash
# From tads/ directory on Expanse
./scripts/run_pipeline.sh all

# Or individual steps:
./scripts/run_pipeline.sh 1  # Extract matrices from .hic files
./scripts/run_pipeline.sh 2  # Run TADCompare differential analysis
./scripts/run_pipeline.sh 3  # ConsensusTADs robustness check
./scripts/run_pipeline.sh 4  # Post-processing (shift distances)
```

### Direct SLURM Submission

```bash
sbatch scripts/01_extract_matrices.sb
sbatch scripts/02_run_tadcompare.sb
sbatch scripts/03_run_consensus.sb
sbatch scripts/04_postprocess.sb
sbatch scripts/05_visualizations.sb
```

### Running Visualizations Locally

```bash
cd /path/to/tads
Rscript scripts/tad_visualizations.R --timepoint late
```

### Running Individual R Scripts

```bash
# TADCompare analysis
Rscript scripts/02_run_tadcompare.R

# ConsensusTADs robustness
Rscript scripts/03_run_consensus.R

# Shift distance analysis
Rscript scripts/04_analyze_shifts.R
```

## Pipeline Architecture: 5-Step Workflow

| Step | Script | Purpose | Key Output |
|------|--------|---------|------------|
| 1 | `01_extract_matrices.sb` | Extract sparse contact matrices from .hic files | `data/extracted/merged/*.txt`, `data/extracted/replicates/*.txt` |
| 2 | `02_run_tadcompare.R` | Differential boundary detection (merged ctrl vs mut) | `results/tadcompare/tadcompare_all_boundaries.tsv` |
| 3 | `03_run_consensus.R` | Replicate robustness assessment | `results/consensus/tadcompare_with_robustness.tsv` |
| 4 | `04_analyze_shifts.R` | Calculate shift distances for shifted boundaries | `results/final/tadcompare_final_annotated.tsv` |
| 5 | `tad_visualizations.R` | Generate 40+ publication plots | `results/visualizations/late/` |

## Key Data Structures

**TADCompare TAD_Frame columns:**
- `Boundary`: Genomic coordinate (bp)
- `Gap_Score`: Differential boundary score (Z-score)
- `TAD_Score1`, `TAD_Score2`: Boundary scores in control/mutant
- `Differential`: "Differential" or "Non-Differential"
- `Enriched_In`: "Matrix 1" (control) or "Matrix 2" (mutant)
- `Type`: Split, Merge, Shifted, Strength Change, Complex, Non-Differential

**Sparse matrix format (straw output):**
```
region_i    region_j    counts
0           0           1023
0           25000       45
```

## File Paths (HPC - Expanse)

**Base directory:** `/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads`

**Input .hic files:**
```
/expanse/lustre/projects/csd940/zalibhai/stripes/StripeCaller/data/hic/250402/
├── ctrl_merged.hic, mut_merged.hic          # Merged matrices
├── ctrl_M1.hic, ctrl_M2.hic, ctrl_M3.hic    # Control replicates
└── mut_M1.hic, mut_M2.hic, mut_M3.hic       # Mutant replicates
```

**Output structure:**
```
results/
├── tadcompare/          # Step 2: TADCompare raw output
├── consensus/           # Step 3: Robustness annotations
├── final/               # Step 4: Final annotated results
└── visualizations/late/ # Step 5: All plots (9 subdirectories)
```

## Code Architecture Patterns

### Configuration at Top of Scripts

All R scripts define paths and parameters at the top:
```r
RESOLUTION <- 25000
BASE_DIR <- "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads"
CHROMOSOMES <- paste0("chr", 1:19)  # Autosomes only, excludes chrX
```

### Per-Chromosome Processing Pattern

TADCompare runs chromosome-by-chromosome, results combined:
```r
results_list <- list()
for (chrom in CHROMOSOMES) {
    result <- run_tadcompare_chr(chrom)
    if (!is.null(result)) results_list[[chrom]] <- result
}
all_boundaries <- bind_rows(results_list)
```

### Error Handling with tryCatch

All matrix operations wrapped:
```r
result <- tryCatch({
    TADCompare(cont_mat1 = ctrl_mat, cont_mat2 = mut_mat, resolution = RESOLUTION)
}, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    return(NULL)
})
```

### Progress Indicators

Scripts use emoji indicators for status:
- `✓` / `✅` - Success
- `✗` / `❌` - Failure
- Extensive `cat()` logging for long operations

## Key Dependencies

**R packages:**
- TADCompare (differential boundary detection)
- dplyr, readr, tidyverse (data wrangling)
- ggplot2, patchwork, RColorBrewer (visualization)
- GenomicRanges, rtracklayer (genomic operations)
- TxDb.Mmusculus.UCSC.mm10.knownGene, org.Mm.eg.db (annotations)
- clusterProfiler, enrichplot (GO/KEGG enrichment)

**System tools:**
- `straw` (Juicer contact matrix extraction)
- `bedtools` (genomic distance calculations)
- conda environment: `mariner_env`

**SLURM configuration:**
- Account: csd940
- Partition: shared
- Typical resources: 4 CPUs, 16-32GB RAM, 2-4 hours

## Important Constraints

1. **Exclude chrX**: All chromosome lists exclude sex chromosomes (autosomes 1-19 only)

2. **Resolution is 25kb**: Not 40kb as mentioned in context doc - implementation uses 25kb

3. **Matrix format**: TADCompare accepts sparse 3-column format (region_i, region_j, counts)

4. **No pre_tads argument**: TADCompare detects boundaries internally; existing Homer TAD files are cross-condition and don't match expected format

5. **Shift distance requires post-processing**: TADCompare classifies "Shifted" boundaries but doesn't calculate distances; use bedtools closest

## Biological Context

**TAD boundaries** are regions where chromatin domains change, often marked by CTCF binding. This pipeline identifies:
- Differential boundaries between conditions (16.7% of all boundaries)
- Boundary type changes: Shifted (3.1%), Split, Merge, Strength Change, Complex
- Enrichment direction: Control-enriched (38.4%) vs Mutant-enriched (61.6%)

**Key findings:**
- syt1/nav3 locus (chr10) shows high impact
- Long-range loops appear lost in mutant
- Compartments weakening (±0.25 range) rather than A/B switching

## Visualization Categories

`tad_visualizations.R` generates plots in 9 subdirectories:
1. **overview/** - Gap Score distributions, TAD Score scatter
2. **classification/** - Boundary type pie charts, type by chromosome
3. **shift_analysis/** - Shift distance histograms/violins
4. **robustness/** - Robustness × differential heatmaps
5. **chromosome/** - Per-chromosome differential percentages
6. **chipseq/** - H3K27ac/H3K27me3/H3K4me1 overlap analysis
7. **enrichment/** - GO BP/CC/MF and KEGG pathway analysis
8. **syt1_nav3_focus/** - Locus-specific chr10 analysis
9. **summary/** - Summary statistics document
