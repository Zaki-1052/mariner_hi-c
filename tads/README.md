# Differential TAD Boundary Analysis Pipeline

## Overview

This pipeline performs **differential TAD (Topologically Associating Domain) boundary analysis** from Hi-C data, comparing **BAP1-KO mutant** and **wildtype control** conditions in mouse cerebellum (mm10 genome) with **biological replication** (n=3 per condition) using the TADCompare R/Bioconductor package.

### Key Features

- **TADCompare-based detection**: Uses TADCompare (Cresswell & Dozmorov, Bioinformatics 2020) for differential boundary classification
- **Replicate robustness assessment**: ConsensusTADs validates boundaries across biological replicates
- **Multi-timepoint analysis**: Processes early (P12) and late developmental stages
- **Boundary type classification**: Shifted, Split, Merge, Strength Change, Complex, Non-Differential
- **Shift distance quantification**: Calculates exact distances for shifted boundaries using bedtools
- **ChIP-seq integration**: Anchor classification using H3K27ac, H3K27me3, H3K4me1 peaks
- **GO/KEGG enrichment**: Functional analysis of genes near differential boundaries
- **Publication-ready outputs**: 40+ visualization plots across 9 categories

### Biological Context

**TAD boundaries** are genomic regions where chromatin domains change, often marked by CTCF binding and insulator elements. They represent:

- **Structural organization**: ~900kb median TAD size in mammals, with boundaries isolating regulatory domains
- **Regulatory insulation**: Prevent enhancer-promoter cross-talk between adjacent TADs
- **CTCF-cohesin mediated**: Most boundaries enriched for CTCF and cohesin complex

**BAP1** is a Polycomb regulator (H2AK119ub1 deubiquitinase). Loss of BAP1 may affect:
- Polycomb-repressed domain boundaries
- Chromatin compaction at TAD boundaries
- Insulation strength at developmental gene loci

---

## Pipeline Architecture

### Workflow Overview

```
Input: .hic contact matrices (merged + individual replicates)
  |
[1] 01_extract_matrices.sb  -> Extract sparse contact matrices via straw
  |
[2] 02_run_tadcompare.R     -> Differential boundary detection (merged ctrl vs mut)
  |
[3] 03_run_consensus.R      -> Replicate robustness assessment via ConsensusTADs
  |
[4] 04_analyze_shifts.R     -> Calculate shift distances using bedtools closest
  |
[5] tad_visualizations.R    -> Generate 40+ publication-quality plots
  |
Output: Differential boundaries with classifications, annotations, and visualizations
```

### Sample Structure

```
Control (wildtype):        Mutant (BAP1-KO):
├── ctrl_M1               ├── mut_M1
├── ctrl_M2               ├── mut_M2
├── ctrl_M3               └── mut_M3
└── ctrl_merged           └── mut_merged

Timepoints:
├── early (250831) - P12 developmental stage
└── late (250402)  - Later developmental stage
```

---

## Quick Start

### Single-Command Execution (Recommended)

```bash
# On SDSC Expanse HPC
cd /expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads
./scripts/run_pipeline.sh all
```

**Expected runtime:** ~2-4 hours (both timepoints, all phases + visualization)

**Output:** Differential TAD boundaries for early and late timepoints with all annotations and visualizations.

### Step-by-Step Execution

For manual control or debugging:

```bash
cd /expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads

# Process each timepoint
for TIMEPOINT in late early; do
  echo "Processing ${TIMEPOINT}..."

  # Step 1: Extract matrices from .hic files
  sbatch scripts/01_extract_matrices.sb  # Adjust for timepoint

  # Step 2: TADCompare differential analysis
  Rscript scripts/02_run_tadcompare.R ${TIMEPOINT}

  # Step 3: ConsensusTADs robustness check
  Rscript scripts/03_run_consensus.R ${TIMEPOINT}

  # Step 4: Post-processing (shift distances via bedtools)
  sbatch scripts/04_postprocess.sb  # Runs bedtools closest + R script

  # Step 5: Visualizations
  Rscript scripts/tad_visualizations.R --timepoint ${TIMEPOINT}
done
```

### Direct SLURM Submission

```bash
sbatch scripts/01_extract_matrices.sb
sbatch scripts/02_run_tadcompare.sb
sbatch scripts/03_run_consensus.sb
sbatch scripts/04_postprocess.sb
sbatch scripts/05_visualizations.sb
```

---

## Detailed Pipeline Steps

### Step 1: Matrix Extraction (`01_extract_matrices.sb`)

**Purpose:** Extract sparse contact matrices from .hic files using Juicer's `straw` tool.

**Input:**
- `.hic` files: Merged (ctrl_merged.hic, mut_merged.hic) + individual replicates (ctrl_M1.hic, etc.)
- Location: `/expanse/lustre/projects/csd940/zalibhai/stripes/StripeCaller/data/hic/{timepoint}/`

**Process:**
1. Extract observed contacts at 25kb resolution
2. Use NONE normalization (TADCompare handles internally)
3. Process all 19 autosomes (chr1-chr19, excludes chrX)
4. Generate sparse 3-column format (region_i, region_j, counts)

**Output:** `data/{timepoint}/extracted/`
- `merged/ctrl_merged.chr{N}.25kb.txt` - Merged control matrices
- `merged/mut_merged.chr{N}.25kb.txt` - Merged mutant matrices
- `replicates/{sample}.chr{N}.25kb.txt` - Individual replicate matrices

**Key Parameters:**
- Resolution: 25kb
- Normalization: NONE (raw observed)
- Chromosomes: chr1-chr19 (autosomes only)

---

### Step 2: TADCompare Differential Analysis (`02_run_tadcompare.R`)

**Purpose:** Identify and classify differential TAD boundaries between control and mutant.

**Input:**
- Sparse contact matrices from Step 1 (merged ctrl vs mut)

**Process:**
1. Load sparse matrices per chromosome
2. TADCompare internally:
   - Detects all TAD boundaries using TopDom-like scoring
   - Calculates Gap_Score (differential boundary score)
   - Classifies boundary types based on score patterns
3. Combine results across all chromosomes
4. Generate summary statistics

**Output:** `results/{timepoint}/tadcompare/`
- `tadcompare_all_boundaries.tsv` - All boundaries with classifications
- `tadcompare_differential_only.tsv` - Differential boundaries only
- `tadcompare_summary.txt` - Summary statistics

**TADCompare Classifications:**

| Type | Description |
|------|-------------|
| Non-Differential | Boundary unchanged between conditions |
| Strength Change | Boundary strength differs (present in both) |
| Shifted | Boundary position moved to adjacent region |
| Split | One boundary became two |
| Merge | Two boundaries became one |
| Complex | Multiple changes, difficult to classify |

**TADCompare Output Columns:**
- `Boundary`: Genomic coordinate (bp)
- `Gap_Score`: Differential boundary score (Z-score)
- `TAD_Score1`: Boundary score in control (Matrix 1)
- `TAD_Score2`: Boundary score in mutant (Matrix 2)
- `Differential`: "Differential" or "Non-Differential"
- `Enriched_In`: "Matrix 1" (control) or "Matrix 2" (mutant)
- `Type`: Classification (Shifted, Split, Merge, etc.)

---

### Step 3: ConsensusTADs Robustness (`03_run_consensus.R`)

**Purpose:** Assess boundary consistency across biological replicates.

**Input:**
- TADCompare results from Step 2
- Individual replicate matrices from Step 1

**Process:**
1. Run ConsensusTADs on 3 control replicates per chromosome
2. Run ConsensusTADs on 3 mutant replicates per chromosome
3. Match consensus boundaries to TADCompare results
4. Annotate with robustness information

**Output:** `results/{timepoint}/consensus/`
- `tadcompare_with_robustness.tsv` - All boundaries with robustness annotations
- `consensus_control.tsv` - Control consensus boundaries
- `consensus_mutant.tsv` - Mutant consensus boundaries
- `high_confidence_differential.tsv` - Robust differential boundaries
- `consensus_summary.txt` - Summary statistics

**Robustness Annotations:**
- `Robustness`: "Both", "Control_only", "Mutant_only", "Neither"
- `n_ctrl_reps`: Number of control replicates with boundary (0-3)
- `n_mut_reps`: Number of mutant replicates with boundary (0-3)

---

### Step 4: Post-Processing (`04_postprocess.sb` + `04_analyze_shifts.R`)

**Purpose:** Calculate exact shift distances for shifted boundaries.

**Input:**
- Robustness-annotated boundaries from Step 3

**Process:**
1. Extract boundaries by condition (control-enriched, mutant-enriched)
2. Use `bedtools closest -d` to find nearest boundary in other condition
3. Calculate shift distances (in bp and kb)
4. Merge distances back to main results
5. Generate final annotated output

**Output:** `results/{timepoint}/final/`
- `tadcompare_final_annotated.tsv` - Complete final results with shift distances
- `analysis_summary.txt` - Summary statistics
- `tmp/` - Intermediate bedtools files

**Shift Distance Interpretation:**
- Distance in bp/kb from control boundary to nearest mutant boundary (or vice versa)
- Relevant for "Shifted" type boundaries
- Median shift distances typically 25-50kb (1-2 resolution bins)

---

### Step 5: Visualizations (`tad_visualizations.R`)

**Purpose:** Generate publication-quality figures, ChIP-seq annotations, and functional enrichment.

**Input:**
- Final annotated boundaries from Step 4
- ChIP-seq peak files (H3K27ac, H3K27me3, H3K4me1)
- Gene annotations (TxDb.Mmusculus.UCSC.mm10.knownGene)

**Analyses:**
1. **Overview plots**: Gap Score distributions, TAD Score scatter
2. **Classification**: Boundary type pie charts, type by chromosome
3. **Shift analysis**: Distance histograms, violin plots
4. **Robustness**: Robustness x differential heatmaps
5. **Chromosome analysis**: Per-chromosome differential percentages
6. **ChIP-seq integration**: Anchor classification by histone marks
7. **Enrichment**: GO BP/CC/MF and KEGG pathway analysis
8. **syt1/nav3 focus**: Locus-specific chr10 analysis
9. **Summary**: Statistics document

**Output:** `results/visualizations/{timepoint}/`

| Subdirectory | Contents |
|--------------|----------|
| `overview/` | gap_score_distribution.pdf, tad_score_scatter.pdf, differential_landscape.pdf |
| `classification/` | boundary_type_pie.pdf, boundary_type_by_chromosome.pdf, enrichment_direction_by_type.pdf |
| `shift_analysis/` | shift_distance_histogram.pdf, shift_distance_violin.pdf, shift_vs_gap_score.pdf |
| `robustness/` | robustness_differential_heatmap.pdf, robustness_by_type.pdf |
| `chromosome/` | per_chromosome_differential.pdf |
| `chipseq/` | anchor_classification.pdf, chipseq_overlap_heatmap.pdf, chipseq_by_enrichment_direction.pdf |
| `enrichment/` | go_bp_dotplot.pdf, go_cc_dotplot.pdf, go_mf_dotplot.pdf, kegg_dotplot.pdf, boundary_genes.tsv |
| `syt1_nav3_focus/` | syt1_nav3_regional_overview.pdf, syt1_nav3_statistics.tsv |
| `summary/` | visualization_summary.txt |

---

## Output Directory Structure

```
tads/
├── data/
│   ├── early/
│   │   └── extracted/
│   │       ├── merged/          # Merged matrices per chromosome
│   │       └── replicates/      # Individual replicate matrices
│   └── late/
│       └── extracted/
│           ├── merged/
│           └── replicates/
│
├── results/
│   ├── early/
│   │   ├── tadcompare/          # Step 2: TADCompare raw output
│   │   │   ├── tadcompare_all_boundaries.tsv
│   │   │   ├── tadcompare_differential_only.tsv
│   │   │   └── tadcompare_summary.txt
│   │   ├── consensus/           # Step 3: Robustness annotations
│   │   │   ├── tadcompare_with_robustness.tsv
│   │   │   ├── consensus_control.tsv
│   │   │   ├── consensus_mutant.tsv
│   │   │   ├── high_confidence_differential.tsv
│   │   │   └── consensus_summary.txt
│   │   └── final/               # Step 4: Final annotated results
│   │       ├── tadcompare_final_annotated.tsv
│   │       └── analysis_summary.txt
│   │
│   ├── late/                    # Same structure as early
│   │
│   └── visualizations/          # Step 5: All plots
│       ├── early/               # 9 subdirectories
│       └── late/                # 9 subdirectories
│
├── scripts/
│   ├── 00_setup_and_verify.sh
│   ├── 01_extract_matrices.sb
│   ├── 02_run_tadcompare.R
│   ├── 02_run_tadcompare.sb
│   ├── 03_run_consensus.R
│   ├── 03_run_consensus.sb
│   ├── 04_analyze_shifts.R
│   ├── 04_postprocess.sb
│   ├── 05_visualizations.sb
│   ├── run_early_pipeline.sb
│   ├── run_pipeline.sh
│   └── tad_visualizations.R
│
├── logs/                        # SLURM output logs
│
├── CLAUDE.md                    # AI assistant context
├── TadCompare-CONTEXT.md        # Detailed analysis context
├── TadCompare.md                # TADCompare package documentation
├── Input_Data.Rmd               # Input format documentation
└── README.md                    # This file
```

---

## Configuration

### Resolution and Chromosome Settings

All scripts use consistent settings defined at the top:

```r
RESOLUTION <- 25000              # 25kb resolution
CHROMOSOMES <- paste0("chr", 1:19)  # Autosomes only, excludes chrX
```

### Base Directories (HPC - Expanse)

```r
BASE_DIR <- "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/tads"
```

### Input .hic Files

```
/expanse/lustre/projects/csd940/zalibhai/stripes/StripeCaller/data/hic/{timepoint}/
├── ctrl_merged.hic, mut_merged.hic          # Merged matrices
├── ctrl_M1.hic, ctrl_M2.hic, ctrl_M3.hic    # Control replicates
└── mut_M1.hic, mut_M2.hic, mut_M3.hic       # Mutant replicates
```

### ChIP-seq Peak Files

Located in `../peaks/`:

| Timepoint | H3K27ac | H3K27me3 | H3K4me1 |
|-----------|---------|----------|---------|
| Early | P12_ctrl_H3K27ac_early_peaks.bed | P12_ctrl_H3K27me3_early_peaks.bed | N/A |
| Late | 220310index25H3K27acLatePeakRegions.bed | 220310index29H3K27me3LatePeakRegions.bed | K4me1_aligned_reads_peaks.broadPeak-filtered.bed |

### SLURM Configuration

```bash
#SBATCH --account=csd940
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
```

---

## Quality Control Checkpoints

### After Step 2 (TADCompare)

- [ ] All 19 chromosomes processed successfully
- [ ] 8,000-10,000 total boundaries detected
- [ ] 14-17% differential boundaries (expected range)
- [ ] Boundary types distributed across all categories

### After Step 3 (ConsensusTADs)

- [ ] Robustness annotations added to all boundaries
- [ ] 20-25% boundaries detected in "Both" control and mutant
- [ ] High-confidence differential subset generated

### After Step 4 (Shift Distances)

- [ ] Shift distances calculated for all "Shifted" boundaries
- [ ] Median shift ~25-50kb (1-2 resolution bins)
- [ ] Final annotated file contains all expected columns

### After Step 5 (Visualizations)

- [ ] All 9 subdirectories contain expected plots
- [ ] GO/KEGG enrichment completed (or noted if insufficient genes)
- [ ] ChIP-seq anchor classification present
- [ ] Summary statistics file generated

---

## Expected Results

### Late Timepoint (250402)

| Metric | Value |
|--------|-------|
| Total boundaries | 8,264 |
| Differential | 1,384 (16.7%) |
| Shifted | 254 (3.1%) |
| Strength Change | 714 (8.7%) |
| Complex | 363 (4.4%) |
| Merge | 193 (2.3%) |
| Split | 112 (1.4%) |
| Non-Differential | 6,607 (79.9%) |
| Control-enriched | 532 (38.4%) |
| Mutant-enriched | 852 (61.6%) |
| Median shift distance | 0 kb (most not shifted) |
| Mean shift distance | 18.5 kb |

### Early Timepoint (250831)

| Metric | Value |
|--------|-------|
| Total boundaries | 8,584 |
| Differential | 1,222 (14.2%) |
| Shifted | 278 (3.2%) |
| Strength Change | 620 (7.2%) |
| Complex | 282 (3.3%) |
| Split | 193 (2.3%) |
| Merge | 123 (1.4%) |
| Non-Differential | 7,065 (82.5%) |
| Control-enriched | 784 (64.2%) |
| Mutant-enriched | 438 (35.8%) |

### Key Biological Observations

1. **Mutant-enriched boundaries dominate in late timepoint**: 61.6% of differential boundaries are stronger in mutant, suggesting increased boundary formation or strengthening with BAP1 loss.

2. **Control-enriched boundaries dominate in early timepoint**: 64.2% of differential boundaries are stronger in control, indicating potential developmental timing effects.

3. **Consistent differential rates**: ~15-17% of all boundaries show differential behavior, indicating measurable but not dramatic chromatin organization changes.

4. **syt1/nav3 locus (chr10)**: Previously identified as high-impact region, shows consistent boundary changes.

---

## Troubleshooting

### "ERROR: Control file not found"

**Possible causes:**
- Matrix extraction incomplete or failed
- Wrong timepoint specified
- Incorrect directory structure

**Solutions:**
- Re-run Step 1 extraction
- Verify files exist in `data/{timepoint}/extracted/merged/`
- Check logs for extraction errors

### "Insufficient contacts" warnings

**Possible causes:**
- Low sequencing depth for specific chromosome
- Sparse matrix regions (telomeres, centromeres)

**Solutions:**
- Usually safe to ignore for small chromosomes
- Check .hic file quality for specific samples
- Verify matrix file sizes are non-zero

### ConsensusTADs fails with < 2 replicates

**Possible causes:**
- Missing replicate matrix files
- Replicate extraction incomplete

**Solutions:**
- Verify all 6 replicate files exist per chromosome
- Re-run extraction for missing samples
- Check for R package version compatibility

### No significant enrichment (GO/KEGG)

**Possible causes:**
- Insufficient genes near differential boundaries
- Conservative FDR thresholds

**Solutions:**
- Check `boundary_genes.tsv` for gene counts
- Reduce minimum gene threshold in visualization script
- Use relaxed FDR for exploratory analysis

### Memory errors in TADCompare

**Possible causes:**
- Large chromosomes with dense contact matrices
- Insufficient allocated memory

**Solutions:**
- Increase SLURM memory request to 64GB
- TADCompare converts sparse to dense internally - this is memory-intensive
- Process chromosomes in smaller batches if needed

---

## Dependencies

### R Packages

```r
# Bioconductor
BiocManager::install(c(
  "TADCompare",          # Differential TAD boundary detection
  "GenomicRanges",       # Genomic interval operations
  "rtracklayer",         # BED file I/O
  "TxDb.Mmusculus.UCSC.mm10.knownGene",  # Gene annotations
  "org.Mm.eg.db",        # Mouse gene IDs
  "clusterProfiler",     # GO/KEGG enrichment
  "enrichplot"           # Enrichment visualization
))

# CRAN
install.packages(c(
  "tidyverse",
  "ggplot2",
  "patchwork",
  "RColorBrewer",
  "scales"
))
```

### System Tools

- `straw` (Juicer tool for .hic matrix extraction)
- `bedtools` (genomic distance calculations)
- conda environment: `mariner_env`

### SLURM Configuration

- Account: csd940
- Partition: shared
- Typical resources: 4 CPUs, 16-32GB RAM, 2-4 hours

---

## Important Constraints

1. **Exclude chrX**: All chromosome lists exclude sex chromosomes (autosomes 1-19 only)

2. **Resolution is 25kb**: Standard for TAD-level analysis, matches TADCompare recommendations

3. **Matrix format**: TADCompare accepts sparse 3-column format (region_i, region_j, counts) and converts internally

4. **No pre_tads argument**: TADCompare detects boundaries internally; existing Homer TAD files are cross-condition and don't match expected format

5. **Shift distance requires post-processing**: TADCompare classifies "Shifted" boundaries but doesn't calculate distances; we use bedtools closest

---

## Logs Reference

SLURM logs are stored in `logs/`:

| Log File | Pipeline Step | Description |
|----------|---------------|-------------|
| `extract_*.out` | Step 1 | Matrix extraction from .hic files |
| `tadcompare_*.out` | Step 2 | TADCompare differential analysis |
| `consensus_*.out` | Step 3 | ConsensusTADs robustness |
| `postprocess_*.out` | Step 4 | Shift distance calculations |
| `visualizations_*.out` | Step 5 | Plot generation |
| `early_pipeline_*.out` | Combined | Early timepoint full pipeline |

Example log locations from actual runs:
- `tadcompare_45828115.out` - Late timepoint TADCompare
- `early_pipeline_45828862.out` - Early timepoint combined run
- `visualizations_45828727.out` - Late visualization generation

---

## Citation

**TADCompare:**
- Cresswell KG, Dozmorov MG (2020). TADCompare: An R Package for Differential and Temporal Analysis of Topologically Associated Domains. *Frontiers in Genetics*, 11:158.

**Key Methods:**
- Rao et al. (2014). A 3D Map of the Human Genome at Kilobase Resolution. *Cell*, 159(7):1665-1680.
- Dixon et al. (2012). Topological domains in mammalian genomes. *Nature*, 485(7398):376-380.

---

## Related Documentation

- [CLAUDE.md](CLAUDE.md) - Complete analytical context for AI assistance
- [TadCompare-CONTEXT.md](TadCompare-CONTEXT.md) - Detailed experimental and analysis context
- [Main README](../README.md) - Loop analysis pipeline documentation
- [Stripes README](../stripes/README.md) - Stripe analysis pipeline documentation

---

**Last Updated:** January 2026
**Project:** BAP1-KO Differential TAD Boundary Analysis
**Organism:** Mouse (mm10)
