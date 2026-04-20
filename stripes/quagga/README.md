# Differential Chromatin Stripe Analysis Pipeline

## Overview

This pipeline performs **replicate-aware differential analysis of chromatin stripes** from Hi-C data, comparing **BAP1-KO mutant** and **wildtype control** conditions in mouse cerebellum (mm10 genome) with **biological replication** (n=3 per condition).

### Key Features

- **Quagga-based detection**: Uses Quagga (Feng et al., Genome Research 2025) for stripe calling
- **Anchor-based quantification**: Samples stripe signal at fixed distance into span direction
- **Multi-timepoint analysis**: Processes early (P12) and late developmental stages
- **Tiered confidence scoring**: Detection + replicate + 10kb validation for confidence levels
- **Direction classification**: Lost, gained, strengthened, weakened, unchanged categories
- **GO/KEGG enrichment**: Functional analysis of genes near differential stripe anchors
- **Publication-ready outputs**: Volcano plots, length distributions, ChIP-seq annotations

### Biological Context

**Chromatin stripes** are asymmetric 3D structures emanating from anchored points (CTCF sites or enhancer-promoter contacts) that extend into neighboring chromatin. They represent:

- **CTCF-stripes**: Result of cohesin-mediated loop extrusion from CTCF-bound anchors (~300kb median length)
- **EP-stripes**: Enhancer-promoter interactions without CTCF (~40-70kb median length)

**BAP1** is a Polycomb regulator (H2AK119ub1 deubiquitinase). Loss of BAP1 may affect:
- Polycomb-repressed domain boundaries
- Enhancer-promoter contacts at derepressed loci
- Emergence of ectopic EP-stripes at newly activated enhancers

---

## Pipeline Architecture

### Workflow Overview

```
Input: Quagga stripe calls (BEDPE) + Hi-C contact matrices (.hic)
  |
[1] phase1_detection.R     -> Union set creation, confidence scoring
  |
[2] phase2_quantification.R -> Extract 5x5 Hi-C matrices at sampling points
  |
[3] phase3_edgeR.R         -> Quasi-likelihood GLM differential analysis
  |
[4] phase4_integration.R   -> Merge annotations + statistics, classify direction
  |
[5] stripe_visualizations.R -> Volcano plots, length analysis, ChIP-seq annotation
  |
Output: Differential stripes with confidence tiers and biological annotations
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
cd /expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes
sbatch scripts/run_stripe_pipeline.sb
```

**Expected runtime:** ~4-6 hours (both timepoints, all phases + visualization)

**Output:** Differential stripes for early and late timepoints with all annotations and visualizations.

### Step-by-Step Execution

For manual control or debugging:

```bash
cd /expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes

# Process each timepoint
for TIMEPOINT in early late; do
  echo "Processing ${TIMEPOINT}..."

  # Phase 1: Detection & Union Set
  Rscript scripts/phase1_detection.R ${TIMEPOINT}

  # Phase 2: Hi-C Quantification
  Rscript scripts/phase2_quantification.R ${TIMEPOINT}

  # Phase 3: edgeR Differential Analysis
  Rscript scripts/phase3_edgeR.R ${TIMEPOINT}

  # Phase 4: Integration & Classification
  Rscript scripts/phase4_integration.R ${TIMEPOINT}
done

# Phase 5: Visualization (processes both timepoints)
Rscript scripts/stripe_visualizations.R
```

---

## Detailed Pipeline Steps

### Phase 1: Detection & Union Set Creation (`phase1_detection.R`)

**Purpose:** Create unified stripe set from merged BEDPE files with confidence annotations.

**Input:**
- Quagga BEDPE files: `data/quagga/{5kb,10kb}/{timepoint}/{sample}_stripes.bedpe`
- Merged files (pooled reads) + individual replicate files

**Process:**
1. Load merged 5kb stripes for control and mutant
2. Match anchors with 50kb tolerance
3. Classify each stripe: `control_only`, `mutant_only`, or `shared`
4. Handle direction switches as separate stripes
5. Annotate with replicate counts and 10kb validation

**Output:** `outputs/{timepoint}/`
- `01_unified_stripes.rds` - GRanges with all annotations
- `01_unified_stripes.tsv` - Tab-delimited for inspection

**Confidence Scoring:**
```
high:   n_reps >= 2 AND validated at 10kb
medium: n_reps >= 1 OR in_10kb OR pval < 1e-10
low:    detection only
```

---

### Phase 2: Hi-C Quantification (`phase2_quantification.R`)

**Purpose:** Extract Hi-C signal at stripe sampling points for differential analysis.

**Input:**
- Unified stripes from Phase 1
- 6 .hic files per timepoint

**Process:**
1. Convert stripes to GInteractions (anchor + sampling point)
2. Sampling point: 100kb into span direction from anchor center
3. Extract 5x5 matrices using mariner's `pullHicMatrices()`
4. KR normalization (fallback: VC)
5. Aggregate to count matrix (sum of 25 pixels per stripe)

**Output:** `outputs/{timepoint}/`
- `02_stripe_counts.rds` - Count matrix object
- `02_stripe_counts.tsv` - Stripes x samples matrix
- `02_extracted*.h5` - HDF5 matrices (if on_disk=true)

**Key Parameters:**
- `sampling_distance_bp: 100000` - 100kb into span
- `buffer_bins: 2` - 5x5 pixel extraction
- `normalization: "KR"` - Knight-Ruiz preferred

---

### Phase 3: edgeR Differential Analysis (`phase3_edgeR.R`)

**Purpose:** Identify stripes with significant differential contact frequency.

**Input:**
- Count matrix from Phase 2
- Sample metadata (condition assignment)

**Process:**
1. Create DGEList with count matrix
2. **Skip filtering** (retain all stripes due to small feature count)
3. TMM normalization
4. Estimate dispersion (robust)
5. Fit quasi-likelihood GLM
6. Test mutant effect with `glmQLFTest()`

**Output:** `outputs/{timepoint}/edgeR_results/`
- `03_all_results.tsv` - Full statistics (logFC, FDR, etc.)
- `03_dge_object.rds` - DGEList object
- `03_fit_object.rds` - QL fit object
- `plots/` - MDS, BCV, volcano, MA plots

**Statistical Framework:**
- **Model:** `~condition` (mutant vs control)
- **Test:** Quasi-likelihood F-test
- **Correction:** Benjamini-Hochberg FDR
- **Primary threshold:** FDR < 0.05

---

### Phase 4: Integration & Classification (`phase4_integration.R`)

**Purpose:** Merge Phase 1 annotations with Phase 3 statistics and classify stripe direction.

**Input:**
- Unified stripes from Phase 1
- edgeR results from Phase 3

**Process:**
1. Join detection annotations with differential statistics
2. Classify direction using tiered confidence approach
3. Flag directional consistency (logFC matches source)

**Direction Classification:**
```
lost:         source == "control_only" (by detection)
gained:       source == "mutant_only" (by detection)
strengthened: source == "shared" AND FDR < 0.05 AND logFC > 0.3
weakened:     source == "shared" AND FDR < 0.05 AND logFC < -0.3
unchanged:    source == "shared" (doesn't meet above)
```

**Confidence Tiers (for lost/gained):**
```
high:   FDR < 0.10 AND logFC in expected direction
medium: |logFC| > 0.2 in expected direction
low:    detection only
```

**Output:** `outputs/{timepoint}/`
- `04_final_differential_stripes.tsv` - Complete results
- `04_final_differential_stripes.bedpe` - For Juicebox
- `04_stripes_lost.tsv`, `04_stripes_gained.tsv` - By category
- `04_summary.txt` - Summary statistics

---

### Phase 5: Visualization (`stripe_visualizations.R`)

**Purpose:** Generate publication-quality figures, ChIP-seq annotations, and functional enrichment.

**Input:**
- Phase 4 results from both timepoints
- ChIP-seq peak files (H3K27ac, H3K27me3, H3K4me1)

**Analyses:**
1. **Volcano plots**: Per-timepoint and combined
2. **Length distribution**: CTCF vs EP classification
3. **ChIP-seq annotation**: Active, Polycomb, or indeterminate anchors
4. **GO/KEGG enrichment**: Functional analysis of genes near stripe anchors
5. **Summary statistics**: Confidence tiers, direction categories

**Output:** `outputs/visualizations/`
- `{timepoint}/volcano_{timepoint}.pdf`
- `{timepoint}/length_distribution_{timepoint}.pdf`
- `{timepoint}/anchor_classification_{timepoint}.pdf`
- `{timepoint}/{timepoint}_annotated_stripes.tsv`
- `enrichment/go_bp_dotplot_{timepoint}.pdf` (GO Biological Process)
- `enrichment/go_cc_dotplot_{timepoint}.pdf` (GO Cellular Component)
- `enrichment/go_mf_dotplot_{timepoint}.pdf` (GO Molecular Function)
- `enrichment/kegg_dotplot_{timepoint}.pdf` (KEGG Pathways)
- `enrichment/stripe_anchor_genes_{timepoint}.tsv` (genes near anchors)
- `combined/summary_statistics.txt`

**GO/KEGG Enrichment Details:**
- Extracts genes within 10kb of stripe anchors
- Groups genes by stripe direction (lost vs gained)
- Uses `clusterProfiler::compareCluster()` for enrichment
- Only includes medium/high confidence differential stripes
- Requires minimum 5 genes per category for analysis

---

## Output Directory Structure

```
stripes/
├── data/
│   └── quagga/                    # Quagga stripe calls
│       ├── 5kb/
│       │   ├── 250831/            # Early timepoint
│       │   └── 250402/            # Late timepoint
│       └── 10kb/                  # Validation resolution
│
├── outputs/
│   ├── early/                     # Early timepoint results
│   │   ├── 01_unified_stripes.*   # Phase 1: Detection
│   │   ├── 02_stripe_counts.*     # Phase 2: Quantification
│   │   ├── edgeR_results/         # Phase 3: Differential
│   │   │   ├── 03_all_results.tsv
│   │   │   ├── 03_dge_object.rds
│   │   │   └── plots/
│   │   ├── 04_final_differential_stripes.*  # Phase 4: Integration
│   │   ├── 04_stripes_lost.tsv
│   │   ├── 04_stripes_gained.tsv
│   │   └── 04_summary.txt
│   │
│   ├── late/                      # Late timepoint (same structure)
│   │
│   └── visualizations/            # Phase 5: Plots
│       ├── early/
│       ├── late/
│       ├── enrichment/            # GO/KEGG enrichment results
│       └── combined/
│
├── config/
│   └── stripe_config.yaml         # All parameters
│
├── scripts/
│   ├── phase1_detection.R
│   ├── phase2_quantification.R
│   ├── phase3_edgeR.R
│   ├── phase4_integration.R
│   ├── stripe_visualizations.R
│   ├── run_stripe_pipeline.sb     # Master SLURM script
│   └── *.sb                       # Individual phase wrappers
│
└── logs/                          # SLURM output logs
```

---

## Configuration

### Main Configuration File

`config/stripe_config.yaml` - Controls all analysis parameters.

**Key Settings:**

```yaml
# Detection parameters
detection:
  anchor_tolerance_bp: 50000    # 50kb anchor matching
  min_stripe_length: 200000     # 200kb minimum (Quagga level)

# Quantification parameters
quantification:
  sampling_distance_bp: 100000  # 100kb into span
  buffer_bins: 2                # 5x5 matrix
  normalization: "KR"           # Knight-Ruiz

# edgeR parameters
edger:
  skip_filtering: true          # Keep all stripes
  normalization_method: "TMM"
  fdr_primary: 0.05

# Classification thresholds
classification:
  logFC_threshold: 0.3          # For strengthened/weakened
  fdr_threshold: 0.05
  tiered:
    high_fdr: 0.10
    medium_logfc: 0.2

# Sample information
samples:
  names: ["ctrl_M1", "ctrl_M2", "ctrl_M3", "mut_M1", "mut_M2", "mut_M3"]
  groups: ["ctrl", "ctrl", "ctrl", "mut", "mut", "mut"]
```

### File Paths

Update these in `stripe_config.yaml` for your environment:
- `stripe_data.base_path` - Quagga BEDPE files
- `hic_files.base_path` - .hic contact matrices
- `outputs.base_dir` - Output directory
- `chipseq_peaks.*` - ChIP-seq BED files

---

## Quality Control Checkpoints

### After Phase 1 (Detection)

- [ ] Reasonable stripe counts (100-300 per merged file at 200kb min)
- [ ] Shared stripes present (indicates reproducible features)
- [ ] Confidence distribution balanced (not all "low")

### After Phase 2 (Quantification)

- [ ] Count matrix has no NA values
- [ ] Library sizes within 2-fold across samples
- [ ] Non-zero counts for majority of stripes

### After Phase 3 (edgeR)

- [ ] BCV (biological coefficient of variation) ~0.2-0.5
- [ ] MDS plot shows condition separation
- [ ] Some significant stripes at FDR < 0.10

### After Phase 4 (Integration)

- [ ] Direction classification distribution reasonable
- [ ] Directional consistency ~50% indicates noise (expected)
- [ ] Medium/high confidence stripes present

**Expected Results:**

| Metric | Early | Late |
|--------|-------|------|
| Total stripes | ~200-300 | ~150-250 |
| Lost (control_only) | ~80-130 | ~60-100 |
| Gained (mutant_only) | ~60-100 | ~50-80 |
| Shared | ~50-80 | ~30-60 |
| Significant (FDR<0.05) | ~0-10 | ~0-5 |

**Note:** Low significance rates are expected with ~300 features and high detection noise.

---

## Troubleshooting

### Low Stripe Counts (< 50)

**Possible causes:**
- Quagga parameters too stringent
- Hi-C data quality issues
- Wrong resolution selected

**Solutions:**
- Check `min_length` parameter (should be 200kb)
- Verify merged BEDPE files exist
- Review Quagga detection logs

### Phase 2 Memory Errors

**Possible causes:**
- Too many stripes for available RAM
- HDF5 backend not working

**Solutions:**
- Enable `on_disk: true` in config
- Reduce `block_size`
- Request more memory in SLURM

### No Significant Stripes

**Possible causes:**
- Weak biological effect
- High detection noise (~50% directional consistency)
- Small feature set limits power

**Solutions:**
- Focus on medium/high confidence tiers
- Use exploratory FDR (0.10)
- Interpret detection-based categories directly

### Visualization Fails

**Possible causes:**
- ChIP-seq peak files not found
- Missing Phase 4 outputs

**Solutions:**
- Verify `chipseq_peaks` paths in config
- Re-run Phase 4 for both timepoints
- Check R package dependencies (EnhancedVolcano, ChIPseeker)

---

## Key Differences from Loop Analysis

| Aspect | Loops | Stripes |
|--------|-------|---------|
| **Geometry** | Point-like (anchor-anchor) | Linear (anchor-span) |
| **Detection** | HiCCUPS | Quagga |
| **Quantification** | Loop center | 100kb into span |
| **Feature count** | ~10,000-30,000 | ~150-300 |
| **Filtering** | filterByExpr() | Skip (too few features) |
| **Classification** | Up/down regulated | Lost/gained + strengthened/weakened |
| **Confidence** | Replicate-based | Detection + replicate + 10kb validation |

### Why Anchor-Based Quantification?

Stripes are linear features extending from an anchor into a span region. Unlike loops (point interactions), stripes require:

1. **Direction awareness**: Span extends 3' or 5' from anchor
2. **Sampling strategy**: Fixed 100kb into span captures consistent signal
3. **Anchor matching**: 50kb tolerance for same-stripe identification

---

## Dependencies

### R Packages

```r
# Bioconductor
BiocManager::install(c(
  "mariner",           # Hi-C matrix extraction
  "edgeR",             # Differential analysis
  "InteractionSet",    # GInteractions
  "GenomicRanges",     # Genomic intervals
  "HDF5Array",         # Out-of-memory arrays
  "strawr",            # .hic file reading
  "ChIPseeker",        # Peak annotation
  "EnhancedVolcano",   # Volcano plots
  "clusterProfiler",   # GO/KEGG enrichment
  "enrichplot",        # Enrichment visualization
  "org.Mm.eg.db",      # Mouse gene annotations
  "TxDb.Mmusculus.UCSC.mm10.knownGene"  # Mouse transcript database
))

# CRAN
install.packages(c(
  "tidyverse",
  "yaml",
  "ggplot2",
  "patchwork"
))
```

### System Requirements

- R >= 4.1.0
- HDF5 library (for HDF5Array)
- 64GB RAM recommended
- 12+ CPUs for parallel extraction

---

## Citation

**Quagga:**
- Feng et al. (2025). *Quagga: A tool for chromatin stripe detection*. Genome Research.

**mariner:**
- Kramer et al. (2022). *mariner: Explore the Hi-Cs*. Bioconductor.

**edgeR:**
- Robinson et al. (2010). *edgeR: a Bioconductor package for differential expression analysis*. Bioinformatics.
- Chen et al. (2016). *From reads to genes to pathways*. F1000Research.

---

## Related Documentation

- [CLAUDE.md](CLAUDE.md) - Complete analytical context and design decisions
- [Main README](../README.md) - Loop analysis pipeline documentation
- [config/stripe_config.yaml](config/stripe_config.yaml) - Full parameter reference

---

**Last Updated:** December 2024
**Project:** BAP1-KO Differential Chromatin Stripe Analysis
**Organism:** Mouse (mm10)
