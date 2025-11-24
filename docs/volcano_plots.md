# Volcano Plots: TAD and Compartment Differential Visualization

## Scripts: `tad_volcano_plot.R` and `compartment_volcano_plot.R`

## Overview

These standalone scripts generate **publication-quality volcano plots** for TAD (Topologically Associating Domain) and chromatin compartment differential analyses from HOMER's `getDiffExpression.pl` output. They complement the main loop analysis by visualizing higher-order chromatin organization changes between BAP1-KO and wildtype conditions.

## Purpose

Visualize differential chromatin organization through:

1. **TAD Analysis** - Changes in TAD inclusion ratios (boundary strength)
2. **Compartment Analysis** - A/B compartment switching (active vs. inactive chromatin)
3. **Dual Threshold Output** - Both exploratory (relaxed) and publication (standard) versions
4. **Gene Annotations** - Link compartment changes to nearby genes

---

## What is a Volcano Plot?

### Concept

A volcano plot simultaneously displays:

**X-axis:** Effect size (Difference between conditions)
**Y-axis:** Statistical significance (-log10 adjusted p-value)

```
                    Volcano Plot Anatomy

    -log10(FDR)
        ↑
      5 │           ·                  ·
        │         ·  ·              ·   ·
      4 │        ·    ·            ·     ·
        │       ·      ·          ·       ·
      3 │────────────────────────────────────── FDR threshold
        │     ···        ·      ·        ···
      2 │    ·····        ····         ·····
        │   ·······      ··  ··       ·······
      1 │  ·········    ·    ·      ·········
        │ ···········  ·      ·    ···········
      0 │─────────────────────────────────────→ Difference
       -1.0    -0.5    0     0.5     1.0

           ← Decreased        Increased →
             in Mutant        in Mutant
```

### Interpretation Quadrants

| Region | Significance | Effect | Biological Meaning |
|--------|--------------|--------|-------------------|
| Top Right | High | Positive | Strong increase in mutant |
| Top Left | High | Negative | Strong decrease in mutant |
| Bottom Right | Low | Positive | Modest/uncertain increase |
| Bottom Left | Low | Negative | Modest/uncertain decrease |
| Center | Low | Near zero | No significant change |

---

## TAD Volcano Plot

### Script: `scripts/tad_volcano_plot.R`

### What is a TAD?

**Topologically Associating Domains (TADs)** are self-interacting genomic regions where chromatin interactions occur more frequently within the domain than across domain boundaries.

```
TAD Structure:
                 TAD Boundary           TAD Boundary
                     ↓                      ↓
Gene A  Gene B  Gene C  │  Gene D  Gene E  │  Gene F
──────────────────────────────────────────────────────
       ▓▓▓▓▓▓▓▓▓▓▓     │     ▓▓▓▓▓▓▓▓▓    │
          TAD 1        │       TAD 2      │   TAD 3
                       │                   │

Inclusion Ratio = intra-TAD contacts / inter-TAD contacts
```

**High Inclusion Ratio:** Strong TAD (well-insulated)
**Low Inclusion Ratio:** Weak TAD (permeable boundary)

### TAD Difference Interpretation

| Difference | Meaning | Biological Implication |
|------------|---------|----------------------|
| Positive | Higher IR in mutant | TAD boundaries strengthened |
| Negative | Lower IR in mutant | TAD boundaries weakened |

### Input Data

**Required file from HOMER getDiffExpression.pl:**
```
tad_analysis/Bap1.diff.tad.txt
```

**File structure:**
```
#TAD name    chr1   start1   end1   ...   ctrl vs. mut Difference   ctrl vs. mut adj. p-value
chr1:1000-5000   chr1   1000   5000   ...   0.148   0.234
chr2:8000-12000  chr2   8000   12000  ...   -0.085  0.609
...
```

**Key columns used:**
- Column 1: TAD name (chr:start-end format)
- `ctrl vs. mut Difference`: Effect size
- `ctrl vs. mut adj. p-value`: FDR-adjusted significance

### Command-Line Usage

```bash
# Basic usage (generates both relaxed and standard threshold plots)
Rscript scripts/tad_volcano_plot.R

# Specify custom input file
Rscript scripts/tad_volcano_plot.R tad_analysis/Bap1.diff.tad.txt

# Custom output directory
Rscript scripts/tad_volcano_plot.R --output outputs/custom_tad/

# Custom plot title
Rscript scripts/tad_volcano_plot.R --title "BAP1-KO TAD Differential Analysis"

# Custom dimensions
Rscript scripts/tad_volcano_plot.R --width 12 --height 10
```

**All options:**
```
Rscript scripts/tad_volcano_plot.R [INPUT_FILE] [OPTIONS]

Arguments:
  INPUT_FILE              Path to TAD differential file
                          (default: tad_analysis/Bap1.diff.tad.txt)
  --output DIR            Output directory (default: outputs/tad_analysis/)
  --title TEXT            Custom plot title
  --width WIDTH           Plot width in inches (default: 10)
  --height HEIGHT         Plot height in inches (default: 8)
```

### Threshold Definitions

The script automatically generates **two versions** of each plot:

| Version | FDR Cutoff | Difference Cutoff | Use Case |
|---------|------------|-------------------|----------|
| Relaxed | < 0.15 | \|Diff\| > 0.15 | Exploratory, hypothesis generation |
| Standard | < 0.05 | \|Diff\| > 0.30 | Publication, stringent filtering |

### Output Files

```
outputs/tad_analysis/
├── tad_volcano_relaxed.pdf         # Volcano plot (relaxed thresholds)
├── tad_volcano_standard.pdf        # Volcano plot (standard thresholds)
├── tad_significant_relaxed.tsv     # Significant TADs (relaxed)
├── tad_significant_standard.tsv    # Significant TADs (standard)
├── tad_volcano_summary.txt         # Statistics summary
└── tad_all_annotated.tsv           # Full dataset with annotations
```

### Visual Elements

**Color coding:**
- **Black:** Not significant (NS)
- **Grey:** Difference threshold only
- **Red:** FDR threshold only
- **Dark Red:** Both thresholds (most significant)

**Annotations:**
- **Top left:** Count of TADs with negative difference
- **Top right:** Count of TADs with positive difference
- **Bottom right:** Total TAD count

### Example Output

```
        TAD Differential Expression: Control vs Mutant (standard)

    -log10(FDR)
        ↑
      2.0│    156               234
         │     ·                  ·
      1.5│   ···                ···
         │   ····              ····
      1.3│─────────────────────────────  ← FDR = 0.05
         │  ·····            ·····
      1.0│ ······          ······
         │ ·······        ·······
      0.5│ ········      ········
         │·········    ·········
      0.0│──────────────────────────────→
        -0.5  -0.3   0    0.3   0.5   Difference
              ↑           ↑
              FC          FC
           threshold   threshold

                          total = 5074 TADs
```

---

## Compartment Volcano Plot

### Script: `scripts/compartment_volcano_plot.R`

### What are Chromatin Compartments?

Chromatin compartments represent the spatial organization of active (A) and inactive (B) chromatin:

```
Chromatin Compartment Organization:

A Compartment (Active)          B Compartment (Inactive)
─────────────────────          ─────────────────────────
• Euchromatin                  • Heterochromatin
• Gene-rich                    • Gene-poor
• Open chromatin               • Closed chromatin
• Active transcription         • Silenced genes
• H3K4me3, H3K27ac marks       • H3K9me3, H3K27me3 marks

PC1 Value Interpretation:
  Positive PC1 → A compartment (active)
  Negative PC1 → B compartment (inactive)
```

### Compartment Switching

**PC1 Difference interpretation:**

| Difference | Direction | Biological Meaning |
|------------|-----------|-------------------|
| Positive | B → A | Region becoming MORE active in mutant |
| Negative | A → B | Region becoming MORE inactive in mutant |

```
Compartment Switching Diagram:

                    ctrl              mut
                    ────              ────
B → A switch:       B(-)      →       A(+)    [Difference > 0]
                 (inactive)        (active)

A → B switch:       A(+)      →       B(-)    [Difference < 0]
                  (active)        (inactive)
```

### Input Data

**Required file from HOMER getDiffExpression.pl with annotatePeaks:**
```
tad_analysis/diffcompartments.txt
```

**File structure:**
```
PeakID    Chr   Start    End    ...   Gene Name   ...   ctrl vs. mut Difference   ctrl vs. mut adj. p-value
chr14-115525000   chr14   115525000   115550000   ...   4930505G20Rik   ...   -0.090   0.041
chr2-75475000     chr2    75475000    75500000    ...   9430019J16Rik   ...   0.052    0.342
...
```

**Key columns used:**
- Column 1: PeakID (chr-position format)
- `Chr`, `Start`, `End`: Genomic coordinates
- `Gene Name`: Nearest gene annotation
- `Annotation`: Genomic feature (promoter, intergenic, etc.)
- `ctrl vs. mut Difference`: PC1 difference (effect size)
- `ctrl vs. mut adj. p-value`: FDR-adjusted significance

### Command-Line Usage

```bash
# Basic usage (generates both threshold versions)
Rscript scripts/compartment_volcano_plot.R

# Specify custom input file
Rscript scripts/compartment_volcano_plot.R tad_analysis/diffcompartments.txt

# Custom output directory
Rscript scripts/compartment_volcano_plot.R --output outputs/custom_compartment/

# Label top genes on the plot
Rscript scripts/compartment_volcano_plot.R --label-genes --n-labels 20

# Custom dimensions
Rscript scripts/compartment_volcano_plot.R --width 14 --height 12

# Combine options
Rscript scripts/compartment_volcano_plot.R --label-genes --n-labels 15 --output outputs/pub_figures/
```

**All options:**
```
Rscript scripts/compartment_volcano_plot.R [INPUT_FILE] [OPTIONS]

Arguments:
  INPUT_FILE              Path to compartment differential file
                          (default: tad_analysis/diffcompartments.txt)
  --output DIR            Output directory (default: outputs/compartment_analysis/)
  --title TEXT            Custom plot title
  --width WIDTH           Plot width in inches (default: 12)
  --height HEIGHT         Plot height in inches (default: 10)
  --label-genes           Label top significant genes on plot
  --n-labels N            Number of genes to label (default: 10)
```

### Output Files

```
outputs/compartment_analysis/
├── compartment_volcano_relaxed.pdf         # Volcano plot (relaxed)
├── compartment_volcano_standard.pdf        # Volcano plot (standard)
├── compartment_significant_relaxed.tsv     # Significant regions (relaxed)
├── compartment_significant_standard.tsv    # Significant regions (standard)
├── compartment_volcano_summary.txt         # Statistics + gene summary
└── compartment_all_annotated.tsv           # Full dataset
```

### Visual Elements

**Color coding:**
- **Grey:** Not significant (NS)
- **Grey50:** Difference threshold only
- **Steel Blue:** FDR threshold only
- **Firebrick:** Both thresholds (most significant)

**Annotations:**
- **Top left:** "A→B: N" (regions becoming inactive)
- **Top right:** "B→A: N" (regions becoming active)
- **Bottom right:** Total region count

### Example Output

```
    Compartment Switching: Control vs Mutant (standard)

    -log10(FDR)
        ↑
      3.5│                          B→A: 1234
         │                              ·
      3.0│  A→B: 987                  ··
         │      ·                    ···
      2.5│     ··                  ····
         │    ···                 ·····
      2.0│   ····              ······
         │  ·····            ·······
      1.5│ ······          ········
         │ ·······        ·········
      1.3│─────────────────────────────  ← FDR = 0.05
         │·········    ···········
      1.0│··········  ·············
         │··········· ··············
      0.5│·············  ···············
         │···············  ···············
      0.0│──────────────────────────────→
        -2.0  -1.0    0     1.0    2.0   Difference

       More inactive          More active
         in mutant            in mutant

                        total = 101684 regions
```

---

## Data Files Reference

### TAD Analysis Directory

Located in `tad_analysis/`:

| File | Description | Source |
|------|-------------|--------|
| `Bap1.diff.tad.txt` | TAD differential expression | `getDiffExpression.pl` on TAD scores |
| `BAP1.tad.scores.txt` | Raw TAD inclusion ratios | `findTADsAndLoops.pl` |
| `diffcompartments.txt` | Compartment differential with gene annotations | `getDiffExpression.pl` with `annotatePeaks.pl` |
| `all_PC1.txt` | Raw PC1 values per sample | Hi-C compartment calling |

### Expected Data Sizes

| File | Approximate Entries |
|------|---------------------|
| TAD differential | ~5,000 TADs |
| Compartment differential | ~100,000 regions (25kb bins) |

---

## Biological Interpretation

### TAD Changes in BAP1-KO

**If TADs strengthen (positive difference):**
- Increased insulation between domains
- More restricted enhancer-promoter interactions
- Potential gene silencing

**If TADs weaken (negative difference):**
- Reduced domain insulation
- Increased inter-domain interactions
- Potential ectopic enhancer activation

### Compartment Changes in BAP1-KO

**B → A switching (positive difference):**
- Chromatin becoming more accessible
- Potential gene activation
- Often at developmental genes

**A → B switching (negative difference):**
- Chromatin becoming more compact
- Potential gene silencing
- May affect active regulatory regions

### Integration with Loop Analysis

TAD and compartment changes should be interpreted together with loop differential analysis:

```
Analysis Integration:

Loop Analysis         TAD Analysis         Compartment Analysis
─────────────         ────────────         ────────────────────
Specific E-P          Domain-level         Global chromatin
interactions          organization         state changes

     ↓                    ↓                      ↓
     └──────────────┬─────────────────────────────┘
                    ↓
            Integrated Model:
    How BAP1 loss affects 3D chromatin organization
```

---

## Troubleshooting

### Script Not Finding Input File

**Error:**
```
ERROR: Input file not found: tad_analysis/Bap1.diff.tad.txt
```

**Solution:**
- Verify file exists: `ls -la tad_analysis/`
- Run from project root directory
- Provide full path if needed

### Missing Required Columns

**Error:**
```
ERROR: Could not find 'ctrl vs. mut Difference' column
```

**Solution:**
- Check column names: `head -1 tad_analysis/Bap1.diff.tad.txt | tr '\t' '\n'`
- Ensure file is from `getDiffExpression.pl`
- Verify tab-delimited format

### EnhancedVolcano Not Installed

**Error:**
```
Error in library(EnhancedVolcano) : there is no package called 'EnhancedVolcano'
```

**Solution:**
```r
BiocManager::install("EnhancedVolcano")
```

### Very Few or No Significant Points

**Observation:** Volcano plot shows almost no colored points

**Possible causes:**
- Thresholds too stringent for your data
- Weak biological effect
- High biological variability

**Solutions:**
- Check the relaxed threshold plot first
- Review `*_volcano_summary.txt` for data ranges
- Consider biological context (is weak effect expected?)

### Memory Issues with Large Compartment Files

**Error:**
```
Error: cannot allocate vector of size X GB
```

**Solution:**
- Compartment files are large (~100k rows)
- Increase available memory
- Close other applications
- Run on HPC if needed

---

## Performance

### Expected Runtime

| Script | Data Size | Runtime |
|--------|-----------|---------|
| TAD volcano | ~5,000 TADs | 30-60 seconds |
| Compartment volcano | ~100,000 regions | 2-5 minutes |

### Memory Usage

| Script | Approximate Memory |
|--------|-------------------|
| TAD volcano | 1-2 GB |
| Compartment volcano | 4-8 GB |

---

## Dependencies

### R Packages Required

```r
# Core visualization
library(ggplot2)
library(tidyverse)
library(EnhancedVolcano)
```

### Installation

```r
# Install if needed
install.packages(c("ggplot2", "tidyverse"))
BiocManager::install("EnhancedVolcano")
```

---

## Customization

### Modify Thresholds

Edit the threshold definitions in the script:

```r
# In tad_volcano_plot.R or compartment_volcano_plot.R
thresholds <- list(
  relaxed = list(fdr = 0.15, fc = 0.15, name = "relaxed"),
  standard = list(fdr = 0.05, fc = 0.30, name = "standard"),
  # Add custom threshold:
  stringent = list(fdr = 0.01, fc = 0.50, name = "stringent")
)
```

### Change Color Scheme

```r
# In generate_volcano_plot() function
col = c('black', 'grey50', 'red', 'darkred')  # Default
# Change to:
col = c('grey80', 'grey50', 'steelblue', 'darkblue')  # Blue theme
```

### Adjust Point Size

```r
# For dense datasets (compartments)
pointSize = 1.5  # Smaller points

# For sparse datasets (TADs)
pointSize = 3.0  # Larger points
```

---

## Publication Checklist

Before using figures in manuscripts:

- [ ] Plots generated at appropriate resolution (PDF, 300+ dpi)
- [ ] Both relaxed and standard thresholds reviewed
- [ ] Threshold values clearly stated in figure legend
- [ ] Total region/TAD count visible on plot
- [ ] Significant counts annotated on plot
- [ ] Color scheme accessible (colorblind-friendly if needed)
- [ ] Biological interpretation documented
- [ ] Statistical methods described in Methods section

---

**Last Updated:** November 24, 2025
