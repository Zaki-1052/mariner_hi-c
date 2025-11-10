# Visualization Analysis: Publication-Quality Figures

## Script: `visualizations.R`

## Overview

This script generates **publication-ready visualizations** from differential chromatin loop analysis, including volcano plots, functional enrichment analyses, genomic feature distributions, and loop type classifications. It integrates results from all resolutions and produces figures suitable for manuscripts and presentations.

## Purpose

Transform statistical results into biological insights through:
1. **Volcano plots** - Visualize significance vs. effect size
2. **Feature annotations** - Genomic distribution of loop anchors
3. **Functional enrichment** - GO/KEGG pathway analysis
4. **Loop classification** - Distribution by regulatory type
5. **Comparative analyses** - Patterns across resolutions

## Input Data

### Required Files

**From Downstream Analysis:**
```
outputs/merged_loops/characterized_loops.tsv    # Main annotated dataset
outputs/merged_loops/non_redundant_loops.rds    # GInteractions object
```

**From edgeR Analysis (all resolutions):**
```
outputs/edgeR_results_res_5kb/primary_analysis/all_results_primary.tsv
outputs/edgeR_results_res_10kb/primary_analysis/all_results_primary.tsv
outputs/edgeR_results_res_25kb/primary_analysis/all_results_primary.tsv
```

**Optional (for complete merged volcano):**
```
outputs/merged_loops/merged_all_results.tsv    # All loops merged across resolutions
```

## Visualization Sections

---

## Section 1: Volcano Plots

### What is a Volcano Plot?

**X-axis:** log2 Fold Change (effect size)
**Y-axis:** -log10(FDR) (significance)

**Purpose:** Simultaneously display both statistical significance and biological magnitude

**Interpretation:**
- **Top corners:** Most significant differential loops
- **Far left/right:** Largest effect sizes
- **Center top:** Significant but modest changes
- **Bottom:** Non-significant loops

### Implementation

Uses **EnhancedVolcano** package for publication-quality plots with:
- Custom color schemes matching edgeR plots
- Automated significance threshold lines
- Count annotations for up/down-regulated loops
- Professional styling with clear legends

### Plots Generated

**Per-Resolution Volcano Plots:**
```
outputs/visualizations/volcano/volcano_5kb.pdf
outputs/visualizations/volcano/volcano_10kb.pdf
outputs/visualizations/volcano/volcano_25kb.pdf
```

Each shows differential loops at that resolution only.

**Multi-Resolution Merged Plot:**
```
outputs/visualizations/volcano/volcano_merged_multiresolution.pdf
```

Shows non-redundant union of all resolutions with proper overlap handling.

### Visual Elements

**Color Coding:**
- **Red/Dark Red:** Significant upregulated (mut > ctrl)
- **Blue/Dark Blue:** Significant downregulated (ctrl > mut)
- **Grey:** Not significant
- **Black:** Baseline

**Threshold Lines:**
- **Vertical:** |logFC| = 0.3 (fold change cutoff)
- **Horizontal:** FDR = 0.05 (significance cutoff)

**Annotations:**
- **Top left:** Count of upregulated loops (blue text)
- **Top right:** Count of downregulated loops (blue text)
- **Bottom right:** Total loops analyzed (black text)

### Example Interpretation

```
                Volcano Plot: WT vs KO Differential Loops

        ↑ Significance
    100 │           ·234        ·456
        │         ··  ··      ··  ··
     50 │       ··      ··  ··      ··
        │    ···          ··          ···
     10 │ ···                            ···
      5 │─────────────────────────────────────────→
        │                  ·····
      1 │         ············
        ├──────┼────────┼────────┼────────┼────────
       -2     -1       0       1        2      logFC

    Down in KO    ←   →    Up in KO
```

**Reading:**
- 234 loops significantly UP in KO (right corner)
- 456 loops significantly DOWN in KO (left corner)
- Most loops cluster near zero (no change)
- Clear separation validates biological effect

---

## Section 2: Genomic Feature Distribution

Uses **ChIPseeker** for annotating loop anchors relative to genomic features.

### Features Analyzed

**Gene-centric:**
- **Promoter** (≤1kb from TSS)
- **Promoter** (1-2kb from TSS)
- **5' UTR**
- **3' UTR**
- **Exon**
- **Intron**
- **Downstream** (≤3kb)

**Intergenic:**
- **Distal Intergenic** (>3kb from any gene)

### Plots Generated

```
outputs/visualizations/features/
├── anchor1_features.pdf       # First anchor distribution
├── anchor2_features.pdf       # Second anchor distribution
├── combined_features.pdf      # Both anchors aggregated
└── feature_enrichment.pdf     # Enrichment vs. genome background
```

### Feature Enrichment Analysis

**Method:** Compare observed vs. expected distribution

**Expected:** Genome-wide background frequency
**Observed:** Frequency in differential loops

**Enrichment score:** log2(observed / expected)

**Interpretation:**
- **Positive enrichment:** Feature over-represented in differential loops
- **Negative enrichment:** Feature under-represented
- **Near zero:** Proportional to genome background

**Example:**
```
Feature              Expected   Observed   Enrichment   p-value
─────────────────────────────────────────────────────────────
Promoter (≤1kb)        3.5%      15.2%      +2.12      < 0.001
Distal Intergenic     45.0%      28.3%      -0.67      < 0.01
Exon                   8.2%       6.8%      -0.27       0.12
```

**Biological insight:** Differential loops are enriched at promoters (regulatory changes) and depleted in intergenic regions (targeted effects).

---

## Section 3: Functional Enrichment Analysis

Uses **clusterProfiler** for Gene Ontology (GO) and KEGG pathway enrichment.

### Gene Set Preparation

**For each anchor:**
1. Extract nearest gene
2. Remove duplicates
3. Separate by direction (up vs. down regulated)

**Gene sets:**
- **Up-regulated loops:** Genes near loops with logFC > 0
- **Down-regulated loops:** Genes near loops with logFC < 0
- **All differential:** Combined set

### GO Enrichment

**Categories:**
- **Biological Process (BP):** Cellular functions
- **Molecular Function (MF):** Biochemical activities
- **Cellular Component (CC):** Subcellular localization

**Method:** Over-representation analysis (hypergeometric test)

**Correction:** Benjamini-Hochberg FDR

### KEGG Pathway Enrichment

**Database:** KEGG (Kyoto Encyclopedia of Genes and Genomes)

**Focus:** Signaling and metabolic pathways

### Plots Generated

```
outputs/visualizations/enrichment/
├── go_bp_enrichment_up.pdf         # GO BP for upregulated
├── go_bp_enrichment_down.pdf       # GO BP for downregulated
├── go_mf_enrichment_up.pdf         # GO MF for upregulated
├── go_mf_enrichment_down.pdf       # GO MF for downregulated
├── go_cc_enrichment_up.pdf         # GO CC for upregulated
├── go_cc_enrichment_down.pdf       # GO CC for downregulated
├── kegg_enrichment_up.pdf          # KEGG for upregulated
├── kegg_enrichment_down.pdf        # KEGG for downregulated
├── go_comparison_dotplot.pdf       # Compare up vs. down
└── enrichment_network.pdf          # Network of related terms
```

### Visualization Types

**Dot Plot:**
- X-axis: Gene ratio (fraction of genes in term)
- Y-axis: GO term / pathway
- Color: Adjusted p-value
- Size: Number of genes

**Bar Plot:**
- X-axis: Gene count
- Y-axis: GO term / pathway
- Color: Adjusted p-value

**Network Plot:**
- Nodes: Enriched terms
- Edges: Shared genes
- Clustering: Related biological processes

### Example Enrichment Results

**Upregulated Loops (BAP1-KO > WT):**
```
Term                                    GeneRatio   p.adjust   Genes
─────────────────────────────────────────────────────────────────────
Chromatin organization                  45/234      1.2e-12    ARID1A, SMARCA4, ...
Histone modification                    32/234      3.4e-10    KMT2A, EZH2, ...
DNA repair                              28/234      5.6e-08    BRCA1, RAD51, ...
```

**Downregulated Loops (WT > BAP1-KO):**
```
Term                                    GeneRatio   p.adjust   Genes
─────────────────────────────────────────────────────────────────────
Cell differentiation                    52/289      2.1e-15    SOX2, NANOG, ...
Transcription factor binding            38/289      4.5e-11    MYC, TP53, ...
Developmental process                   41/289      8.9e-10    NOTCH1, WNT3A, ...
```

**Biological interpretation:** BAP1 loss leads to upregulation of chromatin remodeling loops (compensatory?) and downregulation of differentiation loops (loss of cell identity?).

---

## Section 4: Loop Type Classification

Visualize distribution of loop types from ChIP-seq annotation.

### Loop Type Categories (10 total)

Based on anchor classifications:
1. **Promoter-Promoter (P-P)**
2. **Promoter-Active_Enhancer (P-AE)**
3. **Promoter-Poised_Enhancer (P-PE)**
4. **Promoter-Other (P-O)**
5. **Active_Enhancer-Active_Enhancer (AE-AE)**
6. **Active_Enhancer-Poised_Enhancer (AE-PE)**
7. **Active_Enhancer-Other (AE-O)**
8. **Poised_Enhancer-Poised_Enhancer (PE-PE)**
9. **Poised_Enhancer-Other (PE-O)**
10. **Other-Other (O-O)**

### Plots Generated

```
outputs/visualizations/loop_classification/
├── loop_type_distribution.pdf          # Overall counts
├── loop_type_by_direction.pdf          # Up vs. down
├── loop_type_by_resolution.pdf         # Across 5kb, 10kb, 25kb
├── loop_type_effect_size.pdf           # Box plot of logFC by type
└── loop_type_distance.pdf              # Distance distribution by type
```

### Key Visualizations

**1. Overall Distribution**
```
Bar plot showing count of each loop type:

P-AE     ████████████████████ 456 (34%)
AE-AE    ██████████████ 312 (23%)
P-P      ████████ 178 (13%)
P-PE     ██████ 134 (10%)
...
```

**2. Directional Distribution**
```
Stacked bar plot:
Up in KO  | Down in KO
──────────────────────
P-AE      █████│████
AE-AE     ███│█████
P-P       ██│███
```

Shows which loop types are preferentially up or down regulated.

**3. Effect Size by Type**
```
Box plot of logFC for each loop type:

         │     ┌───┐
P-AE     ├─────┤   ├─────  (median ≈ 0.5)
         │     └───┘

         │  ┌───┐
AE-AE    ├──┤   ├──────    (median ≈ 0.3)
         │  └───┘
```

Shows which types have stronger differential effects.

### Biological Insights from Loop Types

**High P-AE count:**
→ BAP1 affects enhancer-promoter communication

**High AE-AE count:**
→ BAP1 regulates super-enhancer organization

**High P-P count:**
→ BAP1 affects gene regulatory hubs

**Enrichment of P-PE in upregulated:**
→ Developmental genes becoming primed

**Enrichment of P-AE in downregulated:**
→ Loss of active enhancer connections

---

## Section 5: Distance and Spatial Analyses

### Loop Distance Distribution

```
outputs/visualizations/loop_classification/loop_distance_distribution.pdf
```

**Histogram of genomic distances:**
```
Count
  │
500 │  ██
    │  ████
400 │  ████
    │  ██████
300 │  ██████
    │  ████████
200 │  ██████████
    │  ████████████
100 │  ██████████████
    │  ████████████████
  0 │──────────────────────→ Distance (kb)
    0   100  200  300  400  500
```

**Categories:**
- **Short-range** (< 100kb): Local regulation
- **Medium-range** (100-500kb): Typical E-P distance
- **Long-range** (> 500kb): Inter-TAD, structural

**Analysis:**
- Median distance by loop type
- Distance vs. effect size correlation
- Resolution-specific distance preferences

### Distance vs. Effect Size

```
Scatter plot: logFC vs. loop distance

logFC
  2 │         ·    ·
    │       ·  ·   ·
  1 │    ··    ·  ·
    │  ··       ··
  0 │···       ···········
    │ ··        ··
 -1 │   ··    ·
    │     ·  ·
 -2 │       ·
    ├─────────────────────→
    0   200  400  600  800  Distance (kb)
```

**Interpretation:**
- **No correlation:** Effect independent of distance (good)
- **Positive correlation:** Long-range preferentially UP
- **Negative correlation:** Short-range preferentially DOWN

---

## Section 6: Multi-Resolution Comparison

### Resolution-Specific Effects

Identify loops unique to each resolution:
- **5kb-only:** Short-range, fine-scale
- **10kb-only:** Intermediate
- **25kb-only:** Long-range, structural
- **Shared:** Core regulatory loops

### Plots Generated

```
outputs/visualizations/resolution_comparison/
├── loops_by_resolution.pdf             # Venn diagram
├── effect_size_by_resolution.pdf       # Box plot of logFC
├── loop_type_by_resolution.pdf         # Type distribution
└── distance_by_resolution.pdf          # Distance distribution
```

### Concordance Analysis

**Fold Change Correlation:**
- Loops detected at multiple resolutions
- Compare logFC values
- Expected: High correlation (r > 0.8)

**Example:**
```
5kb logFC vs. 10kb logFC:

10kb
  2 │          ·
    │      ·   ·
  1 │    ·   ·
    │  ·   ·
  0 │···············
    │  ·   ·
 -1 │    ·   ·
    │      ·   ·
 -2 │          ·
    ├─────────────────→
   -2 -1  0  1  2  5kb

r = 0.92, p < 0.001
```

High correlation validates robustness across resolutions.

---

## Command-Line Usage

### Run All Visualizations (Default)

```bash
Rscript scripts/visualizations.R
```

Uses 5kb resolution for main volcano plot.

### Specify Volcano Resolution

```bash
Rscript scripts/visualizations.R --resolution 10000
```

Create volcano using 10kb dataset.

### Skip Time-Intensive Analyses

```bash
Rscript scripts/visualizations.R --skip-apa
```

Skips APA-related plots (useful if APA not yet run).

---

## Output Organization

```
outputs/visualizations/
├── volcano/
│   ├── volcano_5kb.pdf                    # 5kb resolution
│   ├── volcano_10kb.pdf                   # 10kb resolution
│   ├── volcano_25kb.pdf                   # 25kb resolution
│   └── volcano_merged_multiresolution.pdf # Non-redundant merged
│
├── features/
│   ├── anchor1_features.pdf               # First anchor annotation
│   ├── anchor2_features.pdf               # Second anchor annotation
│   ├── combined_features.pdf              # Aggregated
│   └── feature_enrichment.pdf             # Enrichment analysis
│
├── enrichment/
│   ├── go_bp_enrichment_up.pdf            # GO Biological Process (up)
│   ├── go_bp_enrichment_down.pdf          # GO Biological Process (down)
│   ├── go_mf_enrichment_up.pdf            # GO Molecular Function (up)
│   ├── go_mf_enrichment_down.pdf          # GO Molecular Function (down)
│   ├── go_cc_enrichment_up.pdf            # GO Cellular Component (up)
│   ├── go_cc_enrichment_down.pdf          # GO Cellular Component (down)
│   ├── kegg_enrichment_up.pdf             # KEGG pathways (up)
│   ├── kegg_enrichment_down.pdf           # KEGG pathways (down)
│   ├── go_comparison_dotplot.pdf          # Comparative analysis
│   └── enrichment_network.pdf             # Term relationships
│
└── loop_classification/
    ├── loop_type_distribution.pdf         # Overall counts
    ├── loop_type_by_direction.pdf         # Up vs. down
    ├── loop_type_by_resolution.pdf        # Across resolutions
    ├── loop_type_effect_size.pdf          # logFC by type
    └── loop_type_distance.pdf             # Distance by type
```

---

## Performance

**Expected Runtime:** 20-30 minutes

Breakdown:
- Volcano plots: 2-3 min
- Feature annotation: 5-8 min
- GO enrichment: 10-15 min (internet dependent)
- Loop classification: 2-3 min
- Multi-resolution: 3-5 min

**Memory Usage:** 4-8 GB RAM

**Disk Space:** ~200 MB (PDFs are high-resolution)

---

## Dependencies

### R Packages

**Core visualization:**
```r
library(ggplot2)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
```

**Genomic annotation:**
```r
library(GenomicRanges)
library(InteractionSet)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
```

**Enrichment analysis:**
```r
library(clusterProfiler)
library(enrichplot)
library(DOSE)
```

**Optional (for APA preview):**
```r
library(mariner)
library(strawr)
```

---

## Troubleshooting

### EnhancedVolcano Not Installed

```r
BiocManager::install("EnhancedVolcano")
```

### Enrichment Analysis Fails

**Issue:** No internet connection

**Solution:** GO/KEGG queries require internet. Run on connected system or skip enrichment section.

**Issue:** No significant terms

**Solution:**
- Check if enough differential loops
- Verify gene annotations are working
- Try less stringent FDR threshold

### ChIPseeker Errors

**Issue:** Different genome versions

**Solution:**
- Ensure using mm10 for all annotations
- Check TxDb version matches genome

### Too Many Loop Types

**Issue:** Hard to visualize 10 categories

**Solution:** Script automatically collapses rare types or uses grouped visualization. Edit color palette if needed:

```r
# In script, modify:
library(RColorBrewer)
my_colors <- colorRampPalette(brewer.pal(9, "Set1"))(10)
```

---

## Customization

### Modify Volcano Plot Appearance

Edit `create_publication_volcano()` function:

```r
# Change point size
pointSize = 2.5   # Make larger: 4.0

# Change colors
col = c('black', 'grey', 'red', 'darkred')
# To: c('grey70', 'grey50', 'blue', 'darkblue')

# Change thresholds
pCutoff = 0.05    # More stringent: 0.01
FCcutoff = 0.3    # More stringent: 0.5
```

### Add Custom Annotations

Add specific genes to volcano plot:

```r
# In create_publication_volcano(), add:
genes_of_interest <- c("Bap1", "Asxl1", "Brca1")

# Add labels
geom_text_repel(
  data = df[df$anchor1_nearest_gene %in% genes_of_interest, ],
  aes(label = anchor1_nearest_gene),
  size = 3
)
```

### Modify Enrichment Thresholds

```r
# Change p-value cutoff
enrichGO(..., pvalueCutoff = 0.05)  # Less stringent: 0.10

# Change q-value cutoff
enrichGO(..., qvalueCutoff = 0.2)   # Standard: 0.05

# Minimum gene set size
enrichGO(..., minGSSize = 10)       # Smaller: 5
```

---

## Interpretation Guidelines

### Volcano Plot

**Good patterns:**
- Clear separation of significant loops
- Balanced up/down (unless strong directional biology)
- Moderate number of significant loops (100-5,000)

**Concerning patterns:**
- All loops in one direction → batch effect
- Very few significant → low power
- Too many significant (>50%) → technical issue

### Enrichment Analysis

**Strong enrichment:**
- p.adjust < 0.01
- GeneRatio > 0.1 (>10% of genes)
- Biological coherence

**Validation:**
- Terms should make biological sense
- Consistent with known BAP1 function
- Supported by literature

### Loop Type Distribution

**Expected patterns:**
- P-AE most common (30-40%)
- AE-AE second (20-30%)
- P-P moderate (10-20%)
- Poised elements variable (0-15%)

**Unexpected patterns:**
- Very high O-O → poor ChIP-seq quality
- No P-AE → missing enhancer marks
- All one type → classification error

---

## Publication Checklist

Before using figures in manuscripts:

- [ ] All plots at 300+ dpi
- [ ] Consistent color schemes across figures
- [ ] Clear axis labels and legends
- [ ] Figure legends written (separate file)
- [ ] Statistics reported (n, p-values, FDR)
- [ ] Gene symbols verified (human: GENE, mouse: Gene)
- [ ] Supplementary tables prepared
- [ ] Data deposition planned (GEO, etc.)

---

**Last Updated:** November 10, 2025
