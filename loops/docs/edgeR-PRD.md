```

Input Data (Mariner) → Primary edgeR Analysis → Results Formatting
                                    ↓
                    Secondary Comparison ← Hiccups Results
                                    ↓
                    Integrated Analysis & Visualization
```

### Critical Design Decisions

1. **Dispersion Strategy**: Use BCV = 0.15 for mouse controlled experiment (between technical and biological replicate expectations)
2. **Filtering Threshold**: Use `filterByExpr()` defaults but with `min.count=5` (not 10) to retain more loops
3. **FDR Cutoff**: Standard 0.05 for primary analysis, explore 0.10 for sensitivity analysis
4. **Comparison Framework**: Coordinate-based matching with 10kb tolerance (matching Mariner merge radius)

## Phase 1: Primary edgeR Analysis Implementation

### 1.1 Data Preparation Module

**Input Requirements**:
- Primary: `06_counts_matrix.rds` (22,108 × 2 matrix)
- Coordinates: `03_binned.rds` (GInteractions object for genomic positions)
- Metadata: `02_merged.rds` (for cluster membership if needed)

**Key Implementation Details**:
```
DGEList Construction:
- Extract genomic coordinates as annotation dataframe
- Format: chr1_start1_end1_chr2_start2_end2 for unique IDs
- Preserve original loop indices for back-mapping
- Manual lib.size setting if equal pairing needed (though data shows only 4% difference)
```

**Coordinate Extraction Strategy**:
- Use both anchor coordinates to create unique loop identifiers
- Maintain mapping between loop_1...loop_22108 generic IDs and actual positions
- Store width information (all should be 5001bp at 5kb resolution)

### 1.2 Statistical Analysis Pipeline

**Filtering Configuration**:
```
filterByExpr parameters:
- min.count = 5 (lowered from default 10)
- min.total.count = 10 (lowered from default 15)
- Keep loops expressed in at least 1 sample (given n=1 per group)
```

**Dispersion Estimation Strategy**:
Given no replicates, implement three-tier approach:
1. Primary analysis: BCV = 0.15 (reasonable for controlled mouse experiment)
2. Sensitivity bounds: BCV = 0.10 (lower bound) and 0.20 (upper bound)
3. Report convergence of results across dispersion values

**Statistical Testing Framework**:
```
Design matrix: ~group (intercept + treatment effect)
Test: coefficient 2 (mutant effect)
Method: exactTest with fixed dispersion (given no replicates)
Alternative: If using GLM, apply glmFit with prior.count=1
```

### 1.3 Results Processing Module

**Output Structure**:
```
Primary results table:
- loop_id (original index)
- genomic_coords (chr1:start1-end1_chr2:start2-end2)
- baseMean (average normalized counts)
- log2FoldChange (mut vs ctrl)
- pvalue (raw p-value)
- padj (FDR-adjusted)
- significant (FDR < 0.05)
- shift_status (from QC analysis - whether loop showed positional shift)
```

**Categorization Scheme**:
- Strong differential: |logFC| > 1 & FDR < 0.05
- Moderate differential: |logFC| > 0.5 & FDR < 0.05
- Weak differential: FDR < 0.05 & |logFC| ≤ 0.5
- Trending: FDR < 0.10 & FDR ≥ 0.05

## Phase 2: Hiccups Comparison Analysis

### 2.1 Data Integration Module

**Hiccups Data Processing**:
```
File mapping:
- file1/postprocessed_pixels_5000.bedpe → control-enriched (19,386 loops)
- file2/postprocessed_pixels_5000.bedpe → mutant-enriched (19,207 loops)

Coordinate standardization:
- Convert to GRanges/GInteractions
- Apply same seqlevels style as Mariner data
- Create unified coordinate system
```

**Matching Algorithm**:
```
Two-tier matching strategy:
1. Exact match: Same bin positions for both anchors
2. Fuzzy match: Within 10kb (2 bins) using mergePairs logic
   - Track match quality (exact vs fuzzy)
   - Record distance of match if fuzzy
```

### 2.2 Comparison Metrics Module

**Overlap Analysis Framework**:
```
Categories to track:
1. Shared differential (both methods agree on direction)
2. Shared but discordant (both significant, opposite directions)
3. Mariner-specific (captured by buffer approach)
4. Hiccups-specific (missed by aggregation)
5. Non-differential in both
```

**Quantitative Metrics**:
- Cohen's kappa for agreement
- Correlation of effect sizes for shared loops
- Enrichment of shifted loops in method-specific calls
- Spatial distribution analysis (are method-specific loops clustered?)

### 2.3 Biological Validation Module

**Shifted Loop Analysis**:
```
Key question: Do the 11,924 loops with positional shifts show different patterns?
- Compare differential rates in shifted vs non-shifted
- Assess if Mariner specifically captures shifted differential loops
- Calculate enrichment statistics
```

**Effect Size Comparison**:
- For matched loops, compare fold changes between methods
- Identify systematic biases (does aggregation dampen effects?)
- Flag loops with large discrepancies for manual review

## Phase 3: Visualization & Reporting

### 3.1 Core Visualizations

**Statistical Plots**:
1. **MA Plot** with overlays:
   - Color by significance status
   - Shape by shift status
   - Highlight Hiccups overlap

2. **Volcano Plot** with annotations:
   - Label top differential loops
   - Indicate method agreement
   - Show dispersion sensitivity bounds

3. **Method Comparison Plots**:
   - Venn diagram of differential calls
   - Scatter plot of fold changes (Mariner vs Hiccups)
   - Heatmap of top differential loops

### 3.2 Summary Statistics Module

**Key Metrics to Report**:
```
Overall statistics:
- Total loops tested (after filtering)
- Significant at FDR < 0.05
- Direction bias (up vs down in mutant)
- Median absolute fold change

Comparison statistics:
- Method agreement rate
- Unique contributions per method
- Shift-specific enrichments
```

### 3.3 Output Files Structure
```
outputs/
├── edgeR_results/
│   ├── differential_loops_fdr05.tsv
│   ├── all_results_annotated.tsv
│   ├── ma_plot.pdf
│   ├── volcano_plot.pdf
│   └── summary_stats.txt
├── comparison_results/
│   ├── method_overlap.tsv
│   ├── concordance_metrics.txt
│   ├── venn_diagram.pdf
│   └── fold_change_correlation.pdf
└── integrated_report.html
```
