```yaml
# config/edgeR_config.yaml
# Configuration for Mariner differential loop analysis with edgeR

## Input/Output Paths
paths:
  base_dir: "/expanse/lustre/projects/csd940/zalibhai/mariner"
  
  # Input files from Mariner pipeline
  input:
    counts_matrix: "outputs/full/06_counts_matrix.rds"
    coordinates: "outputs/full/03_binned.rds"
    merged_loops: "outputs/full/02_merged.rds"
    qc_summary: "outputs/qc_report/qc_report_summary.rds"
  
  # Hiccups comparison files
  hiccups:
    ctrl_enriched: "/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/merged/diffhiccups_results/resorted_ctrl_vs_resorted_mut/file1/postprocessed_pixels_5000.bedpe"
    mut_enriched: "/expanse/lustre/projects/csd940/ctea/nf-hic/juicerpre/merged/diffhiccups_results/resorted_ctrl_vs_resorted_mut/file2/postprocessed_pixels_5000.bedpe"
  
  # Output directories
  output:
    base: "outputs/edgeR_results"
    primary: "outputs/edgeR_results/primary_analysis"
    comparison: "outputs/edgeR_results/hiccups_comparison"
    plots: "outputs/edgeR_results/plots"
    logs: "outputs/edgeR_results/logs"

## Sample Information
samples:
  names: ["ctrl", "mut"]
  group_labels: ["Control", "Mutant"]
  description: "BAP1 mutant vs control Hi-C merged replicates"

## Statistical Parameters
statistics:
  # Dispersion settings (no replicates, using fixed BCV)
  primary_bcv: 0.15      # Primary analysis (mouse controlled experiment)
  sensitivity_bcv:
    lower: 0.10          # Lower bound (near model organism expectation)
    upper: 0.20          # Upper bound (conservative estimate)
  
  # Filtering thresholds
  filtering:
    min_count: 5         # Minimum count threshold (lowered from default 10)
    min_total_count: 10  # Minimum total across samples (lowered from 15)
    min_prop: 0.7        # Minimum proportion for large groups
  
  # Significance thresholds
  fdr_primary: 0.05      # Primary FDR cutoff
  fdr_exploratory: 0.10  # Exploratory threshold for sensitivity
  
  # Effect size categories
  fold_change_thresholds:
    strong: 1.0          # |logFC| > 1 (2-fold)
    moderate: 0.5        # |logFC| > 0.5 (1.4-fold)
    weak: 0.0            # Any significant change
  
  # Normalization
  normalization_method: "TMM"  # Trimmed Mean of M-values

## Comparison Settings
comparison:
  # Hiccups matching parameters
  match_radius_kb: 10    # 10kb tolerance (matches Mariner merge radius)
  require_both_anchors: true
  
  # Overlap categories
  categories:
    - "shared_concordant"     # Both methods, same direction
    - "shared_discordant"     # Both methods, opposite direction
    - "mariner_specific"      # Only Mariner significant
    - "hiccups_specific"      # Only Hiccups significant
    - "non_differential"      # Neither method significant

## Visualization Settings
visualization:
  # Plot parameters
  ma_plot:
    point_size: 0.8
    alpha: 0.6
    fdr_threshold: 0.05
    fc_lines: [-1, 1]    # ±2-fold change reference lines
  
  volcano_plot:
    point_size: 0.8
    alpha: 0.6
    fdr_threshold: 0.05
    fc_threshold: 1.0
    label_top_n: 20      # Label top N significant loops
  
  # Color schemes
  colors:
    significant_up: "#d73027"      # Red for up in mutant
    significant_down: "#4575b4"    # Blue for down in mutant
    non_significant: "#999999"     # Gray for NS
    shifted_loops: "#fdae61"       # Orange for positionally shifted
    hiccups_overlap: "#1b7837"     # Green for Hiccups agreement
  
  # Output formats
  output_formats: ["pdf", "png"]
  dpi: 300
  width: 8
  height: 6

## Quality Control
qc:
  # Thresholds for warnings
  min_library_size: 1e6
  max_bcv: 1.0           # Warn if dispersions exceed 1
  min_filtered_loops: 10000  # Expect at least this many after filtering
  
  # Shifted loop analysis
  track_positional_shifts: true
  shift_enrichment_test: true

## Output Specifications
outputs:
  # File formats and naming
  results_formats:
    primary: ["tsv", "rds"]
    annotated: ["tsv", "xlsx"]
  
  # Results tables to generate
  tables:
    - name: "all_results"
      description: "Complete results for all tested loops"
      include_columns: ["loop_id", "chr1", "start1", "end1", "chr2", "start2", "end2", 
                       "logFC", "logCPM", "PValue", "FDR", "significant", "category", "shift_status"]
    
    - name: "significant_loops"
      description: "FDR < 0.05 loops only"
      filter: "FDR < 0.05"
    
    - name: "top_differential"
      description: "Top 100 by effect size"
      sort_by: "abs(logFC)"
      top_n: 100
  
  # Summary statistics to compute
  summary_stats:
    - "total_loops_tested"
    - "significant_loops_fdr05"
    - "significant_loops_fdr10"
    - "loops_up_in_mutant"
    - "loops_down_in_mutant"
    - "median_abs_logfc"
    - "strong_differential_count"
    - "moderate_differential_count"
    - "shifted_loops_tested"
    - "shifted_loops_significant"

## Runtime Configuration
runtime:
  parallel_cores: 8      # Matches SLURM --cpus-per-task
  memory_limit_gb: 32    # Matches SLURM --mem
  seed: 42              # For reproducibility
  verbose: true
  progress_updates: true

## Validation Checks
validation:
  # Pre-analysis checks
  check_input_files: true
  verify_coordinates: true
  check_count_integrity: true
  
  # Post-analysis checks
  verify_output_files: true
  check_result_validity: true
  generate_qc_report: true

## Metadata
metadata:
  pipeline_version: "1.0.0"
  edgeR_version: "4.x"
  analysis_date: "2025-10-20"
  analyst: "Zaki"
  project: "BAP1_mutant_Hi-C_differential_loops"
  description: |
    Differential chromatin loop analysis using Mariner buffer aggregation approach
    with edgeR statistical testing. Designed for merged replicate data (n=1 per condition)
    with comparison to Hiccups differential results for validation.
```

---

### Implementation Plan: Script 1 - Primary edgeR Analysis

**Script Purpose**: `01_run_edgeR_analysis.R`

**Workflow Structure** (following guidelines.md patterns):

1. **Configuration & Setup Module**
   - Load YAML config
   - Define all paths and constants
   - Set random seed
   - Load required libraries
   - Create output directories

2. **Data Loading & Validation Module**
   - Load count matrix (06_counts_matrix.rds)
   - Load coordinates (03_binned.rds)
   - Load QC summary for shift status
   - Extract genomic coordinates
   - Create annotation dataframe
   - Validate data integrity (no NAs, negatives, zeros)

3. **DGEList Construction Module**
   - Build count matrix with proper rownames
   - Create gene annotation with coordinates
   - Set group factor
   - Create DGEList object
   - Manual lib.size setting (optional based on config)

4. **Filtering Module**
   - Apply filterByExpr with custom thresholds
   - Track number filtered
   - Update library sizes
   - Log filtering statistics

5. **Normalization Module**
   - Apply TMM normalization
   - Calculate effective library sizes
   - Log normalization factors
   - Generate pre/post normalization plots

6. **Dispersion Analysis Module** (3-tier approach)
   - Primary: BCV = 0.15
   - Lower bound: BCV = 0.10
   - Upper bound: BCV = 0.20
   - Run exactTest for each
   - Track result convergence

7. **Results Extraction Module**
   - Extract topTags for each dispersion
   - Merge with coordinate annotations
   - Add shift status from QC
   - Categorize by effect size
   - Format output tables

8. **Statistical Summary Module**
   - Compute all summary statistics from config
   - Generate decideTests summaries
   - Calculate enrichments
   - Create summary report

9. **Visualization Module**
   - MA plot with multiple layers
   - Volcano plot with annotations
   - Dispersion sensitivity comparison
   - Distribution plots (logFC, logCPM)
   - Shifted loop enrichment plot

10. **Output Generation Module**
    - Save all result tables
    - Export plots in multiple formats
    - Generate summary text file
    - Save R objects for downstream use

**Key Features**:
- ✅ Comprehensive error handling at each step
- ✅ Progress logging with emoji milestones (matching clustering.py style)
- ✅ Defensive programming with input validation
- ✅ Modular design with clear function separation
- ✅ Process validation after each major step
- ✅ Complete output specification

**Expected Outputs**:
```
outputs/edgeR_results/
├── primary_analysis/
│   ├── edgeR_dge_object.rds
│   ├── all_results_bcv_0.15.tsv
│   ├── all_results_bcv_0.15.rds
│   ├── significant_loops_fdr05.tsv
│   ├── top100_differential.tsv
│   ├── sensitivity_comparison.tsv
│   └── summary_statistics.txt
├── plots/
│   ├── ma_plot_bcv_0.15.pdf
│   ├── volcano_plot_bcv_0.15.pdf
│   ├── dispersion_sensitivity.pdf
│   ├── logfc_distribution.pdf
│   ├── normalization_effects.pdf
│   └── shifted_loop_enrichment.pdf
└── logs/
    ├── analysis_log.txt
    └── session_info.txt
```

