```
(base) zakiralibhai@Zakirs-MacBook-Air downstream % Rscript /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/downstream/scripts/viz_sections/section_40_h3k36me2_me3_combined.R

================================================================================
BIOMODAL DUET evoC DIFFERENTIAL METHYLATION VISUALIZATION
================================================================================
Date: 2026-04-01 01:10:58 
Base directory: /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/downstream 
Output directory: /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/downstream/plots/visualizations 

Loading required packages...
Loaded multi_format_output.R - Functions available:
  - save_multiformat_ggplot(plot, base_path, width, height, dpi)
  - save_multiformat_base(plot_expr, base_path, width, height, dpi)
  - save_multiformat_pheatmap(pheatmap_call, base_path, width, height, dpi)
Packages loaded successfully.

================================================================================
LOADING DATA
================================================================================

Loading gene body mC DMRs...
  Loaded 20991 genes
Loading gene body hmC DMRs...
  Loaded 20991 genes
Loading Biological QC data...
  Loaded QC data successfully
Loading upstream sequencing metrics...
  Loaded 8 samples

Loading regional DMR data for comparison...
  Gene bodies: 10779 significant / 20991 total (51.4%)
  CpG islands: 446 significant / 8910 total (5.0%)
  CpG shores: 9837 significant / 32581 total (30.2%)
  CpG shelves: 6914 significant / 29094 total (23.8%)
  Promoters: 1700 significant / 20417 total (8.3%)
  TSS regions: 195 significant / 14165 total (1.4%)

Data loading complete.


================================================================================
SECTION 40: H3K36me2/me3 Combined Analysis
================================================================================

  Ratio table: 20915 genes
  me3 gene table loaded: 11952 genes
  me2 gene table loaded: 2909 genes
  Multi-mark data joined (4 marks)
  Genes with both me2 and me3: 1517
  Total combined profile: 20915 genes

--- FIGURE 40a: Expanded Correlation Heatmap ---

Warning messages:
1: In cor.test.default(cor_data[complete, i], cor_data[complete, j],  :
  Cannot compute exact p-value with ties
2: In cor.test.default(cor_data[complete, i], cor_data[complete, j],  :
  Cannot compute exact p-value with ties
3: In cor.test.default(cor_data[complete, i], cor_data[complete, j],  :
  Cannot compute exact p-value with ties
4: In cor.test.default(cor_data[complete, i], cor_data[complete, j],  :
  Cannot compute exact p-value with ties
5: In cor.test.default(cor_data[complete, i], cor_data[complete, j],  :
  Cannot compute exact p-value with ties
6: In cor.test.default(cor_data[complete, i], cor_data[complete, j],  :
  Cannot compute exact p-value with ties
7: In cor.test.default(cor_data[complete, i], cor_data[complete, j],  :
  Cannot compute exact p-value with ties
8: In cor.test.default(cor_data[complete, i], cor_data[complete, j],  :
  Cannot compute exact p-value with ties
  Saved: 40a_expanded_correlation_heatmap/{pdf,svg,jpg}
  Correlation heatmap saved

--- FIGURE 40b: me2 vs me3 Fold-Change Scatter ---

`geom_smooth()` using formula = 'y ~ x'
`geom_smooth()` using formula = 'y ~ x'
`geom_smooth()` using formula = 'y ~ x'
  Saved: 40b_me2_vs_me3_scatter/{pdf,svg,jpg}
  me2 vs me3 rho = -0.082, p = 1.37e-03
  Quadrants: Q1=458, Q2=426, Q3=276, Q4=357

--- FIGURE 40c: me2-me3 Ratio Delta by DMR Status ---

  Saved: 40c_me2_me3_ratio_delta/{pdf,svg,jpg}
Warning messages:
1: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on '(mC↑/hmC↓)' in 'mbcsToSbcs': for ↑ (U+2191)
2: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on '(mC↑/hmC↓)' in 'mbcsToSbcs': for ↓ (U+2193)

--- FIGURE 40d: Three-Way Venn ---

  Saved: 40d_three_way_venn/{pdf,svg,jpg}
  Triple overlap (me2 + me3 + mC): 55 genes
  me2 sig: 2909 | me3 sig: 333 | mC DMR: 10764

--- FIGURE 40e: Direction Flow (Grouped Bar) ---

  Direction patterns (triple-significant genes):
    me2↓ + me3↑ + mC↑: 19 genes (34.5%)
    me2↑ + me3↓ + mC↓: 11 genes (20.0%)
    me2↑ + me3↑ + mC↑: 8 genes (14.5%)
    me2↑ + me3↑ + mC↓: 7 genes (12.7%)
    me2↓ + me3↑ + mC↓: 5 genes (9.1%)
    me2↑ + me3↓ + mC↑: 4 genes (7.3%)
    me2↓ + me3↓ + mC↓: 1 genes (1.8%)
  Saved: 40e_direction_flow/{pdf,svg,jpg}
There were 21 warnings (use warnings() to see them)

--- FIGURE 40f: GO Enrichment Comparison ---

  me2-only: 2827 | me3-only: 251 | shared: 82
  Saved: 40f_go_comparison/{pdf,svg,jpg}

--- Exporting Tables ---

  Saved h3k36_combined_gene_profile.tsv (20915 rows)
  Saved h3k36_expanded_correlations.tsv

================================================================================
SECTION 40 COMPLETE
================================================================================
(base) zakiralibhai@Zakirs-MacBook-Air downstream %    
```