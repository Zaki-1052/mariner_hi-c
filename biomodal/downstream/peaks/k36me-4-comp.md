```
(base) zakiralibhai@Zakirs-MacBook-Air downstream % Rscript /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/downstream/scripts/viz_sections/section_41_dnmt3a_vs_dnmt3b_discrimination.R

================================================================================
BIOMODAL DUET evoC DIFFERENTIAL METHYLATION VISUALIZATION
================================================================================
Date: 2026-04-01 01:15:01 
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
SECTION 41: DNMT3A vs DNMT3B Mechanistic Discrimination
================================================================================

  Combined profile loaded: 20915 genes, 18 columns
  Loaded 4 missing columns from multi-mark table: atac_fdr, k27ac_fdr, k27me3_fdr, k119ub_fdr
  Hypermethylated genes: 7502
  K119ub gained genes: 3016
  me3 gained genes: 237 | me3 lost genes: 41

--- FIGURE 41a: me3 Fold at Hyper Genes with/without K119ub ---

  Group sizes:
    Hyper + K119ub Gained: 1312
    Hyper + No K119ub Gain: 3970
    Non-Hyper: 5019
    Hyper + K119ub Gained vs Hyper + No K119ub Gain: p = 1.23e-03
    Hyper + K119ub Gained vs Non-Hyper: p = 6.83e-01
    Hyper + No K119ub Gain vs Non-Hyper: p = 1.66e-05
  Saved: 41a_me3_at_hyper_k119ub_stratified/{pdf,svg,jpg}

--- FIGURE 41b: me2/me3 at Hyper Genes WITHOUT K119ub ---

  H3K36me2: Hyper(no K119ub) vs Non-Hyper p = 5.68e-20
  H3K36me3: Hyper(no K119ub) vs Non-Hyper p = 1.46e-02
  Saved: 41b_me2_me3_independent_pathway/{pdf,svg,jpg}

--- FIGURE 41c: 6-Mark Logistic Regression ---

  Complete cases for 6-mark model: 483 genes
  4-mark AUC: 0.869 [0.838, 0.900]
  6-mark AUC: 0.872 [0.842, 0.903]
  AUC improvement: +0.003
`height` was translated to `width`.
`height` was translated to `width`.
`height` was translated to `width`.
  Saved: 41c_logistic_regression_forest/{pdf,svg,jpg}
Warning messages:
1: `geom_errorbarh()` was deprecated in ggplot2 4.0.0.
ℹ Please use the `orientation` argument of `geom_errorbar()` instead. 
2: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'N=483 genes | 4-mark AUC=0.869 | 6-mark AUC=0.872 (Δ=+0.003)' in 'mbcsToSbcs': for Δ (U+0394)

--- FIGURE 41d: Hypermethylation Rate by Pathway Group ---

  Pathway group summary:
    K119ub Only: 1480/2483 hyper (59.6%)
    H3K36 Only: 236/817 hyper (28.9%)
    Convergent (K119ub + H3K36): 301/533 hyper (56.5%)
    Neither: 5485/17082 hyper (32.1%)
  Saved: 41d_pathway_hypermethylation_rate/{pdf,svg,jpg}

--- FIGURE 41e: me3 Status x K119ub Status Decision Matrix ---

  Saved: 41e_decision_matrix_heatmap/{pdf,svg,jpg}

--- FIGURE 41f: K119ub vs me3 Scatter ---

  Saved: 41f_k119ub_vs_me3_scatter/{pdf,svg,jpg}
Warning messages:
1: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'K119ub↑ + me3↓' in 'mbcsToSbcs': for ↑ (U+2191)
2: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'K119ub↑ + me3↓' in 'mbcsToSbcs': for ↓ (U+2193)
3: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'K119ub↑ + me3↑' in 'mbcsToSbcs': for ↑ (U+2191)
4: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'K119ub↑ + me3↑' in 'mbcsToSbcs': for ↑ (U+2191)
5: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'K119ub↓ + me3↑' in 'mbcsToSbcs': for ↓ (U+2193)
6: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on 'K119ub↓ + me3↑' in 'mbcsToSbcs': for ↑ (U+2191)

--- FIGURE 41g: Pathway Attribution ---

  Pathway attribution for hypermethylated genes:
    Unknown Mechanism: 5485 (73.1%)
    DNMT3A via H2AK119ub: 1480 (19.7%)
    Convergent (K119ub + me2): 268 (3.6%)
    DNMT3A via H3K36me2: 152 (2.0%)
    DNMT3B via H3K36me3: 84 (1.1%)
    Convergent (K119ub + me3): 33 (0.4%)
  Saved: 41g_pathway_attribution_pie/{pdf,svg,jpg}

--- FIGURE 41h: Top-50 Summary Heatmap ---

  Saved: 41h_top50_summary_heatmap/{pdf,svg,jpg}

--- Exporting Tables ---

  Saved dnmt3a_vs_dnmt3b_pathway_attribution.tsv (7502 rows)

================================================================================
SECTION 41 COMPLETE
================================================================================
(base) zakiralibhai@Zakirs-MacBook-Air downstream % 
```