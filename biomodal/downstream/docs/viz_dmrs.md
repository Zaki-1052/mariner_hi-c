```
(base) zakiralibhai@Zakirs-MacBook-Air modality % Rscript scripts/visualize_dmr_results.R
============================================
DMR Visualization Pipeline
============================================
Output directory: /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/downstream/modality/DMR_processed/summary/plots 

Loading summary data...
  Loaded 36 rows from summary file
Loading detailed comparison data...
  CG/genes: 8401 features
  CG/promoters: 210 features
  CG/cpg_islands: 112 features
  CG/cpg_shores: 2663 features
  CG/cpg_shelves: 3455 features
  CG/tss_region: 23 features
  Total comparison data: 14864 features

Generating Plot 1: Summary bar chart...
  Saved: 01_summary_significant_dmrs.pdf/png
Generating Plot 2: Hyper/Hypo direction chart...
Warning messages:
1: No shared levels found between `names(values)` of the manual scale and the data's fill values. 
2: No shared levels found between `names(values)` of the manual scale and the data's fill values. 
Warning messages:
1: No shared levels found between `names(values)` of the manual scale and the data's fill values. 
2: No shared levels found between `names(values)` of the manual scale and the data's fill values. 
  Saved: 02_hyper_hypo_direction.pdf/png
Generating Plot 3: mC vs hmC correlation scatter...
`geom_smooth()` using formula = 'y ~ x'
`geom_smooth()` using formula = 'y ~ x'
  Saved: 03_mc_hmc_correlation_scatter.pdf/png
  Correlation (all): r = -0.1818 
`geom_smooth()` using formula = 'y ~ x'
`geom_smooth()` using formula = 'y ~ x'
  Saved: 03b_mc_hmc_correlation_by_annotation.pdf/png
Generating Plot 4: Volcano plots...
  Saved: 04_volcano_plots_CG.pdf/png
Generating Plot 5: Heatmap...
  Saved: 05_heatmap_significant_dmrs.pdf/png
Generating Plot 6: Direction concordance...
  Saved: 06_direction_concordance.pdf/png

  Direction Concordance Summary (CG genes):
# A tibble: 4 × 2
  Pattern             Count
  <chr>               <int>
1 Both hyper            116
2 Both hypo             407
3 mC hyper + hmC hypo  6614
4 mC hypo + hmC hyper  1264
Generating Plot 7: Effect size distributions...
  Saved: 07_effect_size_distributions.pdf/png
Generating Plot 8: Context comparison...
  Saved: 08_context_comparison.pdf/png

============================================
Visualization Complete!
============================================

Plots saved to: /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/downstream/modality/DMR_processed/summary/plots 

=== Summary Statistics ===

CG Context (Primary Analysis):
  ModType Total_Significant Total_Hyper Total_Hypo
1     hmc             21477        4470      17007
2      mc             29871       21066       8805

Non-CG Contexts (minimal/no signal expected):
  Context ModType Total_Significant
1     CHG     hmc                24
2     CHG      mc               102
3     CHH     hmc                93
4     CHH      mc               389

Correlation (mC vs hmC) in CG genes:
  All genes: r = -0.6344 
  Significant in both: r = -0.6344 

=== Key Findings ===
1. CG context shows strong DMR signal (~4000+ significant)
2. CHG/CHH contexts show minimal/no significant DMRs (as expected)
3. Strong negative correlation between mC and hmC changes
4. Dominant pattern: mC hypermethylation + hmC hypomethylation in mutant

Done!
(base) zakiralibhai@Zakirs-MacBook-Air modality % 
```