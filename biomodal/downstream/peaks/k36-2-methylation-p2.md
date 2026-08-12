```
(base) zakiralibhai@Zakirs-MacBook-Air downstream % Rscript scripts/viz_sections/section_39_h3k36me2_boundary_analysis.R

================================================================================
BIOMODAL DUET evoC DIFFERENTIAL METHYLATION VISUALIZATION
================================================================================
Date: 2026-04-01 01:06:08 
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
SECTION 39: H3K36me2 Polycomb Boundary Analysis (NSD/DNMT3A Axis)
================================================================================

Loading H3K36me2 DiffBind data...
  H3K36me2 (all peaks): 403123 peaks (4176 sig up, 2450 sig down at FDR<0.05)
  Ratio table: 20915 genes loaded
Loading H3K36me2 differential BEDs...
  H3K36me2 Up: 4176 peaks
  H3K36me2 Down: 2450 peaks
Loading H3K27me3 peaks for boundary analysis...
  H3K27me3: 15809 peaks
Loading H3K27me3 DiffBind data...
  H3K27me3 DiffBind: 18324 peaks (2293 sig up, 4811 sig down at FDR<0.05)

Annotating significant H3K36me2 peaks with ChIPseeker...
  6626 significant peaks to annotate
>> preparing features information...             2026-04-01 01:06:18 
>> identifying nearest features...               2026-04-01 01:06:19 
>> calculating distance from peak to TSS...      2026-04-01 01:06:19 
>> assigning genomic annotation...               2026-04-01 01:06:19 
>> adding gene annotation...                     2026-04-01 01:06:25 
>> assigning chromosome lengths                  2026-04-01 01:06:25 
>> done...                                       2026-04-01 01:06:25 
  Annotated 6626 peaks
  2909 unique genes with significant H3K36me2 peaks
  2178 genes with both methylation and me2 data (10.4%)

--- FIGURE 39a: H3K36me2 DiffBind Overview ---

  Saved: 39a_h3k36me2_volcano_annotation/{pdf,svg,jpg}

--- FIGURE 39b: me2 x mC Direction O/E Heatmap ---

  Genes with both sig mC DMR and sig me2: 1540
  Saved: 39b_me2_x_mc_oe_heatmap/{pdf,svg,jpg}
  Fisher OR = 0.37, p = 2.33e-19

--- FIGURE 39c: me2 Fold at H3K27me3 Boundaries ---

  me2 peaks at boundary: 660 (10.0%)
  Saved: 39c_me2_at_k27me3_boundary/{pdf,svg,jpg}
  Wilcoxon p = 2.97e-47
  Median fold — Boundary: +0.636, Away: +0.524

--- FIGURE 39d: me2 x K27me3 Direction O/E Heatmap ---

>> preparing features information...             2026-04-01 01:06:27 
>> identifying nearest features...               2026-04-01 01:06:27 
>> calculating distance from peak to TSS...      2026-04-01 01:06:27 
>> assigning genomic annotation...               2026-04-01 01:06:27 
>> adding gene annotation...                     2026-04-01 01:06:28 
>> assigning chromosome lengths                  2026-04-01 01:06:28 
>> done...                                       2026-04-01 01:06:28 
  Genes with both sig me2 and sig K27me3: 741
  Saved: 39d_me2_x_k27me3_oe_heatmap/{pdf,svg,jpg}
  Fisher OR = 0.05, p = 1.32e-63

--- FIGURE 39e: me2 Fold vs Methylation Change Scatters ---

`geom_smooth()` using formula = 'y ~ x'
`geom_smooth()` using formula = 'y ~ x'
`geom_smooth()` using formula = 'y ~ x'
`geom_smooth()` using formula = 'y ~ x'
`geom_smooth()` using formula = 'y ~ x'
`geom_smooth()` using formula = 'y ~ x'
  Saved: 39e_me2_vs_methylation_scatter/{pdf,svg,jpg}
  mC correlation:  rho = -0.185, p = 3.75e-18
  hmC correlation: rho = +0.320, p = 5.52e-53

--- FIGURE 39f: me2 Status by Chromatin State ---

  Saved: 39f_me2_by_chromatin_state/{pdf,svg,jpg}

--- FIGURE 39g: me2 Fold by Genomic Compartment ---

  Saved: 39g_me2_genic_vs_intergenic/{pdf,svg,jpg}
  Wilcoxon p = 4.10e-06

--- Exporting Tables ---

  Saved h3k36me2_gene_level_summary.tsv (2909 rows)
  Saved h3k36me2_k27me3_boundary_analysis.tsv

================================================================================
SECTION 39 COMPLETE
================================================================================
(base) zakiralibhai@Zakirs-MacBook-Air downstream % 
```