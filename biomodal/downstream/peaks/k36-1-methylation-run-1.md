```
(base) zakiralibhai@Zakirs-MacBook-Air mariner_hi-c % cd downstream/
cd: no such file or directory: downstream/
(base) zakiralibhai@Zakirs-MacBook-Air mariner_hi-c % cd biomodal/downstream
(base) zakiralibhai@Zakirs-MacBook-Air downstream % ls
CLAUDE.md       docs            modality        plots           scripts
data            logs            peaks           qc              tables
(base) zakiralibhai@Zakirs-MacBook-Air downstream % Rscript /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/downstream/scripts/viz_sections/section_38_h3k36me3_gene_body_analysis.R

================================================================================
BIOMODAL DUET evoC DIFFERENTIAL METHYLATION VISUALIZATION
================================================================================
Date: 2026-04-01 00:59:09 
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
SECTION 38: H3K36me3 Gene Body Analysis (SETD2/DNMT3B Axis)
================================================================================

Loading H3K36me3 DiffBind data...
  H3K36me3 (all peaks): 39636 peaks (1465 sig up, 228 sig down at FDR<0.05)
  Ratio table: 20915 genes loaded
Loading H3K36me3 differential BEDs...
  H3K36me3 Up: 1465 peaks
  H3K36me3 Down: 228 peaks

Aggregating H3K36me3 peaks to gene level...
  11952 unique genes with H3K36me3 peaks
  10301 genes with both methylation and me3 data (49.3%)

--- FIGURE 38a: H3K36me3 DiffBind Overview ---

  Saved: 38a_h3k36me3_volcano_annotation/{pdf,svg,jpg}

--- FIGURE 38b: me3 x mC Direction O/E Heatmap ---

  Genes with both sig mC DMR and sig me3: 207
  Saved: 38b_me3_x_mc_oe_heatmap/{pdf,svg,jpg}
  Fisher OR = 4.30, p = 1.24e-04

--- FIGURE 38c: me3 x hmC Direction O/E Heatmap ---

  Genes with both sig hmC DMR and sig me3: 185
  Saved: 38c_me3_x_hmc_oe_heatmap/{pdf,svg,jpg}
  Fisher OR = 0.16, p = 1.47e-05

--- FIGURE 38d: me3 Fold vs Methylation Change Scatters ---

`geom_smooth()` using formula = 'y ~ x'
`geom_smooth()` using formula = 'y ~ x'
`geom_smooth()` using formula = 'y ~ x'
`geom_smooth()` using formula = 'y ~ x'
`geom_smooth()` using formula = 'y ~ x'
`geom_smooth()` using formula = 'y ~ x'
  Saved: 38d_me3_vs_methylation_scatter/{pdf,svg,jpg}
  mC correlation:  rho = +0.018, p = 6.71e-02
  hmC correlation: rho = -0.034, p = 5.63e-04

--- FIGURE 38e: me3 Fold at Coordinated Genes ---

  Group sizes:
    Coordinated (mC↑/hmC↓): 4786 genes
    Other DMR Genes: 3864 genes
    Non-DMR Genes: 1651 genes
  Saved: 38e_me3_fold_coordinated_violin/{pdf,svg,jpg}
Warning messages:
1: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on '(mC↑/hmC↓)' in 'mbcsToSbcs': for ↑ (U+2191)
2: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  conversion failure on '(mC↑/hmC↓)' in 'mbcsToSbcs': for ↓ (U+2193)
  Wilcoxon p = 3.73e-02 (coord vs non-DMR)
  Median me3 fold — Coordinated: +0.096, Non-DMR: +0.069

--- FIGURE 38f: DMR Overlap with me3 Up/Down BEDs ---

  Saved: 38f_dmr_me3_bed_overlap/{pdf,svg,jpg}
  Fisher OR = 6.10, p = 1.25e-13

--- FIGURE 38g: Chromosome Distribution ---

  Saved: 38g_chromosome_distribution/{pdf,svg,jpg}
  chr13: 87 / 1693 sig peaks (5.1%)

--- FIGURE 38h: me3 Fold by Chromatin State ---

  Saved: 38h_me3_fold_by_chromatin_state/{pdf,svg,jpg}
  Kruskal-Wallis p = 9.79e-02

--- Exporting Tables ---

  Saved h3k36me3_gene_level_summary.tsv (11952 rows)
  Saved h3k36me3_direction_concordance.tsv (2 rows)

================================================================================
SECTION 38 COMPLETE
================================================================================
(base) zakiralibhai@Zakirs-MacBook-Air downstream % Rscript /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/downstream/scripts/viz_sections/section_39_h3k36me2_boundary_analysis.R

================================================================================
BIOMODAL DUET evoC DIFFERENTIAL METHYLATION VISUALIZATION
================================================================================
Date: 2026-04-01 01:01:25 
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
>> preparing features information...             2026-04-01 01:01:36 
>> identifying nearest features...               2026-04-01 01:01:37 
>> calculating distance from peak to TSS...      2026-04-01 01:01:37 
>> assigning genomic annotation...               2026-04-01 01:01:37 
>> adding gene annotation...                     2026-04-01 01:01:44 
>> assigning chromosome lengths                  2026-04-01 01:01:44 
>> done...                                       2026-04-01 01:01:44 
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

>> preparing features information...             2026-04-01 01:01:46 
>> identifying nearest features...               2026-04-01 01:01:46 
>> calculating distance from peak to TSS...      2026-04-01 01:01:46 
>> assigning genomic annotation...               2026-04-01 01:01:46 
>> adding gene annotation...                     2026-04-01 01:01:47 
>> assigning chromosome lengths                  2026-04-01 01:01:47 
>> done...                                       2026-04-01 01:01:47 
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

Error in h(simpleError(msg, call)) : 
  error in evaluating the argument 'x' in selecting a method for function 'slice': Rle of type 'list' is not supported
Calls: geom_text ... Rle -> new_Rle -> .Call2 -> .handleSimpleError -> h
Execution halted
(base) zakiralibhai@Zakirs-MacBook-Air downstream % 
```