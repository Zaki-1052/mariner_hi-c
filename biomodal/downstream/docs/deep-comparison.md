```
(base) zakiralibhai@Zakirs-MacBook-Air biomodal %   cd downstream/ && Rscript scripts/viz_sections/compare_shallow_vs_deep.R

================================================================================
BIOMODAL DUET evoC DIFFERENTIAL METHYLATION VISUALIZATION
================================================================================
Date: 2026-03-30 23:24:27 
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
DEEP-4 (RUN-3, 4 samples) vs DEEP-8 (RUN-4, 8 samples) COMPARISON
================================================================================

Validating deep-4 (run-3) paths...
  All deep-4 paths validated.
Validating deep-8 (run-4) paths...
  All deep-8 paths validated.

Loading deep-4 and deep-8 DMR data...
  Deep-4 gene body: 20957 mC genes, 20957 hmC genes (deduplicated)
  Deep-8 gene body: 20969 mC genes, 20969 hmC genes (deduplicated)
  All regional data loaded.

--- Analysis A: DMR count comparison by region ---
  Saved: 20a_dmr_counts_by_region/{pdf,svg,jpg}
  Saved 20a_dmr_counts_by_region

--- Analysis B: Gene status classification ---

  5mC gene status:
    Deep only            17
    Lost significance    572
    Never significant    9615
    Newly significant    2511
    Retained             8254
    Shallow only         5

  5hmC gene status:
    Deep only            17
    Lost significance    482
    Never significant    9005
    Newly significant    2028
    Retained             9437
    Shallow only         5
  Saved: 20b_gene_status_summary/{pdf,svg,jpg}

  Saved 20b_gene_status_summary

--- Analysis C: Newly significant genes ---
  2511 newly significant mC genes
  2028 newly significant hmC genes
  Saved newly_significant_deep_seq.tsv (4539 rows)

--- Analysis D: Effect size scatter ---
  Saved: 20d_effect_size_scatter/{pdf,svg,jpg}
  5mC correlation: 0.9090 (20952 genes)
  5hmC correlation: 0.9234 (20952 genes)
  Saved 20d_effect_size_scatter

--- Analysis E: Q-value improvement ---
  Saved: 20e_qvalue_improvement/{pdf,svg,jpg}
  74.6% of genes have improved q-values in deep-8 vs deep-4
  Saved 20e_qvalue_improvement

--- Analysis F: Coverage comparison ---
  Saved: 20f_coverage_comparison/{pdf,svg,jpg}
  Median coverage: 18.7 (deep-4) vs 19.4 (deep-8)
  Saved 20f_coverage_comparison

--- Analysis G: Coordinated pattern stability ---
  Deep-4 coordinated genes: 5686 (84.6% of co-significant)
  Deep-8 coordinated genes: 6587 (78.8% of co-significant)
  Shared:       5429
  Deep-4 only:  257
  Deep-8 only:  1158
  Saved: 20g_coordinated_pattern_stability/{pdf,svg,jpg}
  Saved 20g_coordinated_pattern_stability

--- Creating master summary table ---
  Master table: 20974 genes
  Columns: gene, shallow_mc_q, deep_mc_q, shallow_mc_effect, deep_mc_effect, shallow_mc_coverage, deep_mc_coverage, mc_status, shallow_hmc_q, deep_hmc_q, shallow_hmc_effect, deep_hmc_effect, hmc_status, mc_effect_delta, hmc_effect_delta
  Saved shallow_vs_deep_comparison.tsv

================================================================================
COMPARISON COMPLETE
================================================================================

Output plots:   /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/downstream/plots/visualizations/comparison_shallow_vs_deep 
Output tables:  /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/downstream/plots/visualizations/tables 

Key findings:
  mC significant: 8826 (deep-4) -> 10765 (deep-8)
  hmC significant: 9919 (deep-4) -> 11465 (deep-8)
  Newly significant mC genes: 2511
  Newly significant hmC genes: 2028
  Effect size correlation: mC r=0.909, hmC r=0.923
  Coordinated pattern: 84.6% (deep-4) -> 78.8% (deep-8)
  Coverage: 18.7x (deep-4) -> 19.4x (deep-8) median

(base) zakiralibhai@Zakirs-MacBook-Air downstream % 
```