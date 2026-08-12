## Section 19: H3K27ac Peak Analysis at DMR Genes
**Key numbers:**
- Conditional gain among peak-bearing genes: mC Up = 31.5% (1,320/4,194), mC Down = 52.2% (335/642), hmC Up = 52.4%; background = 35.4% (source: h3k27ac_peak_analysis_summary.tsv)
- Fisher gain-vs-background: mC Down p = 3.20e-16 (enriched gain), hmC Down p = 0.317 (NS) (source: h3k27ac_peak_analysis_summary.tsv)
- Cliff's delta vs background: mC Down = +0.214 (small), mC Up = −0.052, hmC Up = +0.164 (source: h3k27ac_peak_analysis_summary.tsv)
- 4-mark O/E (mC-perspective enriched quadrants): ATAC 1.40/1.52, K119ub 1.30/1.55, H3K27ac 1.34/1.11, MeCP2 1.04/1.01 (source: h3k27ac_all_marks_oe_comparison.tsv)

**What this shows:** H3K27ac behaves opposite to the naive expectation — it is hypomethylated (mC Down) genes that significantly gain H3K27ac (OR>1, Cliff's delta +0.21), consistent with hypomethylated loci becoming more active enhancers, while hypermethylated genes show no significant H3K27ac change above background. In the 4-mark ranking, ATAC and K119ub remain the strongest directional couplings, H3K27ac intermediate, MeCP2 weakest.

**Figures:**
- `19a_h3k27ac_status_breakdown` — 100% stacked status bars
- `19b_h3k27ac_conditional_direction` — gain/equal/lost vs background + Fisher
- `19c_h3k27ac_effect_sizes` — per-gene net-peak strip/box
- `19d_methylation_vs_h3k27ac_scatter` — mod_difference vs net-H3K27ac (Spearman)
- `19e_h3k27ac_waterfall` — all genes ranked by net H3K27ac change
- `19f_h3k27ac_condition_overlap` — ctrl-vs-mut H3K27ac overlap at DMRs + Fisher
- `19g_h3k27ac_oe_heatmaps` — 2×2 mC × H3K27ac and hmC × H3K27ac O/E
- `19h_chromatin_mark_oe_comparison` — 4-mark O/E dot plot
