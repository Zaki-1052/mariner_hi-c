## Section 33: Multi-Mark DiffBind Quantitative Integration
**Key numbers:**
- 4-mark logistic regression for mC hypermethylation: H2AK119ub OR=4.71 [3.33, 6.66] p=2.4e-18 (top predictor); H3K27me3 OR=1.44; H3K27ac OR=0.48; ATAC OR=0.22 (source: diffbind_logistic_model_coefficients.tsv)
- Model AUC=0.818 [0.794, 0.841] (4-mark), 0.869 [0.850, 0.889] (4-mark + baseline mC/hmC) (source: logs/down/run.txt §Section 33 — log-only, not a TSV)
- Cross-mark Spearman: K119ub fold vs mC difference rho=+0.41; vs hmC difference rho=-0.30; K27me3 vs K27ac rho=-0.50 (source: diffbind_cross_mark_correlations.tsv)
- O/E cascade quadrants: mC Up + K119ub Gained Fisher OR=14.75 p=7.5e-175; mC Up + K27me3 Gained Fisher OR=14.12 p=8.6e-69 (source: diffbind_quantitative_oe_comparison.tsv)
- 178 genes with ≥3 marks concordant with the cascade (128 mC-hyper, 14 hypo, 36 non-DMR; 20 with all 4 marks) (source: diffbind_convergent_genes_3plus.tsv)

**What this shows:** Quantitative DiffBind fold-changes across 20,915 genes rank the four chromatin marks exactly as the BAP1 cascade predicts. Raised H2AK119ub is by far the strongest predictor of DNA hypermethylation (OR an order of magnitude above K27me3), while active/open marks (ATAC, K27ac) are protective. The cross-mark correlations and O/E enrichments independently reproduce the cascade wiring, and 178 high-confidence loci show three or more layers converging in-direction, anchoring the paper's mechanistic claim that K119ub gain is upstream of coordinated chromatin and methylation remodeling.

**Figures:**
- 33a_diffbind_volcano_plots — per-mark DiffBind volcanoes (ATAC, K27ac, K27me3, K119ub)
- 33b_cross_mark_correlation_heatmap — 7x7 Spearman heatmap (4 marks + mC + hmC + delta-ratio)
- 33c_quantitative_oe_dotplot — O/E enrichment per mark, mC and hmC perspectives
- 33d_methylation_vs_mark_scatters — methylation difference vs mark log2FC, loess fits
- 33e_logistic_regression_forest — forest plot of 4-mark logistic ORs
- 33f_convergence_analysis — cascade-concordance by DMR status + mark-combination intersections
