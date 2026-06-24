## Section 27: Methylation × Hi-C Loop Anchor Integration
**Key numbers:**
- Lost-vs-Gained anchor hypermethylation: Fisher OR = 2.589, p = 1.10e-35 (lost 46.86% vs gained 25.40% hyper) (methylation_direction_by_loop_direction.tsv)
- Logistic regression: loop_dir_binaryLost OR = 2.290 (p = 5.53e-32); Active_Enhancer anchor OR = 1.660, Polycomb anchor OR = 0.232 (logistic_regression_methylation_loop.tsv)
- Pooled anchor genes LESS coordinated than background: GREAT OR = 0.672, p = 5.33e-23 (43.37% vs 53.26%) (methylation_loop_anchor_coordinated_enrichment.tsv)
- K119ub-gain × hyper convergence OR = 1.791 (p = 2.05e-06); 119 triple-convergence genes (k119ub_loop_methylation_convergence.tsv)
- Shared anchors lower coordinated rate 34.32% vs non-shared 44.17% vs bg 53.26%; shared-vs-BG OR = 0.459 (shared_anchor_methylation_profile.tsv)

**What this shows:** Methylation tracks loop *direction*, not loop presence. Lost loops (down_in_mutant) mark TET-impaired hypermethylating genes; gained loops trend hypomethylated, so pooling cancels and anchor genes look depleted of coordination. Loop direction stays the dominant predictor of hypermethylation after adjusting for distance-to-TSS, loop distance, and anchor chromatin state, with active-enhancer anchors vulnerable and Polycomb/repressed anchors protected. K119ub gain converges with hypermethylation at lost-loop anchors, supporting a causal H2AK119ub→TET-block→5mC chain.

**Figures:**
- 27a_coordinated_enrichment_at_loop_anchors — % coordinated at anchors vs background (GREAT + nearest-gene)
- 27b_mc_direction_by_loop_direction — hypermethylation rate Lost/Gained/Background
- 27b_mc_diff_violin_by_loop_direction — mC difference by loop direction
- 27b_delta_ratio_violin_by_loop_direction — delta demethylation ratio by loop direction
- 27c_k119ub_methylation_loop_convergence_heatmap — K119ub × methylation O/E heatmap
- 27c_triple_convergence_summary — multi-layer convergence gene counts
- 27d_logistic_regression_forest_plot — odds ratios for predictors of hypermethylation
- 27d_linear_model_coefficients — delta-ratio linear-model coefficients
- 27e_shared_anchor_coordinated_rate — coordinated rate at shared vs non-shared anchors
- 27e_shared_anchor_delta_ratio_violin — delta-ratio across shared/non-shared/background
