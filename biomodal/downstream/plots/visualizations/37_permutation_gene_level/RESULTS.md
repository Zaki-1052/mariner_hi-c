## Section 37: Gene-Level Label-Shuffle Permutation Tests (Confirmatory)
**Key numbers:**
- 15 / 15 gene-level tests Confirmed (0 Weakened, 0 Strengthened, 0 Concordant-NS) (source: permutation_37_summary.tsv)
- Strongest: mC x expression direction perm z=-31.93 (OR=0.012, n=1,177); mC x ATAC perm z=-30.12 (OR=0.095); hmC x ATAC perm z=+28.89 (OR=10.48); mC x K119ub DiffBind perm z=+28.75 (OR=14.75) (source: permutation_37_summary.tsv)
- DiffBind cascade tests (validating 33c) all Confirmed: mC x ATAC OR=0.059, mC x K27ac OR=0.062, mC x K27me3 OR=14.12, mC x K119ub OR=14.75 (source: permutation_37_summary.tsv)
- Compartment-shift tests (validating 29): B→A x mC OR=0.030 z=-25.77; A→B x mC OR=2.73 z=+5.09 (source: permutation_37_summary.tsv)
- 13/15 tests hit the empirical-p floor (p≈9.999e-05 at 100,000 permutations) (source: permutation_37_summary.tsv)

**What this shows:** Chromosome-stratified label-shuffle permutation of 15 gene-level 2x2 contingency tables — shuffling categorical labels within chromosomes to preserve marginal counts and neighborhood structure — gives honest p-values for the Fisher/O-E enrichments in sections 12e, 15a-c, 19g, 20d, 27c, 29, 31b, and 33c. Every test passes, confirming the mC↔expression and mC↔ATAC anticorrelations, the hmC↔ATAC positive correlation, the four DiffBind cascade directions, the compartment-shift coupling, and the loop/MeCP2 anchor relationships are driven by genuine co-assignment at the same genes, not chromosomal co-residence. This is the only one of the four permutation sections (34-37) that ran to completion, so it carries the formal confirmatory weight.

**Figures:**
- 37a_zscore_forest_plot — forest plot of all 15 permutation z-scores, coloured by source section
- 37b_null_distribution_top4 — null-distribution histograms for the 4 strongest effects, observed OR marked
- 37c_fisher_vs_permutation_table — Fisher-OR vs permutation-z comparison table
- 37d_log2or_vs_zscore_scatter — observed log2(OR) vs permutation z-score scatter with trend line
