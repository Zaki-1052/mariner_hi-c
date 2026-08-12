## Section 26: TET Triple-KO Comparison (Attenuation)
**Key numbers:**
- Absolute attenuation = 0.033 (BAP1-KO reproduces ~3.3% of the TET-KO ratio shift); relative attenuation = 0.087 (~8.7%) (tet_ko_comparison_summary.tsv)
- Median delta: BAP1-KO -0.0085 vs TET-KO -0.2597; KS D = 0.674 (p < 2.2e-16) (tet_ko_comparison_summary.tsv)
- Per-gene Spearman rho = 0.220 but residualized to 0.098; baseline accounts for 99.9% of explained variance; relative rho = -0.012 (n.s.) (tet_ko_comparison_summary.tsv)
- Response decomposition: TET-KO binary (68.5% complete loss + 27.0% no signal) vs BAP1-KO graded (0.8% strong, 45.1% moderate, 54.1% weak) over 20,851 matched genes (tet_ko_comparison_summary.tsv)

**What this shows:** BAP1-KO is a partial, graded attenuation of the TET triple-KO (GSE166423, Lopez-Moyado et al.) signature, recovering only ~3-9% of the effect. Where TET-KO drives all-or-nothing 5hmC collapse, BAP1-KO produces a continuous, mostly moderate dimming. The apparent per-gene correlation is almost entirely baseline-driven, so the two perturbations do not converge on the same genes — consistent with BAP1 impairing TET access/efficiency rather than eliminating TET enzyme.

**Figures:**
- 26a_wt_ratio_comparison — BAP1 WT vs TET-KO WT ratio densities.
- 26b_delta_ratio_density — delta densities, BAP1-KO vs TET-KO.
- 26c_qq_plot — absolute delta QQ (attenuation slope).
- 26d_effect_size_comparison — Cliff's delta / %dec / Cohen's d / relative attenuation.
- 26e_per_gene_scatter — per-gene BAP1 vs TET-KO delta.
- 26f_chromatin_stratified — per-state delta boxplots.
- 26g_relative_delta_density — baseline-normalized densities.
- 26h_relative_qq_plot — relative-delta QQ (step structure).
- 26i_response_decomposition — binary (TET) vs graded (BAP1) bars.
