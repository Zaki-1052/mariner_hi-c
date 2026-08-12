## Section 30: Polycomb Target Gene Enrichment
**Key numbers:**
- Polycomb targets DEPLETED from hypermethylation: Chromatin-State Polycomb × mC-hyper OR = 0.0633 (95% CI 0.0533-0.0748), p = 0, q = 0 (4.41% vs 42.14% hyper) (polycomb_fisher_tests.tsv)
- Polycomb ENRICHED in hypomethylation: × mC-hypo OR = 9.797 (8.996-10.673), p = 0, q = 0 (49.14% vs 8.98%) (polycomb_fisher_tests.tsv)
- Per-state hyper rate: Active_Promoter 71.77% (OR=10.00), Active_Enhancer 65.31% (OR=3.37) are targets; Polycomb 1.79% (OR=0.0323), Repressed_Promoter 2.99% (OR=0.0427) protected (polycomb_per_chromatin_state_enrichment.tsv)
- Robust across definitions: Strict × mC-hyper OR = 0.0422; H3K27me3-overlap OR = 0.1027 (both q=0) (polycomb_fisher_tests.tsv)

**What this shows:** Classic Polycomb/H3K27me3 genes are essentially excluded from de novo methylation and instead drift hypomethylated, exactly as the dual-mechanism model predicts — compact heterochromatin is inaccessible to DNMT3A. The hypermethylation burden falls on normally active genes (Active_Promoter ~72%, Active_Enhancer ~65%) that acquire ectopic H2AK119ub and lose TET activity. The effect holds across three independent Polycomb definitions and is reflected in smaller methylation magnitude at Polycomb genes (universe = 20,915 genes).

**Figures:**
- 30a_polycomb_vs_non_polycomb_stacked_bar — mC DMR status, Polycomb vs non-Polycomb
- 30b_fisher_forest_plot — ORs across Polycomb definitions and directions
- 30c_mc_magnitude_violin — |mc_diff| Polycomb vs non-Polycomb by direction
- 30d_per_state_hypermethylation_rate — % hyper per chromatin state vs genome-wide
- 30e_hmc_magnitude_violin — |hmc_diff| Polycomb vs non-Polycomb by direction
- 30f_composite_polycomb_summary — composite (30a+30b+30d)
