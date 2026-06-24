## Section 72: K119ub Is Constitutively Enriched at Neuronal Genes (Absolute Signal, No MeCP2)
**Key numbers:**
- Neuronal enrichment among ctrl K119ub top decile: OR = 1.701 [1.544, 1.874], p = 3.36e-26 (q = 1.34e-25); top quartile OR = 1.378, p = 7.46e-19 (source: 72_fisher_results.tsv)
- Dose-response: D1 (lowest K119ub) = 20.87% neuronal → D10 (highest) = 32.99% neuronal; genome-wide fraction 23.50% (source: 72_neuronal_decile_summary.tsv)
- GSEA by absolute ctrl K119ub signal: 1,177 sig BP terms; neuron fate specification NES = +2.48 (q = 6.92e-17) (source: 72_gsea_ctrl_signal_go_bp.tsv)
- GO BP of ctrl-high K119ub genes: pattern specification process q = 1.57e-43; neuron fate commitment q = 5.81e-28 (source: 72_k119ub_topq_ctrl_go_bp.tsv)

**What this shows:** Using absolute (not differential) H2AK119ub gene-body signal, neuronal/axon-guidance genes (5,614-gene GO-derived set, 5,077 in the 21,604-gene quantifiable universe) are constitutively over-represented among the most K119ub-marked genes in wildtype, with an overall increasing (non-strict; D6 dips to 19.54%) dose-response (Spearman rho = 0.648, p = 0.049) and a Polycomb developmental-repression GO signature. Enrichment persists but weakens in mutant; the gained (mutant-only-high) set has no significant GO terms, so BAP1-KO amplifies signal at already-marked loci rather than recruiting a new functional class.

**Figures:**
- `72a_k119ub_signal_neuronal_vs_other/` — K119ub signal violins, neuronal vs non-neuronal (ctrl + mut)
- `72b_k119ub_high_neuronal_fisher/` — Fisher log2(OR) forest plot across 4 thresholds
- `72c_k119ub_topq_ctrl_go_dotplot/` — GO BP dotplot of ctrl top-quartile K119ub genes
- `72d_ctrl_vs_mut_k119ub_go_comparison/` — neuronal GO terms: ctrl-high vs mut-high vs gained
- `72e_gsea_ctrl_signal_dotplot/`, `72e_gsea_mut_signal_dotplot/` — GSEA by absolute signal
- `72e_running_score_neuron_fate_specification/` — GSEA running-score for top neuronal term
- `72f_neuronal_fraction_by_k119ub_decile/` — dose-response fraction + log2(OR) by decile
- `72_composite/` — assembled section-72 figure
