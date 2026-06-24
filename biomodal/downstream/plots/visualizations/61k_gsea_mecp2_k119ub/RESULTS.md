## Section 61k: GSEA (MeCP2 and K119ub Rankings) + 61jk Composites
**Key numbers:**
- MeCP2-ranked GSEA: 1 sig GO BP term — "synapse assembly", NES = +1.68, setSize = 299, q = 0.0245 (source: 61k_gsea_mecp2_go_bp.tsv)
- K119ub-ranked GSEA: 115 sig terms (61 positive NES, 54 negative NES) (source: 61k_gsea_k119ub_go_bp.tsv)
- K119ub top positive terms are RNA/tRNA metabolism, not synaptic: tRNA processing NES = +2.17 (q = 9.2e-06), RNA modification NES = +2.15 (source: 61k_gsea_k119ub_go_bp.tsv)
- K119ub neuronal sig terms = 3 (translation at presynapse NES = +1.94, q = 0.016) (source: 61k_gsea_k119ub_go_bp.tsv)

**What this shows:** Threshold-free GSEA gives a weaker neuronal signal than the peak-level ORA: only synapse assembly survives the MeCP2 ranking, and the K119ub ranking is dominated by RNA/tRNA-processing rather than synaptic terms. Neuronal involvement is real at the binding-site level (61j) but is not a fold-change-driven program — directly anticipating the section 78 neuronal-set-bias correction. Section 61jk renders the composite figures (no new statistics).

**Figures:**
- 61k_gsea_mecp2_dotplot, 61k_gsea_k119ub_dotplot — GSEA dotplots
- 61k_mecp2_running_synapse_assembly, 61k_k119ub_running_translation_at_presynapse — running-score plots
- 61k_composite — side-by-side NES bar composite (from section_61jk_composites.R)
