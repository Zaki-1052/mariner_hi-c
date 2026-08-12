## Section 61h: MeCP2-Up + K119ub-Up Quadrant GO Enrichment
**Key numbers:**
- Strict top-right quadrant (MeCP2 FDR < 0.05 + K119ub log2FC > 0) = 72 genes (source: 61h_mecp2_up_k119ub_up_genes.tsv)
- Strict + custom bg: 1 sig GO BP term — "basement membrane assembly", q = 0.021 (source: 61h_quadrant_go_bp.tsv)
- Relaxed (fold > 0 + K119ub > 0) + custom bg: 50 sig terms (source: 61h_relaxed_go_bp.tsv)
- Relaxed + genome-wide bg: 1,117 sig terms; top = synapse structure/organization (q = 6.95e-13), synapse assembly (q = 1.88e-12); 13/25 top terms neuronal (source: 61h_relaxed_genomewide_go_bp.tsv)

**What this shows:** The strict 72-gene quadrant set is too small/background-matched to yield neuronal enrichment (one non-neuronal term), but relaxing the MeCP2 FDR cut and using a genome-wide background recovers a strong synaptic/axon-guidance program. The neuronal signal is therefore real but background- and threshold-dependent — a caveat carried into 61i-61k.

**Figures:**
- 61h_relaxed_go_dotplot — relaxed quadrant GO dotplot (custom bg)
- 61h_relaxed_genomewide_dotplot — relaxed quadrant GO dotplot (genome-wide bg)
- (no strict dotplots: <3 significant terms)
