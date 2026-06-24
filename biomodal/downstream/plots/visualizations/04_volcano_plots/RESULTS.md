## Section 04: Volcano Plots

**Key numbers:**
- 5mC significant genes = 10,775 (q<0.05) (source: run-5 gene-body mC BED)
- 5hmC significant genes = 11,484 (q<0.05) (source: run-5 gene-body hmC BED)
- Top-gene q-value floor ≈ 4.82e-306; -log10(q) display-capped at 300 (source: `top_mc_dmrs.tsv`)
- 5mC cloud right-shifted (hyper); 5hmC cloud left-shifted (hypo); key genes labeled: Syt1, Zbtb20, Trpm3, Cntnap2 (source: run-5 BEDs / `KEY_GENES`)

**What this shows:** The 5mC and 5hmC volcanoes are mirror images — genes that gain 5mC lose 5hmC — visually encoding the coordinated demethylation-block model. The most extreme, machine-precision-significant genes are synaptic/neuronal.

**Figures:**
- `04a_volcano_mc/` — 5mC volcano (red hyper / blue hypo / grey NS).
- `04b_volcano_hmc/` — 5hmC volcano (dominant blue hypo cloud).
- `04_volcano_plots/` — combined side-by-side 5mC | 5hmC panel.

> Note: subtitle/Venn numbers in `FIGURES.md` (8,836 / 9,930) are stale OLD-run values; current run-5 = 10,775 / 11,484.
