## Section 08: GO/KEGG Enrichment

**Key numbers:**
- Input = 7,513 hypermethylated mC genes (source: run-5 gene-body mC BED)
- Top GO BP = RNA splicing (274 genes, FoldEnrichment 2.42, p.adjust 6.27e-54) (source: `enrichment_go_bp.tsv`)
- Top GO CC = nuclear speck (256, 2.55e-49); top GO MF = ubiquitin-like protein transferase activity (269, 1.36e-52) (source: `enrichment_go_cc.tsv`, `enrichment_go_mf.tsv`)
- Top KEGG = Autophagy - animal (113, 8.82e-24); also Ubiquitin-mediated proteolysis, neurodegeneration pathways (source: `enrichment_kegg.tsv`)
- Delta-ratio deciles: D10 (most TET-impaired) = developmental (pattern specification, 124 genes, 6.91e-37); D1 (least impaired) = metabolic (oxidative phosphorylation, 46, 1.31e-14) (source: `enrichment_delta_ratio_top/bottom_decile_go_bp.tsv`)

**What this shows:** Hypermethylated genes are enriched for RNA-processing and ubiquitin-transferase/autophagy/neurodegeneration machinery — the ubiquitin axis thematically linking back to BAP1's deubiquitinase role. The genes most blocked from demethylation (D10) are developmental/patterning/cell-fate genes, tying the lesion to the cerebellar neuronal differentiation program.

**Figures:**
- `08a_enrichment_go_bp/`, `08b_enrichment_go_cc/`, `08c_enrichment_go_mf/`, `08d_enrichment_kegg/` — dotplots for hypermethylated genes.
- `08e_..compare_go_bp/`, `08f_..compare_kegg/` — compareCluster D10 vs D1 decile.
- `08g_..top_decile_go_bp/`, `08h_..bottom_decile_go_bp/` — standalone GO BP for D10 / D1.

> Note: `TODO.md` "RNA splicing 248 genes, q=3.4e-48" is OLD-run; current = 274 genes, p.adjust 6.27e-54.
