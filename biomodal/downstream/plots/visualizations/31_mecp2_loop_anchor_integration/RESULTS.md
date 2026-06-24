## Section 31: MeCP2 × Hi-C Loop Anchor Integration
**Key numbers:**
- MeCP2 significant-gain × loop direction Fisher OR = 0.424, p = 1.32e-03 (built Gained-vs-Lost, so <1 = MeCP2 gain favors Lost) (mecp2_loop_direction_fisher.tsv)
- MeCP2-Sig-Up enriched at Lost loops: Obs 37 vs Exp 24.5, enrichment 1.51; depleted at Gained: Obs 23 vs Exp 35.5, enrichment 0.65 (mecp2_loop_direction_fisher.tsv)
- Gene-level groups with MeCP2 data: Gained n=2,246, Lost n=1,553, Background (no loop) n=17,701 (mecp2_fold_at_loop_anchor_genes.tsv)
- 2,910 loops scored for anchor MeCP2 overlap; 2,645 loops have gene-level MeCP2 data for the loop-logFC scatter (mecp2_loop_anchor_overlap.tsv, mecp2_loop_gene_level_scatter_data.tsv)

**What this shows:** MeCP2 (CUT&RUN, 5mC-preferring reader) gain associates with the lost-loop anchors that hypermethylate in BAP1-KO, consistent with Mellen-et-al reader logic — more 5mC at disrupted anchors draws more MeCP2 (gain enrichment 1.51 at lost loops vs 0.65 at gained loops). This connects the chemical change (5hmC→5mC) to reader occupancy and to 3D architecture (loop loss). The 31a Fisher ORs, 31c Wilcoxon p-values, and 31d Spearman rho are console-only; persisted tables hold per-loop/per-gene data.

**Figures:**
- 31a_mecp2_peak_overlap_at_loop_anchors — % loops with anchor MeCP2 overlap by loop direction
- 31b_mecp2_loop_direction_fisher_heatmap — Obs/Exp/Enrichment, MeCP2-Sig-Up × loop direction
- 31c_mecp2_fold_by_loop_direction — MeCP2 fold change, Gained/Lost/Background
- 31d_loop_logfc_vs_mecp2_fold_scatter — loop logFC vs mean MeCP2 fold with gene labels
