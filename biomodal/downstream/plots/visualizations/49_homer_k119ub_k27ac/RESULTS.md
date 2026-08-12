## Section 49: HOMER Motif Enrichment (K119ub + K27ac) — matrix dot plots + family bars
**Key numbers:**
- Significant motifs (q < 0.05, uncapped): B1 = 192, B2 = 3, B3 = 216, B4 = 15, C1 = 15, C2 = 60 (501 total motif-comparison rows) (homer_all_significant_motifs.tsv)
- K119ub matrix dot plot built from B1–B4 union (top-20 per comparison); K27ac matrix from C1–C2 (homer_top25_per_comparison.tsv)
- bHLH transcription factors dominate H2AK119ub-gained regions (script-confirmed; top hits Atoh1/Atoh7/NeuroD1/BHLHA15)

**What this shows:** Companion to section_49_homer_motif_enrichment, presenting the K119ub (B1–B4) and K27ac (C1–C2) differential-site motif enrichment as compareCluster-style matrix dot plots plus a TF-family bar chart faceted by mark. The gained-vs-lost contrasts (B1, B3, C1) carry the bHLH neuronal motif grammar; the lost-vs-gained contrasts (B2, B4) are nearly empty, confirming the enrichment is directional.

**NOTE:** This file's `comp_label` column contains literal newlines, so naive parsing of homer_all_significant_motifs.tsv splits rows; counts above were obtained by counting lines beginning with a valid comparison ID (B1–B4/C1–C4).

**Figures:**
- homer_k119ub_matrix_dotplot — B1–B4 matrix dot plot (union of top-20 TFs)
- homer_k27ac_matrix_dotplot — C1–C2 matrix dot plot
- homer_family_barchart_all — TF-family counts, all comparisons, faceted by mark
