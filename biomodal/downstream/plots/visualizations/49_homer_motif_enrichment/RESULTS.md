## Section 49: HOMER Motif Enrichment (K119ub + K27ac) — dotplot/family/wordcloud/MOI
**Key numbers:**
- Significant motifs (q < 0.05, uncapped): B1 (K119ub gained-vs-lost) = 192, B3 (DiffBind) = 216, C2 (K27ac lost-vs-gained) = 60; lost-side B2 = 3, B4 = 15, C1 = 15 (49_homer_k119ub_k27ac/homer_all_significant_motifs.tsv)
- Top H2AK119ub-gained (B1) motifs are neuronal bHLH: Atoh1 (31.5% vs 23.7%, FE 1.33, q≈0), Atoh7 (FE 1.49), NeuroD1 (FE 1.49), BHLHA15 (homer_top25_significant_motifs.tsv)
- Top H3K27ac lost (C2) also bHLH: Atoh1 (37.9% vs 33.8%, q≈0), Atoh7, BHLHA15 (homer_top25_significant_motifs.tsv)
- MOI directionality: NeuroD1 enriched B1/B3 (FE 1.49) but FE<1 in lost-side B2/B4; YY1 significant only in B3 (FE 1.85) (homer_motifs_of_interest_all.tsv)

**What this shows:** H2AK119ub sites that GAIN signal in BAP1-KO are enriched for cerebellar neuronal bHLH motifs (Atoh1, NeuroD1, Atoh7), and the same grammar marks LOST H3K27ac enhancers (C2). Because BAP1 normally removes H2AK119ub, its loss spreads Polycomb ubiquitination precisely over neuronal-differentiation TF binding sites. The signal is one-directional (lost-side B2/B4 nearly empty), and YY1 surfacing in the DiffBind K119ub contrast (B3) hints at a ubiquitination–architecture recruitment link.

**Figures:**
- homer_dotplot_top15 — top-15 motifs per key comparison (B1/B3/C1/C2)
- homer_family_barchart — TF-family counts per comparison
- homer_wordclouds — TF word clouds (B1/B3/C2)
- homer_motifs_of_interest_heatmap — curated MOI -log10(q) across B/C comparisons
