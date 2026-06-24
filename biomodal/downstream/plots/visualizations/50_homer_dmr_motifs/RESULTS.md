## Section 50: HOMER Motif Enrichment (DMRs — A5 coordinated vs discordant)
**Key numbers:**
- A5 significant motifs (q < 0.05) = 128 (homer_a5_significant_motifs.tsv)
- Family breakdown: zinc-finger (Zf) = 34 (dominant), bHLH = 17, NR = 16, Homeobox = 12, ETS = 12 (homer_a5_significant_motifs.tsv)
- Top motifs (all Zf/ETS): ZEB2 (p≈1e-25, 98.1% vs 95.4%, FE 1.03), Maz (p≈1e-24), Klf15, Sp5, ETV4, KLF17/14/5/6, ZNF148, ZNF711 (homer_a5_significant_motifs.tsv)

**What this shows:** Regions with the coordinated mC-up/hmC-down signature (TET-block) are distinguished from discordant DMRs by an excess of GC-rich zinc-finger / KLF / Sp-family motifs (ZEB2, Maz, KLF15/14/5/6, Sp5) — the CpG-rich binding context where 5mC↔5hmC turnover by TET enzymes is most consequential and where methylation can directly modulate TF binding. Fold-enrichments are modest (≈1.02–1.07) because both target and background are CpG-rich; significance comes from large n. A1–A4 (genome background) produced no significant motifs and are excluded.

**Figures:**
- homer_a5_coordinated_dotplot — top-25 A5 motifs, size = %target, color = TF family
- homer_a5_family_counts — TF-family counts among significant A5 motifs
