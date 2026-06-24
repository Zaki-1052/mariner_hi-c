## Section 03: DMR Statistics by Genomic Region

**Key numbers:**
- Gene bodies significant mC DMRs = 10,775 / 20,969 (51.4%) — dominant region (source: run-5 gene-body mC BED / `summary_statistics.txt`)
- Other regions (raw, non-deduped): CpG shores 9,842; shelves 6,924; promoters 1,692; CpG islands 442; TSS 192 (source: run-5 region BEDs)
- 5mC direction: 7,513 hyper (69.7%) / 3,262 hypo; 5hmC direction: 1,963 up / 9,521 down (82.9%) (source: run-5 BEDs)
- Non-CG baseline: CHG 0.628%, CHH 0.862% (vs CG 72.22%), 0 significant DMRs (source: upstream CSV)

**What this shows:** Differential methylation is overwhelmingly a gene-body CpG phenomenon (51% of gene bodies vs 1.4% of TSS), and the directionality is reciprocally asymmetric — gene bodies gain 5mC while losing 5hmC. Non-CG methylation is negligible and unchanged.

**Figures:**
- `03a_dmr_by_region/` — significant mC DMRs across 6 region classes (gene bodies lead).
- `03b_direction_comparison/` — % hyper/hypo (5mC) vs % up/down (5hmC).
- `03_dmr_region_statistics/` — combined 4-panel region + context summary.
