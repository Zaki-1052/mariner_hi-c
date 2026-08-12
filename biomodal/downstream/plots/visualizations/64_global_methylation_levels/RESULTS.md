## Section 64: Global CpG Methylation Levels (BAP1-KO vs Wildtype)

**Key numbers:**
- 5mC: ctrl 72.216% -> mut 72.534%, delta +0.317% (paired t-test p=0.0412, Cohen's d=2.79) (64_statistics.tsv)
- 5hmC: ctrl 10.192% -> mut 9.797%, delta -0.395% (paired t-test p=0.0184, Cohen's d=-4.54) (64_statistics.tsv)
- modC (total): delta -0.079%, p=0.492 — NOT significant (64_statistics.tsv)
- 8 run-5 samples confirmed: 4 control + 4 mutant, 2 sexes x 2 batches (64_global_methylation_summary.tsv)

**What this shows:** Genome-wide autosomal CpG methylation shifts by a small but extremely reproducible amount in BAP1-KO: 5mC rises and 5hmC falls by near-equal magnitude while total modified-C stays flat. This is the TET-pathway-blockade signature (5mC->5hmC oxidation blocked), not global hypermethylation. Every individual mutant sample shows the same direction. Paired Wilcoxon is floored at p=0.125 by n=4 pairs, so the t-test is the operative significance test.

**Figures:**
- 64_global_methylation_levels/ — combined 3-panel figure
- 64a_global_methylation_dotbar/ — per-sample dot+bar by modality, paired t-test
- 64b_replicate_delta/ — matched-replicate (mut-ctrl) deltas per sex/batch
- 64c_composition_stacked/ — 5mC+5hmC+unmethylated composition, ctrl vs mut
