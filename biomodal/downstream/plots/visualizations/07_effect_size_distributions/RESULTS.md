## Section 07: Effect Size Distributions

**Key numbers:**
- 5mC net mean change (significant) = +1.72% (n=10,775); hyper-only = +3.45% (n=7,513, 70%) (source: run-5 gene-body mC BED)
- 5hmC net mean change (significant) = -1.66% (n=11,484); hypo-only = -2.29% (n=9,521, 83%) (source: run-5 gene-body hmC BED)
- Both means match `summary_statistics.txt`; 5mC density shifts positive, 5hmC negative (mirror-image violins)

**What this shows:** The opposing net shifts (+1.7% mC vs -1.7% hmC) are the population-level expression of the coordinated pattern. Per-gene effects are modest (~2-3%) but consistent across thousands of genes, making the directional signal statistically overwhelming and ruling out a noise explanation.

**Figures:**
- `07_effect_size_distributions/` — 5mC effect histogram, 5hmC effect histogram, 5mC-vs-5hmC violin comparison.

> Note: `FIGURES.md`/`TODO.md` quote +2.27%/-2.08% (hyper-only +3.96%, n=6635) — OLD-run; current run-5 = +1.72%/-1.66% (+3.45%, n=7,513).
