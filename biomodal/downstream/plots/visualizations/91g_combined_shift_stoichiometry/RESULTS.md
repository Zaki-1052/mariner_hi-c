## Section 91: Compartment Shift Stoichiometry

**Q1: Co-significance overlap by subcompartment (91_cosignificance_by_subcompartment.tsv):**
- A.1: 5,644 / 10,783 (52.3%) significant for BOTH 5mC and 5hmC
- A.2: 1,331 / 3,046 (43.7%)
- B.1: 598 / 1,897 (31.5%)
- B.2: 559 / 4,172 (13.4%)
- Monotonic decline from active to heterochromatic — co-occurrence of mC/hmC changes is compartment-dependent

**Q2: Compartment transition asymmetry (91_shift_methylation_summary.tsv):**
- A→B (becoming inactive, n=202): mean Δ5mC = +2.44%, Δ5hmC = −2.25%, 86.7% hypermethylated
- B→A (becoming active, n=427): mean Δ5mC = −1.04%, Δ5hmC = +1.29%, 93.4% hypomethylated
- A→B effect is 2.4× stronger in magnitude than B→A despite having half the genes
- The asymmetry suggests methylation gain may DRIVE compartment inactivation rather than the reverse

**Q3: Stoichiometry — redistribution, not removal (91_stoichiometric_ratios.tsv):**
- A→B: |Δ5mC|/|Δ5hmC| = 108.7% — near-stoichiometric; 91% of 5hmC loss recaptured as 5mC
- A→B median total modC change = −0.0018%, p = 0.987 (NOT significant) — total methylation perfectly conserved at the gene level
- Stable: ratio 88.1%, median total = −0.0027%, p < 2.2e-16 (slight net loss genome-wide)

**Q3b: Is 5hmC completely removed? No — it is partially depleted (91_residual_hmc_levels.tsv):**
- Median 5hmC retention across all genes: 91.5% of control levels
- 99.3% of genes retain >50% of their 5hmC; only 23 genes (0.1%) show near-complete depletion
- A.1 (most affected): ctrl 14.14% → mut 12.52% (88% retained, Δ = −1.57%)
- B.1 (facultative het): ctrl 7.33% → mut 7.40% (101% retained — actually slightly increases)
- Baseline 5hmC predicts depletion: Spearman rho = −0.589, p < 2.2e-16 (higher baseline → larger loss)

**TET impediment consistency:** The partial depletion pattern (not removal) is diagnostic. TET is the only enzyme that can actively remove 5hmC (via further oxidation to 5fC/5caC/C). If TET is blocked, 5hmC cannot be actively removed — it can only deplete through blocked production and normal turnover. The 91.5% median retention, the baseline-dependent depletion (rho = −0.59), and the near-zero total modC change at A→B genes (p = 0.99) are all consistent with impeded TET production, not active DNMT3A dehydroxymethylase conversion (which would drive 5hmC toward zero).

**Figures:**
- 91a — Co-significance stacked bar by subcompartment (52.3% → 13.4% gradient)
- 91b — Δ5mC and Δ5hmC by compartment shift (A→B vs B→A vs Stable)
- 91c — Stoichiometric ratio bar with 100% reference line
- 91d — Total modC balance violin (Δ5mC + Δ5hmC per gene, with medians/p/n)
- 91e — 5hmC depletion vs baseline scatter (rho = −0.589)
- 91f — Residual 5hmC levels ctrl vs mut by subcompartment (with medians/Δ/p/n)
