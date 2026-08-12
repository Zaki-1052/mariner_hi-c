# Developmental Comparison: Analysis Reference

**Question:** Does adult BAP1-KO 3D genome resemble P13 wildtype, suggesting blocked developmental chromatin remodeling?

**Background:** Between P13 and adulthood, wildtype cerebellum undergoes extensive 3D genome remodeling -- loop acquisition, compartment strengthening, TAD reorganization. If BAP1-KO blocks that remodeling, the adult KO should look more like P13 WT than adult WT. Jesse's puzzle: BAP1 is more expressed early but the phenotype is stronger late. Blocked remodeling would explain this paradox -- the damage isn't from acute loss of BAP1 activity, it's from failing to remodel during a critical developmental window.

**Key distinction:** This is NOT "early KO effects vs late KO effects" (which `tads/scripts/timepoint_comparison.R` already does). This is a cross-condition, cross-timepoint similarity test: adult-KO vs P13-WT, compared to adult-WT vs P13-WT.

---

## Results Summary

| Layer | Key Metric | Supports Blocking? |
|-------|-----------|-------------------|
| Compartments (PC1) | r(late_KO, early_WT) = 0.899 vs r(late_WT, early_WT) = 0.882; Fisher Z=17.3, p=1.7e-67 | **Yes** |
| Loop positions (gained) | 84% of gained adult-KO loops overlap P13-WT vs 65% null; z=27.7, p<1e-4 | **Yes** |
| Loop positions (Jaccard) | J(late_KO, early_WT) = 0.479 vs J(late_WT, early_WT) = 0.456 | Modest yes |
| Loop distance | Late_KO median 225kb closer to late_WT (240kb) than early_WT (300kb) | **No** |
| TAD boundaries | Fisher OR=1.68, p=5.7e-7 but concordance only 45% (below 50%) | Mixed |
| Stripes | 35.6% of late-gained stripes overlap early-WT anchors | Supplementary |

### Maturation Index (% of WT developmental change completed by KO)

- **Compartments: 86%** -- partially blocked
- **Loop positions: 96%** -- nearly complete
- **Loop distance: 125%** -- overshot (KO loops even shorter than WT)

---

## Module-by-Module Explanations

### Module 1: PC1 Compartment Correlation (`01_pc1_correlation/`)

PC1 values are the first principal component of the Hi-C contact matrix at 25kb resolution. Positive = A compartment (active/euchromatin), negative = B compartment (inactive/heterochromatin). Each of ~95K genomic bins gets a PC1 value per condition. Pearson correlation between two conditions tells you how similar their compartment profiles are genome-wide.

**The three key correlations:**

- **r(late_KO, late_WT) = 0.985** -- Adult KO vs adult WT are very similar. BAP1 loss doesn't demolish compartments, it nudges them.
- **r(late_WT, early_WT) = 0.882** -- Adult WT vs P13 WT. This is the *normal developmental distance* -- compartments change substantially between P13 and adulthood as the cerebellum matures.
- **r(late_KO, early_WT) = 0.899** -- Adult KO vs P13 WT. This is higher than 0.882, meaning the adult KO compartment profile is *closer* to the immature P13-WT than the adult WT is.

**Fisher Z-test:** Is the difference between 0.899 and 0.882 statistically significant? With 95K bins, even a small correlation difference is powered -- z=17.3, p=1.7e-67.

**Biological meaning:** Normal development moves compartments from the P13 state (r=1.0 with early_WT) to the adult state (r=0.882). BAP1-KO doesn't make that full journey -- it stops at r=0.899, retaining more of the immature compartment signature. The effect is modest (~0.017 in correlation units) but genome-wide and highly significant, consistent with partial blocking of developmental compartment remodeling.

---

### Module 2: Loop Jaccard Overlap Matrix (`02_loop_jaccard/`)

Four loop sets were constructed from all-loops BEDPEs using the `direction` column:

- `early_WT` = unchanged + down_in_mutant (38,929 loops)
- `early_KO` = unchanged + up_in_mutant (38,958 loops)
- `late_WT` = unchanged + down_in_mutant (35,091 loops)
- `late_KO` = unchanged + up_in_mutant (35,616 loops)

Jaccard similarity (with 10kb anchor tolerance) measures what fraction of loops are shared between two sets.

**Key numbers:**

- Within-timepoint Jaccard is very high (early_WT vs early_KO = 0.993) because at P13, BAP1-KO barely changes the loop landscape (only 361 differential loops).
- Cross-timepoint Jaccard drops to ~0.46-0.48, meaning roughly half the loop repertoire turns over during normal development.
- J(late_KO, early_WT) = 0.479 vs J(late_WT, early_WT) = 0.456 -- the KO retains slightly more P13-WT loops, but the effect is modest because the Jaccard is dominated by the ~31K unchanged loops that are shared across everything.

---

### Module 3: Gained Loop Overlap with P13-WT (`03_loop_overlap_gained/`)

This is the most direct test. It asks: are the loops gained in adult BAP1-KO found at genomic positions that already had loops in P13-WT? If so, the adult KO is "retaining" an immature loop architecture rather than creating novel contacts.

**The permutation null distribution plot:**

- **Grey histogram (null):** 5,000 times, we randomly sampled 4,253 loops from the late-unchanged background and asked what fraction overlap P13-WT loop positions. By chance, ~65% of any random late-timepoint loops will overlap P13-WT simply because both timepoints share a large common loop repertoire. This is the baseline expectation.
- **Solid line at right (~84%, "Gained"):** The actual fraction of adult-KO gained loops (up_in_mutant) that overlap P13-WT loop positions. At 83.9%, this is 28 standard deviations above the null -- not even close to random. These gained loops are overwhelmingly found at positions that already existed in the immature P13-WT genome.
- **Dashed line (~72.5%, "Lost"):** The fraction of adult-KO lost loops (down_in_mutant) that overlap P13-WT. Also above null, but less extreme.

**Biological interpretation:** The loops that BAP1-KO gains in adulthood aren't new -- they're loops that were present in the immature wildtype cerebellum. The adult KO is preferentially strengthening contacts at P13-WT positions, consistent with a failure to complete developmental loop remodeling.

---

### Module 4: Loop Distance eCDF Overlay (`04_distance_ecdf/`)

Empirical CDFs of loop distance for all four groups on a single plot.

**Medians:** early_WT 300kb, early_KO 300kb, late_WT 240kb, late_KO 225kb.

**Interpretation:** During normal development, loops get shorter (300kb -> 240kb). BAP1-KO loops get *even shorter* (225kb), overshooting the WT developmental trend. The KS test confirms late_KO is more similar to late_WT (D=0.018) than to early_WT (D=0.081) in terms of distance distribution.

This does NOT support blocking at the loop distance level. Instead it suggests BAP1-KO accelerates or exaggerates the developmental shift toward shorter loops -- likely driven by Polycomb-mediated gains of short-range compartmental contacts combined with loss of long-range active loops.

---

### Module 5: TAD Boundary Direction Concordance (`05_tad_boundary_concordance/`)

If BAP1-KO blocks developmental remodeling, then boundaries that become KO-enriched in adults should correspond to boundaries that were WT-enriched (normal) at P13 -- i.e., the KO retained the immature boundary architecture.

**Result:** 1,520 boundaries are differential at both timepoints (within 20kb). Of those that are late-KO-enriched, 45% were early-ctrl-enriched. This is below the 50% threshold the simple blocking model predicts, though the Fisher OR of 1.68 (p=5.7e-7) shows the association is real -- it's just not a clean majority.

**Interpretation:** The TAD boundary data shows a significant but incomplete blocking pattern. Many boundaries that are KO-enriched in adults were *also* KO-enriched at P13 (55%), meaning the KO effect is partly constitutive rather than purely a failure to remodel.

---

### Module 6: Stripe Footprint Comparison (`06_stripe_overlap/`)

35.6% of late-gained stripes (1,078 / 3,032) overlap early-WT stripe anchors within 50kb tolerance. This is framed as supplementary given the weak early signal (only 96 significant stripes at P13 vs 2,320 at adult).

---

## Module 7: Intuitive Summary Panels (`07_summary_panels/`)

These are the PI-friendly versions of the above analyses.

### 7A: PCA of 12 Replicates (`pca_12_replicates/`)

Same format as an RNA-seq MDS plot. Each dot is one biological replicate (3 per group). PC1 (91% of variance) separates by developmental timepoint -- P13 on the left, adult on the right. PC2 (5%) separates by genotype.

The grey arrow shows the WT developmental trajectory (P13-WT centroid to Adult-WT centroid). Adult KO triangles sit along the same PC1 axis as Adult WT circles but are shifted upward on PC2, meaning BAP1-KO changes the compartment landscape without fully blocking the P13-to-adult transition. The subtitle reports KO completed 88% of the WT trajectory (projection of KO centroid onto the WT developmental vector).

### 7B: PC1 Scatter vs P13 Baseline (`pc1_scatter_vs_p13/`)

Each of ~95K genomic bins plotted: x = its PC1 at P13-WT, y = its PC1 at adult. The dashed white line is the diagonal (no change from P13). Points off the diagonal are bins that remodeled during development.

Compare the two panels: Adult KO (left, R^2=0.807, RMSD=0.450) hugs the diagonal slightly tighter than Adult WT (right, R^2=0.778, RMSD=0.485). That means the KO genome stayed closer to the P13 state -- less developmental remodeling. The difference is subtle but genome-wide across 95K bins.

### 7C: Gained Loop Origin (`gained_loop_origin/`)

The simplest plot. Of the 4,253 loops that BAP1-KO gains in adulthood, 84% (3,570) already existed in the P13 wildtype genome. Only 16% (683) are truly novel contacts. The dashed line shows that by chance you'd expect 65% overlap.

So gained KO loops are overwhelmingly "retained immature loops" rather than aberrant new contacts. This is the single most compelling piece of evidence for blocked remodeling.

### 7D: Developmental Maturation Index (`maturation_index/`)

The one-number-per-layer summary. The dashed line at 100% = WT level of developmental change. Bars below 100% = KO didn't fully mature (blocked). Bars above 100% = KO overshot.

- **Compartments: 86%** -- partially blocked from maturing
- **Loop positions: 96%** -- nearly complete maturation
- **Loop distance: 125%** -- overshot; KO shifted to even shorter loops than WT

This tells you the blocking is layer-specific: compartment identity is partially retained in an immature state, but cohesin-driven loop shortening (likely from Polycomb spreading impeding extrusion) goes beyond normal development.

---

## Overall Interpretation

The evidence for "blocked developmental remodeling" is **layer-specific**:

1. **Compartments** show the clearest blocking -- the adult KO retains more of the P13 compartment signature than WT does, genome-wide.
2. **Loop positions** show strong blocking in the gained loops specifically -- 84% of what the KO gains were already present at P13 -- but the overall loop Jaccard shows only modest blocking because the background of ~31K unchanged loops dominates.
3. **Loop distance** shows the *opposite* of blocking -- the KO overshoots the developmental shortening trend, consistent with Polycomb-mediated gain of short-range contacts.
4. **TAD boundaries** show a significant but incomplete pattern -- some boundaries are blocked, but many KO boundary effects are constitutive across timepoints.

The narrative framing: BAP1-KO partially blocks developmental compartment maturation while simultaneously driving Polycomb-mediated loop shortening that goes *beyond* normal development. The KO genome is not simply "frozen at P13" -- it's a hybrid state with immature compartment identity but exaggerated short-range Polycomb contacts. This is consistent with BAP1's role as an H2AK119ub deubiquitinase: loss of BAP1 increases K119ub, which promotes Polycomb spreading, which both blocks compartment remodeling *and* drives aberrant short-range chromatin contacts.

---

## File Index

```
output/developmental_comparison/
  01_pc1_correlation/
    pc1_correlation_replicate/   12x12 replicate-level heatmap
    pc1_correlation_summary/     4x4 group-level heatmap (panel-ready)
    pc1_fisher_z_test.txt
  02_loop_jaccard/
    loop_jaccard_matrix/         4x4 Jaccard heatmap
    loop_jaccard_stats.tsv
  03_loop_overlap_gained/
    gained_vs_early_wt_overlap/  Bar: gained vs lost overlap fractions
    permutation_null_distribution/  Null histogram with observed lines
    permutation_stats.tsv
  04_distance_ecdf/
    distance_ecdf_4groups/       4-condition eCDF overlay
    distance_ks_tests.tsv
  05_tad_boundary_concordance/
    boundary_direction_heatmap/  2x2 contingency heatmap
    boundary_concordance_bar/    Stacked bar of direction concordance
    boundary_contingency_table.tsv
    fisher_exact_result.txt
  06_stripe_overlap/
    stripe_anchor_overlap/       Stripe overlap bar
  07_summary_panels/
    pca_12_replicates/           MDS-style PCA (PI-friendly)
    pc1_scatter_vs_p13/          Side-by-side density scatter
    gained_loop_origin/          Single stacked bar (84%/16%)
    maturation_index/            % maturation per structural layer
  tables/
    developmental_comparison_stats.tsv
  developmental_comparison_report.txt
  analysis.md                    This file
```

**Script:** `scripts/developmental_comparison.R`
**Generated:** 2026-04-21
