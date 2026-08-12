# SNIPER + Integration Results: What We Found (Track B and C)

Personal reference for understanding the SNIPER and integration analysis results. Last updated: 2026-06-02.

> **This is the sequel to `docs/subcompartment_results_guide.md`, which covers Track A (CALDER2) results: compartment biology background, what we ran, and Findings 1-3. This document covers Track B (SNIPER neural network) and Track C (integration analyses). It does not repeat the biology or CALDER2 results -- read the CALDER2 guide first.**

---

## 1. What This Document Covers

Track B is a secondary validation of our subcompartment calls using SNIPER, a deep learning model that predicts subcompartments from inter-chromosomal contact patterns -- a completely different input signal than CALDER2's intra-chromosomal eigenvector approach. SNIPER was originally built for human hg19, so we retrained it from scratch on mm10 using our CALDER2 labels as ground truth. Track C is the integration layer that connects subcompartments to the rest of the project: loop clusters, stripes, HOMER A/B compartments, and epigenomic marks. All Track B and C stages completed 2026-06-01/02.

---

## 2. Track B: Retraining SNIPER for mm10

### 2.1 Why SNIPER Needed Retraining

SNIPER's pre-trained models are hardcoded for human hg19 with five hard blockers: (1) chromosome ranges locked to chr1-22, (2) neural network input dimensions baked into the .h5 weights (13,393/13,594 bins), (3) hg19-specific crop maps, (4) hg19-specific coordinate maps for BED output, and (5) a 5-class output layer (A1/A2/B1/B2/B3) that doesn't apply to mouse (no B.3 in mm10).

We retrained anyway because SNIPER provides orthogonal validation: it uses inter-chromosomal contacts (how chromosomes interact with each other in 3D space) while CALDER2 uses intra-chromosomal contact eigenvectors (how a chromosome is organized within itself). If two fundamentally different input signals agree on the subcompartment calls, we can be confident the biology is real.

### 2.2 What the mm10 Adaptation Changed

| Component | hg19 (original) | mm10 (adapted) |
|-----------|-----------------|----------------|
| Autosomes | 22 (11 odd, 11 even) | 19 (10 odd, 9 even) |
| Inter-chr matrix | ~13,393 x 13,594 bins | 12,808 x 11,831 bins |
| Subcompartment classes | 5 (A1/A2/B1/B2/B3) | 4 (A1/A2/B1/B2) |
| Training labels | Human Rao et al. GM12878 | CALDER2 calls on ctrl_merged |
| Training protocol | Cross-cell-type (high-cov GM12878 target) | Self-supervised (input = target, same .hic) |
| DAE epochs | 10 (upstream bug) | 25 (paper-specified) |

The architecture itself is unchanged: a denoising autoencoder (DAE) learns a 128-dimensional latent representation of inter-chromosomal contact patterns, then a small MLP classifier maps those latent features to 4 subcompartment labels. The DAE is trained self-supervised (the ctrl_merged.hic is both input and target), meaning the autoencoder learns to denoise its own contact matrix. The latent bottleneck forces it to compress the contact pattern into features that capture compartment-level structure.

### 2.3 Pipeline Steps

| Stage | What | Result |
|-------|------|--------|
| B0 | Environment setup (Python 3.6, TF-GPU 1.12, CUDA 9.0 via conda) | Verified on Tesla V100 |
| B1 | mm10 crop map generation from ctrl_merged.hic | rowMap ~(10,686, 3), colMap ~(9,462-9,463, 3). ~80-83% bin retention after blacklist + sparse filtering |
| B2 | Convert CALDER2 labels to SNIPER format (.mat) | 250402: rows=10,686 cols=9,463. 250831: rows=10,686 cols=9,462. 50-70 uncallable bins imputed per TP |
| B3 | Train SNIPER on ctrl_merged (1 per TP, GPU) | DAE converged by epoch ~20 (val_loss=0.503 late). All 12 model files (6 per TP) saved |
| B4 | Apply SNIPER to all 4 samples | 20,148-20,149 predictions per sample. Label distributions consistent with CALDER2 |
| B5 | Concordance analysis (SNIPER vs CALDER2) | Kappa 0.64-0.69. See Section 3 |

---

## 3. Track B Results: How Well Does SNIPER Agree with CALDER2?

### 3.1 Training Fidelity and Generalization

| Timepoint | Condition | n_bins | Accuracy | Kappa | F1(A.1) | F1(A.2) | F1(B.1) | F1(B.2) |
|-----------|-----------|--------|----------|-------|---------|---------|---------|---------|
| 250402 (late) | ctrl | 20,076 | 78.0% | 0.695 | 0.858 | 0.570 | 0.580 | 0.899 |
| 250402 (late) | mut | 20,073 | 74.8% | 0.656 | 0.815 | 0.532 | 0.572 | 0.893 |
| 250831 (early) | ctrl | 20,078 | 74.8% | 0.660 | 0.878 | 0.646 | 0.545 | 0.810 |
| 250831 (early) | mut | 20,078 | 73.7% | 0.643 | 0.855 | 0.625 | 0.534 | 0.816 |

Cohen's kappa of 0.64-0.69 is "substantial agreement" -- well above the 0.40 threshold for moderate agreement. SNIPER recovers the extreme compartments well: A.1 and B.2 have F1 scores above 0.81 in all conditions. These are the compartments with the most distinctive inter-chromosomal contact signatures (A.1 = active hubs contacting other active chromosomes; B.2 = constitutive heterochromatin with sparse inter-chromosomal contacts).

The ~3% accuracy drop from ctrl to mut (78.0% -> 74.8% late, 74.8% -> 73.7% early) is the key generalization metric. SNIPER was trained only on ctrl_merged, so the mutant samples are unseen data with BAP1-KO perturbation. A 3% drop means the learned contact-pattern features transfer well -- the perturbation doesn't break the fundamental inter-chromosomal signature of subcompartments.

> **Figures:** `outputs/sniper/plots/250402_confusion_heatmap_ctrl/` and `250402_per_class_concordance/`

### 3.2 Where SNIPER and CALDER2 Disagree

The intermediate classes -- A.2 and B.1 -- are where the tools diverge. F1 scores of 0.53-0.65 mean SNIPER is essentially coin-flipping on these bins. The confusion pattern is systematic: SNIPER over-predicts A.2 and B.1 (high recall, low precision ~0.50) at the expense of under-predicting A.1 and B.2 (high precision >0.94, lower recall). In plain terms, SNIPER fragments the extreme compartments into intermediates -- it calls some genuine A.1 bins as A.2 and some genuine B.2 bins as B.1.

This makes biological sense. A.2 and B.1 sit at the A/B compartment boundary where the inter-chromosomal contact signal is weakest. These bins interact with both active and inactive chromosomes, so the contact-based features that SNIPER learns are inherently ambiguous there. CALDER2, which uses intra-chromosomal eigenvectors (a different signal), can resolve them better.

The practical consequence: when SNIPER confirms a transition, it's high-confidence. When SNIPER disagrees with CALDER2, it's usually on a borderline A.2/B.1 bin where both tools are uncertain.

### 3.3 Transition Concordance

**250402 (late) -- 20,063 bins with all 4 labels (SNIPER ctrl+mut, CALDER2 ctrl+mut):**

| Category | Bins | % | Meaning |
|----------|------|---|---------|
| Both_stable | 15,123 | 75.4% | Same label in both tools, both conditions |
| Both_change_agree | 707 | 3.5% | Both detect same transition |
| Both_change_disagree | 339 | 1.7% | Both detect change but disagree on direction |
| SNIPER_only | 1,754 | 8.7% | SNIPER sees change, CALDER2 says stable |
| CALDER2_only | 2,140 | 10.7% | CALDER2 sees change, SNIPER says stable |

CALDER2 detects ~22% more transitions than SNIPER in the shared bin set (3,186 vs 2,800 changed bins at late). SNIPER is more conservative -- it misses transitions at the A.2/B.1 boundary where its resolution is lowest. The 707 concordant bins are mined in C5 (Finding 7) as the highest-confidence subset of subcompartment changes.

**250831 (early):** Similar pattern. Both_stable=72.2%, concordant=695 bins, CALDER2_only=13.7%. The higher CALDER2_only fraction at early (13.7% vs 10.7%) reflects the fact that early has more CALDER2 transitions (18.3%) while SNIPER detects a similar count (14.1%).

---

## 4. Key Finding 4: Developmental Reorganization Dominates Genotype Effect

### The Numbers

| Comparison | What's being compared | Changed bins | % of genome |
|------------|----------------------|-------------|------------|
| ctrl: early vs late | Normal developmental maturation (P12 to adult) | 7,790 / 23,796 | **32.7%** |
| mut: early vs late | BAP1-KO + developmental maturation | 8,436 / 23,782 | **35.5%** |
| late: ctrl vs mut | Genotype effect in adult cerebellum | 3,659 / 23,853 | 15.3% |
| early: ctrl vs mut | Genotype effect at P12 | 4,371 / 23,823 | 18.3% |

The genome reorganizes roughly twice as many subcompartments during normal P12-to-adult cerebellum maturation (32.7%) as it does from BAP1 knockout at a single timepoint (15.3-18.3%). This is the developmental axis: as the cerebellum matures, neurons permanently silence developmental genes and compact their genome. That process moves ~1 in 3 bins between subcompartments.

### Interpretation

The mutant is significantly more developmentally dynamic than the control (35.5% vs 32.7%, p=3.34e-38 homogeneity test). BAP1 loss doesn't just open the genome at a fixed timepoint -- it destabilizes subcompartment maintenance over developmental time. The mutant genome is less locked into its identity as it matures.

This is a different effect from "BAP1-KO causes compartment opening" (Finding 1 in the CALDER2 guide). Finding 1 says the genome opens. Finding 4 says the opened genome is also more fluid during development. Together: BAP1 loss both shifts the genome toward active compartment and makes that shifted state less stable.

> **Figures:** `outputs/calder2/combined/combined_stability_barplot/` and `outputs/calder2/combined/combined_developmental_heatmap/`

---

## 5. Key Finding 5: Polycomb Loops Form Inside B Compartment Territory

This is the most biologically direct connection between the subcompartment analysis and the rest of the project. It answers the question Jesse raised: where in the compartment hierarchy do the gained Polycomb loops sit?

### Background

The loop clustering pipeline (in `cluster/`) identified 6 loop clusters from 38,948 differential loops. Cluster 5 (clust5, 1,296 anchors at late) is the gained-Polycomb-loop cluster: 97% of its loops are up in mutant, and its anchors are enriched for H3K27me3 and H2AK119ub. Cluster 6 (clust6, 4,658 anchors) contains the pre-existing active loops anchored by CTCF/cohesin. The question: do clust5 loops form in B.1 (facultative heterochromatin, the Polycomb compartment)?

### The Answer

**Clust5 anchor subcompartment composition (250402, ctrl):**

| Subcompartment | % of clust5 anchors | Genome background |
|---------------|--------------------|--------------------|
| A.1 | 13.2% | 30.6% |
| A.2 | 16.8% | 14.1% |
| B.1 | **27.5%** | 13.8% |
| B.2 | **42.4%** | 41.6% |

70% of clust5 anchors fall in B compartment (B.1 + B.2). The B.1 enrichment is striking: OR=3.22 (95% CI [2.83, 3.65], p=2.96e-63). For comparison, clust6 (pre-existing active loops) has the opposite profile: 49.1% A.1, 22.9% A.2, and is depleted from B.1 (OR=0.64, p=1.19e-49).

At the early timepoint, the B.1 enrichment is present but weaker: OR=2.39 (p=7.67e-33). This is consistent with the Polycomb compartment being less established at P12 -- H3K27me3 enrichment in B.1 is 1.41x at early vs 2.25x at late (from Finding 2, CALDER2 guide).

### The Paradox: Loops Forming in a Compartment That's Opening

At clust5 anchors, 37.2% change subcompartment between ctrl and mut (vs ~15% genome-wide), of which 16.8% are B-to-A transitions. BAP1 loss is simultaneously:

1. **Opening the B.1 compartment** -- bins shifting from B.1 to A.2 (the dominant genome-wide transition)
2. **Forming new PRC1-mediated loops within B.1** -- clust5 loops are gained in the mutant

This is the clust5 paradox. One interpretation: a subset of B.1 is becoming more internally organized (gaining Polycomb loops, possibly as a compensatory response to H2AK119ub accumulation) while another subset is derepressing. Another: the loops form transiently in bins that are in the process of opening -- Polycomb attempts to hold them in B.1 before they escape to A.2.

### Stripes Follow the Same Pattern

Gained stripes (from the Stripenn pipeline) show modest B.1 enrichment at late (OR=1.28, p=3.76e-4) but only a trend at early (OR=1.21, p=0.06). Stripe body composition analysis: gained stripes have 15.1% B.1 + 15.5% B.2 in their body, compared to 10.3% B.1 + 7.4% B.2 for lost stripes. Gained 3D structures (both loops and stripes) preferentially occupy repressed compartment territory.

> **Figures:** `outputs/integration/loops_stripes/plots/c3_loop_subcompartment_enrichment/` (balloon plot showing clust5 B.1 enrichment) and `outputs/integration/loops_stripes/plots/c3_clust5_clust6_transitions/` (alluvial showing transitions at clust5/6 anchors)

---

## 6. Key Finding 6: Only 2.8% of the Genome Undergoes True A<->B Flips

This is the definitive version of Finding 3 from the CALDER2 guide, now run at native 100kb resolution (C4) with effect-size validation, the genome-wide "iceberg" framing, and SNIPER overlay.

### The Iceberg

Most of the genome is invisible to compartment analysis. Of the bins that are HOMER-significant, most are not true compartment flips. The framing: it's an iceberg, and the visible tip (true flips) is tiny.

**250402 (late/adult) -- 23,840 callable bins:**

| Category | Bins | % of genome | Meaning |
|----------|------|------------|---------|
| Non-significant | 21,092 | 88.5% | No detectable HOMER PC1 change |
| Significant + Stable | 1,246 | 5.2% | HOMER flags a PC1 shift, but CALDER2 label unchanged |
| Significant + Within-shift | 828 | 3.5% | Subcompartment moves within A or within B (e.g. A.2->A.1, B.2->B.1) |
| Significant + True flip | 674 | **2.8%** | Genuine A-to-B or B-to-A compartment flip |

**250831 (early/P12):** 92.2% non-significant, 3.8% stable, 2.0% within-shift, **1.9% true flip**. The early timepoint has fewer and smaller true flips -- the developing genome is more plastic but its subcompartment structure is also less firmly established, so "flips" are harder to call.

Note: the true flip fraction is lower at 100kb (24.5% of significant bins) than at 25kb in the A4 analysis (28.4%). Aggregating from 25kb to 100kb smooths borderline bins, making the decomposition cleaner.

### Effect Sizes Validate the Classification

This is important: the classification is not circular. HOMER significance is based on continuous PC1 change (from the Hi-C eigenvector). CALDER2 labels come from a separate contact-pattern analysis. Yet the three categories show a clean monotonic gradient in |mean delta PC1|:

**250402 (late):**

| Category | Median |delta PC1| | Wilcoxon vs True flip |
|----------|----------------------|----------------------|
| True flip | 0.414 | -- |
| Within-shift | 0.352 | p = 4.01e-15 |
| Stable | 0.311 | p = 3.33e-53 |

**250831 (early):** Same gradient (0.376 > 0.327 > 0.303), all pairwise p < 10^-11.

True compartment flips involve larger PC1 magnitude changes than within-compartment shifts, which in turn are larger than "stable" bins. CALDER2's categorical labels track the magnitude of the continuous PC1 signal, even though they were derived independently.

### Cross-Timepoint Concordance Is Low

Only 439 bins are HOMER-significant at both timepoints. Of those, only 33.9% agree on their transition category (true flip at both, or within-shift at both). This likely reflects developmental progression rather than noise: a bin showing "within-B shift" at P12 may have progressed to a "true flip" by adulthood as the B compartment crystallized.

> **Figures:** `outputs/integration/homer_decomposition/plots/c4_iceberg_stacked/` (the iceberg visual) and `outputs/integration/homer_decomposition/plots/c4_publication_panel/` (2x2 combined panel: iceberg + violin + decomposition + chromosome heatmap)

---

## 7. Key Finding 7: Concordant Transitions Reveal the Highest-Confidence Changes

C5 mines the bins where both SNIPER and CALDER2 agree on the exact ctrl-to-mut transition. These are the 707 bins (late) where two completely different methods see the same thing.

### What "Concordant" Means

Of the 3,659 bins that CALDER2 calls as changed at the late timepoint, SNIPER confirms 707 of them (19.3% confirmation rate). At early, it's 695 / 4,371 = 15.9%.

The 80% that SNIPER doesn't confirm are not necessarily wrong. They're concentrated in the A.2/B.1 boundary zone where SNIPER has F1 scores of 0.53-0.65. SNIPER can't reliably distinguish A.2 from B.1, so it can't confirm transitions between them. The 707 confirmed bins are those where even the boundary-confused SNIPER agrees -- making them the highest-confidence subset.

### What the Concordant Transitions Show

**Observation 1: All concordant transitions are nearest-neighbor.**

Zero skip-transitions were observed in the concordant set. No bins jumped from A.1 directly to B.2, or from B.2 to A.1. Every confirmed transition moves exactly one step in the subcompartment hierarchy: A.1<->A.2, A.2<->B.1, or B.1<->B.2. This validates the ordered hierarchy -- subcompartment changes follow the ranking, not random jumps.

**Observation 2: B.2->B.1 is the most robustly detected transition.**

| Transition | Concordant bins | Confirmation rate | Meaning |
|-----------|----------------|-------------------|---------|
| B.2 -> B.1 | 261 | **29.4%** | Constitutive het. becomes Polycomb het. |
| B.1 -> A.2 | 180 | 21.1% | Polycomb compartment opening |
| A.2 -> A.1 | 81 | 10.3% | Moderate active strengthening |
| A.1 -> A.2 | 71 | 18.5% | Active weakening |
| B.1 -> B.2 | 52 | 18.3% | Polycomb becoming constitutive |
| A.2 -> B.1 | 50 | 21.5% | Moderate active closing |

B.2->B.1 has the highest confirmation rate (29.4%) and the most concordant bins (261). This is constitutive heterochromatin becoming more Polycomb-like -- the most reliably detected event across both tools. It's the first step in the B.2->B.1->A.2 cascade described in Finding 1 (CALDER2 guide).

**Observation 3: SNIPER's enrichment gradients are steeper than CALDER2's.**

For matched bins, SNIPER produces purer extreme-compartment calls:

| Mark | SNIPER gradient (late ctrl) | CALDER2 gradient (late ctrl) |
|------|---------------------------|----------------------------|
| H3K27ac | 4.25 > 1.55 > 0.80 > 0.31 | 3.89 > 1.58 > 0.89 > 0.36 |
| H3K4me3 | 15.3 > 2.67 > 0.75 > 0.28 | 14.75 > 2.74 > 0.86 > 0.34 |

SNIPER's A.1 enrichments are higher and its B.2 enrichments are lower because SNIPER classifies fewer borderline bins into the extreme classes (it puts them in A.2/B.1 instead). The bins it does assign to A.1 and B.2 are the unambiguous ones with the clearest epigenomic profiles.

> **Figures:** `outputs/integration/sniper_concordant/plots/c5_concordant_heatmap_late/` (definitive transition heatmap) and `outputs/integration/sniper_concordant/plots/c5_three_way_comparison_late/` (CALDER2 | SNIPER | Concordant side-by-side)

---

## 8. Publication Figure: ComplexHeatmap Enrichment (C2)

C2 is the polished version of the A3 working heatmap, designed for publication. It adds several features over the ggplot2-based A3 output:

- **ComplexHeatmap** (Bioconductor) for proper clustering and annotation
- Row annotations with mark functional category (Active / Gene body / Methylation / Repressive)
- Column annotations with subcompartment color bar and n_bins barplot
- Row split by mark category and column split by A/B compartment class
- **Wilcoxon significance stars** overlaid on differential cells (paired rank-sum ctrl vs mut per bin)

Key numbers: 32/36 mark-subcompartment pairs are significant per timepoint (Wilcoxon, BH-adjusted p < 0.05). The non-significant pairs are DNA methylation (relatively uniform across subcompartments) and H3K27me1 in some subcompartments.

The main output is a 4-panel heatmap (2 timepoints x 2 conditions) with a shared color scale, and a 2-panel differential with asterisks marking significant changes. The normalization caveat from Section 4 of the CALDER2 guide still applies to the differential panels.

> **Figure (show to PI):** `outputs/calder2/enrichment/c2_enrichment_panel/`
> **Differential with stars:** `outputs/calder2/enrichment/c2_differential_panel/`

---

## 9. Open Questions

### A.2/B.1 Boundary: The Most Uncertain Subcompartment

Both SNIPER (F1=0.53-0.65 for A.2 and B.1) and the modest concordance at the A/B boundary suggest these intermediate subcompartments are the hardest to distinguish. Are they genuinely overlapping in their contact patterns, or does the model/method lack resolution? Answering this would require Hi-C at finer resolution or additional marks (H4K20me1, LMNB1) that better separate B.1 from B.2.

### What Genes Are in the 707 Concordant Bins?

We know where the highest-confidence transitions are. The next step is gene content analysis: do the 261 B.2->B.1 concordant bins (the largest group) correspond to known Polycomb targets? Do they overlap with the DEG list from the RNA-seq analysis? This would connect the 3D genome reorganization directly to gene expression changes.

### The Clust5 Paradox: Forming Loops in a Compartment That's Opening

B.1 is simultaneously losing bins to A.2 (compartment opening) and gaining Polycomb loop structures (clust5 enriched with OR=3.22). Does this mean a subset of B.1 is becoming more internally organized while another subset is derepressing? Or are these the same bins, where Polycomb loops form transiently before derepression completes? Resolving this requires tracking which specific B.1 bins gain clust5 loops versus which ones transition to A.2 -- a bin-level intersection of the loop and compartment data.

### Developmental vs Genotype: The Relative Magnitudes

32.7% developmental change vs 15.3% genotype change. The developmental reorganization is not the biological question of this project, but it sets the effect-size context. For a reviewer: the genotype effect is real and significant, but it rides on top of a larger developmental program. Should the developmental transition matrix (C1 outputs) appear as supplemental material to make this explicit?

### SNIPER's Conservative Transition Detection

SNIPER detects 13.9% changed bins vs CALDER2's 15.3% (late). But the per-transition confirmation rates are uneven: B.2->B.1 is confirmed at 29% while A.2->A.1 is confirmed at only 10%. Is SNIPER specifically better at detecting B-compartment transitions (where inter-chromosomal contacts are sparse and distinctive) and worse at A-compartment transitions (where inter-chromosomal contacts are denser and more ambiguous)?

---

## 10. Figure-by-Figure Explainer

### 10.1 SNIPER Confusion Heatmaps (B5)

`outputs/sniper/plots/250402_confusion_heatmap_ctrl/` and `*_mut/`

**Format:** 4x4 matrix. Rows = CALDER2 label (ground truth). Columns = SNIPER prediction. Color intensity = number of bins. Diagonal = agreement.

**How to read it:** The diagonal should dominate (SNIPER agrees with CALDER2). Look at the A.2 and B.1 rows/columns -- the off-diagonal mass concentrates there. Specifically: SNIPER calls some CALDER2 A.1 bins as A.2 (row A.1, column A.2), and some CALDER2 B.2 bins as B.1 (row B.2, column B.1). The extreme compartments (A.1 and B.2 diagonal) are dark, confirming strong agreement there.

**Key result to point out:** The confusion is systematic, not random -- it's always SNIPER fragmenting extremes into intermediates, never the reverse. This means SNIPER's confirmed bins are higher-confidence than CALDER2's borderline calls.

### 10.2 SNIPER Per-Class Concordance Bar (B5)

`outputs/sniper/plots/250402_per_class_concordance/`

**Format:** Grouped bar chart. X-axis = subcompartment class. Y-axis = F1 score. Bars grouped by condition (ctrl blue, mut orange).

**How to read it:** A.1 and B.2 bars are tall (~0.81-0.90). A.2 and B.1 bars are shorter (~0.53-0.65). The ctrl bars are slightly taller than mut bars (training data vs unseen data).

**Key result to point out:** The consistent pattern across all 4 conditions -- extremes agree, intermediates diverge -- is the simplest summary of SNIPER-CALDER2 concordance.

### 10.3 C1 Combined Transition Heatmap

`outputs/calder2/combined/combined_transition_heatmap/`

**Format:** Side-by-side 4x4 transition heatmaps (early left, late right). Color = log10(count+1). Rows = ctrl label, columns = mut label.

**How to read it:** Compare the off-diagonal intensities between panels. The early heatmap has more off-diagonal mass (18.3% changed vs 15.3%), especially A.2->A.1 (1,258 bins at early vs 786 at late). This is the developmental plasticity effect: P12 genome is more malleable.

### 10.4 C1 Developmental Stability Barplot

`outputs/calder2/combined/combined_stability_barplot/`

**Format:** Four grouped bars. Groups = genotype effect (ctrl, mut) and developmental effect (ctrl, mut). Bar height = fraction of bins changed.

**How to read it:** The developmental bars (~33-35%) tower over the genotype bars (~15-18%). This single figure makes the Finding 4 point visually: developmental reorganization is ~2x the genotype effect. The mutant developmental bar is slightly taller than the ctrl developmental bar (35.5% vs 32.7%) -- BAP1-KO adds instability.

### 10.5 C2 Publication Enrichment Panel

`outputs/calder2/enrichment/c2_enrichment_panel/`

**Format:** 4-panel heatmap (2 timepoints x 2 conditions). Rows = 9 epigenomic marks, grouped by category (Active, Gene body, Methylation, Repressive). Columns = 4 subcompartments (A.1, A.2, B.1, B.2) with A|B split. Row annotations on left, column annotations on top (color bar + bin count barplot). Color = fold-enrichment over genome median.

**How to read it:** Each cell = how enriched this mark is in this subcompartment relative to the genome-wide median. Dark red = highly enriched. Active marks (H3K27ac, H3K4me3, ATAC) should be dark red in A.1 and fade to blue in B.2. H3K27me3 should light up in B.1 only.

**Key result to point out:** The consistency across all 4 panels (both conditions, both timepoints) validates the subcompartment calls. Even in the mutant, the epigenomic gradients hold. The subtlety is in the column annotations -- the A.1 bin count is larger in mutant panels than ctrl panels, reflecting compartment opening.

### 10.6 C2 Differential Panel with Significance Stars

`outputs/calder2/enrichment/c2_differential_panel/`

**Format:** 2-panel (early, late). Same layout as enrichment but color = log2(mut/ctrl). Asterisks mark Wilcoxon-significant cells (BH-adjusted p < 0.05).

**How to read it:** Red cells = mark increased in mutant relative to ctrl at that subcompartment. Blue = decreased. Asterisks confirm statistical significance. H3K27me1 in B.1/B.2 should be the most intensely red cells with asterisks. **Show with the normalization caveat from CALDER2 guide Section 4: these are redistribution patterns, not absolute changes for ChIP marks.**

### 10.7 C3 Loop Enrichment Balloon Plot

`outputs/integration/loops_stripes/plots/c3_loop_subcompartment_enrichment/`

**Format:** Balloon plot. X-axis = subcompartment. Y-axis = loop cluster (clust1-6). Dot size = number of anchors. Dot color = fold-enrichment (observed/expected). Red = enriched, blue = depleted. Faceted by timepoint.

**How to read it:** The clust5 row is the outlier: deeply blue for A.1 (depleted) and red for B.1/B.2 (enriched). Clust1-3 rows are deeply red for A.1 -- most loops live in active compartment. Clust4 is intermediate. This visual makes the clust5 exception immediately obvious.

**Key result to point out:** Clust5 is the only cluster enriched in B.1. Point at the B.1 column, clust5 row -- that red dot (OR=3.22) is the Finding 5 result.

### 10.8 C3 Clust5/Clust6 Transition Alluvial

`outputs/integration/loops_stripes/plots/c3_clust5_clust6_transitions/`

**Format:** Alluvial/sankey diagram showing ctrl->mut subcompartment transitions at loop anchors, split by cluster. Colored ribbons = stable, gray = changed.

**How to read it:** For clust5, gray ribbons (changed) are much thicker than for clust6 -- 37.2% vs 19% of anchors change. The gray ribbons from B.1 crossing up to A.2 show the compartment-opening-at-loop-anchors effect.

### 10.9 C4 Iceberg Stacked Bar

`outputs/integration/homer_decomposition/plots/c4_iceberg_stacked/`

**Format:** Two stacked bars (early, late). Segments from bottom: Non-significant (gray, largest), Sig+Stable (yellow), Sig+Within-shift (orange), Sig+True-flip (red, smallest). Y-axis = fraction of genome.

**How to read it:** The vast non-significant gray dominates (88-92%). The small red tip at the top is the true flips (2-3%). This is the iceberg: most of the genome is not significantly changing, and of the significant bins, only a fraction are genuine A<->B flips.

### 10.10 C4 Publication Panel

`outputs/integration/homer_decomposition/plots/c4_publication_panel/`

**Format:** 2x2 patchwork. Top-left: iceberg stacked bar (the overview). Top-right: PC1 violin by transition category (validates that true flips have larger |delta PC1|). Bottom-left: decomposition comparison (A->B vs B->A breakdown, early vs late). Bottom-right: per-chromosome true-flip heatmap.

**Key result to point out:** The violin plot (top-right) is the most important subplot -- it shows the monotonic gradient (true flip > within-shift > stable) that validates the CALDER2 classification independently against HOMER's continuous PC1 metric.

### 10.11 C5 Three-Way Comparison Heatmap

`outputs/integration/sniper_concordant/plots/c5_three_way_comparison_late/`

**Format:** Three 4x4 heatmaps side by side. Left = CALDER2 transitions (all changed bins), center = SNIPER transitions (same bins), right = Concordant (bins both tools agree on). Same color scale across all three.

**How to read it:** The concordant matrix (right) is cleaner than either individual tool's matrix. Off-diagonal mass concentrates in the nearest-neighbor cells only (A.1-A.2, A.2-B.1, B.1-B.2). No skip-transitions. The concordant matrix is the definitive, high-confidence transition map.

### 10.12 C5 Confirmation Rates Bar

`outputs/integration/sniper_concordant/plots/c5_confirmation_rates_late/`

**Format:** Bar chart. X-axis = specific transition type (B.2->B.1, B.1->A.2, etc.). Y-axis = confirmation rate (% of CALDER2 transitions confirmed by SNIPER).

**How to read it:** B.2->B.1 has the tallest bar (29.4%). A.2->A.1 has the shortest meaningful bar (10.3%). Skip-transitions (A.1->B.2, B.2->A.1) have 0% confirmation -- they don't appear in the concordant set at all.

---

## Quick Reference: Figure Paths (Track B and C)

### Tier 1 -- Show First

| Figure | Path (from ML/cmpts/) |
|--------|----------------------|
| C2 publication enrichment panel | `outputs/calder2/enrichment/c2_enrichment_panel/` |
| C3 loop enrichment balloon | `outputs/integration/loops_stripes/plots/c3_loop_subcompartment_enrichment/` |
| C4 publication panel (2x2) | `outputs/integration/homer_decomposition/plots/c4_publication_panel/` |
| C5 three-way comparison (late) | `outputs/integration/sniper_concordant/plots/c5_three_way_comparison_late/` |
| C1 stability barplot | `outputs/calder2/combined/combined_stability_barplot/` |

### Tier 2 -- Show If Time

| Figure | Path (from ML/cmpts/) |
|--------|----------------------|
| B5 confusion heatmap (late ctrl) | `outputs/sniper/plots/250402_confusion_heatmap_ctrl/` |
| B5 per-class F1 bars | `outputs/sniper/plots/250402_per_class_concordance/` |
| C1 combined transition heatmap | `outputs/calder2/combined/combined_transition_heatmap/` |
| C3 loop subcompartment stacked bar | `outputs/integration/loops_stripes/plots/c3_loop_subcompartment_stacked/` |
| C4 iceberg stacked bar | `outputs/integration/homer_decomposition/plots/c4_iceberg_stacked/` |
| C4 effect size violin | `outputs/integration/homer_decomposition/plots/c4_effect_size_violin/` |
| C5 confirmation rates (late) | `outputs/integration/sniper_concordant/plots/c5_confirmation_rates_late/` |
| C5 concordant heatmap (late) | `outputs/integration/sniper_concordant/plots/c5_concordant_heatmap_late/` |

### Supplementary / QC

| Figure | Path (from ML/cmpts/) |
|--------|----------------------|
| C2 differential (Wilcoxon stars) | `outputs/calder2/enrichment/c2_differential_panel/` |
| C3 stripe enrichment | `outputs/integration/loops_stripes/plots/c3_stripe_subcompartment_enrichment/` |
| C3 developmental comparison | `outputs/integration/loops_stripes/plots/c3_developmental_comparison/` |
| C4 chromosome heatmap | `outputs/integration/homer_decomposition/plots/c4_chromosome_heatmap/` |
| C4 SNIPER validation bar | `outputs/integration/homer_decomposition/plots/c4_sniper_validation/` |
| C4 cross-TP alluvial | `outputs/integration/homer_decomposition/plots/c4_cross_tp_alluvial/` |
| C5 concordant sankey (late) | `outputs/integration/sniper_concordant/plots/c5_concordant_sankey_late/` |
| C5 concordant vs discordant | `outputs/integration/sniper_concordant/plots/c5_concordant_vs_discordant_late/` |
| C5 SNIPER enrichment heatmaps | `outputs/integration/sniper_concordant/plots/c5_sniper_enrichment_ctrl_late/` |
| B5 discordant alluvial | `outputs/sniper/plots/250402_discordant_alluvial/` |
| B5 transition agreement | `outputs/sniper/plots/250402_transition_agreement/` |
| C1 developmental heatmap | `outputs/calder2/combined/combined_developmental_heatmap/` |
| C1 karyotype track | `outputs/calder2/combined/combined_karyotype_track/` |

All figure directories contain 4 formats: pdf, svg, png, jpg. Use PDF for PI meetings, PNG for Slack.
