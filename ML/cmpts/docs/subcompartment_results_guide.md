# Subcompartment Analysis: What We Found So Far

Personal reference for understanding the CALDER2 results. Last updated: 2026-05-28.

> **IMPORTANT CAVEAT (added 2026-05-28):** The differential enrichment panel (Section 4) uses BigWigs that have been RPKM + 99th-percentile normalized across replicates. This normalization **masks global abundance changes** for ChIP marks. H2AK119ub and DNA methylation differential values in the heatmaps do NOT reflect their known genome-wide increases -- they show relative redistribution only. See Section 4 for full explanation.

---

## 1. Background: What Are We Looking For?

### The BAP1 Story

BAP1 is an enzyme called a deubiquitinase. Its job is to remove a specific chemical tag -- H2AK119ub (a ubiquitin mark on histone H2A at position 119) -- from chromatin. BAP1 does this as part of a complex called PR-DUB.

H2AK119ub is placed there by Polycomb Repressive Complex 1 (PRC1). Polycomb complexes are gene silencers: they mark regions of the genome for repression. So BAP1 essentially counteracts Polycomb silencing by stripping off PRC1's mark.

In our experiment, we knocked out BAP1 in mouse cerebellum. The question is: **what happens to the 3D organization of the genome when you remove this Polycomb antagonist?**

### Compartments and Subcompartments

When you do Hi-C (a technique that maps which parts of the genome physically touch each other in 3D), the genome falls into two broad categories:

- **A compartment** -- "active" chromatin. Gene-rich, transcribed, open, enriched for active histone marks (H3K27ac, H3K4me3).
- **B compartment** -- "inactive" chromatin. Gene-poor, repressed, compact, enriched for repressive marks (H3K27me3, H3K9me3).

This A/B split is what our earlier HOMER analysis detected. We found 8,189 bins (at 25kb resolution) that significantly changed between ctrl and BAP1-KO, with most shifting from B toward A.

But A and B are coarse. Each one contains distinct subtypes:

- **A.1** -- the most active subcompartment. Highest transcription, strongest active marks.
- **A.2** -- moderately active. Lower transcription than A.1 but still open.
- **B.1** -- the Polycomb/facultative heterochromatin compartment. Marked by H3K27me3. These are genes that are silenced but *could* be turned on (hence "facultative"). This is the compartment most directly relevant to BAP1's biology.
- **B.2** -- constitutive heterochromatin. Permanently silenced, gene-poor, late-replicating. Very low signal for everything.

(In human cells there's also a B.3 specific to certain chromosomes, but mouse doesn't have it.)

### Why Subcompartments Matter Here

At the Dixon lab meeting (2026-04-10), Jesse made a critical observation: the HOMER A/B changes we found might not be true compartment "flips." Instead, they could be **quantitative weakening within subcompartments** -- a B.1 region becoming a slightly weaker B.1 (PC1 value shifting toward A) without actually crossing the threshold into A. Or an A.2 region strengthening to A.1.

Subcompartment calling lets us test this directly: we can check whether bins that HOMER flags as "B to A" actually change their subcompartment label, or just shift within the same class.

### Why CALDER2

We used CALDER2 (an R package) because it works natively on mouse mm10 genome data -- no training data needed. The alternative, SNIPER, is a neural network trained only on human hg19 and would need retraining from scratch (which we're doing as secondary validation). CALDER2 classifies every genomic bin into the A.1/A.2/B.1/B.2 hierarchy using a reference-guided approach, correlating the Hi-C contact pattern with a known mm10 reference.

---

## 2. What We Ran

**Samples:** 4 merged Hi-C datasets from mouse cerebellum:

| Sample | Timepoint | Condition | Hi-C Quality |
|--------|-----------|-----------|-------------|
| 250402 ctrl_merged | Late (adult) | Wild-type control | KR normalization fallback on some chromosomes |
| 250402 mut_merged | Late (adult) | BAP1 knockout | Same |
| 250831 ctrl_merged | Early (postnatal day 12) | Wild-type control | Full KR normalization, no issues |
| 250831 mut_merged | Early (postnatal day 12) | BAP1 knockout | Same |

**Resolution:** CALDER2 ran at 50kb (automatically extends to 100kb for multi-resolution optimization). The final subcompartment labels are reported per 100kb bin for downstream analysis. 19 autosomes (chr1-19), ~23,800-23,900 callable bins per sample.

**Quality check -- reference correlation:** CALDER2 scores each chromosome's subcompartment ranking against a reference. All chromosomes passed (>0.59 correlation). The early timepoint has systematically better correlation (median 0.77 vs 0.69 at late), meaning the P12 compartment structure is more cleanly defined -- probably because adult cerebellum has diverged further from the reference cell type through developmental compaction.

**Runtime:** CALDER2 itself ran in ~2-3 minutes per sample. Differential analysis (A2) took 10 seconds. Epigenomic validation (A3) took 2.8 minutes. HOMER integration (A4) took 13 seconds. Total compute: under 15 minutes of SLURM time for the entire Track A pipeline.

---

## 3. Key Finding 1: BAP1-KO Opens Up the Genome

The single clearest result: **BAP1 knockout shifts the genome toward active (A) compartment identity, at the expense of constitutive heterochromatin (B.2).**

### The numbers

| Subcompartment | Late ctrl | Late mut | Late change | Early ctrl | Early mut | Early change |
|---------------|-----------|----------|------------|-----------|----------|-------------|
| A.1 (most active) | 30.6% | 32.5% | **+1.9 pp** | 32.0% | 35.9% | **+3.9 pp** |
| A.2 (moderate active) | 14.1% | 15.3% | +1.2 pp | 18.8% | 19.3% | +0.5 pp |
| B.1 (Polycomb) | 13.8% | 13.5% | -0.3 pp | 15.5% | 14.2% | -1.3 pp |
| B.2 (constitutive het.) | 41.6% | 38.7% | **-2.9 pp** | 33.7% | 30.5% | **-3.2 pp** |
| **Total A** | **44.7%** | **47.8%** | **+3.1 pp** | **50.8%** | **55.2%** | **+4.4 pp** |

> **Figures:** `outputs/calder2/late/250402_subcompartment_genome_pct/` and `outputs/calder2/early/250831_subcompartment_genome_pct/`

### Transition matrices

Looking at individual bins tells a sharper story. At the late timepoint, 15.3% of bins (3,659 / 23,853) changed their subcompartment label between ctrl and mut. At the early timepoint, 18.3% (4,371 / 23,823) changed.

The transitions are heavily directional -- B-to-A outnumbers A-to-B by 3.7:1 at late and 4.9:1 at early.

**Late (adult) -- top transitions:**

| Transition | Bins | % of changed | What it means |
|-----------|------|-------------|--------------|
| B.1 -> A.2 | 852 | 23.3% | Polycomb compartment opening to moderately active |
| B.2 -> B.1 | 887 | 24.2% | Deep heterochromatin becoming slightly less repressed |
| A.2 -> A.1 | 786 | 21.5% | Already-active regions becoming more active |
| A.1 -> A.2 | 384 | 10.5% | Some active regions weakening (the minority direction) |

**Early (P12) -- top transitions:**

| Transition | Bins | % of changed | What it means |
|-----------|------|-------------|--------------|
| A.2 -> A.1 | 1,258 | 28.8% | Massive within-A strengthening |
| B.1 -> A.2 | 1,250 | 28.6% | Polycomb opening (same as late but ~50% more bins) |
| B.2 -> B.1 | 966 | 22.1% | Heterochromatin loosening |

> **Figures (show these to PI):** `outputs/calder2/late/250402_transition_heatmap/` and `outputs/calder2/early/250831_transition_heatmap/`

### Interpretation

The genome-level story is a cascade of opening: B.2 feeds into B.1, B.1 feeds into A.2, and A.2 strengthens to A.1. BAP1 normally counteracts Polycomb silencing at specific loci. Without it, you might expect *more* silencing (since H2AK119ub should accumulate). Instead, we see the opposite -- the genome opens. This could reflect compensatory mechanisms or a more complex relationship between BAP1 and chromatin organization than simple H2AK119ub removal.

Note that B.1 itself barely changes in size (-0.3pp at late, -1.3pp at early). This is because it's simultaneously losing bins to A.2 (the B.1->A.2 transition) and gaining bins from B.2 (the B.2->B.1 transition). It's a throughput compartment, not a stable one.

---

## 4. Key Finding 2: Epigenomic Marks Confirm the Subcompartment Calls

We overlaid 9 epigenomic marks (ChIP-seq, ATAC-seq, RNA-seq, DNA methylation) onto the subcompartment labels to validate them and find differential signals.

### Important caveat: what the differential panel can and cannot show

The BigWig files used in this analysis were generated through a pipeline of RPKM normalization, then 99th-percentile cross-sample scaling, then replicate averaging. This normalization chain **compresses global abundance differences** for ChIP-seq marks. If a mark increases uniformly across the genome (as H2AK119ub does in BAP1-KO), the 99th-percentile step rescales the mutant signal back down toward the control level before we ever compute the differential.

The differential panel therefore shows **relative redistribution between subcompartments**, not absolute changes in abundance. A mark that goes up everywhere will appear flat or even negative in the heatmap because the normalization has already absorbed the global shift.

This caveat applies to all ChIP-seq marks (H3K27ac, H3K4me3, H3K36me3, H3K27me1, H3K27me3, H2AK119ub, ATAC). It does NOT apply to RNA-seq (count-based) or DNA methylation (fraction-based, 0-1 scale -- see below). The ctrl enrichment panel (which compares subcompartments within a single condition) is unaffected -- it correctly validates the calls.

**Known ground-truth from other analyses in this repo:**
- **H2AK119ub increases genome-wide in BAP1-KO.** DiffBind found 21,812 significant peaks, ALL with positive fold change (zero significant decreases). This is the expected biology: BAP1 is a deubiquitinase for H2AK119ub, so its loss causes H2AK119ub accumulation. (Source: `peaks/diffbind/K119ub_diffbind_results_summit_appended_ap.txt`)
- **DNA methylation (5mC) increases globally in BAP1-KO.** Biomodal analysis: 75.1% of significant genes are hypermethylated, mean +2.27%. (Source: `notes/BIOMODAL_RESULTS.md`)

### Validation: the calls are correct

The ctrl enrichment heatmaps show exactly the expected pattern:

| Mark | Expected pattern | Observed (late ctrl) |
|------|-----------------|---------------------|
| H3K27ac (active enhancer) | A.1 >> A.2 > B.1 > B.2 | 3.89 > 1.58 > 0.89 > 0.36 |
| H3K4me3 (active promoter) | A.1 >> A.2 > B.1 > B.2 | 14.75 > 2.74 > 0.86 > 0.34 |
| ATAC (open chromatin) | A.1 >> A.2 > B.1 > B.2 | 3.90 > 1.74 > 0.97 > 0.40 |
| RNA (transcription) | A.1 >> A.2 > B.1 > B.2 | 21.82 > 5.04 > 0.49 > 0.01 |
| H3K27me3 (Polycomb) | B.1 >> others | A.1=0.66, A.2=1.75, **B.1=2.25**, B.2=0.88 |

All gradients are monotonic and in the right direction. The subcompartment calls are biologically real.

### Differential: what the redistribution pattern shows

Given the normalization caveat above, the differential panel is best read as "where does this mark's signal concentrate or deplete *relative to the normalized genome-wide level*" -- not as absolute changes. Large signals here are informative about redistribution even though absolute levels are masked.

**H3K27me1 -- massive redistribution into B compartments.**

| Subcompartment | Late log2FC | Early log2FC |
|---------------|------------|-------------|
| A.1 | -0.24 | -0.18 |
| A.2 | +0.50 | +0.65 |
| B.1 | **+1.18** | **+1.17** |
| B.2 | **+1.35** | **+1.29** |

Even after normalization, H3K27me1 shows a 2-2.5 fold relative increase in B compartments of the mutant. This is the largest differential signal for any mark in the heatmap. H3K27me1 is placed by PRC2 (the same complex that trimethylates to H3K27me3). The B-compartment concentration could indicate PRC2 redistribution -- spreading into or concentrating in heterochromatic regions when BAP1 is lost. Because the normalization compresses global changes, the true absolute increase may be even larger.

**H3K27me3 -- relative enrichment in A compartments.**

| Subcompartment | Late log2FC | Early log2FC |
|---------------|------------|-------------|
| A.1 | **+0.40** | +0.35 |
| A.2 | +0.31 | +0.25 |
| B.1 | +0.13 | +0.22 |
| B.2 | -0.07 | -0.14 |

The redistribution pattern shows H3K27me3 concentrating in A compartments in the mutant. Even if the absolute level doesn't change much genome-wide, the relative shift toward A.1 (+0.40 log2FC) suggests ectopic Polycomb encroachment into active regions. This is biologically coherent: without BAP1 to antagonize Polycomb, H3K27me3 spreads into territory it normally wouldn't occupy.

**H2AK119ub -- redistribution pattern (NOT absolute decrease).**

| Subcompartment | Late log2FC | Early log2FC |
|---------------|------------|-------------|
| A.1 | -0.29 | -0.37 |
| A.2 | -0.09 | -0.13 |
| B.1 | -0.16 | -0.15 |
| B.2 | -0.32 | -0.34 |

The negative values here are a **normalization artifact**, not a real decrease. H2AK119ub is known to increase genome-wide from DiffBind (21,812 peaks, all up). The 99th-percentile normalization rescaled the mutant BigWig downward to match the control, making the signal appear decreased. What the redistribution pattern does tell us: relative to the (elevated) genome-wide level, B.2 bins have relatively less H2AK119ub enrichment than A bins. Combined with the cluster-level analysis (which found K119ub specifically elevated at lost-loop anchors in clust6: ctrl 1.16 -> mut 1.58), the picture is that H2AK119ub spreads broadly but concentrates at specific loci (lost enhancers, Polycomb targets) rather than uniformly.

**DNA methylation -- negative values need investigation.**

| Subcompartment | Late log2FC | Early log2FC |
|---------------|------------|-------------|
| A.1 | -0.24 | -0.24 |
| A.2 | -0.22 | -0.23 |
| B.1 | -0.26 | -0.26 |
| B.2 | -0.34 | -0.37 |

The Biomodal analysis shows global hypermethylation (+2.27% mean, 75% of genes up). DNA methylation BigWigs encode methylation fractions (0-1), not RPKM counts, so the 99th-percentile normalization issue should not apply in the same way. However, these BigWigs may have been generated from a different data type (CpG fraction scores from DUET evoC) and further investigation is needed to confirm whether the apparent decrease here reflects a genuine subcompartment-level pattern, a difference in CpG density between bins, or a BigWig construction artifact. The Biomodal result (gene-level, replicate-controlled) should be treated as the ground truth.

**Other marks (redistribution patterns):**
- H3K36me3 depletes from B compartments (-0.32 to -0.40) -- gene-body mark redistributing toward A
- H3K4me3 depletes from B.2 (-1.05 log2FC) but gains in A.2 (+0.36) -- active promoter mark concentrating in active regions
- ATAC shows very slight uniform increases (+0.06) -- consistent with mild global chromatin opening

> **Figure (show to PI):** `outputs/calder2/late/250402_enrichment_combined/` -- the three-panel heatmap with ctrl, mut, and differential side by side. The ctrl and mut panels are reliable for subcompartment validation. **The differential panel should be shown with the verbal caveat that it shows redistribution, not absolute changes, due to normalization.** Also `outputs/calder2/early/250831_enrichment_combined/` for comparison.

---

## 5. Key Finding 3: Most HOMER "Compartment Switches" Are Not True Flips

This directly tests Jesse's hypothesis from the Dixon meeting.

### Setup

Our earlier HOMER analysis (at 25kb resolution) found 8,189 significant bins where the A/B compartment PC1 value changed between ctrl and mut (FDR<0.05, |Difference|>0.30). We called these "compartment switches." But PC1 is a continuous value -- a bin can shift *within* a compartment (say, from weak-A to strong-A) and still look significant.

By overlaying CALDER2 subcompartment labels, we can classify each HOMER-significant bin:

- **True A-to-B flip:** CALDER2 ctrl = A.x, CALDER2 mut = B.y
- **True B-to-A flip:** CALDER2 ctrl = B.x, CALDER2 mut = A.y
- **Within-A shift:** stays in A but changes subtype (A.1 <-> A.2)
- **Within-B shift:** stays in B but changes subtype (B.1 <-> B.2)
- **Stable:** same CALDER2 label in both conditions despite significant HOMER PC1 change

### Results

**Late (adult) -- 7,575 significant autosomal HOMER bins:**

Of the 5,158 HOMER bins called "B to A":
| Category | Bins | % | Meaning |
|----------|------|---|---------|
| True B-to-A | 1,788 | 34.7% | Genuine compartment flip confirmed by CALDER2 |
| Stable | 1,683 | 32.6% | Same CALDER2 label -- quantitative weakening only |
| Within-B shift | 947 | 18.4% | B.2 -> B.1 (stays in B, just less repressed) |
| Within-A shift | 735 | 14.2% | Already in A, strengthening within A |

Of the 2,417 HOMER bins called "A to B":
| Category | Bins | % |
|----------|------|---|
| Stable | 1,416 | 58.6% |
| Within-A shift | 368 | 15.2% |
| True A-to-B | 339 | 14.0% |
| Within-B shift | 274 | 11.3% |

**Bottom line: only 28.4% of HOMER-significant bins are true compartment flips (2,152 / 7,575). The other ~71% are either stable or shifting within their compartment class.**

Jesse was right. The HOMER analysis overestimates the extent of compartment reorganization. The genome is weakening and strengthening within its existing subcompartment identity -- a subtler effect than wholesale A<->B switching.

That said, the 28% that are true flips are real and biologically meaningful. The dominant true flip is B.1 -> A.2 (Polycomb compartment opening), which is consistent with the transition matrix results.

> **Figures (show to PI):** `outputs/integration/late/250402_homer_calder2_sankey/` -- the 3-axis sankey diagram showing HOMER call -> CALDER2 ctrl label -> CALDER2 mut label. The gray "change" ribbons and the colored "stable" bands tell the story visually. Also `outputs/integration/late/250402_homer_calder2_dotplot/` as a cleaner summary.

---

## 6. Early vs Late: Developmental Context

The two timepoints tell complementary stories.

**The early (P12) genome is more open at baseline.** In ctrl, early has 50.8% A compartment vs 44.7% at late. The adult cerebellum accumulates more B.2 constitutive heterochromatin (41.6% vs 33.7%) -- this is normal developmental compaction as cell types mature and silence unused genes permanently.

**BAP1-KO has a larger effect at P12.** 18.3% of bins change subcompartment at early vs 15.3% at late. The B-to-A ratio is also higher (4.9:1 vs 3.7:1). The developing brain is more plastic -- chromatin states are still being established, so BAP1 loss has more room to alter subcompartment identity.

**But HOMER detects more significant changes at late.** 7,575 significant bins at late vs 5,036 at early. This seems contradictory but makes sense: HOMER measures continuous PC1 shift magnitude, not discrete label changes. The adult genome may have larger quantitative shifts in PC1 value (bigger "weakening" signal) even though fewer bins actually cross a subcompartment boundary. The adult compartment structure is more rigid, so when it does weaken, the PC1 delta is larger.

**Polycomb is weaker at P12.** H3K27me3 enrichment in B.1 is 1.41x at early vs 2.25x at late. The Polycomb compartment is still being established at postnatal day 12 -- the adult-level Polycomb concentration hasn't been reached yet. This is consistent with B.1 being a smaller fraction at early (15.5%) than might be expected.

---

## 7. What's Left: SNIPER + Integration (Track B & C)

### Track B: SNIPER (secondary validation)

SNIPER is a neural network that classifies subcompartments using inter-chromosomal contact patterns. It was originally trained on human hg19, so we're retraining it on mm10 using our CALDER2 labels as ground truth (B1-B5 pipeline steps).

The key deliverable: if SNIPER agrees with CALDER2 on the control samples (Cohen's kappa > 0.6), we can be confident the calls are robust across two completely different methods. On the mutant samples (which SNIPER hasn't seen during training), agreement would further validate the biological reality of the transitions.

**Status:** Environment set up (B0), mm10 adaptation code written (B1). Crop map generation and training (B2-B4) are pending HPC submission.

### Track C: Integration analyses

These are the "so what" analyses that connect subcompartments to the rest of the project:

- **C1 (genome-wide tracks):** Visual chromosome-by-chromosome view showing subcompartment changes along each chromosome. Good for identifying hotspot regions.

- **C2 (publication figure):** Polished version of the enrichment heatmap using ComplexHeatmap. The current ggplot versions are fine for working, but a pub-quality version with proper clustering and annotation is needed.

- **C3 (loops x subcompartments):** The key biological question still unanswered -- **do the gained Polycomb loops (cluster 5 from the loop clustering pipeline, which are 97% up-in-mutant) sit in B.1?** If yes, BAP1-KO is simultaneously opening B.1 at the compartment level while forming new PRC1-mediated loops within it. That's a nuanced and publishable finding.

- **C4 (HOMER decomposition):** The definitive version of the HOMER x CALDER2 analysis from Finding 3. Better figures, more complete breakdown.

---

## 8. Open Questions

### BigWig normalization limits the differential panel

The most important methodological limitation of the A3 analysis: the RPKM + 99th-percentile normalization pipeline used to generate ChIP BigWigs compresses global abundance differences. For marks that change globally in BAP1-KO (H2AK119ub increases, DNA methylation increases), the differential panel shows redistribution patterns, not absolute changes. The ctrl and mut enrichment panels are unaffected -- they correctly validate subcompartment identity within each condition.

**Possible future fix:** Re-run A3 with spike-in-normalized or unnormalized BigWigs if available, or apply a correction factor derived from the DiffBind global fold-change estimate. Alternatively, compute the differential at peak regions rather than genome-wide bins (matching the DiffBind approach).

### H3K27me1 -- what does the massive redistribution mean?

Even after normalization compresses global changes, H3K27me1 shows the largest redistribution signal: +1.2-1.35 log2FC in B compartments. This is placed by PRC2 (the same complex that places H3K27me3). The pattern -- relative increase in me1 at B compartments, relative increase in me3 at A compartments -- could indicate PRC2 redistribution: PRC2 spreading more broadly into heterochromatin (depositing initial me1) while also encroaching on active regions (depositing me3 at new targets). Whether this redistribution reflects a real mechanistic change in PRC2 targeting or is secondary to H2AK119ub accumulation (which recruits PRC2 via PRC1) is an open question.

### B.1 as a throughput compartment

B.1 barely changes in total size, but this masks huge turnover: ~850-1250 bins leave B.1 for A.2, while a similar number arrive from B.2. Is this balanced flow a coincidence, or does it reflect a biological mechanism where BAP1 loss creates a "conveyor belt" of chromatin states moving from constitutive toward active?

### HOMER early timepoint polarity

A polarity fix was applied to the HOMER eigenvector for the early timepoint (the PC1 sign was inverted relative to CALDER2 orientation). The A4 script now includes a gene-density polarity check. The re-run showed correct directionality, but the early sankey/dotplot should be interpreted knowing this correction was applied.

### DNA methylation BigWig source

The `DNAmethylationCtrl.bw` / `DNAmethylationMut.bw` files encode CpG methylation fractions (0-1, from Biomodal DUET evoC), not RPKM counts. In principle this should preserve global methylation changes. But the A3 differential shows global decreases (-0.22 to -0.37), contradicting the Biomodal gene-level result (+2.27% hypermethylation). This discrepancy needs investigation -- it may reflect differences in how the BigWigs were constructed (CpG density weighting, genomic context, or a different normalization applied during BigWig generation) versus the gene-level methylKit analysis.

---

## 9. Figure-by-Figure Explainer

How to read each figure, what the axes mean, and the key result to point out.

### Transition Heatmap (both timepoints)

`250402_transition_heatmap/` and `250831_transition_heatmap/`

**Format:** 4x4 grid (confusion matrix style). Rows = ctrl subcompartment label, columns = mut subcompartment label. Each cell contains the number of 100kb genomic bins with that ctrl-to-mut transition. Color intensity = log10(count+1), so darker blue = more bins. Red borders highlight the diagonal (bins that kept the same label). Subtitle shows the chi-squared statistic and the percentage of bins that changed.

**How to read it:** The diagonal is stability -- bins that stayed the same. Everything off-diagonal is a transition. Look at the *asymmetry*: in the late heatmap, B.1->A.2 has 852 bins but A.2->B.1 has only 233. That 3.7:1 ratio is the directional bias of BAP1-KO.

**Key results to point out:**
- Diagonal dominates (84-82% of bins are stable -- overall architecture is preserved)
- Below the diagonal is mostly empty (A->B transitions are rare)
- Above the diagonal has substantial counts (B->A transitions are common)
- B.1->A.2 is the single largest off-diagonal cell at both timepoints -- the Polycomb compartment opening up
- A.2->A.1 is the second-largest -- already-active regions strengthening
- Early has 18.3% changed vs late 15.3% -- more plasticity at P12

### Genome Fraction Bar Chart (both timepoints)

`250402_subcompartment_genome_pct/` and `250831_subcompartment_genome_pct/`

**Format:** Two stacked bar charts side by side. Left bar = control, right bar = mutant. Each segment is a subcompartment (A.1 red, A.2 orange, B.1 green, B.2 blue), sized by its percentage of the genome. Y-axis = % of 100kb bins.

**How to read it:** Compare the heights of each color segment between ctrl and mut. The red (A.1) section grows, the blue (B.2) section shrinks. This is the simplest visual of the A-compartment expansion.

**Key results to point out:**
- Late: A.1 grows from 30.6% to 32.5% (+1.9pp), B.2 shrinks from 41.6% to 38.7% (-2.9pp)
- Early: A.1 grows from 32.0% to 35.9% (+3.9pp), B.2 shrinks from 33.7% to 30.5% (-3.2pp)
- B.1 (green) barely changes -- it's a throughput compartment (losing bins to A.2 but gaining from B.2)
- Early ctrl starts with more A (51%) than late ctrl (45%) -- developmental difference, P12 brain is more open

### Enrichment Combined Heatmap (both timepoints)

`250402_enrichment_combined/` and `250831_enrichment_combined/`

**Format:** Three panels side by side, each a 9-row x 4-column grid. Rows = epigenomic marks (H3K27ac, H3K4me3, ATAC, H3K36me3, RNA, DNAmethylation, H3K27me1, H2AK119ub, H3K27me3). Columns = subcompartments (A.1, A.2, B.1, B.2). Color = red for high values, blue for low, white = neutral.

- **Left panel (Ctrl Enrichment):** fold-enrichment of each mark over the genome-wide median, using ctrl labels and ctrl signal. This validates the subcompartment calls.
- **Middle panel (Mut Enrichment):** same but using mut labels and mut signal.
- **Right panel (Differential):** log2(mut_median / ctrl_median) per mark per subcompartment, using ctrl labels. Shows redistribution of signal between conditions. **Caveat: BigWig normalization masks global changes for ChIP marks -- these values show relative redistribution, not absolute direction. See Section 4.**

**How to read it:** In the ctrl panel, look for the expected gradient: H3K27ac/H3K4me3/ATAC/RNA should be dark red in A.1 and white/blue in B.2 (active marks concentrated in active compartment). H3K27me3 should be red in B.1 (Polycomb compartment). If these gradients are correct, the subcompartment calls are validated. In the differential panel, look for the largest color blocks -- H3K27me1 in B.1/B.2 is bright red (+1.2-1.35), the biggest redistribution signal.

**Key results to point out:**
- Ctrl panel: all expected gradients are correct and monotonic -- calls are real
- H3K27me3 lights up in B.1 (2.25x at late, 1.41x at early) -- Polycomb compartment confirmed
- RNA is 21.8x enriched in A.1 vs 0.01 in B.2 -- extreme dynamic range validates A/B identity
- Differential: H3K27me1 has the strongest redistribution signal (B.2 = +1.35 log2FC at late)
- Differential: H3K27me3 shifts toward A.1 (+0.40 log2FC) -- Polycomb encroaching on active regions
- **When presenting:** note that the differential panel shows WHERE marks concentrate, not WHETHER they go up or down globally. H2AK119ub goes UP genome-wide (DiffBind) but appears negative here due to normalization

### HOMER -> CALDER2 Sankey (late timepoint)

`250402_homer_calder2_sankey/`

**Format:** Three-axis alluvial/sankey flow diagram. Left axis = HOMER call direction (A->B or B->A). Middle axis = CALDER2 ctrl subcompartment label. Right axis = CALDER2 mut subcompartment label. Colored ribbons connect them, showing the flow. Ribbon width = number of 25kb bins. Colors: red=A.1, orange=A.2, green=B.1, blue=B.2, gray=change (bins that switched subcompartment). Y-axis = number of 25kb bins. Title shows total significant HOMER bins (7,575 at late).

**How to read it:** Follow the ribbons from left to right. The bottom group (HOMER: B->A) is the majority. From there, ribbons flow through B.1 and B.2 in the middle (their ctrl label) and then either stay the same color on the right (stable) or turn gray and jump to a different subcompartment (true flip). The amount of gray vs colored ribbon is the key: gray = true flips, colored = stable despite HOMER flagging them.

**Key results to point out:**
- This directly tests Jesse's hypothesis from the Dixon meeting
- The large mass of colored (non-gray) ribbons staying on the diagonal = bins that HOMER calls significant but CALDER2 says didn't actually change subcompartment
- Only ~28% of significant HOMER bins are true A<->B flips (the gray ribbons that cross from B to A)
- ~33% are "stable" (same CALDER2 label in ctrl and mut -- HOMER sees a quantitative PC1 shift but it doesn't cross a subcompartment boundary)
- ~32% are within-compartment shifts (B.2->B.1 or A.2->A.1 -- moving within the same broad compartment)
- **Punchline: most of what HOMER calls a "compartment switch" is actually quantitative weakening, not a true flip. Jesse was right.**

### HOMER -> CALDER2 Dotplot (late timepoint)

`250402_homer_calder2_dotplot/`

**Format:** Two panels side by side. Left = HOMER A->B bins, Right = HOMER B->A bins. X-axis = CALDER2 ctrl label, Y-axis = CALDER2 mut label. Each dot = one transition type; dot size = number of bins; dot color = transition category (green=True_A_to_B, pink=True_B_to_A, orange=Within_A_shift, blue=Within_B_shift, gray=Stable). Dashed diagonal = no change.

**How to read it:** Dots on the diagonal are stable (same ctrl and mut label). Dots above the diagonal moved "up" (toward more active). Dots below moved "down." The largest dot in the B->A panel (right) is the pink B.1->A.2 bubble -- the dominant true compartment flip. The large gray dots on the diagonal are the "stable" bins that HOMER flagged but CALDER2 says didn't change.

**Key results to point out:**
- Right panel (B->A): pink B.1->A.2 is the largest single transition -- Polycomb opening
- Right panel: large gray dots on the diagonal show many bins are HOMER-significant but CALDER2-stable
- Left panel (A->B): dominated by gray diagonal dots -- most HOMER A->B calls are not confirmed by CALDER2
- The asymmetry between panels (right has more and larger off-diagonal dots than left) shows the directional bias

### Subcompartment Flow Sankey (both timepoints)

`250402_transition_sankey/` and `250831_transition_sankey/`

**Format:** Two-axis alluvial diagram. Left axis = ctrl subcompartment distribution (sized by bin count). Right axis = mut distribution. Colored ribbons show which bins flow where. Gray ribbons = bins that changed subcompartment. Y-axis = number of 100kb bins.

**How to read it:** Wide colored bands on the diagonal = stability. Thin gray threads crossing between levels = transitions. The more gray, the more change. The B.2 block (bottom, blue) is visibly smaller on the right (mut) than the left (ctrl) -- the genome is opening.

**Key results to point out:**
- Mostly stable (thick colored diagonal bands)
- Visible gray threads from B.1 crossing up to A.2 level
- B.2 block shrinks, A.1 block grows
- These are visually busy -- the heatmap and bar chart are cleaner for presenting. Use this as supplementary if someone asks "but what does the flow look like?"

---

## Quick Reference: Figure Paths

### Tier 1 -- Show First

| Figure | Path (from ML/cmpts/) |
|--------|----------------------|
| Late transition heatmap | `outputs/calder2/late/250402_transition_heatmap/` |
| Late enrichment combined | `outputs/calder2/late/250402_enrichment_combined/` |
| Late HOMER x CALDER2 sankey | `outputs/integration/late/250402_homer_calder2_sankey/` |

### Tier 2 -- Show If Time

| Figure | Path (from ML/cmpts/) |
|--------|----------------------|
| Late genome fractions | `outputs/calder2/late/250402_subcompartment_genome_pct/` |
| Early genome fractions | `outputs/calder2/early/250831_subcompartment_genome_pct/` |
| Late HOMER x CALDER2 dotplot | `outputs/integration/late/250402_homer_calder2_dotplot/` |
| Early transition heatmap | `outputs/calder2/early/250831_transition_heatmap/` |
| Early enrichment combined | `outputs/calder2/early/250831_enrichment_combined/` |

### QC / Reference Only

| Figure | Path (from ML/cmpts/) |
|--------|----------------------|
| Late ctrl ref correlation | `outputs/calder2/late/250402/ctrl_merged/sub_compartments/cor_with_ref.pdf` |
| Late mut ref correlation | `outputs/calder2/late/250402/mut_merged/sub_compartments/cor_with_ref.pdf` |
| Early ctrl ref correlation | `outputs/calder2/early/250831/ctrl_merged/sub_compartments/cor_with_ref.pdf` |
| Early mut ref correlation | `outputs/calder2/early/250831/mut_merged/sub_compartments/cor_with_ref.pdf` |
| Late subcompartment flow (sankey) | `outputs/calder2/late/250402_transition_sankey/` |
| Early subcompartment flow (sankey) | `outputs/calder2/early/250831_transition_sankey/` |

All figure directories contain 4 formats: pdf, svg, png, jpg. Use PDF for PI meetings, PNG for Slack.
