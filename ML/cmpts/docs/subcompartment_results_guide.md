# Subcompartment Analysis: What We Found So Far

Personal reference for understanding the CALDER2 results. Last updated: 2026-05-28.

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

## 4. Key Finding 2: Epigenomic Marks Confirm the Story

We overlaid 9 epigenomic marks (ChIP-seq, ATAC-seq, RNA-seq, DNA methylation) onto the subcompartment labels to validate them and find differential signals.

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

### Differential: what changes in BAP1-KO

Three marks stand out in the differential panel (log2 fold-change, mut vs ctrl):

**H3K27me1 -- the biggest mover.**

| Subcompartment | Late log2FC | Early log2FC |
|---------------|------------|-------------|
| A.1 | -0.24 | -0.18 |
| A.2 | +0.50 | +0.65 |
| B.1 | **+1.18** | **+1.17** |
| B.2 | **+1.35** | **+1.29** |

H3K27me1 (monomethylation of H3K27) increases by 2-2.5 fold in B compartments of the mutant. This mark is associated with PRC2 activity (PRC2 places mono-, di-, and trimethylation on H3K27). The massive B-compartment increase could indicate PRC2 redistribution -- PRC2 spreading into or concentrating in heterochromatic regions when BAP1 is lost.

**H3K27me3 -- Polycomb spreading into active regions.**

| Subcompartment | Late log2FC | Early log2FC |
|---------------|------------|-------------|
| A.1 | **+0.40** | +0.35 |
| A.2 | +0.31 | +0.25 |
| B.1 | +0.13 | +0.22 |
| B.2 | -0.07 | -0.14 |

H3K27me3 (the canonical Polycomb repressive mark) increases in A compartments. This is ectopic Polycomb spreading into active regions, which is biologically coherent with BAP1 loss: without the antagonist, Polycomb marks encroach on territory they normally wouldn't occupy.

**H2AK119ub -- BAP1's direct substrate.**

| Subcompartment | Late log2FC | Early log2FC |
|---------------|------------|-------------|
| A.1 | -0.29 | -0.37 |
| A.2 | -0.09 | -0.13 |
| B.1 | -0.16 | -0.15 |
| B.2 | **-0.32** | **-0.34** |

H2AK119ub *decreases* in the mutant. This is surprising because BAP1 removes H2AK119ub, so knocking BAP1 out should cause H2AK119ub to accumulate, not decrease. Possible explanations:
- Compensatory mechanisms (other deubiquitinases upregulated)
- The BAP1-KO samples reflect a steady-state after adaptation, not the acute response
- The effect may be cell-type or developmental-stage specific
- This is worth discussing with the PI

**Other marks:**
- DNA methylation globally decreases (-0.22 to -0.37 log2FC in all compartments)
- H3K36me3 and H3K4me3 decrease in B compartments, increase slightly in A.2
- ATAC shows very slight increases everywhere (open chromatin gaining)

> **Figure (show to PI):** `outputs/calder2/late/250402_enrichment_combined/` -- the three-panel heatmap with ctrl, mut, and differential side by side. Also `outputs/calder2/early/250831_enrichment_combined/` for comparison.

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

### H2AK119ub goes down -- why?

BAP1 removes H2AK119ub. Knocking BAP1 out should cause H2AK119ub to accumulate. But we see it decrease (log2FC -0.15 to -0.37 across compartments). Three possibilities:

1. **Compensation:** Other deubiquitinases (USP16, MYSM1) may be upregulated to compensate for BAP1 loss, overshooting and removing more H2AK119ub than normal.
2. **Steady state vs acute:** Our samples are from stable BAP1-KO tissue, not acute depletion. The chromatin has had time to reorganize. The initial response (H2AK119ub accumulation) may have been followed by adaptation.
3. **Indirect effect:** If BAP1 loss disrupts PRC1 recruitment or stability at its target sites (not just its catalytic activity), PRC1 may leave those sites entirely, taking H2AK119ub deposition with it.

This needs discussion with the PI.

### H3K27me1 -- what does the massive increase mean?

H3K27me1 increases by 2-2.5 fold in B compartments. This mark is placed by PRC2 (the same complex that places H3K27me3). The pattern -- increase in me1, slight increase in me3 in A regions -- could indicate PRC2 is being redistributed: spreading more broadly (hence more me1, the initial methylation step) rather than deeply methylating specific targets (me3).

Alternatively, if BAP1 loss disrupts the normal Polycomb feedback loop (H2AK119ub recruits PRC2 for H3K27me3), PRC2 may be depositing more "shallow" methylation (me1) without progressing to full me3.

### B.1 as a throughput compartment

B.1 barely changes in total size, but this masks huge turnover: ~850-1250 bins leave B.1 for A.2, while a similar number arrive from B.2. Is this balanced flow a coincidence, or does it reflect a biological mechanism where BAP1 loss creates a "conveyor belt" of chromatin states moving from constitutive toward active?

### HOMER early timepoint polarity

A polarity fix was applied to the HOMER eigenvector for the early timepoint (the PC1 sign was inverted relative to CALDER2 orientation). The A4 script now includes a gene-density polarity check. The re-run showed correct directionality, but the early sankey/dotplot should be interpreted knowing this correction was applied.

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
