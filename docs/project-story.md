# BAP1 Hi-C Project: Personal Reference

## What BAP1 does and why we care about 3D chromatin

BAP1 is a deubiquitinase. It removes H2AK119ub -- the monoubiquitin mark that PRC1 deposits on histone H2A at lysine 119. In normal cells, PRC1 puts the mark on and BAP1 takes it off, creating a dynamic cycle. This cycling matters because it's not just about the steady-state level of the mark -- it's about the cell's ability to toggle regulatory elements between states. When BAP1 is gone (our CRISPR knockout in the cerebellum), H2AK119ub accumulates because PRC1 keeps depositing it and nothing removes it.

The lab generated BAP1-KO mice and observed neurodegeneration. The question is: what's actually going wrong at the level of chromatin organization? We have two timepoints -- P12 (early, before severe neurodegeneration) and P60/adult (late, after significant progression) -- with n=3 biological replicates per condition per timepoint. All Hi-C is from cerebellum, processed through both Juicer and Nextflow-HiC pipelines on SDSC Expanse.

The reason to look at 3D chromatin specifically is that we already have CUT&RUN data for the major histone marks (H3K27ac, H3K27me3, H3K4me1, H3K4me3, H2AK119ub), ATAC-seq, RNA-seq, and Biomodal methylation data. The 1D epigenomic picture was partly established. The open question was whether these epigenetic changes actually restructure the 3D genome -- do loops, TADs, and compartments change? And if they do, does that matter for gene expression?

---

## What we found at the loop level

The core pipeline uses mariner (Bioconductor) + edgeR to do replicate-aware differential loop analysis. We run it at three resolutions (5kb, 10kb, 25kb) and merge the results.

At the late timepoint (adult), the numbers are substantial. At 10kb resolution alone: 22,632 loops tested, 3,981 significant at FDR < 0.05. Of those, 2,240 are gained (up in mutant) and 1,741 are lost (down in mutant). After merging across all three resolutions and removing redundant calls, we get 2,910 non-redundant differential loops: 1,723 gained (59%) and 1,187 lost (41%).

At the early timepoint (P12), it's a different story: only 87 significant loops at 10kb resolution (34 lost, 53 gained), and 165 total after merging. So there's roughly an 18-fold amplification from early to late.

### The distance shift

This is the central structural finding. Lost loops and gained loops have very different distance distributions:

- Loops >1 Mb: 413 lost vs 181 gained. That's a 3.3x enrichment for loss at long range.
- Loops 100-500kb: 419 lost vs 1,012 gained. Short-to-mid range is where the gains pile up.
- Lost loops have a median distance of 625 kb; gained loops have a median of 320 kb -- roughly a 2x difference.

All the statistical tests confirm this isn't noise: KS test D=0.279 (p < 2.2e-16), Spearman correlation between distance and logFC is rho = -0.244 (p < 2.2e-16). Longer loops tend to get weaker in BAP1-KO.

We call this "loop rewriting" -- the genome is losing its long-range contacts and replacing them with shorter-range ones.

### Shared anchors

At 212 genomic loci, we found anchors that participate in BOTH a lost loop and a gained loop. The same spot in the genome loses a long-range partner and gains a short-range one. This accounts for 604 of the 2,910 differential loops (about 21%).

The directionality is strong: 83% of shared anchors have the expected pattern (lost loop is longer than gained loop), with a paired Wilcoxon p = 1.17e-20. The median lost distance at these anchors is 1,150 kb and the median gained distance is 340 kb -- a 3.4x ratio.

What's at these shared anchors? H3K27me3 is enriched (OR = 2.04, p = 1.75e-24) and H3K27ac is depleted (OR = 0.68, p = 3.9e-6). These are Polycomb-marked sites, not active regulatory elements. Bivalent promoters are also enriched (OR = 1.74, p = 0.013).

So the switching phenomenon is concentrated at Polycomb-regulated loci. A Polycomb hub loses its long-range partner in another Polycomb domain and instead picks up a shorter-range contact in its local neighborhood.

---

## The mechanism: K119ub in different chromatin contexts

This is where the biology gets interesting and where a lot of the analytical work went.

### The global picture

We have DiffBind results for H2AK119ub: 41,392 peaks tested, with 6,164 significantly up in mutant (K119ub accumulation) and 1,250 down. The asymmetry makes sense -- BAP1 removes K119ub, so when it's gone, the mark accumulates.

The logistic regression for predicting whether a loop is lost gives an OR of 10.70 for K119ub fold-change (p = 1.71e-91), and an OR of 5.36 for log(distance) (p = 1.67e-67). Both K119ub accumulation and distance independently predict loop loss, and the effects are multiplicative.

### Context-dependent effects

Here's where it gets non-obvious. When you split loops by anchor chromatin state and correlate K119ub change with loop strength change, you get different signs:

- **Active anchors** (H3K27ac+): rho = -0.314, p = 7.46e-19 (n=764). More K119ub = weaker loop. This makes sense -- K119ub is a foreign mark at active sites. BAP1 was keeping these elements clean, and now the repressive mark is accumulating where it doesn't belong.

- **Polycomb/Repressive anchors** (H3K27me3+): rho = +0.177, p = 2.82e-09 (n=1,118). More K119ub = *stronger* loop. This seems contradictory at first.

- **Other anchors**: rho = -0.013, p = 0.676 (n=973). No relationship.

### Reconciling the positive rho at Polycomb anchors

The positive correlation at Polycomb sites doesn't mean K119ub is helping loops there. It means K119ub promotes *local* compaction at Polycomb sites (short-range nucleosome-nucleosome interactions) while the *long-range* Polycomb contacts that require higher-order organization are what collapse. The loops being measured at Polycomb anchors with positive rho are the short-range ones that persist or get stronger. The long-range ones are the ones that disappear.

This is a distance-dependent dual effect: K119ub without turnover supports local compaction but disrupts long-range connectivity. The shared anchor data backs this up -- the same Polycomb anchor gains a short contact while losing a long one.

### The enrichment data tells the same story from a different angle

From the differential ChIP enrichment analysis:

- K119ub_up peaks are enriched at long-range lost loops vs short-range gained: OR = 4.87, FDR = 1.86e-90
- K27me3_down peaks are enriched at short-range gained loops vs unchanged: OR = 8.78, FDR = 6.86e-102
- K27me3_up peaks are enriched at long-range lost loops vs unchanged: OR = 3.18, FDR = 3.58e-28

Translation: where K119ub goes up and K27me3 goes up, long-range loops are lost. Where K27me3 goes *down*, short-range loops are gained. The Polycomb landscape is being redistributed, and the loop architecture follows.

### Chromatin state enrichment at loop anchors

The enrichment-by-chromatin-state analysis (14 tests) shows the specificity:

- Active_Enhancer anchors with K119ub_up: OR = 4.72, FDR = 1.81e-09. These are overwhelmingly lost loops -- 73.6% of lost active enhancer anchors overlap K119ub_up peaks vs 37.1% of gained.
- Polycomb anchors with K119ub_down: OR = 2.92, FDR = 2.25e-04. Polycomb sites losing K119ub are associated with lost loops too, but through a different mechanism (loss of PRC1 occupancy).
- Poised_Enhancer anchors with K119ub_up: OR = 2.05, FDR = 2.16e-06.

The unifying principle: BAP1 loss eliminates K119ub turnover. At active sites where BAP1 was keeping things clean, K119ub accumulates and directly weakens enhancer-promoter contacts. At Polycomb domains where K119ub was being dynamically cycled as part of Polycomb body organization, the loss of cycling disrupts the higher-order structure that maintained long-range contacts.

---

## The temporal story

### Early (P12)

Only 165 differential loops. The anchor type distribution is dominated by Repressed_Promoter (~36% of anchors vs ~9-10% at late), with CTCF_Site second (~28%). Active marks are nearly absent (Active_Promoter is 1.8% of Anchor1 at early vs 8.7% at late). The early disruption is concentrated at Polycomb-repressed sites. There are only 3 Active_Promoter anchors at P12 vs 253 at adult.

The direction at P12 is 57% down (lost) -- the dominant early signal is loops weakening. This makes sense as the initial effect of BAP1 loss: existing developmental contacts begin to fail before the secondary reorganization kicks in.

### Late (P60/adult)

2,910 differential loops. The anchor distribution is much more diverse: Poised_Enhancer (~20%), Polycomb (~16%), Active_Enhancer (~11%), CTCF (~27%), Repressed_Promoter (~10%), Active_Promoter (~9%). The disease has spread beyond Polycomb-specific sites to affect the broader regulatory landscape.

The direction reverses: 59% up (gained). This doesn't mean the genome is getting healthier. The gained loops include ectopic contacts, compensatory interactions, and the structural consequence of Polycomb domain reorganization (short-range contacts filling the void left by long-range losses). The 2,240 gained loops at 10kb resolution have weaker individual effect sizes than the 1,741 lost ones -- many small gains vs fewer strong losses.

### The boundary data across timepoints

From the timepoint comparison file: early has 4,349 differential boundaries (18.8%) and late has 4,144 (19.0%). The percentages are remarkably similar -- boundary dynamics don't change much between timepoints even as loop changes amplify 18-fold. But the composition shifts: late has more Strength Change boundaries (1,250 vs 975) and fewer Splits (593 vs 808). The boundaries are responding to the changing loop architecture -- more strengthening as gained loops densify local TADs.

---

## Other structural scales

### TADs

TAD boundaries are moderately affected: about 16-20% are differential at either timepoint. This is notably less dramatic than what happens at the loop or compartment level.

But the TAD boundary changes aren't random with respect to loop changes. Lost loops sit closer to differential boundaries than gained loops (median 45kb vs 75kb, OR = 1.46, p = 4.8e-6). The boundary-loop directional concordance is 69.6% at the late timepoint (chi-squared p = 0.0005) -- when a loop is lost near a boundary, that boundary tends to be control-enriched (stronger in WT), and when a loop is gained, nearby boundaries tend to be mutant-enriched.

The boundary type breakdown is mechanistically informative:

- **Merge** boundaries are 3x enriched near lost loops (OR = 0.32). When a long-range loop collapses, the TAD boundary it spanned becomes unnecessary, and two TADs fuse.
- **Strength Change** boundaries are 2x enriched near gained loops (OR = 2.11). The new shorter-range contacts densify local structure.
- **Split** boundaries are 1.5x enriched near gained loops (OR = 1.48). Some TADs subdivide as new contacts create internal structure.

This is the intermediate-scale link: loop collapse causes TAD merging; loop gain causes TAD strengthening and splitting.

### Compartments

Compartments show dramatic changes at both timepoints. The analysis uses HOMER's `getDiffExpression.pl` on PC1 eigenvectors at 25kb resolution.

**Early (P13):** 101,684 regions analyzed. At standard thresholds (FDR < 0.05, |Diff| > 0.30): 8,154 significant regions covering 203.8 Mb (7.46% of the genome). B->A: 5,282 regions (132.1 Mb); A->B: 2,872 regions (71.8 Mb). At relaxed thresholds: 26,733 regions, 668.3 Mb (24.5% of genome).

**Late (adult):** 104,071 regions analyzed. At standard thresholds: 8,189 significant regions covering 204.7 Mb (7.50% of genome). B->A: 5,485 regions (137.1 Mb); A->B: 2,704 regions (67.6 Mb). At relaxed thresholds: 24,189 regions, 604.7 Mb (22.1% of genome).

The early and late numbers are remarkably similar in total magnitude (~8,150 regions at standard, ~7.5% of genome) but the composition shifts. Early has proportionally more A->B shifts (2,872 vs 2,704) and fewer B->A shifts (5,282 vs 5,485) compared to late. The early 7-category breakdown shows more flipping (A->B flips: 1,694 regions at 42.4 Mb early vs 656 at 16.4 Mb late; B->A flips: 2,426 at 60.6 Mb early vs 1,517 at 37.9 Mb late) and less strengthening. By late, the flipping decreases and strengthening/weakening increases -- the compartment landscape is settling into a remodeled but stable state rather than actively flipping.

The direction is asymmetric at both timepoints: B-to-A shifts outnumber A-to-B roughly 2:1. This makes sense as a consequence of loss of Polycomb repression -- without proper Polycomb-mediated silencing, formerly repressed regions shift toward active.

The 7-category breakdown at standard thresholds (late):
- Flipped B->A: 37.9 Mb (1.39% of genome)
- Strengthened A: 53.8 Mb (1.97%)
- Weakened B: 45.4 Mb (1.66%)
- Flipped A->B: 16.4 Mb (0.60%)
- Strengthened B: 27.3 Mb (1.00%)
- Weakened A: 23.9 Mb (0.87%)

Compartments are defined by the aggregate behavior of many loci, so even modest per-locus changes sum to large compartment shifts. This is why compartments are the most affected scale -- they're sensitive to the widespread Polycomb redistribution that BAP1 loss causes.

### Stripes

Stripes are essentially unaffected. Out of ~200-286 stripes detected, zero reach FDR < 0.05 and only one reaches FDR < 0.10 (stripe_0014 at the Apaf1 locus on chr10, logFC = 0.34, FDR = 0.075). The BCV for stripes is remarkably low (~6-7%), meaning we had good power to detect changes -- there just aren't any.

This is actually informative. Stripes reflect cohesin-mediated loop extrusion, which is mechanistically independent of Polycomb/H2AK119ub. The fact that stripes are preserved confirms that BAP1's effects are channeled through the Polycomb axis, not through the cohesin/extrusion machinery.

### The hierarchy of sensitivity

Putting it together:
1. **Compartments** -- most affected (~7.5% at standard, ~22-24% at relaxed, both timepoints)
2. **Loops** -- primary functional unit (2,910 differential at late, 165 at early)
3. **TAD boundaries** -- moderately stable (16-20% differential), responsive to loop changes
4. **Stripes** -- preserved (0 significant)

This hierarchy makes biological sense for a Polycomb regulator: compartments and Polycomb loops depend on PRC1/PRC2 activity, while CTCF boundaries and cohesin stripes don't.

---

## The enhancer/ABC integration

The Activity-By-Contact model connects enhancers to genes: ABC = (Activity x Contact) / normalization. We computed this for both WT and KO using an identical consensus enhancer universe (75,371 ATAC peaks) to ensure we're measuring real changes in linkage, not artifacts of different enhancer definitions.

### The key finding: unnormalized Delta(AxC) works better

The strongest correlation between enhancer-gene linkage changes and RNA-seq is with the unnormalized sum of Delta(AxC) across all enhancers for each gene: Spearman rho = 0.582, which is very strong for this kind of analysis. The normalized version (standard ABC) gives much weaker correlation because per-gene normalization compresses real activity changes when there's widespread remodeling.

### Enhancer classes tell a mechanistic story

The enhancer subset analysis splits ~55,000 enhancers into four classes based on what changes:

**Activity_Lost** (7,503 enhancers): H3K27ac goes down. These show median loop logFC of -0.088 (loops weakening) and 61.9% RNA-seq concordance. When an enhancer loses its active mark AND its contacts weaken, the target gene is likely to be downregulated. This is the most functionally consequential class.

**Activity_Gain** (2,851 enhancers): H3K27ac goes up. Median loop logFC of +0.066 and 67.2% RNA-seq concordance. Mirror image of Activity_Lost.

**K119ub_Only** (2,479 enhancers): K119ub goes up but active marks don't change. Median loop logFC of -0.054 (loops weakening, but less than Activity_Lost) and 48.7% RNA-seq concordance -- indistinguishable from chance. This is important: K119ub accumulation alone produces a real but sub-functional contact weakening (3.7% reduction, statistically significant) that does NOT translate to detectable gene expression changes. The contact perturbation doesn't cross a functional threshold.

**Stable** (42,864 enhancers): No significant changes. 48.4% concordance (chance level).

### The threshold model

K119ub_Only vs Activity_Lost is the crux. Both involve K119ub accumulation at enhancers. But Activity_Lost also loses H3K27ac, while K119ub_Only keeps its active marks. The contact weakening at K119ub_Only sites (median logFC = -0.054) is real but about 40% weaker than at Activity_Lost sites (median logFC = -0.088). More importantly, the K119ub_Only contact change doesn't produce gene expression effects.

This establishes something like a threshold: K119ub accumulation is upstream of contact disruption, but the contact change must either be large enough or be accompanied by loss of activating marks before the downstream gene is affected. In the adult cerebellum, K119ub alone isn't sufficient.

### Geometric constraint and paired anchors

When you just look at whether loops and ABC connections overlap in the genome, concordance is about 51% -- chance. But when you require the enhancer to be at one loop anchor and the gene's TSS at the other anchor (geometric constraint), concordance jumps to 89.7% (p = 1.67e-48). And three-way concordance (loop + ABC + RNA-seq all agree) is 88.2% (p = 1.69e-45).

At paired enhancers specifically, K119ub correlates with loop loss at rho = -0.401 (p = 5.48e-13) -- 4.5x stronger than the genome-wide K119ub-loop correlation. So the K119ub -> contact disruption -> gene expression chain is strong when you look at geometrically constrained enhancer-gene pairs.

### GO enrichment splits by loop type

The GO terms for genes at long-range lost loops vs short-range gained loops are different:

Long-range lost loops: locomotory behavior, eye development, embryonic organ morphogenesis, sex differentiation, gonad development. These are developmental regulatory programs maintained by long-range contacts.

Short-range gained loops: pattern specification, regionalization, dorsal/ventral patterning, cell fate commitment, kidney development. These are developmental programs associated with shorter-range regulatory architecture.

The functional interpretation: the genome is losing its long-range developmental regulatory connections (which may have been set up during earlier development and maintained by Polycomb) and gaining shorter-range contacts near developmental transcription factor loci.

---

## Reasoning through the tricky parts

### Why the positive rho at Polycomb anchors isn't contradictory

When you first see that K119ub accumulation correlates with loop STRENGTHENING at Polycomb sites (rho = +0.177) and loop WEAKENING at active sites (rho = -0.314), it seems like K119ub is doing opposite things. But the resolution is that these are different populations of loops.

At Polycomb anchors, the loops that survive (or get stronger) are the short-range ones. K119ub promotes local compaction -- nucleosome-nucleosome interactions within ~100kb are strengthened by PRC1-mediated compaction. The long-range Polycomb contacts (>1 Mb) that require higher-order organization (Polycomb bodies, possibly phase separation) are what collapse. But the correlation analysis includes all loops at Polycomb anchors, and the short-range ones that persist or strengthen dominate the positive rho because there are simply more of them after the long-range ones are lost.

The distance shift data at shared Polycomb anchors confirms this: the same anchor has a lost partner at median 1,150 kb and a gained partner at median 340 kb.

At active anchors, the story is simpler: K119ub is a foreign repressive mark. It doesn't belong there. As it accumulates, enhancer-promoter contacts weaken at all distances. The correlation is negative across the board.

### Why compartments change more than TADs

TAD boundaries are primarily maintained by CTCF and cohesin -- structural proteins whose binding isn't directly regulated by H2AK119ub. CTCF binds DNA motifs (we used 114,081 genome-wide CTCF motifs for annotation); cohesin extrudes loops until blocked by CTCF. Neither process depends on Polycomb.

Compartments, by contrast, are defined by the aggregate association preferences of chromatin -- active regions cluster with active regions (A compartment), repressed with repressed (B compartment). Polycomb is a major organizer of the B compartment. When BAP1 loss disrupts Polycomb organization, the association preferences shift. Many formerly B-compartment regions lose their repressive identity and shift toward A. This is why B->A shifts outnumber A->B shifts 2:1.

TADs sit in between: their boundaries are CTCF-dependent (stable), but their internal structure is influenced by the loops within them. When loops collapse and reform at shorter range, TADs merge (lost-loop sites) or densify (gained-loop sites) -- but the boundaries themselves mostly persist because CTCF is still there.

### Why we see more gained than lost loops at the late timepoint

At early (P12), the majority of differential loops are lost (57% down). At late (adult), it reverses: 59% are gained. This reversal isn't because the genome is recovering. There are several things contributing to the gain signal:

1. **Redistribution of contact probability.** When a long-range loop collapses, the interaction "budget" doesn't disappear -- it redistributes locally. The shared anchor data shows this directly: lost contacts at 1,150 kb become gained contacts at 340 kb.

2. **Boundary failure creates ectopic contacts.** As TADs merge, regions that were previously insulated can now interact. These show up as gained loops.

3. **Compensatory interactions.** Some gained loops may represent the cell attempting to maintain gene expression through alternative enhancer-promoter contacts.

4. **Polycomb domain compaction.** At Polycomb sites, K119ub accumulation promotes local compaction even as long-range contacts fail. This creates many short-range gains.

The key insight is that the gained loops at late timepoint are not "healthy" -- they're the structural consequence of architectural collapse. The gained loops have weaker individual effect sizes than the lost ones, are shorter range, and many sit at Polycomb hubs that have lost their long-range partners.

### Why K119ub_Only enhancers don't affect gene expression

This one has real implications for the disease model. K119ub_Only enhancers (2,479 of them) accumulate K119ub but keep their H3K27ac. The contact weakening is real (median logFC = -0.054, p < 2.2e-16) but the RNA-seq concordance is 48.7% -- chance.

Several possible explanations:

1. **Threshold effect.** The 3.7% contact reduction simply isn't enough to change transcription. The enhancer-gene linkage has built-in robustness -- you need to weaken it substantially before the gene notices.

2. **Redundancy.** Each gene has multiple enhancers. If one weakens slightly due to K119ub, others compensate. From the class-level summary: K119ub_Only enhancers have a median of 1 loop and a mean of 2.17, while Activity_Lost has mean 1.79. With multiple enhancers per gene, a slight weakening of one isn't enough.

3. **Active marks still work.** H3K27ac is still there, ATAC accessibility is still there. The enhancer is still "on" -- it's just slightly less connected. The transcriptional machinery can still find and use it.

4. **Developmental context matters.** In steady-state adult tissue, these K119ub_Only perturbations might be tolerable. During developmental transitions (like P12 to adult), where enhancers need to switch states, the inability to remove K119ub could be catastrophic. This would explain the temporal amplification.

### Why stripes don't change but loops do

Stripes and loops both involve chromatin contacts, but they're generated by different mechanisms. Stripes arise from cohesin-mediated loop extrusion -- cohesin loads onto chromatin and extrudes loops until it hits CTCF barriers. This process is mechanical and doesn't depend on histone modifications. Loops, especially enhancer-promoter loops and Polycomb loops, depend on specific protein-protein interactions mediated by chromatin state.

BAP1 loss perturbs chromatin state (H2AK119ub accumulation -> downstream changes in H3K27me3, H3K27ac). This directly affects loops that depend on those marks but doesn't affect the extrusion machinery. The BCV for stripes was ~6-7% with zero significant calls -- we had good power and there was genuinely nothing to find.

This is actually a nice internal control: it tells us the Hi-C data quality is good (stripes are reproducible) and that the loop changes we observe are specifically due to chromatin state perturbation, not some global artifact of BAP1 loss on chromatin compaction or Hi-C library quality.

### The HOMER motif enrichment pattern

The motif analysis at enhancer subsets is informative:

- Activity_Lost enhancers: 369 of 450 motifs enriched (top = Olig2 at p = 1e-120). Massive motif enrichment suggests these enhancers have very specific regulatory grammar.
- K119ub_Only enhancers: 357 of 450 motifs enriched (top = Atoh1 at p = 1e-41). Also highly enriched but less extreme than Activity_Lost.
- Activity_Gain enhancers: only 30 of 450 motifs enriched (top = TCF4 at p = 1e-9). These are less motif-specific.
- K119ub tertile comparison (high vs low): 0 significant motifs. K119ub dose doesn't select for specific TF binding sites.

The Olig2 enrichment at Activity_Lost enhancers is interesting for a cerebellar context -- Olig2 is involved in glial and neuronal development. Atoh1 at K119ub_Only sites is a cerebellar granule cell specification factor. These are tissue-relevant TF binding sites being affected.

### How the different data types connect

The chain runs like this:

1. BAP1 is gone -> K119ub accumulates (DiffBind: 6,164 K119ub_up peaks)
2. K119ub accumulation at active enhancers weakens E-P contacts (rho = -0.314 at active anchors; OR = 10.70 in logistic regression)
3. K119ub accumulation at Polycomb domains disrupts long-range contacts while promoting local compaction (shared anchors: 83% directionality, Polycomb-enriched)
4. Contact changes propagate to TAD boundaries (69.6% concordance, merges at lost sites, strengthening at gained sites)
5. Aggregate effects shift compartments (7.5-22% of genome, net B->A)
6. Where contact weakening is strong enough AND active marks are lost (Activity_Lost class), gene expression changes follow (61.9% concordance)
7. Where contact weakening is weak and active marks persist (K119ub_Only class), gene expression is unaffected (48.7%, chance)
8. Three-way concordance at geometrically constrained enhancer-gene pairs: 88.2%

The temporal dimension adds that this cascade starts small (165 differential loops at P12, concentrated at Repressed_Promoter sites) and amplifies massively by adulthood (2,910 loops, broad chromatin type involvement, 44% compartment shifts), consistent with a progressive regulatory breakdown rather than an acute event.

---

## Key loci

- **Syt1/Nav3** (chr10): The most impacted region in the data. Essentially collapsed by P60.
- **Apaf1** (chr10): The only significant stripe change -- an apoptosis regulator gaining a stripe in mutant. Interesting given the neurodegeneration phenotype.
- **Zbtb20**: Lost loop at Active_Promoter-Active_Enhancer (logFC = -0.389, FDR < 1e-24 for RNA), one of the stronger DEG-loop associations.
- **Nfib**: Second-highest combined structural score, lost loops, -0.57 log2FC RNA.
- **Tgfbr1**: Gained loops + boundary change + ABC change -- a 3-layer hit.
- **Col15a1**: Strongest ABC effect (delta unnorm = 3.86), multi-layer disruption.

---

## Numbers to have in your head

- 2,910 differential loops at adult (1,723 gained, 1,187 lost)
- 165 at P12 (18x amplification)
- Lost median 625kb, gained median 320kb (1.95x ratio)
- >1 Mb loops: 3.3x enriched for loss
- K119ub OR for predicting loop loss: 10.70
- Active anchor K119ub-loop correlation: rho = -0.314
- Polycomb anchor K119ub-loop correlation: rho = +0.177
- 212 shared anchor hubs, 83% directionality
- TAD boundaries: 16-20% differential, 69.6% directional concordance with loops
- Compartments: ~7.5% (standard) to ~22-24% (relaxed) of genome shifted at both timepoints
- Stripes: 0 significant
- ABC 3-way concordance: 88.2% with geometric constraint
- Activity_Lost concordance: 61.9%; K119ub_Only: 48.7% (chance)
- Unnormalized DABC-RNA rho: 0.582
- 457 DEGs at loop anchors
- K119ub at paired enhancers: rho = -0.401
