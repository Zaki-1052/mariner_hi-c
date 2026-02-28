Good TODO document — the structure mirrors the Hi-C one well and the priority ordering makes sense. I spent time thinking through the biology carefully, and there are several dimensions I think are missing or underexplored. I'll go from what I think is most impactful to more speculative.

---

## Analyses I'd add

### 1. 5hmC/(5mC+5hmC) ratio as a demethylation efficiency index

This might be the single cleanest metric you're not computing. Rather than treating mC and hmC as two separate signals (which is what all 21 sections do), compute the ratio per gene in WT vs KO. This ratio directly represents TET conversion efficiency — what fraction of modified cytosines at a locus have been successfully oxidized. A decrease in this ratio in KO is the most direct measure of impaired active demethylation, gene by gene.

Why this matters beyond what you already have: the coordinated pattern (mC↑/hmC↓) is compelling but binary. The ratio gives you a continuous, per-gene "demethylation activity score" that you can then use as the response variable in all your downstream models. It collapses two marks into one biologically interpretable number, and it sidesteps the issue of mC and hmC having different dynamic ranges and variance structures. You could also compare the WT distribution of this ratio to published TET-KO data from the Rao lab — if your BAP1-KO ratio shift phenocopies direct TET loss, that's a strong claim.

### 2. Methylation entropy / read-level heterogeneity

Your DUET evoC data is at 6bp resolution with deep coverage (~400M reads/sample). That means you have multiple reads covering each CpG, which means you can compute methylation entropy — a measure of how heterogeneous the methylation state is across reads at a given locus.

This distinguishes two very different biological scenarios that your current gene-level means can't separate. Scenario A: every cell shifts from 40% to 42% methylation (uniform mild increase, consistent with a graded DNMT3A recruitment signal). Scenario B: half the cells stay at 40% and the other half jump to 44% (bimodal, consistent with stochastic TET failure or cell-type heterogeneity). These look identical at the mean level but have different entropy signatures. Scenario B would show increased entropy in KO, while Scenario A wouldn't.

This is computationally tractable from the Zarr stores — you'd extract read-level methylation calls at individual CpGs and compute Shannon entropy or epipolymorphism per window. It's also directly relevant to the cell-type question: Math1-Cre targets granule cell progenitors specifically, so if the methylation changes are clean (low entropy), that suggests the targeted cells are uniformly affected. If entropy increases, it suggests mosaicism in the Cre-mediated knockout or cell-type mixing effects.

### 3. CTCF binding site methylation → 3D architecture connection

This is a glaring missing link between your methylation and Hi-C stories. CTCF binding is methylation-sensitive — CpG methylation at CTCF motifs blocks binding, and CTCF is the primary insulator/loop anchor protein. You already have 114,081 CTCF motif positions (`ctcf_motifs_mm10.bed`) from the Hi-C pipeline. The question: do CTCF motif sites at differential loop anchors show methylation changes?

If CTCF sites at lost loop anchors gain methylation, that's a direct mechanistic link between your two arms — methylation changes could be *causing* loop loss by disrupting CTCF binding, rather than both being parallel consequences. If CTCF site methylation is unchanged at loop anchors, that rules out this mechanism and supports the "parallel consequences" interpretation from your revised conclusions. Either way, it's informative and directly tests a specific causal hypothesis.

You can do this with your existing data: extract methylation levels at CTCF motif CpGs from the Zarr stores, stratify by whether the motif is at a differential loop anchor vs. a stable loop anchor, and compare.

### 4. Repeat element / transposable element methylation

This is a dimension your entire analysis pipeline currently ignores. Gene bodies and promoters are the focus, but repeats comprise ~40% of the mouse genome, and they're major DNMT3A substrates. Three reasons this matters for your model:

First, if the DNMT3A-UDR recruitment arm is correct, DNMT3A is being pulled toward H2AK119ub-marked regions. That means DNMT3A may be depleted *from* its normal repeat targets. You could see paradoxical hypomethylation at repeats (especially young LINEs and LTRs that depend on active DNMT3A maintenance) simultaneous with gene body hypermethylation. This would parallel the López-Moyado DNMT redistribution model but with different geography.

Second, repeat de-repression has functional consequences — activation of LINE1 elements or endogenous retroviruses can cause genomic instability and neuroinflammation. If you see repeat hypomethylation, that's a disease-relevant finding.

Third, this is quick to compute. Extract methylation at RepeatMasker-annotated elements (available for mm10), stratify by repeat family (LINE, SINE, LTR, DNA transposons), and compare WT vs KO. You could even look at repeat age/divergence to test whether young (actively repressed) repeats are more affected than ancient (epigenetically inert) ones.

### 5. Polycomb target gene enrichment (systematic, not just chromatin state)

Your Section 10 analyzes chromatin states, and Section 14/17 look at K119ub peaks. But you haven't done a systematic Polycomb target gene analysis using published PRC1/PRC2 target gene lists. This is biologically central — BAP1 is part of PR-DUB, which opposes PRC1. The prediction from the dual-mechanism model is nuanced: classic Polycomb target genes (repressed, H3K27me3-marked) should *not* be the primary hypermethylation targets, because they're already in heterochromatin and inaccessible to DNMT3A. Instead, the model predicts that genes which are normally active but gain ectopic H2AK119ub (due to the genome-wide spreading after BAP1 loss) should be the hypermethylated ones.

You can test this by comparing your DMR list against published Polycomb target gene sets (there are well-curated lists for mouse ESCs and neural progenitors). If classic Polycomb targets are *not* preferentially hypermethylated but non-Polycomb genes are, that's strong evidence for the "ectopic H2AK119ub at active genes" arm of the model. If classic Polycomb targets *are* hypermethylated, that would suggest a different mechanism (Polycomb-mediated compaction driving methylation).

### 6. Spatial spreading / chromosomal autocorrelation

Your analysis treats each gene independently. But H2AK119ub is known to spread in cis along chromatin (Conway et al. 2021 showed this). If the methylation changes are downstream of H2AK119ub spreading, you'd expect spatial autocorrelation — hypermethylated genes should cluster in chromosomal neighborhoods rather than being randomly distributed. A runs test or sliding-window autocorrelation along each chromosome would test this. If you see methylation "domains" (contiguous blocks of hypermethylated genes), their boundaries might correspond to topological features like TAD boundaries or CTCF sites, connecting back to the architecture story.

Conversely, if DMRs are randomly scattered across chromosomes with no spatial structure, that argues against a spreading mechanism and more for a gene-intrinsic susceptibility model (consistent with the "baseline 5hmC as predictor" hypothesis in your TODO Section 3).

### 7. H3K36me2/3 connection (conceptual gap)

This is a data availability question, but it's worth flagging as a conceptual gap. H3K36me2 is a key player you're not discussing. Recent work (Weinberg et al. 2019, Nature; Streubel et al. 2018) showed that H3K36me2 antagonizes PRC2 and recruits DNMT3A to intergenic regions, while H3K36me3 (deposited by SETD2 in actively transcribed gene bodies) recruits DNMT3B. The interplay between H2AK119ub and H3K36me2/3 could be important: if BAP1 loss and H2AK119ub accumulation alter H3K36me2 domains, that's another route to DNMT redistribution. If the lab has or can generate H3K36me2/3 CUT&RUN, it would directly test this. Even without the data, you should probably discuss this in the manuscript as an alternative/additional recruitment mechanism.

### 8. Cross-species comparison with human BAP1 tumor methylation

Field et al. (2019) profiled 5mC (not 5hmC) in BAP1-mutant uveal melanoma. Your bibliography notes this. A direct comparison — which genes are hypermethylated in both your mouse cerebellum KO and human BAP1-mutant tumors — would be extremely valuable for the manuscript. Shared hits suggest a conserved BAP1-dependent methylation program; divergent hits suggest tissue-specific responses. This is straightforward to do with the published data, probably just a few hours of work with ortholog mapping (mouse-human) and overlap testing.

### 9. Non-CG methylation — expanding on your suggestion

Your instinct is right. Let me add specificity: the most important non-CG context is mCA (methylation of CA dinucleotides), which in brain:

- Accumulates specifically in postmitotic neurons (Lister et al. 2013)
- Is written almost exclusively by DNMT3A (not DNMT3B)
- Is bound by MeCP2 with higher affinity than mCG in some contexts
- Is associated with gene repression (opposite to gene body mCG)

If DNMT3A is being redistributed via UDR-H2AK119ub, mCA should change at the same loci that show mCG gain. This is a specific prediction of the DNMT3A recruitment arm that's distinct from the TET impediment arm (TET primarily acts on CG context). So if you see correlated mCA and mCG changes, that's evidence for DNMT3A recruitment specifically. If mCG changes but mCA doesn't, that's more consistent with TET impediment alone.

The sensitivity concern is real, but even at <1% levels, if you're comparing KO vs WT, you need detectable *differences*, not high absolute levels. With ~400M reads per sample, you should have reasonable coverage at CA sites in gene bodies even if individual-site methylation is low.

### 10. Cell type deconvolution

Math1-Cre targets cerebellar granule cell progenitors, but cerebellum contains Purkinje neurons, interneurons, Bergmann glia, astrocytes, oligodendrocytes, etc. Your bulk methylation signal is a weighted average across all cell types. Reference-based deconvolution using published single-cell methylation atlases (or even single-cell RNA-seq atlases with methylation proxies) could estimate the cell-type composition of your samples and, more importantly, whether the methylation changes are driven by composition shifts vs. within-cell-type changes. If the KO cerebellum has fewer granule neurons (due to degeneration), some of your "hypomethylation" could actually be a change in cell-type proportions rather than a change in methylation within surviving cells.

This is particularly relevant for interpreting the minority hypomethylated genes (the 15.4% discordant pattern) — are those genuinely hypomethylated in neurons, or is it a deconvolution artifact from cell death/loss?

---

## Things I'd deprioritize or fold into existing sections

- **DNA methylation valleys (DMVs/canyons):** Interesting but probably a sub-analysis of the CpG density work that would emerge from the DNMT3A prediction model (TODO Section 11). Don't need a separate section.
- **Super-enhancer analysis:** Your Active_Enhancer category already captures most of this, and super-enhancers are a small fraction. Could be a panel in Section 1 rather than standalone.
- **Imprinted gene analysis:** Quick sanity check, but unlikely to be a major finding. Could be a supplementary figure.

## Summary of what I'd add

In rough priority order for biological impact:

1. **5hmC/(5mC+5hmC) demethylation ratio** — cleanest single metric, should probably be computed for every gene and used as a primary variable going forward
2. **CTCF site methylation at loop anchors** — directly links your two project arms with a causal mechanism test
3. **Repeat element methylation** — entirely missing dimension, tests DNMT3A redistribution prediction
4. **Spatial autocorrelation** — tests spreading vs. gene-intrinsic susceptibility
5. **Methylation entropy from read-level data** — distinguishes uniform vs. stochastic mechanisms
6. **Polycomb target gene systematic analysis** — central to the biology, currently only indirect
7. **Non-CG (mCA specifically) × MeCP2 × DNMT3A** — expanding your suggestion with the specific DNMT3A prediction
8. **Cross-species comparison with Field et al. 2019** — low effort, high manuscript value
9. **Cell type deconvolution** — important caveat analysis
10. **H3K36me2/3 conceptual discussion** — manuscript framing even if no data

The first four are things I'd genuinely want done before writing the manuscript. The rest are "would strengthen the paper" tier.