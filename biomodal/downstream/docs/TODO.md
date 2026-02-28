Congrats on submitting! Let me pull up your materials and think through next steps.Congrats on getting it submitted, Zara! I've now gone through the abstract, proposal, figures document, and all the methylation context files. Here's my review and then a structured set of further analyses worth considering.

---

## Quick Review of Submitted Materials

**Abstract (249 words):** Reads well. The gap → approach → preliminary finding → summer plan → significance arc is clean. One small note: the phrase "causing locally paired hypermethylation" slightly undersells the dual-mechanism model from your revised conclusions — it reads more like the original TET-blocking model. But for 249 words that's a reasonable simplification, and it's already submitted, so this is just a note for future versions.

**Proposal:** The Background → Goal → Significance flow is strong. I particularly like that you were upfront about the n=2 limitation and sex confounding — reviewers respect that. The methods section is appropriately specific (GLM with batch covariate, Fisher's exact tests, logistic regression for methylation-loop association). The timeline is a bit compressed in the description — "Weeks 1-3 refine pipeline / Weeks 4-6 differential analysis / remaining time downstream" — which is realistic but could have been more granular. Minor issue for a submitted application.

**FIGURES.md:** This is extremely thorough documentation. The figure descriptions are genuinely publication-ready — someone could reconstruct every panel from text alone. This will be valuable for the eventual manuscript.

---

## Further Analyses: What's Worth Doing Next

I'm organizing these by what you can do *now* with existing data vs. what requires the expanded cohort. The "now" analyses are the ones I'd prioritize — they strengthen the story before new replicates arrive and help you design the expanded analysis more effectively.

### Tier 1: High-value analyses doable now

**1. Methylation × Hi-C loop anchor integration**

This is the most obvious cross-modal analysis and it's the core of what your proposal promises. You have 2,910 differential loops and ~6,750 co-significant methylation genes. The question: are genes at differential loop anchors enriched for the coordinated mC↑/hmC↓ pattern compared to non-anchor genes?

Concretely, you'd take your differential loop anchors, extend them to capture nearby gene promoters/bodies (say ±10kb), overlap with DMR gene lists, and run Fisher's exact test for enrichment. Then stratify: do *lost* loops associate with one methylation direction and *gained* loops with another? Your Hi-C analysis already showed H2AK119ub (OR = 10.7) predicts loop loss — does that same set of loci also show hypermethylation? This directly tests the "parallel consequences vs. causal chain" question from your revised conclusions.

**2. A/B compartment-level methylation mapping**

This follows directly from the López-Moyado et al. (2019) framework your revised conclusions highlight. If you have compartment calls from your Hi-C data (eigenvector decomposition), you can ask: does hypermethylation concentrate in A compartment (euchromatin) while the minority hypomethylated genes fall in B compartment? That's the signature of DNMT3A redistribution. If you see it, it's a strong parallel to the TET-KO phenotype and supports the "convergent mechanisms" interpretation. If you *don't* see it, that's also informative — it would suggest the BAP1-mediated mechanism is distinct from direct TET loss.

**3. Baseline 5hmC as a predictor of DMR susceptibility**

Your revised conclusions hypothesize that "the selectivity of which genes are affected likely reflects their baseline 5hmC levels and chromatin accessibility rather than their local H2AK119ub dynamics." You can test this directly: take WT 5hmC levels per gene and ask whether high-baseline-5hmC genes are more likely to show significant hmC loss. A logistic regression with WT 5hmC as predictor and DMR status as outcome would be clean. If baseline 5hmC is a strong predictor (and K119ub is not, which you've already shown), that's direct evidence for the "substrate availability" interpretation.

**4. MeCP2 binding site analysis**

You mention MeCP2 binding affinity in the proposal and your CV shows you've already integrated MeCP2 data. The biological logic is straightforward: 5hmC → 5mC conversion at gene bodies should increase MeCP2 binding (Mellén et al. 2012, 2017 showed MeCP2 binds 5mC with higher affinity than 5hmC). So genes with the coordinated mC↑/hmC↓ pattern should gain MeCP2 occupancy. If you have MeCP2 ChIP/CUT&RUN data from the lab, testing this overlap would add a functional readout to the methylation phenotype — it connects the chemical modification change to an actual downstream reader protein.

**5. Transcription factor motif enrichment at DMRs**

Run a motif enrichment analysis (HOMER `findMotifsGenome.pl` or similar) on your DMR regions. What TF binding sites are enriched in hypermethylated vs. hypomethylated regions? This is relatively quick to run and could reveal whether specific transcription factor families are preferentially affected. If DNMT3A binding motifs or TET-associated cofactor motifs come up, that's mechanistically informative.

### Tier 2: Analyses that deepen the existing story

**6. Single-CpG resolution analysis within gene bodies**

Your current analysis is at gene-level resolution. But within gene bodies, are the methylation changes concentrated at specific sub-features — exon-intron boundaries, splice sites, particular exons? This is relevant because exonic 5hmC has different regulatory roles than intronic 5hmC (the Szulwach et al. 2011 data showed exonic enrichment). If the mC↑/hmC↓ changes cluster at exons specifically, that strengthens the functional interpretation.

**7. RNA splicing analysis**

Your updated BIOMODAL_RESULTS show RNA splicing jumped to the #1 GO term (248 genes, q=3.4e-48) with the deeper sequencing. That's worth following up computationally. If you have RNA-seq data, you could run a differential splicing analysis (rMATS or SUPPA2) to test whether genes with gene body methylation changes also show alternative splicing. Gene body methylation is known to regulate splicing kinetics (slower elongation → more exon inclusion), so mC↑ at gene bodies could directly affect splicing patterns.

**8. Expression-methylation quantitative analysis**

Beyond enrichment, build a quantitative model: does the *magnitude* of mC change or hmC change at a gene predict its expression change? A simple scatter plot of Δ5mC vs. Δexpression (from RNA-seq) per gene, stratified by chromatin state, would show whether the methylation changes are functionally consequential at the transcriptional level. The current ATAC-methylation correlation is weak (ρ = -0.076), but expression might be more tightly coupled.

**9. Developmental timepoint comparison**

If the lab has early-timepoint (P12) methylation data or plans to generate it, the developmental trajectory is important. Your Hi-C shows progressive loop remodeling (165 → 2,910 differential loops from P12 to adult). Does methylation show the same progressive pattern? If the methylation changes are already present at P12, they might be upstream of the loop changes; if they appear later, the loop architecture might be driving methylation rather than vice versa.

### Tier 3: For the expanded cohort (summer)

**10. Sex-stratified analysis**

Once you have n=4 per condition with proper sex balance, run the full analysis with sex as a covariate but *also* run sex-stratified analyses. Sex differences in brain methylation are well-documented, and your current confound means you can't distinguish genotype from sex effects. The expanded cohort lets you quantify how much of the current signal is robust vs. sex-driven.

**11. Batch effect characterization**

With two sequencing batches, you'll need careful normalization. Run PCA on the raw methylation matrix before and after batch correction to quantify the batch effect. If the biological signal (genotype) dominates PC1 and batch is a lower PC, you're in good shape. If batch is PC1, you'll need more aggressive correction (ComBat or similar).

**12. DNMT3A binding prediction**

Without ChIP-seq, you can computationally predict DNMT3A binding sites using known sequence preferences and chromatin features. Train a simple logistic model: does the combination of H2AK119ub level + chromatin accessibility + CpG density predict which gene bodies become hypermethylated? This tests the DNMT3A-UDR recruitment arm of your dual-mechanism model without needing the ChIP-seq data.

---

## What I'd prioritize

If I were choosing a sequence: start with **#1 (methylation × Hi-C)** and **#3 (baseline 5hmC as predictor)** because they directly test the central claims of your proposal and revised model. Then **#2 (compartment mapping)** because it connects to the Rao lab framework and potential collaboration. These three could realistically be done in a few weeks and would substantially strengthen the manuscript.

Want me to sketch out the specific code/analysis plan for any of these?