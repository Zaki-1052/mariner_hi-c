Alright, I've read through everything carefully. Here's the consolidated plan, incorporating your feedback and corrections throughout.

---

## Scoping Decisions That Shape Everything

The biggest structural decision is that the non-CG methylation / MeCP2-as-ubiquitin-reader story may become its own paper. That means this paper's core narrative is about **BAP1 loss → H2AK119ub accumulation → coordinated 5mC↑/5hmC↓ at gene bodies → disrupted neuronal gene regulation**, with MeCP2 present in a supporting capacity (developmental trajectory, redistribution) rather than as the headline claim. The title needs to reflect that shift — centering the dual methylation phenotype rather than MeCP2 restructuring.

This also means the Sections 51–58 non-CG analysis, the Section 62 "MeCP2 reads chromatin not methylation" regression, and the Section 67 linchpin (359 genes, K119ub at unmethylated loci) are strong results but may belong in paper two. They shouldn't drive figure composition here.

What stays central: the methylation phenotype itself, H2AK119ub as the upstream driver, TET-impediment vs DNMT3A-dehydroxymethylase mechanism, chromatin geography (compartments, Polycomb exclusion, chromatin states), neuronal gene specificity and stoichiometry, and the developmental MeCP2 trajectory as supporting evidence for progressive neurodevelopmental disruption.

---

## Figure-by-Figure Plan with New Section Assignments

### Figure 1: The Coordinated Gene-Body Methylation Phenotype

**What it shows:** BAP1-KO produces a coordinated, reciprocal 5mC↑/5hmC↓ shift concentrated at gene bodies, with flat total modified cytosine — the signature of blocked TET-mediated demethylation.

**Panels to consolidate:**

The bulk genome-wide shift from Section 64 (the paired replicate deltas with Cohen's d values) establishes the systemic phenotype. The dual volcanoes from Section 04 show the directional asymmetry across thousands of genes. The quadrant scatter from Section 05 shows the 78.7% coordinated pattern. The region breakdown from Section 03 shows gene-body dominance. The effect-size distributions from Section 07 show the mirror-image magnitude. One genome browser track from Section 46 (Syt1 is the natural choice as the most-affected gene at +17.3% mC / −15.0% hmC) anchors the abstract numbers to a real locus.

**New section needed → Section 79: `figure1_methylation_phenotype`**

This script reads existing section TSVs and Section 46 Gviz output, builds a unified multi-panel patchwork with consistent color coding (standardize on the `COLORS` palette from `_shared_config.R` for mC/hmC direction throughout all figures). Key design choices: the bulk shift panel should show individual replicate pairs connected by lines (not just bars), since the within-pair consistency is the whole point of those huge Cohen's d values. The quadrant scatter should use the same hyper/hypo color scheme as the volcanoes. Gene-body vs other-region comparison should be a simple fraction bar, not a table.

No new computation needed — all data exists in Sections 03, 04, 05, 07, 46, and 64 TSVs.

### Figure 2: H2AK119ub Drives Hypermethylation at Active Chromatin

**What it shows:** K119ub gain is the strongest predictor of hypermethylation, but the effect is restricted to active/euchromatic regions — Polycomb heterochromatin is protected.

**Panels to consolidate:**

The multi-mark logistic model from Section 33 (K119ub OR=4.707, far above any other mark). The chromatin-state direction split from Section 10 (Active_Promoter 93% hyper vs Repressed_Promoter 94% hypo — the mirror image). The A/B compartment enrichment from Section 29 (mC-hyper in A compartment OR=13.64). The Polycomb exclusion from Section 30 (Polycomb targets excluded from hypermethylation, OR=0.063). The subcompartment stratification from Section 66 (A.1 hypermethylation-dominated, B.2 hypomethylation-dominated).

Your point about wanting to "actually show the data" rather than just reporting that you logistically regressed everything is important here. The Section 33 logistic model is analytically strong, but the figure should lead with the data that *motivated* the model — the chromatin-state split, the compartment geography — and use the model coefficients as quantitative confirmation, not as the primary visual. The forest plot of odds ratios is supporting evidence for what the reader can already see in the chromatin-state bars and compartment enrichments.

**New section needed → Section 80: `figure2_k119ub_chromatin_geography`**

Reads from Sections 10, 29, 30, 33, and 66. The expanded chromatin-state heatmap you mentioned would fit here — showing per-state hyper/hypo rates alongside K119ub, K27me3, K27ac signal levels to give the reader the full picture in one panel. The honest correction from Section 17 (58.2% of DMR genes have no K119ub peak; conditional enrichment +14.2 pp) should be incorporated as context, not hidden — but presented as part of the data, not as a separate "correction."

The 6-mark model from Section 41 (K119ub OR=10.29, H3K36me3 OR=1.00 ruling out DNMT3B) could go here or supplemental depending on space. The DNMT3A vs DNMT3B attribution (19:1 ratio) is worth showing since it narrows the enzymatic mechanism.

### Figure 3: TET-Impediment Mechanism and Dehydroxymethylase Activity

**What it shows:** The coordinated phenotype is driven by impaired TET-mediated demethylation, not indiscriminate de novo methylation — and at neuronal genes specifically, the stoichiometry points to DNMT3A acting as a dehydroxymethylase (direct 5hmC→5mC conversion).

**Panels to consolidate:**

The baseline-5hmC predictor from Section 23 (AUC 0.762 — genes with more TET substrate show the biggest changes). The demethylation ratio by chromatin state from Section 22 (strongest impairment at Active_Promoter/Active_Enhancer, reversed at Repressed_Promoter). The stoichiometry scatter from Section 61/78 showing the genome-wide slope (−0.959, distinguishable from −1) versus the neuronal gene slope (−1.0, consistent with stoichiometric conversion). The total methylation conservation test from Section 78 (total methylation rises at coordinated genes but *decreases* at neuronal genes with the unbiased broad set).

On the TET-KO comparison from Section 26 — you're right that it's not one-to-one with your model. The "3.3% attenuation / dimmer not off switch" framing is conceptually useful but the comparison to a triple-KO in a different tissue/system is inherently imperfect. I'd suggest making this supplemental rather than a main figure panel, unless Cole feels strongly about it. The mechanism stands on the baseline-5hmC prediction, the demethylation ratio, and the stoichiometry without needing the external comparison.

The ROC curves from Section 24 (TET-impediment AUC 0.793 vs DNMT3A-recruitment AUC 0.696) are analytically decisive, but your instinct that there might be a cleaner way to show this is worth thinking about. One option: rather than ROCs, show the actual data that the models are fitting — scatter the genes by their baseline 5hmC level (the TET model's top predictor) vs their observed methylation change, and let the reader see the dose-response directly. The model summary statistics can go in the text or supplemental. That way the figure shows data, not model output.

**New section needed → Section 81: `figure3_tet_mechanism_stoichiometry`**

Reads from Sections 22, 23, 24, 61, and 78. The stoichiometry result from Section 78 (neuronal slope = −1.0 vs genome-wide −0.959) is important to show with the *unbiased broad neuronal set* as the final result — not framing it as a correction of Section 61, just presenting it as the data. The coordinated-gene slope of −1.29 (steeper, meaning net mC gain beyond conversion) is the contrast that shows different loci undergo different mechanisms.

New computation: consider adding a panel that directly shows the baseline-5hmC dose-response — bin genes by their wildtype 5hmC level and plot the resulting hmC change per bin. This is essentially what the Section 23 model captures, but as visible data rather than an ROC curve. Could be done within this same script.

### Figure 4: Neuronal Gene Enrichment and Chromatin Remodeling

**What it shows:** Neuronal-identity genes are constitutively H2AK119ub-marked (the natural BAP1 substrate), and synapse/axon genes specifically undergo selective Polycomb de-repression upon BAP1 loss.

**Panels to consolidate:**

The constitutive K119ub enrichment from Section 72 (top-decile OR=1.70, dose-response across deciles). The synapse/axon specificity from Section 76 (selective K27me3 loss, not extra ATAC or K27ac gain). The GO enrichment from Section 08 (hyper genes enriched for neuronal/synaptic/neurodegenerative pathways). The chromatin remodeling interaction from Section 73 (K119ub × neuronal interaction significant for ATAC and K27ac at K119ub-high tier).

The gene-set overlap from Section 74 (Coordinated × MeCP2-Up OR=5.16 >> Neuronal × MeCP2-Up OR=1.73) could go here if MeCP2 stays in this paper, or be deferred to paper two. The key point — that the effect is "preferentially/disproportionately" at neuronal genes, not "exclusively" — should shape the text regardless.

**New section needed → Section 82: `figure4_neuronal_specificity`**

Reads from Sections 08, 72, 73, 76, and potentially 74. The K119ub decile dose-response (Section 72f) paired with the synapse-specific K27me3 loss (Section 76c) tells a clean story: neuronal genes sit in Polycomb-marked chromatin, BAP1 loss amplifies K119ub there, and synapse/axon genes respond with selective de-repression of that Polycomb repression.

One thing to think about: Section 73's result that K119ub-high neuronal genes gain *more* accessibility and K27ac (not less) was noted as being opposite to a naive heterochromatin-shift prediction. The figure should present this as what it is — de-repression of poised neuronal promoters — rather than trying to force-fit it into a compaction narrative.

### Figure 5: MeCP2 Developmental Trajectory (if MeCP2 stays in this paper)

**What it shows:** MeCP2 binding increases progressively over neurodevelopment, and BAP1 loss amplifies this developmental ramp — more aging-UP peaks, higher fold at shared loci, mutant-unique aging genes enriched for neuronal programs.

**Panels to consolidate:**

The aging peak counts from Section 77 (mutant 2.1× more aging-UP peaks). The gene-level Venn from Section 77 (2,620 shared / 1,654 mutant-unique / 288 ctrl-unique). The shared-peak fold comparison (median fold 2.24 vs 1.83, mutant climbs higher). The GO enrichment of mutant-unique aging genes (435 terms, 49 neuronal).

The MeCP2 redistribution from Section 75 (UP peaks 52% distal-intergenic) and annotation split could go here or supplemental, as you mentioned. If MeCP2 is being downweighted in this paper, showing the developmental amplification as a consequence of the methylation phenotype (rather than as the headline) makes sense — it's evidence that the methylation changes have functional consequences that worsen over development, which ties back to the progressive phenotype (165 differential loops at P12 → 2,910 in adulthood from the Hi-C arm).

**New section needed → Section 83: `figure5_mecp2_aging`** (contingent on MeCP2 staying in paper)

Reads from Sections 75 and 77. If the scope narrows further and MeCP2 becomes entirely paper two, this figure gets cut and the developmental progression evidence comes from the Hi-C arm instead.

**New analysis needed → Section 84: `mecp2_aging_methylation_overlay`**

You agreed this should be a new section. The question: do the 1,654 mutant-unique aging MeCP2 genes co-localize with coordinated mC↑/hmC↓ gene bodies? This is a Fisher's exact test of the mutant-unique aging gene set (from `77_aging_overlap.tsv`) against the coordinated gene set (from `_shared_config.R`'s `mc_dmr`/`hmc_dmr` joined data). If they overlap significantly, it directly connects the methylation phenotype to the developmental MeCP2 amplification. If they don't, that's also informative — it would mean MeCP2 is being recruited to loci independent of where methylation changes, which feeds back into the "MeCP2 reads chromatin not methylation" story for paper two. Either way, the result is useful. This is probably 100–150 lines of R plus a Venn/bar panel and a Fisher table.

### Figure 6: QC and Sample Structure

**What it shows:** PCA of samples by condition/sex/batch, correlation structure, and quality metrics confirming the data is clean and the sex covariate is properly handled.

**Panels needed:**

The sample correlation from Section 02. PCA colored by condition, sex, and batch from Section 42. The chrX analysis from Section 43 showing the phenotype is autosomal.

**New section needed → Section 85: `figure6_qc_pca`**

You mentioned needing to show the PCA for sex. Section 42 already produces per-sample methylation PCAs, but they're colored by condition/batch. Running the same PCA with sex as the color variable would directly address whether samples cluster by sex, and whether the sex covariate is doing what it should. This is a small addition to the existing Section 42 framework — recolor the same PCA coordinates. If samples don't cluster by sex (after the covariate), that confirms the Methods statement and resolves the R3 tension. If they do cluster partially, that's worth showing because it validates the decision to include the sex covariate in run-5.

---

## Summary of New Sections

| New Section | Purpose | Reads From | New Computation |
|---|---|---|---|
| **79** | Figure 1: coordinated methylation phenotype composite | Sec 03, 04, 05, 07, 46, 64 TSVs | None — consolidation only |
| **80** | Figure 2: K119ub + chromatin geography composite | Sec 10, 17, 29, 30, 33, 41, 66 TSVs | Expanded chromatin-state heatmap |
| **81** | Figure 3: TET mechanism + stoichiometry composite | Sec 22, 23, 24, 61, 78 TSVs | Baseline-5hmC dose-response binned panel |
| **82** | Figure 4: neuronal specificity composite | Sec 08, 72, 73, 76 TSVs | None — consolidation only |
| **83** | Figure 5: MeCP2 aging trajectory composite (conditional) | Sec 75, 77 TSVs | None — consolidation only |
| **84** | MeCP2 aging × methylation DMR overlay | Sec 77 + coordinated gene set | Fisher's test, Venn, overlap characterization |
| **85** | QC/PCA with sex stratification | Sec 02, 42, 43 + DUET sample data | Sex-colored PCA recomputation |

---

## Phasing

**Phase 1 (before tomorrow's meeting):** Finalize the figure plan and panel assignments so you can present an outline. No code needed — this document plus your notes should be sufficient to discuss with Cole which figures make the cut and whether MeCP2 stays or goes to paper two.

**Phase 2 (this week):** Write Sections 79 and 80 first. Figure 1 (the methylation phenotype) and Figure 2 (K119ub + chromatin geography) are the backbone of the paper regardless of MeCP2 scoping decisions, and they're pure consolidation of existing TSVs with no new computation. Getting these to publication quality first gives you concrete material to build the manuscript around.

**Phase 3 (next week, parallel with writing):** Section 81 (TET mechanism), Section 82 (neuronal specificity), and Section 84 (the new MeCP2 aging × methylation overlay). Section 81 has a small new computation (baseline-5hmC binned dose-response). Section 84 is the only genuinely new analysis.

**Phase 4 (after scoping decision on MeCP2):** Sections 83 and 85, which depend on whether MeCP2 stays in this paper and how the QC figure is positioned.

---

## Things Explicitly Not in This Plan

Per your feedback: no summary model schematic, no cell-type deconvolution (wet lab limitation), no separate `scripts/figures/` directory (continuing with sections), no framing of the Section 78 broad-set result as a "self-correction" (just present it as the data), and no heavy reliance on the TET-KO comparison as a main figure panel (supplemental at most). The non-CG methylation analysis (Sections 51–58) and the MeCP2-reads-chromatin regression suite (Sections 62, 67, 71) are strong results held for paper two.

---

## Addendum: Figure 6 — MeCP2 Responds to Chromatin State, Not Methylation

You're right to keep this in paper one. The distinction is clean: paper two is about non-CG methylation as a binding modality, while this result is about CG methylation being *insufficient* to explain MeCP2 redistribution — K119ub is the real driver. That's a conclusion drawn entirely from data you already have (CG DUET + MeCP2 CUT&RUN + DiffBind histone marks), and it directly extends the paper's main story about H2AK119ub as the upstream effector.

### What this figure shows

MeCP2 binding changes in BAP1-KO track chromatin state (primarily K119ub gain), not DNA methylation change. Methylation is necessary but not sufficient — most hypermethylated genes don't recruit more MeCP2, and genes that gain MeCP2 without any methylation change do so because K119ub accumulates there.

### Panels

The Section 62 regression is the quantitative core: CG methylation explains less than 2% of MeCP2 binding variance while chromatin marks explain about 25%, with K119ub as the top standardized predictor (β=+0.199). But per your point about showing the data rather than just saying "I regressed everything," the figure should lead with what the reader can see directly. The Section 65 scale comparison (51.4% of genes change methylation but only 9.8% change MeCP2 binding, and only 12.6% of hypermethylated genes recruit more MeCP2) immediately shows the dissociation without any modeling. The Section 67 result (359 genes where MeCP2 rises with no methylation change, K119ub gained at 72.8% of them vs 45.9% background, OR=3.15) is visually concrete — you can show the K119ub violin at those loci versus background and the reader sees it. The Section 71 variance partition (K119ub explains 7.3% of MeCP2 fold, methylation ratio explains 0.0%) is a clean bar that quantifies the hierarchy. The Section 60 mirror-image profiles (MeCP2-Up loci show K119ub↑/5mC↑/5hmC↓/K27ac↓/K27me3↑ while MeCP2-Down loci mirror everything except K119ub stays flat) could be the anchor panel showing the actual epigenetic landscape at these genes.

The regression R² comparison from Section 62 then comes in as confirmation — the reader has already seen the dissociation in the raw data, so the model summary is reinforcement rather than the primary evidence.

### New section needed → Section 86: `figure6_mecp2_chromatin_reader`

Reads from Sections 59, 60, 62, 65, 67, and 71. No new computation — this is consolidation of existing TSVs into a unified figure. The design priority is making each panel show data (violins, scatters, bars of actual signal values) with the statistical summaries as annotations rather than as standalone model-output panels.

This slots into **Phase 3** alongside Sections 81 and 82, since it depends on the same underlying data and doesn't require any new analysis.

---

### On the framing points you shared

Those two statements read like the paper's two main conclusions, and they fit together well. Point 1 (H2AK119ub regulates 5mC/5hmC balance broadly, not just at neuronal loci) is supported by the genome-wide coordinated phenotype in Figures 1–3 and by the Section 78 result that the methylation effect is actually broader than neuronal genes — neuronal genes are disproportionately affected because they're constitutively K119ub-rich, but the methylation shift is everywhere in active chromatin. Point 2 (MeCP2 responds to chromatin shift, not methylation) is this new Figure 6.

The relationship between them is that the methylation story is the biology (Figures 1–4), and the MeCP2 finding (Figure 6) is a mechanistic insight about how that biology is *read* by the cell. That framing keeps MeCP2 in the paper without making it the headline — it's a downstream consequence that reveals the chromatin-level mechanism, which is the H2AK119ub axis established in Figure 2.

Whether the MeCP2 aging trajectory (Section 83/Figure 5) stays as its own main figure or folds into supplemental supporting Figure 6 is probably a decision for the meeting tomorrow. It depends on how much real estate the paper can afford and whether the developmental amplification (2.1× more aging-UP peaks) is compelling enough as a standalone result or works better as supporting context for the chromatin-reader claim.