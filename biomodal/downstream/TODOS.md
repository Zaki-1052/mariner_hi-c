The core challenge here is that you have 78 sections of rigorous exploratory analysis, but the paper needs maybe 6–8 consolidated figures where each panel earns its space by advancing the narrative, not just reporting a test. Let me walk through how I'd think about organizing this, what new work would help, and where the current data has tensions that need resolving before you write.

## The Big-Picture Figure Architecture

Based on the paper outline (R1–R5), the synthesis doc, and the section documentation, I think the paper naturally supports about 7 main figures plus a model schematic. The key organizing principle is that each figure should answer one question the reader has at that point in the story, and every panel within it should contribute to that answer without requiring the reader to mentally stitch together unrelated analyses.

Here's how I'd map the existing work onto a figure plan, and where I think new consolidated scripts would help.

### Figure 1: The global methylation phenotype

This is your "what happened" figure. It establishes the BAP1-KO system and the coordinated mC↑/hmC↓ signature before you explain mechanism.

The strongest panels for this come from Sections 03, 04, 05, 07, and 64. The global bulk shift from Section 64 (the +0.32% mC / −0.39% hmC / flat total with those enormous Cohen's d values) is a clean opener because it immediately tells the reader "this is a TET-block signature, not global hypermethylation." Then you layer on the gene-body dominance from Section 03 (51.4% of gene bodies significant vs 8.3% of promoters), the dual volcanoes from Section 04 showing the directional asymmetry (70% hyper for mC, 83% hypo for hmC), the quadrant scatter from Section 05 showing 78.7% of co-significant genes are coordinated mC↑/hmC↓, and the effect-size distributions from Section 07.

What I think is missing here is a single, clean composite panel that shows these together with consistent styling, axis labels, and color coding. Right now each section has its own plot output folder with its own aesthetic. A new script that reads from `_shared_config.R` and builds a unified Figure 1 patchwork, with shared legends and a consistent color palette for "hyper/hypo" and "5mC/5hmC" across all subpanels, would be worth writing. The genome browser tracks from Section 46 (Syt1, Zbtb20, etc.) could anchor one subpanel as a representative locus, but you'd want to pick your single best example and show it at publication resolution rather than including all ten Gviz tracks.

### Figure 2: H2AK119ub as the upstream driver + chromatin geography

This is your "why did it happen" figure. It establishes the causal chain from K119ub accumulation to methylation change.

The spine here is the multi-mark logistic model from Section 33 (K119ub OR=4.707, dominant predictor, AUC 0.818) combined with the A/B compartment geography from Section 29 (mC-hyper enriched in A compartment, OR=13.64) and the decisive Polycomb falsification from Section 30 (Polycomb targets are excluded from hypermethylation, OR=0.063). The chromatin-state bar from Section 10 (Active_Promoter 93% hyper vs Repressed_Promoter 94% hypo) ties this together visually.

The challenge with the current section outputs is that Section 33 is quantitative DiffBind coefficients, Section 29 is compartment enrichments, and Section 30 is Polycomb exclusion tests — they're scattered across different analytical frames. A consolidated figure could present the logistic model forest plot (from Sec 33), a compartment-stratified direction bar (combining Sec 29/30), and the chromatin-state split (Sec 10) as three panels of a single argument: "K119ub drives hypermethylation, but only in accessible, active chromatin — Polycomb domains are protected."

One piece of new analysis that would strengthen this figure: the Section 17 "honest correction" (showing that 58.2% of DMR genes have no K119ub peak at all, and the conditional enrichment is +14.2 pp among peak-bearing genes) is important for intellectual honesty but gets buried as its own section. Folding that nuance into the Figure 2 caption or as an inset would preempt reviewer skepticism about the K119ub effect size.

### Figure 3: TET-impediment mechanism

This is probably your most mechanistically novel figure. The adjudication between TET-impediment and DNMT3A-recruitment models from Section 24 (AUC 0.793 vs 0.696, DeLong p=9.43e-49) is decisive and publishable on its own.

I'd build this around four panels: the baseline-5hmC predictor ROC from Section 23 (AUC 0.762 — genes with more 5hmC substrate show the largest changes), the head-to-head model comparison ROC from Section 24, the TET-KO comparison scatter from Section 26 showing the "dimmer not off switch" (3.3% absolute attenuation), and the demethylation ratio by chromatin state from Section 22 (showing the strongest impairment at Active_Promoter and Active_Enhancer, with a reversal at Repressed_Promoter).

What's important to note: the synthesis doc flags that K119ub is a positive predictor of hypermethylation as a differential DiffBind fold (Sec 14/33/41) but a negative predictor as baseline absolute signal in the TET model (Sec 24, β=−1.05). This is not contradictory — genes that *gain* K119ub hypermethylate, while genes with high *constitutive* K119ub are already compacted and inaccessible to DNMT3A — but the paper needs to explain this clearly, and a figure panel that explicitly visualizes both relationships side by side would help the reader see why they're compatible rather than conflicting.

### Figure 4: MeCP2 redistribution (the title figure)

This is where the paper's title lives, so it needs to be particularly strong. The challenge is that the P12-vs-adult developmental contrast (R1's headline) comes from the CUT&RUN/DiffBind arm, not the methylation pipeline. The synthesis doc is very clear that no P12 DUET methylation data exists.

For the methylation pipeline's contribution to this figure, I'd center on Section 75 (the redistribution reconciliation: 7,686 UP peaks but genome-wide gene-body signal drops, because MeCP2 concentrates on ~2,000 genes while thinning elsewhere) and Section 77 (the aging trajectory: mutant has 2.1× more aging-UP peaks, +22% fold at shared loci). The annotation split from Section 75 — UP peaks are 52% distal-intergenic/2.2% promoter, DOWN peaks are 62% intronic — is visually striking and mechanistically important.

A tension to resolve before writing: the paper outline says MeCP2 is "binding active promoters mostly," but Section 75 shows the opposite for MeCP2 *gain* (2.2% promoter). The *methylation* hypermethylation lands at active promoters (Sec 10), but MeCP2 gain lands distally. You need to decide whether R1's framing should be about where MeCP2 *accumulates* (distal-intergenic) or where the *methylation change* occurs (active promoters), because those are different answers.

A new analysis that would help: Section 77 ran fully now that Jai's files arrived, which is great. But what I don't see in the docs is a direct overlay of MeCP2 aging peaks with the methylation DMRs — do the mutant-unique aging MeCP2 peaks co-localize with hypermethylated gene bodies? That would tie the aging trajectory directly to the methylation phenotype and strengthen the R1→R2 connection. This could be a relatively straightforward Fisher's test or overlap analysis using the 1,654 mutant-unique aging genes against the coordinated gene set.

### Figure 5: MeCP2 reads chromatin, not methylation

This is the quantitative heart of R4, and Sections 62, 67, and 71 provide everything you need.

The regression from Section 62 (CG R²=0.017 vs chromatin R²=0.246; K119ub as top standardized β=+0.199) is the cleanest single result in the paper for this claim. Pair it with the Section 67 linchpin (359 genes where MeCP2 rises with no methylation change, K119ub gained 3.15× over genome, OR=3.15) and the Section 71 variance partition (K119ub explains 7.3% of MeCP2 fold, methylation ratio explains 0.0%).

The mirror-image profiles from Section 60 — MeCP2-Up loci show K119ub↑/5mC↑/5hmC↓/K27ac↓/K27me3↑, while MeCP2-Down loci mirror every mark except K119ub which stays flat — could be a visually compelling heatmap or parallel violin panel. The flat K119ub at MeCP2-Down is the "cleanest single discriminator" between the two mechanisms, per the synthesis, and it's worth highlighting.

What would be new and useful here: a single integrated scatter or ranked-gene plot that shows, for all ~20,000 genes, the K119ub fold on one axis and the MeCP2 fold on the other, with color encoding whether the gene has a significant methylation change. This exists as Section 59a, but in its current form it's one of seven quadrant scatters. A publication version would be more focused — one clean panel with the 359 "MeCP2-up-no-methylation" genes called out explicitly, perhaps with a marginal density showing their K119ub distribution vs background. That turns a statistical table into a visual argument.

### Figure 6: Neuronal specificity and Polycomb de-repression

This consolidates R5's findings. The key results are: neuronal genes are constitutively K119ub-rich (Sec 72, top-decile OR=1.70), the gene-set overlap shows MeCP2-gain tracks methylation (OR=5.16) far more than neuronal identity (OR=1.73) per Section 74, and synapse/axon genes show selective K27me3 loss (Sec 76, the Polycomb de-repression).

The most important framing decision here is about the word "exclusively" in the R5 title. The data don't support it — Section 74 explicitly shows MeCP2 redistribution tracks coordinated methylation five times more strongly than neuronal identity. I'd recommend softening to "preferentially" or "disproportionately" and building the figure around the refined mechanism: neuronal genes are the natural BAP1 substrate (constitutively K119ub-marked), so they're disproportionately affected, but the redistribution signal is fundamentally methylation-driven.

The Section 78 self-correction (narrow neuronal set shows spurious +0.012 total increase; unbiased broad set shows −0.0022 decrease) is methodologically important and shows intellectual rigor. I'd include it as a supplemental panel or an explicit callout, because reviewers will appreciate that you caught and corrected your own selection bias.

### Figure 7 (or Supplemental): Stoichiometry and mechanism model

The stoichiometry results from Sections 61 and 78 are mechanistically rich but complex. The genome-wide slope of −0.959 (distinguishable from −1, total methylation not conserved) favoring TET-inhibition plus continued DNMT3A, versus the neuronal-gene slope of −1.0 (stoichiometric, dehydroxymethylase-like) is a genuinely interesting mechanistic split. Whether this goes in a main figure or supplemental depends on how central the stoichiometry argument is to the paper's narrative. If the main story is "BAP1→K119ub→TET-block→methylation change→MeCP2 redistribution," the stoichiometry is supporting evidence for the TET-block step and might fit better in Figure 3 or as supplemental.

A summary model schematic (probably the final figure) should lay out the full cascade: BAP1 loss → H2AK119ub accumulation → (at active chromatin) TET access/efficiency reduced → 5mC accumulates, 5hmC depleted → MeCP2 reads K119ub + methylation shift → redistributes to distal sites at neuronal genes → developmental MeCP2 ramp amplified. This is a conceptual diagram, not a data figure, but it's essential for the reader.

---

## Priority Areas for New Work

Based on all of the above, here's what I think you should focus on, roughly in order of importance.

**Tier 1 — resolve before writing starts:**

The R1 framing tension needs resolution. Decide which panels come from the CUT&RUN arm versus the methylation pipeline and label them clearly. The "binding active promoters mostly" claim contradicts Section 75 and needs to be reconciled or dropped. This is a writing decision more than an analysis decision, but it shapes figure composition.

The R3 sex-difference PCA panel is unsupported by any verified table. Either run a dedicated sex-stratified PCA from the Section 42 framework and decide if it shows anything meaningful, or drop the panel and lean into the Methods statement that "none of the main findings varied by sex." The current plan to show sex differences while simultaneously claiming no sex variation is internally contradictory.

**Tier 2 — new consolidated figure scripts:**

Write figure-level scripts that pull from `_shared_config.R` and existing section outputs but produce unified multi-panel figures with consistent aesthetics, shared legends, and publication-quality formatting. One script per main figure, each reading from the section TSVs rather than recomputing from scratch. This is probably 5–7 new R scripts in a `scripts/figures/` directory.

The MeCP2 aging × methylation overlap I mentioned (do mutant-unique aging MeCP2 peaks co-localize with coordinated DMRs?) would be a short new analysis, maybe 100–150 lines, and would tie Section 77 to the rest of the methylation story.

**Tier 3 — analyses that would strengthen but aren't blocking:**

The FURTHER_ANALYSIS.md items the synthesis doc flagged as undone — methylation entropy (uniform vs stochastic mechanism), cell-type deconvolution (the Math1-Cre targets granule progenitors concern), and spatial chromosomal autocorrelation — are all things a reviewer might ask about. Having at least a preliminary answer for the cell-type question would be valuable, since the bulk signal is a cell-type-weighted average and the minority hypomethylated genes could plausibly be composition artifacts. Even a brief analysis showing the hypomethylated set doesn't correlate with known cell-type marker genes would help.

---

## What I'd Actually Do This Week

If you're starting writing next week, I'd focus the remaining prep time on three things. First, draft the figure plan on paper (or in a planning doc) with specific panel assignments, so when you sit down to write figure scripts you know exactly what each one needs to produce. Second, write the 2–3 highest-priority consolidated figure scripts (I'd start with Figure 1 and Figure 5, since R2 and R4 are your strongest results and benefit most from clean presentation). Third, resolve the R1 and R3 framing questions with Cole, because those are decisions that cascade into figure composition and you don't want to build panels that get cut.

The sections from ~30 onward that you mentioned getting more complicated — the non-CG exploration (51–58), the stoichiometry (61), the MeCP2 mechanistic suite (59–63), the neuronal characterization (72–78) — are actually in good shape analytically. The issue isn't that the analyses are incomplete; it's that they're presented as 40+ independent exploratory sections when the paper needs them distilled into maybe 3–4 figures. The consolidation work is the bottleneck, not additional computation.


---

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