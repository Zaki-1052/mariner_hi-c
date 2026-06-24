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