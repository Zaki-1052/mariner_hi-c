# Paper Results Synthesis: BAP1-KO Cerebellum Methylation Analysis

> **Header note.** Paper: *"Loss of BAP1 Deubiquitinase Restructures MeCP2 Binding and DNA Methylation Profiles During Cerebellar Neurodevelopment"* (Ramchandra, Alibhai, Wang, Song, Cui, Ferguson). Methylation analysis = **run-5** (deep-seq, **8 samples = 4 control + 4 mutant**, 2 batches × 1 per sex per condition; **sex covariate** in the modality GLM). All quantitative values below are sourced from the **verified `docs/results/section_*.md` group docs**, which recomputed every headline number directly from run-5 DMR BED files and exported TSV tables. Where a number is inherited unverified from a prose doc (FIGURES.md, TODO.md, CLAUDE.md), it is tagged `[UNVERIFIED]` and collected in the Caveats section. Mechanistic statements are tagged `[INTERPRETATION]`. **Do NOT cite FIGURES.md / TODO.md numbers — they are stale prior-run values (see Caveats §1).**

---

## Abstract of Findings

BAP1 is the catalytic deubiquitinase of the PR-DUB complex and the enzyme that removes H2AK119ub. In *Math1-Cre; Bap1^fl/fl^* mouse cerebellum, its loss produces a coherent, multi-layer epigenetic restructuring whose DNA-methylation arm this analysis characterizes across 78 sections.

The central methylation phenotype is a **coordinated, reciprocal gene-body shift: 5mC rises while 5hmC falls at the same loci**. Across 20,969 tested gene bodies, 10,775 (51.4%) carry a significant 5mC DMR (70% hypermethylated) and 11,484 (54.8%) a significant 5hmC DMR (83% hypomethylated); of the 8,371 genes significant in both, **6,589 (78.7%) show the diagnostic mC↑/hmC↓ pattern**. Because 5hmC is the obligate first oxidation product of 5mC by TET, simultaneous mC-gain/hmC-loss at thousands of gene bodies is the signature of a **TET-mediated active-demethylation block**. This is confirmed at the bulk level (control→mutant 5mC +0.32%, 5hmC −0.39%, total modC unchanged) and adjudicated mechanistically: baseline 5hmC (the TET substrate) is the single best predictor of hmC-DMR status (AUC 0.762, Sec 23) and the #1 standardized predictor of which genes hypermethylate (β=+1.25 in the Sec 24 TET-impediment model, AUC 0.793), the TET-impediment model decisively out-discriminates a DNMT3A-recruitment model (AUC 0.793 vs 0.696, DeLong p=9.43e-49), and BAP1-KO reproduces only ~3–9% of a published TET triple-KO's ratio shift in a graded — not binary — fashion (a "dimmer," not an "off switch").

The **upstream driver is H2AK119ub**. In a quantitative 4-mark DiffBind logistic model, K119ub gain is by far the strongest positive predictor of hypermethylation (OR≈4.7; 6-mark OR=10.29), an order of magnitude above any other mark, with active marks (ATAC, K27ac) protective. The methylation switch is spatially organized: hypermethylation lands on **active, A-compartment chromatin** (mC-hyper→A OR=13.64; Active_Promoter 72% hyper) while Polycomb/H3K27me3 domains are *protected* from de-novo methylation (hyper OR=0.063) and drift hypomethylated — exactly the DNMT3A-redistribution geography of TET-KO. At lost CTCF loop anchors, flanking dynamic-region hypermethylation parallels the IDH-glioma insulator-decay model (OR≈3.3), linking methylation to 3D-architecture loss.

The **MeCP2 arm** resolves the paper's title. MeCP2 redistribution tracks **chromatin state (K119ub), not methylation per se**: CG methylation explains <2% of MeCP2 binding variance while chromatin marks explain ~25%, with K119ub the dominant predictor. The linchpin is a set of **359 genes where MeCP2 binding rises with no methylation change** — there, H2AK119ub is gained 3.15× over genome (p=1.8e-24), proving MeCP2 follows the Polycomb/BAP1 axis independently of DNA methylation. At neuronal genes specifically, BAP1 loss relocates MeCP2 from gene bodies to distal-intergenic sites (UP peaks 52% distal-intergenic) and amplifies the normal developmental MeCP2 ramp (mutant: 2× more age-related UP peaks). Neuronal-identity genes are constitutively K119ub-rich (top-decile OR 1.70–2.34), making them the natural BAP1 substrate, and synapse/axon genes show selective Polycomb de-repression (extra K27me3 loss).

A key methodological caution threads through: the title's lead claim — *MeCP2 binding increases in mature (adult) but not young (P12) neurons* — is established by the **MeCP2 CUT&RUN/DiffBind arm**, which is upstream of the methylation visualization pipeline; **no P12 DUET methylation data exists**, so the methylation analysis is adult-only and cannot itself supply the developmental-trajectory contrast. The methylation pipeline's contribution to R1 is the aging-trajectory section (77) and the genomic-redistribution section (75).

---

# R1: BAP1 loss increases MeCP2 binding in mature neurons (not young)

> **Paper framing (Methylation Paper.md, Results §1):** *"Bap1 loss results in increased MeCP2 binding in mature neurons, with no changes in binding in young neurons."* Planned panels: MeCP2 volcano plots (P12 = no DMRs / no PCA clustering; Adult = thousands of DMRs / clear condition clustering); P12-vs-adult mutant; P12-vs-adult control; "seems to be binding active promoters mostly."
>
> **Primary sections:** 11, 31, 60, 62, 63, 65, 67–71, 75, 77. **Supporting:** 73, 74, 76.

## R1.1 Evidence for an age-dependent MeCP2 increase

**The core P12-vs-adult contrast is a CUT&RUN/DiffBind result, NOT a methylation-pipeline result.** The Methylation Paper.md R1 panels (MeCP2 volcano plots at P12 vs adult, condition PCA clustering) are produced by the upstream MeCP2 CUT&RUN DiffBind analysis. None of the 13 methylation group docs contains a P12-vs-adult MeCP2 DMR/volcano table. The methylation viz pipeline operates exclusively on **adult/late** DUET data (run-5), and **no P12 DUET methylation data exists** (TODO.md §12: *"Currently no P12 methylation data exists"*; the Hi-C arm has P12 and adult, but methylation does not). Consequently R1's headline age-dependence is **supported by collaborator data outside this synthesis's verified scope** — flag clearly in the manuscript which panels come from which assay.

What the methylation pipeline *does* contribute to R1 is the **developmental amplification** evidence (§R1.4) and the **redistribution mechanism** (§R1.3), both adult.

## R1.2 Quantitative MeCP2 gain (adult, run-5)

At the DMR/peak level the adult MeCP2 gain is large and well-characterized:

- **MeCP2 differential peaks:** 7,686 significant UP peaks vs 1,200 DOWN peaks (Sec 75; same totals in Sec 56, 61j, 63). UP outnumbers DOWN ~6.4:1.
- **Gene-level:** 2,570 genes have a significant MeCP2 change (12.3% of tested); 2,052 are MeCP2-Up (Sec 65). At stricter gene-body classification, 79 MeCP2-Up / 34 MeCP2-Down genes (Sec 59/60).
- **DMR-level coupling (Sec 11):** hypermethylated DMRs preferentially overlap MeCP2-Up peaks, Fisher **OR=5.13, p=1.27e-33** — the predicted "more 5mC → more MeCP2" enrichment. (FIGURES.md cites OR=5.60 for this panel — `[UNVERIFIED]`, prior run; verified value 5.13.)
- **Caveat — DMR-level vs gene-level divergence (Sec 11):** the OR=5.13 enrichment is driven by *peak co-location*; aggregated to gene level the directional coupling collapses (stratified Spearman rho≈0.02–0.08; binary-gain Fisher OR=1.04, p=0.795, NS; gene-level mC×MeCP2 2×2 OR≈0.92, p≈0.08). This is the motivation for the entire "MeCP2 reads chromatin, not methylation" thesis (R4).

## R1.3 Genomic redistribution (Section 75)

Section 75 resolves the apparent paradox that DiffBind reports ~7,686 MeCP2 UP peaks yet **gene-body MeCP2 BigWig signal drops genome-wide**:

- The 7,686 UP peaks concentrate on **~2,000 genes (≈10% of MeCP2-bound genes)**; the other ~90% lose signal slightly, so the genome-wide gene-body median falls even as a subset gains strongly. **MeCP2 redistributes, it does not uniformly increase.**
- **Annotation split (the mechanistic kicker):** UP peaks are **51.7% Distal Intergenic + 41.0% Intron, only 2.2% Promoter**; DOWN peaks are 61.8% Intron, 19.5% Distal Intergenic, 8.0% Promoter (Sec 75). [INTERPRETATION] MeCP2 is recruited to new distal regulatory sites while vacating gene bodies.
- **Functional footprint contrast:** K119ub GSEA = 115 significant GO BP terms (broad developmental/Polycomb); MeCP2 GSEA = **1** significant term ("synapse assembly") (Sec 75, cross-checked vs 61k). [INTERPRETATION] MeCP2 is the selective, neuronal-specific reader; K119ub blankets many programs.

> **Note on paper claim "binding active promoters mostly" (R1 panel):** Section 75 shows MeCP2 UP peaks are predominantly **distal-intergenic/intronic, NOT promoter** (2.2% promoter). The "active promoters mostly" framing in the outline is in tension with the redistribution result — the *methylation* hypermethylation lands at active promoters (R2), but *MeCP2 gain* lands distally. Reconcile in manuscript.

## R1.4 Developmental trajectory (Section 77)

Section 77 uses **young-vs-adult MeCP2 CUT&RUN DiffBind** (files from collaborator Jai, now present — see Caveats) to separate normal developmental MeCP2 gain from mutant excess:

- **Aging peaks per genotype:** Control 10,930 UP / 2,822 DOWN; Mutant **23,117 UP / 10,646 DOWN** — mutant has **2.1× more aging-UP** and 3.8× more aging-DOWN peaks.
- **Overlap:** 7,305 control aging-UP peaks overlap mutant (66.8% of control preserved); mut-unique = 15,812 peaks. Gene-level: 2,620 shared, **1,654 mutant-unique aging genes**, 288 ctrl-unique.
- **Mutant-unique aging genes** are GO-enriched for 435 BP terms, **49 neuronal**.
- **At 7,305 shared aging loci:** median ctrl fold = 1.829 vs mutant 2.241 (mutant +0.41 log2 higher, ~+22% more fold); 2,266/7,305 shared peaks at neuronal genes.

[INTERPRETATION] MeCP2 normally accumulates as neurons mature; BAP1-KO **amplifies** this developmental ramp (it does not rewire it wholesale — two-thirds of the control aging program is preserved). A plausible driver: excess H2AK119ub creates aberrant chromatin that attracts additional MeCP2 during the maturation window, converting a developmental ramp into an overshoot. This is the methylation-pipeline's strongest direct support for the title's "mature neurons" age-dependence.

## R1.5 Supporting MeCP2-recruitment mechanism (Sec 60, 62, 63, 67–71)

- **Sec 60 (mirror-image profiles):** MeCP2-Up loci show the full cascade (K119ub BigWig +0.63 q=1.8e-19; 5mC +0.049; 5hmC −0.018; K27ac −0.72; K27me3 +1.11), and MeCP2-Down loci mirror every mark **except K119ub, which stays flat** (−0.039, p=0.182, NS) — the cleanest single discriminator that MeCP2 loss is a separate redistribution, not the reverse of K119ub gain.
- **Sec 62 (regression verdict):** chromatin marks explain MeCP2 binding ~15× better than CG methylation (binding R²: chromatin 0.246 vs CG 0.017); K119ub is the top standardized predictor (β=+0.199). (Full detail in R4.)
- **Sec 67 (linchpin):** at 359 MeCP2-up-no-methylation genes, K119ub gained 3.15× over genome. (Full detail in R4/R5.)
- **Sec 71 (dose model rejected):** the 5mC/5hmC skew ratio has ~zero predictive power for MeCP2 (R²=0.001); K119ub alone explains 7.3%.
- **Sec 63 (heterogeneity):** MeCP2 sites split into 4 chromatin-defined clusters (not one Polycomb block); modal state "Unmarked" across clusters.

### R1 — supported vs missing panels

| Planned R1 panel (Methylation Paper.md) | Status |
|---|---|
| MeCP2 volcano P12 (no DMRs) vs Adult (thousands) | **MISSING from methylation docs** — CUT&RUN/DiffBind arm; no P12 DUET data |
| P12-vs-adult mutant; P12-vs-adult control | **MISSING** — CUT&RUN arm |
| PCA clustering by condition (adult) | **MISSING** — CUT&RUN arm (Sec 42 has DUET *methylation* PCA, not MeCP2) |
| "Binding active promoters mostly" | **In tension** — Sec 75 shows UP peaks distal-intergenic, not promoter |
| Adult MeCP2 gain magnitude | **Supported** — Sec 75 (7,686 UP), Sec 11 (OR=5.13) |
| Developmental amplification | **Supported** — Sec 77 (2.1× aging-UP peaks) |

---

# R2: H2AK119ub accumulation → global 5mC increase / 5hmC decrease at gene bodies

> **Paper framing (Results §2):** *"Increase in H2AK119Ub after Bap1 loss leads to a global 5mC increase and 5hmC decrease trend at gene bodies."* Planned panels: volcano with top genes; 5mC-vs-5hmC quadrant (concordant/discordant); asymmetric-direction bar; effect-size plot; methylation-vs-chromatin-state bar (active promoters gain, repressed lose); heatmap showing MeCP2-down specific to Polycomb targets.
>
> **Primary:** 03–08, 14, 17, 18, 22–26, 64. **Supporting:** 09, 12, 13, 15, 16, 29, 30, 38–41, 61, 78.

## R2.1 DMR statistics & directional bias (Sec 03, 04, 06, 64)

**Canonical run-5 headline numbers** (Sec 03/09, recomputed from `DMR_*_20260402_191818.bed`, matching `summary_statistics.txt`):

- **20,969 gene bodies tested.**
- **5mC significant = 10,775 (51.4%):** 7,513 hypermethylated (69.7%) / 3,262 hypomethylated (30.3%).
- **5hmC significant = 11,484 (54.8%):** 1,963 increased (17.1%) / 9,521 decreased (82.9%).
- **Co-significant (both) = 8,371;** 5mC-only = 2,404; 5hmC-only = 3,113 (Venn, Sec 06).
- **Region localization (Sec 03):** gene bodies 51.4% sig >> CpG shores 30.2%, CpG shelves 23.8%, promoters 8.3%, CpG islands 5.0%, TSS 1.4%. The signal is overwhelmingly a **gene-body** phenomenon.

**Bulk genome-wide confirmation (Sec 64):** control→mutant mean 5mC **+0.317%** (paired t p=0.041, Cohen's d=2.79); 5hmC **−0.395%** (p=0.018, d=−4.54); total modC −0.079% (**NS**, p=0.492). Every individual mutant sample has higher 5mC and lower 5hmC than its matched control. [INTERPRETATION] Directionally-opposed, near-equal 5mC↑/5hmC↓ with flat total modC = classic TET-block signature; tiny absolute deltas but huge Cohen's d = exceptionally reproducible systemic effect, NOT global hypermethylation.

**Top genes (Sec 04, 06, 42):** the q-floor (q≈4.82e-306) is dominated by neuronal/synaptic genes — Syt1 (+17.3% mC / −15.0% hmC), Cntnap2, Zbtb20 (+7.8%/−5.6%), Trpm3, Mcu, Dlgap1. Top-20 5mC DMRs are 19/20 hypermethylated; top-20 5hmC are 16/20 hypomethylated.

## R2.2 The coordinated quadrant & effect sizes (Sec 05, 07, 20, 21, 28)

- **Coordinated mC↑/hmC↓ = 6,589 of 8,371 co-significant (78.7%)** (Sec 05/09). Reverse-coordinated (mC↓/hmC↑, "discordant" Q2) = 1,255 (15.0%); same-direction = 527 (6.3%).
- **Four-quadrant breakdown (Sec 20/21/28):** Q1 (mC↓/hmC↓)=411; **Q2 (discordant, mC↓/hmC↑)=1,255**; Q3 (mC↑/hmC↑)=116; **Q4 (coordinated, mC↑/hmC↓)=6,589 (78.7%)**.
- **Effect sizes (Sec 07):** among significant genes, net mean 5mC = **+1.72%** (hyper-only +3.45%, n=7,513); net 5hmC = **−1.66%** (hypo-only −2.29%, n=9,521). Per-gene effects are modest (~2–3%) but mirror-image and consistent across thousands of genes.
- **Q4 vs Q2 are near-mirror regulatory programs (Sec 21):** Q4 = repressed, K119ub-gaining, MeCP2-up, Active-promoter; Q2 = upregulated, K119ub-losing, K27ac-gaining, Repressed/Bivalent-promoter, ~2.6× more loop-anchored. Expression Up/Down Q2-vs-Q4 Fisher OR≈95 (p≈1.2e-104); MeCP2 Up/Down Q2-vs-Q4 OR≈0.23; loop OR≈2.84.

## R2.3 GO / KEGG enrichment (Sec 08)

Of 7,513 hypermethylated genes submitted:
- **GO BP #1 = RNA splicing** (GO:0008380, 274 genes, FoldEnrichment 2.42, p.adjust 6.27e-54).
- **GO CC #1 = nuclear speck** (256 genes, p.adjust 2.55e-49).
- **GO MF #1 = ubiquitin-like protein transferase activity** (269 genes, p.adjust 1.36e-52).
- **KEGG #1 = Autophagy-animal** (113 genes, p.adjust 8.82e-24); then Amyotrophic lateral sclerosis, Endocytosis, Ubiquitin-mediated proteolysis, multiple neurodegeneration pathways.
- **Delta-ratio decile split:** D10 (most TET-impaired) → developmental terms (pattern specification, regionalization, cell-fate commitment); D1 (least impaired) → metabolic (oxidative phosphorylation, respiration).

[INTERPRETATION] The ubiquitin axis (MF#1 ubiquitin-transferase, KEGG ubiquitin-mediated proteolysis) thematically loops back to BAP1's deubiquitinase biology; the decile dichotomy shows the **neuronal differentiation program** is the primary demethylation-block casualty.

## R2.4 Demethylation-ratio / TET mechanism (Sec 22–26, 61, 78)

The 5hmC/(5mC+5hmC) **TET conversion-efficiency ratio** is the group's unifying continuous metric:

- **Sec 22:** WT median ratio 0.1275 → KO 0.1193 (delta −0.0087); **71.4% of genes show a decreased ratio** (impaired TET). Strongest at Active_Promoter (Δ−0.028 ***) / Active_Enhancer (Δ−0.025 ***); *reversed* at Repressed_Promoter (Δ+0.005). Partial magnitude (0.128→0.118, not →0) = indirect, graded defect.
- **Sec 23 (substrate model):** baseline 5hmC predicts hmC-DMR status with **AUC 0.762** vs K119ub AUC 0.573; combined 0.790. Genes with more 5hmC substrate show the largest hmC changes.
- **Sec 24 (mechanism adjudication):** **TET-impediment model AUC 0.793 >> DNMT3A-recruitment AUC 0.696** (DeLong **p=9.43e-49**; exclusive-feature p=6.13e-75). Baseline 5hmC is the #1 predictor (standardized β=+1.25); **K119ub is a *negative* predictor (β=−1.05)** — opposite to direct DNMT3A-UDR recruitment. Robust across non-promoter stratum and all 5hmC tertiles; K119ub×5hmC interaction is significantly *negative*.
- **Sec 25 (continuous response):** delta-ratio linear refits reproduce the same feature hierarchy with exact sign-flips; OLS robust under quantile regression.
- **Sec 26 (TET triple-KO comparison, GSE166423):** BAP1-KO reproduces only **~3.3% (absolute) to ~8.7% (baseline-normalized)** of the TET-KO ratio shift; TET-KO is binary (68.5% complete 5hmC loss), BAP1-KO is graded (54% weak, 45% moderate). Apparent per-gene rho 0.22 collapses to 0.10 after residualizing on baseline (99.9% baseline-driven). [INTERPRETATION] A **dimmer switch, not an off switch** — BAP1 loss impairs TET *access/efficiency* (via PRC1/H2AK119ub), it does not eliminate the enzyme.
- **Sec 61 (stoichiometry):** genome-wide Δ5mC~Δ5hmC slope = **−0.959** (distinguishable from −1), and total methylation **rises** at MeCP2-Up (+0.032) and coordinated (+0.007) loci rather than being conserved → favors **TET inhibition + continued de-novo DNMT3A** over pure dehydroxymethylase. K27me3→ATAC compaction mediates only ~18% of the TET-efficiency drop.
- **Sec 78 (unbiased self-correction):** with an unbiased GO-neuronal set, neuronal genes' total methylation *decreases* (−0.0022) and their slope sits at **−1.0** (stoichiometric, dehydroxymethylase-like), whereas the narrow DMR-selected set spuriously showed +0.012. (See Caveats — this corrects Sec 61's narrow-set bias.)

## R2.5 K119ub as the upstream driver (Sec 14, 17, 18, 33, 38–41)

- **Sec 14 (peak overlap):** differential K119ub Fisher **OR=4.46, p=5.80e-97** — hyper DMRs gain K119ub (19.6% overlap K119ub-Gained vs 10.4% Lost), hypo DMRs lose it. Condition-overlap OR=0.700 reflects hypo DMRs' high baseline K119ub.
- **Sec 17 (honest correction):** including the no-peak majority (58.2% of DMR genes have no K119ub), conditional gain among peak-bearing genes is mC-Up 47.8% vs background 33.6% (**+14.2 pp**); a real but partial, minority effect.
- **Sec 18 (continuous BigWig):** confirms a genuine **global K119ub increase** (background median log2FC +0.0070, p=1.8e-20); mC-Up genes shift positive, mC-Down negative — but gene-specific effect sizes are tiny (Cliff's delta 0.10–0.18, "negligible-to-small"). The mechanism is a diffuse global rise, not sharp locus-specific redistribution.
- **Sec 33 (quantitative DiffBind, the spine):** 4-mark logistic model — **H2AK119ub OR=4.707** (top) vs H3K27me3 1.443, H3K27ac 0.480, ATAC 0.221; 4-mark AUC=0.818 (extended +baseline 0.869). Cross-mark Spearman: K119ub fold vs mC diff rho=+0.410 (strongest positive correlate). 178 genes have ≥3 marks converging in cascade direction (72% mC-hyper).
- **Sec 41 (6-mark model):** **H2AK119ub OR=10.29** dominant; H3K36me3 OR=1.00 (NS, rules out DNMT3B/SETD2); pathway attribution DNMT3A-axis 1,630 vs DNMT3B-axis 84 genes (~19:1) but 73% of hyper genes have *no* recruiting-mark evidence (→ upstream TET-block, not active DNMT recruitment).

> **Reconciling the K119ub sign (important for the manuscript).** K119ub is the **strongest predictor of hypermethylation as a differential DiffBind fold** (Sec 14/33/41: OR 4.5–10.3) but a **negative predictor as a baseline absolute signal in the TET model** (Sec 24: β=−1.05). [INTERPRETATION] These are not contradictory: genes that *gain* K119ub in the mutant hypermethylate (the cascade), but genes with high *constitutive* K119ub (Polycomb, already silenced/inaccessible) are protected from de-novo methylation. The DNMT3A-UDR direct-recruitment model (Gretarsson/Chen 2024) predicts a positive baseline-K119ub effect; its absence (Sec 24) is the evidence against direct recruitment and for the indirect TET-block (Sec 26 magnitude).

## R2.6 A/B-compartment & Polycomb context (Sec 29, 30)

- **Sec 29:** mC-hyper enriched in **A compartment OR=13.64** (48.8% of A genes hyper vs 6.5% of B); hmC-loss → A OR=9.77; reciprocal mC-hypo→B OR=1.73, hmC-gain→B OR=3.03. B→A-shifted bins are *depleted* for hyper (OR=0.063) — de-novo methylation targets stable-A euchromatin, not freshly-opened regions. [INTERPRETATION] cleanest cross-modality confirmation of TET-KO DNMT3A-redistribution geography.
- **Sec 30 (the decisive Polycomb falsification):** Polycomb/H3K27me3 targets are **excluded from hypermethylation (OR=0.063, q=0)** and enriched in hypomethylation (OR=9.80); robust across 3 Polycomb definitions. Per-state hyper rates: Active_Promoter 71.8%, Active_Enhancer 65.3% (the targets) vs Polycomb 1.8%, Repressed_Promoter 3.0% (protected). [INTERPRETATION] Falsifies the naive "K119ub marks Polycomb → Polycomb hypermethylates"; compact heterochromatin is inaccessible to DNMT3A.

### R2 — supported vs missing panels

| Planned R2 panel | Status |
|---|---|
| Volcano with top genes | **Supported** — Sec 04 (10,775/11,484 sig; Syt1, Cntnap2, etc.) |
| 5mC-vs-5hmC quadrant (concordant/discordant) | **Supported** — Sec 05a, 21b (Q4=6,589 / Q2=1,255) |
| Asymmetric-direction bar | **Supported** — Sec 03b (70% hyper / 83% hmC-down) |
| Effect-size plot | **Supported** — Sec 07 (+1.72% / −1.66%) |
| Methylation-vs-chromatin-state bar (active gain, repressed lose) | **Supported** — Sec 10 (Active_Promoter 93% hyper), Sec 30 |
| Heatmap: MeCP2-down specific to Polycomb targets | **Partial** — Sec 60/67 show MeCP2-down mirror profile & Polycomb context, but no single "MeCP2-down = Polycomb" heatmap table; Sec 31 ties MeCP2-down to gained (hypomethylating) loops. Verify the exact intended panel exists |

---

# R3: Most DMRs at gene bodies, CpG shores/shelves; sex differences

> **Paper framing (Results §3):** *"Most 5mC and 5hmC DMRs occur in gene bodies, CpG Shores, and CpG shelves."* Planned panels: DMR heatmap per context; sex differences for both 5mC and 5hmC regardless of genotype (PCA).
>
> **Primary:** 01, 02, 03, 43, 44, 64. **Supporting:** 09, 32, 42.

## R3.1 Region breakdown (Sec 03)

Significant DMR counts by region class (run-5, Sec 03):

| Region | Significant mC DMRs | % of class |
|---|---|---|
| **Gene bodies** | 10,775 / 20,969 | **51.4%** |
| CpG shores | 9,842 / 32,581 | 30.2% |
| CpG shelves | 6,924 / 29,094 | 23.8% |
| Promoters | 1,692 / 20,417 | 8.3% |
| CpG islands | 442 / 8,910 | 5.0% |
| TSS regions | 192 / 14,165 | 1.4% |

[INTERPRETATION] Gene bodies and CpG shores/shelves carry the signal; promoters/islands/TSS are largely spared — the hallmark of dysregulated DNMT3/TET cycling over transcribed bodies and *dynamic* CpG flanks, not promoter switching.

> **Caveat (Sec 03 flag #5):** non-gene-body region counts are *non-deduplicated* (per-feature, multiple per gene), while gene bodies are per-gene-deduped. The region bar chart mixes the two count types; qualitative "gene bodies dominate" holds but exact cross-region ratios should be stated carefully.

**CpG-island specificity (Sec 48):** of 8,910 islands, only 442 are mC-significant (122 hyper / 320 hypo). Hyper islands are overwhelmingly **amplification of pre-existing methylation** (baseline mC 0.566; only ~15% de-novo even at the loosest threshold), and K119ub is **depleted** (not enriched) at hyper islands (17.2% vs 30.2% background). [INTERPRETATION] The K119ub→methylation coupling operates at gene bodies/dynamic flanks, NOT at islands; island hyper is secondary amplification.

## R3.2 CpG-context specificity & non-CG (Sec 03, 32)

- **CG dominates:** baseline modC — CG 72.22%, CHG 0.628%, CHH 0.862% (~100× lower). **0 significant non-CG DMRs at the global summary level.**
- **CHG exploratory (Sec 32):** the modality GLM does call **320 significant gene-body CHG mC DMRs** (157 hyper / 163 hypo), but they are ~100× smaller in effect than CG (permille vs percent), concentrated on chr7+chr8 (76.3%), and **61% direction-discordant** with CG mC where they overlap. [INTERPRETATION] CHG is low-amplitude, mostly CG-independent background, not a corroborating layer. (Reconcile with FIGURES.md's "No Signal" framing — see Caveats; recommend "statistically significant but biologically negligible and CG-discordant.")

## R3.3 Sex differences (Sec 02, 43, 44; QC Sec 01)

This is the weakest-supported R3 subsection in the verified docs. What exists:

- **Sec 02 (correlation):** 5mC within-control mean r=0.886, within-mutant 0.897, between-group 0.884; 5hmC within-control 0.685, within-mutant 0.692, between-group 0.668. The structure shows a detectable but modest genotype effect; **per-sex clustering is not separately tabulated in the verified group docs.**
- **Sec 43 (chrX-removal panel):** the `43m_chrX_direction_comparison` panel exists *specifically to demonstrate the global direction asymmetry is NOT an artifact of X-linked DMRs in a mixed-sex cohort* — i.e., the methylation phenotype is autosomal and robust to sex-chromosome contribution. Per-chromosome Fisher (Sec 09 tables) shows **chrX depletion** of DMRs (mC obs/exp 0.564, p=2.31e-42; hmC 0.712, p=1.98e-21). **Note: Sec 43 quantitative counts were NOT independently verified against the TSV** (`[UNVERIFIED]` per Sec 42–50 doc flag; only panel inventory confirmed).
- **Sec 44 (allele-specific methylation):** BAP1-KO ~doubles ASM sites (control mean 11,295/sample → mutant 22,080, ratio 1.95×); 58,955 unique significant mC ASM loci; 1,408 gene-body DMRs harbor ASM sites. This is allele- not sex-resolved, but bears on the "asymmetry" theme.

> **The paper's planned R3 sex-difference PCA panel ("Sex differences observed for both 5mC and 5hmC, regardless of mutant or control status") is NOT quantified in any verified group doc.** Section 42 produces per-sample PCA of methylation (max-significance + all-gene), colored by condition/batch, which *could* reveal sex structure, but no group doc reports a sex-axis PCA or sex-stratified DMR statistic. **This panel appears unsupported by a verified table and likely needs a dedicated sex-PCA analysis or must be drawn from Sec 42's PCA visuals (visual-only).** Also note the paper's own Methods states *"None of the main findings varied by sex"* — partially in tension with a headline sex-difference panel; reconcile.

### R3 — supported vs missing panels

| Planned R3 panel | Status |
|---|---|
| DMR heatmap per context | **Supported** — Sec 03 region breakdown (gene bodies 51% > shores 30% > shelves 24%); a per-context "heatmap" specifically may need confirming |
| Sex differences PCA (5mC & 5hmC, both genotypes) | **MISSING / unverified** — no sex-stratified table in group docs; Sec 42 PCA is condition/batch-colored; tension with Methods "no sex variation" |

---

# R4: MeCP2 acts as a reader of chromatin-state shift, not solely methylation

> **Paper framing (Results §4):** *"In neurons, MeCP2 acts as a reader of chromatin state shift, rather than solely changes in methylation."* Planned panels: MeCP2 peak overlap at DMRs bar; #sig DMRs in modCs + genomic elements; #sig 5mC vs 5hmC + elements; 5mC-up/5hmC-down at same loci (bigwig example loci); Ub-vs-MeCP2 log2 plot.
>
> **Primary:** 10, 11, 19, 33, 55–58, 59, 60, 62, 63, 71. **Supporting:** 13, 15, 48–54.

## R4.1 Chromatin state determines the MeCP2 response (Sec 10, 62, 68, 69)

**Chromatin state of DMRs (Sec 10):** of 10,775 significant mC DMRs, Active_Promoter = 4,906 (45.5% of all; 60.7% of hyper) and is **93.0% hypermethylated**; Repressed_Promoter = 1,718 and is **94.4% hypomethylated**; Unmarked = 3,952 (36.7%). The mirror-image partitioning (hyper→active, hypo→repressed) means a naive pooled enrichment would cancel — the directional split *is* the signal.

> **Caveat:** TODO.md's "49.9% DMRs at Active_Promoter" is `[UNVERIFIED]`; the table gives 45.5% (all) / 60.7% (hyper).

**The decisive regression (Sec 62 — "Jai's hypothesis" test):**
- **Binding-level R²:** CG-only **0.017** vs Chromatin-only **0.246** vs Full 0.260 (n=202,574 peaks). Chromatin explains MeCP2 binding ~15× better than CG methylation.
- **Differential R²:** CG-only 0.013 vs Chromatin-only 0.170 vs Full 0.173 — MeCP2 *redistribution* tracks chromatin remodeling, not methylation.
- **Variance partition (binding):** Chromatin-unique 24.3%, CG-unique 1.46%, Shared 0.22%, Unexplained 74.0%.
- **Top standardized β (binding):** **K119ub +0.199** (dominant), ATAC +0.114, CG-5mC +0.089, K27me3 +0.061; CG-5hmC −0.018. LASSO retains K119ub first.
- **R² gain (Full−CG) largest at Active_Enhancer (+0.376) and Poised_Enhancer (+0.210)**; smallest at Repressed_Promoter (+0.008), Polycomb (+0.027). The chromatin advantage concentrates where MeCP2 is being redistributed.

**Context-gated coupling (Sec 68, 69):**
- modC-MeCP2 Spearman by context (Sec 68): K27me3+K27ac Bivalent rho=**0.461** (strongest) > K27ac-only 0.183 > K27me3-only 0.189 > Neither 0.070 (barely couples).
- MeCP2 distinguishes 5mC from 5hmC (Sec 69): in Active and Bivalent contexts 5mC rho positive (0.16–0.31) vs 5hmC weak/negative; interaction significant — **except in facultative heterochromatin (K27me3-only), where the interaction is NS** (MeCP2 does not discriminate), pointing to non-methylation recruitment there.

## R4.2 Non-CG methylation's role (Sec 51–58)

The full non-CG (mCA/mCH reader) hypothesis funnel — does evoC CG data miss a non-CG MeCP2 axis?

- **Sec 51/53/54 (negative at every scale):** evoC non-CG is at the detection floor (gene-body ~1e-5; only 0.51% CHG / 1.64% CHH of MeCP2 peaks detectable; per-TAD medians = 0; **0 significant CHG DMRs, 2 CHH**). MeCP2 peaks are modestly non-CG-favored above shuffled control (CHG Fisher OR=1.73; CHH OR=1.60) but absolute levels are negligible.
- **Sec 55/56 (MeCP2 as indirect proxy):** CG explains only ~1–2% of MeCP2 distribution; the CG-corrected residual is **not** TAD-organized (variance ratio ~0.066 vs CG-5mC 1.018) and is governed by *low-CG* context, not a structured non-CG field. 35.5% of MeCP2-Up peaks are "non-CG candidates" (gain MeCP2 with no CG increase; 2,726 peaks), enriched in Unmarked chromatin (77%).
- **Sec 57 (external Ecker WGBS validation — the one positive):** non-CG-candidate peaks sit at genuinely higher Ecker CH (median 0.0159 vs 0.0129, p=2.6e-60) *and* lower Ecker CG (p=1.4e-55) — a double dissociation matching the mCA-reader model categorically.
- **Sec 58 (dose-response rejected):** higher CH does **not** scale CG-unexplained MeCP2 (Jonckheere increasing p≈1.0; decreasing p=2e-09; same in MeCP2-Down control). [INTERPRETATION] non-CG association is real *as a category* but NOT a quantitative causal mCA-dose mechanism; this evoC+Ecker analysis cannot support a non-CG-dose model for BAP1-KO MeCP2 redistribution.

[INTERPRETATION] Throughout 51–58, **CG-5hmC repeatedly emerges as the most TAD-organized, most MeCP2-predictive track** (Sec 54 variance ratio 2.08 vs CG-5mC 1.02; Sec 55 5hmC regression coefficient steeper than 5mC) — quietly redirecting mechanistic weight from non-CG back to the TET/5hmC axis.

## R4.3 Multifeature regression (Sec 62)

Covered in R4.1 — Sec 62 is the quantitative heart of R4: **CG methylation explains <2% of MeCP2 binding; chromatin marks (K119ub foremost) explain ~25%.** This is the direct evidence that "MeCP2 reads chromatin-state shift, not solely methylation."

## R4.4 K119ub at unmethylated MeCP2 loci (Sec 67) — and the Ub-vs-MeCP2 plot (Sec 11, 59)

- **Sec 67 (KEY RESULT, also R5):** at **359 genes with significant MeCP2 gain but NO significant methylation change** (the Sec 65 `mecp2_only` cell), K119ub is gained: our-group median log2FC +0.202 vs background −0.036 (Mann-Whitney p=5.50e-32); **72.8% gain K119ub vs 45.9% background, Fisher OR=3.15, p=1.82e-24**; MeCP2 fold vs K119ub rho=+0.241 (p=4.4e-06). Co-incident K27me3 gain + K27ac loss. [INTERPRETATION] At loci where methylation is flat, MeCP2 still rises and the dominant concurrent change is **K119ub gain** — MeCP2 follows the BAP1→K119ub Polycomb axis independently of DNA methylation.
- **The "Ub-vs-MeCP2 log2 plot" (paper R4 panel):** delivered by **Sec 59 (59a: K119ub log2FC vs MeCP2 fold)** and **Sec 71** (K119ub vs MeCP2 rho=0.187, p=3.5e-96; ratio rho=0.014 NS; K119ub variance-unique 7.3% vs ratio 0.0%). Sec 11 (11f/11g) adds the delta-ratio LM/GLM (delta_ratio coefficient −0.367, p=1.4e-06).
- **Sec 65 (scale):** methylation change is ~5× more pervasive than MeCP2 change (51.4% vs 9.8% of genes); only 12.6% of hyper genes recruit more MeCP2 (Fisher OR=2.63, p=2.2e-61 but dominant mc_only cell of 8,620 genes). Methylation is necessary but not sufficient for a MeCP2 response.

### R4 — supported vs missing panels

| Planned R4 panel | Status |
|---|---|
| MeCP2 peak overlap at DMRs bar | **Supported** — Sec 11a (OR=5.13) |
| #sig DMRs in modCs + genomic elements | **Supported** — Sec 03 region breakdown |
| #sig 5mC vs 5hmC + elements | **Supported** — Sec 03, 06 (10,775 / 11,484) |
| 5mC-up/5hmC-down at same loci | **Supported (stat)** — Sec 05 (6,589 coordinated) |
| **Bigwig example loci** | **VISUAL-ONLY (Sec 46)** — Gviz tracks at 10 KEY_GENES (Syt1, Zbtb20, Trpm3, Epha3, Mcu, Cntnap2, Lpp, Dlgap1, Arhgap26, Cdh8); **NO quantification table** — flag as figure-only |
| Ub-vs-MeCP2 log2 plot | **Supported** — Sec 59a, 71, 11f/g |

---

# R5: At neuronal genes, BAP1 loss relocates MeCP2 from euchromatic to heterochromatic loci

> **Paper framing (Results §5):** *"Exclusively at genes involved in neuronal structure and development, Bap1 loss leads to relocation of MeCP2 from euchromatic to heterochromatic loci."* Planned panel: "the figure with a lot of log2 plots."
>
> **Primary:** 72–78, 61h–61k. **Supporting:** 27–31, 65–67.

## R5.1 Neuronal gene-set definition & constitutive K119ub enrichment (Sec 72)

Two neuronal gene-set definitions are used (carefully distinguished in Sec 78):
- **Broad (unbiased):** 5,614 genes from org.Mm.eg.db GO BP `synap|neuron|axon|dendrit|nervous` (Sec 72). 4,118–5,077 have valid data.
- **Narrow (biased/circular):** 1,149 genes from GO enrichment of significant DMRs (Sec 61) — DMR-selected, must not be used for unbiased estimates.

**Sec 72 (constitutive K119ub):** neuronal genes are **constitutively H2AK119ub-enriched in wildtype** cerebellum — top-quartile OR 1.378 (p=7.5e-19), **top-decile OR 1.701** (p=3.4e-26); dose-response across deciles (Spearman rho=0.648, p=0.049). GO of K119ub-high genes = Polycomb developmental program (pattern specification, cell-fate commitment, regionalization) layered with neuronal terms. The *gained* (mutant-only-high) set yields **0** significant GO terms → BAP1-KO amplifies signal at already-marked loci, it doesn't create a new target class. [INTERPRETATION] neuronal-identity genes sit in constitutively PRC1/K119ub-marked chromatin — the natural BAP1 substrate.

## R5.2 Triple overlap & "MeCP2 tracks methylation, not lineage" (Sec 74, 76)

- **Sec 74 (gene-set overlap):** Coordinated × MeCP2-Up **OR=5.16** (p=2.8e-12) >> Neuronal × MeCP2-Up **OR=1.73** (p=0.034); Neuronal × Coordinated OR=1.05 (NS). **Triple overlap (Neuronal∩Coordinated∩MeCP2-Up) = 16 genes.** [INTERPRETATION] MeCP2 redistribution tracks the *methylation* redistribution far more than generic neuronal identity.
- **The 16 triple-overlap genes (Sec 76):** Ap3b1, Astn2, Cntn6, Epyc, Fgf1, Fut9, Gprin3, Hcn1, Hif1a, Il1rapl1, Lgi1, Micu3, Ntn4, Prom1, Snca, Tspan7. Too few for robust per-mark stats and NOT K119ub-extreme (0/16 in top decile).

## R5.3 Synapse/axon specificity (Sec 76)

- **Synapse/axon vs broader-neuronal, K27me3:** Δmedian **−0.044** (Wilcoxon p=2.95e-3, adj 6.65e-3) — synapse genes lose *more* K27me3. ATAC p=0.654 (NS), K27ac p=0.757 (NS) — NOT special for accessibility/enhancer activation. Broader-neuronal vs non-neuronal K27me3 NS (p=0.289) → the specificity resides in the synapse/axon subset.
- **K119ub top-decile:** Neuronal OR=2.34 (p=1.3e-54), Synapse/axon OR=1.68 (p=2.1e-13).

[INTERPRETATION] Synapse/axon genes are special specifically for **selective Polycomb de-repression** (K27me3 erosion), not active enhancer activation — fitting BAP1-loss disorganizing Polycomb domains. (Sec 73 corroborates with a positive K119ub×neuronal interaction for ATAC p=0.017 and K27ac p=4.6e-07, and shows K119ub-high neuronal genes are disproportionately Repressed_Promoter OR=1.74 / Bivalent_Promoter OR=2.87.)

## R5.4 Stoichiometry at neuronal loci (Sec 78) — and the neuronal log2-plot figure (Sec 61h–k)

- **Sec 78 (self-correction, the cleanest neuronal result):** with the unbiased broad set, neuronal genes' **total methylation decreases** (mean Δ −0.0022) — opposite to the narrow set's spurious +0.012. Neuronal and synapse/axon **slopes sit at −1.0** (broad neuronal −0.995 [−1.039,−0.952]; synapse/axon −1.020) = **stoichiometric mC-for-hmC exchange, the dehydroxymethylase-like DNMT3A signature**, whereas the genome-wide deviation from −1 is driven by non-neuronal genes (−0.949). Neuronal genes resemble TET-KO *least* (rho 0.137 vs non-neuronal 0.246). [INTERPRETATION] neuronal loci undergo direct 5hmC→5mC conversion rather than TET-inhibition-plus-de-novo.
- **The R5 "figure with a lot of log2 plots":** this is the **Sec 59 (quadrant log2 scatters: K119ub vs MeCP2 vs 5mC vs 5hmC vs K27ac vs K27me3)** and the **Sec 60 mirror-image profiles** (MeCP2-Up = K119ub↑/5mC↑/5hmC↓/K27ac↓/K27me3↑; MeCP2-Down mirror but flat K119ub). The "euchromatic→heterochromatic relocation" is encoded as: MeCP2-Up loci gain K27me3 (+1.11) and lose K27ac (−0.72), i.e. MeCP2 concentrates where chromatin shifts toward Polycomb/heterochromatin.
- **Sec 61h–k (neuronal enrichment of the MeCP2-Up+K119ub-Up program):** strict 72-gene quadrant → 1 GO term (power-limited); **relaxed + genome-wide background → 1,117 terms, synaptic/axon-guidance top** (synapse assembly q=1.9e-12, axon guidance q=3.1e-11). Peak-level (Sec 61j): 818 GO terms, top synapse assembly q=5.3e-22 — robust at peak resolution. **GSEA (Sec 61k) is weaker** (MeCP2 ranking = 1 term; K119ub ranking strongest = tRNA/RNA-processing, not synaptic) → the strong ORA neuronal signal is partly a peak-location artifact, not a fold-change-driven program. (Caveat: neuronal enrichment is background- and weighting-dependent.)

> **Tension with the paper's word "Exclusively."** The verified data say MeCP2 redistribution is **predominantly** methylation-/Polycomb-driven and **enriched** at neuronal genes (OR 1.73–2.34), but Sec 74 explicitly shows MeCP2-Up tracks *coordinated methylation* (OR 5.16) **far more than neuronal identity** (OR 1.73), and Sec 75 shows MeCP2-Down peaks *also* land at synaptic genes (redistribution, not exclusive gain). The "exclusively at neuronal genes" framing overstates the verified specificity — recommend softening to "preferentially/disproportionately at neuronal genes."

### R5 — supported vs missing panels

| Planned R5 panel | Status |
|---|---|
| "The figure with a lot of log2 plots" | **Supported** — Sec 59 (7 quadrant scatters), Sec 60 (mirror profiles) |
| Neuronal-gene specificity of MeCP2 relocation | **Supported but qualified** — Sec 74/76/78; "exclusively" overstates (OR 1.73 neuronal vs 5.16 methylation) |
| Euchromatin→heterochromatin relocation | **Supported (indirect)** — Sec 60 (MeCP2-Up gains K27me3/loses K27ac), Sec 75 (gene-body→distal redistribution); no single direct "euchromatic→heterochromatic" track table |

---

# Methodological Validation

## Permutation tests (Sec 34–37)

- **Sec 37 (load-bearing, gene-level label-shuffle, ntimes=100,000):** **15/15 gene-level Fisher/O-E enrichments Confirmed; 0 Weakened.** Strongest |perm z|: mC×expression z=−31.93, mC×ATAC z=−30.12, hmC×ATAC z=+28.89, mC×K119ub-DiffBind z=+28.75. The four DiffBind cascade directions, the compartment-shift coupling (B→A z=−25.77; A→B z=+5.09), and the loop/MeCP2 anchor tests (K119ub×hyper-at-anchor z=+12.45; loop-dir×MeCP2-up z=−4.30) all pass. 13/15 hit the empirical p-floor (~1e-04). [INTERPRETATION] associations survive a null preserving per-chromosome marginals → driven by genuine co-assignment, not chromosomal co-residence.
- **Sec 34–36 (interval-level regioneReloaded):** permutation computations and heatmaps completed (z-scores up to 65 in Sec 34) but **crashed before emitting their Fisher-vs-permutation comparison TSVs** → effectively visual-only/partial in run-5; their confirmation is read from the Sec 37 gene-level analogues. (See Caveats.)

## Field cross-species replication (Sec 45)

Field et al. (2019, human BAP1-KD uveal melanoma) chr8 methylation hotspot: 81/85 genes mapped to mouse orthologs. **Gene-body mC concordance enriched over chance: Fisher p=0.00887; hypergeometric enrichment of Field orthologs among our DMRs p=4.61e-6; binomial vs 50% p=0.0363.** **Promoter level: p=0.40 (NS)** — replicates at gene-body, collapses at promoter. Trisomy-8 ruled out (45k diagnostic). [INTERPRETATION] the conserved BAP1 methylation response is a gene-body phenomenon, evolutionarily conserved.

## CTCF insulator model (Sec 47)

Lost CTCF loop anchors are hypermethylated at flanking **dynamic** CpG regions (shores+shelves): mC-hyper lost-vs-gained **OR=3.28** (p=5.4e-24); hmC-hypo OR=2.08; coordinated OR=2.18. **CpG islands null as predicted (OR=0.84, p=1).** Distance-stratified: OR=11.21 at <200 kb, NS >1 Mb (CMH common OR=2.87). Loop-logFC vs anchor mC Spearman rho=−0.244 (p=4.2e-9). [INTERPRETATION] local hypermethylation at CTCF sites weakens short-range loops — the IDH-glioma insulator-decay logic, linking methylation to 3D loop loss.

## HOMER motif grammar (Sec 49, 50)

H2AK119ub-gained sites enriched for **neuronal bHLH motifs** (Atoh1 FE 1.33, NeuroD1 FE 1.49, Atoh7, BHLHA15); same grammar marks lost-H3K27ac sites; YY1 specific to the DiffBind K119ub gained contrast. Coordinated mC↑/hmC↓ DMRs enriched for **zinc-finger/KLF/Sp** GC-box motifs (ZEB2, Maz, KLF15/14/5/6) — methyl-sensitive TF context. [INTERPRETATION] connects Polycomb-ubiquitination spread to cerebellar neuronal regulatory grammar.

---

# Cross-Reference Index

## A. Analysis section → paper Results section(s)

| Section(s) | Title | Paper R-section(s) |
|---|---|---|
| 01 | QC overview | R3 (QC), R2 (bulk baseline) |
| 02 | Sample correlation | R3 (sex/replicate structure) |
| 03 | DMR statistics by region | **R2, R3** |
| 04 | Volcano plots | **R2** |
| 05 | Coordinated changes | **R2** |
| 06 | Top genes / Venn | R2, R4 |
| 07 | Effect-size distributions | **R2** |
| 08 | GO/KEGG enrichment | **R2** |
| 09 | Summary statistics | R2, R3 |
| 10 / 10f | Chromatin-state classification | **R2, R4** |
| 11 | MeCP2 correlation | **R1, R4** |
| 12 | ATAC correlation | R2 (support) |
| 13 | ATAC chromatin + loops | R2/R4 (support) |
| 14 | H2AK119ub correlation | **R2** |
| 15 | hmC cross-mark | R2/R4 (support) |
| 16 | Raw concordance | R2 (support) |
| 17 | K119ub honest assessment | **R2** |
| 18 | K119ub BigWig signal | **R2** |
| 19 | H3K27ac peak analysis | **R2, R4** |
| 20 | Coordinated RNA-seq | R2 |
| 21 | Discordant Q2 vs Q4 | R2 |
| 22 | Demethylation ratio | **R2** |
| 23 | Baseline-5hmC predictor | **R2** |
| 24 | DNMT3A vs TET prediction | **R2** |
| 25 | Delta-ratio models | **R2** |
| 26 | TET-KO comparison | **R2** |
| 27 | Methylation × loop anchors | R5 (support), R2 |
| 28 | Coordinated Q4 characterization | R2 (support) |
| 29 | A/B compartment mapping | **R2** |
| 30 | Polycomb target enrichment | **R2** |
| 31 | MeCP2 × loop anchors | **R1**, R5 (support) |
| 32 | CHG exploratory | **R3** |
| 33 | Multi-mark DiffBind | **R2** |
| 34–36 | Interval permutation tests | Methods validation |
| 37 | Gene-level permutation | Methods validation |
| 38–41 | H3K36me2/3 + DNMT3A/B | **R2** (support) |
| 42 | Max-significance genes / PCA | R2, R3 (PCA) |
| 43 | CG exploratory / chrX | **R3** |
| 44 | Allele-specific methylation | R3 |
| 45 | Field chr8 cross-species | Methods validation |
| 46 | Genome browser loci (Gviz) | **R4** (bigwig example loci — visual-only) |
| 47 | CTCF anchor methylation | Methods validation / R2 |
| 48 | CpG-island ubiquitination | R3 (support) |
| 49 / 50 | HOMER motifs | Methods validation |
| 51–58 | Non-CG / MeCP2 CG-corrected | **R4** |
| 59 | Quadrant log2 scatters | **R4, R5** |
| 60 | MeCP2-status epigenetic profiles | **R1, R5** |
| 61 | Stoichiometry mechanism | R2, R5 |
| 61h–k | MeCP2-Up+K119ub-Up neuronal enrichment | **R5** |
| 62 | Multifeature chromatin regression | **R1, R4** |
| 63 | MeCP2 master heatmap | R1, R4 |
| 64 | Global methylation levels | **R2** |
| 65 | MeCP2-methylation scale | **R1, R4** |
| 66 | Subcompartment methylation | R2 (support) |
| 67 | MeCP2/K119ub at unmethylated genes | **R1, R4, R5** (KEY) |
| 68–71 | modC/MeCP2 by mark; ratio-vs-K119ub | **R1, R4** |
| 72 / 73 | Neuronal K119ub characterization + remodeling | **R5** |
| 74 | Gene-set overlaps | **R5** |
| 75 | MeCP2 signal reconciliation | **R1, R5** |
| 76 | Triple overlap / synapse specificity | **R5** |
| 77 | MeCP2 aging trajectory | **R1** |
| 78 | Unbiased-neuronal stoichiometry | **R2, R5** |

## B. Master quantitative table (every key OR / p / % / AUC / slope, with source)

| Metric | Value | p / q | Source |
|---|---|---|---|
| Gene bodies tested | 20,969 | — | Sec 03/09 |
| 5mC significant | 10,775 (51.4%) | q<0.05 | Sec 03 |
| — hyper / hypo | 7,513 / 3,262 | — | Sec 03 |
| 5hmC significant | 11,484 (54.8%) | q<0.05 | Sec 03 |
| — up / down | 1,963 / 9,521 | — | Sec 03 |
| Co-significant (both) | 8,371 | — | Sec 05/06 |
| Coordinated mC↑/hmC↓ | 6,589 (**78.7%**) | — | Sec 05/09 |
| Q1/Q2/Q3/Q4 | 411 / 1,255 / 116 / 6,589 | — | Sec 20/21/28 |
| Net mean 5mC change (sig) | +1.72% (hyper-only +3.45%) | — | Sec 07 |
| Net mean 5hmC change (sig) | −1.66% (hypo-only −2.29%) | — | Sec 07 |
| Bulk 5mC ctrl→mut | 72.216% → 72.534% (+0.317%) | t p=0.041 | Sec 64 |
| Bulk 5hmC ctrl→mut | 10.192% → 9.797% (−0.395%) | t p=0.018 | Sec 64 |
| Bulk total modC | 82.565% → 82.487% (NS) | p=0.492 | Sec 64 |
| GO BP #1 (hyper genes) | RNA splicing, 274 genes | p.adj 6.27e-54 | Sec 08 |
| KEGG #1 | Autophagy-animal, 113 genes | p.adj 8.82e-24 | Sec 08 |
| **Demethylation ratio decreased** | 71.36% of genes | — | Sec 22 |
| WT→KO median ratio | 0.1275 → 0.1193 | — | Sec 22 |
| Baseline-5hmC predictor AUC | **0.762** | — | Sec 23 |
| K119ub-signal predictor AUC | 0.573 | — | Sec 23 |
| TET-impediment model AUC | **0.793** | — | Sec 24 |
| DNMT3A-recruitment model AUC | 0.696 | — | Sec 24 |
| DeLong TET vs DNMT3A | — | **9.43e-49** | Sec 24 |
| Baseline-5hmC std β (hyper) | +1.251 (rank 1) | — | Sec 24 |
| K119ub std β (hyper, TET model) | **−1.049** (negative) | — | Sec 24 |
| TET-KO absolute attenuation | ~3.3% (relative ~8.7%) | — | Sec 26 |
| K119ub DiffBind logistic OR (4-mark) | **4.707** | 2.4e-18 | Sec 33 |
| 4-mark logistic AUC | 0.818 (+baseline 0.869) | — | Sec 33 |
| K119ub 6-mark logistic OR | **10.29** | 1.0e-10 | Sec 41 |
| H3K36me3 6-mark OR | 1.00 (NS) | 0.996 | Sec 41 |
| DNMT3A:DNMT3B attribution | 1,630 : 84 (~19:1) | — | Sec 41 |
| K119ub differential overlap OR | 4.46 | 5.8e-97 | Sec 14 |
| K119ub conditional gain enrichment | +14.2 pp (47.8% vs 33.6%) | — | Sec 17 |
| K119ub global BigWig median log2FC | +0.0070 | 1.8e-20 | Sec 18 |
| ATAC DMR-overlap OR | 0.068 | 4.4e-178 | Sec 12 |
| mC-hyper → A compartment OR | **13.642** | ~0 | Sec 29 |
| hmC-loss → A compartment OR | 9.774 | ~0 | Sec 29 |
| Polycomb × mC-hyper OR | **0.0633** | ~0 | Sec 30 |
| Polycomb × mC-hypo OR | 9.797 | ~0 | Sec 30 |
| Active_Promoter hyper rate | 71.8% | — | Sec 30 |
| MeCP2 DMR-overlap OR | 5.13 | 1.27e-33 | Sec 11 |
| MeCP2 gene-level 2×2 OR | 0.92 (NS) | 0.077 | Sec 11 |
| MeCP2 binding R² (CG vs chromatin) | 0.017 vs **0.246** | — | Sec 62 |
| MeCP2 differential R² (CG vs chromatin) | 0.013 vs 0.170 | — | Sec 62 |
| K119ub std β (MeCP2 binding) | +0.199 (top) | <2.2e-16 | Sec 62 |
| MeCP2-no-meth genes K119ub gained | 72.8% vs 45.9%, **OR=3.15** | 1.82e-24 | Sec 67 |
| MeCP2 fold vs K119ub (no-meth) rho | +0.241 | 4.4e-06 | Sec 67 |
| K119ub vs MeCP2 (all) rho | 0.187 | 3.5e-96 | Sec 71 |
| Ratio vs MeCP2 rho | 0.014 (NS) | 0.113 | Sec 71 |
| K119ub variance-unique (MeCP2) | 7.3% (ratio 0.0%) | — | Sec 71 |
| MeCP2-Up / Down peaks | 7,686 / 1,200 | — | Sec 56/75 |
| MeCP2-Up peaks distal-intergenic | 51.7% (promoter 2.2%) | — | Sec 75 |
| Loop lost-vs-gained hyper OR | 2.589 | 1.1e-35 | Sec 27 |
| Loop logistic Lost OR | 2.290 | 5.5e-32 | Sec 27 |
| MeCP2 gain × lost-loop enrichment | 1.51 (vs gained 0.65) | 1.3e-03 | Sec 31 |
| CTCF lost-anchor mC-hyper OR (dynamic) | 3.28 | 5.4e-24 | Sec 47 |
| CTCF islands (null) OR | 0.84 | 1 | Sec 47 |
| Field gene-body concordance | — | 0.00887 | Sec 45 |
| Field promoter concordance | — | 0.40 (NS) | Sec 45 |
| ASM mutant/control ratio | 1.95× | — | Sec 44 |
| Global Δ5mC~Δ5hmC slope | **−0.959** (≠−1) | — | Sec 61 |
| Neuronal (broad) slope | −0.995 (=−1) | — | Sec 78 |
| Neuronal Δtotal (broad / narrow) | −0.0022 / **+0.012** | — | Sec 78 |
| Neuronal K119ub top-decile OR | 1.70 (Sec 72) / 2.34 (Sec 76) | — | Sec 72/76 |
| Synapse/axon extra K27me3 loss | Δ−0.044 | 2.95e-3 | Sec 76 |
| Mutant aging-UP peaks | 23,117 (vs ctrl 10,930; 2.1×) | — | Sec 77 |
| Coordinated × MeCP2-Up OR | 5.16 | 2.8e-12 | Sec 74 |
| Neuronal × MeCP2-Up OR | 1.73 | 0.034 | Sec 74 |
| Permutation tests confirmed | **15/15** | floor ~1e-04 | Sec 37 |
| Ecker CH at non-CG candidates | 0.0159 vs 0.0129 | 2.6e-60 | Sec 57 |
| Non-CG dose-response (Jonckheere inc.) | NS | ~1.0 | Sec 58 |
| CHG significant DMRs | 320 (chr7/8 76%) | — | Sec 32 |

## C. Figure inventory → paper section mapping (key plot folders)

| Plot folder(s) | Content | Paper section |
|---|---|---|
| `04_volcano_plots`, `04a/04b` | 5mC/5hmC volcanoes | R2 |
| `05a_mc_hmc_scatter`, `05b/05c/05d` | coordinated quadrant + top genes | R2 |
| `03a/03b/03_dmr_region_statistics` | region & direction bars | R2/R3 |
| `07_effect_size_distributions` | effect-size violins | R2 |
| `08a–08h` | GO/KEGG dotplots | R2 |
| `10a–10f`, `10f_expanded` | chromatin-state distribution | R2/R4 |
| `11a–11i` | MeCP2 overlap/scatter/LM/GLM | R1/R4 |
| `14a–14e`, `17a–17c`, `18a–18d` | K119ub overlap/honest/BigWig | R2 |
| `22a–22h`, `23a–23e`, `24a–24k`, `26a–26i` | TET ratio & mechanism | R2 |
| `27a–27e`, `29a–29g`, `30a–30f` | loop/compartment/Polycomb | R2/R5 |
| `31a–31d` | MeCP2 × loops | R1/R5 |
| `33a–33f`, `37a–37d` | DiffBind spine + permutation | R2/Methods |
| `45a–45k`, `47a–47e`, `49*`, `50*` | Field/CTCF/HOMER | Methods validation |
| `46/Syt1_locus` etc. (Gviz) | **bigwig example loci** | R4 (visual-only) |
| `59a–59h`, `60a–60e` | log2 quadrant scatters + mirror profiles | R4/R5 |
| `62a–62g`, `67a–67d`, `71a–71d` | regression / K119ub-unmethylated | R1/R4 |
| `64a–64c`, `66a–66d` | global & subcompartment methylation | R2 |
| `72a–72f`, `73*`, `74a–74_composite`, `75a–75_composite`, `76a–76d`, `77a–77_composite`, `78a–78g` | neuronal characterization & aging | R5/R1 |
| `42_pca_*` | per-sample methylation PCA | R3 (sex PCA candidate) |

---

# Data Quality & Caveats

## 1. Run-5 vs earlier-run number discrepancies (RESOLVED — cite run-5)

**The single most important caveat: `FIGURES.md` and `TODO.md` carry OLD-run (run-4-era) numbers throughout and must NOT be cited.** The Sec 01–09 group doc independently recomputed every headline number from the run-5 BED files and confirmed `summary_statistics.txt` IS the canonical run-5 data.

| Quantity | **Run-5 canonical (cite this)** | Stale FIGURES.md/TODO.md (do NOT cite) |
|---|---|---|
| 5mC significant | **10,775** | 8,836 |
| 5hmC significant | **11,484** | 9,930 |
| Co-significant / Venn overlap | **8,371** | 6,722 / 6,750 |
| Coordinated % | **78.7%** | 84.6% / 85% |
| Net mean 5mC / 5hmC | **+1.72% / −1.66%** | +2.27% / −2.08% |
| Hyper-only mean (n) | **+3.45% (n=7,513)** | +3.96% (n=6,635) |
| 5mC sample correlation | **0.87–0.90** | 0.76–0.79 |
| 5hmC sample correlation | **0.64–0.71** | 0.48–0.51 |
| Q4 / Q2 quadrant | **6,589 / 1,255** | 5,708 / 755 |
| Region bars (gene body) | **10,775** | 8,836 |
| Active_Promoter % of DMRs | **45.5% (60.7% of hyper)** | 49.9% |
| Loop lost-vs-gained hyper OR | **2.589** | 2.54 |
| mC-hyper → A compartment OR | **13.642** | 14.71 |
| Polycomb × mC-hypo OR | **9.797** | 8.71 |

> **The `summary_statistics.txt` header reads `Samples: 4 (2 Control, 2 Mutant)` — this is a COSMETIC hardcoded-template bug** (`section_09_summary.R` line 50); ALL numeric fields populate live from the 8-sample run-5 objects. The BioQC JSON is literally named `..._8_samples_...`. Fix the header string; the data are 8-sample run-5.

**The 10,775/11,484 vs 8,836/9,930 conflict named in the task is RESOLVED: the canonical value is 10,775 / 11,484 (run-5).** FIGURES.md's 8,836/9,930 are stale run-4 numbers, NOT reproducible from current tables.

## 2. The "92.3% coordinated" CLAUDE.md figure is [UNVERIFIED]

Both `biomodal/CLAUDE.md` and `downstream/CLAUDE.md` state **"92.3% of co-significant genes show coordinated mC↑/hmC↓."** This is **NOT reproducible** from any run-5 table. Three different values appear across sources: **78.7% (run-5 TSV, canonical)**, 85% (FIGURES.md, stale), 92.3% (CLAUDE.md, unverified — possibly a different denominator excluding Q1/Q3, a large-effect subset, or a prior run). **Cite 78.7%** (6,589/8,371). `[UNVERIFIED: 92.3%]`.

## 3. Sample-size limitations

- **n=4 effective replicates** (4 ctrl + 4 mut = 2 sexes × 2 batches per condition). Paired Wilcoxon in Sec 64 is underpowered (min two-sided p=0.125 at n=4 pairs); the **paired t-test is the operative test** for the global shift.
- **Sex/genotype remain partially confounded at n=2/group per sex** even with the sex covariate (run-5). run-1 (sex covariate, n=2/group) found NO significant DMRs precisely because of this confound; run-5's 8 samples are what made the sex covariate viable.
- ASM (Sec 44) uses uncorrected per-sample testing — absolute counts are exploratory.
- Small-n cells: triple-overlap n=16 (Sec 74/76), MeCP2-Up stoichiometry slope n=79 (wide CI [−1.31,−0.54], Sec 78), several quadrant per-class cells n=4–17 — interpret cautiously.

## 4. Pending / blocked / partial analyses

- **NO P12 (early-timepoint) DUET methylation data exists** (TODO.md §12). The methylation analysis is **adult-only**. R1's developmental P12-vs-adult contrast comes from the MeCP2 CUT&RUN arm and the Hi-C arm, NOT methylation. This is the largest scope gap relative to the paper's title.
- **Sec 77 (MeCP2 aging trajectory) is NO LONGER BLOCKED.** The downstream CLAUDE.md / MEMORY flag it as "blocked on Jai's young-vs-adult DiffBind data" — **that status is stale**; both input files arrived (`MeCP2_{ctrl,mut}_adultvsyoung_diffbind_results.txt`, ~53 MB each, dated Jun 23) and all four Sec 77 tables are populated. The section ran fully.
- **Sec 34/35/36 (interval permutation tests) are partial:** SLURM jobs crashed at the Fisher-vs-permutation table step (`replacement has length zero` / `differing number of rows`); only `.rds` caches + heatmap panels exist, no comparison TSVs. Confirmation is read from the gene-level Sec 37 (which ran fully, 15/15 confirmed). Treat 34–36 as visual-only.
- **Sec 73 table-prefix mismatch:** on-disk tables carry stale `72g_*` prefix (the committed script writes `73_*` but was not re-run after rename commit 8ba6681c). Data are current and valid; only filenames are stale.
- **Sec 65→67 provenance gap:** the 359-gene `mecp2_only` list is read by Sec 67 but the current Sec 65 script does not write it — generated by an external/earlier script. Count is internally consistent (359, matches Fisher cell), but generating step is undocumented.

### Open items from FURTHER_ANALYSIS.md not yet done (manuscript-discussion candidates)

Several analyses proposed in FURTHER_ANALYSIS.md remain **not done** and are candidate gaps to acknowledge: methylation entropy / read-level heterogeneity (uniform vs stochastic mechanism); repeat/transposable-element methylation (DNMT3A-redistribution prediction at repeats); spatial chromosomal autocorrelation (K119ub cis-spreading test); cell-type deconvolution (Math1-Cre targets granule progenitors; bulk signal is a cell-type-weighted average — relevant to interpreting the minority hypomethylated genes as real vs composition artifact). Note: items 1 (demethylation ratio), 2/CTCF (Sec 47), 5 (Polycomb, Sec 30), 8 (Field, Sec 45), 9 (non-CG, Sec 51–58) from that doc **were subsequently done**.

## 5. Console-only / unexported statistics (cannot be re-verified from tables)

These are real results but live only in script stdout / figure annotations, not TSVs — quote with the caveat or re-run to capture:
- Sec 28: nine cross-dimension Fisher/Wilcoxon stats; Sec 31: 31a Fisher OR / 31c Wilcoxon / 31d Spearman.
- Sec 33: 4-mark AUC=0.818 / extended 0.869 (coefficients ARE in TSV; AUC log-only).
- Sec 41: 4-/6-mark AUC; 40b/40d quadrant & Venn counts.
- Sec 76a: neuronal-vs-non-neuronal 2-group violin p-values (only 4-group tests in TSV). `[UNVERIFIED: ATAC 6.95e-15 / K27ac 0.027 / K27me3 5.53e-4 per root prose]`.
- Sec 58: MeCP2-Up-only correlation rho (only OLS-vs-QR concordance 0.906 in TSV).
- Sec 78e: narrow/broad Δtotal means (re-derived: +0.012 / −0.0022 / −0.0077, confirmed).

## 6. Specific [UNVERIFIED] numbers inherited from group docs

| Claim | Source | Status |
|---|---|---|
| 92.3% coordinated | CLAUDE.md | Not in tables; canonical 78.7% |
| MeCP2 DMR OR=5.60 | FIGURES.md Fig 14 | Verified value 5.13 (Sec 11) |
| 11e Obs counts 651/1,830 | FIGURES.md Fig 15 | OR/p match (~0.92/0.08) but per-cell counts don't |
| Sec 23 K119ub AUC=0.592; dose-response rho=−0.586 | TODO.md | Tables give AUC 0.573; rho not exported (direction confirmed) |
| Sec 22 Cliff's delta=0.455; 72.5% decreased | TODO.md | Table gives 71.36% decreased; Cliff's delta not exported |
| Sec 24 hyper rate 46.2%; 11,936 genes | TODO.md | Table gives 51.9% / 11,937 |
| Sec 27 OR=2.54 / logistic 2.42 | TODO.md/FIGURES.md | Tables give 2.589 / 2.290 |
| Sec 20d mC×expression OR=0.03, p=5.90e-118, n=1,330 | FIGURES.md | Not reconstructable from group TSVs (source gene set not exported); synthesis cites Sec 21 expression Q2-vs-Q4 OR≈95 instead, not the 20d figure |
| Sec 29 mC-hyper→A OR=14.71 | TODO.md/FIGURES.md | Table gives 13.642 |
| Sec 30 Polycomb hypo OR=8.71 | TODO.md/FIGURES.md | Table gives 9.797 |
| Sec 43 per-chromosome counts | — | NOT independently read from TSV (panel inventory only) |
| Sec 45 trisomy-8 coverage ratio / MWU p | — | Console/figure-only, not in TSV |
| Sec 76a 2-group p-values | root prose | Not in TSV (4-group only) |
| DeLong "p<2.2e-16" | TODO.md/script | Display floor; true value 9.43e-49 (Sec 24 TSV) |

## 7. Labeling & framing conventions for the manuscript

- **MeCP2 is CUT&RUN, not ChIP** — use "signal"/"mark"/"binding," never "ChIP" (all group docs comply).
- **Bivalent = Repressed Enhancer** per PI (K4me1+K27me3); the 12-state model retains `Bivalent_Enhancer` for backward compatibility only.
- **CHG framing tension:** FIGURES.md says non-CG "No Signal"; Sec 32 reports 320 significant CHG DMRs. Reconcile as "statistically significant but biologically negligible (~100× smaller, chr7/8-skewed, 61% CG-discordant)."
- **"Exclusively at neuronal genes" (R5 title) overstates** the verified specificity (MeCP2-Up OR 5.16 for methylation vs 1.73 for neuronal identity) — soften to "preferentially/disproportionately."
- **"Binding active promoters mostly" (R1)** is in tension with Sec 75 (MeCP2-Up peaks 52% distal-intergenic, 2.2% promoter); the *methylation* hyper is at active promoters, the *MeCP2 gain* is distal.
- **Region counts mix per-gene (gene bodies) and per-feature (other regions)** — state count basis explicitly (Sec 03 flag).

---

## Bottom line for the manuscript

The verified run-5 data robustly support **R2, R4, R5, and the Methodological Validation** with strong, permutation-confirmed effect sizes. **R1's developmental P12-vs-adult MeCP2 contrast is supplied by the CUT&RUN/Hi-C arms, not the methylation pipeline** (no P12 DUET data) — the methylation pipeline's R1 contribution is the adult aging-amplification (Sec 77) and redistribution (Sec 75) evidence. **R3's sex-difference PCA panel is the least-supported claim** — no sex-stratified table exists in the verified docs, and it sits in tension with the Methods statement that no main finding varied by sex. Cite run-5 numbers exclusively; the 78.7% coordinated figure (not 85% or 92.3%) and the 10,775/11,484 DMR counts (not 8,836/9,930) are canonical.
