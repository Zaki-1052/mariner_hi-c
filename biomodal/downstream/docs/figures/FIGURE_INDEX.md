# FIGURE INDEX — BAP1-KO Cerebellum Methylation Paper

Master index of the eight consolidated publication figures (Fig 1–7 + Fig S1). Each figure has a standalone findings doc in this directory (`figureN_findings.md`); this index gives the one-line takeaway, paper R-section mapping, and a one-line-per-panel summary with the key verified numbers. All values are run-5 canonical (deep-seq, 8 samples = 2 ctrl + 2 mut per sex, sex covariate; n = 4 effective replicates/condition). MeCP2 is CUT&RUN, never ChIP.

---

## Narrative flow

The figures trace a single causal cascade — **BAP1 loss → H2AK119ub accumulation → coordinated 5mC↑ / 5hmC↓ at gene bodies (a TET-demethylation block) → MeCP2 redistribution → preferential disruption of neuronal genes**. **Figure 1** documents *what* happens: a reciprocal gene-body 5mC-gain / 5hmC-loss signature with a flat total modified-cytosine pool (the molecular fingerprint of a blocked active-demethylation pathway), reproducible across all four replicate pairs and both sexes. **Figure 2** pins *why*: H2AK119ub gain — the mark whose eraser is BAP1 — is the dominant predictor of where hypermethylation lands (OR ≈ 4.7), and that hypermethylation is restricted to active/euchromatic (A-compartment) chromatin while Polycomb heterochromatin is actively protected. **Figure 3** adjudicates the *enzymatic mechanism*: the signature is caused by impaired (substrate-limited) TET demethylation, not by K119ub directly recruiting DNMT3A (TET-impediment AUC 0.793 vs DNMT3A-recruitment 0.696, DeLong p = 9.43e-49), with neuronal loci undergoing a stoichiometric (slope ≈ −1) 5hmC→5mC exchange. **Figure 4** shows the *downstream consequence*: MeCP2 binding redistributes toward distal-intergenic sites, preferentially overlaps hypermethylated DMRs (OR = 5.13), and its normal developmental ramp is amplified (2.1× more age-increased peaks in mutant) onto the same coordinated-methylation gene set. **Figure 5** delivers the *reader-vs-follower verdict*: MeCP2 redistribution tracks chromatin state — above all H2AK119ub — not CG methylation (chromatin R² = 0.246 vs CG R² = 0.017). **Figure 6** explains the *neuronal specificity*: neuronal genes are preferentially (not exclusively) affected because they sit in constitutively K119ub-rich Polycomb chromatin (the natural BAP1 substrate), with synapse/axon genes undergoing selective Polycomb de-repression. **Figure 7** lines up the single strongest number from each step as the data companion to the model schematic. **Figure S1** is the robustness/generalizability backstop — the phenotype survives a genomic null, replicates across species, propagates into 3D loop anchors and allelic asymmetry, and behaves as predicted on an unbiased gene set.

---

## Figure 1 — The Coordinated Gene-Body Methylation Phenotype
**Paper mapping:** R2 + R3. **Takeaway:** BAP1-KO cerebellum has a coordinated, reciprocal gene-body 5mC-gain / 5hmC-loss signature with a flat total modified-cytosine pool — the fingerprint of a blocked TET-mediated active-demethylation pathway.

- **A** (Sec 64): Paired bulk autosomal CpG shift — 5mC +0.317% (Cohen's d +2.79, p=0.041), 5hmC −0.395% (d −4.54, p=0.018), total modC −0.079% (NS); every mutant up in 5mC / down in 5hmC = redistribution, not global hypermethylation.
- **B** (Sec 03): % significant 5mC DMRs by region — gene bodies 51.4% >> CpG shores 30.2% > shelves 23.8% > promoters 8.3% > islands 5.0% > TSS 1.4%; a gene-body phenomenon.
- **C** (Sec 04): Dual volcano — 10,775 significant 5mC DMRs (70% hyper / 30% hypo) and 11,484 significant 5hmC DMRs (17% up / 83% down); mirror-image clouds.
- **D** (Sec 05): mC-vs-hmC quadrant — 6,589 of 8,371 co-significant genes (78.7%) fall in the coordinated mC↑/hmC↓ quadrant (reverse 15.0%, same-direction 6.3%).
- **E** (Sec 07): Effect-size densities — median change among significant DMRs ≈ +2.03% 5mC / −1.83% 5hmC (computed live); modest per-gene but reciprocal.
- **F** (Sec 42): Side-by-side 5mC | 5hmC PCA — PC1 = condition (≈42.6–48.8% variance), PC2 = sex (≈11.1–11.8%); genotype dominant, sex secondary (R3).
- **G** (Sec 46): Gviz Syt1 browser (chr10:108,747,000–108,812,000) — +17.3% 5mC, −15.0% 5hmC with concurrent H2AK119ub; whole cascade visible in cis.

---

## Figure 2 — H2AK119ub Drives Hypermethylation at Active Chromatin
**Paper mapping:** R2 (mechanism). **Takeaway:** H2AK119ub gain is the single strongest predictor of where hypermethylation lands, and that hypermethylation is restricted to active, euchromatic (A-compartment) chromatin while Polycomb heterochromatin is protected.

- **A** (Sec 10): Stacked direction bar by chromatin state — Active_Promoter 93.0% hyper (4,562/4,906); Repressed_Promoter 94.4% hypo (1,621/1,718); mirror-image split.
- **B** (Sec 33): 4-mark logistic forest — H2AK119ub OR = 4.71 [3.33–6.66], p=2.4e-18 (dominant); H3K27me3 OR = 1.44; H3K27ac OR = 0.48 and ATAC OR = 0.22 (both protective); model AUC = 0.818 (run-log).
- **C** (Sec 29): A/B compartment enrichment — mC-hyper → A OR = 13.6 [11.92–15.67]; hmC-loss → A OR = 9.8; mC-hypo → B OR = 1.7; hmC-gain → B OR = 3.0.
- **D** (Sec 30): Polycomb exclusion — Polycomb × mC-hyper OR = 0.063 [0.053–0.075] (excluded); mC-hypo OR = 9.80; hmC-gain OR = 11.39; hmC-loss OR = 0.131; falsifies the naive "Polycomb gets hypermethylated" prediction.
- **E** (Sec 17): Honest K119ub breakdown — 58.16% of DMR genes carry no K119ub peak; among peak-bearing genes, mC-Up gains K119ub at 47.8% vs 33.6% background (+14.2 pp); coupling real but partial.

---

## Figure 3 — TET-Impediment, Not DNMT3A Recruitment
**Paper mapping:** R2 (mechanistic adjudication). **Takeaway:** The coordinated phenotype is caused by an impaired (substrate-limited) TET demethylation pathway, not by K119ub recruiting DNMT3A; at neuronal genes the mC-for-hmC exchange is stoichiometric (slope ≈ −1, dehydroxymethylase-like).

- **A** (Sec 23): Baseline-5hmC dose-response (deciles) — baseline-5hmC predictor AUC = 0.762 [0.755–0.769] (vs K119ub-signal 0.573), McFadden R² = 0.141; higher-substrate deciles lose the most 5hmC.
- **B** (Sec 24): Model comparison — TET impediment AUC = 0.793 [0.785–0.801] >> DNMT3A recruitment AUC = 0.696 [0.687–0.706], DeLong p = 9.43e-49 (Full 0.861, K119ub-only 0.658); K119ub enters the full model negative (β = −1.05), baseline 5hmC top positive (β = +1.25).
- **C** (Sec 22): Demethylation-ratio lollipop by state — Active_Promoter Δ = −0.0281 (strongest impairment, n=6,356), Active_Enhancer −0.0248, with reversal at Repressed_Promoter +0.0051 and Polycomb +0.0022 (ns); TET throttled at active chromatin.
- **D** (Sec 78): Stoichiometry scatter (Δ5mC vs Δ5hmC) — all genes slope = −0.959 [−0.978,−0.940] (≠ −1, n=20,915); neuronal broad = −0.995 [−1.039,−0.952] (≈ −1, n=4,118); synapse/axon −1.020; coordinated −1.288.
- **E** (Sec 26): Attenuation bar vs published TET triple-KO — BAP1-KO median Δ = −0.008 vs TET-KO −0.260; absolute attenuation 3.3% (relative 8.7%, n matched 20,851); "dimmer, not an off switch" (supplemental-eligible).

---

## Figure 4 — MeCP2 Redistribution and Developmental Amplification
**Paper mapping:** R1 (title figure). **Takeaway:** MeCP2 binding redistributes toward distal-intergenic sites, preferentially overlaps hypermethylated DMRs (OR = 5.13), and the normal developmental MeCP2 ramp is amplified (2.1× more age-increased peaks in mutant), landing on the coordinated mC↑/hmC↓ gene set.

- **A**: Adult MeCP2 volcano (CUT&RUN) — 202,650 total peaks; 7,686 UP / 1,200 DOWN at FDR<0.05 (6.4:1 skew toward gains).
- **B** (Sec 75): Peak annotation by direction — UP peaks 51.7% distal-intergenic / 41.0% intron / 2.2% promoter; DOWN peaks 61.8% intron / 19.5% distal / 8.0% promoter; gains land distally, losses are genic.
- **C** (Sec 11): DMR overlap — hyper DMRs 7.16% overlap MeCP2-Up vs 1.89% MeCP2-Down; Fisher OR = 5.13, p = 1.27e-33.
- **D** (Sec 77): Aging peaks by genotype — control 10,930 UP / 2,822 DOWN; mutant 23,117 UP / 10,646 DOWN; mutant 2.1× more age-increased peaks.
- **E** (Sec 77): Shared-aging fold scatter — 7,305 loci both genotypes age-gain; median aging fold 2.24 (mut) vs 1.83 (ctrl), +0.41 log2 (~+22%); 2,266 at neuronal genes.
- **F** (NEW): Mutant-unique age-increased genes (~1,654) vs coordinated mC↑/hmC↓ set over ~20,915-gene universe; Fisher OR/p computed at runtime → `fig4f_aging_methylation_overlap.tsv`; ties the developmental overshoot to the methylation phenotype.

---

## Figure 5 — MeCP2 Reads Chromatin State, Not Methylation
**Paper mapping:** R4. **Takeaway:** MeCP2 binding change is governed by chromatin state — above all H2AK119ub — not by CG methylation: CG explains <2% of MeCP2 binding while chromatin marks explain ~25%, and where MeCP2 rises without methylation change it is K119ub that is gained (OR = 3.15).

- **A** (Sec 62): R² bars — binding level CG 0.017 vs Chromatin 0.246 vs Full 0.260 (n=202,574 peaks); differential level CG 0.013 vs Chromatin 0.170; chromatin ~15× CG.
- **B** (Sec 62): Variance partition — binding level Chromatin-unique 24.3% vs CG-unique 1.5% (Shared 0.2%, Unexplained 74.0%); differential Chromatin-unique 16.0% vs CG-unique 0.3%.
- **C** (Sec 67): K119ub at MeCP2-up-no-methylation genes — 359 genes (356 tested) gain K119ub at 72.8% vs 45.9% genome; Fisher OR = 3.15, p = 1.8e-24.
- **D** (Sec 60): Mirror-image mark profiles — MeCP2-Up: K119ub +0.63, 5mC +0.049, 5hmC −0.017, K27ac −0.72, K27me3 +1.11; MeCP2-Down mirrors every mark except K119ub, which stays flat (−0.04, ns, wilcox_q=0.18) — the cleanest single discriminator.
- **E** (Sec 71): MeCP2-fold variance partition — K119ub-unique 7.3% vs methylation-ratio-unique 0.0% (Shared 0.1%, Unexplained 92.6%, n=12,149 genes).

---

## Figure 6 — Neuronal Genes Are Preferentially Affected
**Paper mapping:** R5. **Takeaway:** Neuronal-identity genes are preferentially (not exclusively) affected because they sit in constitutively H2AK119ub-rich Polycomb chromatin; the downstream MeCP2 redistribution tracks the methylation shift ~3× more than neuronal lineage, synapse/axon genes show selective Polycomb de-repression, and neuronal methylation change is stoichiometric.

- **A** (Sec 72): Neuronal fraction across constitutive-K119ub deciles — top decile 33.0% (713/2,161), Fisher OR = 1.70 [1.54,1.87], q=3.4e-25 vs genome baseline 23.5% (5,077/21,604); bottom decile 20.9% (OR 0.84).
- **B** (Sec 74): Pairwise Fisher ORs — Coordinated × MeCP2-Up OR = 5.16 [3.19,8.51], p_adj=2.8e-12; Neuronal × MeCP2-Up OR = 1.73 [1.05,2.82], p_adj=0.034; Neuronal × Coordinated OR = 1.05 (NS); methylation > lineage.
- **C** (Sec 76): Synapse/axon vs broader-neuronal DiffBind — H3K27me3 Δmedian = −0.044 (−0.132 vs −0.087), Wilcoxon p=3.0e-3 (p_adj=6.6e-3); ATAC p=0.65 (NS), H3K27ac p=0.76 (NS); selective Polycomb de-repression.
- **D** (Sec 78): Stoichiometry slopes by gene class — neuronal broad −0.995 [−1.039,−0.952] (stoichiometric), synapse/axon −1.020 (stoichiometric); all genes −0.959 and non-neuronal −0.949 (deviate); coordinated −1.29.

---

## Figure 7 — Mechanism Summary Data Panel (cascade)
**Paper mapping:** All R-sections (data companion to the Illustrator model schematic). **Takeaway:** Each panel quantifies one link of the BAP1 → H2AK119ub → TET-block → 5mC↑/5hmC↓ → MeCP2 → neuronal cascade so the whole chain reads as one internally consistent argument.

- **A** (Sec 18): BAP1 → K119ub — background median K119ub log2FC +0.007, p vs 0 = 1.81e-20; mC-Up +0.059 (n=6,993), mC-Down −0.080 (n=2,974); real but diffuse global rise.
- **B** (Sec 33): K119ub → hypermethylation — H2AK119ub OR = 4.71 [3.33,6.66], p=2.4e-18 (dominant); H3K27me3 1.44; ATAC 0.22, H3K27ac 0.48 (<1).
- **C** (Sec 23+24): TET impediment, not DNMT3A — baseline 5hmC AUC 0.762, TET impediment 0.793 [0.785,0.801], DNMT3A recruitment 0.696 [0.687,0.706], DeLong p = 9.43e-49.
- **D** (Sec 62): Methylation → MeCP2 — binding-level chromatin R² = 0.246 vs CG R² = 0.017 (~15×, n=202,574); K119ub top single predictor (β=+0.199).
- **E** (Sec 74+72): Convergence on neuronal genes — Coordinated × MeCP2-Up OR = 5.16 [3.19,8.51]; Broad-neuronal × MeCP2-Up OR = 1.73 [1.05,2.82]; constitutive K119ub top-decile × neuronal OR = 1.70 [1.54,1.87], q=1.3e-25.

---

## Figure S1 — Supplemental Evidence
**Paper mapping:** Methods / Support (supplemental to R2–R5). **Takeaway:** The phenotype is real, robust, and generalizable — it survives a genomic null, replicates across species, propagates into 3D loop anchors and allelic asymmetry, and behaves as the TET-block / stoichiometric-conversion mechanism predicts on an unbiased gene set.

- **A** (Sec 37): Permutation validation — 15/15 confirmed, 0 weakened (chr-stratified, 1e5 perms); strongest mC×expression z=−31.93, mC×ATAC z=−30.12, hmC×ATAC z=+28.89, mC×K119ub-DiffBind z=+28.75; 13/15 hit the p-floor.
- **B** (Sec 45): Field 2019 human BAP1 chr8 orthologs (n=81) — 21 concordant / 38 discordant by gene-body 5mC; Fisher p=0.0089; hypergeometric enrichment p=4.6e-6 (promoter-level NS, p=0.40); conserved, gene-body-centric.
- **C** (Sec 47): CTCF loop-anchor methylation — lost anchors 22.3% (242/1,084) vs gained 8.1% (115/1,427) hypermethylated flanks; OR = 3.28 [2.57–4.20], p=5.4e-24; islands null (OR≈0.84).
- **D** (Sec 48): CpG-island K119ub depletion — mC-hyper islands 17.2% (21/122) vs non-sig 30.2% (2,561/8,468) carry mutant K119ub peak (depletion, OR<1; gained-K119ub ≈0%); island hyper is amplification of pre-existing methylation, not de-novo.
- **E** (Sec 44): Allele-specific methylation doubling — control mean 11,295.5 vs mutant 22,080.25 ASM sites/sample (1.95×); 1,408 gene-body DMR genes harbour ≥1 ASM site (per-sample uncorrected → indicative).
- **F** (Sec 78): Unbiased broad-neuronal stoichiometry — neuronal broad Δ total-meth = −0.0022 (decreases, n=4,118), all genes −0.0014, non-neuronal −0.0012, synapse −0.0013; rises only at methylation-defined sets: coordinated +0.0070, Neur+Coord +0.0076, MeCP2-Up +0.0313 (n=79).

---

## Framing decisions (apply to every figure and caption)

- **"Preferentially / disproportionately," not "exclusively."** Neuronal genes are over-represented because they are constitutively K119ub-rich (the BAP1 substrate), not because the mechanism is neuron-restricted. MeCP2 redistribution tracks the coordinated methylation set (OR = 5.16) ~3× more strongly than neuronal identity (OR = 1.73), and Neuronal × Coordinated is not even significant (OR = 1.05). R5's original "exclusively" is softened throughout.
- **MeCP2 is CUT&RUN, never ChIP.** Use "binding," "signal," or "mark" for MeCP2 in all labels and captions. (H2AK119ub, H3K27me3, H3K27ac are the ChIP-derived histone marks; "ChIPseeker" appears only as a tool name in Figure 4's annotation step.)
- **Sex PCA reinstated (R3).** Figure 1 Panel F shows PC1 = condition (≈42.6–48.8% variance), PC2 = sex (≈11.1–11.8%) for both modalities — genotype is the dominant axis, sex a real but secondary structure, validating the run-5 sex covariate.
- **P12 MeCP2 gap.** The paper outline's "P12 = no DMRs" volcano requires a P12 ctrl-vs-mut MeCP2 DiffBind contrast not currently in the repo (only the young-vs-adult aging files are present). Figure 4 Panel A is the *adult* contrast; the developmental story is carried by the aging trajectory (Panels D–F), not a P12 volcano.
- **K119ub sign reconciliation.** K119ub is a *positive* predictor of hypermethylation as a *differential* DiffBind fold (genes that *gain* K119ub hypermethylate; OR ≈ 4.71) but a *negative* predictor as *baseline constitutive* signal in the TET-impediment model (β ≈ −1.05; high-resting-K119ub loci are already compacted/inaccessible). These are compatible — gaining K119ub drives the switch; constitutively high K119ub marks already-repressed regions.
- **Coordinated = 78.7%, not 92.3%.** The canonical run-5 figure for coordinated mC↑/hmC↓ among co-significant genes is 78.7% (6,589/8,371). The "92.3%" string in several CLAUDE.md files is stale/unverified against run-5 tables — do not cite it.
- **n = 4 effective replicates/condition** (8 samples, 2 ctrl + 2 mut per sex; sex covariate in run-5, sex/genotype partially confounded at n=2/group per sex) — bounds the precision of all effect-size estimates.
