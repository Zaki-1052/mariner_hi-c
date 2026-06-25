# figure7: Mechanism Summary Data Panel (cascade)

**Paper mapping:** All R-sections (cascade summary). Data companion to the Illustrator model schematic.
**Script:** `scripts/figures/figure7_mechanism_summary.R`
**Verification verdict:** PASS — production-ready, every annotated statistic matches its source TSV exactly; the one malformed input (`62_model_comparison_summary.tsv`) is handled with a strict de-spilling guard, no values fabricated.

---

## One-line takeaway

Each panel quantifies one link of the BAP1 → H2AK119ub → TET-block → coordinated 5mC↑/5hmC↓ → MeCP2 → neuronal-gene cascade, so the whole mechanistic chain reads as a single, internally consistent argument rather than five unrelated tests.

---

## Biological story

Figure 7 is the quantitative backbone for the model schematic. It does not introduce new analysis; it lifts the single strongest number from each step of the mechanism and lines them up left-to-right so the reader can follow the causal chain as one continuous argument. The chain is: BAP1 loss removes the deubiquitinase that normally erases H2AK119ub, so K119ub accumulates genome-wide (Panel A); that K119ub gain is the dominant chromatin predictor of which loci hypermethylate (Panel B); the hypermethylation arises because TET-mediated demethylation is impeded — genes with more 5hmC substrate change most, and a TET-impediment model out-discriminates a DNMT3A-recruitment model decisively (Panel C); the resulting methylation shift is read out by MeCP2, but MeCP2 binding is explained far better by chromatin state (K119ub foremost) than by CG methylation itself (Panel D); and all of these layers converge disproportionately — but not exclusively — on constitutively K119ub-marked neuronal genes (Panel E).

The figure is deliberately built to "show the data" at each step (effect-size bars, odds ratios with CIs, AUCs with CIs, R² bars, and an OR forest) rather than reproducing model-internal artifacts. It is the cascade in numbers: a global, diffuse upstream signal (K119ub) that becomes a sharp, direction-specific methylation phenotype at active chromatin, is read by MeCP2 as a chromatin state rather than a methylation level, and lands preferentially on the neuronal Polycomb genes that are BAP1's natural substrate.

---

## Per-panel findings

| Panel | What it shows | Key numbers (verified) | Reader's conclusion |
|-------|---------------|------------------------|---------------------|
| **A — BAP1 → H2AK119ub** (effect-size bar; Sec 18, `k119ub_bigwig_signal_summary.tsv`) | Per-group median K119ub BigWig log2FC (mut/ctrl). Background ("All DMR genes") anchors the global-increase claim; mC-Up / mC-Down directional bars map onto the methylation story. | Background median log2FC = **+0.007** (0.006956), p vs 0 = **1.81e-20**. mC-Up = **+0.059** (0.0586, n=6,993). mC-Down = **−0.080** (−0.0797, n=2,974). | BAP1 loss raises H2AK119ub globally (significant background shift), and the rise is directionally correct at DMRs — hypermethylated genes shift positive, hypomethylated negative. The upstream signal is real but diffuse (effect sizes small; see caveats). |
| **B — H2AK119ub → hypermethylation** (OR waterfall; Sec 33, `diffbind_logistic_model_coefficients.tsv`) | Logistic model of hypermethylation on four differential DiffBind marks, ordered by odds ratio. | **H2AK119ub OR = 4.71** [3.33, 6.66], p = 2.4e-18 (dominant). H3K27me3 OR = 1.44 [1.19, 1.75] (also >1). ATAC OR = 0.22; H3K27ac OR = 0.48 (both <1, suppress hyper). | Of the four differential chromatin marks, a *gain* in K119ub is by far the strongest positive driver of hypermethylation — the quantitative center of the "K119ub drives the methylation switch" claim. |
| **C — TET impediment, not DNMT3A recruitment** (AUC comparison; Sec 23+24, `dnmt3a_model_comparison.tsv` + `baseline_hmc_model_comparison.tsv` + `dnmt3a_exclusive_model_comparison.tsv`) | Predictive AUC of three models for hypermethylation: baseline-5hmC substrate (TET pathway support), TET-impediment, and DNMT3A-recruitment. | Baseline 5hmC AUC = **0.762** [0.755, 0.769]. TET impediment AUC = **0.793** [0.785, 0.801]. DNMT3A recruitment AUC = **0.696** [0.687, 0.706]. DeLong p = **9.43e-49** (TET > DNMT3A). | The methylation gain is best explained by impaired TET-mediated demethylation (substrate-limited: high-5hmC genes change most), not by K119ub-driven de-novo DNMT3A recruitment. The model adjudication is decisive. |
| **D — DNA methylation → MeCP2 binding** (R² cascade; Sec 62, `62_model_comparison_summary.tsv`) | Variance in MeCP2 binding explained by CG methylation alone vs chromatin marks alone vs full model, at binding level and differential level (n = 202,574 peaks). | Binding-level chromatin-only R² = **0.246** vs CG-only R² = **0.017** (~15×). (Differential level: chromatin 0.170 vs CG 0.013.) | MeCP2 reads chromatin state, not CG methylation per se — methylation is necessary but not sufficient. K119ub is the top single predictor of MeCP2 binding (Sec 62 standardized β = +0.199), tying MeCP2 redistribution back to the same upstream signal. |
| **E — Convergence on neuronal genes** (OR forest, log scale; Sec 74+72, `74_pairwise_fisher.tsv` + `72_fisher_results.tsv`) | Three odds ratios: methylation-driven MeCP2 redistribution, lineage-driven MeCP2 redistribution, and constitutive K119ub enrichment at neuronal genes. | Coordinated mC↑/hmC↓ × MeCP2-Up OR = **5.16** [3.19, 8.51], adj p = 2.8e-12. Broad-neuronal × MeCP2-Up OR = **1.73** [1.05, 2.82], adj p = 0.034. Constitutive K119ub top-decile × neuronal OR = **1.70** [1.54, 1.87], q = 1.3e-25. | MeCP2 redistribution tracks the *methylation* phenotype ~3× more strongly than generic neuronal identity (5.16 vs 1.73). Neuronal genes are affected because they are constitutively K119ub-rich (the natural BAP1 substrate, OR 1.70) — preferentially, not exclusively. |

---

## Caveats & framing

- **"Preferentially," not "exclusively."** Panel E is the explicit evidence that neuronal involvement is disproportionate but not the primary axis: MeCP2 redistribution tracks coordinated methylation (OR 5.16) about three times more strongly than neuronal identity (OR 1.73), and Neuronal × Coordinated is not even significant (OR 1.05, Sec 74). Neuronal genes are over-represented because they sit in constitutive K119ub/Polycomb chromatin (OR 1.70, Sec 72), which is the BAP1 substrate — not because the mechanism is neuron-specific. Any text built on this figure should say "preferentially / disproportionately."

- **MeCP2 is CUT&RUN, not an antibody-pulldown ChIP assay.** Panel D and Panel E labels use "binding"/"signal." Do not relabel as ChIP.

- **K119ub sign reconciliation (the key non-obvious point).** K119ub appears with *opposite* signs in two panels, and this is not a contradiction:
  - In Panel B (Sec 33) K119ub is a strong **positive** predictor (OR = 4.71) because the predictor is the *differential* DiffBind fold — genes that *gain* K119ub in the mutant hypermethylate.
  - In the TET adjudication behind Panel C (Sec 24) K119ub enters as a **negative** predictor (standardized β = −1.05) because there the predictor is *baseline absolute* K119ub signal — genes with high *constitutive* K119ub are already compacted/inaccessible and are less available to the de-novo machinery.
  - The reconciliation: gain of K119ub drives hypermethylation; high resting K119ub marks already-repressed loci that change less. Panels B and C therefore use different K119ub quantities and are compatible.

- **Panel A effect sizes are small (diffuse, not locus-sharp).** The global K119ub increase is statistically robust (background p = 1.81e-20) but the per-gene magnitudes are modest (Cliff's delta ~0.10–0.18, "negligible-to-small"; hmC-Down genes are indistinguishable from background, p = 0.30). The honest reading (Sec 17/18) is a genome-wide diffuse K119ub rise rather than sharp DMR-specific redistribution — Panel A should be cited as evidence of a *global* upstream shift, with the locus-specific selectivity coming from Panels B–C.

- **Panel C model comparison is cross-validated and robust.** TET vs DNMT3A is not an overfit artifact: 10-fold CV AUCs are near-identical to in-sample (optimism ≈ 0.0003), the result persists in the non-promoter stratum and across all 5hmC tertiles, and the exclusive-feature contrast is even more extreme (DeLong p = 6.13e-75). The DeLong p stored in the TSV (9.43e-49) is the exact value; the script display floor "p < 2.2e-16" used elsewhere is only a formatting floor.

- **Panel D source file is malformed but handled correctly.** `62_model_comparison_summary.tsv` was written unquoted while its `label` column embeds a literal newline, so each real record spills into a ragged second line. The script reads it with `fill = TRUE`, then selects the six real records by exact `type × model` identity and enforces `nrow == 6` — no values are guessed or imputed. Verified binding-level values: CG-only r_squared = 0.016795 (≈0.017), Chromatin-only = 0.245650 (≈0.246).

- **Writer self-report inaccuracy (does not affect the figure).** The script's own header/comment text characterizes the non-K119ub marks loosely; in fact H3K27me3 has OR = 1.44 (>1, also promotes hyper), while ATAC (0.22) and H3K27ac (0.48) are <1. The plotted Panel B colors marks correctly by `or > 1` vs `or < 1`, so the rendered figure is data-faithful — only the prose comment is imprecise.

- **n = 4 effective replicates.** All numbers are run-5 (deep-seq, 8 samples = 2 ctrl + 2 mut per sex, sex covariate included). Sex and genotype remain partially confounded at n=2/group per sex.
