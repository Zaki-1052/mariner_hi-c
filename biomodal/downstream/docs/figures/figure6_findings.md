# figure6: Neuronal Genes Are Preferentially Affected

**Paper mapping:** Results section R5 (neuronal-gene specificity / Polycomb de-repression).
**Script:** `scripts/figures/figure6_neuronal_specificity.R`
**Output:** `plots/figures/figure6_neuronal_specificity/` (PDF + SVG + JPG) plus standalone panel C `plots/figures/figure6C_synapse_chromatin/`.
**Layout:** 4 panels, `design = "AABB\nCCDD"`, 180 × 160 mm.

---

## One-line takeaway

Neuronal-identity genes are *preferentially* (not exclusively) affected by BAP1 loss because they sit in constitutively H2AK119ub-rich Polycomb chromatin — but the MeCP2 redistribution that follows tracks the methylation shift (~3× more strongly than neuronal lineage), the synapse/axon subset undergoes selective Polycomb de-repression, and neuronal methylation changes are stoichiometric 5hmC→5mC exchange.

---

## Biological story

This figure positions the neuronal phenotype within the paper's central cascade: BAP1 loss → H2AK119ub accumulation → coordinated 5mC↑/5hmC↓ at gene bodies → (TET-demethylation block; MeCP2/chromatin readout) → disrupted neuronal gene regulation. The question it answers is *why neuronal genes specifically*. The answer is a chromatin-substrate argument, not a lineage-specificity argument. Neuronal-identity and axon-guidance genes are constitutively decorated with H2AK119ub in wildtype cerebellum (Polycomb-poised developmental loci), so they are the natural substrate for the deubiquitinase BAP1 normally polices. When BAP1 is lost, the K119ub that accumulates lands disproportionately on these already-marked neuronal genes (Panel A), making them over-represented among the most-affected loci.

The figure then deliberately bounds that claim. The MeCP2 redistribution downstream of the methylation change overlaps the coordinated mC↑/hmC↓ gene set far more than it overlaps generic neuronal identity (Panel B: OR 5.16 vs 1.73; neuronal × coordinated is not even significant), so MeCP2 follows the *methylation/chromatin state*, not cell-type lineage — this is why the paper softens R5's original "exclusively at neuronal genes" to "preferentially/disproportionately." Within the neuronal compartment, the specificity is sharpened further: the synapse/axon subset loses extra H3K27me3 (selective Polycomb de-repression) without any extra accessibility or enhancer-mark gain (Panel C), refining the mechanism toward de-repression of poised synaptic promoters rather than blanket heterochromatinization. Finally, neuronal and synaptic loci show a δ-5mC vs δ-5hmC slope of ≈ −1 (Panel D), the stoichiometric, conserved-total signature of direct 5hmC→5mC conversion (dehydroxymethylase-like DNMT3A activity), distinguishing them from the genome-wide deviation from −1 that is driven by non-neuronal genes (where TET inhibition plus continued de novo methylation prevails).

---

## Per-panel findings

| Panel | What it shows | Key numbers (verified against source TSVs) | Reader's conclusion |
|-------|---------------|--------------------------------------------|---------------------|
| **A** (Sec 72) | Neuronal-gene fraction across deciles of *constitutive* (wildtype) H2AK119ub gene-body signal, with binomial CI ribbon and a genome-wide reference line. `72_neuronal_decile_summary.tsv` + `72_fisher_results.tsv`. | Top decile (D10) neuronal fraction = **33.0%** (713/2,161), Fisher **OR = 1.70** [1.54, 1.87], q = 3.4e-25, vs genome-wide baseline **23.5%** (5,077/21,604). Bottom decile D1 = 20.9% (OR = 0.84). Trend overall increasing (non-strict; D6 dips to 19.5%). | Neuronal genes are constitutively K119ub-enriched, and the enrichment scales with signal — they are Polycomb-poised loci, the natural BAP1 substrate. |
| **B** (Sec 74) | Pairwise Fisher odds ratios on a log-y axis for three gene-set overlaps, colored to highlight the methylation-driven one. `74_pairwise_fisher.tsv`. | **Coordinated × MeCP2-Up OR = 5.16** [3.19, 8.51], p_adj = 2.8e-12 (overlap 51 genes). **Neuronal × MeCP2-Up OR = 1.73** [1.05, 2.82], p_adj = 0.034 (overlap 27). **Neuronal × Coordinated OR = 1.05** [0.98, 1.13], p_adj = 0.141 (**NS**). | MeCP2 redistribution tracks the methylation shift, not neuronal lineage — ~3× stronger. This is the evidence that softens "exclusively" to "preferentially." |
| **C** (Sec 76) | Gene-body DiffBind differential signal (log2 fold) for ATAC, H3K27ac, H3K27me3, faceted by mark and split into synapse/axon vs broader-neuronal vs non-neuronal classes (synapse set is a strict subset of the broader neuronal set). Class partition independently reproduces the source table's n's and medians. `diffbind_gene_level_all_marks.tsv` + `76_synapse_axon_gene_set.tsv` + `72_neuronal_gene_set_go_derived.tsv`; stat from `76_synapse_vs_neuronal_stats.tsv`. | Synapse/axon vs broader-neuronal **H3K27me3 Δmedian = −0.044** (−0.132 vs −0.087), Wilcoxon **p = 3.0e-3** (p_adj = 6.6e-3, significant). Same contrast: ATAC p = 0.65 (NS), H3K27ac p = 0.76 (NS). Broader-neuronal vs non-neuronal K27me3 p = 0.289 (NS) — the K27me3-loss specificity lives in the synapse/axon subset. | Synapse/axon genes undergo *selective* Polycomb de-repression (extra K27me3 loss only), not blanket activation — de-repression of poised synaptic promoters. |
| **D** (Sec 78) | OLS slope of δ-5mC vs δ-5hmC per gene class with 95% CI, against a dashed reference at −1; bars colored by whether the CI overlaps −1 (stoichiometric) or deviates. `78_stoichiometry_slopes.tsv`. | **Neuronal (broad) slope = −0.995** [−1.039, −0.952], CI overlaps −1 (**stoichiometric**). **Synapse/axon = −1.020** [−1.078, −0.962] (**stoichiometric**). **All genes = −0.959** [−0.978, −0.940] (**deviates**). **Non-neuronal = −0.949** [−0.970, −0.928] (**deviates**). (Context: Coordinated = −1.29, steeper, net mC gain.) | Neuronal/synaptic methylation changes are stoichiometric mC-for-hmC exchange (direct 5hmC→5mC, dehydroxymethylase-like), distinct from the genome-wide deviation driven by non-neuronal loci. |

---

## Caveats & framing

- **"Preferentially," not "exclusively."** The paper outline (R5) originally phrased this as "Exclusively at genes involved in neuronal structure and development." Panel B falsifies the strong claim: MeCP2 redistribution overlaps the coordinated methylation set (OR = 5.16) far more than neuronal identity (OR = 1.73), and the neuronal × coordinated overlap is not significant (OR = 1.05). The figure and any caption must use "preferentially / disproportionately." The mechanism is that neuronal genes are *enriched* in the affected (K119ub-rich, methylation-shifted) compartment, not that the effect is restricted to them.

- **MeCP2 is CUT&RUN, not ChIP.** All MeCP2 data feeding Panel B (and the MeCP2-Up gene set) is CUT&RUN. Histone marks in Panel C are DiffBind differential signal. The script labels marks as "BAP1-KO differential signal (log2 fold change)" and never says "ChIP" for MeCP2 — preserve this in captions.

- **K119ub sign reconciliation (two roles, not a contradiction).** Panel A uses *absolute / constitutive* wildtype K119ub signal to define which genes are Polycomb-poised — a baseline-state property. Elsewhere in the paper K119ub appears as a *positive differential* predictor of hypermethylation (genes that *gain* K119ub hypermethylate, e.g. Figure 2 logistic OR ≈ 4.7), while in the TET-impediment model high *constitutive* K119ub is a negative predictor of methylation change (already-compacted, inaccessible loci). Figure 6 Panel A is the baseline-state usage: neuronal genes are constitutively K119ub-high, which is what makes them the disproportionate substrate. These are compatible because they describe different quantities (constitutive level vs differential gain).

- **De-repression, not heterochromatinization.** Section 73 (the K119ub × neuronal interaction analysis underlying this figure's biology) found that K119ub-high neuronal genes *gain* accessibility and H3K27ac under BAP1 loss — the opposite of a naive compaction/heterochromatin-shift prediction. Panel C is consistent with this: the selective signal at synaptic genes is K27me3 *loss* (Polycomb de-repression), with no extra ATAC/K27ac gain beyond other neuronal genes. Frame as de-repression of poised neuronal promoters.

- **Stoichiometry self-correction context (Sec 78).** The broad GO-derived neuronal set (5,614 genes; 4,118 with marks data) is used deliberately to avoid the circularity of the section-61 narrow DMR-selected set (1,149 genes). With the unbiased set, total methylation at neuronal genes slightly *decreases* (mean Δ −0.0022), and the slope is stoichiometric (−1.0); the narrow set spuriously showed a total *increase* (+0.012). Panel D presents the unbiased result as the finding, not as a "correction."

- **Power / replicate caveats.** n = 8 samples (2 ctrl + 2 mut per sex, sex covariate in run-5; sex/genotype partially confounded at n=2/group per sex). Some neighboring Section 76/78 cells are small-n (e.g. the 16 triple-overlap genes, the 72–79-gene MeCP2-Up quadrant); those are not the panels in this figure, but the broad-set slopes in Panel D (neuronal n=4,118; synapse n=2,099) are well-powered while the MeCP2-Up slope shown elsewhere (n=79) is not — Panel D uses the well-powered classes.

- **Section 73 table-prefix caveat (provenance only).** The chromatin-remodeling tables that motivate Panel C's biology are on disk under a stale `72g_*` prefix (the committed script writes `73_*` but was not re-run after a rename). The *data* is current and valid; Panel C itself reads `diffbind_gene_level_all_marks.tsv` and `76_*`/`72_*` gene-set files, which are unaffected.

---

## Verification verdict

PASS (no fixes required). All four panels are fully implemented and composed via patchwork. Every annotated statistic was confirmed by reading the raw TSVs: D10 OR = 1.70 and genome baseline 23.5% (5,077/21,604); Coordinated×MeCP2-Up OR = 5.16, Neuronal×MeCP2-Up OR = 1.73, p_adj 0.034; synapse-vs-broader-neuronal K27me3 Δ = −0.044, Wilcoxon p = 3.0e-3; slopes −0.96 / −0.95 / −1.00 (neuronal) / −1.02 (synapse) with `differs_from_neg1` correctly mapped (deviates / deviates / stoichiometric / stoichiometric). The Panel C gene-class partition was independently reproduced in Python and matches the source table's n's and medians exactly for all three marks. Sourcing order, color/theme consistency, output formats (PDF + SVG + JPG, no PNG), balanced delimiters, and loud-failure (stop) guards all pass; no forbidden patterns (no "ChIP", no `run-N` paths, no `--vanilla`, no fallbacks/tryCatch). One trivial cosmetic note: the on-panel K27me3 p-value formats as `3.0e-03` (`%.1e`) while prose may write `2.95e-3`; both faithfully represent the source `wilcox_p = 0.002954`.
