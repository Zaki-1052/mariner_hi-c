# figure2: H2AK119ub Drives Hypermethylation at Active Chromatin

**Paper mapping:** Results section R2 (the mechanism figure — "Increase in H2AK119Ub after Bap1 loss leads to a global 5mC increase and 5hmC decrease trend at gene bodies").
**Script:** `scripts/figures/figure2_k119ub_chromatin_geography.R`
**Output:** `plots/figures/figure2_k119ub_chromatin_geography/` (PDF + SVG + JPG, 180 × 200 mm)
**Verification verdict:** PASS (every annotated statistic matches the source TSVs exactly; one cleanup fix applied — removed an unused `format_p()` helper; one non-blocking WARN — Panel B's AUC=0.818 is a run-log value, not table-traceable).

---

## One-line takeaway

H2AK119ub gain is the single strongest predictor of where BAP1-KO hypermethylation lands, and that hypermethylation is restricted to **active, euchromatic (A-compartment) chromatin** while Polycomb heterochromatin is actively **protected** — establishing the BAP1 substrate as the upstream driver of the methylation phenotype.

---

## Biological story

Figure 1 establishes *what* happens in BAP1-KO cerebellum — a coordinated, reciprocal 5mC↑/5hmC↓ shift concentrated at gene bodies. Figure 2 is the *why*: it pins the cause on H2AK119ub, the histone mark whose deubiquitinase is BAP1 itself. When BAP1 is lost, H2AK119ub accumulates, and this figure shows that the resulting methylation gain is not random. A quantitative four-mark logistic model ranks H2AK119ub fold-change as the dominant positive predictor of hypermethylation (OR ≈ 4.7 per unit log2FC), an order of magnitude above any other mark, with the active/open marks (ATAC, H3K27ac) coming out protective (OR < 1). Because BAP1 *is* the H2AK119ub eraser, its loss raising K119ub being the top predictor places K119ub upstream in the causal chain rather than as a passenger.

The figure then shows *where* this lands, and the geography is the mechanistic punchline. Hypermethylation is overwhelmingly an **active-chromatin** event — Active_Promoter DMRs are 93% hypermethylated, the A (euchromatin) compartment is enriched ~14-fold for hypermethylation, and the same loci are losing 5hmC. By the mirror logic, Polycomb/H3K27me3-repressed chromatin is excluded from hypermethylation (OR ≈ 0.063) and instead drifts hypomethylated. This is the decisive falsification of the naive "K119ub marks Polycomb, so Polycomb genes get hypermethylated" expectation: compact heterochromatin is inaccessible to DNMT3A, so the de-novo/retained methylation burden falls on normally active genes where TET is normally most active and 5hmC normally accumulates. This is exactly the signature of a **TET-demethylation block** at active loci (developed mechanistically in Figure 3) — and it geographically mirrors the Lopez-Moyado TET-KO DNMT3A-redistribution phenotype. The final panel keeps the K119ub claim honest: the coupling is real but partial (most DMR genes carry no K119ub peak at all), so K119ub is the upstream effector and elevated broadly, not a locus-by-locus gatekeeper. The downstream consequence of this active-chromatin methylation shift — MeCP2 redistribution at neuronal genes — is the subject of later figures.

---

## Per-panel findings

| Panel | What it shows | Key numbers (verified) | Reader's conclusion |
|-------|---------------|------------------------|---------------------|
| **A** (Sec 10) | Stacked fraction bar: per chromatin state, the split of significant mC DMRs into hyper vs hypo. | **Active_Promoter 93.0% hyper** (4,562 hyper / 4,906 total). Repressed_Promoter 94.4% hypo (1,621 / 1,718). Mirror-image partitioning across the 7-state system. | Hypermethylation lands at *active* promoters/enhancers; hypomethylation at *repressed* promoters. The direction split is the signal — pooling would cancel it. |
| **B** (Sec 33) | Forest plot: 4-mark logistic-regression odds ratios (per unit DiffBind log2FC) predicting mC hypermethylation. | **H2AK119ub OR = 4.71** [3.33–6.66], p = 2.4e-18 (dominant). H3K27me3 OR = 1.44. H3K27ac OR = 0.48 (protective). ATAC-seq OR = 0.22 (protective). Model AUC = 0.818 (run-log context, see caveats). | K119ub gain is the strongest quantitative driver of hypermethylation, an order of magnitude above K27me3; open/active marks protect. Places K119ub upstream. |
| **C** (Sec 29) | A/B compartment enrichment OR bars (log x-axis): where each methylation change concentrates. | **mC hyper → A OR = 13.6** [11.92–15.67], p≈0. hmC loss → A OR = 9.8. mC hypo → B OR = 1.7. hmC gain → B OR = 3.0. | Hypermethylation and 5hmC loss are euchromatin (A-compartment) events — the TET-KO DNMT3A-redistribution geography. |
| **D** (Sec 30) | Polycomb exclusion/enrichment bars (log x-axis): the four methylation directions vs chromatin-state Polycomb membership. | **Polycomb × mC hyper OR = 0.063** [0.053–0.075], p≈0 (excluded). mC hypo OR = 9.80 (enriched). hmC gain OR = 11.39 (enriched). hmC loss OR = 0.131 (excluded). | Polycomb heterochromatin is *protected* from hypermethylation and instead drifts hypomethylated — falsifies the naive "Polycomb gets hypermethylated" prediction. |
| **E** (Sec 17) | The "honest" K119ub breakdown shown as data: % of DMR genes with no K119ub peak, and conditional gain rate among peak-bearing genes vs background. | **58% of DMR genes carry no K119ub peak** (58.16%). Among peak-bearing genes, mC-Up gains K119ub at 47.8% vs 33.6% background = **+14.2 pp**. | The K119ub coupling is real but partial — a modest enrichment in a minority of genes, atop a diffuse global K119ub rise. Pre-empts reviewer skepticism without hiding the limitation. |

---

## Caveats & framing

- **"Preferentially," not "exclusively."** The active-chromatin restriction is strong but not absolute. Hypermethylation is *concentrated* at active promoters/enhancers and the A compartment, and Polycomb is *protected* — but Unmarked regions also carry substantial DMRs (2,744 hyper). Frame as "preferentially / disproportionately at active chromatin," matching the project-wide softening of absolute claims.

- **K119ub coupling is partial (Panel E is load-bearing for honesty).** Most DMR genes (58%) have no K119ub peak at all, and the gene-specific enrichment is a modest +14.2 pp over background. The mechanism is best stated as a *global* K119ub increase that is upstream and elevated genome-wide (Section 18: background median K119ub log2FC > 0, p < 1e-19; gene-specific Cliff's delta only 0.10–0.18, "negligible-to-small"), rather than sharp locus-specific redistribution. K119ub is the upstream effector, not a per-locus gatekeeper.

- **K119ub sign reconciliation (important for the manuscript).** K119ub appears with *opposite* signs depending on whether you measure it as a **differential fold-change** or as **baseline absolute signal** — these are compatible, not contradictory:
  - As a *differential* DiffBind log2FC (this figure, Panel B; also Sections 14/33/41), genes that *gain* K119ub hypermethylate → **positive predictor** (OR = 4.71).
  - As *baseline constitutive* signal (the TET-impediment model in Figure 3 / Section 24, β ≈ −1.05), genes with high pre-existing K119ub are already compacted and inaccessible to DNMT3A → **negative predictor**.
  - The reconciling statement: *gaining* K119ub drives hypermethylation in active chromatin, whereas *constitutively high* K119ub marks already-repressed regions that are protected. The manuscript should state this explicitly to avoid an apparent conflict between figures.

- **A/B and Polycomb numbers are run-5 canonical; older docs carry stale values.** TODO.md / FIGURES.md cite mC-hyper→A "OR=14.71" and Polycomb-hyper "OR=0.064 / hypo OR=8.71" — these are prior-run numbers. The verified run-5 TSV values used in this figure are **A OR = 13.6** and **Polycomb hyper OR = 0.063 / hypo OR = 9.80**. Use the figure's values.

- **Panel B AUC = 0.818 is run-log context, not table-traceable.** `diffbind_logistic_model_coefficients.tsv` stores only the odds ratios; the AUC = 0.818 (4-mark model) lives in the run log / Section 32-37 doc, not in any figure-source TSV. It is a true documented value but a reader cannot reproduce it from the panel's table — the caption should cite the run-log source. (Adding baseline mC/hmC raises the extended-model AUC to 0.869.)

- **All H2AK119ub references are the histone mark / BAP1 PR-DUB substrate.** This is a ChIP-derived mark, correctly labeled "mark"/"signal." (MeCP2, which appears in adjacent figures, is **CUT&RUN, not ChIP** — never label MeCP2 as "ChIP.") This figure does not use MeCP2.

- **The two "protective" marks in Panel B are quantitative, not just qualitative.** ATAC (OR=0.22) and H3K27ac (OR=0.48) both protect against hypermethylation, consistent with the broader chromatin profile (Section 19): it is *hypomethylated* genes that gain H3K27ac, not hypermethylated genes that lose it. Do not over-state H3K27ac loss as a hypermethylation correlate.
