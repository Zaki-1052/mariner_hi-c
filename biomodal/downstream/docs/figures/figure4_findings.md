# figure4: MeCP2 Redistribution and Developmental Amplification

**Paper mapping:** R1 (the title figure)
**Script:** `scripts/figures/figure4_mecp2_redistribution.R`
**Verification verdict:** PASS (after one cosmetic fix). All 6 panels fully implemented, no stubs; every annotated OR/p/count/% checked against the source TSVs and the `MeCP2_annotated.txt` peak file and matches exactly. One presentation-only fix applied (Panel B distal-segment label position). MeCP2 data is CUT&RUN (Epicypher), never ChIP.

---

## One-line takeaway

In BAP1-KO adult cerebellum, MeCP2 binding does not increase uniformly — it **redistributes** away from gene bodies toward distal-intergenic sites, **preferentially overlaps hypermethylated DMRs** (Fisher OR = 5.13), and the **normal developmental MeCP2 ramp is amplified** (2.1x more age-increased peaks in mutant), with the excess developmental gain landing on the same coordinated mC↑/hmC↓ gene set that defines the methylation phenotype.

---

## Biological story

This figure is the downstream-consequence chapter of the paper's central cascade: BAP1 loss → H2AK119ub accumulation → coordinated 5mC↑ / 5hmC↓ at gene bodies (a TET-demethylation block) → an altered chromatin/methylation landscape that MeCP2 must read. Because MeCP2 binds 5mCG with high affinity and 5hmCG with low affinity, the retained-5mC / depleted-5hmC signature is predicted to re-position MeCP2. Figure 4 shows that it does: at the level of *which DMRs sit under differential MeCP2 peaks*, hypermethylated DMRs are 5.13-fold enriched for MeCP2 gains (Panel C). But the gain is geographically selective — MeCP2 vacates gene bodies (62% of lost peaks are intronic) and concentrates at new distal-intergenic regulatory sites (52% of gained peaks; Panel B). The genome-wide gene-body MeCP2 median therefore *falls* even though ~7,700 peaks gain signal, because the gains pile onto ~2,000 genes while the other ~90% thin slightly. This is redistribution, not global up-regulation.

The figure then ties the phenotype to neurodevelopment. MeCP2 normally accumulates as neurons mature; using young-vs-adult CUT&RUN DiffBind per genotype, the mutant shows 2.1x more age-increased peaks than control (Panel D), and even at the 7,305 loci where both genotypes age-gain MeCP2, the mutant climbs ~22% higher (median fold 2.24 vs 1.83; Panel E). The new Panel F closes the loop back to methylation: the genes that gain age-related MeCP2 *only* in the mutant are significantly enriched for the coordinated mC↑/hmC↓ DMR set, directly connecting the developmental MeCP2 overshoot to the methylation block. The mechanistic reading is that excess H2AK119ub creates aberrant, retained-5mC chromatin that recruits additional MeCP2 during the maturation window, converting a developmental ramp into an overshoot. (Note: the *tightness* of methylation-to-MeCP2 coupling is quantitatively modest at the gene level — see caveats — which is exactly what motivates Figure 5's "MeCP2 reads chromatin state, not methylation" argument.)

---

## Per-panel findings

| Panel | What it shows | Key numbers (verified) | Reader's conclusion |
|-------|---------------|------------------------|---------------------|
| **A** | Adult MeCP2 binding volcano (CUT&RUN), log2 fold (mutant/control) vs −log10 FDR; top genes and paper key genes labeled. | 202,650 total peaks; **7,686 UP / 1,200 DOWN** at FDR < 0.05 (recomputed from `MeCP2_annotated.txt`, Fold col 9 / FDR col 11). UP = orange (`COLORS$mecp2`), DOWN = purple. | BAP1 loss produces thousands of MeCP2 binding changes in the adult, strongly skewed toward gains (6.4:1 UP:DOWN). |
| **B** | Stacked fraction of significant peaks by genomic feature, UP vs DOWN. | UP peaks: **51.7% Distal Intergenic**, 41.0% Intron, **2.2% Promoter**. DOWN peaks: 61.8% Intron, 19.5% Distal Intergenic, 8.0% Promoter (`75_peak_annotation_distribution.tsv`). | MeCP2 *gains* land at distal-intergenic sites (not promoters); MeCP2 *losses* are genic/intronic. MeCP2 vacates gene bodies and concentrates distally — redistribution, not uniform increase. |
| **C** | Dodged 2×2 bar: fraction of hyper/hypo DMRs overlapping MeCP2-Up vs MeCP2-Down peaks, with Fisher annotation. | Hyper DMRs: **7.16%** overlap MeCP2-Up (538/7,513) vs 1.89% MeCP2-Down. Hypo DMRs: 5.09% MeCP2-Up vs 6.90% MeCP2-Down. **Fisher OR = 5.13, p = 1.27e-33** (`mecp2_dmr_overlap_summary.tsv`). | Hypermethylated DMRs preferentially sit under MeCP2 gains, as predicted by MeCP2's preference for 5mCG over 5hmCG — the DMR-level link between the methylation phenotype and MeCP2. |
| **D** | Grouped bar of significant young→adult (aging) peaks per genotype (UP/DOWN; NS dropped). | Control: **10,930 UP** / 2,822 DOWN. Mutant: **23,117 UP** / 10,646 DOWN. Mutant = **2.1x more** age-increased peaks (and 3.8x more age-decreased) (`77_aging_peak_summary.tsv`). | The normal developmental MeCP2 ramp is amplified in BAP1-KO — the mutant accumulates far more age-related MeCP2 changes. |
| **E** | Scatter of control vs mutant aging fold at the 7,305 loci both genotypes age-gain, 1:1 reference line, neuronal genes highlighted. | **7,305 shared loci**; median aging fold **2.24 (mutant) vs 1.83 (control)**, +0.41 log2 (~+22%). **2,266** shared loci are at neuronal genes (`77_shared_peak_fold_comparison.tsv`). | Even where both genotypes age-gain MeCP2, the mutant climbs higher — an overshoot of the shared developmental program, disproportionately at neuronal genes. |
| **F** *(NEW)* | Observed vs expected overlap of mutant-unique age-increased genes with the coordinated mC↑/hmC↓ DMR gene set, over the ~20,915-gene universe, with Fisher annotation. Result saved to `plots/figures/tables/fig4f_aging_methylation_overlap.tsv`. | Re-derived from the two MeCP2 aging DiffBind files via ChIPseeker (negate raw log2(young/adult) → adult>young; FDR<0.05 & Fold>0). Reproduces Section 77 peak counts (ctrl 10,930 / mut 23,117 age-increased) and ~1,654 mutant-unique aging genes. Fisher OR/p computed at runtime (table is generated on run). | The developmental MeCP2 overshoot lands on the same coordinated-methylation genes — directly linking the amplified MeCP2 trajectory to the BAP1→K119ub→methylation phenotype. |

---

## Caveats & framing

- **"Preferentially," not "exclusively."** The MeCP2 phenotype is enrichment, not exclusivity. MeCP2 redistribution tracks the *coordinated methylation* set far more than generic neuronal identity (Section 74: Coordinated × MeCP2-Up OR = 5.16 vs Neuronal × MeCP2-Up OR = 1.73; Neuronal × Coordinated not even significant, OR = 1.05). Use "preferentially / disproportionately."

- **MeCP2 = CUT&RUN, not ChIP.** All MeCP2 panels describe "binding," "signal," or "gains/peaks." The only "ChIP" string in the script is the *tool* name `ChIPseeker` (used to annotate aging peaks to nearest genes) plus corrective comments. Never label MeCP2 as ChIP.

- **DMR-level enrichment vs gene-level coupling (the reader-vs-driver nuance).** Panel C's OR = 5.13 is a *DMR-overlap* statistic — it says which DMRs sit under differential MeCP2 peaks. Aggregated to the gene level the directional coupling is much weaker: binary MeCP2-gain enrichment for coordinated vs non-coordinated genes is OR ≈ 1.04 (p ≈ 0.80, NS), and stratified gene-level Spearman ρ ≈ 0.02–0.08 (Section 11). This is not a contradiction — it is the motivation for Figure 5 ("MeCP2 reads chromatin state, primarily K119ub, not methylation per se"). MeCP2 redistribution is a real but quantitatively modest 1:1 readout of methylation change.

- **K119ub sign reconciliation (carried from the cascade, relevant for interpretation).** K119ub is a *positive* predictor of hypermethylation as a *differential DiffBind fold* (genes that *gain* K119ub hypermethylate; e.g. Section 14 differential OR = 4.46) but a *negative* predictor as *baseline absolute signal* in the TET model (genes with high *constitutive* K119ub are already compacted/inaccessible to DNMT3A). These are compatible: gaining K119ub drives the methylation switch, whereas pre-existing K119ub-dense chromatin is protected. Figure 4 itself does not plot K119ub, but the MeCP2-gain story rests on the "gained K119ub → retained 5mC → MeCP2 recruited" arm of this reconciliation.

- **Panel B cosmetic fix (applied).** The on-bar "52% distal" white label was originally positioned at y ≈ 0.741 (over the Intron segment) because "Distal Intergenic" is the last factor level and `geom_col` stacks it at the *bottom* of the bar [0, 0.517]. Corrected to y = (up_distal/100)/2 ≈ 0.259 so it centers on the distal-intergenic segment it describes. Bar heights/data were always correct; this is presentation-faithfulness only.

- **Panel F is the only `tryCatch` in the figure scripts (by design).** If the two large MeCP2 aging DiffBind files (≈53 MB each, arrived Jun 23) are absent, Panel F skips with a warning and Panels A–E still render. This is an explicitly permitted exception, not a silent fallback; all *required* inputs (annotated peak file, aging Fold/FDR columns) hard-`stop()` if missing.

- **Panel F gene-count not hard-asserted.** The ~1,654 mutant-unique aging-gene figure depends on runtime ChIPseeker annotation, so it is logged (not `stopifnot`-asserted); the verifiable peak-level inputs (10,930 / 23,117) are asserted exactly. The displayed Fisher p uses 1-decimal mantissa formatting (`%.1e`); full precision is in the saved `fig4f_aging_methylation_overlap.tsv`.

- **Effective replication.** 8 samples (2 ctrl + 2 mut per sex); sex/genotype remain partially confounded at n=2/group per sex. The aging contrast comes from collaborator (Jai) young-vs-adult MeCP2 CUT&RUN, a separate DiffBind from the adult methylation arm.

- **No P12 ctrl-vs-mut MeCP2 volcano.** The paper outline's "P12 = no DMRs" panel requires a P12 DiffBind contrast not currently in the repo (only the young-vs-adult aging files are present). Panel A is the *adult* contrast; the developmental story is carried by the aging trajectory (Panels D–F), not a P12 volcano.
