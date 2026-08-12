# Sections 74-78: Gene-set overlaps (neuronal x coordinated x MeCP2-up), MeCP2 signal-direction reconciliation, synapse/axon specificity, MeCP2 developmental aging trajectory, unbiased-neuronal stoichiometry

## Summary

This group answers a set of follow-up questions about how the BAP1-KO methylation/MeCP2 phenotype maps onto neuronal gene identity. It tests (a) whether the three biologically interesting gene sets — broad GO-neuronal, coordinated mC↑/hmC↓, and MeCP2-gained — overlap more than chance; (b) reconciles the apparent paradox that MeCP2 DiffBind shows ~7,686 UP peaks yet gene-body MeCP2 signal drops genome-wide; (c) asks whether synapse/axon genes are chromatin-special beyond generic neuronal genes; (d) quantifies excess age-related MeCP2 accumulation in mutant cerebellum using young-vs-adult CUT&RUN DiffBind; and (e) re-runs the section-61 stoichiometry analyses with an *unbiased* GO-derived neuronal set. The headline answers: MeCP2 redistribution tracks methylation redistribution (OR=5.16) far more than neuronal identity (OR=1.73); MeCP2 UP peaks are 52% distal-intergenic (it vacates gene bodies and concentrates at distal sites); synapse/axon genes are special specifically for Polycomb de-repression (extra K27me3 loss); mutants accumulate ~2x more age-related MeCP2 UP peaks and gain ~22% more fold at shared aging loci; and — the key self-correction — total methylation at neuronal genes *decreases* with the unbiased broad set (mean −0.0022) whereas the section-61 narrow set spuriously showed an *increase* (+0.012), a DMR-selection artifact.

Gene-set definitions used throughout this group:
- **Neuronal (broad):** 5,614 genes, org.Mm.eg.db GO BP terms matching `synap|neuron|axon|dendrit|nervous` (section 72, unbiased by methylation). 4,118 have valid all-marks data.
- **Neuronal (narrow):** 1,149 genes from GO enrichment of significant DMRs (section 61 — circular/biased).
- **Synapse/axon:** GO BP terms matching `synap|axon` only (section 76).
- **Coordinated:** genes with mC significantly up AND hmC significantly down (`mecp2_coordinated_genes.tsv`).
- **MeCP2-Up:** genes with `mecp2_nearest_fdr < 0.05` and `mecp2_mean_fold > 0`.

---

## Section 74: section_74_geneset_overlap_methylation

### Analysis question
Do the broad-neuronal, coordinated, and MeCP2-Up gene sets overlap more than chance (3-way Venn + pairwise Fisher), and what are the absolute 5mC / 5hmC / total methylation levels at neuronal genes in control vs mutant?

### Key results
- Universe = 23,150 genes (quadrant master); partition "None" = 13,165 genes = 56.87% (source: 74_geneset_overlap_counts.tsv).
- Triple overlap (Neuronal ∩ Coordinated ∩ MeCP2-Up) = 16 genes = 0.07%; Neuronal-only = 3,888 (16.79%); Coordinated-only = 4,592 (19.84%); MeCP2-Up-only = 17 (0.07%) (source: 74_geneset_overlap_counts.tsv).
- Coordinated × MeCP2-Up: overlap = 51, OR = 5.16 [3.19, 8.51], adj p = 2.82e-12 (source: 74_pairwise_fisher.tsv).
- Neuronal × MeCP2-Up: overlap = 27, OR = 1.73 [1.05, 2.82], adj p = 0.034 (source: 74_pairwise_fisher.tsv).
- Neuronal × Coordinated: overlap = 1,442, OR = 1.05 [0.98, 1.13], adj p = 0.141 (NS) (source: 74_pairwise_fisher.tsv).
- Neuronal genes with valid methylation = 4,110; per-gene median total methylation = 0.8235 (ctrl) vs 0.8225 (mut); median 5mC = 0.6782 (ctrl) → 0.6858 (mut, ↑); median 5hmC = 0.1272 (ctrl) → 0.1168 (mut, ↓) (source: 74_neuronal_methylation_levels.tsv, medians computed from per-gene columns).

### [INTERPRETATION] Biological meaning
[INTERPRETATION] The strongest association by far is Coordinated × MeCP2-Up (OR=5.16), an order of magnitude above Neuronal × MeCP2-Up (OR=1.73), and Neuronal × Coordinated is not even significant (OR=1.05). This says MeCP2 redistribution in BAP1-KO tracks the *methylation* redistribution, not generic neuronal identity — methylation state, not cell-type lineage, is what re-positions MeCP2. At neuronal genes the canonical mC↑/hmC↓ signature is present, but because the per-gene hmC drop (median −0.010) slightly exceeds the mC gain (+0.008), *total* methylation edges down — a precursor to the section-78 self-correction.

### Plot inventory
- `74a_geneset_venn/` — 3-way Venn (Neuronal/Coordinated/MeCP2-Up) with pairwise Fisher ORs in caption.
- `74b_neuronal_total_methylation/` — total (mC+hmC) violins ctrl vs mut at neuronal genes.
- `74c_neuronal_5mc_levels/` — 5mC violins ctrl vs mut.
- `74d_neuronal_5hmc_levels/` — 5hmC violins ctrl vs mut.
- `74_composite/` — all panels combined.

---

## Section 75: section_75_mecp2_signal_reconciliation

### Analysis question
Reconcile the paradox that DiffBind reports ~7,686 MeCP2 UP peaks while gene-body MeCP2 BigWig signal drops genome-wide, and contrast K119ub's broad GSEA footprint against MeCP2's narrow (synapse-only) one.

### Key results
- 21,195 genes carry ≥1 MeCP2 peak; only 2,052 have a significant UP peak, 518 have DOWN-only peaks, 18,625 have no significant peaks (source: 75_mecp2_peak_gene_summary.tsv).
- Peak-level totals: 7,686 UP and 1,200 DOWN significant peaks (source: 75_peak_annotation_distribution.tsv, summed; gene-summed UP = 7,683 — three UP peaks lack a gene SYMBOL).
- MeCP2 UP peaks are 51.7% Distal Intergenic and 41.0% Intron, only 2.2% Promoter; DOWN peaks are 61.8% Intron, 19.5% Distal Intergenic, 8.0% Promoter (source: 75_peak_annotation_distribution.tsv).
- K119ub GSEA = 115 significant GO BP terms (3 neuronal); MeCP2 GSEA = 1 significant term (1 neuronal — "synapse assembly") (source: 75_go_term_comparison.tsv; cross-checked against 61k_gsea_k119ub_go_bp.tsv and 61k_gsea_mecp2_go_bp.tsv).
- The MeCP2 sig term contains 0 non-neuronal terms (1/1 = 100% neuronal) vs K119ub 3/115 ≈ 2.6% neuronal (source: 75_go_term_comparison.tsv).

### [INTERPRETATION] Biological meaning
[INTERPRETATION] There is no real paradox: the 7,686 UP peaks concentrate on just ~2,000 genes (≈10% of MeCP2-bound genes) while the other ~90% lose signal slightly, so the genome-wide gene-body median falls even as a subset gains strongly. MeCP2 *redistributes* rather than uniformly increasing. The annotation split is the mechanistic kicker — UP peaks are predominantly distal-intergenic (52%) while DOWN peaks are genic (62% intronic), meaning MeCP2 is being recruited to new distal regulatory sites and vacating gene bodies. The GSEA contrast confirms the division of labor: K119ub blankets ~115 developmental/Polycomb programs broadly, whereas MeCP2's only coherent functional enrichment is a single synaptic-assembly term — MeCP2 is the selective, neuronal-specific reader.

### Plot inventory
- `75a_peak_distribution_per_gene/` — bar chart of genes by MeCP2 peak status (UP / DOWN-only / none).
- `75b_genebody_signal_by_peak_status/` — violins showing UP-peak genes do carry positive gene-body fold.
- `75c_gsea_term_comparison/` — stacked bar, K119ub (115) vs MeCP2 (1) significant terms.
- `75d_peak_annotation_distribution/` — annotation breakdown for UP vs DOWN peaks (chi-squared).
- `75_composite/` — all panels combined.

---

## Section 76: section_76_triple_overlap_synapse_chromatin

### Analysis question
Attach explicit p-values to the neuronal-vs-non-neuronal chromatin violins, characterize the 16 triple-overlap genes, and test whether synapse/axon genes show stronger chromatin remodeling than broader neuronal genes; also test K119ub top-decile enrichment.

### Key results
- 16 triple-overlap genes (Neuronal ∩ Coordinated ∩ MeCP2-Up): Ap3b1, Astn2, Cntn6, Epyc, Fgf1, Fut9, Gprin3, Hcn1, Hif1a, Il1rapl1, Lgi1, Micu3, Ntn4, Prom1, Snca, Tspan7 (source: 76_triple_overlap_genes.tsv, is_triple==TRUE).
- Synapse/axon vs broader neuronal, K27me3: Δ median = −0.0444 (−0.1316 vs −0.0872), Wilcoxon p = 2.95e-3, adj p = 6.65e-3 (significant) (source: 76_synapse_vs_neuronal_stats.tsv).
- Synapse/axon vs broader neuronal, ATAC: p = 0.654 (NS); K27ac: p = 0.757 (NS) — synapse genes are NOT special for accessibility or enhancer activation (source: 76_synapse_vs_neuronal_stats.tsv).
- Broader-neuronal vs non-neuronal K27me3 is NOT significant (p = 0.289), so the K27me3-loss specificity resides in the synapse/axon subset (source: 76_synapse_vs_neuronal_stats.tsv).
- K119ub top-decile (D10) enrichment: Neuronal 677/4,030, OR = 2.34 [2.11, 2.60], p = 1.34e-54; Synapse/axon 306/2,069, OR = 1.68 [1.47, 1.92], p = 2.12e-13; Triple overlap 0/16, OR = 0, p = 0.396 (NS) (source: 76_top_decile_fisher.tsv).
- 4-group chromatin (triple/neuronal-only/coordinated-only/rest) pairwise Wilcoxon: e.g. ATAC neuronal-only vs rest median 0.107 vs 0.116 (p adj 0.102); K27me3 neuronal-only median −0.103 vs coordinated-only +0.237 (adj p 7.0e-46) (source: 76_chromatin_stats.tsv).

### [INTERPRETATION] Biological meaning
[INTERPRETATION] Synapse/axon genes ARE special, but in a specific way: they lose more K27me3 than broader neuronal genes while showing no extra accessibility or enhancer gain. That is selective Polycomb de-repression — H2AK119ub/H3K27me3 erosion at synaptic loci — rather than active enhancer activation, fitting the model where BAP1 loss disorganizes Polycomb domains. The triple-overlap 16 are too few for robust per-mark statistics and are not K119ub-extreme (0/16 in D10), so they are not driven by the strongest K119ub sites; instead it is the *broad* neuronal and synapse sets that are disproportionately K119ub-high (OR 2.34 and 1.68), placing neuronal genes squarely in the K119ub-rich chromatin that BAP1 normally polices.

### Plot inventory
- `76a_neuronal_chromatin_with_pvalues/` — ATAC/K27ac/K27me3 neuronal-vs-non-neuronal violins with annotated p-values and n.
- `76b_four_group_chromatin/` — triple / neuronal-only / coordinated-only / rest comparison across 3 marks.
- `76c_synapse_vs_neuronal_chromatin/` — synapse/axon vs broader-neuronal vs non-neuronal violins.
- `76d_k119ub_decile_enrichment/` — forest plot of D10 odds ratios.
- `76_composite/` — all panels combined.

---

## Section 77: section_77_mecp2_aging_trajectory

### Analysis question
Using young-vs-adult MeCP2 CUT&RUN DiffBind (files supplied by collaborator Jai), separate normal developmental MeCP2 gain from mutant-excess: count aging peaks per genotype, find mutant-unique aging loci/genes, GO-characterize them, and compare fold magnitude at shared aging loci.

### Key results
- Control aging: 10,930 UP / 2,822 DOWN / 383,418 NS peaks; Mutant aging: 23,117 UP / 10,646 DOWN / 361,952 NS peaks — mutant has 2.1× more aging-UP and 3.8× more aging-DOWN peaks (source: 77_aging_peak_summary.tsv).
- Peak overlap: 7,305 ctrl aging-UP peaks overlap mut aging-UP (66.8% of 10,930 ctrl); ctrl-unique = 3,625; mut-unique = 15,812 peaks (source: 77_aging_overlap.tsv).
- Gene-level: ctrl aging-UP = 2,908 genes, mut aging-UP = 4,274 genes, shared = 2,620; ⇒ mut-unique aging genes = 4,274 − 2,620 = 1,654; ctrl-unique = 288 (source: 77_aging_overlap.tsv, mut-unique by arithmetic).
- GO enrichment of mut-unique aging genes = 435 significant BP terms, of which 49 are neuronal (`synap|neuron|axon|dendrit|nervous`) (source: 77_mut_specific_aging_go.tsv, is_neuronal column).
- Fold at 7,305 shared aging loci: median ctrl fold = 1.829, median mut fold = 2.241 (mut +0.41 log2 higher) (source: 77_shared_peak_fold_comparison.tsv, medians computed from columns); 2,266/7,305 shared peaks are at neuronal genes (source: 77_shared_peak_fold_comparison.tsv, is_neuronal column).

### [INTERPRETATION] Biological meaning
[INTERPRETATION] MeCP2 normally accumulates as neurons mature; the data confirm that and then quantify the mutant excess. Two-thirds of the control aging program is preserved in the mutant (7,305/10,930 shared), so the mutant is not rewiring aging wholesale — it is *amplifying* it: 1,654 genes gain age-related MeCP2 only in BAP1-KO, these are GO-enriched for neuronal/synaptic programs (49 of 435 terms), and even at loci both genotypes share the mutant climbs higher (median fold 2.24 vs 1.83). [INTERPRETATION] A plausible driver is that excess H2AK119ub in BAP1-KO creates aberrant chromatin that attracts additional MeCP2 during the normal maturation window, converting a developmental ramp into an overshoot.

### Plot inventory
- `77a_aging_overview/` — dodged bar of UP/DOWN/NS peaks by genotype.
- `77b_aging_overlap_venn/` — gene-level Venn (2,620 shared / 288 ctrl-only / 1,654 mut-only).
- `77c_mut_specific_go_enrichment/` — GO dotplot for the 1,654 mut-unique aging genes.
- `77d_shared_peak_fold_comparison/` — scatter of ctrl vs mut aging fold at 7,305 shared peaks (above diagonal = mut ages more).
- `77_composite/` — all panels combined.

---

## Section 78: section_78_stoichiometry_broad_neuronal

### Analysis question
Re-run the section-61 stoichiometry/mechanism analyses with the *unbiased* GO-derived neuronal set (5,614) instead of the DMR-enriched narrow set (1,149), to test whether section 61's "total methylation increases at neuronal genes" was a selection artifact, and to characterize the MeCP2-Up + K119ub-Up quadrant.

### Key results
- Total methylation Δ (mut − ctrl), mean: All genes −0.00139 (n=20,915), Non-neuronal −0.00119 (n=16,797), Neuronal (broad) −0.00222 (n=4,118), Synapse/axon −0.00128 (n=2,099) — all DOWN, all q < 5e-9; Coordinated +0.00697 (n=6,069, UP), MeCP2-Up +0.03132 (n=79, UP) (source: 78_total_methylation_summary.tsv).
- **Direction self-correction:** narrow neuronal (1,149) mean Δtotal = +0.01200 (INCREASES) vs broad neuronal (4,118) mean Δtotal = −0.00222 (DECREASES); broad-only (in broad not narrow, 2,969 genes) = −0.00772 (source: computed from diffbind_gene_level_all_marks.tsv against 61_neuronal_gene_set.tsv and 72_neuronal_gene_set_go_derived.tsv; narrow/broad in-script values not stored in a TSV).
- Stoichiometry slopes (δ-5mC ~ δ-5hmC): Neuronal (broad) = −0.995 [−1.039, −0.952], consistent with −1; Synapse/axon = −1.020 [−1.078, −0.962], consistent with −1; All genes = −0.959 [−0.978, −0.940] differs; Non-neuronal = −0.949 differs; Coordinated = −1.288 [−1.323, −1.252] differs (steeper) (source: 78_stoichiometry_slopes.tsv).
- BAP1-KO vs TET-KO Spearman rho: All = 0.220 (n=20,851), Non-neuronal = 0.246, Neuronal (broad) = 0.137, Synapse/axon = 0.121, Coordinated = 0.176 — all p < 3e-8; weakest at neuronal/synapse (source: 78_tetko_comparison.tsv).
- Quadrant (72 MeCP2-Up + K119ub-Up genes): broad neuronal = 25/72 (35%; Synapse/axon 17 + Broader-neuronal 8), Synapse/axon = 17/72 (24%), non-neuronal = 47/72 (source: 78_quadrant_neuronal_comparison.tsv, class n's; 61h_mecp2_up_k119ub_up_genes.tsv = 72 genes).

### [INTERPRETATION] Biological meaning
[INTERPRETATION] This section is a deliberate methodological self-correction and it lands cleanly. The section-61 narrow set was built from genes selected *because* they had significant DMRs, so testing "do neuronal genes change methylation?" on that set is circular — and indeed it produces a spurious +0.012 increase. With the unbiased GO set the sign flips to −0.0022: total methylation at neuronal genes actually *decreases*, and total methylation only rises where it should (coordinated genes +0.007, MeCP2-Up +0.031). [INTERPRETATION] The slopes are the mechanistically informative result: neuronal and synapse/axon genes sit at −1.0 (stoichiometric mC-for-hmC exchange, the signature of dehydroxymethylase-like DNMT3A activity), whereas the genome-wide deviation from −1 is driven by non-neuronal genes (−0.949). Coordinated genes overshoot to −1.29 (net mC gain, matching their +0.007 total). And neuronal genes resemble TET-KO *least* (rho 0.137 vs 0.246 non-neuronal), reinforcing that neuronal loci undergo direct 5hmC→5mC conversion rather than TET inhibition plus de-novo methylation.

### Plot inventory
- `78a_total_methylation_forest/` — forest plot of mean Δtotal by gene group with Wilcoxon CIs.
- `78b_stoichiometry_scatter/` — δ-5mC vs δ-5hmC colored by neuronal class, OLS vs −1 reference.
- `78c_stoichiometry_slope_forest/` — slopes with 95% CI vs −1 reference (colored by consistency).
- `78d_bap1_vs_tetko/` — BAP1-KO vs TET-KO δ-ratio scatter.
- `78e_narrow_vs_broad_bias/` — narrow / broad-only / full-broad Δtotal violins (the direction flip).
- `78g_quadrant_characterization/` — 72 quadrant genes by class (MeCP2 / K119ub / 5mC / 5hmC).
- `78_composite/` — all panels combined.

---

## Cross-section synthesis

These five sections converge on one refinement of the paper's thesis: BAP1 loss restructures MeCP2 and methylation, and that restructuring is read out through *methylation state* and *Polycomb chromatin*, not through generic neuronal identity. Section 74 shows MeCP2-gain tracks coordinated mC↑/hmC↓ (OR=5.16) far more than neuronal membership (OR=1.73); section 75 explains that MeCP2 redistributes — concentrating at distal-intergenic sites on ~2,000 genes while thinning gene bodies elsewhere — and that MeCP2's functional footprint is narrowly synaptic while K119ub's is broad; section 76 localizes the neuronal "specialness" to selective K27me3 loss at synapse/axon genes (Polycomb de-repression), and shows neuronal genes are K119ub-rich, the natural BAP1 substrate; section 77 shows the mutant amplifies the normal developmental MeCP2 ramp (2× more aging-UP peaks, +22% fold at shared loci, 1,654 mutant-unique neuronal-enriched genes), tying the phenotype to neurodevelopment; and section 78 corrects an earlier selection-biased claim, establishing that neuronal genes lose total methylation but do so *stoichiometrically* (slope −1.0, dehydroxymethylase-like), distinct from the TET-inhibition pattern at non-neuronal loci. Together they sharpen the mechanism from "neuronal genes are affected" to "methylation-redistributed, K119ub-rich, Polycomb synaptic genes undergo stoichiometric demethylation-block and become MeCP2 sinks during maturation."

## Tables used
- `74_geneset_overlap_counts.tsv` — 8 exclusive Venn partition counts + percentages.
- `74_pairwise_fisher.tsv` — 3 pairwise Fisher tests (overlap, OR, CI, adj p).
- `74_neuronal_methylation_levels.tsv` — per-gene 5mC/5hmC/total ctrl & mut for 4,110 neuronal genes.
- `75_mecp2_peak_gene_summary.tsv` — per-gene UP/DOWN/NS MeCP2 peak counts + gene_peak_class (21,195 genes).
- `75_go_term_comparison.tsv` — K119ub vs MeCP2 neuronal/other significant GSEA term counts.
- `75_peak_annotation_distribution.tsv` — annotation breakdown (n + pct) for UP vs DOWN MeCP2 peaks.
- `76_synapse_axon_gene_set.tsv` — synapse/axon gene list (reusable).
- `76_triple_overlap_genes.tsv` — merged 20,915-gene table with all set flags + K119ub decile.
- `76_chromatin_stats.tsv` — 4-group pairwise Wilcoxon (Kruskal-Wallis context) per mark.
- `76_synapse_vs_neuronal_stats.tsv` — synapse vs broader-neuronal vs non-neuronal pairwise tests.
- `76_top_decile_fisher.tsv` — K119ub D10 enrichment ORs for triple/synapse/neuronal.
- `77_aging_peak_summary.tsv` — UP/DOWN/NS aging peak counts per genotype.
- `77_aging_overlap.tsv` — shared/unique aging peak and gene counts.
- `77_mut_specific_aging_go.tsv` — 435 GO BP terms for the 1,654 mut-unique aging genes (+ is_neuronal).
- `77_shared_peak_fold_comparison.tsv` — paired ctrl/mut fold at 7,305 shared aging loci.
- `78_total_methylation_summary.tsv` — mean/median Δtotal + Wilcoxon per gene group.
- `78_stoichiometry_slopes.tsv` — OLS δ-5mC~δ-5hmC slopes with CIs per group.
- `78_tetko_comparison.tsv` — BAP1-KO vs TET-KO Spearman rho per group.
- `78_quadrant_neuronal_comparison.tsv` — per-class mark means for the 72 quadrant genes.
- Cross-checked source tables (not produced here): `61k_gsea_k119ub_go_bp.tsv`, `61k_gsea_mecp2_go_bp.tsv`, `61_neuronal_gene_set.tsv`, `72_neuronal_gene_set_go_derived.tsv`, `61h_mecp2_up_k119ub_up_genes.tsv`, `diffbind_gene_level_all_marks.tsv`.

## Data-quality flags
- **Section 77 was NOT blocked.** Despite the `data_available` graceful-skip guard for Jai's young-vs-adult DiffBind files, both inputs are present (`peaks/MeCP2_ctrl_adultvsyoung_diffbind_results.txt` and `..._mut_...`, ~53 MB each, dated Jun 23) and all four section-77 output tables are populated. The section ran fully. (Note: the prior root doc and `downstream/CLAUDE.md`/MEMORY flag section 77 as "blocked on external DiffBind data from Jai" — that status is now stale; the files arrived.)
- **Section 76a p-values are NOT in any TSV.** The neuronal-vs-non-neuronal violin p-values quoted in the root prose doc (ATAC p≈6.95e-15, K27ac p≈0.027, K27me3 p≈5.53e-4) are computed in-script and printed to stdout only; `76_chromatin_stats.tsv` holds the *4-group* pairwise tests, not the 2-group 76a tests. [UNVERIFIED: ATAC 6.95e-15 / K27ac 0.027 / K27me3 5.53e-4 per root sections_74-78_results.md, not confirmed in tables.]
- **Section 74 "Δ median" convention.** The root prose lists total Δ median = −0.003, 5mC +0.004, 5hmC −0.009. These are *medians of paired per-gene differences* (mut−ctrl), printed in-script, not stored in a TSV. The difference-of-medians computed directly from the per-gene TSV is total −0.0010, 5mC +0.0076, 5hmC −0.0104. Both describe the same mC↑/hmC↓/total↓ direction; the small numeric gap is the median-of-differences vs difference-of-medians distinction, not an error.
- **Section 78e narrow/broad means are in-script only.** The direction-flip values (narrow +0.012, broad −0.002, broad-only −0.008) are not written to a TSV; they were re-derived here directly from `diffbind_gene_level_all_marks.tsv` against the narrow (`61_neuronal_gene_set.tsv`, 1,149) and broad (`72_...`, 5,614) gene lists and confirmed: +0.01200 / −0.00222 / −0.00772.
- **Two universe sizes coexist (not an error).** Section 74 uses the 23,150-gene quadrant master; sections 76 and 78 use the 20,915-gene `diffbind_gene_level_all_marks.tsv`. Neuronal-with-data counts differ accordingly (4,110 in 74's methylation table vs 4,118 in the diffbind table). This reflects different input tables, not a stale run.
- **Peak vs gene UP totals.** 7,686 UP peaks (annotation table) vs 7,683 summed over gene rows — three UP peaks have no gene SYMBOL and drop out of the per-gene summary. Minor, expected.
- **Small-n caveats.** Triple overlap n=16 and several quadrant per-class cells (n=4–17) are underpowered; MeCP2-Up stoichiometry slope (n=79) has a wide CI [−1.31, −0.54]. Interpret these comparisons cautiously.
- **n=4 effective replicates.** 8 samples (2 ctrl, 2 mut per sex); sex covariate included in run-5 but sex/genotype remain partially confounded at n=2/group per sex.
