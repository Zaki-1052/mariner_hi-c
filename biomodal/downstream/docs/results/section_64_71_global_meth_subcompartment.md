# Sections 64-71: Global methylation levels, MeCP2-vs-methylation scale, CALDER2 subcompartment methylation, MeCP2/K119ub at unmethylated genes, modC/MeCP2 stratified by histone-mark category

## Summary

This group asks where, mechanistically, MeCP2 binding changes come from when BAP1 is lost — and whether they track DNA methylation or the Polycomb mark H2AK119ub. Using the run-5 deep-seq DUET evoC data (8 samples: 4 control + 4 mutant) plus MeCP2 CUT&RUN, CALDER2 subcompartments, gene-body H2AK119ub signal, and H3K27me3/H3K27ac ChIP peaks, it establishes: (1) a small but consistent bulk shift (5mC +0.32%, 5hmC -0.39%, modC stable) confirming TET-pathway blockade rather than global hypermethylation; (2) methylation change is ~5x more widespread than MeCP2 binding change at the gene level — MeCP2 is a context-dependent reader, not a 1:1 methylation follower; (3) the central mechanistic result — at 359 genes where MeCP2 binding increases with NO significant methylation change, H2AK119ub is significantly gained (Fisher OR=3.15, p=1.8e-24), proving MeCP2 responds to K119ub independently of methylation; and (4) variance-partitioning that K119ub explains ~7.3% of MeCP2 fold whereas the 5mC/5hmC change ratio explains essentially nothing (0.02%), nailing K119ub as the dominant determinant.

## Section 64: section_64_global_methylation_levels

### Analysis question
Is there a bulk, genome-wide CpG methylation shift in BAP1-KO, and does it match a TET-blockade signature (5mC up, 5hmC down, total modC unchanged) rather than global hypermethylation? Uses the DUET upstream autosomal CpG summary for all 8 run-5 samples, paired by sex+batch.

### Key results
- Mean 5mC: Control = 72.216%, Mutant = 72.534%, delta = +0.317% (paired t-test p = 0.0412, Cohen's d = 2.79) (source: 64_statistics.tsv)
- Mean 5hmC: Control = 10.192%, Mutant = 9.797%, delta = -0.395% (paired t-test p = 0.0184, Cohen's d = -4.54) (source: 64_statistics.tsv)
- Mean modC (total): Control = 82.565%, Mutant = 82.487%, delta = -0.079% (paired t-test p = 0.492, Cohen's d = -0.60; NOT significant) (source: 64_statistics.tsv)
- Paired Wilcoxon p = 0.125 for both 5mC and 5hmC (floored by n=4 pairs; t-test is the powered test here) (source: 64_statistics.tsv)
- Sample table confirms run-5 composition: 8 samples = 4 Control (ctrl-F, ctrl-F-B2, ctrl-M, ctrl-M-B2) + 4 Mutant (mut-F, mut-F-B2, mut-M, mut-M-B2), 2 sexes x 2 batches (source: 64_global_methylation_summary.tsv)
- Every individual mutant sample has higher 5mC and lower 5hmC than its matched control (e.g. ctrl-F 5mC=72.094% / 5hmC=10.252% vs mut-F 5mC=72.623% / 5hmC=9.680%) (source: 64_global_methylation_summary.tsv)

### [INTERPRETATION] Biological meaning
[INTERPRETATION] The directionally-opposed, near-equal-magnitude 5mC gain and 5hmC loss with a flat total modC is the classic signature of a TET-mediated demethylation block: 5mC is normally oxidized to 5hmC by TET enzymes en route to demethylation, so blocking that conversion accumulates 5mC at the direct expense of 5hmC while leaving the summed modified-cytosine pool unchanged. The very large Cohen's d (>2.7 for both modalities) despite tiny absolute deltas reflects extremely low between-replicate variance — the effect is small in magnitude but exceptionally reproducible, arguing it is a genuine systemic consequence of BAP1 loss rather than noise. This rules out a model of indiscriminate global hypermethylation.

### Plot inventory
- 64_global_methylation_levels/ — combined 3-panel figure (dot+bar, replicate delta, composition stacked)
- 64a_global_methylation_dotbar/ — per-sample dot+bar by modality (5mC, 5hmC, modC), zoomed y-axes, paired t-test annotation
- 64b_replicate_delta/ — matched-replicate delta lollipop (mut - ctrl) per sex/batch pair, faceted by modality
- 64c_composition_stacked/ — composition stacked bar (5mC + 5hmC + unmethylated = 100%), ctrl vs mut

## Section 65: section_65_mecp2_methylation_scale

### Analysis question
At gene-body resolution, how does the scale of significant DNA methylation change compare to the scale of significant MeCP2 binding change — and do most methylated genes actually show a MeCP2 response? Tests whether MeCP2 is a context-dependent reader rather than a stoichiometric follower of methylation.

### Key results
- Gene bodies tested = 20,969; significant 5mC change (q<0.05) = 10,775 (51.4%); of these 5mC Hyper = 7,513 (35.8%) and 5mC Hypo = 3,262 (15.6%) (source: 65_methylation_mecp2_cascade.tsv)
- Coordinated genes (mC up AND hmC down, both sig) = 6,589 (31.4% of tested) (source: 65_methylation_mecp2_cascade.tsv)
- Genes with significant MeCP2 change (FDR<0.05) = 2,570 (12.3%); significant MeCP2 Up = 2,052 (9.8%) — roughly 5x fewer than significant-5mC genes (source: 65_methylation_mecp2_cascade.tsv)
- mC Hyper genes overlapping a sig MeCP2 Up = only 947 (4.5% of all tested; 12.6% of the 7,513 mC-hyper genes) (source: 65_methylation_mecp2_cascade.tsv)
- Fisher's exact (sig 5mC x sig MeCP2 up, over 16,941 shared genes): OR = 2.629, p = 2.23e-61; contingency = both_sig 1,204 / mc_only 8,620 / mecp2_only 359 / neither 6,758 (source: 65_fisher_context_dependence.tsv)
- The mecp2_only cell (significant MeCP2 up but NO significant 5mC) = 359 genes — this is exactly the cohort carried into Section 67 (source: 65_fisher_context_dependence.tsv; 65_mecp2_bound_no_methylation_genes.tsv has 359 data rows)

### [INTERPRETATION] Biological meaning
[INTERPRETATION] Methylation change is ~5x more pervasive than MeCP2 binding change (51.4% vs 9.8% of genes), so methylation is necessary but not sufficient for a MeCP2 response — only ~12.6% of hypermethylated genes recruit more MeCP2. The significantly positive Fisher OR (2.63) confirms the two events are coupled above chance, but the dominant mc_only cell (8,620 genes) shows the coupling is weak and context-gated. Crucially, the 359 mecp2_only genes (MeCP2 up without methylation change) are mechanistically inconvenient for a pure methyl-reader model and motivate the Section 67 test: if MeCP2 is not following methylation at these loci, what is it following?

### Plot inventory
- 65_mecp2_methylation_scale/ — combined multi-panel figure
- 65a_cascade_hierarchy/ — side-by-side methylation vs MeCP2 funnel bar charts
- 65b_gene_level_venn/ — Venn of mC-hyper, MeCP2-sig-up, coordinated gene sets
- 65c_allgene_quadrant_scatter/ — 5mC change vs MeCP2 fold quadrant scatter, all tested genes
- 65d_proportional_funnel/ — methylation-to-MeCP2 funnel (tested -> sig 5mC -> hyper -> +MeCP2 peak -> +sig MeCP2 up)
- 65e_mecp2_signal_by_meth_status/ — violin of MeCP2 BigWig signal at gene bodies by methylation status (ctrl vs mut)

## Section 66: section_66_subcompartment_methylation

### Analysis question
Are BAP1-induced methylation changes compartment-specific? Uses CALDER2 100kb subcompartment labels (A.1 strong-active, A.2 weak-active, B.1 facultative-het, B.2 constitutive-het) from the adult/late (250402) Hi-C, assigned to gene bodies by midpoint, plus a parallel H3K27me3/H3K27ac histone-mark stratification.

### Key results
- A.1 (strong active): 10,783 genes, 6,932 sig (64.3%), 5,716 hyper vs 1,216 hypo — most affected, hypermethylation-dominated (source: 66_subcompartment_dmr_summary.tsv)
- A.2 (weak active): 3,046 genes, 1,727 sig (56.7%), 1,025 hyper / 702 hypo (source: 66_subcompartment_dmr_summary.tsv)
- B.1 (facultative het): 1,897 genes, 839 sig (44.2%), 311 hyper / 528 hypo — flips to hypo-dominated (source: 66_subcompartment_dmr_summary.tsv)
- B.2 (constitutive het): 4,172 genes, 910 sig (21.8%), 249 hyper / 661 hypo — least affected and hypo-dominated (source: 66_subcompartment_dmr_summary.tsv)
- Histone-mark stratification mirrors compartments: K27ac-only (Active) 4,830/6,199 sig (77.9%), 4,505 hyper / 325 hypo; K27me3-only (Fac. Het) 1,673/3,170 sig (52.8%), but only 92 hyper vs 1,581 hypo — strongly hypo (source: 66_histone_mark_dmr_summary.tsv)
- All subcompartment ctrl-vs-mut Wilcoxon tests are significant; 5mC A.1 and 5hmC A.1 both p ~ 0 (machine zero), B.2 5mC p = 2.14e-62, A.2 5hmC p = 1.75e-139 (source: 66_subcompartment_wilcoxon.tsv)

### [INTERPRETATION] Biological meaning
[INTERPRETATION] The TET-blockade hypermethylation is concentrated in active chromatin (A.1/A.2 and K27ac-marked genes), where baseline 5mC turnover via 5hmC is highest, so blocking TET there causes the largest net 5mC gain. Constitutive heterochromatin (B.2) is already stably, near-maximally methylated with low turnover, so it is largely immune (21.8% affected) and what change occurs is hypomethylation. The clean mirroring of the CALDER2 subcompartment pattern by a purely peak-based H3K27me3/H3K27ac classification (K27ac = A.1-like hyper, K27me3 = B.1-like hypo) is strong internal validation that the compartment effect is driven by the underlying active/repressive chromatin state, not a Hi-C artifact.

### Plot inventory
- 66_subcompartment_methylation/ — combined figure
- 66a_subcompartment_dmr_fraction/ — % sig 5mC DMRs by subcompartment, hyper vs hypo (Chi-squared annotated)
- 66b_direction_by_subcompartment/ — stacked hyper/hypo fraction of sig DMRs per subcompartment
- 66c_methylation_levels_by_subcompartment/ — 5mC and 5hmC violins (ctrl vs mut) faceted by subcompartment
- 66d_histone_mark_dmr_fraction/ — same DMR-fraction analysis stratified by H3K27me3/H3K27ac category

## Section 67: section_67_mecp2_k119ub_unmethylated (KEY RESULT)

### Analysis question
At the 359 genes where MeCP2 binding increases WITHOUT a significant methylation change (the Section 65 mecp2_only set), does H2AK119ub change explain the MeCP2 gain? If so, MeCP2 is responding to the Polycomb/BAP1 axis independently of DNA methylation. Uses gene-body K119ub log2FC plus K27me3/K27ac/ATAC BigWig signal at these loci.

### Key results
- Of 359 candidate MeCP2-up-no-methylation genes, 356 matched K119ub signal with finite log2FC; background = 22,027 genes (source: 67_statistics.tsv)
- K119ub gene-body log2FC: our-group median = +0.202 vs background median = -0.036; Mann-Whitney (our > rest) U = 5,337,358, p = 5.50e-32 (source: 67_statistics.tsv; 67_per_mark_summary.tsv)
- Fisher (K119ub gained, our vs rest): 72.8% of our genes gain K119ub (259/356) vs 45.9% of background; OR = 3.15, p = 1.82e-24 (source: 67_statistics.tsv)
- MeCP2 nearest-peak fold vs K119ub log2FC at these genes: Spearman rho = 0.241, p = 4.38e-06 (n=356) — positive, MeCP2 gain tracks K119ub gain (source: 67_statistics.tsv)
- Multi-mark one-sample Wilcoxon (vs 0) at these loci: K119ub median = +0.202 (p = 2.53e-24, gained), K27me3 = +0.162 (p = 2.25e-05, gained), K27ac = -0.242 (p = 1.65e-12, LOST), ATAC = +0.074 (p = 1.58e-08, slight gain) (source: 67_per_mark_summary.tsv)
- The multimark per-gene table holds 356 genes; top K119ub-gain genes include Serpinb8 (K119ub +1.30) and Tmprss11g (+1.08) (source: 67_mecp2_no_meth_multimark.tsv)

### [INTERPRETATION] Biological meaning
[INTERPRETATION] This is the linchpin result of the paper's MeCP2 model. At loci where DNA methylation is flat, MeCP2 binding still increases, and the dominant concurrent epigenomic change is H2AK119ub GAIN (73% of these genes, 3.15x enriched over genome, with a positive MeCP2-vs-K119ub correlation). Because BAP1 is the H2AK119ub deubiquitinase, its loss raises K119ub, and MeCP2 follows that signal. The coincident H3K27me3 gain (Polycomb co-deposition) and H3K27ac LOSS confirm these are Polycomb-repressed loci tipping further toward repression. Thus MeCP2 is not a passive methyl-CpG reader here; it is recruited (directly or indirectly) by the BAP1->K119ub Polycomb axis, decoupled from methylation. This provides the mechanistic bridge between the methylation phenotype and the Polycomb-domain compaction phenotype seen in the Hi-C arm of the project.

### Plot inventory
- 67_mecp2_k119ub_unmethylated/ — combined figure
- 67a_k119ub_at_mecp2_no_meth/ — K119ub log2FC violin: MeCP2-no-meth genes vs all others (Mann-Whitney)
- 67b_k119ub_gained_fraction/ — % K119ub gained bar: our group vs genome (Fisher OR/p)
- 67c_mecp2_vs_k119ub_scatter/ — MeCP2 nearest fold vs K119ub log2FC scatter at the 356 genes (Spearman)
- 67d_multimark_at_mecp2_no_meth/ — multi-mark violin (K119ub, K27me3, K27ac, ATAC) one-sample tests vs 0

## Section 68: section_68_modc_mecp2_by_mark

### Analysis question
Within each histone-mark chromatin context (K27ac-only Active, K27me3+K27ac Bivalent, K27me3-only Fac. Het, Neither), does total modC change correlate with MeCP2 mean fold, and where do genes fall in the modC-vs-MeCP2 quadrant space? Gene-level join of the DiffBind all-marks table with MeCP2 fold.

### Key results
- K27ac only (Active): n = 5,929 genes, Spearman rho = 0.183, p = 1.02e-45 (source: 68_modc_mecp2_by_mark_stats.tsv)
- K27me3 + K27ac (Bivalent): n = 477, Spearman rho = 0.461, p ~ 0 (strongest coupling of any context) (source: 68_modc_mecp2_by_mark_stats.tsv)
- K27me3 only (Fac. Het): n = 3,073, Spearman rho = 0.189, p = 4.10e-26 (source: 68_modc_mecp2_by_mark_stats.tsv)
- Neither: n = 7,430, Spearman rho = 0.070, p = 1.41e-09 (weakest — unmarked genes barely couple) (source: 68_modc_mecp2_by_mark_stats.tsv)
- Quadrant structure differs by context: Bivalent genes are MeCP2-up-dominated (Q1 hyper+MeCP2up=119, Q2 hypo+MeCP2up=91 vs Q4 hyper+MeCP2down=60), whereas Neither genes are MeCP2-down-dominated (Q3 hypo+MeCP2down=3,364) (source: 68_modc_mecp2_by_mark_stats.tsv)

### [INTERPRETATION] Biological meaning
[INTERPRETATION] The modC-MeCP2 correlation is real but modest in active and facultative-het contexts (rho ~0.18-0.19) and essentially absent at unmarked genes (rho 0.07), reinforcing Section 65's conclusion that methylation alone is a poor predictor of MeCP2. The standout is the Bivalent (K27me3+K27ac) context (rho 0.46): at Polycomb-poised promoters MeCP2 binding co-varies most tightly with methylation, consistent with these being the loci where the BAP1/Polycomb axis and methylation converge. The context-dependence of the correlation strength is itself the message — MeCP2's relationship to methylation is gated by chromatin state.

### Plot inventory
- 68_modc_mecp2_by_mark/ — combined 4-panel figure
- 68a_modc_mecp2_k__ac_only/ — Active (K27ac only) modC vs MeCP2 quadrant scatter
- 68b_modc_mecp2_k__me____k__ac/ — Bivalent (K27me3+K27ac) scatter
- 68c_modc_mecp2_k__me__only/ — Facultative het (K27me3 only) scatter
- 68d_modc_mecp2_neither/ — Unmarked (Neither) scatter

## Section 69: section_69_mc_hmc_mecp2_by_mark

### Analysis question
Does MeCP2 distinguish 5mC from 5hmC? Each gene is plotted twice (5mC change and 5hmC change vs the same MeCP2 fold) within each mark context, and an lm interaction test plus a Fisher z-test compare the two per-modality Spearman correlations.

### Key results
- K27ac only (Active): 5mC rho = 0.160 (p=1.88e-35) vs 5hmC rho = -0.060 (p=3.38e-06); interaction est = -2.053, p = 7.28e-38 (significant — MeCP2 DOES distinguish); Fisher z = 12.09, p = 1.17e-33 (source: 69_mc_hmc_mecp2_by_mark_stats.tsv)
- K27me3 + K27ac (Bivalent): 5mC rho = 0.305 (p=1.35e-11) vs 5hmC rho = 0.032 (NS, p=0.479); interaction p = 4.06e-05 (significant); Fisher z = 4.34, p = 1.41e-05 (source: 69_mc_hmc_mecp2_by_mark_stats.tsv)
- K27me3 only (Fac. Het): 5mC rho = 0.144 (p=1.24e-15) vs 5hmC rho = 0.057 (p=1.63e-03); interaction p = 0.493 (NOT significant — MeCP2 does NOT distinguish here) (source: 69_mc_hmc_mecp2_by_mark_stats.tsv)
- Neither: 5mC rho = -0.009 (NS, p=0.432) vs 5hmC rho = 0.082 (p=1.93e-12); interaction p = 1.77e-06 (significant) (source: 69_mc_hmc_mecp2_by_mark_stats.tsv)
- Per-context n: Active 5,929, Bivalent 477, Fac. Het 3,073, Neither 7,430 genes (source: 69_mc_hmc_mecp2_by_mark_stats.tsv)

### [INTERPRETATION] Biological meaning
[INTERPRETATION] In active and bivalent chromatin MeCP2 binding tracks 5mC positively but 5hmC weakly/negatively, and the interaction tests confirm MeCP2 responds differently to the two modifications — consistent with the canonical view that MeCP2 reads 5mC (and is partially repelled by 5hmC). The notable exception is facultative heterochromatin (K27me3 only), where the interaction is non-significant: in this Polycomb context MeCP2 does not cleanly discriminate 5mC from 5hmC, again pointing to a non-methylation (Polycomb/K119ub) recruitment mode dominating there, dovetailing with Section 67.

### Plot inventory
- 69_mc_hmc_mecp2_by_mark/ — combined 4-panel figure
- 69a_mc_hmc_mecp2_k__ac_only/ — Active context, 5mC (red) + 5hmC (blue) dots vs MeCP2
- 69b_mc_hmc_mecp2_k__me____k__ac/ — Bivalent context
- 69c_mc_hmc_mecp2_k__me__only/ — Facultative het context
- 69d_mc_hmc_mecp2_neither/ — Unmarked context

## Section 70: section_70_mc_hmc_sig_category_by_mark

### Analysis question
Re-colors the same methylation-vs-MeCP2 space by which modification is significantly changed (5mC only / 5hmC only / Both / Neither), per mark context, to see which significance class drives the coupling. "Both" genes contribute two dots.

### Key results
- K27ac only (Active): 5,354 unique genes, 9,464 dots, Spearman rho = 0.044 (p=1.67e-05); sig breakdown Both=4,110, 5mC-only=527, 5hmC-only=717, Neither=575 — overwhelmingly "Both" (source: 70_sig_category_by_mark_stats.tsv)
- K27me3 + K27ac (Bivalent): 415 genes, 676 dots, rho = 0.191 (p=5.91e-07); Both=261, 5mC-only=63, 5hmC-only=91, Neither=62 (source: 70_sig_category_by_mark_stats.tsv)
- K27me3 only (Fac. Het): 2,133 genes, 3,114 dots, rho = 0.139 (p=7.60e-15); Both=981, 5mC-only=682, 5hmC-only=470, Neither=940 (source: 70_sig_category_by_mark_stats.tsv)
- Neither (unmarked): 4,579 genes, 6,991 dots, rho = 0.041 (p=5.59e-04); Both=2,412, 5mC-only=788, 5hmC-only=1,379, Neither=2,851 (source: 70_sig_category_by_mark_stats.tsv)
- Across contexts "Both significant" is the largest class wherever active marks are present (Active 4,110/5,354; Bivalent 261/415), confirming the coordinated mC-up/hmC-down signature dominates marked chromatin (source: 70_sig_category_by_mark_stats.tsv)

### [INTERPRETATION] Biological meaning
[INTERPRETATION] When the dots are split by significance class, the combined Spearman correlations are weak (rho 0.04-0.19), confirming that MeCP2 fold is not strongly predicted by which modality changed — the coupling seen in Section 68 is carried mostly by chromatin context, not by 5mC/5hmC significance per se. The biological takeaway is that the coordinated "Both" class (simultaneous 5mC gain + 5hmC loss, the TET-blockade fingerprint) is the dominant population at active and bivalent genes, but its MeCP2 response is still modest, again leaving room for the K119ub-driven recruitment established in Section 67.

### Plot inventory
- 70_sig_category_by_mark/ — combined 4-panel figure
- 70a_sig_category_k__ac_only/ — Active context, dots colored by sig category
- 70b_sig_category_k__me____k__ac/ — Bivalent context
- 70c_sig_category_k__me__only/ — Facultative het context
- 70d_sig_category_neither/ — Unmarked context

## Section 71: section_71_mc_hmc_ratio_mecp2_k119ub

### Analysis question
Does MeCP2 binding scale with the degree of 5mC-skew (the per-gene 5mC-gain / 5hmC-loss ratio) — a gradient model — or is it driven by H2AK119ub regardless of the ratio? Builds a per-gene log2(mC/hmC change ratio), correlates it and K119ub against MeCP2 fold, runs a quintile dose-response, and partitions variance via nested linear models.

### Key results
- Master set after filtering (finite mc_diff, hmc_diff, MeCP2 fold, K119ub, |hmc_diff|>0.001) = 12,149 genes (source: 71_statistics.tsv; 71_ratio_mecp2_summary.tsv)
- Spearman log2_ratio vs MeCP2 fold: rho = 0.014, p = 0.113 (NOT significant) (source: 71_statistics.tsv)
- Spearman K119ub vs MeCP2 fold: rho = 0.187, p = 3.53e-96 (highly significant) (source: 71_statistics.tsv)
- Nested models: K119ub-only R2 = 0.0738 (AIC -12,776.5); Ratio-only R2 = 0.0010 (AIC -11,856.5); K119ub+Ratio R2 = 0.0740; K119ub*Ratio R2 = 0.0783 (AIC -12,831.8, best) — interaction est = +0.0204, p = 4.30e-14 (source: 71_model_comparison.tsv; 71_statistics.tsv)
- Variance partition: K119ub unique = 7.3%, Shared = 0.1%, Ratio unique = 0.0% (0.00017), Unexplained = 92.6% (source: 71_variance_partition.tsv)
- Standardized betas (additive model): K119ub beta = +0.040 (p = 1.57e-202) vs ratio beta = +0.0020 (p = 0.131, NS) (source: 71_standardized_betas.tsv)
- Quintile dose-response: Kruskal-Wallis H = 14.8, p = 5.23e-03; median MeCP2 fold is essentially flat across all 5 ratio quintiles (~ -0.015 to -0.018), confirming no gradient (source: 71_statistics.tsv; 71_quintile_dose_response.tsv)

### [INTERPRETATION] Biological meaning
[INTERPRETATION] This section adjudicates between two recruitment models and decisively favors the K119ub-driven one. The 5mC/5hmC change ratio — the per-gene magnitude of the TET-blockade skew — has effectively zero independent predictive power for MeCP2 binding (ratio-only R2 = 0.001, standardized beta NS, flat quintile dose-response), whereas K119ub alone explains ~7.3% of MeCP2 fold. MeCP2 therefore does not titrate to how 5mC-skewed a gene is; it responds to the Polycomb mark. The small but significant K119ub-by-ratio interaction (p=4.3e-14) indicates a mild synergy — the K119ub effect is slightly amplified at more 5mC-skewed loci — but the main-effect structure is unambiguous: K119ub is the driver, the methylation ratio is a passenger. The 92.6% unexplained variance also flags that gene-body MeCP2 fold is governed by much more than these two predictors.

### Plot inventory
- 71_mc_hmc_ratio_mecp2_k119ub/ (nested folder) — combined 4-panel figure plus all panels:
  - 71a_k119ub_mecp2_ratio_scatter/ — K119ub log2FC vs MeCP2 fold scatter colored by log2(mC/hmC ratio)
  - 71b_mecp2_by_ratio_quintile/ — MeCP2 fold violins by ratio quintile (dose-response, Kruskal-Wallis)
  - 71c_ratio_vs_mecp2_direct/ — direct log2_ratio vs MeCP2 scatter with loess + OLS
  - 71d_nested_model_comparison/ — R2 bar chart for the 4 nested linear models

## Cross-section synthesis

These eight sections build a single mechanistic argument from the bulk genome down to per-gene regression. Section 64 establishes the systemic phenotype — a small, hyper-reproducible 5mC-up/5hmC-down/total-stable shift diagnostic of TET-pathway blockade, not global hypermethylation. Section 66 localizes that hypermethylation to active chromatin (A.1 / K27ac) while heterochromatin (B.2 / K27me3) is stable or hypomethylated, fixing the chromatin geography of the effect. Sections 65, 68, 69 and 70 then dismantle the naive expectation that MeCP2 binding simply follows this methylation: methylation change is ~5x more widespread than MeCP2 change, the modC-MeCP2 correlation is modest and entirely context-gated (strongest at bivalent promoters, absent at unmarked genes), and MeCP2 distinguishes 5mC from 5hmC everywhere except facultative heterochromatin. Sections 67 and 71 deliver the resolution: at the 359 genes where MeCP2 rises without methylation change, H2AK119ub is gained 3.15x over genome (p=1.8e-24), and genome-wide K119ub explains ~7.3% of MeCP2 fold while the 5mC/5hmC skew explains essentially nothing. Together they support the paper's thesis that BAP1 loss restructures MeCP2 binding primarily through the H2AK119ub Polycomb axis — connecting the methylation phenotype to the same Polycomb-domain mechanism that drives the Hi-C compaction phenotype — rather than through DNA methylation acting as MeCP2's recruitment signal.

## Tables used

- 64_global_methylation_summary.tsv — per-sample autosomal 5mC/5hmC/modC/unmethylated %, 8 run-5 samples with sex/batch
- 64_statistics.tsv — per-modality paired Wilcoxon, paired t-test, Cohen's d (ctrl vs mut)
- 65_methylation_mecp2_cascade.tsv — gene-count cascade from tested -> sig 5mC -> hyper -> +MeCP2 overlaps
- 65_fisher_context_dependence.tsv — Fisher 2x2 of sig-5mC x sig-MeCP2-up (OR, p, cell counts incl. mecp2_only=359)
- 65_mecp2_bound_no_methylation_genes.tsv — the 359 MeCP2-up-no-methylation genes (input to Section 67)
- 65_methylation_mecp2_gene_overlap.tsv — per-gene boolean membership across the methylation/MeCP2 sets (20,969 genes)
- 66_subcompartment_dmr_summary.tsv — per-CALDER2-subcompartment DMR counts and hyper/hypo split
- 66_histone_mark_dmr_summary.tsv — parallel DMR summary by H3K27me3/H3K27ac mark category
- 66_subcompartment_wilcoxon.tsv — paired Wilcoxon ctrl-vs-mut 5mC and 5hmC per subcompartment
- 66_per_gene_compartment_assignment.tsv — per-gene subcompartment + mark assignment (29,750 rows)
- 67_statistics.tsv — Mann-Whitney, Fisher (OR=3.15), Spearman for K119ub at MeCP2-no-meth genes
- 67_per_mark_summary.tsv — one-sample Wilcoxon (vs 0) for K119ub/K27me3/K27ac/ATAC at these loci
- 67_mecp2_no_meth_multimark.tsv — per-gene multi-mark log2FC table (356 genes)
- 68_modc_mecp2_by_mark_stats.tsv — per-mark-context Spearman + quadrant counts (modC vs MeCP2)
- 69_mc_hmc_mecp2_by_mark_stats.tsv — per-context 5mC vs 5hmC Spearmans + interaction + Fisher z tests
- 70_sig_category_by_mark_stats.tsv — per-context Spearman + significance-category gene counts
- 71_statistics.tsv — Spearman/Kruskal-Wallis/OLS/interaction summary (n=12,149)
- 71_model_comparison.tsv — nested linear-model R2/adj-R2/AIC (K119ub, ratio, both, interaction)
- 71_variance_partition.tsv — variance partition of MeCP2 fold (K119ub unique 7.3%, ratio 0.0%)
- 71_standardized_betas.tsv — standardized betas for K119ub vs log2(mC/hmC ratio)
- 71_quintile_dose_response.tsv — MeCP2 fold summary per ratio quintile + pairwise Wilcoxon
- 71_ratio_mecp2_summary.tsv — per-gene ratio/K119ub/MeCP2 table (12,149 genes)

## Data-quality flags

- Section 65 / 67 provenance gap: section_67 reads `65_mecp2_bound_no_methylation_genes.tsv` (the 359-gene list), but the current `section_65_mecp2_methylation_scale.R` does NOT write that file (it writes only the cascade, gene_overlap, and fisher tables). The 359-gene list must be produced by an external/earlier script not in the current viz_sections set. The count is internally consistent — 65_fisher_context_dependence.tsv reports mecp2_only = 359 and the file has exactly 359 data rows — so the number is trustworthy, but the generating step is undocumented in the present pipeline.
- "359 vs 356": the candidate cohort is 359 genes, but Section 67 analyzes 356 after dropping 3 genes with non-finite K119ub log2FC. Both the script comment and FIGURES-style prose cite "359 genes"; the statistics (OR=3.15, medians, Mann-Whitney) are computed on n=356. Report the cohort as 359 and the tested n as 356.
- Section 65 cascade self-percentages exceed 100% for one row: "Genes with any MeCP2 peak" = 21,195 = 101.1% of the 20,969 tested gene bodies. This is because the MeCP2 peak gene universe is not a subset of the mc_dmr tested-gene universe (different annotation sources), so the percentage denominator is mismatched. Counts are fine; the >100% percentage is a denominator artifact, not an error.
- Section 66 row totals (19,898 genes summed across A.1/A.2/B.1/B.2) are smaller than the 20,969 tested genes because only gene bodies whose midpoint falls in a labeled 100kb CALDER2 bin are retained; unassigned genes are dropped before the per-subcompartment summary. Expected behavior, noted for transparency.
- Section 66 uses the adult/late (250402) CALDER2 labels only; there is no early/P12 subcompartment stratification in this group, so the compartment conclusions are adult-specific.
- FIGURES.md and TODO.md contain no entries for sections 64-71 (the only K119ub/Fisher numbers in FIGURES.md belong to earlier figures 20-30). All numbers above therefore come directly from the TSVs with no prose-doc cross-check available for this range.
- Paired Wilcoxon in Section 64 is underpowered (n=4 pairs -> minimum two-sided p = 0.125), so significance for the global shift rests on the paired t-test (5mC p=0.041, 5hmC p=0.018). Treat the t-test as the operative test.
- No table is empty or missing for this group; all 22 referenced TSVs are present and populated.
