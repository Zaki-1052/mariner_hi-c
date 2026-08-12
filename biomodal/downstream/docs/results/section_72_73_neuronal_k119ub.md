# Sections 72-73: Constitutive H2AK119ub enrichment at neuronal genes (absolute signal) and disproportionate chromatin remodeling at K119ub-high neuronal genes

## Summary
This group asks two linked questions, both deliberately MeCP2-independent. Section 72 uses **absolute** (not differential) H2AK119ub gene-body signal to test whether neuronal/axon-guidance genes are *constitutively* K119ub-enriched in wildtype cerebellum, defining a GO-derived neuronal gene set (5,614 genes) reused by sections 73-78. Section 73 then tests whether those constitutively K119ub-high neuronal genes undergo *disproportionate* chromatin remodeling upon BAP1-KO (DiffBind ATAC / K27ac / K27me3 log2FC). The answers: (1) neuronal genes carry significantly elevated K119ub in both ctrl and mut (Wilcoxon p < 2.2e-16, script-reported), with an overall increasing dose-response (Spearman rho = 0.648, p = 0.049, decile vs neuronal-fraction) and strong enrichment among the K119ub top decile (OR = 1.70); (2) at the same K119ub-high tier, neuronal genes gain MORE accessibility/K27ac than non-neuronal genes, and a significant K119ub×neuronal interaction term drives ATAC and K27ac change — partly consistent with, but more nuanced than, a pure heterochromatin-shift model.

## Section 72: section_72_k119ub_neuronal_characterization
### Analysis question
Are neuronal genes constitutively (in wildtype) enriched for absolute H2AK119ub gene-body signal, independent of MeCP2, and does the enrichment scale with signal level? (Tests the PMC12652333 framework that axon-guidance / neuronal-identity genes are disproportionately H2AK119ub-associated.)

### Key results
- Working universe = 21,604 quantifiable genes; GO-derived neuronal set = 5,614 genes; neuronal genes in universe = 5,077 = 23.50% genome-wide fraction (source: 72_neuronal_gene_k119ub_signal.tsv; 72_neuronal_gene_set_go_derived.tsv; 72_fisher_results.tsv `frac_neuronal_genome`)
- Neuronal enrichment among ctrl K119ub top quartile (Q3): OR = 1.378 [1.284, 1.479], p = 7.46e-19, q = 1.49e-18, 1,512/5,401 high genes neuronal (27.99%) (source: 72_fisher_results.tsv)
- Neuronal enrichment among ctrl K119ub top decile (D9): OR = 1.701 [1.544, 1.874], p = 3.36e-26, q = 1.34e-25, frac_neuronal_high = 32.99% (source: 72_fisher_results.tsv)
- Mutant top quartile OR = 1.186 [1.103, 1.274], p = 3.28e-06; mutant top decile OR = 1.374 [1.242, 1.517], p = 6.15e-10 — i.e. neuronal enrichment persists in mutant but is weaker than ctrl at every threshold (source: 72_fisher_results.tsv)
- Dose-response across K119ub deciles: D1 (lowest) = 20.87% neuronal (OR = 0.845, q = 7.64e-03) rising overall (non-strict; D6 dips to 19.54%) to D10 (highest) = 32.99% neuronal (OR = 1.701, q = 3.36e-25); D10 signal range 2.595–44.998 (source: 72_neuronal_decile_summary.tsv)
- Spearman dose-response trend rho = 0.648 (p = 0.049; reproduced from the decile table as cor.test(decile, frac_neuronal) — see Data-quality flags) (source: 72_neuronal_decile_summary.tsv)
- GO BP of ctrl top-quartile K119ub genes: top terms are developmental/Polycomb-type — pattern specification process (q = 1.57e-43, 273 genes), cell fate commitment (q = 1.48e-42, 192), regionalization (q = 1.68e-39, 247), with neuronal terms present (neuron fate commitment q = 5.81e-28, 73 genes) (source: 72_k119ub_topq_ctrl_go_bp.tsv)
- GSEA ranked by absolute ctrl K119ub signal: 1,177 significant BP terms (table stores significant terms only); top positively-enriched neuronal terms include neuron fate specification (NES = +2.48, q = 6.92e-17) and neuron fate commitment (NES = +2.47, q = 1.02e-23) (source: 72_gsea_ctrl_signal_go_bp.tsv)
- GSEA ranked by absolute mut K119ub signal: 919 significant BP terms (source: 72_gsea_mut_signal_go_bp.tsv)
- Mut top-quartile GO: 1,172 significant terms, 61 neuronal; Gained (mut-high not ctrl-high) GO: 0 significant terms (min p.adjust = 0.659 over 4,169 tested) (source: 72_mut_high_go_bp.tsv; 72_gained_high_go_bp.tsv)

### [INTERPRETATION] Biological meaning
[INTERPRETATION] Neuronal-identity and axon-guidance genes sit in a constitutively Polycomb/H2AK119ub-marked state in wildtype cerebellum: they are over-represented among the most K119ub-decorated genes and the enrichment scales with signal (dose-response), exactly as predicted for PRC1 target domains. [INTERPRETATION] The dominant GO signature of the K119ub-high set (pattern specification, cell-fate commitment, regionalization) is the canonical Polycomb developmental-repression program, with neuronal differentiation terms layered on top — consistent with K119ub marking poised/repressed lineage genes. [INTERPRETATION] That the *gained* (mutant-only-high) set yields ZERO significant GO terms while the *ctrl-high* set is sharply enriched indicates BAP1-KO does not create a new coherent functional class of K119ub targets; instead it amplifies signal at loci that are already constitutively marked (the mutant Fisher ORs are weaker than ctrl because the spread compresses, not because new neuronal targets appear). This makes section 72 a clean *baseline-state* definition rather than a differential result.

### Plot inventory
- `72a_k119ub_signal_neuronal_vs_other/` — violin/box of K119ub gene-body signal (log10), neuronal vs non-neuronal, faceted ctrl/mut
- `72b_k119ub_high_neuronal_fisher/` — forest plot of log2(OR) for neuronal enrichment at the 4 thresholds (ctrl/mut × Q3/D9)
- `72c_k119ub_topq_ctrl_go_dotplot/` — GO BP dotplot of ctrl top-quartile K119ub genes (no MeCP2 filter)
- `72d_ctrl_vs_mut_k119ub_go_comparison/` — neuronal-only GO terms across ctrl-high / mut-high / gained sets
- `72e_gsea_ctrl_signal_dotplot/` — GSEA dotplot, genes ranked by absolute ctrl K119ub signal
- `72e_gsea_mut_signal_dotplot/` — GSEA dotplot, genes ranked by absolute mut K119ub signal
- `72e_running_score_neuron_fate_specification/` — GSEA running-score (enrichment) plot for the top neuronal term (neuron fate specification)
- `72f_neuronal_fraction_by_k119ub_decile/` — dose-response: neuronal fraction + log2(OR) by K119ub decile (2-panel)
- `72_composite/` — assembled multi-panel figure for section 72

## Section 73: section_73_neuronal_chromatin_remodeling
### Analysis question
Do constitutively K119ub-high neuronal genes show stronger BAP1-KO chromatin remodeling (DiffBind ATAC / K27ac / K27me3 log2FC) than other genes, and is there a statistical K119ub×neuronal interaction? (No MeCP2 filter; merges K119ub signal with DiffBind all-marks and the section-72 neuronal set.)
NOTE: on-disk output tables for this section carry the `72g_*` prefix (stale from the pre-rename script version) — see Data-quality flags. The committed script writes `73_*`.

### Key results
- Marks, neuronal vs non-neuronal (all genes): ATAC Δmedian = +0.0199 (neuronal MORE positive; n = 3,586 vs 10,606; p = 1.10e-14); K27ac Δmedian = +0.0063 (p = 0.0164); K27me3 Δmedian = −0.0423 (neuronal MORE negative; n = 1,388 vs 2,781; p = 6.22e-04) (source: 72g_mark_stats_neuronal_vs_other.tsv)
- 4-group key test (K119ub-High tier, neuronal vs non-neuronal): ATAC neuronal median = +0.1667 vs other +0.1088 (n = 1,314 vs 2,570; p = 2.90e-11); K27ac neuronal +0.0469 vs other −0.0155 (n = 668 vs 1,370; p = 2.15e-05); K27me3 neuronal −0.0876 vs other −0.0747 (n = 809 vs 1,201; p = 0.219, NS) (source: 72g_4group_stats.tsv)
- 4-group, K119ub-Low tier: differences shrink markedly (ATAC p = 5.12e-03; K27ac p = 0.813 NS; K27me3 neuronal −0.1414 vs other −0.0484 p = 2.80e-04), showing the neuronal-vs-other gap is concentrated in the K119ub-High tier (source: 72g_4group_stats.tsv)
- Interaction model ATAC_fold ~ k119ub_ctrl_signal * is_neuronal (n = 14,192): interaction beta = +0.00832, p = 0.0172 (significant); k119ub main beta = +0.0218 (p = 2.16e-20); additive R² = 0.0190 → interaction R² = 0.0194 (ANOVA p = 0.0172) (source: 72g_interaction_models.tsv)
- Interaction model K27ac_fold (n = 8,338): interaction beta = +0.0462, p = 4.55e-07 (significant, strongest); k119ub main beta = −0.0242 (p = 3.13e-05); ANOVA p = 4.55e-07. K27me3_fold (n = 4,169): interaction beta = +0.00304, p = 0.760 (NS); ANOVA p = 0.760 (source: 72g_interaction_models.tsv)
- Chromatin-state composition of K119ub-high genes (neuronal vs non-neuronal Fisher): Repressed_Promoter OR = 1.737, p = 6.72e-17 (679 neuronal vs 1,009 other); Bivalent_Promoter OR = 2.866, p = 3.38e-07 (57 vs 42); "Other"/unmarked OR = 0.407, p = 1.08e-34 (290 vs 1,129, neuronal depleted); Active_Promoter OR = 1.139, p = 0.082 (NS) (source: 72g_chromatin_state_fisher.tsv)
- Dose-response by decile (neuronal genes, D1→D10): ATAC median +0.0897 → +0.2172; K27ac +0.0980 → +0.0767; K27me3 −0.4331 → −0.0940 (source: 72g_decile_mark_summary.tsv)

### [INTERPRETATION] Biological meaning
[INTERPRETATION] At the highest constitutive K119ub tier, neuronal genes do respond differently from other genes to BAP1 loss, and the significant positive K119ub×neuronal interaction terms for ATAC (p = 0.017) and especially K27ac (p = 4.6e-07) show neuronal genes are *more sensitive* to their K119ub level when chromatin is remodeled. [INTERPRETATION] However, the direction is the opposite of a naive heterochromatin-shift prediction: neuronal K119ub-high genes show MORE accessibility gain and MORE K27ac gain than non-neuronal genes (positive medians/deltas), not the ATAC-down/K27ac-down/K27me3-up signature the script set out to test; the K27me3 interaction is non-significant. [INTERPRETATION] The chromatin-state Fisher results explain why: K119ub-high neuronal genes are disproportionately Repressed_Promoter (OR 1.74) and Bivalent_Promoter (OR 2.87) and strongly depleted of unmarked "Other" (OR 0.41) — they are poised/bivalent Polycomb loci, and upon BAP1-KO these specific neuronal Polycomb-poised promoters gain accessibility/acetylation. [INTERPRETATION] So the data support "disproportionate remodeling at neuronal K119ub-high loci" (the group thesis) but argue for de-repression/activation of poised neuronal genes rather than a uniform heterochromatinization — the simple "MeCP2 follows heterochromatin" narrative in the script header is not cleanly confirmed by the mark directions here.

### Plot inventory
- `73_multimark_neuronal_vs_other/` — (panel 73a) violin/box of ATAC/K27ac/K27me3 DiffBind log2FC, neuronal vs non-neuronal
- `73_4group_k119ub_x_neuronal/` — (panel 73b) 4-group violin/box: K119ub level × neuronal status, per mark
- `73_chromatin_state_k119ub_high/` — (panel 73c) stacked chromatin-state composition of K119ub-high neuronal vs non-neuronal genes
- `73_dose_response_decile_marks/` — (panel 73d) median mark change across K119ub deciles, neuronal vs non-neuronal lines
- `73_interaction_coefficients/` — (panel 73e) forest plot of K119ub×neuronal interaction coefficients per mark
- `73_composite/` — assembled multi-panel figure for section 73

## Cross-section synthesis
Section 72 establishes the baseline: neuronal-identity genes are constitutively H2AK119ub-marked in wildtype cerebellum, enriched among the most-K119ub genes (top-decile OR 1.70) with a clean dose-response, and dominated by a Polycomb developmental-repression GO signature. Section 73 then shows that this constitutive K119ub status predicts *how* the same genes remodel under BAP1-KO: at the K119ub-high tier, neuronal genes change chromatin marks significantly more than non-neuronal genes (positive K119ub×neuronal interaction for ATAC and K27ac), and these neuronal loci are specifically the bivalent/repressed-promoter Polycomb-poised class. Together they support the paper's thesis that BAP1 loss restructures the epigenome preferentially at constitutively Polycomb/K119ub-marked neuronal genes during cerebellar neurodevelopment — though the observed gain of accessibility/acetylation at poised neuronal promoters refines the mechanism toward de-repression of poised neuronal genes rather than blanket heterochromatinization. Both sections are MeCP2-independent and thus provide an orthogonal, chromatin-first anchor for the MeCP2-redistribution and methylation findings elsewhere in the paper.

## Tables used
- `72_fisher_results.tsv` — neuronal-gene Fisher enrichment among K119ub-high genes at ctrl/mut × Q3/D9 thresholds
- `72_neuronal_decile_summary.tsv` — per-decile neuronal fraction, binomial CI, Fisher OR/q, signal ranges (ctrl)
- `72_neuronal_gene_set_go_derived.tsv` — GO-derived neuronal gene set (5,614 genes; reused by sections 73-78)
- `72_neuronal_gene_k119ub_signal.tsv` — per-gene K119ub ctrl/mut signal, log2FC, is_neuronal, decile (21,604 genes)
- `72_k119ub_topq_ctrl_go_bp.tsv` — GO BP enrichment of ctrl top-quartile K119ub genes
- `72_mut_high_go_bp.tsv` — GO BP of mut top-quartile K119ub genes (1,172 sig, 61 neuronal)
- `72_gained_high_go_bp.tsv` — GO BP of gained (mut-high not ctrl-high) genes (0 significant)
- `72_gsea_ctrl_signal_go_bp.tsv` — GSEA ranked by absolute ctrl K119ub signal (1,177 sig terms)
- `72_gsea_mut_signal_go_bp.tsv` — GSEA ranked by absolute mut K119ub signal (919 sig terms)
- `72g_mark_stats_neuronal_vs_other.tsv` — (section 73 panel a) DiffBind mark medians, neuronal vs non-neuronal
- `72g_4group_stats.tsv` — (section 73 panel b) per-mark medians for K119ub level × neuronal status
- `72g_chromatin_state_fisher.tsv` — (section 73 panel c) chromatin-state Fisher OR, neuronal vs non-neuronal among K119ub-high
- `72g_decile_mark_summary.tsv` — (section 73 panel d) per-decile median/mean mark change, neuronal vs non-neuronal
- `72g_interaction_models.tsv` — (section 73 panel e) lm interaction models mark_fold ~ k119ub_signal × is_neuronal

## Data-quality flags
- **Table-prefix mismatch / stale-prefix (section 73).** The committed `section_73_neuronal_chromatin_remodeling.R` writes tables with the `73_*` prefix (lines 447/481/485/489/494), but the on-disk tables carry the OLD `72g_*` prefix and there are NO `73_*` tables in `tables/`. Commit 8ba6681c ("rename") renamed the script `72g`→`73` and the plot folders, and changed the table-write prefixes, but the script was NOT re-run afterward, so the `72g_*` files persist. The DATA is current and valid (generated 2026-06-21 15:58, same session as section 72's 15:45-15:46; schemas match the script's write blocks exactly) — only the filename prefix is stale. Re-running section 73 would emit identically-contented `73_*` tables.
- **Plot-panel naming inside section 73.** Plot subfolders use `73_<descriptor>` (e.g. `73_4group_k119ub_x_neuronal`) rather than `73b_…`; some files inside were renamed from `72gb_…`/`72gc_…` by commit 8ba6681c. Panel-letter mapping (a→multimark, b→4group, c→chromatin_state, d→dose_response, e→interaction) is inferred from the script.
- **Spearman rho (section 72f) is not stored in a TSV but is reproducible from the decile table.** `72_neuronal_decile_summary.tsv` contains the per-decile fractions; computing `cor.test(decile, frac_neuronal, method="spearman")` directly yields rho = 0.648, p = 0.049 (CORRECTED from the doc's original 0.68, which did not match). The trend is overall increasing but NOT strictly monotone — D6 dips to 19.54% before D10 rises to 32.99%.
- **GSEA tables store significant terms only.** `72_gsea_ctrl_signal_go_bp.tsv` and `72_gsea_mut_signal_go_bp.tsv` total-rows equal sig-rows (1,177 and 919) because gseGO was run with pvalueCutoff = 0.05; the "n significant" counts are exact, not truncated.
- **Gained-high GO is genuinely empty of hits** (4,169 terms tested, min p.adjust = 0.659). This is a real negative result (no coherent functional class among newly-K119ub-high mutant genes), not a missing/broken table.
- **Universe count off-by-one.** Per-gene table has 21,604 data rows; the script summary text references the quantifiable universe and Fisher `frac_neuronal_genome` = 0.23500 reconciles with 5,077/21,604 = 0.2350. No material discrepancy.
- **No MeCP2 anywhere** in either section (by design); no CUT&RUN/ChIP labeling concerns arise here.
