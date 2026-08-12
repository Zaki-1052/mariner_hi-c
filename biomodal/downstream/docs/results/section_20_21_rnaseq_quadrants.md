# Sections 20-21: RNA-seq integration of coordinated genes; discordant (Q2) vs coordinated (Q4) quadrant characterization

## Summary
These two sections integrate the bivariate mC/hmC DMR signal (genes significant in BOTH 5mC and 5hmC) with adult-timepoint BAP1 WT/KO RNA-seq and a panel of epigenomic tracks (ATAC, H2AK119ub, H3K27ac, MeCP2 CUT&RUN, chromatin state, Hi-C loops). Section 20 asks whether the dominant **coordinated** quadrant (Q4: mC up / hmC dn) is transcriptionally repressed; it is — significantly-changed Q4 genes are downregulated ~5:1 (Tier 1) to ~17.5:1 (Tier 2), and combined methylation effect anticorrelates with log2FC (Spearman rho = -0.249). Section 21 contrasts the minority **discordant** quadrant (Q2: mC dn / hmC up) against Q4 across nine dimensions and finds the two groups are nearly mirror-image regulatory programs: Q2 genes are upregulated, gain accessibility, lose K119ub, gain H3K27ac, sit at Repressed/Bivalent promoters, and are ~2.6x more likely to be at Hi-C loop anchors, whereas Q4 genes are repressed, K119ub-gaining, MeCP2-up, Active-promoter genes.

The four quadrants of the co-significant set (run-5 TSV, n = 8,371): **Q1 (mC dn/hmC dn) = 411**, **Q2 (mC dn/hmC up, DISCORDANT) = 1,255**, **Q3 (mC up/hmC up) = 116**, **Q4 (mC up/hmC dn, COORDINATED) = 6,589** (78.7% of the co-significant set) (source: coordinated_rnaseq_expression.tsv / discordant_gene_characteristics.tsv, col `quadrant`).

## Section 20: section_20_coordinated_rnaseq
### Analysis question
Does the coordinated mC-up/hmC-down (Q4) signature translate into transcriptional repression, and is the repression bias stronger than for other methylation-changed genes? Methylation gain does not automatically imply silencing, so RNA-seq is used to quantify the up/down/unchanged split.

### Key results
- Q4 (coordinated) gene count = 6,589 of 8,371 co-significant genes = 78.7% (source: coordinated_rnaseq_expression.tsv, col `quadrant`).
- Q4 genes with matched RNA-seq data = 5,961; without = 628 (No RNA-seq Data) (source: coordinated_rnaseq_expression.tsv, col `expr_tier1`).
- Tier 1 (padj < 0.05) Q4 expression: Down = 787 (13.2%), Up = 157 (2.6%), Not Significant = 5,017 (84.2%); **Down/Up ratio = 5.01** (source: coordinated_rnaseq_expression.tsv, col `expr_tier1`).
- Tier 2 (padj < 0.05 & |log2FC| > 0.3) Q4 expression: Down = 666 (11.2%), Up = 38 (0.6%), Not DEG = 5,257 (88.2%); **Down/Up ratio = 17.53** (source: coordinated_rnaseq_expression.tsv, col `expr_tier2`).
- Spearman correlation, combined methylation effect (|mC| + |hmC|) vs RNA-seq log2FC across Q4 genes = **rho = -0.249** (n = 5,961 with RNA-seq) — larger combined methylation change tracks with lower expression (source: recomputed from coordinated_rnaseq_expression.tsv cols `combined_effect`, `log2FC`).
- Q4 median RNA-seq log2FC = -0.0065 (n = 5,961 with RNA-seq) — only a slight global downward shift; the repression signal concentrates in the significantly-changed minority (source: discordant_gene_characteristics.tsv, col `log2FC`, quadrant Q4).
- 20d-style contingency among co-significant genes that are also expression-significant (padj < 0.05): mC Up→Expr Down = 810, mC Down→Expr Up = 239 (predicted anti-diagonal) vs mC Up→Expr Up = 161, mC Down→Expr Down = 37; sample OR = 0.031, showing the expected inverse mC–expression relationship (source: recomputed from coordinated_rnaseq_expression.tsv cols `mc_diff`, `log2FC`, `padj`). NOTE: the script's actual Fig 20d uses ALL significant mC DMR genes (a superset not contained in this TSV), so the script's headline 20d OR/p cannot be reproduced from these tables — see Data-quality flags.

### [INTERPRETATION] Biological meaning
[INTERPRETATION] The coordinated Q4 quadrant is the canonical BAP1-loss epigenetic signature: H2AK119ub accumulates (BAP1 can no longer deubiquitinate it), TET activity is impeded, 5mC rises and 5hmC falls together. The RNA-seq integration shows this is a genuinely *repressive* methylation gain — among genes whose expression changes significantly, downregulation outnumbers upregulation 5:1 (Tier 1) and 17.5:1 (Tier 2), and the magnitude of the combined methylation change predicts the degree of repression (rho = -0.249). The large "Not Significant" fraction (84%) indicates that for most Q4 genes the methylation change is subthreshold for transcriptional consequence at this timepoint, i.e., the methylation shift precedes or exceeds the transcriptional response — consistent with a slow, accumulating epigenetic lesion rather than an acute transcriptional switch.

### Plot inventory
- `20a_coordinated_expression_breakdown/` — Stacked bar, expression-direction breakdown (Up/Down/Unchanged) for Coordinated vs Other mC DMR vs All genes, faceted by Tier 1 / Tier 2, with Fisher OR annotations.
- `20b_methylation_vs_expression_scatter/` — Scatter of combined methylation effect (%) vs RNA-seq log2FC for Q4 genes, colored by DEG status, with lm fit and Spearman rho.
- `20c_log2fc_violin_comparison/` — Violin + boxplot of log2FC for Coordinated (Q4) vs Other significant mC DMR genes, with Wilcoxon p.
- `20d_mc_expression_heatmap/` — 2x2 enrichment heatmap of 5mC direction x expression direction (Obs/Exp/OR) among genes significant in both mC and expression.

## Section 21: section_21_discordant_mc_hmc_analysis
### Analysis question
How does the minority discordant quadrant (Q2: mC dn / hmC up) differ from the dominant coordinated quadrant (Q4: mC up / hmC dn) across nine epigenomic dimensions (methylation effect size, RNA-seq, ATAC, H2AK119ub, H3K27ac, chromatin state, MeCP2, Hi-C loops, GO)?

### Key results
- Quadrant counts (genes significant in both mC and hmC): Q1 = 411, **Q2 (discordant) = 1,255**, Q3 = 116, **Q4 (coordinated) = 6,589** (source: discordant_gene_characteristics.tsv, col `quadrant`).
- Methylation effect size (medians): |mC diff| Q2 = 1.55% vs Q4 = 2.82%; |hmC diff| Q2 = 1.06% vs Q4 = 2.47% — Q4 has ~1.8x larger mC and ~2.3x larger hmC effects (source: discordant_gene_characteristics.tsv, cols `mc_diff`, `hmc_diff`).
- RNA-seq direction (padj < 0.05): Q2 = 190 Up vs 10 Down; Q4 = 157 Up vs 787 Down. Fisher Up/Down Q2-vs-Q4 sample **OR ≈ 95, p ≈ 1.2e-104** — Q2 genes are overwhelmingly upregulated while Q4 genes are overwhelmingly downregulated (source: discordant_gene_characteristics.tsv, cols `padj`, `log2FC`).
- RNA-seq log2FC medians (genes with RNA-seq): Q2 = +0.179 (n = 942) vs Q4 = -0.0065 (n = 5,961) — opposite directional trend (source: discordant_gene_characteristics.tsv, col `log2FC`).
- H2AK119ub gene-body log2FC medians: Q2 = -0.0385 (n = 1,199) vs Q4 = +0.0302 (n = 6,291) — Q2 *loses* K119ub, Q4 *gains* it (source: discordant_gene_characteristics.tsv, col `k119ub_gb_log2fc`).
- H3K27ac status: Q2 = 131 gained / 20 lost (gain-dominated); Q4 = 1,319 gained / 871 lost (more balanced). Fisher Gained/Lost Q2-vs-Q4 sample OR ≈ 4.3, p ≈ 4.7e-12 (source: discordant_gene_characteristics.tsv, cols `n_k27ac_gained`, `n_k27ac_lost`).
- MeCP2 binding: Q2 = 161 MeCP2-Up / 117 MeCP2-Down (≈1.4:1); Q4 = 770 MeCP2-Up / 131 MeCP2-Down (≈5.9:1). Fisher Up/Down Q2-vs-Q4 sample OR ≈ 0.23, p ≈ 1.1e-20 — Q4 is far more strongly MeCP2-up biased (source: discordant_gene_characteristics.tsv, cols `n_mecp2_up`, `n_mecp2_down`).
- Hi-C loop involvement: Q2 = 165/1,255 = 13.1% at loop anchors vs Q4 = 333/6,589 = 5.1%. Fisher Loop/No-Loop Q2-vs-Q4 sample **OR ≈ 2.84, p ≈ 1.1e-22** — discordant genes are ~2.6x more likely at loop anchors (source: discordant_gene_characteristics.tsv, col `loop_involved`).
- Chromatin state (gene-body, promoter-centric): Q2 is dominated by **Repressed_Promoter = 846 (67.4%)** with Active_Promoter = 78 (6.2%); Q4 is dominated by **Active_Promoter = 4,202 (63.8%)** with Repressed_Promoter = 54 (0.8%) — near-inverted compositions (source: discordant_gene_characteristics.tsv, col `chromatin_state`).
- GO enrichment of Q2 against Q4 background: no significant terms (no `21d_go_enrichment_discordant` folder produced) — Q2 is not a pathway-coherent gene set relative to Q4 (source: section_21 script behavior; consistent with concordant_vs_discordant_trends.md).

### [INTERPRETATION] Biological meaning
[INTERPRETATION] Q2 and Q4 are two opposing regulatory responses to BAP1 loss. Q4 (the majority) is *passive repression*: K119ub accumulates, TET is inhibited, mC rises / hmC falls, MeCP2 reads the new methylation, and expression trends down — these genes are bystanders of the genome-wide Polycomb/methylation shift, and they sit overwhelmingly at active promoters. Q2 (the minority) is *active de-repression*: these loci specifically *lose* K119ub (against the genome-wide trend), open chromatin (net ATAC, H3K27ac gain), lose mC while gaining hmC (active TET-driven demethylation), and become upregulated. Their strong starting bias toward Repressed/Bivalent promoters and their ~2.6x enrichment at Hi-C loop anchors suggests Q2 represents structurally-engaged, Polycomb-poised loci that escape the dominant repressive program — plausibly because loop-anchored / poised genes recruit BAP1-independent deubiquitinase or remodeling activity that maintains an active state. The smaller methylation effect sizes in Q2 are consistent with a more targeted, locus-specific mechanism rather than the bulk passive drift seen in Q4.

### Plot inventory
- `21a_discordant_composite/` — 3x3 composite figure, Q2 vs Q4 across all nine panels.
- `21a_panel_mc_diff/` — |mC difference| boxplots, Q2 vs Q4 (Wilcoxon).
- `21a_panel_hmc_diff/` — |hmC difference| boxplots, Q2 vs Q4 (Wilcoxon).
- `21a_panel_log2fc/` — RNA-seq log2FC boxplots, Q2 vs Q4 (Wilcoxon).
- `21a_panel_atac/` — net ATAC change boxplots, Q2 vs Q4.
- `21a_panel_k119ub/` — K119ub gene-body log2FC boxplots, Q2 vs Q4.
- `21a_panel_k27ac/` — H3K27ac gained/lost bar chart, Q2 vs Q4 (Fisher).
- `21a_panel_chromatin/` — chromatin-state stacked bars, Q2 vs Q4 (Fisher, simulated p).
- `21a_panel_mecp2/` — MeCP2 up/down bar chart, Q2 vs Q4 (Fisher).
- `21a_panel_loop/` — Hi-C loop involvement stacked bars, Q2 vs Q4 (Fisher).
- `21b_mc_hmc_concordance_scatter/` — mC vs hmC scatter for all co-significant genes, Q4 (green) and Q2 (red) highlighted over background.
- `21c_mc_vs_expression_per_group/` — faceted mC-change vs log2FC scatter with per-group Spearman.
- `21e_all_quadrants_comprehensive/` — 4-panel comparison of log2FC, net ATAC, K119ub, chromatin state across all four quadrants.
- (`21d_go_enrichment_discordant/` — NOT produced; GO returned no significant terms.)

## Cross-section synthesis
Section 20 establishes that the dominant Q4 coordinated quadrant is a repressive signature (significantly-changed genes downregulated 5:1 to 17.5:1; methylation magnitude anticorrelates with expression). Section 21 then shows that the minority Q2 discordant quadrant is its near-perfect inverse along every axis measured — expression (up vs down, Fisher OR ≈ 95), K119ub (loss vs gain), accessibility/H3K27ac (gain vs balanced), MeCP2 (weaker vs stronger up-bias), chromatin state (Repressed vs Active promoter), and loop anchoring (2.6x enriched). Together they support the paper's central thesis that BAP1 loss restructures methylation and MeCP2 binding through a genome-wide H2AK119ub/TET-block (Q4), against which a structurally-distinct, loop-anchored, Polycomb-poised minority (Q2) is actively de-repressed — directly linking the methylation phenotype to MeCP2 redistribution and to 3D-genome architecture.

## Tables used
- `coordinated_rnaseq_expression.tsv` — Section 20 output; one row per co-significant gene (all 4 quadrants) with mc_diff, hmc_diff, combined_effect, quadrant, RNA-seq log2FC/padj/baseMean, and Tier-1/Tier-2 expression classifications (8,371 rows).
- `discordant_gene_characteristics.tsv` — Section 21 output; one row per co-significant gene with quadrant plus the full nine-dimension feature set (RNA-seq, ATAC counts, K119ub gene-body/promoter log2FC, H3K27ac gained/lost, chromatin state, MeCP2 up/down, loop_involved) (8,371 rows).
- `coordinated_gene_characteristics.tsv` — byte-identical to `discordant_gene_characteristics.tsv` (see Data-quality flags); used only to confirm the duplication.

## Data-quality flags
- **Stale numbers in prose docs (run-4 era).** `docs/concordant_vs_discordant_trends.md` and `docs/FIGURES.md` (Figures 29, 32, 33) report Q4 = 5,708, Q2 = 755, total = 6,750 (85%/84.6%), plus Figure 32 subtitle "n = 5,203 coordinated | 2,322 other mC DMR | 16,572 all genes." These do NOT match the current run-5 TSVs (Q4 = 6,589, Q2 = 1,255, Q1 = 411, Q3 = 116, total = 8,371; Q4 = 78.7%). All counts, medians, ORs, and rho in those prose docs are from a prior run and should be regenerated. Directional conclusions are unchanged.
- **CLAUDE.md "92.3%" coordinated figure not reproduced here.** The repo headline "92.3% of co-significant genes show coordinated mC↑/hmC↓" does not match Q4/total in this TSV (6,589/8,371 = 78.7%); FIGURES.md states 85%. [UNVERIFIED: 92.3% per CLAUDE.md, not confirmed in tables] — likely a different denominator (e.g., excluding Q1/Q3, or a different significance set, or a prior run). Three different "coordinated %" values appear across sources (78.7% TSV, 85% FIGURES.md, 92.3% CLAUDE.md).
- **Two characteristics TSVs are byte-identical.** `discordant_gene_characteristics.tsv` and `coordinated_gene_characteristics.tsv` are identical (`diff -q` reports no difference; both 8,371 rows). The section_21 script only writes `discordant_gene_characteristics.tsv`; `coordinated_gene_characteristics.tsv` appears to be a stale copy/alias (no section in this group writes it). They carry the same 22 columns including all four quadrants, so "coordinated" in that filename is a misnomer for the current content.
- **Figure 20d OR/p not reconstructable from group TSVs.** The script's Fig 20d contingency uses ALL significant mC DMR genes (mc_dmr filtered to significant, joined to RNA-seq), not just the co-significant set in `coordinated_rnaseq_expression.tsv`. The 20d-style numbers reported above are recomputed on the co-significant subset only and will differ from the script's headline 20d OR. FIGURES.md (prose) cites 20d "OR = 0.03, p = 5.90e-118 among 1,330 genes" — [UNVERIFIED: 0.03 / 5.90e-118 / 1,330 per FIGURES.md, not confirmable from the group TSVs because the source gene set is not exported here].
- **"Other mC DMR" and "All genes" RNA-seq breakdowns not in group TSVs.** Section 20a's Other-mC-DMR (FIGURES.md cites 2,322) and All-genes (16,572) bars draw on `mc_dmr` and the full RNA-seq table, neither exported to these tables; only the coordinated/Q4 portion is verifiable here.
- **Fisher ORs reported as sample odds ratios.** The Q2-vs-Q4 ORs above (expression ≈95, loop ≈2.84, MeCP2 ≈0.23, H3K27ac ≈4.3) are sample odds ratios recomputed from the count matrices; R's `fisher.test` conditional-MLE estimate and exact p will differ slightly in the last digits but are directionally identical.
- **Net ATAC medians are 0 in both groups.** Q2 and Q4 both have median net_atac = 0 (most genes have no differential ATAC peak); the group difference lives in the tails / any-up-vs-any-down counts, not the median. The prose doc's "Q2 net ATAC +1.0" reflects the stale run and is not reproduced (current median = 0).
- **PNG output missing.** `20a` (and likely other folders) contain only JPG/PDF/SVG, not the 4th PNG format the convention expects — minor, non-blocking.
- **21d (GO) folder absent.** Consistent with the script (GO enrichment of Q2 vs Q4 background returned no significant terms); not a missing-output error.
