# Paper Figure TODO — Analyses to Complete

Items where upstream data exists in this repo but the specific figure/analysis called for in the paper outline is missing or incomplete.

---

## Figure 1

### TODO 1.1 — P13 TAD Volcano Plot ✅ DONE
**Paper panel:** 1B — "Volcano plot - Differential TAD IDs at P13 and adult"
**What was missing:** Adult volcano exists (`tads/tad-pc-analysis/output/tad_analysis/tad_volcano_*.pdf`), but no equivalent P13/early volcano had been generated.
**Resolution:** Added TADCompare format auto-detection to `scripts/tad_volcano_plot.R`. P-values derived from Gap_Score Z-scores (`2*pnorm(-|Z|)`, BH FDR). Thresholds adjusted for Z-score scale (|Z|>2.0 relaxed, |Z|>3.0 standard). Output: `tads/tad-pc-analysis/output/tad_analysis_early/` — 23,104 boundaries analyzed; 1,419 significant at standard, 2,950 at relaxed.
**Upstream data:**
- `tads/results/early/tadcompare/tadcompare_differential_only.tsv`
- `tads/results/early/consensus/high_confidence_differential.tsv`
- `tads/results/early/final/tadcompare_final_annotated.tsv`
**Reference script:** `scripts/tad_volcano_plot.R` (currently uses HOMER `getDiffExpression.pl` format — may need adaptation for TADCompare output columns)
**Notes:** TADCompare output has different columns than HOMER differential (boundary type, gap score vs logFC/FDR). May need a TADCompare-specific volcano or a unified representation that works across both timepoints.

### TODO 1.2 — % of Genome with Differential PC1
**Paper panel:** 1D — "% of genome with differential PC1"
**What's missing:** The compartment volcano and significant region TSVs exist, but no explicit calculation or figure showing the % of the genome affected.
**Upstream data:**
- `tads/tad-pc-analysis/output/compartment_analysis/compartment_significant_relaxed.tsv`
- `tads/tad-pc-analysis/output/compartment_analysis/compartment_significant_standard.tsv`
- `tads/tad-pc-analysis/output/compartment_analysis/compartment_all_annotated.tsv`
**Task:** Sum region widths from significant compartment shifts / mm10 genome size (~2.73 Gb). Produce a summary statistic and possibly a bar/pie chart showing % A-to-B vs B-to-A vs unchanged. Should be computed at both thresholds (relaxed and standard).

---

## Figure 3

### TODO 3.1 — Cross-Timepoint Permutation Test (P12 Marks Predict Adult Structures)
**Paper panel:** 3D — "Permutation tests - P12 marks predict adult structures; P12 and adult loops against adult loops/boundaries/TADs/compartments"
**What's missing:** No permutation test has been run testing whether P12 (early) ChIP-seq peaks are enriched at adult structural features (differential loops, TAD boundaries, compartment shifts). The timepoint comparison in `tads/results/visualizations/comparison/` compares boundary types across timepoints but does NOT test the predictive question.
**Upstream data:**
- P12 ChIP peaks: `peaks/beds/H3K27acCerebellumEarly2.bed`, `H3K27me3CerebellumEarly1.bed`, `H3K4me1CerebellumEarly1.bed`, `H3K4me3CerebellumEarly2.bed`, `Bivalent_Cerebellum_Early.bed`
- P12 differential peaks: `peaks/new/P12_H3K27me3_up.bed`, `P12_H3K27me3_down.bed`
- Adult structural features:
  - Differential loops: `outputs/250402-late_outputs/merged_loops/characterized_loops.tsv` or `peaks/loop_annotation_extended/late/extended_characterized_loops.tsv`
  - Differential boundaries: `tads/results/late/final/differential_boundaries_final.bed`
  - Compartment shifts: `tads/tad-pc-analysis/output/compartment_analysis/compartment_significant_*.tsv`
- Adult ChIP peaks (for comparison): `peaks/beds/H3K27acCerebellumLate2.bed`, etc.
- Early differential loops: `outputs/250831-early_outputs/merged_loops/characterized_loops.tsv`
**Task:** Permutation-based enrichment test: Do P12 epigenetic marks (especially differential marks) overlap adult structural changes more than expected by chance? Test both P12 and adult marks against adult structures to show that early marks are predictive. Consider using `regioneR` or a custom shuffle-based approach.
**Reference scripts:** `tads/results/late/boundary_loop_analysis/permutation_test_results.tsv` shows the pattern for boundary-loop permutation tests; `scripts/diff_chip_polycomb_enrichment.R` shows ChIP enrichment at loops.

---

## Figure 4

### TODO 4.1 — Concordance Pie Chart
**Paper panel:** 4B — "Pie chart of concordant vs discordant and of the 4 categories"
**What's missing:** The existing figure (`abc/results/enhancer_subset_analysis/16_abc_rnaseq_concordance_bar/`) is a bar plot, not a pie chart. The paper specifically requests a pie chart showing concordant vs discordant proportions and the 4 concordance categories.
**Upstream data:**
- `abc/results/enhancer_subset_analysis/abc_rnaseq_concordance_by_class.tsv`
- `abc/results/discordant_gene_characteristics.tsv`
- `abc/results/figures/discordant_analysis/01_discordant_composite/01_discordant_composite.pdf` (may already contain a pie — verify)
**Task:** Generate a pie chart (or dual pie chart) showing: (1) concordant vs discordant proportion out of 957 DEGs, (2) breakdown of the 4 categories (concordant-up, concordant-down, discordant-up-ABC/down-RNA, discordant-down-ABC/up-RNA). May already be in the discordant composite — check before creating new.

### TODO 4.2 — K27ac and ATAC Volcano Plots Within Combined ATAC Peaks
**Paper panel:** 4C — "Volcano plots for K27ac, ATAC, K119ub within combined ATAC peaks"
**What's missing:** K119ub volcano at enhancers exists (`abc/results/figures/k119ub_correlation/09_k119ub_volcano_at_enhancers/`). K27ac and ATAC volcanoes restricted to the combined ATAC peak set do not exist locally.
**Upstream data (partial):**
- Differential peak BEDs: `peaks/atac_seq/ATAC_up.bed`, `ATAC_down.bed`
- Differential K119ub: `peaks/new/H2AK119ub_up.bed`, `H2AK119ub_down.bed`
- Combined ATAC peaks: `peaks/atac_seq/consensus_all.bed`
- K119ub DiffBind full results: `abc/K119ub_allATAC_diffbind_results_summit_appended_ap.txt`
**What may be needed from HPC:** Full DiffBind statistical tables for H3K27ac and ATAC (with logFC + FDR per peak, not just significant BEDs). If these exist on HPC, download them. If not, re-run DiffBind subsetting to combined ATAC regions.
**Notes:** Check if K27ac DiffBind results exist somewhere in the repo or on HPC. The K119ub full table exists locally so that volcano is done. The remaining two marks may need their DiffBind tables transferred.

---

## Figure 5

### TODO 5.1 — Combined Structural Score & Network Analysis
**Paper panel:** 5C — "Network analysis: Nodes = genes with combined structural changes (boundary, loop, ABC); Node size = AxC score change; Node color = Gene expression logFC"
**What's missing:** No network figure exists. This requires integrating three structural layers (TAD boundaries, loops, ABC) per gene, computing a combined disruption score, and visualizing as a network.
**Upstream data:**
- Boundary genes: `tads/results/visualizations/late/boundary_genes.tsv`
- Loop-gene assignments: `abc/results/loops_with_gene_assignments.tsv`
- Loop data: `peaks/loop_annotation_extended/late/extended_characterized_loops.tsv`
- ABC gene summary: `abc/results/gene_level_summary.tsv` (13,588 genes with strongest ΔABC per gene)
- Paired anchor matches: `abc/results/paired_anchor_matches.tsv`
- RNA-seq: `tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx`
**Task:**
1. Build a per-gene table joining: (a) boundary proximity/disruption, (b) differential loop count and max logFC, (c) ΔABC score
2. Compute combined structural score (e.g., sum of z-scored disruption metrics)
3. Filter to genes with changes in 2+ layers
4. Visualize as network using `igraph` or `ggraph` in R, or `networkx` in Python
   - Nodes = genes
   - Node size = AxC score change magnitude
   - Node color = RNA-seq log2FC (blue = down, red = up)
   - Edges = shared regulatory elements or proximity
5. Consider using GO-based clustering for node grouping

### TODO 5.2 — Top 50 Genes Heatmap by Combined Structural Score
**Paper panel:** 5D — "Heatmap of top 50 genes by combined structural score; Columns = logFC genes, logFC AxC"
**What's missing:** No heatmap of top genes ranked by combined structural disruption exists.
**Upstream data:**
- ABC gene summary: `abc/results/gene_level_summary.tsv`
- Paired anchor matches: `abc/results/paired_anchor_matches.tsv`
- TAD boundary data: `tads/results/late/final/tadcompare_final_annotated.tsv`
- Boundary genes: `tads/results/visualizations/late/boundary_genes.tsv`
- Loop data: `peaks/loop_annotation_extended/late/extended_characterized_loops.tsv`
- RNA-seq: `tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx`
**Task:**
1. Reuse the combined structural score from TODO 5.1
2. Rank genes, take top 50
3. Generate a heatmap (pheatmap or ComplexHeatmap) with columns:
   - Gene expression logFC
   - AxC (Activity × Contact) logFC or ΔABC
   - Optionally: boundary shift magnitude, loop logFC, number of structural layers affected
4. Row annotation: gene name, number of disrupted layers, GO category
5. Column clustering should group related metrics

### TODO 5.3 — 3x2 Model Figure Data Summary
**Paper panel:** 5E — "3x2 model of epigenetic layer, structural layer, and functional (DEG) layer at P13 and adult"
**What's missing:** This is a conceptual schematic (Illustrator/BioRender), but it needs quantitative annotations from the data.
**Upstream data:** All data exists across the repo.
**Task:** Compile a summary table of key numbers for each cell of the 3x2 grid:

| | P13 (Early) | Adult (Late) |
|---|---|---|
| **Epigenetic** | # diff H3K27ac peaks, # diff H3K27me3 peaks, # diff K119ub peaks | Same |
| **Structural** | # diff loops, # diff TAD boundaries, % genome compartment shift | Same |
| **Functional** | # DEGs, # DEGs at structural features, concordance rate | Same |

Pull these numbers from existing summary statistics and TSVs. Output a clean table that can be used to annotate the schematic.

---

## Summary

| Priority | TODO | Figure | Complexity | Dependencies |
|----------|------|--------|------------|-------------|
| High | 5.1 | 5C | High | Integrates 3 data layers + visualization |
| High | 5.2 | 5D | Medium | Builds on 5.1 combined score |
| Medium | 3.1 | 3D | Medium | Permutation testing framework needed |
| Medium | 1.1 | 1B | Low-Medium | Adapt volcano script for TADCompare format |
| Low | 1.2 | 1D | Low | Simple arithmetic from existing TSV |
| Low | 4.1 | 4B | Low | May already exist in discordant composite — verify first |
| Low | 4.2 | 4C | Low-Medium | May need HPC data transfer for K27ac/ATAC DiffBind tables |
| Low | 5.3 | 5E | Low | Data compilation for schematic annotation |
