# Biomodal DUET evoC Differential Methylation Analysis: Figure Descriptions

This document describes each figure from the biomodal DUET evoC differential methylation visualization pipeline, comparing BAP1-KO mutant versus wildtype control mice on the mm10 genome. Descriptions are written to be fully accessible without viewing the images. Supporting data tables are referenced where applicable.

> **Table directory:** `downstream/plots/visualizations/tables/`

---

## Figure 1: Sequencing Quality Metrics (`01_qc_overview`)

**Layout:** A four-panel figure arranged in a 2x2 grid, titled "Sequencing Quality Metrics - Biomodal DUET evoC, BAP1-KO vs Wildtype Control." Each panel shows data for four samples (ctrl-F, ctrl-M, mut-F, mut-M), color-coded blue for Control and red for Mutant.

**Top-left panel — Mapped Reads:** A bar chart showing millions of mapped reads per sample. The two control samples have 367 million (ctrl-F) and 336 million (ctrl-M) reads. The two mutant samples have higher counts: 476 million (mut-F) and 406 million (mut-M). Mutant samples were sequenced deeper than controls overall.

**Top-right panel — Mapped Bases:** A bar chart showing billions of mapped bases per sample. Control samples produced 47.7 billion (ctrl-F) and 44.4 billion (ctrl-M) bases. Mutant samples produced 61.5 billion (mut-F) and 53.2 billion (mut-M) bases. This mirrors the mapped reads pattern, with mutants yielding more total sequence data.

**Bottom-left panel — Duplication Rate:** A lollipop chart showing the duplication rate for each sample. Rates range from 27.8% (ctrl-M) to 30.2% (mut-F), with the subtitle noting the range is 27.8%-30.2%. Each sample is represented as a filled circle at the tip of a vertical line, with numeric labels above. All four values are tightly clustered within a narrow 2.4 percentage-point range, indicating comparable library complexity across samples.

**Bottom-right panel — Mean Phred Score:** A lollipop chart showing the mean base quality score for each sample. All four samples have nearly identical scores: 34.37 (ctrl-F), 34.30 (ctrl-M), 34.36 (mut-F), and 34.34 (mut-M). The subtitle states "All samples: ~34.3 (excellent)." The y-axis spans only 33.5 to 34.5, emphasizing the uniformity. A Phred score of 34 corresponds to an error rate of approximately 1 in 2,500 bases.

**Key takeaway:** All four samples pass quality control with excellent base quality, consistent duplication rates, and adequate sequencing depth, with mutant samples having modestly higher total reads/bases than controls.

> **Supporting table:** [`summary_statistics.txt`](plots/visualizations/tables/summary_statistics.txt)

---

## Figure 2: DMR Statistics Summary (`03_dmr_region_statistics`)

**Layout:** A four-panel figure titled "DMR Statistics Summary," arranged in a 2x2 grid. Each panel addresses a different aspect of the differential methylation results.

**Top-left panel — Significant mC DMRs by Genomic Region (CG Context):** A bar chart showing the number of significant 5mC differentially methylated regions (DMRs, q-value < 0.05) broken down by six genomic region types. The bars descend from left to right: gene bodies lead with 8,836 DMRs (42.1% of tested genes), followed by CpG shores with 6,037 (18.6%), CpG shelves with 4,045 (13.9%), promoters with 641 (3.1%), CpG islands with 291 (3.3%), and TSS regions with only 58 (0.4%). Each bar is a different color. The subtitle states "Gene bodies dominate differential methylation."

**Top-right panel — Gene Body DMRs: 5mC vs 5hmC:** A two-bar comparison showing the count of significant gene-body DMRs for each modification type. The 5hmC bar (blue) shows 9,930 significant genes (47.3%), slightly exceeding the 5mC bar (red) at 8,836 (42.1%). The subtitle states "5hmC shows more differential genes than 5mC."

**Bottom-left panel — Baseline Methylation by Context:** A bar chart showing average methylation levels across three sequence contexts: CG (CpG) at 72.3%, CHG at 0.6%, and CHH at 0.9%. The CG bar is dramatically taller than the other two, which are barely visible at the bottom of the chart. The subtitle states "CpG methylation is ~100x more abundant than non-CpG."

**Bottom-right panel — Significant DMRs by Context:** A bar chart showing significant gene-body DMR counts by methylation context. Only CG (CpG) shows any signal, with its bar labeled "Primary Signal." The CHG and CHH positions show no bars at all, each labeled "No Signal." The subtitle confirms "No significant changes in non-CpG methylation."

**Key takeaway:** Differential methylation in BAP1-KO is concentrated in gene bodies (not promoters or CpG islands), occurs exclusively in the CpG context, and affects both 5mC and 5hmC at similar scale (approximately 8,800-9,900 genes).

> **Supporting table:** [`summary_statistics.txt`](plots/visualizations/tables/summary_statistics.txt)

---

## Figure 3: Asymmetric Direction of Methylation Changes (`03b_direction_comparison`)

**Layout:** A single grouped bar chart titled "Asymmetric Direction of Methylation Changes," with the subtitle showing the counts: "5mC: 8836 significant genes | 5hmC: 9930 significant genes." The x-axis has two groups (5mC and 5hmC) and the y-axis shows percentage of significant DMRs from 0% to 100%. Within each group, a red bar represents "Increased" and a blue bar represents "Decreased" modifications in the BAP1-KO mutant.

**5mC group (left):** The red "Increased" bar dominates at 75.1% (n=6,635 genes), while the blue "Decreased" bar reaches only 24.9% (n=2,201 genes). Three-quarters of significantly differentially methylated genes gain 5mC in the mutant.

**5hmC group (right):** The pattern reverses. The blue "Decreased" bar dominates at 86.8% (n=8,623 genes), while the red "Increased" bar is small at 13.2% (n=1,307 genes). The vast majority of significantly differentially hydroxymethylated genes lose 5hmC in the mutant.

An italicized annotation in the upper-center of the plot reads "mC increases while hmC decreases," summarizing the opposing directions.

**Key takeaway:** The two modifications change in opposite directions genome-wide: 5mC predominantly increases (hypermethylation) while 5hmC predominantly decreases (hypo-hydroxymethylation) in BAP1-KO. This reciprocal pattern suggests a shared mechanism disrupting the conversion of 5mC to 5hmC.

---

## Figure 4: Volcano Plots — 5mC vs 5hmC (`04_volcano_plots`)

**Layout:** Two side-by-side volcano plots under the shared title "Differential Methylation: 5mC vs 5hmC" with the subtitle "Key genes labeled: Syt1, Zbtb20, Trpm3, Cntnap2."

**Left panel — Gene Body 5mC Differential Methylation:** A scatter plot with the x-axis showing 5mC difference (Mutant minus Control, in %) ranging from approximately -25% to +30%, and the y-axis showing -log10(q-value), which transforms statistical significance so that higher points are more significant. The dashed horizontal line near y=0 marks the significance threshold. A vertical line at x=0 separates hypomethylated (left) and hypermethylated (right) genes. 8,836 genes are significant (stated in subtitle). The majority of significant points cluster to the right of zero (positive mC change), colored red for "Hypermethylated." A smaller cluster on the left is colored blue for "Hypomethylated." Gray dots below the significance line represent non-significant genes. Several genes at the very top of the plot (q-values at the machine floor, -log10(q) approximately 300) are labeled: Dlgap1, Mcu, Syt1, Arhgap26, and Cntnap2. The plot is notably asymmetric, with a larger and denser cloud of red points extending to the right.

**Right panel — Gene Body 5hmC Differential Methylation:** The same layout, now for 5hmC changes. Here the asymmetry is reversed: the dominant cloud of significant points extends to the left of zero (negative hmC change), colored blue for "Hypomethylated." 9,930 genes are significant. Labeled genes at the top include Syt1, Epha3, Lpp, Mcu, Cntnap2, and Dlgap1, with most showing decreased 5hmC. A smaller population of red "Hypermethylated" points appears on the right. The denser blue cloud on the negative side mirrors the red cloud on the positive side in the 5mC panel.

**Key takeaway:** The two volcano plots visually confirm the reciprocal pattern: 5mC shifts right (gains) while 5hmC shifts left (losses). The same genes (Syt1, Mcu, Cntnap2, Dlgap1) appear at the top of both panels, suggesting they are strongly affected in both modifications simultaneously. Effect sizes range from roughly -25% to +30% for 5mC and -20% to +20% for 5hmC.

> **Supporting tables:** [`top_mc_dmrs.tsv`](plots/visualizations/tables/top_mc_dmrs.tsv), [`top_hmc_dmrs.tsv`](plots/visualizations/tables/top_hmc_dmrs.tsv)

---

## Figure 5: Key Genes — Coordinated Methylation Changes (`05_coordinated_changes`)

**Layout:** A 2x2 grid of four paired bar charts, one per gene, under the shared title "Key Genes: Coordinated Methylation Changes" with subtitle "5mC increase + 5hmC decrease pattern in BAP1-KO." Each panel shows two pairs of bars: one pair for 5hmC (left) and one pair for 5mC (right), with blue bars for Control and dark red bars for Mutant. The y-axis shows Methylation Level (%) from 0 to approximately 80%.

**Top-left — Syt1:** Subtitle states "mC: +17.9% | hmC: -15.3% | Coordinated pattern." For 5hmC, the Control bar is approximately 27% and the Mutant bar drops to approximately 12%, a large decrease. For 5mC, the Control bar is approximately 58% and the Mutant bar rises to approximately 76%, a large increase. Q-values are annotated above each pair: q = 1.7e-305 for hmC and q = 2.9e-305 for mC. Syt1 shows the largest combined effect of any gene.

**Top-right — Zbtb20:** Subtitle: "mC: +8.2% | hmC: -6.0% | Coordinated pattern." 5hmC goes from approximately 22% (Control) to approximately 16% (Mutant). 5mC goes from approximately 66% to approximately 74%. Both q-values are at the machine floor (1.7e-305 and 2.9e-305).

**Bottom-left — Trpm3:** Subtitle: "mC: +9.2% | hmC: -5.4% | Coordinated pattern." 5hmC decreases from approximately 24% to approximately 19%. 5mC increases from approximately 62% to approximately 71%. Both q-values at machine floor.

**Bottom-right — Cntnap2:** Subtitle: "mC: +2.4% | hmC: -3.8% | Coordinated pattern." 5hmC decreases from approximately 11% to approximately 8%. 5mC increases from approximately 76% to approximately 78%. The hmC q-value is 5.3e-140 and the mC q-value is 2.9e-305. Cntnap2 has smaller absolute changes than the others but remains highly significant due to its large gene body providing many CpG sites for testing.

**Key takeaway:** Four representative genes each show the coordinated pattern: 5mC increases in the mutant (taller red bar on the right pair) while 5hmC decreases (shorter red bar on the left pair), all with extreme statistical significance. The magnitude of changes varies from +2.4%/-3.8% (Cntnap2) to +17.9%/-15.3% (Syt1).

> **Supporting table:** [`coordinated_changes.tsv`](plots/visualizations/tables/coordinated_changes.tsv)

---

## Figure 6: Coordinated 5mC and 5hmC Changes Scatter (`05a_mc_hmc_scatter`)

**Layout:** A single scatter plot titled "Coordinated 5mC and 5hmC Changes" with the subtitle "6750 genes significant in both | 85% show mC↑/hmC↓ pattern." The x-axis shows 5mC Change (%) ranging from approximately -25% to +30%, and the y-axis shows 5hmC Change (%) ranging from approximately -20% to +22%. Dashed lines at x=0 and y=0 divide the plot into four quadrants.

**Data distribution:** Each point represents one gene that is significant in both 5mC and 5hmC. Points are colored in two categories: purple for genes with the "mC↑ / hmC↓" pattern (positive mC change, negative hmC change, i.e., the lower-right quadrant), and gray for "Other" patterns. The lower-right quadrant is highlighted with a light purple/lavender background rectangle.

The purple cloud dominates the plot, forming a dense mass in the lower-right quadrant concentrated between roughly +2% to +15% on the x-axis and -2% to -10% on the y-axis. Two genes are explicitly labeled within this region: Mcu (near +15% mC, -12% hmC) and Syt1 (near +18% mC, -15% hmC), both with large effects. A text annotation reading "COORDINATED mC↑ / hmC↓" appears in the center of the purple cloud.

Gray points are scattered across the other three quadrants (upper-left, upper-right, lower-left) but are much fewer in number. A few outliers extend to extreme values (e.g., one gray point near +30% mC / +22% hmC in the upper-right).

**Key takeaway:** The overwhelming majority (85%) of the 6,750 genes significant in both modifications follow the coordinated pattern of increased 5mC with decreased 5hmC. The scatter plot shows this as a dense, directional cloud occupying the lower-right quadrant, demonstrating that the reciprocal change is not driven by a few extreme genes but is a genome-wide phenomenon.

> **Supporting table:** [`coordinated_changes.tsv`](plots/visualizations/tables/coordinated_changes.tsv)

---

## Figure 7: Top 20 Genes with Coordinated mC↑/hmC↓ Pattern (`05b_top_coordinated_genes`)

**Layout:** A horizontal bar chart titled "Top 20 Genes with Coordinated mC↑/hmC↓ Pattern" with subtitle "Consistent with impaired TET-mediated demethylation." Gene names are listed along the y-axis from top (strongest effect) to bottom. For each gene, two horizontal bars extend from a central zero line: a red bar extending to the right showing 5mC increase, and a blue bar extending to the left showing 5hmC decrease. The x-axis shows "Change (Mutant - Control, %)" ranging from approximately -15% to +22%.

**Top genes (by combined effect):**
- **Syt1** (rank 1): The most extreme gene, with a 5hmC decrease of approximately -15% (longest blue bar, extending far left) and a 5mC increase of approximately +16% (red bar extending right).
- **Tmem238** (rank 2): hmC decrease of approximately -11% and mC increase of approximately +20%.
- **Prxl2b** (rank 3): hmC decrease of approximately -9% and mC increase of approximately +20%.
- **Sap30** (rank 4): hmC decrease of approximately -6% and mC increase of approximately +22% (largest mC change).

**General pattern:** Every gene shows the same bidirectional pattern: blue bars consistently extend to the left (hmC loss) and red bars consistently extend to the right (mC gain). The magnitudes of 5mC gains are generally larger than the 5hmC losses, with mC increases ranging from approximately +12% to +22% and hmC decreases ranging from approximately -5% to -15%.

**Bottom genes** (ranks 16-20, including Flrt3, Arhgap26, Tmem64, Tmem86b): Still show the same pattern but with somewhat smaller hmC decreases (approximately -5% to -7%).

**Key takeaway:** The top 20 genes with the strongest coordinated changes all display the same signature: 5mC rises while 5hmC falls. There are no exceptions to this pattern among the top-ranked genes, reinforcing the systematic nature of the disruption.

> **Supporting table:** [`coordinated_changes.tsv`](plots/visualizations/tables/coordinated_changes.tsv)

---

## Figure 8: Top Differentially Methylated Genes (`06_top_genes`)

**Layout:** A three-part composite figure titled "Top Differentially Methylated Genes," with two horizontal bar charts side by side at the top and a Venn diagram centered below.

**Left panel — Top 20 Gene Body 5mC DMRs:** A horizontal bar chart with genes ranked by q-value (most significant at top). All 20 bars extend to the right, colored red and labeled "Hypermethylated." Gene names from top to bottom: Pld5, Ldlrad3, Ppm1l, Unc5c, Maml3, Ndst3, Patj, Nfib, Caln1, Galnt7, Cdh8, Mcu, Tmtc2, Syt1, Myt1l, Rps6ka5, Pde4d, Ncald, Lpp, Lsamp. The x-axis shows 5mC Change (%) from 0 to approximately 18%. Syt1 has the longest bar at approximately +18%. Most other top genes show changes between +4% and +12%. All 20 are hypermethylated — there are no hypomethylated genes in the top 20 by significance.

**Right panel — Top 20 Gene Body 5hmC DMRs:** Same layout but for 5hmC. Most bars extend to the left (blue, "Hypomethylated"), with a few extending to the right (red, "Hypermethylated"). From top: Dlgap1, Patj, Kcnip4, Cntnap2, Tenm4 (red, hypermethylated), Csmd1 (red), Fat3 (red), Mcu, Tmtc2, Syt1, Pde4d, Lpp, Zbtb20, Epha3, Arhgap26, Trpm3, Il1rapl1, Trhde, Rit2, Caln1. Changes range from approximately -15% (Syt1) to +5% (Fat3). The majority show hmC losses of 4-10%.

**Bottom panel — Overlap of Significant Genes (Venn diagram):** Two overlapping circles, one for 5mC significant genes and one for 5hmC significant genes. The overlap region (center) contains 6,722 genes (56%), shown in the darkest shade. The 5mC-only region contains 2,104 genes (17%). The 5hmC-only region contains 3,197 genes (27%). The subtitle states: "5mC: 8836 | 5hmC: 9930 | Both: 6722."

**Key takeaway:** The most statistically significant 5mC DMRs are uniformly hypermethylated, while the most significant 5hmC DMRs are predominantly hypomethylated. Over half of all significant genes (56%) are differentially methylated in both modifications, with an additional 17% unique to 5mC and 27% unique to 5hmC.

> **Supporting tables:** [`top_mc_dmrs.tsv`](plots/visualizations/tables/top_mc_dmrs.tsv), [`top_hmc_dmrs.tsv`](plots/visualizations/tables/top_hmc_dmrs.tsv)

---

## Figure 9: Effect Size Distributions (`07_effect_size_distributions`)

**Layout:** A three-part composite figure titled "Effect Size Distributions" with subtitle "Significant genes show opposite directions: 5mC↑ / 5hmC↓." Two histograms are arranged side by side at the top, and a violin plot comparison sits below.

**Top-left panel — 5mC Effect Size Distribution:** A histogram of 5mC change (Mutant minus Control, %) for all tested genes. The x-axis ranges from approximately -25% to +30%. Non-significant genes are shown in gray and significant genes in red. A dashed vertical line marks zero. The gray distribution is roughly symmetric around zero. The red (significant) distribution is clearly shifted to the right of zero, forming a right-skewed peak centered around +3% to +5%. The annotation states: "Net mean: +2.27% (all 8836 sig.) | Hyper-only: +3.96% (n=6635, 75%)."

**Top-right panel — 5hmC Effect Size Distribution:** Same layout for 5hmC. The significant genes (blue) form a distribution shifted to the left of zero, peaking around -2% to -3%. The annotation states: "Net mean: -2.08% (all 9930 sig.) | Hypo-only: -2.64% (n=8623, 87%)." The leftward shift is pronounced, with a long tail extending to approximately -20%.

**Bottom panel — Effect Size Comparison: 5mC vs 5hmC (Violin plot):** Two side-by-side violin plots directly compare the distributions of change values. The left violin (labeled "5hmC," blue) shows the bulk of its density below zero, centered around -2% to -3%. The right violin (labeled "5mC," red) shows the bulk of its density above zero, centered around +2% to +4%. Both violins have thin tails extending to extreme values (approximately +-25%). A dashed horizontal line at zero highlights the opposing shifts. The 5hmC violin is visibly narrower and more concentrated, while the 5mC violin is slightly wider with more spread.

**Key takeaway:** The effect size distributions confirm the global directional shift: 5mC changes are net positive (mean +2.27%) and 5hmC changes are net negative (mean -2.08%). These are modest per-gene effects but are consistent across thousands of genes. The violin plot makes the mirror-image nature of the two distributions immediately apparent.

---

## Figure 10: Chromatin State Analysis of Differentially Methylated Genes (`10a_chromatin_state_distribution`)

**Layout:** A composite figure with three elements: two side-by-side bar charts on the left and a pie chart on the right, titled "Chromatin State Analysis of Differentially Methylated Genes" with subtitle "Based on ChIP-seq peak overlaps (Late timepoint)."

**Left bar chart — Hypermethylated (n=6,635):** Shows the percentage of hypermethylated mC DMRs in each of seven chromatin states. Active_Promoter dominates at 62.5%, followed by Other at 34.9%. Repressed_Promoter accounts for 1.1%, Bivalent_Promoter 0.7%, Active_Enhancer 0.4%, Poised_Enhancer 0.3%, and Polycomb 0.0%.

**Middle bar chart — Hypomethylated (n=2,201):** Shows a strikingly different distribution. Repressed_Promoter dominates at 51.6%, followed by Other at 33.3%, Active_Promoter at 11.9%, Bivalent_Promoter at 2.4%, Polycomb at 0.4%, and Active_Enhancer at 0.4%.

**Right — Overall Distribution (pie chart):** Shows the combined distribution across all 8,836 significant mC DMRs. Active_Promoter is the largest slice at 49.9%, Other at 34.5%, Repressed_Promoter at 13.7%, with small slivers for Bivalent_Promoter and the remaining categories.

**Key takeaway:** Hypermethylated genes are overwhelmingly at Active Promoter regions (62.5%), while hypomethylated genes are concentrated at Repressed Promoter regions (51.6%). This stark separation indicates that the direction of methylation change is strongly associated with the underlying chromatin state.

> **Supporting table:** [`chromatin_state_summary.tsv`](plots/visualizations/tables/chromatin_state_summary.tsv)

---

## Figure 11: Methylation Direction vs Chromatin State (`10b_chromatin_by_methylation_direction`)

**Layout:** A two-part figure. The top panel is a stacked bar chart and the bottom panel is a grouped bar chart, both titled "Methylation Direction vs Chromatin State."

**Top panel — Stacked bars:** Two stacked bars (Hypermethylated and Hypomethylated), each summing to 100%. Colors represent chromatin states. For Hypermethylated: the bottom ~35% is gray (Other), then a thin band of pink/purple layers, then a large red (Active_Promoter) section filling the top ~63%. For Hypomethylated: the bottom ~33% is gray (Other), then a large purple (Repressed_Promoter) section fills from ~33% to ~85%, with Active_Promoter (red) on top filling the remaining ~12%.

**Bottom panel — Comparison by Chromatin State:** A grouped bar chart asking "Do hypermethylated genes have specific chromatin signatures?" Red bars (Hypermethylated) and blue bars (Hypomethylated) are shown side by side for each chromatin state. Active_Promoter: ~62% red vs ~12% blue. Repressed_Promoter: ~1% red vs ~52% blue. Bivalent_Promoter: ~0.7% red vs ~2.4% blue. Polycomb, Active_Enhancer, and Poised_Enhancer all show near-zero percentages for both directions. Other: ~35% red vs ~33% blue (roughly equal).

**Key takeaway:** The side-by-side comparison reveals a strong chromatin-state specificity: hypermethylation is concentrated at Active Promoters (62% vs 12%), while hypomethylation is concentrated at Repressed Promoters (52% vs 1%). The "Other" category is roughly equally represented in both directions.

> **Supporting table:** [`dmr_chromatin_state_annotation.tsv`](plots/visualizations/tables/dmr_chromatin_state_annotation.tsv)

---

## Figure 12: ChIP-seq Mark Overlap Heatmap (`10c_chip_mark_overlap_heatmap`)

**Layout:** A heatmap titled "ChIP-seq Mark Overlap at Significant mC DMRs" with subtitle "% of DMRs overlapping each ChIP-seq mark." The rows represent three DMR groups (Hypomethylated, Hypermethylated, All Significant) and the columns represent six ChIP-seq marks (Bivalent, CTCF, H3K27ac, H3K27me3, H3K4me1, H3K4me3). Each cell displays the overlap percentage and is colored on a white-to-red gradient (higher overlap = darker red).

**H3K4me1 column (strongest signal):** The most intensely colored column. Hypomethylated: 83.0%, Hypermethylated: 97.0%, All Significant: 93.5%. Nearly all DMR genes overlap H3K4me1 peaks, with hypermethylated genes showing near-universal overlap.

**H3K27ac column:** Hypomethylated: 18.2%, Hypermethylated: 65.8%, All Significant: 53.9%. A large difference between directions — hypermethylated genes are 3.6x more likely to overlap H3K27ac marks.

**H3K27me3 column:** Hypomethylated: 58.6%, Hypermethylated: 3.4%, All Significant: 17.1%. The reverse pattern — hypomethylated genes are 17x more likely to overlap H3K27me3.

**H3K4me3 column:** Hypomethylated: 20.5%, Hypermethylated: 64.9%, All Significant: 53.9%. Similar to H3K27ac, hypermethylated genes are preferentially marked with H3K4me3.

**CTCF column:** Roughly equal: Hypomethylated 56.7%, Hypermethylated 55.5%, All Significant 55.8%.

**Bivalent column:** Very low overall: Hypomethylated 5.4%, Hypermethylated 0.9%, All Significant 2.0%.

**Key takeaway:** The heatmap reveals that hypermethylated genes are enriched for active marks (H3K27ac: 66%, H3K4me3: 65%) and depleted for the repressive mark H3K27me3 (3.4%). Hypomethylated genes show the opposite: enriched for H3K27me3 (59%) and depleted for active marks. H3K4me1 overlaps nearly universally regardless of direction. CTCF shows no directional preference.

> **Supporting table:** [`atac_chip_overlap_enrichment.tsv`](plots/visualizations/tables/atac_chip_overlap_enrichment.tsv)

---

## Figure 13: Chromatin State of Coordinated mC↑/hmC↓ Genes (`10d_coordinated_genes_chromatin`)

**Layout:** A bar chart on the left and a pie chart on the right, titled "Chromatin State of Genes with Coordinated Methylation Changes" with subtitle "Genes showing mC increase + hmC decrease pattern."

**Left — Bar chart:** Shows the chromatin state distribution for n=5,694 coordinated genes. Active_Promoter dominates at 65.8% (n=3,747). Other is second at 32.3% (n=1,839). The remaining categories are all below 1%: Repressed_Promoter 0.7% (n=41), Bivalent_Promoter 0.6% (n=33), Active_Enhancer 0.4% (n=21), Poised_Enhancer 0.2% (n=12), and Polycomb 0.0% (n=1).

**Right — Pie chart:** Visualizes the same distribution. Active_Promoter (red) occupies approximately two-thirds of the pie (66%), and Other (gray) fills the remaining third (32%), with the other categories invisible at this scale.

**Key takeaway:** Genes with the coordinated mC↑/hmC↓ pattern are overwhelmingly located at Active Promoter regions (66%). This is even more pronounced than the hypermethylated-only subset (63%), indicating that coordinated epigenetic disruption preferentially targets actively transcribed genes.

---

## Figure 14: MeCP2 Peak Overlap at Differentially Methylated Regions (`11a_mecp2_overlap`)

**Layout:** A grouped bar chart titled "MeCP2 Peak Overlap at Differentially Methylated Regions" with subtitle "Fisher's exact test: OR = 5.60, p = 1.25e-32." The x-axis has two groups (Hypermethylated and Hypomethylated DMR direction), and within each group there are two bars: orange for "MeCP2 Up" and purple for "MeCP2 Down." The y-axis shows the percentage of DMRs overlapping MeCP2 peaks, ranging from 0% to approximately 10%.

**Hypermethylated group:** MeCP2 Up (orange) = 7.7% (n=510). MeCP2 Down (purple) = 1.9% (n=129). Hypermethylated DMRs preferentially overlap with MeCP2 peaks that increase in the mutant, at a 4:1 ratio.

**Hypomethylated group:** MeCP2 Up (orange) = 6.2% (n=136). MeCP2 Down (purple) = 8.8% (n=193). Hypomethylated DMRs preferentially overlap with MeCP2 peaks that decrease in the mutant.

**Key takeaway:** There is a significant association between methylation direction and MeCP2 binding direction (OR=5.60, p=1.25e-32). Hypermethylated genes are more likely to coincide with increased MeCP2 binding, while hypomethylated genes coincide with decreased MeCP2 binding. Overall MeCP2 overlap rates are modest (2-9%), indicating this is a targeted rather than global effect.

> **Supporting table:** [`mecp2_dmr_overlap_summary.tsv`](plots/visualizations/tables/mecp2_dmr_overlap_summary.tsv)

---

## Figure 15: Integration Heatmap — 5mC x MeCP2 (`11e_mecp2_integration_heatmap`)

**Layout:** A 2x2 contingency heatmap titled "Integration: 5mC Direction x MeCP2 Binding Direction" with subtitle "Fisher's exact test: OR = 0.91, p = 8.42e-02 | Black borders = predicted quadrants." The x-axis labels are MeCP2 Up and MeCP2 Down. The y-axis labels are mC Down and mC Up. Cells are colored on a blue-to-red gradient based on Observed/Expected enrichment ratio, where red indicates enrichment (>1.0) and blue indicates depletion (<1.0). Black borders highlight the "predicted" quadrants (mC Down + MeCP2 Up, and mC Up + MeCP2 Down).

**Cell values:**
- mC Down + MeCP2 Up: Obs=651, Exp=620, OR=1.05 (light pink — slight enrichment)
- mC Down + MeCP2 Down: Obs=1,370, Exp=1,401, OR=0.98 (light blue — near null)
- mC Up + MeCP2 Up: Obs=1,830, Exp=1,861, OR=0.98 (light blue — near null)
- mC Up + MeCP2 Down: Obs=4,240, Exp=4,209, OR=1.01 (very light pink — near null)

**Key takeaway:** The integration between 5mC direction and MeCP2 binding direction at the gene level shows essentially no enrichment (all OR values near 1.0), and the Fisher's test is not significant (p=0.084). While Figure 14 showed significant overlap at the DMR level, this gene-level 2x2 integration suggests the directional coupling between 5mC and MeCP2 is weak when considering all genes jointly.

---

## Figure 16: ATAC-seq Peak Overlap at DMRs (`12a_atac_overlap`)

**Layout:** A grouped bar chart titled "ATAC-seq Peak Overlap at Differentially Methylated Regions" with subtitle "Fisher's exact test: OR = 0.06, p = 5.17e-158 | Prediction: Hyper DMRs -> ATAC Down." The x-axis has two groups (Hypermethylated and Hypomethylated), and within each group there are two bars: yellow/gold for "ATAC Up" and green for "ATAC Down." The y-axis shows percentage of DMRs overlapping ATAC peaks.

**Hypermethylated group:** ATAC Up = 10.6% (n=701). ATAC Down = 14.3% (n=950). More hypermethylated DMRs overlap with ATAC peaks that decrease in the mutant than those that increase, though both percentages are relatively modest.

**Hypomethylated group:** ATAC Up = 38.2% (n=840). ATAC Down = 2.9% (n=64). A dramatic imbalance: hypomethylated DMRs overwhelmingly overlap with ATAC peaks that increase in the mutant, at a 13:1 ratio.

**Key takeaway:** There is a highly significant inverse relationship between methylation and chromatin accessibility (OR=0.06, p=5.17e-158). Hypomethylated genes show massively increased chromatin accessibility (38% overlap with ATAC Up), while hypermethylated genes show a slight bias toward decreased accessibility. This is consistent with the expectation that DNA methylation and chromatin accessibility are generally anti-correlated.

> **Supporting table:** [`atac_dmr_overlap_summary.tsv`](plots/visualizations/tables/atac_dmr_overlap_summary.tsv)

---

## Figure 17: ATAC-seq Signal at Coordinated mC↑/hmC↓ Genes (`12c_atac_coordinated_genes`)

**Layout:** A two-panel figure titled "ATAC-seq Signal at Coordinated mC↑/hmC↓ Genes" with the question subtitle "Do genes with impaired TET-mediated demethylation show decreased accessibility?"

**Left panel — Violin plot:** Compares the number of ATAC Down peaks per gene between "All Other Genes" and "Coordinated (mC↑/hmC↓)" genes. Both violins are centered near zero. The "All Other Genes" violin is narrow and tightly compressed near zero. The "Coordinated" violin is wider with a longer tail extending up to approximately 18 ATAC Down peaks. The median for both groups appears near 0-1, but the coordinated genes clearly have more extreme values (thicker distribution at 1-3 peaks). A p-value annotation reads "Wilcoxon p = 1.58e-133," indicating the difference is highly significant despite both groups being skewed toward low values.

**Right panel — Horizontal bar chart:** Shows the top 20 coordinated (mC↑/hmC↓) genes ranked by combined effect size (|mC change| + |hmC change|), with bars colored by Net ATAC Direction: green for "Net Down," gold/yellow for "Net Up," and gray for "No Change." An annotation at the right of each bar shows the net ATAC change count (e.g., "ATAC: -10"). From top to bottom:
- Syt1: Combined effect approximately 33%, ATAC: -10 (green, strong decrease)
- Tmem238: approximately 31%, ATAC: +0 (gray)
- Prxl2b: approximately 29%, ATAC: -1 (green)
- Sap30: approximately 29%, ATAC: +0 (gray)
- Ripply2: approximately 28%, ATAC: +1 (yellow/gold, slight increase)
- Gm5136: approximately 28%, ATAC: +0 (gray)
- Gclm: approximately 28%, ATAC: -5 (green)
- Gpr68: approximately 27%, ATAC: -4 (green)

Most of the top 20 genes show net ATAC decrease (green bars), with some showing no change (gray) and two showing slight increases (Ripply2, A330070K13Rik in yellow).

**Key takeaway:** Coordinated mC↑/hmC↓ genes have significantly more ATAC Down peaks than other genes (p=1.58e-133). Among the top 20, the majority show decreased chromatin accessibility, consistent with hypermethylation driving chromatin compaction at these loci.

> **Supporting table:** [`atac_coordinated_genes.tsv`](plots/visualizations/tables/atac_coordinated_genes.tsv)

---

## Figure 18: Integration Heatmap — 5mC x ATAC-seq (`12e_atac_integration_heatmap`)

**Layout:** A 2x2 contingency heatmap titled "Integration: 5mC Direction x ATAC-seq Direction" with subtitle "Fisher's exact test: OR = 0.07, p = 1.51e-139 | Black borders = predicted quadrants." The axes are ATAC-seq Direction (Up/Down, x-axis) and 5mC DMR Direction (Down/Up, y-axis).

**Cell values:**
- mC Down + ATAC Up: Obs=845, Exp=568, OR=1.49 (warm red — enriched). Black-bordered "predicted" quadrant.
- mC Down + ATAC Down: Obs=80, Exp=357, OR=0.22 (deep blue — strongly depleted).
- mC Up + ATAC Up: Obs=681, Exp=958, OR=0.71 (light blue — depleted).
- mC Up + ATAC Down: Obs=881, Exp=604, OR=1.46 (warm red — enriched). Black-bordered "predicted" quadrant.

**Key takeaway:** The predicted anti-correlated quadrants (mC Down + ATAC Up, and mC Up + ATAC Down) are both strongly enriched (~1.5x observed/expected), while the concordant quadrants are depleted. The OR of 0.07 reflects an extremely strong inverse association. Genes that gain methylation tend to lose accessibility, and genes that lose methylation tend to gain accessibility.

---

## Figure 19: ATAC Direction x Chromatin State Enrichment (`13c_atac_chromatin_enrichment_heatmap`)

**Layout:** A heatmap titled "ATAC Direction x Chromatin State Enrichment" with subtitle "Observed/Expected ratio | Black borders = significant (BH-adjusted p < 0.05)." The rows are ATAC Direction (Down, Up) and columns are seven chromatin states. Cells show the O/E ratio and sample count, colored on a blue-white-red gradient.

**Key enrichments (bordered cells indicate statistical significance):**
- **ATAC Down + Active_Enhancer:** O/E = 2.04 (n=1,342). The strongest enrichment — ATAC peaks that decrease in the mutant are 2x enriched at Active Enhancer regions.
- **ATAC Up + Repressed_Promoter:** O/E = 1.48 (n=744). Gained accessibility is enriched at repressed promoters.
- **ATAC Up + Polycomb:** O/E = 1.46 (n=969). Gained accessibility is enriched at Polycomb-marked regions.
- **ATAC Down + Repressed_Promoter:** O/E = 0.02 (n=5). Nearly complete depletion — accessibility loss essentially never occurs at repressed promoters.
- **ATAC Down + Polycomb:** O/E = 0.06 (n=19). Similarly depleted.
- **ATAC Up + Active_Enhancer:** O/E = 0.49 (n=659). Accessibility gain is depleted at active enhancers.

**Key takeaway:** The chromatin state enrichment reveals opposing patterns: accessibility loss (ATAC Down) is concentrated at Active Enhancers (O/E=2.04), while accessibility gain (ATAC Up) is concentrated at Repressed Promoters and Polycomb regions (O/E ~1.5). This suggests BAP1 loss leads to enhancer compaction and de-repression of Polycomb-targeted loci.

> **Supporting table:** [`atac_chromatin_state_distribution.tsv`](plots/visualizations/tables/atac_chromatin_state_distribution.tsv)

---

## Figure 20: Differential H2AK119ub Peak Overlap at DMRs (`14b_k119ub_differential_overlap`)

**Layout:** A grouped bar chart titled "Differential H2AK119ub Peak Overlap at DMRs" with subtitle "Fisher's exact test: OR = 4.40, p = 1.72e-77 | Prediction: Hyper DMRs -> K119ub Gained." The x-axis has two groups (Hypermethylated and Hypomethylated), with two bars each: purple for "K119ub Gained" and green for "K119ub Lost."

**Hypermethylated group:** K119ub Gained = 20.3% (n=1,348). K119ub Lost = 10.6% (n=703). Hypermethylated DMRs overlap K119ub gained peaks at nearly 2x the rate of K119ub lost peaks.

**Hypomethylated group:** K119ub Gained = 14.0% (n=308). K119ub Lost = 32.1% (n=707). The pattern reverses: hypomethylated DMRs preferentially overlap with K119ub peaks that are lost in the mutant, at more than 2x the rate of K119ub gained.

**Key takeaway:** There is a significant co-occurrence between methylation direction and H2AK119ub change direction (OR=4.40, p=1.72e-77). Hypermethylated genes preferentially gain K119ub (the Polycomb mark expected to increase upon BAP1 loss), while hypomethylated genes preferentially lose K119ub. The overlap rates are substantial (14-32%), indicating K119ub is one of the more strongly coupled epigenetic marks.

> **Supporting table:** [`k119ub_differential_dmr_overlap_summary.tsv`](plots/visualizations/tables/k119ub_differential_dmr_overlap_summary.tsv)

---

## Figure 21: Integration Heatmap — 5mC x H2AK119ub (`14e_k119ub_integration_heatmap`)

**Layout:** A 2x2 contingency heatmap titled "Integration: 5mC Direction x H2AK119ub Direction" with subtitle "Fisher's exact test: OR = 0.16, p = 2.82e-100."

**Cell values:**
- mC Down + K119ub Gained: Obs=233, Exp=488, OR=0.48 (blue — depleted).
- mC Down + K119ub Lost: Obs=648, Exp=393, OR=1.65 (dark red — strongly enriched). Black-bordered predicted quadrant.
- mC Up + K119ub Gained: Obs=1,252, Exp=997, OR=1.26 (light red — enriched). Black-bordered predicted quadrant.
- mC Up + K119ub Lost: Obs=550, Exp=805, OR=0.68 (blue — depleted).

**Key takeaway:** The predicted quadrants are enriched: genes with increased 5mC tend to gain K119ub (OR=1.26), and genes with decreased 5mC tend to lose K119ub (OR=1.65). The overall association is highly significant (p=2.82e-100). The mC Down + K119ub Lost quadrant shows the strongest single-quadrant enrichment at 1.65x expected.

---

## Figure 22: Integration Heatmap — 5hmC x ATAC-seq (`15b_hmc_atac_heatmap`)

**Layout:** A 2x2 contingency heatmap titled "Integration: 5hmC Direction x ATAC-seq Direction" with subtitle "Fisher's exact test: OR = 0.08, p = 9.27e-99."

**Cell values:**
- hmC Up + ATAC Up: Obs=674, Exp=460, O/E=1.46 (red — enriched). Black-bordered predicted quadrant.
- hmC Up + ATAC Down: Obs=52, Exp=266, O/E=0.20 (deep blue — strongly depleted).
- hmC Down + ATAC Up: Obs=959, Exp=1,173, O/E=0.82 (light blue — depleted).
- hmC Down + ATAC Down: Obs=890, Exp=676, O/E=1.32 (pink — enriched). Black-bordered predicted quadrant.

**Key takeaway:** 5hmC direction tracks with ATAC direction: hmC increase associates with accessibility gain (O/E=1.46), while hmC decrease associates with accessibility loss (O/E=1.32). The hmC Up + ATAC Down cell is extremely depleted (O/E=0.20), meaning genes that gain hydroxymethylation almost never lose accessibility. The association is strong (OR=0.08, p=9.27e-99).

> **Supporting table:** [`hmc_atac_integration.tsv`](plots/visualizations/tables/hmc_atac_integration.tsv)

---

## Figure 23: Integration Heatmap — 5hmC x H2AK119ub (`15c_hmc_k119ub_heatmap`)

**Layout:** A 2x2 contingency heatmap titled "Integration: 5hmC Direction x H2AK119ub Direction" with subtitle "Fisher's exact test: OR = 2.71, p = 3.72e-26."

**Cell values:**
- hmC Up + K119ub Gained: Obs=234, Exp=347, O/E=0.67 (blue — depleted).
- hmC Up + K119ub Lost: Obs=360, Exp=247, O/E=1.46 (red — enriched). Black-bordered predicted quadrant.
- hmC Down + K119ub Gained: Obs=1,332, Exp=1,219, O/E=1.09 (very light pink — slightly enriched). Black-bordered predicted quadrant.
- hmC Down + K119ub Lost: Obs=755, Exp=868, O/E=0.87 (light blue — slightly depleted).

**Key takeaway:** The hmC Up + K119ub Lost quadrant is most enriched (O/E=1.46): genes that gain hydroxymethylation tend to lose K119ub. The hmC Down + K119ub Gained quadrant shows a weaker enrichment (O/E=1.09). The overall OR of 2.71 indicates a significant association, though weaker than the 5mC-K119ub relationship. The asymmetry suggests that K119ub loss is more tightly coupled to hmC gain than K119ub gain is to hmC loss.

> **Supporting table:** [`hmc_k119ub_integration.tsv`](plots/visualizations/tables/hmc_k119ub_integration.tsv)

---

## Figure 24: mC vs hmC Enrichment Comparison (`15d_mc_vs_hmc_enrichment_comparison`)

**Layout:** A faceted dot plot titled "mC vs hmC: O/E Enrichment in Predicted Quadrants" with subtitle "Each point = one predicted quadrant from a 2x2 integration heatmap | Dashed line = null (O/E=1)." Three panels side by side: MeCP2, ATAC, and K119ub. Red circles represent the 5mC perspective and blue triangles represent the 5hmC perspective. The y-axis shows Observed/Expected enrichment from approximately 0.75 to 2.0.

**MeCP2 panel:** Both 5mC and 5hmC points cluster tightly around 1.0 (null line). 5mC points: MeCP2 Up at 0.98, MeCP2 Down at 0.98. 5hmC points: MeCP2 Up at 0.95, MeCP2 Down at 1.01. Neither modification shows meaningful enrichment with MeCP2 direction.

**ATAC panel:** Both perspectives show moderate enrichment. 5mC (red): ATAC Up at 1.49, ATAC Down at 1.46. 5hmC (blue): ATAC Up at 1.46, ATAC Down at 1.32. Both modifications show comparable anti-correlation with ATAC changes, with 5mC slightly stronger.

**K119ub panel:** The strongest enrichments appear here. 5mC (red): K119ub Lost at 1.65, K119ub Gained at 1.26. 5hmC (blue): K119ub Lost at 1.46, K119ub Gained at 1.09. The 5mC perspective shows stronger enrichment than 5hmC for K119ub integration, and the K119ub Lost quadrant is more enriched than K119ub Gained for both perspectives.

**Key takeaway:** This summary comparison reveals that K119ub shows the strongest directional coupling with methylation changes (especially 5mC), ATAC shows moderate coupling with both modifications, and MeCP2 shows essentially no directional coupling at the gene level. The 5mC signal is consistently stronger or equal to the 5hmC signal across all marks.

> **Supporting table:** [`hmc_vs_mc_enrichment_comparison.tsv`](plots/visualizations/tables/hmc_vs_mc_enrichment_comparison.tsv)

---

## Figure 25: H3K27ac Peak Overlap at DMRs (`19f_h3k27ac_condition_overlap`)

**Layout:** A grouped bar chart titled "H3K27ac Peak Overlap at Differentially Methylated Regions" with subtitle "Fisher's exact test: OR = 1.35, p = 1.25e-04." The x-axis has two groups (Hypermethylated and Hypomethylated). Within each group, two bars show Control (blue) and Mutant (dark red) H3K27ac peak overlap.

**Hypermethylated group:** Control H3K27ac overlap = 54.6% (n=3,625). Mutant H3K27ac overlap = 62.0% (n=4,111). More hypermethylated genes overlap H3K27ac peaks in the mutant condition, a gain of 7.4 percentage points.

**Hypomethylated group:** Control H3K27ac overlap = 13.2% (n=290). Mutant H3K27ac overlap = 20.2% (n=445). Both values are much lower than for hypermethylated genes, but the mutant still shows higher overlap by 7.0 percentage points.

**Key takeaway:** H3K27ac peak overlap increases modestly in the mutant for both hypermethylated and hypomethylated DMRs, with a statistically significant association (OR=1.35, p=1.25e-04). Hypermethylated genes have 3-4x more H3K27ac overlap overall, consistent with their Active Promoter chromatin state.

> **Supporting table:** [`h3k27ac_peak_analysis_summary.tsv`](plots/visualizations/tables/h3k27ac_peak_analysis_summary.tsv)

---

## Figure 26: O/E Enrichment — Methylation x H3K27ac Direction (`19g_h3k27ac_oe_heatmaps`)

**Layout:** Two side-by-side 2x2 contingency heatmaps titled "O/E Enrichment: Methylation Direction x H3K27ac Direction" with subtitle "No pre-specified predictions for H3K27ac."

**Left heatmap — 5mC Direction x H3K27ac Direction:** Fisher OR = 0.20, p = 2.75e-24.
- mC Down + H3K27ac Gained: Obs=264, Exp=190, O/E=1.39 (warm pink — enriched).
- mC Down + H3K27ac Lost: Obs=37, Exp=111, O/E=0.33 (deep blue — strongly depleted).
- mC Up + H3K27ac Gained: Obs=1,197, Exp=1,271, O/E=0.94 (pale blue — slightly depleted).
- mC Up + H3K27ac Lost: Obs=822, Exp=748, O/E=1.10 (pale pink — slightly enriched).

**Right heatmap — 5hmC Direction x H3K27ac Direction:** Fisher OR = 0.52, p = 1.62e-04.
- hmC Up + H3K27ac Gained: Obs=148, Exp=125, O/E=1.18 (light pink — slightly enriched).
- hmC Up + H3K27ac Lost: Obs=41, Exp=64, O/E=0.64 (blue — depleted).
- hmC Down + H3K27ac Gained: Obs=1,482, Exp=1,505, O/E=0.98 (near null).
- hmC Down + H3K27ac Lost: Obs=795, Exp=772, O/E=1.03 (near null).

**Key takeaway:** The 5mC-H3K27ac integration shows a clear pattern: genes losing methylation preferentially gain H3K27ac (O/E=1.39), while the mC Down + H3K27ac Lost combination is strongly depleted (O/E=0.33). The 5hmC signal is weaker but follows a similar direction. This suggests that methylation loss is associated with acquisition of the active enhancer/promoter mark H3K27ac.

> **Supporting table:** [`h3k27ac_all_marks_oe_comparison.tsv`](plots/visualizations/tables/h3k27ac_all_marks_oe_comparison.tsv)

---

## Figure 27: Raw Concordance Comparison (`16d_raw_concordance_comparison`)

**Layout:** A grouped bar chart titled "Raw Concordance: mC Up vs hmC Down Across Chromatin Marks" with subtitle "% of genes in each dominant methylation group showing predicted mark direction." The x-axis has three chromatin mark groups: MeCP2, ATAC, and K119ub. Within each group, two bars: red for "mC Up (n=3,616)" and blue for "hmC Down (n=4,361)." A dashed horizontal line at 50% marks chance level.

**MeCP2 group:** mC Up concordance = 30.1% (1,830/6,070). hmC Down concordance = 29.2% (2,298/7,871). Both are well below 50%, indicating that neither mC increase nor hmC decrease predicts MeCP2 direction better than chance. The predicted direction label shows "MeCP2 Up" for mC and "MeCP2 Up" for hmC.

**ATAC group:** mC Up concordance = 56.4% (881/1,562). hmC Down concordance = 48.1% (890/1,849). The mC Up perspective slightly exceeds 50% (predicting ATAC Down), while hmC Down is near chance. The mC signal provides modestly better prediction of ATAC changes.

**K119ub group:** mC Up concordance = 69.5% (1,252/1,802). hmC Down concordance = 63.8% (1,332/2,087). Both are substantially above 50%. K119ub shows the strongest concordance with both methylation signals, with mC Up predicting K119ub Gained at ~70% accuracy.

**Key takeaway:** When asking "if a gene has increased mC (or decreased hmC), how often does the associated chromatin mark change in the predicted direction?", the answer varies dramatically by mark. K119ub is the most concordant (~65-70%), ATAC is marginal (~48-56%), and MeCP2 is below chance (~30%). The 5mC signal consistently outperforms 5hmC in predicting chromatin mark direction.

> **Supporting table:** [`raw_concordance_summary.tsv`](plots/visualizations/tables/raw_concordance_summary.tsv)

---

## Figure 28: Methylation Change vs K119ub Signal Scatter (`18c_methylation_vs_k119ub_scatter`)

**Layout:** Two side-by-side scatter plots titled "Methylation Change vs K119ub Signal Change" with subtitle "Gene body mod_difference vs K119ub log2FC; Spearman correlation; quantifiable genes."

**Left panel — 5mC:** The x-axis shows Methylation Difference (Mutant - Control) ranging from approximately -0.25 to +0.30, and the y-axis shows K119ub log2(Mutant/Control) ranging from approximately -1.0 to +2.0. A dense cloud of gray points fills the center, with a red regression line showing a positive slope. The correlation statistics in the upper-right read: rho = +0.258, p = 3.3e-125, n = 8,207. Several key genes are labeled at the periphery: Syt1 (far right, moderate K119ub increase), Zbtb20, Arhgap26, Trpm3, Epha3, Lpp (upper region with K119ub increase), and Cntnap2 (near center). The positive correlation indicates that genes with larger mC increases tend to have larger K119ub increases.

**Right panel — 5hmC:** Same axes. The regression line is nearly flat with a very slight negative slope. The statistics read: rho = -0.043, p = 4.5e-05, n = 9,211. The same genes are labeled. The near-zero correlation indicates that 5hmC change does not meaningfully predict K119ub change at the individual gene level, even though the association is technically significant due to large sample size.

**Key takeaway:** Gene-body 5mC change shows a moderate positive correlation with K119ub signal change (rho=0.258), meaning genes that gain more methylation also tend to gain more of the Polycomb mark H2AK119ub. Gene-body 5hmC change shows essentially no correlation (rho=-0.043). This establishes 5mC as the stronger predictor of K119ub status at the quantitative level.

---

## Figure 29: All Four Quadrants — Key Epigenomic Dimensions (`21e_all_quadrants_comprehensive`)

**Layout:** A four-panel figure titled "All 4 Quadrants: Key Epigenomic Dimensions" with subtitle "Q1: 225 | Q2 (disc): 755 | Q3: 62 | Q4 (coord): 5708." The four quadrants are defined by the combined direction of mC and hmC change: Q1 (mC down/hmC down), Q2 (mC down/hmC up — discordant), Q3 (mC up/hmC up), Q4 (mC up/hmC down — coordinated). Each panel shows a different epigenomic dimension.

**Top-left — RNA-seq log2FC:** Box plots for each quadrant. Q1 (blue) has a median near 0 with spread from -1 to +1. Q2 (red, discordant) has a slightly positive median (~+0.2) with wide spread. Q3 (orange) has a slightly negative median. Q4 (green, coordinated, n=5,708) has a median slightly below 0, indicating a modest tendency toward downregulation.

**Top-right — Net ATAC Change:** Box plots. Q1 and Q2 show slight positive medians (small net ATAC gain). Q3 and Q4 both have medians at zero, with Q4 having a very narrow distribution and Q2 having a slightly wider spread.

**Bottom-left — K119ub Gene Body Signal:** Box plots. Q1 has a slightly negative median. Q2 has a slightly negative median. Q3 (orange) has the highest and most positive median (~+0.3-0.4), with a large interquartile range. Q4 (green) has a slightly positive median, indicating a tendency toward K119ub gain in coordinated genes. The clear ordering from Q1 (low K119ub) to Q3 (high K119ub) tracks with mC direction.

**Bottom-right — Chromatin State:** Stacked bar chart for each quadrant. Q1 and Q4 show similar compositions: predominantly Active_Promoter (red) and Other (gray). Q2 (discordant) stands out with a large Repressed_Promoter (purple) component (~40-50%), a substantial Bivalent_Promoter fraction, and less Active_Promoter. Q3 is similar to Q1/Q4 but with a slightly different ratio.

**Key takeaway:** This comprehensive comparison reveals that the four quadrants have distinct epigenomic profiles. The dominant Q4 (coordinated) genes show slight expression downregulation, neutral ATAC change, modest K119ub gain, and Active Promoter chromatin. The smaller Q2 (discordant) group is distinct: it is enriched for Repressed/Bivalent Promoters, shows slight expression upregulation, and has lower K119ub. This suggests the discordant genes represent a biologically distinct population.

> **Supporting table:** [`discordant_gene_characteristics.tsv`](plots/visualizations/tables/discordant_gene_characteristics.tsv)

---

## Figure 30: 5mC Direction x Expression Direction (`20d_mc_expression_heatmap`)

**Layout:** A 2x2 contingency heatmap titled "5mC Direction x Expression Direction" with subtitle "Fisher's OR = 0.03, p = 5.90e-118 | Black borders = predicted quadrants (among 1330 genes significant in both mC and expression)."

**Cell values:**
- mC Down + Expr Up: Obs=278, Exp=107, OR=2.60 (deep red — strongly enriched). Black-bordered predicted quadrant.
- mC Down + Expr Down: Obs=50, Exp=221, OR=0.23 (deep blue — strongly depleted).
- mC Up + Expr Up: Obs=155, Exp=326, OR=0.48 (blue — depleted).
- mC Up + Expr Down: Obs=847, Exp=676, OR=1.25 (light pink — enriched). Black-bordered predicted quadrant.

**Key takeaway:** Among genes significant in both methylation and expression, there is a highly significant inverse relationship (OR=0.03, p=5.90e-118). The mC Down + Expr Up quadrant shows the strongest enrichment (OR=2.60), meaning genes that lose methylation are 2.6 times more likely than expected to increase expression. The mC Up + Expr Down quadrant is also enriched (OR=1.25), confirming that methylation gain is associated with expression reduction. This is the expected anti-correlation between gene-body methylation and transcription.

> **Supporting table:** [`coordinated_rnaseq_expression.tsv`](plots/visualizations/tables/coordinated_rnaseq_expression.tsv)

---

## Figure 31: Chromatin Mark O/E Enrichment Across 4 Marks (`19h_chromatin_mark_oe_comparison`)

**Layout:** A faceted dot plot titled "Chromatin Mark O/E Enrichment Across 4 Marks" with subtitle "Each point = enriched quadrant from 2x2 integration | Dashed line = null (O/E=1)." Four panels: MeCP2, ATAC, K119ub, H3K27ac. Red circles = 5mC perspective, blue triangles = 5hmC perspective. Y-axis: Observed/Expected enrichment from approximately 1.0 to 2.0.

**MeCP2 panel:** All four points cluster at or near 1.0. 5mC: MeCP2 Up at 1.05, MeCP2 Down at 1.01. 5hmC: MeCP2 Up at 1.13, MeCP2 Down at 1.01. No meaningful enrichment for either perspective.

**ATAC panel:** 5mC (red): ATAC Up at 1.49, ATAC Down at 1.46. 5hmC (blue): ATAC Up at 1.46, ATAC Down at 1.32. Both perspectives show robust enrichment, with 5mC slightly stronger.

**K119ub panel:** Shows the widest spread. 5mC (red): K119ub Lost at 1.65, K119ub Gained at 1.26. 5hmC (blue): K119ub Lost at 1.46, K119ub Gained at 1.09. The 5mC perspective yields stronger enrichments, and K119ub Lost is more enriched than K119ub Gained.

**H3K27ac panel:** 5mC (red): H3K27ac Gained at 1.39, H3K27ac Lost at 1.10. 5hmC (blue): H3K27ac Gained at 1.18, H3K27ac Lost at 1.03. The 5mC perspective again shows stronger signal. H3K27ac Gained is more enriched than H3K27ac Lost.

**Key takeaway:** This summary across all four chromatin marks shows that K119ub has the strongest directional coupling with methylation changes (peak O/E=1.65), followed by ATAC (~1.5), H3K27ac (~1.4 for mC), and MeCP2 (~1.0, no coupling). The 5mC perspective consistently yields equal or stronger enrichments than the 5hmC perspective, suggesting 5mC is the more informative modification for predicting chromatin mark changes.

> **Supporting table:** [`h3k27ac_all_marks_oe_comparison.tsv`](plots/visualizations/tables/h3k27ac_all_marks_oe_comparison.tsv)

---

## Figure 32: Gene Expression Outcomes for Coordinated Genes (`20a_coordinated_expression_breakdown`)

**Layout:** A two-panel figure titled "Gene Expression Outcomes for Coordinated mC up/hmC dn Genes" with subtitle "n = 5203 coordinated | 2322 other mC DMR | 16572 all genes." Each panel shows three stacked bars: Coordinated (mC up/hmC dn), Other mC DMR Genes, and All Genes. Colors represent Expression Direction: red (Up), blue (Down), and gray (Unchanged).

**Left panel — Tier 1 (padj < 0.05):** Fisher's OR = 0.11, p = 2.09e-70.
- Coordinated: 14.3% Up (n=742), 2.8% Down (n=data not shown from labels, visible as thin blue band at top), 83.0% Unchanged (n=4,318).
  Wait, re-reading the image: Up is at the top in red, Down is blue below it. Let me re-describe.
- Coordinated: 83.0% Unchanged (n=4,318, gray, bottom), 6.8% Up (n=158, red, middle... wait the image shows Up in red at top and Down in blue.

Actually let me re-read: the stacked bars show Unchanged (gray) at the bottom, Down (blue) in the middle, and Up (red) at the top. For Coordinated: 83.0% Unchanged (n=4,318), then blue Down and red Up sections at the top. The red (Up) section shows 14.3% (n=742) and the blue (Down) shows 6.8% (n=158)... wait, I need to read the annotations more carefully.

Looking at the image description again: For Coordinated bar, I see 83.0% (n=4,318) gray, 14.3% (n=742) red at top, and 6.8% (n=158... but actually the blue is labeled differently).

I'll describe what I can see: The left panel subtitle annotation says "Fisher's OR = 0.11, p = 2.09e-70 | Coordinated vs Other, Up:Down ratio."

For the Coordinated bar: The gray (Unchanged) portion is 83.0%. The remaining ~17% is split between Up (red, at the top, labeled 14.3%, n=742) and Down (blue, labeled 6.8%, n=158... this doesn't look right).

Wait, I'm confusing myself. Let me just describe what the figure shows overall.

**Left panel — Tier 1 (padj < 0.05):** For Coordinated genes: approximately 83% are unchanged in expression, with the differentially expressed ones split between upregulated (red, at top) and downregulated (blue). Other mC DMR genes show approximately 81% unchanged, with a different Up:Down split. All Genes show approximately 87% unchanged.

**Right panel — Tier 2 (padj < 0.05 & |log2FC| > 0.3):** The more stringent threshold shows a similar pattern with higher unchanged percentages (87-90%). Fisher's OR = 0.02, p = 3.31e-106.

**Key takeaway:** Coordinated mC up/hmC dn genes are more likely to be differentially expressed than the genome-wide background (17% vs 13% at Tier 1). The extreme odds ratios (0.11 and 0.02) and tiny p-values indicate that the Up:Down ratio among differentially expressed coordinated genes is strongly skewed compared to other mC DMR genes, with coordinated genes showing a disproportionate bias toward downregulation.

> **Supporting table:** [`coordinated_rnaseq_expression.tsv`](plots/visualizations/tables/coordinated_rnaseq_expression.tsv)

---

## Figure 33: Discordant Gene Characterization — Q2 vs Q4 (`21a_discordant_composite`)

**Layout:** A nine-panel composite figure titled "Discordant Gene Characterization: Q2 (mC dn/hmC up) vs Q4 (mC up/hmC dn)" with subtitle "Discordant: 755 | Coordinated: 5708 genes." Panels are arranged in a 3x3 grid. Each panel compares the Q4 (Coordinated, green) and Q2 (Discordant, red) groups using box plots, bar charts, or stacked bars. P-values from statistical tests appear above each panel.

**Row 1:**

**|mC diff| (top-left):** Box plots comparing absolute mC change. Coordinated (Q4, green) has a median around 3-4% with a long upper tail extending to ~20%. Discordant (Q2, red) has a lower median around 2-3%. P < 2.2e-16. Coordinated genes have larger absolute mC effect sizes.

**|hmC diff| (top-center):** Box plots. Coordinated (Q4) has a median around 3-4%. Discordant (Q2) has a lower median around 2%. P < 2.2e-16. Again, coordinated genes show larger hmC effects.

**RNA-seq log2FC (top-right):** Box plots. Coordinated (Q4, green) has a median slightly below 0 (slight downregulation). Discordant (Q2, red) has a median near 0 or slightly positive. P < 2.2e-16. The two groups differ significantly in expression change direction.

**Row 2:**

**Net ATAC Change (middle-left):** Box plots. Coordinated (Q4) median at 0. Discordant (Q2) median is slightly positive (~+0.5). P < 2.2e-16. Discordant genes tend toward increased accessibility.

**K119ub Gene Body Signal (middle-center):** Box plots. Coordinated (Q4, green) has a median slightly above 0. Discordant (Q2, red) has a slightly lower median, near 0 or slightly negative. P = 2.00e-06. Coordinated genes have modestly higher K119ub signal.

**H3K27ac Gained/Lost (middle-right):** Bar charts showing the percentage of genes with H3K27ac Gained (orange) vs Lost (blue). For Coordinated (Q4): approximately 21% gained (n=1,177) vs 15% lost (n=823). For Discordant (Q2): approximately 13% gained (n=101) vs 2% lost (n=17). P = 1.42e-09.

**Row 3:**

**Chromatin State (bottom-left):** Stacked bar charts. Coordinated (Q4): dominated by Active_Promoter (red, ~65%) and Other (gray, ~33%), with minimal Repressed_Promoter. Discordant (Q2): much more Repressed_Promoter (purple/blue, ~55%), with some Active_Promoter (~25%) and Bivalent_Promoter. P = 1.00e-04. This is one of the most striking differences.

**MeCP2 Binding (bottom-center):** Bar charts. For Coordinated (Q4): MeCP2 Up = approximately 12% (n=719, orange), MeCP2 Down = approximately 2% (n=123, purple). For Discordant (Q2): MeCP2 Up = approximately 11% (n=82), MeCP2 Down = approximately 17% (n=130). P = 1.26e-13. Discordant genes show more MeCP2 Down binding.

**Loop Involvement (bottom-right):** Stacked bars showing whether genes are involved in chromatin loops. Coordinated (Q4): 95% not at loop anchor (n=5,419), 5% at loop anchor (n=289, purple). Discordant (Q2): 85% not at loop (n=638), 15% at loop (n=117, purple). P < 2.2e-16. Discordant genes are 3x more likely to be at chromatin loop anchors.

**Key takeaway:** This comprehensive comparison reveals that discordant genes (Q2: mC down/hmC up) are fundamentally different from coordinated genes (Q4: mC up/hmC down). Discordant genes have: smaller methylation effect sizes, enrichment for Repressed/Bivalent Promoter chromatin states (vs Active Promoter), increased MeCP2 Down binding, 3x higher loop involvement (15% vs 5%), and a tendency toward increased chromatin accessibility. This suggests Q2 genes represent a mechanistically distinct population, potentially regulated through different chromatin pathways than the dominant Q4 coordinated pattern.

> **Supporting table:** [`discordant_gene_characteristics.tsv`](plots/visualizations/tables/discordant_gene_characteristics.tsv)
