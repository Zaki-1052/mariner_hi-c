## Section 10f: Expanded 10-Mark Chromatin Overlap Heatmap
**Key numbers:**
- H3K4me1 overlaps 91.0% and ATAC 90.3% of all significant DMRs (highest) (source: 10f_expanded_chip_overlap_percentages.tsv)
- H3K27ac: hyper = 64.0% vs hypo = 15.9%; H3K27me3: hyper = 3.6% vs hypo = 55.4% (source: 10f_expanded_chip_overlap_percentages.tsv)
- H2AK119ub: all = 46.2%, hyper = 38.7%, hypo = 63.4% (source: 10f_expanded_chip_overlap_percentages.tsv)
- H3K36me2 = 11.0% and H3K36me3 = 7.3% overall (lowest) (source: 10f_expanded_chip_overlap_percentages.tsv)

**What this shows:** Extends the 6-mark heatmap to 10 marks (adds H2AK119ub, ATAC, H3K36me2/me3 as condition-union peaks). The mark profile separates the two DMR classes cleanly: hypermethylated DMRs carry active marks (H3K27ac, H3K4me3, H3K4me1), hypomethylated DMRs carry repressive marks (H3K27me3, K119ub).

**Figures:**
- `10f_expanded_chip_heatmap` — 10-mark × 3-row (All/Hyper/Hypo) overlap % heatmap
