## Section 12: ATAC-seq Accessibility Correlation at DMRs
**Key numbers:**
- DMR-overlap Fisher OR = 0.068, p = 4.44e-178 (strong directional coupling) (source: atac_dmr_overlap_summary.tsv)
- Hyper DMRs: 13.4% overlap ATAC-Down (1,006/7,513) vs 10.2% ATAC-Up (source: atac_dmr_overlap_summary.tsv)
- Hypo DMRs: 33.4% overlap ATAC-Up (1,088/3,262) vs 2.9% ATAC-Down (source: atac_dmr_overlap_summary.tsv)
- 6,589 coordinated mC↑/hmC↓ genes annotated (top genes trend ATAC-down, e.g. Syt1 net=−10) (source: atac_coordinated_genes.tsv)

**What this shows:** Accessibility is the most cleanly coupled mark at the DMR level — hypermethylated loci close (ATAC-down) and hypomethylated loci open (ATAC-up), as predicted if K119ub-driven compaction blocks TET. The extreme p-value reflects a large directional asymmetry.

**Figures:**
- `12a_atac_overlap` — ATAC-Up/Down overlap % by DMR direction + Fisher
- `12b_consensus_accessibility` — consensus ctrl-vs-mut accessibility at DMRs
- `12c_atac_coordinated_genes` — ATAC-down box + top-20 coordinated bar
- `12d_mc_vs_atac_scatter` — gene-level mC% vs net-ATAC scatter
- `12e_atac_integration_heatmap` — 2×2 mC × ATAC O/E heatmap
