## Section 13: ATAC Chromatin-State Enrichment + Loop-Anchor Accessibility
**Key numbers:**
- ATAC-Up peaks = 7,620; ATAC-Down = 3,744 (source: atac_chip_overlap_enrichment.tsv)
- H3K27me3 overlap: ATAC-Up 23.6% vs ATAC-Down 1.0% (OR=0.034, padj<0.001); H3K27ac: ATAC-Up 12.3% vs ATAC-Down 40.2% (OR=4.79) (source: atac_chip_overlap_enrichment.tsv)
- Polycomb = 12.7% of ATAC-Up vs 0.5% of ATAC-Down; Active_Enhancer = 35.8% of ATAC-Down vs 8.6% of ATAC-Up (source: atac_chromatin_state_distribution.tsv)
- 2,910 differential loops (1,723 up, 1,187 down); Active_Enhancer–Active_Enhancer loops 79.6% ATAC-concordant (113/142) vs CTCF_Site–CTCF_Site 15.4% (52/338) (source: loop_anchor_atac_overlap.tsv; loop_atac_concordance_by_type.tsv)

**What this shows:** ATAC-up peaks are unexpectedly enriched at Polycomb/H3K27me3 regions. At Hi-C loop anchors, regulatory (enhancer-enhancer) loops change accessibility concordantly with loop direction (~80%) while structural CTCF loops do not (~15%).

**Figures:**
- `13a_atac_chromatin_state_distribution` — state composition of ATAC-Up vs Down
- `13b_atac_chip_overlap_enrichment` — per-mark overlap + Fisher
- `13c_atac_chromatin_enrichment_heatmap` — O/E heatmap ATAC-direction × state
- `13d_loop_anchor_atac_overlap` — ATAC at loop anchors by loop direction
- `13e_anchor_atac_by_chromatin_state` — anchor ATAC by 8-category state
- `13f_loop_atac_concordance_by_type` — concordance % per loop type
