## Section 63: Master Multi-Mark Heatmap at MeCP2 Binding Sites
**Key numbers:**
- Significant MeCP2 peaks = 8,886 (7,686 Up + 1,200 Down); heatmap subsampled to 5,000 (source: 62_mecp2_peak_chromatin_signal.tsv; 63_cluster_assignments.tsv)
- Subsample: 4,306 MeCP2-Up + 694 MeCP2-Down peaks (source: 63_cluster_assignments.tsv)
- 4 ward.D2 clusters: 1,874 / 1,692 / 1,064 / 370 peaks (source: 63_cluster_assignments.tsv)
- Every cluster's modal state is "Unmarked"; only Cluster 3 has appreciable Polycomb (168 peaks) (source: 63_cluster_assignments.tsv)
- Exported matrix = 5,000 peaks x 10 marks (source: 63_peak_signal_matrix.tsv)

**What this shows:** MeCP2 binding sites do not form a single homogeneous Polycomb block; Z-scored multi-mark clustering (k=4) splits them into chromatin-defined clusters, with enhancer-flavored clusters distinct from a more Polycomb-leaning one. The persistence of "Unmarked" as the modal state across clusters (echoing 61i) indicates much MeCP2 occupancy is at regions lacking strong canonical histone marks, where CG mC and K119ub carry the signal.

**Figures:**
- 63a_multimark_heatmap — Z-scored ctrl-signal heatmap (10 marks, k=4)
- 63b_cluster_profiles — per-cluster mean mark bars
- 63c_cluster_composition — chromatin state + MeCP2 direction per cluster
- 63d_log2fc_heatmap — mut/ctrl log2FC heatmap
