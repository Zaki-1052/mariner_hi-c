## Section 53: MeCP2 Non-CG Methylation at Peak Resolution
**Key numbers:**
- CHG detection 45/8,886 peaks (0.51%); CHH 146/8,886 (1.64%) have any mCH in ≥1 of 8 samples (mecp2_peak_noncg_summary.tsv)
- MeCP2 vs shuffled control: CHG 0.51% vs 0.29% (Fisher OR=1.734, p=0.0316); CHH 1.64% vs 1.04% (OR=1.597, p=5.17e-04) (mecp2_peak_noncg_summary.tsv)
- Differential DMRs: CHG 0/8,356 at q<0.05 (best q=0.247); CHH 2/8,886 at q<0.05 (best q=0.00317) (mecp2_peak_noncg_summary.tsv)
- Peak Wilcoxon ctrl-vs-mut: CHG p=0.0286, CHH p=0.1143 (mecp2_peak_noncg_summary.tsv)

**What this shows:** Non-CG methylation is modestly but significantly enriched at MeCP2 peaks over shuffled controls (Fisher OR ~1.6-1.7), yet absolute detection is ≤1.6% of peaks and there is essentially no differential signal between mutant and control (0 CHG DMRs, 2 CHH DMRs in a chr8 cluster). The coverage/context-density panel confirms the negative is not a power artifact, and gene-body vs peak resolution give the same flat answer.

**Figures:**
- `53a_detection_rate/` — detection categories per context × MeCP2 direction
- `53b_nonzero_peak_scatter/` — ctrl vs mut at non-zero peaks
- `53c_mecp2_vs_control_baseline/` — MeCP2 vs shuffled control, Fisher tests
- `53d_dmr_rank_plot/` — DMR q-value ranks (sex-covariate GLM)
- `53e_chr8_cluster_spotlight/` — chr8 35-50 Mb CHH cluster
- `53f_resolution_comparison/` — gene-body vs peak non-CG
- `53g_coverage_quality/` — coverage and context density (power check)
