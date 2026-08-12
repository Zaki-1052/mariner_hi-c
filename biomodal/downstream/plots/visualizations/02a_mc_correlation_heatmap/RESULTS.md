## Section 02: Sample Correlation Analysis

**Key numbers:**
- 5mC sample-sample correlation = 0.868-0.903 (mean 0.887) (source: run-5 BioQC JSON `biological_qc_report_8_samples_...json`)
- 5hmC sample-sample correlation = 0.639-0.709 (mean 0.677) (source: run-5 BioQC JSON)
- 5mC within-mutant 0.897 / within-control 0.886 / between 0.884 (source: run-5 BioQC JSON)
- 5hmC within-mutant 0.692 / within-control 0.685 / between 0.668 (source: run-5 BioQC JSON)

**What this shows:** 5mC is highly reproducible across all 8 samples (r≈0.89), while the sparse 5hmC mark is noisier (r≈0.68) as expected for a low-abundance oxidation intermediate. Within-group correlations slightly exceed between-group, indicating a modest, locus-concentrated genotype effect rather than a global bulk shift.

**Figures:**
- `02a_mc_correlation_heatmap/` — 8x8 pheatmap of 5mC correlations (Blues).
- `02b_hmc_correlation_heatmap/` — 8x8 pheatmap of 5hmC correlations (Greens).
- `02c_correlation_comparison/` — within/between-group bar comparison, 5mC vs 5hmC.
