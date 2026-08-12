## Section 01: QC & Data Overview

**Key numbers:**
- Mapped reads = 336M-489M per sample; mapped bases = 44.4B-62.3B (source: upstream `evoC_Bap1_run_duet-evoC_Summary.csv`)
- Duplication = 27.8%-31.6%; mean Phred = 34.30-34.40 (≈Q34) (source: upstream summary CSV)
- Baseline control CG 5mC = 72.22% mean (source: upstream summary CSV)
- 8 samples (4 control + 4 mutant, both batches, both sexes) (source: run-5 BioQC JSON index)

**What this shows:** All eight run-5 deep-seq libraries pass QC with excellent base quality, comparable duplication, and adequate depth; mutant libraries are marginally deeper than controls. Bulk CpG methylation sits at the expected ~72% mammalian level, confirming assay calibration.

**Figures:**
- `01_qc_overview/` — 2x2 panel: mapped reads, mapped bases, duplication-rate lollipop, mean-Phred lollipop.
