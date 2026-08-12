## Section 16: Raw Concordance Bar Charts
**Key numbers:**
- MeCP2: mC Up = 30.0% MeCP2-Up (2,058/6,862), mC Down = 68.2% MeCP2-Down; hmC Down = 29.2% MeCP2-Up (source: raw_concordance_all_marks.tsv)
- ATAC: mC Up = 54.9% ATAC-Down (938/1,710), mC Down = 89.7% ATAC-Up; hmC Down = 47.9% ATAC-Down (source: raw_concordance_all_marks.tsv)
- K119ub: mC Up = 68.3% Gained (1,358/1,988), mC Down = 73.7% Lost; hmC Down = 62.3% Gained (source: raw_concordance_all_marks.tsv)
- Dominant-group (mC Up vs hmC Down): MeCP2 30.0/29.2, ATAC 54.9/47.9, K119ub 68.3/62.3 (source: raw_concordance_summary.tsv)

**What this shows:** Raw concordance rates (% of genes in each methylation group showing the predicted mark direction) make effect sizes visible. Coupling is real for K119ub (~68%) and ATAC (~55%) in the dominant mC-Up group but barely above chance for MeCP2 (~30%). Motivates the Section 17 honest re-assessment.

**Figures:**
- `16a_mecp2_concordance_bars` — MeCP2 concordance across 4 groups
- `16b_atac_concordance_bars` — ATAC concordance across 4 groups
- `16c_k119ub_concordance_bars` — K119ub concordance across 4 groups
- `16d_raw_concordance_comparison` — mC-Up vs hmC-Down summary, all 3 marks
