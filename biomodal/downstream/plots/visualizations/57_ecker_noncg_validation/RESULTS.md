## Section 57: Ecker WGBS Validation of Non-CG MeCP2 Candidate Peaks
**Key numbers:**
- KEY TEST — Ecker CH (non-CG) at non-CG-candidate vs CG-concordant MeCP2-Up peaks: median 0.0159 vs 0.0129, Wilcoxon p=2.63e-60 (n=2,726 vs 4,960) (57_ecker_noncg_validation_summary.tsv)
- SPECIFICITY CONTROL — Ecker CG goes the other way: median 0.00955 vs 0.0136, p=1.39e-55 (candidates have lower CG) (57_ecker_noncg_validation_summary.tsv)
- MeCP2-Down peaks have the highest Ecker CH (median 0.0219, n=1,200); Not-Significant lowest (0.00523, n=193,764) (57_ecker_noncg_validation_summary.tsv)
- Residual vs Ecker CH (pooled Up+Down): Spearman rho=-0.176, p=1.81e-62, n=8,886 — NEGATIVE (57_ecker_correlation.tsv)
- evoC CHG/CHH at these peaks ~0 (median 0, <1.2% detection) — invisible to evoC (57_ecker_noncg_validation_summary.tsv)

**What this shows:** External wildtype WGBS validates that non-CG-candidate MeCP2-Up peaks sit at higher non-CG (CH) and lower CG methylation regions — a double dissociation matching MeCP2 reading mCA where mCG is absent, and undetectable by evoC. However the residual-vs-CH correlation is negative when pooled and MeCP2-Down peaks have the highest CH, warning that "more CH" does not simply mean "more CG-unexplained MeCP2" (dissected in section 58).

**Figures:**
- `57a_ecker_ch_noncg_vs_concordant/` — Ecker CH, candidate vs concordant (key panel)
- `57b_residual_vs_ecker_ch/` — residual vs CH, pooled (Simpson's-paradox plot)
- `57c_detection_rate/` — CH detection by group (~100%)
- `57d_ecker_cg_control/` — Ecker CG specificity control
- `57e_three_way_ecker_ch/` — CH by MeCP2 Up/Down/NS
- `57f_ecker_vs_evoc_sidebyside/` — Ecker CH vs evoC CHG vs CHH
- `57g_composite/` — CH + residual + CG-control composite
- (`57f_composite/` is a stale duplicate folder; cosmetic)
