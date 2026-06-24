## Section 58: Non-CG Methylation Dose-Response on MeCP2 Residual
**Key numbers:**
- DOSE-RESPONSE FAILS — MeCP2-Up OLS residual by Ecker CH quartile is flat-to-decreasing: Q1 1.456, Q2 1.510, Q3 1.502, Q4 (high CH) 1.338 (58_noncg_dose_response_up.tsv)
- Jonckheere increasing p=0.99999 (n.s.); DECREASING p=2.01e-09; Kruskal-Wallis p=1.96e-15 (58_noncg_trend_tests.tsv)
- Negative control MeCP2-Down also decreasing: Q1 -0.985 → Q4 -1.069; Jonckheere decreasing p=1.95e-05 (58_noncg_dose_response_down.tsv, 58_noncg_trend_tests.tsv)
- Joint CG×CH: residual highest at CG-Low/CH-Low (1.602, n=1,970), lowest at CG-High/CH-High (1.316) — tracks LOW CG, not CH level (58_noncg_joint_stratification.tsv)
- OLS slope 0.1909 vs QR(tau=0.5) slope 0.0145; residuals 91% concordant (Spearman 0.906) (58_noncg_regression_comparison.tsv)

**What this shows:** Higher non-CG methylation does NOT produce more CG-unexplained MeCP2 — the residual significantly DECREASES with CH (Jonckheere decreasing p=2e-09), and the same pattern appears in the MeCP2-Down negative control, indicating confounding with genomic context rather than CH-driven recruitment. The joint stratification shows the residual is governed by low CG methylation, not CH level. Result is robust to OLS-vs-quantile-regression. The dose-response (causal mCA-dose) hypothesis is rejected, even though candidates are categorically CH-enriched (section 57).

**Figures:**
- `58a_corrected_residual_vs_ecker_ch/` — 57b corrected to MeCP2-Up only
- `58b_dose_response_up/` — residual by CH quartile, Up (key panel)
- `58c_dose_response_down/` — residual by CH quartile, Down control
- `58d_scatter_by_quartile/` — residual vs CH scatter with LOESS
- `58e_joint_cg_ch_heatmap/` — joint CG×CH median-residual heatmap (key panel)
- `58f_ols_vs_qr_comparison/` — OLS vs quantile-regression concordance
- `58g_composite/` — corrected correlation + dose-response + joint heatmap
