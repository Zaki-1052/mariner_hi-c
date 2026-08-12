## Section 70: Methylation vs MeCP2 by Significance Category (by mark)

**Key numbers:**
- Active (K27ac): 5,354 genes, rho=0.044, p=1.67e-05; Both=4,110, 5mC-only=527, 5hmC-only=717, Neither=575 (70_sig_category_by_mark_stats.tsv)
- Bivalent (K27me3+K27ac): 415 genes, rho=0.191, p=5.91e-07; Both=261 (70_sig_category_by_mark_stats.tsv)
- Fac. Het (K27me3 only): 2,133 genes, rho=0.139, p=7.60e-15; Both=981, 5mC-only=682, 5hmC-only=470 (70_sig_category_by_mark_stats.tsv)
- Neither: 4,579 genes, rho=0.041, p=5.59e-04; Both=2,412, 5hmC-only=1,379 (70_sig_category_by_mark_stats.tsv)

**What this shows:** Coloring the methylation-vs-MeCP2 space by significance class (5mC-only / 5hmC-only / Both / Neither), the combined Spearman correlations are weak (rho 0.04-0.19), confirming MeCP2 fold is not strongly predicted by which modality changed — the Section 68 coupling is carried by chromatin context, not significance class. The coordinated "Both" class (simultaneous 5mC gain + 5hmC loss, the TET-blockade fingerprint) is the dominant population wherever active marks are present (Active 4,110/5,354; Bivalent 261/415).

**Figures:**
- 70a_sig_category_k__ac_only/ — Active context, dots by sig category
- 70b_sig_category_k__me____k__ac/ — Bivalent context
- 70c_sig_category_k__me__only/ — Facultative het context
- 70d_sig_category_neither/ — Unmarked context
