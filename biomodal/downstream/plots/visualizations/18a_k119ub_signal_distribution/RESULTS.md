## Section 18: K119ub Continuous BigWig Signal Analysis
**Key numbers:**
- mC Up median K119ub log2FC = +0.0586 (n=6,993), p(vs 0) = 4.67e-96; Cliff's delta = +0.099 (negligible) (source: k119ub_bigwig_signal_summary.tsv)
- mC Down median log2FC = −0.0797 (n=2,974); Cliff's delta = −0.181 (small) (source: k119ub_bigwig_signal_summary.tsv)
- hmC Down median log2FC = +0.0141 (n=8,757), p(vs background) = 0.299 (NS), Cliff's delta = +0.008 (source: k119ub_bigwig_signal_summary.tsv)
- Background (all DMR genes) median log2FC = +0.0070 (n=12,721), p(vs 0) = 1.81e-20 — global K119ub increase (source: k119ub_bigwig_signal_summary.tsv)

**What this shows:** Continuous signal confirms a genome-wide K119ub increase in BAP1-KO (background median log2FC > 0). Hypermethylated genes shift positive and hypomethylated genes negative, directionally correct, but the gene-specific effect sizes are tiny (Cliff's delta mostly "negligible"); hmC-Down genes are indistinguishable from background. The mechanism is diffuse global K119ub elevation, not sharp locus-specific redistribution.

**Figures:**
- `18a_k119ub_signal_distribution` — ctrl-vs-mut signal per group (paired Wilcoxon)
- `18b_k119ub_log2fc_distributions` — log2FC per group vs background
- `18c_methylation_vs_k119ub_scatter` — mod_difference vs K119ub log2FC (Spearman)
- `18d_k119ub_waterfall` — all genes ranked by K119ub log2FC, DMR genes colored
