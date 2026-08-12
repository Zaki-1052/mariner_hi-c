## Section 05: Coordinated Changes Analysis (Key Insight)

**Key numbers:**
- Genes significant in both mC and hmC = 8,371 (source: `coordinated_changes.tsv`)
- Coordinated mC↑/hmC↓ = 6,589 (78.7%) (source: `coordinated_changes.tsv` coordinated_pattern==TRUE)
- Reverse (mC↓/hmC↑) = 1,255 (15.0%); same-direction = 527 (6.3%) (source: run-5 BED join)
- Top coordinated genes: Tmem238 (+20.3%/-12.5%), Syt1 (+17.3%/-15.0%), Gclm (+16.4%/-10.9%) (source: `coordinated_changes.tsv` / `summary_statistics.txt`)

**What this shows:** This is the paper's quantitative core: 78.7% of co-significant genes show the diagnostic reciprocal mC-gain/hmC-loss pattern. Because 5hmC is the obligate TET oxidation product of 5mC, simultaneous mC↑ and hmC↓ across thousands of gene bodies indicates a genome-wide TET-mediated active-demethylation block downstream of BAP1 loss.

**Figures:**
- `05a_mc_hmc_scatter/` — mc_diff vs hmc_diff, mC↑/hmC↓ quadrant highlighted.
- `05b_top_coordinated_genes/` — top-20 coordinated genes (paired mC/hmC bars).
- `05c_syt1_detail/`, `05d_zbtb20_detail/` — per-gene control-vs-mutant detail.
- `05_coordinated_changes/` — 4-gene panel (Syt1, Zbtb20, Trpm3, Cntnap2).

> Note: `FIGURES.md`/`CLAUDE.md` cite 6,750 / 84.6% / 92.3%; table-confirmed run-5 value is 6,589 / 78.7%.
