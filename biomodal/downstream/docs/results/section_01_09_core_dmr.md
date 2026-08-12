# Sections 01-09: Core DMR statistics, QC, volcano/effect-size, GO/KEGG enrichment

## Summary

This group establishes the central differential-methylation result of the paper: in BAP1-KO mutant cerebellum, gene-body 5mC rises while 5hmC falls at the same loci, a coordinated reciprocal signature consistent with a TET-mediated active-demethylation block. Sections 01-02 confirm data quality (8 deep-seq samples, excellent base quality, tight within-condition correlations); Section 03 localizes the signal to gene bodies (not promoters/islands); Sections 04-07 quantify the asymmetric directionality (5mC hypermethylation, 5hmC hypomethylation) at the gene level; Sections 05/09 quantify the coordinated mC↑/hmC↓ pattern; Section 08 shows the affected genes are enriched for RNA splicing / ubiquitin-transferase / autophagy / neurodegeneration pathways, with a developmental-vs-metabolic split between the most- and least-TET-impaired deciles. All quantitative values below are recomputed directly from the run-5 (8-sample, sex-covariate) DMR BED files and the exported TSV tables.

**Canonical run-5 headline numbers (recomputed from `modality/outputs/run-5/.../DMR_*_20260402_191818.bed`, matching `summary_statistics.txt`):**
20,969 genes tested; 5mC significant = 10,775 (51.4%; 7,513 hyper / 3,262 hypo); 5hmC significant = 11,484 (54.8%; 1,963 up / 9,521 down); 8,371 co-significant; 6,589 coordinated mC↑/hmC↓ (78.7%).

> **See the prominent Data-quality flag at the bottom.** `FIGURES.md` and `TODO.md` rows 01-09 still carry OLDER-run numbers (8,836 / 9,930 sig; 6,750 co-significant; 84.6% coordinated; mC r=0.76-0.79). Those are NOT reproducible from the current run-5 tables and should be treated as stale. The current canonical numbers are above.

---

## Section 01: section_01_qc_overview

### Analysis question
Do all sequenced libraries pass QC (depth, base quality, duplication) and is sequencing balanced across conditions?

### Key results
- Mapped reads per sample range = 336M-489M (source: upstream `evoC_Bap1_run_duet-evoC_Summary.csv`; 8 samples)
- Mapped bases per sample range = 44.4B-62.3B (source: upstream summary CSV)
- Duplication rate range = 27.8%-31.6% (source: upstream summary CSV)
- Mean Phred score range = 34.30-34.40 (≈Q34, ~1 error in 2,500 bases) (source: upstream summary CSV)
- Baseline CG 5mC (control autosomes) = 72.22% mean (per-sample 72.09-72.43%) (source: upstream summary CSV)
- Sample count = 8 (4 control: ctrl-F, ctrl-F-B2, ctrl-M, ctrl-M-B2; 4 mutant: mut-F, mut-F-B2, mut-M, mut-M-B2) (source: upstream summary CSV + run-5 BioQC JSON index)

### [INTERPRETATION] Biological meaning
[INTERPRETATION] All eight libraries are high quality (Q34, ~30% duplication, hundreds of millions of reads), so downstream DMR differences reflect biology rather than coverage or quality artifacts. Mutant libraries were sequenced marginally deeper than controls, but the sex-covariate GLM and per-gene q-value thresholds control for this. The ~72% bulk CG methylation is the expected mammalian genome-wide level, confirming the assay is calibrated.

### Plot inventory
- `01_qc_overview/` — 2x2 panel: mapped reads, mapped bases, duplication-rate lollipop, mean-Phred lollipop (3 image formats present).

---

## Section 02: section_02_correlation

### Analysis question
Do replicate samples cluster by methylation profile, and is 5mC more reproducible than 5hmC?

### Key results
- 5mC sample-sample correlation (off-diagonal) range = 0.868-0.903, mean = 0.887 (source: run-5 BioQC JSON `biological_qc_report_8_samples_20260402_185215.json`)
- 5hmC sample-sample correlation (off-diagonal) range = 0.639-0.709, mean = 0.677 (source: run-5 BioQC JSON)
- 5mC within-control mean = 0.886; within-mutant mean = 0.897; between-group mean = 0.884 (source: run-5 BioQC JSON)
- 5hmC within-control mean = 0.685; within-mutant mean = 0.692; between-group mean = 0.668 (source: run-5 BioQC JSON)
- Number of off-diagonal sample pairs = 28 per modality (8 samples) (source: run-5 BioQC JSON)

### [INTERPRETATION] Biological meaning
[INTERPRETATION] 5mC is a stable, high-abundance bulk mark (r≈0.89 across all samples), whereas 5hmC is intrinsically sparser and noisier (r≈0.68), which is the expected biology of a low-abundance oxidation intermediate. The slightly higher within-group than between-group correlations (especially for 5hmC: 0.69 within vs 0.67 between) indicate a detectable but modest genotype effect on the global profile, consistent with the differential signal being concentrated at a subset of gene bodies rather than a genome-wide shift in bulk levels.

### Plot inventory
- `02a_mc_correlation_heatmap/` — 8x8 pheatmap of 5mC sample correlations (Blues).
- `02b_hmc_correlation_heatmap/` — 8x8 pheatmap of 5hmC sample correlations (Greens).
- `02c_correlation_comparison/` — faceted bar chart of within-control / within-mutant / between-group correlations for 5mC vs 5hmC.

---

## Section 03: section_03_dmr_statistics

### Analysis question
Which genomic region class carries the differential-methylation signal, and is the mC vs hmC directionality asymmetric?

### Key results
- Gene bodies significant mC DMRs = 10,775 / 20,969 deduped genes (51.4%) (source: run-5 gene-body mC BED, matches `summary_statistics.txt`)
- CpG shores significant mC DMRs = 9,842 / 32,581 (30.2%) (source: run-5 cpg_shores mC BED, raw/non-deduped as loaded by `region_dmrs`)
- CpG shelves = 6,924 / 29,094 (23.8%); Promoters = 1,692 / 20,417 (8.3%); CpG islands = 442 / 8,910 (5.0%); TSS regions = 192 / 14,165 (1.4%) (source: run-5 region BEDs)
- 5mC direction: 7,513 hypermethylated (69.7%) vs 3,262 hypomethylated (30.3%) (source: run-5 gene-body mC BED / `summary_statistics.txt`)
- 5hmC direction: 1,963 increased (17.1%) vs 9,521 decreased (82.9%) (source: run-5 gene-body hmC BED / `summary_statistics.txt`)
- Non-CG baseline (control autosomes): CHG = 0.628%, CHH = 0.862% modC (vs CG 72.22%) — ~100x lower, 0 significant DMRs (source: upstream summary CSV; section_03 hardcodes 0 non-CG DMRs)

### [INTERPRETATION] Biological meaning
[INTERPRETATION] The differential signal is overwhelmingly a gene-body phenomenon (51% of gene bodies vs 1.4% of TSS regions and 5% of CpG islands are significant), which is the hallmark of dysregulated DNMT3/TET cycling over transcribed gene bodies rather than promoter switching. The directional asymmetry is the key mechanistic clue: gene bodies preferentially GAIN 5mC (70% hyper) while LOSING 5hmC (83% down) — the exact signature of a stalled 5mC→5hmC oxidation step. Non-CG methylation is negligible and unchanged, so the mechanism is CpG-specific.

### Plot inventory
- `03a_dmr_by_region/` — bar chart of significant mC DMRs across 6 region classes (gene bodies dominate).
- `03b_direction_comparison/` — grouped bars: % hyper vs hypo for 5mC and % up vs down for 5hmC, showing the reciprocal asymmetry.
- `03_dmr_region_statistics/` — combined 4-panel: region bars, mC-vs-hmC bars, baseline methylation by context, significant DMRs by context.

---

## Section 04: section_04_volcano_plots

### Analysis question
Visualize per-gene effect size vs significance for 5mC and 5hmC, highlighting key neuronal genes.

### Key results
- 5mC volcano significant genes (subtitle) = 10,775 (q<0.05) (source: run-5 gene-body mC BED via `sum(mc_dmr$significant)`)
- 5hmC volcano significant genes (subtitle) = 11,484 (q<0.05) (source: run-5 gene-body hmC BED)
- Floor q-value (machine-precision cap) ≈ 4.82e-306 for mC top genes; -log10(q) capped at 300 for display (source: `top_mc_dmrs.tsv`)
- 5mC cloud is right-shifted (hypermethylation dominant); 5hmC cloud is left-shifted (hypomethylation dominant) (source: direction counts above, run-5 BEDs)
- Labeled key genes: Syt1, Zbtb20, Trpm3, Cntnap2 (and Dlgap1, Mcu, Arhgap26) (source: section_04 `KEY_GENES`)

### [INTERPRETATION] Biological meaning
[INTERPRETATION] The two volcanoes are mirror images — the same genes that gain 5mC lose 5hmC — which visually encodes the coordinated demethylation-block model. The genes pinned at the machine-precision q-floor (Syt1, etc.) are synaptic/neuronal, indicating the strongest dysregulation falls on the neuronal differentiation program of the developing cerebellum.

### Plot inventory
- `04a_volcano_mc/` — 5mC volcano (mod_difference vs -log10 q), red=hyper / blue=hypo / grey=NS.
- `04b_volcano_hmc/` — 5hmC volcano, dominant blue (hypo) cloud.
- `04_volcano_plots/` — side-by-side combined 5mC|5hmC panel.

---

## Section 05: section_05_coordinated_changes

### Analysis question
Among genes significant in BOTH modifications, what fraction follow the coordinated mC↑/hmC↓ pattern (the demethylation-block signature)?

### Key results
- Genes significant in both mC and hmC = 8,371 (source: `coordinated_changes.tsv`, 8,371 data rows; matches run-5 BED inner-join)
- Coordinated mC↑/hmC↓ genes = 6,589 (78.7%) (source: `coordinated_changes.tsv` col `coordinated_pattern`==TRUE = 6,589; FALSE = 1,782)
- Reverse coordinated (mC↓/hmC↑) = 1,255 (15.0%); same-direction (both up or both down) = 527 (6.3%) (source: run-5 BED join of deduped per-gene tables)
- Top combined-effect coordinated genes: Tmem238 (mc +20.3%/hmc -12.5%), Syt1 (+17.3%/-15.0%), Gclm (+16.4%/-10.9%), Sap30 (+21.2%/-5.4%), Prxl2b (+21.6%/-4.5%) (source: `coordinated_changes.tsv` / `summary_statistics.txt`)
- Syt1: mc_diff = +17.3%, mc_q = 4.82e-306, hmc_diff = -15.0%, hmc_q = 7.08e-306, coordinated = TRUE (source: `coordinated_changes.tsv`); Zbtb20: mc_diff = +7.8%, hmc_diff = -5.6%, coordinated = TRUE (source: `coordinated_changes.tsv`)

### [INTERPRETATION] Biological meaning
[INTERPRETATION] This is the quantitative core of the paper: 78.7% of co-significant genes show the diagnostic mC↑/hmC↓ reciprocity. Because 5hmC is the obligate first oxidation product of 5mC by TET, a simultaneous mC gain + hmC loss at thousands of gene bodies is best explained by TET being unable to access/process 5mC — i.e., a demethylation block downstream of BAP1 loss (BAP1 loss → elevated H2AK119ub → restricted TET access). The fact that this is genome-wide (6,589 genes, not a handful of outliers) argues for a systemic chromatin-level mechanism rather than locus-specific regulation.

### Plot inventory
- `05a_mc_hmc_scatter/` — scatter of mc_diff vs hmc_diff, lower-right (mC↑/hmC↓) quadrant highlighted.
- `05b_top_coordinated_genes/` — top-20 coordinated genes, paired mC/hmC bars.
- `05c_syt1_detail/` — Syt1 control-vs-mutant 5mC/5hmC bar detail.
- `05d_zbtb20_detail/` — Zbtb20 control-vs-mutant detail.
- `05_coordinated_changes/` — 4-gene panel (Syt1, Zbtb20, Trpm3, Cntnap2).

---

## Section 06: section_06_top_genes

### Analysis question
Which individual genes are the most statistically significant DMRs, and how much do the two significant gene sets overlap?

### Key results
- Top-20 5mC DMRs (by q-value): 19 hypermethylated / 1 hypomethylated (source: `top_mc_dmrs.tsv` direction column)
- Top-20 5hmC DMRs (by q-value): 16 hypomethylated / 4 hypermethylated (source: `top_hmc_dmrs.tsv` direction column)
- Top mC DMRs all at q-floor ≈ 4.82e-306 (e.g. 1700025G04Rik mc +6.3%, Add3 +9.7%, Afap1l2 +11.5%, Aig1 +8.1%, Ank3 +4.0%) (source: `top_mc_dmrs.tsv`)
- Venn overlap (both significant) = 8,371 genes; 5mC-only = 10,775 − 8,371 = 2,404; 5hmC-only = 11,484 − 8,371 = 3,113 (source: run-5 BED significant sets / `coordinated_changes.tsv`)
- Total significant 5mC = 10,775; total significant 5hmC = 11,484 (source: run-5 BEDs)

### [INTERPRETATION] Biological meaning
[INTERPRETATION] The most significant mC DMRs are uniformly hypermethylated, while the most significant hmC DMRs are predominantly hypomethylated — the top of each ranked list independently recapitulates the global asymmetry. Large overlap (8,371 of ~11-12k significant genes shared) shows the two marks change together at the same loci, again supporting a single coupled mechanism rather than two independent processes.

> NOTE: `FIGURES.md` Figure 8 states the Venn as "5mC: 8836 | 5hmC: 9930 | Both: 6722" — those are OLD-run values; the current run-5 overlap is 8,371 (see Data-quality flag).

### Plot inventory
- `06a_top_dmrs/` — side-by-side top-20 5mC and top-20 5hmC horizontal bar charts.
- `06b_venn_overlap/` — 2-circle Venn of significant 5mC vs 5hmC gene sets.
- `06_top_genes/` — composite of both bar charts plus the Venn.

---

## Section 07: section_07_effect_size

### Analysis question
What are the genome-wide effect-size distributions of 5mC and 5hmC change, and do they shift in opposite directions?

### Key results
- 5mC net mean change among significant = +1.72% (n=10,775) (source: run-5 gene-body mC BED; matches `summary_statistics.txt`)
- 5mC hyper-only mean = +3.45% (n=7,513, 70% of significant) (source: run-5 gene-body mC BED)
- 5hmC net mean change among significant = -1.66% (n=11,484) (source: run-5 gene-body hmC BED; matches `summary_statistics.txt`)
- 5hmC hypo-only mean = -2.29% (n=9,521, 83% of significant) (source: run-5 gene-body hmC BED)
- Effect sizes are modest per gene (~2-3%) but consistent across thousands of genes; violin plot shows 5mC density above zero, 5hmC density below zero (source: run-5 BEDs)

### [INTERPRETATION] Biological meaning
[INTERPRETATION] The opposing net shifts (+1.7% mC vs -1.7% hmC) are the population-level expression of the coordinated pattern. Individual gene effects are small (a few percent), which is expected for bulk gene-body methylation, but their consistency and reciprocity (mirror-image violins) make the directional signal statistically overwhelming. This rules out a stochastic or noise explanation and supports a directed enzymatic perturbation.

> NOTE: `FIGURES.md`/`TODO.md` quote "+2.27% / -2.08%" with "hyper-only +3.96% (n=6635)" — OLD-run values; current run-5 means are +1.72% / -1.66% (hyper-only +3.45%, n=7,513).

### Plot inventory
- `07_effect_size_distributions/` — 3-part: 5mC effect histogram, 5hmC effect histogram, and 5mC-vs-5hmC violin comparison.

---

## Section 08: section_08_enrichment

### Analysis question
What biological pathways are over-represented among hypermethylated genes, and do the most- vs least-TET-impaired deciles differ functionally?

### Key results
- Input: hypermethylated mC genes (mod_difference>0, significant) = 7,513 genes submitted for enrichment (source: run-5 gene-body mC BED; section_08 `hyper_genes`)
- Top GO Biological Process term = RNA splicing (GO:0008380), 274 genes, GeneRatio 274/6992, FoldEnrichment = 2.42, p.adjust = 6.27e-54 (source: `enrichment_go_bp.tsv`); 6,455 GO BP terms total in table
- Top GO Cellular Component = nuclear speck, 256 genes, p.adjust = 2.55e-49 (source: `enrichment_go_cc.tsv`; 770 terms)
- Top GO Molecular Function = ubiquitin-like protein transferase activity, 269 genes, p.adjust = 1.36e-52 (source: `enrichment_go_mf.tsv`; 1,398 terms)
- Top KEGG pathway = Autophagy - animal (mmu04140), 113 genes, p.adjust = 8.82e-24; followed by Amyotrophic lateral sclerosis (195, 7.86e-21), Endocytosis (152, 7.86e-21), Ubiquitin mediated proteolysis (99, 2.91e-18), several neurodegeneration pathways (source: `enrichment_kegg.tsv`; 348 terms)
- Delta-ratio decile split: D10 (most TET-impaired) top GO BP = pattern specification process (124 genes, p.adjust 6.91e-37), regionalization, cell fate commitment — developmental; D1 (least impaired) top GO BP = oxidative phosphorylation (46, 1.31e-14), aerobic respiration, cellular respiration — metabolic (source: `enrichment_delta_ratio_top_decile_go_bp.tsv`, `enrichment_delta_ratio_bottom_decile_go_bp.tsv`)
- compareCluster GO BP rows: 1,994 enriched in D10 vs 804 in D1 (source: `enrichment_delta_ratio_compare_go_bp.tsv`)

### [INTERPRETATION] Biological meaning
[INTERPRETATION] The hypermethylated gene set is enriched for RNA-processing and ubiquitin-transferase/autophagy/neurodegeneration machinery — notably the ubiquitin axis (MF #1 = ubiquitin-like transferase; KEGG = ubiquitin-mediated proteolysis), which thematically links back to BAP1's role as a deubiquitinase and the H2AK119ub upstream signal. The decile dichotomy is mechanistically informative: the genes most blocked from demethylation (D10) are developmental/patterning/cell-fate genes — i.e., the neuronal differentiation program is the primary casualty — whereas the least-impaired genes (D1) are housekeeping metabolic/respiration genes that are normally constitutively low-methylated and thus less dependent on active TET cycling.

> NOTE: `TODO.md` row 08 cites "RNA splicing #1 GO term (248 genes, q=3.4e-48)"; the current table reports 274 genes at p.adjust 6.27e-54 (still RNA splicing as #1, but different counts — OLD-run figure).

### Plot inventory
- `08a_enrichment_go_bp/` — GO Biological Process dotplot (top 15) for hypermethylated genes.
- `08b_enrichment_go_cc/` — GO Cellular Component dotplot.
- `08c_enrichment_go_mf/` — GO Molecular Function dotplot.
- `08d_enrichment_kegg/` — KEGG pathway dotplot.
- `08e_enrichment_delta_ratio_compare_go_bp/` — compareCluster GO BP, D10 vs D1 decile.
- `08f_enrichment_delta_ratio_compare_kegg/` — compareCluster KEGG, D10 vs D1.
- `08g_enrichment_delta_ratio_top_decile_go_bp/` — standalone GO BP for D10 (most TET-impaired).
- `08h_enrichment_delta_ratio_bottom_decile_go_bp/` — standalone GO BP for D1 (least impaired).

---

## Section 09: section_09_summary

### Analysis question
Aggregate all core statistics into one machine-generated summary report.

### Key results
- Total genes tested = 20,969 (both modalities) (source: `summary_statistics.txt`, regenerated from run-5 `mc_dmr`/`hmc_dmr`)
- 5mC significant = 10,775 (51.4%); 7,513 hyper / 3,262 hypo; mean +1.72% (source: `summary_statistics.txt`)
- 5hmC significant = 11,484 (54.8%); 1,963 up / 9,521 down; mean -1.66% (source: `summary_statistics.txt`)
- Co-significant = 8,371; coordinated mC↑/hmC↓ = 6,589 (78.7%) (source: `summary_statistics.txt`)
- Top coordinated genes table: Tmem238, Syt1, Gclm, Sap30, Prxl2b, Ly6g6e, Gpr68 (source: `summary_statistics.txt`)

### [INTERPRETATION] Biological meaning
[INTERPRETATION] Section 09 is a reporting/aggregation step; its scientific content is identical to Sections 03/05/07. Its numbers ARE the canonical run-5 values (the script populates every field live from the run-5 DMR objects). Its only defect is a hardcoded `Samples: 4 (2 Control, 2 Mutant)` header string in the template — a cosmetic bug, not a data error (see Data-quality flag).

### Plot inventory
- (No image folder — Section 09 emits the text file `tables/summary_statistics.txt` only.) Primary RESULTS.md for this section is placed in `09_summary_statistics` per convention; since no plot folder exists, the section's results live in the group doc and in `summary_statistics.txt`.

---

## Cross-section synthesis

These nine sections build the paper's foundational argument in a single logical chain: the data are clean and well-powered (01-02, 8 deep-seq samples, r≈0.89 for 5mC); the differential signal is specifically a gene-body CpG phenomenon (03); it is directionally asymmetric — 5mC up, 5hmC down (03-04, 07); and at the level of individual genes the two marks move reciprocally at the same loci (05-06), with 78.7% of 8,371 co-significant genes showing the diagnostic mC↑/hmC↓ pattern. Pathway enrichment (08) then shows the affected genes are RNA-processing, ubiquitin-transferase, autophagy and — most tellingly — neuronal-developmental/patterning genes, tying the molecular lesion to the cerebellar neurodevelopmental phenotype and to BAP1's ubiquitin biology. Together this is the quantitative case for the thesis that BAP1 loss imposes a TET-mediated active-demethylation block (via elevated H2AK119ub restricting TET access), restructuring the methylome of the developing cerebellum.

## Tables used

- `summary_statistics.txt` — master text summary of all core DMR counts (run-5 numbers; cosmetic "Samples: 4" header bug).
- `top_mc_dmrs.tsv` — top-20 gene-body 5mC DMRs by q-value (17-column annotated BED-derived).
- `top_hmc_dmrs.tsv` — top-20 gene-body 5hmC DMRs by q-value.
- `coordinated_changes.tsv` — 8,371 co-significant genes with mc/hmc diff, q, ctrl/mut means, coordinated flag, combined effect.
- `cg_top50_mc_dmr_genes.tsv` — top-50 mC DMR genes with paired hmC values (used by section 43/CG-exploratory; corroborates Sec 06).
- `cg_mc_chromosome_distribution.tsv` — per-chromosome mC DMR enrichment (Fisher), shows chrX depletion (obs/exp 0.564, p=2.31e-42).
- `cg_hmc_chromosome_distribution.tsv` — per-chromosome hmC DMR enrichment, chrX depletion (obs/exp 0.712, p=1.98e-21).
- `enrichment_go_bp.tsv` — GO Biological Process over-representation of hypermethylated genes (6,455 terms; #1 RNA splicing).
- `enrichment_go_cc.tsv` — GO Cellular Component (770 terms; #1 nuclear speck).
- `enrichment_go_mf.tsv` — GO Molecular Function (1,398 terms; #1 ubiquitin-like transferase activity).
- `enrichment_kegg.tsv` — KEGG pathways (348 terms; #1 Autophagy - animal).
- `enrichment_delta_ratio_compare_go_bp.tsv` — compareCluster GO BP, D10 vs D1 decile (2,798 rows).
- `enrichment_delta_ratio_top_decile_go_bp.tsv` — GO BP for most-TET-impaired decile (developmental terms).
- `enrichment_delta_ratio_bottom_decile_go_bp.tsv` — GO BP for least-impaired decile (metabolic/respiration terms).
- `enrichment_delta_ratio_compare_kegg.tsv` — compareCluster KEGG, D10 vs D1 (87 rows).
- (Upstream, not in `tables/`) `evoC_Bap1_run_duet-evoC_Summary.csv` — 8-sample sequencing QC; `biological_qc_report_8_samples_20260402_185215.json` — run-5 correlation matrices.

## Data-quality flags

1. **RESOLVED — the 4-vs-8 sample / number discrepancy.** `summary_statistics.txt` header reads `Samples: 4 (2 Control, 2 Mutant)` and is dated 2026-05-12, which superficially contradicts run-5 being 8 samples. INVESTIGATION: `section_09_summary.R` line 50 hardcodes the string `Samples: 4 (2 Control, 2 Mutant)` in its template, but populates ALL numeric fields live from the run-5 `mc_dmr`/`hmc_dmr` objects (lines 111-130). I independently recomputed every headline number directly from the run-5 BED files (`DMR_*_20260402_191818.bed`) using the script's exact significance rule (`dmr_qvalue < 0.05`, dedup-by-gene keeping min q): 20,969 genes, 10,775 mC sig (7,513/3,262), 11,484 hmC sig (1,963/9,521), 8,371 co-significant, 6,589 coordinated (78.7%) — an EXACT match to `summary_statistics.txt`. The run-5 BioQC JSON is literally named `..._8_samples_...` and its index lists all 8 sample IDs (ctrl/mut × F/M × B1/B2). **Conclusion: `summary_statistics.txt` IS the canonical run-5 (8-sample) data; only its `Samples: 4` header string is a stale hardcoded template bug.** Recommend fixing line 50 of `section_09_summary.R` to `Samples: 8 (4 Control, 4 Mutant)`.

2. **`FIGURES.md` (Figures 1-9) carries OLD-run numbers throughout — DO NOT cite from it.** Confirmed stale values that are NOT reproducible from current run-5 tables:
   - Fig 2/4: 8,836 mC sig / 9,930 hmC sig (current: 10,775 / 11,484).
   - Fig 3a region bars: gene bodies 8,836, CpG shores 6,037, shelves 4,045, promoters 641, islands 291, TSS 58 (current run-5: 10,775 / 9,842 / 6,924 / 1,692 / 442 / 192).
   - Fig 6: "85% of 6,750 co-significant" (current: 78.7% of 8,371).
   - Fig 8 Venn: "5mC 8836 | 5hmC 9930 | Both 6722" (current overlap 8,371).
   - Fig 9: mean +2.27% / -2.08%, hyper-only +3.96% n=6635 (current: +1.72% / -1.66%, hyper-only +3.45% n=7,513).
   - Fig 1 prose describes only 4 samples (ctrl-F, ctrl-M, mut-F, mut-M) — the OLD 4-sample run; run-5 has 8.
   - Fig 10 references n=6,635 hyper / 2,201 hypo (current: 7,513 / 3,262).

3. **`TODO.md` rows 01-09 carry OLD-run numbers** (row 02 mC r=0.76-0.79 / hmC r=0.48-0.51 → current 0.87-0.90 / 0.64-0.71; row 04 8,836/9,930; row 05 5,708/6,750/84.6%; row 07 +2.27%/-2.08%; row 08 RNA splicing "248 genes, q=3.4e-48" → current 274 genes, p.adjust 6.27e-54). Update rows 01-09 to run-5 values.

4. **`CLAUDE.md` "92.3% coordinated" is UNVERIFIED in tables.** Both biomodal CLAUDE.md and downstream/CLAUDE.md state "92.3% of co-significant genes show coordinated mC↑/hmC↓." The run-5 `coordinated_changes.tsv` and BED recomputation both give 78.7% (6,589/8,371). [UNVERIFIED: 92.3% per CLAUDE.md, not confirmed in tables — possibly a different denominator (e.g., restricting to large-effect genes) or a different/older run.] Canonical, table-confirmed figure = 78.7%.

5. **Region DMR counts for non-gene-body classes are non-deduplicated.** In `_shared_config.R`, `region_dmrs` loads "Gene bodies" as the deduped `mc_dmr` (one row/gene) but loads CpG islands/shores/shelves/promoters/TSS via raw `load_dmr_bed()` (one row per region feature, multiple per gene). The Section 03a region bar chart therefore mixes a per-gene count (gene bodies) with per-feature counts (other regions); this is a methodological inconsistency to note, though it does not change the qualitative "gene bodies dominate" conclusion.

6. **Section 09 has no plot folder** (text-only output). Its RESULTS.md is written to a placeholder folder per the task convention.

7. **All section 01-08 plot folders present and populated** (3 image formats each — JPG/PDF/SVG; PNG not separately present in these folders). No empty/missing plot folders in the 01-09 range.
