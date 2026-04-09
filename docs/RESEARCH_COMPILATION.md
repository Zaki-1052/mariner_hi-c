# Research Compilation: Regulation of Chromatin Conformation by BAP1 in the Brain

**Author:** Zakir Alibhai
**Affiliation:** Ferguson Lab, UCSD
**Organism:** Mouse (mm10 genome)
**Tissue:** Cerebellum
**Conditions:** BAP1-KO (Mutant) vs Wildtype (Control)
**Timepoints:** Early (P13) and Late (P60/Adult)
**Replicates:** n=3 biological replicates per condition per timepoint
**Samples:** ctrl_M1, ctrl_M2, ctrl_M3, mut_M1, mut_M2, mut_M3

---

## 1. Project Overview & Biological Context

### The Biological Question

BAP1 is a histone deubiquitinase that removes H2AK119ub (monoubiquitination of histone H2A at lysine 119), a mark deposited by Polycomb Repressive Complex 1 (PRC1). This project investigates how BAP1 loss reshapes three-dimensional chromatin architecture in the mouse cerebellum across two developmental timepoints.

### Epigenomic Data Available

| Mark | Type | Assay |
|------|------|-------|
| H3K27ac | Active enhancer/promoter | CUT&RUN |
| H3K27me3 | Polycomb repression | CUT&RUN |
| H3K4me1 | Enhancer (poised/active) | CUT&RUN |
| H3K4me3 | Active promoter | CUT&RUN |
| H2AK119ub | PRC1 ubiquitination (BAP1 substrate) | CUT&RUN |
| ATAC-seq | Chromatin accessibility | ATAC-seq |
| RNA-seq | Gene expression | Bulk RNA-seq (both timepoints) |
| Methylation | 5mC/5hmC | Biomodal DUET evoC |

---

## 2. Hi-C Data Processing

Hi-C data was processed on SDSC Expanse (HPC, SLURM). Two complementary pipelines were used:
- **Nextflow-HiC pipeline** for deep sequencing production (cool/mcool files, ValidPairs)
- **Juicer pipeline** for .hic file generation, HiCCUPS loop calling, and APA

Post-processing included juicer_pre for Juicer v2.0 .hic file generation, HiCCUPS loop calling on individual and merged samples, and HOMER tag directory generation for differential TAD/compartment analysis.

---

## 3. Differential Chromatin Loop Analysis (Mariner/edgeR)

### Pipeline

11-step replicate-aware differential loop analysis using the mariner R/Bioconductor package and edgeR quasi-likelihood GLM framework, run at three resolutions (5kb, 10kb, 25kb).

**Statistical framework:**
- Normalization: TMM (Trimmed Mean of M-values)
- Design: `~group` (intercept + mutant effect)
- Testing: `glmQLFTest` for mutant vs control
- Dispersion: Data-driven, `robust=TRUE`
- Multiple testing: Benjamini-Hochberg FDR

### Per-Resolution Results

**Late (P60/Adult):**

| Resolution | Total Tested | Significant (FDR<0.05) | Down (lost) | Up (gained) |
|------------|-------------|------------------------|-------------|-------------|
| 5 kb | 17,982 | 1,766 (9.8%) | 639 | 1,127 |
| 10 kb | 22,632 | 3,981 (17.6%) | 1,741 | 2,240 |
| 25 kb | 20,398 | 4,774 (23.4%) | 2,262 | 2,512 |

*Source: data/tsvs/figure_2_loop_rewiring/2A_late_edgeR_{5,10,25}kb.tsv*

**Early (P13):**

| Resolution | Total Tested | Significant (FDR<0.05) | Down (lost) | Up (gained) |
|------------|-------------|------------------------|-------------|-------------|
| 5 kb | 16,350 | 6 (0.04%) | 1 | 5 |
| 10 kb | 22,122 | 87 (0.4%) | 34 | 53 |
| 25 kb | 20,991 | 295 (1.4%) | 149 | 146 |

*Source: data/tsvs/figure_2_loop_rewiring/2A_early_edgeR_{5,10,25}kb.tsv*

### Merged Results (Cross-Resolution, Non-Redundant)

| Metric | Early (P13) | Late (P60) |
|--------|-------------|------------|
| Total differential loops | 165 | 2,910 |
| Down in BAP1-KO (lost) | 92 (55.8%) | 1,187 (40.8%) |
| Up in BAP1-KO (gained) | 73 (44.2%) | 1,723 (59.2%) |

*Source: data/upstream/loop_calls/{early,late}_characterized_loops.tsv*

---

## 4. Loop Annotation & ChIP-seq Integration

### 8-Category Chromatin State Classification

Priority-based classification using 5 histone marks + CTCF DNA motifs (114,081 genome-wide):

| Priority | Category | Definition |
|----------|----------|-----------|
| 1 | Active_Promoter | H3K4me3+ AND NOT K27me3 AND <=2kb TSS |
| 2 | Repressed_Promoter | H3K27me3+ AND NOT K27ac AND <=2kb TSS |
| 3 | Bivalent_Promoter | K4me3+K27me3 overlap |
| 4 | Polycomb | H3K27me3+ AND >2kb TSS |
| 5 | Active_Enhancer | H3K27ac+ AND >2kb TSS |
| 6 | Poised_Enhancer | H3K4me1+ AND NOT K27ac AND NOT K27me3 AND >2kb |
| 7 | CTCF_Site | CTCF DNA motif+ |
| 8 | Other | No marks |

### Anchor Type Distribution (Late, 2,910 Loops)

| Category | Anchor 1 (%) | Anchor 2 (%) |
|----------|-------------|-------------|
| CTCF_Site | 26.6 | 26.7 |
| Poised_Enhancer | 21.2 | 19.1 |
| Polycomb | 16.0 | 15.8 |
| Active_Enhancer | 10.8 | 12.0 |
| Repressed_Promoter | 9.2 | 10.0 |
| Active_Promoter | 8.7 | 8.8 |
| Bivalent_Promoter | 0.9 | 0.8 |
| Other | 6.6 | 6.7 |

*Source: data/tsvs/figure_2_loop_rewiring/2F_late_anchor_type_summary.tsv*

### Anchor Type Distribution (Early, 165 Loops)

| Category | Anchor 1 (%) | Anchor 2 (%) |
|----------|-------------|-------------|
| Repressed_Promoter | 36.4 | 35.2 |
| CTCF_Site | 28.5 | 27.9 |
| Polycomb | 10.9 | 12.7 |
| Poised_Enhancer | 8.5 | 6.1 |
| Active_Enhancer | 4.2 | 3.0 |
| Active_Promoter | 1.8 | 6.7 |
| Bivalent_Promoter | 1.2 | 0.6 |
| Other | 8.5 | 7.9 |

*Source: data/tsvs/figure_2_loop_rewiring/2F_early_anchor_type_summary.tsv*

### Top Loop Type Combinations (Late)

| Loop Type | Count | % |
|-----------|-------|---|
| CTCF_Site-CTCF_Site | 338 | 11.6 |
| Poised_Enhancer-CTCF_Site | 306 | 10.5 |
| Polycomb-CTCF_Site | 190 | 6.5 |
| Polycomb-Polycomb | 175 | 6.0 |
| Repressed_Promoter-Polycomb | 172 | 5.9 |
| Active_Enhancer-Poised_Enhancer | 162 | 5.6 |
| Poised_Enhancer-Poised_Enhancer | 148 | 5.1 |
| Active_Enhancer-Active_Enhancer | 142 | 4.9 |

*Source: data/tsvs/figure_2_loop_rewiring/2G_late_loop_type_summary.tsv*

---

## 5. Loop Distance Shift ("Loop Rewriting")

BAP1-KO preferentially loses long-range loops and gains shorter-range loops.

### Distance Statistics (Late)

| Metric | Lost Loops | Gained Loops | Ratio |
|--------|-----------|-------------|-------|
| Count | 1,187 | 1,723 | 0.69 |
| Median distance | 625 kb | 320 kb | 1.95x |
| <100 kb | 96 (8.1%) | 206 (12.0%) | 0.47 |
| 100-500 kb | 419 (35.3%) | 1,012 (58.7%) | 0.41 |
| 500 kb - 1 Mb | 259 (21.8%) | 324 (18.8%) | 0.80 |
| >1 Mb | 413 (34.8%) | 181 (10.5%) | 2.28 |

>1 Mb loops: 3.3x fold enrichment for loss.

*Source: output/loops_visualization_extended/late/distance_shift_statistics.txt; data/tsvs/figure_2_loop_rewiring/2B_late_distance_shift_summary.tsv*

### Statistical Tests (Late)

| Test | Statistic | P-value |
|------|-----------|---------|
| Kolmogorov-Smirnov | D = 0.279 | < 2.2e-16 |
| Wilcoxon Rank-Sum | W = 1,356,071 | < 2.2e-16 |
| Spearman (distance vs logFC) | rho = -0.244 | < 2.2e-16 |
| Chi-Square (distance categories) | chi2 = 294.94, df = 3 | < 2.2e-16 |

*Source: output/loops_visualization_extended/late/distance_shift_statistics.txt*

---

## 6. Shared Anchor Analysis (Loop Switching)

Specific genomic loci participate in BOTH a lost loop and a gained loop at the same anchor.

### Results (Late)

| Metric | Value |
|--------|-------|
| Shared anchor regions | 212 |
| Total loops at shared anchors | 604 |
| Direction support (lost distance > gained distance) | 176/212 = 83.0% |
| Median lost loop distance | 1,150 kb |
| Median gained loop distance | 340 kb |
| Distance ratio (lost/gained) | 3.4x |

*Source: data/tsvs/supplemental/shared_anchors.tsv (212 rows); data/tsvs/supplemental/shared_anchor_paired_distance.tsv (212 rows, direction computed from median_lost_distance vs median_gained_distance columns); data/tsvs/supplemental/shared_anchor_loops.tsv (604 rows)*

### ChIP-seq Enrichment at Shared Anchors (vs Non-Shared)

| Mark | Odds Ratio | P-value | Direction |
|------|-----------|---------|-----------|
| H3K27me3 | 2.04 | 1.75e-24 | ENRICHED |
| Bivalent_Promoter | 1.74 | 0.013 | Enriched |
| H3K4me1 | 1.27 | 4.6e-04 | Enriched |
| H3K4me3 | 1.05 | 0.64 | NS |
| H3K27ac | 0.68 | 3.9e-06 | DEPLETED |

*Source: data/tsvs/supplemental/shared_anchor_chip_enrichment.tsv*

### Polycomb-Specific Shared Loops

312 Polycomb-classified shared loops analyzed separately.

*Source: data/tsvs/supplemental/polycomb_shared_loops.tsv (312 rows)*

---

## 7. H2AK119ub Integration (BAP1 Substrate)

### DiffBind Peak Data

| Category | Count |
|----------|-------|
| Total K119ub peaks tested | 41,392 |
| K119ub UP in mutant (FDR<0.05, Fold>0) | 16,715 |
| K119ub DOWN in mutant (FDR<0.05, Fold<0) | 5,097 |

*Source: data/upstream/diffbind/K119ub_diffbind_results_summit_appended_ap.txt (counted with awk)*

### Continuous Signal Correlation (K119ub Change vs Loop Strength Change)

| Anchor Group | Spearman rho | P-value | N |
|-------------|-------------|---------|---|
| Active | -0.314 | 7.46e-19 | 764 |
| Other | -0.013 | 0.676 | 973 |
| Polycomb/Repressive | +0.177 | 2.82e-09 | 1,118 |

*Source: data/tsvs/figure_3_epigenetic_integration/3A_k119ub_correlation_results.tsv*

### Logistic Regression: P(lost) ~ K119ub_fc + log(distance)

| Term | Estimate | P-value | Odds Ratio |
|------|----------|---------|------------|
| mean_anchor_k119ub_fc | 2.370 | 1.71e-91 | 10.70 |
| log_distance | 1.680 | 1.67e-67 | 5.36 |

*Source: data/tsvs/figure_3_epigenetic_integration/3A_k119ub_logistic_regression.tsv*

### K119ub Enrichment by Chromatin State (Selected Significant Results)

| Chromatin State | Overlap | OR | FDR | Significant |
|-----------------|---------|-----|-----|-------------|
| Active_Enhancer | K119ub_up | 4.72 | 1.81e-09 | TRUE |
| Active_Promoter | K119ub_up | 2.94 | 2.40e-06 | TRUE |
| Polycomb | K119ub_down | 2.92 | 2.25e-04 | TRUE |
| Poised_Enhancer | K119ub_up | 2.05 | 2.16e-06 | TRUE |
| Repressed_Promoter | K119ub_down | 3.39 | 2.16e-06 | TRUE |
| Other | K119ub_up | 2.98 | 2.57e-09 | TRUE |

*Source: data/tsvs/figure_3_epigenetic_integration/3A_k119ub_enrichment_by_chromstate.tsv (14 tests, 11 significant)*

---

## 8. Differential ChIP-seq Polycomb Enrichment

### All Loops: Significant Enrichments (FDR < 0.05)

| Mark | Test vs Reference | OR | FDR | Direction |
|------|-------------------|-----|-----|-----------|
| K27me3_down | short_range_gained vs unchanged | 8.78 | 6.86e-102 | enriched |
| H2AK119ub_up | long_range_lost vs short_range_gained | 4.87 | 1.86e-90 | enriched |
| K27me3_up | long_range_lost vs short_range_gained | 5.05 | 2.85e-28 | enriched |
| K27me3_up | long_range_lost vs unchanged | 3.18 | 3.58e-28 | enriched |
| H2AK119ub_down | short_range_gained vs unchanged | 3.20 | 4.46e-41 | enriched |
| H2AK119ub_up | long_range_lost vs unchanged | 2.18 | 4.56e-41 | enriched |

*Source: data/tsvs/figure_3_epigenetic_integration/3C_enrichment_tests_all_loops.tsv (12 tests)*

### Polycomb-Specific Loops

At Polycomb loops specifically, the only significant test: H2AK119ub_down at long_range_lost vs short_range_gained (OR = 3.92, FDR = 0.0024). All other tests NS within Polycomb subset due to small N (121 lost, 111 gained).

*Source: data/tsvs/figure_3_epigenetic_integration/3C_enrichment_tests_polycomb.tsv*

---

## 9. DiffBind Summary (All Marks)

| Mark | Total Peaks | FDR<0.05 Up | FDR<0.05 Down |
|------|------------|-------------|---------------|
| H2AK119ub | 41,392 | 16,715 | 5,097 |
| H3K27ac | 25,669 | 5,077 | 6,706 |
| H3K27me3 | 18,324 | 2,293 | 4,811 |
| ATAC-seq | 75,867 | 9,263 | 4,159 |

*Source: data/upstream/diffbind/ files (counted with awk: FDR column < 0.05 AND Fold column > or < 0)*

---

## 10. DEG-Loop Anchor Integration

457 differentially expressed genes mapped to differential loop anchors using regulatory domain assignment.

*Source: data/tsvs/figure_5_model_functional/5B_deg_anchor_genes.tsv (457 data rows)*

---

## 11. TAD Boundary Analysis (TADCompare)

TADCompare pipeline at 10kb resolution for differential boundary detection.

### Boundary Counts

| Metric | Early (P13) | Late (P60) |
|--------|-------------|------------|
| Total boundaries | 23,104 | 21,755 |
| Differential boundaries | 4,349 (18.8%) | 4,144 (19.0%) |
| Control-enriched | 2,114 | 1,927 |
| Mutant-enriched | 2,235 | 2,217 |
| Merge events | 724 | 781 |
| Split events | 808 | 593 |
| Strength Change | 975 | 1,250 |
| Complex | 936 | 771 |
| Shifted | 906 | 746 |

*Source: data/tsvs/figure_3_epigenetic_integration/3D_timepoint_comparison_stats.tsv*

### Boundary-Loop Cross-Reference (Late)

| Metric | Value |
|--------|-------|
| Permutation fold enrichment (10kb threshold) | 2.08x |
| Permutation p-value | < 0.001 |
| Fisher OR (lost vs gained near boundary, 10kb) | 1.46 |
| Fisher p-value | 4.79e-06 |
| Population-level directional concordance | 69.6% |

*Source: data/tsvs/figure_1_tads_boundaries_compartments/1F_permutation_test_results.tsv; 1F_enrichment_statistics.tsv; data/tsvs/supplemental/shared_boundary_concordance_summary.tsv*

---

## 12. Chromatin Compartment (PC1) Analysis

HOMER `getDiffExpression.pl` on PC1 eigenvectors at 25kb resolution.

### Early (P13) — 101,684 regions analyzed

| Threshold | Significant Regions | B->A (More Active) | A->B (More Inactive) |
|-----------|--------------------|--------------------|---------------------|
| Standard (FDR<0.05, \|Diff\|>0.30) | 8,154 (7.46% genome) | 5,282 (132.1 Mb) | 2,872 (71.8 Mb) |
| Relaxed (FDR<0.15, \|Diff\|>0.15) | 26,733 (24.5% genome) | 14,411 (360.3 Mb) | 12,322 (308.1 Mb) |

*Source: tads/tad-pc-analysis/output/compartment_analysis_early/ (generated 2026-04-09 from inputs/old/diffcompartments.txt)*

### Early 7-Category Breakdown (Standard)

| Category | bp (Mb) | % Genome | Regions |
|----------|---------|----------|---------|
| Flipped A->B | 42.4 | 1.55% | 1,694 |
| Flipped B->A | 60.6 | 2.22% | 2,426 |
| Strengthened A | 56.5 | 2.07% | 2,258 |
| Strengthened B | 12.2 | 0.45% | 489 |
| Weakened A | 17.2 | 0.63% | 689 |
| Weakened B | 14.9 | 0.55% | 598 |

### Late (P60/Adult) — 104,071 regions analyzed

| Threshold | Significant Regions | B->A (More Active) | A->B (More Inactive) |
|-----------|--------------------|--------------------|---------------------|
| Standard (FDR<0.05, \|Diff\|>0.30) | 8,189 (7.50% genome) | 5,485 (137.1 Mb) | 2,704 (67.6 Mb) |
| Relaxed (FDR<0.15, \|Diff\|>0.15) | 24,189 (22.1% genome) | 12,926 (323.2 Mb) | 11,263 (281.6 Mb) |

*Source: data/tsvs/figure_1_tads_boundaries_compartments/1D_compartment_genome_pct_summary.txt*

### Late 7-Category Breakdown (Standard)

| Category | bp (Mb) | % Genome | Regions |
|----------|---------|----------|---------|
| Flipped A->B | 16.4 | 0.60% | 656 |
| Flipped B->A | 37.9 | 1.39% | 1,517 |
| Strengthened A | 53.8 | 1.97% | 2,151 |
| Strengthened B | 27.3 | 1.00% | 1,093 |
| Weakened A | 23.9 | 0.87% | 955 |
| Weakened B | 45.4 | 1.66% | 1,817 |

Both timepoints show ~7.5% of genome with significant compartment shifts at standard thresholds, with B->A (toward active) outnumbering A->B roughly 2:1. Early has more complete flips; late has more strengthening/weakening.

---

## 13. Chromatin Stripe Analysis

Differential stripe analysis using Quagga with mariner/edgeR quantification.

### Results

| Metric | Early (P13) | Late (P60) |
|--------|-------------|------------|
| Total stripes | 286 | 200 |
| Lost (control_only) | 126 (44%) | 83 (42%) |
| Gained (mutant_only) | 86 (30%) | 73 (36%) |
| Shared | 74 (26%) | 44 (22%) |
| FDR < 0.05 | 0 | 0 |
| FDR < 0.10 | 0 | 1 |

No significant widespread stripe changes. BCV ~6-7%, indicating good power.

*Source: stripes/outputs/ (early and late summary files)*

---

## 14. ABC Enhancer-Gene Linkage Model

Activity-By-Contact model: `ABC = (Enhancer Activity x Contact Frequency) / normalization`. Both WT and KO use an identical consensus enhancer universe (ATAC-seq peaks) to ensure Delta-ABC measures real linkage changes.

### Scale

| Metric | Value |
|--------|-------|
| Total E-G pairs analyzed | 179,709 |
| E-G pairs with RNA-seq | 113,943 |

*Source: data/tsvs/figure_4_abc_analysis/4A_delta_abc_all_pairs.tsv (179,709 rows); 4A_delta_abc_with_rnaseq.tsv (113,943 rows)*

### Enhancer Subset Stratified Analysis

| Class | N Enhancers | Median Loop logFC | RNA Concordance |
|-------|------------|-------------------|----------------|
| Activity_Lost | 7,503 | -0.088 | 62.0% |
| K119ub_Only | 2,479 | -0.054 | 48.7% |
| Activity_Gain | 2,851 | +0.066 | 67.2% |
| Stable | 42,864 | -0.013 | 48.4% |

K119ub_Only shows real but sub-functional contact weakening (concordance at chance level).

*Source: data/tsvs/figure_4_abc_analysis/4D_class_level_summary.tsv; 4B_concordance_by_class.tsv*

### K119ub-ABC Correlation

| Comparison | Spearman rho | P-value | N |
|-----------|-------------|---------|---|
| delta_activity vs K119ub_Fold | -0.208 | 0 | 35,777 |
| mean_delta_unnorm vs K119ub_Fold | -0.248 | 0 | 35,777 |
| mean_delta_ABC vs K119ub_Fold | -0.155 | 2.91e-190 | 35,777 |

*Source: data/tsvs/figure_4_abc_analysis/4F_k119ub_abc_correlation_summary.tsv*

### Activity vs Contact Decomposition (by Enhancer Class)

| Class | Median Activity (raw delta) | Median Contact (raw delta) | Spearman rho (raw) |
|-------|---------------------------|---------------------------|-------------------|
| Activity_Lost | -4.354 | -0.00096 | 0.127 |
| K119ub_Only | -0.958 | 0.000 | 0.273 |
| Activity_Gain | +2.385 | +0.00173 | 0.082 |
| Stable | +0.076 | +0.00016 | 0.205 |
| Unclassified | +0.622 | +0.000075 | 0.192 |

*Source: data/tsvs/figure_4_abc_analysis/4A_activity_contact_summary.tsv*

### Paired-Anchor Analysis Results

The following numbers are from abc/docs/ analysis documentation generated during the pipeline run. They were not independently re-computed from raw TSVs for this compilation:

- Unnormalized Delta(AxC) sum vs RNA-seq: Spearman rho = 0.582
- Paired-anchor dABC concordance (with geometric constraint): 89.7% (p = 1.67e-48)
- 3-way concordance (loop + ABC + RNA-seq): 88.2% (p = 1.69e-45)
- K119ub at paired enhancers: rho = -0.401 (p = 5.48e-13, n = 490)
- Genes with both significant dABC and DE: 940 (553 concordant = 58.8%)

*Source: abc/docs/step9b-paired-anchor-analysis.md; abc/docs/steps5-7-9-abc-pipeline-results.md*

---

## 15. HOMER Motif Analysis at Enhancer Subsets

| Comparison | Sig Motifs (q<0.05) | Top Motif | Top P-value |
|-----------|-------------------|-----------|------------|
| Activity_Lost vs Stable | 369/450 | Olig2 | 1e-120 |
| K119ub_Only vs Stable | 357/450 | Atoh1 | 1e-41 |
| Activity_Gain vs Stable | 30/450 | TCF4 | 1e-9 |
| T3_high vs T1_low (K119ub dose) | 0/450 | Rfx2 | 0.1 (NS) |

*Source: data/tsvs/supplemental/homer_motif_summary_stats.tsv*

---

## 16. Gene Ontology Enrichment

### Long-Range Lost Loops (Top 5 GO BP Terms)

| Term | Gene Ratio | Fold Enrichment | Adj. P-value |
|------|-----------|----------------|-------------|
| Locomotory behavior | 29/684 | 4.19 | 4.69e-07 |
| Camera-type eye development | 32/684 | 3.17 | 2.67e-05 |
| Embryonic organ morphogenesis | 28/684 | 3.39 | 3.60e-05 |
| Sex differentiation | 29/684 | 3.26 | 3.63e-05 |
| Gonad development | 25/684 | 3.52 | 4.68e-05 |

*Source: data/tsvs/figure_5_model_functional/5A_go_long_lost_loops.tsv*

### Short-Range Gained Loops (Top 5 GO BP Terms)

| Term | Gene Ratio | Fold Enrichment | Adj. P-value |
|------|-----------|----------------|-------------|
| Pattern specification process | 77/1311 | 3.39 | 2.62e-17 |
| Regionalization | 71/1311 | 3.45 | 1.49e-16 |
| Sensory organ morphogenesis | 50/1311 | 3.21 | 6.41e-10 |
| Cell fate commitment | 44/1311 | 3.12 | 3.02e-08 |
| Embryonic organ morphogenesis | 47/1311 | 2.97 | 3.02e-08 |

*Source: data/tsvs/figure_5_model_functional/5A_go_short_gained_loops.tsv*

---

## 17. K119ub Enrichment at Loop Categories (Global)

| K119ub Category | Test | Reference | OR | FDR | Direction |
|----------------|------|-----------|-----|-----|-----------|
| K119ub_up | long_range_lost | short_range_gained | 4.87 | 3.73e-90 | enriched |
| K119ub_up | long_range_lost | unchanged | 2.18 | 4.56e-41 | enriched |
| K119ub_up | short_range_gained | unchanged | 0.45 | 2.60e-45 | depleted |
| K119ub_down | short_range_gained | unchanged | 3.20 | 5.57e-41 | enriched |
| K119ub_down | long_range_lost | unchanged | 2.08 | 6.09e-08 | enriched |
| K119ub_down | long_range_lost | short_range_gained | 0.65 | 0.0012 | depleted |

*Source: data/tsvs/figure_3_epigenetic_integration/3A_k119ub_enrichment_global.tsv (12 tests)*

---

## 18. Complete Results Summary

| Analysis | Key Metric | Value | Source |
|----------|-----------|-------|--------|
| **Differential Loops (Late)** | Non-redundant | 2,910 | characterized_loops.tsv |
| | Up in BAP1-KO | 1,723 (59.2%) | " |
| | Down in BAP1-KO | 1,187 (40.8%) | " |
| **Differential Loops (Early)** | Non-redundant | 165 | " |
| | Up | 73 (44.2%) | " |
| | Down | 92 (55.8%) | " |
| **Loop Rewriting** | Median distance lost/gained | 625 kb / 320 kb (1.95x) | distance_shift_statistics.txt |
| | >1Mb lost enrichment | 3.31x | 2B_late_distance_shift_summary.tsv |
| **Shared Anchors** | Hub regions | 212 | shared_anchors.tsv |
| | Direction support | 83.0% | shared_anchor_paired_distance.tsv |
| | H3K27me3 enrichment OR | 2.04 (p=1.75e-24) | shared_anchor_chip_enrichment.tsv |
| **K119ub Integration** | Logistic regression OR | 10.70 (p=1.71e-91) | 3A_k119ub_logistic_regression.tsv |
| | Active anchor rho | -0.314 | 3A_k119ub_correlation_results.tsv |
| **TAD Boundaries** | Differential (Late) | 4,144 (19.0%) | 3D_timepoint_comparison_stats.tsv |
| | Differential (Early) | 4,349 (18.8%) | " |
| | Boundary-loop concordance | 69.6% | shared_boundary_concordance_summary.tsv |
| **Compartments (Late)** | Significant (standard) | 8,189 (7.50% genome) | 1D_compartment_genome_pct_summary.txt |
| **Compartments (Early)** | Significant (standard) | 8,154 (7.46% genome) | compartment_analysis_early/ |
| **Stripes** | Significant (FDR<0.10) | 1 (late), 0 (early) | stripes/outputs/ |
| **ABC Model** | E-G pairs | 179,709 | 4A_delta_abc_all_pairs.tsv |
| | Activity_Lost concordance | 62.0% | 4B_concordance_by_class.tsv |
| | K119ub_Only concordance | 48.7% (chance) | " |
| **DiffBind: K119ub** | Up in mutant (FDR<0.05) | 16,715 | DiffBind file (awk count) |
| | Down in mutant | 5,097 | " |
| **DiffBind: K27ac** | Up / Down | 5,077 / 6,706 | " |
| **DiffBind: K27me3** | Up / Down | 2,293 / 4,811 | " |
| **DiffBind: ATAC** | Up / Down | 9,263 / 4,159 | " |
| **DEG-Loop** | Genes at diff anchors | 457 | 5B_deg_anchor_genes.tsv |
| **HOMER Motifs** | Activity_Lost sig motifs | 369/450 | homer_motif_summary_stats.tsv |

---

## 19. Data Provenance Notes

- All per-resolution edgeR counts were computed from `data/tsvs/figure_2_loop_rewiring/2A_*.tsv` files by counting rows with FDR < 0.05
- All DiffBind counts were computed from `data/upstream/diffbind/` files by counting rows with FDR < 0.05 and Fold > 0 or < 0
- Merged loop counts are row counts from `data/upstream/loop_calls/{early,late}_characterized_loops.tsv`
- Shared anchor direction support was computed from `shared_anchor_paired_distance.tsv` column comparison (median_lost_distance > median_gained_distance)
- Compartment early analysis was regenerated on 2026-04-09 from `tads/tad-pc-analysis/inputs/old/diffcompartments.txt`
- Distance statistics (medians, KS, Wilcoxon, Spearman, chi-square) are from R script output `distance_shift_statistics.txt` generated 2026-01-28
- ABC paired-anchor statistics (rho=0.582, 89.7%, 88.2%, rho=-0.401, 940 genes) are from `abc/docs/` analysis documentation and were not independently re-derived from raw TSVs for this compilation
