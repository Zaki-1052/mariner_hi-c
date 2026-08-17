# Statistical methodology reference for grant Aims 1.1 and 1.2

This document records the statistical analysis frameworks for the two methylation pipelines in the BAP1 project. Each number comes from a verified pipeline output file. Use this document as a reference for adding analytical detail to the grant text at three specific locations.

## Grant locations that need statistical detail

**Aim 1.1 (GRANT.md line 54).** Two sentences. States that H2Aub (CUT&RUN), DNMT3A (ChIP-seq), and 5mC/5hmC (Biomodal evoC) will be profiled at P7, P12, P21, and adult. States that cross-correlation and lag analysis will test whether DNMT3A occupancy follows H2Aub redistribution. Does not specify the DMR calling framework, the statistical model, the covariates, the significance threshold, or how 5mC and 5hmC are jointly interpreted.

**Aim 1.2 (GRANT.md line 55).** Two sentences. States that Enzymatic Methyl-seq will profile mCH across the same time course. Does not specify the differential testing framework, the gene-body aggregation rationale, the edgeR model design, or how conversion noise is controlled.

**Rigor and Reproducibility (GRANT.md line 44).** Names SEACR for peak calling, DiffBind for differential binding, and DESeq2 for transcript abundance. For methylation, states only "standard bisulfite/EM-seq aligners." Does not name the modality XPLR or edgeR frameworks. Mentions "spike-in normalization strategies" without specifying how these are implemented.

## evoC statistical pipeline (Aim 1.1)

### DMR calling: modality XPLR

The Biomodal modality XPLR CLI (v1.3, core-workflow-v1.3.sh) calls differentially methylated regions from the Zarr store produced by the upstream DUET evoC pipeline (v1.5.0). The tool runs a GLM-based test separately for each genomic region file and each methylation context.

Configuration for the primary CG run (config_CG.txt):
- `Condition_Order=control mutant` (control is the reference level)
- `Covariates=sex`
- `Depth_Filter=10` (minimum per-CpG read coverage)
- `Overdispersion=False`

The GLM error family (binomial, quasi-binomial, or other) is not documented. The vendor materials, the platform paper (Fullgrabe et al., 2023, Nat Biotechnol 41:1457-1464), and the modality XPLR workflow documentation do not specify it. The tool is proprietary. Biomodal confirmed that the internals are not publicly available.

### Region files

Six GENCODE vM25 annotation sets define the genomic regions tested:
1. Gene bodies (gencode.vM25.mouse.genes.annotation.bed.gz)
2. Promoters (gencode.vM25.mouse.promoters.annotation.bed.gz)
3. CpG islands (gencode.vM25.mouse.cpg_islands.annotation.bed.gz)
4. CpG shores (gencode.vM25.mouse.cpg_shores.annotation.bed.gz)
5. CpG shelves (gencode.vM25.mouse.cpg_shelves.annotation.bed.gz)
6. TSS regions (gencode.vM25.mouse.tss_region.annotation.bed.gz)

### Independent 5mC and 5hmC testing

The modality XPLR tests 5mC and 5hmC as separate GLM fits per region. Each produces a BED file with 13-14 columns. The columns are: chr, start, end, num_contexts, mean_coverage, mean_mod_group1, mean_mod_group2, mod_fold_change, mod_difference, test_statistic, dmr_pvalue, dmr_qvalue, annotation. Gene-body regions include a 14th column for gene name.

### Significance threshold

Benjamini-Hochberg FDR correction at q < 0.05. This threshold is set in the downstream R pipeline (_shared_config.R, line 176: `Q_THRESHOLD <- 0.05`).

### Coordinated-change classification

The downstream R pipeline joins the 5mC and 5hmC DMR tables by gene name (inner join). Each table is first deduplicated to one row per gene by keeping the row with the lowest q-value. Each co-significant gene is then classified into one of four quadrants by the sign of mod_difference:
- Q4: 5mC increased, 5hmC decreased (the coordinated pattern)
- Q2: 5mC decreased, 5hmC increased (the discordant pattern)
- Q1: both decreased
- Q3: both increased

### Canonical numbers (run-5, verified against PAPER_RESULTS_SYNTHESIS.md)

- 20,969 gene bodies tested
- 10,775 carry a significant 5mC DMR (51.4%): 7,513 hypermethylated (69.7%), 3,262 hypomethylated (30.3%)
- 11,484 carry a significant 5hmC DMR (54.8%): 1,963 increased (17.1%), 9,521 decreased (82.9%)
- 8,371 genes are significant in both contexts
- 6,589 of 8,371 co-significant genes (78.7%) show the coordinated mC-up/hmC-down pattern (Q4)
- 1,255 (15.0%) show the discordant pattern (Q2)

### Bulk genome-wide shift (Section 64)

Paired autosomal analysis across 4 matched sample pairs:
- 5mC: +0.317%, paired t-test p = 0.041, Cohen d = 2.79
- 5hmC: -0.395%, paired t-test p = 0.018, Cohen d = -4.54
- Total modC: unchanged, p = 0.492

The large Cohen d values and flat total modC are consistent with a TET-mediated demethylation block. A de novo methylation surge alone would increase total modC.

## EM-seq statistical pipeline (Aim 1.2)

### Gene-body aggregation

Site-level differential mCH testing is underpowered. The signal (about 0.6-1.3% mCH at gene bodies) is close to the conversion noise floor (about 0.5-0.9% at lambda phage). Gene-body aggregation pools all CH sites within each protein-coding gene body (GENCODE vM25) per sample. This produces a total methylated count and a total coverage count per gene. The aggregation uses GenomicRanges `findOverlaps` in R (script: 02_mch_aggregate_sample.R).

### edgeR quasi-likelihood F-test

The primary statistical engine is the edgeR quasi-likelihood framework (script: 03_mch_differential.R, line 364-374).

Design matrix: `~ genotype + sex + lambda_ch_rate`

- `genotype`: the treatment variable (control vs mutant)
- `sex`: biological covariate (male vs female)
- `lambda_ch_rate`: per-sample spike-in conversion noise rate, extracted from lambda phage CH sites before chromosome filtering

Offset: `log(total_coverage)` per gene per sample. This accounts for differences in sequencing depth across genes and samples without library-size normalization.

Dispersion estimation: `estimateDisp(y, design, robust = TRUE)`

Model fit: `glmQLFit(y, design)`

Test: `glmQLFTest(fit, coef = "genotypemut")`

Significance: Benjamini-Hochberg FDR at q < 0.05.

### Lambda spike-in as covariate

The lambda phage spike-in conversion rate enters the edgeR model as a numeric covariate. Pre-subtraction of counts zeros out too many genes at these low signal levels. The covariate approach adjusts the genotype coefficient for between-sample variation in conversion efficiency.

### All-CH vs CA-only

The primary analysis combines CHH and CHG methylKit files per sample (step 01). This gives about 3 times more sites per gene body than CA dinucleotides alone. A parallel CA-only pipeline exists for comparability to published mCA atlases (Luo et al., 2017). The CA pipeline reuses the same R scripts (02 and 03) with different input paths.

### Verified numbers (from mch_analysis_summary.txt)

- 21,097 testable genes (21,846 total, 749 filtered for low coverage)
- 5,557 genes at FDR < 0.05
- Direction by model-adjusted logFC: 4,915 lose mCH in mutant, 642 gain (7.7:1 ratio)
- Common dispersion: 0.0299 (BCV: 0.1728)
- Tagwise median dispersion: 0.0139
- Genomic inflation factor (lambda_GC): 4.962
- Top gene: Auts2 (logFC = -0.59, F = 1406.6, FDR = 6.49e-20)

The genomic inflation factor of 4.962 reflects the genome-wide nature of BAP1-dependent mCH redistribution. Loss of a chromatin modifier alters methylation at tens of thousands of loci. This is expected biology, confirmed by the directional consistency and the biological coherence of the top genes.

### Lambda conversion rates per sample

| Sample   | Lambda CH rate |
|----------|---------------|
| ctrl_M1  | 1.33%         |
| ctrl_M2  | 0.99%         |
| ctrl_F1  | 1.02%         |
| ctrl_F2  | 1.12%         |
| mut_M1   | 0.87%         |
| mut_M2   | 1.62%         |
| mut_F1   | 1.18%         |
| mut_F2   | 1.81%         |

Genotype-lambda correlation: r = +0.417. The model covariate absorbs this correlation so that it does not inflate the genotype effect.

## Downstream statistical toolkit

Both the evoC and EM-seq pipelines use additional statistical tests in their downstream integration analyses. These tests appear in the 93-section evoC R pipeline and the 25-section EM-seq R pipeline.

**Fisher's exact test.** Used for categorical enrichment: chromatin state by DMR direction, A/B compartment by direction, peak overlap with DiffBind marks. Reported as odds ratio with 95% CI and p-value. Examples: K119ub differential-overlap OR = 4.46 (evoC section 14). mC-hyper in A-compartment OR = 13.64 (section 29). Polycomb exclusion from hypermethylation OR = 0.063 (section 30).

**Multi-mark logistic regression.** GLM with binomial family. Predicts hypermethylation status from multiple DiffBind marks. The 4-mark model gives K119ub OR = 4.71, AUC = 0.818 (section 33). The 6-mark model gives K119ub OR = 10.29 (section 41).

**ROC/AUC with DeLong test.** Compares discriminative ability of competing models. TET-impediment model AUC = 0.793 vs DNMT3A-recruitment model AUC = 0.696, DeLong p = 9.43e-49 (section 24).

**Permutation testing.** Interval-level via regioneReloaded (circular per-chromosome randomization, 5,000 iterations). Gene-level via chromosome-stratified label-shuffle (100,000 iterations, section 37). The gene-level test confirmed 15 of 15 prior Fisher associations (z-scores up to 32).

**GO/KEGG enrichment.** clusterProfiler (enrichGO, gseGO) with BH correction (section 08). The EM-seq pipeline also uses goseq for gene-length correction.

**Correlation and effect size.** Spearman rho for continuous gene-level relationships. Cliff's delta for nonparametric effect size (K119ub BigWig shifts). Wilcoxon rank-sum for two-group comparisons. Cohen's d for the paired bulk methylation shift.

## Sample structure

### evoC (Aim 1.1)

8 samples. 4 control (Bap1-ff), 4 mutant (Bap1-ff-cre, Math1-Cre). 2 batches. 1 male and 1 female per condition per batch. Adult cerebellum (P60). Metadata file: upstream/metadata.csv. The sex covariate is properly powered at n = 4 per genotype. At n = 2 per genotype, sex was confounded and the pipeline found zero significant DMRs.

### EM-seq (Aim 1.2)

8 samples. 4 control, 4 mutant. 2 male and 2 female per genotype. Single batch. Adult cerebellum. Coverage range: 38x to 118x (two deep outliers at 101x and 118x, handled by the log-coverage offset).

### Sample size note

The Rigor and Reproducibility section (line 44) states: "we will maintain n >= 6 per genotype per timepoint." The current evoC and EM-seq datasets both have n = 4 per genotype. This discrepancy needs verification with the PI.

## Cross-modality integration

EM-seq step 04d and downstream section 80_01 test the correlation between mCH (from EM-seq) and CG 5mC/5hmC (from evoC) at the gene level. The analysis uses Spearman rho and quadrant concordance to ask whether genes that gain CG 5mC also gain or lose mCH.

evoC non-CG methylation (CHG/CHH contexts) is not used. The Biomodal DUET upstream pipeline does not faithfully capture non-CG methylation. This limitation is confirmed by the vendor. EM-seq is the sole platform for non-CG profiling in this project.

## Open questions

**Temporal lag analysis.** The grant proposes cross-correlation and lag analysis across the P7/P12/P21/adult time course (Aim 1.1, line 54). This analysis design is not yet specified. The existing pipelines run two-condition comparisons at a single timepoint. The PI will determine the temporal modeling strategy.

**Modality XPLR GLM internals.** The error family and dispersion model are proprietary. The grant can describe the inputs (Zarr store, metadata with covariates, region files) and outputs (per-region BED files with test statistic, p-value, and BH-corrected q-value).

**Sample size.** The grant states n >= 6 per genotype per timepoint. Current methylation data has n = 4 per genotype. Whether the time-course experiments will scale to n = 6 needs confirmation.

## Key citations

- Fullgrabe J et al. (2023). Simultaneous sequencing of genetic and epigenetic bases in DNA. Nature Biotechnology 41:1457-1464. [The evoC platform paper. Describes the sequencing chemistry. Does not describe the modality XPLR analysis software.]
- Robinson MD, McCarthy DJ, Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26:139-140. [The edgeR framework used for mCH differential testing.]
- Lun ATL, Chen Y, Smyth GK (2016). It's DE-licious: A recipe for differential expression analyses of RNA-seq experiments using quasi-likelihood methods in edgeR. Methods in Molecular Biology 1418:391-416. [The quasi-likelihood F-test methodology.]
- Vaisvila R et al. (2021). Enzymatic methyl-seq: a next generation methyltransferase inhibitor-free and bisulfite-free methylation analysis tool. Genome Research. [The EM-seq chemistry.]
