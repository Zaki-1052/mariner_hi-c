# Hi-C statistical methodology reference

This document records the analysis parameters for the chromatin loop and stripe pipelines in the BAP1 Hi-C paper. Each number comes from a verified script, config file, conda environment YAML, or pipeline output log. The source file for each number is named in the text. Use this document as the single source of truth when writing the methods section.

Challana's request: "Can you write a section up for the loop analysis? It'll need details on the loop caller we used, what resolutions we used, our mariner & edgeR parameters, then how we did our loop annotation. Try to have all the software versions included b/c some journals are real sticklers about it. Also for stripe detection and differential analysis."

## Upstream Hi-C processing

Deep-sequenced Hi-C libraries were trimmed with fastp (v1.0.1) and aligned to the mm10 mouse genome with the nf-core/hic pipeline (v2.1.0) run on Nextflow (v25.04.6). The pipeline uses Bowtie2 (v2.5.4) for alignment with the following end-to-end parameters: `--very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder`. The digestion protocol was set to Arima. The pipeline produced valid-pair files, BAM files, and cooler-format contact maps at bin sizes 1000000, 500000, 250000, 100000, 40000, 25000, 10000, 5000, 2000, and 1000 bp. TADs were called with insulation and hicExplorer at 40, 20, and 10 kb. Compartments were computed at 1000, 250, and 100 kb. The pipeline ran in Singularity containers on SDSC Expanse.

BAM files were processed with SAMtools (v1.22). Valid pairs were converted to Juicer v2.0 format .hic files with juicerpre, which computes VC, VC_SQRT, GW_SCALE, SCALE, and INTER_SCALE normalization vectors.

Source: `Hi-C_Pipeline.md`, `envs/mariner_env.yml` (Bowtie2 2.5.4, fastp 1.0.1, SAMtools 1.22). Nextflow log (nf-core/hic v2.1.0, Nextflow v25.04.6, revision fe4ac65631).

The upstream pipeline was developed by Challana Tea. The loop and stripe analysis pipelines described below were developed by Zakir Alibhai.

## Loop detection: HiCCUPS

Chromatin loops were called with HiCCUPS from Juicer Tools v2.0.1, run inside a Singularity container (`juicer_2.0.1.sif`) on GPU nodes. The HiCCUPS command used default parameters with no explicit flags. Default HiCCUPS settings produce loop calls at three resolutions: 5,000 bp (5 kb), 10,000 bp (10 kb), and 25,000 bp (25 kb). Each of the 6 individual replicate .hic files was processed separately.

Source: `lab/pipeline-scripts/juicer.sb` (line 30: container path, line 50-51: command with no additional flags).

## Mariner loop processing pipeline

Loop calls from 6 replicates were processed with the mariner R/Bioconductor package (v1.8.1) at each resolution independently.

### Step 1: Loop merging (prep_loops.R)

HiCCUPS BEDPE files from 6 replicates were loaded as GInteractions objects. Loops were merged across replicates with `mergePairs()` using a 10,000 bp (10 kb) radius, Manhattan distance metric, and `selectMax = TRUE` on the observed contact count. This union strategy retains a loop if any replicate detected it, and selects the representative position by the highest observed count.

Merged loops were assigned to resolution-specific genomic bins with `assignToBins(binSize = RESOLUTION, pos1 = "center", pos2 = "center")`. Each loop was then buffered with `pixelsToMatrices(buffer = 2)` to define a 5x5 pixel extraction window around each loop center.

chrX was excluded from analysis (`edgeR_config.yaml` filtering.exclude_chromosomes: `["chrX"]`).

Source: `loops/scripts/prep_loops.R` (lines 134-140: mergePairs call, lines 181-186: assignToBins, lines 198-201: pixelsToMatrices). `loops/config/edgeR_config.yaml` (line 91 in paths_config: chrX exclusion).

### Step 2: Hi-C count extraction (extract_counts.R)

5x5 observed contact matrices were extracted from 6 .hic files at each buffered loop position with `pullHicMatrices()`. The extraction used VC (Vanilla Coverage) normalization, observed matrix values, on-disk HDF5 storage (`onDisk = TRUE`, `compressionLevel = 1`), and a block size of 1,000,000.

The script explicitly sets `norm_to_use <- "VC"` with the comment "Use VC normalization (KR not available)."

Source: `loops/scripts/extract_counts.R` (lines 80-81: VC normalization, lines 121-131: pullHicMatrices call).

### Step 3: Aggregation (aggregate.R)

Each 5x5 matrix was collapsed to a single count per loop per sample by summing all 25 pixels. This simple sum was the primary aggregation method. A center-weighted alternative (Gaussian-like kernel, center weight 0.16) was computed for comparison but not used for downstream analysis.

A DGEList object was constructed with `group = factor(c("ctrl","ctrl","ctrl","mut","mut","mut"))`.

Source: `loops/scripts/aggregate.R` (line 75: sum, lines 233-234: primary method selection, lines 353-357: DGEList).

### Step 4: Quality control (qc-val.R)

Nine QC analyses were run per resolution. Key thresholds: within-group Pearson correlation above 0.90 (good) or 0.95 (excellent). Library size ratio below 2.0-fold. No negative count values. Biological signal detected when within-group correlation exceeds between-group correlation by more than 0.05. Positional shift detection flagged loops where the peak-signal pixel differed between conditions.

Source: `loops/scripts/qc-val.R`, `loops/config/edgeR_config.yaml` (lines 76-79: QC thresholds).

## edgeR differential loop analysis

### Model specification

The differential analysis used the edgeR quasi-likelihood GLM (v4.4.0 on HPC, v4.6.3 on local Mac for visualization).

**Filtering.** Low-count loops were removed with `filterByExpr(y, group, min.count = 5, min.total.count = 15, min.prop = 0.5)`.

**Normalization.** TMM (Trimmed Mean of M-values) normalization was applied with `normLibSizes(y, method = "TMM")`. This is the edgeR >= 4.0 function name for the operation previously called `calcNormFactors()`.

**Design matrix.** `design <- model.matrix(~group, data = y$samples)` with `group = factor(c("ctrl","mut"), levels = c("ctrl","mut"))`. The design has two columns: Intercept and MutantEffect. No additional covariates (no sex, no batch).

**Dispersion estimation.** `estimateDisp(y, design, robust = TRUE)`. This estimates common, trended, and tagwise dispersions from the biological replicates.

**Model fit.** `glmQLFit(y, design, robust = TRUE)`. Quasi-likelihood GLM with robust estimation.

**Test.** `glmQLFTest(fit, coef = 2)`. Tests coefficient 2 ("MutantEffect"), the mutant-vs-control contrast.

**FDR correction.** Benjamini-Hochberg, applied by `topTags()` (edgeR default).

**Primary significance threshold.** FDR < 0.05.

**Exploratory threshold.** FDR < 0.10.

**Final stringent filter.** |logFC| > 0.3 AND FDR < 0.03. Applied in both `edgeR.R` (line 783) and `generate_final_results.py`. This produces the non-redundant publication loop set.

**No offset term.** The model uses edgeR's built-in effective library size (lib.size x norm.factors) from TMM normalization. No additional log-offset was specified.

Source: `loops/scripts/edgeR.R` (lines 213-226: DGEList, lines 253-259: filterByExpr, lines 278-289: normLibSizes, lines 330-336: design, line 340: estimateDisp, line 369: glmQLFit, line 391: glmQLFTest, lines 782-791: stringent filter). `loops/config/edgeR_config.yaml` (lines 36-59: all statistical parameters).

### Effect-size categories

| Category | Threshold | Linear fold change |
|----------|-----------|-------------------|
| Strong | |logFC| > 1.0 | > 2-fold |
| Moderate | |logFC| > 0.5 | > 1.4-fold |
| Weak | |logFC| > 0.0 | any change |

Source: `loops/config/edgeR_config.yaml` (lines 54-56).

## Multi-resolution comparison

Loops were matched across the three resolutions (5 kb, 10 kb, 25 kb) by requiring the same chromosome pair and both anchor start positions within 10,000 bp (10 kb) of each other. The closest match by summed anchor distance was selected. Both anchors had to match within tolerance.

Source: `loops/scripts/compare_resolutions.R` (line 239: tolerance_kb = 10).

## Loop anchor annotation

### 8-category priority system

Each loop anchor was classified into one of 8 chromatin-state categories based on overlap with CUT&RUN peak files. The classification uses a priority system: the first matching category is assigned. Five histone marks and CTCF motifs define the categories.

**Peak files used (late/adult timepoint):**
- H3K27ac: `H3K27acCerebellumLate2.bed` (15,105 regions)
- H3K27me3: `H3K27me3CerebellumLate1.bed` (15,809 regions)
- H3K4me1: `H3K4me1CerebellumLate1.bed` (113,781 regions)
- H3K4me3: `H3K4me3CerebellumLate2.bed` (6,581 regions)
- Bivalent: `Bivalent_Cerebellum_Late.bed` (318 regions, pre-computed K4me3+K27me3 intersection)
- CTCF: `ctcf_motifs_mm10.bed` (114,081 motif sites)

Peak region counts are from `loops/output/master_loop_table/late/build_master.log`.

**Priority order:**

| Priority | Category | Rule |
|----------|----------|------|
| 1 | Active_Promoter | H3K4me3+ AND NOT H3K27me3 AND distance to TSS <= 2,000 bp |
| 2 | Repressed_Promoter | H3K27me3+ AND NOT H3K27ac AND distance to TSS <= 2,000 bp |
| 3 | Bivalent_Promoter | Overlaps pre-computed K4me3+K27me3 intersection region |
| 4 | Polycomb | H3K27me3+ AND distance to TSS > 2,000 bp |
| 5 | Active_Enhancer | H3K27ac+ AND distance to TSS > 2,000 bp |
| 6 | Poised_Enhancer | H3K4me1+ AND NOT H3K27ac AND NOT H3K27me3 AND distance to TSS > 2,000 bp |
| 7 | CTCF_Site | Overlaps CTCF DNA motif AND not classified above |
| 8 | Other | No marks |

**TSS distance threshold:** 2,000 bp. TSS positions from `TxDb.Mmusculus.UCSC.mm10.knownGene` (v3.10.0).

**Gene annotation:** GENCODE vM25 (`gencode.vM25.mouse.genes.annotation.bed.gz`). Nearest gene assigned per anchor via `distanceToNearest()`.

**Loop type:** Each loop receives a pairwise label from its two anchor types (for example, Active_Promoter-Active_Enhancer). This produces up to 36 combinations.

Source: `loops/scripts/annotate_loops_extended.R` (lines 233-305: classify_anchor_type_extended function), `loops/scripts/build_master_loop_table.R` (lines 138-209: reproduced classification logic). Paper figures 2F and 2G were generated by `peaks/scripts/annotate_loops_extended.R` (same logic, confirmed by `data/INDEX.md`).

### Loop class (structural/CRE/mixed)

Loops were also classified by anchor composition into 4 classes based on CTCF and enhancer/promoter overlap:
- Structural: both anchors overlap CTCF AND fewer than 2 anchors overlap H3K27ac or a promoter
- CRE: both anchors overlap H3K27ac or a promoter AND fewer than 2 overlap CTCF
- Mixed: both anchors overlap both CTCF and H3K27ac/promoter
- Unclassified: all other loops

Source: `loops/scripts/build_master_loop_table.R` (step 5), reproducing `cluster/scripts/05_grouped_analyses.py` section 4.2.

## APA analysis

Aggregate Peak Analysis (APA) was computed with GENOVA. GENOVA operates on cooler-balanced cool files or Juicer v2.0 .hic files and computes pixel enrichment around loop coordinates from BEDPE input. For the paper, merged ctrl and merged mutant contact maps were used (single per-condition maps, not per-replicate). GENOVA produces side-by-side pyramid plots, differential pyramid plots, and enrichment heatmaps.

GENOVA is not replicate-aware. It performs single ctrl-vs-mutant comparisons on the merged contact maps.

Source: `Hi-C_Pipeline.md` (GENOVA APA section).

Note: A separate mariner-based APA script (`loops/scripts/apa_analysis.R`) also exists in the repo. It uses `pullHicMatrices()` with a 50 kb window, mean aggregation, P2LL enrichment (center/corner ratio), and Wilcoxon + Cliff's delta statistics across individual replicates. The GENOVA results are the ones used for the paper figures.

## GO/KEGG enrichment

Gene Ontology and KEGG pathway enrichment were computed with clusterProfiler (v4.16.0). Entrez gene IDs were identified within 10,000 bp of up-regulated or down-regulated loop anchors. Enrichment tests used `enrichGO()` for three GO ontologies (BP, CC, MF) and `enrichKEGG()` for KEGG pathways, with `org.Mm.eg.db` (v3.21.0) and BH correction at p < 0.05 and q < 0.05.

Source: `loops/scripts/visualizations.R` (lines 613-647: gene universe, lines 660-746: enrichment calls).

## TAD and compartment analysis

Two independent systems produced the paper's TAD and compartment results. HOMER handled consensus TAD calling, differential inclusion ratios, and compartment (PC1) analysis (Figure 1B TAD volcanos, Figure 1D compartment volcanos). TADCompare handled differential boundary detection, boundary shift analysis, and boundary-loop cross-reference (Figure 1C differential boundaries, Figure 1E Syt1/Nav3 locus, Figure 1F boundary-loop cross-ref and ChIP classification).

### HOMER TAD and compartment analysis

TAD boundaries and inclusion ratios were computed with HOMER (v5.1). Tag directories were generated from Hi-C BAM files. The HOMER TAD pipeline had 4 steps:

1. `findtad` per sample to call TADs from individual tag directories
2. `mergetad2` to merge per-sample TADs into a consensus set (`merged.tad.2D.bed`)
3. `scoretad` to compute inclusion ratios (interaction density) per TAD per sample
4. `getDiffExpression.pl` for differential TAD analysis

Compartment (PC1) values were computed per sample at 25x50kb resolution, stored as bedGraph files. Differential compartment analysis used `getDiffExpression.pl` on the PC1 values.

TAD volcano plots were generated with `scripts/tad_volcano_plot.R` at two threshold levels: relaxed (FDR < 0.15, |Difference| > 0.15) and standard (FDR < 0.05, |Difference| > 0.30). Compartment volcano plots used the same dual-threshold approach (`scripts/compartment_volcano_plot.R`).

The compartment volcano script includes PC1 polarity auto-correction: it validates that genic bins have higher PC1 than intergenic bins and negates values if flipped. Positive Difference = shift toward A compartment (more active) in mutant. Negative Difference = shift toward B compartment (more repressed).

Genome percentage of compartment shifts was computed by `scripts/compartment_genome_percentage.R` against the mm10 assembled genome (2,730,871,774 bp).

Source: `Hi-C_Pipeline.md` (HOMER section, lines 198-231), `loops/CLAUDE.md` (TAD analysis section).

### TADCompare differential boundary analysis

Differential TAD boundaries were detected with TADCompare (Cresswell & Dozmorov, Frontiers in Genetics 11:158, 2020). TADCompare version was not pinned in conda environments and needs to be confirmed from HPC.

**Pipeline: 6 steps.**

**Step 1: Matrix extraction (`tads/scripts/01_extract_matrices.sb`).** Sparse contact matrices were extracted from .hic files with Juicer `straw` at 10,000 bp (10 kb) resolution. The extraction used NONE normalization (`straw observed NONE`). Matrices were extracted for merged ctrl and merged mutant files, and for all 6 individual replicates. Chromosomes: chr1-chr19 autosomes only (chrX excluded).

**Step 2: TADCompare (`tads/scripts/02_run_tadcompare.R`).** Differential boundary detection on merged ctrl vs merged mutant sparse matrices per chromosome. The `TADCompare()` function was called with `cont_mat1 = ctrl_mat, cont_mat2 = mut_mat, resolution = 10000`. No `pre_tads` argument was used; TADCompare detected boundaries internally. The output `TAD_Frame` contains per-boundary columns: `Boundary` (genomic coordinate), `Gap_Score` (differential boundary score, Z-score), `TAD_Score1` and `TAD_Score2` (boundary scores in control/mutant), `Differential` ("Differential" or "Non-Differential"), `Enriched_In` ("Matrix 1" = control or "Matrix 2" = mutant), and `Type` (one of 6 classifications: Non-Differential, Strength Change, Shifted, Split, Merge, Complex).

Source: `tads/scripts/02_run_tadcompare.R` (lines 90-95: TADCompare call, line 25: resolution, line 30: chromosomes).

**Step 3: ConsensusTADs (`tads/scripts/03_run_consensus.R`).** Replicate robustness assessment across 3 biological replicates per condition. Resolution 10 kb, chromosomes chr1-chr19.

Source: `tads/scripts/03_run_consensus.R` (line 25: resolution, line 31: chromosomes).

**Step 4: Shift analysis (`tads/scripts/04_analyze_shifts.R`).** Shift distances were computed with `bedtools closest -d` for boundaries classified as "Shifted." Only "Shifted" boundaries received a shift_distance value; all other types were set to NA. Binned distribution: 0-25 kb, 25-50 kb, 50-100 kb, 100-200 kb, 200-500 kb, >500 kb.

A filtered shift analysis (`tads/scripts/shift_analysis_filtered.R`) used only shifts >= 10 kb (1 bin at 10 kb resolution). Input: `results/{timepoint}/final/shifted_boundaries_min1bin.bed`.

Source: `tads/scripts/04_analyze_shifts.R` (lines 97-100: Shifted-only assignment, lines 126-130: bins), `tads/scripts/shift_analysis_filtered.R` (line 4: >= 10 kb threshold).

**Step 5: Blacklist filtering (`tads/scripts/05_filter_blacklist.R`).** Two blacklist files were applied as a union:
- Lab-specific blacklist: `tads/250123blacklist.bed` (254 regions)
- ENCODE mm10 blacklist: `tads/mm10-blacklist.v2.bed` (3,435 regions)

The two lists were merged with `reduce(c(blacklist_lab, blacklist_encode))`. Any boundary whose 10 kb bin overlapped any blacklisted region was removed.

Source: `tads/scripts/05_filter_blacklist.R` (lines 24-26: blacklist files, lines 55-57: union, lines 64-78: filtering logic).

**Step 6: Visualization (`tads/scripts/tad_visualizations.R`).** 40+ publication plots across 9 categories.

### TAD volcano plot dual-format support

The TAD volcano script (`scripts/tad_volcano_plot.R`) auto-detects two input formats: HOMER `getDiffExpression.pl` output and TADCompare output. For TADCompare input, p-values were computed from Gap_Score as Z-scores: `2 * pnorm(-abs(Gap_Score))`, then BH FDR correction. TADCompare thresholds: relaxed (FDR < 0.15, |Gap_Score Z| > 2.0), standard (FDR < 0.05, |Gap_Score Z| > 3.0).

Source: `scripts/tad_volcano_plot.R` (lines 146-156: format detection, lines 199-202: Z-score conversion, lines 221-224: TADCompare thresholds).

### Boundary-loop cross-reference

Differential loops were tested for preferential localization near differential TAD boundaries (`tads/scripts/boundary_loop_crossref.R`). Distance thresholds for "near boundary" classification: 10,000 bp, 25,000 bp, 50,000 bp, and 100,000 bp. Statistical tests: Fisher's exact test, bootstrap (1,000 iterations), and permutation (1,000 iterations). FDR threshold: 0.05.

Source: `tads/scripts/boundary_loop_crossref.R` (lines 55-56: distance thresholds, lines 58-60: statistical parameters).

### Boundary ChIP-seq classification

TAD boundaries were classified with the same 8-category chromatin state system as loop anchor annotation (`tads/scripts/tad_chip_classification.R`). The system uses the same 5 histone marks and the same 2,000 bp TSS distance threshold.

Source: `tads/scripts/tad_chip_classification.R` (lines 1-9: description).

### Verified TADCompare numbers

**Late timepoint (250402):**

| Metric | Count |
|--------|-------|
| Total boundaries | 8,264 |
| Differential | 1,384 (16.7%) |
| Shifted | 254 (3.1%) |
| Strength Change | 714 (8.7%) |
| Complex | 363 (4.4%) |
| Merge | 193 (2.3%) |
| Split | 112 (1.4%) |
| Non-Differential | 6,607 (79.9%) |
| Control-enriched | 532 (38.4% of differential) |
| Mutant-enriched | 852 (61.6% of differential) |

**Early timepoint (250831):**

| Metric | Count |
|--------|-------|
| Total boundaries | 8,584 |
| Differential | 1,222 (14.2%) |
| Shifted | 278 (3.2%) |
| Strength Change | 620 (7.2%) |
| Complex | 282 (3.3%) |
| Split | 193 (2.3%) |
| Merge | 123 (1.4%) |
| Non-Differential | 7,065 (82.5%) |
| Control-enriched | 784 (64.2% of differential) |
| Mutant-enriched | 438 (35.8% of differential) |

**HOMER TAD counts:**

| Metric | Late (Adult) | Early (P13) |
|--------|-------------|-------------|
| Consensus TADs | 4,366 | 5,074 |
| Significant relaxed (FDR < 0.15, |Diff| > 0.15) | 1,283 | 45 |
| Significant standard (FDR < 0.05, |Diff| > 0.30) | 355 | 0 |

Source for TADCompare numbers: `tads/README.md` (lines 535-566). Source for HOMER numbers: `data/INDEX.md` (lines 63-68).

### TADCompare version

TADCompare v1.16.0. Not pinned in conda environment YAMLs; installed via Bioconductor. Version confirmed by user (`packageVersion("TADCompare")`).

---

## Stripe detection: Stripenn

### Format conversion (Stage 0)

.hic files were converted to multi-resolution .mcool format with hictk (v2.1.4) at resolutions 5,000 bp and 10,000 bp. One set of merged .hic files (250402 timepoint, ctrl_merged and mut_merged) lacked native KR normalization vectors. ICE balancing weights were computed on those two .mcool files with `cooler balance --name KR` and stored under the "KR" column name. All other .mcool files retained native KR vectors from the source .hic files.

Source: `stripes/stripenn/scripts/00_convert_hic_to_cool.sb`, `envs/hictk.yml` (hictk 2.1.4), `stripes/stripenn/README.md` (section 5.2: 250402 KR fix).

### Stripe calling (Stage 1)

Chromatin stripes were called with Stripenn (v1.1.65.22) on merged ctrl and merged mutant .mcool files at both 5 kb and 10 kb resolution. Stripenn uses Canny edge detection and pixel saturation thresholds to identify architectural stripes.

**Parameters:**

| Parameter | Value | Stripenn flag |
|-----------|-------|---------------|
| Normalization | KR | `--norm KR` |
| Chromosomes | all (sex/mito filtered post-hoc) | `-k all` |
| Canny sigma | 2.0 | (default) |
| Minimum length | 10 bins | `-l 10` |
| Maximum width | 8 bins | `-w 8` |
| Maxpixel percentiles | 0.95, 0.96, 0.97, 0.98, 0.99 | `-m 0.95,0.96,0.97,0.98,0.99` |
| P-value threshold | 0.1 | `-p 0.1` |
| Background filter | 3 | `-b 3` |
| Seed | 42 | `--seed 42` |

The default Stripenn seed is 123456789. The pipeline uses 42.

chrX, chrY, and chrM stripes were removed post-hoc in R (Stage 2).

Source: `stripes/stripenn/scripts/01_call_stripes.sb` (lines 64-75: exact command), `stripes/stripenn/config/stripenn_config.yaml` (lines 59-68: parameters).

### Union set construction (Stage 2)

Ctrl and mutant stripe calls were merged into a union set per timepoint per resolution. Anchor matching used a 50,000 bp (50 kb) tolerance with direction-aware overlap detection. A ctrl stripe and a mutant stripe were considered "shared" if their anchors overlapped within 50 kb and they had the same directionality (3-prime or 5-prime). Unmatched stripes were labeled control_only or mutant_only.

Source: `stripes/stripenn/scripts/02_build_union.R`, `stripes/stripenn/config/stripenn_config.yaml` (line 74: anchor_tolerance_bp = 50000).

### Replicate scoring (Stage 3)

Each of 6 individual replicate .mcool files was scored against the union BEDPE with `stripenn score`. The key output column is `O_Sum_added`: the sum of observed contacts per stripe per sample. Row counts were verified to match the union BEDPE exactly.

Source: `stripes/stripenn/scripts/03_score_replicates.sb` (lines 72-78: exact command).

## Stripe differential analysis (Stage 4)

### edgeR model for stripes

The differential analysis used the same edgeR quasi-likelihood GLM as the loop pipeline, with two differences.

**Difference 1:** `filterByExpr()` was skipped. The stripe feature sets were small (1,800 to 7,400 stripes per timepoint-resolution combination). The config sets `skip_filtering: true`.

**Difference 2:** `calcNormFactors()` was used (not the newer `normLibSizes()` alias).

**Design matrix.** `design <- model.matrix(~group, data = dge$samples)` with `group = factor(c("ctrl","mut"), levels = c("ctrl","mut"))`. Two columns: Intercept and MutantEffect.

**Dispersion.** `estimateDisp(dge, design, robust = TRUE)`.

**Fit.** `glmQLFit(dge, design, robust = TRUE)`.

**Test.** `glmQLFTest(fit, coef = 2)`. BH FDR correction by `topTags()` default.

**Primary threshold.** FDR < 0.05.

Source: `stripes/stripenn/scripts/04_edgeR.R` (line 163: DGEList, lines 186-193: skip filtering, line 200: calcNormFactors, lines 234-235: design, line 237: estimateDisp, line 256: glmQLFit, line 267: glmQLFTest). `stripes/stripenn/config/stripenn_config.yaml` (lines 79-84: edgeR parameters).

## Stripe classification (Stage 5)

Each stripe was classified by direction:
- control_only stripes were labeled "lost"
- mutant_only stripes were labeled "gained"
- shared stripes with FDR < 0.05 and logFC > 0.3 were labeled "strengthened"
- shared stripes with FDR < 0.05 and logFC < -0.3 were labeled "weakened"
- all other shared stripes were labeled "unchanged"

Confidence tiers were assigned based on FDR and logFC agreement:
- High: control_only with FDR < 0.10 and logFC < 0, or mutant_only with FDR < 0.10 and logFC > 0, or shared with FDR < 0.05
- Medium: control_only with logFC < -0.2, or mutant_only with logFC > 0.2, or shared with FDR < 0.10
- Low: all other stripes

A directional consistency flag recorded whether the edgeR logFC agreed with the detection-based source label (for example, control_only stripes with logFC < 0 are consistent).

Source: `stripes/stripenn/scripts/05_integration.R` (lines 95-101: direction, lines 103-115: confidence, lines 117-123: consistency).

## Stripe cross-resolution comparison (Stage 6)

5 kb and 10 kb stripe results were matched per timepoint. Anchors were extended by 100,000 bp (2 x 50 kb tolerance) and matched with `findOverlaps()`. One-to-one matching was enforced by deduplication. Concordant stripes were significant at FDR < 0.05 at both resolutions with the same direction label.

Pearson correlation of logFC values was computed on matched stripes.

Source: `stripes/stripenn/scripts/06_compare_resolutions.R`.

## Stripe anchor annotation (Stage 7)

Stripe anchors were classified with a 7-category system using the same 5 histone marks as the loop annotation (H3K27ac, H3K27me3, H3K4me1, H3K4me3, Bivalent) and the same 2,000 bp TSS threshold.

**Priority order (stripes):**

| Priority | Category | Rule |
|----------|----------|------|
| 1 | Bivalent_Promoter | Overlaps Bivalent region AND distance to TSS <= 2,000 bp |
| 2 | Active_Promoter | H3K4me3+ AND H3K27ac+ AND distance to TSS <= 2,000 bp |
| 3 | Repressed_Promoter | H3K27me3+ AND NOT H3K27ac AND distance to TSS <= 2,000 bp |
| 4 | Active_Enhancer | H3K27ac+ AND distance to TSS > 2,000 bp |
| 5 | Poised_Enhancer | H3K4me1+ AND NOT H3K27ac AND distance to TSS > 2,000 bp |
| 6 | Polycomb | H3K27me3+ AND distance to TSS > 2,000 bp |
| 7 | Other | No marks |

The stripe annotation does not include a CTCF_Site category (7 categories total vs 8 for loops). The Active_Promoter rule in the stripe pipeline requires both H3K4me3+ and H3K27ac+, while the loop pipeline requires only H3K4me3+ AND NOT H3K27me3. Bivalent_Promoter is checked first in the stripe pipeline (priority 1) versus third in the loop pipeline (priority 3).

Source: `stripes/stripenn/scripts/stripenn_visualizations.R` (lines 738-746: classification logic).

---

## Sample structure

### Loops

6 samples: ctrl_M1, ctrl_M2, ctrl_M3, mut_M1, mut_M2, mut_M3. 3 biological replicates per condition. BAP1-KO (Bap1f/f,Math1-cre) mutant vs wildtype control. Mouse cerebellum, mm10 genome.

Two timepoints:
- Early: P12/P13 developmental stage (batch 250831)
- Late: adult (batch 250402)

Sex composition: to be filled in (not recorded in analysis configs).

### Stripes

Same 6 replicates as loops, scored individually. Stripe calling was done on merged samples (ctrl_merged, mut_merged) per condition per timepoint. Two resolutions: 5 kb and 10 kb.

---

## Verified numbers

### Loop counts (late/adult timepoint)

| Stage | Count | Source |
|-------|-------|--------|
| Total loops tested (all resolutions merged) | 39,344 | build_master.log |
| Significant at FDR < 0.05 | 7,981 | build_master.log |
| Up in mutant (of 7,981) | 4,253 | build_master.log |
| Down in mutant (of 7,981) | 3,728 | build_master.log |
| Non-redundant final set (|logFC| > 0.3 + FDR < 0.03) | 2,910 | non_redundant_loops.tsv row count |
| Structural class | 3,537 | build_master.log |
| CRE class | 1,294 | build_master.log |
| Mixed class | 599 | build_master.log |
| Unclassified | 2,551 | build_master.log |
| Super-enhancer anchor overlap | 688 (8.6%) | build_master.log |
| Shared/switched anchor | 2,250 (28.2%) | build_master.log |
| DEG overlap | 2,297 (28.8%) | build_master.log |

Source: `loops/output/master_loop_table/late/build_master.log`.

### Loop counts (early/P13 timepoint)

| Stage | Count | Source |
|-------|-------|--------|
| Total loops tested | 39,124 | merged_all_results.tsv row count |
| Non-redundant final set | 165 | non_redundant_loops.tsv row count |

Source: `loops/outputs/250831-early_outputs/merged_loops/`.

### Stripe counts (late/adult, 5 kb)

| Stage | Count | Source |
|-------|-------|--------|
| Ctrl stripes called | 5,327 | stripenn README.md section 6.2 |
| Mut stripes called | 5,850 | stripenn README.md section 6.2 |
| Union set | 7,371 | stripenn README.md section 6.2 |
| Significant at FDR < 0.05 | 2,320 (31.5%) | stripenn README.md section 6.2 |
| Lost | 1,528 | stripenn README.md section 6.2 |
| Gained | 2,052 | stripenn README.md section 6.2 |
| High-confidence (lost + gained) | 1,005 | stripenn README.md section 6.3 |
| Max |logFC| | 0.389 | stripenn README.md section 6.3 |
| BCV | 0.012 | stripenn README.md section 6.3 |
| Cross-res logFC correlation (5 kb vs 10 kb) | r = 0.850 | stripenn README.md section 6.3 |

### Stripe counts (early/P12, 5 kb)

| Stage | Count | Source |
|-------|-------|--------|
| Union set | 4,008 | stripenn README.md section 6.2 |
| Significant at FDR < 0.05 | 96 (2.4%) | stripenn README.md section 6.2 |
| Max |logFC| | 0.321 | stripenn README.md section 6.3 |
| BCV | 0.020 | stripenn README.md section 6.3 |
| Cross-res logFC correlation | r = 0.808 | stripenn README.md section 6.3 |

Source for all stripe numbers: `stripes/stripenn/README.md` sections 6.2-6.3.

---

## Software versions

### Upstream Hi-C processing

| Software | Version | Source |
|----------|---------|--------|
| nf-core/hic | 2.1.0 | Nextflow log (revision fe4ac65631) |
| Nextflow | 25.04.6 | Nextflow log |
| fastp | 1.0.1 | `envs/mariner_env.yml` line 51 |
| Bowtie2 | 2.5.4 | `envs/mariner_env.yml` line 43 |
| SAMtools | 1.22 | `envs/mariner_env.yml` line 192 |
| BEDTools | 2.31.1 | `envs/mariner_env.yml` line 12 |
| Juicer Tools (HiCCUPS) | 2.0.1 | `lab/pipeline-scripts/juicer.sb` line 30 |
| hictk | 2.1.4 | `envs/hictk.yml` line 13 |
| HOMER | v5.1 | User-confirmed (2024-07-16 release) |
| GENOVA | 1.0.1 | HPC `genova_env` (R 4.4.2, `packageVersion('GENOVA')`) |

### R packages (edgeR models ran on HPC: R 4.4.2; visualization/annotation ran locally: R 4.5.2)

| Package | Version (local Mac) | Version (HPC mariner_env) | Source |
|---------|--------------------|-----------------------------|--------|
| R | 4.5.2 | 4.4.2 | `software_versions.txt`, `envs/mariner_env.yml` |
| Bioconductor | 3.21 | 3.20 (inferred from R 4.4.2) | `software_versions.txt` |
| mariner | 1.8.1 | (not in HPC env) | `software_versions.txt` |
| edgeR | 4.6.3 | 4.4.0 | `software_versions.txt`, `envs/mariner_env.yml` line 23 |
| limma | 3.64.3 | 3.62.1 | `software_versions.txt`, `envs/mariner_env.yml` line 29 |
| InteractionSet | 1.36.1 | — | `software_versions.txt` |
| GenomicRanges | 1.62.1 | 1.58.0 | `software_versions.txt`, `envs/mariner_env.yml` line 27 |
| strawr | 0.0.92 | — | `software_versions.txt` |
| HDF5Array | 1.36.0 | — | `software_versions.txt` |
| ChIPseeker | 1.44.0 | — | `software_versions.txt` |
| clusterProfiler | 4.16.0 | — | `software_versions.txt` |
| EnhancedVolcano | 1.26.0 | — | `software_versions.txt` |
| TxDb.Mmusculus.UCSC.mm10.knownGene | 3.10.0 | — | `software_versions.txt` |
| org.Mm.eg.db | 3.21.0 | — | `software_versions.txt` |
| TADCompare | 1.16.0 | 1.16.0 | User-confirmed (`packageVersion("TADCompare")`) |
| rtracklayer | 1.70.1 | — | `software_versions.txt` |

### Python / Stripenn (HPC mariner_env)

| Package | Version | Source |
|---------|---------|--------|
| Python | 3.13.5 | `hpc_versions.txt` |
| Stripenn | 1.1.65.22 | `hpc_versions.txt`, `envs/mariner_env.yml` line 269 |
| cooler | 0.10.4 | `hpc_versions.txt`, `envs/mariner_env.yml` line 232 |
| NumPy | 2.4.1 | `hpc_versions.txt` |
| SciPy | 1.17.1 | `hpc_versions.txt` |
| pandas | 2.3.3 | `hpc_versions.txt` |
| OpenCV | 4.13.0 | `hpc_versions.txt` |
| scikit-image | 0.26.0 | `hpc_versions.txt` |
| hic-straw (Python) | 1.3.1 | `envs/mariner_env.yml` line 238 |

### Genome and annotation

| Resource | Version |
|----------|---------|
| Genome assembly | mm10 (GRCm38) |
| Gene annotation (R) | TxDb.Mmusculus.UCSC.mm10.knownGene v3.10.0 |
| Gene annotation (regions) | GENCODE vM25 |
| Organism database | org.Mm.eg.db v3.21.0 |

---

## Notes for Challana

**Sex composition of Hi-C samples.** The sample names (ctrl_M1, ctrl_M2, ctrl_M3, mut_M1, mut_M2, mut_M3) do not encode sex. The methods text does not need to state sex composition. If Challana wants to include it, the information is in the lab records.

**HPC vs Mac R versions.** The edgeR differential models ran on HPC with R 4.4.2 and edgeR 4.4.0. Visualization and annotation scripts ran locally with R 4.5.2 and edgeR 4.6.3. The methods section can report the HPC versions for the statistical analysis. The core edgeR QL-GLM algorithm did not change between these versions.

---

## Key citations

| Tool | Citation |
|------|----------|
| HiCCUPS | Rao SS et al. (2014) Cell 159(7):1665-1680 |
| Juicer | Durand NC et al. (2016) Cell Systems 3(1):95-98 |
| mariner | Kramer NE et al. (2022) mariner: Explore the Hi-Cs. Bioconductor |
| edgeR | Robinson MD, McCarthy DJ, Smyth GK (2010) Bioinformatics 26(1):139-140 |
| edgeR QL method | Lun ATL, Chen Y, Smyth GK (2016) Methods in Molecular Biology 1418:391-416 |
| Bowtie2 | Langmead B, Salzberg SL (2012) Nature Methods 9:357-359 |
| SAMtools | Li H et al. (2009) Bioinformatics 25(16):2078-2079 |
| Stripenn | Yoon S, Chandra R, Vahedi G (2022) Genome Biology 23:153 |
| cooler | Abdennur N, Mirny LA (2020) Bioinformatics 36(1):311-316 |
| TADCompare | Cresswell KG, Dozmorov MG (2020) Frontiers in Genetics 11:158 |
| HOMER | Heinz S et al. (2010) Molecular Cell 38(4):576-589 |
| ChIPseeker | Yu G, Wang LG, He QY (2015) Bioinformatics 31(14):2382-2383 |
| clusterProfiler | Yu G et al. (2012) OMICS 16(5):284-287 |
| fastp | Chen S et al. (2018) Bioinformatics 34(17):i884-i890 |
| hictk | Rossini R et al. (2024) Bioinformatics 40(7):btae408 |
| deepTools | Ramirez F et al. (2016) Nucleic Acids Research 44(W1):W160-W165 |
| DiffBind | Stark R, Brown G. DiffBind: Differential binding analysis of ChIP-Seq peak data |
| DESeq2 | Love MI, Huber W, Anders S (2014) Genome Biology 15:550 |
| EnhancedVolcano | Blighe K, Rana S, Lewis M (2018) EnhancedVolcano. Bioconductor |
