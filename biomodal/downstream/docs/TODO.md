# Biomodal Methylation Remaining Analyses TODO

Reference: `URS_Proposal.md`, `FIGURES.md`, `docs/urs/methylation-bio-revised-conclusions_2026-02-09.md`

---

## Status Legend

- `[x]` = Done (script exists and has been run)
- `[~]` = Partially done (infrastructure exists, analysis incomplete)
- `[ ]` = Not started

---

## Completed Work (Sections 1-22 of Visualization Pipeline)

The following are **done** and produce outputs in `plots/visualizations/`:

| Section | Analysis | Key Finding |
|---------|----------|-------------|
| 01 | QC & Sequencing Metrics | 336-476M mapped reads per sample |
| 02 | Sample Correlations | mC r=0.76-0.79, hmC r=0.48-0.51 |
| 03 | DMR Statistics by Region | Gene bodies primary (42% mC, 47% hmC); promoters <0.2% |
| 04 | Volcano Plots | mC: 8,836 sig; hmC: 9,930 sig |
| 05 | Coordinated Changes | 5,708/6,750 (84.6%) show mC up/hmC down |
| 06 | Top Genes | Syt1 top hit (+18% mC, -15% hmC) |
| 07 | Effect Size Distributions | Mean mC change +2.27%, hmC change -2.08% |
| 08 | GO/KEGG Enrichment | RNA splicing #1 GO term (248 genes, q=3.4e-48); delta-ratio decile GO/KEGG (08e-08h) |
| 09 | Summary Statistics | Integrated tables |
| 10 | Chromatin State Analysis | 49.9% DMRs at Active_Promoter; 94% of those hypermethylated |
| 11 | MeCP2 Integration | MeCP2 peak overlap with DMR direction; delta-ratio lm/glm (11f-11g) |
| 12 | ATAC-seq Correlation | 14.3% hypermethylated DMRs overlap ATAC-down (weak coupling) |
| 13 | ATAC + Chromatin + Loops | Loop-ATAC concordance 39.5%; Active_Enhancer-Active_Enhancer 80% |
| 14 | K119ub Peak Integration | K119ub-up at hypermethylated: OR=4.40 |
| 15 | hmC Cross-Mark Correlations | hmC direction O/E heatmaps for MeCP2/ATAC/K119ub |
| 16 | Raw Concordance Bars | Direct concordance rates per mark |
| 17 | Honest K119ub Assessment | Including "No Peaks" category; conditional direction |
| 18 | K119ub Bigwig Signal | Genome-wide increase confirmed (median log2FC +0.012, p<1e-24); weak gene-specific enrichment (Cliff's delta 0.089) |
| 19 | H3K27ac Peak Analysis | Status breakdown, waterfall, O/E, 4-mark comparison |
| 20 | Expression Integration | mC vs log2FC scatter, expression outcome bars |
| 21 | Discordant Gene Characterization | 4-quadrant analysis, composite panels |
| 22 | Demethylation Efficiency Ratio | 72.5% genes show decreased 5hmC/(5mC+5hmC); Cliff's delta=0.455 (medium); Active_Promoter most affected (med=-0.030) |
| 23 | Baseline 5hmC Predictor | 5hmC AUC=0.762 >> K119ub AUC=0.592; substrate availability confirmed (rho=-0.586) |
| 24 | DNMT3A Binding Prediction | TET impediment AUC=0.793 >> DNMT3A recruitment AUC=0.693 (DeLong p<2.2e-16); baseline 5hmC #1 predictor; K119ub negative direction argues against direct UDR recruitment |
| 25 | Delta-Ratio Linear Model Refits (14d) | Refits Section 23/24 logistic models as linear with delta_ratio response; feature importance rank correlation between logistic and linear frameworks |
| 26 | TET Triple-KO Comparison (14c) | GSE166423 BS/OxBS-seq; absolute attenuation 3.9%, relative 9.8%; QQ slope=0.106; rho=0.217 (57% baseline-driven, residualized=0.092); TET-KO binary (68.6% complete loss) vs BAP1-KO graded (47.9% moderate); 9 figures (26a-26i) + 2 tables |
| 27 | Methylation x Hi-C Loop Anchor Integration | Lost-loop anchors 2.5x more likely hypermethylated than gained (OR=2.54, p<2.2e-16); K119ub-gained x hyper at anchors OR=1.84; 113 triple-convergence genes; loop direction strongest predictor in logistic model (OR=2.42); GREAT-style mapping; 10 figures (27a-27e) + 7 tables |
| 28 | Coordinated Q4 Gene Characterization | Q4 (5,708) vs non-Q4 (1,042) across 9 dimensions; all comparisons p<1e-4; Q4 has larger effect sizes, lower expression, less ATAC-up, more K119ub gain, more MeCP2-up, less loop involvement; 1,438 GO terms enriched; 5 figures (28a-28e) + 1 table |
| 29 | A/B Compartment Methylation Mapping | mC hyper enriched in A (OR=14.71, p<2.2e-16); hmC loss enriched in A (OR=9.35); Spearman rho=0.348 (PC1 vs mC diff); B->A shifted bins 3.67x enriched for mC hypo (O/E); supports convergent TET-KO mechanism (DNMT3A redistribution); 7 figures (29a-29g) + 3 tables |
| 30 | Polycomb Target Gene Enrichment (17) | Polycomb targets DEPLETED from hyper (OR=0.064, p<2.2e-16) and ENRICHED in hypo (OR=8.71); Active_Promoter 65.2% hyper vs Repressed_Promoter 2.3%; 3 Polycomb definitions all consistent; magnitude also smaller at Polycomb; confirms dual-mechanism model; 6 figures (30a-30f) + 4 tables |

---

## 1. Methylation x Hi-C Loop Anchor Integration

**Source:** TODO Tier 1 #1, URS Proposal Methods

**Core question:** Are genes at differential loop anchors enriched for the coordinated mC up/hmC down pattern compared to non-anchor genes? Does methylation direction associate with loop direction (lost vs gained)?

### Tasks

- [x] **1a. Overlap differential loop anchors with DMR gene lists.** Take 2,910 differential loop anchors, extend to capture nearby gene promoters/bodies (+/-10kb), overlap with DMR gene lists. Fisher's exact test for enrichment of coordinated mC up/hmC down genes at loop anchors vs genome-wide. *(Done: section_27, panel 27a. GREAT-style regulatory domains (5kb up/1kb down/100kb max) mapped 2,897 anchor genes in DMR universe. Pooled anchors show OR=0.68 — anchor genes are LESS coordinated than background, because lost and gained loops have opposite methylation signatures that cancel. Both GREAT and nearest-gene methods agree.)*
- [x] **1b. Stratify by loop direction.** Do lost loops associate with hypermethylation and gained loops with hypomethylation, or vice versa? Split differential loops into lost (n=1,187) and gained (n=1,723), test whether anchor-proximal genes show directional methylation bias. *(Done: section_27, panel 27b. Lost-loop anchors 43.0% hypermethylated vs gained 22.9% vs background 31.6%. Lost vs Gained OR=2.54, p<2.2e-16. Delta demethylation ratio: Lost median=-0.015, Gained median~0, confirming greater TET impairment at lost-loop genes. THE main finding — direction-specific association.)*
- [x] **1c. Test K119ub-loop-methylation convergence.** The Hi-C analysis showed K119ub predicts loop loss (OR=10.7). Do the same loci also show hypermethylation? Overlap the K119ub-high loop anchors with hypermethylated DMRs to test "parallel consequences vs causal chain" (revised conclusions). *(Done: section_27, panel 27c. K119ub-gained at anchor x hypermethylated: OR=1.84, p=8.4e-7. O/E=1.40 for the convergent cell. 113 genes show triple convergence (K119ub gained + hypermethylated + lost loop). Supports causal chain model.)*
- [x] **1d. Logistic regression: methylation direction ~ loop direction + distance + chromatin state.** Quantitative model testing whether loop loss/gain at an anchor predicts methylation change direction at the associated gene, controlling for confounders. *(Done: section_27, panel 27d. Loop direction is dominant predictor: Lost OR=2.42, p=4.2e-35. Distance features NS. Active_Enhancer anchors OR=1.58 (vulnerable), Repressed_Promoter OR=0.17, Polycomb OR=0.21 (protected). Linear model R²=0.149. Confirms direction effect survives multivariate adjustment.)*
- [x] **1e. Shared anchor methylation profile.** At the 212 shared anchor hubs (where loops switch partners), do the associated genes show the coordinated pattern at higher rates than non-shared anchors? *(Done: section_27, panel 27e. Opposite to prediction — shared anchors show LOWER coordinated rates (35.2%) than non-shared (45.7%) or background (54.6%), OR=0.45. Delta-ratio also less negative (median=-0.004 vs -0.011). Shared anchors may be structural hubs resistant to methylation perturbation.)*

### Existing resources

- **Loop data (Late):** `../../25042-late_outputs/merged_loops/characterized_loops.tsv` (2,910 loops, 57 cols: coordinates, logFC, FDR, direction, distance, 5-mark ChIP overlaps, 7-category anchor_type, gene annotations)
- **Loop data (Early):** `../../250831-early_outputs/merged_loops/characterized_loops.tsv` (165 loops)
- **DMR data:** Loaded via `_shared_config.R` (`mc_dmr`, `hmc_dmr`)
- **Coordinated gene list:** `plots/visualizations/tables/` (from Section 05)
- **Shared anchor data:** `../../output/shared_anchor_analysis/late/` (212 anchors)
- **K119ub-loop integration:** `../../scripts/h2ak119ub_loop_integration.R` (logistic regression, OR=10.7)
- **Gene association pattern:** `../../tads/scripts/deg_tad_violin.R` (GREAT-style regulatory domains, reusable)

### Notes

- This is the most direct test of the URS proposal's central promise: connecting methylation to 3D architecture.
- The GREAT-style approach (5kb upstream, 1kb downstream, 100kb max extension) is already implemented for TAD boundaries and can be adapted.
- Consider both gene-level (nearest gene to anchor) and region-level (DMR overlapping anchor +/-10kb) approaches.

---

## 2. A/B Compartment Methylation Mapping

**Source:** TODO Tier 1 #2, Lopez-Moyado et al. (2019) framework

**Core question:** Does hypermethylation concentrate in A compartment (euchromatin) while hypomethylation falls in B compartment? This is the signature of DNMT3A redistribution seen in TET-KO.

### Tasks

- [x] **2a. Assign DMR genes to A/B compartments.** Use PC1 eigenvector from Hi-C HOMER analysis (late timepoint: `tads/tad-pc-analysis/inputs/late/diffPC/all_PC1.txt`) to classify genomic bins as A (positive PC1) or B (negative PC1). Map each DMR gene to its compartment. *(Done: section_29. Used `diffcompartments.txt` (104,071 bins, LATE timepoint). Control mean PC1 for A/B classification. Deduplicated by closest-to-TSS bin per gene. 22,143 unique genes -> 16,086 matched to DMRs (76.8% match rate). A: 12,347 genes (76.8%), B: 3,739 genes (23.2%).)*
- [x] **2b. Test compartment enrichment of methylation direction.** Fisher's exact: are hypermethylated genes enriched in A compartment? Are hypomethylated genes enriched in B compartment? This tests the Lopez-Moyado TET-KO phenotype parallel. *(Done: section_29, Step 2. mC hyper -> A enriched: OR=14.71 (12.65-17.20), p<2.2e-16 (43.8% in A vs 5.0% in B). mC hypo -> B enriched: OR=1.67, p<2.2e-16. hmC loss -> A enriched: OR=9.35, p<2.2e-16. All four directional tests significant — strong TET-KO parallel confirmed.)*
- [x] **2c. Differential compartment x DMR overlap.** The Hi-C analysis found 44% of the genome shows significant compartment shifts (44,703 regions). Do genes in B-to-A shifted regions show hypermethylation (gained euchromatin + DNMT3A access)? *(Done: section_29, Step 3. 629 genes at shifted compartments (427 B->A, 202 A->B). B->A bins show 3.67x enrichment for mC hypomethylation (O/E), opposite to static prediction — these bins are actively losing heterochromatic character. A->B x hyper O/E=1.17. Shift Fisher's tests also computed.)*
- [x] **2d. Visualization: compartment-stratified methylation effect sizes.** Violin or box plots of mC/hmC change magnitude in A vs B compartment genes, and in stable vs shifted compartment regions. *(Done: section_29, Step 4. 7 figures: 29a/29b violins by A/B, 29c/29d violins by shift, 29e stacked bar DMR direction proportions, 29f PC1 vs mC scatter (Spearman rho=0.348, p<2.2e-16), 29g composite panel. Wilcoxon A vs B: mC p<2.2e-16 (median A: +0.012, B: -0.004), hmC p<2.2e-16.)*

### Existing resources

- **PC1/compartment data (LATE timepoint, bwt2_bam deep-seq, 104,071 regions):**
  - `../../tads/tad-pc-analysis/inputs/late/diffPC/all_PC1.txt` — Raw PC1 eigenvector, 6 samples (3 ctrl, 3 mut), gene annotations. Positive PC1 = A compartment, negative = B compartment. 25x50kb resolution.
  - `../../tads/tad-pc-analysis/inputs/late/diffPC/diffcompartments.txt` — Differential PC1 from HOMER `getDiffExpression.pl` with Difference, p-value, adj. p-value columns. Positive Difference = B→A shift (more active in mutant).
  - `../../tads/tad-pc-analysis/inputs/late/diffPC/regions.Up_mut_vs_ctrl.txt` — B→A shifted regions (full annotations)
  - `../../tads/tad-pc-analysis/inputs/late/diffPC/regions.Down_mut_vs_ctrl.txt` — A→B shifted regions (full annotations)
  - `../../tads/tad-pc-analysis/inputs/late/diffPC/corrdiff/` — Per-replicate compartment correlation differences (ctrl_M{1,2,3}_vs_mut_M{1,2,3})
- **Per-sample PC1 bedGraphs:** `../../tads/tad-pc-analysis/inputs/late/tags/{ctrl,mut}_M{1,2,3}/tagdirs.25x50kb.PC1.{bedGraph,txt}` (12 files)
- **Differential TAD data (LATE):**
  - `../../tads/tad-pc-analysis/inputs/late/diffTAD/Bap1.diff.tad.txt` — Differential TAD inclusion ratios
  - `../../tads/tad-pc-analysis/inputs/late/diffTAD/BAP1.tad.scores.txt` — TAD scores per sample
  - `../../tads/tad-pc-analysis/inputs/late/diffTAD/merged.tad.2D.bed` — TAD boundaries
- **Differential compaction (LATE):** `../../tads/tad-pc-analysis/inputs/late/diffcompaction/` — Per-replicate DLR, ICF, PC1 bedGraphs + region-specific (CTCF, H3K27ac, H3K27me3, TSS) histogram/table files
- **PC-associated regions:** `../../tads/tad-pc-analysis/inputs/late/PC_regions/` — Per-replicate compartment-mark association tables (H3K27ac, H3K27me3, Superenhancers, TSS)
- **Early timepoint (250831, for comparison):** `../../tads/tad-pc-analysis/inputs/new/` — Same file types, 101,684 regions
- **Differential compartment analysis script:** `../../scripts/compartment_volcano_plot.R` (ran on early data: 44,703 significant, 6,462 B->A, 4,184 A->B — re-run needed on late data)
- **Lopez-Moyado framework:** Hypermethylation in A compartment + hypomethylation in B compartment = DNMT3A redistribution signature (docs/urs/methylation-bio-revised-conclusions)
- **HPC source:** `/expanse/lustre/projects/csd940/ctea/homer/bwt2_bam/` (late timepoint, bowtie2-aligned deep-seq BAMs → HOMER tag dirs → PC1/TAD analysis by ctea)

### Notes

- The early timepoint PC1 data (250831_bams, in `inputs/new/`) was previously used for compartment volcano plots. The late timepoint data (bwt2_bam, in `inputs/late/`) has ~2.4K more genomic regions (104K vs 101K) due to deeper sequencing.
- The compartment volcano plot script needs to be re-run on the late timepoint `diffcompartments.txt` to get updated significance counts.
- If the A/B pattern matches TET-KO, it supports the "convergent mechanisms" interpretation from the revised conclusions.
- If the pattern does NOT match, it suggests BAP1-mediated mechanism is distinct from direct TET loss, which is also informative.

---

## 3. Baseline 5hmC as Predictor of DMR Susceptibility

**Source:** TODO Tier 1 #3, Revised conclusions hypothesis

**Core question:** Does baseline (wildtype) 5hmC level predict which genes become DMRs? This tests the "substrate availability" interpretation: genes with more 5hmC have more substrate for TET impediment to affect.

### Tasks

- [x] **3a. Extract per-gene WT 5hmC levels.** From the Zarr stores or modality output, obtain mean 5hmC fraction per gene body in control samples. *(Done: section_23, averaged ctrl-M + ctrl-F regional-frac, 20,947 genes)*
- [x] **3b. Logistic regression: DMR status ~ WT 5hmC level.** Binary outcome (significant hmC DMR vs not), predictor is baseline 5hmC. If baseline 5hmC is a strong predictor, supports substrate availability model. *(Done: Model A AUC=0.762, OR=663K, p<2.2e-16)*
- [x] **3c. Compare predictive power: WT 5hmC vs K119ub.** Run competing models: (a) DMR ~ WT 5hmC, (b) DMR ~ K119ub signal, (c) DMR ~ both. AIC comparison. Section 18 already showed K119ub is a weak gene-specific predictor (Cliff's delta 0.089); if 5hmC is stronger, it confirms the revised conclusions. *(Done: 5hmC AUC=0.762 >> K119ub AUC=0.592; combined AUC=0.800; R² 6x larger for 5hmC)*
- [x] **3d. Dose-response visualization.** Scatter plot of WT 5hmC level (x) vs magnitude of hmC change (y), colored by DMR significance. Expect negative correlation (high baseline = more loss). *(Done: Spearman rho=-0.586, p<2.2e-16, confirms substrate availability)*

### Existing resources

- **Zarr stores:** `/expanse/lustre/projects/csd940/zalibhai/biomodal/evoC-run/output/duet-1.5.0_evoC_Bap1_run_6bp/sample_outputs/zarr_store/` (CG context)
- **Modality feature extraction:** `modality/outputs/run-3/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/Extract_*/` (count, mean, regional-frac per gene)
- **K119ub signal:** Bigwig-derived per-gene signal (from Section 18 preprocessing)

### Notes

- The modality `get mean` output may already contain per-sample mean 5hmC at gene bodies. Check `Extract_*/hmC_mean_*.tsv` files before extracting from Zarr.
- Cerebellar neurons have ~40% of modified cytosines as 5hmC (Kriaucionis & Heintz 2009), so there is a wide range of baseline levels to test.

---

## 4. MeCP2 Functional Readout

**Source:** TODO Tier 1 #4, URS Proposal Methods

**Core question:** Does the coordinated mC up/hmC down pattern at gene bodies predict increased MeCP2 binding? MeCP2 binds 5mC with higher affinity than 5hmC (Mellen et al. 2012, 2017), so 5hmC-to-5mC conversion should increase MeCP2 occupancy.

### Tasks

- [~] **4a. Basic MeCP2-DMR overlap.** Section 11 already computes MeCP2 peak overlap with DMR direction. -> `plots/visualizations/11a_mecp2_overlap/` through `11e_mecp2_integration_heatmap/`. **Done but may need extension.**
- [~] **4b. Predict MeCP2 gain from coordinated mC up/hmC down.** Among coordinated genes (n=5,708), test whether those with MeCP2 gain (from CUT&RUN differential) are enriched compared to non-coordinated genes. Fisher's exact test. *Partially done: Section 11d Wilcoxon test (p=8.64e-08) shows coordinated genes have significantly different MeCP2 fold change vs all other genes. Fisher's exact on binary gain not yet done.*
- [~] **4c. Quantitative MeCP2-methylation model.** Scatter plot: delta_mC (x) vs delta_MeCP2 (y), per gene. Expect positive correlation (more mC gain = more MeCP2 gain). Stratify by coordinated vs discordant genes. *Partially done: Section 11c scatter (rho=0.015, p=0.18) + 11f-11g delta-ratio regression models (lm coeff=-0.319, p=6.5e-05; glm OR<1, p=5.2e-04). Stratification by coordinated vs discordant not yet done.*
- [ ] **4d. MeCP2 at loop anchors.** Cross-reference MeCP2 differential binding with loop anchor positions. Do anchors where MeCP2 increases also show loop changes? This connects the methylation reader to 3D architecture.

### Existing resources

- **Section 11 outputs:** `plots/visualizations/11{a-e}_*/` (overlap, fold change, scatter, heatmap)
- **MeCP2 peak files:** Loaded via `_shared_config.R` (`MECP2_FILES`: annotated, up, down)
- **Coordinated gene list:** From Section 05 tables

### Notes

- Section 11 provides the foundation but tests only overlap/direction. Tasks 4b-4d add the functional prediction (if 5hmC becomes 5mC, MeCP2 should bind more).
- This connects the chemical modification change to a downstream reader protein, strengthening the functional interpretation.

---

## 5. Transcription Factor Motif Enrichment at DMRs

**Source:** TODO Tier 1 #5

**Core question:** What TF binding sites are enriched in hypermethylated vs hypomethylated DMR regions? If DNMT3A-binding motifs or TET-cofactor motifs appear, that is mechanistically informative.

### Tasks

- [ ] **5a. Generate DMR BED files for HOMER input.** Extract hypermethylated DMR coordinates and hypomethylated DMR coordinates as separate BED files (gene body regions).
- [ ] **5b. Run HOMER findMotifsGenome.pl.** On hypermethylated DMRs (foreground) vs genome background. Repeat for hypomethylated DMRs.
- [ ] **5c. Compare motif enrichment between directions.** Which motifs are specific to hypermethylated vs hypomethylated regions? Are neuron-specific TFs (e.g., NeuroD, MEF2, CREB) enriched?
- [ ] **5d. Test DNMT3A/TET-associated motifs specifically.** Look for known DNMT3A sequence preferences and TET-associated cofactor binding sites (e.g., IDAX/CXXC motifs) in the results.

### Existing resources

- **DMR data:** `mc_dmr`, `hmc_dmr` from `_shared_config.R` (contain genomic coordinates)
- **HOMER installation:** Available on Expanse HPC
- **CTCF motif reference:** `../../peaks/ctcf_motifs_mm10.bed` (114,081 motifs, already used in Hi-C pipeline)
- **Hi-C HOMER motif analysis:** `../../scripts/` (pattern exists for HOMER motif enrichment at enhancer subsets)

### Notes

- Quick to run on HPC (~30 min per set). Low effort, potentially high insight.
- Consider running on the coordinated (mC up/hmC down) gene bodies vs the discordant (mC down/hmC up) gene bodies as foreground/background pair.

---

## 6. Single-CpG Resolution Gene Body Analysis

**Source:** TODO Tier 2 #6

**Core question:** Within gene bodies, are mC/hmC changes concentrated at specific sub-features (exon-intron boundaries, splice sites, specific exons)? Exonic 5hmC has different regulatory roles than intronic 5hmC (Szulwach et al. 2011).

### Tasks

- [ ] **6a. Extract sub-gene-body methylation from Zarr stores.** Use biomodal DUET 6bp-resolution data to get CpG-level methylation across gene bodies for top affected genes (e.g., Syt1, Tmem238, Mcu).
- [ ] **6b. Annotate CpGs by gene feature.** Classify each CpG as exonic, intronic, exon-intron boundary (+/-50bp of splice site), 5'UTR, 3'UTR using GENCODE vM25 annotations.
- [ ] **6c. Compare methylation change by gene feature.** Are mC up/hmC down changes concentrated at exons specifically? Box plots of delta_mC and delta_hmC stratified by feature type.
- [ ] **6d. Metagene profile.** Aggregate CpG-level changes across all significant gene bodies, scaled to gene length, showing positional bias (5' vs 3', exons vs introns).

### Existing resources

- **Zarr stores:** CG context, 6bp resolution on Expanse
- **GENCODE annotations:** `modality/mm10/gencode.vM25.mouse.genes.annotation.bed.gz`
- **Top genes:** From Section 06 tables

### Notes

- This requires working with the raw Zarr stores on HPC, not the gene-level DMR summaries.
- Consider the modality `get regional-frac` output which may provide sub-gene-body resolution.
- If exonic changes dominate, this strengthens the functional interpretation (exonic 5hmC affects transcription elongation and splicing).

---

## 7. RNA Splicing Analysis

**Source:** TODO Tier 2 #7

**Core question:** Do genes with gene body methylation changes also show altered splicing? Gene body 5mC regulates splicing kinetics (slower elongation -> more exon inclusion), so mC gain could directly affect splicing.

### Tasks

- [ ] **7a. Determine RNA-seq BAM availability.** Check if aligned BAM files from the adult timepoint RNA-seq are accessible (needed for rMATS/SUPPA2). Files may be on Expanse or need to be requested from the lab.
- [ ] **7b. Run differential splicing analysis.** rMATS or SUPPA2 on BAP1-KO vs WT RNA-seq BAMs. Identify differentially spliced events (skipped exons, alternative 5'/3' splice sites, retained introns, mutually exclusive exons).
- [ ] **7c. Overlap spliced genes with DMR genes.** Fisher's exact: are genes with mC-up DMRs enriched for differential splicing events? Do the directions make biological sense (mC gain -> more exon inclusion)?
- [ ] **7d. Correlate splicing magnitude with methylation magnitude.** For genes with both DMRs and splicing changes, scatter plot of delta_mC vs delta_PSI (percent spliced in).

### Existing resources

- **RNA-seq DEG results:** `../../tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx` (gene-level, already used in Section 20)
- **GO enrichment:** RNA splicing is #1 GO term (248 genes, q=3.4e-48) from Section 08

### Notes

- RNA splicing being the top GO term (Section 08) makes this a high-priority follow-up despite being Tier 2.
- Requires BAM files, not just DEG tables. Check data availability before planning.
- rMATS needs aligned BAMs + GTF annotation. SUPPA2 can work with transcript quantifications (Salmon/kallisto).

---

## 8. Expression-Methylation Quantitative Model

**Source:** TODO Tier 2 #8

**Core question:** Does the magnitude of methylation change predict expression change? Section 20 shows basic overlap but not quantitative modeling.

### Tasks

- [~] **8a. Basic scatter: delta_mC vs delta_expression per gene.** Section 20 produces `20b_methylation_vs_expression_scatter/`. **Partially done.** Needs review for stratification and quantitative fit.
- [ ] **8b. Stratify by chromatin state.** Separate scatter plots for Active_Promoter, Polycomb, Active_Enhancer, etc. The current Section 20 may not stratify by chromatin state.
- [ ] **8c. Multivariate model: delta_expression ~ delta_mC + delta_hmC + chromatin_state + K119ub.** Linear model quantifying the relative contribution of each epigenetic change to expression change. Report R-squared and coefficient significance.
- [ ] **8d. Concordance by methylation magnitude.** Bin genes by magnitude of mC change (tertiles or quartiles). Show that genes with larger mC gain show stronger expression downregulation. This tests the dose-response relationship.

### Existing resources

- **Section 20 outputs:** `plots/visualizations/20{a-d}_*/` (expression breakdown, scatter, violin, heatmap)
- **RNA-seq data:** `../../tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx`
- **Chromatin state classification:** Via `_shared_config.R` `classify_chromatin_state()`
- **K119ub bigwig signal:** From Section 18 preprocessing

### Notes

- Section 20 exists but is primarily visualization. The quantitative model (8c) and dose-response (8d) are the new contributions.
- ATAC-methylation correlation is weak (rho=-0.076 from Section 12), but expression may be more tightly coupled.

---

## 9. ATAC-seq Follow-Up Analyses

**Source:** Hi-C TODO Section 9, biomodal Section 13

**Core question:** Section 13 established baseline ATAC-chromatin-loop relationships. These follow-ups test specific mechanistic predictions.

### Tasks

- [ ] **9a. Anchor-level (unpooled) ATAC overlap rates.** Show ATAC overlap separately for anchor1 vs anchor2 instead of pooling. Confirms whether the concordance signal is symmetric or driven by one anchor position.
- [ ] **9b. Polycomb ATAC-up at hypermethylated DMR overlap.** Section 13a showed ATAC-up peaks are unexpectedly enriched at Polycomb regions (12.7% vs 0.5%). Section 12 showed hypomethylated DMRs overlap ATAC-up (50.5%). Overlap the Polycomb-classified ATAC-up peaks with hypermethylated DMRs to test whether Polycomb accessibility gain is linked to methylation changes.
- [ ] **9c. Permutation test for loop anchor ATAC enrichment.** Current ATAC overlap rates at loop anchors (32.6% ATAC-up, 14.0% ATAC-down) lack a null comparison. Shuffle ATAC peaks across the genome (regioneR) and compare observed vs expected overlap at loop anchors.

### Existing resources

- **Section 13 script:** `scripts/viz_sections/section_13_atac_chromatin_and_loops.R`
- **Section 13 outputs:** `plots/visualizations/13{a-f}_*/`, tables in `tables/atac_chromatin_*.tsv`, `tables/loop_anchor_*.tsv`
- **ATAC peaks:** `../../peaks/atac_seq/ATAC_up.bed` (7,620), `../../peaks/atac_seq/ATAC_down.bed` (3,744)
- **Loop annotations:** `../../peaks/loop_annotation_extended/late/extended_characterized_loops.tsv`
- **Permutation framework:** regioneR/regioneReloaded (Bioconductor)

---

## 10. Expanded Cohort Analyses (Summer)

**Source:** TODO Tier 3 #10-12, URS Proposal

**Core question:** With n=4 per condition (properly sex-balanced), which findings are robust and which were driven by sex confounding?

### Tasks

- [ ] **10a. Process new samples through DUET pipeline.** Two additional biological replicates per condition through biomodal DUET v1.5.0, doubled read depth. `upstream/evoc_run.sb`
- [ ] **10b. Run modality XPLR with batch covariate.** DMR calling with `Covariates=batch` in config. Compare to run-2 (no sex covariate) and run-1 (sex covariate, no results).
- [ ] **10c. Sex-stratified analysis.** With proper balance, run full analysis with sex covariate AND sex-stratified subgroup analyses. Quantify how much of the current signal is robust vs sex-driven.
- [ ] **10d. Batch effect characterization (PCA).** PCA on raw methylation matrix before and after batch correction. If genotype dominates PC1 and batch is lower PC, proceed normally. If batch is PC1, apply ComBat or similar.
- [ ] **10e. Power analysis update.** With n=4 vs n=2, report change in number of significant DMRs, effect size distributions, and false discovery proportion estimates.

### Existing resources

- **Current pipeline:** `upstream/evoc_run.sb` (DUET pipeline), `downstream/modality/run_modality.sb` (DMR calling)
- **Metadata template:** `downstream/modality/metadata.tsv` (needs expansion for new samples)
- **Config template:** `downstream/modality/config_CG.txt` (add batch covariate)

### Notes

- This is gated on new samples arriving from the lab and sequencing turnaround.
- The URS proposal timeline: Weeks 1-3 pipeline refinement + upstream processing, Weeks 4-6 differential analysis with full cohort.
- Current limitation: n=2 with sex=genotype confound. Everything in Sections 1-21 needs re-evaluation once the expanded cohort is processed.

---

## 11. DNMT3A Binding Prediction

**Source:** TODO Tier 3 #12, URS Proposal dual-mechanism model

**Core question:** Can we computationally predict which gene bodies become hypermethylated based on the DNMT3A-UDR recruitment model? This tests the dual-mechanism hypothesis without ChIP-seq.

### Tasks

- [x] **11a. Build feature matrix per gene.** For each gene body: H2AK119ub signal (from bigwig), ATAC accessibility, CpG density, baseline 5mC, baseline 5hmC, gene length, expression level. *(Done: section_24, 11,936 genes with all 7 features after inner join cascade; 46.2% hyper-DMR)*
- [x] **11b. Train logistic regression: hypermethylated DMR ~ features.** Binary outcome (significant mC-up DMR vs not), predictors from 11a. Report feature importance and model performance (AUC). *(Done: 5 models — Full AUC=0.857, DNMT3A recruitment AUC=0.693, TET impediment AUC=0.793, K119ub only AUC=0.645, Stepwise AUC=0.857. Random forest also fitted (OOB error 20.8%))*
- [x] **11c. Test DNMT3A-UDR predictions.** If the model works: H2AK119ub + accessibility should be the strongest predictors (DNMT3A recruited by ubiquitin to accessible chromatin). If H2AK119ub alone is sufficient, supports direct recruitment. If accessibility alone suffices, suggests DNMT3A access rather than recruitment. *(Done: K119ub is a NEGATIVE predictor (beta=-1.035), meaning genes with MORE K119ub are LESS likely to be hypermethylated. Baseline 5hmC is the strongest positive predictor (beta=+1.267). RF confirms: baseline_hmc #1, k119ub #2. This argues against direct DNMT3A-UDR recruitment as primary mechanism.)*
- [x] **11d. Cross-validate with TET impediment model.** Compare predictive power of "DNMT3A recruitment features" (K119ub, accessibility, CpG density) vs "TET impediment features" (baseline 5hmC, chromatin compaction) for hypermethylation. *(Done: DeLong test p<2.2e-16 — TET impediment (AUC=0.793) significantly outperforms DNMT3A recruitment (AUC=0.693). Supports TET blockade as primary mechanism.)*

### Existing resources

- **K119ub bigwig signal:** From Section 18 preprocessing
- **ATAC consensus peaks:** `../../peaks/atac_seq/` (union ~71,000 peaks)
- **CpG annotations:** GENCODE vM25 in `modality/mm10/`
- **Baseline methylation:** From modality feature extraction output

### Notes

- This is feasible even with n=2 because it is a cross-sectional feature prediction (not differential testing).
- The Chen et al. (2024) cryo-EM structures showing DNMT3A UDR bidentate interaction with H2AK119ub and nucleosome acidic patch provide the structural basis.
- Consider using a random forest or gradient boosting for feature importance ranking in addition to logistic regression.

### Potential Enhancements

- [x] Add k-fold cross-validation for more robust AUC estimates (current AUCs are in-sample) *(Done: 10-fold stratified CV shows near-zero optimism across all models — Full CV=0.856+/-0.011, TET=0.793+/-0.014, DNMT3A=0.692+/-0.017. Model ranking robust under CV.)*
- [x] Stratified analysis excluding Active_Promoter genes (which dominate at 66.9% hyper-DMR rate) to test whether predictions hold in non-promoter chromatin contexts *(Done: Non-promoter subset (N=6,596, 29.5% hyper-DMR) — TET still outperforms DNMT3A (p=6.4e-15). K119ub beta becomes even more negative (-1.509 vs -1.035), confirming it is NOT driven by Active_Promoter confounding.)*
- [x] Interaction terms (K119ub x baseline_hmc) to test whether K119ub modulates the TET impediment effect — the negative K119ub direction may reverse in high-5hmC genes if dual mechanism applies *(Done: Interaction NS (p=0.425, beta=-0.032). K119ub is negative across ALL 5hmC tertiles (OR=0.25-0.46). No evidence for dual mechanism — TET impediment is the dominant explanation.)*

---

## 12. Developmental Timepoint Comparison

**Source:** TODO Tier 2 #9

**Core question:** If early-timepoint (P12) methylation data exists or is generated, does methylation show the same progressive pattern as Hi-C loops (165 at P12 -> 2,910 at P60, 18x amplification)?

### Tasks

- [ ] **12a. Determine P12 methylation data availability.** Check whether the lab has or plans to generate early-timepoint DUET data from the same BAP1-KO model.
- [ ] **12b. (If available) Run parallel DMR analysis at P12.** Same pipeline, same parameters. Compare number of DMRs, effect sizes, direction bias.
- [ ] **12c. Temporal trajectory comparison.** If methylation changes are already present at P12, they may be upstream of loop changes. If they appear later, loop architecture may drive methylation. This is the key causal ordering question.

### Existing resources

- **Hi-C temporal data:** 165 differential loops at P12 (57% DOWN) vs 2,910 at P60 (59% UP) -- 18x amplification with direction reversal
- **Early Hi-C loop data:** `../../250831-early_outputs/merged_loops/characterized_loops.tsv`
- **Early RNA-seq:** `../../tads/young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx`

### Notes

- Currently no P12 methylation data exists. This is contingent on lab generating it.
- If no P12 methylation, consider using the early-timepoint RNA-seq + Hi-C to at least triangulate the temporal model.

---

## 13. Non-CG Methylation Investigation (mCA Focus)

**Source:** URS Proposal Methods (Weeks 1-3), FURTHER_ANALYSIS #9

**Core question:** Can non-CG methylation (CHG/CHH) be detected with deeper sequencing or bisulfite-based validation? Current data shows <1% methylation at CHG/CHH (below reliable detection). The most important non-CG context is **mCA** (methylation of CA dinucleotides), which in brain accumulates specifically in postmitotic neurons (Lister et al. 2013), is written almost exclusively by DNMT3A (not DNMT3B), is bound by MeCP2 with higher affinity than mCG in some contexts, and is associated with gene repression (opposite to gene body mCG).

### Tasks

- [ ] **13a. Evaluate CHG/CHH DMR calling from run-3 (deep-seq).** The deeper sequencing may improve power. Check `modality/outputs/run-3/outputs_CHG/` and `outputs_CHH/` if they exist.
- [ ] **13b. Compare DUET evoC vs bisulfite for non-CG detection.** Literature comparison of detection limits. Biomodal's enzymatic conversion may have different sensitivity than bisulfite at low methylation levels.
- [ ] **13c. (If detectable) Characterize non-CG methylation changes.** Non-CG methylation in brain is associated with MeCP2 binding (Lister et al. 2013) and neuronal maturation. If BAP1-KO affects it, this adds another layer.
- [ ] **13d. Test DNMT3A-specific mCA prediction.** If DNMT3A is being redistributed via UDR-H2AK119ub, mCA should change at the same loci that show mCG gain. This is a specific prediction of the DNMT3A recruitment arm that is distinct from the TET impediment arm (TET primarily acts on CG context). If mCA and mCG changes are correlated at the same gene bodies, that is evidence for DNMT3A recruitment specifically. If mCG changes but mCA does not, that favors TET impediment alone.
- [ ] **13e. mCA x MeCP2 cross-reference.** MeCP2 binds mCA with high affinity in neurons. If mCA changes are detectable, test whether affected loci overlap with MeCP2 differential binding (from Section 11). This adds a neuron-specific functional readout.

### Existing resources

- **CHG/CHH Zarr stores:** Generated by DUET pipeline with `--chg-chh-contexts` flag
- **CHG config:** `downstream/modality/config_CHG.txt`
- **Previous CHG results:** Run-2 showed no significant DMRs at CHG context
- **MeCP2 data:** Via `_shared_config.R` MECP2_FILES

### Notes

- Even at <1% absolute levels, you need detectable *differences*, not high absolute levels. With ~400M reads per sample, there should be reasonable coverage at CA sites in gene bodies.
- The DNMT3A-specific prediction (13d) makes this higher value than a simple "can we detect CHG" question.
- If no signal, a negative result is still informative and worth documenting.

---

## 14. 5hmC/(5mC+5hmC) Demethylation Efficiency Ratio

**Source:** FURTHER_ANALYSIS #1

**Core question:** Rather than treating mC and hmC as two separate signals, compute the ratio 5hmC/(5mC+5hmC) per gene in WT vs KO. This ratio directly represents TET conversion efficiency -- what fraction of modified cytosines at a locus have been successfully oxidized. A decrease in this ratio in KO is the most direct per-gene measure of impaired active demethylation. It collapses two marks into one biologically interpretable number and sidesteps the issue of mC and hmC having different dynamic ranges and variance structures.

### Tasks

- [x] **14a. Compute per-gene demethylation ratio in WT and KO.** For each gene body, calculate 5hmC/(5mC+5hmC) in control and mutant samples. Use modality feature extraction output (mean mC and hmC per gene) or Zarr stores. *Done: Section 22 computes ratio_ctrl and ratio_mut for 20,898 genes. WT median=0.1284, KO median=0.1182. Per-sample ratios also computed from feature extraction TSVs (22h).*
- [x] **14b. Compute delta-ratio per gene.** KO ratio minus WT ratio. This continuous "demethylation activity score" can replace or supplement the binary coordinated/discordant classification in downstream models (Sections 1, 3, 8, 11). *Done: Section 22 exports `demethylation_ratio_all_genes.tsv` (20,898 genes with delta_ratio). 72.5% genes decreased. Cliff's delta=0.455 (medium). Figures 22a-22h produced.*
- [x] **14c. Compare WT ratio distribution to published TET-KO data.** If the BAP1-KO ratio shift phenocopies direct TET loss (e.g., Rao lab data, Lopez-Moyado et al. 2019), that is a strong claim for convergent mechanisms. If the shift pattern differs, BAP1 is working through a distinct pathway. *Done: Section 26 implements formal quantitative comparison using GSE166423 TET triple-KO BS/OxBS-seq bigWigs. Pipeline: download_tet_ko_data.sb -> preprocess_tet_ko_bigwig.R -> section_26_tet_ko_comparison.R. Key results: absolute attenuation 3.9% (QQ slope=0.106), relative attenuation 9.8% after baseline normalization; per-gene Spearman rho=0.217 but variance decomposition shows 57% was baseline-driven (residualized rho=0.092); TET-KO shows binary response (68.6% complete loss, 27.0% no signal) vs BAP1-KO graded (1.4% strong, 47.9% moderate, 50.7% weak). Supports indirect TET blockade model: BAP1-KO impairs TET access (graded dimmer) rather than eliminating TET activity (binary switch). 9 figures (26a-26i) + 2 tables (37-row summary TSV).*
- [x] **14d. Use delta-ratio as primary response variable.** Refit key models (baseline 5hmC predictor in Section 3, expression-methylation model in Section 8, DNMT3A prediction in Section 11) using the ratio instead of separate mC/hmC metrics. *Done: Section 25 refits all key models from Sections 23 and 24 using delta_ratio as continuous linear response. Section 23 refits: 3 linear models (baseline 5hmC, K119ub, combined) with R² vs AUC comparison. Section 24 refits: 4 linear models (Full, DNMT3A recruitment, TET impediment, K119ub only) with feature importance rank correlation between logistic and linear betas. Includes residual diagnostics, predicted vs observed scatter. 6 figures (25a-25f) + 4 tables. Previously done: Section 8 decile GO/KEGG (08e-08h), Section 11 MeCP2 regression (11f-11g).*

### Existing resources

- **Modality feature extraction:** `modality/outputs/run-3/outputs_CG/.../Extract_*/` (per-gene mean mC and hmC)
- **Zarr stores:** 6bp resolution on Expanse
- **TET-KO reference:** Lopez-Moyado et al. (2019), Rao lab published datasets

### Notes

- High priority despite being a new suggestion -- should be computed early and used as a primary variable going forward.
- The ratio is more robust than raw mC/hmC values because it normalizes for total modification level per gene.

---

## 15. CTCF Binding Site Methylation at Loop Anchors

**Source:** FURTHER_ANALYSIS #3

**Core question:** CpG methylation at CTCF motifs blocks CTCF binding, and CTCF is the primary insulator/loop anchor protein. Do CTCF motif sites at differential loop anchors show methylation changes? If CTCF sites at lost loop anchors gain methylation, that is a direct causal link between the methylation and Hi-C arms (methylation blocks CTCF, loop lost). If CTCF site methylation is unchanged at loop anchors, that rules out this mechanism and supports the "parallel consequences" interpretation.

### Tasks

- [ ] **15a. Extract methylation at CTCF motif CpGs.** From the Zarr stores, pull mC and hmC levels at CpGs within CTCF motif sites (114,081 genome-wide from `ctcf_motifs_mm10.bed`). Compute WT and KO levels and delta.
- [ ] **15b. Stratify by loop anchor status.** Compare CTCF-motif methylation at: (a) differential loop anchors (lost loops), (b) differential loop anchors (gained loops), (c) stable loop anchors, (d) non-anchor CTCF sites (background).
- [ ] **15c. Test causal vs parallel model.** If CTCF sites at lost loop anchors gain methylation, direct causal link. If methylation is unchanged at CTCF sites at loop anchors, rules out this mechanism and supports "parallel consequences" from revised conclusions.
- [ ] **15d. CTCF ChIP-seq signal at methylated CTCF sites.** Using CTCF ChIP data (32,487 peaks in `peaks/CTCF.bed`), test whether CTCF sites that gain methylation show reduced ChIP signal.

### Existing resources

- **CTCF motifs:** `../../peaks/ctcf_motifs_mm10.bed` (114,081 genome-wide)
- **CTCF ChIP peaks:** `../../peaks/CTCF.bed` (32,487 peaks, Late)
- **Loop anchors:** `../../25042-late_outputs/merged_loops/characterized_loops.tsv`
- **Zarr stores:** 6bp resolution on Expanse

### Notes

- This is a glaring missing link between the methylation and Hi-C arms. Either result (causal or parallel) is informative and publishable.
- Requires HPC access to extract CpG-level methylation from Zarr stores at specific coordinates.
- Pre-manuscript priority per FURTHER_ANALYSIS.

---

## 16. Repeat Element / Transposable Element Methylation

**Source:** FURTHER_ANALYSIS #4

**Core question:** The entire analysis pipeline currently ignores repeat elements (~40% of the mouse genome, major DNMT3A substrates). If DNMT3A is being pulled toward H2AK119ub-marked regions (UDR recruitment), DNMT3A may be depleted *from* its normal repeat targets, causing paradoxical repeat hypomethylation simultaneous with gene body hypermethylation. This parallels the Lopez-Moyado DNMT redistribution model but with different geography.

### Tasks

- [ ] **16a. Extract methylation at RepeatMasker elements.** Get mC/hmC levels at LINE, SINE, LTR, and DNA transposon families from Zarr stores or by overlapping DMR coordinates with RepeatMasker annotations (available for mm10).
- [ ] **16b. Compare WT vs KO by repeat family.** Stratify by LINE1, SINE/B1, SINE/B2, LTR/ERV, DNA transposons. Test for hypomethylation (DNMT3A depletion prediction) or hypermethylation.
- [ ] **16c. Stratify by repeat age/divergence.** Young repeats (low divergence from consensus, actively repressed by DNMT3A) should be more affected than ancient repeats (epigenetically inert). RepeatMasker provides divergence scores.
- [ ] **16d. Test DNMT3A redistribution model.** If repeats show hypomethylation while gene bodies show hypermethylation, that parallels Lopez-Moyado. If repeats are unaffected, DNMT3A is not being depleted from its normal targets.

### Existing resources

- **RepeatMasker annotations:** Available for mm10 from UCSC (`rmsk.txt.gz`)
- **Zarr stores:** 6bp resolution on Expanse
- **Gene body DMR data:** From `_shared_config.R`

### Notes

- Repeat de-repression (LINE1/ERV activation) can cause genomic instability and neuroinflammation -- a disease-relevant finding if present.
- Quick to compute if using pre-annotated RepeatMasker BED files overlapped with modality output.
- Pre-manuscript priority per FURTHER_ANALYSIS.

---

## 17. Polycomb Target Gene Systematic Enrichment

**Source:** FURTHER_ANALYSIS #5

**Core question:** The dual-mechanism model predicts that classic Polycomb target genes (repressed, H3K27me3-marked) should NOT be the primary hypermethylation targets because they are in heterochromatin and inaccessible to DNMT3A. Instead, genes that are normally active but gain ectopic H2AK119ub should be hypermethylated. Section 10 uses chromatin state proxies but no systematic test against published Polycomb target gene lists has been done.

### Tasks

- [x] **17a. Obtain published Polycomb target gene lists.** Curated PRC1/PRC2 target gene lists for mouse ESCs and neural progenitors (Boyer et al. 2006, Ku et al. 2008, or more recent mouse brain datasets). *(Done: section_30. Used 3 lab-data-based definitions instead of published lists: Chromatin State Polycomb (Repressed+Polycomb+Bivalent, n=3,445), Strict Polycomb (no Bivalent, n=3,294), Broad H3K27me3 gene body overlap (n=3,804). Published list infrastructure exists — will auto-load from `data/polycomb_gene_lists/*.txt` if placed there. All 3 definitions give consistent results.)*
- [x] **17b. Test Polycomb target enrichment in DMR list.** Fisher's exact: are classic Polycomb targets over- or under-represented among hypermethylated genes? Among hypomethylated genes? *(Done: section_30, Step 4. 19 Fisher's tests with BH correction. Polycomb massively DEPLETED from hypermethylation: OR=0.064 (Chromatin State), 0.040 (Strict), 0.105 (H3K27me3), all p<2.2e-16. Polycomb ENRICHED in hypomethylation: OR=8.71, 8.33, 9.08. hmC reciprocal: Polycomb enriched in hmC gain (OR=10.8) and depleted from hmC loss (OR=0.14).)*
- [x] **17c. Compare Polycomb vs non-Polycomb gene methylation.** If classic Polycomb targets are NOT preferentially hypermethylated but non-Polycomb genes are, supports "ectopic H2AK119ub at active genes." If classic Polycomb targets ARE hypermethylated, suggests Polycomb-mediated compaction driving methylation. *(Done: section_30, Step 5. Polycomb targets show significantly SMALLER effect sizes in all 4 comparisons (mC hyper, mC hypo, hmC hyper, hmC hypo; all Wilcoxon p<1e-6). Confirms heterochromatin is inaccessible to DNMT3A — non-Polycomb genes are the primary targets.)*
- [x] **17d. Cross-reference with Section 10 chromatin states.** Validate consistency with the Active_Promoter-dominated hypermethylation pattern (49.9% of DMRs at Active_Promoter; 94% hypermethylated). *(Done: section_30, Step 4 per-state tests + Step 6 figure 30d. Per-state hyper rates: Active_Promoter 65.2% (OR=9.10***), Active_Enhancer 49.0% (OR=2.07*), Bivalent_Promoter 32.5% (ns), Other 21.1% (OR=0.35***), Poised_Enhancer 19.6% (OR=0.52*), Polycomb 3.6% (OR=0.08***), Repressed_Promoter 2.3% (OR=0.04***). Genome-wide: 31.7%. Fully consistent with Section 10.)*

### Existing resources

- **Section 10 outputs:** `plots/visualizations/10{a-e}_*/` (chromatin state distributions)
- **Chromatin state classification:** `_shared_config.R` `classify_chromatin_state()`
- **H3K27me3 peaks:** `../../peaks/beds/H3K27me3CerebellumLate1.bed` (15,809 peaks)

### Notes

- Biologically central to the BAP1 story and currently only addressed indirectly.
- Mouse ESC Polycomb targets are well-curated but may not perfectly represent cerebellar Polycomb targets. Use H3K27me3 peaks from the lab's own CUT&RUN as a complementary approach.

---

## 18. Spatial Spreading / Chromosomal Autocorrelation

**Source:** FURTHER_ANALYSIS #6

**Core question:** H2AK119ub is known to spread in cis along chromatin (Conway et al. 2021). If methylation changes are downstream of H2AK119ub spreading, hypermethylated genes should cluster in chromosomal neighborhoods rather than being randomly distributed. If DMRs are randomly scattered, that argues against spreading and for gene-intrinsic susceptibility (consistent with Section 3, baseline 5hmC as predictor).

### Tasks

- [ ] **18a. Runs test for DMR clustering.** Along each chromosome, order genes by position and mark as DMR or not. Apply runs test (Wald-Wolfowitz) to detect non-random clustering.
- [ ] **18b. Sliding-window autocorrelation.** Compute methylation change (delta_mC) in sliding windows along each chromosome (e.g., 500kb windows, 100kb steps). Test for spatial autocorrelation (Moran's I or similar).
- [ ] **18c. DMR domain boundary analysis.** If methylation "domains" exist (contiguous blocks of hypermethylated genes), do their boundaries correspond to TAD boundaries, CTCF sites, or compartment boundaries from Hi-C?
- [ ] **18d. Compare spatial pattern with K119ub spreading.** Overlay K119ub bigwig signal with DMR positions along chromosomes. Correlated spatial domains support H2AK119ub spreading driving methylation. Uncorrelated patterns support gene-autonomous susceptibility.

### Existing resources

- **DMR coordinates:** From `_shared_config.R` (`mc_dmr`, `hmc_dmr`)
- **K119ub bigwig signal:** From Section 18 preprocessing
- **TAD boundaries:** `../../tads/results/late/final/tadcompare_final_filtered.tsv`
- **CTCF motifs:** `../../peaks/ctcf_motifs_mm10.bed`
- **Compartment boundaries:** `../../tad_analysis/diffcompartments.txt`

### Notes

- A positive result (spatial clustering with boundaries matching TAD/compartment boundaries) would be a strong integrative finding.
- A negative result (random scattering) supports the gene-intrinsic model and is also informative.
- Pre-manuscript priority per FURTHER_ANALYSIS.

---

## 19. Methylation Entropy / Read-Level Heterogeneity

**Source:** FURTHER_ANALYSIS #2

**Core question:** With ~400M reads/sample at 6bp resolution, multiple reads cover each CpG. Methylation entropy distinguishes two scenarios that gene-level means cannot separate: (A) every cell shifts uniformly from 40% to 42% methylation (graded DNMT3A recruitment), vs (B) half the cells stay at 40% and half jump to 44% (stochastic TET failure or cell-type heterogeneity). These look identical at the mean level but have different entropy signatures.

### Tasks

- [ ] **19a. Extract read-level methylation calls from Zarr stores.** For each CpG in top affected gene bodies, obtain per-read methylation status.
- [ ] **19b. Compute Shannon entropy per CpG window.** Sliding window (4-8 CpGs) across gene bodies. Shannon entropy or epipolymorphism metric. Compare WT vs KO.
- [ ] **19c. Test entropy change at DMR vs non-DMR loci.** Entropy increase at DMR gene bodies in KO suggests mosaicism (stochastic TET failure or cell-type mixing). Constant entropy suggests uniform methylation gain (graded recruitment).
- [ ] **19d. Interpret for Cre mosaicism.** Math1-Cre targets granule cell progenitors. Clean methylation changes (low entropy) mean targeted cells are uniformly affected. Entropy increase suggests mosaicism in Cre-mediated knockout or cell-type mixing.

### Existing resources

- **Zarr stores:** 6bp resolution on Expanse (contain read-level methylation calls)
- **Top genes:** From Section 06 tables (Syt1, Tmem238, Mcu, etc.)

### Notes

- Computationally tractable from Zarr stores but requires extracting read-level data, not just summary statistics.
- Also relevant to interpreting the 15.4% discordant genes (Section 21) -- genuinely discordant in all cells, or cell-type heterogeneity artifact?
- May distinguish dual-mechanism arms: uniform shift (DNMT3A recruitment) vs stochastic failure (TET impediment).

---

## 20. Cross-Species Comparison (Human BAP1 Tumors)

**Source:** FURTHER_ANALYSIS #8

**Core question:** Field et al. (2019) profiled 5mC (not 5hmC) in BAP1-mutant uveal melanoma. Which genes are hypermethylated in both mouse cerebellum KO and human BAP1-mutant tumors? Shared hits suggest a conserved BAP1-dependent methylation program; divergent hits suggest tissue-specific responses.

### Tasks

- [ ] **20a. Obtain Field et al. (2019) DMR gene list.** Download published supplementary tables of hypermethylated genes in BAP1-mutant uveal melanoma.
- [ ] **20b. Map mouse-human orthologs.** Use biomaRt or NCBI HomoloGene to map mouse DMR genes to human orthologs.
- [ ] **20c. Overlap analysis.** Fisher's exact: are BAP1-KO mouse cerebellum DMRs enriched for genes also hypermethylated in human BAP1-mutant tumors?
- [ ] **20d. Characterize shared vs tissue-specific hits.** For shared genes: what pathways? For divergent genes: do tissue-specific programs explain the difference (neuronal vs melanocyte)?

### Existing resources

- **Mouse DMR data:** From `_shared_config.R` (`mc_dmr`)
- **Field et al. (2019):** Published supplementary tables (need to download)
- **Ortholog mapping:** biomaRt R package (Ensembl)

### Notes

- Low effort (~a few hours), high manuscript value.
- Field et al. only profiled 5mC (not 5hmC), so the comparison is limited to the mC arm.

---

## 21. Cell Type Deconvolution

**Source:** FURTHER_ANALYSIS #10

**Core question:** Math1-Cre targets granule cell progenitors, but cerebellum contains Purkinje neurons, interneurons, Bergmann glia, astrocytes, oligodendrocytes, etc. Bulk methylation is a weighted average across cell types. Are observed changes driven by within-cell-type methylation shifts or by cell-type composition changes (e.g., granule cell loss due to degeneration)?

### Tasks

- [ ] **21a. Obtain reference methylation signatures.** Find published single-cell methylation atlases or snRNA-seq-derived methylation proxies for cerebellar cell types (Luo et al. 2017 Nature, Liu et al. 2021 Nature).
- [ ] **21b. Run reference-based deconvolution.** Estimate cell-type proportions per sample. Test whether KO samples have altered composition vs WT.
- [ ] **21c. Adjust DMR analysis for composition.** If composition differs, re-run DMR calling with estimated cell-type fractions as covariates. Compare to unadjusted analysis.
- [ ] **21d. Interpret discordant genes.** The minority hypomethylated genes (15.4% discordant, Section 21 viz) could be genuinely hypomethylated in neurons OR a deconvolution artifact from cell death/loss. If granule cells are lost in KO, their methylation signature drops out, appearing as hypomethylation at granule-cell-specific genes.

### Existing resources

- **Bulk methylation data:** Zarr stores, modality output
- **Section 21 viz outputs:** Discordant gene characterization
- **Published references:** Luo et al. (2017) Nature; Liu et al. (2021) Nature
- **Deconvolution tools:** MethylCIBERSORT, BisqueRNA, CIBERSORTx

### Notes

- Important caveat analysis for the manuscript even if the conclusion is "composition effects are minimal."
- If KO cerebellum has significant cell death, compositional shifts could confound many findings.

---

## 22. H3K36me2/3 Conceptual Discussion

**Source:** FURTHER_ANALYSIS #7

**Core question:** H3K36me2 antagonizes PRC2 and recruits DNMT3A to intergenic regions, while H3K36me3 (deposited by SETD2 in actively transcribed gene bodies) recruits DNMT3B. If BAP1 loss and H2AK119ub accumulation alter H3K36me2 domains, that is another route to DNMT redistribution. This is primarily a conceptual gap -- even without data, it should be discussed in the manuscript (Weinberg et al. 2019 Nature, Streubel et al. 2018).

### Tasks

- [ ] **22a. Check H3K36me2/3 data availability.** Does the lab have or plan to generate H3K36me2/3 CUT&RUN?
- [ ] **22b. (If available) Overlap H3K36me2/3 with DMRs.** Test whether hypermethylated gene bodies gain or lose H3K36me3. If H3K36me3 is lost (less transcription) while mC increases, rules out SETD2/DNMT3B-mediated methylation and points to DNMT3A.
- [ ] **22c. (If unavailable) Manuscript discussion section.** Discuss H3K36me2/3 as an alternative/additional DNMT recruitment mechanism. The interplay between H2AK119ub and H3K36me2/3 could explain why gene bodies (H3K36me3-marked in active genes) are preferentially affected.

### Existing resources

- **Literature:** Weinberg et al. (2019) Nature, Streubel et al. (2018)
- **Gene body DMR data:** From `_shared_config.R`

### Notes

- If no H3K36me2/3 data, this becomes a manuscript-writing task rather than computational analysis.

---

## 23. Q4 Coordinated Gene Sub-Stratification

**Source:** Section 28 follow-up

**Core question:** Q4 (mC up/hmC dn) contains 5,708 genes — too heterogeneous to treat as a single block. Sub-stratifying by effect size, chromatin state, or expression response could separate primary targets (direct BAP1/TET mechanism) from secondary/downstream effects, and identify which subgroups drive the aggregate statistical signals in Section 28.

### Tasks

- [ ] **23a. Stratify by effect size magnitude.** Split Q4 into tertiles/quartiles by combined_effect (|mc_diff| + |hmc_diff|). Compare top vs bottom quartile across all 9 Section 28 dimensions. Strong-effect genes are more likely direct targets.
- [ ] **23b. Stratify by chromatin state.** Q4 genes at Bivalent/Repressed promoters (n~686) have a more direct mechanistic link to BAP1 (PRC1 deubiquitinase) than Q4 genes at Active promoters (n~3,950). Compare epigenomic profiles between chromatin state subgroups.
- [ ] **23c. Stratify by expression response.** Q4 genes with concordant mC up + expression down (functional responders) vs Q4 genes with mC up but no expression change (non-responders). Tests whether methylation gain is functionally consequential or compensated.
- [ ] **23d. Multivariate clustering within Q4.** Unsupervised clustering (k-means or hierarchical) on the Section 28 master table features to identify natural subgroups without pre-defining the stratification axis.

### Existing resources

- **Section 28 master table:** `plots/visualizations/tables/coordinated_gene_characteristics.tsv` (5,708 Q4 genes with all 9 dimensions)
- **Section 28 script:** `scripts/viz_sections/section_28_coordinated_mc_hmc_analysis.R`

### Notes

- **Requires expanded cohort (n>=4 per condition)** before sub-group statistical tests are meaningful. With n=2, splitting further reduces power below useful thresholds.
- Descriptive/exploratory analysis is feasible now but inferential tests should wait for the expanded cohort (TODO Section 10).

---

## Summary: Priority Order

Based on URS proposal commitments, mechanistic importance, data availability, and FURTHER_ANALYSIS recommendations:

**Pre-manuscript priority** (per FURTHER_ANALYSIS: "genuinely want done before writing the manuscript"):

| Priority | Analysis | Data Available? | Depends On |
|----------|----------|----------------|------------|
| **1** | Methylation x Hi-C loop anchors (Section 1) | Yes | Nothing |
| **2** | 5hmC/(5mC+5hmC) demethylation ratio (Section 14) | **14a,14b done; 14d partial** | Nothing |
| **3** | CTCF site methylation at loop anchors (Section 15) | Yes (HPC) | Nothing |
| **4** | A/B compartment methylation (Section 2) | Yes | Nothing |
| **5** | Repeat element methylation (Section 16) | Yes (HPC) | Nothing |
| **6** | Spatial autocorrelation (Section 18) | Yes | Nothing |
| **7** | Baseline 5hmC as predictor (Section 3) | Likely (check modality output) | Nothing |
| **8** | MeCP2 functional readout (Section 4) | **4a~, 4b~, 4c~ done** | Nothing |

**Would strengthen the paper:**

| Priority | Analysis | Data Available? | Depends On |
|----------|----------|----------------|------------|
| **9** | TF motif enrichment at DMRs (Section 5) | Yes (HPC) | Nothing |
| **~~10~~** | ~~Polycomb target gene enrichment (Section 17)~~ | **Done (Section 30)** | — |
| **11** | Expression-methylation model (Section 8) | Yes (Section 20 foundation) | Nothing |
| **12** | Non-CG / mCA methylation (Section 13) | Uncertain | Check run-3 output |
| **13** | Cross-species comparison (Section 20) | Yes (download paper data) | Nothing |
| **14** | Methylation entropy (Section 19) | Yes (HPC) | HPC access |
| **15** | Cell type deconvolution (Section 21) | Yes (reference data needed) | Nothing |
| **16** | ATAC-seq follow-ups (Section 9) | Yes (Section 13 foundation) | Nothing |
| **17** | DNMT3A binding prediction (Section 11) | Yes | Nothing |

**Requires external data or expanded cohort:**

| Priority | Analysis | Data Available? | Depends On |
|----------|----------|----------------|------------|
| **18** | Expanded cohort (Section 10) | No (waiting on new samples) | Lab |
| **19** | RNA splicing analysis (Section 7) | Needs BAM files | Lab data access |
| **20** | Single-CpG resolution (Section 6) | Yes (HPC) | HPC access |
| **21** | Developmental timepoint (Section 12) | No (no P12 DUET data) | Lab |
| **22** | H3K36me2/3 discussion (Section 22) | No data (manuscript task) | Lab / writing |
| **23** | Q4 sub-stratification (Section 23) | Yes (Section 28 table) | Expanded cohort (Section 10) |

---

## Data File Reference

### Methylation Data (primary inputs)

| Data | File / Location | Key Columns |
|------|-----------------|-------------|
| Gene body mC DMRs | `modality/outputs/run-3/outputs_CG/.../DMR_*/num_mc_dmr_results.bed.gz` | chr, start, end, gene, q-value, effect_size, direction |
| Gene body hmC DMRs | Same location, `num_hmc_dmr_results.bed.gz` | Same schema |
| Feature extraction | `modality/outputs/run-3/outputs_CG/.../Extract_*/` | Per-gene count, mean, regional-frac |
| Zarr stores (6bp) | Expanse: `evoC-run/output/.../zarr_store/` | CpG-level 5mC/5hmC fractions |
| BioQC JSON | `modality/outputs/run-3/outputs_CG/.../BioQC_*/biological_qc_data.json` | Sample correlations, PCA |

### Cross-Modal Data (for integration analyses)

| Data | File | Used In |
|------|------|---------|
| Differential loops (Late) | `../../25042-late_outputs/merged_loops/characterized_loops.tsv` | Sections 1, 4d, 9, 15 |
| Differential loops (Early) | `../../250831-early_outputs/merged_loops/characterized_loops.tsv` | Section 1 |
| Shared anchors | `../../output/shared_anchor_analysis/late/` | Section 1e |
| PC1/Compartments | `../../tad_analysis/all_PC1.txt`, `diffcompartments.txt` | Section 2 |
| K119ub bigwig signal | Via Section 18 preprocessing | Sections 1c, 3c, 8c, 11a, 18d |
| ATAC peaks | `../../peaks/atac_seq/ATAC_{up,down}.bed` | Section 9 |
| MeCP2 peaks | Via `_shared_config.R` MECP2_FILES | Sections 4, 13e |
| RNA-seq DEGs | `../../tads/adult_timepoint_rna-seq-*.xlsx` | Sections 7, 8 |
| ChIP-seq peaks | `../../peaks/beds/` (5 marks + CTCF) | Sections 5, 8b, 17 |
| CTCF motifs | `../../peaks/ctcf_motifs_mm10.bed` (114,081) | Section 15 |
| CTCF ChIP peaks | `../../peaks/CTCF.bed` (32,487) | Section 15d |
| TAD boundaries | `../../tads/results/late/final/tadcompare_final_filtered.tsv` | Section 18c |
| RepeatMasker annotations | UCSC mm10 `rmsk.txt.gz` (need to download) | Section 16 |
| Polycomb target gene lists | Published (Boyer et al. 2006, Ku et al. 2008) | Section 17 |
| Field et al. (2019) DMRs | Published supplementary tables | Section 20 |

### Template Scripts (reusable patterns)

| Pattern | Script | Reusable For |
|---------|--------|-------------|
| DMR loading + chromatin annotation | `scripts/viz_sections/_shared_config.R` | All sections |
| Multi-format plot output | `_shared_config.R` `save_multiformat_ggplot()` | All sections |
| GREAT-style gene-anchor association | `../../tads/scripts/deg_tad_violin.R` | Section 1 |
| Fisher's exact enrichment | Sections 11, 12, 13, 14 | Sections 1, 2, 4, 7, 17, 20 |
| Logistic regression | Section 18 / `../../scripts/h2ak119ub_loop_integration.R` | Sections 1d, 3b, 11b |
| Permutation testing | regioneR/regioneReloaded (Bioconductor) | Sections 9c, 18a |
| Ortholog mapping | biomaRt R package | Section 20 |
| Deconvolution | MethylCIBERSORT / BisqueRNA / CIBERSORTx | Section 21 |
