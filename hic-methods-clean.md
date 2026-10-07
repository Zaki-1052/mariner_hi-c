# Hi-C Methods

## Hi-C sample structure

Hi-C was performed on cerebella from Bap1fl/fl;Math1-Cre conditional knockout and littermate control mice using the Arima Hi-C protocol. Three biological replicates per genotype (six libraries total) were deep-sequenced at two developmental timepoints: early (P12/P13) and late (adult). Merged per-condition contact maps were generated for loop calling, stripe detection, and APA. Individual replicate libraries were used for all differential analyses.

## Hi-C data processing

### Full version

Deep-sequenced Hi-C libraries were trimmed with fastp (v1.0.1) and aligned to the mm10 mouse genome with the nf-core/hic pipeline (v2.1.0) run on Nextflow (v25.04.6). The pipeline used Bowtie2 (v2.5.4) for alignment with the Arima digestion protocol. The pipeline produced valid-pair files, BAM files, and cooler-format contact maps at bin sizes from 1,000 bp to 1,000,000 bp. TADs were called with insulation scoring and hicExplorer at 40, 20, and 10 kb resolutions. Compartments were computed at 1,000, 250, and 100 kb. The pipeline ran in Singularity containers on SDSC Expanse. BAM files were processed with SAMtools (v1.22). Valid pairs were converted to Juicer v2.0 format .hic files with juicerpre, which computed VC, VC_SQRT, GW_SCALE, SCALE, and INTER_SCALE normalization vectors.

### Brief version

Deep-sequenced Hi-C libraries were trimmed with fastp (v1.0.1) and aligned to the mm10 genome with the nf-core/hic pipeline (v2.1.0, Bowtie2 v2.5.4). Valid pairs were converted to Juicer v2.0 format .hic files with juicerpre.

## Chromatin loop detection

Chromatin loops were called with HiCCUPS from Juicer Tools (v2.0.1) on GPU nodes at 5, 10, and 25 kb resolutions. Default parameters were used. Each of the six individual replicate .hic files was processed separately.

## Differential loop analysis

Loop calls from six replicates were processed with the mariner R/Bioconductor package (v1.8.1) at each resolution independently. HiCCUPS BEDPE files were merged across replicates with mergePairs using a 10 kb radius, Manhattan distance, and selectMax on observed contact count. Merged loops were assigned to resolution-specific bins (assignToBins, center positioning) and buffered to a 5x5 pixel window (pixelsToMatrices, buffer = 2). Observed 5x5 contact matrices were extracted from the six .hic files at each loop position with pullHicMatrices, using VC normalization and on-disk HDF5 storage. Each matrix was collapsed to a single count per loop per sample by summing all 25 pixels. chrX was excluded from the analysis.

The differential analysis used the edgeR quasi-likelihood GLM (v4.4.0) with a ~group design (control and mutant). Low-count loops were removed with filterByExpr (min.count = 5, min.total.count = 15, min.prop = 0.5). Library sizes were normalized with TMM, and dispersions were estimated with estimateDisp (robust = TRUE). The model was fitted with glmQLFit (robust = TRUE) and tested with glmQLFTest on the mutant-vs-control coefficient. P-values were corrected by Benjamini-Hochberg FDR at a primary threshold of FDR < 0.05. A second filter (|logFC| > 0.3, FDR < 0.03) defined the stringent loop set used for downstream analysis.

## Multi-resolution loop comparison

Loops called at 5, 10, and 25 kb were matched across resolutions. Matches required the same chromosome pair and both anchor start positions within 10 kb. The match with the smallest summed anchor distance was selected.

## Loop anchor annotation

Each loop anchor was classified into one of 8 chromatin-state categories based on overlap with CUT&RUN peak files. The classification used five histone marks (H3K27ac, H3K27me3, H3K4me1, H3K4me3, and a pre-computed bivalent set) and CTCF DNA motifs. The TSS distance threshold was 2 kb (TxDb.Mmusculus.UCSC.mm10.knownGene v3.10.0). Categories were assigned by priority. Anchors within 2 kb of a TSS were classified as Active_Promoter (H3K4me3+ without H3K27me3), Repressed_Promoter (H3K27me3+ without H3K27ac), or Bivalent_Promoter (K4me3+K27me3 overlap). Anchors beyond 2 kb were classified as Polycomb (H3K27me3+), Active_Enhancer (H3K27ac+), or Poised_Enhancer (H3K4me1+ without H3K27ac or H3K27me3). Remaining anchors were classified as CTCF_Site (CTCF motif overlap) or Other (no marks). The nearest gene was assigned per anchor from GENCODE vM25 annotations. Each loop received a pairwise type label from its two anchor categories. Loops were further classified by anchor composition as structural, CRE, mixed, or unclassified based on CTCF and H3K27ac/promoter overlap at both anchors.

## Aggregate peak analysis

Aggregate Peak Analysis (APA) was computed with GENOVA (v1.0.1). Merged contact maps (three replicates combined per condition) were used as input with loop coordinates in BEDPE format.

## Chromatin stripe detection

.hic files were converted to multi-resolution .mcool format with hictk (v2.1.4) at 5 and 10 kb resolutions. Chromatin stripes were called with Stripenn (v1.1.65.22) on merged per-condition .mcool files using KR normalization. Detection parameters were: minimum length 10 bins, maximum width 8 bins, each maxpixel percentile from 0.95 to 0.99, p-value threshold 0.1, background filter 3, and seed 42. Stripes on chrX, chrY, and chrM were removed post-hoc. Control and mutant stripe calls were merged into a union set per resolution with a 50 kb anchor tolerance and direction-aware overlap matching. Each of the six individual replicate .mcool files was scored against the union set with stripenn score, producing per-stripe observed contact sums.

## Differential stripe analysis

The differential stripe analysis used the same edgeR quasi-likelihood GLM as the loop pipeline. Two parameters were different: filterByExpr was not applied, and TMM normalization used calcNormFactors. The primary significance threshold was FDR < 0.05. Stripes detected in only one condition were labeled lost (control-only) or gained (mutant-only). Shared stripes were labeled strengthened or weakened (FDR < 0.05, |logFC| > 0.3) or unchanged. Stripe anchors at 5 and 10 kb were extended by 50 kb on each side and matched per timepoint with findOverlaps. Concordant stripes were significant at both resolutions with the same direction label. Stripe anchors were classified with a 7-category system (CTCF_Site excluded) using the same five histone marks and 2 kb TSS threshold as the loop annotation. Two rules differed: Bivalent_Promoter was checked at priority 1, and Active_Promoter required both H3K4me3+ and H3K27ac+.

## TAD and compartment analysis

TAD boundaries and inclusion ratios were computed with HOMER (v5.1). Tag directories were generated from Hi-C BAM files. TADs were called per sample and merged into a consensus set across all replicates. Inclusion ratios (interaction density) were scored per TAD per sample from the consensus set. Differential TAD analysis was performed with HOMER's getDiffExpression.pl. Significance was assessed at two levels: relaxed (FDR < 0.15, |Difference| > 0.15) and standard (FDR < 0.05, |Difference| > 0.30).

Differential TAD boundaries were detected with TADCompare (v1.16.0) at 10 kb resolution on merged control and mutant sparse contact matrices. Contact matrices were extracted from .hic files with Juicer straw using NONE normalization for chr1-chr19. TADCompare classified each boundary into one of six types: Non-Differential, Strength Change, Shifted, Split, Merge, or Complex. Boundaries overlapping a lab-specific blacklist (254 regions) or the ENCODE mm10 blacklist (3,435 regions) were removed. Shift distances for Shifted boundaries were computed with bedtools closest, with a minimum of 10 kb (one resolution bin). Differential loops were tested for localization near differential boundaries by Fisher's exact test, bootstrap, and permutation (1,000 iterations each). Distance thresholds of 10, 25, 50, and 100 kb were used. Boundary chromatin states were classified with the same 8-category system used for loop anchor annotation.

Compartment (PC1) values were computed per sample at 25x50 kb resolution with HOMER and stored as bedGraph files. Differential compartment analysis used getDiffExpression.pl with the same dual-threshold approach as TADs. PC1 polarity was corrected per chromosome to ensure positive values corresponded to the A (active) compartment. The percentage of the genome with compartment shifts was computed against the mm10 assembled size (2,730,871,774 bp).

## Gene ontology and pathway enrichment

Gene Ontology and KEGG pathway enrichment were computed with clusterProfiler (v4.16.0). Entrez gene IDs were identified within 10 kb of up-regulated or down-regulated loop anchors. Enrichment was tested with enrichGO for three ontologies (BP, CC, MF) and enrichKEGG for KEGG pathways, using org.Mm.eg.db (v3.21.0). Significance thresholds were p < 0.05 and q < 0.05 with Benjamini-Hochberg correction.

---

## Software versions

| Software                           | Version       |
| ---------------------------------- | ------------- |
| **Upstream processing**            |               |
| nf-core/hic                        | 2.1.0         |
| Nextflow                           | 25.04.6       |
| fastp                              | 1.0.1         |
| Bowtie2                            | 2.5.4         |
| SAMtools                           | 1.22          |
| BEDTools                           | 2.31.1        |
| Juicer Tools                       | 2.0.1         |
| hictk                              | 2.1.4         |
| HOMER                              | 5.1           |
| GENOVA                             | 1.0.1         |
| **R packages**                     |               |
| R                                  | 4.4.2         |
| Bioconductor                       | 3.20          |
| mariner                            | 1.8.1         |
| edgeR                              | 4.4.0         |
| limma                              | 3.62.1        |
| InteractionSet                     | 1.36.1        |
| GenomicRanges                      | 1.58.0        |
| strawr                             | 0.0.92        |
| HDF5Array                          | 1.36.0        |
| TADCompare                         | 1.16.0        |
| rtracklayer                        | 1.70.1        |
| ChIPseeker                         | 1.44.0        |
| clusterProfiler                    | 4.16.0        |
| EnhancedVolcano                    | 1.26.0        |
| TxDb.Mmusculus.UCSC.mm10.knownGene | 3.10.0        |
| org.Mm.eg.db                       | 3.21.0        |
| **Python packages**                |               |
| Python                             | 3.13.5        |
| Stripenn                           | 1.1.65.22     |
| cooler                             | 0.10.4        |
| NumPy                              | 2.4.1         |
| SciPy                              | 1.17.1        |
| pandas                             | 2.3.3         |
| OpenCV                             | 4.13.0        |
| scikit-image                       | 0.26.0        |
| **Genome and annotation**          |               |
| Genome assembly                    | mm10 (GRCm38) |
| Gene annotation                    | GENCODE vM25  |

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
