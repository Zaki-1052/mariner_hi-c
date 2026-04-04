## Summary: BAP1 Loss Drives Epigenetic Reprogramming in Aggressive Uveal Melanoma

This study by Field et al. (2019, *Clin Cancer Res*) investigates **why loss of an entire copy of chromosome 3 is required** for the highly metastatic Class 2 phenotype in uveal melanoma (UM), and whether BAP1 loss specifically drives the epigenetic changes that define it.

### Background & Core Question

UM divides into two prognostic subtypes: Class 1 (low metastatic risk, BAP1-wildtype) and Class 2 (high metastatic risk, BAP1-mutant). Class 2 tumors lose one chromosome 3 copy entirely, then acquire a BAP1 mutation on the remaining copy — but it was unclear *why* whole-chromosome loss is needed, since no other commonly mutated tumor suppressors exist on chromosome 3 besides BAP1. The authors hypothesized that the retained chromosome 3 carries **epigenetic alterations** (specifically DNA methylation changes) that cooperate with BAP1 loss to produce the Class 2 phenotype.

### Approach

They performed genome-wide DNA methylation profiling (450K array) on 92 primary UMs plus 80 TCGA samples, integrated this with RNA-Seq data, then functionally validated findings using inducible BAP1 knockdown (BAP1KD) in two UM cell lines (Mel202, 92.1). Orthogonal validation used bisulfite sequencing on an independent 14-tumor cohort.

### Key Findings

**1. Class 2 UMs show extensive, non-random methylomic repatterning.** PCA clearly separated Class 1 and Class 2 tumors by methylation profile alone. Chromosome 3 showed the most significant promoter hypermethylation enrichment — particularly at loci 3p21–23, 3p25–26, 3q12–21, and 3q27.

**2. Methylation changes are functionally linked to gene silencing.** Integrating methylation with expression identified 508 hypermethylated/downregulated genes and 923 hypomethylated/upregulated genes in Class 2 tumors. Hypermethylation clustered in CpG shore and shelf regions near promoters, which are known to inversely correlate with expression.

**3. Axon guidance and melanogenesis are the top deregulated pathways.** This makes biological sense: Class 1 UMs resemble differentiated melanocytes, while Class 2 UMs resemble undifferentiated neural crest-like progenitors. Key silenced genes on chromosome 3 include ROBO1, PLXNB1, SEMA3B, CHL1, SATB1, MITF, DVL3, and RAF1 — all involved in neural crest migration or melanocyte differentiation.

**4. BAP1 itself is epigenetically regulated.** A novel hypermethylated CpG site (cg16871520) within the BAP1 gene body showed a strong inverse correlation with BAP1 expression (Spearman R = −0.79) and correctly classified 79/80 tumors by class.

**5. BAP1 knockdown directly causes methylomic repatterning.** Inducible BAP1KD in cell lines produced methylation changes enriched for the same pathways (axon guidance, melanogenesis, development) and the same chromosome 3p21 locus. The overlap wasn't perfect — expected, since cell lines lack whole-chromosome-3 loss, microenvironment effects, etc. — but the pathway-level concordance was striking.

**6. Methylation probes can distinguish tumor classes with 100% accuracy.** Bisulfite sequencing of IL12RB2, SATB1, SESN1, and ENPP2 probes correctly classified all 14 independent tumors, establishing proof-of-concept for a potential methylation-based liquid biopsy.

### Proposed Evolutionary Model

The authors propose a sequential model: **(1)** Loss of one chromosome 3 copy → **(2)** BAP1 mutation on the remaining copy → **(3)** BAP1 loss triggers methylomic repatterning that silences melanocyte differentiation genes and activates neural crest migration programs → **(4)** Class 2 phenotype with high metastatic potential. This explains why complete chromosome 3 loss (not just regional BAP1 deletion) is necessary — the retained chromosome must carry both the BAP1 mutation *and* the resulting epigenetic alterations for the full Class 2 program to emerge.

### Translational Implications

The identified methylation signatures could serve as liquid biopsy biomarkers (avoiding invasive tissue biopsies), and the centrality of epigenetic reprogramming to Class 2 progression suggests DNA methylation modulators (e.g., DNMT inhibitors) as potential therapeutic strategies.

### Limitations Worth Noting

The BAP1KD cell model only partially recapitulates tumor methylation patterns, which the authors acknowledge transparently — cell lines are heterozygous for chr3, lack the tumor microenvironment, and BAP1 knockdown ≠ whole-chromosome loss. The liquid biopsy concept remains proof-of-concept with significant technical challenges ahead.

---

Integrated analysis of DNA methylation and RNA expression
To further explore the potential functional relevance of these findings, we then focused on integrating methylation with gene expression. As hypermethylation in promoter regions is usually associated with gene silencing (37), we identified genes with hypermethylation associated with decreased gene expression (hypermethylated/downregulated) or hypomethylation associated with increased gene expression (hypomethylated/upregulated) in Class 2 UMs relative to Class 1 UMs (Fig. 2A). This analysis revealed 1621 hypermethylated probes associated with 508 downregulated genes (FDR < 0.05) and 3876 hypomethylated probes associated with 923 upregulated genes (FDR < 0.05) (Supplementary Table S3). Out of the twelve genes in the GEP test, six of the eight downregulated genes were hypermethylated (FXR1, ID2, ROBO1, LMCD1, SATB1, and MTUS1) and all four of the upregulated genes were hypomethylated (HTR2B, ECM1, RAB31, CDH1). Similar to our global analysis, hypermethylated/downregulated genes in Class 2 tumors were enriched within the promoter and 5’UTR “shore” and “shelf” regions (Fig. 2B), and hypomethylated/upregulated genes in Class 2 tumors occurred mostly in open sea regions (38). Chromosomal regions that were significantly enriched (FDR < 0.05) for hypermethylated/downregulated genes included 6p21, 6p24, 19q13, 10q24, 4p14, as well as multiple regions on chromosome 3 (3p21–23, 3p25–26, 3q12–21, and 3q27). Regions enriched (FDR < 0.05) for hypomethylated/upregulated genes included 6p21, 17q21, 15q21, and 12p13, 22q13, 1q31, and several regions on chromosome 8 (8p21, 8q13, 8q21–22)(Fig. 2A–C and Supplementary Table S4). Interestingly, 6p21 contains regions of both hypomethylated/upregulated and hypermethylated/downregulated genes in Class 2 tumors, but only the former is enriched for HLA genes, which are known to be expressed as part of the altered immune microenvironment in Class 2 tumors.

Integration of DNA methylation with RNA expression data. (A) Inner quadrant graph demonstrates methylation of individual probe sites (x-axis, Class 2 Beta – Class 1 Beta) plotted against gene expression (y-axis, log2 fold change Class 2/Class 1) in 80 TCGA samples. Red, hypermethylated/downregulated genes; blue, hypomethylated/upregulated genes; black, other genes. Outer Circos plot demonstrates hypermethylated/downregulated genes (red bars) and hypomethylated/upregulated genes (blue bars) that met filtering cutoffs (Delta Beta ± 0.05 and log2 fold change ± 1) with respect to chromosomal location. (B) Global CpG-island feature analysis of the filtered hypermethylated/downregulated genes and hypomethylated/upregulated genes showing the normalized percent of hypermethylated (red) and hypomethylated (blue) probe sites. (C) Gene set enrichment analysis of chromosomal location for hypermethylated/downregulated genes (red lines) and hypomethylated/upregulated genes (blue lines) in Class 2 UMs compared to Class 1 UMs from the TCGA dataset. Copy number gains (blue circles) and losses (red circles) that commonly occur in UM on chromosomes 3, 6, and 8 are indicated next to the enriched methylated chromosomal regions. CpG-islands are defined as regions > 500 base pairs with > 55% GC content and an expected/observed CpG ratio of > 0.65. CpG shores are ~2Kb from islands and CpG shelves are ~4Kb from islands. *p < 0.05, **p < 0.01, ***p < 0.001 with binomial testing.

