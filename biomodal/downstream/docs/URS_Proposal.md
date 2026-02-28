
# Proposal
## Background

BAP1 encodes a tumor-suppressive deubiquitinase that removes H2AK119ub, a Polycomb-deposited histone modification that silences genes. Its loss causes severe neurodevelopmental defects including epilepsy, but how it shapes DNA methylation in the brain is wholly unknown. When BAP1 is absent, H2AK119ub accumulates genome-wide and chromatin compacts, reducing overall accessibility (Conway et al., 2021). In the brain, DNA carries two functionally distinct modifications: 5-methylcytosine (5mC), which represses transcription, and 5-hydroxymethylcytosine (5hmC), which is generated when TET enzymes oxidize 5mC in the first step of active demethylation. Cerebellar neurons contain the highest 5hmC of any mammalian cell type, approximately 40% of modified cytosines (Kriaucionis and Heintz, 2009). This makes the cerebellum especially dependent on the 5mC/5hmC balance, and especially vulnerable when that balance is disrupted. Despite this, no study has simultaneously profiled both modifications in any BAP1-KO model, and no similar methylation study exists in brain tissue (Parmar et al., 2025). This project will ask how BAP1 loss disrupts active demethylation in neurons, and whether the resulting changes converge with the three-dimensional chromatin remodeling we have independently identified in the same system.

## Goal

This project will primarily validate and extend preliminary methylation findings using an expanded sample size, with biological replicates and deeper sequencing that resolve current statistical limitations. Secondarily, we will integrate methylation data with chromatin loop architecture and histone modification peaks to determine whether epigenetic and structural changes can be causally linked to the chromatin conformation consequences of BAP1 loss. Lastly, we intend to computationally test predictions of a dual-mechanistic hypothesis, in which H2AK119ub accumulation simultaneously recruits DNA methyltransferases and impedes TET-mediated demethylation at the same developmental loci.

## Significance

No previous BAP1 study has profiled both 5mC and 5hmC, and no BAP1 methylation study has been conducted in brain tissue. This is despite 5hmC being unusually abundant in neurons and critical for gene regulation. TET-mediated 5hmC production is required for cerebellar neuron differentiation (Stoyanova et al., 2021), and 5hmC loss in Purkinje cells drives neurodegeneration (Jiang et al., 2015), so our observed disruption in DNA methylation is likely to be the cause of functional dysregulation. Affected genes are enriched for neurodegenerative disease pathways including spinocerebellar ataxia, amyotrophic lateral sclerosis, and Alzheimer's disease, which suggests that altered methylation from BAP1 loss may contribute to neuronal vulnerability, in addition to developmental and synaptic changes. By linking opposing DNA methylation trends to three-dimensional genome architecture in the same genetic model, this work connects Polycomb chromatin biology in the field of Hi-C to dual-epigenetic methylation dynamics in neurons. This is a gap that current literature does not address, and has novel implications for studying the role of BAP1-associated neurodegeneration in mice.


## Work

Since joining the Ferguson Lab in September 2025, I have built the computational infrastructure for both arms of this project, beginning with Hi-C chromatin conformation, and in December 2025, for dual-epigenetic methylation sequencing. For the methylation side, I developed a 22-section analysis pipeline for Biomodal DUET evoC data. The technology quantifies both 5mC and 5hmC simultaneously at 6-base-pair resolution, calling differentially methylated regions. My pipeline takes these base-level calls through chromatin state annotation using matched CUT&RUN peaks, and integrates with ATAC-seq (accessibility), RNA-seq (gene expression), and further histone modification data. In the initial sample set, we identified a coordinated pattern, in which 5mC increases while 5hmC decreases at the same gene body loci. This trend currently appears to affect 6,750 loci with 85% concordance among co-significant genes. Gene bodies are the primary affected compartment (42% of genes significant for 5mC, 47% for 5hmC) while promoters show almost no change (<0.2%), consistent with the known enrichment of 5hmC in neuronal gene bodies (Szulwach et al., 2011; Mellén et al., 2017). Synaptotagmin-1, the primary calcium sensor for fast neurotransmitter release, is the most affected gene (+18% mC, −15% hmC), and has been noted for its clear epigenetic changes in past experiments.

Furthermore, when we examine the presence of ubiquitinated histone, we see that H2AK119ub gain co-occurs with hypermethylation (Fisher's OR = 4.40), while these loci are biased toward transcriptional downregulation, and are inversely correlated with chromatin accessibility. Together, these results are consistent with a dual-mechanism model, in which H2AK119ub accumulation both recruits DNMT3A/B to accessible gene bodies via its ubiquitin-recognition domain (Gretarsson et al., 2024), and globally reduces TET enzyme efficiency or access (López-Moyado et al., 2019). I also developed a parallel differential chromatin loop pipeline for Hi-C data from the same BAP1-KO model, which identified 2,910 dysregulated loops, and showed progressive architectural remodeling from early postnatal development to adulthood (18x increase). Because the initial cohort has the limitation of two samples per condition, and sex is confounded with genotype, we will sequence two additional biological replicates each, and intend to double the read count in order to reduce noise from a lower sequencing depth.

## Methods

The full dataset (n=4 per genotype) will be processed through the established pipeline with additional analyses planned dependent on exploratory results. Differential methylation calling will also incorporate batch as a covariate, and determine significance of MeCP2 binding affinity. To connect methylation changes to three-dimensional genome structure, I will test whether genes with the coordinated hyper-methylation and hypo-hydroxymethylation pattern are enriched at the differential loop anchors identified in our Hi-C analysis, using continued Fisher’s exact tests and logistic regression to quantify any association between methylation direction and chromatin remodeling. To test predictions of the dual-mechanism model, I will evaluate whether gene-body hypermethylation continues to track with local H2AK119ub levels, consistent with DNMT3A/B recruitment via its ubiquitin-dependent domain (Conway et al., 2021; Gretarsson et al., 2024); or whether it occurs independently, consistent with global TET impediment. 
All computation will use SDSC Expanse through existing ACCESS allocations, implemented in R/Bioconductor, Python, and Bash. Weeks 1-3 will refine the expanded pipeline, investigate the detection of non-CG methylation with bisulfite sequencing, and perform upstream processing of the new batch. Weeks 4-6 will proceed to run differential methylation analysis with a full cohort, and update current results. Remaining time will be dedicated to further downstream analysis.

## References

Conway, E. et al. (2021). BAP1 enhances Polycomb repression by counteracting widespread H2AK119ub1 deposition and chromatin condensation. Molecular Cell, 81(17), 3525–3541.

Gretarsson, K.H. et al. (2024). Cancer-associated DNA hypermethylation of Polycomb targets requires DNMT3A dual recognition of histone H2AK119 ubiquitination and the nucleosome acidic patch. bioRxiv. doi:10.1101/2024.03.18.585588

Jiang, D. et al. (2015). Alteration in 5-hydroxymethylcytosine-mediated epigenetic regulation leads to Purkinje cell vulnerability in ATM deficiency. Brain, 138(12), 3520–3536.

Kriaucionis, S. and Heintz, N. (2009). The nuclear DNA base 5-hydroxymethylcytosine is present in Purkinje neurons and the brain. Science, 324(5929), 929–930.

López-Moyado, I.F. et al. (2019). Paradoxical association of TET loss of function with genome-wide DNA hypomethylation. PNAS, 116(33), 16933–16942.

Mellén, M. et al. (2017). 5-Hydroxymethylcytosine accumulation in postmitotic neurons results in functional demethylation of expressed genes. PNAS, 114(37), E7812–E7821.

Parmar, A. et al. (2025). Polycomb repressive complexes 1 and 2 independently and dynamically regulate euchromatin during cerebellar neurodevelopment. PLOS Genetics, 21(9), e1011843.

Stoyanova, E. et al. (2021). 5-Hydroxymethylcytosine-mediated active demethylation is required for mammalian neuronal differentiation and function. eLife, 10, e66973.

Szulwach, K.E. et al. (2011). 5-hmC-mediated epigenetic dynamics during postnatal neurodevelopment and aging. Nature Neuroscience, 14(12), 1607–1616.

