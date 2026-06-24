**Loss of BAP1 Deubiquitinase Restructures MeCP2 binding and DNA Methylation Profiles During Cerebellar Neurodevelopment**

Jai Ramchandra, Zakir Alibhai, Kelly C Wang, Zhendong Song, Rachael Cui, Cole J Ferguson 

**ABSTRACT**

**INTRODUCTION**

* Bap1 is the main deubiquitinase part of the PR-DUB complex. It actively removes Ub, which generally is thought to remodel chromatin to euchrom. Bap1 functions both in the nucleus and the cytoplasm. UBE20 monoubs bap1 at several points in its NLS, inactivating it and localizing it to cytoplasm. To enter the nucleus, Bap1 can deub itself and reactivate the NLS. The PR-DUB subunits (BAP1 and ASXL1) are among the most frequently mutated epigenetic factors in human cancers. Bap1 has been implicated in cancer and thought to mostly be a tumor suppressor, but it has a wide range of activity in the cell.  
  * [https://pmc.ncbi.nlm.nih.gov/articles/PMC7862696/](https://pmc.ncbi.nlm.nih.gov/articles/PMC7862696/)  
* BAP1 participates in a wide range of biological processes by directly targeting proteins as substrates or indirectly via transcription, and these BAP1-regulated processes include cell cycle control, cell survival and proliferation, cell death, the DNA damage response, repair, DNA replication, metabolism, and cell differentiation and development \+ embryogenesis.  
  * [https://www.nature.com/articles/s12276-023-00979-1](https://www.nature.com/articles/s12276-023-00979-1)  
* BAP1 and ca2+ signaling (maybe related to syt1?) → BAP1 promotes apoptosis by deubiquitinating and stabilizing IP3R3 at the endoplasmic reticulum (ER), which stimulates Ca2+ release from the ER into the cytosol and thereby increases the mitochondrial Ca2+ concentration and cytochrome C release. BAP1 also indirectly regulates apoptosis and ferroptosis, a recently identified nonapoptotic form of cell death, by regulating the transcription of the genes critical for these processes.  
  * [https://www.nature.com/articles/s12276-023-00979-1](https://www.nature.com/articles/s12276-023-00979-1)  
* MECP2 binds to methylated DNA (especially methylated CpG sites) and helps regulate transcription. It modulates gene expression, particularly in neurons. It is Critical for neuronal maturation, synaptic function, and plasticity. MECP2 helps condense chromatin and influence the accessibility of genes. Mutations in MECP2 cause Rett syndrome, a severe neurological disorder that primarily affects girls.  
* MeCP2 is very particular with how it binds. It binds to single CG dinucleotides primarily and has an MBD domain \+ lot of disordered regions that people propose help to anchor and recognize binding sequences. It scans genome for chromatin that has exposed linker dna between nucleosomes and can access methylation on the surface and also core. It also prefers to bind to dsDNA over ssDNA \+ methylated over unmethylated DNA. It strongly recognizes mCA, then mCG, then the other non-CGs. It also is affected by sequence context, CG content, \+ gene length.  
  * [https://www.sciencedirect.com/science/article/pii/S0021925826020995](https://www.sciencedirect.com/science/article/pii/S0021925826020995)  
  * [https://www.nature.com/articles/s41593-025-01947-w](https://www.nature.com/articles/s41593-025-01947-w)  
  * [https://academic.oup.com/bfg/article/15/6/420/2555368](https://academic.oup.com/bfg/article/15/6/420/2555368)  
  * https://pmc.ncbi.nlm.nih.gov/articles/PMC9779294/  
* DNMT3a is a de novo methylase that can recognize H2AK119Ub and add 5mC (something about acid patch of histones..?). It also has unique dehydroxymethylase activity. DNMT1 does NOT have this activity.  
  * [https://pmc.ncbi.nlm.nih.gov/articles/PMC11352909/](https://pmc.ncbi.nlm.nih.gov/articles/PMC11352909/)  
  * [https://www.nature.com/articles/s41467-024-50526-3\#Sec11](https://www.nature.com/articles/s41467-024-50526-3#Sec11)  
  * [https://pubmed.ncbi.nlm.nih.gov/39528729/](https://pubmed.ncbi.nlm.nih.gov/39528729/)  
  * [https://pmc.ncbi.nlm.nih.gov/articles/PMC3460417/](https://pmc.ncbi.nlm.nih.gov/articles/PMC3460417/) (dehydroxy paper)  
* TET proteins convert 5mC to 5hmC.

**RESULTS**

1. **Bap1 loss results in increased MeCP2 binding in mature neurons, with no changes in binding in young neurons.**  
* MeCP2 volcano plots  
  * P12 \- no DMRs, no apparent clustering in PCA  
  * Adult \- thousands of DMRs, clear clustering by condition in PCA  
  * P12 vs adult mut \- adult is more  
  * P12 vs adult ctrl \- adult is more  
* Seems to be binding active promoters mostly.  
2. **Increase in H2AK119Ub after Bap1 loss leads to a global 5mC increase and 5hmC decrease trend at gene bodies.**  
   * Volcano plot with top genes  
   * 5mC vs 5hmC quadrant plot showing concordant and discordant genes  
   * Asymmetric direction of methylation changes bar plot  
   * Effect size plot  
   * Active promoters are increasingly methylated while repressed promoters are losing methylation.  
     * Methylation vs chromatin state bar plot  
   * Heatmap shows that MeCP2 down sites are specific to polycomb targets.  
3. **Most 5mC and 5hmC DMRs occur in gene bodies, CpG Shores, and CpG shelves.**  
   * DMR heatmap per context  
   * Sex differences are observed for both 5mC and 5hmC, regardless of mutant or control status.  
     * PCA  
4. **In neurons, MeCP2 acts as a reader of chromatin state shift, rather than solely changes in methylation.**  
* Mecp2 peak overlap at differentially methylated regions bar graph  
* \# sig of DMRs in modCs overall and what genomic elements they are found at.  
* \# sig DMRS 5mC vs 5hmC \+ what genomic elements they are found at.  
* 5mC up and 5hmC down at the same loci.  
  * Show bigwig example (s) of top loci where this is true.  
* Ub vs mecp2 log2 plot  
5. **Exclusively at genes involved in neuronal structure and development, Bap1 loss leads to relocation of MeCP2 from euchromatic to heterochromatic loci.**   
   * The figure with a lot of log2 plots  
6. 

**DISCUSSION**

**MATERIALS & METHODS**

**RESOURCE AVAILABILITY**  
***Lead contact:*** Further information and requests for resources and/or reagents should be  
directed to the Lead Contact, Cole Ferguson (cferguson@health.ucsd.edu).  
***Materials availability:*** Any other non-commercially available reagent or input biomaterial  
generated or used during this study will be made freely available to interested researchers.  
***Data and code availability:***  
Raw and analyzed data is available under the NCBI BioProject PRJNA1150596. All original code regarding the CUT\&RUN pipeline is available at: [https://github.com/FergusonLab/CUT-RUN-Pipeline](https://github.com/FergusonLab/CUT-RUN-Pipeline). All original code regarding downstream analyses for methylation can be found at \*ZAKIR’S GITHUB\*. Requests to access code and pipeline information regarding modality XPLR can be done by contacting Biomodal. Any additional information required to reanalyze the data reported in this paper is  
available from the lead contact upon request.

**EXPERIMENTAL MODELS**  
***Mice***

* Mouse model:  
  * Strains obtained from Jackson Labs:  
    * Bap1 f/f \#1031874: Bap1 double floxed on both alleles  
    * Bap1 \+/+; Math1-cre \#011104: hemizygous cre expression under control of Math1 enhancer element  
  * Bap1 f/f crossed with Bap1 f/f; Math1 cre  
  * Genotyped mice using primers designed to detect cre  
    

Animals were cared for in accordance with NIH guidelines. All experimental methods were  
approved by the UCSD Institutional Committee on the Use and Care of Animals under the  
protocol number S20121. Sex matched-littermate controls on the Jackson Labs C57BL/6  
background were used in all experiments, and the ages at which animals were used is reported  
in figure legends. Animals were housed in a 12:12 light:dark cycle. Biological replicates were  
sex-matched littermates, and both male and female animals were examined in the course of the  
major experiments. None of the main findings in this work varied by sex in the tissues we  
examined.

**METHOD DETAILS**  
**Cerebellar Nuclei Isolation**  
Mice were deeply anesthetized with isofluorane before decapitation and extraction of the  
while cerebellum, liver and/or kidney. The liver and kidney were further dissected, before  
choosing the identical lobe in the case of the liver and pole in the case of the kidney, inputting \~100 mg of tissue. Tissue was finely minced in lysis buffer (Nuclei EZ Lysis Buffer supplemented with 1x Halt combined protease and phosphatase Inhibitor cocktail, 10 mM sodium butyrate and 1.5 mM iodoacetamide) prior to mechanical homogenization with Dounce A pestle in 2 mL of Lysis Buffer. This suspension was incubated on ice for 5 minutes before centrifugation at 500 x g for 5 minutes, with acceleration and deceleration on the lowest setting, at 4 degrees. The supernatant was aspirated and nuclei were resuspended in EZ lysis buffer and incubated on ice for an additional 5 minutes prior to centrifugation. Examination of nuclei under DIC demonstrated abundant intact nuclei for all the tissues examined. Nuclei were resuspended in CUT\&RUN wash buffer and filtered through a 40 µm mesh. This process routinely yielded \~6 million nuclei from a single mouse cerebellum, liver and kidney samples.

**MeCP2 CUT\&RUN**  
For CUT\&RUN performed under native conditions, nuclei were resuspended in wash  
buffer (20 mM HEPES pH 7.5, 150 mM NaCl, 0.5 mM spermidine supplemented with Roche  
complete EDTA-free protease inhibitor tablet). Pelleted nuclei were resuspended twice in wash  
buffer and subjected to centrifugation at 500 x g for 5 minutes. Nuclei were resuspended in 1  
mL, filtered by gravity through a 40 µm mesh, and counted with a Countess automated cell  
counter (Invitrogen). In general, between 250,000-500,000 nuclei were used in each CUT\&RUN  
experiment, and the cellular input was kept the same across all replicates. Nuclei were bound to  
Concanavalin A-coated magnetic beads. Following bead binding, the supernatant was discarded, and the bead-nuclei slurry was resuspended in 75 uL of antibody buffer (wash buffer supplemented with 0.005% digitonin and 2 mM EDTA) containing primary antibodies at a concentration of 1:50. Antibody incubation was carried out overnight at 4°C on a nutator.  
Samples were washed with cold cell permeabilization buffer (Wash buffer with 0.005% digitonin)  
to remove unbound antibody and then incubated with 1x Protein A/G-MNase in cell  
permeabilization buffer for 15 minutes at room temperature on a shaker. Samples were washed  
with cold cell permeabilization buffer and then resuspended in 50 µL cell permeabilization buffer  
containing 2 mM CaCl2. Samples were digested for 2 hours at 4°C on a nutator. The digestion  
reactions were quenched with 33 µL of STOP buffer (340 mM NaCl, 20 mM EDTA pH 8, 4 mM  
EGTA pH 7.7, RNase A 0.05 mg/mL, Glycogen 0.05 mg/mL). Samples were then incubated for  
10 minutes at 37°C to release digested chromatin fragments. The supernatant containing  
CUT\&RUN fragments was transferred to fresh tubes and DNA was purified via column  
purification.

**Library preparation and sequencing**  
Samples were pooled for sequencing 100 bp reads in paired-end  
configuration on an Illumina NovaSeqX platform at the UCSD Institute for Genomic Medicine.  
See further details below regarding specific library preparation.  
***MeCP2 CUT\&RUN***  
Library preparation for CUT\&RUN experiments was performed according to the  
manufacturer’s specifications using SPRIselect beads. The concentration of DNA yielded from  
CUT\&RUN was measured on a Qubit device and between 25 and 50 ng of input DNA was used  
for library preparation, with the amount of input DNA kept constant for all samples using the  
same antibody. 14 PCR cycles were used in the final amplification with Illumina-compatible  
adaptors, universal i5 primer and barcoded i7 primers. After measuring the concentration of the  
amplified product on Nanodrop, all samples were diluted to the same concentration for a given  
antibody, ranging from 10 ng/uL to 40 ng/uL. Tapestation was performed to determine the  
precise concentration and fragment distribution range and to ensure there were not significant  
amount adaptor dimers.   
\>25 mil reads  
***DNA Methylation (5mC and 5hmC)***  
400-500 mil reads

**CUT\&RUN Data Processing and Analysis**  
***Quality Assessment and Sequencing Alignment***  
We used FastQC to evaluate read quality, including base quality scores and adapter  
contamination. Compressed sequenced reads in fastq.gz format were first trimmed using  
Trimmomatic v.0.39111. We then aligned samples to the mm10 genome using Bowtie2 v.2.3.4.1 and sorted BAM files, removing unmapped fragments using SAMtools v.1.14. During alignment, we set a minimum Phred score of 33 and employed the \--dovetail option. We then converted alignment files to .bigwig files using the bamCoverge tool in deepTools v.3.5.062 with 50 bp bins. Duplicate reads were not removed.   
***Peak Calling***  
Peak calling was performed using SEACR (Sparse Enrichment Analysis for CUT\&RUN) v.1.154 and MACS2 (Model-Based Analysis of ChIP-seq) in broad and narrow formats. Samples of different replicates collected in the same experiment were normalized to each other to correct for different loading between samples. .bigwig and .bed files were examined and visualized in Integrated Genome Viewer. We obtained regional overlap between different samples using the ‘intersect’ functionality of BEDTools v.2.29.2.  
***Normalization of Data***  
After parsing the genome into 50 bp bins, we identified local maxima and generated a  
blacklist of 271 intergenic peaks representing spurious alignment artifacts. We quantified the  
height of the leftover authentic peaks and recorded the value of the height for the 99th percentile peak recorded. The ratio of the values of the 99th percentile peaks between samples served as a scaling factor by which the height of every bin in one sample could be multiplied to normalize samples to one another.  
***DiffBind Analysis***  
Differential peak analysis was performed on SEACR-called peaks using DiffBind (v.  
3.14) with default parameters (summit width \= 400 bp). Significantly differential peaks were  
identified using DESeq2 within Diffbind and only peaks with FDR\<0.05 were used for further  
analysis. Functional annotation of significant peaks was performed using the “annotatePeak”  
function of ChIPseeker (v.1.4) using the TxDb.Mmusculus.UCSC.mm10.knownGene (v. 3.10)  
Database.  
***Gene Ontology (GO) Analysis***  
A list of the chromosomal coordinates of all known genes in the mm10 genome along with  
their Refseq IDs was obtained from the UCSC Genome Browser. The closest RefSeq genes  
within 10 kb of each .bed file locus was obtained using BEDTools v.2.29.2. The gene IDs left  
after this filtering was the gene list used as input to ShinyGO v.0.771 for GO analysis with an  
FDR cutoff of 0.05. No background was inputted.

**Generation of Figures**  
All plots were generated in R 4.5.2 and modified in Adobe Illustrator 2021\.

* Wet lab:  
  * Nuclei isolation from cerebellum  
  * MeCP2 C\&R \- epicypher  
  * Biomodal evoC kit (cite their paper)  
    * Experimental set up for both P12 and adult  
    * NovaSeqX Paired-end sequencing to 400-500 million reads  
* Computational:  
  * C\&R pipeline (S1b in AP’s paper)  
  * Roman normalization (S1c in AP’s paper)  
  * MACS2 Broad \- peak calling  
    * Chosen due to broad distribution of MeCP2 binding to CpG sites  
    * Peak curation method (mention this at all?)  
  * DiffBind for differential MeCP2 binding analysis  
  * Biomodal downstream pipeline “Modality” for DMR calling and quality metrics  
  * Zakir’s pipelines

**ACKNOWLEDGMENTS**  
We would like to thank the following individuals and teams for their support on this project:

* Dr. Vicky L. Brandt for many insightful discussions and a commitment to substance and clarity.   
* Ferguson lab members for their commitment and work on this project.  
* Dr. Carmen Guarco, Thao Huynh, Craig Fishman, and Dr. Padmanabhan Srinivasan from the Biomodal team for their thoughtful feedback and support on the use of the evoC kit in this study.  
* Dr. Anjana Rao and Dr. Zhen Dong at La Jolla Institute for Immunology for their guidance and feedback on this project.  
* Dr. Kristen Jepsen and the UC San Diego IGM Genomics Center team for data generated utilizing an Illumina NovaSeq X Plus that was purchased with funding from a National Institutes of Health SIG grant (\#S10 OD026929).

