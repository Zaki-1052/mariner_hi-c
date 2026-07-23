MeCP2 controls neurological function by selectively repressing long genes through binding to methylated CA sites, and when this repression is lost in Rett syndrome, the resulting overexpression of these long, neuron-specific genes drives disease pathology.

## The Core Problem

Rett syndrome (RTT) is a severe neurological disorder caused by mutations in the _MECP2_ gene, which encodes a protein that binds methylated DNA. Despite years of study, no one had a clear model for _how_ MeCP2 actually regulates transcription in neurons. Previous studies found only subtle, inconsistent gene expression changes in MeCP2 mutant mice, and no unifying explanation had emerged for what ties those misregulated genes together.

## The Key Discovery: Gene Length Is the Unifying Feature

Gabel et al. systematically asked what upregulated genes in MeCP2 knockout brains have in common. The answer wasn't a shared pathway, histone mark, or promoter feature — it was **gene length**. Genes that go up when MeCP2 is lost are dramatically longer than the genome average. This relationship is continuous: the longer the gene, the more it's upregulated in the knockout. The effect is small per gene but widespread across the genome, and it was reproducible across multiple brain regions, labs, and measurement platforms (microarray, RNA-seq, qPCR, nCounter). Conversely, when MeCP2 is overexpressed (as in MeCP2 duplication syndrome), long genes are specifically _downregulated_, confirming that MeCP2 directly represses transcription in a length-dependent way.

Critically, this pattern is specific to MeCP2 disruption. Sixteen other mouse models of neurological disease did not show the same long-gene misregulation, ruling out the possibility that it's a generic consequence of neuronal dysfunction.

## The Mechanism: mCA Methylation as the Binding Platform

The paper then asks _why_ long genes are specifically affected. MeCP2 was originally identified for binding methylated CG dinucleotides (mCG), but neurons also have a unique type of DNA methylation: methylated cytosine followed by adenine (mCA), deposited by the enzyme Dnmt3a during postnatal brain development.

Using electrophoretic mobility shift assays, the authors showed MeCP2 binds mCA and mCG with high affinity, but does _not_ bind well to mCT, mCC, or hydroxymethylated CG (hmCG). ChIP-seq confirmed that in living neurons, MeCP2 binding density within gene bodies correlates with mCA levels, especially in genes longer than 100 kb.

The mechanistic logic connects as follows: longer genes physically contain more mCA sites because mCA density increases with gene length. More mCA means more MeCP2 binding across the gene body. More MeCP2 binding means more transcriptional repression. So when MeCP2 is removed, the longest, most mCA-rich genes lose the most repression and become the most overexpressed. Confirming this, long genes with _low_ mCA levels are largely unaffected in the knockout, while long genes with _high_ mCA are strongly upregulated.

The Dnmt3a conditional knockout experiment provides an elegant orthogonal validation: deleting the enzyme that deposits mCA (without affecting MeCP2 itself) produces the same length- and mCA-dependent gene upregulation pattern as deleting MeCP2. This confirms mCA is the critical substrate for MeCP2-mediated repression.

## Why This Causes Neurological Disease Specifically

The paper addresses why disrupting a broadly expressed protein causes specifically _neurological_ dysfunction. Long genes as a population are enriched for neuronal functions (synaptic proteins, ion channels, axon guidance molecules) and are preferentially expressed in the brain relative to other tissues. The 466 genes consistently repressed by MeCP2 are enriched for annotations like post-synaptic density, axonogenesis, and voltage-gated cation channel activity. So neurons are uniquely vulnerable because they disproportionately rely on long genes for their specialized functions, and they use the mCA/MeCP2 system to keep those genes in check.

## Disease Correlation and Therapeutic Implications

The severity of long-gene misregulation tracks with disease progression (worse at 9 weeks than 4 weeks in mice) and with mutation severity (more severe mutations cause greater misregulation). The pattern is also present in human RTT brain tissue and human MECP2-null neurons derived from ES cells.

As a proof of concept for therapeutic intervention, the authors showed that topotecan (a topoisomerase inhibitor known to preferentially suppress long gene transcription) can dose-dependently reverse the overexpression of long genes in MeCP2-knockdown neurons. The concentration that best normalizes gene expression (50 nM) also partially rescues the reduced ribosomal RNA content that serves as a marker of overall cell health in MeCP2-deficient neurons.

## Connection to Broader Autism Biology

The paper also draws a link to fragile X syndrome (FXS), another autism-spectrum disorder. FMRP, the protein lost in FXS, represses _translation_ of mRNAs, and its target mRNAs are also encoded by exceptionally long genes that significantly overlap with MeCP2-repressed genes. This suggests that overactivation of long gene function — whether at the level of transcription (RTT) or translation (FXS) — may be a shared mechanism across neurodevelopmental disorders.