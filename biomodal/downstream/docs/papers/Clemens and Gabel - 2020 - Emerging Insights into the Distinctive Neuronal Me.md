## Summary

This review establishes that neurons possess a fundamentally different epigenetic landscape from other cell types, defined by massive accumulation of non-CG DNA methylation (mCH) that is read by the protein MeCP2 to control gene expression. Disruption of any component in the mCH → MeCP2 → NCoR repression axis causes neurodevelopmental disorders, including Rett syndrome, Tatton-Brown–Rahman syndrome, and autism.

---

## The Neuronal Methylome Is Unique

In most mammalian cells, DNA methylation occurs almost exclusively at CG dinucleotides. Neurons break this rule dramatically. After birth, the enzyme DNMT3A ramps up in neurons and deposits methylation at non-CG sites (called mCH, where H = A, C, or T — predominantly mCA). Because CH sites outnumber CG sites in the genome by roughly 50:1, even though each individual CH site is methylated at a much lower rate (around 1.5–3%), the total number of methylated CH sites in a neuron ends up rivaling or exceeding the total number of methylated CG sites (roughly 16–30 million mCH versus ~17 million mCG). Neurons also accumulate unusually high levels of hydroxymethylation (hmCG), created by TET enzymes oxidizing existing mCG marks. Together, these modifications create an epigenetic environment that exists nowhere else in the body.

## How mCH Patterns Are Established

mCH accumulates postnatally — peaking around 4–6 weeks in mice and taking up to 16 years to fully accumulate in humans — during the same window that neural circuits are being wired and refined. Two major forces shape where mCH lands across the genome.

**Gene expression acts as a shield.** Genes that are highly transcribed during early postnatal life physically block DNMT3A from binding and depositing mCH in their gene bodies. Lowly expressed genes accumulate high mCH, and moderately expressed genes get intermediate levels. Critically, these patterns set during early development persist into adulthood, meaning that a snapshot of early gene expression gets "frozen" into the methylation landscape.

**Chromosome architecture sets regional baselines.** At a larger scale, mCH levels correlate with topologically associating domains (TADs) — megabase-sized chunks of chromatin that fold together in 3D space. Each TAD tends to have a consistent "set-point" of mCH, so genes and enhancers within a high-mCH TAD carry more methylation than equivalent elements in a low-mCH TAD. The molecular mechanism likely involves DNMT3A reading histone marks: its ADD domain binds unmethylated H3K4 (enriched in open euchromatin, absent at active promoters), while its PWWP domain binds H3K36me2 in broad euchromatic regions. Active gene bodies, marked by H3K36me3 instead, may recruit DNMT3A less efficiently.

## mCH as a Cell Identity Marker

One of the most striking features of mCH is its cell type specificity. Global mCH levels vary up to twofold between brain regions and ~1.5-fold between neuron subtypes within the same region. For instance, PV+ and SST+ inhibitory interneurons carry about 30% more mCH than VIP+ neurons, and deep-layer cortical excitatory neurons have 30–50% more than upper-layer ones. At individual gene loci, mCH patterns are even more cell type-distinctive than canonical mCG or chromatin accessibility, making single-neuron methylomes powerful enough to predict a neuron's brain region and laminar position. These patterns arise from cell type-specific gene expression during the postnatal window and then serve to maintain those identity programs in adulthood.

## MeCP2: The Reader

MeCP2 is the primary functional reader of mCH. It accumulates in neurons during postnatal development in parallel with mCH, reaching near-histone levels (~16 million molecules per nucleus). While originally identified as a CG methylation reader, MeCP2 binds mCA with high affinity — particularly mCAC, which happens to be the most common non-CG methylation trinucleotide. Engineered mice carrying MeCP2 that can bind mCG but not mCA still develop Rett-like phenotypes, demonstrating that CG binding alone is insufficient.

An important nuance: the massive conversion of mCG to hmCG in neurons (by TET enzymes) reduces MeCP2's affinity for CG sites, effectively shifting its functional binding toward mCH and hmCH sites. MeCP2 binds essentially everywhere in the genome, but shows modest enrichment at highly methylated regions and modest depletion at unmethylated regions — a low dynamic range that makes traditional "target gene" identification by ChIP-seq alone unreliable. Instead, the field has had to integrate methylation maps with transcriptomic data to identify functional relationships.

## Mechanism of Gene Repression

MeCP2's best-characterized partner is the NCoR co-repressor complex, which contains the histone deacetylase HDAC3. The current model works as follows: MeCP2 binds mCH and mCG throughout gene bodies and at enhancers, recruiting NCoR. This broadly suppresses histone acetylation, particularly at enhancer elements, which reduces their ability to activate target genes. Genes repressed by MeCP2 ("MeCP2-repressed genes") are enriched for mCH in their gene bodies, at associated enhancers, and throughout their surrounding TADs.

Importantly, the mechanism is _not_ a "speed bump" model where MeCP2 physically blocks RNA polymerase elongation — recent GRO-seq and intronic RNA-seq data show that what actually changes is transcription _initiation_, not processivity. The current explanation invokes 3D chromatin contacts: MeCP2 bound to mCH in gene bodies and enhancers can loop to promoters and repress initiation through NCoR-mediated deacetylation. Intragenic enhancers are more susceptible to this repression than extragenic ones, which explains why gene body methylation originally appeared so important.

However, one study found that disabling NCoR's HDAC3 enzymatic activity didn't rescue MeCP2 overexpression toxicity, suggesting NCoR may contribute to repression through mechanisms beyond deacetylation alone. MeCP2 also has emerging biophysical properties — it undergoes liquid-liquid phase separation with chromatin in vitro, and Rett-causing mutations reduce this condensate-forming ability. How phase separation connects to the epigenomic and transcriptomic consequences of MeCP2 disruption remains an open question.

## Disease Implications

The mCH–MeCP2–NCoR axis is exquisitely dose-sensitive: too much or too little MeCP2 both cause severe neurological dysfunction. This sensitivity means the pathway is vulnerable to disruption at multiple points.

**Rett syndrome** results from loss-of-function MECP2 mutations (the reader). **MeCP2 duplication syndrome** results from MECP2 overexpression (too much reader). **Tatton-Brown–Rahman syndrome** results from heterozygous DNMT3A mutations (reduced writer) — and mouse models show that even 50% DNMT3A loss causes ~50% global mCH reduction and partially recapitulates MeCP2 loss-of-function effects on enhancer acetylation and gene expression. Mutations in TBL1XR1, a component of the NCoR complex, cause additional neurodevelopmental phenotypes — and some missense mutations specifically disrupt NCoR's interaction with MeCP2.

These disorders share overlapping features (intellectual disability, autism-spectrum behaviors, seizures) but also have unique phenotypes reflecting the non-overlapping roles of each protein. A particularly hopeful finding is that reintroducing MeCP2 in adult knockout mice dramatically reverses symptoms, because the mCH patterns it needs to read were laid down correctly during development. This has motivated gene therapy approaches for Rett syndrome, though disorders caused by DNMT3A mutations may be harder to reverse if critical methylation patterns were never properly established.

## Functional Model

The review proposes an elegant developmental model: during the postnatal critical period, actively expressed genes resist mCH accumulation and remain accessible, while inactive genes accumulate mCH and become locked into repression by MeCP2. This effectively converts a transient gene expression state into a durable epigenetic memory, closing an "epigenomic critical period" that stabilizes circuit function — analogous to how extracellular matrix buildup closes critical periods of synaptic plasticity. Activity-dependent MeCP2 phosphorylation and stimulus-driven methylation changes then allow more limited, dynamic gene regulation in adulthood on top of this stable baseline.