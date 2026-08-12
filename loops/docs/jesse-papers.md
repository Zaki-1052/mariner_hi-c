# atlas

OK, I've read the core of the paper (Zhou, Wu, Liu et al. 2025 - "Human Body Single-Cell Atlas of 3D Genome Organization and DNA Methylation") and cross-referenced it with your biomodal and Hi-C work. Here's what's relevant:

---

## High relevance

### 1. Methylation and 3D genome are partially independent (Fig 5I, Fig 6)

This is probably the most important finding for your project. The paper shows that **DNA methylation and chromatin looping don't always co-vary across cell types** -- in some lineages (especially those with PMDs), there's weak correlation between DMRs and differential loops (Fig 5I). They go further in Fig 6 to show that mCG-based clustering and chromatin contact-based clustering can assign cells to *different subtypes* (skeletal muscle, placenta examples).

**Why this matters for you:** Your key-takeaways.md question #5 asks "How does the methylation story connect to the loop story?" This paper provides a framework: they may be **parallel consequences of the same upstream perturbation** (H2AK119ub accumulation) rather than one causing the other. Your discordant gene analysis already showed that loop-anchored genes behave differently (3x enrichment at loops, opposite K119ub direction). The paper's finding that methylation may reflect **maturation history** while 3D genome reflects **current cell state** is an interesting lens for your BAP1-KO model.

### 2. Non-CG methylation is pervasive, not neuron-specific (Fig 3)

Major finding: **mCH exists at low levels across nearly all cell types**, not just neurons. They show:
- Enrichment in specific trinucleotide contexts (CAC, CTC most prominent in non-neuronal cells)
- mCH patterns are **non-random and cell-type-specific** -- can distinguish cell types even in non-neuronal tissues
- mCH is **depleted at cCREs** across all lineages, suggesting functional relevance

**Why this matters for biomodal/non_cg_analysis/:** Your non-CG analysis comparing EM-seq vs evoC in cerebellum is well-placed. The paper validates that your CHG/CHH signals aren't just noise. But it also suggests you should look at **trinucleotide context specificity** -- your `non_cg_methylation.py` treats CHG and CHH as aggregate categories, but the paper shows that within those, specific trinucleotides (especially CA-containing) carry the most signal. Your `ca_filter.py` seems to already be going in this direction.

### 3. DMRs enriched at loop anchors + motif analysis (Fig 5H)

They explicitly tested motif enrichment at loop-associated DMRs versus non-loop DMRs and found lineage-specific TFs. **CTCF** was enriched at loop anchors in all major types. They also found **IRF2 and ISRE** in immune lineages and **GATA factors** in trophoblasts.

**Why this matters:** You could directly replicate this analysis -- overlay your biomodal DMRs onto your mariner differential loop anchors. Your concordance analysis already found that discordant genes (mC down/hmC up) are 3x enriched at loops. The paper's approach of k-means clustering loop-DMR pairs (Fig 5H) and then running motif enrichment per cluster could reveal which TFs are driving the differential behavior at loop-anchored vs non-loop loci in your BAP1-KO.

### 4. Methylation compartments are finer-grained than Hi-C compartments (Fig 5C-E)

They define "DNA methylation compartments" from 10kb-resolution mCG data that map to, but subdivide, the traditional A/B compartments from Hi-C (100kb). The partially methylated compartment consistently maps to B compartment, but provides more granularity.

**Why this matters:** Your compartment/TAD analysis (`tad_analysis/`) uses HOMER's PC1 at relatively coarse resolution. The paper argues that methylation-based compartmentalization is **higher resolution** and can reveal fine-grained compartment transitions that Hi-C misses. Your biomodal DUET data, with its 6bp resolution methylation, could in principle provide this kind of fine-grained compartment annotation -- a complementary view to your Hi-C compartment analysis.

---

## Moderate relevance

### 5. 283,606 differential chromatin loops identified (Fig 4E)

They found differential loops are **shared among related cell types** in clusters (fibroblasts, epithelial, neurons, hematopoietic). Loops in A-compartment regions are more lineage-specific. They used ANOVA-based F-statistics with FDR<0.05 for differential loop calling -- different from your edgeR quasi-likelihood GLM but conceptually similar.

**Useful as:** A reference for loop variability baseline. If any of the 35 cell types they profiled are relevant to cerebellum (neurons especially), you could compare your BAP1-KO differential loops against their normal neuronal loop catalog.

### 6. Compartment-dominant vs domain-dominant phenotypes (Fig 4A-C)

Neurons are the most domain-dominant (stronger short-range interactions, weaker compartmentalization). Hematopoietic lineages are most compartment-dominant. This is consistent with your observation of Polycomb loop collapse in BAP1-KO -- if BAP1 loss disrupts long-range Polycomb loops, it could shift cells from domain-dominant toward a more mixed phenotype.

### 7. mCH at gene bodies anticorrelated with expression in neurons (Fig 3H)

In neurons specifically, mCH in gene bodies is negatively correlated with expression, and this relationship **varies by lineage**. In endocrine pancreas, the relationship is different (anticorrelated downstream of TSS but positively correlated upstream).

**Relevant because:** Your biomodal data is from cerebellum (high 5hmC, high mCH). The paper's nuanced treatment of how mCH-expression relationships vary by tissue context is important for interpreting your non-CG findings. The paper also notes that **MeCP2** binding is influenced by mCH levels -- and your concordance analysis already tracks MeCP2.

---

## Lower relevance but worth noting

- **PMDs on a continuum** (Fig 2C-E): PMDs aren't binary but vary in degree across cell types. The "partially methylated compartment" concept could help interpret your biomodal global methylation patterns.
- **APA methodology** (Methods, Fig S20D): They do APA on single-cell differential loops at 10kb resolution with +/-100kb flanking -- validates your apa_analysis.R approach.
- **Human, not mouse**: The paper is all human (hg38). Your data is mouse (mm10). Direct comparisons of specific loci won't work, but the conceptual framework and methodological approaches transfer.

---

## Bottom line

The paper's central insight -- that **methylation and 3D genome organization can diverge** and may reflect different timescales of cellular identity -- is directly relevant to your core question of how BAP1 loss affects both modalities simultaneously. It gives you a framework to argue that your coordinated mC/hmC changes and your differential loop changes may stem from the same upstream cause (H2AK119ub accumulation, chromatin compaction) while operating through partially independent downstream mechanisms.

The non-CG findings also validate your `non_cg_analysis/` work and suggest looking at trinucleotide-level context specificity rather than just aggregate CHG/CHH.

Worth citing in your eventual writeup, especially Fig 5I (loop-DMR correlation variability), Fig 6 (methylation vs 3D genome discrepancies), and the non-CG discussion.

---
## p2

Here's what the paper has on neurons specifically:

**No cerebellum, no granule neurons.** The only brain tissue is **primary motor cortex** (2 donors, from the BRAIN Initiative Cell Census Network). The neuronal subtypes are all cortical:
- **Neu Exc**: IT L2/3 (5 subtypes), IT L4, IT L5, IT L6, CT L6, ET L5/6, NP
- **Neu Inh**: PVALB (2 subtypes), SST (3 subtypes), VIP (2 subtypes), LAMP5

So no granule neurons, no Purkinje neurons, no cerebellum at all. That said, there are several neuron-specific findings worth knowing:

### 3D genome: Neurons are uniquely domain-dominant

- Neurons have the **strongest short-range interactions and weakest compartmentalization** of any cell type (Fig 4A-C, Fig S19). They're at one extreme of a "compartment-dominant ↔ domain-dominant" spectrum.
- They have the **most and longest chromatin loops** (Fig 4D, Fig S20A).
- The paper speculates this could relate to **cohesin-mediated loop extrusion** -- neurons may have more active loop extrusion and weaker heterochromatin, reducing compartment segregation.

**Relevance to your Hi-C:** If BAP1 loss causes chromatin compaction and H2AK119ub spreading, it could shift neurons *away* from this domain-dominant phenotype toward something more compartment-dominant (stronger A/B segregation, weaker loops). That would be consistent with your observation of Polycomb loop collapse -- losing long-range loops while strengthening local compartment interactions.

### mCH: Neurons are the outlier in every analysis

- Neurons have **10-50x more mCH** than any other cell type (Fig 3A). The paper confirms this is real (not bisulfite artifact) via lambda phage spike-in controls (Fig S14A).
- **CAA and CTC are the dominant contexts** in neurons (Fig S14B), but non-neuronal cells also show low-level mCH enriched in these same contexts.
- **mCH distinguishes neuronal subtypes** cleanly -- IT L2/3, SST, PVALB, LAMP5 all separate by mCH alone (Fig S17).
- **mCH in gene bodies is negatively correlated with expression** in neurons (Fig 3H, Fig S18B), consistent with the MeCP2-mediated repression model. But the paper notes this relationship **varies across subtypes** even within excitatory neurons.

**Relevance to biomodal:** Your cerebellum data should show even higher mCH than these cortical neurons (Purkinje neurons have the highest documented mCH in any mammalian cell type). The paper's finding that mCH-expression relationships vary by neuronal subtype is important -- you can't assume cortical patterns extrapolate to cerebellar granule neurons.

### The critical neuron exception: Methylation ≠ compartments

This is the most important neuron finding for your project:

- In **most cell types**, mCG positively correlates with compartment score (high mCG = A compartment = active). In neurons, the correlation is **negative** (Fig S21A-B): Neu Inh r = -0.286, Neu Exc r = -0.249.
- Similarly, **mCH negatively correlates with compartment score only in neurons** (Fig 5B) -- neurons have enriched mCH in B-compartment regions.
- The paper explicitly states that neurons and glia are **excluded** from the compartment score vs PMD density comparison (Fig 5E) because the relationship breaks down.

**Why this matters:** In your BAP1-KO cerebellum, you're trying to connect methylation changes to 3D genome changes. This paper is saying that **in neurons specifically**, the normal rules about methylation tracking with compartments don't apply. The methylation-3D genome relationship in neurons operates differently from every other cell type. This makes it even harder to predict how BAP1-mediated chromatin compaction would affect methylation in your system -- and potentially even more interesting that you see coordinated changes at all.

### DMRs at loops: Neurons are split

From Fig S22:
- **Neu Inh** shows moderate-strong enrichment of DMRs at differential loop anchors
- **Neu Exc** shows lower enrichment
- In the loop-DMR correlation analysis (Fig S23), Neu Inh shows one of the **strongest negative correlations** -- loop gains correspond to methylation loss at those loci

This is consistent with your finding that discordant genes (mC down) are 3x enriched at loop anchors. In neurons, the loop-methylation relationship appears to be: **more looping = less methylation** at those sites.

### Bottom line

The paper has rich cortical neuron data but nothing on cerebellum. The most citable finding for your work is that **neurons are exceptional** -- they break the methylation-compartment relationship seen in every other cell type, they're the most domain-dominant lineage, and they show the strongest anticorrelation between loops and DMRs. This reinforces that your BAP1-KO cerebellum system, where both modalities change coordinately, is operating in a unique biological context where the normal rules may not apply.


---


# paper 3

This paper is **much** more directly relevant than the human atlas. Liu, Zeng, Zhou et al. (2023) *Nature* 624, 366-377 -- "Single-cell DNA methylome and 3D multi-omic atlas of the adult mouse brain." Same lab (Ecker/Dixon), same species (mouse, mm10), and **cerebellum is included**.

Here's what matters for your work:

---

## Cerebellum is explicitly covered

From Extended Data Fig. 1a: cerebellum dissection regions include **CB-1 through CB-5** (slices 15-18). From Extended Data Fig. 4: cerebellum (CB) gets its own panel with **30,424 nuclei** profiled, showing cell subclasses and dissection regions. The cell taxonomy includes:

- **10. CB Glut** -- cerebellar granule neurons (glutamatergic)
- **23. CB GABA** -- cerebellar GABAergic neurons (Purkinje cells, interneurons)

These are your cell types. The paper has single-cell methylome (mCG + mCH) and 3D chromatin conformation data for them.

---

## Directly relevant findings

### 1. mCH as the primary gene body regulatory mark in neurons

The entire paper's cell taxonomy is built on **gene body mCH fraction**, not mCG. They show (Extended Data Fig. 5a) that normalized gene body mCH is anticorrelated with RNA expression at the cell-group level -- this is the primary epigenetic readout for neuronal gene activity across the whole brain. Immediate early genes (*Fos*, *Egr1*, *Arc*, *Bdnf*, *Nr4a2*) and neurotransmitter genes (*Slc17a7*, *Gad1*) all show this pattern.

**Why this matters:** Your biomodal data captures both mC and hmC. The paper's framework says mCH in gene bodies is *the* mark for neuronal cell identity. Your finding that gene body methylation changes (not promoter) dominate in BAP1-KO fits exactly -- you're perturbing the methylation mark that defines neuronal identity in this atlas.

### 2. 2.6 million DMRs as a reference set

They identified **2.6 million cell-type-specific CG-DMRs** across the whole brain (Fig 2). These DMRs:
- Are enriched near TSS (distal cCREs represent a large repertoire)
- 91.5% are >2kb from TSS -- matching your finding that gene bodies (not promoters) are where the action is
- Hypo-methylated DMRs predict cell-type specific regulatory element activity when cross-referenced with snATAC-seq (Extended Data Fig. 6c-d)

**Why this matters:** You could directly compare your BAP1-KO DMRs against their normal cerebellar granule cell DMR catalog to see whether BAP1 loss preferentially disrupts cerebellar-specific vs brain-wide regulatory elements.

### 3. Spatial methylation gradients across brain regions

Fig 3 and Extended Data Fig. 7a-c show that **DMRs and mCH levels vary spatially** across brain axes (anterior-posterior, dorsal-ventral, medial-lateral). Different brain regions show distinct DMR methylation landscapes even within the same cell type.

**Why this matters:** Your BAP1-KO is in whole cerebellum. This paper provides the normal spatial baseline for cerebellar methylation patterns.

### 4. TAD boundaries correlate with mCH (Fig 4)

This is one of the most directly relevant findings:

- **TAD boundaries in neuronal genes are strongly correlated with the transcript body mCH fraction** (Fig 4). They calculated the probability of boundary at TSS and TTS for each gene.
- **Genes within boundary regions are important for neuronal function**: they found genes linked to **essential tremor** (*Lingo1*/*Lingo2*) and **Parkinson's disease** at the boundary between high-mCH gene body domains and low-mCH intergenic regions.
- The boundary probability is **negatively correlated with mCH** (PCC < -0.65, FDR < 0.001, permutation test).

**Why this matters for your Hi-C:** Your TAD analysis shows differential TADs in BAP1-KO. This paper provides the mechanistic connection -- TAD boundaries in neurons are defined by mCH gradients. If BAP1-KO disrupts methylation at gene bodies, you'd expect TAD boundary shifts. The paper explicitly cites **Tan et al. 2023 (ref 38) "Lifelong restructuring of 3D genome architecture in cerebellar granule cells"** -- showing they recognized cerebellar granule cells as a key system for this.

### 5. Chromatin conformation diversity at neuronal genes (Fig 4)

They analyzed variable chromatin interactions around neuronal genes:
- **ANOVA across cell subclasses** at 25kb resolution to identify highly variable interactions
- Four types: upstream, downstream, upstream-intragenic, downstream-intragenic
- Upstream-intragenic and downstream-intragenic interactions were **predominantly negatively correlated with gene body mCH** and positively correlated with gene body mCH fraction
- Upstream-downstream interactions were positively correlated with gene body methylation
- **Gene body domain boundaries** are a key structural feature -- "a potential explanation for the higher gene domain boundary probability observed"

**Why this matters:** This directly informs how to interpret your BAP1-KO loop changes. If BAP1 loss increases gene body methylation (your mC-up finding), this paper predicts it should:
- Strengthen intragenic TAD boundaries
- Weaken upstream-intragenic and downstream-intragenic interactions
- Potentially explain the "loop collapse inward to nearer TAD checkpoint" your mentor described

### 6. Gene Regulatory Network (GRN) connecting TFs, DMRs, and genes (Fig 5)

They built a multi-omic GRN for the whole mouse brain:
- **TF → DMR → target gene** connections using correlation between mCH fraction, DMR mCG, and gene expression
- Cell-type-specific TF importance via **PageRank analysis**
- Numerous TFs with mCH fractions correlated positively or negatively with the TF mCH fraction
- **CTCF motif enriched** at loop summit anchors across all cell types

**Why this matters:** You could query their GRN for cerebellar granule cell-specific TFs and check whether the genes disrupted in your BAP1-KO are downstream targets. Their data is publicly available.

### 7. Alternative splicing predicted by intragenic methylation + chromatin (Fig 6)

They found that intragenic DNA methylation and chromatin conformation patterns **predict alternative gene isoform expression** observed in whole-brain SMART-seq2 data. This is neuron-specific alternative splicing driven by epigenetic state.

**Why this matters:** If BAP1-KO changes gene body methylation, it could alter splicing patterns. This hasn't been explored in your data yet but is a testable prediction.

---

## Data availability -- this is a usable resource

Their data is available at:
- **NeMO Archive** (snmC-seq and snm3C-seq)
- **NCBI GEO** (snATAC-seq)
- **Interactive browser**: mousebrain.salk.edu
- **Code**: github.com/lhqing/wmb2023

You could download their cerebellar granule cell (CB Glut) methylome and chromatin conformation profiles and directly compare them to your BAP1-KO data as a "normal cerebellum" reference.

---

## Bottom line

This paper is arguably the **most important reference** for your entire project. It provides:
1. The normal cerebellar granule cell methylome and 3D genome at single-cell resolution (same species, same genome)
2. The mechanistic connection between gene body mCH and TAD boundaries in neurons
3. A GRN framework for interpreting regulatory consequences of methylation perturbation
4. Direct evidence that chromatin conformation in neuronal genes tracks with methylation state

The key citation for your writeup: their finding that TAD boundaries are defined by mCH gradients at gene bodies (Fig 4) directly supports the model that BAP1-KO methylation changes could drive the chromatin architecture changes you observe. This is stronger evidence for the methylation→3D genome direction of causality than anything in the Zhou et al. 2025 human atlas.