# Annotated Bibliography: The Literature Story

This bibliography is organized to tell the story of the field chronologically and conceptually, starting with foundational discoveries about 5hmC in brain, moving through TET enzyme mechanisms, then BAP1/Polycomb biology, and finally the limited neural-specific work. This order reflects how you might explain the background to your grad mentor before presenting your findings.

---

## Part 1: Discovery of 5hmC and Its Unique Abundance in Brain

### 1. Kriaucionis & Heintz (2009) — *Science*
**"The nuclear DNA base 5-hydroxymethylcytosine is present in Purkinje neurons and the brain"**

https://pubmed.ncbi.nlm.nih.gov/19372393/

*The foundational discovery.* First demonstration that 5hmC exists as a stable modification in mammalian DNA, with highest levels in cerebellar Purkinje neurons (~0.6% of all nucleotides, equivalent to ~40% of 5mC levels). Established the brain as a uniquely 5hmC-enriched tissue.

**Why it matters for your work:** Sets up why cerebellum is the ideal tissue to study 5mC/5hmC dynamics—any perturbation will be maximally detectable here.

---

### 2. Tahiliani et al. (2009) — *Science* [Rao Lab]
**"Conversion of 5-methylcytosine to 5-hydroxymethylcytosine in mammalian DNA by MLL partner TET1"**

https://pubmed.ncbi.nlm.nih.gov/19372391/

*Discovery of TET enzymes.* Identified TET1 as the enzyme that converts 5mC to 5hmC, establishing the enzymatic basis for active DNA demethylation. This paper launched the field of TET biology.

**Why it matters for your work:** Establishes the molecular machinery (TET enzymes) whose function may be impaired in your BAP1-KO model.

---

### 3. Szulwach et al. (2011) — *Nature Neuroscience*
**"5-hmC–mediated epigenetic dynamics during postnatal neurodevelopment and aging"**

https://www.nature.com/articles/nn.2959 | PMID: 21804537

*Developmental dynamics of brain 5hmC.* Genome-wide mapping in mouse hippocampus and cerebellum showing 4.22-fold increase in 5hmC from P7 to 6 weeks, with continued accumulation through aging. Demonstrated 5hmC enrichment in gene bodies of active genes and depletion at promoter CpG islands.

**Why it matters for your work:** Establishes the developmental acquisition pattern of 5hmC in cerebellum and the gene body distribution—exactly where your DMRs occur.

---

### 4. Mellén et al. (2012) — *Cell*
**"MeCP2 binds to 5hmC enriched within active genes and accessible chromatin in the nervous system"**

https://www.sciencedirect.com/science/article/pii/S0092867412014079 | https://pubmed.ncbi.nlm.nih.gov/23260135/

*Functional significance of brain 5hmC.* Showed that MeCP2 (mutated in Rett syndrome) binds both 5mC and 5hmC with similar affinity, but 5hmC marks accessible chromatin in neurons. Established that ~40% of modified cytosines in brain are 5hmC.

**Why it matters for your work:** Demonstrates that 5hmC has specific functions in neurons distinct from 5mC, making the 5mC↑/5hmC↓ shift in your data functionally meaningful.

---

### 5. Mellén et al. (2017) — *PNAS*
**"5-hydroxymethylcytosine accumulation in postmitotic neurons results in functional demethylation of expressed genes"**

https://www.pnas.org/doi/10.1073/pnas.1708044114

*5hmC as "functional demethylation."* Showed that 5hmC accumulation in gene bodies reduces MeCP2 occupancy by replacing high-affinity 5mCG sites with low-affinity 5hmCG sites. This explains why 5hmC in gene bodies promotes transcription even though the cytosine is still modified.

**Why it matters for your work:** Provides the mechanistic logic for why losing 5hmC (and gaining 5mC) at gene bodies would have transcriptional consequences.

---

## Part 2: TET Enzyme Mechanisms and Regulation

### 6. Ito et al. (2011) — *Nature* [Rao Lab]
**"Tet proteins can convert 5-methylcytosine to 5-formylcytosine and 5-carboxylcytosine"**

https://pubmed.ncbi.nlm.nih.gov/21778364/

*Complete oxidation pathway.* Demonstrated that TET enzymes catalyze sequential oxidation: 5mC → 5hmC → 5fC → 5caC, establishing the full active demethylation pathway.

**Why it matters for your work:** If TET access is blocked, only the first step fails—explaining why 5mC accumulates while 5hmC depletes (existing 5hmC still progresses downstream but isn't replenished).

---

### 7. Ko et al. (2015) — *Nature* [Rao Lab]
**"TET proteins and 5-methylcytosine oxidation in hematological cancers"**

https://pubmed.ncbi.nlm.nih.gov/26607761/

*TET loss causes malignancy.* Showed that deletion of multiple TET proteins leads to aggressive blood cancers with increased DNA damage. Established TET enzymes as tumor suppressors.

**Why it matters for your work:** Provides cancer context for why TET dysfunction matters; parallels between cancer and neurodegeneration pathways.

---

### 8. López-Moyado et al. (2019) — *PNAS* [Rao Lab]
**"Paradoxical association of TET loss of function with genome-wide DNA hypomethylation"**

https://pubmed.ncbi.nlm.nih.gov/31371502/

*The TET paradox—critical paper.* Showed TET deficiency causes dual methylation changes: hypermethylation in euchromatin (A compartment) but hypomethylation in heterochromatin (B compartment). Explained by DNMT3A redistribution from heterochromatin to euchromatin.

**Why it matters for your work:** Your hypermethylation at Active Promoter chromatin states mirrors the euchromatin hypermethylation pattern. This paper provides framework for understanding how methylation machinery redistributes when normal patterns are disrupted.

---

### 9. Stoyanova et al. (2021) — *eLife*
**"5-Hydroxymethylcytosine-mediated active demethylation is required for mammalian neuronal differentiation and function"**

https://elifesciences.org/articles/66973

*TET in Purkinje cell differentiation.* Deletion of Tet1/2/3 from postmitotic Purkinje cells prevents proper neuronal differentiation, directly demonstrating TET requirement in cerebellar neurons.

**Why it matters for your work:** Most directly relevant TET-brain paper. Shows TET function is essential in the exact cell type you're studying.

---

### 10. Joshi et al. (2022) — *Cellular and Molecular Life Sciences*
**"Mechanisms that regulate the activities of TET proteins"**

https://link.springer.com/article/10.1007/s00018-022-04396-x

*TET regulation review.* Comprehensive review of TET interaction partners, post-translational modifications, and recruitment mechanisms. Discusses how chromatin context affects TET access.

**Why it matters for your work:** Background on what factors control TET activity—relevant for hypothesizing how chromatin compaction might restrict TET function.

---

## Part 3: BAP1/PR-DUB and Chromatin Architecture

### 11. LaFave et al. (2015) — *Nature Medicine*
**"Loss of BAP1 function leads to EZH2-dependent transformation"**

https://www.nature.com/articles/nm.3947

*BAP1-EZH2 connection.* Showed BAP1 loss increases H3K27me3 and EZH2 expression; mesothelioma cells lacking BAP1 are sensitive to EZH2 inhibitors. First major paper connecting BAP1 to Polycomb repressive marks.

**Why it matters for your work:** Establishes BAP1's role in restraining Polycomb activity; relevant for understanding why H3K27me3 patterns might change.

---

### 12. Conway et al. (2021) — *Molecular Cell*
**"BAP1 enhances Polycomb repression by counteracting widespread H2AK119ub1 deposition and chromatin condensation"**

https://www.sciencedirect.com/science/article/pii/S1097276521005001 | https://pubmed.ncbi.nlm.nih.gov/34186021/

*Key mechanistic paper.* Demonstrated BAP1 constrains H2AK119ub throughout the genome; loss causes >75% of genome to become more compact, titrates PRC2 away from targets, and globally reduces transcription initiation. Identified PCGF3/5-PRC1 as the complexes spreading H2AK119ub.

**Why it matters for your work:** Directly establishes that BAP1 loss → chromatin compaction. This is the foundation for hypothesizing reduced TET access.

---

### 13. Bonnet et al. (2022) — *Genes & Development*
**"PR-DUB preserves Polycomb repression by preventing excessive H2AK118 mono-ubiquitylation"**

https://pubmed.ncbi.nlm.nih.gov/35961776/

*PR-DUB as rheostat—Drosophila.* Showed excessive H2Aub actually interferes with nucleosome stacking and chromatin fiber folding. Creates paradox: high H2Aub increases DNA accessibility despite Polycomb binding.

**Why it matters for your work:** Adds nuance—chromatin effects of BAP1 loss may not be simple compaction everywhere. The relationship is more like a rheostat than an on/off switch.

---

### 14. Thomas et al. (2023) — *Science Advances*
**"Structural basis of histone H2A lysine 119 deubiquitination by Polycomb repressive deubiquitinase BAP1/ASXL1"**

https://www.science.org/doi/10.1126/sciadv.adg9832

*PR-DUB structure.* Cryo-EM structure of human BAP1-ASXL1 on nucleosomes, showing how the complex is positioned for H2AK119Ub removal. Mapped >50 cancer-associated mutations to functional domains.

**Why it matters for your work:** Structural understanding of how BAP1 works; useful for understanding why loss has such widespread effects.

---

## Part 4: H2AK119ub and DNA Methylation Crosstalk

### 15. Weinberg et al. (2021) — *Nature Structural & Molecular Biology*
**"The histone mark H3K36me2 recruits DNMT3A and shapes the intergenic DNA methylation landscape"**

https://pubmed.ncbi.nlm.nih.gov/33432233/

*DNMT3A targeting mechanisms.* Showed how histone marks direct DNA methyltransferase localization. Relevant for understanding how chromatin state controls methylation patterns.

**Why it matters for your work:** Background on how DNMTs are recruited to specific genomic regions.

---

### 16. Gu et al. (2022) — *Nature*
**"DNMT3A and TET2 compete and cooperate to repress lineage-specific transcription factors in hematopoietic stem cells"**

https://pubmed.ncbi.nlm.nih.gov/35508662/

*DNMT3A-TET competition.* Demonstrated that DNMT3A and TET2 directly compete at the same loci; loss of one affects the other's function.

**Why it matters for your work:** Shows these pathways are interconnected—changes in one affect the other.

---

### 17. Chen et al. (2024) — *bioRxiv* (preprint)
**"Cancer-associated DNA Hypermethylation of Polycomb Targets Requires DNMT3A Dual Recognition of Histone H2AK119 Ubiquitination and the Nucleosome Acidic Patch"**

https://www.biorxiv.org/content/10.1101/2024.03.18.585588v1.full

*DNMT3A-H2AK119ub structural basis.* Cryo-EM showing DNMT3A1's UDR domain binds H2AK119ub, recruiting DNA methyltransferase to Polycomb-marked regions. Explains hypermethylation at Polycomb targets.

**Why it matters for your work:** If H2AK119ub accumulates (from BAP1 loss), DNMT3A could be recruited—potentially explaining hypermethylation independent of TET effects.

---

### 18. Zou et al. (2024) — *Nature*
**"RNA m5C oxidation by TET2 regulates chromatin state and leukaemogenesis"**

https://www.nature.com/articles/s41586-024-07969-x

*TET2→MBD6→BAP1 axis.* Discovered that TET2 oxidizes m5C on chromatin-associated RNA, which MBD6 recognizes to recruit BAP1. Direct molecular link between TET activity and H2AK119ub regulation.

**Why it matters for your work:** Shows TET and BAP1 are directly connected through an RNA-mediated pathway—these systems talk to each other.

---

## Part 5: BAP1 in Neural Contexts (Limited Literature)

### 19. Field et al. (2019) — *Clinical Cancer Research*
**"BAP1 loss is associated with DNA methylomic repatterning in highly aggressive Class 2 uveal melanomas"**

https://pmc.ncbi.nlm.nih.gov/articles/PMC6744995/

*BAP1 and methylation in cancer.* BAP1 knockdown in uveal melanoma induced methylomic repatterning with predominant hypermethylation at gene promoters.

**Why it matters for your work:** Only existing study of DNA methylation changes in BAP1 loss—but in cancer, not brain, and only 5mC (not 5hmC), and focused on promoters (not gene bodies).

---

### 20. Küry et al. (2022) — *American Journal of Human Genetics*
**"De novo mutations in BAP1 cause a severe neurodevelopmental disorder"** (Küry-Isidor Syndrome)

https://pubmed.ncbi.nlm.nih.gov/35120633/

*BAP1 in human neurodevelopment.* Described 11 individuals with de novo heterozygous missense BAP1 variants showing developmental delay, intellectual disability, seizures, and autism features.

**Why it matters for your work:** Only human genetic evidence linking BAP1 to brain function. Shows BAP1 is required for normal neurodevelopment.

---

### 21. Liao et al. (2024) — *Journal of Clinical Investigation*
**"BAP1 is required prenatally for differentiation and maintenance of postnatal murine enteric nervous system"**

https://pmc.ncbi.nlm.nih.gov/articles/PMC11060734/

*BAP1 in enteric neurons.* First evidence linking prenatal BAP1 loss to postnatal neurodegeneration (enteric nervous system). Shows BAP1 is required for neuronal differentiation and survival.

**Why it matters for your work:** Closest existing neural BAP1-KO study, though in ENS not CNS. Supports that BAP1 has critical neuronal functions.

---

### 22. Jiang et al. (2015) — *Brain*
**"Alteration in 5-hydroxymethylcytosine-mediated epigenetic regulation leads to Purkinje cell vulnerability in ATM deficiency"**

https://academic.oup.com/brain/article/138/12/3520/413106

*5hmC loss causes Purkinje degeneration.* ATM deficiency reduces 5hmC in cerebellar Purkinje cells, leading to neurodegeneration. Shows 5hmC is specifically protective in Purkinje neurons.

**Why it matters for your work:** Precedent for 5hmC loss causing cerebellar pathology—supports that your 5hmC decrease could have functional consequences.

---

## Part 6: Additional Context Papers

### 23. Ma et al. (2018) — *Epigenetics & Chromatin*
**"Distal regulatory elements identified by methylation and hydroxymethylation haplotype blocks from mouse brain"**

https://link.springer.com/article/10.1186/s13072-018-0248-3

*5hmC in gene bodies as distal regulators.* 5hmC CpG sites frequently have simultaneous 5mC modification; enriched in gene body regions of neuron development genes.

**Why it matters for your work:** Supports gene body as the relevant compartment for 5hmC biology in brain.

---

### 24. Lardenoije et al. (2015) — *Neurobiology of Aging*
**"Age-related epigenetic changes in hippocampal subregions of four animal models of Alzheimer's disease"**

https://www.sciencedirect.com/science/article/pii/S0197458015004121

*Aging effects on cerebellar methylation.* Both 5mC and 5hmC increase with aging in cerebellar Purkinje cells; ratio of 5mC to 5hmC decreases with age.

**Why it matters for your work:** Shows cerebellar methylation is dynamic and functionally relevant throughout life.

---

## Summary: Key Gaps Your Work Addresses

| Gap in Literature | Your Data |
|-------------------|-----------|
| No simultaneous 5mC + 5hmC profiling in BAP1 models | ✓ DUET evoC provides both |
| No BAP1-KO methylation studies in brain | ✓ Mouse cerebellum |
| No gene body methylation analysis in BAP1 loss | ✓ Gene bodies are primary site |
| No coordinated mC↑/hmC↓ pattern described | ✓ 92.3% concordance |

---

# Key Points to Address in Meeting

## Foundational Context to Establish
- [ ] Cerebellum has highest 5hmC levels of any tissue (~40% of modified C)
- [ ] 5hmC is enriched in gene bodies, not promoters—matches your DMR distribution
- [ ] TET enzymes mediate 5mC → 5hmC conversion; blocking this step causes 5mC accumulation + 5hmC depletion
- [ ] BAP1 loss causes H2AK119ub accumulation and chromatin compaction (Conway et al. 2021)

## Your Key Findings to Present
- [ ] 92.3% of genes significant in both modifications show coordinated mC↑/hmC↓
- [ ] Gene bodies (not promoters) are primary affected regions
- [ ] Active Promoter chromatin states preferentially affected (97% of hypermethylated DMRs)
- [ ] Top hit is Syt1 (synaptic transmission) with +18.4% mC / -15.8% hmC

## Mechanistic Interpretation to Discuss
- [ ] Pattern is consistent with impaired TET-mediated oxidation
- [ ] Proposed chain: BAP1 loss → H2AK119ub accumulation → chromatin compaction → TET access restricted
- [ ] Alternative/additional mechanism: DNMT3A recruitment via UDR domain to H2AK119ub regions
- [ ] Connection to loop collapse story: both may stem from same chromatin changes

## What You Can vs. Cannot Claim
- [ ] CAN claim: Coordinated methylation pattern exists and is highly reproducible (92%)
- [ ] CAN claim: First brain BAP1-KO methylation study, first simultaneous 5mC/5hmC in any BAP1 model
- [ ] CANNOT claim: TET access is actually reduced (would need TET ChIP or activity assay)
- [ ] CANNOT claim: H2AK119ub is elevated (would need H2AK119ub ChIP in your samples)

## Caveats to Acknowledge
- [ ] Sample size: n=2 per condition with sex confounding
- [ ] Mechanism is hypothesis—pattern is the discovery
- [ ] Cancer literature may not translate directly to neurons

## Questions to Raise
- [ ] Should we do H2AK119ub ChIP-seq to validate mechanistic model?
- [ ] How does methylation story connect to loop collapse/TAD density story?
- [ ] Rao lab collaboration potential for TET expertise?
- [ ] What experiments would distinguish TET-block vs. DNMT3A-recruitment models?

## Grad Mentor's Framework
> "Long polycomb loops are collapsing inwards to a closer TAD checkpoint, making a higher density TAD"

- [ ] How do Polycomb anchor annotations relate to the DMR locations?
- [ ] Are the genes showing methylation changes also at loop anchors?
- [ ] Is there overlap between loop-associated genes and DMR genes?
