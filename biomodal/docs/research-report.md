# DNA methylation and hydroxymethylation in the brain: BAP1 connections and literature landscape

The brain exhibits exceptionally high levels of 5-hydroxymethylcytosine (5hmC), with Purkinje neurons showing **~40% 5hmC relative to 5mC**—the highest documented in any mammalian cell type. While the mechanisms regulating TET-mediated 5mC-to-5hmC conversion are well characterized, and BAP1/PR-DUB's role in chromatin regulation is established in cancer contexts, **no existing studies have examined whether BAP1 loss affects TET-mediated hydroxymethylation in neural tissue**. The observation of coordinated 5mC increase and 5hmC decrease at gene body loci in BAP1-KO cerebellum would represent a novel mechanistic connection between PR-DUB chromatin function and DNA methylation dynamics—a gap that current literature cannot directly address.

## 5hmC reaches unique abundance in brain and cerebellum

The brain contains approximately **10-fold higher 5hmC levels** than peripheral tissues, a finding consistent across multiple quantification methods. Within the brain, cerebellar neurons are particularly enriched: Kriaucionis and Heintz's 2009 *Science* discovery found that 5hmC accounts for ~0.6% of all nucleotides in Purkinje neurons, equivalent to 40% of total 5mC. Human cerebellum shows substantially higher 5hmC than prefrontal cortex, with **65,563 versus 37,145 detectable probes** in oxBS-450K array studies.

5hmC accumulates specifically in postmitotic neurons rather than proliferating progenitors. Szulwach et al. (2011, *Nature Neuroscience*, PMID: 21804537) demonstrated a **4.22-fold increase** in cerebellar 5hmC from postnatal day 7 to 6 weeks in mice, with continued increase (~21%) through one year of age. This developmental acquisition involves both neurons and glia: recent Nanopore sequencing revealed that astrocytes contribute substantially to brain 5hmC levels, showing lower 5mC and higher 5hmC across specific genes compared to neurons.

Genomically, 5hmC is **enriched in gene bodies** of actively expressed genes (60% of differentially hydroxymethylated regions are intragenic, with 5.7-fold enrichment in exons) while being depleted at promoter CpG islands. This distribution pattern has functional significance: Mellén et al. (2017, *PNAS*) demonstrated that 5hmC accumulation in gene bodies constitutes "functional demethylation" by replacing high-affinity 5mCG binding sites for MeCP2 with low-affinity 5hmCG sites, thereby reducing MeCP2 occupancy and its repressive effects.

## TET enzymes mediate active demethylation through sequential oxidation

TET1, TET2, and TET3 are Fe(II)- and α-ketoglutarate-dependent dioxygenases that catalyze sequential oxidation: 5mC → 5hmC → 5-formylcytosine (5fC) → 5-carboxylcytosine (5caC). The rate constant for 5mC-to-5hmC conversion is relatively fast, whereas further oxidation proceeds considerably more slowly. Two routes complete demethylation: **active (replication-independent)** demethylation via thymine DNA glycosylase (TDG) excision of 5fC/5caC followed by base excision repair, and **passive (replication-dependent)** dilution since DNMT1/UHRF1 does not recognize hemi-modified CpGs containing oxidized bases.

TET enzyme recruitment to chromatin occurs through multiple mechanisms. TET1 and TET3 contain CXXC domains that directly bind CpG-rich DNA, while TET2 (which lacks CXXC) requires recruitment by transcription factors including EGR1, WT1, and lineage-specific factors. All three TET enzymes interact with O-GlcNAc transferase (OGT), which restrains TET activity genome-wide—a 2025 *Nature Structural & Molecular Biology* study demonstrated that OGT prevents inappropriate DNA demethylation in heterochromatin. Post-translational modifications including acetylation (p300-mediated) and monoubiquitination (CRL4^VprBP-mediated) further regulate TET activity and chromatin binding.

In brain, TET3 is the most abundantly expressed TET enzyme. TET-mediated oxidation is essential for proper neuronal maturation: Stoyanova et al. (2021, *eLife*) showed that deletion of Tet1/2/3 from postmitotic Purkinje cells prevents proper differentiation. Individual TET knockout studies reveal distinct functions—Tet1 KO mice show impaired memory extinction and LTD abnormalities, while Tet2 loss paradoxically enhances hippocampal cognitive function.

## TET loss produces paradoxical methylation patterns across chromatin compartments

A critical finding from the Rao lab (López-Moyado et al., 2019, *PNAS*, PMID 31371502) reveals that TET deficiency causes **dual, seemingly contradictory methylation changes**: hypermethylation in euchromatin (Hi-C A compartment) as expected, but also hypomethylation in heterochromatin (Hi-C B compartment). This paradox is explained by DNMT3A redistribution: in TET-deficient cells, DNMT3A relocalizes from heterochromatin to euchromatin sites previously occupied by TET, creating hypermethylation where TETs normally act while depleting methylation maintenance in heterochromatin.

This pattern—global hypomethylation combined with focal hypermethylation—mirrors the epigenetic landscape of cancer genomes. The finding has mechanistic implications: TET enzymes act as "guardians" protecting CpG islands from aberrant methylation. In normal cells, DNMT3A1 can localize to Polycomb-marked CpG islands via its UDR domain recognizing H2AK119ub, but TET enzymes actively counter this methyltransferase activity. When TET activity is insufficient or absent, cancer-associated CGI hypermethylation at Polycomb targets can proceed unopposed.

## BAP1/PR-DUB constrains H2AK119ub to enable proper chromatin architecture

BAP1 is the catalytic subunit of PR-DUB, a Polycomb Repressive-Deubiquitinase complex that removes H2AK119 monoubiquitination. Cryo-EM structures (Thomas et al., 2023, *Science Advances*) reveal that BAP1 binds ASXL1's DEUBAD domain and contacts nucleosomes through multiple anchor points including a "DNA clamp" near the nucleosome dyad. This architecture positions the complex optimally for H2AK119Ub deubiquitination.

The key mechanistic insight from Conway et al. (2021, *Molecular Cell*) is that BAP1 functions as an epigenomic "bookend" constraining H2AK119ub throughout the genome—not just at classic Polycomb targets. BAP1 loss triggers several cascading consequences:

- **H2AK119ub spreads intergenically**, primarily via PCGF3/5-containing non-canonical PRC1 complexes
- This spreading **titrates PRC2 away** from target promoters to intergenic sites
- **H3K27me3 accumulates aberrantly** at intergenic regions while depleting from proper targets
- **Over 75% of the genome** acquires more compact configuration upon BAP1 loss
- **Global reduction in RNA Pol II Ser5 phosphorylation** indicates compromised transcription initiation

Work in *Drosophila* (Bonnet et al., 2022, *Genes & Development*) revealed an important nuance: PR-DUB acts as a "rheostat" where excessive H2Aub1 actually **interferes with nucleosome stacking** and chromatin fiber folding. High H2Aub1 increases DNA accessibility despite maintaining Polycomb protein binding—the mark antagonizes canonical PRC1's compaction activity. This creates apparent paradox: loss of PR-DUB can cause repression defects not from losing deubiquitination per se, but from H2Aub1 accumulation interfering with proper chromatin architecture.

## H2AK119ub directly influences DNA methyltransferase targeting

Multiple 2024 cryo-EM structures have defined how H2AK119ub serves as a recruitment signal for DNMT3A1. The DNMT3A1 N-terminal UDR (ubiquitin-dependent recruitment) domain binds both H2AK119Ub and the nucleosome acidic patch through a bidentate interaction. This mechanism recruits DNMT3A1 to bivalent promoters (H3K4me3+/H3K27me3+) and is required for postnatal neuronal development.

The interplay between Polycomb marks and DNA methylation is complex and context-dependent:

- H2AK119ub depletion leads predominantly to **DNA hypomethylation**
- H3K27me3 depletion paradoxically causes **CpG island hypermethylation** at co-regulated loci
- At certain promoters, H3K27me3 appears to inhibit DNMT recruitment, such that its loss leaves H2AK119ub "unrestrained" to promote aberrant methylation

A major mechanistic discovery (Zou et al., 2024, *Nature*) identified a TET2 → RNA m5C → MBD6 → BAP1 axis: TET2 oxidizes m5C on chromatin-associated retrotransposon RNA, which MBD6 recognizes to recruit BAP1 for H2AK119ub deubiquitination. This pathway represents a direct molecular link between TET activity and H2AK119ub regulation, though it operates through RNA rather than DNA methylation.

## BAP1 studies in brain contexts remain extremely limited

The literature gap for BAP1 in brain is substantial. Three areas have received any attention:

**Küry-Isidor syndrome** (KURIS): Küry et al. (2022, *American Journal of Human Genetics*) described 11 individuals with de novo heterozygous **missense** BAP1 variants exhibiting developmental delay, intellectual disability, speech delay, hypotonia, seizures, and autism spectrum features. These are partial loss-of-function or dominant-negative variants—mechanistically distinct from the complete loss-of-function causing tumor predisposition.

**Enteric nervous system**: A 2024 *JCI* study demonstrated BAP1 is critical for postnatal enteric neuron differentiation and survival, with prenatal Bap1 loss causing severe bowel dysfunction and enteric neuron loss by P15—representing first evidence linking prenatal PR-DUB defects to postnatal neurodegeneration.

**Neuroblastoma**: BAP1 induces cell death via interaction with 14-3-3 proteins, releasing Bax to promote apoptosis, with high expression correlating with better outcomes.

Virtually nothing is known about BAP1 function in CNS neurons, its role in neuronal differentiation beyond ENS, H2AK119ub dynamics in brain development with BAP1 focus, or chromatin compaction effects in neural tissue. The mechanistic insights from cancer and ESC studies have not been translated to neural contexts.

## The mechanistic chain from BAP1 to DNA methylation: what literature supports

The proposed model—BAP1 loss → H2AK119ub accumulation → chromatin compaction → reduced TET access → 5mC increase + 5hmC decrease—has variable support across its components:

**Strongly supported**: BAP1 loss causing H2AK119ub accumulation is well-established. BAP1 loss causing DNA methylation changes is documented in uveal melanoma (Field et al., 2019, *Clinical Cancer Research*), where knockdown induced methylomic repatterning with predominant hypermethylation at gene promoters. DNMT3A1 recruitment to H2AK119ub-marked regions via UDR domain is structurally characterized.

**Mechanistically plausible but not directly demonstrated**: Chromatin compaction restricting TET enzyme access to substrate DNA has not been explicitly shown. TET activity correlates with chromatin accessibility and 5fC/5caC are enriched in open chromatin, but whether H2AK119ub-dependent compaction specifically limits TET access is unstudied.

**Not addressed in existing literature**: Simultaneous profiling of 5mC and 5hmC changes in BAP1-deficient models has not been performed. The specific pattern of coordinated 5mC increase + 5hmC decrease at the same gene body loci in BAP1-KO is not described anywhere. Gene body methylation dynamics in BAP1 loss are largely unexplored (existing work focuses on promoters). Neural tissue/cerebellum BAP1-KO methylation studies do not exist.

## Conclusion: significant novelty in the observed methylation phenotype

The literature landscape reveals robust understanding of individual components—brain 5hmC biology, TET enzyme mechanisms, BAP1/PR-DUB function, H2AK119ub effects—but a **critical absence of integrative studies** connecting these systems in neural contexts. The Rao lab's work on TET and the recent structural/mechanistic characterization of PR-DUB provide strong foundations for understanding how chromatin state might regulate DNA methylation dynamics.

An experimental observation of coordinated 5mC increase and 5hmC decrease at gene body loci in BAP1-KO mouse cerebellum would be genuinely novel on multiple fronts: **first characterization of DNA methylation changes in brain BAP1 loss**, **first simultaneous 5mC/5hmC profiling in any BAP1-KO model**, and **first direct evidence that BAP1-dependent chromatin changes affect TET-mediated hydroxymethylation**. The gene body localization adds further novelty since existing BAP1-methylation work has focused on promoter changes.

The mechanistic interpretation—that H2AK119ub accumulation creates chromatin environment restricting TET access—is consistent with existing knowledge but would require additional experiments to establish causality. Such a finding would fill an important conceptual gap linking PR-DUB chromatin regulation to the DNA methylation/demethylation balance that is particularly relevant in the 5hmC-enriched brain environment.