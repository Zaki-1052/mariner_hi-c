# The BAP1 Story: How a Broken Eraser Rewrites the 3D Genome

**A narrative synthesis of the BAP1 Hi-C project**
*For: anyone who wants to understand what this paper is about, what we found, and why it matters*

---

## The One-Paragraph Version

BAP1 is an enzyme that erases a specific Polycomb mark (H2AK119ub) from chromatin. When it's lost in the mouse cerebellum, the 3D folding of the genome gets progressively rewritten: long-range chromatin loops that connected distant developmental regulators collapse, and shorter-range contacts form in their place. This isn't random degradation -- it's a structured reorganization driven by the accumulation of that unerased mark at two types of sites, through two distinct mechanisms, producing a single convergent phenotype: the spatial regulatory grammar of the neuron is being overwritten, one loop at a time. The effect is progressive (18-fold amplification from early postnatal life to adulthood), scale-dependent (compartments massively affected, TAD boundaries mostly spared, cohesin stripes untouched), and functionally consequential (88% concordance between architectural disruption and gene expression changes at the same loci). The affected genes are the ones you'd expect if cerebellar neuronal identity were unraveling.

---

## Part 1: The Setup — What BAP1 Does and Why We Looked at 3D Chromatin

### The Enzyme

BAP1 (BRCA1-Associated Protein 1) is a deubiquitinase. Its job is to remove monoubiquitination from histone H2A at lysine 119 -- a mark we call H2AK119ub, deposited by Polycomb Repressive Complex 1 (PRC1). In the classical textbook view, PRC1 puts this mark down to silence genes, and BAP1 takes it off to allow activation. But the reality in the brain is more nuanced than that.

The Ferguson lab (where this work was done) had previously shown something unexpected: H2AK119ub isn't just sitting at silenced heterochromatin in the cerebellum. It's dynamically present at *active* enhancers. This means BAP1 isn't just a passive clean-up enzyme -- it's part of an active cycling process where PRC1 constantly deposits the mark and BAP1 constantly removes it. The *turnover* matters, not just the steady-state level.

### The Mouse Model

BAP1 was knocked out via CRISPR in the mouse cerebellum. These mice develop neurodegeneration -- seizures, cerebellar atrophy, progressive loss of neuronal identity. The question was: *how does the loss of this one enzyme affect the physical folding of chromosomes?*

### Why 3D Chromatin?

Genes don't work in isolation along a linear chromosome. Enhancers can sit hundreds of kilobases or even megabases away from the genes they regulate, and they're brought into contact through chromatin loops -- physical bridges of folded DNA. If those loops change, the regulatory wiring of the cell changes. We used Hi-C (a genome-wide chromosome conformation capture method) at two developmental timepoints -- postnatal day 12 (P12, early) and adult (P60, late) -- with three biological replicates per condition, to ask whether BAP1 loss disrupts this wiring.

---

## Part 2: The Central Discovery — Loop Rewriting

### The Headline

In the adult BAP1-KO cerebellum, we identified 2,910 differential chromatin loops. But the story isn't just "loops change." It's *how* they change.

Lost loops are long. Gained loops are short.

Lost loops have a median genomic distance of 625 kb. Gained loops have a median of 320 kb. Loops spanning more than a megabase are over three-fold enriched for being lost. This is not subtle -- the Kolmogorov-Smirnov test for the distance distributions has a p-value below the limit of floating-point arithmetic.

### What "Loop Rewriting" Means Biologically

Imagine a neuron's genome as a city with a transit system. Long-range loops are like subway lines connecting distant neighborhoods -- they bring a silenced developmental gene in one district into the same regulatory space as a Polycomb hub in another. Short-range contacts are like walking between adjacent buildings.

What BAP1 loss does is systematically shut down the subway while making the sidewalks more crowded. The long-distance connections that linked developmental regulatory hubs -- the ones maintained by Polycomb-mediated higher-order organization -- preferentially fail. The interaction "budget" that was spent maintaining those long-range contacts redistributes locally, creating denser short-range contact neighborhoods.

This isn't cells randomly losing loops and gaining others. At 212 specific genomic loci, the *same anchor site* loses a distant contact partner and gains a closer one. These "shared anchor" hubs are enriched for H3K27me3 (the Polycomb mark) and depleted for H3K27ac (the active mark), confirming the switching is concentrated at Polycomb-regulated sites. The directionality is remarkably consistent -- at 83% of these hubs, the lost loop is longer than the gained one.

---

## Part 3: The Mechanism — Two Contexts, One Broken Switch

This is where the biology gets interesting, and where the apparent paradox in the data resolves into a coherent model.

### The Paradox

H2AK119ub accumulates in the mutant. At active enhancers, this accumulation correlates with loop *weakening*. At Polycomb domains, the story is different -- long-range contacts weaken but short-range ones strengthen. How can the same mark have seemingly opposite effects?

### The Resolution: Context-Dependent Consequences of Lost Turnover

The answer is that BAP1's job isn't just "removal" -- it's *dynamic cycling*. In the wildtype, PRC1 deposits H2AK119ub and BAP1 removes it, over and over. This turnover is what the cell actually depends on. Different chromatin contexts use that cycling differently, and when it stops, different things break.

**At active enhancers and promoters**, BAP1 was keeping these sites "clean" -- removing H2AK119ub that PRC1 deposits as part of its normal scanning activity. Without BAP1, the mark accumulates like rust on an unlocked gate. The enhancer's ability to maintain contact with its target gene degrades. This is the direct, enzymatic consequence of losing BAP1: substrate accumulates at the site where the enzyme was needed, and the functional output (the enhancer-promoter loop) weakens. The data show this clearly -- higher K119ub accumulation at an anchor makes that loop ten times more likely to be lost (odds ratio 10.7 in logistic regression).

**At Polycomb-repressed domains**, the mechanism is different. These are regions already decorated with H3K27me3, where long-range contacts connect distant silenced regions into nuclear compartments through higher-order Polycomb body organization. These megabase-scale contacts depend on dynamic reorganization of Polycomb bodies -- possibly involving phase separation -- and that dynamism requires turnover. When H2AK119ub gets "stuck," the Polycomb landscape can't reorganize properly. Long-range contacts (the ones requiring the most organizational energy to maintain) weaken, while local compaction (driven by nucleosome-nucleosome interactions that H2AK119ub actually promotes at short range) persists or even increases. This is why shared anchors gain short-range contacts as they lose long-range ones -- the interaction potential just redistributes locally.

**The unifying principle**: BAP1 loss eliminates H2AK119ub turnover, causing the same mark to act as a "stuck switch" in two contexts. At active elements, it introduces a foreign repressive signal that blocks enhancer-promoter contacts. At Polycomb domains, it breaks the dynamic organization needed for long-range connectivity while preserving local compaction. Both effects push interactions toward shorter range. The result, at the genome-wide level, is loop rewriting.

---

## Part 4: The Temporal Story — From a Whisper to a Catastrophe

### The Numbers

At postnatal day 12, only 165 differential loops are detected. By adulthood, 2,910. That's an 18-fold amplification.

### What This Means

This isn't just "more of the same thing happening over time." The *character* of the disruption changes.

**Early (P12)**: The disease is just starting. The changes are subtle, concentrated at poised and repressed regulatory elements -- Repressed_Promoter anchors dominate at 36% of all differential anchors, compared to a much more mixed distribution later. The majority of differential loops are weakened (57% down). What you're seeing is the initial destabilization -- developmental contacts that should be maintained during cerebellar maturation are beginning to fail. The chromatin landscape is still largely intact, and the system is trying to hold together.

**Late (P60)**: By adulthood, the Polycomb landscape has been fundamentally remodeled. 44% of the genome shows compartment shifts. Now the majority of differential loops are *gained* (59% up) -- but these aren't healthy new connections. They're the structural consequence of architectural collapse: ectopic contacts forming as boundaries fail, compensatory interactions as the cell tries to maintain gene expression through alternative routes, and the passive redistribution of contact probability into local neighborhoods as long-range connections dissolve.

The direction reversal from early to late is key. Early, the dominant signal is *loss of what should be maintained*. Late, it's *gain of what shouldn't exist*. The cell starts by losing its developmental wiring, then fills in the gaps with noise.

### Why Progressive?

The temporal cascade likely reflects a self-amplifying cycle:

1. Initial H2AK119ub accumulation destabilizes a small number of contacts
2. Those contact changes alter local chromatin environment
3. The altered environment permits further Polycomb redistribution
4. Which destabilizes more contacts
5. Which further alters chromatin state...

By P60, 44% of the genome has shifted compartment identity. That's not the direct footprint of BAP1 loss at individual loci -- it's the downstream consequence of a cascade that started small and propagated through the interconnected chromatin regulatory network.

---

## Part 5: The Multi-Scale Picture — What Changes and What Doesn't

We analyzed four scales of chromatin organization, and they form a clear hierarchy of sensitivity:

### Compartments: Massively Affected (44% of genome)

Compartments are the broadest scale of 3D organization -- the division of the genome into transcriptionally active (A) and inactive (B) regions. Nearly half the genome shifts compartment identity, with a net trend toward more A-like (active) states -- about 137 Mb shifts B-to-A versus 68 Mb shifting A-to-B, a 2:1 ratio. This makes sense: loss of Polycomb repression (due to H2AK119ub accumulation disrupting normal Polycomb function, paradoxically) shifts chromatin toward a more permissive state. But "more permissive" doesn't mean "properly regulated." These aren't regions becoming healthily active -- they're regions losing the repressive identity that was keeping them in the B compartment, without gaining proper activation.

### Loops: The Primary Functional Unit (2,910 differential)

Loops are where the regulatory consequence is most directly measurable. Each loop connects a specific regulatory element to a specific gene. The loop rewriting phenomenon is the core finding.

### TAD Boundaries: Moderately Stable (16-20% differential)

Topologically associating domains -- the megabase-scale neighborhoods that organize the genome -- are relatively preserved. Their boundaries are maintained primarily by CTCF and cohesin, structural proteins whose binding isn't directly affected by H2AK119ub. BAP1 loss reorganizes the regulatory landscape *within* TADs more than the boundaries *between* them.

That said, TAD boundaries aren't completely immune. Where they do change, the changes are coordinated with loop changes: TAD merges are enriched near lost loops (boundaries become unnecessary when the long-range contacts that defined them dissolve), and boundary strengthening is enriched near gained loops (increased local contact density makes existing boundaries more prominent).

### Stripes: Preserved (no significant changes)

Stripes reflect cohesin-mediated loop extrusion -- a structural process mechanistically independent of Polycomb/H2AK119ub. Their preservation is the "dog that didn't bark": it confirms that BAP1's effects are channeled specifically through the Polycomb axis, not through the cohesin/extrusion machinery.

### Why This Hierarchy Matters

This pattern is exactly what you'd predict if BAP1 were specifically a Polycomb regulator. Polycomb-dependent structures (compartments driven by Polycomb-associated heterochromatin identity, Polycomb-mediated long-range loops) are heavily affected. CTCF/cohesin-dependent structures (TAD boundaries, stripes) are largely spared. The hierarchy isn't just a description of what changed -- it's mechanistic evidence that the disruption flows through the Polycomb pathway.

---

## Part 6: The Enhancer Story — When Contact Loss Actually Matters for Genes

### The ABC Model

To connect architectural changes to gene regulation, we implemented an Activity-by-Contact (ABC) model that predicts enhancer-gene linkages based on two inputs: enhancer activity (from H3K27ac and ATAC-seq) and 3D contact frequency (from Hi-C). By comparing WT and KO ABC scores, we could ask: at sites where loops change, do enhancer-gene connections actually break? And do genes actually change expression?

### The Concordance

**88% three-way concordance.** When we required geometric co-localization (enhancer at one loop anchor, gene TSS at the other), loop changes, enhancer-gene score changes, and differential gene expression agreed at the same loci 88% of the time. This is the strongest evidence that the architectural changes are not epiphenomenal -- they directly track with transcriptional disruption.

### The Threshold: K119ub Alone Is Necessary But Not Sufficient

One of the most informative findings is what happens at enhancers where *only* K119ub changes but the active marks (H3K27ac, ATAC accessibility) stay the same. These "K119ub_Only" enhancers show real contact weakening -- about 3.7% reduction, statistically significant. But this contact perturbation doesn't cross a functional threshold: the concordance between contact changes and gene expression changes at these sites is 49.4%, which is indistinguishable from chance.

By contrast, enhancers that actually lose their H3K27ac mark ("Activity_Lost" class) show 59.8% concordance with gene expression.

This establishes a dose-response model: H2AK119ub accumulation is mechanistically upstream of contact disruption, but the contact change must reach a critical magnitude -- or be accompanied by loss of activating marks -- before the downstream gene is affected. K119ub alone perturbs the system; it takes additional cooperative disruption to break gene expression. This has implications for understanding why BAP1-associated diseases develop progressively: the initial epigenetic perturbation may be tolerated until it crosses a threshold.

---

## Part 7: What This Means — The Model and Its Implications

### The Unified Model

BAP1 loss eliminates the dynamic cycling of H2AK119ub. This has cascading consequences across the 3D genome:

1. **Direct enzymatic effect**: H2AK119ub accumulates at sites where BAP1 was needed, weakening enhancer-promoter contacts (particularly at active regulatory elements)

2. **Structural reorganization**: The Polycomb landscape can no longer maintain long-range contacts, leading to collapse of developmental loops and passive gain of short-range contacts

3. **Progressive amplification**: Initial perturbations cascade through the interconnected chromatin network, amplifying from 165 to 2,910 differential loops over postnatal development

4. **Functional consequence**: Disrupted enhancer-gene connections lead to dysregulation of synaptic and developmental genes (BDNF, Mef2c, Pax6 at lost loops; Sox2/6/9, Foxg1 at gained loops)

5. **Scale-dependent sensitivity**: Effects propagate primarily through Polycomb-dependent structures (compartments, loops) while sparing Polycomb-independent ones (CTCF boundaries, cohesin stripes)

### Why This Matters for the Field

**For chromatin biology**: This is among the first demonstrations that a single histone modification enzyme can cause structured, genome-wide 3D reorganization -- not just at its direct target sites, but through cascading effects across the Polycomb regulatory network. The "loop rewriting" phenomenon (systematic replacement of long-range with short-range contacts) hasn't been described before in this context.

**For BAP1 biology**: BAP1 mutations cause cancer (mesothelioma, uveal melanoma, renal cell carcinoma) and neurodegeneration. This work suggests a unifying mechanism: progressive 3D genome reorganization that gradually overwrites the cell's regulatory identity. The disease isn't about losing individual gene contacts -- it's about losing *control* over which contacts form and when.

**For neurodegeneration**: The temporal progression -- subtle at P12, catastrophic by P60 -- mirrors the clinical trajectory of neurodegenerative diseases. The self-amplifying cascade model suggests there may be a window early in disease where intervention could prevent the architectural collapse from propagating.

### The Syt1/Nav3 Locus

The most dramatically affected region is the Syt1/Nav3 locus on chromosome 10. Syt1 (Synaptotagmin 1) is a calcium sensor essential for synaptic vesicle fusion -- about as core a neuronal gene as you can get. By adulthood, the regulatory architecture around this locus has essentially collapsed. This isn't just a statistical finding; it's a concrete example of how 3D genome disruption connects directly to neuronal dysfunction.

---

## Part 8: Key Loci and Genes

| Gene | What It Does | What Happens |
|------|-------------|-------------|
| **Syt1** | Synaptic vesicle fusion | Regulatory architecture collapses -- most impacted locus |
| **Nav3** | Neuronal migration | Co-impacted with Syt1 |
| **BDNF** | Neuronal survival/plasticity | At lost loops -- connections to enhancers severed |
| **Mef2c** | Neuronal differentiation TF | At lost loops |
| **Pax6** | Brain development TF | At lost loops |
| **Sox2/6/9** | Developmental TFs | At gained loops -- ectopic contacts |
| **Foxg1** | Forebrain development | At gained loops |
| **Shh** | Morphogen signaling | At gained loops |
| **Apaf1** | Apoptosis regulator | Only significant stripe change (gained) |

The pattern: neuronal identity and function genes at lost loops (developmental wiring severed), developmental transcription factors at gained loops (ectopic contacts forming).

---

## Analytical Reasoning & Extended Notes

*This section is a reference for understanding the logic behind the interpretations above.*

### Why "Loop Rewriting" and Not Just "Loop Loss"

It would be simpler to say "BAP1 loss weakens chromatin loops." But that misses something fundamental: 59% of differential loops are *gained*, not lost. If this were simple degradation, you'd expect predominantly loss. The gain of short-range contacts is not noise -- it's the structural consequence of the redistribution of interaction probability when long-range contacts dissolve. The genome has a finite interaction budget, and when long-range contacts fail, that budget flows into the local neighborhood.

The shared anchor analysis proves this isn't independent gain and loss at different sites. At 212 hubs, the same genomic position is involved in both a lost long-range loop and a gained short-range one. The switching is directionally consistent (83%), distance-consistent (lost loops 3.4x longer than gained at the same anchor), and enriched at Polycomb sites. This is structured reorganization, not random shuffling.

### The Positive Rho at Polycomb Anchors

One data point that initially seems contradictory: at Polycomb/repressive anchors, higher K119ub correlates with *more* loop strength (positive Spearman rho of +0.177), while at active anchors, higher K119ub correlates with *less* loop strength (negative rho of -0.314).

This actually confirms the model. At Polycomb sites, H2AK119ub promotes local nucleosome compaction -- that's what the mark does at the molecular level. So more K119ub at a Polycomb anchor = stronger local compaction = stronger short-range contacts. But the long-range contacts that require higher-order Polycomb body organization (not just nucleosome-level compaction) still weaken. K119ub has a distance-dependent dual effect at Polycomb domains: promotes short-range compaction, disrupts long-range connectivity.

At active anchors, K119ub is entirely foreign -- it's the wrong mark in the wrong place. There's no distance-dependent nuance; it simply interferes with enhancer-promoter contact at any range.

### Why Compartments Change More Than TADs

This seems paradoxical: compartments (the largest-scale organization) change the most, while TADs (an intermediate scale) change the least of any affected structure. But compartment identity is defined by the *aggregate* epigenetic state of a region -- the sum of all the histone marks, binding proteins, and transcriptional activity across hundreds of kilobases. Even modest per-locus Polycomb redistribution, when summed across many sites, produces massive compartment-level shifts. TAD boundaries, by contrast, are point features maintained by specific protein complexes (CTCF/cohesin) that are not directly affected by H2AK119ub.

### The 88% Concordance Is Geometrically Constrained

Without the geometric constraint (requiring enhancer at one anchor, TSS at the other), the concordance between loop changes and ABC model predictions is 51% -- essentially chance. This matters: it means you can't just overlap two genome-wide datasets and expect meaningful biology. The physical requirement that the enhancer and gene sit at opposite ends of the same loop is what makes the integration meaningful. When you respect the geometry, concordance jumps to 88%.

### Why the Early Timepoint Has So Few Loops

165 vs 2,910 is a dramatic difference, and one might worry this is a power issue (fewer sequencing reads, less signal). But the statistical framework is the same, the sequencing depth is comparable, and the BCV (biological coefficient of variation) is actually lower at early timepoint. The difference is biological: at P12, the Polycomb landscape hasn't had time to fully reorganize. The 165 loops represent the *leading edge* of the disruption -- the most vulnerable contacts that fail first. By P60, the cascade has propagated through the network.

### The Stripe Null Finding

Finding zero significant stripe changes (out of ~200-286 tested) is as informative as finding 2,910 differential loops. Stripes are generated by cohesin-mediated loop extrusion -- a mechanical process that doesn't depend on Polycomb marks. If BAP1 loss were causing generalized chromatin destabilization (rather than specific Polycomb pathway disruption), stripes would be affected too. Their preservation is strong evidence for pathway specificity.

### On the Direction Reversal (Early DOWN-dominant, Late UP-dominant)

This reversal has a biological interpretation beyond "more stuff breaks over time." At P12, the dominant signal is failure of existing contacts -- loops that should persist during cerebellar maturation are weakening. This is the *direct* consequence of H2AK119ub accumulation. By P60, the system has had time to reorganize: boundaries have merged, Polycomb domains have restructured, and the interaction probability that was held in long-range contacts has redistributed into short-range neighborhoods. The "gained" loops aren't functional new connections -- they're the structural debris of architectural collapse, plus whatever compensatory contacts the cell has managed to establish.

The analogy: if you dam a river (P12 = initial blockage, downstream flow reduced), eventually the water backs up and flows over the banks in new channels (P60 = redirected flow, more total water movement, but in the wrong places).

---

*Document synthesized April 2026*
*Source: mariner_hi-c repository, Ferguson Lab, UCSD*
