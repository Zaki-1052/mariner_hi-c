# BAP1 Hi-C: What's going on and why it matters

## What BAP1 actually does and why we're looking at chromatin folding

So the starting point is that BAP1 is a deubiquitinase — it removes H2AK119ub, the ubiquitin mark that PRC1 puts on histone H2A. In a normal cell, PRC1 is constantly depositing this mark and BAP1 is constantly removing it. The important thing is that this isn't a static system where the mark is either "on" or "off" — it's a cycle. The mark is being written and erased continuously, and that dynamic turnover is itself functionally important. When you lose BAP1, PRC1 keeps writing but nothing erases, so H2AK119ub accumulates. But the problem isn't just "too much ubiquitin" — it's that the cycling stops.

The reason this matters for chromatin folding specifically comes from a key finding in the Ferguson lab: H2AK119ub isn't just sitting at silenced Polycomb domains in the cerebellum. It's actually enriched at active enhancers too. That's not what you'd naively expect from a "repressive" mark. It means BAP1 is doing real work at active regulatory elements — constantly removing K119ub to keep those enhancers "clean" and functional. And it's also doing work at Polycomb domains, where the dynamic cycling of K119ub is part of how those long-range repressive contacts are organized and maintained.

So when BAP1 is gone, you'd expect problems at both types of sites, but for different reasons. At active enhancers, K119ub accumulates where it doesn't belong — it's a foreign mark that interferes with enhancer-promoter interactions. At Polycomb domains, the mark was already there but was being dynamically cycled; without cycling, the higher-order organization that depends on that dynamism breaks down. The question was whether this actually manifests as changes in 3D chromatin structure, and if so, whether those structural changes actually matter for gene expression.

That's why we did Hi-C.

## What happens to loops

The most immediate finding is that BAP1 loss causes widespread loop dysregulation in the adult cerebellum. Roughly one in five detected loops are significantly different between mutant and control. But it's not a uniform loss or gain — there's a very clear directional pattern that depends on distance.

Long-range loops are preferentially lost. Loops spanning more than a megabase are heavily enriched for weakening in the mutant. Short-to-mid-range loops, on the other hand, are preferentially gained. The median distance of lost loops is about twice the median distance of gained loops. This is highly significant statistically and holds up across every test you throw at it.

We call this "loop rewriting" because the genome isn't just losing contacts or gaining contacts — it's replacing one type with another. Long developmental contacts are being swapped out for shorter, more local ones. The loop distance shift is really the central structural observation of the whole project.

### The shared anchor phenomenon

The strongest evidence for this replacement pattern comes from shared anchors. At a couple hundred genomic loci, the exact same anchor site participates in both a lost loop and a gained loop. It loses a long-range partner and gains a shorter-range one. The directionality is very consistent — at the vast majority of these shared anchors, the lost contact is longer than the gained contact.

What's at these shared anchors? They're enriched for H3K27me3 (Polycomb) and depleted for H3K27ac (active marks). So the switching phenomenon is concentrated at Polycomb-regulated sites. This makes sense with the model: Polycomb hubs that were maintaining long-range contacts through dynamic K119ub cycling can no longer sustain those contacts when BAP1 is gone. The long-range partner is lost, and the interaction probability redistributes locally, creating shorter-range contacts that show up as "gained" loops.

This is important because it means the gained loops at these sites are not new functional enhancer-promoter connections. They're the structural consequence of long-range contact collapse — the interaction "budget" getting redirected into the local neighborhood.

## How K119ub drives the changes depending on context

This is where the biology gets interesting and where a lot of the analytical reasoning lives.

When you correlate K119ub change at loop anchors with loop strength change, you get different answers depending on what kind of chromatin you're looking at. At active enhancer anchors (H3K27ac-positive), more K119ub means weaker loops — negative correlation. At Polycomb/repressive anchors (H3K27me3-positive), more K119ub means stronger loops — positive correlation. At first glance that looks contradictory. It's not.

The key is that these are different populations of loops being measured, and K119ub is doing different things in different contexts.

At active enhancers, BAP1 was keeping K119ub away. These are sites where the mark is foreign — it doesn't belong there. When BAP1 is lost and K119ub accumulates at active enhancers, it directly interferes with the enhancer-promoter contacts those elements maintain. The correlation is negative at every distance. This is the simplest, most direct consequence of BAP1 loss: the substrate accumulates at sites where the enzyme was needed, and the functional output (E-P contact) degrades.

At Polycomb sites, the situation is more nuanced. K119ub was already present — it's a normal part of Polycomb chromatin. What matters here is the loss of dynamic cycling, not the accumulation per se. The positive correlation between K119ub and loop strength at Polycomb anchors is driven by short-range contacts. PRC1-mediated nucleosome compaction is a local phenomenon — K119ub promotes compaction within roughly a hundred kilobases. So at Polycomb sites that retain or gain K119ub, local short-range compaction persists or even increases. But the long-range Polycomb contacts — the ones that span megabases and depend on higher-order organization like Polycomb body formation — those collapse because they need dynamic turnover to be maintained.

So the positive rho at Polycomb anchors doesn't mean K119ub is helping there. It means the surviving loops at those anchors are the short-range ones, and there are more of them after the long-range ones die. The shared anchor data confirms this directly: same anchor, lost partner at a megabase, gained partner at a few hundred kilobases.

The logistic regression ties it together globally: K119ub fold-change at an anchor is by far the strongest predictor of whether a loop is lost, with an odds ratio over ten. Distance is the second strongest predictor. Both effects are independent and multiplicative. A loop at a K119ub-accumulating anchor that also happens to be long-range is extremely likely to be lost.

## The temporal story

The early timepoint (P13) and the late timepoint (adult) tell very different parts of the same story.

At P13, there's barely any signal — the number of differential loops is tiny compared to adult. But the composition of those early changes is informative. The anchor types are dominated by Repressed_Promoter, which accounts for over a third of anchors at P13 versus under ten percent at adult. Active marks are almost absent. The early disruption is concentrated at Polycomb-repressed sites — which makes sense, because those are the sites where K119ub cycling is most critical for maintaining chromatin organization, and they'd be the first to feel the loss of BAP1.

The direction is also different. At P13, the majority of differential loops are weakened (lost). At adult, it reverses — the majority are gained. This reversal doesn't mean things are getting better. The late-stage gains represent the downstream structural consequences of the early losses: redistribution of contacts into shorter range, boundary failure creating ectopic interactions, compensatory contacts, and Polycomb domain compaction filling the void left by collapsed long-range loops. The gained loops have weaker individual effect sizes than the lost ones — many small gains versus fewer strong losses.

There's also a progression in which types of chromatin are affected. Early: mostly Polycomb/repressive. Late: everything — active enhancers, poised enhancers, Polycomb, CTCF sites, all involved. The disease starts as a targeted disruption at Polycomb-regulated loci and expands into a system-wide architectural breakdown.

The mechanistic interpretation — which is speculative but consistent with the data — is that there are two phases. First, BAP1 loss directly disrupts Polycomb-dependent long-range contacts because those contacts depend on dynamic K119ub cycling. Second, the Polycomb landscape globally reorganizes in response, redistributing marks and association preferences across the genome, which then indirectly affects non-Polycomb structures. By adulthood, the cascade has propagated to the point where the chromatin architecture is fundamentally remodeled.

The abstract frames this as: as increasing ubiquitinated histone collapses developmental loops and replaces them with more proximal contacts, the result is dysregulation of synaptic and developmental genes.

## What happens at other structural scales

### Compartments

Compartments are the most dramatically affected scale. About seven to eight percent of the genome shows significant compartment shifts at standard thresholds at both timepoints, with roughly twice as many regions shifting toward the A compartment (active) as toward B (inactive). This net B-to-A shift is consistent with loss of Polycomb repression — without proper Polycomb-mediated silencing, formerly repressed regions shift toward a more active identity.

An interesting detail is that compartment changes are already at full magnitude at P13, even though loop changes are minimal at that age. This seems counterintuitive until you remember that compartments are defined by the aggregate behavior of many loci. They're the coarsest scale of chromatin organization — the sum of many small association-preference changes across millions of base pairs. Even modest per-locus changes in Polycomb mark distribution can sum to detectable compartment shifts. Loops, by contrast, are individual contacts that each need to cross a statistical threshold to be called differential. So compartments can change subtly but measurably while individual loops haven't yet reached significance.

The composition of compartment changes does shift between timepoints though. Early has more complete A-to-B and B-to-A flips (sign changes in PC1), while late has more strengthening and weakening of existing compartment identity. This suggests the compartment landscape is actively flipping at P13 and then settling into a remodeled but more stable state by adulthood.

### TADs

TAD boundaries are moderately affected — about a fifth are differential at either timepoint. This is notably less dramatic than compartments, and the percentages barely change between early and late even as loops amplify eighteen-fold. TAD boundaries are primarily maintained by CTCF and cohesin, structural proteins whose binding is not directly regulated by H2AK119ub. So you'd expect them to be relatively resilient to BAP1 loss, and they are.

But they're not just passive bystanders either. The TAD boundary changes that do occur are spatially and directionally coordinated with loop changes. Lost loops tend to sit closer to differential boundaries than gained loops. There's about seventy percent directional concordance — when a loop is lost near a boundary, that boundary tends to be control-enriched, and when a loop is gained, nearby boundaries tend to be mutant-enriched.

The specific types of boundary changes are mechanistically informative. Merge boundaries (two TADs fusing) are heavily enriched near lost loops — when a long-range loop collapses, the TAD boundary it spanned becomes unnecessary, and the two domains merge. Strength-change boundaries are enriched near gained loops — the new shorter-range contacts densify local structure. Split boundaries are also enriched near gained loops. So the TAD-level story is: loops collapse, TADs merge where the loops were; loops gain locally, TADs densify or subdivide.

This provides the intermediate-scale link between loop rewriting and compartment shifts. The loops change, the TADs respond by reorganizing their domain structure, and the aggregate effect of all of this shows up as compartment shifts.

### Stripes

Stripes are essentially unaffected. No significant changes at either timepoint. This is actually informative because stripes reflect cohesin-mediated loop extrusion — a mechanical process that doesn't depend on Polycomb or H2AK119ub. The fact that stripes are preserved confirms that BAP1's effects are channeled through the Polycomb axis specifically, not through some general disruption of chromatin compaction or extrusion machinery.

### The hierarchy

So you get a hierarchy of sensitivity: compartments (most affected) > loops > TAD boundaries (moderate) > stripes (unaffected). This maps cleanly onto the expected biology: compartments and Polycomb loops depend on PRC1/PRC2 activity; CTCF boundaries and cohesin stripes don't. A Polycomb regulator should primarily affect Polycomb-dependent structures, and it does, at every scale.

## When do contact changes actually matter for genes?

This is where the ABC (Activity-By-Contact) model comes in. The core question is: of all these loop changes, which ones actually translate into gene expression changes?

The answer turns out to depend on what else is happening at the enhancer besides the contact change.

The analysis classifies enhancers into groups based on what changes in the mutant. The critical comparison is between two groups: enhancers where K119ub accumulates AND the active mark H3K27ac is lost ("Activity_Lost"), versus enhancers where K119ub accumulates BUT active marks remain unchanged ("K119ub_Only").

Both groups show real loop weakening. K119ub accumulation at an enhancer weakens its contacts whether or not the enhancer loses its active marks. But the functional consequence is different. Activity_Lost enhancers show clear concordance between contact loss and gene downregulation — the loop weakens, the enhancer-gene connection degrades, and the target gene actually goes down. K119ub_Only enhancers show loop weakening that is statistically real but doesn't translate into detectable gene expression changes. The concordance between contact change and gene expression at these sites is at chance level.

This establishes something like a threshold model. K119ub accumulation is mechanistically upstream of contact disruption — it weakens loops. But the contact change alone doesn't cross a functional threshold. The enhancer needs to also lose its activity (H3K27ac) before the downstream gene is affected. In the steady-state adult cerebellum, K119ub alone is necessary but not sufficient for transcriptional consequence.

This has implications for understanding the progressive nature of the disease. Early on, K119ub accumulates and starts weakening contacts, but many enhancers can tolerate this because their active marks are still intact. Over time, as the Polycomb landscape continues to remodel and active marks are secondarily affected, more enhancers cross the threshold into functional disruption. The initial epigenetic perturbation is tolerable; the cascade isn't.

When you do enforce a geometric constraint — requiring the enhancer to sit at one loop anchor and the gene's TSS at the other — the concordance between loop changes, ABC score changes, and RNA-seq changes jumps dramatically. This tells you that the structural contact and the regulatory connection are mechanistically linked, not just coincidentally overlapping. The loops aren't epiphenomenal; they're carrying the regulatory signal.

## How to think about the paradoxes and connections

### Why gained loops aren't good news

At the late timepoint, there are more gained loops than lost. This seems counterintuitive if BAP1 loss is supposed to be disruptive. But the gains aren't functional recovery. There are several things contributing to the gain signal, and none of them represent healthy chromatin organization.

When a long-range loop collapses, the interaction probability doesn't vanish — it redistributes locally. The shared anchor data shows this literally: the same anchor loses a megabase-scale partner and gains a few-hundred-kilobase partner. The "gained" loop is just where the lost interaction's budget went.

TAD merging creates ectopic contacts — regions that were previously insulated by a boundary can now interact. Those show up as gains. The cell may also attempt to compensate by forming alternative enhancer-promoter contacts. And at Polycomb domains, K119ub accumulation promotes local nucleosome compaction even as long-range contacts fail, creating many weak short-range gains.

The key tell is that gained loops have weaker individual effect sizes than lost loops. They're not robust new functional contacts — they're noise and redistribution.

### Why compartments change before loops

This seems backward — shouldn't the individual contacts (loops) change before the aggregate measure (compartments)? But compartments are defined by the average behavior of chromatin across very large regions. Even small, individually sub-significant changes in association preferences at many loci can sum to a significant compartment shift. Each individual locus might not have a significant loop change yet, but the collective tendency of Polycomb-marked regions to associate differently is already detectable. Loops need to cross a per-loop statistical threshold; compartments measure the mean of the entire distribution.

### The positive rho at Polycomb anchors

Already discussed above, but worth repeating because it comes up: the positive correlation between K119ub and loop strength at Polycomb sites is about which loops survive, not about K119ub helping. K119ub promotes local compaction (short-range, positive), but the long-range contacts that depend on higher-order Polycomb organization collapse (negative). Since the correlation includes all loops at Polycomb anchors and the short-range survivors outnumber the long-range losses in the remaining data, the aggregate rho is positive. It's a composition effect, not a mechanistic reversal.

### Why stripes don't change but loops do

Both involve physical contacts, but they're generated by fundamentally different mechanisms. Stripes arise from cohesin extruding chromatin until it hits a CTCF barrier — a mechanical, motor-driven process that doesn't care about histone modifications. Loops that depend on specific protein-protein interactions mediated by chromatin state (enhancer-promoter loops, Polycomb loops) are sensitive to changes in those marks. BAP1 loss changes chromatin state; it doesn't change the extrusion machinery. The stripes being unaffected is a nice internal control — it confirms the loop changes are specifically about chromatin state perturbation, not some artifact of library quality or general compaction changes.

### The K119ub_Only threshold

This is probably the most interesting biological finding from the ABC analysis for thinking about disease mechanism. The fact that K119ub accumulation alone produces real but sub-functional contact weakening means there's a buffer. Enhancers can tolerate some amount of K119ub accumulation without their target genes being affected, as long as the active marks stay. This predicts that the disease should be progressive — early perturbations are tolerated, and the phenotype only manifests when enough enhancers cross the threshold through secondary loss of active marks. That's consistent with the temporal progression from subtle early changes to massive late-stage dysregulation.

It also explains why the early timepoint is dominated by Polycomb-repressed sites rather than active enhancers. The active enhancers have this buffer; the Polycomb sites, which depend on K119ub cycling for their structural organization, don't. The buffered sites fail later.

### Connecting it all

The chain, from cause to consequence, goes roughly like this: BAP1 is lost, so K119ub cycling stops. At active enhancers, K119ub accumulates and weakens E-P contacts; where the weakening is severe enough and active marks are secondarily lost, target genes go down. At Polycomb domains, loss of cycling destabilizes long-range contacts; the interaction budget redistributes locally, creating short-range gains. The net effect is loop rewriting. The loop changes propagate to TAD boundaries (merges and strengthening) and sum to compartment shifts. Stripes, which depend on a different structural mechanism entirely, are preserved.

The temporal dimension adds that this starts small and focused (P13: few loops, Polycomb-specific) and amplifies into a system-wide remodeling by adulthood (thousands of loops, all chromatin types, massive compartment shifts), consistent with a progressive regulatory cascade rather than an acute structural collapse.

The abstract puts it well: BAP1 loss collapses long-range developmental loops and replaces them with more proximal contacts, resulting in dysregulation of synaptic and developmental genes. That's the one-sentence version.
