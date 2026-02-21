# Abstract Rough Draft

**Title:** Regulation of chromatin conformation by the histone deubiquitinase BAP1 in the brain

---

Gene expression is not solely determined by DNA sequence -- it is also regulated by how DNA is physically organized within the nucleus. Chromosomes fold into three-dimensional (3D) structures, forming loops that bring distant regulatory elements called enhancers into contact with the genes they control. Disruption of this 3D architecture is increasingly linked to developmental disorders and neurodegeneration, yet the molecular mechanisms governing chromatin folding remain poorly understood.

BAP1 is an enzyme that removes a specific chemical modification -- monoubiquitination of histone H2A at lysine 119 (H2AK119ub) -- from the proteins that package DNA. This modification, deposited by the Polycomb repressive complex, is associated with gene silencing. Loss of BAP1 in the mouse brain leads to neurodegeneration, but how BAP1 contributes to 3D chromatin organization has not been characterized.

To address this question, we performed Hi-C chromosome conformation capture sequencing in BAP1-knockout and wildtype mouse cerebellum at two developmental timepoints (postnatal day 12 and adult), with three biological replicates per condition. We developed computational pipelines to analyze chromatin architecture at multiple scales -- loops, topologically associating domains (TADs), stripes, and compartments -- integrating data from five histone modifications, chromatin accessibility, gene expression, and enhancer-gene linkage predictions.

We identified 2,910 differential chromatin loops in the adult cerebellum, revealing a striking "loop rewriting" phenomenon: long-range loops spanning over one megabase were preferentially lost (3.3-fold enrichment for loss), while shorter-range contacts were gained. Accumulation of H2AK119ub at loop anchors strongly predicted loop loss (odds ratio = 10.7, p < 10^-91). Integration with an enhancer-gene linkage model showed 88% concordance between loop changes, altered enhancer-gene connections, and differential gene expression at the same loci. At the broadest scale, 44% of the genome exhibited shifts in compartmentalization, while TAD boundaries remained relatively stable, and stripe architecture was preserved. Critically, this effect was progressive -- only 165 differential loops were detected at postnatal day 12, expanding to 2,910 by adulthood, representing an 18-fold amplification.

These findings establish a multi-scale model of BAP1-mediated chromatin regulation in which the progressive accumulation of H2AK119ub destabilizes long-range developmental contacts, fundamentally reorganizing the regulatory landscape of the brain. This work provides new insight into how epigenetic dysregulation translates into 3D genome disruption, with implications for understanding neurodegenerative disease.

---

**Word count:** ~365 words
