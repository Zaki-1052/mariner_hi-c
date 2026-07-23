## Summary

This 2014 paper demonstrates that non-CpG methylation (called "CpH" methylation, where H = A, C, or T) is a pervasive and functionally significant epigenetic mark in mammalian neurons, challenging the long-standing assumption that DNA methylation in animal genomes is essentially restricted to CpG dinucleotides.

---

### The Core Discovery

The authors generated a single base–resolution DNA methylome from mouse dentate gyrus granule neurons and found that roughly 25% of all methylated cytosines occur in a non-CpG context. This breaks down into about 21% mCHH and 4% mCHG. Importantly, this isn't some low-level background noise — it's a substantial fraction of the total methylation landscape, and it's neuron-specific. The same loci showed virtually no CpH methylation in spleen tissue, establishing this as a brain-enriched phenomenon.

### Where It Lives in the Genome

CpH methylation has a distinct genomic signature compared to CpG methylation. It shows a strong sequence preference for the CAC motif, tends to cluster in regions of low CpG density, and displays intermediate methylation levels (around 25%) rather than the all-or-nothing bimodal pattern typical of CpG methylation. This intermediate level suggests that at any given CpH site, only a fraction of alleles in the cell population carry the methyl mark at steady state — a reflection of its dynamic turnover rather than stable maintenance.

The spacing pattern of neighboring mCpH sites shows an 8-bp periodicity (reminiscent of the DNMT3A-DNMT3L complex footprint) and a broader ~180-bp periodicity tied to nucleosome positioning, suggesting that the physical architecture of chromatin directly shapes where CpH methylation lands.

### Evolutionary Conservation

The pattern isn't mouse-specific. Sanger bisulfite sequencing of orthologous regions in human brain tissue confirmed conserved CpH methylation at the same loci, and about 83% of human CpH-methylated genes had orthologs that were also CpH-methylated in mouse. This conservation across tens of millions of years of evolutionary divergence strongly implies functional importance.

### Functional Consequences

Three lines of evidence establish that CpH methylation represses gene transcription:

First, genome-wide, CpH methylation levels are anticorrelated with expression of nearby genes — and this holds even for mCpH sites that have no neighboring CpGs within 500 bp, ruling out the possibility that CpG methylation is doing all the work while CpH is just a bystander.

Second, in a reporter assay, placing CpH methylation on a GFP-expressing plasmid repressed transcription in cultured neurons just as effectively as CpG methylation at the same density. The repression wasn't indirect (i.e., CpH methylation didn't trigger secondary CpG methylation on the plasmid).

Third, knocking down DNMT3A in vivo preferentially reduced CpH methylation while leaving CpG methylation largely intact, and this was accompanied by derepression (increased expression) of the affected genes.

### The MeCP2 Connection

MeCP2, the well-known methylated-DNA reader protein whose mutations cause Rett syndrome, was shown to bind CpH-methylated DNA both in vitro (via EMSA) and in vivo (via ChIP-bisulfite sequencing). It binds mCpH with lower affinity than mCpG alone, but when both marks coexist, binding is dramatically enhanced. This is particularly interesting because the postnatal onset of Rett syndrome coincides temporally with the emergence of CpH methylation during neuronal maturation — not with CpG methylation, which is established much earlier in development.

### Developmental Timing and Maintenance

Unlike CpG methylation, which is largely set up during early embryonic development and then maintained by DNMT1 through cell division, CpH methylation is established _de novo_ postnatally as neurons mature. It accumulates gradually during postnatal brain development and requires continuous DNMT3A activity for its maintenance. This makes mechanistic sense: CpH methylation is inherently asymmetric (a methylated CAC on one strand doesn't have a complementary methylated site on the opposite strand), so DNMT1's maintenance mechanism — which relies on recognizing hemimethylated CpG palindromes after replication — simply can't maintain it. Every CpH methyl mark that gets removed must be re-established from scratch by DNMT3A.

This dependency on ongoing _de novo_ methylation also explains the higher turnover and intermediate methylation levels observed — the marks are in a constant state of being removed and re-added, and at any snapshot in time, only a subset of alleles will be caught in the methylated state.

### Broader Significance

Since CpG dinucleotides represent only about 4% of the genome, CpH methylation massively expands the proportion of the neuronal genome under cytosine methylation regulation. It's particularly enriched in CpG-depleted regions, suggesting it may compensate for the absence of CpG methylation in those areas. The authors frame this as an entirely new layer of epigenetic control in neurons, with implications for understanding neuronal identity, synaptic plasticity, and neurological disorders — especially Rett syndrome and other MeCP2-related conditions.