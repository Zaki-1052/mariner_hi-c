# Jesse Dixon Meeting Summary — 04/10/26

## Overview

We presented our BAP1 Hi-C paper to Jesse Dixon at the Salk. The meeting covered the full scope of the second paper: epigenetic intro (cut-and-run, EvoC, MeCP2), then the bulk of the Hi-C analysis — TAD interaction density, boundaries, loops, loop annotations, loop lengths, pileups, compartments, and the activity-by-contact stuff. Jesse was generally positive — he thinks the magnitude of 3D genome changes is surprisingly large given the histology looks normal, and that this is enough to tell a story even without a full mechanism. He validated most of our analytical choices and gave us a bunch of concrete suggestions.

## Epigenetics Recap (for Jesse's context)

We walked him through the first paper's findings: BAP1-KO causes genome-wide increase in H2AK119ub, decrease in K27ac, increase in K27me3 at euchromatic loci (heterochromatinization), and the inverse at heterochromatic loci (euchromatinization). The MYC locus is a good example of the heterochromatin → euchromatin switch. Transcriptional changes track with histone modifications — overexpressed genes have more K27ac, downregulated genes have less K27ac and more K27me3. Dysregulated genes are mostly neuronal/synaptic, especially the upregulated ones.

The MeCP2 cut-and-run data is interesting — Jesse pointed out that the peaks we're seeing probably reflect non-CG methylation rather than CG methylation (which is too uniform to produce peaks). We confirmed this lines up almost perfectly with non-CG methylation, which is the next thing to test. EvoC shows methyl-C and hydroxymethyl-C are inversely correlated, again breaking down along euchromatin/heterochromatin lines.

## Hi-C: TADs & Boundaries

No differential TAD interaction densities at P13 (Homer), but ~20–25% differential TADs in adult (166 decreased, 777 increased out of ~4,366). Jesse thought this looked normal and suggested correlating differential TADs with ubiquitination status — the relationship between interaction density and compaction vs. activity isn't straightforward.

TADcompare identifies differential boundaries at both timepoints (~20%), which is weird because Homer finds no differential TAD IDs at P13. Jesse wasn't sure what to make of this — could be statistical differences between tools. He noted that some of the insulation score profiles look unusual: instead of high flanking → strong dip (normal boundary), some show low flanking → rise → dip, which might be sub-TAD structure rather than real boundaries.

He suggested the boundary phenotype might actually be stripes/flares rather than true boundary changes — directional cohesin loading from CTCF/gene/super-enhancer sites. This would explain the association with gene expression changes (gained boundaries = upregulated DEGs, lost boundaries = downregulated DEGs, significant in adult p < 0.00001 but not P13).

## Hi-C: Loops

Our cutoffs (logFC > 0.3, FDR < 0.05) are reasonable — Jesse said fold changes of 3–4 are rare in Hi-C, and 1.3–1.5 fold is standard. He emphasized visual inspection as the gold standard ("look at the data before you do 20 other analyses"). P13: 195 lost, 166 gained (out of ~39K). Adult: 4,253 lost, 3,728 gained (out of ~39K). APA plots look clean and believable.

**HICCUPS caveat:** If you call loops on one condition, that condition will artifactually look stronger due to peak pixel selection bias. Mariner should correct for this, which Jesse confirmed is fine.

**Loop annotations (using control peaks):** Gained loops are enriched for Polycomb/CTCF. Lost loops are enriched for active enhancers. Jesse acknowledged the challenge of annotating with control chromatin state when the state itself is changing — no clear right answer. He thinks post-translational modifications are primary, driving the conformational changes (Erez agreed with this interpretation).

**Key finding — loop type distribution (slide 13):** Most lost loops are active enhancer–active enhancer. Most gained are CTCF–CTCF, poised enhancer–CTCF, Polycomb–Polycomb. The stats back this up: 92.4% of K27ac-anchored loops are lost, 74.6% of K27me3-anchored loops are gained, 82.7% of K4me3-anchored loops are lost.

## Loop Length — The Big Finding

Lost loops are systematically longer than gained loops (median 625 kb vs. 320 kb, KS p = 5.49e-48). This is even more dramatic for K27me3-anchored loops: one+ anchor median lost = 1,560 kb vs. gained = 325 kb; both anchors (Polycomb–Polycomb) median lost = 3,238 kb vs. gained = 325 kb. As lost loops get longer, they become more K27me3-associated and less K27ac-associated.

Jesse's interpretation: lost loops may break into two classes — shorter CRE-type/enhancer–promoter loops (which are naturally shorter and K27ac-associated) and longer architectural/Polycomb loops. Gained loops are mostly Polycomb, which are all shorter-range. He couldn't fully explain why K27me3 trends opposite in lost loops as a function of distance, or why it's not a clean reciprocal in gained loops.

**CTCF–CTCF loops specifically** (slide 20): Lost CTCF loops are longer than gained (median 810 kb vs. 325 kb, p = 8.92e-37), but K27ac-anchored loops show no length difference (p = 0.41). This suggests a possible extrusion phenotype — Polycomb may be impeding loop extrusion, preventing longer CTCF loops from forming.

Jesse connected this to the **Brad Bernstein IDH mutant glioma work** — repressive chromatin gains knock out CTCF binding sites. His question: is Polycomb repression happening at the anchors specifically, or across the loop body? If body → impediment to extrusion. If anchors → those CTCF sites are specifically sensitive to repressive state changes. He also wondered about CpG island / GC content differences at affected anchors.

## Tessa Popay — Anchor vs. Span Analysis

Jesse's biggest concrete suggestion: look at the body/span of loops, not just anchors. His postdoc **Tessa Popay** (P-O-P-A-Y) did exactly this for cohesin depletion — metagene-style ChIP-seq plots from loop anchors showed that the real signal was K27me3 enrichment *flanking* CTCF sites, suggesting the CTCF loops are flanking repressive regions. This would explain why CTCF–CTCF gained loops are enriched: they're weak loops flanking heterochromatin that opens up when BAP1 is lost. Jesse will email Tessa and CC us — she's interested in Polycomb stability and would probably share code.

## Pileups, Permutations, Compartments

Pileup analyses (slides 21–24) show the expected: increased accessibility/K27ac → increased contacts, decreased → decreased. K27me3 is opposite. K119ub effect is modest. Permutation tests (slides 25–28) confirm: K27ac and ATAC correlate with boundary strengthening, K119ub is opposite, but K27me3 is enriched at both strengthened and weakened boundaries (matching the two-class loop loss phenomenon).

**Differential TADs in A/B compartment (slide 26):** Increased TAD ID → B compartment (heterochromatin gaining density). Decreased TAD ID → A compartment (euchromatin losing density). Fits the overall model.

**Compartment switches (slide 29):** We have differential PC1 at many 25 kb bins, but not many contiguous flips. Jesse was surprised — expected more full switches given the magnitude of epigenetic changes. Suggested it's quantitative weakening rather than full flips. For plotting: do PC1 mutant vs. PC1 WT scatter, not volcano — a change of ±1 in PC1 doesn't necessarily mean a flip.

**Subcompartment tools:** Jesse recommended **SNIPER** (deep learning) over DCIC for calling subcompartments (A1, A2, B1, B2 etc.). We tried DCIC and results were iffy. SNIPER needs ~1B reads — Zakir confirmed we have that. The idea: maybe we're not switching A↔B but rather B1↔B2 or A1↔A2, which have different histone modification associations (K79me, H4K20me, etc.). Compare subcompartment calls against our histone mod data in controls to validate.

## Activity-by-Contact

Our ABC analysis (slide 30) shows a decent correlation (Pearson r = 0.537, Spearman ρ = 0.660 for DE genes) between summed ΔActivity×Contact and gene expression changes. Per-enhancer K119ub fold change shows a negative trend with ΔA×C (ρ = -0.248). Jesse thinks this is real signal — our long-term chromatin state change is probably more predictive than acute cohesin depletion experiments (where ABC-expression correlation is ~0). But he noted 90% of the per-enhancer data is a tight blob near zero — would look better if we filter to only differentially ubiquitinated enhancers first.

Slide 32 shows that at the per-enhancer level, contact changes are minimal (BAP1 loss doesn't cause drastic enhancer–TSS contact changes detectable by standard Hi-C), but activity changes are clear. Slide 33: loop strength is decreased at K119ub-only enhancers but not ΔA×C — so ub alone affects loop strength without necessarily changing the ABC score.

## Big Picture — Is This a Story?

Jesse thinks yes. The magnitude of Hi-C changes is large and surprising given normal histology. The chromatin state changes predict the 3D changes. The fact that BAP1 contributes this much to higher-order structure is novel. He mentioned that even mediator knockouts don't show this level of change acutely. The progressive nature (P13 → adult) adds developmental interest. He doesn't think the burden for mechanism is that high for this type of paper.

His one remaining biological puzzle: BAP1 is more highly expressed early and less later, but the phenotype is stronger later. We explained that ubiquitination starts in heterochromatin and shifts to euchromatin over development, and BAP1 tracks with that. Still not fully satisfying as an explanation, but Jesse didn't seem to think it's a dealbreaker.

---

# Action Items

## High Priority (for paper)

- **Loop body/span analysis (Tessa Popay style)**
  - Make metagene-style ChIP-seq plots from loop anchors — look at K27me3, K27ac, K119ub signal across anchor *and flanking regions*
  - Key question: is Polycomb enrichment at the anchors or across the body of gained CTCF–CTCF loops?
  - Body enrichment → extrusion impediment model. Anchor-specific → sensitivity model.
  - Email Tessa Popay (Jesse is CC'ing us) to get her code/approach
  - Relevant paper from Tessa's work on cohesin depletion + Polycomb loop persistence

- **StripENN analysis**
  - Re-run stripe calling with StripENN (not the other stripe caller we used)
  - https://www.nature.com/articles/s41467-022-29258-9
  - Jesse thinks our "differential boundaries" might actually be stripes/flares — directional cohesin loading
  - Check: are differential boundary-proximal DEGs enriched at stripe-like features?
  - If stripes explain the boundary–expression association, that's more interpretable than boundary strength alone

- **Insulation score investigation**
  - Jesse flagged unusual insulation profiles — some boundaries show low flanking → rise → dip instead of normal high → dip
  - Use **Cooltools** pileup analysis at differential boundaries (we did Genova, but Cooltools is more standard)
  - Browse specific loci in the genome browser to visually check what's happening at called "boundaries"
  - Determine if TADcompare is calling sub-TAD structure rather than real boundaries

- **ABC analysis cleanup**
  - Filter per-enhancer ABC plot to only sites with significantly differential K119ub
  - Also plot activity vs. contact separately (not multiplied) to see if effect is mostly activity-driven
  - Will make the scatter plot much cleaner (currently 90% blob at zero)

## Medium Priority

- **Compartment analysis improvements**
  - Plot PC1 mut vs. PC1 WT as scatter (not volcano) to properly visualize flips vs. quantitative weakening
  - Try **SNIPER** for subcompartment calling (we have ~1B reads, should be enough)
    - https://www.nature.com/articles/s41467-019-12954-4
    - https://github.com/ma-compbio/SNIPER
  - Validate subcompartment calls against control histone modification data
  - Look for B1↔B2 or A1↔A2 type switches rather than full A↔B flips

- **Correlate differential TADs with K119ub status**
  - Jesse suggested this directly — which TADs gaining density are at ub-enriched regions vs. not?
  - Unclear whether increased density = compaction (more repressive) or more loops (more active)

- **GC content / CpG island analysis at CTCF anchors**
  - Are the CTCF sites that lose long-range loops in CpG island regions?
  - Connects to Brad Bernstein IDH mutant work on CG methylation disrupting CTCF
  - **Check this paper** (Brad Bernstein repressor phenotype / IDH mutants / CTCF)

- **Developmental comparison**
  - Compare adult mutant Hi-C to P12 wild-type Hi-C — does the mutant adult look more like immature neurons?
  - Could support a "blocked developmental remodeling" narrative

## Lower Priority / Future

- **CTCF/Rad21 CUT&RUN** — we're missing this (noted on slide 20). Would directly test whether CTCF binding changes at affected loop anchors.
- **Non-CG methylation** — MeCP2 binding data predicts we'll see phenotypes. EvoC can't measure this. Need different method.
- **P12 histology** — Jesse asked, we haven't done it. Not urgent if adult cerebellum looks normal.
- **EZH2 inhibitor comparison** — Jesse mentioned they did EZH2 inhibitor experiments and saw less Hi-C change than we're seeing with BAP1-KO. Interesting comparison point for discussion section.
