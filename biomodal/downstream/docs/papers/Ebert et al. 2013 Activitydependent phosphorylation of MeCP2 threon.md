## Summary

This 2013 *Nature* paper demonstrates that activity-dependent phosphorylation of MeCP2 at threonine 308 (T308) acts as a molecular switch controlling MeCP2's interaction with the NCoR co-repressor complex, and that disruption of this mechanism contributes to Rett syndrome (RTT).

---

### Background and Motivation

Rett syndrome is an X-linked neurodevelopmental disorder caused by mutations in MeCP2, a protein that binds methylated DNA across the genome and regulates transcription in neurons. Of the four most common RTT missense mutations, three (R106, R133, T158) sit in the DNA-binding domain and impair MeCP2's ability to latch onto methylated DNA. The fourth, R306C, falls within the repressor domain, but its mechanism of action was poorly understood. Prior work had identified one activity-dependent phosphorylation site on MeCP2 (S421), but knocking that site out in mice produced no detectable transcriptional changes, leaving the field without a clear link between neuronal activity, MeCP2 phosphorylation, and RTT pathology.

### Discovery of New Phosphorylation Sites

Using phosphotryptic mapping on radiolabeled neurons (a technique that catches phosphorylation events mass spectrometry can miss), the authors identified three new activity-dependent phosphorylation sites on MeCP2: S86, S274, and T308. All three sites, along with the previously known S421, become phosphorylated both in cultured neurons stimulated with KCl or bicuculline and in the brains of mice treated with kainic acid to induce seizures. Importantly, these sites respond differently to different upstream signals. S86 and S274 are preferentially triggered by BDNF and forskolin (a PKA activator), while T308 and S421 respond most strongly to membrane depolarization. This differential sensitivity positions MeCP2 as a convergence point for multiple signaling pathways, somewhat analogous to how histone modifications integrate diverse environmental cues to regulate chromatin state.

### T308 Phosphorylation Controls the MeCP2–NCoR Interaction

The authors zeroed in on T308 because of its physical proximity to the R306C RTT mutation. NCoR is a co-repressor complex that includes HDAC3 (a histone deacetylase) and is thought to mediate MeCP2's gene-silencing function. Using peptide pull-down assays, they showed that an unphosphorylated MeCP2 peptide spanning the T308 region successfully pulls down NCoR complex components (HDAC3, TBL1, TBLR1, GPS2) from brain lysates, but a phosphorylated version of the same peptide does not bind NCoR at all. Phosphomimetic mutations (T308D and T308E, which use acidic amino acids to partially imitate a phosphate group) showed intermediate, reduced binding. This establishes a clean mechanistic relationship: when T308 is phosphorylated, MeCP2 releases the NCoR repressor complex.

In a luciferase reporter assay, wild-type MeCP2 fused to a GAL4 DNA-binding domain effectively silenced a reporter gene, while the R306C mutant (which can't bind NCoR) completely lost repressive ability. The phosphomimetic T308D/E mutants showed partial loss of repression, matching their partial loss of NCoR binding. The T308A mutant (which can't be phosphorylated but still binds NCoR normally) retained full repressive capacity. This confirms that T308 phosphorylation weakens MeCP2-mediated transcriptional repression specifically by disrupting the NCoR interaction.

### Functional Consequences in T308A Knock-In Mice

To test what happens when neurons can never release NCoR from MeCP2 in an activity-dependent manner, the authors generated T308A knock-in mice. These mice express normal levels of MeCP2 protein, bind DNA normally, and maintain baseline NCoR interaction, so any phenotype can be attributed specifically to the loss of phosphorylation-regulated NCoR release.

In cultured neurons from these mice, membrane depolarization induced most activity-dependent genes (Arc, Fos, Nptx2, Adcyap1) normally, but *Npas4* induction was significantly reduced. *Bdnf* showed a similar trend. This same selective deficit appeared *in vivo* when dark-reared mice were exposed to light: T308A knock-in mice had attenuated *Npas4* and *Bdnf* induction in their visual cortex. This is particularly relevant because NPAS4 promotes inhibitory synapse development on excitatory neurons, and disrupted excitatory/inhibitory balance is a hallmark of RTT.

The T308A knock-in mice also displayed RTT-like neurological phenotypes: reduced brain weight, hindlimb clasping, impaired rotarod performance (indicating motor deficits), and a lower seizure threshold when challenged with PTZ.

### Connection to RTT Through R306C

A critical finding ties everything together. When R306C knock-in mice were treated with kainic acid, MeCP2 T308 phosphorylation failed to occur, even though the pT308 antibody can still recognize phosphorylated T308 in the presence of the R306C mutation. This makes mechanistic sense: the kinases that phosphorylate T308 (CaMKIV and PKA) require a basic residue two to three amino acids upstream of the phosphorylation site, and R306 is that residue. Mutating it to cysteine eliminates the kinase recognition motif.

This means the R306C RTT mutation has a dual effect: it both abolishes baseline NCoR binding *and* prevents activity-dependent phosphorylation at T308. The authors propose that R306C knock-in mice (which can't bind NCoR at all) phenocopy MeCP2 loss-of-function, while T308A knock-in mice (which bind NCoR constitutively and can never release it) phenocopy a gain-of-function state where MeCP2 is a permanently active repressor. Both extremes produce RTT-like symptoms, paralleling the established observation that both MeCP2 deletion and MeCP2 overexpression cause neurological disease. Proper MeCP2 function requires not just NCoR binding, but *regulated* NCoR binding that can be dynamically toggled by neuronal activity.

### Significance

The paper establishes a model where MeCP2 functions as an activity-responsive epigenetic regulator: neurons fire, calcium-dependent kinases phosphorylate T308, the NCoR repressor complex disengages from MeCP2, and specific activity-dependent genes (particularly *Npas4* and *Bdnf*) become derepressed. When this dynamic regulation is disrupted, whether by preventing phosphorylation (T308A) or by destroying NCoR binding entirely (R306C), the result is aberrant gene expression programs and RTT-like neurological dysfunction.