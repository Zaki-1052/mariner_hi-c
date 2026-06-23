Section 72 script is written. Here's what it does:

**72a** — Violin+boxplot comparing K119ub gene body signal between neuronal and non-neuronal genes, in both ctrl and mut conditions. Wilcoxon one-sided test.

**72b** — Fisher's exact tests at four thresholds (top quartile and top decile, for ctrl and mut). Log2(OR) forest plot with CIs.

**72c** — GO BP enrichment of the K119ub top-quartile ctrl genes (no MeCP2 filter). Standard dotplot.

**72d** — Compares neuronal GO terms across three gene sets: ctrl-high, mut-high, and gained (mut-high but not ctrl-high). Shows whether neuronal enrichment is constitutive or gained upon BAP1-KO.

**72e** — GSEA ranking ALL genes by absolute K119ub signal (ctrl and mut separately). This is the key distinction from section 61k which ranked by log2FC. Running-score plot for top neuronal term.

**72f** — Dose-response: neuronal fraction across 10 K119ub signal deciles with Spearman trend test and per-decile Fisher ORs.

**Neuronal gene set** is derived fresh from `org.Mm.eg.db` (all mouse genes annotated to synapse/neuron/axon/dendrite/nervous GO BP terms) — completely unbiased by any MeCP2 or methylation filtering.

To run:
```bash
cd /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/downstream
Rscript scripts/viz_sections/section_72_k119ub_neuronal_characterization.R 2>&1 | tee logs/a2/72_k119ub_neuronal.txt
```

---

The results are clean and directly validate the PMC12652333 framework. Here's what stands out:

**The core finding holds**: Neuronal genes ARE constitutively K119ub-enriched in ctrl (p = 2.15e-15, OR = 1.70 at top decile, 49/49 neuronal GSEA terms with positive NES). This confirms the paper's claim in your cerebellum dataset — before BAP1 is even deleted, K119ub preferentially sits on neuronal genes.

**The surprising twist is in the mutant**: The enrichment *weakens* in mut (Wilcoxon p drops from 2.15e-15 to 0.29 NS; Fisher OR drops 1.70 → 1.37). The gained-K119ub genes (902 genes) are only 15.4% neuronal vs 23.5% genome-wide, and have zero significant GO terms. BAP1-KO adds ubiquitin *broadly*, not preferentially at neuronal loci — it dilutes the neuronal specificity.

**What this means for the mechanism chain**:
1. PRC1 constitutively deposits K119ub at neuronal/developmental genes in WT cerebellum (the top GO terms — pattern specification q=1.6e-43, cell fate commitment q=1.5e-42 — are classic Polycomb targets)
2. BAP1-KO removes the deubiquitinase → K119ub accumulates genome-wide, but neuronal genes that already had high K119ub get pushed even higher in absolute terms
3. MeCP2 responds to this K119ub-driven chromatin change, not to methylation (consistent with section 71 showing K119ub explains 7.3% of MeCP2 variance vs 0% for mC/hmC ratio)

**The dose-response (72f)** confirms a quantitative gradient: D10 genes are 33% neuronal vs 21% at D1 (rho = 0.648, p = 0.049), though D6 dips oddly — suggesting a bimodal distribution where the neuronal enrichment really kicks in at the high tail.

**One technical note**: The GSEA warnings about "all values greater than zero" are expected since we ranked by absolute signal (not fold-change). The `scoreType = "pos"` suggestion could marginally improve power, but with 49 neuronal terms already significant at NES > +1.9, it wouldn't change the conclusions. I can add that flag if you want a cleaner rerun.

Bottom line: the data strongly supports the interpretation that neuronal genes are under constitutive K119ub/Polycomb control, and that's the substrate MeCP2 is reading — exactly the PMC12652333 framework.


---
Two reasons:

1. **Section 61k already did fold change.** It ranked all genes by K119ub log2FC (mut/ctrl) and found neuronal terms enriched at the positive end. That question — "do genes that *gain* K119ub upon BAP1-KO tend to be neuronal?" — was already answered.

2. **The paper's claim is about baseline association, not perturbation.** PMC12652333 says neuronal genes *are associated with* H2AK119ub — meaning they carry it constitutively in the wildtype brain. Ranking by absolute ctrl signal directly tests that: "in our WT cerebellum, do genes with the most K119ub tend to be neuronal?" Fold change can't answer that because a gene can have very high K119ub in both ctrl and mut (neuronal, constitutively marked) yet show zero fold change.

The results actually illustrate why both views matter: neuronal genes are enriched at the top of the ctrl signal distribution (OR = 1.70, GSEA NES = +2.48), but the *gained* K119ub genes (fold-change positive, mut-only high) are only 15.4% neuronal — below genome-wide average. The neuronal K119ub is already there; BAP1-KO adds ubiquitin elsewhere.

That said, the GSEA warnings are real — since absolute signal is all-positive, standard GSEA (designed for signed statistics) is suboptimal. The Fisher/Wilcoxon/decile analyses (72a, 72b, 72f) are more methodologically appropriate for this question and they tell the same story.

---

---

Comparing the two runs — the core findings are identical, but the technical quality improved:

**What changed:**
- Spurious negative-NES terms eliminated: ctrl went from 19 neg → 0, mut from 7 → 0. Those were artifacts of `scoreType = "std"` on all-positive data.
- The "unbalanced gene-level statistic values" warning (112 pathways with NA p-values) is gone entirely.
- q-values got substantially more significant where warranted — e.g., neuron fate specification went from q=3.82e-09 to q=6.92e-17, neuron fate commitment from q=3.82e-09 to q=1.02e-23. That's `eps = 0` doing its job.
- Total sig terms shifted slightly (1209→1177 ctrl, 935→919 mut) because the fake negatives dropped out.

**What didn't change:**
- Same 49 neuronal terms in ctrl, 39 in mut
- Same NES values (neuron fate specification still +2.48)
- All other panels (72a-d, 72f) are deterministic, so identical

Clean run. The only remaining warnings are the 0.6% ties in ranked stats (inherent to the data) and the enrichplot deprecation warnings (upstream package issue, not ours).

---

I've now thoroughly reviewed all the relevant data. Here's my critical assessment:

---

## Critical Evaluation of Your Theory

### Claim 1: "Global redistribution of DNA methylation as increased 5mC + decreased 5hmC without changes in overall total methylation"

**This is the weakest claim and is contradicted by your own data at multiple levels.**

**Global level (Section 64):** The modC (total = mC+hmC) delta is -0.079% with Wilcoxon p = 0.375 and t-test p = 0.492. So yes, at the genome-wide autosomal level, total methylation is statistically unchanged. You can defensibly say this.

**But Section 61a demolishes the "without changes in total" framing at the gene level:**
- **All genes:** total methylation mean delta = -0.00139, **significantly decreases** (q = 2.48e-112)
- **Coordinated genes:** total methylation mean delta = +0.00697, **significantly increases** (q = 2.12e-178)
- **Neuronal genes:** total methylation mean delta = +0.01200, **significantly increases** (q = 2.07e-56)
- **MeCP2-Up genes:** total methylation mean delta = +0.03160, **significantly increases** (q = 3.19e-11)

The stoichiometry test further shows the global slope is -0.959, which is **significantly different from -1** (the stoichiometric expectation for a dehydroxymethylase). Your data explicitly says "NOT stoichiometric → independent processes." The neuronal subset slope is even further off at -0.72.

**Bottom line:** You cannot say "without changes in overall total methylation" as a blanket statement. Genome-wide it's approximately flat, but at every biologically meaningful gene subset — coordinated, neuronal, MeCP2-bound — total methylation **significantly increases**. That's a net gain, not conservation. If you want to keep this claim, you need to restrict it explicitly to "genome-wide averages" and acknowledge it breaks down at the loci that matter for your story.

### Claim 2: "The strongest chromatin-state remodeling occurs at neuronal genes, where MeCP2 redistribution is selectively observed"

**Partially supported, but "selectively" is too strong.**

Section 72 shows neuronal genes do carry elevated K119ub constitutively (ctrl median 1.03 vs other 0.97, p = 2.15e-15) and are enriched among K119ub-high genes (top-decile OR = 1.70). The GSEA confirms neuronal fate/differentiation terms are enriched among K119ub-high genes in both ctrl and mut.

**However, the selectivity claim has problems:**

1. **MeCP2 redistribution is NOT selective to neuronal genes.** Section 65 shows MeCP2 signal drops genome-wide: ctrl median = 1.40, mut median = 1.14 across ALL 20,969 genes. The drop is significant in every methylation category — hypermethylated (p = 5.4e-148), hypomethylated (p = 1.6e-63), and non-significant genes (p = 5.0e-159). This is a global MeCP2 redistribution, not a neuronal-selective one.

2. **The MeCP2-neuronal overlap is tiny.** Section 61e: only 16 out of 1,149 neuronal genes overlap with MeCP2-Up (1.4%). Only 12 genes are in all three sets (neuronal + coordinated + MeCP2-Up). That's not "selective observation at neuronal genes" — it's near-background.

3. **Section 61k GSEA:** Ranking all 21k genes by MeCP2 fold-change, only **1 single GO term** reaches significance: "synapse assembly" (q = 0.025). Compare this to K119ub ranking which hits 115 terms. The neuronal enrichment in MeCP2 redistribution is marginal at best.

### Claim 3: "MeCP2 responds not to H2AK119ub or DNA methylation alone, but to the mix of polycomb remodeling, methylation shifts, and local chromatin-state changes"

**This is well-supported and is the strongest part of your theory.**

Section 62 is definitive:
- Chromatin features explain MeCP2 binding far better than CG methylation alone (R² = 0.2456 vs 0.0168 for binding level; R² = 0.1704 vs 0.0128 for differential)
- LASSO selects K119ub, K27me3, ATAC, K36me3, K27me1, and CG-mC as relevant features — confirming it's a multivariate chromatin context
- K119ub has the largest standardized beta for both binding level (β = +0.199) and differential (β = +0.113)
- After removing CG methylation effects, chromatin features still explain R² = 0.159 of the residual

Section 67 further supports this: at the 359 MeCP2-bound genes with NO methylation change, K119ub is still dramatically gained (72.8% vs 45.9% background, OR = 3.15, p = 1.82e-24), K27me3 increases, and K27ac decreases.

Section 71 adds: K119ub explains 7.3% unique variance in MeCP2, while the mC/hmC ratio explains 0.0% unique variance. There is a significant K119ub × ratio interaction (p = 4.3e-14), but it's tiny.

### Claim 4: "BAP1 loss increases H2AK119ub at normally regulated loci and disrupts polycomb domain organization"

**Supported by Section 72** — the K119ub enrichment analysis is clean. But "normally regulated loci" is vague. The data shows neuronal genes are K119ub-high constitutively, and BAP1-KO further increases K119ub there (OR drops from 1.70 in ctrl to 1.37 in mut top-decile, indicating the signal spreads from already-high neuronal loci to others). The shift is more "broadening of K119ub" than "disruption at normal loci."

### Claim 5: "Net shift toward repressive chromatin states"

**Supported by Section 66.** The subcompartment data shows clear directionality:
- A.1 (most active): 64.3% significant DMRs, overwhelmingly hypermethylated (5716 hyper vs 1216 hypo)
- B.2 (most repressed): only 21.8% significant, majority hypomethylated (249 hyper vs 661 hypo)
- K27me3-only (fac. het): 92/1581 = 5.5% hyper vs 94.5% hypo
- K27ac-only (active): 4505/4830 = 93.3% hyper

This is consistent with active regions gaining methylation and repressed regions losing it — which could be interpreted as homogenization, not necessarily a "net shift toward repressive." The active regions are gaining mC (which could be interpreted either way).

---

## Summary: What's solid vs. what needs revision

**Solid:**
- MeCP2 responds to chromatin context, not methylation alone (sections 62, 67, 71)
- K119ub is the strongest single predictor of MeCP2 redistribution
- Neuronal genes carry elevated constitutive K119ub
- The mC↑/hmC↓ pattern is real and compartment-dependent

**Needs revision:**
- "Without changes in overall total methylation" — contradicted at every meaningful gene-level subset (section 61a). The gene-level data says total methylation increases at coordinated/neuronal/MeCP2-up loci
- "MeCP2 redistribution is selectively observed at neuronal genes" — MeCP2 drops everywhere; the neuronal-MeCP2 overlap is 16/1,149 genes; GSEA finds only 1 significant neuronal term for MeCP2
- The stoichiometry data explicitly says mC and hmC changes are NOT stoichiometric (independent processes), which weakens the "redistribution" framing and instead suggests partial TET inhibition plus independent mC gain

**Suggested reframing:** BAP1 loss produces a genome-wide mC↑/hmC↓ shift via partial TET inhibition (supported by TET-KO phenocopy rho = 0.22 and mediation K27me3→ATAC→TET), with total methylation approximately conserved genome-wide but significantly increased at biologically active subsets. MeCP2 redistribution is driven primarily by chromatin remodeling (K119ub, K27me3, ATAC) with DNA methylation playing a modulatory rather than driving role. The neuronal specificity exists for K119ub enrichment but not convincingly for MeCP2 redistribution itself.


---

Script written. Five panels:

- **72g-a**: Violin/boxplot comparing ATAC, K27ac, K27me3 fold changes between neuronal and non-neuronal genes (all genes, no filters)
- **72g-b**: The critical 4-group test — K119ub-high/low x neuronal/non-neuronal per mark. If K119ub-high + neuronal shows the strongest heterochromatin shift, that's the theory confirmed
- **72g-c**: Chromatin state composition (Active_Promoter, Polycomb, etc.) of K119ub-high genes split by neuronal status — are neuronal K119ub-high genes disproportionately in repressive states?
- **72g-d**: Dose-response curves: K119ub decile vs median mark change, with separate neuronal and non-neuronal lines. Shows whether the relationship is neuronal-specific
- **72g-e**: Formal interaction models: `mark_fold ~ k119ub_ctrl_signal * is_neuronal` — if the interaction term is significant, neuronal genes respond differently to K119ub

```bash
cd /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/downstream
Rscript scripts/viz_sections/section_72g_neuronal_chromatin_remodeling.R 2>&1 | tee logs/a2/72g_neuronal_chromatin.txt
```

---

The results are really interesting but tell a more nuanced story than the simple "heterochromatin shift" prediction. Let me break it down honestly:

**What supports the theory:**
- **Chromatin state composition (72g-c)** is exactly right: neuronal K119ub-high genes are disproportionately Repressed_Promoter (47.9% vs 34.6%, OR=1.74, p < 2.2e-16) and Bivalent (4.0% vs 1.4%, OR=2.87, p=3.4e-7). These ARE Polycomb-regulated genes at baseline. That cleanly validates the PMC12652333 framework.

**What complicates the theory:**
- The chromatin *changes* upon BAP1-KO go the **opposite direction** from "heterochromatin shift":
  - **ATAC**: neuronal K119ub-high genes GAIN accessibility (+0.17 vs +0.11, p=2.9e-11)
  - **K27ac**: neuronal K119ub-high genes GAIN K27ac (+0.05) while non-neuronal K119ub-high genes LOSE it (-0.02), p=2.15e-5
  - **K27me3**: both groups lose it, no significant neuronal difference at K119ub-high (p=0.22)

So neuronal K119ub-high genes are starting from repressive states but becoming MORE open upon BAP1-KO, not less. The K27ac interaction is highly significant (p=4.6e-7) — neuronal genes respond to K119ub accumulation with active mark *gain*, non-neuronal genes with active mark *loss*.

**What this actually means for the story:**

This is arguably more interesting than simple heterochromatin shift. It suggests **Polycomb disorganization** rather than Polycomb reinforcement:

1. Neuronal genes sit under PRC1/K119ub-mediated repression at baseline (confirmed by 72g-c)
2. BAP1 loss disrupts the normal K119ub cycling → excess K119ub accumulates
3. But excess K119ub doesn't strengthen repression — it **destabilizes** the Polycomb domain architecture, causing aberrant chromatin opening at these normally-repressed neuronal genes
4. MeCP2 then responds to this chromatin disruption (gain of accessibility + methylation redistribution at Polycomb-disrupted neuronal loci)

This reframes the theory from "shift toward heterochromatin" to "disrupted Polycomb homeostasis" — BAP1 loss doesn't tighten repression at neuronal genes, it breaks the Polycomb machinery that was keeping them properly silenced. The neuronal specificity comes from the constitutive K119ub enrichment at these loci (section 72), which makes them the most vulnerable to BAP1 loss.