# Revised Conclusions: Integrating K119ub Signal, ATAC-seq, and Methylation Data

**Date**: February 2026
**Context**: Reassessment of the BAP1 → methylation mechanistic model after Section 18 (K119ub continuous bigwig signal analysis) and Sections 12–13 (ATAC-seq accessibility correlation)

---

## 1. What Remains Solid

The core methylation phenotype is unchanged and well-supported:

- **84.6% concordance** of coordinated mC↑/hmC↓ at co-significant genes (5,708 / 6,750)
- **Gene body specificity** — 42.1% of gene bodies significant for mC, 47.3% for hmC
- **Active Promoter enrichment** — 49.9% of DMRs at Active_Promoter chromatin states
- **Synaptic/neuronal gene enrichment** — Syt1 (+17.9% mC, −15.3% hmC) as top hit

These observations are robust and novel. No existing literature describes simultaneous 5mC/5hmC profiling in any BAP1-KO model (Field et al. 2019 profiled only 5mC in uveal melanoma), and no study has reported this coordinated pattern in neural tissue.

---

## 2. What the Original Model Predicted

The initial mechanistic chain proposed:

```
BAP1 loss → H2AK119ub accumulation → chromatin compaction → TET access blocked
→ 5mC cannot be oxidized to 5hmC → coordinated mC↑ / hmC↓
```

This model — primarily a TET-blocking hypothesis — predicted:

1. Gene-specific H2AK119ub accumulation at hypermethylated loci
2. Reduced chromatin accessibility (ATAC-down) at hypermethylated genes
3. Strong correlation between K119ub gain and methylation gain per gene

---

## 3. What the Data Actually Shows

### 3.1 K119ub signal (Section 18)

**Genome-wide H2AK119ub increase is confirmed.** All DMR genes show a positive shift in K119ub signal (median log2FC = +0.012, one-sample Wilcoxon p = 9.83 × 10⁻²⁵). This validates the expected BAP1-KO phenotype — without the deubiquitinase, H2AK119ub accumulates — and is consistent with Conway et al. (2021), who showed >75% of the genome becomes more compact upon BAP1 loss.

**Gene-specific enrichment is weak.** Hypermethylated (mC Up) genes show slightly more K119ub gain than background (median log2FC +0.059 vs +0.012, Wilcoxon rank-sum p = 3.73 × 10⁻²²), but the effect size is small (Cliff's delta = +0.089).

**hmC groups do not separate from background.** Despite hmC Down being the largest group (n = 8,623 genes) and the hallmark of the coordinated pattern, it is not significantly different from all DMR genes in K119ub change (p = 0.666). hmC Up is also non-significant (p = 0.492).

**Interpretation**: H2AK119ub accumulation is genome-wide, not gene-specific. The methylation phenotype does not track local K119ub levels on a per-gene basis.

### 3.2 Chromatin accessibility (Sections 12–13)

**ATAC-down peaks are enriched at Polycomb/H3K27me3 regions** (Section 13a–b, Fisher p < 0.001). Compaction is real and preferentially affects repressive chromatin states. This is consistent with Conway et al. (2021).

**But hypermethylation and accessibility loss are poorly coupled.** Only 14.3% of hypermethylated DMRs overlap ATAC-down peaks. If chromatin compaction were the primary per-gene mechanism driving hypermethylation, much tighter coupling would be expected.

**The inverse relationship is stronger.** Hypomethylated DMRs show 38.2% overlap with ATAC-up peaks vs 2.9% ATAC-down — a strong directional association. Genes that lose methylation tend to gain accessibility, consistent with the known relationship between H2AK119ub depletion and DNA hypomethylation (Chen et al. 2024 preprint).

**Interpretation**: Compaction occurs at Polycomb regions, but it is not tightly coupled to the hypermethylation phenotype at the gene level.

---

## 4. Revised Mechanistic Model

The data supports a **multi-mechanism model** rather than a simple TET-blocking chain. Three parallel pathways, each with different support from the data:

### 4.1 DNMT3A-UDR recruitment (likely primary driver of mC gain)

BAP1 loss → genome-wide H2AK119ub increase → DNMT3A1 recruited broadly via UDR domain → de novo 5mC deposition at accessible gene bodies

**Literature basis**: 2024 cryo-EM structures (Chen et al. 2024 preprint) showed DNMT3A1's N-terminal UDR domain binds H2AK119ub and the nucleosome acidic patch through a bidentate interaction, recruiting the methyltransferase to ubiquitinated regions. This mechanism was shown to be required for postnatal neuronal development.

**What the data supports**:
- H2AK119ub is elevated genome-wide (Section 18: background log2FC > 0, p < 10⁻³⁷)
- Hypermethylation occurs preferentially at Active_Promoter regions (49.9%) — accessible chromatin where DNMT3A can act on substrate
- The weak positive K119ub-mC enrichment (Cliff's delta = +0.089) is consistent with dose-dependent DNMT3A recruitment — more ubiquitin → slightly more methyltransferase → slightly more mC
- Hypermethylation does NOT require local accessibility loss (only 14.3% of hypermethylated DMRs show ATAC-down), which fits a recruitment-based model rather than a compaction-based model

**What it doesn't explain**: DNMT3A recruitment alone does not account for the 5hmC decrease. The coordinated hmC↓ component requires some TET impediment (see 4.2).

### 4.2 Global TET impediment (explains hmC loss)

BAP1 loss → global chromatin environment shift → reduced TET enzyme residence time or processivity → 5mC → 5hmC conversion rate decreases → existing 5hmC continues through the pathway but is not replenished

**Literature basis**: TET enzymes catalyze sequential oxidation: 5mC → 5hmC → 5fC → 5caC → C (Ito et al. 2011). The first step is relatively fast; if blocked, 5mC accumulates while 5hmC depletes as existing 5hmC progresses downstream but is not replenished (Tahiliani et al. 2009). TET activity correlates with chromatin accessibility and open chromatin marks (Joshi et al. 2022). In cerebellar Purkinje neurons, TET function is essential for differentiation (Stoyanova et al. 2021), and 5hmC loss in these cells leads to neurodegeneration (Jiang et al. 2015).

**What the data supports**:
- The coordinated mC↑/hmC↓ pattern (84.6% concordance) is the signature of impaired 5mC → 5hmC conversion — expected when TET access is reduced
- Genome-wide K119ub increase and some ATAC-down at Polycomb regions are consistent with a global environment less favorable to TET activity
- The 5hmC Spearman correlation with K119ub (rho = −0.061) has the predicted sign, though the effect is near zero

**Critical nuance**: The TET impediment appears to operate as a **global environment shift**, not a gene-by-gene targeting mechanism. hmC Down genes do NOT separate from background in K119ub levels (p = 0.666). This is consistent with Conway et al. (2021) describing genome-wide (>75%) compaction rather than locus-specific effects.

**Why gene bodies are selectively affected**: Cerebellar gene bodies have the highest 5hmC of any genomic compartment in any tissue (~40% of modified cytosines; Kriaucionis & Heintz 2009; Mellén et al. 2012, 2017). Genes with the most baseline 5hmC have the most to lose when TET efficiency drops. The genomic selectivity of the methylation phenotype may reflect **baseline 5hmC vulnerability** rather than locus-specific TET blocking.

### 4.3 Polycomb-region compaction (secondary, contributes at specific loci)

BAP1 loss → H2AK119ub spreads beyond normal Polycomb targets → ATAC-down at H3K27me3+ regions → some genes lose accessibility and gain methylation

**Literature basis**: Conway et al. (2021) demonstrated that H2AK119ub spreads primarily via PCGF3/5-containing non-canonical PRC1 complexes, titrating PRC2 away from normal targets and causing widespread compaction. However, Bonnet et al. (2022) showed in *Drosophila* that excessive H2Aub can actually interfere with nucleosome stacking — the mark acts as a "rheostat" rather than a simple compaction signal.

**What the data supports**:
- ATAC-down peaks are enriched at Polycomb/H3K27me3 regions (Section 13, Fisher p < 0.001)
- Individual genes like Syt1 show both ATAC-down peaks and strong hypermethylation (+17.9% mC)
- 14.3% of hypermethylated DMRs do overlap ATAC-down — a real minority but not zero

**What limits this mechanism**: Only ~1 in 7 hypermethylated DMRs shows accessibility loss. Compaction contributes at specific loci but is not the dominant mechanism.

---

## 5. The "Double Whammy" Synthesis

The revised model is a **convergent dual mechanism**:

```
BAP1 loss
    ↓
H2AK119ub accumulation (genome-wide)
    ↓
┌───────────────────────────────────┐
│                                   │
│   GAIN pathway                    │   LOSS pathway
│   (active mC deposition)          │   (passive hmC depletion)
│                                   │
│   H2AK119ub recruits DNMT3A1     │   Global chromatin shift
│   via UDR domain                  │   reduces TET residence time
│          ↓                        │          ↓
│   De novo 5mC at accessible       │   5mC → 5hmC conversion
│   gene bodies                     │   rate decreases
│          ↓                        │          ↓
│   Hypermethylation                │   5hmC depletion
│                                   │
│   + some locus-specific           │
│     compaction (ATAC-down at      │
│     Polycomb regions) further     │
│     reinforces mC gain            │
│                                   │
└───────────────────────────────────┘
                    ↓
        Coordinated mC↑ / hmC↓
        at gene bodies (84.6%)
                    ↓
        Gene selectivity determined by:
        - Baseline 5hmC level (high in cerebellum gene bodies)
        - Chromatin state (Active_Promoter most accessible to DNMT3A)
        - Local accessibility changes (contributes at ~14% of loci)
```

This model explains several features the simple TET-block hypothesis could not:

1. **Why hypermethylation doesn't require accessibility loss**: DNMT3A is recruited by H2AK119ub, not repelled by compaction. It can act at accessible gene bodies.
2. **Why K119ub enrichment at DMR genes is weak but positive**: Genome-wide ubiquitination means DNMT3A is recruited broadly, with slight dose-dependence.
3. **Why the hmC groups don't separate in K119ub levels**: hmC depletion is a global TET efficiency problem, not a gene-specific K119ub targeting event.
4. **Why Active_Promoter states are preferentially affected**: These regions are accessible substrates for DNMT3A — the mark is there (H2AK119ub), and the DNA is reachable.
5. **Why the hypomethylated minority (2,201 genes) strongly overlaps ATAC-up**: Where H2AK119ub is redistributed away from normal targets (Conway et al. 2021), both DNMT3A recruitment and Polycomb repression are lost → accessibility gain + methylation loss.

---

## 6. Implications for the Rao Lab Paradox

López-Moyado et al. (2019) showed TET loss causes hypermethylation in euchromatin (A compartment) and paradoxical hypomethylation in heterochromatin (B compartment), explained by DNMT3A redistribution. Your data has a structural parallel:

- **Hypermethylation at Active_Promoter regions** (euchromatic, accessible) — mirrors their A-compartment hypermethylation
- **Hypomethylation at a minority of genes with ATAC-up** — could reflect DNMT3A redistribution away from normal heterochromatic targets toward newly ubiquitinated regions

The key difference: in López-Moyado et al., TET is directly lost. In your model, TET is not deleted — its access is globally reduced while DNMT3A is simultaneously recruited. The phenotypic outcome is similar (coordinated mC↑/hmC↓) but the upstream cause is distinct (chromatin state change rather than enzyme loss).

---

## 7. Connection to Loop Architecture

The chromatin loop analysis shows Polycomb loop collapse and TAD densification in BAP1-KO (grad mentor interpretation). The revised methylation model connects to this:

- **Shared upstream cause**: H2AK119ub accumulation drives both loop collapse (Polycomb domain destabilization) and methylation changes (DNMT3A recruitment + TET impediment)
- **Loop anchor accessibility**: Active_Enhancer–Active_Enhancer loops show 79.6% ATAC concordance (Section 13f), but overall concordance is ~50% — again suggesting global rather than locus-specific effects
- **Not a linear chain**: The loop changes and methylation changes are likely **parallel consequences** of the same upstream perturbation, not one causing the other

---

## 8. Revised Claims

### Can claim (supported by data)

- BAP1-KO cerebellum shows coordinated 5mC↑ / 5hmC↓ at 84.6% of co-significant genes — **novel finding**
- Gene bodies are the primary affected region, consistent with known 5hmC biology in brain (Mellén et al. 2012, 2017; Szulwach et al. 2011)
- H2AK119ub signal increases genome-wide in mutant, consistent with BAP1 deubiquitinase loss (Conway et al. 2021)
- Hypermethylation does not require locus-specific accessibility loss — only 14.3% of hypermethylated DMRs overlap ATAC-down peaks
- The data is consistent with DNMT3A-UDR recruitment to H2AK119ub-marked chromatin (Chen et al. 2024) as a contributor to hypermethylation
- This is the first simultaneous 5mC/5hmC profiling in any BAP1-KO model and the first BAP1-KO methylation study in brain tissue

### Plausible but not directly demonstrated

- DNMT3A1 is actively recruited to gene bodies via UDR–H2AK119ub interaction (structural basis exists but not shown in this tissue)
- TET enzyme access is globally reduced (no TET ChIP or activity assay performed)
- The two mechanisms (DNMT3A recruitment + TET impediment) operate in parallel — consistent with data but not individually dissected
- Baseline 5hmC levels determine gene vulnerability (would need to correlate with pre-KO 5hmC profiles)

### Cannot claim

- That locus-specific H2AK119ub accumulation targets individual genes for methylation change — the K119ub enrichment at DMR genes is small (Cliff's delta = 0.089)
- That chromatin compaction is the primary driver — ATAC-mC coupling is weak
- Causality in any direction — all associations are correlational
- That methylation changes cause the loop architecture changes or vice versa

---

## 9. Experiments That Would Distinguish Models

| Experiment                                      | What it tests                                               | Expected result if model is correct                                                                     |
| ----------------------------------------------- | ----------------------------------------------------------- | ------------------------------------------------------------------------------------------------------- |
| H2AK119ub CUT&RUN in BAP1-KO cerebellum         | Whether H2AK119ub actually accumulates in these samples     | Genome-wide increase, with no gene-specific enrichment at DMRs above background                         |
| DNMT3A ChIP-seq in BAP1-KO vs WT                | Whether DNMT3A redistributes to newly ubiquitinated regions | DNMT3A gain at Active_Promoter gene bodies that overlap DMRs                                            |
| TET1/2/3 expression (RNA-seq or western)        | Whether TET enzymes are still expressed                     | Normal expression expected — the model predicts access restriction, not expression loss                 |
| TET activity assay (5hmC ELISA on purified DNA) | Global TET conversion rate                                  | Reduced 5mC → 5hmC conversion in mutant                                                                 |
| DNMT3A-UDR mutant rescue                        | Whether UDR-mediated recruitment is required                | UDR-dead DNMT3A fails to rescue hypermethylation at gene bodies                                         |
| EZH2 inhibitor treatment                        | Whether PRC2-mediated compaction contributes                | Partial rescue of accessibility but incomplete rescue of methylation (if DNMT3A-UDR is the main driver) |

---

## 10. Summary

The original TET-blocking hypothesis was a reasonable starting point, but the K119ub bigwig and ATAC-seq data reveal a more complex picture. The methylation phenotype is not driven by gene-specific H2AK119ub targeting or locus-specific compaction. Instead, the data points to **convergent global mechanisms**: genome-wide H2AK119ub accumulation recruits DNMT3A to accessible gene bodies (active mC gain), while a global chromatin environment shift reduces TET efficiency (passive hmC loss). Locus-specific compaction at Polycomb regions contributes at a minority (~14%) of hypermethylated loci. The selectivity of which genes are affected likely reflects their baseline 5hmC levels and chromatin accessibility rather than their local H2AK119ub dynamics.

This revised interpretation is arguably more interesting than the simple model: it identifies **parallel pathways converging on the same phenotype**, explains the Active_Promoter enrichment pattern, and connects to the broader DNMT3A-UDR structural biology that is an active area of discovery (Chen et al. 2024; Zou et al. 2024).

---

## References

1. Chen et al. (2024). Cancer-associated DNA hypermethylation of Polycomb targets requires DNMT3A dual recognition of histone H2AK119 ubiquitination and the nucleosome acidic patch. *bioRxiv*. doi:10.1101/2024.03.18.585588
2. Conway et al. (2021). BAP1 enhances Polycomb repression by counteracting widespread H2AK119ub1 deposition and chromatin condensation. *Molecular Cell*. PMID: 34186021
3. Bonnet et al. (2022). PR-DUB preserves Polycomb repression by preventing excessive H2AK118 mono-ubiquitylation. *Genes & Development*. PMID: 35961776
4. Field et al. (2019). BAP1 loss is associated with DNA methylomic repatterning in highly aggressive Class 2 uveal melanomas. *Clinical Cancer Research*. PMC6744995
5. Ito et al. (2011). Tet proteins can convert 5-methylcytosine to 5-formylcytosine and 5-carboxylcytosine. *Nature*. PMID: 21778364
6. Jiang et al. (2015). Alteration in 5-hydroxymethylcytosine-mediated epigenetic regulation leads to Purkinje cell vulnerability in ATM deficiency. *Brain*. doi:10.1093/brain/awv284
7. Joshi et al. (2022). Mechanisms that regulate the activities of TET proteins. *Cellular and Molecular Life Sciences*. doi:10.1007/s00018-022-04396-x
8. Kriaucionis & Heintz (2009). The nuclear DNA base 5-hydroxymethylcytosine is present in Purkinje neurons and the brain. *Science*. PMID: 19372393
9. López-Moyado et al. (2019). Paradoxical association of TET loss of function with genome-wide DNA hypomethylation. *PNAS*. PMID: 31371502
10. Mellén et al. (2012). MeCP2 binds to 5hmC enriched within active genes and accessible chromatin in the nervous system. *Cell*. PMID: 23260135
11. Mellén et al. (2017). 5-hydroxymethylcytosine accumulation in postmitotic neurons results in functional demethylation of expressed genes. *PNAS*. doi:10.1073/pnas.1708044114
12. Stoyanova et al. (2021). 5-Hydroxymethylcytosine-mediated active demethylation is required for mammalian neuronal differentiation and function. *eLife*. doi:10.7554/eLife.66973
13. Szulwach et al. (2011). 5-hmC-mediated epigenetic dynamics during postnatal neurodevelopment and aging. *Nature Neuroscience*. PMID: 21804537
14. Tahiliani et al. (2009). Conversion of 5-methylcytosine to 5-hydroxymethylcytosine in mammalian DNA by MLL partner TET1. *Science*. PMID: 19372391
15. Thomas et al. (2023). Structural basis of histone H2A lysine 119 deubiquitination by Polycomb repressive deubiquitinase BAP1/ASXL1. *Science Advances*. doi:10.1126/sciadv.adg9832
16. Zou et al. (2024). RNA m5C oxidation by TET2 regulates chromatin state and leukaemogenesis. *Nature*. doi:10.1038/s41586-024-07969-x

---

As the most well-studied epigenetic modification in cancer, DNA methylation is the process of 5-methylcytosine (5-mC) formation, which is dynamically regulated by DNA methyltransferases (DNMTs) and ten-eleven translocation (TET) family members that are responsible for the transfer of methyl groups to cytosine residues and DNA demethylation via the oxidative catalysis of 5-mC to form 5-hmC (58). The downregulation of 5-hmC induced by the reduction of isocitrate dehydrogenase 2 (IDH2) and TET family members is a valuable epigenetic marker with relatively high diagnostic and prognostic significance in melanoma (59). Quite a few genes implicated in cell differentiation, epithelial–mesenchymal transformation, PI3K/mTOR signaling, and cell metabolism are differentially methylated in the promoter regions between melanomas and benign nevus (60, 61). Genome-wide DNA methylation analysis also unveils the epigenetic heterogeneity of melanoma metastases in specific organs (62). Moreover, several studies have demonstrated the engagement of DNA methylation in the regulation of specific oncogene or tumor suppressors in melanoma. For example, Field et al. showed that canonical tumor suppressor BAP1 is hypermethylated in uveal melanoma. In addition, integrative analysis reveals that SOX9 is methylated and lowly expressed in the highly proliferative group, and SOX9 overexpression results in decreased proliferation but increased invasion of melanoma cells (15). Therefore, DNA methylation plays a versatile role in mediating various aspects of melanoma cell biology.



https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.868786/full