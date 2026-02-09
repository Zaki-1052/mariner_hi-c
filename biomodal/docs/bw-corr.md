# BigWig Correlation Analysis: MeCP2 and H2AK119ub at DMRs

## Purpose

Examine MeCP2 occupancy and H2AK119ub levels at differentially methylated regions to test whether the observed methylation changes correlate with predicted chromatin changes from BAP1 loss.

## Background Mechanisms

### MeCP2 Binding Affinity (Mellén et al. 2017)

MeCP2 binding affinity depends on the specific modified cytosine context:

| Modification | Binding Affinity | Functional Effect |
|--------------|------------------|-------------------|
| **5mCG** | High | Repressive (MeCP2 recruited) |
| **5hmCG** | Low (similar to unmodified) | Permissive (MeCP2 dissociates) |
| **5mCA** | High | Repressive |
| **5hmCA** | High | Repressive (no change from 5mCA) |

**Key insight:** Converting 5mCG → 5hmCG reduces MeCP2 binding = "functional demethylation"

**In control cerebellum:**
- Gene bodies accumulate 5hmCG through TET-mediated oxidation
- High 5hmCG / low 5mCG → Low MeCP2 occupancy → Transcription facilitated

**Predicted in BAP1-KO (if TET blocked):**
- 5mC accumulates, 5hmC depletes at gene bodies
- High 5mCG / low 5hmCG → **Increased MeCP2 occupancy** → More repression

### H2AK119ub Accumulation (Conway et al. 2021)

BAP1 (via PR-DUB complex) removes H2AK119ub from chromatin. When BAP1 is lost:

1. **H2AK119ub accumulates** and spreads beyond normal Polycomb targets
2. **>75% of genome** becomes more compact
3. **Chromatin compaction** may restrict TET enzyme access to DNA
4. PRC2 is titrated away from normal targets to intergenic regions
5. Transcription initiation is globally reduced

**Predicted at DMRs:**
- Regions with increased 5mC / decreased 5hmC may show elevated H2AK119ub
- Elevated H2AK119ub may correlate with increased MeCP2 occupancy
- Combined effect: chromatin compaction + increased MeCP2 binding = enhanced repression

### DNMT3A Recruitment (Chen et al. 2024)

DNMT3A1 contains a UDR domain that directly binds H2AK119ub, recruiting DNA methyltransferase to ubiquitinated regions.

**Potential mechanism:**
- H2AK119ub accumulation → DNMT3A recruitment → 5mC deposition
- Combined with reduced TET access → enhanced hypermethylation

## Analysis Goals

### 1. MeCP2 Correlation at DMRs
- Compare MeCP2 ChIP-seq signal at hypermethylated vs hypomethylated DMRs
- Test prediction: Hypermethylated DMRs (5mC↑/5hmC↓) show **increased MeCP2**

### 2. H2AK119ub Correlation at DMRs
- Compare H2AK119ub ChIP-seq signal at hypermethylated DMRs
- Test prediction: Hypermethylated DMRs show **increased H2AK119ub** in mutant

### 3. MeCP2 vs H2AK119ub Co-occurrence
- Examine overlap between H2AK119ub-enriched and MeCP2-enriched DMRs
- Test whether chromatin compaction (H2AK119ub) and functional repression (MeCP2) co-localize

## Expected Patterns (Hypotheses)

If the TET-block model is correct:

| DMR Type | 5mC Change | 5hmC Change | Expected MeCP2 | Expected H2AK119ub |
|----------|------------|-------------|----------------|---------------------|
| Hypermethylated | ↑ | ↓ | Increased (mutant > control) | Increased (mutant > control) |
| Hypomethylated | ↓ | ↑ | Decreased (mutant < control) | Variable |
| Not significant | No change | No change | No change | Variable |

## Data Requirements

### ChIP-seq BigWig Files Needed
1. **MeCP2 ChIP-seq** (control and mutant cerebellum)
2. **H2AK119ub ChIP-seq** (control and mutant cerebellum)

### DMR Regions (Already Available)
- Gene body DMRs from `modality/outputs/run-2/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_*/`
- Significant DMRs (q < 0.05) stratified by direction (hyper vs hypo)

## Analysis Steps

1. **Load DMR coordinates** (BED format with direction annotation)
2. **Compute mean ChIP signal** over DMR windows (±500bp from DMR center)
3. **Compare signal distributions:**
   - Hypermethylated DMRs vs control regions
   - Hypomethylated DMRs vs control regions
   - Mutant vs control at same DMR coordinates
4. **Visualize:**
   - Heatmaps of ChIP signal at DMRs
   - Box plots comparing signal by DMR type
   - Scatterplots: Methylation change (x) vs ChIP signal change (y)

## Interpretation Guidelines

**Strong support for model:**
- Hypermethylated DMRs show increased MeCP2 AND increased H2AK119ub in mutant

**Partial support:**
- Either MeCP2 OR H2AK119ub correlates, but not both (suggests one mechanism dominates)

**Against model:**
- No correlation between methylation changes and chromatin marks
- Would suggest methylation changes are secondary or indirect

## References

- Mellén et al. (2017) PNAS — MeCP2 binding affinity differences (5hmCG = low affinity)
- Conway et al. (2021) Molecular Cell — BAP1 loss causes H2AK119ub spread and compaction
- Chen et al. (2024) bioRxiv — DNMT3A recruited by H2AK119ub via UDR domain
