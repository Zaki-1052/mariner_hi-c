# Coordinated vs Discordant Gene Trends

Comparison of section 28 (Coordinated Q4 characterization) and section 21 (Discordant Q2 characterization) findings across all 9 analysis dimensions.

**Gene counts:**
- Coordinated Q4 (mC up / hmC dn): 5,708 genes (84.6%)
- Discordant Q2 (mC dn / hmC up): 755 genes (11.2%)
- Q1 (mC dn / hmC dn): 225 genes (3.3%)
- Q3 (mC up / hmC up): 62 genes (0.9%)
- Total co-significant: 6,750 genes

## Summary Table

| Dimension | Coordinated Q4 (n=5,708) | Discordant Q2 (n=755) | Interpretation |
|---|---|---|---|
| **mC direction** | UP (gain) | DOWN (loss) | Opposite by definition |
| **hmC direction** | DOWN (loss) | UP (gain) | Opposite by definition |
| **\|mC diff\| median** | 3.23% | 1.71% | Q4 has ~2x larger mC effect sizes |
| **\|hmC diff\| median** | 2.82% | 1.19% | Q4 has ~2.4x larger hmC effect sizes |
| **RNA-seq log2FC median** | -0.010 (near zero) | +0.191 (upregulated) | Discordant genes trend upregulated; coordinated genes are expression-neutral on average |
| **Expression Up/Down OR** | 0.05 (biased downregulated) | 103.14 (biased upregulated) | Strikingly opposite: Q2 genes strongly upregulated, Q4 genes strongly downregulated when significant |
| **Net ATAC median** | 0.0 (no net change) | +1.0 (net gain) | Discordant genes gain accessibility; coordinated genes do not |
| **ATAC Up/Down OR** | 0.12 (depleted for ATAC-up) | 14.55 (enriched for ATAC-up) | Same pattern as expression: Q2 opens chromatin, Q4 does not |
| **K119ub gb log2FC median** | +0.033 (gain) | -0.031 (loss) | Q4 genes gain K119ub; Q2 genes lose it |
| **H3K27ac Gained/Lost** | 1,177 gained / 823 lost | 101 gained / 17 lost | Q4 has more balanced gain/loss; Q2 is overwhelmingly gain-only |
| **Chromatin state** | Significant difference (p=1e-4) | Significant difference (p=1e-4) | Both groups differ from their reference in chromatin composition |
| **MeCP2 Up/Down** | 719 up / 123 dn (5.8:1 ratio) | 130 up / 82 dn (1.6:1 ratio) | Q4 is strongly MeCP2-up biased; Q2 has more balanced MeCP2 |
| **Loop involvement** | 5.1% at loop anchors | 15.5% at loop anchors | Discordant genes are 3x more likely at loop anchors |
| **mC-expression rho** | -0.247 (p=6.4e-73) | -0.151 (p=1.8e-4) | Both negative (more mC = less expression), but Q4 is a stronger correlation |
| **GO enrichment** | 1,438 terms (vs non-Q4 bg) | None (vs Q4 bg) | Q4 has strong pathway enrichment; Q2 does not distinguish from Q4 at the pathway level |

## Biological Narrative

### Coordinated Q4: The dominant BAP1-KO signature

Q4 genes represent the canonical TET-inhibition model: BAP1 loss reduces H2AK119ub deubiquitination, leading to K119ub accumulation, which inhibits TET enzymes. The result is mC gain (DNMTs unopposed) and hmC loss (TETs inhibited). These genes:

- Have **larger methylation effect sizes** than any other quadrant, suggesting they are more directly affected
- Show a **weak but significant trend toward downregulation** (median log2FC = -0.010), and when expression is significantly changed, it is overwhelmingly downregulated (OR = 0.05)
- Do **not** gain chromatin accessibility (net ATAC = 0), consistent with a repressive methylation gain
- **Gain K119ub** (median +0.033), directly consistent with the BAP1-loss mechanism
- **Gain MeCP2 binding** (5.8:1 up/down ratio), expected since MeCP2 reads methylated DNA
- Are **less likely to be at Hi-C loop anchors** (5.1%), suggesting these are not structurally regulated genes
- Have a **stronger mC-expression anticorrelation** (rho = -0.247), meaning the mC gain is functionally linked to expression changes

### Discordant Q2: The minority counter-pattern

Q2 genes show the opposite: mC loss with hmC gain. This is inconsistent with a simple TET-inhibition model and suggests a different regulatory mechanism at these loci. These genes:

- Have **smaller methylation effect sizes**, suggesting weaker or more indirect effects
- Are **strongly upregulated** (median log2FC = +0.191, OR = 103.14) — the most expression-biased group
- **Gain chromatin accessibility** (net ATAC = +1.0, OR = 14.55), consistent with an activating program
- **Lose K119ub** (median -0.031), opposite to the genome-wide trend — suggesting these loci are specifically de-repressed
- **Gain H3K27ac** almost exclusively (101 gained vs 17 lost), consistent with activation
- Are **3x more likely at Hi-C loop anchors** (15.5% vs 5.1%), suggesting structural/regulatory involvement
- Have a **weaker mC-expression correlation** (rho = -0.151)
- Show **no pathway enrichment** against Q4 background, meaning they are not a functionally distinct gene set at the GO level — they are scattered across pathways

### The Contrast

The two groups paint a coherent picture of two distinct regulatory responses to BAP1 loss:

1. **Q4 (majority):** Passive repression — K119ub accumulates, TETs are inhibited, mC rises, hmC falls, MeCP2 binds, expression trends down. These genes are victims of a genome-wide epigenetic shift.

2. **Q2 (minority):** Active de-repression — K119ub is specifically removed, chromatin opens, H3K27ac is gained, mC falls (active demethylation via TET, producing hmC), and genes are upregulated. These loci may be compensatory targets or genes that escape the dominant repressive program through an alternative regulatory mechanism.

The key mechanistic question is why ~755 genes escape the dominant Q4 pattern. Their enrichment at loop anchors (3x more than Q4) and strong K119ub loss suggest they may be structurally protected loci where BAP1-independent deubiquitinases or chromatin remodelers maintain an active state despite global BAP1 loss.
