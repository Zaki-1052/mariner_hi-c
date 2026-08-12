# Sections 74-78: Responses to Jai's Follow-Up Questions

Results from sections 74-78, addressing each point Jai raised on 2026-06-23.

Gene set definitions used throughout:
- **Neuronal (broad)**: 5,614 genes from org.Mm.eg.db GO BP terms matching `synap|neuron|axon|dendrit|nervous` (section 72, unbiased by methylation data)
- **Neuronal (narrow)**: 1,149 genes from GO enrichment of significant DMRs matching `synap|neuron|axon` (section 61, biased by selection)
- **Synapse/axon**: 2,321 genes from GO BP terms matching `synap|axon` only (section 76)
- **Coordinated**: 6,069 genes with mC significantly up AND hmC significantly down (both q<0.05)
- **MeCP2-Up**: 79 genes with mecp2_nearest_fdr < 0.05 and mecp2_mean_fold > 0

---

## Point 1: Gene Set Overlap + Neuronal Methylation Levels

> Jai: "could we create a venn diagram or counts table for coordinated, neuronal, mecp2 up, overlap? ... make a plot for neuronal genes showing total methylation ctrl vs mut, and 5mc and 5hmc levels"

### Overlap (Section 74a)

Universe: 23,150 genes (quadrant master).

| Partition | Count | % |
|-----------|------:|----:|
| Neuronal only | 3,888 | 16.8% |
| Coordinated only | 4,592 | 19.8% |
| MeCP2-Up only | 17 | 0.07% |
| Neuronal + Coordinated | 1,426 | 6.2% |
| Neuronal + MeCP2-Up | 11 | 0.05% |
| Coordinated + MeCP2-Up | 35 | 0.15% |
| **All three** | **16** | **0.07%** |
| None | 13,165 | 56.9% |

Pairwise Fisher enrichment (BH-corrected):

| Pair | Overlap | OR | adj p |
|------|--------:|---:|------:|
| Neuronal x Coordinated | 1,442 | 1.05 | 0.14 (NS) |
| Neuronal x MeCP2-Up | 27 | 1.73 | 0.034* |
| **Coordinated x MeCP2-Up** | **51** | **5.16** | **2.82e-12*** |

The strongest association is Coordinated x MeCP2-Up (OR=5.16): MeCP2 redistribution tracks methylation redistribution far more than neuronal identity. Neuronal genes are NOT preferentially coordinated (OR=1.05, NS).

### Methylation at Neuronal Genes (Sections 74b-d)

n = 4,110 neuronal genes with valid methylation data.

| Measure | Ctrl median | Mut median | Delta median | p |
|---------|------------|------------|-------------|---|
| Total (mC+hmC) | 0.8235 | 0.8225 | **-0.003** | < 2.2e-16 |
| 5mC | 0.6782 | 0.6858 | **+0.004** | < 2.2e-16 |
| 5hmC | 0.1272 | 0.1168 | **-0.009** | < 2.2e-16 |

The mC up / hmC down pattern is present at neuronal genes, but total methylation **decreases** slightly because the hmC loss (-0.009) exceeds the mC gain (+0.004). This contradicts the section 61a finding that total methylation increases at neuronal genes — see Point 6 below for explanation.

### Key plots
- `74a_geneset_venn` — 3-way Venn with Fisher ORs in caption
- `74b_neuronal_total_methylation` — total methylation violins
- `74c_neuronal_5mc_levels` — 5mC violins
- `74d_neuronal_5hmc_levels` — 5hmC violins
- `74_composite` — all panels combined

---

## Point 2: MeCP2 Signal Direction Paradox

> Jai: "how does mecp2 signal drop in the mutant? bc in the diffbind volcano plot, there were like 8000 loci for increased mecp2 binding in the mutant"

### Resolution (Section 75)

The 7,686 DiffBind UP peaks concentrate at only **2,052 genes** (~10% of genes with MeCP2 peaks). The remaining 90% of genes have no significant MeCP2 peaks.

| Gene category | n genes | Median gene-body fold |
|--------------|--------:|---------------------:|
| Has UP peak(s) | 2,052 | **+0.224** |
| Has DOWN peak(s) only | 518 | -0.128 |
| No significant peaks | 18,625 | -0.021 |

Genes with UP peaks genuinely gain MeCP2 at the gene-body level (median +0.224). The genome-wide median drops because 88% of genes have no sig peaks and show slight negative fold. MeCP2 **redistributes**: it concentrates at ~2,000 genes while thinning everywhere else.

### Why UP peaks don't register as gene-body signal

MeCP2 UP peaks are disproportionately **distal intergenic** — far from gene bodies:

| Annotation | UP peaks | DOWN peaks |
|-----------|--------:|----------:|
| Distal Intergenic | 51.7% | 19.5% |
| Intron | 41.0% | 61.8% |
| Promoter | 2.2% | 8.0% |
| Exon/UTR | 5.0% | 10.4% |

Chi-squared X^2 = 503.5, p = 1e-4. UP peaks are being recruited to new distal/intergenic sites; DOWN peaks are more genic.

### K119ub vs MeCP2 neuronal specificity

> Jai: "K119ub controls a broad range of regular + neuronal loci, but mecp2 was selective to synapse/axon genes right?"

Correct. GSEA GO BP terms at q < 0.05:
- **K119ub**: 115 significant terms (3 matching neuronal pattern = 2.6%)
- **MeCP2**: 1 significant term (**synapse assembly**, NES=1.68, q=0.025 = 100% neuronal)

K119ub broadly targets developmental/Polycomb genes. MeCP2's only GSEA hit is neuronal.

### Key plots
- `75a_peak_distribution_per_gene` — bar chart of gene categories by MeCP2 peak status
- `75b_genebody_signal_by_peak_status` — violins showing UP-peak genes DO have positive fold
- `75c_gsea_term_comparison` — stacked bar comparing K119ub (115) vs MeCP2 (1) sig terms
- `75d_peak_annotation_distribution` — annotation breakdown for UP vs DOWN peaks
- `75_composite` — all panels combined

---

## Point 3: P-Values for 72g Violins + Triple-Overlap Chromatin

> Jai: "do we have p values for the neuronal vs non-neuronal violin plots? ... also make this same thing for the overlap set of mecp2/ub/coordinated genes"

### Chromatin mark p-values (Section 76a)

| Mark | Neuronal median | Non-neuronal median | p (Wilcoxon) |
|------|:-:|:-:|---:|
| ATAC | +0.107 | +0.087 | 6.95e-15 |
| K27ac | +0.046 | +0.041 | 0.027 |
| K27me3 | -0.103 | -0.061 | 5.53e-4 |

All three marks show significant differences. Neuronal genes gain MORE accessibility (ATAC), gain slightly more K27ac, and lose MORE K27me3 than non-neuronal genes. This is chromatin **opening**, not closing — consistent with Polycomb disorganization, not heterochromatin reinforcement.

### Triple-overlap genes (Section 76b)

The 16 triple-overlap genes (Neuronal AND Coordinated AND MeCP2-Up): Ap3b1, Astn2, Cntn6, Epyc, Fgf1, Fut9, Gprin3, Hcn1, Hif1a, Il1rapl1, Lgi1, Micu3, Ntn4, Prom1, Snca, Tspan7.

n = 16 is too small for robust violin comparisons, but Kruskal-Wallis across 4 groups (triple/neuronal-only/coordinated-only/rest) is significant for all three marks (all p < 2.2e-16).

### Are triple-overlap genes in top K119ub deciles? (Section 76d)

No. 0/16 triple-overlap genes are in the top K119ub decile (D10). They're scattered across D1-D9. The triple-overlap genes are not K119ub-extreme.

However, neuronal and synapse/axon genes **are** enriched in D10:

| Gene set | In D10 | Total | OR | p |
|----------|-------:|------:|---:|--:|
| Triple overlap | 0 | 16 | 0.00 | 0.40 (NS) |
| Synapse/axon | 306 | 2,069 | 1.68 | 2.1e-13 |
| Neuronal | 677 | 4,030 | 2.34 | < 2.2e-16 |

### Key plots
- `76a_neuronal_chromatin_with_pvalues` — violins with annotated p-values and n
- `76b_four_group_chromatin` — triple/neuronal/coordinated/rest comparison
- `76d_k119ub_decile_enrichment` — forest plot of D10 ORs

---

## Point 4: Synapse/Axon Genes vs Broader Neuronal

> Jai: "maybe there is something special about the synapse/axon genes rather than just neuronal genes overall"

### Synapse/axon chromatin comparison (Section 76c)

| Mark | Synapse vs Broader neuronal | adj p |
|------|:--:|------:|
| ATAC | +0.003 difference | 0.74 (NS) |
| K27ac | -0.004 difference | 0.76 (NS) |
| **K27me3** | **-0.044 difference** | **0.007** |

Synapse/axon genes are NOT different from broader neuronal for ATAC or K27ac. But they lose **significantly more K27me3** (p = 0.003, adj p = 0.007). The broader neuronal vs non-neuronal K27me3 difference is actually NS (p = 0.29), meaning the K27me3 loss specificity resides almost entirely in the synapse/axon subset.

This partially confirms Jai's intuition: synapse/axon genes ARE special, but specifically for **Polycomb de-repression** (K27me3 loss), not for accessibility or enhancer activation.

### Key plots
- `76c_synapse_vs_neuronal_chromatin` — 3-class violins (synapse / broader neuronal / non-neuronal)
- `76d_k119ub_decile_enrichment` — both sets enriched in D10 but neuronal more so

---

## Point 5: Young vs Adult MeCP2 Developmental Trajectory

> Jai: "MeCP2 binding increases as neurons mature, but super increases in mutants. Need to 1) account for age-related increases in both, 2) find loci uniquely increased in mut aging"

### Overview (Section 77a)

After negating Fold (raw files had young/adult, we converted to adult/young):

| Genotype | Aging-UP | Aging-DOWN | NS | Total peaks |
|----------|--------:|-----------:|---:|----------:|
| Control | 10,930 | 2,822 | 383,418 | 397,170 |
| Mutant | 23,117 | 10,646 | 361,952 | 395,715 |

Mutants have **2.1x more aging-UP peaks** and **3.8x more aging-DOWN peaks** than controls.

### Overlap analysis (Section 77b)

| Category | Peaks | Genes |
|----------|------:|------:|
| Ctrl aging-UP | 10,930 | 2,908 |
| Mut aging-UP | 23,117 | 4,274 |
| Shared (ctrl peaks overlapping mut) | 7,305 (66.8% of ctrl) | 2,620 |
| Ctrl-unique | 3,625 | 288 |
| **Mut-unique** | **15,812** | **1,654** |

66.8% of ctrl aging-UP peaks also appear in mut, meaning most normal aging changes are captured in both genotypes. But the mutant has **15,812 additional unique peaks** at **1,654 unique genes** that don't appear in ctrl aging.

### Mut-unique aging genes: enriched for neuronal (Section 77c)

- 404 / 1,654 mut-unique genes are neuronal (24.4% — slightly above the 19.7% genome-wide rate)
- GO enrichment: **435 significant BP terms**, of which **49 are neuronal** (11.3%)
- The mut-specific aging genes have strong neuronal enrichment

### Fold comparison at shared loci (Section 77d)

At the 7,305 peaks shared between ctrl and mut aging:
- Ctrl aging fold median: **1.829**
- Mut aging fold median: **2.241**
- Paired Wilcoxon: **p < 2.2e-16**

At loci where BOTH genotypes gain MeCP2 with age, the mutant gains significantly MORE (22% higher fold on average). This confirms Jai's "super-increase" observation quantitatively.

### Key plots
- `77a_aging_overview` — bar chart of UP/DOWN/NS by genotype
- `77b_aging_overlap_venn` — gene-level Venn (2,620 shared / 288 ctrl-only / 1,654 mut-only)
- `77c_mut_specific_go_enrichment` — GO dotplot for the 1,654 mut-unique genes
- `77d_shared_peak_fold_comparison` — scatter of ctrl vs mut aging fold (above diagonal = mut ages more)
- `77_composite` — all panels combined

---

## Point 6: Neuronal Gene Set Bias (Section 78)

The section 61 neuronal gene set (1,149 genes from DMR enrichment) was derived from GO enrichment of genes with **significant methylation changes** — circular reasoning when asking whether neuronal genes have different methylation. Section 78 redoes the key section 61 analyses with the unbiased broad set.

### Total methylation: direction flips (Section 78a, 78e)

| Gene set | n | Mean delta total | Direction | p |
|----------|--:|:--:|:-:|--:|
| All genes | 20,915 | -0.00139 | DOWN | < 2.2e-16 |
| Non-neuronal | 16,797 | -0.00119 | DOWN | < 2.2e-16 |
| **Neuronal (broad)** | **4,118** | **-0.00222** | **DOWN** | **< 2.2e-16** |
| Synapse/axon | 2,099 | -0.00128 | DOWN | 4.23e-9 |
| Coordinated | 6,069 | +0.00697 | UP | < 2.2e-16 |
| Neur + Coord | 1,442 | +0.00757 | UP | < 2.2e-16 |
| MeCP2-Up | 79 | +0.03132 | UP | 2.23e-11 |

The direction flips:
- **Narrow neuronal (1,149)**: mean delta = **+0.012 (INCREASES)**
- **Broad neuronal (4,118)**: mean delta = **-0.002 (DECREASES)**
- Broad-only (in broad but not narrow, 2,969 genes): mean delta = **-0.008**

The narrow set's total methylation increase was entirely driven by DMR selection bias. Total methylation increases only at **coordinated** genes (those with significant mC+hmC changes) and **MeCP2-Up** genes — not at neuronal genes per se.

### Stoichiometry slopes: neuronal genes ARE stoichiometric (Section 78c)

| Gene group | Slope | 95% CI | Consistent with -1? |
|-----------|------:|--------|:---:|
| All genes | -0.959 | [-0.978, -0.940] | No |
| Non-neuronal | -0.949 | [-0.970, -0.928] | No |
| **Neuronal (broad)** | **-0.996** | **[-1.039, -0.952]** | **Yes** |
| **Synapse/axon** | **-1.020** | **[-1.078, -0.962]** | **Yes** |
| Coordinated | -1.288 | [-1.323, -1.252] | No (steeper) |
| Neur + Coord | -1.303 | [-1.376, -1.229] | No (steeper) |

This is a key finding: neuronal genes have a stoichiometry slope of -0.996, **consistent with -1** (dehydroxymethylase conversion). The genome-wide deviation from stoichiometry (-0.959) is driven by **non-neuronal genes** (-0.949). Coordinated genes overshoot at -1.29, meaning they gain more mC than they lose hmC (net methylation gain, consistent with their +0.007 delta total).

### TET-KO phenocopy: weaker at neuronal genes (Section 78d)

| Gene group | Spearman rho | p |
|-----------|:---:|--:|
| All genes | 0.220 | < 2.2e-16 |
| Non-neuronal | 0.246 | < 2.2e-16 |
| Neuronal (broad) | 0.137 | < 2.2e-16 |
| Synapse/axon | 0.121 | 2.98e-8 |
| Coordinated | 0.176 | < 2.2e-16 |

The BAP1-KO vs TET-KO correlation is significant everywhere but **weakest at neuronal/synapse genes**. Neuronal genes resemble TET-KO less than non-neuronal genes do, consistent with the stoichiometric (dehydroxymethylase-like) conversion at neuronal loci rather than TET inhibition.

### Quadrant characterization with broad set (Section 78g)

The 72 MeCP2-Up + K119ub-Up quadrant genes:
- **Broad neuronal**: 25/72 (35%) — up from 15/72 (21%) with narrow set
- **Synapse/axon**: 17/72 (24%)

The broad neuronal set captures 10 additional genes in this quadrant that the narrow set missed, increasing the neuronal fraction from 21% to 35%.

### Key plots
- `78a_total_methylation_forest` — forest plot of delta total by gene group
- `78b_stoichiometry_scatter` — delta-5mC vs delta-5hmC colored by neuronal class
- `78c_stoichiometry_slope_forest` — slopes with CIs vs reference -1
- `78d_bap1_vs_tetko` — BAP1-KO vs TET-KO scatter
- `78e_narrow_vs_broad_bias` — side-by-side violins showing the direction flip
- `78g_quadrant_characterization` — quadrant genes with broad neuronal classification
- `78_composite` — all panels combined

---

## Open Questions and Future Directions

1. **Why are neuronal genes stoichiometric while non-neuronal deviate?** The slope of -1.0 at neuronal genes suggests direct 5hmC-to-5mC conversion (DNMT3A dehydroxymethylase) rather than independent TET inhibition + de novo methylation. Is there a structural reason (chromatin context, Polycomb state) that channels the DNMT3A reaction differently at neuronal loci?

2. **The 1,654 mut-unique aging genes**: These are genes where MeCP2 increases with age ONLY in mutants. 24.4% are neuronal, with 435 significant GO terms (49 neuronal). What drives the additional age-related MeCP2 accumulation at these loci specifically in the BAP1-KO context? Is it the excess K119ub creating aberrant chromatin states that attract MeCP2 during neuronal maturation?

3. **Synapse/axon K27me3 specificity**: Synapse/axon genes lose more K27me3 than broader neuronal genes (p=0.007) but don't differ on ATAC or K27ac. This selective Polycomb de-repression without corresponding accessibility gain is unusual — are these genes transitioning from Polycomb-repressed to a different silencing mechanism? Or is the K27me3 loss accompanied by K119ub gain that maintains repression?

4. **MeCP2 intergenic redistribution**: 52% of MeCP2 UP peaks are distal intergenic. What are these sites? Enhancers? Repetitive elements? The genic DOWN peaks (62% intronic, 8% promoter) suggest MeCP2 is vacating gene bodies and concentrating at distal regulatory elements.

5. **The coordinated gene set as the true MeCP2 substrate**: The OR=5.16 for Coordinated x MeCP2-Up (vs OR=1.73 for Neuronal x MeCP2-Up) suggests that methylation redistribution, not neuronal identity, is what MeCP2 reads. Are the 51 Coordinated+MeCP2-Up genes (minus the 16 triple-overlap) a distinct biological group from the 16 triple-overlap genes?

6. **Aging fold magnitude**: At shared aging-UP loci, the mutant fold is 22% higher than ctrl (2.24 vs 1.83). Is this excess proportional to the K119ub gain at those loci? A correlation of (mut_aging_fold - ctrl_aging_fold) vs K119ub_fold would test whether K119ub accumulation amplifies age-related MeCP2 gain.

---

## Limitations

- **n = 4 samples** (2 ctrl, 2 mut, 1 per sex per condition): sex and genotype effects are partially confounded despite the sex covariate in run-5. All methylation results should be interpreted with this caveat.
- **Broad neuronal set (5,614 genes, 23% of genome)** is large and heterogeneous. The GO pattern `synap|neuron|axon|dendrit|nervous` captures everything from ion channels to transcription factors. The synapse/axon subset (2,321) is more specific but still broad.
- **Triple overlap (n=16)** is too small for robust statistical comparisons. Individual gene-level interpretations should be cautious.
- **ChIPseeker annotation for section 77**: 137-276 peaks per file (~0.03%) couldn't be annotated (fell outside known mm10 chromosomes). These were matched back by coordinates; ~0.07% of peaks remain unannotated.
- **Young vs adult DiffBind Fold convention**: The raw files had young/adult (Fold < 0 = adult-higher). We negated to adult/young for biological interpretability. Verified against mecp2.png expected sig counts (13,752 ctrl, 33,763 mut).
- **Section 78 stoichiometry slopes at small n**: The MeCP2-Up group (n=79) has wide CIs [-1.31, -0.54], technically consistent with -1 but underpowered to distinguish mechanism.

---

## Plot Index

All plots saved in 4 formats (PDF/PNG/SVG/JPG) under `plots/visualizations/`.

### Section 74 — Gene Set Overlap and Neuronal Methylation
| Plot | File | Best for |
|------|------|----------|
| 3-way Venn | `74_geneset_overlap_methylation/74a_geneset_venn/` | Showing overlap counts + Fisher ORs |
| Total methylation violins | `74_geneset_overlap_methylation/74b_neuronal_total_methylation/` | Total meth DECREASES at neuronal |
| 5mC violins | `74_geneset_overlap_methylation/74c_neuronal_5mc_levels/` | 5mC increases |
| 5hmC violins | `74_geneset_overlap_methylation/74d_neuronal_5hmc_levels/` | 5hmC decreases (more than mC gain) |
| Composite | `74_geneset_overlap_methylation/74_composite/` | Full panel |

### Section 75 — MeCP2 Signal Reconciliation
| Plot | File | Best for |
|------|------|----------|
| Peak-per-gene bar | `75_mecp2_signal_reconciliation/75a_peak_distribution_per_gene/` | 7,686 UP peaks at only 2,052 genes |
| Gene-body signal by peak status | `75_mecp2_signal_reconciliation/75b_genebody_signal_by_peak_status/` | UP-peak genes DO have positive fold |
| GSEA term comparison | `75_mecp2_signal_reconciliation/75c_gsea_term_comparison/` | K119ub 115 terms vs MeCP2 1 term |
| Peak annotation | `75_mecp2_signal_reconciliation/75d_peak_annotation_distribution/` | UP peaks are 52% intergenic |
| Composite | `75_mecp2_signal_reconciliation/75_composite/` | Full panel |

### Section 76 — Triple-Overlap and Synapse-Specific
| Plot | File | Best for |
|------|------|----------|
| Violins with p-values | `76_triple_overlap_synapse_chromatin/76a_neuronal_chromatin_with_pvalues/` | Direct answer to "do we have p-values" |
| 4-group comparison | `76_triple_overlap_synapse_chromatin/76b_four_group_chromatin/` | Triple overlap in context |
| Synapse vs neuronal | `76_triple_overlap_synapse_chromatin/76c_synapse_vs_neuronal_chromatin/` | K27me3 specificity at synapse genes |
| K119ub decile enrichment | `76_triple_overlap_synapse_chromatin/76d_k119ub_decile_enrichment/` | Neuronal/synapse enriched in D10 |
| Composite | `76_triple_overlap_synapse_chromatin/76_composite/` | Full panel |

### Section 77 — MeCP2 Developmental Trajectory
| Plot | File | Best for |
|------|------|----------|
| Aging overview bar | `77_mecp2_aging_trajectory/77a_aging_overview/` | 2.1x more aging-UP peaks in mut |
| Overlap Venn | `77_mecp2_aging_trajectory/77b_aging_overlap_venn/` | 1,654 mut-unique aging genes |
| Mut-specific GO | `77_mecp2_aging_trajectory/77c_mut_specific_go_enrichment/` | 435 sig terms, 49 neuronal |
| Shared peak fold scatter | `77_mecp2_aging_trajectory/77d_shared_peak_fold_comparison/` | Mut ages 22% more than ctrl at shared loci |
| Composite | `77_mecp2_aging_trajectory/77_composite/` | Full panel |

### Section 78 — Stoichiometry with Unbiased Neuronal Set
| Plot | File | Best for |
|------|------|----------|
| Total methylation forest | `78_stoichiometry_broad_neuronal/78a_total_methylation_forest/` | Which groups gain vs lose total meth |
| Stoichiometry scatter | `78_stoichiometry_broad_neuronal/78b_stoichiometry_scatter/` | Neuronal colored in mc vs hmc scatter |
| Slope forest | `78_stoichiometry_broad_neuronal/78c_stoichiometry_slope_forest/` | Neuronal slope = -1.0 (stoichiometric!) |
| BAP1 vs TET-KO | `78_stoichiometry_broad_neuronal/78d_bap1_vs_tetko/` | TET-KO phenocopy weaker at neuronal |
| Bias demonstration | `78_stoichiometry_broad_neuronal/78e_narrow_vs_broad_bias/` | Direction flip: narrow +0.012 vs broad -0.002 |
| Quadrant characterization | `78_stoichiometry_broad_neuronal/78g_quadrant_characterization/` | 35% neuronal with broad set (was 21%) |
| Composite | `78_stoichiometry_broad_neuronal/78_composite/` | Full panel |

### Saved Tables (in `plots/visualizations/tables/`)
| File | Contents |
|------|----------|
| `74_geneset_overlap_counts.tsv` | 8 exclusive partition counts |
| `74_pairwise_fisher.tsv` | 3 pairwise Fisher tests |
| `74_neuronal_methylation_levels.tsv` | Per-gene methylation for neuronal set |
| `75_mecp2_peak_gene_summary.tsv` | Per-gene UP/DOWN/NS MeCP2 peak counts |
| `75_go_term_comparison.tsv` | K119ub vs MeCP2 neuronal term counts |
| `75_peak_annotation_distribution.tsv` | Annotation breakdown for UP vs DOWN peaks |
| `76_synapse_axon_gene_set.tsv` | 2,321 synapse/axon genes (reusable) |
| `76_triple_overlap_genes.tsv` | 20,915 genes with all flags |
| `76_chromatin_stats.tsv` | Kruskal-Wallis + pairwise Wilcoxon results |
| `76_synapse_vs_neuronal_stats.tsv` | Synapse vs broader neuronal test results |
| `76_top_decile_fisher.tsv` | D10 enrichment ORs |
| `77_aging_peak_summary.tsv` | UP/DOWN/NS counts per genotype |
| `77_aging_overlap.tsv` | Shared/unique peak and gene counts |
| `77_mut_specific_aging_go.tsv` | GO enrichment of mut-unique aging genes |
| `77_shared_peak_fold_comparison.tsv` | Paired fold values at shared peaks |
| `78_total_methylation_summary.tsv` | Delta total per gene group |
| `78_stoichiometry_slopes.tsv` | OLS slopes per group with CIs |
| `78_tetko_comparison.tsv` | BAP1-KO vs TET-KO correlations |
| `78_quadrant_neuronal_comparison.tsv` | Quadrant gene mark means by class |
