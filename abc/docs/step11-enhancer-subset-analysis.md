# Step 11: Enhancer Subset Stratified Analysis

## Position in the Pipeline

Step 11 sits downstream of the ABC enhancer-gene linkage pipeline (Steps 0-10) and the differential Hi-C loop analysis (mariner/edgeR). It integrates three independent data modalities:

1. **Differential Hi-C loops** (mariner pipeline): 39,344 non-redundant loops with logFC and FDR from edgeR quasi-likelihood GLM (n=3 biological replicates per condition)
2. **ABC enhancer-gene pairs** (Steps 6-7): 180,423 distal E-G pairs with delta-ABC scores from the Broad Institute ABC model
3. **Enhancer epigenomic classification** (DiffBind): 55,697 ATAC-defined enhancers classified by H3K27ac, H2AK119ub, and ATAC-seq differential signal

The central biological question: **Does H2AK119ub accumulation at enhancers in BAP1-KO alter chromatin looping and enhancer-gene linkage independently of H3K27ac/ATAC activity changes?** BAP1 is a deubiquitinase that removes H2AK119ub. If K119ub itself modulates 3D chromatin contacts, then enhancers that gain K119ub *without* losing activity marks should still show contact changes.

---

## Enhancer Classification

Four mutually exclusive phenotypic classes are defined from DiffBind analysis of the control ATAC enhancer universe (55,697 peaks with log-concentration K27ac >= 2.0). Binary significance flags (FDR < 0.05) determine membership:

| Class | N | Definition | Biological Meaning |
|-------|---|-----------|-------------------|
| **Activity_Lost** | 7,503 | K27ac down; ATAC/K119ub any | Classic activity loss — enhancer decommissioning in KO |
| **K119ub_Only** | 2,479 | K119ub up; K27ac and ATAC both unchanged | Key test class — K119ub accumulation without activity change |
| **Activity_Gain** | 2,851 | K27ac up; ATAC/K119ub any | Enhancer activation in KO — de-repression or compensatory |
| **Stable** | 42,864 | No significant changes in any mark | Background/control — unaffected by BAP1 loss |

The **K119ub_Only** class (4.5% of enhancers) is the central test: these enhancers accumulate H2AK119ub without measurable changes in chromatin accessibility (ATAC) or active enhancer marks (H3K27ac). Any contact-level phenotype in this class would support a direct role for K119ub in 3D genome organization, independent of transcriptional activity.

---

## Results and Interpretation

### Parts A-B: Loop Anchoring and Promoter Connectivity

**Loop anchoring rates differ by class** (KW chi2 = 244.1, p < 2.2e-16):

| Class | Median loops/enhancer | % with >= 1 loop | % with promoter loop |
|-------|----------------------|-------------------|---------------------|
| Activity_Lost | 1 | 60.0% | 35.6% |
| K119ub_Only | 1 | 63.5% | 43.3% |
| Activity_Gain | 0 | 47.8% | 25.0% |
| Stable | 1 | 57.7% | 37.8% |

K119ub_Only enhancers have the **highest** loop anchoring rate (63.5%) and the **highest** promoter-loop fraction (43.3%). This is not trivially expected from their definition — these are enhancers defined by K119ub gain, not by loop status. The elevated promoter connectivity suggests K119ub_Only enhancers are disproportionately located at functional regulatory loci that actively contact gene promoters.

Activity_Gain enhancers have the lowest anchoring (47.8%), consistent with many being distal or newly activated elements not yet engaged in stable loops.

**Gene expression at promoter-connected targets** shows a clear gradient: Activity_Lost target genes have median log2FC = -0.112 (downregulated), K119ub_Only targets show median log2FC = -0.040 (weak downregulation), Activity_Gain targets show median log2FC = +0.089 (upregulated), and Stable targets are near zero (+0.009). The K119ub_Only vs Stable difference is highly significant (Wilcoxon p < 2.2e-16), but the effect size is small.

### Part C: Loop Strength Changes — The Core Finding

**K119ub_Only enhancers show significant loop weakening** (one-sample Wilcoxon vs 0: p < 2.2e-16, median logFC = -0.054, rank-biserial r = 0.346).

| Class | N loop pairs | Median logFC | % FDR < 0.05 |
|-------|-------------|-------------|--------------|
| Activity_Lost | 13,410 | **-0.088** | 26.3% |
| K119ub_Only | 5,372 | **-0.054** | 19.2% |
| Activity_Gain | 3,492 | **+0.066** | 28.7% |
| Stable | 70,395 | **-0.013** | 14.8% |

The ordering is internally consistent: Activity_Lost shows the strongest weakening, K119ub_Only is intermediate, Stable is near zero, and Activity_Gain shows strengthening. All pairwise comparisons are significant (Holm-adjusted Wilcoxon p < 1e-62 for K119ub_Only vs Activity_Lost, p < 1e-99 for K119ub_Only vs Stable).

**Critical assessment:** The median logFC of -0.054 for K119ub_Only is statistically unambiguous (p < 2.2e-16) but biologically modest. This corresponds to a ~3.7% reduction in loop strength (2^-0.054 = 0.963). For comparison, Activity_Lost enhancers show ~5.9% reduction (2^-0.088 = 0.941). Whether a 3.7% shift has functional consequences is not established by this analysis alone.

The effect size (r = 0.346) indicates a moderate shift in the distribution — the entire distribution of loop strengths at K119ub_Only enhancers is shifted negative, not just a few outliers. This supports a broad, mild effect rather than a few strongly affected loci.

### Part D: ABC Score Changes

**ABC scores show expected class-level patterns** (KW chi2 = 2349, p effectively 0):

| Class | N enhancers with ABC | Median delta-ABC | Median delta-unnorm |
|-------|---------------------|-----------------|-------------------|
| Activity_Lost | 3,289 | **-0.0085** | **-0.077** |
| K119ub_Only | 1,064 | **+0.0018** | **+0.008** |
| Activity_Gain | 890 | **+0.0097** | **+0.076** |
| Stable | 12,887 | **+0.0010** | **+0.013** |

Activity_Lost and Activity_Gain show large, directional ABC changes consistent with their definitions. K119ub_Only has a weakly positive median delta-ABC (+0.0018) and delta-unnorm (+0.008). Both are significantly different from Activity_Lost (p < 1e-120) but only weakly different from Stable (p = 0.0002 for delta-ABC, p = 0.0003 for delta-unnorm).

**An apparent paradox:** K119ub_Only enhancers show **negative** loop logFC (weakening) but **positive** delta-ABC (strengthening). This is not a contradiction. Loops and ABC scores measure different things:

- **Loop logFC**: Measures strength change at specific called loop interactions (HiCCUPS point interactions). Captures disruption of discrete 3D contacts.
- **Delta-ABC**: Measures change in the *normalized* enhancer-gene linkage score, which incorporates activity (ATAC signal) and contact across all E-G pairs within 5 Mb. For K119ub_Only, activity is unchanged by definition, so delta-ABC reflects contact changes as seen by the broader Hi-C matrix — including powerlaw-estimated contacts where real Hi-C data is sparse.

The positive delta-ABC likely reflects (1) the ABC normalization compressing real changes (documented in abc-analysis-context.md: unnormalized Δ(A×C) outperforms normalized ΔABC, 65.3% vs 58.8% concordance), and (2) the small positive shift in short-range contact (<50 kb) seen in the contact decay profile (Part E). Loop weakening occurs at specific structured interactions, while background contact may slightly increase as loops dissolve and contacts redistribute.

### Part E: Contact Decay Profile

The contact decay data reveals distance-dependent behavior:

| Distance | Activity_Lost delta | K119ub_Only delta | Stable delta |
|----------|-------------------|------------------|-------------|
| < 50 kb | -4.0e-4 | **+5.8e-4** | +4.2e-4 |
| 50-100 kb | -4.8e-4 | +8.3e-6 | +9.0e-5 |
| 100-200 kb | -4.6e-4 | +1.7e-5 | +7.0e-5 |
| 200 kb-1 Mb | -2.0e-4 to -1.4e-4 | +6.8e-5 to -3.6e-5 | +1.4e-5 to +1.1e-6 |
| 1-5 Mb | -7.2e-5 | +2.9e-5 | -4.5e-6 |

K119ub_Only shows a **short-range contact gain** (<50 kb, delta = +5.8e-4) that is larger than Stable (+4.2e-4) and opposite to Activity_Lost (-4.0e-4). At longer ranges, K119ub_Only contact changes are near zero. This is consistent with a model where loop weakening at K119ub_Only loci redistributes contact into the local neighborhood without large-scale contact loss.

Activity_Lost shows **consistent contact loss at all distances**, confirming these enhancers undergo genuine 3D reorganization.

### Part G: K119ub Dose-Response (Tertile Analysis)

2,479 K119ub_Only enhancers were split into tertiles by log2FC_K119ub magnitude:

| Tertile | K119ub log2FC range | N loops | Median loop logFC | Median delta-ABC |
|---------|-------------------|---------|------------------|-----------------|
| T1_low | [0.32, 0.55] | 1,868 | -0.054 | +0.0012 |
| T2_mid | (0.55, 0.68] | 1,732 | -0.045 | +0.0021 |
| T3_high | (0.68, 1.59] | 1,772 | **-0.064** | +0.0029 |

**Loop logFC trend:** KW significant (p = 3.9e-5), Wilcoxon T1 vs T3 significant (p = 0.028). But the trend is **non-monotonic**: T2_mid (-0.045) is *less* negative than T1_low (-0.054), with only T3_high (-0.064) showing the expected stronger weakening. Spearman rho for continuous K119ub vs loop logFC = -0.035 (p = 0.010) — statistically significant but **extremely weak** (R^2 ~ 0.1%).

**Delta-ABC trend:** KW p = 0.054 (borderline, **not significant** at alpha = 0.05). The direction is **positive** — higher K119ub correlates with *more positive* delta-ABC — which is the opposite of what a simple K119ub-weakens-contacts model would predict.

**Honest assessment:** The dose-response evidence is weak. The loop logFC gradient, while significant by KW, explains virtually none of the variance (rho = -0.035). The non-monotonic pattern (T2 less negative than T1) and the discordant ABC direction (positive, not negative) indicate that K119ub magnitude is, at best, a very minor determinant of contact changes among K119ub_Only enhancers. Other factors — genomic context, distance to target, local chromatin architecture — likely dominate.

### Part H: RNA-seq Concordance

This tests whether ABC changes at enhancers predict target gene expression changes.

| Class | N pairs | Concordant | Rate | vs 50% (binom p) |
|-------|---------|-----------|------|-------------------|
| Activity_Lost | 7,066 | 4,226 | **59.8%** | — |
| K119ub_Only | 2,851 | 1,409 | **49.4%** | **p = 0.55** |
| Activity_Gain | 1,463 | 906 | **61.9%** | — |
| Stable | 31,750 | 15,489 | **48.8%** | — |

**K119ub_Only concordance is indistinguishable from chance** (49.4%, binomial p = 0.55). It is also indistinguishable from Stable (48.8%, Fisher's p = 0.52, OR = 1.03). In contrast, Activity_Lost (59.8%) and Activity_Gain (61.9%) show significant concordance, and K119ub_Only is significantly *worse* than both (Fisher's p < 2.2e-16 and p = 5.5e-15, respectively).

Only 231/2,136 (10.8%) of K119ub_Only target genes are differentially expressed (padj < 0.05, |log2FC| > 0.5).

**This is the most important negative result in the analysis.** The contact-level changes at K119ub_Only enhancers — while statistically detectable — do **not** propagate to measurable gene expression changes. The 49.4% concordance is what you would expect if delta-ABC and gene log2FC were independent random variables. This contrasts sharply with Activity_Lost and Activity_Gain, where the ABC model's E-G predictions are functionally validated by RNA-seq.

**Possible explanations:**
1. **The contact changes are too small to be functional.** A 3.7% loop weakening may be below the threshold for transcriptional consequences.
2. **K119ub acts through mechanisms not captured by ABC.** The ABC model measures activity x contact. If K119ub affects chromatin compaction, nuclear compartmentalization, or transcription factor binding without altering ATAC-measurable accessibility, the ABC framework is the wrong tool.
3. **K119ub effects on expression may require co-occurring chromatin changes.** The K119ub_Only class explicitly excludes enhancers with K27ac/ATAC changes — perhaps K119ub's transcriptional impact requires concurrent activity mark changes.
4. **Statistical power.** With 2,851 concordance-testable pairs, the analysis has reasonable power to detect concordance rates of ~53%+. A very subtle effect (51-52%) could be missed, but would be biologically negligible anyway.

### Part I: HOMER BED Export

Seven BED files were exported for motif enrichment analysis on HPC:

- 4 class-level: Activity_Lost (7,503), K119ub_Only (2,479), Activity_Gain (2,851), Stable (42,864)
- 3 K119ub_Only tertiles: T1_low (827), T2_mid (826), T3_high (826)

Planned HOMER analyses:
- Each class vs Stable background (identifies class-specific motifs)
- T3_high vs T1_low background (tests whether high-K119ub enhancers have distinct TF binding motifs)

These have not yet been run. The motif analysis could reveal whether K119ub_Only enhancers are enriched for Polycomb response elements (PREs) or specific TF motifs that distinguish them from Stable enhancers.

---

## Summary of Claims and Evidence Strength

| Claim | Evidence | Strength |
|-------|----------|----------|
| K119ub_Only enhancers show loop weakening vs Stable | Wilcoxon p < 2.2e-16, median logFC = -0.054 | **Strong statistically**, modest effect size |
| K119ub_Only loop weakening is intermediate between Activity_Lost and Stable | Pairwise tests all p < 1e-62 | **Strong** — consistent ordering |
| Higher K119ub gain causes more loop weakening (dose-response) | Spearman rho = -0.035 (p = 0.01), non-monotonic tertiles | **Weak** — barely above noise, non-monotonic |
| K119ub_Only ABC changes predict gene expression | Concordance 49.4%, binomial p = 0.55 | **Not supported** — concordance at chance |
| K119ub_Only differs from Stable in functional output | Concordance: OR = 1.03, p = 0.52; DE fraction ~10.8% | **Not supported** |
| Contact changes are distance-dependent | Short-range gain, long-range neutral | **Moderate** — consistent with redistribution |

---

## Limitations

1. **ATAC-only ABC mode.** The ABC pipeline used ATAC signal as the sole activity measure because H3K27ac BAMs were unavailable. The ABC model benchmarks best with DNase + H3K27ac. ATAC-only mode may underestimate enhancer activity, particularly for elements where H3K27ac and ATAC diverge.

2. **Hi-C sparsity.** Approximately 47-49% of Hi-C bins at 5 kb resolution are zero. The ABC model uses powerlaw extrapolation for these bins. For K119ub_Only enhancers specifically, contact changes are small, and powerlaw-filled values could mask or distort real signal.

3. **Consensus enhancer universe.** Both WT and KO used the same 75,371 ATAC peaks. This is correct for ΔABC computation but means condition-specific peaks (present only in KO or WT) are excluded.

4. **K119ub_Only class size.** At 2,479 enhancers (4.5% of total), this class is small. Only 1,064 match to ABC pairs after self-promoter removal, and only 2,851 pairs merge successfully with RNA-seq data. The tertile analysis splits ~826 enhancers three ways, limiting power for dose-response detection.

5. **Single timepoint.** RNA-seq is from adult cerebellum. If K119ub effects on gene expression are developmental (manifesting during differentiation rather than in steady-state tissue), this analysis would miss them entirely.

6. **Loop-ABC methodological mismatch.** Loop logFC and delta-ABC measure fundamentally different quantities (discrete point interactions vs integrated E-G linkage). Their discordance for K119ub_Only enhancers may reflect this measurement difference rather than biology.

---

## Recommended Figures for Presentation

Step 11 tells a story with two acts: (1) K119ub_Only enhancers show a contact-level phenotype, and (2) that phenotype does not reach gene expression. Both acts are important — the positive finding establishes the biology, and the negative finding defines its limits. A good presentation shows both honestly.

### If You Have One Slide: Panel 12 (Summary Patchwork)

The `12_summary_patchwork/` figure is a 4-panel composite (violin of loop logFC, density overlay, unnormalized ABC boxplot, contact decay by distance). It puts the core Part C result and supporting evidence on a single slide. The top-left violin is the punchline: Activity_Lost and K119ub_Only are both shifted negative, Activity_Gain is shifted positive, Stable is near zero. The bottom-left density overlay shows the same data as continuous distributions, making the K119ub_Only-vs-Stable shift visible as a slight leftward displacement of the red curve relative to the light blue curve. The top-right boxplot confirms the ABC-level pattern, and the bottom-right contact decay shows the distance-dependent behavior (Activity_Lost loses contact at all distances, K119ub_Only has a short-range gain that decays quickly).

### If You Have Two Slides

**Slide 1 — Panel 05 (Loop logFC Violin) + Panel 16 (ABC-RNA Concordance Bar):**

Panel 05 is the single most important positive result: the violin plot of loop logFC by enhancer class shows the clean ordering Activity_Lost (-0.088) > K119ub_Only (-0.054) > Stable (-0.013) > Activity_Gain (+0.066), with K119ub_Only intermediate and clearly shifted from Stable (p < 2.2e-16). The four violins are color-coded and the medians are visually distinct. This is the figure that establishes K119ub as having a contact-level effect independent of activity mark changes.

Then immediately pair it with Panel 16, the ABC-RNA concordance stacked bar. This is the most important *negative* result: K119ub_Only concordance is 49.4% (the green/red split sits right at the 50% dashed line), indistinguishable from chance, while Activity_Lost (59.8%) and Activity_Gain (61.9%) show clear concordance. The visual contrast between the K119ub_Only bar (50/50 split) and the Activity_Lost/Activity_Gain bars (clear green majority) is striking and communicates the finding without any statistics needed.

Showing both panels together is essential for honest presentation. Panel 05 alone would imply K119ub_Only enhancers are functionally important. Panel 16 tempers that: the contact perturbation is real but does not propagate to measurable transcriptional changes.

**Slide 2 — Panel 06 (Loop logFC Density Overlay):**

Panel 06 is the cleanest visualization of the distributional shift. The overlaid density curves show that K119ub_Only (red) tracks Activity_Lost (blue) on the left tail but peaks slightly to the right, while Activity_Gain (orange) peaks clearly to the right of zero. The visual makes the "intermediate phenotype" argument intuitive. This is the figure best suited for a paper main figure or a slide where you want to show the full distribution shapes rather than just medians.

### Why Not the Others?

- **Panels 01-04** (Parts A-B: loop anchoring, distance, promoter loops, gene logFC) are characterization panels that set the stage but don't carry the central finding. Panel 04 (gene logFC by class) is useful if asked about expression effects, but shows small effect sizes (median log2FC = -0.040 for K119ub_Only)
- **Panels 07-08** (delta-ABC, delta-unnorm boxplots) show the paradoxical positive delta-ABC for K119ub_Only, which requires careful explanation of the loop-vs-ABC methodological mismatch. Not ideal for a talk unless you plan to discuss this explicitly
- **Panel 09** (ABC directionality stacked bar) is informative but secondary — shows the gained/lost/unchanged proportions per class
- **Panels 10-11** (contact decay) are important for the paper (distance-dependent behavior) but are secondary for a talk. Panel 11 is included in the summary patchwork (Panel 12) already
- **Panels 13-15** (Part G: K119ub tertile analysis) show the dose-response, which is the *weakest* result (rho = -0.035, non-monotonic tertiles). Panel 13 (tertile violin) visually shows T2_mid is *less* negative than T1_low, which undermines the dose-response narrative. Panel 15 (continuous scatter) is a cloud with an essentially flat regression line. These are important for completeness and the paper, but showing them in a talk invites criticism of the dose-response claim without adding to the story
- **Panel 17** (delta-ABC vs gene logFC scatter) is visually noisy — the four classes overlap heavily and the K119ub_Only points (pink) are indistinguishable from the Stable background (grey). Better as a supplementary figure
- **Panel 18** (K119ub_Only target gene logFC histogram) shows that most target genes are near zero log2FC, reinforcing the negative finding from Panel 16 but less efficiently

### The Two-Slide Narrative

> **Slide 1:** "K119ub_Only enhancers — which gain H2AK119ub without losing activity marks — show significant loop weakening (Panel 05: median logFC = -0.054, p < 2.2e-16), intermediate between Activity_Lost and Stable. This supports a direct contact-level role for K119ub independent of transcriptional activity."

> **Slide 2:** "However, this contact perturbation does not propagate to detectable gene expression changes (Panel 16: 49.4% concordance = chance, p = 0.55). The K119ub-driven loop weakening is real but below the functional threshold needed to alter transcription in the adult cerebellum."

---

## Output File Descriptions

### Plots (18 panels, each in PDF + SVG + JPG subfolders)

| # | File | Part | Content |
|---|------|------|---------|
| 01 | `01_loops_per_enhancer_boxplot/` | A | Boxplot of loop count per enhancer by class |
| 02 | `02_loop_distance_violin/` | A | Violin plot of loop distance by class |
| 03 | `03_promoter_loop_proportion/` | B | Stacked bar: promoter vs non-promoter loops |
| 04 | `04_gene_logfc_by_class/` | B | Boxplot of promoter-connected gene log2FC |
| 05 | `05_loop_logfc_violin/` | C | Violin of loop logFC by class (core result) |
| 06 | `06_loop_logfc_density/` | C | Overlaid density of loop logFC distributions |
| 07 | `07_delta_abc_boxplot/` | D | Boxplot of mean delta-ABC per enhancer |
| 08 | `08_delta_unnorm_boxplot/` | D | Boxplot of unnormalized Δ(Activity x Contact) |
| 09 | `09_abc_direction_stacked_bar/` | D | Gained/lost/unchanged ABC pair proportions |
| 10 | `10_contact_decay_wt/` | E | WT contact frequency vs distance by class |
| 11 | `11_contact_decay_delta/` | E | Delta contact vs distance by class |
| 12 | `12_summary_patchwork/` | F | Multi-panel summary (plots 05, 06, 08, 11) |
| 13 | `13_k119ub_tertile_loop_logfc/` | G | Loop logFC by K119ub tertile (violin) |
| 14 | `14_k119ub_tertile_delta_abc/` | G | Delta-ABC by K119ub tertile (boxplot) |
| 15 | `15_k119ub_vs_loop_logfc_scatter/` | G | Continuous K119ub vs loop logFC scatter |
| 16 | `16_abc_rnaseq_concordance_bar/` | H | Concordance rate stacked bar by class |
| 17 | `17_delta_abc_vs_gene_logfc_scatter/` | H | Delta-ABC vs gene log2FC scatter |
| 18 | `18_k119ub_target_gene_logfc_hist/` | H | K119ub_Only target gene log2FC histogram |

### TSV Files (10)

| File | Rows | Description |
|------|------|-------------|
| `class_level_summary.tsv` | 4 | One row per class: all metrics aggregated |
| `statistical_tests.tsv` | 30 | All KW, Wilcoxon, Spearman, Fisher's, binomial tests |
| `enhancer_class_loop_metrics.tsv` | 55,697 | Per-enhancer: loop count, distance, promoter loops |
| `enhancer_class_abc_metrics.tsv` | 18,130 | Per-enhancer: ABC pair count, delta-ABC, gained/lost |
| `contact_decay_by_class.tsv` | 24 | Mean contact by distance bin and class |
| `promoter_loop_gene_logfc.tsv` | 26,377 | Enhancer-gene pairs via promoter loops with DE |
| `k119ub_tertile_assignments.tsv` | 2,479 | K119ub_Only enhancers with tertile labels |
| `k119ub_tertile_loop_summary.tsv` | 3 | Per-tertile: loop and ABC metrics |
| `abc_rnaseq_concordance_by_class.tsv` | 4 | Per-class concordance rates |
| `k119ub_only_target_gene_de.tsv` | 2,136 | K119ub_Only target genes with DE status |

### BED Files for HOMER (7)

| File | Regions | Purpose |
|------|---------|---------|
| `Activity_Lost.bed` | 7,503 | Class-level motif enrichment |
| `K119ub_Only.bed` | 2,479 | Class-level motif enrichment |
| `Activity_Gain.bed` | 2,851 | Class-level motif enrichment |
| `Stable.bed` | 42,864 | Background for motif enrichment |
| `K119ub_Only_T1_low.bed` | 827 | Tertile-level motif analysis |
| `K119ub_Only_T2_mid.bed` | 826 | Tertile-level motif analysis |
| `K119ub_Only_T3_high.bed` | 826 | Tertile-level motif analysis |

---

## Bottom Line

**What the data shows:** K119ub_Only enhancers exhibit a statistically significant but biologically modest loop weakening phenotype (median logFC = -0.054, ~3.7% reduction). This effect is intermediate between Activity_Lost and Stable, consistent with K119ub contributing to contact changes independently of activity mark loss.

**What the data does not show:** This contact-level phenotype does not translate to detectable gene expression changes. ABC-RNA concordance for K119ub_Only (49.4%) is indistinguishable from chance and from the Stable control class. The dose-response of K119ub magnitude to contact changes is present but extremely weak (rho ~ -0.035) and non-monotonic. The ABC score direction is paradoxically positive, likely reflecting normalization artifacts and contact redistribution.

**Interpretation:** K119ub accumulation at enhancers produces a detectable perturbation in chromatin looping that is insufficient — in the adult cerebellum steady-state — to drive measurable transcriptional changes. The effect may be (a) below a functional threshold, (b) compensated by other regulatory mechanisms, (c) relevant only during developmental transitions when chromatin states are being established rather than maintained, or (d) acting through mechanisms (compartmentalization, phase separation) not captured by loop-level or ABC-level metrics. The HOMER motif analysis (Part I, pending) may help distinguish these possibilities by revealing whether K119ub_Only enhancers have distinct regulatory grammar.
