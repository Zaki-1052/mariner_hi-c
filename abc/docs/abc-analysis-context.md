# ABC Model Analysis — Context Document

> **Date:** 2026-02-18
> **Purpose:** Reference document summarizing QC findings, pipeline decisions, and results from Steps 5–9 of the ABC enhancer-gene linkage analysis (BAP1-KO vs WT mouse cerebellum).
> **Companion documents:** `abc-prd-v2.md` (requirements), `abc-execution-plan.md` (implementation)

---

## 1. Pipeline Configuration Summary

- **Organism:** Mouse mm10, GENCODE vM25
- **Tool:** Broad Institute ABC v1.1.2 (Snakemake pipeline)
- **Conditions:** WT_cerebellum (ctrl) and KO_cerebellum (BAP1-KO mutant)
- **Mode:** ATAC-only (no H3K27ac BAMs available; only bigWigs exist)
- **Hi-C:** Condition-specific `.hic` files at 5kb resolution
- **Consensus peaks:** 75,371 peaks from `consensus_all.bed` (union of Batch 1 ctrl/mut consensus ATAC peaks, ≥2/3 replicates within genotype, merged across genotypes)
- **Key config changes from default:** `genome_size: mm`, `use_qnorm: False`, `threshold: 0.02`
- **QNorm disabled because:** Reference is human K562 DHS — cross-species AND cross-assay for mouse ATAC
- **Threshold set explicitly to 0.02** (not auto) because disabling qnorm invalidates auto-calibration

---

## 2. Reference File Statistics

| File | Count | Notes |
|------|-------|-------|
| mm10_genes.bed | 28,230 genes | 21,833 protein_coding, 5,620 lincRNA, 777 processed_transcript |
| mm10_tss.bed | 28,230 TSS regions | All 500bp, matches gene count |
| mm10_merged_blacklist.bed | 3,581 regions | Lab (254) + ENCODE (3,435), merged |
| mm10_ubiquitous_genes.tsv | 730 genes | Top 5% by baseMean from RNA-seq |
| consensus_all.narrowPeak | 75,371 peaks | 21 standard chromosomes |

### Ubiquitous gene list notes
- Actb: present (rank 575/16,572, baseMean=2987.2)
- Gapdh: **absent** (rank 1,054, baseMean=1797.2) — just below 95th percentile cutoff, not a bug
- Rpl13a: **absent** (rank 11,925, baseMean=23.5) — unexpectedly low expression in this cerebellum RNA-seq dataset, possibly a mapping issue with ribosomal genes
- The ubiquitous list is used only for flagging, not filtering — does not affect ABC scores

---

## 3. ABC Pipeline Output Summary (Step 5 QC)

### Per-condition thresholded predictions (ABC ≥ 0.02)

| Metric | WT_cerebellum | KO_cerebellum |
|--------|---------------|---------------|
| Total E-G pairs | 122,212 | 119,965 |
| Self-promoter entries | 29,152 | 29,152 |
| Distal E-G pairs | 93,060 | 90,813 |
| Unique genes | 27,492 | 27,492 |
| Mean enhancers/gene (all) | 4.44 | 4.36 |
| Mean distal enhancers/gene | 3.38 | 3.30 |
| Median ABC score (all) | 0.0451 | 0.0457 |
| Median ABC score (distal) | 0.0343 | 0.0344 |
| Median distance (distal) | 127,108 bp | 115,795 bp |
| Chromosomes represented | 20 (chr1-19, chrX) | 20 |

### Key observations
- **100% gene overlap** between conditions (27,492 genes in both) — confirms consensus peak injection worked correctly
- **2,247 fewer E-G pairs in KO** (1.8% reduction) — modest net loss of enhancer-gene connections crossing threshold
- Element class distribution (KO): 60,311 genic / 30,502 intergenic / 29,152 promoter (~2:1 genic:intergenic)
- WT class distribution confirmed identical pattern: 61,673 genic / 31,387 intergenic / 29,152 promoter

### AllPutative (unthresholded) statistics

| Metric | WT | KO |
|--------|----|----|
| Total E-G pairs | 8,318,699 | 8,220,686 |
| Genes in AllPutative | 19,536 | 19,477 |

- AllPutative contains fewer genes (~19.5k) than thresholded (~27.5k) because AllPutative only includes expressed genes, while thresholded includes forced self-promoters for all genes

### Per-gene ABC score sums (AllPutative)

| Metric | WT | KO |
|--------|----|----|
| Mean sum (all entries) | 1.888 | 1.877 |
| Mean sum (excluding self-promoters) | 0.806 | 0.795 |
| Genes outside [0.90, 1.05] excl. self-promoter | 57.7% | 59.9% |

- Self-promoters are forced to ABC=1.0, so all-entry sums of ~1.89 (not 2.0) indicate distal components average ~0.80
- Substantial gene-to-gene variation: some genes are promoter-dominated (sum_distal ~0.5), others are enhancer-driven (sum_distal ~0.99)
- 1,536 genes have 2 self-promoter entries (~8%), 36 have 3+ — caused by overlapping consensus ATAC peaks at promoter regions

### Hi-C sparsity

| Distance bin | Zero hic_contact (WT) |
|--------------|-----------------------|
| <50kb | 0.1% |
| 50–200kb | 1.1% |
| 200–500kb | 7.0% |
| 500kb–1Mb | 22.0% |
| 1–2Mb | 43.4% |
| 2–5Mb | 65.1% |
| **Overall** | **47.3% (WT) / 49.1% (KO)** |

- **Critical finding:** overall sparsity is ~47-49%, but this is dominated by long-range pairs
- **Thresholded distal pairs (ABC ≥ 0.02):** only 3.1–3.4% have zero Hi-C contact
- 98.7% of thresholded pairs have real Hi-C data in at least one condition
- Powerlaw-only pairs have lower mean ABC (0.031-0.032) vs Hi-C-backed pairs (0.053-0.054)
- Conclusion: Hi-C sparsity is not a concern for the differential analysis — the pairs that matter are well-covered

---

## 4. Delta-ABC Results (Step 6)

### Join statistics
- Both conditions: 7,942,829 pairs (92.9% overlap)
- WT only: 354,726 / KO only: 256,770
- After threshold (ABC ≥ 0.02 in either condition): **180,423 distal E-G pairs**
- Unique genes: 20,265

### ΔABC distribution

| Statistic | Value |
|-----------|-------|
| Mean | -0.000459 |
| Median | -0.000042 |
| Std | 0.031158 |
| Q5 / Q95 | -0.03694 / +0.03488 |
| Gained (ΔABC > 0.01) | 44,060 |
| Lost (ΔABC < -0.01) | 47,499 |
| Unchanged (|Δ| ≤ 0.01) | 88,864 |

- Slightly negative mean — modest net loss of ABC scores in KO
- 7.8% excess of lost over gained links

### ΔABC by distance bin

| Distance | n | mean Δ | gained | lost | lost/gained |
|----------|---|--------|--------|------|-------------|
| <50kb | 44,234 | -0.00194 | 7,825 | 8,002 | 1.02 |
| 50–200kb | 47,044 | -0.00250 | 7,375 | 9,792 | 1.33 |
| 200–500kb | 32,060 | -0.00302 | 6,522 | 9,181 | 1.41 |
| 500kb–1Mb | 18,721 | -0.00353 | 4,851 | 6,791 | 1.40 |
| 1–5Mb | 31,720 | -0.00318 | 10,843 | 13,733 | 1.27 |

- Loss of E-G connections preferentially affects mid-to-long range interactions (50kb–1Mb)
- Short-range (<50kb) is near-symmetric (ratio 1.02)

### Unnormalized Δ(Activity × Contact)
- Mean: -0.003230
- Median: +0.003745
- Sign disagreement between mean and median indicates skewed distribution

### Hi-C coverage in thresholded pairs
- With Hi-C in at least one condition: 178,026/180,423 (98.7%)
- Powerlaw-only: 2,397

---

## 5. RNA-seq Integration Results (Step 7)

### Gene matching
- ABC genes: 20,265
- RNA-seq genes: 16,572
- Matched: 13,588 (67.1% of ABC genes)
- Unmatched ABC genes are likely lincRNAs/processed transcripts not in DESeq2 output
- E-G pairs with RNA-seq: 113,453

### Gene-level concordance (normalized ΔABC)

| Metric | Value |
|--------|-------|
| Genes with \|ΔABC\| > 0.01 | 12,665 |
| Genes with DE (padj<0.05, \|LFC\|>0.5) | 978 |
| **Genes with both** | **940** |
| **Concordant** | **553 (58.8%)** |
| Discordant | 387 (41.2%) |
| **Binomial p-value** | **6.84e-08** |

### 2×2 directional breakdown

| | Upregulated | Downregulated |
|---|---|---|
| **Gained enhancer** | 154 | 252 |
| **Lost enhancer** | 135 | 399 |

- Lost+downregulated (399) is the dominant quadrant (42% of all)
- Total downregulated: 651. Total upregulated: 289. Ratio 2.25:1
- Gained+downregulated (252) is the largest discordant group

### Unnormalized concordance (Δ(A×C))
- Concordance: **595/911 (65.3%)**
- Binomial p-value: **1.64e-20**
- Substantially better than normalized (58.8%) — unnormalized score better predicts expression direction
- Among the 252 "gained(norm)+downregulated" discordant genes, only 31 (12.3%) have negative unnormalized delta — most discordance is genuine, not a normalization artifact

### NaN values in top dysregulated genes
- Several top genes (Cbr1, Dnah1, Tigit, etc.) had NaN for `top_enh_distance` and `top_enh_class`
- Cause: KO-only pairs in the outer join lacked shared columns (distance, class) because Step 6 only included SHARED_COLS in the WT slim, not KO slim
- Fix applied: added `SHARED_COLS` to `ko_cols` in Step 6, then re-ran Steps 6 and 7

---

## 6. Cross-Reference Results (Step 9)

### 9A: Differential loop overlap
- Total loops: 2,910 (from characterized_loops.tsv)
- Significant loops: (subset used for analysis)
- ABC enhancers overlapping anchor1: 1,394/33,004 (4.2%)
- ABC enhancers overlapping anchor2: 1,330/33,004 (4.0%)
- ABC target genes in loop annotations: 2,433/20,265 (12.0%)
- **Directional concordance (loop vs ΔABC): 1,209/2,352 (51.4%)** — essentially at chance
- Loop annotations use nearest-gene assignment, which is inherently noisy. ABC's activity×contact scoring is more principled, so weak concordance reflects methodology mismatch, not model failure

### 9B: H3K27ac peak overlap
- ABC enhancers overlapping H3K27ac: 12,993/33,004 (39.4%)
- ABC enhancers without H3K27ac: 20,011/33,004 (60.6%)
- E-G pairs at H3K27ac+ enhancers: 98,440
- E-G pairs at H3K27ac- enhancers: 81,983
- Mean |ΔABC| at H3K27ac+ enhancers: 0.01805
- Mean |ΔABC| at H3K27ac- enhancers: 0.01598

### H3K27ac-stratified gain/loss asymmetry

| | Gained | Lost | Lost/Gained | Mean ΔABC |
|---|---|---|---|---|
| **H3K27ac+** | 24,346 | 28,409 | **1.17** | **-0.00102** |
| **H3K27ac-** | 19,714 | 19,089 | **0.97** | +0.00021 |

- The gain/loss asymmetry is concentrated at H3K27ac-marked enhancers
- H3K27ac- enhancers are essentially symmetric (no net directional change)

---

## 7. Output File Inventory

All under `/expanse/lustre/projects/csd940/zalibhai/abc/results/`:

| File | Rows | Description |
|------|------|-------------|
| `delta_abc_all_pairs.tsv` | 180,423 | All thresholded distal E-G pairs with ΔABC, components, Hi-C flags |
| `delta_abc_significant.tsv` | 91,559 | Subset with \|ΔABC\| > 0.01 |
| `delta_abc_with_rnaseq.tsv` | 113,453 | E-G pairs merged with RNA-seq (log2FC, padj, baseMean) |
| `gene_level_summary.tsv` | 13,588 | One row per gene: strongest ΔABC, aggregate stats, DE results, dysregulation flag |
| `delta_abc_annotated.tsv` | 180,423 | E-G pairs annotated with `has_H3K27ac` flag |
| `loops_with_gene_assignments.tsv` | — | Loops with nearest-gene annotations |
| `figures/delta_abc_vs_log2fc.pdf` | — | Scatter: ΔABC vs log2FC (normalized and unnormalized panels) |

### Column inventory for `delta_abc_all_pairs.tsv`
1. chr, start, end — enhancer coordinates
2. TargetGene — gene symbol
3. distance — enhancer-to-TSS distance
4. class — genic/intergenic/promoter
5. ABC.Score_WT, ABC.Score_KO — normalized scores per condition
6. ABC.Score.Numerator_WT, ABC.Score.Numerator_KO — unnormalized A×C
7. activity_base_WT, activity_base_KO — ATAC activity component
8. hic_contact_pl_scaled_adj_WT, hic_contact_pl_scaled_adj_KO — adjusted Hi-C contact
9. hic_contact_WT, hic_contact_KO — raw Hi-C contact
10. delta_ABC, delta_unnorm — KO minus WT deltas
11. hic_data_wt, hic_data_ko, has_hic_either — Hi-C coverage flags

---

## 8. Known Limitations and Caveats

1. **ATAC-only mode**: Activity defined purely from ATAC signal, not H3K27ac (BAMs unavailable). H3K27ac peaks used only for downstream annotation, not ABC scoring.
2. **Hi-C sparsity at long range**: ~65% of 2-5Mb pairs have zero raw Hi-C contact. For these, ABC falls back to powerlaw pseudocount (identical between conditions), so long-range ΔABC reflects only ΔACTIVITY, not Δcontact. This primarily affects the AllPutative universe; 98.7% of thresholded pairs have real Hi-C data.
3. **~8% of genes have 2+ self-promoter entries** due to overlapping consensus ATAC peaks at promoter regions. This slightly inflates their ABC denominators, deflating distal enhancer relative scores. Effect cancels in delta-ABC (same peaks both conditions).
4. **Unnormalized Δ(A×C) outperforms normalized ΔABC** for predicting DE direction (65.3% vs 58.8% concordance). The normalization compresses real differences when BAP1-KO causes widespread activity changes across many enhancers for the same gene.
5. **Loop-ABC directional concordance is at chance (51.4%)** due to nearest-gene assignment in loop annotations vs ABC's principled E-G linking. Not a failure of either dataset, but the comparison methodology is limited.
6. **H3K27ac peaks are from a single merged dataset** (H3K27acCerebellumLate2.bed, adult timepoint), not condition-specific. The overlap annotation tells you which enhancers had active marks in the reference, not whether those marks changed between WT and KO.

---

## 9. Pending Work

- **Goal 2 (deferred):** Correlate ΔABC / Δ(A×C) with H2AK119ub levels at enhancers using H2AK119ub bigWigs (4 WT + 4 KO replicates in `heatmaps/`). Prediction from current findings: enhancers with lost ABC score in KO should show increased H2AK119ub signal.
- **Scatter plot inspection:** The `delta_abc_vs_log2fc.pdf` has not been visually reviewed in this conversation. Contains Pearson and Spearman correlations for DE genes.
- **Paired-anchor loop analysis:** A more rigorous loop-ABC comparison would check whether enhancer overlaps anchor1 AND gene TSS overlaps anchor2 of the *same* loop, rather than independent anchor overlap.
