# Plan: Sections 74-77 — Addressing Jai's Follow-Up Analysis Requests

## Context

Jai raised 5 groups of follow-up questions about the MeCP2/BAP1-KO/K119ub theory after reviewing sections 61-73 results and the 72g chromatin remodeling figure. These questions probe: (1) gene set overlaps and neuronal methylation levels, (2) the MeCP2 signal direction paradox, (3) p-values and chromatin analysis for triple-overlap and synapse-specific genes, (4) synapse/axon vs broader neuronal gene distinction, and (5) young vs adult MeCP2 binding trajectory. Four new section scripts address all points.

## Data Availability

**Available now (sections 74-76):**
- `tables/59_quadrant_master.tsv` — 23,150 genes (K119ub + MeCP2 + methylation + chromatin)
- `tables/diffbind_gene_level_all_marks.tsv` — 20,916 genes (ATAC/K27ac/K27me3/K119ub fold + FDR)
- `tables/mecp2_coordinated_genes.tsv` — 5,357 coordinated genes (mC↑ AND hmC↓)
- `tables/72_neuronal_gene_set_go_derived.tsv` — 5,614 neuronal genes (GO-derived)
- `peaks/mecp2/MeCP2_annotated.txt` — 202,650 peaks (Fold/FDR/SYMBOL)
- `tables/61k_gsea_mecp2_go_bp.tsv` + `61k_gsea_k119ub_go_bp.tsv` — pre-computed GSEA
- `data/k119ub_gene_signal.tsv` — K119ub gene body signal

**Blocked on Jai (section 77):**
- MeCP2 DiffBind: ctrl young vs adult (13,752 FDR<0.05 per mecp2.png)
- MeCP2 DiffBind: mut young vs adult (33,763 FDR<0.05 per mecp2.png)

---

## Section 74: Gene Set Overlap & Neuronal Methylation Levels

**File:** `section_74_geneset_overlap_methylation.R`
**Addresses:** Jai's point 1 — "create a venn diagram / counts table... then make a plot for neuronal genes showing total methylation ctrl vs mut, and 5mC and 5hmC levels"

**Inputs:** `59_quadrant_master.tsv`, `mecp2_coordinated_genes.tsv`, `72_neuronal_gene_set_go_derived.tsv`

### Panels

| Panel | Content | Stats |
|-------|---------|-------|
| 74a | UpSet-style intersection plot — all 7 exclusive partitions of Neuronal ∩ Coordinated ∩ MeCP2-Up, with gene counts annotated | 3 pairwise Fisher tests (N×C, N×M, C×M), BH-corrected |
| 74b | Total methylation (mC+hmC) ctrl vs mut for neuronal genes — paired violin+box | Wilcoxon signed-rank (paired by gene) |
| 74c | 5mC levels ctrl vs mut for neuronal genes — paired violin+box | Wilcoxon signed-rank |
| 74d | 5hmC levels ctrl vs mut for neuronal genes — paired violin+box | Wilcoxon signed-rank |
| 74e | Patchwork composite: `74a / (74b | 74c | 74d)` | — |

### Key Logic
- Universe: 59_quadrant_master (has mc_ctrl/mc_mut/hmc_ctrl/hmc_mut)
- MeCP2-Up: `mecp2_nearest_fdr < 0.05 & mecp2_mean_fold > 0` (yields ~79 genes)
- Neuronal: gene ∈ 72_neuronal_gene_set_go_derived.tsv (5,614 genes)
- Coordinated: gene ∈ mecp2_coordinated_genes.tsv (5,357 genes)
- total_ctrl = mc_ctrl + hmc_ctrl; total_mut = mc_mut + hmc_mut
- Build UpSet manually with ggplot (geom_col + geom_point matrix) — no UpSetR dependency

### Saved Tables
- `74_geneset_overlap_counts.tsv` — 8-row exclusive partition counts
- `74_pairwise_fisher.tsv` — Fisher test results
- `74_neuronal_methylation_levels.tsv` — per-gene methylation values for neuronal set

---

## Section 75: MeCP2 Signal Direction Reconciliation

**File:** `section_75_mecp2_signal_reconciliation.R`
**Addresses:** Jai's point 2 — "how does MeCP2 signal drop in the mutant when DiffBind showed ~8000 loci with increased binding?" + K119ub vs MeCP2 neuronal GO term comparison

**Inputs:** `MeCP2_annotated.txt`, `59_quadrant_master.tsv`, `61k_gsea_*.tsv`, `72_neuronal_gene_set_go_derived.tsv`

### Panels

| Panel | Content | Stats |
|-------|---------|-------|
| 75a | Peak-to-gene aggregation — per-gene counts of UP/DOWN/NS MeCP2 peaks; bar chart showing how 7,686 UP peaks distribute across genes vs 1,200 DOWN peaks | Wilcoxon: gene-body signal at UP-peak genes still drops? |
| 75b | Gene-body MeCP2 signal (from quadrant master `mecp2_mean_fold`) partitioned by peak status: genes with UP peaks vs DOWN peaks vs no sig peaks — violin+box | Kruskal-Wallis + pairwise Wilcoxon |
| 75c | GO term count comparison: K119ub sig terms vs MeCP2 sig terms; neuronal subset highlighted; specific term listing | Fisher: neuronal fraction in K119ub vs MeCP2 sig terms |
| 75d | MeCP2 UP vs DOWN peak annotation distribution (Promoter/Intron/Intergenic from annotation column) — stacked bar | Chi-squared |
| 75e | Patchwork composite: `(75a | 75b) / (75c | 75d)` | — |

### Key Logic
- Load MeCP2_annotated.txt (202k peaks). Aggregate to gene level: per SYMBOL, count peaks with FDR<0.05 & Fold>0 (UP), FDR<0.05 & Fold<0 (DOWN), rest (NS)
- Join with 59_quadrant_master for mecp2_mean_fold (gene-body proxy)
- Narrative: MeCP2 redistributes — concentrates at ~X thousand loci (UP peaks) while dropping at the remaining ~Y thousand. Gene-body signal reflects the net effect.
- For 75c: load pre-computed GSEA tables, filter q<0.05, classify neuronal terms by pattern `synap|neuron|axon|dendrit|nervous`

### Saved Tables
- `75_mecp2_peak_gene_summary.tsv` — per-gene UP/DOWN/NS peak counts + mecp2_mean_fold
- `75_go_term_comparison.tsv` — K119ub vs MeCP2 neuronal term counts
- `75_peak_annotation_distribution.tsv` — annotation category counts for UP vs DOWN

---

## Section 76: Triple-Overlap & Synapse-Specific Chromatin Analysis

**File:** `section_76_triple_overlap_synapse_chromatin.R`
**Addresses:** Jai's points 3+4 — "p-values for the violins... same analysis for MeCP2/K119ub/coordinated overlap... synapse/axon genes might be special vs neuronal overall"

**Inputs:** `diffbind_gene_level_all_marks.tsv`, `59_quadrant_master.tsv`, `72_neuronal_gene_set_go_derived.tsv`, `mecp2_coordinated_genes.tsv`, `data/k119ub_gene_signal.tsv`, `org.Mm.eg.db`

### Panels

| Panel | Content | Stats |
|-------|---------|-------|
| 76a | Reproduce 73a violins (ATAC/K27ac/K27me3 neuronal vs non-neuronal) with explicit p-value annotations and group sizes | Wilcoxon rank-sum per mark (3 tests) |
| 76b | 4-group chromatin comparison: triple-overlap vs neuronal-only vs coordinated-only vs rest, for ATAC/K27ac/K27me3 | Kruskal-Wallis + pairwise Wilcoxon (BH) |
| 76c | Synapse/axon (GO pattern `synap|axon`) vs broader neuronal vs non-neuronal — chromatin fold changes | Pairwise Wilcoxon (BH) |
| 76d | K119ub decile membership + top-decile enrichment — where do triple-overlap and synapse/axon genes sit? Forest plot of ORs | Fisher exact test per gene set × mark |
| 76e | Patchwork composite: `(76a) / (76b | 76c) / 76d` | — |

### Key Logic
- Triple-overlap: neuronal AND coordinated AND MeCP2-Up (~12-16 genes based on section 61e)
- Synapse/axon subset: fresh derivation from org.Mm.eg.db with pattern `synap|axon` only (stricter than the neuronal set's `synap|neuron|axon|dendrit|nervous`)
- For 76d: define top decile per mark by |fold_change|, then Fisher test for enrichment of each gene set
- Merge diffbind_gene_level_all_marks + 59_quadrant_master (for MeCP2 status) + K119ub signal + neuronal/coordinated flags
- Pairwise Wilcoxon with `p.adjust(method = "BH")` for multi-group comparisons (no extra packages needed)

### Saved Tables
- `76_synapse_axon_gene_set.tsv` — derived synapse/axon gene set (reusable downstream)
- `76_triple_overlap_genes.tsv` — gene-level table with all flags
- `76_chromatin_stats.tsv` — all statistical test results
- `76_top_decile_fisher.tsv` — Fisher OR/p for top-decile enrichment

---

## Section 77: MeCP2 Developmental Trajectory (Young vs Adult)

**File:** `section_77_mecp2_aging_trajectory.R`
**Addresses:** Jai's point 5 — "account for age-related increases... find loci uniquely increased in mut aging... what % overlap"

**BLOCKED:** Requires data from Jai. Script will include `stopifnot()` with clear error message.

### Data Request for Jai

Two annotated DiffBind result files, placed at `peaks/mecp2/`:

| File | Contrast | Expected Sig |
|------|----------|-------------|
| `MeCP2_ctrl_young_vs_adult_diffbind.txt` | Control: young vs adult | ~13,752 |
| `MeCP2_mut_young_vs_adult_diffbind.txt` | Mutant: young vs adult | ~33,763 |

**Required columns** (matching existing `MeCP2_annotated.txt` schema):
```
seqnames  start  end  width  strand  Conc  Conc_adult  Conc_young  Fold  p.value  FDR
annotation  geneChr  geneStart  geneEnd  geneLength  geneStrand  geneId  transcriptId
distanceToTSS  ENSEMBL  SYMBOL  GENENAME
```

**Minimum viable:** If ChIPseeker annotation isn't available, a 6-column TSV is sufficient:
`seqnames  start  end  Fold  p.value  FDR` — the script can annotate using TxDb.Mmusculus.UCSC.mm10.knownGene.

**Critical:** Must include ALL tested peaks (not just significant) to enable volcano plots, GSEA, and proper background calculations.

### Message for Jai

> Hey Jai, for the young vs adult MeCP2 analysis — could you export the two DiffBind result tables from your young-vs-adult comparisons? Specifically:
>
> 1. **Control young vs adult**: all tested peaks (not just significant), as a TSV
> 2. **Mutant young vs adult**: same format
>
> Ideally in the same annotated format as the existing `MeCP2_annotated.txt` (seqnames/start/end/Fold/p.value/FDR/SYMBOL/distanceToTSS/annotation etc). But if you don't have the ChIPseeker annotation, just the core columns work too: `seqnames, start, end, Fold, p.value, FDR` — I can annotate on my end.
>
> Name them `MeCP2_ctrl_young_vs_adult_diffbind.txt` and `MeCP2_mut_young_vs_adult_diffbind.txt` and I'll drop them in `peaks/mecp2/`.

### Panels

| Panel | Content | Stats |
|-------|---------|-------|
| 77a | Overview — stacked bar of ctrl vs mut aging (UP/DOWN/NS counts) | — |
| 77b | Overlap Venn: ctrl aging-UP peaks vs mut aging-UP peaks (via GenomicRanges::findOverlaps) | Fisher 2×2 |
| 77c | Mut-specific aging genes — GO enrichment dotplot (top 20); test for neuronal term enrichment | enrichGO (ORA) |
| 77d | Fold comparison at shared peaks — scatter (ctrl aging fold vs mut aging fold) with diagonal reference | Wilcoxon signed-rank on (mut_fold - ctrl_fold) |
| 77e | Patchwork composite: `(77a | 77b) / (77c | 77d)` | — |

### Key Logic
- Aging-UP: FDR<0.05 & Fold>0 in each contrast
- Shared peaks: findOverlaps with 1bp minimum (standard for DiffBind peaks)
- Mut-specific: aging-UP in mut but NOT overlapping any ctrl aging-UP peak
- 77d tests Jai's core question: at loci where BOTH genotypes gain MeCP2 with age, does the mutant gain MORE?

### Saved Tables
- `77_aging_peak_summary.tsv` — UP/DOWN/NS counts per genotype
- `77_aging_overlap.tsv` — shared/ctrl-only/mut-only counts + gene lists
- `77_mut_specific_aging_go.tsv` — GO enrichment results
- `77_shared_peak_fold_comparison.tsv` — paired fold values at shared peaks

---

## Execution Order

```
Section 72 (exists — provides neuronal gene set)
    ↓
Section 74 ──── can run immediately
Section 75 ──── can run immediately
Section 76 ──── can run immediately (also derives synapse/axon set)
Section 77 ──── BLOCKED until Jai provides young-vs-adult DiffBind files
```

Sections 74-76 are independent of each other and can run in any order.

## Verification

For each section:
```bash
cd /Users/zakiralibhai/Documents/GitHub/mariner_hi-c/biomodal/downstream
Rscript scripts/viz_sections/section_7X_*.R 2>&1 | tee logs/a2/7X_*.txt
```

Check:
1. All panels render without error
2. Saved tables have expected row counts
3. Statistical tests produce finite p-values
4. Composite plot assembles correctly
5. Multi-format output (PNG/PDF/SVG/JPG) in per-panel subfolders
