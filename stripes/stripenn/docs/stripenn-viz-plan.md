# Stripenn Differential Stripe Visualization & Results Documentation Plan

## Context

The Stripenn differential chromatin stripe pipeline (7 stages) has completed successfully for both timepoints (250831=early/P12, 250402=late/adult) at both resolutions (5kb, 10kb). Jesse Dixon suggested Stripenn as an alternative stripe caller during the April meeting — our Quagga pipeline found very few significant differential calls. Stripenn found dramatically more stripes (~3,000-7,000 per timepoint) but effect sizes remain small (all |logFC| < 0.39). The late timepoint (250402) has strong differential signal (31.5% significant at FDR<0.05) while the early (250831) has weak signal (2.4%).

We need to: (1) visualize results comprehensively, (2) create JuiceBox BEDPE files, (3) write results documents, following the Quagga `stripe_visualizations.R` pattern but adapted for Stripenn's larger scale and different column schema.

**User decisions:**
- **BEDPEs**: Create simple BEDPEs NOW from local cross_res data (for immediate JuiceBox use), then the R viz script produces full 28-column annotated versions later on HPC
- **ChIP marks**: Use ALL 5 marks (H3K27ac, H3K27me3, H3K4me1, H3K4me3, Bivalent) for BOTH timepoints
- **Results documents**: Minimal speculation — report actual results only

---

## Deliverables (8 files + generated outputs)

### 1. `stripes/stripenn/scripts/stripenn_visualizations.R` (~800-1000 lines)

Main visualization script. Processes both timepoints, both resolutions, generates all plots + BEDPE files + summary stats. Runs as Stage 7 of the pipeline.

**Section structure:**

| Section | Description | Plots Generated |
|---------|-------------|-----------------|
| SETUP | Load packages, config, define helpers, color constants | — |
| 1 | Load & validate data (all tp x res combinations) | — |
| 2 | Volcano plots (logFC vs -log10 FDR) | 4 per-tp/res + 1 combined 2x2 |
| 3 | Stripiness score analysis (NEW — Stripenn-specific) | 2 per tp: scatter + boxplot |
| 4 | Length & width distribution by direction | 3 per tp + 1 combined early/late |
| 5 | Source & direction distribution bars | 2 combined summary panels |
| 6 | Cross-resolution concordance | 2 per tp: logFC scatter + direction heatmap |
| 7 | Replicate correlation heatmap | 1 per tp/res (from count_matrix.tsv) |
| 8 | ChIP-seq anchor annotation | 2 per tp + annotated TSVs |
| 9 | BEDPE export (3 tiers per tp) | 6 BEDPE + 4 JuiceBox format files |
| 10 | GO/KEGG enrichment | Up to 8 dotplots (4 ontologies x 2 tp) |
| 11 | Summary statistics export | combined_summary.txt |

**Key adaptations from Quagga `stripe_visualizations.R`:**
- Column translation: `pos1/pos2` = anchor, `pos3/pos4` = span (not `anchor_x1/x2`, `span_y1/y2`)
- Width column: `stripe_width` (not `anchor_width`)
- No `detection_confidence` — use `direction_confidence` for both slots
- All 3 ChIP marks available for both timepoints (improvement over Quagga which only had 2 per tp)
- Scale handling: `alpha=0.3, size=1` for scatter plots; skip `geom_jitter` for n>2000; use `trim=TRUE` for violins
- New Stripenn-specific plots: stripiness ctrl vs mut scatter, direction_type (3prime/5prime) distribution, count matrix heatmap

**Libraries required:**
```r
ggplot2, patchwork, RColorBrewer, GenomicRanges, IRanges, 
TxDb.Mmusculus.UCSC.mm10.knownGene, org.Mm.eg.db,
dplyr, tidyr, yaml, pheatmap, svglite, clusterProfiler, enrichplot
```

**Reused utilities:**
- `scripts/utils/multi_format_output.R` — `save_multiformat_ggplot()`, `save_multiformat_pheatmap()`
- Color scheme matches Quagga: lost=blues (#00008B/#4169E1/#87CEEB), gained=reds (#8B0000/#DC143C/#FFA07A)

**Path pattern (same as all Stripenn scripts):**
```r
CODE_DIR <- "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/stripes/stripenn"
DATA_DIR <- "/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn"
PEAKS_DIR <- "/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/peaks/beds"
```

**Data loading:** From `05_final_differential.tsv` (5kb primary, 30 cols) and `cross_res_merged.tsv` (35 cols) per timepoint.

**BEDPE generation (Section 9) — 3 tiers per timepoint:**

| Tier | Filter | Expected counts (250402 / 250831) | Filename |
|------|--------|-----------------------------------|----------|
| High-confidence | FDR<0.05 AND direction_confidence=="high" AND lost/gained | ~1,005 / ~22 | `{tp}_stripes_highconf.bedpe` |
| All significant | FDR<0.05 AND lost/gained | ~2,320 / ~96 | `{tp}_stripes_allsig.bedpe` |  
| Cross-res concordant | resolution_support=="both_concordant" AND lost/gained | ~133 / ~8 | `{tp}_stripes_concordant.bedpe` |

Plus JuiceBox diagonal + rectangle BEDPE from the highconf tier:
- `{tp}_stripes_diagonal.bedpe` — both x,y span full extent (stripe along diagonal)
- `{tp}_stripes_rectangle.bedpe` — narrow anchor x, full extent y (off-diagonal body)

**BEDPE column schema (28 columns, matching Quagga format):**
```
chr1 x1 x2 chr2 y1 y2 name score strand1 strand2 color
direction direction_confidence logFC FDR source
pval_ctrl pval_mut detection_confidence in_10kb
nearest_gene distance_to_tss anchor_type
h3k27ac h3k27me3 h3k4me1
stripe_length_kb anchor_width_kb
```
Column mapping: x1=pos1, x2=pos2, y1=pos3, y2=pos4, detection_confidence=direction_confidence, anchor_width_kb=stripe_width/1000. ChIP columns populated by Section 8 annotation.

---

### 2. `stripes/stripenn/scripts/stripenn_visualizations.sb` (~50 lines)

SLURM wrapper. Adapted from `stripes/quagga/scripts/stripe_visualizations.sb`.

```
#SBATCH --cpus-per-task=8   (up from Quagga's 4 — larger datasets)
#SBATCH --mem=64G           (up from 32G — GO enrichment on thousands of genes)
#SBATCH --time=04:00:00     (up from 2h)
#SBATCH --account=csd940
#SBATCH --partition=shared
```

Checks for `cross_res_merged.tsv` existence for both timepoints before running. Uses `mariner_env` conda environment.

---

### 3. Config update: `stripes/stripenn/config/stripenn_config.yaml`

Add `chipseq_peaks` section with paths to all available marks. Files confirmed at `peaks/beds/`:

```yaml
chipseq_peaks:
  early:
    h3k27ac:   "peaks/beds/H3K27acCerebellumEarly2.bed"
    h3k27me3:  "peaks/beds/H3K27me3CerebellumEarly1.bed"
    h3k4me1:   "peaks/beds/H3K4me1CerebellumEarly1.bed"
    h3k4me3:   "peaks/beds/H3K4me3CerebellumEarly2.bed"
    bivalent:  "peaks/beds/Bivalent_Cerebellum_Early.bed"
  late:
    h3k27ac:   "peaks/beds/H3K27acCerebellumLate2.bed"
    h3k27me3:  "peaks/beds/H3K27me3CerebellumLate1.bed"
    h3k4me1:   "peaks/beds/H3K4me1CerebellumLate1.bed"
    h3k4me3:   "peaks/beds/H3K4me3CerebellumLate2.bed"
    bivalent:  "peaks/beds/Bivalent_Cerebellum_Late.bed"
```

All 5 marks for both timepoints (confirmed files exist in `peaks/beds/`). Paths are relative to the mariner_hi-c repo root. Resolved in script as `file.path(CODE_DIR, "../..", peak_path)`.

**Anchor classification (6 categories using all 5 marks):**
1. Bivalent_Promoter: Bivalent peak overlap AND <=2kb from TSS
2. Active_Promoter: H3K4me3+ AND H3K27ac+ AND <=2kb from TSS
3. Repressed_Promoter: H3K27me3+ AND <=2kb from TSS (no H3K27ac)
4. Active_Enhancer: H3K27ac+ AND >2kb from TSS
5. Poised_Enhancer: H3K4me1+ AND NOT H3K27ac AND >2kb from TSS
6. Polycomb: H3K27me3+ AND >2kb from TSS
7. Other: no marks

---

### 4. Results documents (written after viz script runs)

- `stripes/stripenn/outputs/250402/250402_results.md` — Late timepoint analysis
- `stripes/stripenn/outputs/250831/250831_results.md` — Early timepoint analysis

**Template following Quagga's `early-stripes.md` / `late-stripes.md` format:**

1. **Summary Statistics** — Tables for overall counts, direction breakdown, confidence tiers, effect sizes, cross-resolution concordance, anchor type distribution, ChIP-seq breakdown
2. **High-Priority Stripes** — Top 8-10 stripes for JuiceBox validation, prioritized by: resolution_support=="both_concordant" > direction_confidence=="high" > FDR > anchor_type relevance > logFC magnitude
3. **JuiceBox Coordinate Blocks** — Copy-paste `chr:start-end` with 50kb padding, grouped by lost/gained
4. **Biological Interpretation** — Minimal speculation per user request. Report what the data shows: direction bias, anchor type enrichment, ChIP mark overlaps, cross-timepoint contrast
5. **Caveats** — Effect size limitations, directional consistency rates, 250831 weak signal caveat

---

### 5. `stripes/stripenn/docs/stripenn-analysis-prompt.md`

Updated analysis prompt template with Stripenn-specific column definitions, prioritization rubric (cross-res > confidence > FDR > biology), and context notes on stripiness scores and direction_type.

---

## Output Directory Structure

After running the visualization script:

```
outputs/
├── 250402/                            # Late timepoint
│   ├── visualizations/                # All viz outputs for this tp
│   │   ├── volcano_250402_5kb.{pdf,svg,jpg}
│   │   ├── volcano_250402_10kb.{pdf,svg,jpg}
│   │   ├── stripiness_250402.{pdf,svg,jpg}
│   │   ├── length_distribution_250402.{pdf,svg,jpg}
│   │   ├── direction_breakdown_250402.{pdf,svg,jpg}
│   │   ├── cross_res_250402.{pdf,svg,jpg}
│   │   ├── replicate_correlation_250402_5kb.{pdf,svg,jpg}
│   │   ├── anchor_annotation_250402.{pdf,svg,jpg}
│   │   ├── 250402_annotated_stripes.tsv
│   │   ├── 250402_stripes_highconf.bedpe     # Tier 1
│   │   ├── 250402_stripes_allsig.bedpe       # Tier 2
│   │   ├── 250402_stripes_concordant.bedpe   # Tier 3
│   │   ├── 250402_stripes_diagonal.bedpe     # JuiceBox diagonal
│   │   └── 250402_stripes_rectangle.bedpe    # JuiceBox rectangle
│   └── 250402_results.md
├── 250831/                            # Early timepoint (same structure)
│   └── ...
└── combined/                          # Cross-timepoint comparison
    ├── volcano_combined.{pdf,svg,jpg}
    ├── length_comparison.{pdf,svg,jpg}
    ├── source_direction_summary.{pdf,svg,jpg}
    ├── cross_res_comparison.{pdf,svg,jpg}
    ├── comparison_summary.{pdf,svg,jpg}
    ├── enrichment/                    # GO/KEGG per tp
    │   ├── go_bp_dotplot_{tp}.{pdf,svg,jpg}
    │   ├── go_cc_dotplot_{tp}.{pdf,svg,jpg}
    │   ├── go_mf_dotplot_{tp}.{pdf,svg,jpg}
    │   └── kegg_dotplot_{tp}.{pdf,svg,jpg}
    └── combined_summary.txt
```

---

## Implementation Order

### Phase A: Immediate (local, no HPC needed)
1. **Add ChIP peak paths to `stripenn_config.yaml`** (small edit — add all 5 marks for both tp)
2. **Generate simple BEDPEs now** from `cross_res_merged.tsv` (locally available) — 3 tiers per timepoint, 15-column format sufficient for JuiceBox (chr1,x1,x2,chr2,y1,y2,name,score,strand1,strand2,color,direction,logFC,FDR,source). These go into `stripes/stripenn/outputs/{tp}/` directly.
3. **Write results documents** (`250402_results.md`, `250831_results.md`) — filled with actual numbers from pipeline output files (all synced locally). ChIP annotation section marked as "pending visualization script" since that requires HPC.
4. **Write `stripenn-analysis-prompt.md`** (analysis prompt template for future use)

### Phase B: Scripts (to run on HPC)
5. **Write `stripenn_visualizations.R`** (main script, ~800-1000 lines)
6. **Write `stripenn_visualizations.sb`** (SLURM wrapper, ~50 lines)

### Phase C: After HPC run
7. **Update results documents** with ChIP annotation findings, enrichment results
8. **Full 28-column annotated BEDPEs** replace the simple ones

---

## Key Data Points for Results Documents (from pipeline outputs)

### 250402 (Late/Adult) — Strong signal

| Metric | 5kb | 10kb | Cross-res |
|--------|-----|------|-----------|
| Total union stripes | 7,371 | 3,566 | 8,415 |
| Significant FDR<0.05 | 2,320 (31.5%) | 1,325 (37.2%) | — |
| Lost | 1,528 | 766 | 1,757 |
| Gained | 2,052 | 967 | 2,435 |
| Unchanged | 3,784 | 1,833 | 4,216 |
| Strengthened/Weakened | 3 / 4 | 0 / 0 | 3 / 4 |
| High-confidence lost/gained | 367 / 638 | 202 / 353 | — |
| Directional consistency | 60.7% | 61.5% | — |
| Common BCV | 0.012 | 0.015 | — |
| Cross-res matched | — | — | 2,522 |
| Both significant concordant | — | — | 377/672 |
| logFC correlation (5kb vs 10kb) | — | — | r=0.850 |

### 250831 (Early/P12) — Weak signal

| Metric | 5kb | 10kb | Cross-res |
|--------|-----|------|-----------|
| Total union stripes | 4,008 | 1,879 | 4,483 |
| Significant FDR<0.05 | 96 (2.4%) | 53 (2.8%) | — |
| Lost | 949 | 483 | 1,127 |
| Gained | 776 | 401 | 887 |
| Unchanged | 2,283 | 995 | 2,469 |
| High-confidence lost/gained | 12 / 10 | 14 / 10 | — |
| Directional consistency | 40.9% | 39.9% | — |
| Common BCV | 0.020 | 0.021 | — |
| Cross-res matched | — | — | 1,404 |
| Both significant concordant | — | — | 15/22 |
| logFC correlation (5kb vs 10kb) | — | — | r=0.808 |

---

## Verification

After implementation, verify:

1. **BEDPE files load in JuiceBox** — Column count matches 28, coordinates are valid genomic positions, colors render correctly
2. **Plot file sizes** — PDFs should be 50KB-5MB (not empty, not bloated)
3. **Row counts in BEDPE** — Match expected filter counts from the data tables above
4. **ChIP annotation** — `anchor_type` column is populated, `nearest_gene` is not all NA
5. **Enrichment** — At least 250402 should produce GO BP results given thousands of significant stripes
6. **Results documents** — All {N} placeholders filled with actual numbers, JuiceBox coordinates are valid
7. **No hardcoded local paths** — Script uses CODE_DIR/DATA_DIR for HPC, not `/Users/zakiralibhai/...`

---

## Critical Files to Modify/Create

| Action | File |
|--------|------|
| EDIT | `stripes/stripenn/config/stripenn_config.yaml` (add chipseq_peaks section) |
| CREATE | `stripes/stripenn/scripts/stripenn_visualizations.R` |
| CREATE | `stripes/stripenn/scripts/stripenn_visualizations.sb` |
| CREATE | `stripes/stripenn/docs/stripenn-analysis-prompt.md` |
| CREATE | `stripes/stripenn/outputs/250402/250402_results.md` |
| CREATE | `stripes/stripenn/outputs/250831/250831_results.md` |

## Reused Existing Code

| Source File | What's Reused |
|-------------|---------------|
| `stripes/quagga/scripts/stripe_visualizations.R` | Volcano plot function, length distribution, ChIP annotation logic, BEDPE construction, GO/KEGG enrichment, color scheme |
| `scripts/utils/multi_format_output.R` | `save_multiformat_ggplot()`, `save_multiformat_pheatmap()` |
| `stripes/stripenn/scripts/05_integration.R` | BEDPE color assignment logic, direction classification |
| `stripes/stripenn/scripts/06_compare_resolutions.R` | Cross-resolution matching pattern |
