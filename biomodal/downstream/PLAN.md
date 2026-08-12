# Plan: Consolidated Publication Figure Scripts for BAP1 Methylation Paper

## Context

The BAP1 methylation analysis has 78 exploratory sections producing 255 pre-computed TSV tables and 256 figure directories. The paper needs ~7 consolidated figures where each panel advances the narrative. The bottleneck is not additional computation — it's consolidation. This plan creates figure-level R scripts that read from existing section TSVs and produce unified multi-panel patchwork figures with consistent aesthetics, publication-quality formatting, and shared color palettes.

The figure order follows a narrative logic (methylation phenotype → mechanism → MeCP2 → neuronal specificity), which differs from the paper's R-section numbering but tells a more coherent story.

**Priority:** The repo `PLAN.md` (sections 74-77 analysis plan) is the authoritative reference for panel definitions, table names, and analytical decisions for those sections. Where `TODOS.md` (a discussion/brainstorm document) suggests something different, `PLAN.md` takes precedence. Figures that draw from sections 74-78 use the exact table names and panel logic specified in `PLAN.md`.

**Framing decisions resolved:**
- R3 sex-difference PCA: **reinstated** as Figure 1 panel G. Section 42 PCA exists (`42_pca_all_genes`, `42_pca_max_significance_genes`) and clearly shows PC1=condition (42.6-48.8% variance), PC2=sex (11.1-11.8% variance) for both 5mC and 5hmC. Effect-size density (panel E) also kept — both panels earn their space.
- R1 "binding active promoters mostly": **clarified** — methylation hyper lands at active promoters; MeCP2 gain lands at distal-intergenic sites. Different assays, different answers, coherent story.
- R5 "exclusively": **softened** to "preferentially/disproportionately" (OR=5.16 for methylation vs 1.73 for neuronal identity).
- R1 CUT&RUN panels: **included** using available data. Adult MeCP2 volcano from `MeCP2_annotated.txt` (202,651 peaks), aging trajectory from the two aging DiffBind files (both 51MB, arrived June 23). Only missing piece: P12 ctrl-vs-mut DiffBind file (would show "P12 = no DMRs" as in paper outline) — not currently in the repo. Per-sample normalized counts at P13/adult ARE available for PCA if needed.
- Section 77 (MeCP2 aging trajectory): **no longer blocked** — DiffBind files arrived June 23. All 4 tables populated per PLAN.md specification.

---

## New Files

### Directory Structure

```
downstream/
  scripts/
    figures/
      _figure_config.R                         # Publication settings
      figure1_methylation_phenotype.R           # R2/R3: What happened (6 panels, incl. Gviz)
      figure2_k119ub_chromatin_geography.R      # R2: Why it happened
      figure3_tet_mechanism.R                   # R2: Mechanism adjudication
      figure4_mecp2_redistribution.R            # R1: MeCP2 consequence + NEW aging overlap
      figure5_mecp2_chromatin_reader.R          # R4: Chromatin, not methylation
      figure6_neuronal_specificity.R            # R5: Neuronal genes disproportionately affected
      figure7_mechanism_summary.R              # Data panel for model schematic
      figureS1_supplemental.R                   # Supporting evidence
      run_all_figures.sh                        # Batch runner
  plots/
    figures/                                    # Output directory (auto-created)
      tables/                                   # New analysis results
```

### Critical Existing Files (reused, not modified)

- `scripts/viz_sections/_shared_config.R` — sources first; provides `mc_dmr`, `hmc_dmr`, `region_dmrs`, `COLORS`, `theme_biomodal()`, all path constants, all helpers
- `../../scripts/utils/multi_format_output.R` — auto-sourced by `_shared_config.R`; provides `save_multiformat_ggplot()`
- `plots/visualizations/tables/` — 255 pre-computed TSVs (the data source for all figure scripts)

---

## Shared Figure Utilities: `_figure_config.R`

Sourced by every figure script AFTER `_shared_config.R`. Provides:

1. **`theme_pub(base_size = 8)`** — extends `theme_biomodal()` for publication:
   - Helvetica/Arial font family
   - `panel.border = element_rect(colour = "black", linewidth = 0.5)`
   - `axis.ticks = element_line(linewidth = 0.3)`
   - `legend.key.size = unit(3, "mm")`
   - `plot.tag = element_text(face = "bold", size = 12)` (patchwork panel letters)

2. **`FIGURE_DIR`** — `file.path(BASE_DIR, "plots/figures")`, auto-created

3. **`save_figure(plot, name, width_mm, height_mm)`** — wraps `save_multiformat_ggplot()` with mm-to-inches conversion and 300 DPI

4. **`read_table_tsv(filename)`** — reads from `TABLES_DIR` with standard options

5. **`add_panel_labels(plot)`** — wraps `plot_annotation(tag_levels = "A")`

6. **`stat_label(name, value, fmt)`** — helper for annotating OR, AUC, p-values on panels

---

## Figure 1: The Coordinated Gene-Body Methylation Phenotype

**File:** `figure1_methylation_phenotype.R`
**Maps to:** R2 + R3
**Question:** What is the methylation phenotype in BAP1-KO?

| Panel | Sec | Type | Table/Data | Key Stat | Notes |
|-------|-----|------|-----------|----------|-------|
| A | 64 | Paired dot-line | `64_global_methylation_summary.tsv` | +0.32% mC, -0.39% hmC, total flat | 3 facets (mC/hmC/modC), lines connect matched samples, bar for means. Cohen's d annotated. |
| B | 03 | Horizontal bar | Pre-loaded `region_dmrs` | Gene body 51.4% >> Promoter 8.3% | 6 region bars showing % significant. Fill by region class. |
| C | 04 | Dual volcano | Pre-loaded `mc_dmr`, `hmc_dmr` | 10,775 mC / 11,484 hmC sig | Side-by-side. `geom_point(size=0.3, alpha=0.4)`, KEY_GENES labeled via `geom_text_repel`. Dashed FDR line. |
| D | 05 | Quadrant scatter | `coordinated_changes.tsv` | 78.7% coordinated (Q4=6,589) | mc_diff vs hmc_diff, colored by quadrant. Q4 count annotated. KEY_GENES labeled. |
| E | 07 | Density overlay | Pre-loaded `mc_dmr`, `hmc_dmr` (sig only) | Median mC +3.07% / hmC -2.14% | Mirror-image `geom_density` in mC red / hmC blue. Vertical dashed at medians. |
| F | 42 | PCA (5mC + 5hmC) | Section 42 per-sample methylation fractions (re-derive from `EXTRACT_PATHS` or read existing PCA coordinates) | PC1=condition (42.6%), PC2=sex (11.8%) | Side-by-side 5mC and 5hmC PCA. Points colored by condition, shaped by batch. Sample labels. Shows sex as secondary axis (R3 claim). |
| G | 46 | Genome browser (Gviz) | BigWig tracks via `METHYLATION_BIGWIGS` + `HISTONE_BIGWIGS` | Syt1: +17.3% mC, -15.0% hmC | Gviz `DataTrack` for mC ctrl/mut + hmC ctrl/mut + K119ub ctrl/mut at the Syt1 locus (chr10:108,747,000-108,812,000). Uses `Gviz::plotTracks()` saved via `save_multiformat_base()`. Gene annotation track from TxDb. |

**Layout:** `design = "AABB\nCCCC\nDDEE\nFFFF\nGGGG"`, `heights = c(1, 1.3, 1, 1, 1.4)`
**Dimensions:** 180 × 290 mm (tall — 7-panel figure with PCA + browser track; may split into Fig 1a/1b for journal)
**Computation:** Panels A-E: minimal (pre-loaded objects + TSVs). Panel F: PCA via `prcomp()` on per-sample gene-body methylation fractions from `EXTRACT_PATHS` (~20 lines). Panel G: Gviz rendering from BigWig files (~30s). Requires `library(Gviz)`.
**Note:** Panel G renders separately via `save_multiformat_base()` (Gviz uses base R graphics), then composed into the final patchwork via `patchwork::wrap_elements()` or assembled in the output PDF by stacking below the ggplot composite. If the figure is too tall for one page, panels A-F can be Figure 1 and panel G can be a supplemental genome browser panel.

---

## Figure 2: H2AK119ub Drives Hypermethylation at Active Chromatin

**File:** `figure2_k119ub_chromatin_geography.R`
**Maps to:** R2 (mechanism)
**Question:** Why does methylation change? K119ub is the upstream driver, targeting active chromatin.

| Panel | Sec | Type | Table/Data | Key Stat | Notes |
|-------|-----|------|-----------|----------|-------|
| A | 10 | Stacked fraction bar | `chromatin_state_summary.tsv` | Active_Promoter 93% hyper | `geom_col(position="fill")`. States on y-axis (CHROMATIN_STATE_ORDER). Fill by direction. |
| B | 33 | Forest plot | `diffbind_logistic_model_coefficients.tsv` | K119ub OR=4.71 | `geom_pointrange`: y=display_name, x=log2(or), xmin/xmax from or_lower/or_upper. Vertical line at OR=1. Color by mark. |
| C | 29 | OR bar chart | `compartment_fisher_tests.tsv` | A-compartment mC-hyper OR=13.64 | `geom_col` showing mC-hyper→A, hmC-loss→A, mC-hypo→B enrichments. Log-scale x. Annotate ORs. |
| D | 30 | Polycomb exclusion bar | `polycomb_fisher_tests.tsv` | Polycomb × hyper OR=0.063 | Mirror of panel C: exclusion from hyper, enrichment in hypo. Per-state hyper rates as inset or secondary axis. |
| E | 17 | Honest correction | `k119ub_honest_breakdown.tsv` | 58.2% no K119ub peak; +14.2pp gain | Small `geom_col`: met_group × pct_gained. Intellectual honesty panel. |

**Layout:** `design = "AABB\nCCDD\n#EE#"`, `heights = c(1.2, 1, 0.7)`
**Dimensions:** 180 × 200 mm
**Computation:** All from pre-computed TSVs. Bar fractions computed from count columns.

---

## Figure 3: TET-Impediment, Not DNMT3A Recruitment

**File:** `figure3_tet_mechanism.R`
**Maps to:** R2 (mechanistic adjudication)
**Question:** Is the methylation change caused by TET impediment or DNMT3A recruitment?

| Panel | Sec | Type | Table/Data | Key Stat | Notes |
|-------|-----|------|-----------|----------|-------|
| A | 23 | Dose-response bar | `demethylation_ratio_all_genes.tsv` | Baseline-5hmC AUC=0.762 | NEW computation: bin genes into deciles by WT 5hmC (mean_mod_group1 from hmc_dmr), plot median hmC delta per bin. Shows substrate-availability thesis as raw data. Annotate AUC from `baseline_hmc_model_comparison.tsv`. |
| B | 24 | Model comparison dot | `dnmt3a_model_comparison.tsv` | TET AUC=0.793 vs DNMT3A 0.696 | `geom_pointrange`: y=model, x=auc, xmin/xmax from CI columns. DeLong p=9.43e-49 annotated. |
| C | 22 | Chromatin-state lollipop | `demethylation_ratio_by_chromatin_state.tsv` | Active_Promoter strongest impairment (Δ-0.028) | `geom_segment` + `geom_point` of median_delta per state. Color by CHROMATIN_STATE_COLORS. P-value annotations from p_adj column. |
| D | 78 | Stoichiometry scatter | `78_stoichiometry_slopes.tsv` + `coordinated_changes.tsv` | Genome slope=-0.959; Neuronal=-1.0 | `geom_point(size=0.3)` of mc_diff vs hmc_diff, two overlaid regression lines (all genes, neuronal). Annotate slopes + whether ≠ -1. |
| E | 26 | Attenuation bar | `tet_ko_comparison_summary.tsv` | 3.3% attenuation ("dimmer not switch") | `geom_col` comparing BAP1-KO vs TET-KO delta. Small panel. |

**Layout:** `design = "AABB\nCCDD\n#EE#"`, `heights = c(1, 1, 0.6)`
**Dimensions:** 180 × 200 mm
**Computation:** Panel A requires new binning (~15 lines): read `demethylation_ratio_all_genes.tsv`, add WT 5hmC from `hmc_dmr$mean_mod_group1`, bin into deciles, compute median hmC change per bin. Panel D overlays raw scatter with regression lines from `78_stoichiometry_slopes.tsv`.

---

## Figure 4: MeCP2 Redistribution and Developmental Amplification

**File:** `figure4_mecp2_redistribution.R`
**Maps to:** R1 (the title figure)
**Question:** How does MeCP2 binding change? It redistributes to distal sites and the developmental ramp is amplified.

| Panel | Sec | Type | Table/Data | Key Stat | Notes |
|-------|-----|------|-----------|----------|-------|
| A | 11 | MeCP2 adult volcano | `MeCP2_annotated.txt` (202,651 peaks via `MECP2_FILES$annotated`) | 7,686 UP / 1,200 DOWN peaks | `geom_point(size=0.3, alpha=0.4)` of Fold vs -log10(FDR). Color by direction (UP orange / DOWN purple from `COLORS$mecp2`). Top genes labeled. Dashed FDR line. |
| B | 75 | Peak annotation bar | `75_peak_annotation_distribution.tsv` | UP: 52% distal-intergenic, 2.2% promoter | `geom_col(position="fill")` for UP vs DOWN peaks by annotation category. Shows MeCP2 gains at distal sites, not promoters. |
| C | 11 | DMR overlap enrichment | `mecp2_dmr_overlap_summary.tsv` | OR=5.13, p=1.27e-33 | 2×2 bar: Hyper/Hypo × MeCP2-Up/Down. Annotate Fisher OR. |
| D | 77 | Aging peak count bar | `77_aging_peak_summary.tsv` | Mut 2.1× more aging-UP | Grouped `geom_col`: Control vs Mutant × UP/DOWN. Annotate fold difference. |
| E | 77 | Shared-peak fold scatter | `77_shared_peak_fold_comparison.tsv` | +22% fold at shared loci | `geom_point(size=0.3, alpha=0.3)` of ctrl_fold vs mut_fold. `geom_abline` at 1:1. Annotate medians. Color by neuronal status. |
| F | NEW | Aging × methylation overlap | Re-derived from MeCP2 aging DiffBind + `mecp2_coordinated_genes.tsv` | Fisher OR for mut-unique aging genes × coordinated DMRs | Small panel: observed vs expected overlap + Fisher annotation. See New Analysis below. |

**Layout:** `design = "AABB\nCCDD\nEEFF"`, `heights = c(1.2, 1, 0.8)`
**Dimensions:** 180 × 220 mm
**Computation:** Panel A reads `MeCP2_annotated.txt` (202k rows) via `load_diffbind_flex()`. Panel F is the new aging × methylation Fisher test (~50 lines). Other panels from pre-computed TSVs.

**Note on P12 MeCP2 data:** The paper outline calls for a P12 MeCP2 volcano showing "no DMRs." A P12 ctrl-vs-mut DiffBind result file is not currently in the repo — only the aging trajectory files (young-vs-adult within each genotype). Per-sample normalized counts at P13 ARE available in `MeCP2_ctrl_adultvsyoung_normalized_counts.csv` and could be used for a PCA panel if needed, but a proper P12 DiffBind contrast requires either Jai's file or re-running DiffBind on the P12 samples.

---

## Figure 5: MeCP2 Reads Chromatin State, Not Methylation

**File:** `figure5_mecp2_chromatin_reader.R`
**Maps to:** R4
**Question:** Does MeCP2 redistribution track DNA methylation or chromatin state?

| Panel | Sec | Type | Table/Data | Key Stat | Notes |
|-------|-----|------|-----------|----------|-------|
| A | 62 | R² comparison bar | `62_model_comparison_summary.tsv` | CG R²=0.017 vs Chromatin R²=0.246 | `geom_col` grouped by analysis level (Binding/Differential), fill by model type (CG/Chromatin/Full). |
| B | 62 | Variance partition | `62_variance_partition.tsv` | CG unique 1.5%, Chromatin unique 24.3% | `geom_col(position="stack")`: 4 components (CG unique, Shared, Chromatin unique, Unexplained). |
| C | 67 | K119ub at MeCP2-no-meth genes | `67_statistics.tsv` + `60_mecp2_status_stats.tsv` | 359 genes, K119ub OR=3.15 | `geom_violin` or `geom_boxplot` of K119ub log2FC partitioned by MeCP2 status (Up/Down/NS). P-value annotated. |
| D | 60 | Mirror-image profiles | `60_mecp2_status_stats.tsv` | MeCP2-Up: K119ub↑; Down: K119ub flat | `geom_point` + `geom_linerange` per mark. Marks on x-axis, mean change on y. Three groups (Up/Down/NS) as facets or colors. The flat K119ub at MeCP2-Down is the visual punchline. |
| E | 71 | K119ub vs ratio variance | `71_variance_partition.tsv` | K119ub 7.3% vs ratio 0.0% | `geom_col(position="stack")`: 4 components. Drives home K119ub dominance. |

**Layout:** `design = "AABB\nCCDD\nEEEE"`, `heights = c(1, 1, 0.8)`
**Dimensions:** 180 × 200 mm
**Computation:** All from pre-computed TSVs. Panel D requires reshaping `60_mecp2_status_stats.tsv` from wide to long for the multi-mark comparison.

---

## Figure 6: Neuronal Genes Are Preferentially Affected

**File:** `figure6_neuronal_specificity.R`
**Maps to:** R5
**Question:** Why are neuronal genes disproportionately affected? Constitutive K119ub enrichment + Polycomb de-repression.

| Panel | Sec | Type | Table/Data | Key Stat | Notes |
|-------|-----|------|-----------|----------|-------|
| A | 72 | Dose-response line | `72_neuronal_decile_summary.tsv` | Top-decile OR=1.70 | `geom_line` + `geom_point`: neuronal fraction per K119ub decile. Annotate OR at top quartile/decile thresholds. |
| B | 74 | OR comparison bar | `74_pairwise_fisher.tsv` | Coord×MeCP2-Up OR=5.16 >> Neuro×MeCP2-Up OR=1.73 | `geom_col` with 3 pairwise Fisher ORs. Log-scale y. P-values annotated. Key message: MeCP2 redistribution tracks methylation 5× more than lineage. |
| C | 76 | Synapse K27me3 violin | `76_chromatin_stats.tsv` | Synapse K27me3 Δ=-0.044, p=2.95e-3 | `geom_boxplot` faceted by mark (K27me3/K27ac/ATAC), split by gene class (synapse/broader-neuronal/non-neuronal). Synapse-specific K27me3 loss is the punchline. |
| D | 78 | Stoichiometry by gene class | `78_stoichiometry_slopes.tsv` | Neuronal slope=-1.0 (stoichiometric) | `geom_col` with error bars (CI): slopes for all/non-neuronal/neuronal/synapse. Dashed line at -1. Annotate which groups differ from -1. |

**Layout:** `design = "AABB\nCCDD"`, `heights = c(1, 1)`
**Dimensions:** 180 × 160 mm
**Computation:** All from pre-computed TSVs. Panel C requires reshaping `76_chromatin_stats.tsv` for the faceted comparison.

---

## Figure 7: Mechanism Summary Data Panel

**File:** `figure7_mechanism_summary.R`
**Maps to:** All R-sections (cascade summary)
**Question:** What is the quantitative evidence for each step of the BAP1 → K119ub → TET-block → methylation → MeCP2 cascade?

This is the data companion to the Illustrator model schematic. Each panel quantifies one link in the chain.

| Panel | Step | Type | Table/Data | Key Stat | Notes |
|-------|------|------|-----------|----------|-------|
| A | BAP1→K119ub | Effect-size bar | `k119ub_bigwig_signal_summary.tsv` | Global K119ub increase (median log2FC +0.007, p=1.8e-20) | `geom_col` showing background vs mC-Up vs mC-Down K119ub signal. |
| B | K119ub→mC | OR cascade waterfall | `diffbind_logistic_model_coefficients.tsv` | K119ub OR=4.71 (dominant) | Horizontal `geom_col` with all 4 marks, ordered by OR magnitude. Reuses Figure 2B data but in waterfall/cascade style. |
| C | TET-block evidence | AUC comparison | `dnmt3a_model_comparison.tsv` + `baseline_hmc_model_comparison.tsv` | TET 0.793 >> DNMT3A 0.696 | `geom_col` with model AUCs side by side. Annotate DeLong p. |
| D | Methylation→MeCP2 | R² cascade | `62_model_comparison_summary.tsv` | Chromatin R²=0.246 >> CG R²=0.017 | `geom_col` showing R² at each step: CG alone, Chromatin alone, K119ub unique variance. |
| E | Neuronal convergence | OR summary | `74_pairwise_fisher.tsv` + `72_fisher_results.tsv` | Coordinated OR=5.16; K119ub top-decile OR=1.70 | `geom_col` with key ORs from each analysis layer. Summary of how all pathways converge on neuronal genes. |

**Layout:** Single-column cascade: `design = "A\nB\nC\nD\nE"`, `heights = c(0.6, 0.8, 0.6, 0.8, 0.6)` — or 2-column: `design = "AABB\nCCDD\n#EE#"`
**Dimensions:** 85 × 220 mm (half-width, tall — designed to sit next to the Illustrator model schematic) or 180 × 140 mm (full-width, short)
**Computation:** All from pre-computed TSVs. Minimal reshaping.

---

## Figure S1: Supplemental Evidence

**File:** `figureS1_supplemental.R`
**Panels supporting the main story but too detailed for main figures:**

| Panel | Sec | Content | Table |
|-------|-----|---------|-------|
| A | 37 | Permutation validation (15/15 confirmed) | `permutation_37_summary.tsv` |
| B | 45 | Field chr8 cross-species concordance | `field_chr8_comparison_full.tsv` + `field_chr8_statistical_tests.tsv` |
| C | 47 | CTCF anchor methylation (lost→hyper OR=3.28) | `47a_dynamic_ctcf_anchor_methylation.tsv` |
| D | 48 | CpG island K119ub depletion | `48_cpg_island_ubiquitination_summary.tsv` |
| E | 44 | ASM doubling (1.95× in mutant) | `asm_dmr_overlap_summary.tsv` |
| F | 78 | Narrow vs broad neuronal self-correction | `78_total_methylation_summary.tsv` |

**Layout:** 3×2 grid, `design = "AABB\nCCDD\nEEFF"`
**Dimensions:** 180 × 240 mm

---

## New Analysis: MeCP2 Aging × Methylation Overlap

**Location:** Inside `figure4_mecp2_redistribution.R` (panel E), ~50 lines.

**What it computes:**
1. Re-derive mut-unique aging genes from MeCP2 DiffBind files (Section 77 saves summary counts but NOT the gene list itself):
   - Read `MECP2_FILES$ctrl_aging` and `MECP2_FILES$mut_aging` via `load_diffbind_flex()`
   - Filter to FDR < 0.05 & Fold > 0 (aging-UP)
   - Annotate with ChIPseeker → extract SYMBOL
   - `setdiff(mut_up_genes, ctrl_up_genes)` → ~1,654 mut-unique genes
2. Load coordinated genes from `mecp2_coordinated_genes.tsv`
3. Fisher 2×2: (mut-unique? y/n) × (coordinated? y/n) over the ~20,915-gene universe
4. Save to `plots/figures/tables/fig4e_aging_methylation_overlap.tsv`
5. Plot: `geom_col` of observed vs expected overlap + Fisher OR annotation

**Dependencies:** ChIPseeker and TxDb already loaded by `_shared_config.R`. MeCP2 aging DiffBind files must be present (confirmed: both arrived June 23).

**Fallback:** If DiffBind files are missing, skip panel E with `tryCatch` warning (not a hard `stopifnot` — the rest of Figure 4 should still render).

---

## Execution Plan

### Batch Runner: `run_all_figures.sh`

Modeled on existing `run_sections_74_78.sh`:
```bash
#!/bin/bash
SCRIPT_DIR="scripts/figures"
LOG_DIR="logs/figures"
mkdir -p "$LOG_DIR"
# Sequential execution: figure1 through figureS1
# Each: Rscript $SCRIPT_DIR/figureN.R 2>&1 | tee $LOG_DIR/figureN.txt
```

### Execution Order

All scripts are independent (each sources `_shared_config.R`). For first run, execute sequentially to catch errors early:

1. `figure1` — exercises pre-loaded `mc_dmr`/`hmc_dmr` + basic TSVs + Gviz (BigWig I/O)
2. `figure2` — exercises logistic model + compartment tables
3. `figure3` — exercises model comparison + new dose-response computation
4. `figure4` — exercises MeCP2 aging files + NEW overlap analysis (most likely error source)
5. `figure5` — exercises regression/variance tables
6. `figure6` — exercises neuronal characterization tables
7. `figure7` — exercises cross-figure summary tables (lightweight)
8. `figureS1` — catchall supplemental

**Estimated runtime:** ~7 minutes total (20s config loading + 10-30s plotting per script; figure1 ~60s due to Gviz BigWig rendering).

---

## Verification

```bash
cd downstream/

# 1. All outputs exist
ls plots/figures/figure*/

# 2. No hardcoded run paths
grep -rn 'run-[0-9]' scripts/figures/*.R  # should return 0

# 3. No "ChIP" references (MeCP2 is CUT&RUN)
grep -i "chip" scripts/figures/*.R  # should return 0

# 4. New analysis table saved
cat plots/figures/tables/fig4e_aging_methylation_overlap.tsv

# 5. Panel labels present (spot-check PDFs)
# 6. Color consistency (5mC always #E41A1C, 5hmC always #377EB8)
# 7. Key stats match source TSVs (K119ub OR=4.71, coordinated 78.7%, CG R²=0.017)
```

---

## Implementation Order

Write scripts in this order:

1. **`_figure_config.R`** — shared utilities, ~60 lines
2. **`figure1_methylation_phenotype.R`** — validates template pattern; Gviz panel exercises BigWig paths
3. **`figure5_mecp2_chromatin_reader.R`** — strongest result (R4), cleanest tables
4. **`figure2_k119ub_chromatin_geography.R`** — forest plot + compartment
5. **`figure3_tet_mechanism.R`** — model comparison + new dose-response binning
6. **`figure6_neuronal_specificity.R`** — neuronal characterization
7. **`figure4_mecp2_redistribution.R`** — includes new aging analysis (do last, most risk)
8. **`figure7_mechanism_summary.R`** — cascade summary data panel
9. **`figureS1_supplemental.R`** — supplemental panels
10. **`run_all_figures.sh`** — batch runner

Total: 10 files. Each figure script ~150-300 lines (figure1 longest due to Gviz). `_figure_config.R` ~60 lines. `run_all_figures.sh` ~50 lines.

---

## Paper Figure ↔ R-Section Mapping Summary

| Figure | Paper R-Section | Primary Source Sections | TSV Count | Special |
|--------|----------------|----------------------|-----------|---------|
| Fig 1 | R2, R3 | 03, 04, 05, 07, 42, 46, 64 | ~6 + BigWigs | PCA + Gviz browser track |
| Fig 2 | R2 | 10, 17, 29, 30, 33 | ~5 | — |
| Fig 3 | R2 | 22, 23, 24, 26, 78 | ~6 | New dose-response binning |
| Fig 4 | R1 | 11, 75, 77, MeCP2 DiffBind, NEW | ~7 + MeCP2 peaks | Adult MeCP2 volcano + new aging Fisher |
| Fig 5 | R4 | 59, 60, 62, 67, 71 | ~7 | — |
| Fig 6 | R5 | 72, 74, 76, 78 | ~6 | — |
| Fig 7 | All (cascade) | 33, 24, 62, 72, 74 | ~5 | Data companion to model schematic |
| Fig S1 | Methods/Support | 37, 44, 45, 47, 48, 78 | ~7 | — |

---

## Session Phases

This plan is designed to be executed across **4 Claude sessions**, each producing 2-3 working figure scripts. Each session starts by reading this master plan, writing the targeted scripts, running them, and verifying output. Sessions are independent — if one runs long, the others aren't blocked.

### Session 1: Infrastructure + Figure 1 + Figure 2
**Scope:** Build the shared utilities and the two most foundational figures.
**Scripts:**
- `_figure_config.R` — shared publication settings (~60 lines; write first, validates the sourcing pattern)
- `figure1_methylation_phenotype.R` — 7 panels (A-G) including PCA and Gviz; the most complex script (~300 lines)
- `figure2_k119ub_chromatin_geography.R` — 5 panels; forest plot + compartment/Polycomb (~200 lines)

**Why grouped:** Figure 1 is the template — getting it right establishes the pattern all other scripts follow. Figure 2 is a natural pair (both R2). `_figure_config.R` must exist before anything else runs.

**Verify:** Both figures render, PCA shows condition/sex separation, Gviz Syt1 track renders, forest plot OR values match TSV, all 4 output formats generated.

### Session 2: Figure 3 + Figure 5
**Scope:** The mechanistic adjudication (TET vs DNMT3A) and the MeCP2-reads-chromatin verdict.
**Scripts:**
- `figure3_tet_mechanism.R` — 5 panels; includes new dose-response binning computation (~220 lines)
- `figure5_mecp2_chromatin_reader.R` — 5 panels; regression/variance tables (~200 lines)

**Why grouped:** Both are the quantitative analytical core of the paper — model comparisons, variance partitions, regression results. Similar plot types (bars, forest plots, scatter). Figure 5 is the cleanest data and good for validating the pattern after the heavier Figure 3.

**Verify:** Dose-response binning produces monotonic trend, ROC AUCs match TSV values, R² bars render correctly, mirror-image profiles show flat K119ub at MeCP2-Down.

### Session 3: Figure 4 + Figure 6
**Scope:** MeCP2 redistribution (the title figure) and neuronal specificity.
**Scripts:**
- `figure4_mecp2_redistribution.R` — 5 panels; includes the NEW aging × methylation overlap analysis (~250 lines)
- `figure6_neuronal_specificity.R` — 4 panels; neuronal characterization (~180 lines)

**Why grouped:** Both draw from sections 74-78 (the most recently completed analyses per PLAN.md). Figure 4 has the only genuinely new computation (aging × methylation Fisher test). Figure 6 uses the neuronal gene sets established in sections 72/74/76.

**Verify:** New aging overlap Fisher test produces finite OR/p, table saved to `plots/figures/tables/`, synapse K27me3 violin shows selective loss, stoichiometry slopes match `78_stoichiometry_slopes.tsv`.

### Session 4: Figure 7 + Figure S1 + Runner + Final Verification
**Scope:** Wrap up with the summary data panel, supplemental, batch runner, and end-to-end verification.
**Scripts:**
- `figure7_mechanism_summary.R` — 5 panels; cascade summary (~150 lines)
- `figureS1_supplemental.R` — 6 panels; permutation/Field/CTCF/ASM/CpGI/correction (~200 lines)
- `run_all_figures.sh` — batch runner with logging (~50 lines)

**Why last:** Figure 7 and S1 are the least complex and most forgiving. The batch runner needs all scripts to exist. This session ends with a full end-to-end run of all 8 figure scripts, verifying consistency across the complete set.

**Verify:** Full `run_all_figures.sh` completes without error, all `plots/figures/` directories populated, `grep -rn 'run-[0-9]' scripts/figures/*.R` returns 0, `grep -i 'chip' scripts/figures/*.R` returns 0, spot-check 5 key stats across figures match source TSVs.

### Session Overview

| Session | Scripts | Lines (est.) | New Computation | Risk |
|---------|---------|-------------|-----------------|------|
| 1 | `_figure_config.R`, Fig 1, Fig 2 | ~560 | PCA derivation, Gviz rendering | Medium (Gviz BigWig I/O) |
| 2 | Fig 3, Fig 5 | ~420 | Dose-response binning | Low (clean tables) |
| 3 | Fig 4, Fig 6 | ~430 | Aging × methylation Fisher | Medium (new analysis) |
| 4 | Fig 7, Fig S1, runner | ~400 | None | Low (wrap-up) |
| **Total** | **10 files** | **~1,810** | | |
