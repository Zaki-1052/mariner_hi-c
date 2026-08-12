# Section 46: Multi-Omic Genome Browser Locus Views — Continuation Plan

## Current State (2026-04-05)

### What Works

The script `scripts/viz_sections/section_46_genome_browser_loci.R` is functional and generates
multi-omic Gviz genome browser figures for key gene loci. It sources `_shared_config.R` and
follows the same conventions as all other section scripts.

**Successfully generated (full + compact) for 7/10 KEY_GENES:**
Syt1, Epha3, Mcu, Cntnap2, Lpp, Arhgap26, Cdh8

**Failed with "Too many stacks" (3 genes):**
Zbtb20, Trpm3, Dlgap1 — the GeneRegionTrack has too many overlapping transcripts in
large regions. Fix is already in place (`collapseTranscripts = "meta"` added to
GeneRegionTrack constructor around line 382) but needs re-run to verify.

**Output location:** `plots/visualizations/46_genome_browser/{Gene}_locus/` and
`{Gene}_locus_compact/` — each contains PDF, SVG, JPG.

### Tracks Included

**Full variant (~15 tracks):**

| # | Track | Source | Type |
|---|-------|--------|------|
| 0 | Genome axis | Gviz built-in | GenomeAxisTrack |
| 1 | Gene model | TxDb.Mmusculus.UCSC.mm10.knownGene | GeneRegionTrack |
| 2 | 5mC% Control | Mean of 4 ctrl BigWigs | DataTrack histogram |
| 3 | 5mC% Mutant | Mean of 4 mut BigWigs | DataTrack histogram |
| 4 | 5mC% Difference | Mut - Ctrl (x100 for %) | DataTrack histogram |
| 5 | 5hmC% Control | Mean of 4 ctrl BigWigs | DataTrack histogram |
| 6 | 5hmC% Mutant | Mean of 4 mut BigWigs | DataTrack histogram |
| 7 | 5hmC% Difference | Mut - Ctrl (x100 for %) | DataTrack histogram |
| 8 | H2AK119ub | Ctrl/Mut condition BEDs | AnnotationTrack (dense) |
| 9 | H3K27me3 | Late timepoint BED | AnnotationTrack (dense) |
| 10 | H3K27ac | Ctrl/Mut condition BEDs | AnnotationTrack (dense) |
| 11 | ATAC-seq | Differential up/down BEDs | AnnotationTrack (dense) |
| 12 | MeCP2 | Differential up/down BEDs | AnnotationTrack (dense) |
| 13 | CTCF | Peak BED | AnnotationTrack (dense) |
| 14 | CpG features | Islands/shores/shelves BEDs | AnnotationTrack (dense) |
| 15-16 | Hi-C loops | Lost/Gained InteractionTracks | InteractionTrack arcs |

**Compact variant (~8 tracks):** Tracks 0-4, 5 (ctrl only), 8, 13, 14, 15-16.

### Key Implementation Details

**BigWig averaging** (`average_bigwigs_in_region`, ~line 125):
- Imports each of 4 BigWigs for the target region using `rtracklayer::import.bw(path, which = region_gr)`
- Auto-selects bin size: `max(50, round(region_width / 1000))` (~1000 bins per view)
- Computes **mean methylation at CpG sites per bin** using `viewSums(weighted_cov) / viewSums(count_cov)`
  - NOT `viewMeans` — that divides by bin width, diluting signal with inter-CpG zeros
- BigWig values are 0-1 fractions; script multiplies by 100 for percentage display

**GenomicInteractions seqinfo fix** (line 23):
- `GenomicInteractions` class (v1.42.0) lacks a `seqinfo` method that Gviz requires
- Script defines it: `setMethod("seqinfo", "GenomicInteractions", function(x) seqinfo(regions(x)))`
- Anchors must have `genome(a1) <- "mm10"` set before creating GenomicInteractions

**Chromosome naming:**
- BigWigs, ChIP BEDs, ATAC BEDs, MeCP2 BEDs, CTCF BED, Hi-C loops: all use `chr` prefix
- CpG annotation BEDs: bare numbers — `load_cpg_annotation()` adds `chr` prefix
- DMR BEDs: bare numbers — `load_dmr_bed()` in `_shared_config.R` adds `chr` prefix
- The `import.bw` warning "levels in 'seqnames' with no entries in 'seqinfo' were dropped"
  is harmless — BigWig has all chromosomes, we query one, unused levels are dropped

**RNA-seq integration:**
- Loads DESeq2 results from `tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx`
- Column: `ensembl_gene_id` (gene symbol), `log2FoldChange`, `padj`
- Currently shown as text in the gene track title: `[RNA-seq: log2FC=X, padj=Y]`
- `RNASEQ_BW_CTRL` / `RNASEQ_BW_MUT` config vars (default NULL) are ready for
  BigWig coverage tracks when the user locates the files

---

## Cosmetic Fixes Needed

These are all in `plot_locus_browser()`:

1. **Track labels truncated on the left.** Increase `title.width` from 1.2 to ~1.6 in the
   `plotTracks()` call (around line 490).

2. **Y-axis labels too small.** Bump `fontsize` from 10 to 11 or 12, or set `cex.axis = 0.8`
   on individual DataTracks.

3. **Gene model transcript clutter.** `collapseTranscripts = "meta"` is now set (line ~382)
   but the 3 failed genes need re-running to verify the fix works.

4. **Difference track single-color.** Currently the difference track uses one color for all
   bars (hyper color). Ideally it should be bipolar: red above 0, blue below 0. Two approaches:
   - Split into two DataTracks (positive/negative) and use OverlayTrack
   - Or use `type = "horizon"` for built-in bipolar display

5. **`load_peak_as_granges` defensiveness.** If a BED file has a header row or NA entries,
   the function can produce `'seqnames' cannot contain NAs`. Add filtering of rows where
   `V2`/`V3` are non-numeric (header rows) before creating GRanges.

---

## Composite Multi-Panel Figure (Not Yet Started)

### Goal

A single publication-quality composite figure combining:
- **Panel A:** Syt1 genome browser view (compact variant)
- **Panel B:** Coordinated mC/hmC scatter plot (from section_05 data)
- **Panel C:** Cross-gene summary heatmap (KEY_GENES x multi-omic metrics)

### Panel A: Browser View

Use the compact Syt1 figure already generated. To embed in a composite with ggplot panels,
capture the Gviz output as a grob:

```r
library(gridGraphics)  # or use grid::grid.grab()

# Render Gviz to current device, then grab as grob
grid.newpage()
plotTracks(track_list, chromosome = chr, from = view_start, to = view_end, ...)
panel_a <- grid.grab()
```

Alternatively, `ggplotify::as.grob()` or `cowplot::as_grob()` can wrap base graphics.
If these don't work cleanly, save Panel A as a separate SVG and assemble in Illustrator.

### Panel B: Coordinated mC/hmC Scatter

Recompute or load from section_05:

```r
# Data source: plots/visualizations/tables/coordinated_changes.tsv
# Or recompute from mc_dmr/hmc_dmr (already loaded by _shared_config.R):
mc_sig <- mc_dmr %>% filter(significant) %>% select(gene, mc_diff = mod_difference, mc_q = dmr_qvalue)
hmc_sig <- hmc_dmr %>% filter(significant) %>% select(gene, hmc_diff = mod_difference, hmc_q = dmr_qvalue)
coordinated <- inner_join(mc_sig, hmc_sig, by = "gene") %>%
  mutate(coordinated = (mc_diff > 0 & hmc_diff < 0),
         combined_effect = abs(mc_diff) + abs(hmc_diff))

p_scatter <- ggplot(coordinated, aes(x = mc_diff * 100, y = hmc_diff * 100)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(aes(color = coordinated), alpha = 0.4, size = 1.5) +
  scale_color_manual(values = c("TRUE" = "#D7191C", "FALSE" = "grey60")) +
  # Highlight Syt1
  geom_point(data = . %>% filter(gene == "Syt1"), size = 3, shape = 21,
             fill = "gold", color = "black", stroke = 1) +
  geom_text_repel(data = . %>% filter(gene == "Syt1"), aes(label = gene),
                  fontface = "italic", size = 4) +
  labs(x = "5mC change (%)", y = "5hmC change (%)",
       title = "Coordinated mC/hmC Changes") +
  theme_biomodal()
```

### Panel C: Cross-Gene Summary Heatmap

Show KEY_GENES as rows, multi-omic metrics as columns:

```r
# Data sources to merge per gene:
# 1. mC difference: mc_dmr$mod_difference (from _shared_config.R)
# 2. hmC difference: hmc_dmr$mod_difference
# 3. Delta ratio: computed as hmc/(mc+hmc) change
# 4. RNA-seq log2FC: from rnaseq_df$log2FoldChange (loaded in section_46)
# 5. Loop status: from loops_df — count lost/gained loops at each gene's locus
# 6. K119ub status: overlap K119ub_mut peaks with gene body -> binary gain/loss
# 7. Chromatin state: from plots/visualizations/tables/dmr_chromatin_state_annotation.tsv

# Build matrix:
summary_df <- data.frame(gene = KEY_GENES) %>%
  left_join(mc_dmr %>% filter(gene %in% KEY_GENES) %>% select(gene, mC_diff = mod_difference)) %>%
  left_join(hmc_dmr %>% filter(gene %in% KEY_GENES) %>% select(gene, hmC_diff = mod_difference)) %>%
  left_join(rnaseq_df %>% select(gene = ensembl_gene_id, log2FC = log2FoldChange, padj))
  # ... add loop counts, K119ub overlap, chromatin state

# Render as pheatmap with column clustering off, row order by combined_effect
```

### Assembly

```r
library(patchwork)
library(cowplot)

# Option 1: patchwork (if Panel A can be converted to grob)
composite <- wrap_elements(panel_a) / (panel_b | panel_c) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(tag_levels = "A")

save_multiformat_ggplot(composite, file.path(SECTION_OUTPUT_DIR, "composite_syt1_panel"),
                        width = 14, height = 16)

# Option 2: cowplot::plot_grid (more flexible with mixed grob/ggplot)
composite <- plot_grid(
  panel_a_grob,
  plot_grid(panel_b, panel_c, nrow = 1, labels = c("B", "C")),
  nrow = 2, rel_heights = c(2, 1), labels = c("A", "")
)

# Option 3: Save panels separately, assemble in Illustrator
# (most reliable for mixed Gviz + ggplot)
```

### Output

`plots/visualizations/46_genome_browser/composite_syt1_panel/{pdf,svg,jpg}`

---

## Data Source Reference (for future agents)

All paths are relative to `downstream/` (the working directory).

| Data | Path | Format | Loaded By |
|------|------|--------|-----------|
| 5mC BigWigs (8) | `/Users/zakiralibhai/Documents/BIO_LAB/methylation-tracks/bigwigs/mc/*.bw` | BigWig, 0-1 fractions, 1bp CpG intervals | `average_bigwigs_in_region()` |
| 5hmC BigWigs (8) | `/Users/zakiralibhai/Documents/BIO_LAB/methylation-tracks/bigwigs/hmc/*.bw` | BigWig, same format | same |
| Gene body DMRs | via `_shared_config.R` -> `mc_dmr`, `hmc_dmr` | data.frame, 14 cols, chr prefix added | `load_dmr_bed()` |
| RNA-seq DESeq2 | `../../tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx` | Excel, cols: ensembl_gene_id, log2FoldChange, padj | `readxl::read_excel()` |
| Hi-C loops | `../../peaks/loop_annotation_extended/late/extended_characterized_loops.tsv` | TSV, 57 cols: chr1/start1/end1/chr2/start2/end2/logFC/FDR/direction | `read.table()` |
| H2AK119ub peaks | `../../peaks/intersect/P51_K119ub_{ctrl,mut}_intersect.bed` | BED 6-col, chr prefix | `load_peak_as_granges()` |
| H3K27ac peaks | `../../peaks/intersect/P60_K27ac_{ctrl,mut}_intersect.bed` | BED 6-col, chr prefix | same |
| H3K27me3 peaks | `../../peaks/beds/H3K27me3CerebellumLate1.bed` | BED 6-col, chr prefix | same |
| CTCF peaks | `../../peaks/CTCF.bed` | BED 10-col (narrowPeak), chr prefix | same |
| ATAC-seq diff | `../../peaks/atac_seq/ATAC_{up,down}.bed` | BED 3-col, chr prefix | same |
| MeCP2 diff | `../../peaks/mecp2/MeCP2_{up,down}.bed` | BED 3-col, chr prefix | same |
| CpG islands | `modality/mm10/gencode.vM25.mouse.cpg_islands.annotation.bed` | BED 4-col with header, NO chr prefix | `load_cpg_annotation()` |
| CpG shores | `modality/mm10/gencode.vM25.mouse.cpg_shores.annotation.bed.gz` | gzipped, same format | same |
| CpG shelves | `modality/mm10/gencode.vM25.mouse.cpg_shelves.annotation.bed.gz` | gzipped, same format | same |
| Gene coordinates | TxDb.Mmusculus.UCSC.mm10.knownGene + org.Mm.eg.db | Bioconductor packages | `get_gene_region()` |
| Coordinated genes | `plots/visualizations/tables/coordinated_changes.tsv` | TSV from section_05 | direct read |
| Chromatin states | `plots/visualizations/tables/dmr_chromatin_state_annotation.tsv` | TSV from section_10 | direct read |

## R Package Dependencies (beyond _shared_config.R)

- `Gviz` 1.51.0 — genome browser track rendering
- `GenomicInteractions` 1.42.0 — InteractionTrack for Hi-C arcs (needs seqinfo method patch)
- `readxl` — RNA-seq Excel file
- `AnnotationDbi` — gene symbol to Entrez ID mapping

## Key Gene Coordinates (mm10, from TxDb)

| Gene | Chr | Start | End | Width | Strand |
|------|-----|-------|-----|-------|--------|
| Syt1 | chr10 | 108,497,650 | 109,010,982 | 513 kb | - |

Other genes: use `get_gene_region()` with `extend_bp = 50000`.

## Syt1 Hi-C Loops (13 differential, within +/- 50kb)

Notable loops:
- loop_2627: chr10:108,450,000-109,500,000 (logFC=-0.70, **lost**, spans gene body)
- loop_4181: chr10:108,450,000-109,875,000 (logFC=-0.71, **lost**, extends past gene)
- loop_20413: chr10:108,730,000-109,010,000 (logFC=-0.44, **lost**, within gene body)
- loop_20167: chr10:108,725,000-109,875,000 (logFC=+0.65, **gained**, extends past gene)

Mix of 5 lost + 8 gained loops. Lost loops have larger effect sizes.
