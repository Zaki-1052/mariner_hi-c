# pub_browser.R — Publication-Quality Genome Browser Locus Tool

## Why this script exists

Our lab produces genome browser locus views (stacked signal tracks showing ChIP-seq, ATAC-seq, RNA-seq, etc. at a specific gene) for publication figures. The standard bioinformatics tool for this in R is **Gviz**, which renders functional but visually crude output: thick colored title bars, histogram-style tracks, cramped axis labels, and generally a lot of "chrome" that doesn't match our publication aesthetic.

Until now, every locus figure went through a manual **Adobe Illustrator cleanup pipeline**. The PI recorded a tutorial for undergraduates (see Reference Files below) documenting this workflow — roughly 18 minutes of removing clipping masks, deleting extraneous elements, adjusting line weights, redistributing track spacing, recoloring tracks, resizing gene models, and adding labels by hand. This process is:

- **Time-consuming**: 15-30 minutes per gene, per figure
- **Error-prone**: easy to accidentally scale ctrl and mut tracks independently, or to compress the gene model out of alignment with the signal tracks
- **Not reproducible**: each figure is a one-off Illustrator artifact
- **A bottleneck**: when the PI asks for "the same view but for 10 more genes," the manual pipeline doesn't scale

`pub_browser.R` eliminates this bottleneck. It produces SVG/PDF/JPEG output that matches the lab's publication style directly from the command line — no Illustrator post-processing required. The SVGs it produces are also Illustrator-friendly (editable text, clean structure via `svglite`), so manual touch-ups remain possible when needed.

## The target aesthetic

The script replicates the visual style seen in published figures from the lab (Ferguson et al., PLOS Genetics, 2024). The defining features are:

| Feature | Description |
|---|---|
| **Signal rendering** | Filled area curves (`geom_area`) with a thin matching-color outline. Not histogram bars. |
| **Color scheme** | One color per epigenetic mark. Control and mutant tracks use the *same color*; they're distinguished by vertical position, not hue. Colors carry semantic meaning: blue = H2AK119ub, green = H3K27ac (activating, "go"), red = H3K27me3 (repressive, "stop"). |
| **Left mark labels** | The mark name (e.g., "H2AK119ub") is rotated 90 degrees and placed to the left of its track pair, in the mark's color, bold. |
| **Condition labels** | "control" and the mutant genotype (e.g., "Bap1-KO") appear as small black horizontal text inside each track, anchored near the top-left. |
| **Track scaling** | A small "0-X" annotation in the mark's color appears at the top of the control track (e.g., "0-50", "0-100"). The number is always a "nice" round value ({1, 2, 5} x 10^n), chosen to give ~30% headroom above the data maximum. |
| **Gene model** | At the bottom: thin black backbone line for introns, filled black rectangles for exons (union of all isoforms), arrowheads on intron segments indicating strand direction, and the gene symbol in *italics* above the gene body. |
| **Scale bar** | A short black horizontal line with "X kb" label, positioned at the top-right. Auto-sized to ~8% of the view width if not specified. |
| **Highlights** | Optional semi-transparent gray rectangles behind the signal, spanning the full height of all tracks, for marking regions of interest. Can have text labels above them. |
| **Background** | Clean white everywhere. No colored panels, no axis ticks, no borders. |
| **Track spacing** | Uniform distribution — the PI specifically instructs to "distribute by the bottom" for even visual spacing between tracks (Illustrator transcript, ~5:38). |

## Reference files

These files informed the design of the script and should be examined by anyone iterating on it:

| File | What it shows |
|---|---|
| `FIG.png` (repo root) | PI's polished locus view of *Slc16a11* — 3 marks (H2AK119ub, H3K27ac, H3K27me3), ctrl/mut, gray highlight box, gene model at bottom. **This is the primary aesthetic target.** |
| `FIGK.png` (repo root) | PI's polished locus view of *Anxa2r1* — same 3 marks, different gene. Panel label "K" visible. Shows the style at a different locus with different signal characteristics. |
| `figs/journal.pgen.1011843.g002.PNG` | Published Figure 2. Panels F and G show locus views with 5 marks stacked, highlight boxes labeled "repressed enhancer" / "active enhancer", and multi-gene display at two different zoom levels. |
| `figs/journal.pgen.1011843.g003.PNG` | Published Figure 3. Panel A shows 4-track locus views comparing early/late timepoints (H2AK119ub early, H2AK119ub late, H3K27me3 early, H3K27ac early). |
| `figs/journal.pgen.1011843.g005.PNG` | Published Figure 5. Panel A shows locus views at euchromatic gene loci. Bottom-right panel K shows a 4-mark stacked view with a custom color for each mark. |
| `figs/journal.pgen.1011843.g007.PNG` | Published Figure 7. Panels A-C show cross-tissue comparisons (cerebellum, liver, kidney) — same mark, three tissues. Demonstrates the label flexibility needed (tissue names instead of genotypes). Also has colored (orange/salmon) highlight regions instead of gray. |
| `figs/journal.pgen.1011843.g008.PNG` | Published Figure 8. Panel A shows a wide-view (several Mb) with 5 marks, two gene rows, heterochromatin/euchromatin colored background blocks (pink/green), and "50 kb" scale bar. Panel B shows zoomed views with gray highlights. |
| `figs/journal.pgen.1011843.g009.PNG` | Published Figure 9. Panel A shows developmental-timepoint comparisons (early/late) with gray highlight boxes labeled "developmental gain" / "developmental reduction". |
| `transcript.txt` (repo root) | Full transcript of PI's Illustrator tutorial for undergraduates. Documents every manual step: removing clipping masks, scaling tracks together, distributing by bottom, gene model line weight (0.1pt) and height (0.05in), color choices (brighter red, olive green, royal blue), label placement, layer organization, and highlight rectangle workflow. **Read this to understand the PI's reasoning.** |
| `figures/ai/` | (To be populated) Illustrator `.ai` source files from the unpublished manuscript. These contain the polished versions of the figures above and additional locus views. Useful for pixel-level comparison when auditing script output. |

## Architecture decisions

### Why karyoploteR (previously ggplot2 + patchwork + cowplot)

The original implementation used ggplot2 + patchwork + cowplot to assemble panels manually. This worked but required complex multi-panel layout logic (cowplot::plot_grid to combine two patchwork columns, manual height arithmetic, separate left-label and right-data columns).

The current implementation uses **karyoploteR**, a Bioconductor package designed specifically for genomic locus visualizations. Key advantages:

- **Native genomic coordinates**: `plotKaryotype(zoom=...)` handles coordinate axes, clipping, and chromosome-aware layout automatically.
- **Built-in track stacking**: `r0/r1` parameters divide the vertical data area into tracks without manual panel arithmetic.
- **`kpArea` / `kpRect` / `kpArrows` / `kpPlotLinks`**: purpose-built primitives for signal curves, highlights, strand arrows, and Hi-C arcs.
- **`kpAddLabels`**: rotated left-margin labels are a one-liner instead of a custom cowplot panel.
- **Base R graphics**: output uses standard R devices (pdf/png/svg/jpeg) directly instead of ggsave, and SVGs are Illustrator-friendly via svglite.
- **Simpler code**: the rendering layer is ~200 lines instead of ~500, with no multi-column patchwork assembly.

**Previous alternatives considered (still rejected):**
- **Gviz**: colored title panels, histogram-style tracks, and excessive chrome don't match the PI's publication aesthetic.
- **plotgardener**: absolute inch-based page coordinates are inflexible for parameterized layouts.

**Trade-off vs the old ggplot2 version:**
- `ggtext` richtext labels (`*italic*^super^`) are not supported by karyoploteR (base R graphics). Use `--genotype-italic` for italic condition labels, or the font parameter in base R. For full richtext control, Illustrator touch-up on the exported SVG remains an option.

### Why standalone CLI (not a biomodal/downstream section script)

The existing section_46 script is tightly coupled to the BAP1 project:
- Sources `_shared_config.R` (which defines project-specific paths, KEY_GENES, etc.)
- Hardcodes BigWig paths for specific methylation/histone/RNA-seq samples
- Lives inside the biomodal/downstream pipeline

`pub_browser.R` is designed as a **general-purpose lab tool** that works for any project:
- All mark definitions come from CLI flags or YAML config
- No hardcoded file paths
- Works with any genome (specify `--txdb` and `--orgdb`)
- Lives at `scripts/pub_browser.R` (top-level, not inside any pipeline)

### Key rendering decisions

**kpArea for signal tracks**: karyoploteR's `kpArea()` draws filled area curves from binned BigWig data, matching the PI's aesthetic. The `base.y=0` parameter fills from zero upward, with `ymin/ymax` controlling the data range per track.

**r0/r1 proportional layout**: Each track occupies a vertical slice of the karyoploteR data panel. `compute_track_layout()` calculates proportional r0/r1 values for all tracks (signal pairs, diff panels, gene model, scale bar, Hi-C arcs) with appropriate gaps between mark groups.

**kpAddLabels for mark names**: Rotated colored labels on the left margin use karyoploteR's built-in `kpAddLabels(srt=90)`, replacing the custom cowplot panel approach.

**Manual gene model rendering**: Rather than using `kpPlotGenes()` (which renders all TxDb transcripts with default styling), the script draws gene models manually with `kpSegments` (backbone), `kpRect` (exons), `kpArrows` (strand direction), and `kpText` (italic symbol). This preserves the PI's thin-line aesthetic and focal-gene filtering.

**nice_ceiling for y-axis**: Track scaling like "0-50" or "0-100" uses round numbers from {1, 2, 5} x 10^n, chosen by rounding the data maximum up by ~30%. This gives headroom for labels above the peaks and matches the PI's manually-chosen scale values.

## Usage guide

### Prerequisites

R packages (all available from CRAN/Bioconductor):

```r
# Core rendering
BiocManager::install("karyoploteR")
install.packages(c("svglite", "dplyr"))

# Genomics
BiocManager::install(c("rtracklayer", "GenomicRanges", "GenomicFeatures"))

# For mm10 genes (install the appropriate TxDb for your genome)
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("org.Mm.eg.db")

# Optional: YAML config support
install.packages("yaml")

# Optional: colorspace (for desaturated diff-track loss color)
install.packages("colorspace")
```

### Basic usage

The minimum invocation needs a gene (or region), at least one mark, and an output path:

```bash
Rscript scripts/pub_browser.R \
  --gene Syt1 \
  --mark 'H2AK119ub:#2171B5:H2AK119ubCtrl.bw:H2AK119ubMut.bw' \
  --output out/syt1_browser
```

### Multi-mark view (matching FIG.png / FIGK.png)

The three marks used in most of our published locus figures, with the PI's richtext genotype label (italic *Bap1* + superscript f/f + plain Math1-cre):

```bash
Rscript scripts/pub_browser.R \
  --gene Syt1 \
  --mark 'H2AK119ub:#2171B5:H2AK119ubCtrl.bw:H2AK119ubMut.bw' \
  --mark 'H3K27ac:#41AB5D:H3K27acCtrl.bw:H3K27acMut.bw' \
  --mark 'H3K27me3:#CB181D:H3K27me3Ctrl.bw:H3K27me3Mut.bw' \
  --label 'control' \
  --label '*Bap1*^f/f^,Math1-cre' \
  --output figures/syt1_pub
```

Note the **repeatable `--label` flag** (passed exactly twice). The legacy `--labels 'ctrl,mut'` shortcut still works for simple cases but cannot contain commas inside a label — use the repeatable form whenever a genotype contains punctuation.

**Note on richtext**: The karyoploteR version uses base R graphics, so ggtext-style markdown (`*italic*`, `^super^`) is not rendered automatically. Use `--genotype-italic` to italicize the second condition label. For complex richtext formatting (italic gene names + superscript alleles), edit the exported SVG in Illustrator.

### Color reference (PI's semantic palette)

| Mark | Hex | Meaning |
|---|---|---|
| H2AK119ub | `#2171B5` | Royal blue |
| H3K27ac | `#41AB5D` | Green (activating = "go") |
| H3K27me3 | `#CB181D` | Red (repressive = "stop") |
| H3K4me3 | `#00CED1` | Cyan/teal (activating) |
| H3K27me1 | `#2171B5` | Blue (same family as K119ub) |
| RNA-seq | `#984EA3` | Purple |
| ATAC-seq | `#E6AB02` | Amber |
| MeCP2 | `#D95F02` | Orange |
| 5mC | `#D95F02` | Orange (distinct from Ub blue) |
| 5hmC | `#7570B3` | Purple (distinct from 5mC orange and Ub blue) |

These are defaults seen in published figures. Use any hex color — the script doesn't enforce a palette.

### Multi-replicate averaging

For marks with multiple replicate BigWigs per condition, pass comma-separated paths in the ctrl and mut slots:

```bash
Rscript scripts/pub_browser.R \
  --gene Syt1 \
  --mark '5mC:#2166AC:rep1_ctrl.bw,rep2_ctrl.bw,rep3_ctrl.bw,rep4_ctrl.bw:rep1_mut.bw,rep2_mut.bw,rep3_mut.bw,rep4_mut.bw' \
  --output out/syt1_mc
```

The script averages across replicates using `rowMeans` of binned coverage.

### Sparse binning for CpG methylation tracks

Standard BigWig binning uses `viewMeans` (sum of signal / bin width), which works for continuous ChIP-seq coverage. CpG methylation BigWigs are **sparse point data** — values exist only at individual CpG positions, with zero signal between them. When `viewMeans` averages a 780bp bin, it divides by the full bin width, diluting a 0.8 methylation fraction down to ~0.005.

The `sparse: true` YAML flag switches to **CpG-aware binning**: `sum(scores) / count(CpGs)` per bin. The result is the mean methylation fraction per CpG within each bin — the biologically meaningful quantity.

```yaml
marks:
  - name: '5mC'
    color: '#D95F02'
    sparse: true          # mean over CpG positions, not full bin width
    ctrl: [rep1_ctrl.bw, rep2_ctrl.bw, rep3_ctrl.bw, rep4_ctrl.bw]
    mut:  [rep1_mut.bw,  rep2_mut.bw,  rep3_mut.bw,  rep4_mut.bw]
    average: true
```

Use `sparse: true` for any BigWig containing per-site point data (CpG methylation, hydroxymethylation, etc.). Do **not** use it for continuous coverage tracks (ChIP-seq, ATAC-seq, RNA-seq) — those are correctly handled by the default `viewMeans` binning.

### Percent difference tracks

When ctrl and mut signals are visually similar (common with CpG methylation), a difference track reveals the actual changes. The `diff: true` YAML flag adds a third panel below each mark pair showing `(mut - ctrl)` as a diverging area chart:

- Gains (mut > ctrl): filled in the mark's color above the zero line
- Losses (mut < ctrl): filled in a desaturated shade below the zero line
- Symmetric y-axis with `nice_ceiling` rounding

For sparse/fraction marks (`sparse: true`), the difference is automatically multiplied by 100 and displayed as percentage points (e.g., ±5% means 5 percentage point change in methylation).

```yaml
marks:
  - name: '5mC'
    color: '#D95F02'
    sparse: true
    diff: true            # add (mut - ctrl) diverging panel
    ctrl: [...]
    mut:  [...]
```

Requires `colorspace` (installed with ggplot2). The diff panel uses `colorspace::desaturate(colorspace::lighten(...))` for the loss color.

### Multiple focal genes

By default, the gene model panel shows either a single focal gene (from `--gene` or `--gene-symbol`) or all overlapping genes (`--show-all-genes`). The `focal_genes` option provides a middle ground: specify exactly which genes to display.

```yaml
focal_genes: [Syt1, Pawr]       # show only these two genes
```

Or via CLI (repeatable flag):

```bash
--focal-gene Syt1 --focal-gene Pawr
```

Each symbol is resolved via OrgDb → entrez → TxDb, so strand and exons are correct. Genes not overlapping the view region are silently skipped. If none of the specified genes overlap, falls back to showing all overlapping genes.

`focal_genes` takes priority over `--show-all-genes` and `--gene-symbol` for gene model filtering. The `--gene-symbol` flag still controls region resolution metadata (printed in the startup echo).

### Region workflow (iterate from gene to coordinates)

The PI's preferred workflow is: render with `--gene`, look at the output, decide the view needs to be extended (e.g., to capture an entire H3K27me3 domain that exceeds the gene body), then re-render with `--region` using explicit coordinates. The script makes this seamless.

**Step 1 — render with a gene name** (script picks the gene's coordinates + the `--extend` flank):

```bash
Rscript scripts/pub_browser.R \
  --gene Anxa2r1 \
  --extend 10000 \
  --mark 'H3K27me3:#CB181D:ctrl.bw:mut.bw' \
  --output out/anxa2r1_gene
```

The script prints the resolved region so you know what coordinates to adjust:

```
Resolving region...
  chr13:120,024,605-120,062,194 (37.6 kb)
  gene model: focal gene Anxa2r1 (entrez 100861753)
```

**Step 2 — re-render with explicit coordinates** to extend the view as needed:

```bash
Rscript scripts/pub_browser.R \
  --region 'chr13:119990000-120070000' \
  --gene-symbol Anxa2r1 \
  --mark 'H3K27me3:#CB181D:ctrl.bw:mut.bw' \
  --output out/anxa2r1_extended
```

`--gene-symbol` is the key piece: when paired with `--region`, the script resolves the symbol through OrgDb→TxDb to fetch the gene's strand and exons. The gene model panel then shows just that focal gene (with correct arrow direction and italic name) instead of every gene in the window.

**Three behaviors when using `--region`:**

| Flags | Gene model shows |
|---|---|
| `--region` alone | All genes overlapping the view (no focal selection) |
| `--region --gene-symbol Foo` | Only Foo, with strand/exons resolved from TxDb |
| `--region --show-all-genes` | All genes overlapping the view (explicit) |

A startup message confirms the resolved mode so you can sanity-check before opening the figure.

### Region-based view without a focal gene

```bash
Rscript scripts/pub_browser.R \
  --region 'chr11:108400000-109100000' \
  --gene-symbol Syt1 \
  --show-all-genes \
  --mark 'H2AK119ub:#2171B5:ctrl.bw:mut.bw' \
  --output out/syt1_region
```

### Highlight regions

```bash
Rscript scripts/pub_browser.R \
  --gene Slc16a11 \
  --mark 'H2AK119ub:#2171B5:ctrl.bw:mut.bw' \
  --highlight 'chr11:70135000-70140000' \
  --highlight-label 'differential region' \
  --output out/slc16a11_highlight
```

Multiple `--highlight` flags are supported (paired with `--highlight-label` flags).

### YAML config (recommended for complex setups)

For figures with many marks or batch processing, use a YAML config file:

```yaml
# configs/syt1_locus.yaml
genome: mm10
gene: Syt1
extend_bp: 50000
# YAML arrays sidestep the comma-split issue, so richtext is safe here too
labels: ['control', '*Bap1*^f/f^,Math1-cre']
output: figures/syt1_pub
scale_bar_kb: 50

marks:
  - name: H2AK119ub
    color: '#2171B5'
    ctrl: /path/to/H2AK119ubCtrl.bw
    mut:  /path/to/H2AK119ubMut.bw
  - name: H3K27ac
    color: '#41AB5D'
    ctrl: /path/to/H3K27acCtrl.bw
    mut:  /path/to/H3K27acMut.bw
  - name: H3K27me3
    color: '#CB181D'
    ctrl: /path/to/H3K27me3Ctrl.bw
    mut:  /path/to/H3K27me3Mut.bw

highlights:
  - region: 'chr10:108700000-108800000'
    label: 'K119ub gain region'
```

Run with: `Rscript scripts/pub_browser.R --config configs/syt1_locus.yaml`

CLI flags override YAML values when both are provided.

### Hi-C loop arcs (optional overlay)

```bash
Rscript scripts/pub_browser.R \
  --gene Syt1 \
  --mark 'H2AK119ub:#2171B5:ctrl.bw:mut.bw' \
  --hic-loops peaks/loop_annotation_extended/late/extended_characterized_loops.tsv \
  --output out/syt1_with_loops
```

The TSV needs columns: `chr1, start1, end1, chr2, start2, end2, direction` where direction is `up_in_mutant` or `down_in_mutant`.

### Layout tuning

```bash
# Wider figure with taller tracks
Rscript scripts/pub_browser.R \
  --gene Syt1 \
  --mark 'H2AK119ub:#2171B5:ctrl.bw:mut.bw' \
  --width 10 \
  --track-height 0.6 \
  --output out/syt1_wide

# Smaller flanking region (tighter zoom)
Rscript scripts/pub_browser.R \
  --gene Anxa2r1 \
  --extend 5000 \
  --mark 'H2AK119ub:#2171B5:ctrl.bw:mut.bw' \
  --output out/anxa2r1_tight

# Explicit scale bar size
Rscript scripts/pub_browser.R \
  --gene Syt1 \
  --mark 'H2AK119ub:#2171B5:ctrl.bw:mut.bw' \
  --scale-bar 100 \
  --output out/syt1_100kb_bar
```

### Using a different genome

```bash
Rscript scripts/pub_browser.R \
  --gene TP53 \
  --genome hg38 \
  --txdb TxDb.Hsapiens.UCSC.hg38.knownGene \
  --orgdb org.Hs.eg.db \
  --mark 'H3K27me3:#CB181D:ctrl.bw:mut.bw' \
  --output out/tp53_human
```

### Output

The script creates a subfolder at the output path containing four files:
```
out/syt1_pub/
  syt1_pub.pdf   # Vector (publication submission)
  syt1_pub.svg   # Vector (Illustrator-editable, text as text not paths)
  syt1_pub.png   # Raster 300dpi (web, markdown, transparent background)
  syt1_pub.jpg   # Raster 300dpi (slides, previews)
```

## All CLI flags

| Flag | Type | Default | Description |
|---|---|---|---|
| `--gene` | string | — | Gene symbol (resolved via TxDb) |
| `--region` | string | — | Explicit region as `chr:start-end` |
| `--genome` | string | `mm10` | Genome assembly name |
| `--mark` | string | — | `'name:color:ctrl.bw:mut.bw'` (repeatable) |
| `--config` | path | — | YAML config file |
| `--label` | string | — | Condition label (repeatable; pass exactly twice). Supports richtext markdown via ggtext: `*italic*`, `**bold**`, `^super^`, `~sub~`, `<sup>...</sup>`, `<sub>...</sub>`. **Preferred over `--labels`** when a genotype contains commas. |
| `--labels` | string | `control,mutant` | Legacy comma-separated shortcut. Cannot contain commas in either label; use `--label` (repeatable) for richtext or genotypes with commas. |
| `--output` | path | `out/locus` | Output path prefix |
| `--extend` | int | `50000` | Flanking bp on each side |
| `--scale-bar` | float | auto | Scale bar length in kb |
| `--txdb` | string | `TxDb.Mmusculus.UCSC.mm10.knownGene` | TxDb package |
| `--orgdb` | string | `org.Mm.eg.db` | OrgDb package |
| `--bin-size` | int | auto | BigWig binning resolution (bp) |
| `--width` | float | `7.0` | Figure width (inches) |
| `--track-height` | float | `0.40` | Per-track height (inches) |
| `--highlight` | string | — | `chr:start-end` highlight box (repeatable) |
| `--highlight-label` | string | — | Label above highlight (repeatable, paired) |
| `--hic-loops` | path | — | Hi-C loops TSV |
| `--show-all-genes` | flag | `false` | Show all overlapping genes (vs. just the focal gene) |
| `--focal-gene` | string | — | Gene symbol to include in gene model (repeatable). Shows exactly these genes, resolved via OrgDb→TxDb. Takes priority over `--show-all-genes`. |
| `--gene-symbol` | string | — | Displayed gene symbol. With `--region`, this also resolves the focal gene from TxDb so the gene model shows just that gene with correct strand/exons. |
| `--genotype-italic` | flag | `false` | Legacy: italicize second condition label as a whole. Superseded by `--label '*genotype*'`. |

## For AI sessions: visual refinement notes

If you are an AI assistant working on this script in a future session, read this section.

**Goal**: The output should match the published figures in `figs/` and the reference figures `FIG.png`/`FIGK.png` as closely as possible. The PI's Illustrator workflow is documented in `transcript.txt` — read it to understand the reasoning behind every aesthetic choice.

**What to compare against**: Open the script output side-by-side with `FIG.png` and `FIGK.png`. These are the primary targets. The `figs/journal.pgen.*.PNG` files show the full range of published variants. The `figs/ai/` directory contains PNGs of unpublished-manuscript Illustrator files with additional locus views.

### karyoploteR migration (2026-06)

Replaced the ggplot2 + patchwork + cowplot rendering engine with karyoploteR. The CLI interface, YAML config, BigWig processing, and output format convention are unchanged. Key differences:

- **Rendering engine**: all signal tracks, gene models, highlights, and Hi-C arcs now rendered via karyoploteR primitives (`kpArea`, `kpRect`, `kpArrows`, `kpPlotLinks`, etc.) instead of ggplot2 geom layers.
- **Layout**: proportional `r0/r1` track positions calculated by `compute_track_layout()` instead of patchwork height vectors and cowplot grid assembly.
- **Left-margin labels**: `kpAddLabels(srt=90)` replaces custom cowplot label panels.
- **Output**: base R devices (`pdf()`, `png()`, `jpeg()`, `svglite()`) replace `ggsave()`.
- **Richtext dropped**: ggtext-style markdown labels (`*italic*^super^`) are not supported by base R graphics. Use `--genotype-italic` for italic condition labels, or edit SVG output in Illustrator for complex formatting.
- **Dependencies removed**: ggplot2, patchwork, cowplot, grid, ggtext no longer required.
- **Dependencies added**: karyoploteR (Bioconductor).

### Previous audit notes (2026-05)

These changes from the ggplot2 era are preserved in the karyoploteR version:

- **Gene model**: thin backbone/arrows (lwd 0.3), slim exons, gene panel height ~55% of a signal track.
- **Scale bar**: bold (lwd 3.0), auto-sized to ~8% of view width.
- **Repeatable `--label` flag**: fixes comma-in-genotype parsing bug.
- **`--region` + `--gene-symbol`**: resolves focal gene from TxDb for correct strand/exons.
- **Runtime echo**: prints resolved region and gene-model mode.
- **Robustness**: numeric arg validation, missing-file checks, TxDb/OrgDb install guidance.

### Methylation and multi-gene enhancements (2026-06)

This pass added support for CpG methylation tracks and multi-gene views. Reference render: `figures/syt1_pawr_methylation/`.

Changes made in this pass:

- **Sparse binning mode** (`sparse: true`): per-mark YAML flag that switches `import_bigwig_binned` from `viewMeans` (sum/bin_width) to `viewSums(scores) / viewSums(occupancy)` — mean over CpG positions only. Fixes the signal dilution bug where methylation fractions (0-1 per CpG) were averaged over full bin widths (mostly zero between CpGs), producing near-zero values (e.g., 0.005 instead of 0.65). Backward-compatible: default is `FALSE`, existing ChIP-seq marks unaffected.
- **Percent difference tracks** (`diff: true`): per-mark YAML flag that adds a third panel below each ctrl/mut pair showing `(mut - ctrl)` as a diverging area chart. For sparse marks, values are automatically converted to percentage points (×100). Uses `colorspace::desaturate(colorspace::lighten(...))` for the loss color. Assembly and left-column labels adjusted to span 3 rows when a diff panel is present.
- **Multiple focal genes** (`focal_genes: [Syt1, Pawr]` / `--focal-gene`): YAML array or repeatable CLI flag that specifies exactly which genes to show in the gene model panel. Each symbol is resolved via OrgDb→TxDb. Takes priority over `--show-all-genes` and single `--gene-symbol`. Solves the problem where `--show-all-genes` pulled in pseudogenes (Gm36283, Brcc3dc, 1r12a) alongside the two genes of interest.
- **BigWig path update**: per-sample methylation BigWigs moved from the old `/Users/zakiralibhai/Documents/BIO_LAB/methylation-tracks/` path to `/Users/zakiralibhai/sdsc/bigwigs/methylation/{mc,hmc}/`.

### Known areas for potential refinement (remaining)

- Condition label font size and positioning for very long richtext strings (multi-line wrap not yet supported)
- Highlight rectangle colors — currently gray-only; published figures (g007, g008 panel A) use pink/green for functional annotations. Adding `--highlight-color` per region would close this gap.
- Panel letters (the "K" in `FIGK.png`, the A/B/C/... in `figs/journal.pgen.*.PNG`) — would need a `--panel-label` flag.
- Signal curve smoothness — `--bin-size` controls this; defaults are reasonable but extreme zooms may need explicit overrides.

### What NOT to change

- The same-color-per-mark principle (this is non-negotiable in the PI's style)
- The rotated left labels (this is the defining feature)
- The general layout (left labels | signal tracks | gene model at bottom)
- The `nice_ceiling` rounding logic (produces y-axis values matching PI's manual choices)
- The repeatable `--label` flag — comma-split parsing was a bug, do not reintroduce it.
- The sparse binning algorithm (`viewSums / cpg_counts`) — `viewMeans` produces biologically wrong values for CpG methylation.
- The karyoploteR rendering engine — don't revert to ggplot2 unless there's a specific feature requirement that karyoploteR cannot support.

### How to test

Run the script against the BigWig files in `/Users/zakiralibhai/sdsc/bigwigs/` (H2AK119ub, H3K27ac, H3K27me3, H3K4me3, ATAC, RNA, DNA methylation — ctrl/mut pairs are available). The canonical command used in the audit is:

```bash
Rscript scripts/pub_browser.R \
  --gene Syt1 \
  --mark 'H2AK119ub:#2171B5:/Users/zakiralibhai/sdsc/bigwigs/H2AK119ubCtrl.bw:/Users/zakiralibhai/sdsc/bigwigs/H2AK119ubMut.bw' \
  --mark 'H3K27ac:#41AB5D:/Users/zakiralibhai/sdsc/bigwigs/H3K27acCtrl.bw:/Users/zakiralibhai/sdsc/bigwigs/H3K27acMut.bw' \
  --mark 'H3K27me3:#CB181D:/Users/zakiralibhai/sdsc/bigwigs/H3K27me3Ctrl.bw:/Users/zakiralibhai/sdsc/bigwigs/H3K27me3Mut.bw' \
  --label 'control' --label '*Bap1*^f/f^,Math1-cre' \
  --output test_pub_browser/syt1_v2
```

Compare the resulting JPEG/PNG to `FIG.png` (PI's polished Slc16a11 view). The most important visual checks: italic+superscript genotype label, thin gene model, bold scale bar, no Illustrator post-processing needed.
