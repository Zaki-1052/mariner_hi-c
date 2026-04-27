# CONTEXT-CLUSTER.md

Comprehensive context document for the Popay-style Hi-C loop clustering analysis adapted for BAP1-KO cerebellum. This file plus `PLAN-CLUSTER.md` provide full context for AI sessions -- no need to @ individual project files.

**Last updated:** 2026-04-26

---

## 1. Project Abstract (Verbatim)

> BAP1 is a tumor-suppressor gene that removes the histone modification H2AK119ub, traditionally understood to mediate gene repression. Its loss causes severe neurodevelopmental defects including epilepsy, but how it shapes the three-dimensional folding of the genome has never been characterized in the brain. Given that H2AK119ub has been found to modulate active enhancer levels, we used Hi-C to map the physical contacts of DNA regions. We found genome-wide chromatin interactions to be perturbed in the BAP1 mutant, where 1 in 5 detected chromatin loops were dysregulated in the adult mouse cerebellum. At over 200 loci, long range loops associated with repression were preferentially lost in exchange for shorter range contacts. Additionally, the presence of H2AK119ub strongly predicted this distance-dependent loop-loss: if the histone modification was present at an anchor, contacts were ten times more likely to become disrupted. Analysis of differentially expressed genes revealed repression of developmental and synaptic genes where connections to enhancers were broken. This effect was progressive over neurodevelopment, expanding from under 200 differential loops in early development to almost three thousand in adulthood. These findings indicate that BAP1 is required for organization of the 3D genome in the developing brain, and its loss leads to architectural changes resulting in dysregulated synaptic gene expression. We propose that elevation of H2AK119ub, as a result of BAP1 loss, collapses long-range developmental loops, replacing them with proximal repressive contacts.

**Organism:** Mouse (mm10). **Conditions:** BAP1-KO mutant vs wildtype control, n=3 replicates per condition. **Timepoints:** 250831 (early/P12), 250402 (late/adult).

---

## 2. Why This Analysis: The Central Biological Question

The Dixon meeting (2026-04-10) and CTEA April meeting identified a critical unresolved question:

**Is Polycomb (H3K27me3) enrichment at the anchors of differential loops, or across the loop body/span?**

- **Anchor-specific enrichment** = sensitivity model: CTCF binding sites are directly disrupted by repressive chromatin state changes
- **Body/span enrichment** = extrusion impediment model: Polycomb spreading blocks cohesin loop extrusion, preventing formation of longer loops

This distinction is the mechanistic heart of the paper. The Popay anchor-vs-span ChromHMM analysis (Fig 2f of Popay et al.) directly tests this by comparing chromatin state enrichment at loop anchors vs. loop spans across k-means clusters of loops grouped by their sensitivity to perturbation.

### Key findings that motivate this analysis

- Lost loops are systematically longer than gained loops (median 625 kb vs 320 kb, KS p = 5.49e-48)
- 92.4% of K27ac-anchored loops are lost; 74.6% of K27me3-anchored loops are gained
- For CTCF-CTCF loops: lost are longer than gained (median 810 kb vs 325 kb, p = 8.92e-37)
- K27ac-anchored loops show NO length difference (p = 0.41) -- the length effect is specific to structural loops
- Both Popay et al. and Yu et al. show the same pattern: long active loops lost, short repressive contacts gained
- Our K119ub data provides a unique mechanistic link (H2AK119ub elevation -> Polycomb spreading -> extrusion impediment)

---

## 3. Dixon Meeting Summary (2026-04-10)

### What Jesse Dixon validated
- Our cutoffs (logFC > 0.3, FDR < 0.05) are reasonable; fold changes of 3-4 are rare in Hi-C, 1.3-1.5 fold is standard
- The magnitude of Hi-C changes is large and surprising given normal histology
- HICCUPS peak pixel selection bias is corrected by mariner
- This is enough to tell a story even without a full mechanism

### Key finding -- loop type distribution
- Most lost loops: active enhancer-active enhancer
- Most gained loops: CTCF-CTCF, poised enhancer-CTCF, Polycomb-Polycomb
- 92.4% of K27ac-anchored loops are lost, 74.6% of K27me3-anchored loops are gained, 82.7% of K4me3-anchored loops are lost

### Jesse's interpretation of loop length
Lost loops may break into two classes: shorter CRE-type/enhancer-promoter loops (K27ac-associated) and longer architectural/Polycomb loops. He connected this to the Brad Bernstein IDH mutant glioma work -- repressive chromatin gains knock out CTCF binding sites.

### Tessa Popay suggestion (highest priority)
Jesse's biggest concrete suggestion: look at the body/span of loops, not just anchors. His postdoc Tessa Popay did exactly this for cohesin depletion -- metagene-style ChIP-seq plots from loop anchors showed K27me3 enrichment *flanking* CTCF sites. Jesse emailed Tessa and CC'd us.

### Other suggestions
- Stripe analysis with Stripenn (DONE -- see Section 8)
- Insulation score investigation with Cooltools
- ABC analysis cleanup (DONE)
- PC1 scatter plot not volcano (DONE)
- SNIPER for subcompartments
- Developmental comparison (DONE -- partial blocking at compartment level)

---

## 4. CTEA April Meeting Notes

Key directives from the lab meeting:

- **Superenhancers within 2MB window of DEGs:** How frequently does gene contact SE or enhancers (APA)? Sub-classify contact frequency by K27ac change.
- **Center on Rad21/CTCF peak/motif within differential anchor:** Are changes happening at the actual motif? 10kb anchors are large; center annotation around the CTCF region.
- **CTCF motif filter:** At least one anchor must have a CTCF motif. Re-run annotation with this filter.
- **DEG-centric analysis:** Pull DEGs and invariants as comparison. Contact strength between just those genes and enhancers. May clear up ABC plot.
- **Stripes as SE hubs:** Stripes can be super-enhancer hubs (enhancers modifying multiple promoters).
- **Polycomb flanking anchors:** "extending along euchromatin stretch/span, but then hits polycomb K27me3 and stops" -- test at differential stripe anchors.
- **Shorter CTCF loops:** Why loops vs stripes? Cohesin trying to extrude/grow longer but more polycomb regions blocking it.

---

## 5. Popay Collaboration Context

### Email chain (received 2026-04-21)
Tessa Popay (tpopay@salk.edu) is a postdoc in Jesse Dixon's lab at the Salk Institute. She did her PhD on HCF-1, which associates with BAP1 in the PR-DUB complex. She shared her analysis code and offered to help.

### What Popay shared
- GitHub repo: `https://github.com/tpopay/HiC-clustering/tree/main` (cloned to `cluster/`)
- Key notebook: `grouped_loops_figures.ipynb` > "ChromHMM relationship" section
- Format guidance: BEDPE with first 6 columns as coordinates
- Conda environment specs: `cluster/environment_mac.yml` / `environment_linux.yml`
- ChromHMM OverlapEnrichment used as subprocess from the notebook
- Example data link for formatting reference

### What we need to adapt
Popay's pipeline was built for hTERT RPE-1 cells with NIPBL depletion (dTAG system). We need to adapt for:
- mm10 genome (her data is hg38)
- BAP1-KO vs wildtype (not DMSO vs dTAG time course)
- 5 histone marks instead of her 12-state ChromHMM from RPE-1
- No RAD21 ChIP data (she normalizes to RAD21)
- 2 conditions (ctrl_merged, mut_merged) instead of 3 (DMSO, 4hr, 24hr)

---

## 6. Popay et al. (Nat Genet 2026) -- Key Findings

**Paper:** "Acute NIPBL depletion reveals in vivo dynamics of loop extrusion and its role in transcription activation" -- Popay, Pant, Munting, Tastemel, Black, Haghani & Dixon. Nature Genetics 58, 869-880 (April 2026).

### Experimental system
- Acute NIPBL depletion via dTAG in hTERT RPE-1 cells
- Also RAD21 depletion and CTCF degron for comparison
- Hi-C at DMSO, dTAG 4h, dTAG 24h
- SLAM-seq for nascent transcription

### Six k-means clusters of chromatin loops (Fig 1c-d)
Loops classified by cohesin/RAD21 dependence into 6 clusters:
1. **Minimum-dependency (min-dep):** n=1,171 -- least affected by NIPBL loss
2. **Low-dependency:** n=3,164
3. **Medium-dependency (med-dep):** n=3,931
4. **Mixed-dependency:** n=1,757 -- unique class, stronger reduction at 24h than 4h
5. **High-dependency:** n=4,124
6. **Maximum-dependency (max-dep):** n=2,662 -- most sensitive to NIPBL loss

### The critical Figure 2f -- ChromHMM anchor vs span enrichment
- Used ChromHMM 12-state model to classify chromatin at anchors and spans separately
- Less-dependent clusters: enhancer chromatin states at anchors AND spans
- More-dependent clusters: promoter states at anchors, enhancer states within loops
- **Mixed-dependency cluster: anchors devoid of K27me3 but STRONG K27me3 signal flanking the anchor** (Fig 2g)
- This is the "extrusion impediment" signature we want to test in our data
- Repressive chromatin is the majority state for ~65% of loop anchors within the mixed-dependency cluster

### Other relevant findings
- NIPBL-dependent loops are longer; less-dependent are shorter (median loop size ~450kb for mixed-dep)
- Structural loops (cohesin/CTCF at both anchors): 38.9% for low-dep to 53.6% for mixed-dep
- CRE loops (enhancer/promoter at both anchors): 25.2% for mixed-dep, 9.7% for max-dep
- STAG1 stabilizes persistent loops; STAG2 at less-dependent loops
- SMC3 acetylation higher at less NIPBL-dependent loops
- NIPBL-dependent genes have stripe-like features emanating from TSS

### Methodological details from the paper
- Clustering: Cluster 3.0, k-means with k=6, normalized balanced counts to DMSO control
- ChromHMM: 12-state model, OverlapEnrichment run separately for anchors and spans
- Loop classification: structural (CTCF+RAD21 at both anchors), CRE (enhancer/promoter), mixed, unclassified
- ChIP signal: RPKM-normalized BigWig, mean at RAD21 peaks overlapping anchors, normalized to RAD21
- deepTools: computeMatrix reference-point for metagene plots at anchors +/-5kb
- Cooltools: off-diagonal pileup for aggregate Hi-C signal per cluster
- Gene annotation: bedtools closest with promoter-proximal BED

---

## 7. Cluster 3.0 Reference

### Installation
Installed locally. Available on PATH as `cluster` (the compiled executable).

### CLI usage (what Popay's pipeline uses)
```bash
cluster -f input_matrix.txt -r 100 -g 7 -k 6
```

**Flags used by the pipeline:**
- `-f filename` -- input data file (tab-delimited, first column = row IDs, first row = column headers)
- `-r number` -- number of runs/trials for k-means (100 in Popay's pipeline)
- `-g [0..9]` -- distance measure for gene clustering: 7 = Euclidean distance
- `-k number` -- number of clusters for k-means

**Distance measure codes:** 0=none, 1=uncentered corr, 2=Pearson corr, 3=uncentered abs, 4=Pearson abs, 5=Spearman, 6=Kendall, 7=Euclidean, 8=city-block

### Input format
Tab-delimited text file. First column = row identifier (ID). First row = column labels. Data cells = numeric values. Missing values = empty cells.

```
UNIQID	col1	col2	col3
row1	1.0	0.5	0.8
row2	0.3	0.7	0.2
```

### Output files
- `JobName_K_G{k}.kgg` -- gene clustering assignments (two columns: ID and GROUP)
- `JobName_K_G{k}.cdt` -- clustered data table (reordered)

The `.kgg` file is what Popay's pipeline reads back to merge with the original coordinates.

---

## 8. ChromHMM v1.27 Reference

### Installation
Installed at `cluster/ChromHMM/ChromHMM.jar`. Wrapper script at `cluster/ChromHMM/chromhmm`:
```bash
#!/usr/bin/env bash
CHROMHMM_JAR="/Users/zakiralibhai/Documents/GitHub/mariner_hi-c/cluster/ChromHMM/ChromHMM.jar"
exec java -mx8G -jar "$CHROMHMM_JAR" "$@"
```
The wrapper is on PATH via `~/.zshrc`. Invocation: `chromhmm <Command> [options]`.

Chromosome sizes for mm10 are at `cluster/ChromHMM/CHROMSIZES/mm10.txt`.

### Commands used by this pipeline

**BinarizeBed** -- converts peak BED files to binarized data for model learning:
```bash
chromhmm BinarizeBed -peaks \
    cluster/ChromHMM/CHROMSIZES/mm10.txt \
    <inputbeddir> \
    <cellmarkfiletable> \
    <outputbinarydir>
```
- `-peaks` flag: treat BED files as peak calls (give '1' to any bin overlapping a peak)
- `cellmarkfiletable`: tab-delimited, columns = `cell  mark  bedfile.bed`
- Output: `CELL_CHROM_binary.txt` files in outputbinarydir

**LearnModel** -- learn chromatin state model from binarized data:
```bash
chromhmm LearnModel -p 4 \
    -l cluster/ChromHMM/CHROMSIZES/mm10.txt \
    <inputdir> <outputdir> <numstates> <assembly>
```
- `-p 4`: use 4 processors
- `-l`: chromosome length file (required for proper segmentation)
- Outputs: `emissions_*.txt`, `transitions_*.txt`, `*_segments.bed`, webpage summary
- Default binsize: 200bp. Default max iterations: 200.

**OverlapEnrichment** -- compute fold enrichment of each state for external coordinates:
```bash
chromhmm OverlapEnrichment \
    -noimage -uniformscale \
    -m <labelmappingfile> \
    -f <coordlistfile> \
    -colfields 0,1,2 \
    <inputsegment> <inputcoorddir> <outfileprefix>
```
- `-noimage`: suppress image output (we make our own heatmaps)
- `-uniformscale`: all columns on same color scale in heatmap
- `-m <file>`: two-column tab-delimited mapping state IDs (e.g. `U1`) to names (e.g. `Active_Promoter`)
- `-f <file>`: lists BED filenames to compute enrichment for, controls display order
- `-colfields 0,1,2`: specify chromosome, start, end column indices
- Fold enrichment = (C/A)/(B/D) where A=bases in state, B=bases in annotation, C=bases in both, D=bases in genome
- Outputs: `outfileprefix.txt` (enrichment table) and optionally `.png`/`.svg` heatmaps

---

## 9. Popay Python Modules (in `cluster/`)

### `cluster_tools.py`
- `elbow(out_dir, count_matrix, data_cols)` -- produces elbow plot for k=1..14
- `sort_clusters(cluster_sort_df, data_cols)` -- sorts cluster names by descending mean signal
- `comparison_type(data_cols)` -- determines if data is 'pairwise' (cols have `_` separator like `DMSO_4hr`) or 'multiple'
- `sort_by_strings(df, sort_col, order)` -- sort dataframe by custom string order

### `plotting.py`
- `heat(count_matrix, data_cols, ...)` -- clustered heatmap with optional proportional subplots
- `line(melted_df, xcol, ycol, ...)` -- line plots per cluster
- `box(melted_df, xcol, ycol, ...)` -- box plots with optional subplots
- `strip(melted_df, xcol, ycol, ...)` -- strip plots
- `stacked(count_table, ...)` -- stacked bar charts (used for ChromHMM proportions, loop classification)
- `bar(melted_df, xcol, ycol, ...)` -- grouped bar charts
- `joint(melted_df, xcol, ycol, ...)` -- joint distribution plots (scatter + marginals)
- `randomize_hex()` -- random color generation

**Bug (line 24-25):** hardcoded path `/Users/tessapopay/example_data/custom_params.json`. Fix: use `os.path.join(os.path.dirname(os.path.abspath(__file__)), 'custom_params.json')`.

### `chromHMM_heatmap.py`
- `heatmap_plot(path, normalize=False)` -- reads ChromHMM OverlapEnrichment output `.txt` and renders heatmap

**Bug (line 9-10):** same hardcoded `custom_params.json` path.

### `deeptools_plotting.py`
- `heatmap_plot(matrix_file, ...)` -- renders deeptools computeMatrix output as heatmap

**Bug (line 8-9):** same hardcoded path. Also line 108: hardcoded `if bam == 'HA-NIPBL': line.set_ylim([0,2])`.

### `deepTools_pipeline.py`
- `bam_coverage(bam_dict, ...)` -- runs `bamCoverage` to generate BigWig files
- `bed_pileup(out_dir, bigWig_dict, bam_dict, bed_dict, up_down, ...)` -- runs `computeMatrix reference-point` + `plotHeatmap` for metagene analysis at anchor BEDs

**Bug (line 15):** hardcoded `--effectiveGenomeSize 2913022398` (hg38, should be mm10: 2494787188).

### `cooltools_called.py`
- `mcool_pileup(mcool_dict, bedpe_dict, out_dir, ...)` -- cooltools off-diagonal pileup analysis

**Bug (line 16):** `plt.style.use('seaborn-poster')` removed in matplotlib 3.6+. **Bug (line 18-19):** hardcoded custom_params path. **Bug (line 27):** default `genome='hg38'`, should be `'mm10'`.

### `bedpe_analysis.py`
- `loop_size(out_dir, bedpe_dict, logY=False)` -- compute and plot loop sizes per group
- `bedtools_annotation(out_dir, bedpe_dict, FPKM_df=None, temp_dir=None)` -- annotate loop anchors with nearest genes using bedtools

**Bug (line 67):** hardcoded gene annotation path. **Bug (lines 142-159):** crashes when `FPKM_df=None`.

### `statistics_functions.py`
- `kruskal_wilcoxon(data_names, data_list)` -- Kruskal-Wallis test + pairwise Wilcoxon rank-sum

### `custom_params.json`
Matplotlib rcParams for consistent styling (tick sizes, label sizes, font). Located at `cluster/custom_params.json`.

---

## 10. Stripenn Pipeline Results Summary

All stages 0-7 complete. Results are final.

### 250402 (late/adult) -- strong differential signal
- 7,371 union stripes at 5kb, 2,320 significant (31.5% FDR<0.05)
- More gained than lost (2,052 vs 1,528)
- High-confidence: 367 lost + 638 gained = 1,005 stripes
- All effect sizes minimal (max |logFC| = 0.389), BCV very low (0.011-0.021)
- Cross-res logFC correlation: r=0.850
- H3K27me3+ anchors 1.75x enriched in gained vs lost (9.1% vs 5.2%)
- Active_Enhancer enriched in lost (19.5% vs 14.0%)
- GO enrichment: gained stripes -- developmental TFs (Hox, Shh, Dlx5/6), synaptic genes, ion channels

### 250831 (early/P12) -- weak differential signal
- 4,008 union stripes at 5kb, 96 significant (2.4% FDR<0.05)
- More lost than gained (949 vs 776)
- Bivalent_Promoter anchors 3x more prevalent than adult (3.6% vs 1.2%)
- GO/KEGG: lost stripes -- stem cell pluripotency signaling (Wnt5a, Fgfr2, Klf4)

### Cross-timepoint pattern
Directional bias reverses between timepoints: more lost at P12, more gained in adult. Consistent with progressive Polycomb spreading model.

---

## 11. Main Loop Analysis Results Summary

### Late/adult timepoint (250402)
- ~39,344 nonredundant merged loops across 3 resolutions (5kb: 7,901; 10kb: 14,553; 25kb: 16,890)
- 4,253 lost, 3,728 gained (FDR<0.05, |logFC|>0.3)
- Lost loops median 625 kb, gained median 320 kb (KS p = 5.49e-48)
- Loop type: most lost are Active_Enhancer-Active_Enhancer; most gained are CTCF-CTCF, Polycomb-Polycomb
- K119ub at anchor: 10x more likely to become disrupted

### Early/P12 timepoint (250831)
- 195 lost, 166 gained
- Progressive effect: expands from ~200 differential loops early to ~3,000 in adult

### Key data files for clustering
- `outputs/250402-late_outputs/bedpe_final/merged_all_loops_nonredundant.bedpe` -- 39,344 loops
- Per-resolution count matrices: `outputs/250402-late_outputs/res_{5,10,25}kb/06_counts_matrix.tsv`
- edgeR results: `outputs/250402-late_outputs/edgeR_results_res_{RES}kb/primary_analysis/all_results_primary.tsv`
- Characterized loops: `outputs/250402-late_outputs/merged_loops/characterized_loops.tsv`

---

## 12. Available Data Files

### ChIP-seq peak BEDs (for ChromHMM and annotation)

| Mark | Late file | Peaks |
|------|-----------|-------|
| H3K27ac | `peaks/beds/H3K27acCerebellumLate2.bed` | 15,105 |
| H3K27me3 | `peaks/beds/H3K27me3CerebellumLate1.bed` | 15,809 |
| H3K4me1 | `peaks/beds/H3K4me1CerebellumLate1.bed` | 113,781 |
| H3K4me3 | `peaks/beds/H3K4me3CerebellumLate2.bed` | 6,581 |
| CTCF | `peaks/CTCF.bed` | 32,487 |

### BigWig files (for deepTools metagene)

**On disk at `peaks/bigwigs/macs2.narrow.aug18.dedup/`:**
- ctrl H3K27ac: `index_13_ctrl_1_H3K27ac_S25_L001_aligned_reads.sorted.bw`
- ctrl H3K27me3: `index_25_ctrl_1_H3K27me3_S37_L001_aligned_reads.sorted.bw`
- mut H3K27ac: `index_19_mut_1_H3K27ac_S46_L002_R2_001.fastq.gz_aligned_reads.sorted.bw`
- mut H3K27me3: `index_23_mut_1_H3K27me3_S50_L002_R1_001_aligned_reads.sorted.bw`

**Present locally on Mac (`/Users/zakiralibhai/sdsc/bigwigs/`):**
- H2AK119ubCtrl.bw, H2AK119ubMut.bw
- H3K27me1Ctrl.bw, H3K27me1Mut.bw

### DiffBind results (for cluster x differential ChIP)
- `peaks/diffbind/K27ac_diffbind_results_summit_appended_ap.txt`
- `peaks/diffbind/K27me3_diffbind_results_summit_appended_ap.txt`
- `peaks/diffbind/K119ub_diffbind_results_summit_appended_ap.txt`

Column mapping vs Popay: our files use `Peak_Chr/Peak_Start/Peak_End` (cols 4-6) and `Fold` (col 12), `FDR` (col 14); Popay uses `Chr/Start/End/Fold/sig`.

### Other data
- Superenhancers: `peaks/Superenhancers_P60.bed` (1,046 SEs), `peaks/Superenhancers_encode.bed` (52 SEs)
- DEG lists: `tads/adult_timepoint_rna-seq-BAP1_WT_KO_v2_Results.xlsx`, `tads/young_timepoint_rna-seq-Bap1Math1paired_ctrl_mut_Results.xlsx`
- CTCF motifs: `peaks/ctcf_motifs_mm10.bed`
- Blacklist: `tads/mm10-blacklist.v2.bed`
- Compartment PC1: `tads/tad-pc-analysis/`

---

## 13. Completed Analyses (Don't Redo)

| Analysis | Status | Output location |
|----------|--------|-----------------|
| PC1 scatter plot (not volcano) | DONE | `output/compartment_analysis/compartment_pc1_scatter/` |
| TAD-K119ub correlation | DONE | `output/tad_k119ub_analysis/` |
| ABC cleanup (filtered scatter) | DONE | `abc/results/figures/k119ub_filtered_scatter/` |
| Activity vs Contact decomposition | DONE | `abc/results/figures/activity_contact_scatter/` |
| Stripenn full pipeline (Stages 0-7) | DONE | `stripes/stripenn/outputs/` |
| Developmental comparison | DONE | `output/developmental_comparison/` |

---

## 14. TODO Priority for This Clustering Analysis

From `TODO.md` Section 4, this is priority #1 for the paper:

1. **Popay anchor-vs-span analysis** -- Tessa just shared code; highest scientific value for distinguishing extrusion impediment vs. sensitivity model
2. **CTCF motif re-centering** (Section 2 of TODO) -- refines existing annotation
3. **SE-DEG proximity** (Section 1a) -- direct meeting request
4. **DEG-centric enhancer contacts** (Section 3) -- fills gap in ABC analysis

The clustering analysis (this work) addresses #1 directly.

---

## 15. Environment & Tools

### Conda environment
Created from `cluster/cluster.yml`. Key packages: pandas, numpy, seaborn, matplotlib, scipy, scikit-learn, bioframe, pyBigWig, pybedtools, cooltools, deeptools.

### Cluster 3.0
Installed and available on PATH as `cluster`. Location: `~/apps/cluster-1.59/src/cluster`.

### ChromHMM v1.27
JAR at `cluster/ChromHMM/ChromHMM.jar`. Wrapper on PATH as `chromhmm`. CHROMSIZES for mm10 included. Mac has 24GB RAM (relevant for `-mx8G` allocation).

### System
- macOS, 24GB RAM
- Java required for ChromHMM (verified working)
- Python 3.x in conda env
- R available for gene annotation script (TxDb.Mmusculus.UCSC.mm10.knownGene, org.Mm.eg.db installed in mariner_env)

---

## 16. Key Adaptations from Popay to BAP1-KO

| Aspect | Popay (hTERT RPE-1) | Our adaptation (BAP1-KO cerebellum) |
|--------|---------------------|-------------------------------------|
| Genome | hg38 | mm10 |
| Conditions | DMSO, dTAG 4h, dTAG 24h | ctrl_merged, mut_merged |
| Perturbation | Acute NIPBL depletion (hours) | Genetic BAP1 KO (developmental) |
| Loop source | RAD21 ChIP-based, 16,860 loops | HiCCUPS + mariner, 39,344 merged loops |
| Count normalization | Normalize to DMSO_merge column | Normalize to ctrl_merged column |
| Cluster count | k=6 (elbow-determined) | Start with k=6, adjust by elbow |
| ChromHMM states | 12-state from RPE-1 | 12-state learned from our 5 marks |
| ChIP normalization | Normalized to RAD21 RPKM | Raw RPKM (no RAD21 data) |
| Loop classification | CTCF+RAD21 = structural | CTCF alone = structural (no RAD21) |
| Gene annotation | gencode v25 (hg38) | TxDb mm10 knownGene |
| DiffBind marks | RAD21 only | K27ac, K27me3, K119ub |
| Pairwise comparison | comparison_type returns 'pairwise' | Returns 'multiple' (no `_` in col names) |

---

## 17. Scope and Focus

**This clustering analysis focuses on the late/adult timepoint (250402) only.** Rationale: strong differential signal (4,253 lost + 3,728 gained loops), clear directional patterns, most ChIP-seq data available for this timepoint. Early timepoint has too few differential loops (~200) for meaningful clustering.

**Input:** All 39,344 nonredundant merged loops (multi-resolution). Replicates averaged to 2 columns (ctrl_merged, mut_merged).

**Target output:** ChromHMM anchor vs span heatmaps (Fig 2f equivalent) showing whether Polycomb enrichment is at anchors or across loop body, broken down by loop clusters with varying sensitivity to BAP1 loss. Plus: loop size distributions, loop classification, ChIP signal, DiffBind relationship, gene annotation, and deepTools metagene plots per cluster.
