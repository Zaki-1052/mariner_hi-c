# Differential Loop Analysis - Methods and Visualizations Index

## Analysis Pipeline Overview

**Normalization:** Faster normalization method from multiHiCompare

**Differential Analysis:** multiHiCompare on Hi-C data comparing WT vs KO conditions

**Annotation:** chIPseeker (R) to assign each loop anchor to genes and genomic regions

**Enrichment Analysis:** clusterProfiler (R) for Gene Ontology analysis

**Visualization:** EnhancedVolcano package (R) for volcano plots

---

## Visualization Types Present

### 1. Volcano Plots (EnhancedVolcano)
**Input:** _primary.tsv files from differential analysis

**Resolutions analyzed:**
- 5kb resolution (17,982 total variables; 458 up, 903 down)
- 10kb resolution (22,632 total variables; 640 up, 988 down)  
- 25kb resolution (20,398 total variables; 668 up, 556 down)
- Merged multi-resolution (61,012 total variables; 1,449 up, 2,392 down)

**Axes:**
- X-axis: LogFC Pixel Enrichment
- Y-axis: FDR (False Discovery Rate)

**Color coding:**
- Black: Non-significant (NS)
- Grey: Log₂FC only
- Red: p-value only
- Dark red: p-value and log₂FC

---

### 2. Feature Distribution Analysis

**Method:** chIPseeker annotation of loop anchors

**Format:** Stacked horizontal bar chart showing percentage distribution

**Categories analyzed:** up1, up2, down1, down2

**Genomic features annotated:**
- Promoter
- 5' UTR
- 3' UTR
- 1st Exon
- Other Exon
- 1st Intron
- Other Intron
- Downstream (<=300bp)
- Distal Intergenic

---

### 3. Gene Ontology Enrichment Dot Plots

**Method:** clusterProfiler

**Ontologies analyzed:**
- KEGG pathways
- GO Biological Process (GO BP)
- GO Cellular Component (GO CC)
- GO Molecular Function (GO MF)

**Format:** Dot plot matrices

**Visual encoding:**
- Dot color: p.adjust value (gradient scale)
- Dot size: GeneRatio
- Rows: Enriched terms/pathways
- Columns: Clusters (up1, up2, down1, down2)
- Sample counts shown below each cluster

---

### 4. Loop Type Classification Pie Charts

**Categories:** Two separate pie charts for "Up loops" and "Down loops"

**Loop type classifications:**
- promoter-promoter
- enhancer-promoter
- enhancer-enhancer
- distal-promoter
- distal-enhancer
- distal-distal

**Additional info:** Lists of genes contained within each loop set, including:
- Up loops: Smad2/5/6, Myc, Lhx9, Zbtb20, Wnt5a, Sox2/6/9, Slc7a11, Shh, Ptprt, Nav1/3, Mki67, Lsamp, Hox locus, Foxg1, Cntnap5a, Auts2, Apc, Epha5, Tgfbr1, Satb2, Pvt1, Prune2, Nrcam, Map2, En1, Ep300, Dscam, Dcx, Ano1
- Down loops: Aebp2, Auts2, Bdnf, Bmp5, C9orf72, Cbx3, Cgas, Crbn, Dgkb, Dlg2, Eed, Hox locus, Lsamp, Mef2c, Nfib, Nkx2-3, Patj, Pax6, Phc3, Prkn, Rit2, Zdbf2, Syt1, Nav3, Snca, Zic1/2/4, Prune2, Zbtb20, Lmnb1, Eomes, Calb1

---

### 5. Aggregate Peak Analysis (APA)

**Conditions compared:** HI-C ctrl vs BAP1-KO

**Panels shown for both up and down loops:**

1. **Contact frequency heatmaps (2x2 grid):**
   - Top left: KO condition (µ Contacts scale)
   - Top right: WT condition (µ Contacts scale)
   - Bottom left: Individual KO
   - Bottom right: Individual WT
   - Center bottom: Difference heatmap (Difference scale: -5.0 to +5.0)

2. **Pixel enrichment box plots:**
   - Comparison between KO (red) and WT (grey)
   - Y-axis: pixel enrichment_final
   - Shows distribution with median, quartiles, and outliers

**Spatial features:** 50kb windows centered on loops with annotation at -50kb, 3', and 50kb positions

---

### 6. Loop Length Distribution

**Format:** Strip plot (jitter/dot plot)

**Categories:** 
- Up (blue dots)
- Down (red dots)

**X-axis:** Loop length (0 to 100k+ bp range)

**Observation:** Majority of loops cluster at short distances (<20kb), with sparse longer-range loops extending to 100kb+

---

### 7. Genome Browser Views - Chromosome-wide Hi-C Maps

**Example: Chromosome 1**

**Layout:**
- Mutant (KO) on bottom left triangle
- Control (WT) on top right triangle
- RefSeq genes track at top

**Loop annotations overlaid:**
- Blue dots: Downregulated loops
- Black dots: Upregulated loops

**Resolution:** Full chromosome view (0 MB to ~180 MB)

**Note:** Shows patches of spatially clustered downregulated loops

---

### 8. Locus-Specific Hi-C Contact Maps

**Loci examined:**
- Hox locus
- Syt1/Nav3 locus

**Format:** 
- High-resolution contact maps (heatmaps)
- Side-by-side comparison panels (typically 2-3 panels per locus)
- KO condition on left, WT on right (or comparisons)

**Annotations:**
- Blue dots/squares: Downregulated loops
- Black dots/squares: Upregulated loops

**Color scale:** Red intensity indicating contact frequency (µ Contacts)

---

## Data Processing Notes

- Normalization performed using faster method from multiHiCompare
- Analysis run on _primary.tsv output files
- Gene assignment and region annotation via chIPseeker prior to GO analysis
- Clustering performed to separate up/down regulated loops into subclusters (up1, up2, down1, down2)