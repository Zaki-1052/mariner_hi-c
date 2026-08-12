## Section 46: Multi-Omic Genome Browser Locus Views (VISUAL ONLY — no quantification)
**Key numbers:**
- No TSV output. All values are rendered live onto Gviz tracks (per-track RNA-seq log2FC/padj annotations from the RNA-seq xlsx; mC/hmC from condition-averaged BigWigs).
- 10 KEY_GENES rendered (full + compact each): Syt1, Zbtb20, Trpm3, Epha3, Mcu, Cntnap2, Lpp, Dlgap1, Arhgap26, Cdh8 (_shared_config.R)
- Each view extends ±50 kb (EXTEND_BP = 50,000) around the gene body

**What this shows:** Per-locus integrative browser views stacking condition-averaged 5mC/5hmC, mC/hmC difference, RNA-seq coverage, H2AK119ub/H3K27me3/H3K4me3/H3K27ac/ATAC/H3K27me1/MeCP2 signal (MeCP2 is CUT&RUN), CpG island/shore/shelf annotation, CTCF peaks, and Hi-C loop arcs (lost/gained). They let a reader see the full BAP1 mechanism — gained 5mC, lost 5hmC, RNA change, H2AK119ub accumulation, and loop disruption — co-localized in cis at single loci (Syt1 is the showcase). This RESULTS.md lives in the alphabetically-first panel folder; the section has no plain base folder.

**Figures:**
- Syt1_locus (+_compact, +_poster, +_poster_v2, +_poster_v2_aligned) — Syt1 browser variants
- Zbtb20/Trpm3/Epha3/Mcu/Cntnap2/Lpp/Dlgap1/Arhgap26/Cdh8 _locus (+_compact) — per-gene full + compact views
- composite_syt1_panel — Syt1 browser + coordinated mC/hmC scatter + multi-omic z-score heatmap
