## Section 67: MeCP2 at Unmethylated Genes Tracks H2AK119ub (KEY RESULT)

**Key numbers:**
- 359 candidate MeCP2-up-no-methylation genes; 356 with finite K119ub log2FC analyzed; background n=22,027 (67_statistics.tsv)
- K119ub gained in 72.8% of our genes vs 45.9% genome-wide: Fisher OR=3.15, p=1.82e-24 (67_statistics.tsv)
- K119ub log2FC our median +0.202 vs rest -0.036: Mann-Whitney p=5.50e-32 (67_statistics.tsv)
- Multi-mark at these loci: K119ub +0.202 (gained), K27me3 +0.162 (gained), K27ac -0.242 (LOST), ATAC +0.074 (67_per_mark_summary.tsv)

**What this shows:** At genes where MeCP2 binding increases without any methylation change, the dominant concurrent epigenomic event is H2AK119ub gain (73% of genes, 3.15x over genome). MeCP2 fold positively correlates with K119ub change (Spearman rho=0.241, p=4.4e-06). Because BAP1 is the K119ub deubiquitinase, its loss raises K119ub and MeCP2 follows that Polycomb signal — decoupled from DNA methylation. Coincident K27me3 gain and K27ac loss confirm these are Polycomb-repressed loci. This is the mechanistic bridge to the Hi-C Polycomb-compaction phenotype.

**Figures:**
- 67a_k119ub_at_mecp2_no_meth/ — K119ub log2FC violin: our genes vs all others
- 67b_k119ub_gained_fraction/ — % K119ub gained bar (Fisher OR/p)
- 67c_mecp2_vs_k119ub_scatter/ — MeCP2 fold vs K119ub log2FC (Spearman)
- 67d_multimark_at_mecp2_no_meth/ — K119ub/K27me3/K27ac/ATAC violins (vs 0)
