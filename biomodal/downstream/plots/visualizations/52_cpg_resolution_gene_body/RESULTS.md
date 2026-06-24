## Section 52: Single-CpG-Resolution Gene Body Analysis
**Key numbers:**
- 6,196 sub-gene feature intervals across 123 genes (1,742 Intron, 1,567 Exon, 1,341 SpliceSite_Acceptor, 1,330 SpliceSite_Donor, 124 3UTR, 92 5UTR) (section52_feature_methylation.tsv)
- delta-5mC: 3/15 feature pairs significant (q<0.05), all vs Intron — Exon-Intron Z=-3.864 q=4.18e-04; SS_Donor-Intron Z=-3.898 q=7.26e-04; SS_Acceptor-Intron Z=-3.424 q=1.55e-03 (section52_dunn_posthoc.tsv)
- delta-5hmC: 8/15 pairs significant — incl 5UTR-Intron Z=3.845 q=9.06e-04, SS_Acceptor-Intron q=6.40e-04, 5UTR-3UTR q=2.44e-03 (section52_dunn_posthoc.tsv)

**What this shows:** The BAP1-KO CG methylation response is structured by sub-gene feature: introns differ from exons and splice sites for both 5mC and 5hmC, and 5hmC additionally separates the 5' end (5'UTR) from the rest of the body. 5hmC carries more spatially-resolved signal (8 vs 3 significant pairs), consistent with the dataset's TET-block thesis.

**Figures:**
- `52a_feature_distribution/` — feature-interval counts per type
- `52b_delta_mc_by_feature/` / `52c_delta_hmc_by_feature/` — delta-methylation boxplots by feature
- `52d_key_gene_loci/` — per-KEY_GENE lollipop panels
- `52e_metagene_profile/` — 5'-to-3' metagene delta-5mC/5hmC
- `52f_dunn_posthoc/` — Dunn post-hoc q-value heatmaps
- `52g_chromatin_overlay/` — ChIP-mark overlap and H3K27ac splits per feature
