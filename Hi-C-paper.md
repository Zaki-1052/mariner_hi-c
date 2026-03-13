Figure 1: Progressive intensification of 3D genome phenotype \- TADs, boundaries, representative compartments

- Contact map of similar site in P13 vs adult where changes are present in both, but more subtle at P13, differential insulation track below each map  
- Volcano plot \- Differential TAD IDs at P13 and adult  
- GENOVA insulation plot \- Differential TAD boundaries  
- % of genome with differential PC1  
  - Smaller contact map of obvious differential region  
- Representative locus for differential TAD, triangular (probably Cntnap5a)  
- Permutation test of differential TADs and differential boundaries in K27ac/K27me3/K119ub

Figure 2: Chromatin loop rewiring \- differential loop analysis, APA, annotation w/ ctrl peaks (might split into 2 figures down the road)

- Volcano plots \- P13 and adult differential loop  
- Lost vs gained loop eCDF at both P13 and adult  
- APA analysis at both timepoints  
- Example loci of gained and lost loops  
- Lost loop length distribution annotation for K27me3 & K27ac → lost loops generally divided b/tw active enhancer loops and long-range polycomb loops  
- Genomic anchor distribution (combined anchor 1&2)  
- Bar plot \- Number of loops by category  
- Loop strength changes by category  
- Lost vs gained loop eCDF for K27ac-anchored, K27me3-anchored, CTCF-CTCF loops

Figure 3: Integration with epigenetic changes

- ARA for ATAC, H3K27me3, H3K27ac, K119ub  
- Histone marks at differential TAD boundaries  
- Permutation tests \- adult differential peaks in adult structures → changes to the B compartment are largely indirect of epigenetic changes, which are occurring largely in the A compartment  
- Permutation tests \- P12 marks predict adult structures  
  - P12 and adult loops against adult loops/boundaries/TADs/compartments  
- H2Az variant dynamics at structural features (if published or on biorxiv by summer)  
- One integrated locus \- probably Syt1

Figure 4: Enhancer ABC analysis

- ABC-RNA correlation using unnormalized score  
- Concordance analysis for 957 DEGs  
  - Pie chart of concordant vs discordant and of the 4 categories  
  - Characteristics of discordant in supplemental  
- Volcano plots for K27ac, ATAC, K119ub within combined ATAC peaks  
- Model for subsetting enhancers  
- Loop logFC, delta activity, delta contacts, gene logFC  
- ABC category vs K119ub significance

Figure 5: Model and functional implications

- GO analysis \- diff boundaries, loops, AxC connections  
- DEG logFC for genes near diff boundaries, at diff loops, with lost AxC connections  
- Network analysis  
  - Nodes \- genes with combined structural changes (boundary, loop, ABC)  
  - Node size \- AxC score change  
  - Node color \- Gene expression logFC  
- Heatmap of top 50 genes by combined structural score  
  - Columns \- logFC genes, logFC AxC  
- 3x2 model of epigenetic layer, structural layer, and functional (DEG) layer at P13 and adult