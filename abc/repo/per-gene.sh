echo "=== Per-gene ABC.Score sums (WT, single-pass) ===" && \
zcat results/WT_cerebellum/Predictions/EnhancerPredictionsAllPutative.tsv.gz \
  | awk -F'\t' 'NR>1{sum[$10]+=$26; n[$10]++} END{
      # Print 10 random-ish genes spanning different sum ranges
      i=0;
      for(g in sum){
        i++;
        if(i%4800==1 || i%4800==1000 || i%4800==2400 || i%4800==3600 || i%4800==4000){
          printf "  %s: sum(ABC)=%.6f from %d enhancers\n", g, sum[g], n[g]
        }
      }
      # Also compute global statistics on per-gene sums
      n_genes=0; sum_sums=0; n_bad=0;
      for(g in sum){
        n_genes++;
        sum_sums+=sum[g];
        if(sum[g] < 0.95 || sum[g] > 1.05) n_bad++;
      }
      printf "\n  Total genes: %d\n", n_genes;
      printf "  Mean per-gene sum: %.6f\n", sum_sums/n_genes;
      printf "  Genes with sum outside [0.95, 1.05]: %d (%.1f%%)\n", n_bad, 100*n_bad/n_genes;
    }'

echo "=== Per-gene ABC.Score sums (KO, single-pass) ===" && \
zcat results/KO_cerebellum/Predictions/EnhancerPredictionsAllPutative.tsv.gz \
  | awk -F'\t' 'NR>1{sum[$10]+=$26; n[$10]++} END{
      i=0;
      for(g in sum){
        i++;
        if(i%4800==1 || i%4800==1000 || i%4800==2400 || i%4800==3600 || i%4800==4000){
          printf "  %s: sum(ABC)=%.6f from %d enhancers\n", g, sum[g], n[g]
        }
      }
      n_genes=0; sum_sums=0; n_bad=0;
      for(g in sum){
        n_genes++;
        sum_sums+=sum[g];
        if(sum[g] < 0.95 || sum[g] > 1.05) n_bad++;
      }
      printf "\n  Total genes: %d\n", n_genes;
      printf "  Mean per-gene sum: %.6f\n", sum_sums/n_genes;
      printf "  Genes with sum outside [0.95, 1.05]: %d (%.1f%%)\n", n_bad, 100*n_bad/n_genes;
    }'
