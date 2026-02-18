cd /expanse/lustre/projects/csd940/zalibhai/abc/ABC-Enhancer-Gene-Prediction

# ---- Diagnostic 1: Count isSelfPromoter entries per gene ----
# If hypothesis is correct, most genes should have exactly 1 isSelfPromoter=True
# But may also have a nearby "enhancer" class element at the same promoter location
echo "=== Diagnostic 1: isSelfPromoter counts per gene ==="
zcat results/WT_cerebellum/Predictions/EnhancerPredictionsAllPutative.tsv.gz \
  | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="isSelfPromoter")sp=i; if($i=="TargetGene")tg=i}; next}
    $sp=="True"{count[$tg]++}
    END{
      for(g in count) freq[count[g]]++;
      for(f in freq) printf "  Genes with %d self-promoter entries: %d\n", f, freq[f];
      total=0; for(g in count) total++;
      printf "  Total genes with any self-promoter: %d\n", total
    }'

# ---- Diagnostic 2: Per-gene sum WITH vs WITHOUT self-promoter entries ----
# If sums drop to ~1.0 when excluding isSelfPromoter=True, the duplicate is confirmed
echo ""
echo "=== Diagnostic 2: Sum with vs without self-promoter ==="
zcat results/WT_cerebellum/Predictions/EnhancerPredictionsAllPutative.tsv.gz \
  | awk -F'\t' 'NR==1{
      for(i=1;i<=NF;i++){
        if($i=="isSelfPromoter")sp=i;
        if($i=="TargetGene")tg=i;
        if($i=="ABC.Score")abc=i
      }; next}
    {
      sum_all[$tg]+=$abc; n_all[$tg]++;
      if($sp!="True"){sum_nosp[$tg]+=$abc; n_nosp[$tg]++}
    }
    END{
      ng=0; s_all=0; s_nosp=0;
      for(g in sum_all){
        ng++;
        s_all+=sum_all[g];
        s_nosp+=(g in sum_nosp ? sum_nosp[g] : 0)
      }
      printf "  Total genes: %d\n", ng;
      printf "  Mean per-gene sum (all entries): %.4f\n", s_all/ng;
      printf "  Mean per-gene sum (excluding isSelfPromoter=True): %.4f\n", s_nosp/ng;
    }'

# ---- Diagnostic 3: Inspect a specific gene's promoter-region entries ----
# Pick Slc2a1 (sum=1.60, so it has a moderate effect) and show
# all entries within 2kb of its TSS
echo ""
echo "=== Diagnostic 3: Slc2a1 entries near TSS ==="
echo "(columns: chr start end name class isSelfPromoter distance ABC.Score)"
zcat results/WT_cerebellum/Predictions/EnhancerPredictionsAllPutative.tsv.gz \
  | awk -F'\t' 'NR==1{
      for(i=1;i<=NF;i++){
        if($i=="TargetGene")tg=i;
        if($i=="isSelfPromoter")sp=i;
        if($i=="distance")dist=i;
        if($i=="ABC.Score")abc=i;
        if($i=="class")cl=i;
        if($i=="name")nm=i
      }; next}
    $tg=="Slc2a1" && $dist+0 < 2000 {
      printf "%s\t%s\t%s\t%s\t%s\tisSP=%s\tdist=%s\tABC=%.6f\n", $1,$2,$3,$nm,$cl,$sp,$dist,$abc
    }'
