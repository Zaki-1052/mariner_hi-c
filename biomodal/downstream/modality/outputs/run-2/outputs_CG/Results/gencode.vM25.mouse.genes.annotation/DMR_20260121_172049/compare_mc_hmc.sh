#!/bin/bash
# compare_mc_hmc.sh

MC_FILE="DMR_mc_control__mutant_20260121_172049.bed"
HMC_FILE="DMR_hmc_control__mutant_20260121_172049.bed"
OUTPUT="mc_hmc_stats.txt"

{
echo "=== FILE SUMMARY ==="
echo "mC file: $(wc -l < "$MC_FILE") total entries"
echo "hmC file: $(wc -l < "$HMC_FILE") total entries"
echo ""

echo "=== SIGNIFICANT DMRs (q < 0.05) ==="
echo "mC significant: $(awk -F'\t' 'NR>1 && $12 < 0.05' "$MC_FILE" | wc -l)"
echo "hmC significant: $(awk -F'\t' 'NR>1 && $12 < 0.05' "$HMC_FILE" | wc -l)"
echo ""

echo "=== DIRECTION OF CHANGE ==="
echo "mC hyper (mutant > control): $(awk -F'\t' 'NR>1 && $12 < 0.05 && $9 > 0' "$MC_FILE" | wc -l)"
echo "mC hypo (mutant < control): $(awk -F'\t' 'NR>1 && $12 < 0.05 && $9 < 0' "$MC_FILE" | wc -l)"
echo "hmC hyper (mutant > control): $(awk -F'\t' 'NR>1 && $12 < 0.05 && $9 > 0' "$HMC_FILE" | wc -l)"
echo "hmC hypo (mutant < control): $(awk -F'\t' 'NR>1 && $12 < 0.05 && $9 < 0' "$HMC_FILE" | wc -l)"
echo ""

echo "=== GENE OVERLAP (q < 0.05) ==="
awk -F'\t' 'NR>1 && $12 < 0.05 {print $14}' "$MC_FILE" | sort -u > /tmp/mc_sig_genes.txt
awk -F'\t' 'NR>1 && $12 < 0.05 {print $14}' "$HMC_FILE" | sort -u > /tmp/hmc_sig_genes.txt
echo "Unique significant genes in mC: $(wc -l < /tmp/mc_sig_genes.txt)"
echo "Unique significant genes in hmC: $(wc -l < /tmp/hmc_sig_genes.txt)"
echo "Genes significant in BOTH: $(comm -12 /tmp/mc_sig_genes.txt /tmp/hmc_sig_genes.txt | wc -l)"
echo "Genes significant in mC ONLY: $(comm -23 /tmp/mc_sig_genes.txt /tmp/hmc_sig_genes.txt | wc -l)"
echo "Genes significant in hmC ONLY: $(comm -13 /tmp/mc_sig_genes.txt /tmp/hmc_sig_genes.txt | wc -l)"
echo ""

echo "=== DIRECTION COMPARISON (q < 0.05 in BOTH) ==="
awk -F'\t' 'NR>1 && $12 < 0.05 {
    if ($9 > 0) dir="hyper"; else dir="hypo"
    print $14"\t"dir
}' "$MC_FILE" | sort -u > /tmp/mc_directions.txt

awk -F'\t' 'NR>1 && $12 < 0.05 {
    if ($9 > 0) dir="hyper"; else dir="hypo"
    print $14"\t"dir
}' "$HMC_FILE" | sort -u > /tmp/hmc_directions.txt

join -t'	' /tmp/mc_directions.txt /tmp/hmc_directions.txt | \
awk -F'\t' '
BEGIN {same=0; opposite=0; mc_hyper_hmc_hypo=0; mc_hypo_hmc_hyper=0}
{
    if ($2 == $3) same++
    else {
        opposite++
        if ($2 == "hyper" && $3 == "hypo") mc_hyper_hmc_hypo++
        if ($2 == "hypo" && $3 == "hyper") mc_hypo_hmc_hyper++
    }
}
END {
    print "Same direction: " same
    print "Opposite direction: " opposite
    print "  mC hyper + hmC hypo: " mc_hyper_hmc_hypo
    print "  mC hypo + hmC hyper: " mc_hypo_hmc_hyper
}'
echo ""

echo "=== SUMMARY STATISTICS ==="
echo "mC difference (mutant - control):"
awk -F'\t' 'NR>1 {sum+=$9; sumsq+=$9*$9; n++} END {
    mean=sum/n; sd=sqrt(sumsq/n - mean*mean)
    printf "  All genes: mean=%.4f, sd=%.4f\n", mean, sd
}' "$MC_FILE"
awk -F'\t' 'NR>1 && $12 < 0.05 {sum+=$9; sumsq+=$9*$9; n++} END {
    mean=sum/n; sd=sqrt(sumsq/n - mean*mean)
    printf "  Significant: mean=%.4f, sd=%.4f\n", mean, sd
}' "$MC_FILE"

echo "hmC difference (mutant - control):"
awk -F'\t' 'NR>1 {sum+=$9; sumsq+=$9*$9; n++} END {
    mean=sum/n; sd=sqrt(sumsq/n - mean*mean)
    printf "  All genes: mean=%.4f, sd=%.4f\n", mean, sd
}' "$HMC_FILE"
awk -F'\t' 'NR>1 && $12 < 0.05 {sum+=$9; sumsq+=$9*$9; n++} END {
    mean=sum/n; sd=sqrt(sumsq/n - mean*mean)
    printf "  Significant: mean=%.4f, sd=%.4f\n", mean, sd
}' "$HMC_FILE"
echo ""

echo "=== PEARSON CORRELATION (mC_diff vs hmC_diff) ==="
awk -F'\t' 'NR>1 {print $14"\t"$9}' "$MC_FILE" | sort -k1,1 > /tmp/mc_diff.txt
awk -F'\t' 'NR>1 {print $14"\t"$9}' "$HMC_FILE" | sort -k1,1 > /tmp/hmc_diff.txt

join -t'	' /tmp/mc_diff.txt /tmp/hmc_diff.txt | awk -F'\t' '{
    n++; sum_x+=$2; sum_y+=$3; sum_xy+=$2*$3; sum_x2+=$2*$2; sum_y2+=$3*$3
} END {
    mean_x=sum_x/n; mean_y=sum_y/n
    cov=(sum_xy/n)-(mean_x*mean_y)
    sd_x=sqrt((sum_x2/n)-(mean_x*mean_x))
    sd_y=sqrt((sum_y2/n)-(mean_y*mean_y))
    r=cov/(sd_x*sd_y)
    printf "All genes (n=%d): r = %.4f\n", n, r
}'

awk -F'\t' 'NR>1 && $12 < 0.05 {print $14"\t"$9"\t"$12}' "$MC_FILE" | sort -k1,1 > /tmp/mc_sig.txt
awk -F'\t' 'NR>1 && $12 < 0.05 {print $14"\t"$9"\t"$12}' "$HMC_FILE" | sort -k1,1 > /tmp/hmc_sig.txt

join -t'	' /tmp/mc_sig.txt /tmp/hmc_sig.txt | awk -F'\t' '{
    n++; sum_x+=$2; sum_y+=$4; sum_xy+=$2*$4; sum_x2+=$2*$2; sum_y2+=$4*$4
} END {
    mean_x=sum_x/n; mean_y=sum_y/n
    cov=(sum_xy/n)-(mean_x*mean_y)
    sd_x=sqrt((sum_x2/n)-(mean_x*mean_x))
    sd_y=sqrt((sum_y2/n)-(mean_y*mean_y))
    r=cov/(sd_x*sd_y)
    printf "Significant in BOTH (n=%d): r = %.4f\n", n, r
}'
echo ""

echo "=== TOP 15 GENES: mC hyper + hmC hypo ==="
join -t'	' /tmp/mc_sig.txt /tmp/hmc_sig.txt | \
awk -F'\t' '$2 > 0 && $4 < 0' | sort -t'	' -k2 -rn | head -15 | \
awk -F'\t' '{printf "%-15s mC_diff: %+.4f (q=%.2e)  hmC_diff: %+.4f (q=%.2e)\n", $1, $2, $3, $4, $5}'
echo ""

echo "=== TOP 15 GENES: mC hypo + hmC hyper ==="
join -t'	' /tmp/mc_sig.txt /tmp/hmc_sig.txt | \
awk -F'\t' '$2 < 0 && $4 > 0' | sort -t'	' -k2 -n | head -15 | \
awk -F'\t' '{printf "%-15s mC_diff: %+.4f (q=%.2e)  hmC_diff: %+.4f (q=%.2e)\n", $1, $2, $3, $4, $5}'

} > "$OUTPUT"

echo "Results written to $OUTPUT"
