#!/bin/bash
# significant_to_bigbed.sh
# Convert significant DMR BED files to BigBed format

echo "========================================="
echo "Significant DMR BED to BigBed Conversion"
echo "========================================="
echo "Start time: $(date)"
echo ""

# Output directory
BIGBED_DIR="bigbed"
mkdir -p "${BIGBED_DIR}"

# Get mm10 chromosome sizes if not present
CHROM_SIZES="mm10.chrom.sizes"
if [ ! -f "${CHROM_SIZES}" ]; then
    echo "Downloading mm10 chromosome sizes..."
    curl -s "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes" > "${CHROM_SIZES}"
fi

# Create AutoSql file describing DMR BED format
cat > "${BIGBED_DIR}/dmr.as" << 'EOF'
table dmr
"Differentially Methylated Region from modality XPLR"
    (
    string chrom;        "Chromosome"
    uint chromStart;     "Start position"
    uint chromEnd;       "End position"
    string name;         "Gene name"
    uint score;          "Score (scaled -log10 qvalue, 0-1000)"
    char[1] strand;      "Strand"
    float modDiff;       "Methylation difference (mutant - control)"
    float qvalue;        "Adjusted p-value (FDR)"
    float foldChange;    "Fold change"
    )
EOF

# List of BED files to convert
BED_FILES=(
    "mc_significant_q0.05.bed"
    "mc_significant_hyper.bed"
    "mc_significant_hypo.bed"
    "hmc_significant_q0.05.bed"
    "hmc_significant_hyper.bed"
    "hmc_significant_hypo.bed"
)

echo "Converting ${#BED_FILES[@]} BED files..."
echo ""

for BED in "${BED_FILES[@]}"; do
    if [ ! -f "${BED}" ]; then
        echo "WARNING: ${BED} not found, skipping"
        continue
    fi

    BASENAME=$(basename "${BED}" .bed)
    TMPBED="${BIGBED_DIR}/${BASENAME}.tmp.bed"
    BIGBED="${BIGBED_DIR}/${BASENAME}.bb"

    echo "Processing: ${BED}"

    # Convert DMR BED to BigBed-compatible format
    # Input columns: Chromosome Start End num_contexts mean_coverage mean_mod_group_1 mean_mod_group_2 mod_fold_change mod_difference test_statistic dmr_pvalue dmr_qvalue Annotation Name
    # Output columns: chr, start, end, name, score, strand, modDiff, qvalue, foldChange
    tail -n +2 "${BED}" | \
        awk -F'\t' '
        $9 != "" && $9 != "nan" && $12 != "" && $12 != "nan" {
            # Scale -log10(qvalue) to 0-1000 score
            qval = $12 + 0
            if (qval <= 0) qval = 1e-10
            score = int(-log(qval)/log(10) * 100)
            if (score > 1000) score = 1000
            if (score < 0) score = 0

            # Handle missing/invalid fold change
            fc = ($8 == "" || $8 == "nan" || $8 == "inf" || $8 == "-inf") ? 1 : $8

            # Use gene name ($14) as name field
            print "chr"$1, $2, $3, $14, score, ".", $9, $12, fc
        }' OFS='\t' | \
        sort -k1,1 -k2,2n > "${TMPBED}"

    # Check if temp BED has content
    if [ -s "${TMPBED}" ]; then
        # Try with AutoSql first
        bedToBigBed -type=bed6+3 -as="${BIGBED_DIR}/dmr.as" "${TMPBED}" "${CHROM_SIZES}" "${BIGBED}" 2>/dev/null
        if [ $? -eq 0 ]; then
            echo "  -> ${BIGBED} ($(wc -l < "${TMPBED}") entries)"
        else
            # Fallback: try without AutoSql (simple BED6)
            awk -F'\t' '{print $1,$2,$3,$4,$5,$6}' OFS='\t' "${TMPBED}" > "${TMPBED}.bed6"
            bedToBigBed "${TMPBED}.bed6" "${CHROM_SIZES}" "${BIGBED}"
            rm -f "${TMPBED}.bed6"
            echo "  -> ${BIGBED} (BED6 fallback, $(wc -l < "${TMPBED}") entries)"
        fi
        rm -f "${TMPBED}"
    else
        echo "  WARNING: Empty after filtering, skipping"
        rm -f "${TMPBED}"
    fi
done

echo ""
echo "========================================="
echo "Conversion complete"
echo "BigBed files saved to: ${BIGBED_DIR}/"
ls -lh "${BIGBED_DIR}"/*.bb 2>/dev/null || echo "No BigBed files created"
echo "End time: $(date)"
echo "========================================="
