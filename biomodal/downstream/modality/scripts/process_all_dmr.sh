#!/bin/bash
# biomodal/downstream/modality/scripts/process_all_dmr.sh
# Comprehensive DMR processing pipeline for all contexts and annotations

set -e

# =============================================================================
# CONFIGURATION
# =============================================================================

# Base directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="$BASE_DIR/outputs/run-5"
OUTPUT_DIR="$BASE_DIR/DMR_processed"
CHROM_SIZES="$OUTPUT_DIR/mm10.chrom.sizes"
BEDTOBIGBED="$SCRIPT_DIR/bedToBigBed"

# Contexts and annotations to process
CONTEXTS="CG CHG CHH"
ANNOTATIONS="genes promoters cpg_islands cpg_shores cpg_shelves tss_region"

# Significance threshold
Q_THRESHOLD=0.05

# Function to get annotation folder name
get_annotation_folder() {
    local annotation="$1"
    case "$annotation" in
        genes) echo "gencode.vM25.mouse.genes.annotation" ;;
        promoters) echo "gencode.vM25.mouse.promoters.annotation" ;;
        cpg_islands) echo "gencode.vM25.mouse.cpg_islands.annotation" ;;
        cpg_shores) echo "gencode.vM25.mouse.cpg_shores.annotation" ;;
        cpg_shelves) echo "gencode.vM25.mouse.cpg_shelves.annotation" ;;
        tss_region) echo "gencode.vM25.mouse.tss_region.annotation" ;;
        *) echo "" ;;
    esac
}

# =============================================================================
# FUNCTIONS
# =============================================================================

log_msg() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

find_dmr_file() {
    # Find the DMR BED file for a given context, annotation, and modification type
    local context="$1"
    local annotation="$2"
    local mod_type="$3"  # mc or hmc

    local annotation_folder
    annotation_folder=$(get_annotation_folder "$annotation")
    local results_dir="$INPUT_DIR/outputs_${context}/Results/${annotation_folder}"

    # Find the most recent DMR directory (excluding DMR_Report_* directories)
    local dmr_dir=$(find "$results_dir" -maxdepth 1 -type d -name "DMR_*" ! -name "DMR_Report_*" 2>/dev/null | sort -r | head -1)

    if [[ -z "$dmr_dir" ]]; then
        echo ""
        return
    fi

    # Find the BED file matching the modification type
    local bed_file=$(find "$dmr_dir" -name "DMR_${mod_type}_control__mutant_*.bed" 2>/dev/null | head -1)
    echo "$bed_file"
}

process_dmr_file() {
    # Process a single DMR file: filter significant, split by direction, create BigBed
    local input_file="$1"
    local output_dir="$2"
    local mod_type="$3"  # mc or hmc
    local context="$4"
    local annotation="$5"

    local bed_dir="$output_dir/bed"
    local bigbed_dir="$output_dir/bigbed"
    local stats_dir="$output_dir/stats"

    mkdir -p "$bed_dir" "$bigbed_dir" "$stats_dir"

    if [[ ! -f "$input_file" ]]; then
        log_msg "  WARNING: Input file not found: $input_file"
        return
    fi

    # Count total entries (excluding header)
    local total_count=$(tail -n +2 "$input_file" | wc -l | tr -d ' ')

    # Extract significant DMRs (q < 0.05)
    local sig_file="$bed_dir/${mod_type}_significant_q0.05.bed"
    awk -F'\t' -v q="$Q_THRESHOLD" 'NR==1 || $12 < q' "$input_file" > "$sig_file"
    local sig_count=$(tail -n +2 "$sig_file" | wc -l | tr -d ' ')

    # Extract hyper (mutant > control, mod_difference > 0)
    local hyper_file="$bed_dir/${mod_type}_significant_hyper.bed"
    awk -F'\t' -v q="$Q_THRESHOLD" 'NR==1 || ($12 < q && $9 > 0)' "$input_file" > "$hyper_file"
    local hyper_count=$(tail -n +2 "$hyper_file" | wc -l | tr -d ' ')

    # Extract hypo (mutant < control, mod_difference < 0)
    local hypo_file="$bed_dir/${mod_type}_significant_hypo.bed"
    awk -F'\t' -v q="$Q_THRESHOLD" 'NR==1 || ($12 < q && $9 < 0)' "$input_file" > "$hypo_file"
    local hypo_count=$(tail -n +2 "$hypo_file" | wc -l | tr -d ' ')

    log_msg "  ${mod_type}: total=$total_count, sig=$sig_count (hyper=$hyper_count, hypo=$hypo_count)" >&2

    # Create BigBed files (only if there are entries)
    for bed_file in "$sig_file" "$hyper_file" "$hypo_file"; do
        local base_name=$(basename "$bed_file" .bed)
        local bb_file="$bigbed_dir/${base_name}.bb"

        # Count entries (excluding header)
        local entry_count=$(tail -n +2 "$bed_file" | wc -l | tr -d ' ')

        if [[ "$entry_count" -gt 0 ]]; then
            # Create BED3 format for BigBed (chr, start, end) with chr prefix
            local tmp_bed="/tmp/${base_name}_bed3.bed"
            tail -n +2 "$bed_file" | awk -F'\t' '{
                chr = $1
                if (chr !~ /^chr/) chr = "chr" chr
                print chr"\t"$2"\t"$3
            }' | sort -k1,1 -k2,2n > "$tmp_bed"

            # Convert to BigBed
            if "$BEDTOBIGBED" "$tmp_bed" "$CHROM_SIZES" "$bb_file" 2>/dev/null; then
                :  # Success, do nothing
            else
                log_msg "  WARNING: BigBed conversion failed for $bed_file"
            fi
            rm -f "$tmp_bed"
        fi
    done

    # Return counts as tab-separated values
    echo -e "${total_count}\t${sig_count}\t${hyper_count}\t${hypo_count}"
}

generate_comparison_stats() {
    # Generate combined mc/hmc statistics for an annotation
    local output_dir="$1"
    local context="$2"
    local annotation="$3"

    local bed_dir="$output_dir/bed"
    local stats_dir="$output_dir/stats"

    local mc_sig="$bed_dir/mc_significant_q0.05.bed"
    local hmc_sig="$bed_dir/hmc_significant_q0.05.bed"

    # Skip if files don't exist or are empty
    if [[ ! -s "$mc_sig" ]] || [[ ! -s "$hmc_sig" ]]; then
        log_msg "  Skipping comparison stats - insufficient data"
        return
    fi

    local stats_file="$stats_dir/mc_hmc_stats.txt"
    local comparison_file="$stats_dir/mc_hmc_comparison.tsv"

    {
        echo "=== DMR Statistics: ${context} - ${annotation} ==="
        echo "Generated: $(date)"
        echo ""

        echo "=== SIGNIFICANT DMRs (q < $Q_THRESHOLD) ==="
        echo "mC significant: $(tail -n +2 "$mc_sig" | wc -l | tr -d ' ')"
        echo "hmC significant: $(tail -n +2 "$hmc_sig" | wc -l | tr -d ' ')"
        echo ""

        echo "=== DIRECTION OF CHANGE ==="
        echo "mC hyper (mutant > control): $(tail -n +2 "$bed_dir/mc_significant_hyper.bed" 2>/dev/null | wc -l | tr -d ' ')"
        echo "mC hypo (mutant < control): $(tail -n +2 "$bed_dir/mc_significant_hypo.bed" 2>/dev/null | wc -l | tr -d ' ')"
        echo "hmC hyper (mutant > control): $(tail -n +2 "$bed_dir/hmc_significant_hyper.bed" 2>/dev/null | wc -l | tr -d ' ')"
        echo "hmC hypo (mutant < control): $(tail -n +2 "$bed_dir/hmc_significant_hypo.bed" 2>/dev/null | wc -l | tr -d ' ')"
        echo ""

        # Gene/feature overlap analysis
        # Use Name (col 14) if it exists and is non-empty, otherwise use coordinates
        echo "=== FEATURE OVERLAP (q < $Q_THRESHOLD) ==="
        tail -n +2 "$mc_sig" | awk -F'\t' '{key = ($14 != "" ? $14 : $1":"$2"-"$3); print key}' | sort -u > /tmp/mc_sig_features.txt
        tail -n +2 "$hmc_sig" | awk -F'\t' '{key = ($14 != "" ? $14 : $1":"$2"-"$3); print key}' | sort -u > /tmp/hmc_sig_features.txt
        echo "Unique significant in mC: $(wc -l < /tmp/mc_sig_features.txt | tr -d ' ')"
        echo "Unique significant in hmC: $(wc -l < /tmp/hmc_sig_features.txt | tr -d ' ')"
        echo "Features significant in BOTH: $(comm -12 /tmp/mc_sig_features.txt /tmp/hmc_sig_features.txt | wc -l | tr -d ' ')"
        echo "Features significant in mC ONLY: $(comm -23 /tmp/mc_sig_features.txt /tmp/hmc_sig_features.txt | wc -l | tr -d ' ')"
        echo "Features significant in hmC ONLY: $(comm -13 /tmp/mc_sig_features.txt /tmp/hmc_sig_features.txt | wc -l | tr -d ' ')"
        echo ""

        # Direction comparison for features significant in both
        echo "=== DIRECTION COMPARISON (significant in BOTH) ==="
        tail -n +2 "$mc_sig" | awk -F'\t' '{key = ($14 != "" ? $14 : $1":"$2"-"$3); if ($9 > 0) dir="hyper"; else dir="hypo"; print key"\t"dir}' | sort -u > /tmp/mc_dir.txt
        tail -n +2 "$hmc_sig" | awk -F'\t' '{key = ($14 != "" ? $14 : $1":"$2"-"$3); if ($9 > 0) dir="hyper"; else dir="hypo"; print key"\t"dir}' | sort -u > /tmp/hmc_dir.txt

        join -t'	' /tmp/mc_dir.txt /tmp/hmc_dir.txt | awk -F'\t' '
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

        # Summary statistics
        echo "=== SUMMARY STATISTICS ==="
        echo "mC difference (mutant - control):"
        tail -n +2 "$mc_sig" | awk -F'\t' '{sum+=$9; sumsq+=$9*$9; n++} END {
            if (n > 0) {mean=sum/n; sd=sqrt(sumsq/n - mean*mean); printf "  Significant: mean=%.4f, sd=%.4f, n=%d\n", mean, sd, n}
        }'
        echo "hmC difference (mutant - control):"
        tail -n +2 "$hmc_sig" | awk -F'\t' '{sum+=$9; sumsq+=$9*$9; n++} END {
            if (n > 0) {mean=sum/n; sd=sqrt(sumsq/n - mean*mean); printf "  Significant: mean=%.4f, sd=%.4f, n=%d\n", mean, sd, n}
        }'

    } > "$stats_file"

    # Create merged comparison TSV (for R visualization)
    # Use Name if available, otherwise use coordinate key
    {
        echo -e "Name\tChromosome\tStart\tEnd\tmc_diff\tmc_qvalue\thmc_diff\thmc_qvalue\tAnnotation"
        tail -n +2 "$mc_sig" | awk -F'\t' '{
            key = ($14 != "" ? $14 : $1":"$2"-"$3)
            print key"\t"$1"\t"$2"\t"$3"\t"$9"\t"$12
        }' | sort -k1,1 > /tmp/mc_for_merge.txt
        tail -n +2 "$hmc_sig" | awk -F'\t' '{
            key = ($14 != "" ? $14 : $1":"$2"-"$3)
            print key"\t"$9"\t"$12
        }' | sort -k1,1 > /tmp/hmc_for_merge.txt
        join -t'	' /tmp/mc_for_merge.txt /tmp/hmc_for_merge.txt | awk -F'\t' -v ann="$annotation" '{
            print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"ann
        }'
    } > "$comparison_file"

    log_msg "  Stats written to: $stats_file"
}

# =============================================================================
# MAIN PIPELINE
# =============================================================================

log_msg "=============================================="
log_msg "DMR Processing Pipeline - All Contexts"
log_msg "=============================================="
log_msg "Input directory: $INPUT_DIR"
log_msg "Output directory: $OUTPUT_DIR"
log_msg ""

# Verify input directory exists
if [[ ! -d "$INPUT_DIR" ]]; then
    log_msg "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

# Verify bedToBigBed exists
if [[ ! -x "$BEDTOBIGBED" ]]; then
    log_msg "ERROR: bedToBigBed not found or not executable: $BEDTOBIGBED"
    exit 1
fi

# Create output directories
mkdir -p "$OUTPUT_DIR/summary"

# Initialize summary file
SUMMARY_FILE="$OUTPUT_DIR/summary/all_contexts_summary.tsv"
echo -e "Context\tAnnotation\tModType\tTotal\tSignificant\tHyper\tHypo" > "$SUMMARY_FILE"

# Process each context and annotation
for context in $CONTEXTS; do
    log_msg ""
    log_msg "=== Processing Context: $context ==="

    for annotation in $ANNOTATIONS; do
        log_msg ""
        log_msg "  Annotation: $annotation"

        output_subdir="$OUTPUT_DIR/$context/$annotation"

        # Process mC
        mc_file=$(find_dmr_file "$context" "$annotation" "mc")
        if [[ -n "$mc_file" ]] && [[ -f "$mc_file" ]]; then
            mc_stats=$(process_dmr_file "$mc_file" "$output_subdir" "mc" "$context" "$annotation")
            echo -e "${context}\t${annotation}\tmc\t${mc_stats}" >> "$SUMMARY_FILE"
        else
            log_msg "  WARNING: No mC DMR file found for $context/$annotation"
            echo -e "${context}\t${annotation}\tmc\t0\t0\t0\t0" >> "$SUMMARY_FILE"
        fi

        # Process hmC
        hmc_file=$(find_dmr_file "$context" "$annotation" "hmc")
        if [[ -n "$hmc_file" ]] && [[ -f "$hmc_file" ]]; then
            hmc_stats=$(process_dmr_file "$hmc_file" "$output_subdir" "hmc" "$context" "$annotation")
            echo -e "${context}\t${annotation}\thmc\t${hmc_stats}" >> "$SUMMARY_FILE"
        else
            log_msg "  WARNING: No hmC DMR file found for $context/$annotation"
            echo -e "${context}\t${annotation}\thmc\t0\t0\t0\t0" >> "$SUMMARY_FILE"
        fi

        # Generate comparison stats
        generate_comparison_stats "$output_subdir" "$context" "$annotation"
    done
done

log_msg ""
log_msg "=============================================="
log_msg "Processing Complete!"
log_msg "=============================================="
log_msg ""
log_msg "Summary file: $SUMMARY_FILE"
log_msg ""

# Display summary
log_msg "=== SUMMARY ==="
column -t -s$'\t' "$SUMMARY_FILE"

log_msg ""
log_msg "Next step: Run R visualization script"
log_msg "  Rscript $SCRIPT_DIR/visualize_dmr_results.R"
