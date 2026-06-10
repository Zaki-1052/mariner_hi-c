#!/bin/bash
# poster/pub_browser/run_mecp2_noncg_browser.sh
#
# Batch genome browser views of non-CG methylation at MeCP2 binding sites.
# Reads the peaks TSV, generates a per-peak YAML config, and runs pub_browser.R.
#
# Prerequisites:
#   1. CHG/CHH bigwigs rsynced from HPC to exports/bigwig/{CHG,CHH}_mc/
#   2. MeCP2 bigwigs at /Users/zakiralibhai/sdsc/bigwigs/MeCP2{Ctrl,Mut}.bw
#
# Usage:
#   bash poster/pub_browser/run_mecp2_noncg_browser.sh          # all 191 peaks
#   bash poster/pub_browser/run_mecp2_noncg_browser.sh --top 20  # top 20 by signal
#   bash poster/pub_browser/run_mecp2_noncg_browser.sh --top-yamls  # just the 2 curated YAMLs

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
PUB_BROWSER="${SCRIPT_DIR}/pub_browser.R"
PEAKS_TSV="${REPO_ROOT}/biomodal/downstream/plots/visualizations/tables/mecp2_noncg_peaks_for_browser.tsv"
OUTPUT_BASE="${SCRIPT_DIR}/output/mecp2_noncg"

# BigWig paths
BW_DIR="/Users/zakiralibhai/sdsc/bigwigs"
MECP2_CTRL="${BW_DIR}/MeCP2Ctrl.bw"
MECP2_MUT="${BW_DIR}/MeCP2Mut.bw"

CHG_BW_DIR="${REPO_ROOT}/biomodal/downstream/modality/exports/bigwig/CHG_mc"
CHH_BW_DIR="${REPO_ROOT}/biomodal/downstream/modality/exports/bigwig/CHH_mc"

# Per-sample bigwig names (same naming as CG mc bigwigs)
CTRL_SAMPLES=("evoC-Bap1-ctrl-F" "evoC-Bap1-ctrl-M" "evoC-Bap1-ctrl-F-B2" "evoC-Bap1-ctrl-M-B2")
MUT_SAMPLES=("evoC-Bap1-mut-F" "evoC-Bap1-mut-M" "evoC-Bap1-mut-F-B2" "evoC-Bap1-mut-M-B2")

EXTEND=10000       # 10kb flank around each 400bp peak
MAX_PEAKS=""        # empty = all peaks

# ---- Parse arguments --------------------------------------------------------
if [[ "$1" == "--top-yamls" ]]; then
    echo "Running curated top-peak YAMLs only..."
    Rscript "${PUB_BROWSER}" --config "${SCRIPT_DIR}/configs/mecp2_noncg_top_chg.yaml"
    Rscript "${PUB_BROWSER}" --config "${SCRIPT_DIR}/configs/mecp2_noncg_top_chh.yaml"
    echo ""
    echo "Done. Outputs in: ${OUTPUT_BASE}/"
    exit 0
fi

if [[ "$1" == "--top" && -n "$2" ]]; then
    MAX_PEAKS="$2"
    echo "Running top ${MAX_PEAKS} peaks by mean signal..."
fi

# ---- Validate bigwig availability -------------------------------------------
missing=0
if [[ ! -f "$MECP2_CTRL" ]]; then echo "ERROR: Missing $MECP2_CTRL"; missing=1; fi
if [[ ! -f "$MECP2_MUT" ]]; then echo "ERROR: Missing $MECP2_MUT"; missing=1; fi

chg_count=$(ls "${CHG_BW_DIR}"/*.bw 2>/dev/null | wc -l)
chh_count=$(ls "${CHH_BW_DIR}"/*.bw 2>/dev/null | wc -l)
if [[ "$chg_count" -eq 0 ]]; then
    echo "WARNING: No CHG bigwigs found in ${CHG_BW_DIR}/"
    echo "  Run the HPC export + rsync first (see mecp2-noncg.md)"
    missing=1
fi
if [[ "$chh_count" -eq 0 ]]; then
    echo "WARNING: No CHH bigwigs found in ${CHH_BW_DIR}/"
    echo "  Run the HPC export + rsync first (see mecp2-noncg.md)"
    missing=1
fi
if [[ "$missing" -eq 1 ]]; then
    echo ""
    echo "Aborting. Fix missing bigwigs and re-run."
    exit 1
fi

echo "BigWig check passed: MeCP2 (2), CHG (${chg_count}), CHH (${chh_count})"

# ---- Build bigwig path lists ------------------------------------------------
build_bw_list() {
    local bw_dir="$1"
    shift
    local samples=("$@")
    local paths=""
    for s in "${samples[@]}"; do
        local bw="${bw_dir}/${s}.genome.mc_bedgraph.bw"
        if [[ -f "$bw" ]]; then
            if [[ -n "$paths" ]]; then paths="${paths}, "; fi
            paths="${paths}${bw}"
        fi
    done
    echo "$paths"
}

CHG_CTRL_BWS=$(build_bw_list "$CHG_BW_DIR" "${CTRL_SAMPLES[@]}")
CHG_MUT_BWS=$(build_bw_list "$CHG_BW_DIR" "${MUT_SAMPLES[@]}")
CHH_CTRL_BWS=$(build_bw_list "$CHH_BW_DIR" "${CTRL_SAMPLES[@]}")
CHH_MUT_BWS=$(build_bw_list "$CHH_BW_DIR" "${MUT_SAMPLES[@]}")

# ---- Process peaks ----------------------------------------------------------
mkdir -p "${OUTPUT_BASE}"
TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

total=0
succeeded=0
failed=0

# TSV columns: chr start end context mecp2_direction ctrl_mean mut_mean all_mean max_value peak_id
# Header is at the LAST line (line 192); data lines are 1-191, sorted by all_mean descending
while IFS=$'\t' read -r chr start end context direction ctrl_mean mut_mean all_mean max_value peak_id; do
    # Skip header
    [[ "$chr" == "chr" ]] && continue
    # Skip empty lines
    [[ -z "$chr" ]] && continue

    total=$((total + 1))

    if [[ -n "$MAX_PEAKS" && "$total" -gt "$MAX_PEAKS" ]]; then
        break
    fi

    # Region with flanking
    region_start=$((start - EXTEND))
    if [[ "$region_start" -lt 1 ]]; then region_start=1; fi
    region_end=$((end + EXTEND))
    region="${chr}:${region_start}-${region_end}"

    # Output name
    safe_id=$(echo "${peak_id}" | tr ':' '_' | tr '-' '_')
    out_name="${context}_${safe_id}"
    out_path="${OUTPUT_BASE}/${out_name}"

    # Select bigwigs based on context
    if [[ "$context" == "CHG" ]]; then
        meth_name="mCHG"
        meth_color="#7570B3"
        meth_ctrl_bws="$CHG_CTRL_BWS"
        meth_mut_bws="$CHG_MUT_BWS"
    else
        meth_name="mCHH"
        meth_color="#1B9E77"
        meth_ctrl_bws="$CHH_CTRL_BWS"
        meth_mut_bws="$CHH_MUT_BWS"
    fi

    # Format all_mean as percentage for highlight label
    pct=$(awk "BEGIN {printf \"%.2f\", ${all_mean} * 100}")

    # Generate per-peak YAML
    yaml_file="${TMPDIR}/${out_name}.yaml"
    cat > "$yaml_file" <<YAML
genome: mm10
region: '${region}'
extend_bp: 0
labels: ['control', '*Bap1*^f/f^,Math1-cre']
output: ${out_path}

marks:
  - name: MeCP2
    color: '#D95F02'
    ctrl: ${MECP2_CTRL}
    mut:  ${MECP2_MUT}

  - name: '${meth_name}'
    color: '${meth_color}'
    sparse: true
    diff: true
    ctrl: [${meth_ctrl_bws}]
    mut:  [${meth_mut_bws}]
    average: true

highlights:
  - region: '${chr}:${start}-${end}'
    label: '${direction} (${meth_name} ${pct}%)'
YAML

    echo "[${total}] ${context} ${peak_id} (${pct}% mean) -> ${out_name}"
    if Rscript "${PUB_BROWSER}" --config "$yaml_file" 2>&1; then
        succeeded=$((succeeded + 1))
    else
        echo "  FAILED: ${peak_id}"
        failed=$((failed + 1))
    fi

done < "$PEAKS_TSV"

echo ""
echo "========================================"
echo "Batch complete"
echo "Total: ${total}, Succeeded: ${succeeded}, Failed: ${failed}"
echo "Output: ${OUTPUT_BASE}/"
echo "========================================"
