












#!/usr/bin/env bash

###############################################################################
# modality XPLR Core Workflow Pipeline
#
# This script runs a comprehensive analysis pipeline including:
# - Biological QC
# - Feature extraction
# - DMR calling and visualisation
#
# Input requirements:
# - config.txt file in current directory with required parameters
# - Zarr store(s) containing methylation data
# - Sample metadata CSV or TSV file
# - Directory with region BED/TSV files
#
# This script is intended for running on modality XPLR version 1.0.0 or later.
#
###############################################################################

#------------------------------------------------------------------------------
# CONFIGURATION
#------------------------------------------------------------------------------
# Set script behaviour
set -euo pipefail

# Record start time for performance tracking
start_time=$(date +%s)

# Log file setup - uncomment to enable logging to file
# LOG_FILE="./modality_pipeline_$(date +%Y%m%d_%H%M%S).log"
# exec > >(tee -a "$LOG_FILE") 2>&1


#------------------------------------------------------------------------------
# HELPER FUNCTIONS
#------------------------------------------------------------------------------
# Logging functions
log_info() {
  echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_warning() {
  echo "[WARNING] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

log_error() {
  echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2
}

log_debug() {
  if [[ "${DEBUG:-false}" == "true" ]]; then
    echo "[DEBUG] $(date '+%Y-%m-%d %H:%M:%S') - $1"
  fi
}


 # Helper to trim whitespace
trim() {
  awk '{$1=$1};1'
}

# Function to check if a header exists in metadata
header_exists() {
  local col="$1"
  local headers=("${@:2}")

  for h in "${headers[@]}"; do
    [[ "$h" == "$col" ]] && return 0
  done
  return 1
}

# Function to create directories with validation
create_directory() {
  local dir="$1"

  if [[ ! -d "$dir" ]]; then
    mkdir -p "$dir" 2>/dev/null || { log_error "Cannot create directory $dir"; exit 1; }
  fi
  if [[ ! -w "$dir" ]]; then
    log_error "No write access to directory $dir"
    exit 1
  fi
  log_info "Directory ready: $dir"
}

# Function to copy reports to destination
copy_reports() {
  local source_dir="$1"
  local pattern="$2"
  local dest_dir="$3"
  local description="$4"

  log_info "Copying $description HTML reports to $dest_dir..."
  local files=$(find "$source_dir" -type f -name "$pattern")

  if [[ -z "$files" ]]; then
    log_error "No $description HTML reports found in $source_dir"
    return 1
  fi

  local count=0
  while IFS= read -r file; do
    cp "$file" "$dest_dir/"
    ((count++))
  done <<< "$files"

  log_info "Copied $count $description HTML report(s)"
  return 0
}


#------------------------------------------------------------------------------
# INITIALISATION
#------------------------------------------------------------------------------
printf "\n=======================================\n"
printf "STEP 0 - CONFIGURATION & INITIALISATION\n"
printf "=======================================\n"
log_info "Starting modality XPLR Core Workflow Pipeline"

# Default values for variables that could be customised
CONFIG_FILE="./config.txt"
DEBUG=false  # Set to true to enable debug output


#------------------------------------------------------------------------------
# LOAD CONFIGURATION
#------------------------------------------------------------------------------
log_info "Reading configuration from $CONFIG_FILE"

# Validate config file exists
if [[ ! -f "$CONFIG_FILE" ]]; then
  log_error "config.txt not found in the current working directory: $(pwd)"
  log_error "Please ensure config.txt is present before running this script."
  exit 1
fi

# Read config.txt variables
ZARR=""                 # Path(s) to zarr store(s), comma-separated
METADATA=""             # Path to sample metadata CSV or TSV file
GROUP_COLUMN=""         # Column name in metadata for grouping samples
COVARIATES=""           # Space-separated list of covariate column names (optional)
CONDITION_ORDER=""      # Comma-separated order of conditions for DMR analysis (optional)
QC_REGION=""            # Region for Biological QC analysis (optional)
REGIONS_DIR=""          # Directory containing region BED/TSV files
MAIN_OUTPUT=""          # Base output directory
DEPTH_FILTER=""         # Minimum depth filter for analyses, e.g. 10 (optional)
OVERDISPERSION=""       # Apply overdispersion correction: "True" or "False" (optional)

# Parse config file
while IFS='=' read -r key value; do
  key=$(echo "$key" | trim)
  value=$(echo "$value" | trim)
  case "$key" in
    Zarr) ZARR="$value" ;;
    Metadata) METADATA="$value" ;;
    Group_Column) GROUP_COLUMN="$value" ;;
    Covariates) COVARIATES="$value" ;;
    Condition_Order) CONDITION_ORDER="$value" ;;
    QC_Region) QC_REGION="$value" ;;
    Regions_Directory) REGIONS_DIR="$value" ;;
    Main_Output) MAIN_OUTPUT="$value" ;;
    Depth_Filter) DEPTH_FILTER="$value" ;;
    Overdispersion) OVERDISPERSION="$value" ;;
    # You could add additional parameters here as needed
  esac
done < "$CONFIG_FILE"

# Log parsed configuration
printf "\n\n>>Provided inputs\n"
log_info "ZARR path(s): $ZARR"
log_info "Metadata file: $METADATA"
log_info "Group column: $GROUP_COLUMN"
log_info "Covariates: $COVARIATES"
log_info "Condition order: $CONDITION_ORDER"
log_info "QC region: $QC_REGION"
log_info "Regions directory: $REGIONS_DIR"
log_info "Main output directory: $MAIN_OUTPUT"
log_info "Depth filter: $DEPTH_FILTER"
log_info "Overdispersion correction: $OVERDISPERSION"


#------------------------------------------------------------------------------
# VALIDATE INPUTS
#------------------------------------------------------------------------------
printf "\n\n>>Validation & directory creation\n"
log_info "Validating inputs..."

# Parse and validate zarr paths
if [[ -n "$ZARR" ]]; then
  # Split comma-separated zarr paths into array
  IFS=',' read -r -a ZARR_RAW_ARRAY <<< "$ZARR"

  # Process each zarr path to remove quotes and trim whitespace
  ZARR_ARRAY=()
  for zarr_raw in "${ZARR_RAW_ARRAY[@]}"; do
    # Trim leading/trailing whitespace
    zarr_clean=$(echo "$zarr_raw" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

    # Remove surrounding single or double quotes if present
    zarr_clean=$(echo "$zarr_clean" | sed 's/^["'"'"']\(.*\)["'"'"']$/\1/')

    # Add to processed array if not empty
    if [[ -n "$zarr_clean" ]]; then
      ZARR_ARRAY+=("$zarr_clean")
    fi
  done

  # Log the number of zarr paths identified and their strings
  log_info "Identified ${#ZARR_ARRAY[@]} zarr path(s):"
  for i in "${!ZARR_ARRAY[@]}"; do
    log_info "  [$((i+1))] ${ZARR_ARRAY[i]}"
  done

  # Validate each zarr path
  for zarr in "${ZARR_ARRAY[@]}"; do
    if [[ ! -s "$zarr" ]]; then
      log_error "Zarr file '$zarr' doesn't contain any data or is not readable"
      exit 1
    fi
  done
else
  log_error "No zarr paths specified in config.txt"
  exit 1
fi

# Process and validate metadata file path
if [[ -n "$METADATA" ]]; then
  # Trim leading/trailing whitespace
  METADATA=$(echo "$METADATA" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

  # Remove surrounding single or double quotes if present
  METADATA=$(echo "$METADATA" | sed 's/^["'"'"']\(.*\)["'"'"']$/\1/')

  log_info "Processed metadata file path: $METADATA"
else
  log_error "No metadata file specified in config.txt"
  exit 1
fi

# Process and validate main output directory path
if [[ -n "$MAIN_OUTPUT" ]]; then
  # Trim leading/trailing whitespace
  MAIN_OUTPUT=$(echo "$MAIN_OUTPUT" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

  # Remove surrounding single or double quotes if present
  MAIN_OUTPUT=$(echo "$MAIN_OUTPUT" | sed 's/^["'"'"']\(.*\)["'"'"']$/\1/')

  log_info "Processed main output directory: $MAIN_OUTPUT"
else
  log_error "No main output directory specified in config.txt"
  exit 1
fi

# Process and validate regions directory path
if [[ -n "$REGIONS_DIR" ]]; then
  # Trim leading/trailing whitespace
  REGIONS_DIR=$(echo "$REGIONS_DIR" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

  # Remove surrounding single or double quotes if present
  REGIONS_DIR=$(echo "$REGIONS_DIR" | sed 's/^["'"'"']\(.*\)["'"'"']$/\1/')

  log_info "Processed regions directory: $REGIONS_DIR"
else
  log_error "No regions directory specified in config.txt"
  exit 1
fi

# Check metadata file
log_info "Validating metadata file headers..."
if [[ ! -r "$METADATA" ]]; then
  log_error "Metadata file not found or not readable: $METADATA"
  exit 1
fi

# Check group column and covariates in metadata
HEADER_LINE=$(head -n 1 "$METADATA")

# Remove any line ending characters (Windows \r\n, Mac \r, Unix \n)
HEADER_LINE=$(echo "$HEADER_LINE" | tr -d '\r\n')

# Detect delimiter: use tab if present, otherwise comma
if [[ "$HEADER_LINE" == *$'\t'* ]]; then
  DELIM=$'\t'
else
  DELIM=','
fi

# Split header into array using the detected delimiter
IFS="$DELIM" read -r -a HEADERS_RAW <<< "$HEADER_LINE"

# Trim whitespace from each header individually
HEADERS=()
for header in "${HEADERS_RAW[@]}"; do
  # Trim leading and trailing whitespace from each header
  header_clean=$(echo "$header" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
  if [[ -n "$header_clean" ]]; then
    HEADERS+=("$header_clean")
  fi
done

if [[ -n "$GROUP_COLUMN" ]]; then
  if ! header_exists "$GROUP_COLUMN" "${HEADERS[@]}"; then
    log_error "Group column '$GROUP_COLUMN' not found in metadata file header. Please review the metadata file and provide a valid group column."
    exit 1
  fi
fi

if [[ -n "$COVARIATES" ]]; then
  # Convert space-separated list to array
  read -r -a COVS <<< "$COVARIATES"
  for cov in "${COVS[@]}"; do
    cov=$(echo "$cov" | trim)
    if ! header_exists "$cov" "${HEADERS[@]}"; then
      log_error "Covariate column '$cov' not found in metadata file header"
      exit 1
    fi
  done
fi

log_info "Metadata validation passed"

# Check regions directory
if [[ ! -d "$REGIONS_DIR" ]]; then
  log_error "Regions directory not found: $REGIONS_DIR"
  exit 1
fi


#------------------------------------------------------------------------------
# SETUP OUTPUT DIRECTORIES
#------------------------------------------------------------------------------
log_info "Creating output directories under $MAIN_OUTPUT..."

# Create output directories
for sub in "" "Results" "Reports"; do
  dir="$MAIN_OUTPUT"
  [[ -n "$sub" ]] && dir="$MAIN_OUTPUT/$sub"
  create_directory "$dir"
done

# Define output subdirectories for easier reference
RESULTS_DIR="$MAIN_OUTPUT/Results"
REPORTS_DIR="$MAIN_OUTPUT/Reports"

# Create subdirectories for each region file found in the regions directory, to organise outputs by region file input.
log_info "Creating region-specific subdirectories..."
for bedfile in "$REGIONS_DIR"/*.bed "$REGIONS_DIR"/*.bed.gz "$REGIONS_DIR"/*.tsv "$REGIONS_DIR"/*.tsv.gz; do
  # Skip if no files match pattern
  [[ -e "$bedfile" ]] || continue
  [[ "${METADATA##*.}" == "tsv" && $(realpath "$bedfile") == $(realpath "$METADATA") ]] && continue

  # Extract region file name without path and extension
  region_name=$(basename "$bedfile")
  region_name="${region_name%%.bed.gz}"
  region_name="${region_name%%.bed}"
  region_name="${region_name%%.tsv.gz}"
  region_name="${region_name%%.tsv}"

  # Create region-specific subdirectories in both Results and Reports
  region_results_dir="$RESULTS_DIR/$region_name"
  region_reports_dir="$REPORTS_DIR/$region_name"
  create_directory "$region_results_dir"
  create_directory "$region_reports_dir"
  log_info "Created region directories: $region_results_dir and $region_reports_dir"
done


#------------------------------------------------------------------------------
# STEP 1: BIOLOGICAL QC
#------------------------------------------------------------------------------
printf "\n\n======================================\n"
printf "STEP 1: RUNNING BIOLOGICAL QC ANALYSES\n"
printf "======================================\n"

# Build region argument if QC_REGION is set
REGION_ARG=""
[[ -n "$QC_REGION" ]] && REGION_ARG="--region $QC_REGION"

# Run QC without aggregation
log_info "Running biological QC without aggregation..."
modality -v biological-qc \
  --zarr-path "${ZARR_ARRAY[@]}" \
  --sample-sheet "$METADATA" \
  $REGION_ARG \
  --order-by-group "$GROUP_COLUMN" \
  --output-dir "$RESULTS_DIR"

# Copy QC reports to Reports directory
copy_reports "$RESULTS_DIR" "biological_qc_report*.html" "$REPORTS_DIR" "Biological QC" || exit 1


#------------------------------------------------------------------------------
# STEP 2: FEATURE EXTRACTION
#------------------------------------------------------------------------------
printf "\n\n==================================================\n"
printf "STEP 2: RUNNING FEATURE EXTRACTION ON REGION FILES\n"
printf "==================================================\n"

# Process each region file for feature extraction
for bedfile in "$REGIONS_DIR"/*.bed "$REGIONS_DIR"/*.bed.gz "$REGIONS_DIR"/*.tsv "$REGIONS_DIR"/*.tsv.gz; do
  # Skip if no files match pattern
  [[ -e "$bedfile" ]] || continue
  [[ "${METADATA##*.}" == "tsv" && $(realpath "$bedfile") == $(realpath "$METADATA") ]] && continue

  log_info "Processing region file: $bedfile"

  # Extract region file name without path and extension
  region_name=$(basename "$bedfile")
  region_name="${region_name%%.bed.gz}"
  region_name="${region_name%%.bed}"
  region_name="${region_name%%.tsv.gz}"
  region_name="${region_name%%.tsv}"

  # Set region-specific output directories
  region_output_dir="$RESULTS_DIR/$region_name"
  region_reports_dir="$REPORTS_DIR/$region_name"

  # 1. Extract count for total_c
  log_info "Running feature extraction: count total_c..."
  modality -v get count \
    --zarr-path "${ZARR_ARRAY[@]}" \
    --sample-sheet "$METADATA" \
    --bedfile "$bedfile" \
    --fields total_c \
    --output-dir "$region_output_dir"

  # 2. Extract mean for total_c with plots generated automatically
  log_info "Running feature extraction: mean total_c (plots generated automatically)..."
  modality -v get mean \
    --zarr-path "${ZARR_ARRAY[@]}" \
    --sample-sheet "$METADATA" \
    --bedfile "$bedfile" \
    --fields total_c \
    --order-by-group "$GROUP_COLUMN" \
    --output-dir "$region_output_dir"

  # Copy mean HTML reports immediately
  for extract_dir in "$region_output_dir"/Extract_*; do
    if [[ -d "$extract_dir" ]]; then
      find "$extract_dir" -type f -name "*mean_extract_report*.html" -exec cp {} "$region_reports_dir/" \; 2>/dev/null || true
    fi
  done

  # 3. Extract regional-frac for mc and hmc with depth filter (plots generated automatically)
  log_info "Running feature extraction: regional-frac for mc and hmc with minimum depth $DEPTH_FILTER (plots generated automatically)..."
  # Build the optional argument if DEPTH_FILTER is set
  MCV_ARG=""
  [[ -n "$DEPTH_FILTER" ]] && MCV_ARG="--min-coverage $DEPTH_FILTER"
  # --min-coverage 30, if DEPTH_FILTER=30 in config.txt

  modality -v get regional-frac \
    --zarr-path "${ZARR_ARRAY[@]}" \
    --sample-sheet "$METADATA" \
    --bedfile "$bedfile" \
    --fields mc hmc \
    $MCV_ARG \
    --order-by-group "$GROUP_COLUMN" \
    --output-dir "$region_output_dir"

  # Copy regional-frac HTML reports immediately
  for extract_dir in "$region_output_dir"/Extract_*; do
    if [[ -d "$extract_dir" ]]; then
      find "$extract_dir" -type f -name "*regional-frac_extract_report*.html" -exec cp {} "$region_reports_dir/" \; 2>/dev/null || true
    fi
  done

  log_info "Feature extraction completed for $region_name"
done

log_info "Feature extraction HTML reports copied to region-specific directories"


#------------------------------------------------------------------------------
# STEP 3: DMR CALLING
#------------------------------------------------------------------------------
printf "\n\n===========================\n"
printf "STEP 3: RUNNING DMR CALLING\n"
printf "===========================\n"

# Build covariate / depth_filter / condition_order / overdispersion flags if present
COV_FLAGS=""
DEPTH_ARG=""
CONDITION_ARG=""
OVERDISPERSION_ARG=""

[[ -n "$DEPTH_FILTER" ]] && DEPTH_ARG="--filter-context-depth $DEPTH_FILTER"
[[ -n "$CONDITION_ORDER" ]] && CONDITION_ARG="--condition-order $CONDITION_ORDER"
[[ "$OVERDISPERSION" == "True" ]] && OVERDISPERSION_ARG="--overdispersion"

if [[ -n "$COVARIATES" ]]; then
  # Use space-separated covariates directly for --covariates flag
  COV_FLAGS="--covariates $COVARIATES"
fi

# Run DMR calling on each region file individually
for bedfile in "$REGIONS_DIR"/*.bed "$REGIONS_DIR"/*.bed.gz "$REGIONS_DIR"/*.tsv "$REGIONS_DIR"/*.tsv.gz; do
  # Skip if no files match pattern
  [[ -e "$bedfile" ]] || continue
  [[ "${METADATA##*.}" == "tsv" && $(realpath "$bedfile") == $(realpath "$METADATA") ]] && continue

  log_info "Running DMR calling on region file: $bedfile"

  # Extract region file name without path and extension
  region_name=$(basename "$bedfile")
  region_name="${region_name%%.bed.gz}"
  region_name="${region_name%%.bed}"
  region_name="${region_name%%.tsv.gz}"
  region_name="${region_name%%.tsv}"

  # Set region-specific output directory
  region_output_dir="$RESULTS_DIR/$region_name"

  # Run DMR calling
  # --covariates cov1 cov2, if COVARIATES="cov1 cov2" in config.txt
  # --filter-context-depth 15, if DEPTH_FILTER=15 in config.txt
  # --condition-order condition1,condition2, if CONDITION_ORDER="condition1,condition2" in config.txt
  # --overdispersion is applied if OVERDISPERSION="True" in config.txt

  modality -v dmr call \
    --zarr-path "${ZARR_ARRAY[@]}" \
    --sample-sheet "$METADATA" \
    --condition-array-name "$GROUP_COLUMN" \
    --methylation-contexts num_mc num_hmc \
    --bedfile "$bedfile" \
    $COV_FLAGS \
    $DEPTH_ARG \
    $CONDITION_ARG \
    $OVERDISPERSION_ARG \
    --output-dir "$region_output_dir"

  log_info "DMR calling completed for $bedfile"
done

log_info "All DMR calling completed"


#------------------------------------------------------------------------------
# STEP 4: DMR VISUALISATION
#------------------------------------------------------------------------------
printf "\n\n=====================================\n"
printf "STEP 4: GENERATING DMR VISUALISATIONS\n"
printf "=====================================\n"

# Loop through each region directory
for region_dir in "$RESULTS_DIR"/*; do
  [[ -d "$region_dir" ]] || continue

  region_name=$(basename "$region_dir")
  # Skip if this is not a region directory (e.g., if it's a BioQC directory)
  [[ "$region_name" == "BioQC_"* ]] && continue

  log_info "Processing DMR visualisation for region: $region_name"

  # Set region-specific reports directory
  region_reports_dir="$REPORTS_DIR/$region_name"

  # Loop through each DMR results directory within this region
  for dmr_dir in "$region_dir"/DMR_*; do
    [[ -d "$dmr_dir" ]] || continue

    # For each .bed file in the DMR directory
    for dmrresults in "$dmr_dir"/*.bed; do
      [[ -e "$dmrresults" ]] || continue

      log_info "Generating DMR plots for $dmrresults..."

      # Generate default plots - output to the same region directory
      modality -v dmr plot \
        --dmr-results "$dmrresults" \
        --output-dir "$region_dir"

      # Copy unfiltered DMR reports immediately
      for dmr_report_dir in "$region_dir"/DMR_Report_*; do
        if [[ -d "$dmr_report_dir" ]]; then
          find "$dmr_report_dir" -type f -name "*.html" -not -name "*max_q_*" -exec cp {} "$region_reports_dir/" \; 2>/dev/null || true
        fi
      done

      # Generate filtered plots (q-value < 0.05) - output to the same region directory
      modality dmr plot \
        --dmr-results "$dmrresults" \
        --filter-qvalue 0.05 \
        --output-dir "$region_dir"

      # Copy filtered DMR reports immediately
      for dmr_report_dir in "$region_dir"/DMR_Report_*; do
        if [[ -d "$dmr_report_dir" ]]; then
          find "$dmr_report_dir" -type f -name "*max_q_*.html" -exec cp {} "$region_reports_dir/" \; 2>/dev/null || true
        fi
      done

      log_info "DMR plots generated for $dmrresults"
    done
  done

  log_info "DMR visualisation completed for region: $region_name"
done

log_info "DMR HTML reports copied to region-specific directories"


#------------------------------------------------------------------------------
# COMPLETION
#------------------------------------------------------------------------------
log_info "Pipeline completed successfully. Reports are available in $REPORTS_DIR and region-specific subdirectories"

# Calculate and report duration
end_time=$(date +%s)
duration=$((end_time - start_time))
log_info "Script completed in $duration seconds"

# Exit cleanly
exit 0


