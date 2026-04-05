# Updating the Downstream Pipeline for a New Modality XPLR Run

This document describes how to update the biomodal downstream analysis pipeline when
a new modality XPLR run is added. It covers the four files that contain hardcoded
run paths, the order in which to update them, and how to verify that no stale
references remain.

## Prerequisites

- A completed modality XPLR run on Expanse, with results in `outputs/run-N/`
- Access to rsync the results locally
- Familiarity with the timestamped directory structure inside each run

## Background: Directory Structure

Each run lives in `downstream/modality/outputs/run-N/` and contains context
subdirectories (`outputs_CG`, `outputs_CHG`, `outputs_CHH`). Inside each context's
`Results/` directory, there are region folders (one per genomic annotation), and
within those, timestamped subdirectories:

```
outputs/run-N/outputs_CG/Results/
  gencode.vM25.mouse.genes.annotation/
    BioQC_20260402_185215/
      biological_qc_report_8_samples_20260402_185215.json
    Extract_20260402_190326/
      Extract_mc_regional-frac_20260402_190326.tsv.gz
      Extract_hmc_regional-frac_20260402_190326.tsv.gz
    DMR_20260402_191818/
      DMR_mc_control__mutant_20260402_191818.bed
      DMR_hmc_control__mutant_20260402_191818.bed
    DMR_Report_20260402_191818/
      (volcano plots, filtered BED files)
  gencode.vM25.mouse.promoters.annotation/
    DMR_20260402_192045/
      ...
  (... 4 more regions ...)
```

Key observations:

- Timestamps are in `YYYYMMDD_HHMMSS` format.
- BED filenames embed the same timestamp as their parent directory.
- The BioQC JSON filename includes the sample count (e.g., `8_samples`).
- Some regions may have multiple `DMR_*` directories from re-runs. Always use the
  **latest** timestamp.
- `DMR_Report_*` directories are separate from `DMR_*` directories. The `find`
  commands below exclude report directories.

## Run History

| Run | Samples | Sex Covariate | Notes |
|-----|---------|---------------|-------|
| run-1 | 4 (shallow-seq) | Yes | No significant DMRs |
| run-2 | 4 (shallow-seq) | No | First significant results |
| run-3 | 4 (deep-seq) | No | Improved power from deeper sequencing |
| run-4 | 8 (deep-seq) | No | Added batch 2 replicates |
| run-5 | 8 (deep-seq) | Yes | Current primary: sex covariate included |

## Step-by-Step Update Procedure

### Step 1: Transfer Results from Expanse

```bash
rsync -avz --progress \
  expanse:/expanse/lustre/projects/csd940/zalibhai/biomodal/modality/outputs/run-N/ \
  downstream/modality/outputs/run-N/
```

Replace `N` with the new run number throughout this document.

### Step 2: Discover Timestamps

Run these commands from the repository root to identify the timestamps you will need.

**BioQC timestamp and filename:**

```bash
find downstream/modality/outputs/run-N/outputs_CG/Results \
  -type d -name "BioQC_*" | sort
```

Note the full JSON filename inside (includes sample count).

**Extract timestamp (gene body only):**

```bash
find downstream/modality/outputs/run-N/outputs_CG/Results/gencode.vM25.mouse.genes.annotation \
  -type d -name "Extract_*" | sort
```

**DMR timestamps per region (excluding report directories):**

```bash
find downstream/modality/outputs/run-N/outputs_CG/Results \
  -type d -name "DMR_*" ! -name "DMR_Report_*" | sort
```

For regions with multiple `DMR_*` directories, record only the latest timestamp.

**CHG context DMR timestamp (gene body only):**

```bash
find downstream/modality/outputs/run-N/outputs_CHG/Results/gencode.vM25.mouse.genes.annotation \
  -type d -name "DMR_*" ! -name "DMR_Report_*" | sort
```

Collect the timestamps into a reference like this before editing any files:

```
BioQC:       20260402_185215  (filename: biological_qc_report_8_samples_20260402_185215.json)
Extract:     20260402_190326
genes:       20260402_191818
cpg_islands: 20260402_191006
cpg_shores:  20260402_191531
cpg_shelves: 20260402_191227
promoters:   20260402_192045
tss_region:  20260402_192302
CHG genes:   20260402_222845
```

### Step 3: Update `_shared_config.R` (Central Config)

**File:** `downstream/scripts/viz_sections/_shared_config.R`

This is the most important file. Every visualization section script sources it.
Three path blocks must be updated:

#### 3a. `DATA_PATHS`

Update the run number and all timestamps in the `DATA_PATHS` list (lines 14-32).
There are 13 entries total:

| Key | Region | Context |
|-----|--------|---------|
| `mc_dmr` | Gene bodies | CG mC |
| `hmc_dmr` | Gene bodies | CG hmC |
| `bioqc_json` | BioQC | CG |
| `cpg_islands_mc` | CpG islands | CG mC |
| `cpg_shores_mc` | CpG shores | CG mC |
| `cpg_shelves_mc` | CpG shelves | CG mC |
| `promoters_mc` | Promoters | CG mC |
| `tss_mc` | TSS regions | CG mC |
| `cpg_islands_hmc` | CpG islands | CG hmC |
| `cpg_shores_hmc` | CpG shores | CG hmC |
| `cpg_shelves_hmc` | CpG shelves | CG hmC |
| `promoters_hmc` | Promoters | CG hmC |
| `tss_hmc` | TSS regions | CG hmC |

For each entry, the run number and timestamp appear in both the directory path and
the filename. For example:

```r
# Before (run-5 example)
mc_dmr = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260402_191818/DMR_mc_control__mutant_20260402_191818.bed"),

# After (run-6 example)
mc_dmr = file.path(BASE_DIR, "modality/outputs/run-6/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_YYYYMMDD_HHMMSS/DMR_mc_control__mutant_YYYYMMDD_HHMMSS.bed"),
```

The `upstream_csv` path does NOT change -- it points to the upstream DUET pipeline
output, which is independent of downstream runs.

The `bioqc_json` entry has a unique filename pattern that includes the sample count:

```r
bioqc_json = file.path(BASE_DIR, "modality/outputs/run-N/outputs_CG/Results/BioQC_TIMESTAMP/biological_qc_report_X_samples_TIMESTAMP.json"),
```

Verify the exact filename after rsync, as the sample count may change between runs.

#### 3b. `EXTRACT_PATHS`

Update the run number and Extract timestamp in the `EXTRACT_PATHS` list (lines 35-39).
Three entries share the same timestamp:

```r
EXTRACT_PATHS <- list(
  extract_dir       = file.path(BASE_DIR, "modality/outputs/run-N/.../Extract_TIMESTAMP"),
  mc_regional_frac  = file.path(BASE_DIR, "modality/outputs/run-N/.../Extract_TIMESTAMP/Extract_mc_regional-frac_TIMESTAMP.tsv.gz"),
  hmc_regional_frac = file.path(BASE_DIR, "modality/outputs/run-N/.../Extract_TIMESTAMP/Extract_hmc_regional-frac_TIMESTAMP.tsv.gz")
)
```

#### 3c. `CHG_DATA_PATHS`

Update the run number and CHG gene body DMR timestamp in the `CHG_DATA_PATHS` list
(lines 42-45). Two entries share the same timestamp:

```r
CHG_DATA_PATHS <- list(
  mc_dmr  = file.path(BASE_DIR, "modality/outputs/run-N/.../DMR_TIMESTAMP/DMR_mc_control__mutant_TIMESTAMP.bed"),
  hmc_dmr = file.path(BASE_DIR, "modality/outputs/run-N/.../DMR_TIMESTAMP/DMR_hmc_control__mutant_TIMESTAMP.bed")
)
```

#### What NOT to change in `_shared_config.R`

The following path blocks are external data and do not change between runs:

- `upstream_csv` (inside `DATA_PATHS`)
- `CHIP_PEAK_FILES`
- `MECP2_FILES`
- `ATAC_FILES`
- `K119UB_FILES`
- `H3K27AC_FILES`
- `DIFFBIND_FILES`
- `H3K36ME2_FILES`
- `H3K36ME3_FILES`
- `LOOP_FILES`

### Step 4: Update `process_all_dmr.sh`

**File:** `downstream/modality/scripts/process_all_dmr.sh`

One variable to change, on line 13:

```bash
# Before
INPUT_DIR="$BASE_DIR/outputs/run-5"

# After
INPUT_DIR="$BASE_DIR/outputs/run-N"
```

This script auto-discovers DMR directories and uses the latest timestamp per region,
so no timestamp changes are needed here.

### Step 5: Update `visualize_dmr_results.R`

**File:** `downstream/modality/scripts/visualize_dmr_results.R`

One hardcoded run reference in the `load_dmr_data()` function, on line 262:

```r
# Before
results_dir <- file.path(BASE_DIR, "outputs", "run-5", paste0("outputs_", context), "Results",

# After
results_dir <- file.path(BASE_DIR, "outputs", "run-N", paste0("outputs_", context), "Results",
```

This function auto-discovers the latest `DMR_*` directory within each region, so no
timestamp changes are needed.

### Step 6: Update `compare_shallow_vs_deep.R`

**File:** `downstream/scripts/viz_sections/compare_shallow_vs_deep.R`

This script compares two runs side-by-side. When a new run becomes the primary, the
previous primary run becomes the comparison baseline. Three areas need updating:

#### 6a. Rotate `SHALLOW_PATHS`

The `SHALLOW_PATHS` block (lines 17-30) should now point to the previous primary run.
Copy the timestamps that were in `DATA_PATHS` (from `_shared_config.R`) before you
updated them in Step 3, and paste them here with the old run number.

For example, if run-5 was primary and run-6 is the new primary:

```r
# SHALLOW_PATHS should now point to run-5 (the old primary)
SHALLOW_PATHS <- list(
  mc_dmr = file.path(BASE_DIR, "modality/outputs/run-5/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260402_191818/DMR_mc_control__mutant_20260402_191818.bed"),
  # ... (all 12 entries with run-5 timestamps)
)
```

`DEEP_PATHS` is set to `DATA_PATHS` from `_shared_config.R` and does not need a
manual change -- it picks up the new run automatically after Step 3.

#### 6b. Update Labels and Colors

Search the file for user-visible text that references specific run numbers or
comparison descriptions. At minimum, update:

- **Script header comment** (line 2): Description of which runs are being compared
- **Log banner** (line 10): The `cat()` statement identifying the comparison
- **`RUN_COLORS`** (line 61): Named vector keys used in plot legends
- **Validation log messages** (lines 39-58): `cat()` statements naming run numbers
- **Data loading messages** (lines 71, 91-94): `cat()` and `sprintf()` run labels
- **`dmr_counts$Run` factor levels** (line 156): The factor level strings
- **Plot titles and legends**: Any `ggtitle()`, `labs()`, or annotation text

#### 6c. Verification

After editing, confirm that `SHALLOW_PATHS` references the OLD run and `DEEP_PATHS`
(via `DATA_PATHS`) references the NEW run.

### Step 7: Verify No Stale References

Run a grep across all downstream scripts to check for leftover references to old run
numbers:

```bash
grep -rn 'run-[0-9]' \
  downstream/scripts/viz_sections/ \
  downstream/modality/scripts/ \
  --include='*.R' --include='*.sh'
```

Every match should either reference the new primary run or be an intentional reference
to a historical run (e.g., `SHALLOW_PATHS` in `compare_shallow_vs_deep.R`).

Also verify that no section scripts have hardcoded run paths. Section scripts should
only access data through the centralized config:

```bash
grep -rn 'run-[0-9]' downstream/scripts/viz_sections/section_*.R
```

This command should return zero results. If any section script contains a hardcoded
run path, move it to `_shared_config.R`.

## Update Checklist

Use this checklist when performing the update:

- [ ] rsync results from Expanse to `downstream/modality/outputs/run-N/`
- [ ] Record all timestamps (BioQC, Extract, DMR per region, CHG)
- [ ] Note the BioQC JSON filename (sample count may differ)
- [ ] **`_shared_config.R`**: Update `DATA_PATHS` (13 entries: run number + timestamps)
- [ ] **`_shared_config.R`**: Update `EXTRACT_PATHS` (3 entries: run number + timestamps)
- [ ] **`_shared_config.R`**: Update `CHG_DATA_PATHS` (2 entries: run number + timestamps)
- [ ] **`_shared_config.R`**: Confirm `upstream_csv` and all ChIP/ATAC/DiffBind paths are unchanged
- [ ] **`process_all_dmr.sh`**: Update `INPUT_DIR` (line 13)
- [ ] **`visualize_dmr_results.R`**: Update `run-N` in `load_dmr_data()` (line 262)
- [ ] **`compare_shallow_vs_deep.R`**: Rotate old `DATA_PATHS` timestamps into `SHALLOW_PATHS`
- [ ] **`compare_shallow_vs_deep.R`**: Update all labels, colors, and log messages
- [ ] Run `grep -rn 'run-[0-9]'` verification (see Step 7)
- [ ] Run one section script to confirm paths resolve: `cd downstream/ && Rscript scripts/viz_sections/section_04_volcano_plots.R`
- [ ] Update `CLAUDE.md` run history and "current primary" references if needed

## Frequent Mistakes

**Mismatched timestamps between directory and filename.** The DMR directory name and
the BED filename inside it share the same timestamp. If a region was re-run, the
directory and file timestamps will differ from the first run. Always verify by listing
the actual files:

```bash
ls downstream/modality/outputs/run-N/outputs_CG/Results/gencode.vM25.mouse.genes.annotation/DMR_*/
```

**Forgetting the BioQC sample count.** The BioQC JSON filename includes the number of
samples (e.g., `8_samples`). If the new run has a different sample count, the filename
pattern changes.

**Not rotating `compare_shallow_vs_deep.R` correctly.** Save the old `DATA_PATHS`
timestamps before overwriting them in `_shared_config.R`, or retrieve them from git
history.

**CHG paths forgotten.** The CHG context has its own timestamp (different from CG) and
its own path block in `_shared_config.R`. It is easy to overlook since it produces
minimal signal.

**`CLAUDE.md` not updated.** The `downstream/CLAUDE.md` and `biomodal/CLAUDE.md` files
reference the current primary run in their documentation. Update these to avoid
confusing future contributors (human or AI).
