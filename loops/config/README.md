# Configuration Files

This directory contains configuration files for the Mariner Hi-C differential loop analysis pipeline.

## Files

### `paths_config.yaml`

**Purpose:** Centralized configuration for all external data file paths (BEDPE files, .hic files, ChIP-seq peaks).

**Why it exists:** When you receive new deep sequencing data, you only need to update this one file instead of modifying multiple scripts.

### `edgeR_config.yaml`

**Purpose:** Statistical analysis parameters and edgeR settings.

**When to modify:** When changing analysis parameters, significance thresholds, or visualization settings.

---

## Updating Paths for New Data

When you receive new Hi-C data, follow these steps to update the pipeline:

### Step 1: Update `paths_config.yaml`

Open `config/paths_config.yaml` and update the following sections:

#### 1. HiCCUPS Loop Calls (BEDPE files)

```yaml
hiccups_loops:
  base_path: "/path/to/your/hiccups_results"  # ← Update this path
  samples:  # ← Update sample names if different
    - ctrl_M1
    - ctrl_M2
    - ctrl_M3
    - mut_M1
    - mut_M2
    - mut_M3
```

**Directory structure expected:**
```
/path/to/your/hiccups_results/
  ├── ctrl_M1/
  │   └── postprocessed_pixels_5000.bedpe
  │   └── postprocessed_pixels_10000.bedpe
  │   └── postprocessed_pixels_25000.bedpe
  ├── ctrl_M2/
  │   └── postprocessed_pixels_5000.bedpe
  │   └── ...
  └── ...
```

#### 2. Hi-C Contact Matrices (.hic files)

```yaml
hic_files:
  base_path: "/path/to/your/hic_files"  # ← Update this path if using pattern

  # OR update individual paths:
  ctrl_M1: "/path/to/ctrl_M1.hic"  # ← Update these
  ctrl_M2: "/path/to/ctrl_M2.hic"
  ctrl_M3: "/path/to/ctrl_M3.hic"
  mut_M1: "/path/to/mut_M1.hic"
  mut_M2: "/path/to/mut_M2.hic"
  mut_M3: "/path/to/mut_M3.hic"
```

#### 3. ChIP-seq Peak Files (optional - for anchor classification)

```yaml
chipseq_peaks:
  h3k27ac: "path/to/H3K27ac_peaks.bed"  # ← Update if you have new ChIP-seq
  h3k4me1: "path/to/H3K4me1_peaks.bed"  # ← Update if you have new ChIP-seq
```

**Note:** These can be relative paths (from the project base directory) or absolute paths.

#### 4. Project Base Directory

```yaml
project:
  base_dir: "/path/to/your/project"  # ← Update if working directory changes
```

#### 5. Sample Metadata

```yaml
samples:
  names: ["ctrl_M1", "ctrl_M2", "ctrl_M3", "mut_M1", "mut_M2", "mut_M3"]
  groups: ["ctrl", "ctrl", "ctrl", "mut", "mut", "mut"]

  descriptions:  # ← Update these descriptions
    ctrl_M1: "Control replicate 1 (wildtype)"
    # ... etc
```

### Step 2: Validate Your Configuration

After updating `paths_config.yaml`, verify that all paths are correct:

```bash
# Check that BEDPE files exist
ls /path/to/your/hiccups_results/ctrl_M1/postprocessed_pixels_5000.bedpe

# Check that .hic files exist
ls /path/to/ctrl_M1.hic

# Check that ChIP-seq files exist (if using)
ls path/to/H3K27ac_peaks.bed
```

### Step 3: Run the Pipeline

Once configured, run the pipeline as usual:

```bash
# Full pipeline
sbatch scripts/run_full_pipeline.sb

# Or step-by-step
Rscript scripts/prep_loops.R 5000
Rscript scripts/extract_counts.R 5000
# ... etc
```

The scripts will automatically read the paths from `config/paths_config.yaml`.

---

## Example: Updating for New Dataset

### Scenario

You received new Hi-C data:
- Dataset name: `250601_Bap1_newdata`
- Location: `/expanse/lustre/projects/csd940/ctea/nf-hic/250601_Bap1_newdata`
- Same sample names (ctrl_M1, ctrl_M2, etc.)

### Changes Needed

Edit `config/paths_config.yaml`:

```yaml
# Old:
hic_files:
  base_path: "/expanse/lustre/projects/csd940/ctea/nf-hic/250402_Bap1_deepseq"
  ctrl_M1: "/expanse/lustre/projects/csd940/ctea/nf-hic/250402_Bap1_deepseq/trimmed_ctrl_M1/juicer/ctrl_M1.hic"
  # ... etc

# New:
hic_files:
  base_path: "/expanse/lustre/projects/csd940/ctea/nf-hic/250601_Bap1_newdata"
  ctrl_M1: "/expanse/lustre/projects/csd940/ctea/nf-hic/250601_Bap1_newdata/trimmed_ctrl_M1/juicer/ctrl_M1.hic"
  # ... etc
```

Then update the metadata:

```yaml
metadata:
  last_updated: "2025-11-12"  # ← Update date
  dataset: "BAP1_mutant_newdata"  # ← Update dataset name
```

That's it! No need to modify any scripts.

---

## Troubleshooting

### Error: "File not found"

**Cause:** Path in config is incorrect.

**Fix:**
1. Check the exact error message to see which file is missing
2. Verify the path exists: `ls /path/to/file`
3. Update `config/paths_config.yaml` with the correct path

### Error: "Cannot read YAML"

**Cause:** Syntax error in YAML file.

**Fix:**
1. Check indentation (YAML is sensitive to spaces)
2. Ensure strings with special characters are quoted
3. Validate YAML syntax online: http://www.yamllint.com/

### Scripts still use old paths

**Cause:** Script not updated to use config.

**Fix:** Verify these scripts are using the config:
- prep_loops.R
- extract_counts.R
- aggregate.R
- qc-val.R
- downstream_analysis.R
- visualizations.R
- apa_analysis.R
- compare_resolutions.R

Each should have:
```R
library(yaml)
config <- yaml::read_yaml("config/paths_config.yaml")
```

---

## Which Scripts Use Which Config Settings?

| Script | Uses |
|--------|------|
| `prep_loops.R` | `hiccups_loops.base_path`, `hiccups_loops.samples` |
| `extract_counts.R` | `hic_files.*`, `project.base_dir` |
| `aggregate.R` | `project.base_dir` |
| `qc-val.R` | `project.base_dir` |
| `edgeR.R` | `edgeR_config.yaml` (different file) |
| `compare_resolutions.R` | `project.base_dir` |
| `downstream_analysis.R` | `chipseq_peaks.*`, `project.base_dir` |
| `visualizations.R` | `project.base_dir` |
| `apa_analysis.R` | `hic_files.*` |

---

## Questions?

If you have questions about updating paths or encounter issues:

1. Check this README
2. Look at comments in `config/paths_config.yaml`
3. Refer to the main project documentation in `CLAUDE.md`
4. Check the error messages carefully - they usually indicate which file is missing
