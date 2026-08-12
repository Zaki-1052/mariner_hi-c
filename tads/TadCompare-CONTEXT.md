# Context Document for TADCompare Analysis: Control vs BAP1-KO Mutant Hi-C Comparison

---

## 1. Project Overview

### Biological Context
This is a Hi-C analysis project comparing **control (wildtype) vs BAP1 KO mutant** mouse cerebellum samples. BAP1 is a histone deubiquitinase, and the goal is to characterize 3D genome organization changes at the TAD (Topologically Associating Domain) level.

**Key biological findings from prior analyses:**
- Compartments are weakening in mutant (±0.25 range), not undergoing overt A/B switching
- The syt1/nav3 locus on chromosome 10 shows high insulation score changes and is one of the most impacted regions
- Long-range loops appear to be lost in mutant
- Developmental genes show evidence of rewiring
- Approximately 900 TAD-level differences were previously identified using a consensus TAD set (via Homer)

### Experimental Design
| Parameter | Value |
|-----------|-------|
| Organism | Mouse (mm10 genome) |
| Protocol | Arima Hi-C |
| Conditions | Control (wildtype) vs BAP1-KO (mutant) |
| Replicates | 3 biological replicates per condition |
| Sequencing depth | ~300-600 million valid pairs per sample |
| Timepoints | Late (250402) and Early (250831) |

**Sample naming convention:**
- `ctrl_M1`, `ctrl_M2`, `ctrl_M3`: Control/wildtype replicates
- `mut_M1`, `mut_M2`, `mut_M3`: BAP1-KO mutant replicates

---

## 2. Analysis Scope and Goals

### Current Task
Run TADCompare to identify and classify differential TAD boundaries between control and BAP1-KO mutant conditions, starting with the **late timepoint (250402)**.

### Specific Objectives
1. Compare TAD boundaries between merged control and merged mutant contact matrices
2. Classify boundary changes using TADCompare's schema (Shifted, Split, Merged, Strength Change, Complex, Non-Differential)
3. Calculate exact shift distances (in bp/kb) for boundaries that moved
4. Generate summary statistics of boundary changes
5. Export results in BED format for IGV visualization and downstream integration

### Downstream Integration (Context Only)
The TAD boundary analysis will eventually be paired with:
- CTCF and RAD21 ChIP-seq data
- Histone modifications: H2AK119ub, H3K27ac, H3K27me3
- Loop-level analyses (already being done separately with Mariner)

---

## 3. Key Decisions Made

These decisions were explicitly confirmed during planning discussions:

| Decision | Choice | Rationale |
|----------|--------|-----------|
| **Resolution** | 40kb | Matches existing HiCExplorer TAD calls; standard resolution for TAD-level analysis |
| **TAD boundary detection** | Let TADCompare detect boundaries (do NOT use `pre_tads` argument) | Existing Homer consensus TAD file is cross-condition (not per-condition as `pre_tads` requires); internal consistency with TADCompare's classification schema |
| **Primary analysis input** | Merged .hic files (`ctrl_merged.hic`, `mut_merged.hic`) | Aggregates replicates for higher signal |
| **Robustness check** | Run ConsensusTADs on individual replicates | Identifies which boundaries are consistently detected across replicates |
| **Shift distance calculation** | Yes, required | Use `bedtools closest` to calculate exact distances for boundaries classified as "Shifted" |
| **Chromosomes** | All autosomes | Exclude chrX (per existing filtering conventions) |
| **Output format** | BED format | For loading into IGV; retain all metadata |

---

## 4. File Paths and Data Locations

### Computing Environment
- **HPC System**: SDSC Expanse
- **Account**: csd940
- **User directory**: `/expanse/lustre/projects/csd940/zalibhai/`
- **Job scheduler**: SLURM

### Input Files: .hic Contact Matrices

**Location**: `/expanse/lustre/projects/csd940/zalibhai/stripes/StripeCaller/data/hic/250402/`

**Merged files (for primary TADCompare analysis):**
- `ctrl_merged.hic`
- `mut_merged.hic`

**Individual replicate files (for ConsensusTADs robustness check):**
- `ctrl_M1.hic`, `ctrl_M2.hic`, `ctrl_M3.hic`
- `mut_M1.hic`, `mut_M2.hic`, `mut_M3.hic`

**File format**: Juicer v2.0 .hic format (created via juicerpre from nf-core/hic pipeline output)

### Existing TAD Calls (For Reference/Validation Only - NOT for pre_tads input)

**HiCExplorer calls** (per-replicate, at 10kb/20kb/40kb):
```
/expanse/lustre/projects/csd940/ctea/nf-hic/250402_Bap1_deepseq/trimmed_{sample}/tads/hicExplorer/
├── {sample}.10000_balanced_hicfindtads_boundaries.bed
├── {sample}.10000_balanced_hicfindtads_domains.bed
├── {sample}.20000_balanced_hicfindtads_boundaries.bed
├── {sample}.20000_balanced_hicfindtads_domains.bed
├── {sample}.40000_balanced_hicfindtads_boundaries.bed
├── {sample}.40000_balanced_hicfindtads_domains.bed
└── ... (score files, etc.)
```

**Homer calls** (consensus across all samples):
```
/expanse/lustre/projects/csd940/ctea/homer/juicer_bam/
├── merged.tad.2D.bed          # Consensus TAD set across all replicates
├── Bap1.diff.tad.txt          # Differential TAD analysis output
└── BAP1.tad.scores.txt        # TAD inclusion ratio scores per replicate
```

**Homer file format example** (`merged.tad.2D.bed`):
```
#merged=5074 chr1    start1  end1    chr2    start2  end2    info    info    info
chrX    38060999    38183998    chrX    38060999    38183998    255,255,0    1.525    1.525
chr5    68006999    68120998    chr5    68006999    68120998    255,255,0    1.502    1.502
...
```

---

## 5. TADCompare: What It Does and Does Not Provide

### What TADCompare Outputs

Based on the TADCompare vignette (TADCompare.Rmd), the main output is a data frame (`TAD_Frame`) with these columns:

| Column | Description |
|--------|-------------|
| `Boundary` | Genomic coordinate of the boundary |
| `Gap_Score` | Differential boundary score (Z-score of difference between boundary scores) |
| `TAD_Score1` | Boundary score in matrix 1 (control) |
| `TAD_Score2` | Boundary score in matrix 2 (mutant) |
| `Differential` | "Differential" or "Non-Differential" |
| `Enriched_In` | Which matrix contains the stronger boundary ("Matrix 1" or "Matrix 2") |
| `Type` | Classification: Split, Merge, Shifted, Strength Change, Complex, or Non-Differential |

**Additional outputs:**
- `Boundary_Scores`: Same info but for ALL regions (not just boundaries)
- `Count_Plot`: ggplot2 stacked barplot of boundary type prevalence

### What TADCompare Does NOT Provide
- **Exact shift distances**: TADCompare classifies a boundary as "Shifted" but does not report how many bp/kb it shifted
- **Distribution statistics**: e.g., "X% of boundaries shifted by Y kb"
- **TAD size change metrics**: Not directly calculated

These metrics require post-processing after TADCompare runs.

### Shift Distance Calculation Approach
To obtain exact shift distances for boundaries classified as "Shifted":
1. Extract control boundaries and mutant boundaries from TADCompare output
2. Use `bedtools closest -d` to find nearest boundary in other condition and report distance
3. This is a standard genomic operation, not a custom analysis

---

## 6. Input Data Format for TADCompare

TADCompare accepts three input formats (see Input_Data.Rmd vignette):
1. **n × n contact matrices** (fastest)
2. **n × (n+3) matrices** (TopDom style)
3. **Sparse 3-column matrices** (region_i, region_j, contacts)

**To extract from .hic files**, use Juicer's `straw` tool:
```bash
straw NONE sample.hic chr chr BP 40000 > contacts_sample.40kb.txt
```
This produces sparse 3-column format that TADCompare can convert internally.

**Resolution note**: When extracting at 40kb, specify `40000` as the resolution parameter to straw.

---

## 7. Analysis Workflow Components

The following components were discussed and confirmed. Implementation details are to be determined by the implementing party.

### Component 1: TADCompare on Merged Matrices
- **Input**: Contact matrices extracted from `ctrl_merged.hic` and `mut_merged.hic` at 40kb resolution
- **Output**: Differential boundary classifications, boundary scores

### Component 2: ConsensusTADs on Individual Replicates
- **Input**: Contact matrices from all 6 individual .hic files at 40kb resolution
- **Output**: Consensus boundary scores identifying robustly-detected boundaries across replicates
- **Purpose**: Robustness/confidence filter for which boundaries to trust

### Component 3: Shift Distance Quantification
- **Input**: Boundary coordinates from TADCompare output
- **Method**: `bedtools closest` to calculate distances
- **Output**: Shift distances for "Shifted" boundaries

### Component 4: Summary and Export
- **Output format**: BED files suitable for IGV
- **Metadata to retain**: Boundary type, shift distance, boundary scores, enriched condition

---

## 8. Output Requirements

### BED File Contents
Each boundary should include:
- Genomic coordinates (chr, start, end)
- Boundary type classification (Shifted/Split/Merged/etc.)
- Shift distance (where applicable)
- Boundary scores from both conditions
- Which condition the boundary is enriched in
- Differential status

### Summary Statistics Needed
- Count and percentage of boundaries in each classification category
- Distribution of shift distances for "Shifted" boundaries
- Overall summary of boundary changes

---

## 9. Reference Documents

The following vignette documents are provided separately and contain detailed usage examples and function documentation:

1. **TADCompare.Rmd** - Main vignette covering:
   - `TADCompare()` function for differential boundary detection
   - `TimeCompare()` function for time-course analysis
   - `ConsensusTADs()` function for consensus boundary calling
   - `DiffPlot()` function for visualization
   - Pre-specification of TADs option
   - Output interpretation

2. **Input_Data.Rmd** - Input format specifications:
   - n × n matrix format
   - n × (n+3) matrix format
   - Sparse 3-column format
   - Working with .hic files (straw extraction)
   - Working with .cool files
   - Runtime benchmarks by format

3. **README.md** - Package overview and quick-start examples

---

## 10. Important Notes and Constraints

1. **Timepoint clarification**: "250402" = late timepoint (doing first); "250831" = early timepoint (doing second after late is confirmed working)

2. **Do NOT use pre_tads**: The existing Homer consensus TAD file (`merged.tad.2D.bed`) is a single consensus across ALL samples (ctrl + mut combined), not separate per-condition files. This does not match the `pre_tads` expected format (list of two BED files, one per condition).

3. **Chromosome filtering**: Exclude chrX from analysis (consistent with other pipelines in this project)

4. **The .hic files in the stripes directory** are the same as those generated by the main pipeline, with KR normalization applied. Either location can be used.

5. **Validation opportunity**: Results can be compared against existing Homer differential TAD calls (`Bap1.diff.tad.txt`) and HiCExplorer boundaries as a sanity check.

---

## 11. Prompt for Implementation

You have been provided with:
1. This context document (single source of truth for project decisions and file locations)
2. TADCompare.Rmd - the main package vignette
3. Input_Data.Rmd - input format specifications
4. README.md - package overview

**Your task**: Analyze all provided context and plan the implementation to accomplish the TADCompare analysis for comparing control vs BAP1-KO mutant TAD boundaries using the late timepoint (250402) data.

Consider:
- How to extract contact matrices from the .hic files at 40kb resolution
- Running TADCompare on merged matrices
- Running ConsensusTADs on individual replicates
- Calculating shift distances for shifted boundaries
- Generating summary statistics and BED output files
- Any data preprocessing or format conversions needed

The analysis will run on SDSC Expanse (SLURM environment). R is available.