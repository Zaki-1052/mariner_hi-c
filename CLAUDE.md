# edgeR Pipeline Modification: Biological Replicate Integration

## Context
You are modifying an existing edgeR differential loop analysis pipeline that currently uses merged/consensus Hi-C samples (n=1 per condition) to instead use the actual biological replicates (n=3 per condition). The current approach only identifies 2 significant differential loops while Hiccups found 38,593, likely due to the lack of statistical power without replicates.

## Required Documents
Please review these documents in order:
1. **edgeR-replicates.md** - Current issue documentation and data structure
2. **edgeR-prep.md** - Original pipeline context and data locations
3. **prepare_loops.R** and **05_extract_counts.R** scripts from **mariner.md** - Current implementation
4. **edgeR-doc.md** - edgeR documentation for proper replicate handling
5. **config.md** - Current configuration that needs updating
6. **guidelines.md** - Coding standards to follow

Reference:
1. @edgeR-replicates.md
2. @edgeR-prep.md
3. @edgeR-PRD.md
4. @mariner.md
5. @mariner.Rmd
6. @config.md
7. @edgeR-doc.md
8. @guidelines.md

## Objective
Modify the Mariner + edgeR pipeline to use individual biological replicates (ctrl_M1, ctrl_M2, ctrl_M3, mut_M1, mut_M2, mut_M3) instead of merged samples, enabling proper variance estimation and increasing statistical power for differential loop detection.

## Key Requirements

### 1. Data Structure Changes
- Input: 6 BEDPE files from individual replicates (paths in edgeR-replicates.md)
- Output: 22,108 loops × 6 samples count matrix
- Maintain the 5×5 buffer approach for handling positional shifts

### 2. Statistical Improvements
- Enable edgeR to estimate dispersions from biological replicates
- Remove all fixed BCV assumptions (0.10, 0.15, 0.20)
- Implement the recommended QL pipeline with robust=TRUE
- Add replicate quality checks (correlation, MDS plots)

### 3. Implementation Plan
Create or modify these scripts in order:

**Script 1: prepare_loops_replicates.R**
- NOTE: This is: @scripts/prep_loops.R
- Read all 6 individual replicate BEDPE files
- Implement consensus loop calling across all samples
- Consider both union (any sample) and intersection (all samples) approaches
- Output merged GInteractions with all 6 samples tracked

**Script 2: 05_extract_counts_replicates.R**
- NOTE: This is: @scripts/extract_counts.R
- Extract Hi-C matrices for all 6 samples
- Handle case where individual .hic files might not exist (check first!)
- If no individual .hic files, propose alternative approach

**Script 3: 01_run_edgeR_replicates.R**
- NOTE: This is: @scripts/edgeR.R
- Implement proper replicate-aware edgeR analysis
- Use estimateDisp() for data-driven dispersion estimation
- Add MDS plot for sample relationships
- Add BCV plot to visualize dispersion estimates
- Compare results to original 2-sample analysis

**Script 4: config/edgeR_replicates_config.yaml**
- NOTE: This is: @config/edgeR_config.yaml
- Update configuration for 6-sample structure
- Remove fixed BCV parameters
- Add replicate-specific QC thresholds

## Expected Outcomes
- Substantially more differential loops detected (target: hundreds to thousands)
- Better concordance with Hiccups results
- Proper statistical confidence from biological replication
- Reproducible variance estimates

## Critical Considerations
1. **Check if individual .hic files exist** before assuming extraction strategy
   - Not needed, I've verified for you
2. **Preserve the buffer approach** - still needed for positional shift handling
3. **Maintain backward compatibility** - keep ability to compare with merged approach
4. **Follow modular design** from guidelines.md with clear function separation
5. **Add comprehensive error handling** for missing files or failed extractions

## Deliverables
Please provide:
1. Modified scripts with clear comments explaining changes
2. Updated configuration file
3. Brief explanation of key architectural decisions
4. Suggested SLURM job script for running the full pipeline
5. Expected computational requirements (time/memory) for 6-sample analysis

Focus on getting the core replicate analysis working first. We can refine the loop merging strategy and comparison analyses in subsequent iterations.