# Plan: Detect Non-CG (CA-context) Methylation from EM-seq and Compare with evoC

## Context

We want to compare non-CG methylation (specifically CpA / "CA context") between EM-seq and evoC to determine whether EM-seq detects more non-CG signal. The EM-seq pipeline already extracts CHH and CHG methylation via MethylDackel (`--CHH --CHG`), producing methylKit files. The evoC pipeline reports genome-wide CHH ~0.86-0.89% modC and CHG ~0.63-0.66% modC, with zero significant DMRs in CHG and only noise-level in CHH. EM-seq M-bias suggests CHH background of ~0.9-1.6%. The question is whether EM-seq captures more real mCA signal (biologically relevant in P60 cerebellum) or simply has a higher false-positive rate.

**Key challenge**: The methylKit format doesn't include trinucleotide context, so filtering for CA requires the reference genome. CA filtering will run on Expanse via SLURM.

**Scope**: Genome-wide aggregate rates AND per-gene comparisons (with coverage caveats for EM-seq at ~1x).

---

## Step 1: Create analysis directory and scaffold

Create `non_cg_analysis/` under `/Users/zakiralibhai/sdsc/emseq/`:
- `non_cg_methylation.py` — Main local analysis script (Steps 2, 3, 5, 6)
- `ca_filter.py` — Cluster-side CA-context filtering script (Step 4)
- `run_ca_filter.sb` — SLURM wrapper for `ca_filter.py`

## Step 2: Genome-wide non-CG rates from EM-seq (local)

Parse existing methylKit files to compute genome-wide CHH and CHG methylation.

**Input files:**
- `em-seq_output/methylDackelExtracts/250504_A1_Bap1_P60_ctrl1_5mC_S44_L002_CHH.methylKit.gz`
- `em-seq_output/methylDackelExtracts/250504_A1_Bap1_P60_ctrl1_5mC_S44_L002_CHG.methylKit.gz`
- `em-seq_output/methylDackelExtracts/250504_B1_Bap1_P60_ctrl2_5mC_S45_L002_CHH.methylKit.gz`
- `em-seq_output/methylDackelExtracts/250504_B1_Bap1_P60_ctrl2_5mC_S45_L002_CHG.methylKit.gz`

**Format**: `chrBase  chr  base  strand  coverage  freqC  freqT`

**Logic:**
```
For each file (CHH, CHG x ctrl1, ctrl2):
  Stream gzipped TSV, skip header
  Separate by chromosome category:
    - Autosomes (chr1-chr19)
    - Sex chromosomes (chrX, chrY)
    - Spike-ins (phage_lambda, plasmid_puc19c, phage_T4, phage_Xp12)
  Accumulate per-category:
    - total_meth_reads = sum(coverage * freqC / 100)
    - total_coverage = sum(coverage)
    - avg_methylation = total_meth_reads / total_coverage * 100
  Also compute per-chromosome breakdown
```

**Output**: Summary table with genome-wide CHH/CHG methylation rates per sample, plus spike-in false-positive rates.

## Step 3: evoC baseline extraction (local)

Extract evoC non-CG methylation from existing local data.

**Source 1 — Genome-wide rates from evoC summary CSV:**
- Path: `biomodal/upstream/duet-1.5.0_evoC_Bap1_run_6bp/reports/evoC_Bap1_run_duet-evoC_Summary.csv`
- Columns: `modality_summary_chh_autosomes_modc` (~0.86-0.89%), `modality_summary_chg_autosomes_modc` (~0.63-0.66%)

**Source 2 — Per-gene means from evoC DMR BED files:**
- CHH: `biomodal/downstream/modality/outputs/run-2/outputs_CHH/Results/gencode.vM25.mouse.genes.annotation/DMR_20260121_184551/DMR_mc_control__mutant_20260121_184551.bed`
- CHG: `biomodal/downstream/modality/outputs/run-2/outputs_CHG/Results/gencode.vM25.mouse.genes.annotation/DMR_20260121_174629/DMR_mc_control__mutant_20260121_174629.bed`
- Format: 14-col BED with `mean_mod_group_1` (control mean), `mean_mod_group_2` (mutant mean), `Name` (gene name)
- Use `mean_mod_group_1` as the control-sample evoC non-CG rate per gene

## Step 4: CA-context filtering (Expanse cluster)

**`ca_filter.py`** — Reads methylKit files, looks up each position in the reference genome, writes CA-only output + summary stats.

**Input:**
- EM-seq methylKit CHH and CHG files (copy to cluster or access from work dir)
- Reference: `/expanse/lustre/projects/csd940/zalibhai/refs/mm10+controls.fa`

**Logic:**
```python
import pysam
fasta = pysam.FastaFile(ref_path)

for each row in methylKit:
    chr, pos (1-based), strand = row['chr'], row['base'], row['strand']
    if strand == 'F':
        # C is at pos (1-based), next base is at pos (0-based in pysam)
        next_base = fasta.fetch(chr, pos, pos + 1).upper()
        is_CA = (next_base == 'A')
    elif strand == 'R':
        # C on reverse strand at pos means G on forward at pos
        # The base before G on forward strand is at pos-2 (0-based)
        prev_base = fasta.fetch(chr, pos - 2, pos - 1).upper()
        is_CA = (prev_base == 'T')  # T on forward = A on reverse (the base after C)
```

**Output:**
- `*_CA_CHH.methylKit.gz` and `*_CA_CHG.methylKit.gz` — CA-only sites
- `ca_summary_stats.tsv` — Genome-wide CA methylation rates, site counts, per-chromosome breakdown

**`run_ca_filter.sb`** — SLURM script:
```bash
#!/bin/bash
#SBATCH --job-name=ca_filter
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --account=csd940

module load cpu/0.17.3b
conda activate emseq  # or appropriate env with pysam

python ca_filter.py \
  --ref /expanse/lustre/projects/csd940/zalibhai/refs/mm10+controls.fa \
  --input-dir /path/to/methylKit/files \
  --output-dir ./ca_filtered/
```

## Step 5: Comparison analysis (local)

### 5a. Genome-wide rates table
Compare EM-seq (ctrl1, ctrl2) vs evoC (ctrl-F, ctrl-M) for:
- Overall CHH and CHG methylation (%)
- CA-specific methylation (once Step 4 complete)
- Lambda false-positive rate

### 5b. Per-chromosome breakdown
- EM-seq per-autosome CHH/CHG rates
- Check if non-CG methylation is uniform or shows chromosome-specific enrichment

### 5c. Coverage-stratified analysis
- EM-seq sites at coverage >=1, >=2, >=3, >=5
- If methylation rate changes with coverage threshold, low-coverage sites may be noisy

### 5d. Per-gene comparison
- Aggregate EM-seq CHH/CHG methylation per gene body using gencode.vM25 gene BED
  - Gene BED: `biomodal/downstream/modality/mm10/gencode.vM25.mouse.genes.annotation.bed.gz`
  - For each gene: intersect EM-seq sites, compute weighted mean methylation
- Compare with evoC per-gene means from DMR BED `mean_mod_group_1` column
- Scatter plot: EM-seq per-gene % vs evoC per-gene % (expect both near 0, look for outliers)
- Caveat: EM-seq at ~1x will have very sparse per-gene coverage; include number of sites per gene

## Step 6: Spike-in control analysis (local, part of main script)

**Lambda phage (unmethylated DNA):**
- Extract CHH/CHG methylation from `phage_lambda`-mapping rows in methylKit files
- Any methylation detected = false positive (incomplete conversion)
- QC report reference: lambda nonconverted reads = 1.07% (ctrl1), 2.27% (ctrl2)
- Note: read-level nonconversion rate != site-level methylation rate

**pUC19 (CpG-methylated, non-CG unmethylated):**
- CHH/CHG methylation on `plasmid_puc19c` should be ~0%
- Second independent false-positive estimate

**evoC comparison:**
- evoC DQS reports: `dqs_lambda_mc_sensitivity` ~0.987 (controls detect 98.7% of true methylation)
- Different metric (sensitivity vs false-positive), but both speak to conversion accuracy

## Step 7: Output

1. **Console summary table**: Genome-wide non-CG rates (EM-seq vs evoC) with spike-in false-positive estimates
2. **`non_cg_comparison.png`**: Bar chart of CHH, CHG, and CA methylation by method/sample
3. **`per_gene_scatter.png`**: Per-gene EM-seq vs evoC non-CG methylation (with coverage >= threshold filter)
4. **`per_chromosome_rates.png`**: Per-autosome non-CG rates for EM-seq samples
5. **`non_cg_summary.tsv`**: Machine-readable output of all computed rates

---

## Files to create

| File | Location | Runs on |
|------|----------|---------|
| `non_cg_methylation.py` | `non_cg_analysis/` | Local |
| `ca_filter.py` | `non_cg_analysis/` | Expanse |
| `run_ca_filter.sb` | `non_cg_analysis/` | Expanse |

## Existing files to read

| File | Purpose |
|------|---------|
| `em-seq_output/methylDackelExtracts/*_CHH.methylKit.gz` | EM-seq CHH per-site calls |
| `em-seq_output/methylDackelExtracts/*_CHG.methylKit.gz` | EM-seq CHG per-site calls |
| `biomodal/upstream/.../evoC_Bap1_run_duet-evoC_Summary.csv` | evoC genome-wide stats |
| `biomodal/downstream/modality/outputs/run-2/outputs_CHH/.../DMR_mc_*.bed` | evoC CHH per-gene |
| `biomodal/downstream/modality/outputs/run-2/outputs_CHG/.../DMR_mc_*.bed` | evoC CHG per-gene |
| `biomodal/downstream/modality/mm10/gencode.vM25.mouse.genes.annotation.bed.gz` | Gene body coordinates |

## Dependencies

- **Local**: Python 3.8+, `pandas`, `numpy`, `matplotlib`
- **Cluster**: Above + `pysam` (for reference genome lookup)

## Verification

1. Run `python non_cg_methylation.py` locally — produces genome-wide rates + comparison table + plots
2. Lambda false-positive rate should be consistent with QC report (~1-2% nonconversion)
3. evoC values should match summary CSV (~0.86% CHH, ~0.63% CHG)
4. Per-gene scatter should show most genes clustered near 0% for both methods
5. After running `ca_filter.py` on Expanse, re-run local script to include CA-specific comparisons
