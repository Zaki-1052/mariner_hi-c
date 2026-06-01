# B0/B1 Checkpoint: Adding KR Normalization to 250402 .hic Files

## Problem

B1 (`scripts/B1_generate_cropmap.sb`) failed on timepoint 250402. All 90 juicer_tools `dump observed KR` calls produced empty files, resulting in an all-zero matrix and 0 retained bins.

## Root Cause

The 250402 merged .hic files (`ctrl_merged.hic`, `mut_merged.hic`) were created without the `-k KR` flag. They have GW_SCALE, INTER_SCALE, SCALE, VC, VC_SQRT — but no KR normalization at any resolution. The 250831 files have KR and work fine (verified: 2M records for chr1×chr2 at 100kb).

This is a known issue documented in `lab/hic-file-auditor/header.md` Section 3 and `stripes/stripenn/README.md` Section 5.2.

## Prior Art

For the Stripenn pipeline, `juicer_tools addNorm` was tried but hictk could not read the KR vectors it wrote (Class 3 cross-tool incompatibility). However, SNIPER's B1 uses juicer_tools for both writing and reading, so the same-tool path should work.

## What Was Done

Wrote `scripts/B0_addnorm_kr.sb` — a SLURM job that:

1. Copies both 250402 merged .hic files to working copies in DATA_DIR
2. Runs `juicer_tools addNorm -k KR -r 100000 -j 8` on each copy (KR at 100kb only)
3. Verifies each copy with `juicer_tools dump norm KR` and `dump observed KR` (prints PASS/FAIL)

Submit: `sbatch scripts/B0_addnorm_kr.sb`

## Key Paths

| What | Path |
|------|------|
| Original ctrl .hic | `/expanse/lustre/projects/csd940/zalibhai/stripes/StripeCaller/data/hic/250402/ctrl_merged.hic` |
| Original mut .hic | `/expanse/lustre/projects/csd940/zalibhai/stripes/StripeCaller/data/hic/250402/mut_merged.hic` |
| KR'd ctrl copy | `/expanse/lustre/projects/csd940/zalibhai/sniper/sniper_mm10/ctrl_merged_addnorm.hic` |
| KR'd mut copy | `/expanse/lustre/projects/csd940/zalibhai/sniper/sniper_mm10/mut_merged_addnorm.hic` |
| B0 script | `ML/cmpts/scripts/B0_addnorm_kr.sb` |
| B1 script | `ML/cmpts/scripts/B1_generate_cropmap.sb` |
| B1 Python | `ML/cmpts/scripts/B1_adapt_sniper_mm10.py` |

## Next Steps (After B0 Completes)

1. Check the B0 log (`b0_addnorm_*.out`) — both files should show PASS
2. ~~Back up the originals~~ and move the KR'd copies into the original HIC_ROOT paths so B1 picks them up without code changes
3. Resubmit B1: `sbatch scripts/B1_generate_cropmap.sb 250402`

~~Steps 2 and 3 can be combined into a single SLURM job that moves the files and then runs B1.~~
User Note: No, I'd prefer to move them manually and then run myself (since I dont want to modify the B1 script).
Also no need to back up the originals if the new ones just have more normalization.