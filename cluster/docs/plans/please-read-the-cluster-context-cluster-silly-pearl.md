# Phase 5 (deepTools metagene, local) + Phase 6 (cooltools pileup, HPC)

## Context

Phases 1–4 of the Popay-style Hi-C loop clustering analysis for BAP1-KO cerebellum (late timepoint, 250402) are complete. Phase 4.4 produced the headline result: clust5 (gain) shows Polycomb at both anchor (6.59×) and span (3.03×) — Polycomb-domain expansion / extrusion impediment; clust6 (loss) shows Polycomb only at anchor (2.09×, span 0.94×) — anchor-disruption / sensitivity model. This is a stronger split than Popay's mixed-dependency cluster.

Phase 5 (deepTools metagene at anchors ±5kb) and Phase 6 (cooltools off-diagonal pileup over loop coordinates) are the remaining downstream visualizations that ground the chromatin-state result in two complementary spatial views:
- **Phase 5** answers: *what does the histone-mark signal look like centered on each cluster's anchors, locally?* — a 1D profile that complements 4.3's anchor-only mean by showing flanking shape (e.g., does K27me3 form a peak at the anchor or a plateau extending into the loop?).
- **Phase 6** answers: *what does the aggregate Hi-C contact frequency look like per cluster in ctrl vs mut?* — a 2D pileup that visually confirms loop strength changes between conditions.

User decisions (made via AskUserQuestion):
- **Phase 5 marks:** 4 marks × ctrl/mut = 8 BigWigs (K27ac, K27me3, K119ub, K27me1) sourced from `/Users/zakiralibhai/sdsc/bigwigs/`.
- **Phase 6 HPC env:** create a fresh `cluster` conda env on Expanse from `cluster/environment_linux.yml` (matches local-Mac parity, isolated from `mariner_env`).
- **Phase 6 mcools:** merged only (`ctrl_merged.mcool`, `mut_merged.mcool`); skip per-replicate (variance work is not the goal here).

No changes to existing Phase 1–4 outputs. All new artifacts go under `cluster/bap1_late/figures/deeptools_input/`, `cluster/bap1_late/figures/deeptools/`, and `cluster/bap1_late/cooltools/`.

---

## Design Overview

### Pipeline placement

```
Phase 4 outputs (combined-clusters.txt, k=6, 38,948 loops)
  │
  ├─► Phase 5 [LOCAL — Mac, cluster env]
  │     06_deeptools_metagene.py
  │       ├─ 5.1 Per-cluster anchor BEDs (3-col, deduplicated)
  │       └─ 5.2 deepTools bed_pileup → 1 heatmap (6 clusters × 8 BigWigs)
  │     run_phase5.sh (driver, tee → cluster/phase5.txt)
  │
  └─► Phase 6 [HPC — Expanse, new `cluster` env]
        07_cooltools_pileup.sb     (SLURM wrapper, account=csd940)
        07_cooltools_pileup.py     (Python orchestrator)
          ├─ Load combined-clusters.txt → 6 BEDPEs (cooler-convention cols)
          └─ cooltools_called.mcool_pileup → 2 heatmaps (ctrl + mut, 6 clusters each)
        Output writes back to repo (synced via git)
```

### Numbering convention

- Script numbers continue the existing sequence: `01–05` already used by Phases 1–4. Phase 5 → `06_deeptools_metagene.py`. Phase 6 → `07_cooltools_pileup.{py,sb}`. (plan-p2 originally reserved `06_` for cooltools; renumber to keep numbers monotonic with phase progression.)
- Drivers: `run_phase5.sh` (local). Phase 6 has no local driver — it is `sbatch`-submitted from Expanse.

### Reuse audit (no new utility code where existing modules suffice)

| Need | Reuse from |
|------|-----------|
| Cluster loading + group_dict | `cluster/scripts/05_grouped_analyses.py:135–148` (`load_clusters`) |
| Multi-format figure output | `cluster/scripts/utils/multi_format_output.py` (patches `Figure.savefig`; covers both `plt.savefig` and `fig.savefig` — proven in Phase 4) |
| deepTools bed_pileup | `cluster/deepTools_pipeline.py:25` (`bed_pileup`) — Phase 0 fixed `genome='mm10'`; pass `bigWig_dict` directly (`bam_dict=None`); shells out to `computeMatrix reference-point` |
| deepTools heatmap render | `cluster/deeptools_plotting.py:15` (`heatmap_plot`) — supports `color_dict`, `vmax_groups`, `up_down`; reads `<out_name>_values` text matrix |
| cooltools pileup | `cluster/cooltools_called.py:28` (`mcool_pileup`) — Phase 0 fixed `genome='mm10'`; supports `bedpe_dict` (cooler col convention) |
| sbatch template | `stripes/stripenn/scripts/01_call_stripes.sb` (gold standard: `csd940`, `shared`, `source ~/.bashrc; conda activate`, two-root paths, manual `STATUS=$?` checks, no `set -euo pipefail`) |
| Driver pattern | `cluster/scripts/run_phase4.sh` (positional `ANALYSES`, `LOG` env var, banner + `tee`, post-run `ls -lh` audit) |

Nothing new to invent — the work is wiring proven helpers to the right inputs.

---

## Phase 5: deepTools Metagene (LOCAL)

### 5.1 — Per-cluster anchor BED preparation

For each of the 6 clusters in `combined-clusters.txt`, concatenate both anchor halves and deduplicate to a 3-column BED:

```python
sub = bedpe_df[bedpe_df['GROUP'] == cluster]
a1 = sub[['chr1', 'x1', 'x2']].rename(columns={'chr1':'chrom','x1':'start','x2':'end'})
a2 = sub[['chr2', 'y1', 'y2']].rename(columns={'chr2':'chrom','y1':'start','y2':'end'})
anchors = pd.concat([a1, a2]).drop_duplicates().sort_values(['chrom','start','end'])
anchors.to_csv(out_path, sep='\t', header=False, index=False)
```

Outputs (one file per cluster): `cluster/bap1_late/figures/deeptools_input/{clust1..clust6}_anchors.bed`

Expected unique-anchor counts (loops × 2 minus hub-anchor dedup, ~88–92% unique):
| Cluster | Loops | Expected unique anchors |
|---------|------:|------------------------:|
| clust1 | 12,298 | ~22,000 |
| clust2 | 10,970 | ~19,500 |
| clust3 |  8,738 | ~15,500 |
| clust4 |  3,916 |  ~7,000 |
| clust5 |    667 |  ~1,200 |
| clust6 |  2,359 |  ~4,200 |

### 5.2 — Run deepTools bed_pileup

Single call to `bed_pileup` produces one combined heatmap (rows = 6 clusters, cols = 8 BigWigs):

```python
BIGWIG_BASE = Path('/Users/zakiralibhai/sdsc/bigwigs')   # also session working dir
bigwig_dict = {
    'H3K27ac_ctrl':   BIGWIG_BASE / 'H3K27acCtrl.bw',
    'H3K27ac_mut':    BIGWIG_BASE / 'H3K27acMut.bw',
    'H3K27me3_ctrl':  BIGWIG_BASE / 'H3K27me3Ctrl.bw',
    'H3K27me3_mut':   BIGWIG_BASE / 'H3K27me3Mut.bw',
    'H2AK119ub_ctrl': BIGWIG_BASE / 'H2AK119ubCtrl.bw',
    'H2AK119ub_mut':  BIGWIG_BASE / 'H2AK119ubMut.bw',
    'H3K27me1_ctrl':  BIGWIG_BASE / 'H3K27me1Ctrl.bw',
    'H3K27me1_mut':   BIGWIG_BASE / 'H3K27me1Mut.bw',
}
vmax_groups = [
    ['H3K27ac_ctrl',   'H3K27ac_mut'],
    ['H3K27me3_ctrl',  'H3K27me3_mut'],
    ['H2AK119ub_ctrl', 'H2AK119ub_mut'],
    ['H3K27me1_ctrl',  'H3K27me1_mut'],
]
color_dict = {
    'H3K27ac_ctrl':'Blues','H3K27ac_mut':'Blues',
    'H3K27me3_ctrl':'Reds','H3K27me3_mut':'Reds',
    'H2AK119ub_ctrl':'Greens','H2AK119ub_mut':'Greens',
    'H3K27me1_ctrl':'Purples','H3K27me1_mut':'Purples',
}
sub = figure_subfolder(REPO_ROOT/'cluster/bap1_late/figures/deeptools', 'histone_anchors')
with multi_format_savefig():
    bed_pileup(
        bed_dict=bed_dict,            # 6-key dict from 5.1
        bigWig_dict=bigwig_dict,
        out_dir=str(sub),
        out_name='histone_anchors',
        blacklisted_regions=str(REPO_ROOT/'tads/mm10-blacklist.v2.bed'),
        up_down=5000,                  # ±5kb (Popay default)
        color_dict=color_dict,
        vmax_groups=vmax_groups,       # pair ctrl/mut → fair visual comparison
    )
```

`bed_pileup` shells out to `computeMatrix reference-point --referencePoint center --missingDataAsZero -bl <blacklist> -p 4 ...` (verbatim flags from `deepTools_pipeline.py:56`), then calls `heatmap_plot()` which renders a custom seaborn/matplotlib figure (rows = clusters, cols = bigwigs, top row = mean profile lines, lower rows = sorted heatmaps, last row = colorbars). With `multi_format_savefig()` wrapping, the .png save fires the patch and emits `.svg`/`.pdf`/`.jpg` siblings.

**Output files** in `cluster/bap1_late/figures/deeptools/histone_anchors/`:
- `histone_anchors` (binary computeMatrix matrix, gzipped)
- `histone_anchors_values` (tab-text values for `heatmap_plot` to read back)
- `histone_anchors_values.{png,pdf,svg,jpg}` (the heatmap)

### 5.3 — Driver script

`cluster/scripts/run_phase5.sh` mirrors `run_phase4.sh` (`set -e`, banner, `tee` to `$LOG`, post-run `ls -lh`) but with one critical difference: **must put cluster env's `bin/` on `$PATH`** so the `computeMatrix` subprocess call inside `bed_pileup` resolves. Verified locally: `/opt/homebrew/anaconda3/envs/cluster/bin/computeMatrix` exists (deepTools 3.5.5) but is not on the default shell PATH.

```bash
#!/usr/bin/env bash
# cluster/scripts/run_phase5.sh
set -e
cd "$(dirname "$0")/../.."   # repo root

CLUSTER_ENV_BIN=/opt/homebrew/anaconda3/envs/cluster/bin
export PATH="${CLUSTER_ENV_BIN}:${PATH}"   # exposes computeMatrix to bed_pileup subprocess
PYTHON="${CLUSTER_ENV_BIN}/python3"
SCRIPT=cluster/scripts/06_deeptools_metagene.py
LOG=${LOG:-cluster/phase5.txt}

{
  echo "============================================================"
  echo "Phase 5: deepTools metagene at loop anchors"
  echo "Repo root: $(pwd)"
  echo "Started:   $(date)"
  echo "log_path:  ${LOG}"
  echo "python:    ${PYTHON}"
  echo "computeMatrix: $(which computeMatrix 2>/dev/null || echo NOT_ON_PATH)"
  echo "============================================================"
  echo
  echo "[1/1] Running 06_deeptools_metagene.py..."
  "$PYTHON" "$SCRIPT"
  echo
  echo "============================================================"
  echo "Phase 5 outputs"
  echo "============================================================"
  echo "--- per-cluster anchor BEDs ---"
  ls -lh cluster/bap1_late/figures/deeptools_input/clust*_anchors.bed 2>/dev/null
  echo "--- metagene heatmap ---"
  ls -lh cluster/bap1_late/figures/deeptools/histone_anchors/ 2>/dev/null
  echo
  echo "Phase 5 finished: $(date)"
} 2>&1 | tee "$LOG"
```

### 5.4 — Phase 5 verification

| Check | Expected |
|-------|---------|
| Anchor BED count per cluster | 6 files, each ≥600 lines (clust5 floor) |
| Anchor BED dedup correctness | `sort -u` of full file = full file |
| Heatmap values text dimensions | rows = sum(unique anchors per cluster) ≈ 70k; cols = 8 BigWigs × (5000+5000)/100bp + 1 group label = ~801 |
| Heatmap PNG file size | 100KB–500KB (vector points × ~70k regions × 8 cols) |
| 4 file formats produced | `histone_anchors_values.{png,pdf,svg,jpg}` all present and non-empty |
| Visual signature (clust5) | K27me3 + K119ub + K27me1 columns show pronounced central peak in BOTH ctrl and mut (Polycomb-rich anchors); K27ac near zero |
| Visual signature (clust1) | Mixed/active marks (K27ac high, K27me3 low) — bulk unchanged loops |
| ctrl vs mut paired panels | Height differences within each `vmax_groups` pair indicate true mut effect (not vmax-scaling artifact) |

---

## Phase 6: Cooltools Pileup (HPC sbatch)

### 6.1 — One-time HPC env setup

User-provided manual step (NOT automated by the pipeline):

```bash
# On Expanse login node (one-time)
ssh <expanse-alias>
cd /expanse/lustre/projects/csd940/zalibhai/mariner_hi-c
git pull   # pull the new Phase 5/6 scripts
conda env create -f cluster/environment_linux.yml -n cluster
conda activate cluster
python -c "import cooltools, bioframe, cooler; print('cooltools', cooltools.__version__)"
```

Document this in the `.sb` header. If `cluster/environment_linux.yml` is missing cooltools/bioframe/pybbi (verify before submission), follow up with:

```bash
conda activate cluster
pip install cooltools bioframe pybbi
```

### 6.2 — Python orchestrator

`cluster/scripts/07_cooltools_pileup.py`:

```python
#!/usr/bin/env python3
# cluster/scripts/07_cooltools_pileup.py
# Phase 6: cooltools off-diagonal pileup of cluster BEDPEs against ctrl/mut mcools.

import argparse, sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')           # disable GUI; headless on HPC compute nodes
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
CLUSTER_DIR = SCRIPT_DIR.parent
REPO_ROOT = CLUSTER_DIR.parent
sys.path.insert(0, str(CLUSTER_DIR))            # cooltools_called
sys.path.insert(0, str(SCRIPT_DIR / 'utils'))   # multi_format_output

from cooltools_called import mcool_pileup
from multi_format_output import multi_format_savefig, figure_subfolder

CLUSTER_FILE = REPO_ROOT / 'cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt'
OUT_BASE     = REPO_ROOT / 'cluster/bap1_late/cooltools'
CLUSTER_ORDER = ['clust1','clust2','clust3','clust4','clust5','clust6']

def load_bedpe_dict(path):
    df = pd.read_csv(path, sep='\t')
    bedpe_dict = {}
    for c in CLUSTER_ORDER:
        sub = df[df['GROUP'] == c].copy()
        # cooler / cooltools convention: chrom1/start1/end1/chrom2/start2/end2
        sub = sub.rename(columns={
            'chr1':'chrom1','x1':'start1','x2':'end1',
            'chr2':'chrom2','y1':'start2','y2':'end2',
        })
        bedpe_dict[c] = sub[['chrom1','start1','end1','chrom2','start2','end2']].reset_index(drop=True)
    return bedpe_dict

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mcool-ctrl', required=True)
    ap.add_argument('--mcool-mut',  required=True)
    ap.add_argument('--resolution', type=int, default=10000, help='10kb pileup binsize')
    ap.add_argument('--flank',      type=int, default=500000, help='±500kb window per Popay')
    ap.add_argument('--out-dir',    default=str(OUT_BASE))
    args = ap.parse_args()

    print(f'[1/2] Loading clusters from {CLUSTER_FILE}...', flush=True)
    bedpe_dict = load_bedpe_dict(CLUSTER_FILE)
    for c, df in bedpe_dict.items():
        print(f'  {c}: {len(df)} loops', flush=True)

    mcool_dict = {'ctrl': args.mcool_ctrl, 'mut': args.mcool_mut}
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sub = figure_subfolder(out_dir, 'obs_exp_contacts')

    print(f'\n[2/2] cooltools pileup at {args.resolution}bp / flank ±{args.flank}bp...', flush=True)
    print('  NOTE: bioframe.fetch_chromsizes(\"mm10\") + fetch_centromeres requires '
          'outbound HTTPS to UCSC. If this fails on a compute node, pre-fetch '
          '`arms` on the login node and patch make_arms to load from cache.', flush=True)
    with multi_format_savefig():
        mcool_pileup(
            mcool_dict=mcool_dict,
            bedpe_dict=bedpe_dict,
            out_dir=str(sub),
            out_name='obs_exp_contacts',
            flank=args.flank,
            resolution=args.resolution,
            over_expected=True,
            split_diagonal=False,    # latent bug if True (Plan-p1 §3 corrections)
            v_range=[-1, 2],         # Popay notebook value
            genome='mm10',           # Phase 0 fix; redundant with default
        )
    print('\nPhase 6 complete. Outputs:', flush=True)
    print(f'  {sub}/obs_exp_contacts_ctrl.{{png,pdf,svg,jpg}}')
    print(f'  {sub}/obs_exp_contacts_mut.{{png,pdf,svg,jpg}}')

if __name__ == '__main__':
    main()
```

`mcool_pileup` writes via `plt.savefig(..., dpi=300)` (lines 158–159 of `cooltools_called.py`), which routes through `Figure.savefig` → `multi_format_savefig()` patch → emits 4 sibling formats per condition.

### 6.3 — sbatch wrapper

`cluster/scripts/07_cooltools_pileup.sb` (mirrors `stripes/stripenn/scripts/01_call_stripes.sb` template):

```bash
#!/bin/bash
# cluster/scripts/07_cooltools_pileup.sb
# Phase 6: cooltools off-diagonal pileup over loop-cluster BEDPEs (ctrl + mut, late tp).
# One-time setup on Expanse login node before first submission:
#   conda env create -f cluster/environment_linux.yml -n cluster
#
# Usage (from CODE_DIR on Expanse login node):
#   sbatch cluster/scripts/07_cooltools_pileup.sb
#
# Optional override (default uses 250402 merged mcools):
#   sbatch cluster/scripts/07_cooltools_pileup.sb <ctrl_mcool> <mut_mcool>
#SBATCH --job-name=cooltools_pileup
#SBATCH --output=cluster/logs/phase6_pileup_%j.out
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --account=csd940

CODE_DIR="/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c"
DATA_DIR="/expanse/lustre/projects/csd940/zalibhai/stripes/stripenn"

MCOOL_CTRL="${1:-${DATA_DIR}/data/cool/250402/ctrl_merged.mcool}"
MCOOL_MUT="${2:-${DATA_DIR}/data/cool/250402/mut_merged.mcool}"

mkdir -p "${CODE_DIR}/cluster/logs"
mkdir -p "${CODE_DIR}/cluster/bap1_late/cooltools"

echo "==========================================="
echo "Phase 6: cooltools off-diagonal pileup"
echo "==========================================="
echo "Job ID:      ${SLURM_JOB_ID}"
echo "Node:        ${HOSTNAME}"
echo "CPUs:        ${SLURM_CPUS_PER_TASK}"
echo "Memory:      ${SLURM_MEM_PER_NODE} MB"
echo "CODE_DIR:    ${CODE_DIR}"
echo "MCOOL_CTRL:  ${MCOOL_CTRL}"
echo "MCOOL_MUT:   ${MCOOL_MUT}"
echo "Start:       $(date)"
echo "==========================================="

if [ ! -s "${MCOOL_CTRL}" ]; then
  echo "ERROR: missing or empty ctrl mcool: ${MCOOL_CTRL}" >&2
  exit 2
fi
if [ ! -s "${MCOOL_MUT}" ]; then
  echo "ERROR: missing or empty mut mcool: ${MCOOL_MUT}" >&2
  exit 2
fi

source ~/.bashrc
conda activate cluster

which python
python -c "import cooltools, bioframe, cooler; print('cooltools', cooltools.__version__, 'bioframe', bioframe.__version__, 'cooler', cooler.__version__)"
echo "env OK"

export PYTHONWARNINGS="ignore::FutureWarning"

cd "${CODE_DIR}"

python cluster/scripts/07_cooltools_pileup.py \
  --mcool-ctrl "${MCOOL_CTRL}" \
  --mcool-mut  "${MCOOL_MUT}" \
  --resolution 10000 \
  --flank      500000
STATUS=$?

OUT_DIR="${CODE_DIR}/cluster/bap1_late/cooltools/obs_exp_contacts"

if [ ${STATUS} -ne 0 ]; then
  echo "ERROR: cooltools_pileup failed (exit=${STATUS})." >&2
  exit 1
fi

if [ ! -s "${OUT_DIR}/obs_exp_contacts_ctrl.png" ]; then
  echo "ERROR: ctrl heatmap missing or empty after run." >&2
  exit 1
fi
if [ ! -s "${OUT_DIR}/obs_exp_contacts_mut.png" ]; then
  echo "ERROR: mut heatmap missing or empty after run." >&2
  exit 1
fi

echo ""
echo "Output files:"
ls -lh "${OUT_DIR}/" 2>/dev/null

echo ""
echo "End: $(date)"
echo "==========================================="
```

Key choices vs the stripenn template:
- Partition / account / no-`set -euo pipefail` / `source ~/.bashrc; conda activate <env>` / two-root paths / manual `STATUS=$?` checks / `if [ ! -s ... ]` validation — all preserved.
- 16 CPUs (cooltools `expected_cis` uses `nproc=2` internally; remaining cores for I/O & balance reads).
- 64G RAM is generous: cooler `expected_cis` with chunksize=1M and ~38k loops × 2 conditions × 6 clusters at 10kb / ±500kb is ~16-32GB peak.
- 4h time budget: ~30–60min for `expected_cis` per mcool + ~10–30min for `pileup` = ~2h total; 4h is safe.
- Log to `cluster/logs/phase6_pileup_${SLURM_JOB_ID}.out` so subsequent re-runs don't clobber.

### 6.4 — Submission

From the Expanse login node, after one-time env setup:

```bash
cd /expanse/lustre/projects/csd940/zalibhai/mariner_hi-c
sbatch cluster/scripts/07_cooltools_pileup.sb
squeue -u $USER
# When complete:
ls -lh cluster/bap1_late/cooltools/obs_exp_contacts/
```

To pull outputs back to local Mac:
```bash
# from local mariner_hi-c
git pull   # if outputs were committed on HPC, OR
rsync -avP <expanse>:/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/cluster/bap1_late/cooltools/ cluster/bap1_late/cooltools/
```

### 6.5 — Phase 6 verification

| Check | Expected |
|-------|---------|
| sbatch submission | Returns Job ID; `squeue` shows R/PD until complete |
| SLURM exit code | 0 (no `STATUS != 0` errors in `phase6_pileup_*.out`) |
| Output files | 8 files total: `obs_exp_contacts_{ctrl,mut}.{png,pdf,svg,jpg}` in `cluster/bap1_late/cooltools/obs_exp_contacts/` |
| File sizes | PNG/JPG 50–250KB; PDF 30–100KB; SVG 50–500KB |
| ctrl vs mut visual difference (clust5) | Mut heatmap should show STRONGER central red dot (loop strengthening) — gain cluster |
| ctrl vs mut visual difference (clust6) | Mut heatmap should show WEAKER central red dot (loop weakening) — loss cluster |
| ctrl vs mut visual difference (clust1) | Heatmaps near-identical (unchanged bulk) |
| Background quadrants | Equal/random pattern outside the central pixel — expected log2(obs/exp) scatter ≈ 0 |
| Internet connectivity | If `bioframe.fetch_chromsizes` fails, error contains `urllib`/`HTTPError`; fallback noted in §6.6 |

### 6.6 — Internet caveat & fallback

`cooltools_called.make_arms()` at `cooltools_called.py:185–190` calls `bioframe.fetch_chromsizes('mm10')` and `bioframe.fetch_centromeres('mm10')`, both of which hit `https://hgdownload.soe.ucsc.edu/`. SDSC Expanse compute nodes typically allow outbound HTTPS, but if firewalled the script fails immediately. **Fallback** (only apply if HTTP fails):

```bash
# On Expanse login node (where internet works):
conda activate cluster
python -c "
import bioframe, pickle
chromsizes = bioframe.fetch_chromsizes('mm10')
cens = bioframe.fetch_centromeres('mm10')
arms = bioframe.make_chromarms(chromsizes, cens)
arms.to_csv('/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/cluster/data/mm10_arms.tsv', sep='\t', index=False)
print('cached', len(arms), 'arms')
"
```

Then patch `07_cooltools_pileup.py` to read from cache via a `make_arms_cached(cache_path)` shim that returns the saved DataFrame. Documenting only — do not apply unless the first sbatch run errors with an HTTP exception.

---

## Files to Create

| Path | Phase | Purpose |
|------|------:|---------|
| `cluster/scripts/06_deeptools_metagene.py` | 5 | Anchor BED prep + bed_pileup orchestrator (single combined heatmap, 6 clusters × 8 BigWigs) |
| `cluster/scripts/run_phase5.sh` | 5 | Local driver: prepends cluster env bin to PATH, tees log to `cluster/phase5.txt` |
| `cluster/scripts/07_cooltools_pileup.py` | 6 | Loads combined-clusters.txt → cooler-convention bedpe_dict; calls `mcool_pileup` |
| `cluster/scripts/07_cooltools_pileup.sb` | 6 | SLURM wrapper: csd940 / shared / 16 cpus / 64G / 4h, two-root paths, manual exit checks |

## Files to Modify

None. Phase 0 corrections to `deepTools_pipeline.py`, `deeptools_plotting.py`, and `cooltools_called.py` are already in place. `multi_format_output.py` (Phase 4 patch) handles `Figure.savefig` correctly for both deepTools and cooltools.

## Outputs Created (under `cluster/bap1_late/`)

| Path | Phase | Format |
|------|------:|--------|
| `figures/deeptools_input/clust{1..6}_anchors.bed` | 5 | 3-col BEDs, deduplicated |
| `figures/deeptools/histone_anchors/histone_anchors` | 5 | Gzip computeMatrix matrix |
| `figures/deeptools/histone_anchors/histone_anchors_values` | 5 | Tab-text values for heatmap_plot |
| `figures/deeptools/histone_anchors/histone_anchors_values.{png,pdf,svg,jpg}` | 5 | Combined heatmap |
| `cooltools/obs_exp_contacts/obs_exp_contacts_ctrl.{png,pdf,svg,jpg}` | 6 | Per-cluster pileup, ctrl |
| `cooltools/obs_exp_contacts/obs_exp_contacts_mut.{png,pdf,svg,jpg}` | 6 | Per-cluster pileup, mut |
| `cluster/phase5.txt` | 5 | Local run log (tee'd) |
| `cluster/logs/phase6_pileup_<jobid>.out` | 6 | SLURM stdout (on HPC, synced back) |

## Critical-File Reference

| File | Why critical |
|------|--------------|
| `cluster/bap1_late/cluster3/k-6/data/combined-clusters.txt` | Phase 3 v2 canonical clustering (38,948 loops × 6 clusters); both Phase 5 & 6 read this |
| `/Users/zakiralibhai/sdsc/bigwigs/{H3K27ac,H3K27me3,H2AK119ub,H3K27me1}{Ctrl,Mut}.bw` | Phase 5 BigWigs (8 files, all verified present, 19MB–127MB each) |
| `tads/mm10-blacklist.v2.bed` | Phase 5 deeptools `-bl` blacklist (verified: 3,435 lines) |
| `/expanse/.../stripes/stripenn/data/cool/250402/{ctrl,mut}_merged.mcool` | Phase 6 inputs; KR weights via cooler balance per stripenn 5.2 history |
| `cluster/deepTools_pipeline.py:25` | `bed_pileup` — Phase 5 entrypoint |
| `cluster/cooltools_called.py:28` | `mcool_pileup` — Phase 6 entrypoint |
| `cluster/scripts/utils/multi_format_output.py` | Patches `Figure.savefig` so `.png` writes emit `.svg/.pdf/.jpg` siblings (works for both deepTools `fig.savefig` and cooltools `plt.savefig`) |
| `stripes/stripenn/scripts/01_call_stripes.sb` | sbatch gold-standard template (account, partition, conda activation, two-root paths) |

## Verification Plan (end-to-end)

1. **Phase 5 anchor BEDs** — `wc -l cluster/bap1_late/figures/deeptools_input/clust*_anchors.bed` shows 6 files, each ≥600 lines, line counts roughly proportional to loop counts.
2. **Phase 5 metagene** — `bash cluster/scripts/run_phase5.sh` completes without error in ~5–15 min; `histone_anchors_values.png` exists at >100KB; `cluster/phase5.txt` is non-empty and ends with "Phase 5 finished".
3. **Phase 5 visual sanity** — open `histone_anchors_values.pdf`. clust5 row should show K27me3+K119ub+K27me1 columns with strong central peaks; K27ac should be flat in clust5. clust1 should be K27ac-dominant. Mut vs ctrl panels (paired vmax) should reveal mark-specific changes.
4. **Phase 6 env on HPC** — `conda env create -f cluster/environment_linux.yml -n cluster` succeeds; `python -c "import cooltools, bioframe, cooler; print(cooltools.__version__)"` returns a version.
5. **Phase 6 sbatch** — `sbatch cluster/scripts/07_cooltools_pileup.sb` returns a JobID; `squeue -u $USER` shows the job; on completion `sacct` shows `COMPLETED` state.
6. **Phase 6 outputs** — 8 files (ctrl+mut × 4 formats) exist in `obs_exp_contacts/`; SLURM log ends with `End: <date>` and no `ERROR` lines.
7. **Phase 6 visual sanity** — clust5 mut heatmap shows stronger central red pixel than clust5 ctrl (gain); clust6 mut weaker than clust6 ctrl (loss); clust1 nearly identical between conditions (unchanged bulk).
8. **No regressions** — Phase 4 outputs untouched: `cluster/bap1_late/chromHMM/{anchor,span}.txt` and `figures/...` from Phase 4 unchanged on disk.

## Post-Completion Task (NOT in this plan's scope)

After Phase 5 + 6 are run successfully, append two new sections to `cluster/plan-p2.md` mirroring the format used for Phases 1–4 in `plan-p1.md` / `plan-p2.md`:

- **Phase 5: deepTools Metagene Analysis — DONE (date)** with Status, Corrections applied during execution (real per-cluster anchor counts, any deepTools quirks discovered, vmax pairing efficacy, etc.), and a result-summary table.
- **Phase 6: Cooltools Pileup — DONE (date)** with Status, HPC-specific corrections (e.g., env setup time, internet connectivity outcome), and visual-result summary (clust5 strengthening / clust6 weakening / clust1 invariance — quantified by central-pixel obs/exp where possible).

Never delete content from `plan-p2.md`; append-only with explicit completion notes per the user's "never deleting things but properly noting" instruction. This update is post-implementation and can be made via Edit tool once plan mode is exited and Phases 5/6 have completed.
