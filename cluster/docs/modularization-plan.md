# Plan: Modularize the Hi-C Loop Clustering Pipeline into a Reusable Tool

## Context

The Popay-style Hi-C loop clustering pipeline (`cluster/`) was built as a single-project analysis for BAP1-KO mouse cerebellum. It works end-to-end for that dataset but is not reusable by another researcher on a different organism, condition, mark set, or machine. The user wants an analysis of the current state and a plan to make it a proper generalizable tool.

This plan covers: (1) a diagnostic of what's hardcoded and where, (2) how Popay's original code is used, and (3) the minimum-viable refactoring to make it reusable — without a full rewrite.

---

## Part 1: Current State Analysis

### 1A. Hardcoding Severity by Category

**Category 1 — Machine-specific absolute paths (blocks ANY other user)**

| Path | Where it appears | Count |
|------|-----------------|-------|
| `/Users/zakiralibhai/sdsc/bigwigs` | `05`, `06`, `08`, `09`, `10`, `12` + both `.conf` files | 8 files |
| `/opt/homebrew/anaconda3/envs/cluster/bin/python3` | `run_phase{1,3,4,5}`, `run_summary`, `run_oriented_metagene`, `run_clust6_subgroups` + error messages in `06`, `09`, `10` | 10 files |
| `/expanse/lustre/projects/csd940/zalibhai/...` | `03b_chromhmm_9mark.sb`, `07_cooltools_pileup.sb`, `12_comprehensive_asymmetry.sb` + config | 5 files |
| `/usr/local/bin/cluster` (Cluster 3.0 binary) | `04_clustering.py` + both config files | 3 files |
| `/usr/local/bin/Rscript` | `run_phase1.sh`, `run_go_kegg.sh` + config | 4 files |

**Category 2 — Project-specific constants (fails silently on new data)**

| Constant | Where | Impact |
|----------|-------|--------|
| `EXPECTED_TOTAL = 39344`, `EXPECTED_PER_RES = {5: 7901, ...}` | `01_build_loop_count_file.py` | Hard assertion blocks any dataset with different loop counts |
| `CTRL_REPS = ['ctrl_M1','ctrl_M2','ctrl_M3']`, `MUT_REPS` | `01_build_loop_count_file.py` | Breaks if different replicate names or number |
| `DATA_COLS = ['ctrl_merge','mut_merge']`, `NORMALIZE_COL = 'ctrl_merge'` | `04_clustering.py` | Breaks if different condition column names |
| `CLUSTER_ORDER = ['clust1',...,'clust6']` | 6+ scripts | Breaks if k != 6 |
| `RESOLUTIONS = [5, 10, 25]` | `01_build_loop_count_file.py` | Breaks if data uses other resolutions |
| `SIZE_THRESHOLD = 800_000` | `10_clust6_subgroup_asymmetry.py`, `12_comprehensive_asymmetry.py` | Duplicated, project-specific split point |

**Category 3 — Biological labels baked into figure code**

| What | Where | Problem |
|------|-------|---------|
| `BIO_ORDER = ['clust6','clust3','clust1','clust2','clust4','clust5']` | `08`, `11`, `13`, `visualize_orientation_asymmetry.py` + 2 others | Biological reordering specific to BAP1-KO results |
| `BIO_NAMES` with `'Strong gain\n(Polycomb domain)'` etc. | `08_summary_figures.py` | Interpretive biology in code |
| `BIO_LABELS` with `'loss (78%)'`, `'gain (97%)'` | `11_histone_anchors_metagene.py`, `visualize_orientation_asymmetry.py` | Specific percentages from one dataset |
| `'11_Polycomb'`, `'12_Bivalent_Enhancer'` row selectors | `08_summary_figures.py` L304-305, L499 | Assumes ChromHMM state E11=Polycomb |
| Panel titles: `'Gained loops (clust5, n=667)'` | `08_summary_figures.py`, `13_mechanism_9mark.py` | Hardcoded n-values |

**Category 4 — Code duplication across scripts**

| Duplicated block | Scripts | Lines each |
|-----------------|---------|-----------|
| `BIGWIG_DICT` (8-entry filename dict) | `05`, `06`, `09`, `10` + partial in `08`, `12` | ~10 lines x 5 |
| `parse_header()` (deepTools matrix parser) | `quantify_orientation_asymmetry`, `10`, `11`, `12`, `visualize_orientation_asymmetry` | ~25 lines x 5 |
| `BIO_ORDER` / `BIO_NAMES` / `BIO_LABELS` | `08`, `11`, `13`, `visualize_...`, `quantify_...` | ~15 lines x 6 |
| `VMAX_GROUPS` + `COLOR_DICT` | `06`, `09`, `10` | ~12 lines x 3 |
| `STATE_COLORS` palette dict | `05`, `08`, `13` | ~25 lines x 3 |

---

### 1B. How Popay's Original Code Is Used

**Architecture:** 7 Python modules live in `cluster/modules/`. They are NOT a Python package (no `__init__.py`). Scripts access them via `sys.path.insert(0, str(CLUSTER_DIR / 'modules'))` at the top of each file.

**Module inventory and modification status:**

| Module | Functions | Modified from Popay? | Used by scripts |
|--------|-----------|---------------------|-----------------|
| `plotting.py` | `heat()`, `line()`, `box()`, `strip()`, `stacked()`, `bar()`, `joint()` | Phase 0 only: fixed custom_params path, removed gene-specific ylim overrides, fixed logX bug | `04`, `05` |
| `cluster_tools.py` | `elbow()`, `sort_clusters()`, `comparison_type()`, `sort_by_strings()` | Phase 0 only: fixed savefig/show order | `04`, `05` |
| `statistics_functions.py` | `kruskal_wilcoxon()` | Untouched | `05` |
| `chromHMM_heatmap.py` | `heatmap_plot()` | Phase 0 only: fixed custom_params path | `05` |
| `deeptools_plotting.py` | `heatmap_plot()` | Phase 0 only: fixed custom_params path, removed HA-NIPBL override | `06`, `09`, `10` (via `bed_pileup`) |
| `deepTools_pipeline.py` | `bam_coverage()`, `bed_pileup()` | Phase 0: added genome param with mm10 size | `06`, `09`, `10` |
| `cooltools_called.py` | `mcool_pileup()` | Added `_viewframe_from_cooler()` (new), changed default genome to mm10 | `07` |
| `bedpe_analysis.py` | `loop_size()`, `bedtools_annotation()` | Phase 0: fixed gene path default, added FPKM None guard | `05` |

**Two-tier coupling pattern:**

- **Tier 1 — Popay-coupled (core pipeline):** Scripts `04`, `05`, `06`, `07`, `09` directly call Popay module functions. These scripts adapt Popay's Jupyter notebook cells 1-25.
- **Tier 2 — Self-contained (later extensions):** Scripts `01`, `08`, `10`, `11`, `12`, `13`, `quantify_*`, `visualize_*` have zero Popay module imports. They only use `utils/multi_format_output.py` (original code). The development pattern moved away from Popay modules as the analysis matured past the core clustering phase.

**Key takeaway:** Popay's modules are a thin library (~800 lines total) used primarily in Phases 3-6. The later half of the pipeline (Phases 7-11) is entirely self-contained. The modules were stabilized in Phase 0 with minimal changes (path fixes, deprecated API updates, hg38->mm10 defaults) and have been untouched since.

---

### 1C. What's Already Parameterized (partial credit)

The pipeline is not starting from zero. A config system already exists:

- **`scripts/config/late.conf` and `early.conf`** — shell-sourceable key=value files covering timepoint identity, clustering params, input path templates, peak BED paths, and some machine-specific paths. Already used by `run_full_pipeline.sb`.
- **`05_grouped_analyses.py`** has the richest env-var override surface: 10+ variables readable via `os.environ.get(VAR, default)`. The `CLUSTER_NSTATES`/`CLUSTER_CHROMHMM_SUBDIR`/`CLUSTER_ENRICH_SUFFIX` system works — proven by `run_phase4_9mark.sh` which sets env vars to rerun Phase 4.4/4.5 for the 9-mark model.
- **`01_build_loop_count_file.py`** reads `CLUSTER_TIMEPOINT_LABEL`, `CLUSTER_TIMEPOINT_ID`, and path templates from env vars.
- **`04_clustering.py`** reads `CLUSTER_OUT_DIR`, `CLUSTER_COUNT_FILE`, `CLUSTER_BIN` from env vars.

What's NOT parameterized (the gaps):
- Replicate names/count — completely absent from any env-var surface
- Condition column names — hardcoded, no override
- BigWig filenames — only the directory is parameterized; individual basenames are hardcoded
- Binary interpreter paths in runner scripts — hardcoded, not env-var gated
- Resolution list — hardcoded `[5, 10, 25]`
- Biological ordering/labels — hardcoded in 6 scripts

---

## Part 2: Modularization Plan

### Design Principles

1. **Surgical, not rewrite.** Extend the existing env-var + `.conf` pattern (proven working) rather than introducing new config formats.
2. **Backward compatible.** All new env vars have defaults matching current BAP1-KO values. Existing analysis produces identical output without setting anything new.
3. **Don't touch Popay modules.** They're stabilized and working. Wrap or work around them.
4. **One new shared module** (`scripts/utils/pipeline_config.py`) to eliminate all duplication.

### Step 1: Create `scripts/utils/pipeline_config.py` (~120 lines)

This is the cornerstone. It provides:

```python
# Deduplicated constants
COORD_COLS = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2']

# Functions replacing copy-pasted code
def parse_header(filepath)          # deepTools matrix parser (removes from 5 scripts)
def get_cluster_order(k)            # ['clust1',...,'clustK'] from CLUSTER_K
def get_bio_order(cluster_order)    # from CLUSTER_BIO_ORDER env or fallback to cluster_order
def get_bio_names(bio_order)        # from CLUSTER_BIO_NAME_{c} env vars or fallback to cluster ID
def get_condition_cols()            # from CLUSTER_COND1_COL / COND2_COL env vars
def get_size_threshold()            # from CLUSTER_SIZE_THRESHOLD env var, default 800_000
def build_bigwig_dict(bigwig_dir)   # from CLUSTER_BW_* env vars, with current filenames as defaults
def build_vmax_groups(bigwig_dict)  # derive ctrl/mut pairs by mark prefix
def build_color_dict(bigwig_dict)   # derive colormaps by mark name
def get_key_state_id(default)       # from CLUSTER_KEY_STATE_ID env var
```

All defaults match current BAP1-KO behavior. A new researcher only sets what differs.

### Step 2: Upgrade `.conf` files (both `late.conf` and `early.conf`)

Add three new blocks (new variables commented-out by default for backward compat):

**Block A — Experimental design:**
```bash
CLUSTER_CTRL_REPS="ctrl_M1,ctrl_M2,ctrl_M3"
CLUSTER_MUT_REPS="mut_M1,mut_M2,mut_M3"
CLUSTER_COND1_COL="ctrl_merge"
CLUSTER_COND2_COL="mut_merge"
CLUSTER_RESOLUTIONS="5,10,25"
# CLUSTER_EXPECTED_TOTAL=39344  # optional validation; omit to skip
```

**Block B — Machine-specific binaries:**
```bash
# CLUSTER_PYTHON="/opt/homebrew/anaconda3/envs/cluster/bin/python3"
# CLUSTER_RSCRIPT="/usr/local/bin/Rscript"
# CLUSTER_BIN="/usr/local/bin/cluster"
```

**Block C — Biological annotation (fill in AFTER inspecting clustering results):**
```bash
# CLUSTER_BIO_ORDER="clust6,clust3,clust1,clust2,clust4,clust5"
# CLUSTER_KEY_STATE_ID="11_Polycomb"
# CLUSTER_BIO_NAME_clust5="Strong gain\n(Polycomb domain)"
# CLUSTER_BIO_NAME_clust6="Strong loss\n(anchor-disrupted)"
```

### Step 3: Parameterize `01_build_loop_count_file.py` (highest value)

4 changes:
- `CTRL_REPS` / `MUT_REPS` read from `CLUSTER_CTRL_REPS` / `CLUSTER_MUT_REPS` env vars (comma-split)
- Resolution list read from `CLUSTER_RESOLUTIONS` env var
- `EXPECTED_TOTAL` check becomes optional: if `CLUSTER_EXPECTED_TOTAL` not set, warn and skip (don't fail)
- Column names `'ctrl_merge'` / `'mut_merge'` read from `CLUSTER_COND1_COL` / `CLUSTER_COND2_COL`

### Step 4: Parameterize `04_clustering.py` (small)

2 changes:
- `DATA_COLS` reads from `get_condition_cols()` instead of hardcoded list
- `NORMALIZE_COL` becomes `COND1_COL` from the same source

### Step 5: Replace `BIGWIG_DICT` copies in 5 scripts

In `05_grouped_analyses.py`, `06_deeptools_metagene.py`, `09_oriented_anchor_metagene.py`, `10_clust6_subgroup_asymmetry.py`, `12_comprehensive_asymmetry.py`:

Replace the ~10-line hardcoded dict + `VMAX_GROUPS` + `COLOR_DICT` with:
```python
from pipeline_config import build_bigwig_dict, build_vmax_groups, build_color_dict
BIGWIG_DICT = build_bigwig_dict(BIGWIG_BASE)
```

### Step 6: Fix `08_summary_figures.py` biological hardcoding

4 changes:
- `BIO_ORDER` -> `get_bio_order()`
- `BIO_NAMES` -> `get_bio_names()`
- `'11_Polycomb'` state ID -> `get_key_state_id('11_Polycomb')` with fallback substring search
- Panel titles compute `n=` from data instead of hardcoding `n=667` / `n=2,359`

### Step 7: Replace `parse_header()` in 5 scripts

In `quantify_orientation_asymmetry.py`, `10_clust6_subgroup_asymmetry.py`, `11_histone_anchors_metagene.py`, `12_comprehensive_asymmetry.py`, `visualize_orientation_asymmetry.py`:

Replace local `def parse_header()` with `from pipeline_config import parse_header`.

### Step 8: Replace duplicated constants in remaining scripts

- `SIZE_THRESHOLD` in `10` and `12` -> `from pipeline_config import get_size_threshold`
- `CLUSTER_ORDER` in `05`, `06`, `07`, `09` -> `from pipeline_config import get_cluster_order`
- `BIO_ORDER` / `BIO_LABELS` in `11`, `13`, `visualize_*` -> `from pipeline_config import get_bio_order, get_bio_names`

### Step 9: Fix binary paths in all 8 runner scripts

Replace hardcoded `/opt/homebrew/.../python3` with:
```bash
[ -n "${CLUSTER_CONF}" ] && source "${CLUSTER_CONF}"
PYTHON="${CLUSTER_PYTHON:-$(command -v python3)}"
RSCRIPT="${CLUSTER_RSCRIPT:-$(command -v Rscript)}"
```

Usage becomes: `CLUSTER_CONF=scripts/config/late.conf bash scripts/run_phase3.sh 6`

### Step 10: Add a template config for new projects

Create `scripts/config/template.conf` with all variables documented and example values, serving as the onboarding entry point for new researchers.

---

## Part 3: What NOT to change

- **Popay modules (`modules/`)** — leave untouched. They're stabilized, working, and the coupling surface is narrow.
- **The phase-by-phase execution model** — users need to run phases interactively and inspect results between them. Don't force a monolithic CLI.
- **`STATE_COLORS` / `DIRECTION_COLORS` dicts** — these are rendering constants, not project-specific data. Different scripts use different subsets for different figures; deduplicating would create a confusing god-object.
- **`DIFFBIND_FILES` dict in `05`** — correctly gated behind sub-analysis 4.7; a new researcher without DiffBind data simply skips that analysis.
- **The `run_phase4_9mark.sh` pattern** — this is already the target architecture (env vars set in runner, Python reads them). Leave it as the exemplar.

---

## Summary

| Metric | Current | After refactoring |
|--------|---------|-------------------|
| Files with hardcoded absolute paths | 18 | 0 (all env-var gated) |
| Copy-pasted `BIGWIG_DICT` blocks | 5 | 0 (single source in `pipeline_config`) |
| Copy-pasted `parse_header()` | 5 | 0 (single source in `pipeline_config`) |
| Hardcoded `CLUSTER_ORDER` for k=6 | 6 | 0 (computed from `CLUSTER_K`) |
| Hardcoded biological labels | 6 scripts | 0 (env-var driven, auto-fallback) |
| Files to create | — | 2 (`pipeline_config.py`, `template.conf`) |
| Files to edit | — | 18 (16 scripts + 2 config files) |
| Popay modules touched | — | 0 |

The refactoring is ~300 lines of new code (`pipeline_config.py` + `template.conf`) and ~400 lines of edits across 18 files (mostly replacing hardcoded blocks with 1-2 line imports). Every change is independently testable and backward compatible.

---

## Verification Plan

1. After creating `pipeline_config.py`: import it from a Python shell and verify all functions return current BAP1-KO defaults when no env vars are set.
2. After each script edit: run the corresponding phase with no env vars set and diff outputs against existing outputs — they should be identical.
3. Full integration test: `CLUSTER_CONF=scripts/config/late.conf bash scripts/run_phase3.sh 6` through Phase 4 — verify identical clustering and ChromHMM heatmaps.
4. Portability test: set `CLUSTER_PYTHON=$(which python3)` and `CLUSTER_BIN=$(which cluster)` to verify env-var discovery works.
