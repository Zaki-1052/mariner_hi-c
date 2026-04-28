# Plan: Phase 1 — Data Preparation for Popay Clustering Pipeline

## Context

We are reproducing Popay et al. (Nat Genet 2026) Figure 2f — ChromHMM anchor-vs-span enrichment across k-means clusters of Hi-C loops — adapted to BAP1-KO mouse cerebellum (mm10) at the late/adult timepoint (250402). The Dixon meeting (2026-04-10) flagged this as the highest-priority outstanding item for the paper because it directly tests the "extrusion impediment" model (Polycomb spreading flanking anchors blocks cohesin) vs. the "sensitivity" model (Polycomb at anchors disrupts CTCF binding).

**This plan covers Phase 1 only**: build a Popay-format loop count file plus its differential-stats sidecar, generate an mm10 promoter BED matching her gene annotation format, and create the project's output directory tree. Phases 0 (bug fixes) is already applied; Phases 2–6 will follow once Phase 1 is verified.

The full project context lives in `cluster/CONTEXT-CLUSTER.md`; the multi-phase pipeline plan in `cluster/PLAN-CLUSTER.md`. This plan is the executable spec for Phase 1, with corrections from data-file inspection.

---

## Critical findings from data inspection (correct the parent plan)

1. **`merged_all_loops_nonredundant.bedpe` schema** (32 columns, 39,344 data rows, header present)
   - Columns 1–6: `chr1, start1, end1, chr2, start2, end2` (verified)
   - Column 30: `kept_from_resolution` ∈ {5, 10, 25} (kb units)
   - Distribution: 5kb=7,901 / 10kb=14,553 / 25kb=16,890 (sums to 39,344 ✓)
   - Other useful columns: `coord_string` (col 7), `logFC` (col 9), `FDR` (col 13), `direction` (col 17), `resolution` (col 19), `is_multi_resolution` (col 23)
   - Row 2 of the merged BEDPE matches `loop_2` of `res_5kb/all_results_primary.tsv` exactly (chr8 110165000–110170000 / 110475000–110480000), confirming coordinates are preserved verbatim — joins are clean.

2. **Per-resolution count matrices** (`res_{5,10,25}kb/06_counts_matrix.tsv`)
   - Header row: empty index col, then `ctrl_M1\tctrl_M2\tctrl_M3\tmut_M1\tmut_M2\tmut_M3`
   - Row counts: 17,983 (5kb), 22,639 (10kb), 20,421 (25kb)
   - Index format: `loop_1`, `loop_2`, …
   - Values are float (mariner-aggregated counts from 5×5 pixel windows), not raw integers.

3. **Per-resolution edgeR results** (`edgeR_results_res_{5,10,25}kb/primary_analysis/all_results_primary.tsv`)
   - 19 columns starting with `loop_id, chr1, start1, end1, chr2, start2, end2, coord_string, …`
   - Row counts: 17,982 / 22,632 / 20,398 (slightly fewer than the count matrix because `filterByExpr` removed low-count loops). All loops in the merged BEDPE are present here — the edgeR pass is upstream of merging.

4. **Popay's INPUT format** (verified from `cluster/clustering_example_data/RAD21_dependent.txt` and `cluster/HiC_cluster3.ipynb`)
   ```
   chr1   x1   x2   chr2   y1   y2   DMSO_merge   4hr_merge   24hr_merge
   ```
   - Tab-separated, header present
   - `data_cols = [c for c in cols if c not in ['chr1','x1','x2','chr2','y1','y2']]`
   - Loaded by `pd.read_csv(loop_count_file, sep='\t', header=0)`
   - Notebook downstream code does `treatment_group.str.replace('_merge','')` — **the suffix MUST be exactly `_merge`, not `_merged`**, otherwise display labels become `ctrld`/`mutd`.

5. **Count-scale heterogeneity across resolutions**
   - 5kb median ctrl_avg = 171, 5%-tile = 59
   - 10kb median ctrl_avg = 403, 5%-tile = 111
   - 25kb median ctrl_avg = 984, 5%-tile = 176
   - This means a single absolute-value filter (`ctrl_merge > X`) would disproportionately remove 5kb loops. The Popay normalization step (each row divided by its own `ctrl_merge`) handles per-loop scale, but the upstream filter does not. **Decision for Phase 1**: write raw averaged counts faithfully and defer the filter-threshold decision to Phase 3 where we can profile the distribution — possibly per-resolution. Phase 1 must NOT silently apply any filter or scaling.

6. **Popay's 7-column promoter BED** (`gencode.v25.annotation_pp.bed`)
   ```
   chr1   10869   12369   ENSG00000223972.5   0   +   DDX11L1
   ```
   - 1500-bp width = TSS ± 750 bp
   - Columns: chr, start, end, gene_id, score (=0), strand, gene_name
   - `bedpe_analysis.py` line 70 confirms: `B_bed_cols = ['chr','start','end','gene_id','idk','strand','gene_name']`. The default path at line 67 is already fixed to `<module_dir>/data/mm10_knownGene_pp.bed`.

7. **Environment status**
   - System R 4.5.2 at `/usr/local/bin/R` has `TxDb.Mmusculus.UCSC.mm10.knownGene`, `org.Mm.eg.db`, `GenomicRanges`, `GenomicFeatures` — all four `requireNamespace()` calls returned TRUE. The CONTEXT-CLUSTER claim that they live in `mariner_env` is incorrect on this Mac (no `mariner_env` exists in `conda env list`); they are in the system R library at `/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library`. **Use system R for the annotation script** — no conda activation needed.
   - `cluster` conda env at `/opt/homebrew/anaconda3/envs/cluster/` has pandas 1.5.3, bioframe 0.6.1, pyBigWig, numpy, sklearn. Direct invocation: `/opt/homebrew/anaconda3/envs/cluster/bin/python3`. (`conda run -n cluster` does NOT activate properly here — use the absolute path or `source <hook> && conda activate cluster`.)
   - Phase 0 module fixes are applied (verified by reading `cluster/bedpe_analysis.py` — line 67 default path already points to repo-relative `data/mm10_knownGene_pp.bed`; CONTEXT-CLUSTER §0 records "Verified 2026-04-26").

---

## What Phase 1 produces

| Artifact | Path | Rows | Purpose |
|---|---|---|---|
| Loop count file (Popay format) | `cluster/data/late_merged_loop_counts.txt` | 39,344 + 1 header | Direct input to `HiC_cluster3` notebook / `04_clustering.py` |
| Differential-stats sidecar | `cluster/data/late_merged_loop_metadata.tsv` | 39,344 + 1 header | Phase 4.8 cluster × differential-status crosstab |
| mm10 promoter BED | `cluster/data/mm10_knownGene_pp.bed` | ~24K rows (whatever TxDb yields) | Used by `bedpe_analysis.bedtools_annotation` (Phase 4.6) and CRE classifier (Phase 4.2) |
| Output directory skeleton | `cluster/bap1_late/{cluster3, figures/{annotation,ChIP_intersect,deeptools,deeptools_input}, chromHMM/{binarized,learned_model,anchor_input,span_input,peak_beds}, cooltools}` | — | Holds all downstream outputs |

---

## Step 1.1 — Build the loop count file

**Script:** `cluster/scripts/01_build_loop_count_file.py`
**Run with:** `/opt/homebrew/anaconda3/envs/cluster/bin/python3 cluster/scripts/01_build_loop_count_file.py`

### Algorithm

1. **Load merged BEDPE** (`outputs/250402-late_outputs/bedpe_final/merged_all_loops_nonredundant.bedpe`, `sep='\t'`, `header=0`). Verify row count == 39,344 — fail loudly otherwise.

2. **For each resolution X in {5, 10, 25}:**
   a. Load `outputs/250402-late_outputs/edgeR_results_res_{X}kb/primary_analysis/all_results_primary.tsv` → keep `[loop_id, chr1, start1, end1, chr2, start2, end2]`.
   b. Load `outputs/250402-late_outputs/res_{X}kb/06_counts_matrix.tsv` (`index_col=0`) → keep `[ctrl_M1, ctrl_M2, ctrl_M3, mut_M1, mut_M2, mut_M3]`.
   c. Compute `ctrl_merge = mean(ctrl_M1, ctrl_M2, ctrl_M3)`, `mut_merge = mean(mut_M1, mut_M2, mut_M3)` per row using `df[[...]].mean(axis=1)`.
   d. Inner-join the two on `loop_id` (count-matrix index ↔ edgeR `loop_id` column). Result: a per-resolution DataFrame keyed by `(chr1, start1, end1, chr2, start2, end2)` with `ctrl_merge` and `mut_merge`.

3. **Build the join key on the merged BEDPE side.** Key columns: `(chr1, start1, end1, chr2, start2, end2)` — must coerce numeric coords to int explicitly to avoid float-vs-int mismatches across resolutions. Use the merged BEDPE's `kept_from_resolution` value to decide which per-resolution table to look up. Implementation: split the merged BEDPE into three resolution-specific frames, do per-resolution coord-keyed merges, then concatenate. Avoids ambiguity vs. a single mass merge.

4. **Verify completeness.** Every merged loop must produce a non-null `ctrl_merge` and `mut_merge`. Any unmatched row indicates a coordinate-encoding mismatch (e.g., int vs. float, or a row that exists in one source but not the other) — **fail with a diagnostic listing the first 10 unmatched rows**, do not silently fill with NaN. This is a fail-fast invariant per the project's "no fallbacks in scientific code" rule.

5. **Rename output columns to Popay's convention:** `start1→x1, end1→x2, start2→y1, end2→y2`. The data-column suffix is `_merge` (not `_merged`) so the notebook's `.str.replace('_merge','')` produces clean `ctrl`/`mut` labels downstream.

6. **Write the count file:** `cluster/data/late_merged_loop_counts.txt`, columns in this order:
   ```
   chr1   x1   x2   chr2   y1   y2   ctrl_merge   mut_merge
   ```
   Tab-separated, header on, no index. Use Unix LF line endings (per repo CLAUDE.md and the user's encoding feedback). In pandas: `to_csv(..., sep='\t', index=False, lineterminator='\n')`.

7. **Write the metadata sidecar** (`cluster/data/late_merged_loop_metadata.tsv`) with one row per loop, joinable to the count file by `(chr1, x1, x2, chr2, y1, y2)`. Columns:
   ```
   chr1  x1  x2  chr2  y1  y2  logFC  FDR  PValue  direction  significant  category  resolution  resolution_kb  kept_from_resolution  is_multi_resolution
   ```
   Pulled directly from the merged BEDPE columns. Used by Phase 4.8 to cross-tabulate cluster assignments against differential status.

### Sanity checks the script must run before exiting (and print to stdout)

- Output row count == 39,344
- `ctrl_merge.min() > 0` and `mut_merge.min() > 0` (Hi-C aggregated counts are non-negative; zeros warrant investigation)
- `ctrl_merge.isna().sum() == 0` and same for `mut_merge`
- Per-resolution row counts match expectation: 7,901 / 14,553 / 16,890 → from a left-join check on `kept_from_resolution`
- Print summary: `ctrl_merge.describe(percentiles=[0.01,0.05,0.50,0.95,0.99])` so we can choose the Phase 3 filter threshold from real data later.

---

## Step 1.2 — mm10 promoter BED

**Script:** `cluster/scripts/02_build_mm10_gene_annotation.R`
**Run with:** `/usr/local/bin/Rscript cluster/scripts/02_build_mm10_gene_annotation.R`

### Algorithm

1. Load `TxDb.Mmusculus.UCSC.mm10.knownGene` and `org.Mm.eg.db`.
2. `genes(txdb)` → `GRanges` of gene bodies (Entrez IDs as names).
3. `promoters(genes_gr, upstream=750, downstream=750)` → 1500-bp windows centered on each gene's TSS, strand-aware (matching the width of Popay's gencode-v25 file).
4. Map Entrez IDs → MGI symbols via `mapIds(org.Mm.eg.db, keys=names(promoters), keytype='ENTREZID', column='SYMBOL', multiVals='first')`. Genes with no symbol mapping: keep the Entrez ID in `gene_id`, leave `gene_name` blank — do not drop them (`bedpe_analysis.py` tolerates blanks; dropping silently would bias annotation).
5. Filter to standard chromosomes (chr1..19, chrX, chrY, chrM) using `seqlevels()` — drops random/alt contigs that won't appear in our Hi-C data. Coerce `seqlevelsStyle()` to `"UCSC"` so we keep `chr` prefixes.
6. Sort by `(chrom, start)` (bedtools sortable).
7. Output 7-column BED, no header, tab-separated:
   ```
   chr  start  end  gene_id  0  strand  gene_name
   ```
   Score column is the literal integer 0 (matching Popay).
   Path: `cluster/data/mm10_knownGene_pp.bed`.

### Sanity checks

- Row count is in the 20–60K range (typical for `TxDb.Mmusculus.UCSC.mm10.knownGene`)
- All `end - start == 1500`
- All `chrom` values match `^chr(\d+|X|Y|M)$`
- No row has both `gene_id` and `gene_name` blank
- File loads without error via `pd.read_csv(..., sep='\t', header=None, names=['chr','start','end','gene_id','score','strand','gene_name'])` (smoke test from Python)

---

## Step 1.3 — Output directory tree

Create only the empty directories Phase 1 doesn't fill itself; Phase 2+ scripts create their own subdirs.

```bash
mkdir -p cluster/data
mkdir -p cluster/bap1_late/cluster3
mkdir -p cluster/bap1_late/figures/{annotation,ChIP_intersect,deeptools,deeptools_input}
mkdir -p cluster/bap1_late/chromHMM/{binarized,learned_model,anchor_input,span_input,peak_beds}
mkdir -p cluster/bap1_late/cooltools
mkdir -p cluster/scripts
```

(The `cluster/scripts/` mkdir is technically the home of the scripts being created in Phase 1 — include it for safety.)

---

## Files to create

| File | Purpose | Lang |
|---|---|---|
| `cluster/scripts/01_build_loop_count_file.py` | Build Popay-format count file + metadata sidecar | Python (cluster env) |
| `cluster/scripts/02_build_mm10_gene_annotation.R` | mm10 promoter BED, 7-col, TSS ±750 bp | R (system R 4.5.2) |
| `cluster/scripts/run_phase1.sh` | Driver: mkdir + run both scripts; tail with `wc -l` and `head` of outputs for at-a-glance verification | Bash (no `set -euo pipefail` per user feedback) |

## Files NOT to modify

- All seven `cluster/*.py` modules (Phase 0 fixes already applied per CONTEXT-CLUSTER §0 verification — do not re-touch)
- `cluster/CONTEXT-CLUSTER.md` and `cluster/PLAN-CLUSTER.md` (parent docs; update separately if execution surfaces new issues)
- Any pipeline-level `outputs/250402-late_outputs/...` files (read-only inputs)

## Existing utilities to reuse (do not reimplement)

- `cluster/cluster_tools.py` `comparison_type()` — will be invoked in Phase 3 with our `data_cols=['ctrl_merge','mut_merge']`. We verified it returns `'multiple'` for that input (treatments={'ctrl','mut'}, len(treatments)==len(data_cols)==2 → 'multiple'). No code changes needed for our 2-condition design.
- `cluster/bedpe_analysis.py` `bedtools_annotation()` — already points at our future `cluster/data/mm10_knownGene_pp.bed`. Phase 1.2's BED schema is engineered to match its `B_bed_cols` exactly.
- pandas / bioframe / pyBigWig from the existing `cluster` conda env — no new packages required for Phase 1.

---

## Verification plan

After running Phase 1:

```bash
# 1. Loop count file: shape, sample, distribution
wc -l cluster/data/late_merged_loop_counts.txt   # expect 39345 (header + 39344 rows)
head -3 cluster/data/late_merged_loop_counts.txt # expect: chr1 x1 x2 chr2 y1 y2 ctrl_merge mut_merge
/opt/homebrew/anaconda3/envs/cluster/bin/python3 -c "
import pandas as pd
df = pd.read_csv('cluster/data/late_merged_loop_counts.txt', sep='\t')
print('cols:', list(df.columns))
print('shape:', df.shape)
print('NaN counts:', df.isna().sum().to_dict())
print(df[['ctrl_merge','mut_merge']].describe(percentiles=[0.01,0.05,0.50,0.95,0.99]))
"

# 2. Metadata sidecar joins cleanly
/opt/homebrew/anaconda3/envs/cluster/bin/python3 -c "
import pandas as pd
counts = pd.read_csv('cluster/data/late_merged_loop_counts.txt', sep='\t')
meta = pd.read_csv('cluster/data/late_merged_loop_metadata.tsv', sep='\t')
joined = counts.merge(meta, on=['chr1','x1','x2','chr2','y1','y2'], how='inner')
assert len(joined) == 39344, f'Inner join lost rows: {len(joined)}'
print('Sidecar joins on all', len(joined), 'rows.')
print('Direction counts:', meta['direction'].value_counts().to_dict())
"

# 3. Promoter BED: shape, schema, smoke test through bedpe_analysis loader
wc -l cluster/data/mm10_knownGene_pp.bed         # expect ~24-60K
awk -F'\t' 'NF!=7{print "row", NR, "has", NF, "fields"; exit 1} {if ($3-$2 != 1500) {print "row", NR, "width =", $3-$2; exit 1}} END{print "all rows: 7 cols, 1500bp width"}' cluster/data/mm10_knownGene_pp.bed
/opt/homebrew/anaconda3/envs/cluster/bin/python3 -c "
import pandas as pd
bed = pd.read_csv('cluster/data/mm10_knownGene_pp.bed', sep='\t', header=None,
                  names=['chr','start','end','gene_id','score','strand','gene_name'])
print('shape:', bed.shape, 'unique chroms:', sorted(bed['chr'].unique()))
print('strand counts:', bed['strand'].value_counts().to_dict())
"

# 4. Notebook smoke test: confirm Popay's loading code accepts our file unchanged
/opt/homebrew/anaconda3/envs/cluster/bin/python3 -c "
import pandas as pd
df = pd.read_csv('cluster/data/late_merged_loop_counts.txt', sep='\t', header=0)
data_cols = [c for c in df.columns if c not in ['chr1','x1','x2','chr2','y1','y2']]
print('data_cols:', data_cols)              # expect ['ctrl_merge', 'mut_merge']
import sys; sys.path.insert(0, 'cluster')
from cluster_tools import comparison_type
print('comparison_type:', comparison_type(data_cols))   # expect ('multiple', ['ctrl','mut'] in some order)
"
```

The Phase 1 deliverables are correct iff:
- Count file has 39,344 rows × 8 columns, no NaN, columns named `ctrl_merge`/`mut_merge` (not `_merged`)
- Metadata sidecar inner-joins to 39,344 rows on coordinates
- Promoter BED has 7 fields per row, all rows 1500 bp wide, only standard chromosomes
- `comparison_type(['ctrl_merge','mut_merge'])` returns `'multiple'`

If any check fails, the script is wrong — do not patch outputs to make checks pass.

---

## Out of scope for Phase 1 (deferred)

- **Filter threshold selection** — defer to Phase 3 after we see the empirical `ctrl_merge` distribution in the actual count file. The cross-resolution count-scale heterogeneity (5kb median 171, 25kb median 984) likely needs either a percentile-based or per-resolution filter; that decision belongs to clustering, not data prep.
- **Any normalization or scaling** of the count values themselves. Popay's normalize-by-`ctrl_merge` step happens inside the clustering notebook. Phase 1 outputs raw averaged counts.
- **ChromHMM peak BED staging** — Phase 2.
- **k-means clustering** — Phase 3.
- **All Phase 4–6 downstream analyses.**
