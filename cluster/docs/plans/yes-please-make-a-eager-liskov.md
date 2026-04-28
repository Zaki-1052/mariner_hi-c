# Plan: Oriented Anchor Metagene — K27me3 Exterior/Interior Asymmetry Test

## Context

The grad mentor's framework: if H3K27me3 is enriched on the loop-exterior side of an anchor but euchromatic on the loop-interior side (the span), the anchor is a Polycomb extrusion-stopping boundary. If K27me3 is on both sides, the anchor is just embedded in a Polycomb domain. The current Phase 5 deepTools metagene pools both anchor halves into unoriented 3-col BEDs and uses a symmetric ±5kb window — it cannot distinguish exterior from interior.

deepTools `computeMatrix reference-point` natively flips the signal for BED6 minus-strand regions. Since all loops are cis with anchor1 always left of anchor2 (`x1 < y1` invariant), we can assign strand `+` to anchor1 (right/downstream = interior) and strand `-` to anchor2 (deepTools flips so downstream = interior). After the flip, the left side of every plot = **exterior** and right side = **interior** for both anchor types.

**Expected results:**
- clust5 (97% gained, Polycomb domain loops): symmetric K27me3 on both sides
- clust6 (78% lost, euchromatic span): asymmetric K27me3, enriched on exterior only

## Changes

### 1. `cluster/deeptools_plotting.py` — add `xticklabels` param

**Line 15** — add `xticklabels=None` to `heatmap_plot()` signature:
```python
def heatmap_plot(values_path,use_height,pileup_type,color_dict=None,vmax_groups=None,line_measure='mean',bed_color_dict=None,up_down=None,body_length=None,xticklabels=None):
```

**Lines 99-100** — add conditional before the existing label line:
```python
if pileup_type == 'referencePoint':
    line.set_xticks([-n_bins/2 * bin_size,0,n_bins/2 * bin_size])
    if xticklabels is not None:
        line.set_xticklabels(xticklabels)
    else:
        line.set_xticklabels([str((-n_bins/2 * bin_size)/1000) + 'kb',0,'+' + str((n_bins/2 * bin_size)/1000) + 'kb'])
```

Backward-compatible — all existing callers pass no `xticklabels` and get existing behavior.

### 2. `cluster/deepTools_pipeline.py` — pass `xticklabels` through

**Line 25** — add `xticklabels=None` to `bed_pileup()` signature.

**Lines 63-64** — add `xticklabels=xticklabels` to the `heatmap_plot()` call.

### 3. `cluster/scripts/09_oriented_anchor_metagene.py` — new script

Modeled on `06_deeptools_metagene.py`. Key difference: `build_oriented_anchor_beds()` produces BED6 files with strand:
- anchor1 (`chr1/x1/x2`) → strand `+` (downstream = loop interior)
- anchor2 (`chr2/y1/y2`) → strand `-` (deepTools flips → downstream = loop interior)
- Dedup on `(chrom, start, end, strand)` — hub anchors with both orientations keep both entries
- Same 8 BigWigs, same ±5kb window, same `VMAX_GROUPS`/`COLOR_DICT`

Passes `xticklabels=['-5kb (exterior)', 'anchor', '+5kb (interior)']` to `bed_pileup()`.

Output: `cluster/bap1_late/figures/deeptools/oriented_anchors/oriented_anchors_values.{png,pdf,svg,jpg}`
Anchor BEDs: `cluster/bap1_late/figures/deeptools_input/clust{1..6}_oriented_anchors.bed`

### 4. `cluster/scripts/run_oriented_metagene.sh` — phase driver

Same pattern as `run_phase5.sh`: `set -e`, `cd` to repo root, prepend cluster env bin to PATH, `PYTHONUNBUFFERED=1`, tee to `cluster/oriented_metagene.txt`.

## Implementation Order

1. Edit `deeptools_plotting.py` (line 15 + lines 99-100)
2. Edit `deepTools_pipeline.py` (line 25 + lines 63-64)
3. Verify imports: `python3 -c "from deepTools_pipeline import bed_pileup; print('OK')"`
4. Write `scripts/09_oriented_anchor_metagene.py`
5. Write `scripts/run_oriented_metagene.sh`
6. Dry-run BED builder: verify `clust6_oriented_anchors.bed` has 6 columns with `+`/`-` strand
7. Full run: `bash cluster/scripts/run_oriented_metagene.sh` (~90-120 min)

## Verification

1. Anchor BEDs: 6-col, last column is `+` or `-` only. clust5 ≈ 1,198 anchors, clust6 ≈ 3,964
2. Output heatmap x-axis reads `-5kb (exterior)`, `anchor`, `+5kb (interior)`
3. **Biological signal**: K27me3 columns — clust6 should show left-skewed enrichment (exterior > interior); clust5 should show bilateral enrichment (both sides high)
4. Log file `cluster/oriented_metagene.txt` has no Python exceptions
