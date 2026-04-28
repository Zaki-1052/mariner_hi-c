Both diffs are clean. Here's what changed and why:

**`cooltools_called.py`** (the library — fix at the source):
- Added `_viewframe_from_cooler(clr)`: builds a 4-column viewframe (`chrom`, `start`, `end`, `name`) directly from the cooler's own `clr.chromsizes`. This is the single source of truth — the cooler knows exactly which chromosomes it contains and their sizes.
- `mcool_get_matrix`: now calls `_viewframe_from_cooler(clr)` right after opening the cooler, and passes that to both `expected_cis` and `pileup`. The `arms` parameter is kept in the signature (default `None`) for backward compatibility but is ignored.
- `mcool_pileup`: removed the `make_arms(genome=genome)` call and the `arms=arms` plumbing to all 4 `mcool_get_matrix` call sites.
- `make_arms` is still defined at the bottom — dead code now, but removing it could break any external script that imports it.
- Also removed Popay's stale commented-out code block (the `use_cols` and "max intensity" notes).

**`07_cooltools_pileup.py`** (the driver — stripped clean):
- Removed all monkey-patch machinery (`make_arms_robust`, `cooltools_called.make_arms = ...`), the `import bioframe` and `import cooltools_called` that only existed to support the patches.
- Updated docstring and runtime NOTE to reflect that no UCSC fetch is needed.
- The driver is now just what it should be: parse args, load BEDPEs, call `mcool_pileup`.

**Why whole-chromosome arms are correct for mm10:** Mouse chromosomes are acrocentric — centromeres sit at the very start (~0–3 Mb). `make_chromarms` would produce a tiny p-arm and a q-arm that's essentially the whole chromosome. The expected distance-decay at your +/-500 kb pileup scale is dominated by local polymer physics, not centromere effects. The Popay arm-splitting was designed for hg38 (human, metacentric chromosomes) where it actually matters.

---

Phase 6 is working. The pileups look solid:

- **All 6 clusters** show the characteristic loop enrichment (central red dot) with anchor stripe patterns (cross shape) — exactly what you expect from off-diagonal Hi-C pileups.
- **clust1/clust2** (unchanged): strong, symmetric enrichment in both ctrl and mut, consistent with high-signal unchanged loops.
- **clust5** (strong gain, n=667): intense central enrichment — fits the Polycomb domain compaction biology.
- **clust6** (strong loss, n=2,359): enrichment with visible stripes, consistent with anchor disruption.
- **Ctrl vs mut are visually similar** at this scale, which makes sense — the majority of loops are unchanged. The differential signal is in the cluster assignments, not in dramatic pileup differences.
- Color scale log2(obs/exp) [-1, 2] captures the dynamic range well.

Three fixes to get here: viewframe from cooler chromsizes, `clr_weight_name` plumbing, and ICE balancing the mcools. Want me to commit these changes?


rsync -avz expanse:/expanse/lustre/projects/csd940/zalibhai/mariner_hi-c/cluster/logs/ cluster/logs/  