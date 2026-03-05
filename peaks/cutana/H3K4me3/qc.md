## CUT&RUN H3K4me3 FastQC Summary

**Samples:** 4 H3K4me3 (2 ctrl, 2 mut, paired-end) + 1 IgG control (paired-end). All post-trimming (20-101 bp reads).

### Read Counts

| Sample | Reads | %GC | Duplication |
|---|---|---|---|
| IgG (control) | 23.5M | 48% | **62.6%** |
| ctrl_1 (idx25) | 9.5M | 53% | 27.5% |
| mut_1 (idx26) | 6.5M | 51% | 20.3% |
| ctrl_2 (idx27) | 7.7M | 52% | 23.7% |
| mut_2 (idx28) | 7.5M | 55% | 28.1% |

### Module Summary

| Module | Status |
|---|---|
| Base quality, tile quality, seq quality scores | **PASS** all 10 files |
| Per base sequence content | FAIL all 10 -- **expected** (Tn5 insertion bias) |
| Sequence length distribution | WARN all 10 -- **expected** (variable after trimming) |
| GC content | FAIL IgG; WARN ctrl2/mut2; PASS ctrl1/mut1 |
| Duplication | FAIL IgG (62.6%); PASS all H3K4me3 |
| Adapter content | PASS all -- trimming was effective |
| Overrepresented sequences | PASS all -- no contamination |

### Interpretation

**The H3K4me3 libraries look good.** All the flags are expected CUT&RUN artifacts:

- **Per base content FAIL** -- universal Tn5 transposase sequence preference at the 5' end. Not actionable.
- **GC shift in H3K4me3 samples (51-55%)** -- biologically expected since H3K4me3 marks CpG island-rich active promoters. The fact that H3K4me3 samples are GC-shifted relative to IgG (48%) actually indicates good target enrichment.
- **H3K4me3 duplication (20-28%)** -- normal range for CUT&RUN, and tracks with sequencing depth (mut_1 has fewest reads and lowest duplication).

### Flags

1. **IgG duplication at 62.6% -- moderate concern.** IgG CUT&RUN controls are notoriously low-complexity (no specific target, Tn5 cuts only at accessible chromatin). After alignment + deduplication, check how many unique reads remain. If the unique count is very low, the IgG may have limited utility as peak-calling background. This is a common issue with CUT&RUN IgG controls.

2. **mut_2 has the highest %GC (55%)** -- slightly above the other H3K4me3 samples (51-53%). Likely reflects stronger enrichment at GC-rich loci. Worth keeping an eye on in downstream analysis (e.g., does it correlate with peak calls) but not alarming.

3. **mut_1 has the fewest reads (6.5M)** -- about 2/3 of ctrl_1's depth. However, its low duplication (20%) indicates good library complexity, so the library itself is fine -- it just got less sequencing. Depending on your peak calling needs, this may or may not be sufficient depth.

---

## IgG R2 FastQC Results

| Metric | Value |
|---|---|
| Total Sequences | 23,538,581 (matches R1) |
| %GC | 48% |
| Seq Length | 20-101 bp |

| Module | Status |
|---|---|
| Per base seq quality | PASS |
| Per tile seq quality | **WARN** (R1 was PASS) |
| Per seq quality scores | PASS |
| Per base seq content | FAIL (expected, Tn5 bias) |
| Per seq GC content | **FAIL** (same as R1) |
| Per base N content | PASS |
| Seq length distribution | WARN (expected, trimming) |
| Seq duplication levels | **FAIL** -- 60.0% duplication (R1 was 62.6%) |
| Overrepresented sequences | PASS |
| Adapter content | PASS |

### What's different from R1

- **Per tile sequence quality: WARN** -- R1 was PASS. This means some flowcell tiles had slightly lower quality on R2. This is common since R2 is sequenced second and quality degrades slightly. Not concerning unless it's a FAIL.

- **Duplication: 60.0%** -- Consistent with R1's 62.6%. Confirms the IgG library itself is low-complexity (not a read-direction artifact).

- **Notable duplication breakdown:** 18.7% of sequences are duplicated >100x, and 21.8% are >500x. This is extreme and means a large fraction of IgG reads come from a small number of genomic loci — consistent with Tn5 cutting the same accessible sites repeatedly in the absence of antibody targeting.

### Updated assessment

No change to the overall picture. The IgG R2 matches R1 closely, which is expected. The IgG high-duplication flag still stands — after alignment and deduplication, you'll want to check how many unique fragments remain to assess whether this IgG control is usable for background normalization in peak calling.