# .hic Binary Header Format: Reverse-Engineered Reference

## Document Purpose

This document records the actual binary layout of .hic file headers as determined by manual hex-dumping and struct probing of real v9 files from this project. There is no formal .hic specification — the canonical reference is the Juicer Java source code, which changes without notice. This document captures what the bytes actually contain, with worked examples from our BAP1-KO Hi-C data.

**Key finding:** The v9 header layout differs from v7/v8 in two ways that existing open-source parsers (hic2cool, mustache) do not handle. These differences cause complete parse failures that cascade into garbage output for every field after byte 16.

---

## 1. Byte-Level Layout

### v7/v8 Layout (what most open-source parsers implement)

```
Offset  Bytes  Type        Field
──────  ─────  ──────────  ─────────────────────────────
0       4      char[3]+\0  Magic string "HIC\0"
4       4      int32 LE    Version (7 or 8)
8       8      int64 LE    Master index offset (footer position)
16      var    cstring     Genome ID (null-terminated)
?       4      int32 LE    nAttributes
?       var    cstring[]   Key-value pairs (nAttributes × 2 strings)
?       4      int32 LE    nChromosomes
?       var    (cstring + int32)[]   Chromosome name + length pairs
?       4      int32 LE    nBpResolutions
?       var    int32[]     BP resolution values
?       4      int32 LE    nFragResolutions
?       var    int32[]     Fragment resolution values
```

### v9 Layout (what the files actually contain)

```
Offset  Bytes  Type        Field
──────  ─────  ──────────  ─────────────────────────────
0       4      char[3]+\0  Magic string "HIC\0"
4       4      int32 LE    Version (9)
8       8      int64 LE    Master index offset (footer position)
16      var    cstring     Genome ID (null-terminated)
?       8      int64 LE    nviPosition  ←── NEW: normalization vector index offset
?       8      int64 LE    nviLength    ←── NEW: normalization vector index size
?       4      int32 LE    nAttributes
?       var    cstring[]   Key-value pairs (nAttributes × 2 strings)
?       4      int32 LE    nChromosomes
?       var    (cstring + int64)[]   Chromosome name + length pairs  ←── int64, not int32
?       4      int32 LE    nBpResolutions
?       var    int32[]     BP resolution values
?       4      int32 LE    nFragResolutions
?       var    int32[]     Fragment resolution values
```

**Two differences from v7/v8:**

1. **nviPosition + nviLength (2 × int64)** inserted between the genome string and the attribute count. "nvi" = Normalization Vector Index — these point to where the normalization vector metadata is stored in the file footer, separate from the main master index.

2. **Chromosome lengths are int64** instead of int32. The maximum int32 value (~2.1 billion) is sufficient for all current genomes, but v9 switched to int64 for future-proofing.

---

## 2. Worked Example: late-ctrl_merged.hic (250402)

### File metadata

| Property | Value |
|---|---|
| File | `late-ctrl_merged.hic` (250402 timepoint, ctrl merged) |
| Size | 6,754,817,391 bytes (6.29 GB) |
| Created by | Juicer Tools Version 2.20.00 |
| Genome ref | mm10 (mouse) |
| Normalizations | VC, VC_SQRT, SCALE only (**no KR, no GW_SCALE, no INTER_SCALE**) |

This is the file that caused the Stripenn pipeline 250402 normalization incident (see `stripes/stripenn/README.md` Section 5.2).

### Raw hex dump (first 216 bytes)

```
Offset  Hex                                              ASCII
──────  ───────────────────────────────────────────────  ────────────────
    0:  48 49 43 00 09 00 00 00 97 5d 90 8f 01 00 00 00  HIC.......].....
   16:  2f 65 78 70 61 6e 73 65 2f 6c 75 73 74 72 65 2f  /expanse/lustre/
   32:  70 72 6f 6a 65 63 74 73 2f 63 73 64 39 34 30 2f  projects/csd940/
   48:  63 74 65 61 2f 48 69 43 2f 72 65 73 74 72 69 63  ctea/HiC/restric
   64:  74 69 6f 6e 5f 73 69 74 65 73 2f 6d 6d 31 30 2e  tion_sites/mm10.
   80:  66 69 6c 74 65 72 65 64 2e 63 68 72 6f 6d 2e 73  filtered.chrom.s
   96:  69 7a 65 73 00 da 29 d4 8f 01 00 00 00 79 7b 00  izes..)......y{.
  112:  00 00 00 00 00 02 00 00 00 73 6f 66 74 77 61 72  .........softwar
  128:  65 00 4a 75 69 63 65 72 20 54 6f 6f 6c 73 20 56  e.Juicer Tools V
  144:  65 72 73 69 6f 6e 20 32 2e 32 30 2e 30 30 00 68  ersion 2.20.00.h
  160:  69 63 46 69 6c 65 53 63 61 6c 69 6e 67 46 61 63  icFileScalingFac
  176:  74 6f 72 00 31 2e 30 00 16 00 00 00 41 6c 6c 00  tor.1.0.....All.
  192:  91 96 29 00 00 00 00 00 63 68 72 31 00 63 aa a6  ..).....chr1.c..
  208:  0b 00 00 00 00 63 68 72                           .....chr
```

### Annotated field-by-field parse

```
Offset    0 [  4 bytes]: Magic = "HIC\0"
Offset    4 [  4 bytes]: Version = 9
Offset    8 [  8 bytes]: Master Index = 6,703,570,327 (0x18f905d97)
Offset   16 [ 85 bytes]: Genome = "/expanse/lustre/projects/csd940/ctea/HiC/restriction_sites/mm10.filtered.chrom.sizes"
Offset  101 [  8 bytes]: nviPosition = 6,708,013,530 (0x18fd429da)    ← v9 only
Offset  109 [  8 bytes]: nviLength = 31,609                           ← v9 only
Offset  117 [  4 bytes]: nAttributes = 2
                          "software" = "Juicer Tools Version 2.20.00"
                          "hicFileScalingFactor" = "1.0"
Offset  184 [  4 bytes]: nChromosomes = 22
                          All    =     2,725,521 bp    ← pseudo-chromosome
                          chr1   =   195,471,971 bp
                          chr2   =   182,113,224 bp
                          chrX   =   171,031,299 bp    ← 3rd by size, not alphabetical
                          chr3   =   160,039,680 bp
                          chr4   =   156,508,116 bp
                          chr5   =   151,834,684 bp
                          chr6   =   149,736,546 bp
                          chr7   =   145,441,459 bp
                          chr10  =   130,694,993 bp
                          chr8   =   129,401,213 bp
                          chr14  =   124,902,244 bp
                          chr9   =   124,595,110 bp
                          chr11  =   122,082,543 bp
                          chr13  =   120,421,639 bp
                          chr12  =   120,129,022 bp
                          chr15  =   104,043,685 bp
                          chr16  =    98,207,768 bp
                          chr17  =    94,987,271 bp
                          chrY   =    91,744,698 bp
                          chr18  =    90,702,639 bp
                          chr19  =    61,431,566 bp
                          Total (excl All): 2,725,521,370 bp
Offset  483 [  4 bytes]: nBpResolutions = 10
                          [1000, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000, 2500000]
Offset  527 [  4 bytes]: nFragResolutions = 0
Offset  531:              ← header ends, data region begins
```

---

## 3. Comparison: 250402 (late) vs 250831 (early)

Both files are v9, created by Juicer Tools 2.20.00 against the same mm10 chrom.sizes, with identical chromosome lists, ordering, and resolutions. The differences are in size and normalization:

| Field | late-ctrl (250402) | early-ctrl (250831) |
|---|---|---|
| File size | 6.29 GB | 8.11 GB |
| Master index | 6,703,570,327 | 8,651,922,151 |
| nviPosition | 6,708,013,530 | 8,654,083,789 |
| **nviLength** | **31,609** | **6,304** |
| Genome string | identical | identical |
| Chromosomes | identical (22, same order) | identical |
| Resolutions | identical (10 BP, 0 FRAG) | identical |
| Available norms | VC, VC_SQRT, SCALE | **KR only** |

The **nviLength** difference is the most informative: 31,609 bytes for the file with 3 normalization types vs 6,304 bytes for the file with 1. The normalization vector index (nvi) stores metadata about where each norm vector lives in the file — more norms = larger index. The ~5:1 ratio (not exactly 3:1) makes sense because each norm type has vectors at multiple resolutions.

### Normalization deep audit results

Tested with hicstraw 1.3.1 via `hic_audit.py --deep`:

| Norm | 250402 ctrl | 250831 ctrl | Interpretation |
|---|---|---|---|
| NONE (raw) | PASS | PASS | Raw contact data always accessible |
| KR | **FAIL** | PASS | 250402 was created without `-k KR` flag |
| VC | PASS | **FAIL** | 250831 was created WITH `-k KR` (KR-only, no VC) |
| VC_SQRT | PASS | **FAIL** | Same — excluded by `-k KR` |
| SCALE | PASS | **FAIL** | Same |
| GW_SCALE | **FAIL** | **FAIL** | Missing from both |
| INTER_SCALE | **FAIL** | **FAIL** | Missing from both |

**Root cause:** The 250402 merged .hic files were created by a different lab member using `juicer pre` without the `-k KR` flag, producing the default normalization set (VC, VC_SQRT, SCALE) but omitting KR. The 250831 merged files were created with `-k KR`, producing KR only. Neither set has GW_SCALE or INTER_SCALE — these are apparently not computed by `juicer pre` for merged files (only for initial .hic creation from the full pipeline).

---

## 4. Surprising Findings

### 4.1 The genome field is a filesystem path, not an assembly name

Expected: `"mm10"` or `"GRCm38"`.

Actual: `"/expanse/lustre/projects/csd940/ctea/HiC/restriction_sites/mm10.filtered.chrom.sizes"`

The genome field stores whatever string was passed to `juicer pre` or `juicer_tools pre` as the genome/chrom.sizes argument. For our lab's pipeline, that's the full HPC path to the chrom.sizes file. This means:
- You cannot reliably determine the genome assembly from this field alone.
- Different labs' files for the same genome will have different genome strings.
- The field is useful for provenance ("which chrom.sizes was used?") but not for automated genome identification. You'd need to parse the filename or match chromosome sizes against known assemblies.

### 4.2 The "All" pseudo-chromosome encodes genome size in kilobases

The chromosome list always includes a pseudo-entry called "All" as the first element:

```
All = 2,725,521 bp
Total genome (excl All) = 2,725,521,370 bp
Ratio: All / Total = 0.001000 (exactly 1/1000)
```

The "All" length is exactly the total genome size divided by 1000 — the genome size in kilobases. This is used internally by Juicer for genome-wide matrix calculations. It should be filtered out when listing actual chromosomes (hicstraw's `getChromosomes()` includes it; tools should skip entries named "All").

### 4.3 Chromosomes are ordered by size, not naturally

```
chr1 (195M) → chr2 (182M) → chrX (171M) → chr3 (160M) → chr4 (156M) → ...
```

chrX appears between chr2 and chr3 because it's the 3rd longest mouse chromosome. chrY appears between chr17 and chr18. This is the order Juicer stores them internally (descending by size). Tools that expect natural chromosome ordering (chr1, chr2, chr3, ..., chrX, chrY) need to sort after reading.

### 4.4 Parsing v9 with a v7/v8 parser produces catastrophic garbage

When we first ran the mustache-ported `read_header()` (which handles v7/v8) against these v9 files, the output was:

```
Resolutions (BP): -2131067586, -2124399734, -2122678964, ...
```

Hundreds of wildly negative numbers instead of the expected 10 resolutions. The cascade:

1. Parser reads bytes 101-108 (nviPosition) as `nAttributes` (int32) → gets `0x8fd429da` = -1,882,051,110 (a negative "attribute count")
2. Actually, struct interprets this as a large positive number or negative depending on signed/unsigned, but the loop tries to read that many key-value pairs
3. Every subsequent field read is misaligned, reading chromosome data as resolutions, chromosome names as attribute keys, etc.
4. The parser doesn't crash — it happily produces thousands of garbage "resolutions" from misinterpreted chromosome length bytes

This is why tools that only support v7/v8 fail on v9: it's not just a version check — the binary layout genuinely shifts by 16 bytes (the two nvi int64 fields), and chromosome lengths double in width from 4 to 8 bytes. Every field after offset 101 is misaligned.

### 4.5 hicstraw silently falls back to raw data when norms are missing

When requesting KR-normalized data from a file without KR vectors:

```python
mzd = hic.getMatrixZoomData("chr1", "chr1", "observed", "KR", "BP", 5000)
records = mzd.getRecords(0, 1_000_000, 0, 1_000_000)
# records is non-empty! No exception raised!
```

hicstraw prints a warning to C-level stderr:
```
File did not contain KR normalization vectors for one or both chromosomes at 5000 BP
```

But the Python call succeeds and returns data — raw (unnormalized) contacts. This is a silent failure: code that expects KR-normalized data will instead process unnormalized data without any Python-level error. This is how a lab could unknowingly analyze unnormalized Hi-C data.

The `hic_audit.py --deep` mode catches this by capturing C-level stderr during each `getMatrixZoomData()` call and checking for "did not contain" in the warning text.

### 4.6 nviLength correlates with normalization count

| File | Norms available | nviLength |
|---|---|---|
| 250402 ctrl | 3 (VC, VC_SQRT, SCALE) | 31,609 |
| 250831 ctrl | 1 (KR only) | 6,304 |

Ratio: 31,609 / 6,304 = 5.01×. The normalization vector index stores per-resolution metadata for each norm type. With 10 resolutions and 3 norms vs 1 norm, the expected ratio would be ~3× just from norm count, but each resolution also has per-chromosome entries which adds overhead non-linearly.

This field could be used as a quick heuristic: nviLength < 10,000 likely means only 1 normalization type is present.

---

## 5. Footer and Normalization Vectors (Not Fully Mapped)

The header points to two locations in the file:

1. **Master index** (offset from header field at byte 8): Points to the main footer which contains the matrix block index — where the contact data for each chromosome pair at each resolution is stored.

2. **nviPosition** (v9 only, offset from header): Points to the normalization vector index — where the normalization vectors (KR, VC, etc.) for each resolution are stored.

The footer and normalization vector formats are where v9 diverges most from v7/v8, and where the cross-tool incompatibilities (Class 3 in HIC-PLAN.md) originate. Full mapping of the footer/nvi structure is out of scope for this document — it would require comparing the Juicer Java source code across versions.

What we know:
- The footer is near the end of the file (master index is ~99% of file size in our files).
- The nvi is after the master index (nviPosition > master index in all test files).
- The gap between master index and nvi (nviPosition - masterIndex) is small: ~4.4 MB for 250402, ~2.2 MB for 250831.
- The file ends shortly after the nvi (nviPosition + nviLength ≈ file size).

---

## 6. Field Reference

### Magic String (4 bytes)

Always `HIC\0` (0x48 0x49 0x43 0x00). All versions.

### Version (int32)

Known values: 6, 7, 8, 9. Our files are all 9 (Juicer Tools 2.20.00).

### Master Index (int64)

Byte offset in the file where the main footer begins. Used by readers to seek to the block index for finding matrix data. Always a large number near the end of the file.

### Genome ID (cstring)

Null-terminated string. In practice, this is whatever was passed as the genome argument to `juicer pre` or `juicer_tools pre`. May be:
- A genome assembly name: `"mm10"`, `"hg38"`
- A path to a chrom.sizes file: `"/path/to/mm10.chrom.sizes"`
- Something else entirely, depending on the tool that created the file

### nviPosition / nviLength (v9 only, 2 × int64)

Normalization Vector Index location. Points to a section of the file that describes where each normalization vector (KR, VC, etc.) is stored for each resolution. The "nvi" name is inferred from the Juicer source variable naming.

nviLength correlates with the number of normalization types × resolutions stored.

### nAttributes (int32)

Number of metadata key-value pairs. Observed values: 2 (in all our files).

Known attributes:
- `software`: The tool and version that created the file (e.g., "Juicer Tools Version 2.20.00")
- `hicFileScalingFactor`: Scaling factor, typically "1.0"

### Chromosomes (nChromosomes × (cstring + int32|int64))

List of chromosome entries. **int32 lengths in v7/v8, int64 in v9.**

Always includes a pseudo-chromosome "All" as the first entry, whose length equals the total genome size in kilobases (total bp / 1000). Real chromosomes follow, ordered by descending length (not alphabetically).

### BP Resolutions (nBpResolutions × int32)

Base-pair resolutions available in the file. Standard Juicer output includes: 1000, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000, 2500000.

### Fragment Resolutions (nFragResolutions × int32)

Restriction fragment resolutions. Typically 0 in modern files (fragment-level resolution is rarely used). Legacy field from early Juicer versions.

---

**Created:** 2026-04-22
**Test files:** BAP1-KO Hi-C merged .hic files (mm10, Juicer Tools 2.20.00, v9)
**Parser:** `hic/hic_audit.py` → `parse_hic_header()`
