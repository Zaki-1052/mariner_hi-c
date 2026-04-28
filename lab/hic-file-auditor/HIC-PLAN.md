# HIC-PLAN.md

# .hic Format Compatibility: Root Cause Analysis and Tooling Plan

## Document Purpose

This document captures the full analysis of .hic file format incompatibilities encountered across the mariner Hi-C project, the Stripenn stripe pipeline, and related lab workflows. It serves as ground truth for building reusable tooling to diagnose and fix these issues.

**Created:** 2026-04-21
**Context:** BAP1-KO Hi-C analysis on SDSC Expanse. Primary tools: Juicer 2.0, hictk, mariner (strawr), Stripenn (cooler/mcool), Juicebox.

---

## 1. Background: How We Got Here

### 1.1 The Triggering Incident (Stripenn Pipeline, 250402 Timepoint)

The Stripenn pipeline (`stripes/stripenn/`) requires .mcool input. During Stage 0 (conversion), the 250402 merged .hic files exposed a chain of failures:

1. **hic2cool could not parse the files at all** — they are Juicer v9 format; hic2cool only handles v7/v8.
2. **hictk converted matrix data but failed on GW_SCALE** — hictk added `unable to read "GW_SCALE" weights`. Adding `--normalization-methods KR` told hictk to skip GW_SCALE and only transfer KR, but then hictk warned `KR normalization vector is missing. SKIPPING!`. The resulting mcool files had pixel data but no KR weights.
3. **juicer_tools addNorm appeared to succeed but didn't** — ran via `abc/scripts/addnorm.sb`, exited 0, printed "Finished writing norms", but hictk still could not read the KR vectors from the .hic file. Likely a v9 sub-variant incompatibility between juicer_tools 2.0.1's writer and hictk v2.1.4's reader.
4. **cooler balance was the workaround** — `cooler balance --name KR` computed ICE balancing weights directly on the mcool files and named the column "KR" so Stripenn's `--norm KR` found it. Applied via `stripes/stripenn/scripts/fix_250402_balance.sb`.

The affected files were `data/cool/250402/ctrl_merged.mcool` and `data/cool/250402/mut_merged.mcool` only. All other 14 mcool files had native KR from hictk conversion of .hic files that contained KR.

### 1.2 Broader Lab Context

This is not an isolated incident. From October 2025 meeting notes:

```
- juicer: 1.6 vs 2 vs v9
    - using 2
    - - hictk
```

The lab has encountered .hic version incompatibilities across multiple projects. Graduate mentor has confirmed encountering additional v9-related issues beyond what's documented here.

### 1.3 Lab's .hic File Generation Workflow

From the lab Google Doc (juicerpre & hiccups instructions):

> If you want to run hiccups (loops) or arrowhead (TADs), your hic files need to be in juicer v2.0 format. In order to generate this, you need to run juicerpre on your pairs output from the nextflow pipeline.
>
> To do this, you will need to pre-format your pairs files from hicpro/valid_pairs/pairix:
> - `cd /expanse/lustre/projects/csd940/ctea/HiC/juicer_scripts`
> - Vim `format_4DN.sb` and `juicer_pre.sb` and change `Sample_folder` to the directory that contains your processed nextflow folders
> - Run `sbatch format_4DN.sb [sample_name]`
> - This produces a juicer folder with `[sample_name]_short.pairs` for juicerpre
> - Run `sbatch juicer_pre.sb [sample_name]`
> - This produces a juicer 2.0 .hic file
>
> Note: if you want these files to be KR normalized and only KR-normalized (needed for some tools), then add `-k KR` as a flag. Otherwise, it will run the standard VC, VC_SQRT, GW_SCALE, SCALE, INTER_SCALE.

This means the lab's standard pipeline produces v9 .hic files via `juicer pre` from Juicer 2.0. The files contain the full set of v9 normalization types (VC, VC_SQRT, GW_SCALE, SCALE, INTER_SCALE) unless explicitly restricted with `-k KR`.

Additional note: GenOVA requires the KR weight column from juicer KR-normalized .hic files specifically, which is a separate compatibility concern.

### 1.4 Format Usage Across the Project

Most of the lab's work uses .hic format directly. The Stripenn pipeline is one of the only cases in two quarters where .mcool was required. The mariner loop pipeline uses .hic via strawr, Juicebox uses .hic, HiCCUPS uses .hic. Any solution should prioritize staying in .hic format where possible — mcool as an intermediary is acceptable but not the preferred endpoint.

---

## 2. Root Cause: .hic v9 is a Spec-Less Binary Format

### 2.1 The Core Problem

The .hic format is a binary container developed by the Aiden Lab (Baylor College of Medicine / Rice University) as part of the Juicer ecosystem. There is **no formal specification document** — the "spec" is the Juicer Java source code. When the Aiden Lab bumped from v8 to v9 (Juicer 2.0), three things changed in the binary layout:

1. **Header/footer index structure** — v9 reorganized how the file footer indexes matrix blocks and normalization vectors. The offsets and compression scheme for locating data within the file changed.

2. **Normalization vector storage** — In v8, KR/VC/VC_SQRT vectors lived at predictable offsets described by the footer. In v9, the normalization entries use a different encoding layout. This is why hic2cool fails entirely (written against v7/v8, never updated) and why hictk can read the matrix data from v9 but sometimes can't find normalization vectors written by a different tool.

3. **No sub-version signaling** — The version field just says "9", but different builds of juicer_tools (2.0.x, 2.13.x, etc.) write slightly different internal layouts. There's no way for a reader to know which v9 variant it's dealing with.

### 2.2 The Version Lineage

| Version | Era | Writer | Notes |
|---------|-----|--------|-------|
| v6 | Early Juicer | Juicer ~1.0 | Basic format |
| v7 | Juicer 1.x | Juicer 1.x | Added normalization vectors (KR, VC, VC_SQRT). Became de facto standard. |
| v8 | Juicer 1.6.x | Juicer 1.6 | Minor changes from v7. Most v7 readers handle v8 fine. |
| v9 | Juicer 2.0+ | Juicer 2.0, hictk | Major format revision. New norm types. Binary layout overhaul. |

Most tools in the Hi-C ecosystem were built against v7/v8. When v9 shipped, each tool was updated independently (or not at all) with its own interpretation of the changes.

### 2.3 Why There's No "samtools for .hic"

For BAM files, the community converged on the SAM spec (a formal document maintained by the GA4GH consortium), and samtools became the universal Swiss army knife. For .hic files, neither happened:

1. **No formal spec** — The format is defined by Juicer's Java source. When the source changes, the "spec" changes. Other tools reverse-engineer from code snapshots.

2. **No neutral governance** — The Aiden Lab controls the format. The 4DN/open2C community responded by developing cooler/mcool as an alternative built on HDF5 (which does have a spec), rather than trying to standardize .hic.

3. **hictk is the closest thing to a universal tool** — Roberto Rossini's hictk (C++ library + CLI) reads v7/v8/v9 .hic and mcool, writes both formats. But it's still catching up to the full surface area of v9 edge cases, and it can only write v9 .hic (not v8).

4. **The community is split:**
   - **4DN / open2C ecosystem** prefers .mcool (cooler format). Tools: cooler, cooltools, pairtools, hictk.
   - **Aiden Lab ecosystem** prefers .hic. Tools: Juicer, Juicebox, HiCCUPS, Arrowhead.
   - **Other tools**: HiCExplorer (prefers .h5 or .mcool), FAN-C (supports both), Mustache (supports both), Stripenn (.mcool).
   - Labs like ours need both formats depending on the tool.

---

## 3. Six Classes of .hic Incompatibility

### Class 1: Total Read Failure (Version Gate)

The first bytes of a .hic file are `HIC\0` + an int32 version number. Many tools read that int, see `9`, and immediately bail because they only have codepaths for v7/v8.

**Affected tools:**
- **hic2cool** — Dead on v9. Never updated. This is usually the first thing that breaks.
- **Older straw versions** (C++/Python/R bindings) — The original straw by Neva Durand only handled through v8. Newer forks added v9 but adoption is uneven.
- **HiCExplorer's hicConvertFormat** — Older versions can't parse v9.
- **FAN-C** — Depends on version; older releases fail on v9.

**Symptom:** Immediate crash or "invalid format" error. The tool doesn't even attempt to read data.

### Class 2: Matrix Data Readable, Normalization Vectors Not

The v9 format changed two things independently: the matrix block index structure, and the footer/normalization vector layout. Some tools updated their matrix reader for v9 but not their normalization reader, or vice versa.

**Sub-issues:**
- **v9 introduced new normalization types** (`GW_SCALE`, `INTER_SCALE`, `SCALE`) that v7/v8 readers don't know how to skip gracefully — they hit an unknown type and crash instead of ignoring it.
- **KR vector encoding changed** — Even for the "same" normalization (KR), the byte layout (compression, offset indexing) differs between v8 and v9. A reader that finds the KR entry in the footer may still misparse the actual vector bytes.
- **Resolution-specific vs genome-wide norms** — v9 stores some norms at a genome-wide level rather than per-chromosome. Tools expecting per-chromosome layout get confused.

**Symptom:** Tool reads raw contact data fine but throws "normalization vector not found" or returns corrupt normalization values.

**Our example:** hictk reading 250402 merged files — could extract matrix data at both resolutions, but GW_SCALE caused an error and KR was missing entirely.

### Class 3: Cross-Tool Writer Incompatibility (Silent Success)

This is the most dangerous class because **nothing errors out**. Tool A writes data that Tool B claims to support but parses incorrectly.

Both tools "support v9." But juicer_tools writes normalization entries using Juicer's internal v9 layout, while hictk's v9 reader expects a slightly different byte arrangement for the same logical data. The footer index says "KR is at offset X, length Y" — but the actual encoding at that offset differs between what juicer_tools wrote and what hictk expects to find.

**This class also appears when:**
- Creating .hic files with `juicer pre` from one Juicer version, then processing with `juicer_tools` from a different minor version.
- Using third-party .hic writers (e.g., pairtools or custom scripts) that implemented v9 from reading the source code at a point in time, then the Juicer codebase evolved.

**Symptom:** Write operation succeeds (exit code 0, "Finished writing norms"), but downstream reads fail or return wrong values.

**Our example:** `juicer_tools addNorm` on the 250402 merged .hic files — exited 0, printed success, but hictk could not read the KR vectors it wrote.

### Class 4: Merge/Aggregation Normalization Loss

When merging multiple .hic files (via `juicer mega`, `juicer pre` on combined pairs, or other tools), normalization vectors are **not carried over** — they must be recomputed on the merged matrix. Different tools handle this differently:

- Some compute KR/VC automatically after merge.
- Some compute only a subset (VC but not KR).
- Some compute nothing and leave it to the user.
- The 250402 merged files fell into this last category — created by a different lab member using a merging tool/version that did not compute KR.

This isn't strictly a v9 problem (can happen with any version), but v9 makes it worse because the "add norms after the fact" pathway (Class 3) is also broken across tools.

**Symptom:** File is structurally valid and contains contact data, but expected normalization vectors are absent. The file only has a subset of norms (e.g., VC, VC_SQRT, GW_SCALE, SCALE, INTER_SCALE but not KR).

**Our example:** The 250402 merged .hic files had VC/VC_SQRT/GW_SCALE/SCALE/INTER_SCALE but no KR. The 250831 merged files (created separately by a different process) did have KR.

### Class 5: Partial/Corrupt Files from Interrupted Operations

Conversions and writes can crash mid-operation, leaving a file that is **partially valid** — it passes basic existence checks and even header parsing, but operations at certain resolutions or chromosomes fail.

**Manifestations:**
- `hictk convert` crashes after writing 5kb but before writing 10kb resolution.
- `juicer pre` crashes mid-write — file has valid header, some resolutions, but missing data blocks for others.
- The file passes `-s file` checks (non-empty) and even header-level tool queries, but fails when accessing specific resolutions or chromosomes.

**Symptom:** Tool works for some operations but fails for others. Stage 0 idempotent check sees the file exists and skips reconversion. Downstream stages fail with resolution-specific errors.

**Our example:** Initial hictk conversion of 250402 merged files crashed on GW_SCALE after writing the 5kb resolution but before writing 10kb. Stage 1 later failed with `KeyError: No cooler found at ... Coolers found in ['/resolutions/5000']`.

### Class 6: Language Binding Fragmentation

The .hic reader exists independently in multiple languages:

| Binding | Language | Used By |
|---------|----------|---------|
| juicer_tools / Juicebox | Java | Juicer pipeline, Juicebox desktop, HiCCUPS, Arrowhead |
| hictk | C++ | hictk CLI, growing ecosystem |
| hic-straw / hicstraw | Python | Various Python Hi-C tools |
| strawr | R (Rcpp) | mariner, our loop pipeline |
| Juicebox.js | JavaScript | Juicebox web viewer |

Each is an independent binary parser. When v9 shipped, each binding was updated on its own timeline with its own interpretation of the format changes. The same .hic file can work in Java juicer_tools but fail in strawr, or work in hictk but fail in Python hic-straw.

**Why this matters for our lab:** We use multiple language ecosystems simultaneously — R/strawr for mariner loops, Python/cooler for Stripenn stripes, Java for Juicebox visualization, C++ via hictk for format conversion. A file that's valid in one ecosystem may not be in another.

---

## 4. The .hic ↔ .mcool Relationship

### 4.1 Why cooler balance Worked as a Workaround

The mcool format (HDF5-based, from the cooler ecosystem) **does** have a well-defined spec. Once matrix data is converted to mcool via hictk (which reads v9 matrices fine), you're in a clean ecosystem. `cooler balance` computes ICE weights and stores them as a named column ("KR") in the HDF5 bins table — a format that every tool in the cooler ecosystem (and Stripenn) reads identically. The name "KR" is just a label; the algorithm is ICE, but ICE and KR are functionally equivalent iterative matrix balancing algorithms.

### 4.2 The Reliable Conversion Path (Current Known-Good)

```
.hic (any version) → hictk convert → .mcool → cooler balance --name KR → normalized mcool
```

This path works because:
- hictk reads v7/v8/v9 matrix data reliably
- The mcool format is spec'd and consistent
- cooler balance operates entirely within the mcool ecosystem
- No cross-tool normalization writing is involved

### 4.3 Limitations of mcool as a Solution

- Most of the lab's tools expect .hic: mariner (strawr), Juicebox, HiCCUPS, Arrowhead, GenOVA.
- The Stripenn pipeline is one of the few cases in two quarters where .mcool was needed.
- Round-tripping .hic → .mcool → .hic may lose metadata.
- The goal should be to fix/validate .hic files directly where possible, with mcool as a fallback.

---

## 5. Solution Plan: Two-Phase Approach

### Phase 1: Auditor (`hic-audit`) — COMPLETE (2026-04-22)

Build a diagnostic tool that inspects a .hic file and reports everything needed to understand its compatibility state. This is immediately useful and doesn't require solving the hard writer problems.

**Implemented:** `hic/hic_audit.py` — single-file Python script, zero required dependencies. See `hic/README.md` for usage guide.

**What it reports:**
- File version (v7/v8/v9)
- All resolutions present
- All normalization vectors present, per resolution
- Chromosome list and naming convention (chr1 vs 1)
- Whether the file can be read by common tools/bindings (hictk, straw, hic2cool)
- File integrity: all declared resolutions actually contain data
- Specific incompatibilities detected and recommendations

**Design decisions (resolved):**
- Language: Python 3.8+ (binary parsing via struct, hicstraw/hictk optional at runtime)
- Scope: Tiered — fast header-only by default, `--deep` for functional data validation
- Output: Human-readable by default, `--json` for machine output

**Validated against real files:** Correctly identifies missing KR in 250402 merged .hic files and present KR in 250831 files. Also discovered that v9 headers differ from v7/v8 in two previously undocumented ways — see `hic/header.md` for the reverse-engineered format reference.

**Additional deliverables:**
- `hic/header.md` — Reverse-engineered .hic v9 binary header format with annotated hex dumps from project files
- `hic/README.md` — Usage guide, CLI reference, common scenarios

### Phase 2: Fixers (per-class)

Build targeted fix scripts for each class of issue, informed by what the auditor reveals about the actual files encountered. These should be modular — each class gets its own fixer, composed as needed.

**Class 1 (version gate):** Convert to target format via hictk. This is essentially solved — hictk is the bridge.

**Class 2 (missing norms):** Compute missing normalization vectors. Two sub-paths:
- Stay in .hic: `juicer_tools addNorm` (works when writer and reader are same tool)
- Via mcool: hictk convert → cooler balance → optionally convert back

**Class 3 (cross-tool writer):** Re-create the file from scratch using a single tool's writer. Dump contacts → write with chosen tool. Nuclear option but guaranteed consistent.

**Class 4 (merge norm loss):** Detect post-merge and compute norms. Essentially a specialized case of Class 2 but triggered by merge detection.

**Class 5 (partial/corrupt):** Detect incomplete resolutions → reconvert the missing ones or full reconvert.

**Class 6 (binding fragmentation):** The auditor itself addresses this by testing readability across bindings. The fix is "use the conversion path that targets your specific downstream tool."

**Target version question:**
- hictk writes v9 .hic — this is probably the best "canonical" v9 writer since it's maintained by the open2C community and tested against the broadest tool set.
- No maintained tool currently writes v8 from arbitrary input (v9→v8 downgrade is not supported).
- For maximum compatibility: produce both .hic (v9 via hictk) and .mcool, let the user choose.

### Phase 2 Priority Order

1. Class 2 (missing norms) — Most common issue, directly caused the 250402 incident
2. Class 5 (partial/corrupt) — Caused the incomplete mcool issue, silent and dangerous
3. Class 4 (merge norm loss) — Variant of Class 2 but important for merged samples
4. Class 1 (version gate) — Already solved by hictk, just needs wrapping
5. Class 3 (cross-tool writer) — Hardest to solve, requires full re-creation
6. Class 6 (binding fragmentation) — Addressed by auditor + documentation

---

## 6. Tool Ecosystem Reference

### 6.1 Tools and Their .hic Version Support

| Tool | Reads v7/v8 | Reads v9 | Writes .hic | Writes .mcool | Notes |
|------|-------------|----------|-------------|---------------|-------|
| hic2cool | Yes | **No** | No | Yes | Dead on v9. Unmaintained. |
| hictk | Yes | Yes (matrix); partial (norms) | Yes (v9 only) | Yes | Best universal reader. Can't write v8. |
| juicer_tools | Yes | Yes | Yes (v9) | No | Reference implementation but cross-tool norm issues |
| strawr (R) | Yes | Partial | No | No | Used by mariner. v9 support depends on version. |
| hic-straw (Python) | Yes | Partial | No | No | Multiple forks with different v9 support levels. |
| cooler | No | No | No | Yes | mcool ecosystem only. Well-spec'd. |
| Juicebox desktop | Yes | Yes | No | No | Java reader, most complete v9 support. |
| GenOVA | Yes | ? | No | No | Requires KR weight column specifically. |

### 6.2 Normalization Types

| Type | Versions | Algorithm | Notes |
|------|----------|-----------|-------|
| KR | v7+ | Knight-Ruiz matrix balancing | The standard. Required by most tools. |
| VC | v7+ | Vanilla Coverage | Simple read-depth normalization. |
| VC_SQRT | v7+ | Square root of VC | Less aggressive than VC. |
| GW_SCALE | v9 | Genome-wide scaling | New in v9. Causes hictk errors. |
| SCALE | v9 | General scaling | New in v9. |
| INTER_SCALE | v9 | Inter-chromosomal scaling | New in v9. |
| ICE | mcool | Iterative Correction and Eigenvector decomposition | cooler balance default. Functionally equivalent to KR. |

### 6.3 Our Lab's Tool Chain

| Pipeline | Format | Reader | Downstream Tools |
|----------|--------|--------|------------------|
| Mariner loops | .hic | strawr (R) | edgeR, Juicebox |
| Stripenn stripes | .mcool | cooler (Python) | edgeR, Juicebox |
| HiCCUPS (loop calling) | .hic | juicer_tools (Java) | — |
| Arrowhead (TADs) | .hic | juicer_tools (Java) | — |
| Juicebox visualization | .hic | Juicebox (Java) | — |
| GenOVA | .hic | GenOVA reader | — |
| HOMER (TADs/compartments) | .hic | HOMER reader | — |

---

## 7. Files Related to This Issue

### One-Off Fix Scripts (Stripenn pipeline)

- `stripes/stripenn/scripts/fix_250402_merged.sb` — Reconverted incomplete mcool files with `--normalization-methods KR` to skip GW_SCALE.
- `stripes/stripenn/scripts/fix_250402_balance.sb` — Computed ICE weights via `cooler balance --name KR` on the two affected mcool files.
- Both marked "delete after use" in their comments.

### Lab Infrastructure

- `/expanse/lustre/projects/csd940/ctea/HiC/juicer_scripts/format_4DN.sb` — Formats pairs files for juicerpre.
- `/expanse/lustre/projects/csd940/ctea/HiC/juicer_scripts/juicer_pre.sb` — Runs juicer pre to create v2.0 .hic files.
- `abc/scripts/addnorm.sb` — Attempted (and failed) juicer_tools addNorm on 250402 merged files.

### Documentation

- `stripes/stripenn/README.md` Section 5.2 — Full history of the 250402 KR normalization issue.
- `stripes/stripenn/CLAUDE.md` — References the KR issue and the two conda environment requirement (hictk env for Stage 0).

---

## 8. Open Questions

1. **What .hic version does strawr (mariner's reader) actually support?** — Need to test v9 .hic files directly with strawr to determine if the main loop pipeline is affected.

2. **What other v9 issues has the graduate mentor encountered?** — Known to exist but not documented. Should be captured when encountered.

3. **Can hictk round-trip .hic → .mcool → .hic without metadata loss?** — Relevant if the mcool intermediary path is used for .hic repair.

4. **Is there a way to write v8 .hic from v9 input?** — Would solve Class 1 universally. Currently no maintained tool does this. `juicer pre` from Juicer 1.6 writes v8 but is old and unmaintained.

5. **Which v9 sub-variants exist in the wild?** — Different juicer_tools builds (2.0.x, 2.13.x, etc.) may produce different internal layouts. Exhaustive testing needed.

6. **Does `-k KR` in juicer pre avoid the GW_SCALE/SCALE/INTER_SCALE issues entirely?** — If so, the lab could adopt KR-only normalization as a standard practice to avoid Class 2 issues with v9-specific norm types.

---

## 9. Decision Log

| Date | Decision | Rationale |
|------|----------|-----------|
| 2025-10 | Use Juicer 2.0 (v9) as standard | Required for HiCCUPS/Arrowhead on Expanse |
| 2025-10 | Use hictk for .hic → .mcool conversion | hic2cool can't parse v9 |
| 2026-04 | Use cooler balance --name KR for missing KR | juicer_tools addNorm wrote norms that hictk couldn't read (Class 3) |
| 2026-04-21 | Build auditor first, fixers second | Biggest gap is diagnosis, not conversion — currently debugging by trial and error |
| 2026-04-21 | Stay in .hic where possible | Lab primarily uses .hic; mcool intermediary is acceptable but not preferred endpoint |

---

**Last Updated:** 2026-04-21
**Project:** BAP1-KO Differential Chromatin Analysis (mariner Hi-C)
**Related:** `stripes/stripenn/README.md` Section 5.2, `stripes/stripenn/CLAUDE.md`
