Great question — this is a genuinely confusing ecosystem issue and understanding the root cause makes future workarounds much less mysterious.

## The Core Problem: .hic v9 is a Spec-Less Binary Format

The .hic format is a binary container developed by the Aiden Lab as part of Juicer. There is no formal specification document — the "spec" is the Juicer Java source code. When they bumped from v8 to v9 (Juicer 2.0), three things changed in the binary layout:

1. **Header/footer index structure** — v9 reorganized how the file footer indexes matrix blocks and normalization vectors. The offsets and compression scheme for locating data within the file changed.

2. **Normalization vector storage** — In v8, KR/VC/VC_SQRT vectors lived at predictable offsets described by the footer. In v9, the normalization entries use a different encoding layout. This is why `hic2cool` fails entirely (it was written against v7/v8 and never updated) and why `hictk` can read the *matrix data* from v9 but sometimes can't find normalization vectors written by a different tool.

3. **No sub-version signaling** — The version field just says "9", but different builds of juicer_tools (2.0.x, 2.13.x, etc.) write slightly different internal layouts. There's no way for a reader to know *which* v9 variant it's dealing with.

## Why Your Specific Chain Failed

Your 250402 merged files hit three separate manifestations of this:

**hic2cool couldn't parse them at all** — hic2cool only understands v7/v8 headers. It reads the version byte, sees "9", and has no codepath for it.

**hictk could read matrix data but not GW_SCALE** — hictk v2.x added v9 matrix support but doesn't handle every normalization variant. The `GW_SCALE` normalization type (genome-wide scaling) uses a format hictk doesn't expect, so it errors unless you explicitly skip it with `--normalization-methods KR`.

**juicer_tools addNorm "succeeded" but hictk couldn't read the result** — This is the most insidious one. `juicer_tools addNorm` wrote KR vectors into the v9 file using Juicer's own internal v9 layout. It exited 0, printed "Finished writing norms." But hictk's v9 normalization reader expects vectors at offsets/encodings that don't match what juicer_tools actually wrote. Both tools "support v9" but they disagree on the binary details. Since there's no formal spec, neither is technically wrong.

## Why cooler balance Worked

The mcool format (HDF5-based, from the cooler ecosystem) *does* have a well-defined spec. Once you convert the matrix data to mcool via hictk (which reads v9 matrices fine), you're in a clean ecosystem. `cooler balance` computes ICE weights and stores them as a named column ("KR") in the HDF5 bins table — a format that every tool in the cooler ecosystem (and stripenn) reads identically. The name "KR" is just a label; the algorithm is ICE, but functionally equivalent for downstream normalization.

## Could You Build Something Reusable?

Yes, and the approach is essentially what you already did, generalized. The reliable path through the v9 minefield is:

```
.hic (any version) → hictk convert → .mcool → cooler balance --name KR → normalized mcool
```

A reusable script would:
1. Detect whether the .hic is v8 or v9 (first 8 bytes: magic "HIC\0" + int32 version)
2. Attempt hictk convert (works for both v8 and v9 matrix data)
3. Check whether KR normalization survived the conversion (`cooler dump -t bins` and check for the KR column)
4. If KR is missing, run `cooler balance --name KR` to compute it
5. Verify the result

Want me to write that as a general-purpose script? It would replace the need for one-off fixes like `fix_250402_balance.sb` and `fix_250402_merged.sb` — any .hic file regardless of version or normalization state goes in, a fully-normalized mcool comes out.


---

You're right, I was scoping the problem way too narrowly. Let me think about this properly.

## The Full Taxonomy of .hic v9 Breakage

The KR normalization issue is just one of at least **six distinct classes** of incompatibility in the .hic ecosystem. They all stem from the same root: the .hic format has no formal specification — the spec *is* the Juicer Java source code, and when Aiden Lab shipped Juicer 2.0 with format v9, they changed the binary layout in ways that no other tool was consulted on.

### Class 1: Total Read Failure (version gate)

The first bytes of a .hic file are `HIC\0` + an int32 version number. Many tools read that int, see `9`, and immediately bail because they only have codepaths for v7/v8.

- **hic2cool** — Dead on v9. Never updated. This was probably the first thing that broke for your lab.
- **Older straw versions** (the C++/Python/R bindings) — The original straw by Neva Durand only handled through v8. Newer forks added v9 but adoption is uneven.
- **HiCExplorer's hicConvertFormat** — Older versions can't parse v9.
- **FAN-C** — Depends on which version; older releases fail on v9.

This is the most common class. The tool doesn't even try to read the data — it stops at the header.

### Class 2: Matrix Data Readable, Normalization Vectors Not

This is what hit you with hictk. The v9 format changed **two things independently**: the matrix block index structure, and the footer/normalization vector layout. Some tools updated their matrix reader for v9 but not their normalization reader, or vice versa.

Specific sub-issues:
- **v9 introduced new normalization types** (`GW_SCALE`, `INTER_SCALE`, `SCALE`) that v7/v8 readers don't know how to skip gracefully — they hit an unknown type and crash instead of ignoring it.
- **KR vector encoding changed** — Even for the "same" normalization, the byte layout (compression, offset indexing) differs between v8 and v9. A reader that finds the KR entry in the footer may still misparse the actual vector bytes.
- **Resolution-specific vs genome-wide norms** — v9 stores some norms at a genome-wide level rather than per-chromosome. Tools expecting per-chromosome layout get confused.

### Class 3: Cross-Tool Writer Incompatibility (the "silent success" problem)

This is the most dangerous class because **nothing errors out**. Tool A writes data that Tool B claims to support but parses incorrectly.

Your juicer_tools addNorm → hictk failure is the canonical example. Both tools "support v9." But juicer_tools 2.0.x writes normalization entries using Juicer's internal v9 layout, while hictk's v9 reader expects a slightly different byte arrangement for the same logical data. The footer index says "KR is at offset X, length Y" — but the actual encoding at that offset differs between what juicer_tools wrote and what hictk expects to find.

This class also shows up when:
- Creating .hic files with `juicer pre` from one version, then processing with `juicer_tools` from a different minor version
- Using third-party .hic writers (e.g., pairtools or custom scripts) that implemented v9 from reading the source code at a point in time, then the Juicer codebase evolved

### Class 4: Merge/Aggregation Normalization Loss

When you merge multiple .hic files (via `juicer mega`, `juicer pre` on combined pairs, or other merge tools), normalization vectors are **not carried over** — they have to be recomputed on the merged matrix. Different tools handle this differently:

- Some compute KR/VC automatically after merge
- Some only compute a subset (VC but not KR)
- Some compute nothing and leave it to the user
- The 250402 merged files fell into this last category

This isn't strictly a v9 problem — it can happen with any version — but v9 makes it worse because the "add norms after the fact" pathway (Class 3) is also broken.

### Class 5: Partial/Corrupt Files from Interrupted Operations

This hit you with the incomplete mcool conversions where 5kb was written but 10kb wasn't. But it happens in .hic land too:

- `juicer pre` crashes mid-write → file has valid header, some resolutions, but missing data blocks for others
- The file passes basic existence checks (`-s file`) and even basic header parsing, but operations at certain resolutions fail
- hictk convert crashes after writing some resolutions → same problem downstream

The pernicious part: these files are **partially valid**. They don't fail fast — they work for some operations and silently break for others depending on which resolution or chromosome you access.

### Class 6: Language Binding Fragmentation

The .hic reader exists independently in:
- **Java** (Juicer/juicer_tools — the "reference" implementation)
- **C++** (straw C++, hictk)
- **Python** (hic-straw, hicstraw, cooler's hic reader)
- **R** (strawr — what mariner uses)

Each is an independent binary parser. When v9 shipped, each binding was updated on its own timeline with its own interpretation of the format changes. So the same .hic file can work in Java juicer_tools but fail in strawr, or work in hictk but fail in Python hic-straw.

Your lab hits this because you use **multiple language ecosystems simultaneously**: R/strawr for mariner, Python for stripenn, Java for Juicebox, C++ via hictk for conversion. A file that's valid in one ecosystem may not be in another.

---

## Why There's No "samtools for .hic"

For BAM files, the community converged on the SAM spec (a formal document maintained by the GA4GH consortium), and `samtools` became the universal Swiss army knife. For .hic files, neither of those things happened:

1. **No formal spec** — The format is defined by Juicer's Java source. When the source changes, the "spec" changes. Other tools reverse-engineer from code snapshots.

2. **No neutral governance** — The Aiden Lab controls the format. The 4DN/open2C community responded by developing cooler/mcool as an alternative built on HDF5 (which *does* have a spec), rather than trying to standardize .hic.

3. **hictk is the closest thing to a universal tool** — Roberto Rossini's hictk (C++ library + CLI) reads v7/v8/v9 .hic and mcool, writes both formats. But it's still catching up to the full surface area of v9 edge cases, and it can only write v9 .hic (not v8).

4. **The community is split** — 4DN/open2C prefers mcool. Aiden Lab ecosystem prefers .hic. Neither is going away. Labs like yours need both.

## What a Reusable Tool Would Actually Need to Do

Given all six classes, the tool needs to handle a pipeline like:

```
INPUT (.hic, any version/state)
  │
  ├─ AUDIT: version, resolutions, norms present, per-binding readability
  │
  ├─ FIX PATH A (.hic → .hic):
  │    Extract contacts → juicer pre (target version) → addNorm
  │    Pro: stays in .hic, JuiceBox-ready
  │    Con: juicer pre is Java/heavy, version targeting is tricky
  │
  ├─ FIX PATH B (.hic → .mcool → .hic round-trip):
  │    hictk convert → cooler balance → hictk convert back
  │    Pro: normalization guaranteed to work
  │    Con: round-trip may lose metadata, two conversions
  │
  ├─ FIX PATH C (.hic → .mcool, keep both):
  │    hictk convert → cooler balance
  │    Pro: simplest, most reliable
  │    Con: not .hic output (fine for cooler tools, not for JuiceBox)
  │
  └─ VALIDATE: confirm output readable by target tool set
```

The hard question is which version to target for .hic output. The honest answer:

- **v8** maximizes compatibility with existing tools but `juicer pre` from Juicer 1.6 is old and unmaintained
- **v9** is the forward-looking choice but currently has the interop problems above
- **hictk's v9** is probably the most "standard" v9 writer at this point, since the open2C community is testing against it

## The Honest Feasibility Assessment

**What's buildable now with existing tools:**
- Audit/diagnosis (read header bytes + check norms + test multi-tool readability)
- .hic → .mcool with normalization guaranteed (hictk + cooler balance)
- .hic → .hic re-creation via contact dump + juicer pre (picks a consistent writer)
- Validation layer (try reading with straw/hictk/cooler and report)

**What's genuinely hard:**
- v9→v8 downgrade (no maintained tool writes v8 from arbitrary input)
- In-place norm addition that works cross-tool (the Class 3 problem — would need a new writer)
- Handling every v9 sub-variant (would need exhaustive testing against files from different juicer_tools builds)

**What would actually move the needle for the community:**
The biggest gap isn't the conversion — it's the **audit**. Right now when a .hic file fails, you get a cryptic error and have to debug by trial and error (which is exactly what you went through). A tool that reports "this is v9, has KR at 5kb but not 10kb, GW_SCALE present but non-standard, readable by hictk and juicer_tools but NOT by strawr or hic2cool" would save enormous time.

---

So I'd frame the project in two phases: (1) build the auditor, because that's immediately useful and doesn't require solving the hard writer problems; (2) build the converter/fixer informed by what the auditor reveals about your actual files. Does that framing match what you're thinking, or are you more focused on one particular class?