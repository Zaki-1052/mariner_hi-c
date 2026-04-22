# hic-audit: .hic File Compatibility Auditor

A diagnostic tool that inspects `.hic` files and reports their version, normalization vectors, structural integrity, and compatibility with common Hi-C tools. Designed to replace trial-and-error debugging of .hic format issues.

## Why This Exists

The .hic format (Aiden Lab / Juicer) has no formal specification. Version 9 changed the binary layout in ways that break many tools — silently in some cases. When a .hic file fails downstream, the error messages are cryptic and the root cause is rarely obvious. This tool reports everything you need to know about a file's compatibility in one command.

See `HIC-PLAN.md` for the full root cause analysis (six classes of .hic v9 incompatibility) and `header.md` for the reverse-engineered binary format reference.

## Quick Start

```bash
# Basic audit (binary header only, zero dependencies)
python hic_audit.py sample.hic

# Deep audit with functional normalization testing (needs hicstraw)
conda activate hic  # or any env with hicstraw
python hic_audit.py --deep sample.hic

# Audit specific resolutions only (faster for large files)
python hic_audit.py --deep --resolution 5000 --resolution 10000 sample.hic

# Batch audit multiple files
python hic_audit.py *.hic

# Machine-readable output
python hic_audit.py --json sample.hic | python -m json.tool

# Only show files with problems
python hic_audit.py --quiet *.hic
```

## Dependencies

**Required:** Python 3.8+ (stdlib only — `struct`, `argparse`, `json`, `subprocess`)

**Optional (detected at runtime):**

| Dependency | What it enables | How to get it |
|---|---|---|
| `hictk` CLI | Normalization listing, structural validation | `conda activate hictk` (or install via conda/pip) |
| `hicstraw` Python package | Functional data readability testing (`--deep`) | `conda activate hic` (or `pip install hicstraw`) |

The tool works with zero external dependencies. When optional tools are available, it uses them for deeper checks and reports what's available in the `[TOOLS]` section.

## Understanding the Output

### Example: file with missing KR normalization

```
=== .hic File Audit: late-ctrl_merged.hic ===

  [HEADER]  Version: 9 | Genome: /.../mm10.filtered.chrom.sizes | Size: 6.3 GB
  [HEADER]  Chromosomes: 21 (chr1..chr19) | Naming: UCSC (chr-prefix)
  [HEADER]  Resolutions (BP): 1000, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000, 2500000
  [HEADER]  Metadata: 2 attribute(s)
            software: Juicer Tools Version 2.20.00
            hicFileScalingFactor: 1.0

  [TOOLS]   hictk: not found
  [TOOLS]   hicstraw: 1.3.1

  [DATA]    Functional readability (chr1 spot-check):
    5000bp:  NONE:PASS  KR:FAIL  VC:PASS  VC_SQRT:PASS  SCALE:PASS  GW_SCALE:FAIL  INTER_SCALE:FAIL
    10000bp: NONE:PASS  KR:FAIL  VC:PASS  VC_SQRT:PASS  SCALE:PASS  GW_SCALE:FAIL  INTER_SCALE:FAIL

  [COMPAT]
    PASS  Juicebox / juicer_tools — v9 native
    FAIL  hic2cool — requires v7/v8
    PASS  hictk — v7-v9 supported
    WARN  strawr (R/mariner) — v9 support is version-dependent
    PASS  hicstraw (Python) — v8/v9 supported
    N/A   cooler ecosystem — requires mcool conversion (hictk convert)

  [ACTION]
    * KR normalization missing at: 5000bp, 10000bp. Fix: hictk balance ice <file> or convert to mcool + cooler balance --name KR.
    * File is v9. Tools that only support v7/v8 (hic2cool, older straw) will fail. Use hictk for conversions.

  [SUMMARY] 1 incompatible, 1 warning(s), 3 compatible, 2 recommendation(s)
```

### Output sections

| Section | Source | What it tells you |
|---|---|---|
| `[HEADER]` | Binary header parse (always available) | Version, genome, chromosomes, resolutions, metadata. Works on any .hic file regardless of installed tools. |
| `[TOOLS]` | Runtime detection | Which optional tools are available. Determines what deeper checks can run. |
| `[NORMS]` | hictk CLI (if available) | Lists normalization vectors declared in the file, per resolution. Shown when hictk is on PATH. |
| `[VALID]` | hictk CLI (if available) | Structural validation result — whether the file is internally consistent. |
| `[DATA]` | hicstraw (if available + `--deep`) | Functional test: actually tries to read data with each normalization type. Catches silent failures where norms are declared but unreadable. |
| `[COMPAT]` | Synthesis of above | Per-tool compatibility verdict based on version and norms. |
| `[ACTION]` | Synthesis of above | Actionable recommendations for fixing detected issues. |
| `[SUMMARY]` | Counts | Quick overview of how many tools are compatible/incompatible. |

### Status meanings

| Status | Color | Meaning |
|---|---|---|
| PASS | Green | Tool can read this file / norm is accessible |
| WARN | Yellow | May work depending on tool version or configuration |
| FAIL | Red | Will not work / norm is missing or unreadable |
| N/A | Plain | Not applicable (e.g., cooler can't read .hic directly) |

## CLI Reference

```
usage: hic_audit.py [-h] [--deep] [--json] [--resolution RES] [--quiet] file [file ...]

Audit .hic files for version, normalization, and compatibility issues.

positional arguments:
  file                  .hic file(s) to audit

options:
  -h, --help            show this help message and exit
  --deep                Functional readability tests (requires hicstraw)
  --json                Machine-parseable JSON output
  --resolution RES      Only check specific resolution(s). Can be repeated.
  --quiet               Only print files with warnings or errors
```

### Flags in detail

**`--deep`**: Tests actual data readability by attempting to fetch contacts from chr1 at each resolution with each known normalization type (NONE, KR, VC, VC_SQRT, SCALE, GW_SCALE, INTER_SCALE). Requires `hicstraw` to be importable. Takes ~5-30 seconds per file depending on size and number of resolutions.

This flag is critical for catching **Class 3 issues** (cross-tool writer incompatibility) where normalization vectors are declared in the file but unreadable by hicstraw. Without `--deep`, the tool can only report what's in the header.

**`--resolution RES`**: Limits which resolutions are checked (for both hictk norms and `--deep` testing). Can be specified multiple times. Useful for speeding up audits when you only care about specific resolutions:

```bash
python hic_audit.py --deep --resolution 5000 --resolution 10000 sample.hic
```

**`--json`**: Outputs a JSON object (single file) or JSON array (multiple files). Useful for piping to other tools or scripting:

```bash
# Check if KR is available at 5kb via jq
python hic_audit.py --json --deep --resolution 5000 sample.hic \
  | jq '.deep_readability.norm_tests["5000"].KR'
```

**`--quiet`**: Only prints output for files that have warnings, errors, or recommendations. Useful for batch scanning a directory of .hic files to find the problematic ones.

## How It Works

### Step 1: Binary header parse (always runs)

Reads the first ~500 bytes of the file using Python's `struct` module. Extracts:
- Magic bytes and version (first 8 bytes)
- Master index offset
- Genome ID string
- v9-specific nvi fields (nviPosition, nviLength)
- Metadata key-value pairs
- Chromosome names and lengths (int32 for v7/v8, int64 for v9)
- BP and fragment resolutions

This is a direct port of the header parser from [hic2cool](https://github.com/4dn-dcic/hic2cool) (via mustache's `diff_mustache.py`), extended to handle v9's additional fields. See `header.md` for the full byte-level layout.

### Step 2: hictk integration (if available)

If `hictk` is on PATH, the tool runs:
- `hictk metadata <file> -f json` — full metadata
- `hictk dump -t normalizations <file> --resolution <res>` — lists available norms per resolution
- `hictk validate <file> -f json` — structural validation

Each command runs in a subprocess with a 30-second timeout. If hictk is not available, these sections are skipped with a note in the output.

### Step 3: hicstraw functional testing (--deep only)

If `hicstraw` is importable and `--deep` is set:
1. Opens the file with `hicstraw.HiCFile(path)`
2. Compares hicstraw's reported resolutions/chromosomes against the header
3. For each resolution: tries `getMatrixZoomData()` with each of 7 known normalization types
4. **Captures C-level stderr** during each call — hicstraw prints warnings like `"File did not contain KR normalization vectors"` to stderr but does NOT raise Python exceptions. Without stderr capture, missing norms appear to succeed (returning raw unnormalized data).

### Step 4: Compatibility assessment

Synthesizes all findings into:
- Per-tool verdicts (based on version and norm availability)
- Actionable recommendations (what to fix and how)

The compatibility logic:

| Tool | PASS when | FAIL when | WARN when |
|---|---|---|---|
| Juicebox / juicer_tools | version >= 7 | — | version < 7 |
| hic2cool | version <= 8 | version 9 | — |
| hictk | always | — | — |
| strawr (R/mariner) | version <= 8 | — | version 9 (support varies) |
| hicstraw (Python) | version >= 8 | — | version 7 |
| cooler ecosystem | — | ��� | always N/A (needs mcool) |

## Common Scenarios

### "My tool says normalization not found"

Run `--deep` to see which norms are actually readable:

```bash
python hic_audit.py --deep --resolution 5000 problem_file.hic
```

Look at the `[DATA]` section. If KR shows FAIL, the file was either:
- Created without `-k KR` in `juicer pre` (has VC/VC_SQRT/SCALE instead)
- Merged without recomputing norms
- Written by a tool whose norm format another tool can't read (Class 3)

### "hic2cool / hic2cool can't open my file"

Almost certainly a v9 file. The `[HEADER]` section confirms the version. Use `hictk convert` instead of hic2cool — hictk handles v7-v9.

### "hictk convert fails on GW_SCALE"

The file has GW_SCALE normalization vectors that hictk can't read. Pass `--normalization-methods KR` to hictk to skip problematic norms and only transfer KR.

### "I need to audit all my .hic files at once"

```bash
# Find all problems in a directory tree
find /path/to/hic/files -name "*.hic" -exec python hic_audit.py --quiet {} +

# JSON summary of everything
find /path/to/hic/files -name "*.hic" | xargs python hic_audit.py --json > audit_results.json
```

### "Which of my files are missing KR?"

```bash
# Deep scan at 5kb resolution, JSON output, filter for KR failures
for f in *.hic; do
  kr=$(python hic_audit.py --json --deep --resolution 5000 "$f" 2>/dev/null \
    | python -c "import sys,json; d=json.load(sys.stdin); print(d.get('deep_readability',{}).get('norm_tests',{}).get('5000',{}).get('KR','N/A'))" 2>/dev/null)
  echo "$f: KR=$kr"
done
```

## Files in This Directory

| File | Purpose |
|---|---|
| `hic_audit.py` | The auditor script |
| `HIC-PLAN.md` | Root cause analysis of .hic v9 incompatibilities, six classes, and the two-phase tooling plan |
| `header.md` | Reverse-engineered .hic binary header format with worked examples from project files |
| `README.md` | This file |

## Related

- `stripes/stripenn/README.md` Section 5.2 — The 250402 KR normalization incident that motivated this tool
- `stripes/stripenn/scripts/fix_250402_balance.sb` — The workaround script (cooler balance) for the missing KR issue
