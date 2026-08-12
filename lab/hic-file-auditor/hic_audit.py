# hic/hic_audit.py
"""
.hic file compatibility auditor.

Reports version, normalizations, integrity, and tool compatibility for .hic
files. Designed around the six classes of v9 incompatibility documented in
hic/HIC-PLAN.md.

Zero required dependencies beyond Python 3.8+ stdlib. Optionally uses:
  - hictk CLI: normalization listing, structural validation
  - hicstraw (Python): functional data readability testing (--deep)

Usage:
  python hic_audit.py sample.hic
  python hic_audit.py --deep sample.hic
  python hic_audit.py --json *.hic
"""
import argparse
import importlib.util
import json
import os
import shutil
import struct
import subprocess
import sys


# ---------------------------------------------------------------------------
# ANSI color helpers (disabled when not a TTY or on Windows without VT)
# ---------------------------------------------------------------------------

_USE_COLOR = sys.stdout.isatty() and os.name != "nt"


def _green(s):
    return f"\033[32m{s}\033[0m" if _USE_COLOR else s


def _yellow(s):
    return f"\033[33m{s}\033[0m" if _USE_COLOR else s


def _red(s):
    return f"\033[31m{s}\033[0m" if _USE_COLOR else s


def _bold(s):
    return f"\033[1m{s}\033[0m" if _USE_COLOR else s


def _status(level):
    if level == "PASS":
        return _green("PASS")
    if level == "WARN":
        return _yellow("WARN")
    if level == "FAIL":
        return _red("FAIL")
    return level


# ---------------------------------------------------------------------------
# Step 1: Binary header parser
# Ported from hic2cool via mustache (diff_mustache.py:186-249).
# Header format is identical across .hic v7/v8/v9.
# ---------------------------------------------------------------------------

def _readcstr(f):
    buf = b""
    while True:
        b = f.read(1)
        if b is None or b == b"\0":
            return buf.decode("utf-8", errors="replace")
        if b == b"":
            raise EOFError("Unexpected EOF reading null-terminated string")
        buf += b


def parse_hic_header(path):
    file_size = os.path.getsize(path)
    if file_size < 8:
        return {"error": "File too small to be a valid .hic file", "file_size": file_size}

    result = {"file_path": path, "file_size": file_size}

    with open(path, "rb") as f:
        magic = struct.unpack("<3s", f.read(3))[0]
        f.read(1)  # null terminator
        if magic != b"HIC":
            return {**result, "error": f"Bad magic bytes: {magic!r} (expected b'HIC')"}

        version = struct.unpack("<i", f.read(4))[0]
        result["version"] = version
        master_index = struct.unpack("<q", f.read(8))[0]
        result["master_index"] = master_index

        if master_index > file_size:
            result["truncated"] = True

        # Genome ID (null-terminated string)
        result["genome"] = _readcstr(f)

        # v9 adds nviPosition + nviLength (2x int64) between genome and attributes
        if version >= 9:
            result["nvi_position"] = struct.unpack("<q", f.read(8))[0]
            result["nvi_length"] = struct.unpack("<q", f.read(8))[0]

        # Metadata key-value pairs
        n_attrs = struct.unpack("<i", f.read(4))[0]
        metadata = {}
        for _ in range(n_attrs):
            key = _readcstr(f)
            value = _readcstr(f)
            metadata[key] = value
        result["metadata"] = metadata

        # Chromosomes — v9 uses int64 for lengths, v7/v8 use int32
        chr_size_fmt = "<q" if version >= 9 else "<i"
        chr_size_bytes = 8 if version >= 9 else 4
        n_chrs = struct.unpack("<i", f.read(4))[0]
        chromosomes = []
        for _ in range(n_chrs):
            name = _readcstr(f)
            length = struct.unpack(chr_size_fmt, f.read(chr_size_bytes))[0]
            if name and length:
                chromosomes.append({"name": name, "length": length})
        result["chromosomes"] = chromosomes

        # BP resolutions
        n_bp_res = struct.unpack("<i", f.read(4))[0]
        bp_resolutions = []
        for _ in range(n_bp_res):
            bp_resolutions.append(struct.unpack("<i", f.read(4))[0])
        result["bp_resolutions"] = sorted(bp_resolutions)

        # Fragment resolutions (same format, immediately after BP resolutions)
        try:
            n_frag_res = struct.unpack("<i", f.read(4))[0]
            frag_resolutions = []
            for _ in range(n_frag_res):
                frag_resolutions.append(struct.unpack("<i", f.read(4))[0])
            result["frag_resolutions"] = sorted(frag_resolutions)
        except struct.error:
            result["frag_resolutions"] = []

    # Derived fields
    chr_names = [c["name"] for c in chromosomes if c["name"] != "All"]
    if chr_names:
        result["chr_naming"] = "UCSC (chr-prefix)" if chr_names[0].startswith("chr") else "Ensembl (no prefix)"
        result["chr_count"] = len(chr_names)
        result["chr_range"] = f"{chr_names[0]}..{chr_names[-1]}"

    return result


# ---------------------------------------------------------------------------
# Step 2: Tool detection + hictk integration
# ---------------------------------------------------------------------------

def detect_tools():
    tools = {"hictk": None, "hicstraw": None}

    hictk_path = shutil.which("hictk")
    if hictk_path:
        try:
            proc = subprocess.run(
                ["hictk", "--version"], capture_output=True, text=True, timeout=10
            )
            version_str = proc.stdout.strip() or proc.stderr.strip()
            tools["hictk"] = version_str
        except (subprocess.TimeoutExpired, OSError):
            tools["hictk"] = "found (version unknown)"

    if importlib.util.find_spec("hicstraw"):
        try:
            import hicstraw
            tools["hicstraw"] = getattr(hicstraw, "__version__", "installed")
        except ImportError:
            pass

    return tools


def _run_hictk(args, timeout=30):
    try:
        proc = subprocess.run(
            ["hictk"] + args, capture_output=True, text=True, timeout=timeout
        )
        return proc.returncode, proc.stdout, proc.stderr
    except subprocess.TimeoutExpired:
        return -1, "", "hictk timed out"
    except OSError as e:
        return -1, "", str(e)


def query_hictk_metadata(path):
    rc, stdout, stderr = _run_hictk(["metadata", path, "-f", "json"])
    if rc != 0:
        return {"error": stderr.strip()}
    try:
        return json.loads(stdout)
    except json.JSONDecodeError:
        return {"error": f"Invalid JSON from hictk metadata: {stdout[:200]}"}


def query_hictk_norms(path, resolutions):
    norms_by_res = {}
    for res in resolutions:
        rc, stdout, _stderr = _run_hictk(
            ["dump", "-t", "normalizations", path, "--resolution", str(res)]
        )
        if rc == 0 and stdout.strip():
            norms_by_res[res] = stdout.strip().splitlines()
        else:
            norms_by_res[res] = []
    return norms_by_res


def query_hictk_validate(path):
    rc, stdout, stderr = _run_hictk(["validate", path, "-f", "json"])
    try:
        result = json.loads(stdout)
    except (json.JSONDecodeError, ValueError):
        result = {}
    result["exit_code"] = rc
    if rc != 0 and stderr.strip():
        result["stderr"] = stderr.strip()
    return result


# ---------------------------------------------------------------------------
# Step 3: hicstraw functional testing (--deep only)
# ---------------------------------------------------------------------------

KNOWN_NORMS = ["NONE", "KR", "VC", "VC_SQRT", "SCALE", "GW_SCALE", "INTER_SCALE"]


def _capture_stderr(func):
    """Call func() while capturing C-level stderr (hicstraw prints warnings there)."""
    old_fd = os.dup(2)
    pipe_r, pipe_w = os.pipe()
    os.dup2(pipe_w, 2)
    os.close(pipe_w)
    try:
        result = func()
    finally:
        os.dup2(old_fd, 2)
        os.close(old_fd)
    captured = b""
    while True:
        chunk = os.read(pipe_r, 4096)
        if not chunk:
            break
        captured += chunk
    os.close(pipe_r)
    return result, captured.decode("utf-8", errors="replace")


def test_hicstraw_readability(path, resolutions, filter_res=None):
    try:
        import hicstraw
    except ImportError:
        return {"error": "hicstraw not installed"}

    results = {"open": "FAIL", "resolutions_match": None, "chromosomes_match": None, "norm_tests": {}}

    try:
        hic = hicstraw.HiCFile(path)
        results["open"] = "PASS"
    except Exception as e:
        results["open_error"] = str(e)
        return results

    straw_res = sorted(hic.getResolutions())
    results["straw_resolutions"] = straw_res
    results["resolutions_match"] = straw_res == sorted(resolutions)

    straw_chroms = [c.name for c in hic.getChromosomes() if c.name != "All"]
    results["straw_chromosomes"] = straw_chroms

    test_chrom = None
    for candidate in ["chr1", "1", straw_chroms[0] if straw_chroms else None]:
        if candidate and candidate in straw_chroms:
            test_chrom = candidate
            break

    if not test_chrom:
        results["norm_tests"] = {"error": "No suitable test chromosome found"}
        return results

    test_resolutions = filter_res if filter_res else resolutions
    for res in test_resolutions:
        if res not in straw_res:
            results["norm_tests"][res] = {"error": "resolution not in file"}
            continue

        norm_results = {}
        for norm in KNOWN_NORMS:
            try:
                def _fetch(n=norm, r=res):
                    mzd = hic.getMatrixZoomData(
                        test_chrom, test_chrom, "observed", n, "BP", r
                    )
                    return mzd.getRecords(0, 1_000_000, 0, 1_000_000)

                records, stderr_text = _capture_stderr(_fetch)

                if "did not contain" in stderr_text or "not found" in stderr_text.lower():
                    norm_results[norm] = "FAIL (not found)"
                elif records:
                    norm_results[norm] = "PASS"
                else:
                    norm_results[norm] = "PASS (empty)"
            except Exception as e:
                err = str(e)
                if "not found" in err.lower() or "normalization" in err.lower():
                    norm_results[norm] = "FAIL (not found)"
                else:
                    norm_results[norm] = f"FAIL ({err[:60]})"
        results["norm_tests"][res] = norm_results

    results["test_chromosome"] = test_chrom
    return results


# ---------------------------------------------------------------------------
# Step 4: Compatibility assessment
# ---------------------------------------------------------------------------

def assess_compatibility(header, norms_by_res, _validation, deep_results):
    version = header.get("version", 0)
    verdicts = []

    verdicts.append({
        "tool": "Juicebox / juicer_tools",
        "status": "PASS" if version >= 7 else "WARN",
        "note": "v9 native" if version == 9 else f"v{version}",
    })

    verdicts.append({
        "tool": "hic2cool",
        "status": "PASS" if version <= 8 else "FAIL",
        "note": "requires v7/v8" if version > 8 else f"v{version} supported",
    })

    verdicts.append({
        "tool": "hictk",
        "status": "PASS",
        "note": "v7-v9 supported",
    })

    if version <= 8:
        strawr_status, strawr_note = "PASS", f"v{version} supported"
    else:
        strawr_status, strawr_note = "WARN", "v9 support is version-dependent"
    verdicts.append({"tool": "strawr (R/mariner)", "status": strawr_status, "note": strawr_note})

    verdicts.append({
        "tool": "hicstraw (Python)",
        "status": "PASS" if version >= 8 else "WARN",
        "note": "v8/v9 supported" if version >= 8 else "v7 may have limited support",
    })

    verdicts.append({
        "tool": "cooler ecosystem",
        "status": "N/A",
        "note": "requires mcool conversion (hictk convert)",
    })

    # Recommendations
    recommendations = []

    if header.get("truncated"):
        recommendations.append("File appears truncated (master index offset > file size). Reconvert from source.")

    has_gw_scale = False
    missing_kr_res = []

    # Check norms from hictk (if available)
    if norms_by_res:
        for res, norms in norms_by_res.items():
            if "KR" not in norms:
                missing_kr_res.append(res)
            if "GW_SCALE" in norms:
                has_gw_scale = True

    # Also check norms from deep hicstraw results (if available and hictk wasn't)
    if not norms_by_res and deep_results and "norm_tests" in deep_results:
        for res, norm_map in deep_results["norm_tests"].items():
            if isinstance(norm_map, dict) and "error" not in norm_map:
                if norm_map.get("KR", "").startswith("FAIL"):
                    missing_kr_res.append(res)
                if not norm_map.get("GW_SCALE", "").startswith("FAIL"):
                    has_gw_scale = True

    if missing_kr_res:
        res_str = ", ".join(f"{r}bp" for r in sorted(missing_kr_res))
        recommendations.append(
            f"KR normalization missing at: {res_str}. "
            "Fix: hictk balance ice <file> or convert to mcool + cooler balance --name KR."
        )

    if has_gw_scale:
        recommendations.append(
            "GW_SCALE normalization present. Some tools (older hictk, hic2cool) may choke on this norm type."
        )

    if version == 9:
        recommendations.append(
            "File is v9. Tools that only support v7/v8 (hic2cool, older straw) will fail. "
            "Use hictk for conversions."
        )

    # Check deep results for Class 3 (declared but unreadable norms)
    if deep_results and "norm_tests" in deep_results:
        for res, norm_map in deep_results["norm_tests"].items():
            if isinstance(norm_map, dict) and "error" not in norm_map:
                for norm, status in norm_map.items():
                    if norm != "NONE" and "FAIL" in status:
                        if norms_by_res and res in norms_by_res and norm in norms_by_res[res]:
                            recommendations.append(
                                f"{norm} at {res}bp: declared in header but unreadable by hicstraw "
                                "(Class 3: cross-tool writer incompatibility)."
                            )

    return {"verdicts": verdicts, "recommendations": recommendations}


# ---------------------------------------------------------------------------
# Step 5: Output formatting
# ---------------------------------------------------------------------------

def _fmt_size(n_bytes):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n_bytes) < 1024.0:
            return f"{n_bytes:.1f} {unit}"
        n_bytes /= 1024.0
    return f"{n_bytes:.1f} PB"


def format_report(audit):
    lines = []
    header = audit["header"]

    lines.append("")
    lines.append(_bold(f"=== .hic File Audit: {os.path.basename(header['file_path'])} ==="))
    lines.append("")

    if "error" in header:
        lines.append(f"  {_red('ERROR')}  {header['error']}")
        return "\n".join(lines)

    # Header section
    v = header["version"]
    size = _fmt_size(header["file_size"])
    genome = header.get("genome", "unknown")
    lines.append(f"  [HEADER]  Version: {v} | Genome: {genome} | Size: {size}")

    if "chr_count" in header:
        naming = header.get("chr_naming", "unknown")
        lines.append(f"  [HEADER]  Chromosomes: {header['chr_count']} ({header['chr_range']}) | Naming: {naming}")

    bp_res = header.get("bp_resolutions", [])
    if bp_res:
        res_strs = [str(r) for r in bp_res]
        lines.append(f"  [HEADER]  Resolutions (BP): {', '.join(res_strs)}")

    frag_res = header.get("frag_resolutions", [])
    if frag_res:
        lines.append(f"  [HEADER]  Resolutions (FRAG): {', '.join(str(r) for r in frag_res)}")

    meta = header.get("metadata", {})
    if meta:
        lines.append(f"  [HEADER]  Metadata: {len(meta)} attribute(s)")
        for k, v_val in meta.items():
            display = v_val if len(v_val) <= 60 else v_val[:57] + "..."
            lines.append(f"            {k}: {display}")

    if header.get("truncated"):
        lines.append(f"  [HEADER]  {_red('TRUNCATED')} — master index offset ({header['master_index']}) > file size ({header['file_size']})")

    lines.append("")

    # Tools section
    tools = audit.get("tools", {})
    hictk_ver = tools.get("hictk")
    straw_ver = tools.get("hicstraw")
    lines.append(f"  [TOOLS]   hictk: {hictk_ver or 'not found'}")
    lines.append(f"  [TOOLS]   hicstraw: {straw_ver or 'not found'}")

    if not hictk_ver:
        lines.append(f"            {_yellow('(normalization listing and validation unavailable)')}")
    if not straw_ver:
        lines.append(f"            {_yellow('(functional readability tests unavailable, use --deep with hicstraw)')}")

    lines.append("")

    # Normalizations section
    norms = audit.get("norms_by_res")
    if norms is not None:
        if norms:
            lines.append("  [NORMS]   via hictk dump -t normalizations:")
            for res in sorted(norms.keys()):
                norm_list = norms[res]
                if norm_list:
                    lines.append(f"    {res}bp:  {_bold('  '.join(norm_list))}")
                else:
                    lines.append(f"    {res}bp:  {_red('(none)')}")
        else:
            lines.append(f"  [NORMS]   {_yellow('No resolutions to check')}")
        lines.append("")

    # Validation section
    validation = audit.get("validation")
    if validation is not None:
        rc = validation.get("exit_code", -1)
        if rc == 0:
            lines.append(f"  [VALID]   hictk validate: {_green('PASS')}")
        else:
            stderr = validation.get("stderr", "unknown error")
            lines.append(f"  [VALID]   hictk validate: {_red('FAIL')} — {stderr}")
        lines.append("")

    # Deep readability section
    deep = audit.get("deep_results")
    if deep:
        if "error" in deep:
            lines.append(f"  [DATA]    {_yellow(deep['error'])}")
        else:
            open_status = deep.get("open", "FAIL")
            if open_status != "PASS":
                lines.append(f"  [DATA]    hicstraw open: {_red('FAIL')} — {deep.get('open_error', 'unknown')}")
            else:
                test_chr = deep.get("test_chromosome", "?")
                res_match = deep.get("resolutions_match")
                if res_match is False:
                    lines.append(f"  [DATA]    {_yellow('Resolution mismatch')}: header={bp_res}, straw={deep.get('straw_resolutions')}")

                lines.append(f"  [DATA]    Functional readability ({test_chr} spot-check):")
                for res in sorted(deep.get("norm_tests", {}).keys()):
                    norm_map = deep["norm_tests"][res]
                    if isinstance(norm_map, dict) and "error" not in norm_map:
                        parts = []
                        for norm in KNOWN_NORMS:
                            if norm in norm_map:
                                s = norm_map[norm]
                                tag = _green("PASS") if "PASS" in s else _red("FAIL")
                                parts.append(f"{norm}:{tag}")
                        lines.append(f"    {res}bp:  {'  '.join(parts)}")
                    elif isinstance(norm_map, dict):
                        lines.append(f"    {res}bp:  {_yellow(norm_map.get('error', 'error'))}")
        lines.append("")

    # Compatibility section
    compat = audit.get("compatibility", {})
    verdicts = compat.get("verdicts", [])
    if verdicts:
        lines.append("  [COMPAT]")
        for v in verdicts:
            lines.append(f"    {_status(v['status'])}  {v['tool']} — {v['note']}")
        lines.append("")

    # Recommendations
    recs = compat.get("recommendations", [])
    if recs:
        lines.append("  [ACTION]")
        for r in recs:
            lines.append(f"    * {r}")
        lines.append("")

    # Summary
    n_fail = sum(1 for v in verdicts if v["status"] == "FAIL")
    n_warn = sum(1 for v in verdicts if v["status"] == "WARN")
    n_pass = sum(1 for v in verdicts if v["status"] == "PASS")
    summary_parts = []
    if n_fail:
        summary_parts.append(_red(f"{n_fail} incompatible"))
    if n_warn:
        summary_parts.append(_yellow(f"{n_warn} warning(s)"))
    if n_pass:
        summary_parts.append(_green(f"{n_pass} compatible"))
    if header.get("truncated"):
        summary_parts.append(_red("TRUNCATED"))
    if recs:
        summary_parts.append(f"{len(recs)} recommendation(s)")
    lines.append(f"  [SUMMARY] {', '.join(summary_parts)}")
    lines.append("")

    return "\n".join(lines)


def build_json_result(audit):
    result = {
        "file": audit["header"].get("file_path"),
        "file_size": audit["header"].get("file_size"),
        "version": audit["header"].get("version"),
        "genome": audit["header"].get("genome"),
        "chromosomes": audit["header"].get("chromosomes"),
        "chr_naming": audit["header"].get("chr_naming"),
        "bp_resolutions": audit["header"].get("bp_resolutions"),
        "frag_resolutions": audit["header"].get("frag_resolutions"),
        "metadata": audit["header"].get("metadata"),
        "master_index": audit["header"].get("master_index"),
        "truncated": audit["header"].get("truncated", False),
    }
    if "error" in audit["header"]:
        result["error"] = audit["header"]["error"]

    result["tools_available"] = audit.get("tools", {})

    if audit.get("norms_by_res") is not None:
        result["normalizations"] = {str(k): v for k, v in audit["norms_by_res"].items()}

    if audit.get("validation") is not None:
        result["validation"] = audit["validation"]

    if audit.get("deep_results"):
        result["deep_readability"] = audit["deep_results"]

    compat = audit.get("compatibility", {})
    result["compatibility"] = compat.get("verdicts", [])
    result["recommendations"] = compat.get("recommendations", [])

    return result


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def audit_file(path, tools, deep=False, filter_res=None):
    audit = {}

    # Step 1: Binary header
    header = parse_hic_header(path)
    audit["header"] = header
    audit["tools"] = tools

    if "error" in header:
        audit["compatibility"] = {"verdicts": [], "recommendations": [header["error"]]}
        return audit

    bp_res = header.get("bp_resolutions", [])
    check_res = filter_res if filter_res else bp_res

    # Step 2: hictk queries
    norms_by_res = None
    validation = None
    if tools.get("hictk"):
        norms_by_res = query_hictk_norms(path, check_res)
        validation = query_hictk_validate(path)
    audit["norms_by_res"] = norms_by_res
    audit["validation"] = validation

    # Step 3: Deep functional testing
    deep_results = None
    if deep:
        if tools.get("hicstraw"):
            deep_results = test_hicstraw_readability(path, bp_res, filter_res=check_res)
        else:
            deep_results = {"error": "--deep requires hicstraw (pip install hicstraw)"}
    audit["deep_results"] = deep_results

    # Step 4: Compatibility assessment
    audit["compatibility"] = assess_compatibility(header, norms_by_res or {}, validation, deep_results)

    return audit


def main():
    parser = argparse.ArgumentParser(
        description="Audit .hic files for version, normalization, and compatibility issues."
    )
    parser.add_argument("files", nargs="+", metavar="file", help=".hic file(s) to audit")
    parser.add_argument("--deep", action="store_true", help="Functional readability tests (requires hicstraw)")
    parser.add_argument("--json", dest="json_out", action="store_true", help="Machine-parseable JSON output")
    parser.add_argument("--resolution", type=int, action="append", dest="resolutions", help="Only check specific resolution(s)")
    parser.add_argument("--quiet", action="store_true", help="Only print warnings and errors")
    args = parser.parse_args()

    tools = detect_tools()
    json_results = []

    for i, path in enumerate(args.files):
        if not os.path.exists(path) or not os.path.isfile(path):
            err = {"header": {"file_path": path, "error": f"File not found: {path}"}, "tools": tools}
            err["compatibility"] = {"verdicts": [], "recommendations": []}
            if args.json_out:
                json_results.append(build_json_result(err))
            else:
                print(format_report(err))
            continue

        audit = audit_file(path, tools, deep=args.deep, filter_res=args.resolutions)

        if args.json_out:
            json_results.append(build_json_result(audit))
        elif args.quiet:
            compat = audit.get("compatibility", {})
            recs = compat.get("recommendations", [])
            fails = [v for v in compat.get("verdicts", []) if v["status"] == "FAIL"]
            warns = [v for v in compat.get("verdicts", []) if v["status"] == "WARN"]
            if recs or fails or warns or audit["header"].get("error"):
                print(format_report(audit))
        else:
            if i > 0:
                print("-" * 60)
            print(format_report(audit))

    if args.json_out:
        output = json_results[0] if len(json_results) == 1 else json_results
        print(json.dumps(output, indent=2))

    any_error = any(
        "error" in (r if isinstance(r, dict) else {})
        for r in (json_results if args.json_out else [])
    )
    sys.exit(1 if any_error else 0)


if __name__ == "__main__":
    main()
