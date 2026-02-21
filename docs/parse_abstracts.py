import re
import sys
import json


def parse_abstracts(text):
    raw = [l.rstrip() for l in text.splitlines()]

    # Strip page numbers and blank/whitespace-only lines so structure is flat
    lines = [
        l for l in raw
        if l.strip() and not re.fullmatch(r'\s*\d+\s*', l)
    ]

    mentor_re = re.compile(r'mentored\s+by', re.IGNORECASE)
    mentor_idx = [i for i, l in enumerate(lines) if mentor_re.match(l.strip())]

    abstracts = []
    for n, mi in enumerate(mentor_idx):
        # Name and major are guaranteed to be the 2 lines immediately before mentor
        name  = lines[mi - 2].strip() if mi >= 2 else ""
        major = lines[mi - 1].strip() if mi >= 1 else ""
        mentor = lines[mi].strip()

        # Body runs from after mentor line up to (but not including)
        # the name/major of the next abstract, i.e. stop at next_mi - 2
        if n + 1 < len(mentor_idx):
            next_mi = mentor_idx[n + 1]
            body_lines = lines[mi + 1 : next_mi - 2]
        else:
            body_lines = lines[mi + 1 :]

        # Collapse everything into one clean paragraph
        body = " ".join(l.strip() for l in body_lines)

        abstracts.append({
            "name":   name,
            "major":  major,
            "mentor": mentor,
            "body":   body,
        })

    return abstracts


def format_text(abstracts):
    sep = "-" * 80
    out = []
    for a in abstracts:
        out.append(f"NAME:   {a['name']}")
        out.append(f"MAJOR:  {a['major']}")
        out.append(f"MENTOR: {a['mentor']}")
        out.append(f"BODY:\n{a['body']}")
        out.append(sep)
    return "\n".join(out)


def main():
    if len(sys.argv) < 2:
        print("Usage: python parse_abstracts_fixed.py input.txt [output.txt] [--json]")
        sys.exit(1)

    in_file  = sys.argv[1]
    out_file = next((a for a in sys.argv[2:] if not a.startswith("--")), None)
    as_json  = "--json" in sys.argv

    with open(in_file, encoding="utf-8") as f:
        text = f.read()

    abstracts = parse_abstracts(text)
    result = (
        "\n".join(json.dumps(a, indent=2, ensure_ascii=False) for a in abstracts)
        if as_json
        else format_text(abstracts)
    )

    if out_file:
        with open(out_file, "w", encoding="utf-8") as f:
            f.write(result)
        print(f"Parsed {len(abstracts)} abstracts → {out_file}")
    else:
        print(result)


if __name__ == "__main__":
    main()