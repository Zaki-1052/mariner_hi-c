import json, re

REMOVE_NAMES = {
    # Redacted
    "Nicole Castro", "Alessandro Cirulli", "Gabriella Ching", "Griffin Hurst",
    "Elysa Loraine Lebig", "Itzel Gonzalez Velazquez", "Michaela Goodman",
    "Andrea Melendez", "Omar Mokhashi", "Arya Lalezarzadeh", "Andrea Mota",
    "Leilani Rivera", "Matthew Kim", "Pooja Parthasarathy", "Elijah Khalil Rosales",
    "Kelsey Schilling", "Risha Sharma", "Madison Wong", "Katelyn Wong",
    "Grace Smith", "Gabriel Soberón Nelson", "Nomy Xin", "Dyllan Mead", "Julia Lee",
    # Structurally broken / not abstracts
    "Amanda Salatino", "Sophie Zhang", "Eamon Lee", "Lan Gao",
    # Too short / missing major components
    "Breanna Fraire", "Kayanne Tran", "Ainsley Gibson", "Cara Chan",
    "Emily Huang", "Stephen Huang", "Natalee Chin", "Richard Madsen",
}

INPUT_FILE = "abstracts.json"
OUTPUT_FILE = "abstracts_filtered.json"

# The file is concatenated JSON objects, not a JSON array.
# Parse by finding each top-level { ... } block.
with open(INPUT_FILE, "r", encoding="utf-8") as f:
    raw = f.read()

entries = []
decoder = json.JSONDecoder()
idx = 0
while idx < len(raw):
    match = re.match(r'\s*', raw[idx:])
    if match:
        idx += match.end()
    if idx >= len(raw):
        break
    try:
        obj, end = decoder.raw_decode(raw, idx)
        if isinstance(obj, dict):
            entries.append(obj)
        idx = end
    except json.JSONDecodeError:
        idx += 1

before = len(entries)
filtered = [e for e in entries if e.get("name", "") not in REMOVE_NAMES]
after = len(filtered)

with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
    for entry in filtered:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")

print(f"Removed {before - after} entries ({before} -> {after})")
print(f"Written to {OUTPUT_FILE}")

# Sanity check: flag any names that weren't found in the data
found = {e.get("name", "") for e in entries}
missed = REMOVE_NAMES - found
if missed:
    print(f"\n⚠️  These names were in REMOVE_NAMES but not found in the file:")
    for m in sorted(missed):
        print(f"   - {m}")