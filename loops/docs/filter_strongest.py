import json, re

# Margaret Jones has two abstracts; match on mentor to get the bioinformatics one
KEEP_NAMES = {
    # Life Sciences / Biomedical
    "Kathleen Bai", "Sebastian Gastelum", "Ethan Lu", "Kamryn Conway",
    "Andrew Dallape", "Daniella Bandari", "Saloni Dangre", "Cadence Seymour",
    "Brandon Saiki", "Amber Lawrence", "Thu Nguyen", "Sophia Trujillo",
    "Julianna Vega Perez", "Nehme Lahoud", "Vivian Chen", "Ruman Das",
    "Celeste Morales", "Amelia Orgill", "Alana Tamayo", "Irum Hasan",
    "Matthew Nunes",
    # Computational / Engineering
    "Olimpia Carrioli", "Manu Bhat", "Alisha Foster",
    "Ifunanya Okoroma", "Marfred Barrera", "Juliette Hamid", "Yang Han",
    "Brina Nguyen",
    # Public Health
    "Jackie Aviles", "Karla Garcia", "Melanie Gallegos", "Sarah Plummer",
    "Michelle Griffith",
    # Social Science / Psychology / Neuroscience
    "Avery Charneski", "Simon Roberts", "Nathalie Gider", "Moumen Gabir",
    "David Ngan", "Julie Qian",
    # Added — originally missed
    "Rebecca Tseng", "Ashley Becker", "Leena Kang", "Connor Stratman",
}

# Margaret Jones appears twice; keep only the bioinformatics abstract (Mesirov mentor)
JONES_MENTOR_KEEP = "Jill Mesirov"

INPUT_FILE = "abstracts.json"
OUTPUT_FILE = "abstracts_strongest.json"

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

def should_keep(entry):
    name = entry.get("name", "")
    mentor = entry.get("mentor", "")
    # Special case: Margaret Jones has two entries
    if name == "Margaret Jones":
        return JONES_MENTOR_KEEP in mentor
    return name in KEEP_NAMES

filtered = [e for e in entries if should_keep(e)]

with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
    for entry in filtered:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")

print(f"Kept {len(filtered)} of {len(entries)} entries")
print(f"Written to {OUTPUT_FILE}")

# Sanity check
found = {e.get("name", "") for e in filtered}
expected = KEEP_NAMES | {"Margaret Jones"}
missed = expected - found
if missed:
    print(f"\n⚠️  Expected but not found in output:")
    for m in sorted(missed):
        print(f"   - {m}")