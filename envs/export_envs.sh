# envs/export_envs.sh
#!/usr/bin/env bash
set -euo pipefail

OUTDIR="$(cd "$(dirname "$0")" && pwd)"

conda env list --json | python3 -c "
import json, sys, os
data = json.load(sys.stdin)
for p in data['envs']:
    name = os.path.basename(p)
    if name.startswith('env-') and len(name) == 36:
        continue
    print(name + '\t' + p)
" | while IFS=$'\t' read -r name prefix; do
    outfile="${OUTDIR}/${name}.yml"
    echo "Exporting ${name} -> ${outfile}"
    conda env export -p "${prefix}" > "${outfile}"
done

echo "Done. Exported to ${OUTDIR}"
