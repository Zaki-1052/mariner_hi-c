# cluster/scripts/inspect_emissions.py
# Print ChromHMM emission matrices with ON/mid/OFF annotations.
# Usage: python3 scripts/inspect_emissions.py [chromHMM_subdir] [nstates...]
#   python3 scripts/inspect_emissions.py chromHMM_9mark_intersect 15 18
#   python3 scripts/inspect_emissions.py  # defaults: chromHMM_9mark_intersect, 15 18

import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
CLUSTER_DIR = SCRIPT_DIR.parent
OUT_BASE = CLUSTER_DIR / 'outputs' / 'bap1_late'


def load_rename(path):
    rename = {}
    if path.exists():
        with open(path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    rename[parts[0]] = parts[1]
    return rename


def print_table(chm_subdir, nstates):
    chm_dir = OUT_BASE / chm_subdir
    emi_path = chm_dir / 'learned_model_{}'.format(nstates) / 'emissions_{}.txt'.format(nstates)
    rename_path = chm_dir / '{}state_rename_cerebellum.txt'.format(nstates)

    if not emi_path.exists():
        print('  NOT FOUND: {}'.format(emi_path))
        return

    emi = pd.read_csv(emi_path, sep='\t', index_col=0)
    emi.index = ['E{}'.format(i) for i in emi.index]
    rename = load_rename(rename_path)

    marks = list(emi.columns)
    short_marks = [m.replace('H3K', 'K').replace('H2AK', 'K') for m in marks]

    print('\n' + '=' * 120)
    print('{} k={}: Emission Probabilities'.format(chm_subdir, nstates))
    if rename:
        print('  (rename file: {})'.format(rename_path.name))
    else:
        print('  (no rename file found)')
    print('=' * 120)

    hdr = '{:5s} {:28s}'.format('State', 'Name')
    for sm in short_marks:
        hdr += ' {:>8s}'.format(sm)
    print(hdr)
    print('-' * len(hdr))

    for state in emi.index:
        name = rename.get(state, '???')
        row = emi.loc[state]
        line = '{:5s} {:28s}'.format(state, name)
        for m in marks:
            v = row[m]
            if v > 0.5:
                tag = '\033[1;32m{:7.3f}+\033[0m'.format(v)
            elif v > 0.2:
                tag = '\033[1;33m{:7.3f}~\033[0m'.format(v)
            else:
                tag = '{:7.3f} '.format(v)
            line += ' ' + tag
        print(line)

    print()
    print('Legend: + = ON (>0.5, green)  ~ = mid (0.2-0.5, yellow)  plain = OFF (<0.2)')


def main():
    args = sys.argv[1:]
    if not args:
        chm_subdir = 'chromHMM_9mark_intersect'
        nstates_list = ['15', '18']
    elif len(args) == 1:
        chm_subdir = args[0]
        nstates_list = ['15', '18']
    else:
        chm_subdir = args[0]
        nstates_list = args[1:]

    for ns in nstates_list:
        print_table(chm_subdir, ns)


if __name__ == '__main__':
    main()
