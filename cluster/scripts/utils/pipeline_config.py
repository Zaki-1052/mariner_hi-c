# cluster/scripts/utils/pipeline_config.py
"""Shared configuration and utility functions for the Hi-C loop clustering pipeline.

Centralizes constants and functions that were previously copy-pasted across
multiple scripts. All env-var defaults match the BAP1-KO cerebellum (late)
analysis so existing runs produce identical output without setting anything new.

A new researcher overrides only the variables that differ for their project
via shell .conf files (sourced before running) or direct env-var exports.
"""
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Shared coordinate column names
# ---------------------------------------------------------------------------
COORD_COLS = ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2']


# ---------------------------------------------------------------------------
# Condition / replicate configuration
# ---------------------------------------------------------------------------

def get_condition_cols():
    # type: () -> Tuple[str, str]
    """Return (cond1_col, cond2_col) from env vars, default ctrl_merge/mut_merge."""
    c1 = os.environ.get('CLUSTER_COND1_COL', 'ctrl_merge')
    c2 = os.environ.get('CLUSTER_COND2_COL', 'mut_merge')
    return c1, c2


def get_data_cols():
    # type: () -> List[str]
    """Return [cond1_col, cond2_col] as a list."""
    return list(get_condition_cols())


def get_normalize_col():
    # type: () -> str
    """Return the column to normalize against (defaults to cond1)."""
    return get_condition_cols()[0]


def get_replicate_cols():
    # type: () -> Tuple[List[str], List[str]]
    """Return (ctrl_reps, mut_reps) from env vars."""
    ctrl = os.environ.get('CLUSTER_CTRL_REPS', 'ctrl_M1,ctrl_M2,ctrl_M3').split(',')
    mut = os.environ.get('CLUSTER_MUT_REPS', 'mut_M1,mut_M2,mut_M3').split(',')
    return ctrl, mut


def get_resolutions():
    # type: () -> List[int]
    """Return list of resolution values (kb) from env var."""
    raw = os.environ.get('CLUSTER_RESOLUTIONS', '5,10,25')
    return [int(x.strip()) for x in raw.split(',')]


# ---------------------------------------------------------------------------
# Cluster ordering and biological labels
# ---------------------------------------------------------------------------

def get_cluster_order(k=None):
    # type: (Optional[int]) -> List[str]
    """Return ['clust1', ..., 'clustK'] computed from CLUSTER_K env var."""
    if k is None:
        k = int(os.environ.get('CLUSTER_K', '6'))
    return ['clust{}'.format(i) for i in range(1, k + 1)]


def get_bio_order(cluster_order=None):
    # type: (Optional[List[str]]) -> List[str]
    """Return biological display ordering from CLUSTER_BIO_ORDER env var.
    Falls back to computational order if unset."""
    raw = os.environ.get('CLUSTER_BIO_ORDER', '')
    if raw.strip():
        return [x.strip() for x in raw.split(',')]
    if cluster_order is not None:
        return list(cluster_order)
    return get_cluster_order()


def get_bio_names(bio_order=None):
    # type: (Optional[List[str]]) -> Dict[str, str]
    """Return {cluster_id: display_name} from CLUSTER_BIO_NAME_* env vars.
    Falls back to the cluster ID itself if no env var is set for that cluster."""
    if bio_order is None:
        bio_order = get_bio_order()
    names = {}
    for c in bio_order:
        raw = os.environ.get('CLUSTER_BIO_NAME_{}'.format(c), c)
        names[c] = raw.replace('\\n', '\n')
    return names


def get_cluster_direction():
    # type: () -> Dict[str, str]
    """Return {cluster_id: direction_label} from CLUSTER_DIRECTION_* env vars.
    Falls back to cluster ID if unset."""
    order = get_bio_order()
    d = {}
    for c in order:
        d[c] = os.environ.get('CLUSTER_DIRECTION_{}'.format(c), c)
    return d


# ---------------------------------------------------------------------------
# Size threshold for subgroup splitting
# ---------------------------------------------------------------------------

def get_size_threshold():
    # type: () -> int
    """Return loop size threshold (bp) for short/long subgroup splitting."""
    return int(os.environ.get('CLUSTER_SIZE_THRESHOLD', '800000'))


# ---------------------------------------------------------------------------
# ChromHMM key state lookup
# ---------------------------------------------------------------------------

def get_key_state_id(default='11_Polycomb'):
    # type: (str) -> str
    """Return the enrichment-table row label for the 'key' state (typically Polycomb).
    Reads from CLUSTER_KEY_STATE_ID env var; defaults to the BAP1-KO 12-state label."""
    return os.environ.get('CLUSTER_KEY_STATE_ID', default)


def get_key_state2_id(default='12_Bivalent_Enhancer'):
    # type: (str) -> str
    """Return the enrichment-table row label for the secondary key state."""
    return os.environ.get('CLUSTER_KEY_STATE2_ID', default)


def find_state_row(df_index, pattern, fallback=''):
    # type: (List, str, str) -> str
    """Search a DataFrame index for a row containing `pattern` (case-insensitive).
    Returns the first match, or `fallback` if none found."""
    pat_lower = pattern.lower()
    for label in df_index:
        if pat_lower in str(label).lower():
            return label
    return fallback


# ---------------------------------------------------------------------------
# BigWig dictionary construction
# ---------------------------------------------------------------------------

_DEFAULT_BIGWIGS = {
    'H3K27ac_ctrl':   'H3K27acCtrl.bw',
    'H3K27ac_mut':    'H3K27acMut.bw',
    'H3K27me3_ctrl':  'H3K27me3Ctrl.bw',
    'H3K27me3_mut':   'H3K27me3Mut.bw',
    'H2AK119ub_ctrl': 'H2AK119ubCtrl.bw',
    'H2AK119ub_mut':  'H2AK119ubMut.bw',
    'H3K27me1_ctrl':  'H3K27me1Ctrl.bw',
    'H3K27me1_mut':   'H3K27me1Mut.bw',
}

_DEFAULT_MARK_COLORS = {
    'H3K27ac':   'Blues',
    'H3K27me3':  'Reds',
    'H2AK119ub': 'Greens',
    'H3K27me1':  'Purples',
}


def build_bigwig_dict(bigwig_dir):
    # type: (Path) -> Dict[str, Path]
    """Build {label: Path} dict from CLUSTER_BW_* env vars or defaults.

    Env vars follow the pattern CLUSTER_BW_<KEY>=<filename>, where <KEY>
    is the uppercased, dot-stripped version of the dict key.
    E.g. CLUSTER_BW_H3K27AC_CTRL="my_K27ac_ctrl.bw"
    """
    bigwig_dir = Path(bigwig_dir)
    result = {}
    for key, default_fname in _DEFAULT_BIGWIGS.items():
        env_key = 'CLUSTER_BW_{}'.format(key.upper().replace('.', ''))
        fname = os.environ.get(env_key, default_fname)
        result[key] = bigwig_dir / fname
    return result


def build_vmax_groups(bigwig_dict):
    # type: (Dict[str, Path]) -> List[List[str]]
    """Derive ctrl/mut vmax-sharing groups from a BigWig dict.

    Groups keys that share the same mark prefix (everything before the last _).
    E.g. 'H3K27ac_ctrl' and 'H3K27ac_mut' form one group.
    """
    from collections import OrderedDict
    groups = OrderedDict()  # type: Dict[str, List[str]]
    for key in bigwig_dict:
        mark = key.rsplit('_', 1)[0]
        groups.setdefault(mark, []).append(key)
    return [v for v in groups.values() if len(v) > 1]


def build_color_dict(bigwig_dict):
    # type: (Dict[str, Path]) -> Dict[str, str]
    """Derive {label: colormap_name} from mark prefixes.

    Reads CLUSTER_MARK_COLORS as a comma-separated 'mark=cmap' list,
    e.g. "H3K27ac=Blues,H3K27me3=Reds". Falls back to built-in defaults.
    """
    raw = os.environ.get('CLUSTER_MARK_COLORS', '')
    overrides = {}
    if raw.strip():
        for pair in raw.split(','):
            if '=' in pair:
                m, c = pair.split('=', 1)
                overrides[m.strip()] = c.strip()

    mark_colors = dict(_DEFAULT_MARK_COLORS)
    mark_colors.update(overrides)

    result = {}
    for key in bigwig_dict:
        mark = key.rsplit('_', 1)[0]
        result[key] = mark_colors.get(mark, 'Greys')
    return result


# ---------------------------------------------------------------------------
# deepTools header parser (previously copy-pasted in 5 scripts)
# ---------------------------------------------------------------------------

def parse_header(filepath):
    # type: (Path) -> Dict
    """Parse the 3-line header of a deepTools computeMatrix --outFileNameMatrix file.

    Returns dict with keys:
      cluster_names  -- list of cluster names (stripped of BED filename suffixes)
      cluster_sizes  -- list of ints (anchors per cluster)
      bin_size       -- int (bp per bin)
      n_bins         -- int (total bins per BigWig per region)
      bigwig_order   -- list of BigWig labels in column order
    """
    with open(str(filepath)) as f:
        line1 = f.readline().strip()
        line2 = f.readline().strip()
        line3 = f.readline().strip()

    tokens1 = [t.replace('#', '') for t in line1.split('\t') if ':' in t]
    cluster_names = []
    cluster_sizes = []
    for t in tokens1:
        name_part, size_part = t.split(':')
        cleaned = name_part
        for suffix in ('_oriented_anchors.bed', '_anchors.bed', '.bed'):
            cleaned = cleaned.replace(suffix, '')
        cluster_names.append(cleaned)
        cluster_sizes.append(int(size_part))

    params = {}
    for p in line2.replace('#', '').split('\t'):
        if ':' in p:
            k, v = p.split(':', 1)
            params[k.strip()] = int(v.strip())
    bin_size = params['bin size']
    n_bins = (params['upstream'] + params['downstream']) // bin_size

    group_set = set(tokens1)
    bigwig_labels = [t for t in line3.split('\t') if t not in group_set]
    bigwig_order = list(dict.fromkeys(bigwig_labels))

    return {
        'cluster_names': cluster_names,
        'cluster_sizes': cluster_sizes,
        'bin_size': bin_size,
        'n_bins': n_bins,
        'bigwig_order': bigwig_order,
    }


# ---------------------------------------------------------------------------
# Expected count validation (optional)
# ---------------------------------------------------------------------------

def get_expected_total():
    # type: () -> Optional[int]
    """Return expected total loop count from CLUSTER_EXPECTED_TOTAL, or None."""
    raw = os.environ.get('CLUSTER_EXPECTED_TOTAL', '')
    if raw.strip():
        return int(raw.strip())
    return None


def get_expected_per_res():
    # type: () -> Optional[Dict[int, int]]
    """Return {res_kb: expected_count} from CLUSTER_EXPECTED_PER_RES, or None.
    Format: "5=7901,10=14553,25=16890"
    """
    raw = os.environ.get('CLUSTER_EXPECTED_PER_RES', '')
    if not raw.strip():
        return None
    result = {}
    for pair in raw.split(','):
        if '=' in pair:
            k, v = pair.split('=', 1)
            result[int(k.strip())] = int(v.strip())
    return result


# ---------------------------------------------------------------------------
# Palette for condition columns (used by clustering visualizations)
# ---------------------------------------------------------------------------

def get_palette():
    # type: () -> Dict[str, str]
    """Return color palette keyed by both full column names and stripped labels."""
    c1, c2 = get_condition_cols()
    c1_color = os.environ.get('CLUSTER_COND1_COLOR', 'darkgrey')
    c2_color = os.environ.get('CLUSTER_COND2_COLOR', 'forestgreen')
    palette = {c1: c1_color, c2: c2_color}
    c1_short = c1.replace('_merge', '')
    c2_short = c2.replace('_merge', '')
    if c1_short != c1:
        palette[c1_short] = c1_color
    if c2_short != c2:
        palette[c2_short] = c2_color
    return palette
