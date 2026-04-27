# cluster/scripts/utils/multi_format_output.py
"""Multi-format plot output utility (Python port of scripts/utils/multi_format_output.R).

Augments matplotlib's plt.savefig() to emit PDF + SVG + JPG (alongside the
caller-requested PNG) when used inside the multi_format_savefig() context
manager, and pre-creates a per-figure subfolder so output matches the R
utility's use_subfolders=TRUE layout:
    {parent_dir}/{name}/{name}.{png,pdf,svg,jpg}

This lets us reuse Popay's plotting.py and cluster_tools.elbow() unchanged
(they call plt.savefig internally) while still delivering Illustrator-editable
SVG + presentation-ready JPG for downstream use.

Usage:
    from multi_format_output import multi_format_savefig, figure_subfolder

    fig_dir = figure_subfolder('cluster/bap1_late/figures', 'heatmap')
    with multi_format_savefig():
        plotting.heat(out_dir=str(fig_dir), out_name='heatmap', ...)
    # produces: cluster/bap1_late/figures/heatmap/heatmap.{png,pdf,svg,jpg}
"""
import os
from contextlib import contextmanager
from pathlib import Path

import matplotlib.pyplot as plt


@contextmanager
def multi_format_savefig():
    # Patches plt.savefig within scope so any .png write also emits .svg/.pdf/.jpg
    # of the same active figure. .pdf/.jpg/.svg writes pass through unchanged.
    # Popay's heat/line/box/strip write .png then .pdf; the .pdf overwrites our
    # vector .pdf with identical content (harmless). cluster_tools.elbow() only
    # writes .png, so this wrapper is what gives it pdf/svg/jpg siblings.
    original = plt.savefig

    def patched(fname, *args, **kwargs):
        base, ext = os.path.splitext(str(fname))
        result = original(fname, *args, **kwargs)
        if ext.lower() == '.png':
            vector_kwargs = {k: v for k, v in kwargs.items() if k != 'dpi'}
            original(base + '.svg', *args, **vector_kwargs)
            original(base + '.pdf', *args, **vector_kwargs)
            jpg_kwargs = dict(kwargs)
            jpg_kwargs.setdefault('dpi', 300)
            original(base + '.jpg', *args, **jpg_kwargs)
        return result

    plt.savefig = patched
    try:
        yield
    finally:
        plt.savefig = original


def figure_subfolder(parent_dir, name) -> Path:
    # Mirror multi_format_output.R use_subfolders=TRUE: {parent_dir}/{name}/
    # so the plotting function (passed this as out_dir) writes
    # {parent_dir}/{name}/{name}.{ext}
    sub = Path(parent_dir) / name
    sub.mkdir(parents=True, exist_ok=True)
    return sub
