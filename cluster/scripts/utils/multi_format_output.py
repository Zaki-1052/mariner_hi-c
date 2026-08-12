# cluster/scripts/utils/multi_format_output.py
"""Multi-format plot output utility (Python port of scripts/utils/multi_format_output.R).

Patches matplotlib.figure.Figure.savefig within scope so that .png writes also emit
.svg/.pdf/.jpg siblings of the same figure. Other extensions pass through unchanged.

Patching Figure.savefig (rather than plt.savefig) covers both calling patterns:
- pyplot.savefig() internally calls gcf().savefig() — i.e. Figure.savefig — so
  Popay's box/strip/heat/line that call plt.savefig still produce 4 formats.
- Code that calls fig.savefig() directly (e.g. plotting.stacked) is also intercepted.

Pre-creates a per-figure subfolder so output matches the R utility's
use_subfolders=TRUE layout: {parent_dir}/{name}/{name}.{png,pdf,svg,jpg}

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

import matplotlib.figure


@contextmanager
def multi_format_savefig():
    # Patches Figure.savefig within scope. Triggered on .png writes only;
    # .pdf/.jpg/.svg writes pass through unchanged. Popay's box/strip write
    # both .png then .pdf — only the .png triggers 4-format emission, so the
    # subsequent .pdf write overwrites our emitted PDF with identical content
    # (harmless). plotting.stacked() uses fig.savefig(.png) then
    # fig.savefig(.pdf); same behavior — single round of 4-format emission.
    original = matplotlib.figure.Figure.savefig

    def patched(self, fname, *args, **kwargs):
        base, ext = os.path.splitext(str(fname))
        result = original(self, fname, *args, **kwargs)
        if ext.lower() == '.png':
            vector_kwargs = {k: v for k, v in kwargs.items() if k != 'dpi'}
            original(self, base + '.svg', *args, **vector_kwargs)
            original(self, base + '.pdf', *args, **vector_kwargs)
            jpg_kwargs = dict(kwargs)
            jpg_kwargs.setdefault('dpi', 300)
            original(self, base + '.jpg', *args, **jpg_kwargs)
        return result

    matplotlib.figure.Figure.savefig = patched
    try:
        yield
    finally:
        matplotlib.figure.Figure.savefig = original


def figure_subfolder(parent_dir, name) -> Path:
    # Mirror multi_format_output.R use_subfolders=TRUE: {parent_dir}/{name}/
    # so the plotting function (passed this as out_dir) writes
    # {parent_dir}/{name}/{name}.{ext}
    sub = Path(parent_dir) / name
    sub.mkdir(parents=True, exist_ok=True)
    return sub
