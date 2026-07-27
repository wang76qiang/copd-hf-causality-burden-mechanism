# ============================================================================
# Module 09 - Figures: shared design system (palette, panel tags, save helper)
# Imported by all make_fig*.py scripts.
# ============================================================================

import os as _os
# --- repository paths -----------------------------------------------------------
# REPO resolves to the repository root (code/<module>/this_file.py -> ../..).
# Override with the REPRO_ROOT environment variable if you run from elsewhere.
REPO = _os.environ.get("REPRO_ROOT", _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..")))
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR = _os.environ.get("GBD_RAW_DIR", _os.path.join(REPO, "data", "raw", "gbd_2021"))
"""Shared design system for T1.1 figures — COPD-HF top-journal suite."""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# ---------- palette ----------
INK     = "#14202E"   # near-black blue
PAPER   = "#FFFFFF"
TEAL    = "#0F7B7C"   # COPD / lung
CRIMSON = "#B53737"   # HF / heart
AMBER   = "#E8A33D"   # highlights / warnings
SLATE   = "#6C7A89"   # neutral
SKY     = "#3D7EAA"   # secondary cool
GREEN   = "#4E8D5B"
LIGHT   = "#E9EDF1"
DIVERGE = {True: CRIMSON, False: TEAL}

mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.edgecolor": INK,
    "axes.labelcolor": INK,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "axes.linewidth": 0.8,
    "xtick.color": INK, "ytick.color": INK,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.frameon": False,
    "figure.facecolor": PAPER,
    "axes.facecolor": PAPER,
    "savefig.facecolor": PAPER,
    "pdf.fonttype": 42,
})

def panel_tag(ax, tag, color=INK, x=-0.14, y=1.06):
    ax.text(x, y, tag, transform=ax.transAxes, fontsize=15, fontweight="bold",
            color=INK, ha="center", va="center")

# --- journal submission mode: suppress all axes (panel) titles ---
_Axes = mpl.axes.Axes
_Axes.set_title = lambda self, *a, **k: None

def strip_spines(ax, keep=("left", "bottom")):
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(side in keep)

def grid(ax, axis="y", alpha=0.35):
    ax.grid(axis=axis, color=LIGHT, lw=0.8, zorder=0)
    ax.set_axisbelow(True)

def supertitle(fig, title, subtitle, x=0.012, ty=0.985, sy=0.955):
    # figure-level titles removed per journal style
    pass

def source_note(fig, note, x=0.012, y=0.008):
    # footnotes removed per minimal-text mode
    pass

def save(fig, name, outdir):
    import os
    fig.savefig(os.path.join(outdir, name + ".png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(outdir, name + ".pdf"), bbox_inches="tight")
    plt.close(fig)
    print("saved", name)

def forest(ax, rows, labels, ors, lo, hi, colors=None, ref=1.0, title="",
           xlab="Odds ratio (log scale)", diamond_idx=None, xlim=None, ann=None):
    """Generic forest plot. diamond_idx: indices drawn as diamonds."""
    y = np.arange(len(rows))[::-1]
    colors = colors or [INK]*len(rows)
    for yi, o, l, h, c, i in zip(y, ors, lo, hi, colors, range(len(rows))):
        if diamond_idx and i in diamond_idx:
            ax.fill([l, o, h, o], [yi, yi+0.28, yi, yi-0.28], color=c, zorder=3)
        else:
            ax.plot([l, h], [yi, yi], color=c, lw=1.6, zorder=2)
            ax.scatter([o], [yi], s=46, color=c, zorder=3, edgecolor="white", lw=0.6)
    ax.axvline(ref, color=SLATE, lw=0.9, ls="--", zorder=1)
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.set_xscale("log")
    if xlim: ax.set_xlim(*xlim)
    ax.set_xlabel(xlab); ax.set_title(title, loc="left")
    if ann:
        for yi, txt in zip(y, ann):
            ax.text(1.02, yi, txt, transform=ax.get_yaxis_transform(),
                    fontsize=7.5, va="center", color=INK)
    strip_spines(ax); grid(ax, "x")
