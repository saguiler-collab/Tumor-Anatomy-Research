"""One journal figure style, applied by every figure script.

WHY A SHARED MODULE. Three scripts draw the paper's figures. Styling each one separately
guarantees they drift -- different fonts, different weights, different greys -- and a reader
notices that before they notice the data. The style also has to be changeable in one place,
because "make the figures look like a journal" is a request that arrives more than once.

THE TARGET. Printed figures in a journal, and the GraphPad-style panels common in
experimental papers: Arial, heavy black axes, ticks pointing outward, saturated fills with
black edges, italic for gene and cell-type symbols, no gridlines, no background tint. The
default matplotlib look -- thin grey spines, muted palette, inward ticks -- reads as a
software default rather than a figure someone made.
"""
from __future__ import annotations

import matplotlib as mpl

#: Saturated, print-safe, and distinguishable in greyscale by lightness as well as hue.
#: Black edges on every patch, which is what makes a fill read as a journal figure rather
#: than a dashboard.
INK = "#000000"
GREY = "#7f7f7f"
LIGHT = "#d9d9d9"
BLUE = "#2e5fa3"
RED = "#c0392b"
YELLOW = "#d4c020"
ORANGE = "#e08214"
GREEN = "#4d8a52"

#: Ordered palette for categorical series.
CYCLE = [BLUE, RED, YELLOW, ORANGE, GREEN, GREY]


def apply() -> None:
    """Set the rcParams. Call once, before any figure is created."""
    mpl.rcParams.update({
        # ---- type -------------------------------------------------------------------
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 12.5,
        "axes.titlesize": 13,
        "axes.labelsize": 12.5,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 11.5,
        "figure.titlesize": 14.5,
        # Journals set panel letters and titles in bold; everything else stays regular.
        "axes.titleweight": "bold",
        "axes.labelweight": "regular",
        # ---- axes -------------------------------------------------------------------
        # Heavy black spines, and only the two that carry information.
        "axes.linewidth": 1.8,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.facecolor": "white",
        "figure.facecolor": "white",
        "axes.grid": False,
        "axes.prop_cycle": mpl.cycler(color=CYCLE),
        # ---- ticks ------------------------------------------------------------------
        # Outward and substantial. Inward ticks sit on top of the data.
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 6.5,
        "ytick.major.size": 6.5,
        "xtick.major.width": 1.8,
        "ytick.major.width": 1.8,
        "xtick.color": INK,
        "ytick.color": INK,
        # ---- marks ------------------------------------------------------------------
        "lines.linewidth": 2.0,
        "lines.markersize": 8,
        "lines.markeredgewidth": 1.2,
        "patch.edgecolor": INK,
        "patch.linewidth": 1.0,
        "patch.force_edgecolor": True,
        "errorbar.capsize": 3,
        # ---- legend -----------------------------------------------------------------
        "legend.frameon": False,
        "legend.handlelength": 1.6,
        "legend.borderpad": 0.3,
        "legend.labelspacing": 0.35,
        # ---- output -----------------------------------------------------------------
        "figure.dpi": 110,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
        # Keep text as text in the PDF so a typesetter can restyle it.
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def italicise(ax, axis: str = "y") -> None:
    """Italicise tick labels. **For GENE SYMBOLS ONLY.**

    Convention: gene symbols are italic (*PITX1*, *CXCR4*); cell-type names, protein names and
    software names are NOT. None of this project's figures label genes -- they label cell types
    (T cell, B cell) and methods (`dwls`, `music`) -- so nothing here calls this, and applying
    it to match the look of a figure that happens to plot genes would be a convention error
    rather than a style choice.
    """
    labels = ax.get_yticklabels() if axis == "y" else ax.get_xticklabels()
    for t in labels:
        t.set_fontstyle("italic")
