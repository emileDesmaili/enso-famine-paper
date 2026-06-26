"""
33_fig_v2_onset_panels.py
==========================

Combined two-panel ED figure:
  Panel a: Central Europe coefficients on lag/contemp/lead NINO3.4 in the
           famine-onset model (from 29_v2_onset_lag_lead.R).
  Panel b: CE coefficient on contemporaneous NINO3.4 with vs without
           Swiss-driven famine years (from 30_v2_onset_excl_swiss.R).

Reads:
  analysis/output/data/figED_v2_onset_lag_lead.csv
  analysis/output/data/figED_v2_onset_excl_swiss.csv

Saves:
  analysis/output/figures/extended data/figED_v2_onset_panels.pdf
"""
from __future__ import annotations
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))

from assemble_figures import (  # noqa: E402
    _label, _title, _polish, _pct_fmt, FS,
)

OUT_DATA = SCRIPT_DIR / "output" / "data"
OUT_ED   = SCRIPT_DIR / "output" / "figures" / "extended data"
OUT_ED.mkdir(parents=True, exist_ok=True)

TIMING_ORDER = ["Lag (t-1)", "Contemporaneous (t)", "Lead (t+1)"]
MODEL_ORDER  = ["No FEs", "FEs", "FEs + Controls"]
COLORS_MODEL = {"No FEs": "gray",
                "FEs":    "black",
                "FEs + Controls": "firebrick"}
SAMPLE_ORDER  = ["Original", "Excl Swiss famines"]
COLORS_SAMPLE = {"Original": "black", "Excl Swiss famines": "firebrick"}


def _panel_a(ax):
    df = pd.read_csv(OUT_DATA / "figED_v2_onset_lag_lead.csv")
    x = np.arange(len(TIMING_ORDER))
    offsets = np.linspace(-0.22, 0.22, len(MODEL_ORDER))
    for mod, off in zip(MODEL_ORDER, offsets):
        sub = df[df["model"] == mod].set_index("timing").reindex(TIMING_ORDER)
        ax.errorbar(x + off, sub["estimate"],
                    yerr=[sub["estimate"] - sub["conf.low"],
                          sub["conf.high"] - sub["estimate"]],
                    fmt="o", color=COLORS_MODEL[mod], elinewidth=2.4,
                    capsize=0, ms=10, label=mod, zorder=3)
    ax.axhline(0, color="gray", lw=0.8, linestyle="--")
    ax.set_xticks(x); ax.set_xticklabels(TIMING_ORDER)
    ax.set_ylabel("Δ Famine-start probability\nper +1°C NINO3.4")
    _pct_fmt(ax)
    _title(ax, "Central Europe: ENSO timing")
    ax.legend(loc="best", frameon=False, fontsize=FS - 2)
    _polish(ax)


def _panel_b(ax):
    df = pd.read_csv(OUT_DATA / "figED_v2_onset_excl_swiss.csv")
    ce = df[df["Region"] == "Central Europe"]
    x = np.arange(len(MODEL_ORDER))
    offsets = np.linspace(-0.18, 0.18, len(SAMPLE_ORDER))
    for sample, off in zip(SAMPLE_ORDER, offsets):
        sub = ce[ce["sample"] == sample].set_index("model").reindex(MODEL_ORDER)
        ax.errorbar(x + off, sub["estimate"],
                    yerr=[sub["estimate"] - sub["conf.low"],
                          sub["conf.high"] - sub["estimate"]],
                    fmt="o", color=COLORS_SAMPLE[sample], elinewidth=2.4,
                    capsize=0, ms=10, label=sample, zorder=3)
    ax.axhline(0, color="gray", lw=0.8, linestyle="--")
    ax.set_xticks(x); ax.set_xticklabels(MODEL_ORDER)
    ax.set_ylabel("Δ CE famine-start probability\nper +1°C NINO3.4")
    _pct_fmt(ax)
    _title(ax, "Central Europe: original vs Swiss-excluded chronology")
    ax.legend(loc="best", frameon=False, fontsize=FS - 2)
    _polish(ax)


def main():
    fig = plt.figure(figsize=(22, 8))
    gs  = gridspec.GridSpec(1, 2, figure=fig, wspace=0.30)
    ax_a = fig.add_subplot(gs[0, 0]); _panel_a(ax_a)
    _label(ax_a, "a", x=-0.12, y=1.06)
    ax_b = fig.add_subplot(gs[0, 1]); _panel_b(ax_b)
    _label(ax_b, "b", x=-0.12, y=1.06)
    out = OUT_ED / "figED_v2_onset_panels.pdf"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
