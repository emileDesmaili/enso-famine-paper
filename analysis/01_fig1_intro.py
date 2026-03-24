"""
Figure 1: ENSO and famines in early modern Europe  (Nature-style combined panel)
Panels:
  A – NINO3.4 timeseries with shaded famine events
  B – Regional famine chronology (Gantt) + count of regions with famine
  C – Geographic bar-map of famine years and periods per region

Output → analysis/output/figures/main/fig1_combined.pdf
Individual panels also saved as fig1A/B/C_*.pdf for flexibility.
"""

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from matplotlib.patches import Rectangle, Patch
from matplotlib.offsetbox import AnchoredOffsetbox, DrawingArea
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pathlib import Path

# ── paths ──────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).parent.parent
DATA = ROOT / "processed data" / "famine_region_data.csv"
OUT  = Path(__file__).parent / "output" / "figures" / "main"
OUT.mkdir(parents=True, exist_ok=True)

# ── Nature-style global rcParams ───────────────────────────────────────────────
FONT_SIZE   = 14   # tick / legend labels
LABEL_SIZE  = 16   # axis labels
PANEL_SIZE  = 20   # panel letter (A, B, C …)

mpl.rcParams.update({
    "font.family":         "sans-serif",
    "font.sans-serif":     ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":           FONT_SIZE,
    "axes.titlesize":      FONT_SIZE,
    "axes.labelsize":      LABEL_SIZE,
    "xtick.labelsize":     FONT_SIZE,
    "ytick.labelsize":     FONT_SIZE,
    "legend.fontsize":     FONT_SIZE - 1,
    "axes.linewidth":      1.4,
    "xtick.major.width":   1.4,
    "ytick.major.width":   1.4,
    "xtick.major.size":    5,
    "ytick.major.size":    5,
    "xtick.direction":     "out",
    "ytick.direction":     "out",
    "axes.spines.top":     False,
    "axes.spines.right":   False,
    "pdf.fonttype":        42,   # embed fonts
    "ps.fonttype":         42,
})

PANEL_KW = dict(fontsize=PANEL_SIZE, fontweight="bold", loc="left",
                x=-0.08, y=1.05)


def load_data():
    return pd.read_csv(DATA)


# ── Panel A: NINO3.4 timeseries ────────────────────────────────────────────────
def draw_panel_A(ax, data):
    enso = data.groupby("Year")["nino34"].mean().reset_index()

    ax.plot(enso["Year"], enso["nino34"], color="black", linewidth=1.6, alpha=0.8)
    ax.axhline(0, color="black", linestyle="--", linewidth=1.0)

    spans = [
        (1590, 1600, "red",    0.25, "1590s/1690s Super Famines"),
        (1690, 1700, "red",    0.25, "_"),
        (1635, 1637, "blue",   0.90, "Famine in Central Europe"),
        (1648, 1652, "gray",   0.45, "Europe-wide Famine"),
        (1788, 1793, "orange", 0.30, "1788–1793 El Niño"),
    ]
    for x0, x1, col, alpha, label in spans:
        ax.axvspan(x0, x1, color=col, alpha=alpha, label=label)

    top7 = enso.nlargest(7, "nino34")
    top7 = top7[top7["nino34"] != top7.nlargest(4, "nino34").iloc[-1]["nino34"]]
    ax.scatter(top7["Year"], top7["nino34"],
               color="red", edgecolor="black", s=70, zorder=5)

    ax.set_xlabel("Year")
    ax.set_ylabel("NINO3.4 (°C)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22),
              ncol=2, frameon=False, fontsize=FONT_SIZE - 1)


# ── Panel B: Famine chronology ─────────────────────────────────────────────────
def draw_panel_B(ax_gantt, ax_count, data):
    data = data.sort_values(["Region", "Year"])

    famine_periods = []
    for region, group in data.groupby("Region"):
        group = group.sort_values("Year")
        active_starts = []
        for _, row in group.iterrows():
            if row["Famine_start"] == 1:
                active_starts.append(row["Year"])
            if row["Famine"] == 0 and active_starts:
                for sy in active_starts:
                    famine_periods.append((region, sy, row["Year"]))
                active_starts = []
            if row["Year"] == group["Year"].max() and active_starts:
                for sy in active_starts:
                    famine_periods.append((region, sy, row["Year"] + 1))
                active_starts = []

    famine_df  = pd.DataFrame(famine_periods, columns=["Region", "Start", "End"])
    all_years  = pd.Series(range(data["Year"].min(), data["Year"].max() + 1))
    yearly_counts = (
        data[data["Famine"] == 1]
        .groupby("Year").size()
        .reindex(all_years, fill_value=0)
        .reset_index(name="Num_Famines")
    )
    yearly_counts.columns = ["Year", "Num_Famines"]

    region_order  = sorted(data["Region"].unique(), reverse=True)
    region_to_y   = {r: i for i, r in enumerate(region_order)}
    colors        = sns.color_palette("tab20", len(region_order))
    region_colors = {r: colors[i] for i, r in enumerate(region_order)}

    for _, row in famine_df.iterrows():
        ax_gantt.barh(
            y=region_to_y[row["Region"]],
            width=row["End"] - row["Start"],
            left=row["Start"],
            height=0.65,
            color=region_colors[row["Region"]],
            edgecolor="black", linewidth=0.5, alpha=0.9,
        )

    region_labels = [r.replace("Russia/Ukraine", "Russia/Ukr.") for r in region_order]
    ax_gantt.set_yticks(list(region_to_y.values()))
    ax_gantt.set_yticklabels(region_labels, fontsize=FONT_SIZE - 1)
    ax_gantt.set_xlim(data["Year"].min(), data["Year"].max())
    ax_gantt.tick_params(labelbottom=False)

    ax_count.plot(yearly_counts["Year"], yearly_counts["Num_Famines"],
                  color="black", linewidth=2)
    ax_count.fill_between(yearly_counts["Year"], yearly_counts["Num_Famines"],
                          color="gray", alpha=0.3)
    ax_count.set_ylabel("No. regions\nw. famines")
    ax_count.set_xlabel("Year")
    ax_count.set_xlim(data["Year"].min(), data["Year"].max())
    ax_count.set_yticks(range(0, int(yearly_counts["Num_Famines"].max()) + 1, 2))


# ── Panel C: Geographic bar-map ────────────────────────────────────────────────
def draw_panel_C(ax, data):
    region_stats = data.groupby("Region").agg(
        famine_years=("Famine",       "sum"),
        famine_periods=("Famine_start", "sum"),
    ).reset_index()

    region_centers = {
        "France":           (5,    48),
        "Central Europe":   (18,   50),
        "Italy":            (13,   42),
        "Spain":            (-4,   40),
        "Nordic Countries": (15,   61),
        "Great Britain":    (-0,   53),
        "Ireland":          (-8.25, 55),
        "Russia/Ukraine":   (30,   50),
        "Low Countries":    (8.5,  52),
    }

    max_val = max(region_stats["famine_years"].max(),
                  region_stats["famine_periods"].max())
    scale = 50 / max_val

    ax.coastlines(linewidth=0.8)
    ax.add_feature(cfeature.LAND,  facecolor="lightgray", alpha=0.45)
    ax.add_feature(cfeature.OCEAN, facecolor="lightblue", zorder=0, alpha=0.5)
    ax.set_extent([-10, 35, 36, 70], crs=ccrs.PlateCarree())

    for _, row in region_stats.iterrows():
        region = row["Region"]
        if region not in region_centers:
            continue
        x, y = region_centers[region]
        h1, h2 = row["famine_years"] * scale, row["famine_periods"] * scale
        v1, v2 = int(row["famine_years"]), int(row["famine_periods"])

        da = DrawingArea(45, 60, 0, 0)
        da.add_artist(Rectangle((0,  5), 14, h1, color="firebrick"))
        da.add_artist(Rectangle((20, 5), 14, h2, color="darksalmon"))
        da.add_artist(plt.Text(7,  h1 + 7, str(v1), ha="center",
                               fontsize=15, fontweight="bold"))
        da.add_artist(plt.Text(27, h2 + 7, str(v2), ha="center",
                               fontsize=15, fontweight="bold"))
        ab = AnchoredOffsetbox(
            loc="center", child=da, frameon=False,
            bbox_to_anchor=(x, y),
            bbox_transform=ccrs.PlateCarree()._as_mpl_transform(ax),
            pad=0,
        )
        ax.add_artist(ab)

    bar_handles = [
        Patch(color="firebrick",   label="Famine Years"),
        Patch(color="darksalmon",  label="Famine Periods"),
    ]
    ax.legend(handles=bar_handles, loc="upper right", fontsize=FONT_SIZE - 1,
              title="Indicators", title_fontsize=FONT_SIZE - 1)


# ── Combined Nature-style Figure 1 ─────────────────────────────────────────────
def make_figure1(data):
    # Layout: A top-left, B top-right (tall), C bottom spanning full width
    fig = plt.figure(figsize=(16, 18))
    gs_outer = gridspec.GridSpec(
        2, 2,
        figure=fig,
        height_ratios=[1, 1.2],
        hspace=0.55,
        wspace=0.35,
    )

    # --- Panel A (top-left) ---
    ax_a = fig.add_subplot(gs_outer[0, 0])
    draw_panel_A(ax_a, data)
    ax_a.set_title("A", **PANEL_KW)

    # --- Panel B (top-right): Gantt + count stacked ---
    gs_b = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=gs_outer[0, 1],
        height_ratios=[3, 1.8], hspace=0.08,
    )
    ax_b1 = fig.add_subplot(gs_b[0])
    ax_b2 = fig.add_subplot(gs_b[1], sharex=ax_b1)
    draw_panel_B(ax_b1, ax_b2, data)
    ax_b1.set_title("B", **PANEL_KW)

    # --- Panel C (bottom, spanning both columns) ---
    proj = ccrs.LambertConformal(central_longitude=10, central_latitude=50,
                                 standard_parallels=(45, 55))
    ax_c = fig.add_subplot(gs_outer[1, :], projection=proj)
    draw_panel_C(ax_c, data)
    # Panel label via text on the figure (cartopy axes don't support set_title loc easily)
    ax_c.text(-0.04, 1.02, "C", transform=ax_c.transAxes,
              fontsize=PANEL_SIZE, fontweight="bold", va="bottom", ha="left")

    fig.savefig(OUT / "fig1_combined.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved fig1_combined.pdf")


# ── Individual panels (for flexibility / appendix reuse) ──────────────────────
def save_individual_panels(data):
    # Panel A
    fig, ax = plt.subplots(figsize=(8, 4))
    draw_panel_A(ax, data)
    fig.tight_layout(rect=[0, 0.1, 1, 1])
    fig.savefig(OUT / "fig1A_nino34_timeseries.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved fig1A_nino34_timeseries.pdf")

    # Panel B
    fig_b, gs_b = plt.subplots(2, 1, figsize=(8, 7),
                                gridspec_kw={"height_ratios": [3, 1.8], "hspace": 0.08})
    draw_panel_B(gs_b[0], gs_b[1], data)
    fig_b.tight_layout()
    fig_b.savefig(OUT / "fig1B_famine_timeline.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig_b)
    print("Saved fig1B_famine_timeline.pdf")

    # Panel C
    proj = ccrs.LambertConformal(central_longitude=10, central_latitude=50,
                                 standard_parallels=(45, 55))
    fig_c = plt.figure(figsize=(8, 8))
    ax_c  = fig_c.add_subplot(1, 1, 1, projection=proj)
    draw_panel_C(ax_c, data)
    fig_c.tight_layout()
    fig_c.savefig(OUT / "fig1C_famine_map.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig_c)
    print("Saved fig1C_famine_map.pdf")


if __name__ == "__main__":
    data = load_data()
    print("Building Figure 1 (combined) …")
    make_figure1(data)
    print("Saving individual panels …")
    save_individual_panels(data)
    print(f"Done. Output → {OUT}")
