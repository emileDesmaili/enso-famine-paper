"""
assemble_figures.py  –  produce all 5 main Nature-style combined figures.

Reads:
  • CSVs in  analysis/output/data/        (from R scripts 02–05)
  • Raw data in  processed data/           (for Figs 1, 2C-E, 3A-B, 4C-D)

Outputs:
  analysis/output/figures/main/
    fig1_combined.pdf   – Fig 1 (A: NINO3.4 series, B: famine Gantt+count, C: geo bar-map)
    fig2_combined.pdf   – Fig 2 (A: boxplot, B: LPM coefs, C: ML chronology,
                                  D: ML importances, E: ML skill)
    fig3_combined.pdf   – Fig 3 (A: PDSI teleco map, B: precip teleco map,
                                  C: WR harvest IRF, D: other grain IRF)
    fig4_combined.pdf   – Fig 4 (A: Cox HRs, B: survival curves,
                                  C: ML concordance, D: ML importances 2×2)
    fig5_combined.pdf   – Fig 5 (A: grain price IRF, B: fish price IRF)

Run after all R scripts:
  Rscript analysis/run_all.R
  python  analysis/assemble_figures.py
"""

from __future__ import annotations
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import pickle
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import matplotlib.font_manager as fm
import seaborn as sns

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.patches import Rectangle, Patch
from matplotlib.offsetbox import AnchoredOffsetbox, DrawingArea

# ── paths ──────────────────────────────────────────────────────────────────────
ROOT      = Path(__file__).parent.parent
DATA_PROC = ROOT / "processed data"
OUT_DATA  = Path(__file__).parent / "output" / "data"
OUT_MAIN  = Path(__file__).parent / "output" / "figures" / "main"
OUT_MAIN.mkdir(parents=True, exist_ok=True)

# ── Nature-style global rcParams ───────────────────────────────────────────────
FS   = 22   # tick / body
LAB  = 24   # axis label
PAN  = 42   # panel letter

mpl.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":         FS,
    "axes.titlesize":    FS,
    "axes.labelsize":    LAB,
    "xtick.labelsize":   FS,
    "ytick.labelsize":   FS,
    "legend.fontsize":   FS,
    "legend.title_fontsize": FS,
    "axes.linewidth":    1.4,
    "xtick.major.width": 1.4,
    "ytick.major.width": 1.4,
    "xtick.major.size":  5,
    "ytick.major.size":  5,
    "xtick.direction":   "out",
    "ytick.direction":   "out",
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "pdf.fonttype":      42,
    "ps.fonttype":       42,
})

# FontProperties for panel letters
_PAN_FP = fm.FontProperties(size=PAN, weight="bold")


def _polish(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for sp in ("bottom", "left"):
        ax.spines[sp].set_linewidth(1.4)
    ax.tick_params(axis="both", which="major", length=5, width=1.4,
                   direction="out")


def _label(ax, letter, x=-0.16, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes,
            font_properties=_PAN_FP, va="bottom", ha="left")


def _title(ax, text):
    """Short bold panel title, placed above the axes."""
    ax.set_title(text, fontsize=FS, fontweight="bold", loc="left", pad=6)


def _pct_fmt(ax, axis="y", decimals=1):
    fmt = f"{{:.{decimals}%}}"
    fn  = lambda v, _: fmt.format(v)
    if axis == "y":
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(fn))
    else:
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(fn))


# colour / linestyle palettes shared across figures
IRF_COL  = {"All regions":        "black",
             "Teleconnected":      "firebrick",
             "Weakly affected":    "cornflowerblue",
             "Teleco w. controls": "orange",
             "All w. controls":    "orange",
             "All":                "black",
             "Herring":            "#d95f02",
             "Cod":                "steelblue"}
IRF_FILL = {**IRF_COL,
             "All regions":  "cornflowerblue",
             "All":          "cornflowerblue"}
IRF_LS   = {"All regions":        "-",
             "Teleconnected":      "--",
             "Weakly affected":    "--",
             "Teleco w. controls": "-",
             "All w. controls":    "-",
             "All":                "-",
             "Herring":            "--",
             "Cod":                "--"}


# ══════════════════════════════════════════════════════════════════════════════
# ① IRF panel helper (used by Figs 3, 5)
# ══════════════════════════════════════════════════════════════════════════════
def _draw_irf(ax, df, group_col, y_label,
              ribbon_group=None, ylim=None,
              col_map=IRF_COL, fill_map=IRF_FILL, ls_map=IRF_LS,
              legend_kw=None):
    """Generic IRF ribbon+line plot from a tidy dataframe."""
    groups      = df[group_col].unique()
    has_nloc    = "n_loc" in df.columns
    ribbon_set  = {ribbon_group} if isinstance(ribbon_group, str) else set(ribbon_group or [])
    for g in groups:
        sub  = df[df[group_col] == g].sort_values("horizon")
        col  = col_map.get(g, "black")
        fill = fill_map.get(g, col)
        ls   = ls_map.get(g, "-")
        if has_nloc and g in ("Teleconnected", "Weakly affected"):
            leg_label = f"{g} (N={int(sub['n_loc'].iloc[0])})"
        else:
            leg_label = g
        if g in ribbon_set:
            ax.fill_between(sub["horizon"], sub["irf_down"], sub["irf_up"],
                            color=fill, alpha=0.20)
        ax.plot(sub["horizon"], sub["irf_mean"], color=col, lw=2.4,
                linestyle=ls, label=leg_label)
    ax.axhline(0, color="gray", lw=0.8, linestyle="--")
    ax.set_xlabel("Horizon (years)")
    ax.set_ylabel(y_label)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    _pct_fmt(ax)
    if ylim:
        ax.set_ylim(ylim)
    kw = dict(loc="upper right", frameon=False, ncol=1, fontsize=FS)
    if legend_kw:
        kw.update(legend_kw)
    ax.legend(**kw)
    _polish(ax)


# ══════════════════════════════════════════════════════════════════════════════
# FIG 1 – ENSO timeseries / famine Gantt+count / geographic bar-map
# ══════════════════════════════════════════════════════════════════════════════
def _draw_1A(ax, data):
    enso = data.groupby("Year")["nino34"].mean().reset_index()
    ax.plot(enso["Year"], enso["nino34"], color="black", lw=1.6, alpha=0.8)
    ax.axhline(0, color="black", lw=1.0, linestyle="--")

    spans = [
        (1590, 1600, "red",    0.25, "1590s/1690s Catastrophic Famines"),
        (1690, 1700, "red",    0.25, "_"),
        (1635, 1637, "blue",   0.85, "Famine in Central Europe"),
        (1648, 1652, "gray",   0.45, "Europe-wide Famine"),
        (1788, 1793, "orange", 0.30, "1788–1793 El Niño"),
    ]
    for x0, x1, c, a, lbl in spans:
        ax.axvspan(x0, x1, color=c, alpha=a, label=lbl)

    top7 = enso.nlargest(7, "nino34").iloc[:6]
    ax.scatter(top7["Year"], top7["nino34"],
               color="red", edgecolor="black", s=60, zorder=5)

    ax.set_xlabel("Year")
    ax.set_ylabel("NINO3.4 (°C)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22),
              ncol=2, frameon=False, fontsize=FS)
    _title(ax, "ENSO reconstruction, 1500–1800")
    _polish(ax)


def _build_region_palette(data):
    """Return ordered list of regions and their color dict (shared by B and C)."""
    regions = sorted(data["Region"].unique(), reverse=True)
    colors  = sns.color_palette("tab20", len(regions))
    return regions, {r: colors[i] for i, r in enumerate(regions)}


def _draw_1B(ax_gantt, ax_count, data, regions, region_colors):
    data = data.sort_values(["Region", "Year"])
    famine_periods = []
    for region, grp in data.groupby("Region"):
        grp = grp.sort_values("Year")
        starts = []
        for _, row in grp.iterrows():
            if row["Famine_start"] == 1:
                starts.append(row["Year"])
            if row["Famine"] == 0 and starts:
                for sy in starts:
                    famine_periods.append((region, sy, row["Year"]))
                starts = []
            if row["Year"] == grp["Year"].max() and starts:
                for sy in starts:
                    famine_periods.append((region, sy, row["Year"] + 1))
                starts = []

    fdf = pd.DataFrame(famine_periods, columns=["Region", "Start", "End"])
    all_yrs = pd.Series(range(int(data["Year"].min()), int(data["Year"].max()) + 1))
    yearly  = (data[data["Famine"] == 1].groupby("Year").size()
               .reindex(all_yrs, fill_value=0)
               .reset_index(name="N"))
    yearly.columns = ["Year", "N"]

    r_to_y = {r: i for i, r in enumerate(regions)}

    for _, row in fdf.iterrows():
        ax_gantt.barh(y=r_to_y[row["Region"]],
                      width=row["End"] - row["Start"],
                      left=row["Start"], height=0.65,
                      color=region_colors[row["Region"]],
                      edgecolor="black", linewidth=0.5, alpha=0.9)

    ax_gantt.set_yticks(list(r_to_y.values()))
    ax_gantt.set_yticklabels(
        [r.replace("Russia/Ukraine", "Russia/Ukr.") for r in regions],
        fontsize=FS)
    ax_gantt.set_xlim(data["Year"].min(), data["Year"].max())
    ax_gantt.tick_params(labelbottom=False)
    _title(ax_gantt, "Famine chronology by region")
    _polish(ax_gantt)

    ax_count.fill_between(yearly["Year"], yearly["N"], color="gray", alpha=0.3)
    ax_count.plot(yearly["Year"], yearly["N"], color="black", lw=1.8)
    ax_count.set_ylabel("No. regions\nw. famines")
    ax_count.set_xlabel("Year")
    ax_count.set_xlim(data["Year"].min(), data["Year"].max())
    _polish(ax_count)


def _draw_1C(ax, data, region_colors):
    region_stats = data.groupby("Region").agg(
        famine_years  =("Famine",       "sum"),
        famine_periods=("Famine_start", "sum"),
    ).reset_index()

    region_centers = {
        "France":           (5,     48),
        "Central Europe":   (18,    50),
        "Italy":            (13,    42),
        "Spain":            (-4,    40),
        "Nordic Countries": (15,    61),
        "Great Britain":    (0,     53),
        "Ireland":          (-8.25, 55),
        "Russia/Ukraine":   (30,    50),
        "Low Countries":    (8.5,   52),
    }
    max_val = max(region_stats["famine_years"].max(),
                  region_stats["famine_periods"].max())
    scale   = 50 / max_val

    ax.coastlines(linewidth=0.8)
    ax.add_feature(cfeature.LAND,  facecolor="lightgray", alpha=0.45)
    ax.add_feature(cfeature.OCEAN, facecolor="lightblue", zorder=0, alpha=0.5)
    ax.set_extent([-10, 35, 36, 70], crs=ccrs.PlateCarree())

    for _, row in region_stats.iterrows():
        r = row["Region"]
        if r not in region_centers:
            continue
        col = region_colors.get(r, "steelblue")
        # Lighter tint for famine periods bar (same hue, higher lightness)
        col_light = tuple(min(1.0, c + 0.35) for c in col) if isinstance(col, tuple) \
                    else col
        x, y = region_centers[r]
        h1 = row["famine_years"]   * scale
        h2 = row["famine_periods"] * scale
        v1, v2 = int(row["famine_years"]), int(row["famine_periods"])
        da = DrawingArea(45, 60, 0, 0)
        da.add_artist(Rectangle((0,  5), 14, h1, color=col))
        da.add_artist(Rectangle((20, 5), 14, h2, color=col_light))
        da.add_artist(plt.Text(7,  h1 + 7, str(v1), ha="center",
                               fontsize=FS, fontweight="bold"))
        da.add_artist(plt.Text(27, h2 + 7, str(v2), ha="center",
                               fontsize=FS, fontweight="bold"))
        ab = AnchoredOffsetbox(
            loc="center", child=da, frameon=False,
            bbox_to_anchor=(x, y),
            bbox_transform=ccrs.PlateCarree()._as_mpl_transform(ax),
            pad=0,
        )
        ax.add_artist(ab)

    ax.legend(handles=[Patch(color="dimgray",  label="Famine Years (left bar)"),
                        Patch(color="lightgray", label="Famine Periods (right bar)")],
              loc="upper left", fontsize=FS, frameon=False,
              prop={"weight": "bold"})


def make_fig1():
    data = pd.read_csv(DATA_PROC / "famine_region_data.csv")
    regions, region_colors = _build_region_palette(data)

    fig = plt.figure(figsize=(18, 22))
    # Row 0: A full width
    # Row 1: B (left, narrower) + C (right, wider map)
    gs = gridspec.GridSpec(2, 2, figure=fig,
                           height_ratios=[1, 1.5],
                           width_ratios=[0.38, 0.62],
                           hspace=0.52, wspace=0.38)

    # A – full-width ENSO timeseries
    ax_a = fig.add_subplot(gs[0, :])
    _draw_1A(ax_a, data)
    _label(ax_a, "a")

    # B – Gantt + count stacked, bottom-left
    gs_b = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=gs[1, 0],
        height_ratios=[3, 1.6], hspace=0.08,
    )
    ax_b1 = fig.add_subplot(gs_b[0])
    ax_b2 = fig.add_subplot(gs_b[1], sharex=ax_b1)
    _draw_1B(ax_b1, ax_b2, data, regions, region_colors)
    _label(ax_b1, "b")

    # C – geo bar-map, bottom-right
    proj = ccrs.LambertConformal(central_longitude=10, central_latitude=50,
                                 standard_parallels=(45, 55))
    ax_c = fig.add_subplot(gs[1, 1], projection=proj)
    _draw_1C(ax_c, data, region_colors)
    _label(ax_c, "c", x=-0.10)
    ax_c.set_title("Famine incidence by region", fontsize=FS, fontweight="bold",
                   loc="left", pad=6)

    fig.savefig(OUT_MAIN / "fig1_combined.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved fig1_combined.pdf")


# ══════════════════════════════════════════════════════════════════════════════
# FIG 2 – Famine onset (A+B from R CSVs, C+D+E from ML)
# ══════════════════════════════════════════════════════════════════════════════
def _draw_2A(ax):
    from scipy.stats import ttest_ind
    df = pd.read_csv(OUT_DATA / "fig2A_onset_box.csv")
    colors = {"Famine Onset": "darksalmon", "No Famine": "lightblue"}
    grp_data = {}
    for i, g in enumerate(["No Famine", "Famine Onset"]):
        sub = df[df["Group"] == g]["nino34"]
        grp_data[g] = sub
        bp  = ax.boxplot(sub, positions=[i], widths=0.45,
                         patch_artist=True,
                         boxprops=dict(facecolor=colors[g], linewidth=1.3),
                         medianprops=dict(color="black", linewidth=2),
                         whiskerprops=dict(linewidth=1.3),
                         capprops=dict(linewidth=1.3),
                         flierprops=dict(marker="o", markersize=0,
                                         markerfacecolor="black"))
        ax.scatter(np.random.normal(i, 0.07, len(sub)), sub,
                   color="black", alpha=0.45, s=18, zorder=3)

    # Welch t-test (one-sided: famine onset > no famine)
    _, pval_two = ttest_ind(grp_data["Famine Onset"], grp_data["No Famine"],
                            equal_var=False)
    pval = pval_two / 2
    y_top = max(grp_data["Famine Onset"].max(), grp_data["No Famine"].max())
    y_range = y_top - min(grp_data["Famine Onset"].min(), grp_data["No Famine"].min())
    y_bar = y_top + 0.08 * y_range
    tick_h = 0.03 * y_range
    ax.plot([0, 0, 1, 1], [y_bar - tick_h, y_bar, y_bar, y_bar - tick_h],
            color="black", lw=1.5)
    stars = "***" if pval < 0.01 else ("**" if pval < 0.05 else
            ("*" if pval < 0.10 else f"p={pval:.2f}"))
    ax.text(0.5, y_bar + 0.01 * y_range, stars, ha="center", va="bottom",
            fontsize=FS + 4, fontweight="bold", color="black")
    ax.set_ylim(top=y_bar + 0.18 * y_range)

    ax.axhline(0, color="black", lw=1.0, linestyle="--")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["No Famine", "Famine Onset"], fontsize=FS + 2)
    ax.set_ylabel("NINO3.4 (°C)", fontsize=LAB + 2)
    ax.tick_params(axis="both", labelsize=FS + 2)
    _title(ax, "ENSO at famine onset\n(Central Europe)")
    _polish(ax)


def _draw_2B(ax):
    df = pd.read_csv(OUT_DATA / "fig2B_onset_coef.csv")
    model_order  = ["No FEs", "FEs", "FEs + Controls"]
    model_colors = {"No FEs": "cornflowerblue", "FEs": "orange",
                    "FEs + Controls": "firebrick"}
    regions = df["Region"].unique()
    x       = np.arange(len(regions))
    width   = 0.18
    offsets = np.linspace(-(len(model_order) - 1) / 2,
                          (len(model_order) - 1) / 2,
                          len(model_order)) * width

    for i, model in enumerate(model_order):
        sub = df[df["model"] == model].set_index("Region").reindex(regions)
        ax.errorbar(x + offsets[i], sub["estimate"],
                    yerr=[sub["estimate"] - sub["conf.low"],
                          sub["conf.high"] - sub["estimate"]],
                    fmt="o", color=model_colors[model],
                    elinewidth=1.4, capsize=0, ms=6,
                    label=model, zorder=3)

    ax.axhline(0, color="black", lw=0.9, linestyle="--")
    ax.set_xticks(x)
    ax.set_xticklabels(regions, rotation=30, ha="right", fontsize=FS + 1)
    ax.set_ylabel("Famine Onset Probability", fontsize=LAB + 2)
    ax.tick_params(axis="y", labelsize=FS + 1)
    ax.legend(loc="best", frameon=False, fontsize=FS + 1, ncol=3)
    _title(ax, "Probability of Famine Onset")
    _polish(ax)


def _draw_2C_chron(ax_top, ax_bot, df_chron, onset_raw):
    """Chronology panel: bars use actual famine duration (Famine_dur), matching notebook."""
    years      = sorted(df_chron["Year"].unique())
    ymin, ymax = min(years), max(years)
    colors     = {"Observed": "firebrick", "Predicted": "darksalmon",
                  "Counterfactual": "cornflowerblue"}
    # Build duration lookup: Year → Famine_dur (from raw onset data)
    dur_lookup = (onset_raw[onset_raw["Famine_start"] == 1]
                  .set_index("Year")["Famine_dur"].to_dict()
                  if "Famine_dur" in onset_raw.columns else {})

    y_pos = {}
    for i, reg in enumerate(sorted(df_chron["Region"].unique())):
        base  = i * 6
        y_pos[reg] = {"Observed": base + 3,
                      "Predicted": base + 1.5,
                      "Counterfactual": base}
        rd = df_chron[df_chron["Region"] == reg].sort_values("Year")
        for sc in ["Observed", "Predicted", "Counterfactual"]:
            for yr in rd.loc[rd[sc] == 1, "Year"]:
                dur = dur_lookup.get(yr, 1)
                ax_top.barh(y=y_pos[reg][sc], width=dur, left=yr, height=1.2,
                            color=colors[sc], edgecolor="black",
                            linewidth=0.9, alpha=0.8)

    ticks, labels = [], []
    for reg in sorted(df_chron["Region"].unique()):
        for sc in ["Counterfactual", "Predicted", "Observed"]:
            n = df_chron[df_chron["Region"] == reg][sc].sum()
            ticks.append(y_pos[reg][sc])
            labels.append(f"{sc} ({n})")
    ax_top.set_yticks(ticks)
    ax_top.set_yticklabels(labels, fontsize=FS - 1)
    ax_top.set_xlim(ymin - 1, ymax + 1)
    ax_top.tick_params(labelbottom=False)
    _title(ax_top, "Predicted vs. observed famine onsets")
    _polish(ax_top)

    nino_a = df_chron.groupby("Year")["NINO34_actual"].first().reindex(years, fill_value=0).values
    nino_c = df_chron.groupby("Year")["NINO34_counterfactual"].first().reindex(years, fill_value=0).values
    ax_bot.plot(years, nino_a, color="crimson",    lw=1.4, label="Observed",       alpha=0.85)
    ax_bot.plot(years, nino_c, color="dodgerblue", lw=1.4, label="Counterfactual", alpha=0.9)
    ax_bot.axhline(0, color="black", lw=0.8, linestyle="-", alpha=0.3)
    ax_bot.set_xlim(ymin - 1, ymax + 1)
    ax_bot.set_xlabel("Year")
    ax_bot.set_ylabel("NINO3.4 (°C)")
    leg = ax_bot.legend(loc="lower right", frameon=True, fontsize=FS,
                        framealpha=0.9, edgecolor="gray",
                        handlelength=2.0, labelspacing=0.4)
    leg.get_frame().set_linewidth(0.8)
    _polish(ax_bot)


def _draw_2D(ax, df_imp):
    med_order = df_imp.groupby("Feature")["Importance"].median().sort_values(ascending=False).index
    sns.boxplot(x="Importance", y="Feature", data=df_imp,
                order=med_order, palette="Reds_r", ax=ax)
    ax.set_xlabel("Permutation Importance", fontsize=LAB + 2)
    ax.set_ylabel("", fontsize=LAB + 2)
    ax.tick_params(axis="both", labelsize=FS + 2)
    _title(ax, "Feature importances")
    _polish(ax)


def _draw_2E(ax, in_sample_acc, cv_scores):
    skill_df = pd.DataFrame({
        "Skill":    ["In-sample", "Cross-Validation"],
        "Accuracy": [in_sample_acc, cv_scores.mean()],
    })
    sns.barplot(data=skill_df, y="Skill", x="Accuracy",
                palette={"In-sample": "firebrick", "Cross-Validation": "darksalmon"},
                ax=ax)
    ax.errorbar(x=cv_scores.mean(), y=1, xerr=cv_scores.std(),
                color="black", fmt="none", capsize=0, elinewidth=1.5)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Accuracy", fontsize=LAB + 2)
    ax.tick_params(axis="both", labelsize=FS + 2)
    _title(ax, "Prediction skill")
    _polish(ax)


def make_fig2():
    _cache2 = OUT_DATA / "_ml_onset_cache.pkl"
    if _cache2.exists():
        print("  Loading ML onset classifier from cache …")
        with open(_cache2, "rb") as f:
            ml = pickle.load(f)
    else:
        print("  Running ML onset classifier …")
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "ml07", Path(__file__).parent / "07_ml_onset_survival.py")
        ml07 = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(ml07)
        ml = ml07.run_onset_classifier()
        with open(_cache2, "wb") as f:
            pickle.dump(ml, f)

    fig = plt.figure(figsize=(22, 20))
    # Row 0: A (left) + B (right)
    # Row 1: C chronology (left ~54%) + D+E stacked (right ~45%)
    gs = gridspec.GridSpec(2, 1, figure=fig,
                           height_ratios=[1.0, 2.2],
                           hspace=0.46)

    # ── Row 0: A + B ──────────────────────────────────────────────────────────
    gs_top = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0],
                                              width_ratios=[0.35, 0.65],
                                              wspace=0.40)
    ax_a = fig.add_subplot(gs_top[0])
    _draw_2A(ax_a)
    _label(ax_a, "a")

    ax_b = fig.add_subplot(gs_top[1])
    _draw_2B(ax_b)
    _label(ax_b, "b")

    # ── Row 1: C (left 54%) + D/E stacked (right 45%) ─────────────────────────
    gs_bot = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1],
                                              width_ratios=[0.54, 0.45],
                                              wspace=0.42)

    # C – chronology (left half, split top/bottom)
    gs_c = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=gs_bot[0],
        height_ratios=[1.6, 1], hspace=0.08,
    )
    ax_c1 = fig.add_subplot(gs_c[0])
    ax_c2 = fig.add_subplot(gs_c[1])
    _draw_2C_chron(ax_c1, ax_c2, ml["df_chron"], ml["onset"])
    _label(ax_c1, "c")

    # D + E stacked (right half)
    gs_de = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=gs_bot[1],
        height_ratios=[1.4, 0.8], hspace=0.52,
    )
    ax_d = fig.add_subplot(gs_de[0])
    _draw_2D(ax_d, ml["df_imp"])
    _label(ax_d, "d")

    ax_e = fig.add_subplot(gs_de[1])
    _draw_2E(ax_e, ml["in_sample_acc"], ml["cv_scores"])
    _label(ax_e, "e")

    fig.savefig(OUT_MAIN / "fig2_combined.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved fig2_combined.pdf")


_IRF_ORDER = ["All regions", "Weakly affected", "Teleconnected", "Teleco w. controls"]

def _draw_irf_errorbar(ax, df, group_col, y_label,
                       ylim=None,
                       col_map=IRF_COL, ls_map=IRF_LS,
                       add_legend=True):
    """IRF plot using points + whiskers (errorbar) – no ribbon lines."""
    all_groups = df[group_col].unique().tolist()
    groups = [g for g in _IRF_ORDER if g in all_groups] + \
             [g for g in all_groups if g not in _IRF_ORDER]
    n       = len(groups)
    offsets = np.linspace(-0.15, 0.15, n) if n > 1 else [0]
    has_nloc = "n_loc" in df.columns
    for g, offset in zip(groups, offsets):
        sub = df[df[group_col] == g].sort_values("horizon")
        col = col_map.get(g, "black")
        # Add N= to legend for Teleconnected/Weakly affected groups
        if has_nloc and g in ("Teleconnected", "Weakly affected"):
            nloc = int(sub["n_loc"].iloc[0])
            leg_label = f"{g} (N={nloc})"
        else:
            leg_label = g
        ax.errorbar(sub["horizon"] + offset,
                    sub["irf_mean"],
                    yerr=[sub["irf_mean"] - sub["irf_down"],
                          sub["irf_up"]   - sub["irf_mean"]],
                    fmt="o", color=col, elinewidth=2.0, capsize=0, ms=8,
                    label=leg_label, zorder=3)
    ax.axhline(0, color="gray", lw=0.8, linestyle="--")
    ax.set_xlabel("Horizon (years)")
    ax.set_ylabel(y_label)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    _pct_fmt(ax)
    if ylim:
        ax.set_ylim(ylim)
    if add_legend:
        ax.legend(loc="upper right", frameon=False, ncol=1, fontsize=FS)
    _polish(ax)


# ══════════════════════════════════════════════════════════════════════════════
# FIG 3 – Teleconnection maps (A+B) + yield IRFs (C+D)
# ══════════════════════════════════════════════════════════════════════════════
def _draw_teleco_map(ax, corr_vals, p_vals, lat2d, lon2d,
                     title, cmap, alpha=0.10):
    """Draw a single teleconnection map panel (used for A and B)."""
    import numpy as np
    from matplotlib.colors import BoundaryNorm
    levels  = np.linspace(-0.2, 0.2, 17)
    cmap_   = plt.get_cmap(cmap, len(levels) + 1)
    norm    = BoundaryNorm(levels, ncolors=cmap_.N, extend="both")

    im = ax.pcolormesh(lon2d, lat2d, corr_vals,
                       cmap=cmap_, norm=norm, shading="auto",
                       transform=ccrs.PlateCarree())
    sig = p_vals <= alpha
    ax.scatter(lon2d[sig], lat2d[sig], s=25, color="black",
               transform=ccrs.PlateCarree(), zorder=5, alpha=0.7)
    ax.coastlines(linewidth=0.9)
    ax.add_feature(cfeature.BORDERS, linestyle=":", linewidth=0.5)
    ax.add_feature(cfeature.LAND, facecolor="whitesmoke", zorder=0)
    ax.add_feature(cfeature.OCEAN, facecolor="aliceblue", zorder=0)
    ax.set_extent([-10, 35, 36, 70], crs=ccrs.PlateCarree())

    # Gridlines with lat/lon labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.6, color="gray", alpha=0.5, linestyle="--")
    gl.top_labels   = False
    gl.right_labels = False
    gl.xlocator = mticker.FixedLocator([-10, 0, 10, 20, 30])
    gl.ylocator = mticker.FixedLocator([40, 50, 60, 70])
    gl.xlabel_style = {"size": FS - 3}
    gl.ylabel_style = {"size": FS - 3}

    cb = plt.colorbar(im, ax=ax, orientation="vertical",
                      fraction=0.046, pad=0.06, label="Pearson r",
                      extend="both")
    cb.ax.tick_params(labelsize=FS - 3)
    cb.set_label("Pearson r", fontsize=FS - 2)
    ax.set_title(title, fontsize=FS, fontweight="bold", loc="left")

    # Central Europe region outline
    import cartopy.io.shapereader as shpreader
    import geopandas as gpd
    from shapely.ops import unary_union
    _CE_COUNTRIES = {
        "Switzerland", "Germany", "Austria", "Czechia",
        "Hungary", "Slovenia", "Bosnia and Herz.",
        "Slovakia", "Croatia", "Serbia", "Poland",
    }
    _ne_path = shpreader.natural_earth(resolution="10m", category="cultural",
                                       name="admin_0_countries")
    _gdf = gpd.read_file(_ne_path)
    _ce_geom = unary_union(_gdf[_gdf["NAME"].isin(_CE_COUNTRIES)]["geometry"].values)
    ax.add_geometries(
        [_ce_geom], crs=ccrs.PlateCarree(),
        facecolor="none", edgecolor="red", linewidth=4.0, zorder=10
    )

    # Yield record locations – purple = teleconnected (PDSI p<0.10), white = not
    import pandas as _pd
    _yield = _pd.read_csv(DATA_PROC / "yield_ljungqvist_2025.csv")
    _ylocs = (
        _yield.groupby(["lat", "lon"])
        .agg(teleco=("teleco_PDSI_10", "max"))
        .reset_index()
    )
    _tele = _ylocs[_ylocs["teleco"] == 1]
    _non  = _ylocs[_ylocs["teleco"] == 0]
    ax.scatter(_non["lon"],  _non["lat"],  s=180, color="white",
               edgecolors="#333333", linewidths=1.1,
               transform=ccrs.PlateCarree(), zorder=11, alpha=0.9)
    ax.scatter(_tele["lon"], _tele["lat"], s=200, color="#7b2d8b",
               edgecolors="white", linewidths=1.1,
               transform=ccrs.PlateCarree(), zorder=12)


def _compute_teleconnections():
    """
    Load ENSO, OWDA, ModE-RA and compute correlation maps.
    Returns (corr_pdsi, p_pdsi, corr_prec, p_prec, lat2d, lon2d_p,
             lat2d_p, lon2d_prec).
    """
    import xarray as xr
    import cftime
    from scipy.stats import pearsonr

    YEAR_RANGE = slice(1500, 1800)

    # Cook 2024 ENSO reconstruction
    enso3   = pd.read_csv(ROOT / "data" / "cook2024-R15-ENSO-Rec-1500-2000.txt",
                          delimiter="\t", comment="#", na_values="NA")
    lat_lon = pd.read_csv(ROOT / "data" / "cook2024-ENSO-latlon.txt",
                          delimiter="\t", comment="#", na_values="NA")
    em      = enso3.melt(id_vars=["Year"], var_name="gridpoint", value_name="enso")
    em["gridpoint"] = em["gridpoint"].astype(int)
    em      = em.merge(lat_lon, on="gridpoint")
    enso_xr = em.set_index(["Year", "lat", "lon"])["enso"].to_xarray()
    enso_xr = enso_xr - enso_xr.rolling(Year=50, min_periods=1).mean()
    nino34  = enso_xr.sel(lat=slice(-5, 5), lon=slice(-170, -120)).mean(dim=["lat", "lon"])
    nino_df = nino34.to_dataframe().reset_index().rename(columns={"enso": "nino34"})
    nino_df = nino_df[nino_df["Year"] <= 1800]

    def _corr_grid(field_xr, nino_series):
        lat = field_xr[list(field_xr.dims)[1]].values
        lon = field_xr[list(field_xr.dims)[2]].values
        lon2d, lat2d = np.meshgrid(lon, lat)
        nv    = nino_series.values
        corrs = np.full((len(lat), len(lon)), np.nan)
        pvals = np.full_like(corrs, np.nan)
        for i in range(len(lat)):
            for j in range(len(lon)):
                fts  = field_xr[:, i, j].values
                mask = np.isfinite(fts) & np.isfinite(nv)
                if mask.sum() > 10:
                    r, p = pearsonr(fts[mask], nv[mask])
                    corrs[i, j] = r
                    pvals[i, j] = p
        return corrs, pvals, lat2d, lon2d

    DATA_RAW = ROOT / "data"

    # Panel A – PDSI (OWDA)
    owda_path = DATA_RAW / "owda.nc"
    if owda_path.exists():
        owda = xr.open_dataset(owda_path)
        owda["PDSI_w"] = (np.sqrt(np.cos(np.deg2rad(owda["lat"])) + 1e-6)
                          * owda["pdsi"])
        owda = (owda.assign_coords(Year=owda.time)
                .swap_dims({"time": "Year"})
                .sel(Year=YEAR_RANGE)
                .sel(lat=slice(35, 70), lon=slice(-15, 35))
                .rename({"lat": "latitude", "lon": "longitude"}))
        nino_xr = xr.Dataset.from_dataframe(nino_df.set_index("Year"))
        merged  = xr.merge([owda, nino_xr])
        pdsi    = merged["PDSI_w"].transpose("Year", "latitude", "longitude")
        nino_s  = merged["nino34"]
        corr_pdsi, p_pdsi, lat2d_p, lon2d_p = _corr_grid(pdsi, nino_s)
    else:
        corr_pdsi = p_pdsi = lat2d_p = lon2d_p = None
        print("  owda.nc not found – skipping PDSI teleconnection")

    # Panel B – Summer precipitation (ModE-RA)
    prec_path = DATA_RAW / "ModE-RA_ensmean_totprec_anom_wrt_1901-2000_1421-2008_mon.nc"
    if prec_path.exists():
        def shift_dec(t_arr):
            return [cftime.DatetimeGregorian(t.year + 1, t.month, t.day)
                    if t.month in [11, 12] else t for t in t_arr]

        xds_p = xr.open_dataset(prec_path, use_cftime=True)
        xds_p = xds_p.sel(time=slice("1500-01-01", "1800-12-31"))
        xds_p = xds_p.assign_coords(time=("time", shift_dec(xds_p.time.values)))
        xds_p["precip"] = (np.sqrt(np.cos(np.deg2rad(xds_p["latitude"])) + 1e-6)
                           * xds_p["totprec"] * 86400 * 30)
        xds_eu = xds_p.sel(latitude=slice(70, 35), longitude=slice(-15, 35))
        xds_eu_amjj = xds_eu.sel(time=xds_eu.time.dt.month.isin([4, 5, 6, 7]))
        summer = xds_eu_amjj.groupby(xds_eu_amjj.time.dt.year).mean(dim="time")
        # rebuild Year coordinate (dim name varies across xarray versions)
        year_dim = [d for d in summer.dims if d not in ("latitude", "longitude")][0]
        summer = summer.rename({year_dim: "Year"})
        nino_xr2 = xr.Dataset.from_dataframe(nino_df.set_index("Year"))
        merged2  = xr.merge([summer, nino_xr2])
        precip   = merged2["precip"].transpose("Year", "latitude", "longitude")
        nino_s2  = merged2["nino34"].sel(Year=YEAR_RANGE)
        precip   = precip.sel(Year=YEAR_RANGE)
        corr_prec, p_prec, lat2d_pr, lon2d_pr = _corr_grid(precip, nino_s2)
    else:
        corr_prec = p_prec = lat2d_pr = lon2d_pr = None
        print("  ModE-RA precip file not found – skipping precip teleconnection")

    return (corr_pdsi, p_pdsi, lat2d_p, lon2d_p,
            corr_prec, p_prec, lat2d_pr, lon2d_pr)


def _draw_pdsi_composite_map(ax):
    """
    Panel B replacement: PDSI composite anomaly for Counterfactual-ML famine years,
    replicating mechanisms.ipynb Fig 1a with black stippling (p < 0.05).
    """
    import xarray as xr
    import cftime
    from scipy.stats import ttest_ind
    from matplotlib.colors import BoundaryNorm

    DATA_RAW  = ROOT / "data"
    DATA_PROC = ROOT / "processed data"

    # Load OWDA
    _coder  = xr.coders.CFDatetimeCoder(use_cftime=True)
    owda    = xr.open_dataset(DATA_RAW / "owda.nc", decode_times=_coder)
    owda    = owda.assign_coords(Year=owda.time).swap_dims({"time": "Year"})
    owda    = owda.sel(Year=slice(1500, 1800))
    owda_eu = owda.sel(lat=slice(35, 70), lon=slice(-15, 35))

    # Load famine-year list (Counterfactual ML: predicted=1, counterfactual=0)
    chron = pd.read_csv(OUT_DATA / "fig2C_chronology_onsets.csv")
    ef_years = sorted(
        chron[(chron["Predicted_ML"] == 1) & (chron["Counterfactual_ML"] == 0)]["Year"].tolist()
    )

    # Compute composite anomaly
    da       = owda_eu["pdsi"].transpose("Year", "lat", "lon")
    all_yrs  = da["Year"].values.astype(int)
    matched  = [y for y in ef_years if y in all_yrs]
    baseline = [y for y in all_yrs  if y not in ef_years]

    ef_3d  = da.sel(Year=matched).values
    bas_3d = da.sel(Year=baseline).values
    anom   = ef_3d.mean(axis=0) - bas_3d.mean(axis=0)

    nlat, nlon = ef_3d.shape[1], ef_3d.shape[2]
    pvals = np.full((nlat, nlon), np.nan)
    for i in range(nlat):
        for j in range(nlon):
            e, b = ef_3d[:, i, j], bas_3d[:, i, j]
            me, mb = np.isfinite(e), np.isfinite(b)
            if me.sum() >= 3 and mb.sum() >= 3:
                _, pvals[i, j] = ttest_ind(e[me], b[mb], equal_var=False)

    lon    = owda_eu["lon"].values
    lat    = owda_eu["lat"].values
    lon2d, lat2d = np.meshgrid(lon, lat)

    levels = np.linspace(-1.0, 1.0, 11)
    cmap_  = plt.get_cmap("BrBG", len(levels) + 1)
    norm   = BoundaryNorm(levels, ncolors=cmap_.N, extend="both")

    im = ax.pcolormesh(lon2d, lat2d, anom,
                       transform=ccrs.PlateCarree(),
                       cmap=cmap_, norm=norm, shading="auto")
    sig = pvals <= 0.10
    ax.scatter(lon2d[sig], lat2d[sig], s=25, color="black",
               transform=ccrs.PlateCarree(), zorder=5, alpha=0.7)

    ax.coastlines(linewidth=0.9)
    ax.add_feature(cfeature.BORDERS, linestyle=":", linewidth=0.5)
    ax.add_feature(cfeature.LAND, facecolor="whitesmoke", zorder=0)
    ax.add_feature(cfeature.OCEAN, facecolor="aliceblue", zorder=0)
    ax.set_extent([-12, 36, 35, 71], crs=ccrs.PlateCarree())

    # Central Europe region outline
    import cartopy.io.shapereader as shpreader
    import geopandas as gpd
    from shapely.ops import unary_union
    _CE_COUNTRIES = {
        "Switzerland", "Germany", "Austria", "Czechia",
        "Hungary", "Slovenia", "Bosnia and Herz.",
        "Slovakia", "Croatia", "Serbia", "Poland",
    }
    _ne_path = shpreader.natural_earth(resolution="10m", category="cultural",
                                       name="admin_0_countries")
    _gdf = gpd.read_file(_ne_path)
    _ce_geom = unary_union(_gdf[_gdf["NAME"].isin(_CE_COUNTRIES)]["geometry"].values)
    ax.add_geometries(
        [_ce_geom], crs=ccrs.PlateCarree(),
        facecolor="none", edgecolor="red", linewidth=4.0, zorder=10
    )

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.6, color="gray", alpha=0.5, linestyle="--")
    gl.top_labels   = False
    gl.right_labels = False
    gl.xlocator = mticker.FixedLocator([-10, 0, 10, 20, 30])
    gl.ylocator = mticker.FixedLocator([40, 50, 60, 70])
    gl.xlabel_style = {"size": FS - 3}
    gl.ylabel_style = {"size": FS - 3}

    cb = plt.colorbar(im, ax=ax, orientation="vertical",
                      fraction=0.046, pad=0.06, extend="both")
    cb.ax.tick_params(labelsize=FS - 3)
    cb.set_label(r"$\Delta$scPDSI", fontsize=FS - 2)
    n = len(matched)
    ax.set_title(f"PDSI anomaly – ENSO Famines (n={n})",
                 fontsize=FS, fontweight="bold", loc="left")


def make_fig3():
    tele = _compute_teleconnections()
    corr_pdsi, p_pdsi, lat2d_p, lon2d_p = tele[:4]
    # tele[4:] (precip teleconnection) no longer used for panel B

    df_c = pd.read_csv(OUT_DATA / "fig3C_irf_yield_WR.csv")
    df_d = pd.read_csv(OUT_DATA / "fig3D_irf_yield_noWR.csv")

    proj = ccrs.EuroPP()
    fig  = plt.figure(figsize=(22, 20))
    gs   = gridspec.GridSpec(2, 2, figure=fig,
                             hspace=0.55, wspace=0.36,
                             height_ratios=[1.5, 1.0])

    # A – PDSI map
    ax_a = fig.add_subplot(gs[0, 0], projection=proj)
    if corr_pdsi is not None:
        _draw_teleco_map(ax_a, corr_pdsi, p_pdsi, lat2d_p, lon2d_p,
                         "ENSO / PDSI correlation", "BrBG")
    else:
        ax_a.text(0.5, 0.5, "Data not found", transform=ax_a.transAxes,
                  ha="center", va="center")
    _label(ax_a, "a")

    # B – PDSI composite anomaly map (mechanisms Fig 1a)
    ax_b = fig.add_subplot(gs[0, 1], projection=proj)
    _draw_pdsi_composite_map(ax_b)
    _label(ax_b, "b")

    # C – Wheat/Rye IRF (points + whiskers)
    ax_c = fig.add_subplot(gs[1, 0])
    _draw_irf_errorbar(ax_c, df_c, "label", "% Harvest response",
                       ylim=(-0.13, 0.06), add_legend=False)
    ax_c.legend(loc="upper center", bbox_to_anchor=(0.5, 1.32),
                ncol=2, frameon=False, fontsize=FS)
    _title(ax_c, "Wheat & rye response")
    _label(ax_c, "c", x=-0.16, y=1.42)

    # D – Other grain IRF (points + whiskers)
    ax_d = fig.add_subplot(gs[1, 1])
    _draw_irf_errorbar(ax_d, df_d, "label", "% Harvest response",
                       ylim=(-0.13, 0.06), add_legend=False)
    ax_d.legend(loc="upper center", bbox_to_anchor=(0.5, 1.32),
                ncol=2, frameon=False, fontsize=FS)
    _title(ax_d, "Other grains response")
    _label(ax_d, "d", x=-0.16, y=1.42)

    fig.savefig(OUT_MAIN / "fig3_combined.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved fig3_combined.pdf")


# ══════════════════════════════════════════════════════════════════════════════
# FIG 4 – Famine duration (A+B Cox, C+D ML survival)
# ══════════════════════════════════════════════════════════════════════════════
def _draw_4A(ax):
    df = pd.read_csv(OUT_DATA / "fig4A_cox_coefs.csv")
    df = df.sort_values("coef")
    y  = np.arange(len(df))
    # 95% CI (thin outer whisker)
    ax.errorbar(df["coef"], y,
                xerr=[df["coef"] - df["lower"], df["upper"] - df["coef"]],
                fmt="none", color="firebrick", elinewidth=2.5, capsize=0, zorder=2)
    # 90% CI (thick inner whisker)
    lo90 = df["coef"] - 1.645 * df["se"]
    hi90 = df["coef"] + 1.645 * df["se"]
    ax.errorbar(df["coef"], y,
                xerr=[df["coef"] - lo90, hi90 - df["coef"]],
                fmt="none", color="firebrick", elinewidth=5.0, capsize=0, zorder=3)
    # Dot on top
    ax.scatter(df["coef"], y, color="firebrick", s=120, zorder=4)
    ax.axvline(0, color="black", lw=0.9, linestyle="--")
    ax.set_yticks(y)
    ax.set_yticklabels(df["Model"], fontsize=FS + 1)
    ax.set_xlabel("NINO3.4 log hazard ratio", fontsize=LAB + 2)
    _title(ax, "Cox hazard ratios")
    _polish(ax)


def _draw_4B(ax):
    df = pd.read_csv(OUT_DATA / "fig4B_survcurves.csv")
    enso_order  = ["La Niña (-1°)", "Neutral", "El Niño (+1°)", "Strong El Niño (+2°)"]
    enso_colors = {"La Niña (-1°)": "cornflowerblue",
                   "Neutral":       "black",
                   "El Niño (+1°)": "darksalmon",
                   "Strong El Niño (+2°)": "firebrick"}
    enso_markers = {"La Niña (-1°)": "o", "Neutral": "s",
                    "El Niño (+1°)": "^", "Strong El Niño (+2°)": "D"}

    for label in enso_order:
        sub = df[df["ENSO"] == label].sort_values("time")
        ax.plot(sub["time"], sub["surv"],
                color=enso_colors[label],
                marker=enso_markers[label], ms=5, lw=1.6, label=label)

    # Delta annotation: El Niño (+2°) − La Niña (−1°) at year 3 (matches Rmd)
    year_of_interest = 3
    surv_elni = (df[df["ENSO"] == "Strong El Niño (+2°)"]
                 .sort_values("time")
                 .set_index("time")["surv"])
    surv_lani = (df[df["ENSO"] == "La Niña (-1°)"]
                 .sort_values("time")
                 .set_index("time")["surv"])
    try:
        idx_en = surv_elni.index.get_indexer([year_of_interest], method="nearest")[0]
        idx_ln = surv_lani.index.get_indexer([year_of_interest], method="nearest")[0]
        s_en   = surv_elni.iloc[idx_en]
        s_ln   = surv_lani.iloc[idx_ln]
        delta_pp = (s_en - s_ln) * 100
        ax.axvline(year_of_interest, color="black", lw=1.0, linestyle="--")
        delta_text = (f"$\\Delta$ El Niño $-$ La Niña = {round(delta_pp):+d} p.p.")
        ax.annotate(
            delta_text,
            xy=(8, 0.40),
            xycoords="data",
            ha="left", va="center", fontsize=FS,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.85,
                      linewidth=0.5),
        )
    except Exception:
        pass

    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:.0%}"))
    ax.set_xlabel("Famine duration (years)")
    ax.set_ylabel("Survival probability")
    ax.legend(loc="best", frameon=False, fontsize=FS, ncol=2)
    _title(ax, "Duration under ENSO scenarios")
    _polish(ax)


def make_fig4():
    _cache4 = OUT_DATA / "_ml_survival_cache.pkl"
    if _cache4.exists():
        print("  Loading ML survival models from cache …")
        with open(_cache4, "rb") as f:
            methods, means, errors, perm_dfs = pickle.load(f)
    else:
        print("  Running ML survival models …")
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "ml07", Path(__file__).parent / "07_ml_onset_survival.py")
        ml07 = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(ml07)
        methods, means, errors, perm_dfs = ml07.run_survival_ml()
        with open(_cache4, "wb") as f:
            pickle.dump((methods, means, errors, perm_dfs), f)

    colors = ["firebrick", "darksalmon", "steelblue", "lightblue"]
    titles = ["Random Survival Forest", "Gradient Boosting Cox",
              "Linear Cox", "Elastic-Net Cox"]

    fig = plt.figure(figsize=(24, 20))
    gs  = gridspec.GridSpec(2, 2, figure=fig,
                            hspace=0.52, wspace=0.48,
                            width_ratios=[0.30, 0.70],
                            height_ratios=[1.0, 1.55])

    ax_a = fig.add_subplot(gs[0, 0])
    _draw_4A(ax_a)
    ax_a.tick_params(axis="both", labelsize=FS + 2)
    ax_a.set_xlabel(ax_a.get_xlabel(), fontsize=LAB + 2)
    _label(ax_a, "a")

    ax_b = fig.add_subplot(gs[0, 1])
    _draw_4B(ax_b)
    ax_b.tick_params(axis="both", labelsize=FS + 2)
    ax_b.set_xlabel(ax_b.get_xlabel(), fontsize=LAB + 2)
    ax_b.set_ylabel(ax_b.get_ylabel(), fontsize=LAB + 2)
    _label(ax_b, "b")

    # C – concordance bar chart (narrower column)
    ax_c = fig.add_subplot(gs[1, 0])
    bars = ax_c.bar(methods, means, yerr=errors, capsize=5,
                    color=colors, edgecolor=None, width=0.55)
    for bar, m, e in zip(bars, means, errors):
        ax_c.text(bar.get_x() + bar.get_width() / 2, m + e + 0.02,
                  f"{m:.3f}", ha="center", va="bottom", fontsize=FS + 2)
    ax_c.set_ylabel("Concordance Index", fontsize=LAB + 2)
    ax_c.set_ylim(0, 0.9)
    ax_c.tick_params(axis="both", labelsize=FS + 2)
    _title(ax_c, "Model concordance (C-index)")
    _polish(ax_c)
    _label(ax_c, "c")

    # D – permutation importances 2×2 (wider column)
    gs_d = gridspec.GridSpecFromSubplotSpec(
        2, 2, subplot_spec=gs[1, 1], hspace=0.55, wspace=0.50)
    for i, (m, title) in enumerate(zip(methods, titles)):
        ax_sub = fig.add_subplot(gs_d[i // 2, i % 2])
        df  = perm_dfs[m]
        med = df.median().sort_values(ascending=False)
        sns.boxplot(data=df[med.index], orient="h",
                    color=colors[i], fliersize=2, ax=ax_sub)
        ax_sub.set_title(title, fontsize=FS + 1, fontweight="bold")
        ax_sub.set_xlabel("ΔC-index", fontsize=LAB)
        ax_sub.set_ylabel("")
        ax_sub.tick_params(axis="both", labelsize=FS)
        _polish(ax_sub)
        if i == 0:
            _label(ax_sub, "d", x=-0.16)

    fig.savefig(OUT_MAIN / "fig4_combined.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved fig4_combined.pdf")


# ══════════════════════════════════════════════════════════════════════════════
# FIG 5 – Grain price & fish price IRFs
# ══════════════════════════════════════════════════════════════════════════════
def make_fig5():
    df_a = pd.read_csv(OUT_DATA / "fig5A_irf_grain_price.csv")
    df_b = pd.read_csv(OUT_DATA / "fig5B_irf_fish_price.csv")

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(18, 7))

    _draw_irf(ax_a, df_a, "label", "Price variation",
              ribbon_group="All regions",
              legend_kw=dict(loc="lower left", ncol=2, fontsize=FS - 2,
                             columnspacing=0.8, handlelength=1.5))
    _title(ax_a, "Grain price response to ENSO")
    _label(ax_a, "a")
    ax_a.tick_params(axis="both", labelsize=FS + 1)
    ax_a.set_xlabel(ax_a.get_xlabel(), fontsize=LAB + 1)
    ax_a.set_ylabel(ax_a.get_ylabel(), fontsize=LAB + 1)

    _draw_irf(ax_b, df_b, "species", "Price variation",
              ribbon_group="All")
    _title(ax_b, "Fish price response to ENSO")
    _label(ax_b, "b")
    ax_b.tick_params(axis="both", labelsize=FS + 1)
    ax_b.set_xlabel(ax_b.get_xlabel(), fontsize=LAB + 1)
    ax_b.set_ylabel(ax_b.get_ylabel(), fontsize=LAB + 1)

    fig.tight_layout(w_pad=4)
    fig.savefig(OUT_MAIN / "fig5_combined.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved fig5_combined.pdf")


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("=== Figure 1 ===")
    make_fig1()

    print("\n=== Figure 2 ===")
    make_fig2()

    print("\n=== Figure 3 ===")
    make_fig3()

    print("\n=== Figure 4 ===")
    make_fig4()

    print("\n=== Figure 5 ===")
    make_fig5()

    print(f"\nAll figures saved to {OUT_MAIN}")
