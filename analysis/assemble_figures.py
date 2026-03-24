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
FS   = 18   # tick / body
LAB  = 20   # axis label
PAN  = 26   # panel letter

mpl.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         FS,
    "axes.titlesize":    FS,
    "axes.labelsize":    LAB,
    "xtick.labelsize":   FS,
    "ytick.labelsize":   FS,
    "legend.fontsize":   FS - 2,
    "legend.title_fontsize": FS - 2,
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

# FontProperties for small-caps panel letters
_PAN_FP = fm.FontProperties(size=PAN, weight="bold", variant="small-caps")


def _polish(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for sp in ("bottom", "left"):
        ax.spines[sp].set_linewidth(1.4)
    ax.tick_params(axis="both", which="major", length=5, width=1.4,
                   direction="out")


def _label(ax, letter, x=-0.10, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes,
            font_properties=_PAN_FP, va="bottom", ha="left")


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
              col_map=IRF_COL, fill_map=IRF_FILL, ls_map=IRF_LS):
    """Generic IRF ribbon+line plot from a tidy dataframe."""
    groups = df[group_col].unique()
    for g in groups:
        sub  = df[df[group_col] == g].sort_values("horizon")
        col  = col_map.get(g, "black")
        fill = fill_map.get(g, col)
        ls   = ls_map.get(g, "-")
        if g == ribbon_group:
            ax.fill_between(sub["horizon"], sub["irf_down"], sub["irf_up"],
                            color=fill, alpha=0.20)
        ax.plot(sub["horizon"], sub["irf_mean"], color=col, lw=1.6,
                linestyle=ls, label=g)
    ax.axhline(0, color="gray", lw=0.8, linestyle="--")
    ax.set_xlabel("Horizon (years)")
    ax.set_ylabel(y_label)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    _pct_fmt(ax)
    if ylim:
        ax.set_ylim(ylim)
    leg = ax.legend(frameon=False, ncol=2, fontsize=FS - 2)
    _polish(ax)


# ══════════════════════════════════════════════════════════════════════════════
# FIG 1 – ENSO timeseries / famine Gantt+count / geographic bar-map
# ══════════════════════════════════════════════════════════════════════════════
def _draw_1A(ax, data):
    enso = data.groupby("Year")["nino34"].mean().reset_index()
    ax.plot(enso["Year"], enso["nino34"], color="black", lw=1.6, alpha=0.8)
    ax.axhline(0, color="black", lw=1.0, linestyle="--")

    spans = [
        (1590, 1600, "red",    0.25, "1590s/1690s Super Famines"),
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
              ncol=2, frameon=False, fontsize=FS - 2)
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
        fontsize=FS - 2)
    ax_gantt.set_xlim(data["Year"].min(), data["Year"].max())
    ax_gantt.tick_params(labelbottom=False)
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
              loc="lower left", fontsize=FS - 2, frameon=False)


def make_fig1():
    data = pd.read_csv(DATA_PROC / "famine_region_data.csv")
    regions, region_colors = _build_region_palette(data)

    fig = plt.figure(figsize=(18, 22))
    # Row 0: A full width
    # Row 1: B (left) + C (right), equal width
    gs = gridspec.GridSpec(2, 2, figure=fig,
                           height_ratios=[1, 1.5],
                           hspace=0.52, wspace=0.38)

    # A – full-width ENSO timeseries
    ax_a = fig.add_subplot(gs[0, :])
    _draw_1A(ax_a, data)
    _label(ax_a, "A")

    # B – Gantt + count stacked, bottom-left
    gs_b = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=gs[1, 0],
        height_ratios=[3, 1.6], hspace=0.08,
    )
    ax_b1 = fig.add_subplot(gs_b[0])
    ax_b2 = fig.add_subplot(gs_b[1], sharex=ax_b1)
    _draw_1B(ax_b1, ax_b2, data, regions, region_colors)
    _label(ax_b1, "B")

    # C – geo bar-map, bottom-right
    proj = ccrs.LambertConformal(central_longitude=10, central_latitude=50,
                                 standard_parallels=(45, 55))
    ax_c = fig.add_subplot(gs[1, 1], projection=proj)
    _draw_1C(ax_c, data, region_colors)
    ax_c.text(-0.05, 1.03, "C", transform=ax_c.transAxes,
              font_properties=_PAN_FP, va="bottom", ha="left")

    fig.savefig(OUT_MAIN / "fig1_combined.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved fig1_combined.pdf")


# ══════════════════════════════════════════════════════════════════════════════
# FIG 2 – Famine onset (A+B from R CSVs, C+D+E from ML)
# ══════════════════════════════════════════════════════════════════════════════
def _draw_2A(ax):
    df = pd.read_csv(OUT_DATA / "fig2A_onset_box.csv")
    groups = df["Group"].unique()
    colors = {"Famine Onset": "darksalmon", "No Famine": "lightblue"}
    for i, g in enumerate(["No Famine", "Famine Onset"]):
        sub = df[df["Group"] == g]["nino34"]
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

    ax.axhline(0, color="black", lw=1.0, linestyle="--")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["No Famine", "Famine Onset"])
    ax.set_ylabel("NINO3.4 (°C)")
    _polish(ax)


def _draw_2B(ax):
    df = pd.read_csv(OUT_DATA / "fig2B_onset_coef.csv")
    model_order  = ["No FEs", "FEs", "FEs + Controls"]
    model_colors = {"No FEs": "cornflowerblue", "FEs": "orange",
                    "FEs + Controls": "firebrick"}
    regions = df["Region"].unique()
    x       = np.arange(len(regions))
    width   = 0.22
    offsets = np.linspace(-(len(model_order) - 1) / 2,
                          (len(model_order) - 1) / 2,
                          len(model_order)) * width * 1.5

    for i, model in enumerate(model_order):
        sub = df[df["model"] == model].set_index("Region").reindex(regions)
        ax.errorbar(x + offsets[i], sub["estimate"],
                    yerr=[sub["estimate"] - sub["conf.low"],
                          sub["conf.high"] - sub["estimate"]],
                    fmt="o", color=model_colors[model],
                    elinewidth=1.2, capsize=3, ms=5,
                    label=model, zorder=3)

    ax.axhline(0, color="black", lw=0.9, linestyle="--")
    ax.set_xticks(x)
    ax.set_xticklabels(regions, rotation=30, ha="right", fontsize=FS - 1)
    ax.set_ylabel("Famine Onset Probability")
    ax.legend(frameon=False, fontsize=FS - 2, ncol=3,
              loc="upper right")
    _polish(ax)


def _run_ml_onset():
    """Run the GBC for 'All_features' and return result dict."""
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.model_selection import GridSearchCV, cross_val_score
    from sklearn.metrics import accuracy_score
    from sklearn.inspection import permutation_importance

    FEAT_LABELS = {
        "nino34": "NINO3.4", "ongoing_wars": "Ongoing Wars",
        "Deaths": "Conflict Deaths", "temp_winter": "NDJF Temp",
        "temp_summer": "AMJJ Temp", "precip_winter": "NDJF Precip",
        "precip_summer": "AMJJ Precip", "PDSI": "JJA scPDSI",
    }

    onset = pd.read_csv(DATA_PROC / "famine_region_data.csv")
    onset["Decade"] = (onset["Year"] // 10 * 10).astype(str)
    onset = onset[onset["Region"] == "Central Europe"].copy()

    feat_set = ["nino34", "ongoing_wars", "Deaths", "temp_winter",
                "temp_summer", "PDSI", "precip_winter"]
    X = onset[feat_set]
    y = onset["Famine_start"]

    param_grid = {"max_depth": [1, 3, 5], "n_estimators": [10, 50],
                  "learning_rate": [0.05, 0.1, 0.3],
                  "subsample": [0.7, 0.9], "max_features": ["sqrt"],
                  "min_samples_leaf": [2, 5], "min_samples_split": [5]}
    clf  = GradientBoostingClassifier(random_state=42)
    grid = GridSearchCV(clf, param_grid, cv=10, scoring="f1", n_jobs=-1)
    grid.fit(X, y)
    best = grid.best_estimator_
    print(f"  GBC best params: {grid.best_params_}")

    pred = best.predict(X)
    X_cf = X.copy()
    mask = (onset["Famine_start"].values == 1) & (onset["nino34"].values > 0)
    X_cf.loc[mask, "nino34"] = 0
    pred_cf = best.predict(X_cf)

    df_chron = pd.DataFrame({
        "Region":                onset["Region"].values,
        "Year":                  onset["Year"].values,
        "Observed":              y.values,
        "Predicted":             ((y.values == 1) & (pred == 1)).astype(int),
        "Counterfactual":        ((onset["Famine_start"].values == 1) &
                                  (pred_cf == 1)).astype(int),
        "NINO34_actual":         onset["nino34"].values,
        "NINO34_counterfactual": X_cf["nino34"].values,
    })

    perm = permutation_importance(best, X, y, n_repeats=50,
                                  random_state=42, scoring="f1", n_jobs=-1)
    df_imp = pd.concat([
        pd.DataFrame({"Feature": FEAT_LABELS.get(f, f),
                      "Importance": perm["importances"][i]})
        for i, f in enumerate(feat_set)
    ])

    cv_scores     = cross_val_score(best, X, y, cv=10, scoring="accuracy")
    in_sample_acc = accuracy_score(y, pred)

    return dict(df_chron=df_chron, df_imp=df_imp,
                cv_scores=cv_scores, in_sample_acc=in_sample_acc)


def _draw_2C_chron(ax_top, ax_bot, df_chron):
    years     = sorted(df_chron["Year"].unique())
    ymin, ymax = min(years), max(years)
    colors    = {"Observed": "firebrick", "Predicted": "darksalmon",
                 "Counterfactual": "cornflowerblue"}

    y_pos = {}
    for i, reg in enumerate(sorted(df_chron["Region"].unique())):
        base  = i * 6
        y_pos[reg] = {"Observed": base + 3,
                      "Predicted": base + 1.5,
                      "Counterfactual": base}
        rd = df_chron[df_chron["Region"] == reg].sort_values("Year")
        for sc in ["Observed", "Predicted", "Counterfactual"]:
            for yr in rd.loc[rd[sc] == 1, "Year"]:
                ax_top.barh(y=y_pos[reg][sc], width=1, left=yr, height=1.2,
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
    _polish(ax_top)

    nino_a = df_chron.groupby("Year")["NINO34_actual"].first().reindex(years, fill_value=0).values
    nino_c = df_chron.groupby("Year")["NINO34_counterfactual"].first().reindex(years, fill_value=0).values
    ax_bot.plot(years, nino_a, color="crimson",    lw=1.4, label="Observed",       alpha=0.85)
    ax_bot.plot(years, nino_c, color="dodgerblue", lw=1.4, label="Counterfactual", alpha=0.9)
    dmask = (nino_a != nino_c) & (nino_a >= 0.5)
    ax_bot.scatter(np.array(years)[dmask], nino_a[dmask],
                   color="red", s=25, zorder=5, label="> 0.5°C")
    ax_bot.axhline(0,   color="black", lw=0.8, linestyle="-",  alpha=0.3)
    ax_bot.axhline(0.5, color="black", lw=0.8, linestyle="--", alpha=0.5)
    ax_bot.set_xlim(ymin - 1, ymax + 1)
    ax_bot.set_xlabel("Year")
    ax_bot.set_ylabel("NINO3.4 (°C)")
    ax_bot.legend(frameon=False, fontsize=FS - 2)
    _polish(ax_bot)


def _draw_2D(ax, df_imp):
    med_order = df_imp.groupby("Feature")["Importance"].median().sort_values(ascending=False).index
    sns.boxplot(x="Importance", y="Feature", data=df_imp,
                order=med_order, palette="Reds_r", ax=ax)
    ax.set_xlabel("Permutation Importance")
    ax.set_ylabel("")
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
                color="black", fmt="none", capsize=4, elinewidth=1.5)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Accuracy")
    _polish(ax)


def make_fig2():
    print("  Running ML onset classifier …")
    ml = _run_ml_onset()

    fig = plt.figure(figsize=(18, 22))
    # Row 0: A (left) + B (right)
    # Row 1: C chronology (left ~54%) + D+E stacked (right ~45%)
    gs = gridspec.GridSpec(2, 1, figure=fig,
                           height_ratios=[1.0, 2.2],
                           hspace=0.52)

    # ── Row 0: A + B ──────────────────────────────────────────────────────────
    gs_top = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0],
                                              wspace=0.40)
    ax_a = fig.add_subplot(gs_top[0])
    _draw_2A(ax_a)
    _label(ax_a, "A")

    ax_b = fig.add_subplot(gs_top[1])
    _draw_2B(ax_b)
    _label(ax_b, "B")

    # ── Row 1: C (left 54%) + D/E stacked (right 45%) ─────────────────────────
    gs_bot = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1],
                                              width_ratios=[0.54, 0.45],
                                              wspace=0.42)

    # C – chronology (left half, split top/bottom)
    gs_c = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=gs_bot[0],
        height_ratios=[2, 1], hspace=0.08,
    )
    ax_c1 = fig.add_subplot(gs_c[0])
    ax_c2 = fig.add_subplot(gs_c[1])
    _draw_2C_chron(ax_c1, ax_c2, ml["df_chron"])
    _label(ax_c1, "C")

    # D + E stacked (right half)
    gs_de = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=gs_bot[1],
        height_ratios=[1.4, 0.8], hspace=0.52,
    )
    ax_d = fig.add_subplot(gs_de[0])
    _draw_2D(ax_d, ml["df_imp"])
    _label(ax_d, "D")

    ax_e = fig.add_subplot(gs_de[1])
    _draw_2E(ax_e, ml["in_sample_acc"], ml["cv_scores"])
    _label(ax_e, "E")

    fig.savefig(OUT_MAIN / "fig2_combined.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved fig2_combined.pdf")


def _draw_irf_errorbar(ax, df, group_col, y_label,
                       ylim=None,
                       col_map=IRF_COL, ls_map=IRF_LS):
    """IRF plot using points + whiskers (errorbar) – no ribbon lines."""
    groups = df[group_col].unique()
    n      = len(groups)
    offsets = np.linspace(-0.15, 0.15, n) if n > 1 else [0]
    for g, offset in zip(groups, offsets):
        sub = df[df[group_col] == g].sort_values("horizon")
        col = col_map.get(g, "black")
        ax.errorbar(sub["horizon"] + offset,
                    sub["irf_mean"],
                    yerr=[sub["irf_mean"] - sub["irf_down"],
                          sub["irf_up"]   - sub["irf_mean"]],
                    fmt="o", color=col, elinewidth=1.2, capsize=3, ms=5,
                    label=g, zorder=3)
    ax.axhline(0, color="gray", lw=0.8, linestyle="--")
    ax.set_xlabel("Horizon (years)")
    ax.set_ylabel(y_label)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    _pct_fmt(ax)
    if ylim:
        ax.set_ylim(ylim)
    ax.legend(frameon=False, ncol=2, fontsize=FS - 2)
    _polish(ax)


# ══════════════════════════════════════════════════════════════════════════════
# FIG 3 – Teleconnection maps (A+B) + yield IRFs (C+D)
# ══════════════════════════════════════════════════════════════════════════════
def _draw_teleco_map(ax, corr_vals, p_vals, lat2d, lon2d,
                     title, cmap, alpha=0.10):
    """Draw a single teleconnection map panel (used for A and B)."""
    import numpy as np
    from matplotlib.colors import BoundaryNorm

    levels  = np.linspace(-0.4, 0.4, 17)
    cmap_   = plt.get_cmap(cmap, len(levels) - 1)
    norm    = BoundaryNorm(levels, ncolors=cmap_.N)

    im = ax.pcolormesh(lon2d, lat2d, corr_vals,
                       cmap=cmap_, norm=norm, shading="auto",
                       transform=ccrs.PlateCarree())
    sig = p_vals <= alpha
    ax.scatter(lon2d[sig], lat2d[sig], s=4, color="black",
               transform=ccrs.PlateCarree(), zorder=5)
    ax.coastlines(linewidth=0.7)
    ax.add_feature(cfeature.BORDERS, linestyle=":", linewidth=0.5)
    ax.set_extent([-10, 35, 36, 70], crs=ccrs.PlateCarree())
    plt.colorbar(im, ax=ax, orientation="vertical",
                 fraction=0.046, pad=0.04, label="Pearson r")
    ax.set_title(title, fontsize=FS)


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
        summer = (xds_eu.sel(time=xds_eu.time.dt.month.isin([4, 5, 6, 7]))
                  .groupby(xds_eu.sel(time=xds_eu.time.dt.month.isin([4, 5, 6, 7]))
                            .time.dt.year)
                  .mean(dim="time"))
        # rebuild Year coordinate
        summer = summer.rename({"group": "Year"})
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


def make_fig3():
    tele = _compute_teleconnections()
    corr_pdsi, p_pdsi, lat2d_p, lon2d_p = tele[:4]
    corr_prec, p_prec, lat2d_pr, lon2d_pr = tele[4:]

    df_c = pd.read_csv(OUT_DATA / "fig3C_irf_yield_WR.csv")
    df_d = pd.read_csv(OUT_DATA / "fig3D_irf_yield_noWR.csv")

    proj = ccrs.EuroPP()
    fig  = plt.figure(figsize=(18, 18))
    gs   = gridspec.GridSpec(2, 2, figure=fig,
                             hspace=0.48, wspace=0.36)

    # A – PDSI map
    ax_a = fig.add_subplot(gs[0, 0], projection=proj)
    if corr_pdsi is not None:
        _draw_teleco_map(ax_a, corr_pdsi, p_pdsi, lat2d_p, lon2d_p,
                         "NINO3.4 × scPDSI (1500–1800)", "BrBG")
    else:
        ax_a.text(0.5, 0.5, "Data not found", transform=ax_a.transAxes,
                  ha="center", va="center")
    ax_a.text(-0.05, 1.04, "A", transform=ax_a.transAxes,
              font_properties=_PAN_FP, va="bottom", ha="left")

    # B – precipitation map
    ax_b = fig.add_subplot(gs[0, 1], projection=proj)
    if corr_prec is not None:
        _draw_teleco_map(ax_b, corr_prec, p_prec, lat2d_pr, lon2d_pr,
                         "NINO3.4 × Summer Precip (1500–1800)", "BrBG")
    else:
        ax_b.text(0.5, 0.5, "Data not found", transform=ax_b.transAxes,
                  ha="center", va="center")
    ax_b.text(-0.05, 1.04, "B", transform=ax_b.transAxes,
              font_properties=_PAN_FP, va="bottom", ha="left")

    # C – Wheat/Rye IRF (points + whiskers)
    ax_c = fig.add_subplot(gs[1, 0])
    _draw_irf_errorbar(ax_c, df_c, "label", "Wheat/Rye Harvest Variation",
                       ylim=(-0.13, 0.06))
    _label(ax_c, "C")

    # D – Other grain IRF (points + whiskers)
    ax_d = fig.add_subplot(gs[1, 1])
    _draw_irf_errorbar(ax_d, df_d, "label", "Other Grain Harvest Variation",
                       ylim=(-0.13, 0.06))
    _label(ax_d, "D")

    fig.savefig(OUT_MAIN / "fig3_combined.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved fig3_combined.pdf")


# ══════════════════════════════════════════════════════════════════════════════
# FIG 4 – Famine duration (A+B Cox, C+D ML survival)
# ══════════════════════════════════════════════════════════════════════════════
def _run_ml_survival():
    """Run RSF/GBSA/CoxPH/CoxNet and return concordance + importance data."""
    from sklearn.model_selection import GridSearchCV, cross_val_score, KFold
    from sklearn.inspection import permutation_importance
    from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
    from sksurv.ensemble import (RandomSurvivalForest,
                                 GradientBoostingSurvivalAnalysis)
    from sksurv.util import Surv
    from sksurv.metrics import concordance_index_censored

    FEAT_LABELS = {
        "avg_nino34": "NINO3.4", "avg_temp_summer": "AMJJ Temp",
        "avg_temp_winter": "NDJF Temp", "avg_precip_winter": "NDJF Precip",
        "avg_PDSI": "JJA scPDSI", "avg_Deaths": "Conflict Deaths",
        "avg_ongoing_wars": "Ongoing Wars",
    }

    famine = pd.read_csv(DATA_PROC / "famine_survival.csv")
    famine["event_observed"] = True
    features = list(FEAT_LABELS.keys())
    X = famine[features]
    X = (X - X.mean()) / X.std()
    y = Surv.from_dataframe("event_observed", "duration", famine.loc[X.index])

    def c_scorer(est, X_, y_):
        return concordance_index_censored(
            y_["event_observed"], y_["duration"], est.predict(X_))[0]

    icv = KFold(n_splits=10, shuffle=True, random_state=42)
    ocv = KFold(n_splits=10, shuffle=True, random_state=42)

    results = {}
    for name, model, param_grid in [
        ("RSF",
         RandomSurvivalForest(random_state=42, n_jobs=-1),
         {"n_estimators": [20, 50], "min_samples_split": [5, 10],
          "min_samples_leaf": [5, 10], "max_features": ["sqrt"]}),
        ("CoxGB",
         GradientBoostingSurvivalAnalysis(random_state=42),
         {"n_estimators": [20, 30], "learning_rate": [0.1, 0.5],
          "max_depth": [1, 3]}),
        ("CoxPH",
         CoxPHSurvivalAnalysis(),
         {}),
        ("CoxNet",
         CoxnetSurvivalAnalysis(max_iter=10000),
         {"alphas": [[0.01],[0.1],[1]], "l1_ratio": [0.5, 1.0]}),
    ]:
        print(f"  Running {name} …")
        if param_grid:
            gs  = GridSearchCV(model, param_grid, scoring=c_scorer,
                               cv=icv, n_jobs=-1)
            scores = cross_val_score(gs, X, y, cv=ocv, scoring=c_scorer)
            gs.fit(X, y)
            fitted = gs.best_estimator_
        else:
            scores = cross_val_score(model, X, y, cv=ocv, scoring=c_scorer)
            model.fit(X, y)
            fitted = model

        perm = permutation_importance(fitted, X, y, n_repeats=30,
                                      random_state=42, scoring=c_scorer)
        results[name] = {
            "scores": scores,
            "perm":   pd.DataFrame(perm.importances.T,
                                   columns=[FEAT_LABELS[f] for f in features])
        }

    methods = list(results.keys())
    means   = [results[m]["scores"].mean() for m in methods]
    errors  = [results[m]["scores"].std()  for m in methods]
    perm_dfs = {m: results[m]["perm"] for m in methods}
    return methods, means, errors, perm_dfs


def _draw_4A(ax):
    df = pd.read_csv(OUT_DATA / "fig4A_cox_coefs.csv")
    df = df.sort_values("coef")
    ax.errorbar(df["coef"], range(len(df)),
                xerr=[df["coef"] - df["lower"], df["upper"] - df["coef"]],
                fmt="o", color="firebrick", elinewidth=1.4, capsize=4, ms=6)
    ax.axvline(0, color="black", lw=0.9, linestyle="--")
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["Model"], fontsize=FS - 1)
    ax.set_xlabel("NINO3.4 log hazard ratio")
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

    # Delta annotation: El Niño (+2°) − La Niña (−1°) at year 3
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
        delta_pp = (surv_elni.iloc[idx_en] - surv_lani.iloc[idx_ln]) * 100
        ax.annotate(
            f"$\\Delta$ El Niño (+2°) $-$ La Niña ($-$1°) = {delta_pp:+.1f} pp\n"
            f"(at year {year_of_interest})",
            xy=(0.97, 0.95), xycoords="axes fraction",
            ha="right", va="top", fontsize=FS - 2,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.7),
        )
    except Exception:
        pass

    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:.0%}"))
    ax.set_xlabel("Famine duration (years)")
    ax.set_ylabel("Survival probability")
    ax.legend(frameon=False, fontsize=FS - 2, ncol=2)
    _polish(ax)


def make_fig4():
    print("  Running ML survival models …")
    methods, means, errors, perm_dfs = _run_ml_survival()

    colors = ["firebrick", "darksalmon", "steelblue", "lightblue"]
    titles = ["Random Survival Forest", "Gradient Boosting Cox",
              "Linear Cox", "Elastic-Net Cox"]

    fig = plt.figure(figsize=(18, 18))
    gs  = gridspec.GridSpec(2, 2, figure=fig,
                            hspace=0.52, wspace=0.40)

    ax_a = fig.add_subplot(gs[0, 0])
    _draw_4A(ax_a)
    _label(ax_a, "A")

    ax_b = fig.add_subplot(gs[0, 1])
    _draw_4B(ax_b)
    _label(ax_b, "B")

    # C – concordance bar chart
    ax_c = fig.add_subplot(gs[1, 0])
    bars = ax_c.bar(methods, means, yerr=errors, capsize=5,
                    color=colors, edgecolor=None, width=0.55)
    for bar, m, e in zip(bars, means, errors):
        ax_c.text(bar.get_x() + bar.get_width() / 2, m + e + 0.02,
                  f"{m:.3f}", ha="center", va="bottom", fontsize=FS)
    ax_c.set_ylabel("Concordance Index")
    ax_c.set_ylim(0, 0.9)
    _polish(ax_c)
    _label(ax_c, "C")

    # D – permutation importances 2×2
    gs_d = gridspec.GridSpecFromSubplotSpec(
        2, 2, subplot_spec=gs[1, 1], hspace=0.5, wspace=0.45)
    for i, (m, title) in enumerate(zip(methods, titles)):
        ax_sub = fig.add_subplot(gs_d[i // 2, i % 2])
        df  = perm_dfs[m]
        med = df.median().sort_values(ascending=False)
        sns.boxplot(data=df[med.index], orient="h",
                    color=colors[i], fliersize=2, ax=ax_sub)
        ax_sub.set_title(title, fontsize=FS - 1)
        ax_sub.set_xlabel("ΔC-index")
        ax_sub.set_ylabel("")
        _polish(ax_sub)
        if i == 0:
            _label(ax_sub, "D", x=-0.22)

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

    _draw_irf(ax_a, df_a, "label", "Log Grain Price",
              ribbon_group="All regions")
    _label(ax_a, "A")

    _draw_irf(ax_b, df_b, "species", "Log Fish Price",
              ribbon_group="All")
    _label(ax_b, "B")

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
