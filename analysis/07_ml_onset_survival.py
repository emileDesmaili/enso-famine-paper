"""
ML analyses for Figures 2 and 4  (Nature-style combined panels)

Figure 2 (5 panels, assembled here after 02_fig2_onset.R saves A & B PNGs):
  A – NINO3.4 boxplot famine onset vs no-famine    [from R: _fig2A_panel.png]
  B – LPM coefficients by region                   [from R: _fig2B_panel.png]
  C – ML classifier chronology & counterfactual
  D – ML permutation feature importances
  E – ML accuracy / skill

Figure 4 (2 ML panels, assembled here; R script provides A+B PNGs):
  C – Nested-CV concordance bar chart
  D – Permutation importances 2×2

Outputs → analysis/output/figures/main/
  fig2_combined.pdf
  fig4_combined.pdf   (if R panels present, else fig4CD_combined.pdf)

Outputs → analysis/output/figures/appendix/
  figA_ML_CI_fullsample.pdf
  figA_ML_featimp_fullsample.pdf
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
import seaborn as sns

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import GridSearchCV, KFold, cross_val_score, cross_val_predict
from sklearn.metrics import accuracy_score, confusion_matrix, roc_curve, auc
from sklearn.inspection import permutation_importance

from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sksurv.ensemble import RandomSurvivalForest, GradientBoostingSurvivalAnalysis
from sksurv.util import Surv
from sksurv.metrics import concordance_index_censored

# ── paths ──────────────────────────────────────────────────────────────────────
ROOT      = Path(__file__).parent.parent
DATA_PROC = ROOT / "processed data"
OUT_MAIN  = Path(__file__).parent / "output" / "figures" / "main"
OUT_APP   = Path(__file__).parent / "output" / "figures" / "appendix"
OUT_MAIN.mkdir(parents=True, exist_ok=True)
OUT_APP.mkdir(parents=True, exist_ok=True)

# ── Nature-style global rcParams ───────────────────────────────────────────────
FONT_SIZE  = 14
LABEL_SIZE = 18
PANEL_SIZE = 42

mpl.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
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
    "pdf.fonttype":        42,
    "ps.fonttype":         42,
})

FEAT_LABELS = {
    "nino34":        "NINO3.4",
    "ongoing_wars":  "Ongoing Wars",
    "Deaths":        "Conflict Deaths",
    "temp_winter":   "NDJF Temp",
    "temp_summer":   "AMJJ Temp",
    "precip_winter": "NDJF Precip",
    "precip_summer": "AMJJ Precip",
    "PDSI":          "JJA scPDSI",
    "avg_nino34":        "NINO3.4",
    "avg_temp_summer":   "AMJJ Temp",
    "avg_temp_winter":   "NDJF Temp",
    "avg_precip_winter": "NDJF Precip",
    "avg_PDSI":          "JJA scPDSI",
    "avg_Deaths":        "Conflict Deaths",
    "avg_ongoing_wars":  "Ongoing Wars",
}

PANEL_KW = dict(fontsize=PANEL_SIZE, fontweight="bold", transform=None)


def _polish(ax):
    """Apply clean Nature-style spine / tick styling."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for sp in ("bottom", "left"):
        ax.spines[sp].set_linewidth(1.4)
    ax.tick_params(axis="both", which="major", length=5, width=1.4,
                   direction="out", labelsize=FONT_SIZE)


def _panel_label(ax, letter, x=-0.16, y=1.05):
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=PANEL_SIZE, fontweight="bold",
            va="bottom", ha="left")


# ══════════════════════════════════════════════════════════════════════════════
# ① FAMINE ONSET – ML CLASSIFIER  (Fig 2 C, D, E)
# ══════════════════════════════════════════════════════════════════════════════

def run_onset_classifier():
    """
    Runs GBC for all 5 feature sets.  Returns results for the 'All_features'
    set (used in the main figure) and saves all individual-panel PDFs.
    Returns dict with keys: df_chron, df_imp, in_sample_acc, cv_scores.
    """
    onset = pd.read_csv(DATA_PROC / "famine_region_data.csv")
    onset["Decade"] = (onset["Year"] // 10 * 10).astype(str)
    onset = onset[onset["Region"] == "Central Europe"].copy()

    feature_combinations = [
        ["nino34"],
        ["nino34", "ongoing_wars", "Deaths"],
        ["nino34", "ongoing_wars", "Deaths", "temp_winter", "temp_summer"],
        ["nino34", "ongoing_wars", "Deaths", "temp_winter", "temp_summer",
         "precip_winter", "precip_summer"],
        ["nino34", "ongoing_wars", "Deaths", "temp_winter", "temp_summer",
         "PDSI", "precip_winter"],
    ]
    save_suffixes = [
        "NINO34", "NINO_conflict", "Nino_conflict_temp",
        "nino_conflict_precip", "All_features",
    ]

    param_grid = {
        "max_depth":          [1, 3, 5, 10],
        "n_estimators":       [3, 5, 10, 50],
        "learning_rate":      [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5],
        "subsample":          [0.7, 0.9],
        "max_features":       ["sqrt", 0.5],
        "min_samples_leaf":   [2, 5],
        "min_samples_split":  [5, 10],
    }

    from sklearn.metrics import f1_score as _f1

    main_result = None
    all_models  = []  # (suffix, best, threshold, X, y, oos_f1)

    for feat_set, suffix in zip(feature_combinations, save_suffixes):
        X = onset[feat_set]
        y = onset["Famine_start"]

        clf  = GradientBoostingClassifier(random_state=42)
        grid = GridSearchCV(clf, param_grid, cv=10,
                            scoring="f1", n_jobs=-1)
        grid.fit(X, y)
        best    = grid.best_estimator_
        best_thr = 0.5
        print(f"[{suffix}] best params: {grid.best_params_}")

        # In-sample predictions
        prob_is = best.predict_proba(X)[:, 1]
        pred    = (prob_is >= best_thr).astype(int)

        # Counterfactual at same threshold
        X_cf = X.copy()
        if "nino34" in feat_set:
            mask = (onset["Famine_start"].values == 1) & (onset["nino34"].values > 0)
            X_cf.loc[mask, "nino34"] = 0
        pred_cf = (best.predict_proba(X_cf)[:, 1] >= best_thr).astype(int)

        df_chron = pd.DataFrame({
            "Region":                onset["Region"].values,
            "Year":                  onset["Year"].values,
            "Observed":              y.values,
            "Predicted":             pred,
            "Counterfactual":        pred_cf,
            "NINO34_actual":         onset["nino34"].values,
            "NINO34_counterfactual": X_cf["nino34"].values if "nino34" in feat_set else 0,
            "Famine_dur":            onset["Famine_dur"].values if "Famine_dur" in onset.columns else 1,
        })
        df_chron["Predicted"]      = ((df_chron["Observed"] == 1) &
                                      (df_chron["Predicted"] == 1)).astype(int)
        df_chron["Counterfactual"] = ((onset["Famine_start"].values == 1) &
                                      (df_chron["Counterfactual"] == 1)).astype(int)

        perm = permutation_importance(best, X, y, n_repeats=50,
                                      random_state=42, scoring="accuracy",
                                      n_jobs=-1)
        df_imp = pd.concat([
            pd.DataFrame({"Feature": f, "Importance": perm["importances"][i]})
            for i, f in enumerate(feat_set)
        ])
        df_imp["Feature"] = df_imp["Feature"].replace(FEAT_LABELS)

        cv_scores     = cross_val_score(best, X, y, cv=10, scoring="accuracy")
        in_sample_acc = accuracy_score(y, pred)

        # OOS F1 via cross_val_predict at default threshold
        prob_oos_cv  = cross_val_predict(best, X, y, cv=10, method="predict_proba")[:, 1]
        pred_oos_thr = (prob_oos_cv >= best_thr).astype(int)
        oos_f1       = _f1(y, pred_oos_thr, zero_division=0)
        print(f"[{suffix}] OOS F1: {oos_f1:.3f}")
        all_models.append((suffix, best, best_thr, X, y, oos_f1))

        _save_chronology_pdf(df_chron, onset, suffix)
        _save_importances_pdf(df_imp, suffix)
        _save_skill_pdf(in_sample_acc, cv_scores, suffix)
        _save_counts_pdf(df_chron, suffix)

        if suffix == "All_features":
            _save_confusion_roc_pdf(best, X, y, threshold=best_thr,
                                    cv=10, suffix="All_features")
            main_result = dict(
                df_chron=df_chron,
                onset=onset,
                df_imp=df_imp,
                in_sample_acc=in_sample_acc,
                cv_scores=cv_scores,
            )

    # Save confusion/metrics + reliability figures for model with highest OOS F1
    best_entry  = max(all_models, key=lambda t: t[5])
    best_suffix, best_best, best_thr, best_X, best_y, best_f1 = best_entry
    print(f"Best OOS F1 model: {best_suffix} (F1={best_f1:.3f})")
    if best_suffix != "All_features":
        _save_confusion_roc_pdf(best_best, best_X, best_y, threshold=best_thr,
                                cv=10, suffix=best_suffix)

    return main_result


def _save_chronology_pdf(df_chron, onset_df, suffix):
    years      = sorted(df_chron["Year"].unique())
    year_min, year_max = min(years), max(years)
    colors     = {"Observed": "firebrick", "Predicted": "darksalmon",
                  "Counterfactual": "cornflowerblue"}

    # Use actual famine duration for bar widths (matches assemble_figures logic)
    dur_lookup = (onset_df[onset_df["Famine_start"] == 1]
                  .set_index("Year")["Famine_dur"].to_dict()
                  if "Famine_dur" in onset_df.columns else {})

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6),
                                   height_ratios=[1, 1], sharex=False)
    y_positions = {}
    for i, region in enumerate(sorted(df_chron["Region"].unique())):
        base_y = i * 6
        y_positions[region] = {"Observed": base_y + 3,
                                "Predicted": base_y + 1.5,
                                "Counterfactual": base_y}
        rd = df_chron[df_chron["Region"] == region].sort_values("Year")
        for scenario in ["Observed", "Predicted", "Counterfactual"]:
            for yr in rd.loc[rd[scenario] == 1, "Year"].values:
                dur = dur_lookup.get(yr, 1)
                ax1.barh(y=y_positions[region][scenario], width=dur, left=yr,
                         height=1.2, color=colors[scenario],
                         edgecolor="black", linewidth=1.2, alpha=0.7)

    y_ticks, y_labels = [], []
    for region in sorted(df_chron["Region"].unique()):
        for scenario in ["Counterfactual", "Predicted", "Observed"]:
            cnt = df_chron[df_chron["Region"] == region][scenario].sum()
            y_ticks.append(y_positions[region][scenario])
            y_labels.append(f"{scenario} ({cnt})")
    ax1.set_yticks(y_ticks)
    ax1.set_yticklabels(y_labels)
    ax1.set_xlim(year_min - 1, year_max + 1)
    ax1.set_xlabel("Year")
    _polish(ax1)

    nino_actual = (df_chron.groupby("Year")["NINO34_actual"]
                   .first().reindex(years, fill_value=0).values)
    nino_cf     = (df_chron.groupby("Year")["NINO34_counterfactual"]
                   .first().reindex(years, fill_value=0).values)
    ax2.plot(years, nino_actual, color="crimson",    linewidth=1.4,
             label="Observed", alpha=0.8)
    ax2.plot(years, nino_cf,     color="dodgerblue", linewidth=1.4,
             label="Counterfactual", alpha=0.9)
    diff_mask = (np.array(nino_actual) != np.array(nino_cf)) & \
                (np.array(nino_actual) >= 0.5)
    ax2.scatter(np.array(years)[diff_mask], np.array(nino_actual)[diff_mask],
                color="red", s=30, zorder=5, label="> 0.5°C")
    ax2.axhline(0,   color="black", linestyle="-",  alpha=0.3, linewidth=0.8)
    ax2.axhline(0.5, color="black", linestyle="--", alpha=0.5, linewidth=0.8)
    ax2.set_xlim(year_min - 1, year_max + 1)
    ax2.set_xlabel("Year")
    ax2.set_ylabel("NINO3.4 (°C)")
    ax2.legend(loc="lower left", frameon=False)
    _polish(ax2)

    fig.tight_layout()
    fig.savefig(OUT_APP / f"MLOnset_{suffix}_chronology.pdf",
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved chronology for {suffix}")


def _save_importances_pdf(df_imp, suffix):
    medians   = df_imp.groupby("Feature")["Importance"].median().sort_values(ascending=False)
    med_order = medians.index
    norm      = plt.Normalize(medians.min(), medians.max())
    palette   = {f: plt.cm.Reds(norm(v)) for f, v in medians.items()}
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.boxplot(x="Importance", y="Feature", hue="Feature", data=df_imp,
                order=med_order, palette=palette, legend=False, ax=ax)
    ax.set_xlabel("Permutation Importance")
    ax.set_ylabel("")
    _polish(ax)
    fig.tight_layout()
    fig.savefig(OUT_APP / f"MLOnset_{suffix}_importances.pdf",
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved importances for {suffix}")


def _save_skill_pdf(in_sample_acc, cv_scores, suffix):
    skill_df = pd.DataFrame({
        "Skill":    ["In-sample", "Cross-Validation"],
        "Accuracy": [in_sample_acc, cv_scores.mean()],
        "CI":       [0, cv_scores.std()],
    })
    fig, ax = plt.subplots(figsize=(5, 2.5))
    sns.barplot(data=skill_df, y="Skill", x="Accuracy", hue="Skill",
                palette={"In-sample": "firebrick", "Cross-Validation": "darksalmon"},
                legend=False, ax=ax)
    ax.errorbar(x=cv_scores.mean(), y=1, xerr=cv_scores.std(),
                color="black", fmt="none", capsize=4, elinewidth=1.5)
    ax.set_xlim(0, 1)
    #remove y axis label that says skill
    ax.set_ylabel("")
    ax.set_xlabel("Accuracy")
    _polish(ax)
    fig.tight_layout()
    fig.savefig(OUT_APP / f"MLOnset_{suffix}_skillHor.pdf",
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved skill for {suffix}")


def _save_counts_pdf(df_chron, suffix):
    """Bar plot of total Observed / Predicted / Counterfactual famine counts."""
    counts = (df_chron[["Observed", "Predicted", "Counterfactual"]]
              .sum()
              .reindex(["Observed", "Predicted", "Counterfactual"]))
    palette = {"Observed": "firebrick", "Predicted": "darksalmon",
               "Counterfactual": "cornflowerblue"}
    fig, ax = plt.subplots(figsize=(4, 5))
    bars = ax.bar(counts.index, counts.values,
                  color=[palette[k] for k in counts.index], edgecolor=None)
    for bar, v in zip(bars, counts.values):
        ax.text(bar.get_x() + bar.get_width() / 2, v + 0.05,
                str(int(v)), ha="center", va="bottom",
                fontsize=FONT_SIZE, fontweight="bold")
    ax.set_ylabel("Number of famines")
    ax.set_xlabel("")
    ax.set_yticks([])
    ax.set_ylim(0, max(counts.values) * 1.25)
    _polish(ax)
    fig.tight_layout()
    fig.savefig(OUT_APP / f"MLOnset_{suffix}_counts.pdf",
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved counts for {suffix}")


def _save_confusion_roc_pdf(best, X, y, threshold=0.5, cv=10, suffix="All_features"):
    """Confusion matrices (A, B) + metrics bar chart (C) + Brier scores (D)."""
    from sklearn.metrics import f1_score, brier_score_loss

    # In-sample
    prob_is  = best.predict_proba(X)[:, 1]
    pred_is  = (prob_is >= threshold).astype(int)
    cm_is    = confusion_matrix(y, pred_is)
    fpr_is, tpr_is, _ = roc_curve(y, prob_is)
    metrics_is = {
        "Accuracy": accuracy_score(y, pred_is),
        "F1":       f1_score(y, pred_is, zero_division=0),
        "AUROC":    auc(fpr_is, tpr_is),
    }
    brier_is = brier_score_loss(y, prob_is)

    # Out-of-sample (stratified CV)
    prob_oos = cross_val_predict(best, X, y, cv=cv, method="predict_proba")[:, 1]
    pred_oos = (prob_oos >= threshold).astype(int)
    cm_oos   = confusion_matrix(y, pred_oos)
    fpr_oos, tpr_oos, _ = roc_curve(y, prob_oos)
    metrics_oos = {
        "Accuracy": accuracy_score(y, pred_oos),
        "F1":       f1_score(y, pred_oos, zero_division=0),
        "AUROC":    auc(fpr_oos, tpr_oos),
    }
    brier_oos = brier_score_loss(y, prob_oos)

    fig = plt.figure(figsize=(18, 10))
    gs  = gridspec.GridSpec(2, 2, figure=fig,
                            height_ratios=[1.4, 1.0],
                            hspace=0.52, wspace=0.42)

    cm_labels = ["No famine", "Famine"]

    # ── Top row: confusion matrices (A, B) ──────────────────────────────────
    for col, (cm, title) in enumerate([
        (cm_is,  "Confusion matrix\n(in-sample)"),
        (cm_oos, "Confusion matrix\n(OOS, stratified 10-fold CV)"),
    ]):
        ax = fig.add_subplot(gs[0, col])
        cm_norm = cm.astype(float) / cm.sum(axis=1, keepdims=True)
        im = ax.imshow(cm_norm, cmap="Reds", vmin=0, vmax=1)
        ax.set_xticks([0, 1]); ax.set_xticklabels(cm_labels, fontsize=FONT_SIZE - 3)
        ax.set_yticks([0, 1]); ax.set_yticklabels(cm_labels, fontsize=FONT_SIZE - 3)
        ax.set_xlabel("Predicted", fontsize=FONT_SIZE - 1)
        ax.set_ylabel("Actual",    fontsize=FONT_SIZE - 1)
        ax.set_title(title, fontsize=FONT_SIZE - 1, pad=8)
        for i in range(2):
            for j in range(2):
                ax.text(j, i, f"{cm[i,j]}\n({cm_norm[i,j]:.0%})",
                        ha="center", va="center",
                        color="white" if cm_norm[i, j] > 0.6 else "black",
                        fontsize=FONT_SIZE - 2)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        letter = "A" if col == 0 else "B"
        ax.text(-0.18, 1.08, letter, transform=ax.transAxes,
                fontsize=PANEL_SIZE, fontweight="bold", va="bottom", ha="left")

    # ── Bottom-left: grouped metric bars (C) ────────────────────────────────
    metric_names = list(metrics_is.keys())
    x     = np.arange(len(metric_names))
    width = 0.35

    ax_m = fig.add_subplot(gs[1, 0])
    bars_is  = ax_m.bar(x - width / 2,
                        [metrics_is[m]  for m in metric_names],
                        width, label="In-sample",        color="firebrick",  alpha=0.85)
    bars_oos = ax_m.bar(x + width / 2,
                        [metrics_oos[m] for m in metric_names],
                        width, label="OOS (10-fold CV)", color="steelblue", alpha=0.85)
    for bar in list(bars_is) + list(bars_oos):
        h = bar.get_height()
        ax_m.text(bar.get_x() + bar.get_width() / 2, h + 0.01,
                  f"{h:.2f}", ha="center", va="bottom",
                  fontsize=FONT_SIZE - 4)
    ax_m.set_xticks(x)
    ax_m.set_xticklabels(metric_names, fontsize=FONT_SIZE - 2)
    ax_m.set_ylim(0, 1.35)
    ax_m.set_ylabel("Score", fontsize=FONT_SIZE)
    ax_m.legend(frameon=False, fontsize=FONT_SIZE - 2,
                loc="upper right", bbox_to_anchor=(1.0, 1.32))
    _polish(ax_m)
    ax_m.text(-0.14, 1.06, "C", transform=ax_m.transAxes,
              fontsize=PANEL_SIZE, fontweight="bold", va="bottom", ha="left")

    # ── Bottom-right: Brier score bar chart (D) ──────────────────────────────
    ax_b = fig.add_subplot(gs[1, 1])
    brier_labels = ["In-sample", "OOS (10-fold CV)"]
    brier_vals   = [brier_is, brier_oos]
    bars_b = ax_b.bar(brier_labels, brier_vals,
                      color=["firebrick", "steelblue"], alpha=0.85, edgecolor=None)
    for bar, v in zip(bars_b, brier_vals):
        ax_b.text(bar.get_x() + bar.get_width() / 2, v + 0.003,
                  f"{v:.3f}", ha="center", va="bottom",
                  fontsize=FONT_SIZE - 1, fontweight="bold")
    ax_b.set_ylabel("Brier score", fontsize=FONT_SIZE)
    ax_b.set_ylim(0, 1)
    ax_b.text(0.97, 0.97, "$\\downarrow$ lower is better", transform=ax_b.transAxes,
              ha="right", va="top", fontsize=FONT_SIZE + 2, color="dimgray")
    _polish(ax_b)
    ax_b.text(-0.14, 1.06, "D", transform=ax_b.transAxes,
              fontsize=PANEL_SIZE, fontweight="bold", va="bottom", ha="left")

    out_path = OUT_APP / f"figA_ML_onset_confusion_roc_{suffix}.pdf"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_path.name}")


# ── Draw individual subpanels into an Axes (for combined figure) ───────────────

def _draw_C_chronology(ax_chron, ax_nino, df_chron):
    years      = sorted(df_chron["Year"].unique())
    year_min, year_max = min(years), max(years)
    colors     = {"Observed": "firebrick", "Predicted": "darksalmon",
                  "Counterfactual": "cornflowerblue"}
    dur_lookup = (df_chron[df_chron["Observed"] == 1]
                  .set_index("Year")["Famine_dur"].to_dict()
                  if "Famine_dur" in df_chron.columns else {})

    y_positions = {}
    for i, region in enumerate(sorted(df_chron["Region"].unique())):
        base_y = i * 6
        y_positions[region] = {"Observed": base_y + 3,
                                "Predicted": base_y + 1.5,
                                "Counterfactual": base_y}
        rd = df_chron[df_chron["Region"] == region].sort_values("Year")
        for scenario in ["Observed", "Predicted", "Counterfactual"]:
            for yr in rd.loc[rd[scenario] == 1, "Year"].values:
                dur = dur_lookup.get(yr, 1)
                ax_chron.barh(y=y_positions[region][scenario], width=dur, left=yr,
                              height=1.2, color=colors[scenario],
                              edgecolor="black", linewidth=1.0, alpha=0.75)

    y_ticks, y_labels = [], []
    for region in sorted(df_chron["Region"].unique()):
        for scenario in ["Counterfactual", "Predicted", "Observed"]:
            cnt = df_chron[df_chron["Region"] == region][scenario].sum()
            y_ticks.append(y_positions[region][scenario])
            y_labels.append(f"{scenario} ({cnt})")
    ax_chron.set_yticks(y_ticks)
    ax_chron.set_yticklabels(y_labels, fontsize=FONT_SIZE - 1)
    ax_chron.set_xlim(year_min - 1, year_max + 1)
    ax_chron.tick_params(labelbottom=False)
    _polish(ax_chron)

    nino_actual = (df_chron.groupby("Year")["NINO34_actual"]
                   .first().reindex(years, fill_value=0).values)
    nino_cf     = (df_chron.groupby("Year")["NINO34_counterfactual"]
                   .first().reindex(years, fill_value=0).values)
    ax_nino.plot(years, nino_actual, color="crimson",    linewidth=1.4,
                 label="Observed", alpha=0.85)
    ax_nino.plot(years, nino_cf,     color="dodgerblue", linewidth=1.4,
                 label="Counterfactual", alpha=0.9)
    diff_mask = (np.array(nino_actual) != np.array(nino_cf)) & \
                (np.array(nino_actual) >= 0.5)
    ax_nino.scatter(np.array(years)[diff_mask], np.array(nino_actual)[diff_mask],
                    color="red", s=25, zorder=5, label="> 0.5°C")
    ax_nino.axhline(0,   color="black", linestyle="-",  alpha=0.3, linewidth=0.8)
    ax_nino.axhline(0.5, color="black", linestyle="--", alpha=0.5, linewidth=0.8)
    ax_nino.set_xlim(year_min - 1, year_max + 1)
    ax_nino.set_xlabel("Year")
    ax_nino.set_ylabel("NINO3.4 (°C)")
    ax_nino.legend(loc="lower left", frameon=False, fontsize=FONT_SIZE - 2)
    _polish(ax_nino)


def _draw_D_importances(ax, df_imp):
    medians   = df_imp.groupby("Feature")["Importance"].median().sort_values(ascending=False)
    med_order = medians.index
    norm      = plt.Normalize(medians.min(), medians.max())
    palette   = {f: plt.cm.Reds(norm(v)) for f, v in medians.items()}
    sns.boxplot(x="Importance", y="Feature", hue="Feature", data=df_imp,
                order=med_order, palette=palette, legend=False, ax=ax)
    ax.set_xlabel("Permutation Importance")
    ax.set_ylabel("")
    _polish(ax)


def _draw_E_skill(ax, in_sample_acc, cv_scores):
    skill_df = pd.DataFrame({
        "Skill":    ["In-sample", "Cross-Validation"],
        "Accuracy": [in_sample_acc, cv_scores.mean()],
    })
    sns.barplot(data=skill_df, y="Skill", x="Accuracy", hue="Skill",
                palette={"In-sample": "firebrick", "Cross-Validation": "darksalmon"},
                legend=False, ax=ax)
    ax.errorbar(x=cv_scores.mean(), y=1, xerr=cv_scores.std(),
                color="black", fmt="none", capsize=4, elinewidth=1.5)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Accuracy")
    _polish(ax)


# ── Assemble Figure 2  (A+B from R PNGs, C+D+E from Python) ──────────────────

def assemble_fig2(main_result):
    png_a = OUT_APP / "_fig2A_panel.png"
    png_b = OUT_APP / "_fig2B_panel.png"
    have_r_panels = png_a.exists() and png_b.exists()

    df_chron      = main_result["df_chron"]
    df_imp        = main_result["df_imp"]
    in_sample_acc = main_result["in_sample_acc"]
    cv_scores     = main_result["cv_scores"]

    if have_r_panels:
        # 5-panel layout: top row (A, B), bottom 3 rows (C, D, E)
        fig = plt.figure(figsize=(18, 22))
        gs  = gridspec.GridSpec(
            3, 2,
            figure=fig,
            height_ratios=[1.0, 1.4, 0.5],
            hspace=0.55,
            wspace=0.38,
        )

        # A  (image from R)
        ax_a = fig.add_subplot(gs[0, 0])
        ax_a.imshow(mpimg.imread(str(png_a)))
        ax_a.axis("off")
        ax_a.text(-0.05, 1.04, "A", transform=ax_a.transAxes,
                  fontsize=PANEL_SIZE, fontweight="bold", va="bottom", ha="left")

        # B  (image from R)
        ax_b = fig.add_subplot(gs[0, 1])
        ax_b.imshow(mpimg.imread(str(png_b)))
        ax_b.axis("off")
        ax_b.text(-0.05, 1.04, "B", transform=ax_b.transAxes,
                  fontsize=PANEL_SIZE, fontweight="bold", va="bottom", ha="left")
    else:
        print("  Note: R panel PNGs not found; building Fig2 with C/D/E only.")
        fig = plt.figure(figsize=(18, 14))
        gs  = gridspec.GridSpec(
            2, 2,
            figure=fig,
            height_ratios=[1.4, 0.5],
            hspace=0.55, wspace=0.38,
        )

    # C  (chronology spans full width, two internal rows)
    if have_r_panels:
        gs_c = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=gs[1, :],
            height_ratios=[2, 1], hspace=0.08,
        )
    else:
        gs_c = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=gs[0, :],
            height_ratios=[2, 1], hspace=0.08,
        )
    ax_c1 = fig.add_subplot(gs_c[0])
    ax_c2 = fig.add_subplot(gs_c[1])
    _draw_C_chronology(ax_c1, ax_c2, df_chron)
    ax_c1.text(-0.05, 1.05, "C", transform=ax_c1.transAxes,
               fontsize=PANEL_SIZE, fontweight="bold", va="bottom", ha="left")

    # D  (importances)
    row_de = 2 if have_r_panels else 1
    ax_d = fig.add_subplot(gs[row_de, 0])
    _draw_D_importances(ax_d, df_imp)
    ax_d.text(-0.12, 1.05, "D", transform=ax_d.transAxes,
              fontsize=PANEL_SIZE, fontweight="bold", va="bottom", ha="left")

    # E  (skill)
    ax_e = fig.add_subplot(gs[row_de, 1])
    _draw_E_skill(ax_e, in_sample_acc, cv_scores)
    ax_e.text(-0.12, 1.05, "E", transform=ax_e.transAxes,
              fontsize=PANEL_SIZE, fontweight="bold", va="bottom", ha="left")

    fig.savefig(OUT_APP / "MLOnset_combined.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved fig2_combined.pdf")


# ══════════════════════════════════════════════════════════════════════════════
# ② SURVIVAL ML  (Fig 4 C, D)
# ══════════════════════════════════════════════════════════════════════════════

def concordance_scorer(estimator, X, y):
    pred = estimator.predict(X)
    return concordance_index_censored(y["event_observed"], y["duration"], pred)[0]


def run_survival_ml():
    famine = pd.read_csv(DATA_PROC / "famine_survival.csv")
    famine["event_observed"] = True

    features = ["avg_nino34", "avg_temp_summer", "avg_temp_winter",
                "avg_precip_winter", "avg_PDSI", "avg_Deaths", "avg_ongoing_wars"]
    X = famine[features]
    X = (X - X.mean()) / X.std()
    y = Surv.from_dataframe("event_observed", "duration", famine.loc[X.index])

    inner_cv = KFold(n_splits=10, shuffle=True, random_state=42)
    outer_cv = KFold(n_splits=10, shuffle=True, random_state=42)

    print("Running RSF …")
    rsf    = RandomSurvivalForest(random_state=42, n_jobs=-1)
    gs_rsf = GridSearchCV(rsf,
                          {"n_estimators": [10, 20, 50],
                           "min_samples_split": [5, 10, 15],
                           "min_samples_leaf": [5, 10, 15],
                           "max_features": ["sqrt", 0.5]},
                          scoring=concordance_scorer, cv=inner_cv, n_jobs=-1)
    rsf_scores = cross_val_score(gs_rsf, X, y, cv=outer_cv,
                                 scoring=concordance_scorer)
    gs_rsf.fit(X, y)
    perm_rsf = permutation_importance(gs_rsf.best_estimator_, X, y,
                                      n_repeats=50, random_state=42,
                                      scoring=concordance_scorer)

    print("Running GBSA …")
    gbsa    = GradientBoostingSurvivalAnalysis(random_state=42)
    gs_gbsa = GridSearchCV(gbsa,
                           {"n_estimators": [10, 20, 30],
                            "learning_rate": [0.1, 0.5, 0.9],
                            "max_depth": [1, 3, 5]},
                           scoring=concordance_scorer, cv=inner_cv, n_jobs=-1)
    gbsa_scores = cross_val_score(gs_gbsa, X, y, cv=outer_cv,
                                  scoring=concordance_scorer)
    gs_gbsa.fit(X, y)
    perm_gbsa = permutation_importance(gs_gbsa.best_estimator_, X, y,
                                       n_repeats=50, random_state=42,
                                       scoring=concordance_scorer)

    print("Running CoxPH …")
    coxph        = CoxPHSurvivalAnalysis()
    coxph_scores = cross_val_score(coxph, X, y, cv=outer_cv,
                                   scoring=concordance_scorer)
    coxph.fit(X, y)
    perm_coxph = permutation_importance(coxph, X, y, n_repeats=50,
                                        random_state=42,
                                        scoring=concordance_scorer)

    print("Running CoxNet …")
    cox_lasso  = CoxnetSurvivalAnalysis(max_iter=10000)
    gs_lasso   = GridSearchCV(cox_lasso,
                              {"alphas":   [[0.001],[0.01],[0.1],[1],[10],[100]],
                               "l1_ratio": [0.01, 0.5, 1.0]},
                              scoring=concordance_scorer, cv=inner_cv, n_jobs=-1)
    lasso_scores = cross_val_score(gs_lasso, X, y, cv=outer_cv,
                                   scoring=concordance_scorer)
    gs_lasso.fit(X, y)
    perm_lasso = permutation_importance(gs_lasso.best_estimator_, X, y,
                                        n_repeats=50, random_state=42,
                                        scoring=concordance_scorer)

    means   = [rsf_scores.mean(),   gbsa_scores.mean(),
               coxph_scores.mean(), lasso_scores.mean()]
    errors  = [rsf_scores.std(),    gbsa_scores.std(),
               coxph_scores.std(),  lasso_scores.std()]
    methods = ["RSF", "CoxGB", "CoxPH", "CoxNet"]
    colors  = ["firebrick", "darksalmon", "steelblue", "lightblue"]
    titles  = ["Random Survival Forest", "Gradient Boosting Cox",
               "Linear Cox", "Elastic-Net Cox"]

    print("Nested CV C-indices:", dict(zip(methods, [f"{m:.3f}" for m in means])))

    cols_rename = {f: FEAT_LABELS.get(f, f) for f in features}
    perm_data   = {
        "RSF":    pd.DataFrame(perm_rsf.importances.T,
                               columns=list(cols_rename.values())),
        "CoxGB":  pd.DataFrame(perm_gbsa.importances.T,
                               columns=list(cols_rename.values())),
        "CoxPH":  pd.DataFrame(perm_coxph.importances.T,
                               columns=list(cols_rename.values())),
        "CoxNet": pd.DataFrame(perm_lasso.importances.T,
                               columns=list(cols_rename.values())),
    }

    # ── Nested-CV concordance bar plot (appendix standalone) ─────────────────
    colors  = ["firebrick", "darksalmon", "steelblue", "lightblue"]
    fig_ci, ax_ci = plt.subplots(figsize=(4, 6))
    bars = ax_ci.bar(methods, means, yerr=errors, capsize=4,
                     color=colors, edgecolor=None)
    for bar, m, e in zip(bars, means, errors):
        ax_ci.text(bar.get_x() + bar.get_width() / 2, m + e + 0.02,
                   f"{m:.3f}", ha="center", va="bottom", fontsize=FONT_SIZE)
    ax_ci.set_ylabel("Concordance Index")
    ax_ci.set_ylim(0, 0.85)
    _polish(ax_ci)
    fig_ci.tight_layout()
    fig_ci.savefig(OUT_APP / "MLSurvival_CI_nestedCV.pdf", dpi=300)
    plt.close(fig_ci)
    print("Saved MLSurvival_CI_nestedCV.pdf")

    # ── Nested-CV feature importance boxplots (appendix standalone) ───────────
    titles  = ["Random Survival Forest", "Gradient Boosting Cox",
               "Linear Cox", "Elastic-Net Cox"]
    fig_fi, axes_fi = plt.subplots(2, 2, figsize=(14, 8), sharey=False)
    axes_fi = axes_fi.flatten()
    for i, (model_name, df) in enumerate(perm_data.items()):
        medians   = df.median().sort_values(ascending=False)
        df_sorted = df[medians.index]
        sns.boxplot(data=df_sorted, orient="h", ax=axes_fi[i],
                    color=colors[i], fliersize=3)
        axes_fi[i].set_title(titles[i], fontsize=FONT_SIZE)
        axes_fi[i].set_xlabel("Permutation Importance (Δ C-index)")
        axes_fi[i].set_ylabel("")
        axes_fi[i].xaxis.set_major_formatter(
            plt.FuncFormatter(lambda val, pos: f"{val:.2f}"))
        _polish(axes_fi[i])
    fig_fi.tight_layout()
    fig_fi.savefig(OUT_APP / "MLSurvival_featimp_nestedCV.pdf", dpi=300)
    plt.close(fig_fi)
    print("Saved MLSurvival_featimp_nestedCV.pdf")

    # ── Assemble Figure 4 C + D ──────────────────────────────────────────────
    png_a4 = OUT_APP / "_fig4A_panel.png"
    png_b4 = OUT_APP / "_fig4B_panel.png"
    have_r4 = png_a4.exists() and png_b4.exists()

    if have_r4:
        fig4 = plt.figure(figsize=(18, 18))
        gs4  = gridspec.GridSpec(2, 2, figure=fig4,
                                 height_ratios=[1, 1.3],
                                 hspace=0.55, wspace=0.38)

        ax_a4 = fig4.add_subplot(gs4[0, 0])
        ax_a4.imshow(mpimg.imread(str(png_a4)))
        ax_a4.axis("off")
        ax_a4.text(-0.05, 1.04, "A", transform=ax_a4.transAxes,
                   fontsize=PANEL_SIZE, fontweight="bold", va="bottom", ha="left")

        ax_b4 = fig4.add_subplot(gs4[0, 1])
        ax_b4.imshow(mpimg.imread(str(png_b4)))
        ax_b4.axis("off")
        ax_b4.text(-0.05, 1.04, "B", transform=ax_b4.transAxes,
                   fontsize=PANEL_SIZE, fontweight="bold", va="bottom", ha="left")

        ax_c4 = fig4.add_subplot(gs4[1, 0])
        ax_d4 = fig4.add_subplot(gs4[1, 1])
    else:
        fig4 = plt.figure(figsize=(18, 10))
        gs4  = gridspec.GridSpec(1, 2, figure=fig4, wspace=0.38)
        ax_c4 = fig4.add_subplot(gs4[0, 0])
        ax_d4 = fig4.add_subplot(gs4[0, 1])

    # Panel C: concordance bar
    bars = ax_c4.bar(methods, means, yerr=errors, capsize=5,
                     color=colors, edgecolor=None, width=0.55)
    for bar, m, e in zip(bars, means, errors):
        ax_c4.text(bar.get_x() + bar.get_width() / 2, m + e + 0.02,
                   f"{m:.3f}", ha="center", va="bottom", fontsize=FONT_SIZE)
    ax_c4.set_ylabel("Concordance Index")
    ax_c4.set_ylim(0, 0.85)
    _polish(ax_c4)
    _panel_label(ax_c4, "c")

    # Panel D: permutation importances 2×2
    gs_d = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=ax_d4.get_subplotspec(),
                                            hspace=0.45, wspace=0.4)
    ax_d4.set_visible(False)
    for i, (model_name, df) in enumerate(perm_data.items()):
        ax_sub = fig4.add_subplot(gs_d[i // 2, i % 2])
        medians   = df.median().sort_values(ascending=False)
        df_sorted = df[medians.index]
        sns.boxplot(data=df_sorted, orient="h", ax=ax_sub,
                    color=colors[i], fliersize=3)
        ax_sub.set_title(titles[i], fontsize=FONT_SIZE)
        ax_sub.set_xlabel("Permutation Importance (\u0394 C-index)")
        ax_sub.set_ylabel("")
        _polish(ax_sub)
        if i == 0:
            _panel_label(ax_sub, "d", x=-0.20)

    out_name = "fig4_combined.pdf" if have_r4 else "MLSurvival_concordance_importances.pdf"
    fig4.savefig(OUT_APP / out_name, dpi=300, bbox_inches="tight")
    plt.close(fig4)
    print(f"Saved {out_name}")

    # ── Appendix: full-sample concordance ────────────────────────────────────
    rsf_fs = RandomSurvivalForest(random_state=42, n_jobs=-1,
                                  max_depth=20, n_estimators=200)
    rsf_fs.fit(X, y)
    gbsa_fs = GradientBoostingSurvivalAnalysis(random_state=42,
                                               n_estimators=200,
                                               learning_rate=0.1, max_depth=3)
    gbsa_fs.fit(X, y)
    cox_fs = CoxPHSurvivalAnalysis(); cox_fs.fit(X, y)
    lasso_fs = CoxnetSurvivalAnalysis(max_iter=10000, alphas=[0.1], l1_ratio=0.5)
    lasso_fs.fit(X, y)

    means_fs = [rsf_fs.score(X, y), gbsa_fs.score(X, y),
                cox_fs.score(X, y), lasso_fs.score(X, y)]
    perm_rsf_fs   = permutation_importance(rsf_fs,  X, y, n_repeats=10,
                                           random_state=42, scoring=concordance_scorer)
    perm_gbsa_fs  = permutation_importance(gbsa_fs, X, y, n_repeats=10,
                                           random_state=42, scoring=concordance_scorer)
    perm_coxph_fs = permutation_importance(cox_fs,  X, y, n_repeats=10,
                                           random_state=42, scoring=concordance_scorer)
    perm_lasso_fs = permutation_importance(lasso_fs,X, y, n_repeats=10,
                                           random_state=42, scoring=concordance_scorer)

    fig_a, ax_a = plt.subplots(figsize=(5, 6))
    bars = ax_a.bar(methods, means_fs, capsize=4, color=colors, edgecolor=None)
    for bar, m in zip(bars, means_fs):
        ax_a.text(bar.get_x() + bar.get_width() / 2, m + 0.01,
                  f"{m:.3f}", ha="center", va="bottom", fontsize=FONT_SIZE)
    ax_a.set_ylabel("Concordance Index")
    ax_a.set_ylim(0, 1.0)
    ax_a.axhline(0.5, color="black", linestyle="--", linewidth=1,
                 label="Random Guessing")
    ax_a.legend(frameon=False)
    _polish(ax_a)
    fig_a.tight_layout()
    fig_a.savefig(OUT_APP / "MLSurvival_CI_fullsample.pdf", dpi=300)
    plt.close(fig_a)

    perm_data_fs = {
        "RSF":      pd.DataFrame(perm_rsf_fs.importances.T,
                                 columns=list(cols_rename.values())),
        "GBSA":     pd.DataFrame(perm_gbsa_fs.importances.T,
                                 columns=list(cols_rename.values())),
        "CoxPH":    pd.DataFrame(perm_coxph_fs.importances.T,
                                 columns=list(cols_rename.values())),
        "LassoCox": pd.DataFrame(perm_lasso_fs.importances.T,
                                 columns=list(cols_rename.values())),
    }
    fig_b, axes_b = plt.subplots(2, 2, figsize=(14, 8), sharey=False)
    axes_b = axes_b.flatten()
    for i, (model_name, df) in enumerate(perm_data_fs.items()):
        medians   = df.median().sort_values(ascending=False)
        df_sorted = df[medians.index]
        sns.boxplot(data=df_sorted, orient="h", ax=axes_b[i],
                    color=colors[i], fliersize=3)
        axes_b[i].set_title(titles[i])
        axes_b[i].set_xlabel("Permutation Importance (\u0394 C-index)")
        axes_b[i].set_ylabel("")
        _polish(axes_b[i])
    fig_b.tight_layout()
    fig_b.savefig(OUT_APP / "MLSurvival_featimp_fullsample.pdf", dpi=300)
    plt.close(fig_b)
    print("Saved appendix ML figures")

    return methods, means, errors, perm_data


if __name__ == "__main__":
    print("=== ML Onset Classifier (Fig 2 C–E) ===")
    main_result = run_onset_classifier()

    print("\nAssembling Figure 2 …")
    assemble_fig2(main_result)

    print("\n=== ML Survival Analysis (Fig 4 C–D) ===")
    run_survival_ml()

    print("\nAll ML figures complete.")
