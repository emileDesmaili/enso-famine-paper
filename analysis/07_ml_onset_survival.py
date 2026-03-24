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
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
import seaborn as sns

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import GridSearchCV, cross_val_score, KFold
from sklearn.metrics import accuracy_score
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
LABEL_SIZE = 16
PANEL_SIZE = 22

mpl.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
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


def _panel_label(ax, letter, x=-0.12, y=1.05):
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
        "max_depth":         [1, 3, 5, 10],
        "n_estimators":      [3, 5, 10, 50],
        "learning_rate":     [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5],
        "subsample":         [0.7, 0.9],
        "max_features":      ["sqrt", 0.5],
        "min_samples_leaf":  [2, 5],
        "min_samples_split": [5, 10],
    }

    main_result = None   # returned for the "All_features" set

    for feat_set, suffix in zip(feature_combinations, save_suffixes):
        X = onset[feat_set]
        y = onset["Famine_start"]

        clf  = GradientBoostingClassifier(random_state=42)
        grid = GridSearchCV(clf, param_grid, cv=10, scoring="f1", n_jobs=-1)
        grid.fit(X, y)
        best = grid.best_estimator_
        print(f"[{suffix}] best params: {grid.best_params_}")

        pred  = best.predict(X)
        X_cf  = X.copy()
        if "nino34" in feat_set:
            mask = (onset["Famine_start"].values == 1) & (onset["nino34"].values > 0)
            X_cf.loc[mask, "nino34"] = 0
        pred_cf = best.predict(X_cf)

        df_chron = pd.DataFrame({
            "Region":                onset["Region"].values,
            "Year":                  onset["Year"].values,
            "Observed":              y.values,
            "Predicted":             pred,
            "Counterfactual":        pred_cf,
            "NINO34_actual":         onset["nino34"].values,
            "NINO34_counterfactual": X_cf["nino34"].values if "nino34" in feat_set else 0,
        })
        df_chron["Predicted"]      = ((df_chron["Observed"] == 1) &
                                      (df_chron["Predicted"] == 1)).astype(int)
        df_chron["Counterfactual"] = ((onset["Famine_start"].values == 1) &
                                      (df_chron["Counterfactual"] == 1)).astype(int)

        perm = permutation_importance(best, X, y, n_repeats=50,
                                      random_state=42, scoring="f1", n_jobs=-1)
        df_imp = pd.concat([
            pd.DataFrame({"Feature": f, "Importance": perm["importances"][i]})
            for i, f in enumerate(feat_set)
        ])
        df_imp["Feature"] = df_imp["Feature"].replace(FEAT_LABELS)

        cv_scores     = cross_val_score(best, X, y, cv=10, scoring="accuracy")
        in_sample_acc = accuracy_score(y, pred)

        # Save individual panel PDFs (for appendix reuse)
        _save_chronology_pdf(df_chron, onset, suffix)
        _save_importances_pdf(df_imp, suffix)
        _save_skill_pdf(in_sample_acc, cv_scores, suffix)

        if suffix == "All_features":
            main_result = dict(
                df_chron=df_chron,
                onset=onset,
                df_imp=df_imp,
                in_sample_acc=in_sample_acc,
                cv_scores=cv_scores,
            )

    return main_result


def _save_chronology_pdf(df_chron, onset_df, suffix):
    years      = sorted(df_chron["Year"].unique())
    year_min, year_max = min(years), max(years)
    colors     = {"Observed": "firebrick", "Predicted": "darksalmon",
                  "Counterfactual": "cornflowerblue"}

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
                ax1.barh(y=y_positions[region][scenario], width=1, left=yr,
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
    fig.savefig(OUT_MAIN / f"fig2C_MLOnset_{suffix}_chronology.pdf",
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved chronology for {suffix}")


def _save_importances_pdf(df_imp, suffix):
    med_order = (df_imp.groupby("Feature")["Importance"]
                 .median().sort_values(ascending=False).index)
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.boxplot(x="Importance", y="Feature", data=df_imp,
                order=med_order, palette="Reds_r", ax=ax)
    ax.set_xlabel("Permutation Importance")
    ax.set_ylabel("")
    _polish(ax)
    fig.tight_layout()
    fig.savefig(OUT_MAIN / f"fig2D_MLOnset_{suffix}_importances.pdf",
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
    sns.barplot(data=skill_df, y="Skill", x="Accuracy",
                palette={"In-sample": "firebrick", "Cross-Validation": "darksalmon"},
                ax=ax)
    ax.errorbar(x=cv_scores.mean(), y=1, xerr=cv_scores.std(),
                color="black", fmt="none", capsize=4, elinewidth=1.5)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Accuracy")
    _polish(ax)
    fig.tight_layout()
    fig.savefig(OUT_MAIN / f"fig2E_MLOnset_{suffix}_skillHor.pdf",
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved skill for {suffix}")


# ── Draw individual subpanels into an Axes (for combined figure) ───────────────

def _draw_C_chronology(ax_chron, ax_nino, df_chron):
    years      = sorted(df_chron["Year"].unique())
    year_min, year_max = min(years), max(years)
    colors     = {"Observed": "firebrick", "Predicted": "darksalmon",
                  "Counterfactual": "cornflowerblue"}

    y_positions = {}
    for i, region in enumerate(sorted(df_chron["Region"].unique())):
        base_y = i * 6
        y_positions[region] = {"Observed": base_y + 3,
                                "Predicted": base_y + 1.5,
                                "Counterfactual": base_y}
        rd = df_chron[df_chron["Region"] == region].sort_values("Year")
        for scenario in ["Observed", "Predicted", "Counterfactual"]:
            for yr in rd.loc[rd[scenario] == 1, "Year"].values:
                ax_chron.barh(y=y_positions[region][scenario], width=1, left=yr,
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
    med_order = (df_imp.groupby("Feature")["Importance"]
                 .median().sort_values(ascending=False).index)
    sns.boxplot(x="Importance", y="Feature", data=df_imp,
                order=med_order, palette="Reds_r", ax=ax)
    ax.set_xlabel("Permutation Importance")
    ax.set_ylabel("")
    _polish(ax)


def _draw_E_skill(ax, in_sample_acc, cv_scores):
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


# ── Assemble Figure 2  (A+B from R PNGs, C+D+E from Python) ──────────────────

def assemble_fig2(main_result):
    png_a = OUT_MAIN / "_fig2A_panel.png"
    png_b = OUT_MAIN / "_fig2B_panel.png"
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

    fig.savefig(OUT_MAIN / "fig2_combined.pdf", dpi=300, bbox_inches="tight")
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
        "RSF":      pd.DataFrame(perm_rsf.importances.T,
                                 columns=list(cols_rename.values())),
        "GBSA":     pd.DataFrame(perm_gbsa.importances.T,
                                 columns=list(cols_rename.values())),
        "CoxPH":    pd.DataFrame(perm_coxph.importances.T,
                                 columns=list(cols_rename.values())),
        "LassoCox": pd.DataFrame(perm_lasso.importances.T,
                                 columns=list(cols_rename.values())),
    }

    # ── Assemble Figure 4 C + D ──────────────────────────────────────────────
    png_a4 = OUT_MAIN / "_fig4A_panel.png"
    png_b4 = OUT_MAIN / "_fig4B_panel.png"
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
    _panel_label(ax_c4, "C")

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
            _panel_label(ax_sub, "D", x=-0.20)

    out_name = "fig4_combined.pdf" if have_r4 else "fig4CD_combined.pdf"
    fig4.savefig(OUT_MAIN / out_name, dpi=300, bbox_inches="tight")
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
    fig_a.savefig(OUT_APP / "figA_ML_CI_fullsample.pdf", dpi=300)
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
    fig_b.savefig(OUT_APP / "figA_ML_featimp_fullsample.pdf", dpi=300)
    plt.close(fig_b)
    print("Saved appendix ML figures")


if __name__ == "__main__":
    print("=== ML Onset Classifier (Fig 2 C–E) ===")
    main_result = run_onset_classifier()

    print("\nAssembling Figure 2 …")
    assemble_fig2(main_result)

    print("\n=== ML Survival Analysis (Fig 4 C–D) ===")
    run_survival_ml()

    print("\nAll ML figures complete.")
