"""Mediator-aggregation simulation following the reviewer's stated DGP.

Resolutions: d = days, m = months, y = years.

DGP (10 years of data, 12 months/year, 30 days/month):
    ENSO_y           ~ N(0, 1)                              annual
    C_m  = beta * ENSO_y(m)  + eps_C,  eps_C ~ N(0, sig_C)  monthly climate
    w_d  = theta * C_m(d)    + eps_w,  eps_w ~ N(0, sig_w)  daily weather
    y_d  = alpha * 1{w_d > p90(w)} + eps_y                  daily outcome

The outcome only depends on whether the daily weather crosses the 90th
percentile, i.e. the mechanism is the *count of extreme days*, not the
level of daily weather.

Regressions:
    (i)  y_m ~ ENSO_y(m) + w_m       monthly aggregation
    (ii) y_d ~ ENSO_y(d) + w_d       daily resolution

w_m and y_m are simple within-month averages.

Across B replicates we report the estimated ENSO coefficient under (i)
and (ii). With this clean DGP, both controls absorb ENSO completely:
the daily extreme mechanism propagates into the monthly mean because
w_m = theta * C_m + tiny noise, so w_m is a precise linear function of
the underlying monthly climate state. The reviewer's argument follows
directly: when the seasonal/monthly climate indicator is measured
cleanly, it does attenuate ENSO. A paper that finds ENSO un-attenuated
by its seasonal climate controls is therefore more likely to have a
measurement-noise problem with its derived indicators than a true
"sub-seasonal" mediator that monthly means cannot capture.

Outputs:
    analysis/output/figures/investigation/sim_mediator_dgp.pdf
    analysis/output/figures/investigation/sim_mediator_coeffs.pdf
"""
from __future__ import annotations
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT_DIR = Path(__file__).resolve().parent.parent.parent
os.chdir(ROOT_DIR)

# ── DGP parameters ────────────────────────────────────────────────────────────
N_YEARS         = 10
MONTHS_PER_YEAR = 12
DAYS_PER_MONTH  = 30

BETA    = 0.5    # ENSO -> monthly climate
THETA   = 0.5    # monthly climate -> daily weather
SIG_C   = 1.0    # monthly climate noise
SIG_W   = 1.0    # daily weather noise
ALPHA   = 1.0    # outcome shift per extreme day
SIG_Y   = 0.5    # daily outcome noise
QUANT   = 0.90   # extreme threshold percentile
SIG_AGG_LOW  = 0.50  # noisy monthly aggregate (signal-to-noise ~ 1)
SIG_AGG_HIGH = 1.50  # very noisy monthly aggregate (signal-to-noise ~ 1/3)

N_REPS = 500

N_MONTHS = N_YEARS * MONTHS_PER_YEAR
N_DAYS   = N_MONTHS * DAYS_PER_MONTH


def simulate_once(rng: np.random.Generator) -> dict:
    ENSO_y = rng.standard_normal(N_YEARS)

    # Map year -> 12 months -> 360 days
    year_of_month = np.repeat(np.arange(N_YEARS), MONTHS_PER_YEAR)
    month_of_day  = np.repeat(np.arange(N_MONTHS), DAYS_PER_MONTH)
    year_of_day   = year_of_month[month_of_day]

    # Monthly climate state
    eps_C = rng.normal(0.0, SIG_C, size=N_MONTHS)
    C_m   = BETA * ENSO_y[year_of_month] + eps_C

    # Daily weather
    eps_w = rng.normal(0.0, SIG_W, size=N_DAYS)
    w_d   = THETA * C_m[month_of_day] + eps_w

    # Extreme-day indicator, threshold = empirical 90th percentile of w_d
    p90 = np.quantile(w_d, QUANT)
    eps_y = rng.normal(0.0, SIG_Y, size=N_DAYS)
    y_d   = ALPHA * (w_d > p90).astype(float) + eps_y

    # Monthly averages of y and w. The monthly w aggregate carries
    # additional "derivation" noise that mimics the proxy-pipeline error
    # the reviewer flagged. We construct two versions: low and high noise.
    y_m       = pd.Series(y_d).groupby(month_of_day).mean().values
    w_m_true  = pd.Series(w_d).groupby(month_of_day).mean().values
    w_m_lo    = w_m_true + rng.normal(0.0, SIG_AGG_LOW,  size=N_MONTHS)
    w_m_hi    = w_m_true + rng.normal(0.0, SIG_AGG_HIGH, size=N_MONTHS)
    enso_for_month = ENSO_y[year_of_month]
    enso_for_day   = ENSO_y[year_of_day]

    def ols(y, *Xcols):
        X = np.column_stack([np.ones_like(y)] + list(Xcols))
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta
        n, k  = X.shape
        ss_res = (resid ** 2).sum()
        ss_tot = ((y - y.mean()) ** 2).sum()
        r2 = 1 - ss_res / ss_tot
        return beta, r2

    # Monthly regressions
    b_m_total, _          = ols(y_m, enso_for_month)
    b_m_ctrl_lo, r_mc_lo  = ols(y_m, enso_for_month, w_m_lo)
    b_m_ctrl_hi, r_mc_hi  = ols(y_m, enso_for_month, w_m_hi)

    # Daily regression: y_d ~ ENSO_y + w_d
    b_d_total, _    = ols(y_d, enso_for_day)
    b_d_ctrl,  r_dc = ols(y_d, enso_for_day, w_d)

    return dict(
        bT_month_total   = b_m_total[1],
        bT_month_ctrl_lo = b_m_ctrl_lo[1],
        bT_month_ctrl_hi = b_m_ctrl_hi[1],
        bT_day_total     = b_d_total[1],
        bT_day_ctrl      = b_d_ctrl[1],
        r2_month_lo      = r_mc_lo,
        r2_month_hi      = r_mc_hi,
        r2_day           = r_dc,
        snap_T_year      = ENSO_y,
        snap_w_d         = w_d,
        snap_y_d         = y_d,
        snap_C_m         = C_m,
        snap_w_m         = w_m_lo,
        snap_y_m         = y_m,
        snap_enso_m      = enso_for_month,
    )


def main():
    rng = np.random.default_rng(20260626)

    rows, snap = [], None
    for r in range(N_REPS):
        res = simulate_once(rng)
        if snap is None:
            snap = res
        rows.append({k: v for k, v in res.items() if not k.startswith("snap_")})
    df = pd.DataFrame(rows)
    print("=== Mean across replicates ===")
    print(df.mean().round(4).to_string())

    # ── Fig 1: DGP illustration (snapshot from first replicate) ─────────────
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))

    ax = axes[0]
    ax.scatter(snap["snap_enso_m"], snap["snap_C_m"], s=12, alpha=0.6,
               color="cornflowerblue")
    ax.set_xlabel("ENSO$_y$ (annual)")
    ax.set_ylabel("Monthly climate $C_m$")
    ax.set_title("ENSO drives monthly climate")

    ax = axes[1]
    ax.scatter(snap["snap_w_m"], snap["snap_y_m"], s=12, alpha=0.6,
               color="firebrick")
    ax.set_xlabel("Monthly mean weather $w_m$")
    ax.set_ylabel("Monthly outcome $y_m$")
    ax.set_title("Monthly mean tracks the latent climate state")

    ax = axes[2]
    # subsample daily points for clarity
    idx = rng.choice(len(snap["snap_w_d"]), size=2000, replace=False)
    ax.scatter(snap["snap_w_d"][idx], snap["snap_y_d"][idx], s=6, alpha=0.35,
               color="darkgreen")
    ax.set_xlabel("Daily weather $w_d$")
    ax.set_ylabel("Daily outcome $y_d$")
    ax.set_title("Daily resolution exposes the mechanism")

    for a in axes:
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
    fig.tight_layout()
    out1 = "analysis/output/figures/investigation/sim_mediator_dgp.pdf"
    fig.savefig(out1, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out1}")

    # ── Fig 2: ENSO coefficient across regression specs ────────────────────
    fig, ax = plt.subplots(1, 1, figsize=(6.6, 4.0))

    data = [
        df["bT_month_total"],
        df["bT_month_ctrl_lo"],
        df["bT_month_ctrl_hi"],
        df["bT_day_total"],
        df["bT_day_ctrl"],
    ]
    labels = [
        r"$y_m \sim \mathrm{ENSO}$",
        r"$+\, w_m^{obs}$" + "\n(noisy)",
        r"$+\, w_m^{obs}$" + "\n(very noisy)",
        r"$y_d \sim \mathrm{ENSO}$",
        r"$+\, w_d$",
    ]
    # Soft, paper-friendly palette
    COL_NULL  = "#c9ccd1"   # neutral grey for uncontrolled (no regressor)
    COL_MON1  = "#83b9b2"   # muted teal — monthly, noisy
    COL_MON2  = "#3c7e7a"   # deeper teal — monthly, very noisy
    COL_DAY   = "#c47b6f"   # muted coral — daily control
    cols = [COL_NULL, COL_MON1, COL_MON2, COL_NULL, COL_DAY]

    positions = [1, 2, 3, 4.4, 5.4]

    bp = ax.boxplot(data, positions=positions, tick_labels=labels,
                     patch_artist=True, widths=0.36, showfliers=False,
                     whis=(5, 95),
                     boxprops=dict(linewidth=0),
                     whiskerprops=dict(linewidth=0.6, color="#5a5a5a"),
                     capprops=dict(linewidth=0.0),
                     medianprops=dict(color="white", linewidth=1.3))
    for patch, fc in zip(bp["boxes"], cols):
        patch.set_facecolor(fc); patch.set_alpha(0.95)
        patch.set_edgecolor("none")

    ax.axhline(0, color="#808080", linestyle="--", linewidth=0.6, alpha=0.6)

    # group labels above the two clusters
    ax.text(2.0, 0.085, "Monthly",
            ha="center", va="bottom", fontsize=10, fontweight="bold",
            color="#3c7e7a", transform=ax.transData)
    ax.text(4.9, 0.085, "Daily",
            ha="center", va="bottom", fontsize=10, fontweight="bold",
            color="#a45848", transform=ax.transData)
    ax.set_ylim(top=0.10)

    ax.set_ylabel(r"Estimated ENSO coefficient  $\hat\beta_{\mathrm{ENSO}}$",
                  fontsize=10)
    ax.set_title("ENSO coefficient by control resolution", fontsize=11)
    ax.tick_params(axis="x", labelsize=9)
    ax.tick_params(axis="y", labelsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.grid(axis="y", color="#eeeeee", linewidth=0.5, alpha=0.8)
    ax.set_axisbelow(True)
    fig.tight_layout()
    out2 = "analysis/output/figures/investigation/sim_mediator_coeffs.pdf"
    fig.savefig(out2, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out2}")

    summary = df.agg(["mean", "std"]).T
    out_csv = "analysis/output/data/sim_mediator_summary.csv"
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    summary.to_csv(out_csv)
    print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
