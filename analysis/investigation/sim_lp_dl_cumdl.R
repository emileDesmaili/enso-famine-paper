###############################################################################
# Simulation: compare LP, individual DL, and cumulative DL estimators of the
# impulse response on a known DGP.
#
# DGP
# ---
#   x_t = ρ · x_{t-1} + ε_t,          ε_t  ~ N(0, 1)     (AR(1) shock)
#   y_t = Σ_{k=0..K} β_k · x_{t-k} + η_t,   η_t ~ N(0, σ_y²)
#   β   = (0.5, 1.0, 0.7, 0.3),  K = 3
#
# Population impulse responses a shock at t = 0 can trace out:
#   • Direct-x  β_h   : if x_0 = 1 and past/future x = 0, y_h = β_h.
#   • Innovation IRF  : if ε_0 = 1, x propagates via AR, so
#                       y_h = Σ_{k=0..min(h,K)} β_k · ρ^{h-k}.
#   • Cumulative-x    : Σ_{k=0..h} β_k  ("long-run" multiplier).
#   (For ρ = 0 the first two collapse onto the same curve.)
#
# Estimators (all fit by OLS with Newey-West SEs)
# ------------------------------------------------
#   LP  (per horizon h):
#       y_{t+h} − y_{t-1} = δ_h · x_t
#                           + γ_{1,h} x_{t-1} + γ_{2,h} x_{t-2}
#                           + γ_{3,h} x_{t-3} + ε_{t,h}
#     → plot δ_h vs h.  Under AR(ρ) shock, δ_h ≈ innovation IRF at h.
#
#   DL  (single fit, lags 0..H):
#       y_t = Σ_{k=0..H} β̂_k · x_{t-k} + u_t
#     → plot individual β̂_h vs h.                        Recovers direct β_h.
#     → cumulative sum  Σ_{k=0..h} β̂_k  vs h.            Recovers Σ β_k.
#
# Output: analysis/output/figures/investigation/sim_lp_dl_cumdl.pdf
###############################################################################
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(fixest); library(ggplot2)
})
setFixest_notes(FALSE)

.script_dir <- function() {
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("--file=", "", file_arg)))
  } else {
    tryCatch(dirname(normalizePath(sys.frame(1)$ofile)),
             error = function(e) getwd())
  }
}
SCRIPT_DIR   <- if (exists("SCRIPT_PATH_HINT")) dirname(SCRIPT_PATH_HINT) else .script_dir()
ANALYSIS_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
ROOT_DIR     <- normalizePath(file.path(ANALYSIS_DIR, ".."), mustWork = FALSE)
setwd(ROOT_DIR)

OUT_FIG <- "analysis/output/figures/investigation"
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# ── DGP settings ──────────────────────────────────────────────────────────────
set.seed(2026)
Tn        <- 1000                                # length of the time series
beta_true <- c(0.5, 1.0, 0.7, 0.3)                # β_0..β_K
K_true    <- length(beta_true) - 1
rho       <- 0                                    # IID shock (no persistence)
sigma_y   <- 1.5                                  # noise sd on y

HOR       <- 10                                   # estimation horizon

# ── Simulate ──────────────────────────────────────────────────────────────────
eps <- rnorm(Tn)
x   <- numeric(Tn)
for (t in 2:Tn) x[t] <- rho * x[t - 1] + eps[t]

eta <- rnorm(Tn, sd = sigma_y)
y   <- rep(NA_real_, Tn)
for (t in (K_true + 1):Tn) {
  y[t] <- sum(beta_true * x[t - (0:K_true)]) + eta[t]
}

df <- data.frame(t = seq_len(Tn), x = x, y = y, id = 1L) %>% drop_na()

# ── Population impulse responses ──────────────────────────────────────────────
# 1) Direct-x β_h (shock only at t = 0)
true_direct <- c(beta_true, rep(0, HOR - K_true))
# 2) Innovation impulse: shock ε_0 = 1, propagate via AR
x_ar <- rho ^ (0:HOR)                                     # x_h under ε_0 = 1
true_innov <- vapply(0:HOR, function(h) {
  sum(beta_true[0:min(h, K_true) + 1] *
        x_ar[h - (0:min(h, K_true)) + 1])
}, numeric(1))
# 3) Cumulative-x Σ β_k
true_cum <- cumsum(true_direct)

truth_df <- bind_rows(
  data.frame(horizon = 0:HOR, mean = true_direct, series = "True β_h (direct)"),
  data.frame(horizon = 0:HOR, mean = true_innov,  series = "True innovation IRF"),
  data.frame(horizon = 0:HOR, mean = true_cum,    series = "True Σ β_k (cumulative)")
)

# ── LP with 3 past-shock controls ─────────────────────────────────────────────
lp_df <- bind_rows(lapply(0:HOR, function(h) {
  fml <- as.formula(paste0("f(y,", h, ") - l(y,1) ~ x + l(x, 1:3)"))
  m <- feols(fml, data = df, panel.id = c("id", "t"),
             vcov = NW(3) ~ t)
  b  <- as.numeric(coef(m)["x"])
  se <- as.numeric(sqrt(vcov(m)["x", "x"]))
  data.frame(horizon = h, mean = b,
             lo = b - 1.96 * se, hi = b + 1.96 * se)
})) %>% mutate(estimator = "Local projection")

# ── DL on levels: single fit with lags 0..HOR ───────────────────────────────
dl_fml <- as.formula(paste0("y ~ l(x, 0:", HOR, ")"))
m_dl   <- feols(dl_fml, data = df, panel.id = c("id", "t"),
                vcov = NW(HOR) ~ t)
lag_names <- paste0("l(x, ", 0:HOR, ")")
b_dl   <- as.numeric(coef(m_dl)[lag_names])
V_dl   <- vcov(m_dl)[lag_names, lag_names]
se_dl  <- sqrt(diag(V_dl))

dl_ind_df <- data.frame(horizon = 0:HOR,
                        mean = b_dl,
                        lo   = b_dl - 1.96 * se_dl,
                        hi   = b_dl + 1.96 * se_dl,
                        estimator = "DL levels (individual β_h)")

dl_cum_df <- bind_rows(lapply(0:HOR, function(h) {
  L <- c(rep(1, h + 1), rep(0, HOR - h))
  est <- sum(L * b_dl)
  se  <- sqrt(as.numeric(t(L) %*% V_dl %*% L))
  data.frame(horizon = h, mean = est,
             lo = est - 1.96 * se, hi = est + 1.96 * se)
})) %>% mutate(estimator = "DL levels (cumulative Σ β_k)")

# ── DL on first differences: d(y) ~ l(x, 0:H) ────────────────────────────────
# Cumulative sum of the differenced γ_k recovers β_h (see notes).
dld_fml <- as.formula(paste0("d(y) ~ l(x, 0:", HOR, ")"))
m_dld   <- feols(dld_fml, data = df, panel.id = c("id", "t"),
                 vcov = NW(HOR) ~ t)
b_dld   <- as.numeric(coef(m_dld)[lag_names])
V_dld   <- vcov(m_dld)[lag_names, lag_names]
se_dld  <- sqrt(diag(V_dld))

dld_ind_df <- data.frame(horizon = 0:HOR,
                         mean = b_dld,
                         lo   = b_dld - 1.96 * se_dld,
                         hi   = b_dld + 1.96 * se_dld,
                         estimator = "DL d(y) (individual g_h)")

dld_cum_df <- bind_rows(lapply(0:HOR, function(h) {
  L <- c(rep(1, h + 1), rep(0, HOR - h))
  est <- sum(L * b_dld)
  se  <- sqrt(as.numeric(t(L) %*% V_dld %*% L))
  data.frame(horizon = h, mean = est,
             lo = est - 1.96 * se, hi = est + 1.96 * se)
})) %>% mutate(estimator = "DL d(y) (cumulative sum g_k)")

# ── Plot ──────────────────────────────────────────────────────────────────────
cols <- c("Local projection"               = "#0072B2",
          "DL levels (individual β_h)"     = "#009E73",
          "DL levels (cumulative Σ β_k)"   = "#D55E00",
          "DL d(y) (individual g_h)"         = "#56B4E9",
          "DL d(y) (cumulative sum g_k)"       = "#CC79A7")

truth_cols <- c("True β_h (direct)"       = "black",
                "True innovation IRF"     = "grey35",
                "True Σ β_k (cumulative)" = "grey55")

p <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey40",
             linetype = "dashed", linewidth = 0.3) +
  # LP: line + ribbon
  geom_ribbon(data = lp_df,
              aes(horizon, ymin = lo, ymax = hi, fill = estimator),
              alpha = 0.18) +
  geom_line(data = lp_df,
            aes(horizon, mean, colour = estimator), linewidth = 0.9) +
  # DL on levels — individual (points+error bars) and cumulative (dashed)
  geom_errorbar(data = dl_ind_df,
                aes(x = horizon, ymin = lo, ymax = hi, colour = estimator),
                width = 0.18, linewidth = 0.5) +
  geom_point(data = dl_ind_df,
             aes(horizon, mean, colour = estimator), size = 2.2) +
  geom_line(data = dl_cum_df,
            aes(horizon, mean, colour = estimator),
            linetype = "dashed", linewidth = 0.9) +
  # DL on Δy — individual (open squares) and cumulative (dotdash)
  geom_errorbar(data = dld_ind_df,
                aes(x = horizon, ymin = lo, ymax = hi, colour = estimator),
                width = 0.18, linewidth = 0.5) +
  geom_point(data = dld_ind_df,
             aes(horizon, mean, colour = estimator),
             shape = 22, size = 2.2, fill = "white") +
  geom_line(data = dld_cum_df,
            aes(horizon, mean, colour = estimator),
            linetype = "dotdash", linewidth = 0.9) +
  # Truth references
  geom_line(data = truth_df,
            aes(horizon, mean, group = series, linetype = series),
            colour = "black", linewidth = 0.45) +
  scale_x_continuous(breaks = 0:HOR) +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_linetype_manual(values = c(
    "True β_h (direct)"       = "dotted",
    "True innovation IRF"     = "twodash",
    "True Σ β_k (cumulative)" = "longdash"
  )) +
  labs(x = "Horizon (years)",
       y = "Response to unit shock in x",
       colour = NULL, fill = NULL, linetype = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.box       = "vertical",
        legend.spacing.y = unit(-2, "pt"))

out_path <- file.path(OUT_FIG, "sim_lp_dl_cumdl.pdf")
ggsave(out_path, p, width = 9, height = 5)
cat("Saved", out_path, "\n")
