###############################################################################
# Investigation figure: LP IRF vs distributed-lag cumulative response
# for the three main outcomes — Wheat/Rye grain harvests, grain prices,
# fish prices. LP specs mirror the paper's main figures:
#   Grain harvest → fig3_harvest.R  run_lp_GT3
#   Grain price   → fig5_prices.R   run_lp_price   ("All regions")
#   Fish price    → fig5_prices.R   run_lp_fish_devcoded  (Cod+Herring avg)
#
# LP uses raw NINO3.4 with l(nino34, 1:3) as controls (paper-consistent).
# DL uses raw NINO3.4 as the shock with a first-differenced LHS
#   d(y) = y_t − y_{t-1}  ~ l(nino34, 0:H) | FE
# so cumulative sums Σ_{k=0}^{h} γ_k of the DL lag coefficients recover the
# level impulse response β_h. Plotted as a dashed line.
#
# All estimators use Driscoll-Kraay SEs (fixest DK ~ Year).
#
# Output: analysis/output/figures/investigation/lp_vs_dl_irf.pdf
###############################################################################
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(fixest)
  library(ggplot2); library(patchwork); library(stringr)
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

source("emileRegs.R")  # lp_panel()

OUT_FIG <- "analysis/output/figures/investigation"
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

HOR_LP_WR  <- 3    # horizons for grain harvest
HOR_LP_P   <- 10   # horizons for grain price / fish price

# DL lag grids per outcome — each entry is a *separate* DL fit
LAGS_WR <- c(3)                        # 1 fit for grain harvest
LAGS_P  <- c(3, 5, 8, 9)               # 4 fits for the two price panels

# ── Data ──────────────────────────────────────────────────────────────────────
yield_v2   <- read.csv("processed data/yield_ljungqvist_v2.csv") %>%
  mutate(decade = floor(Year / 10) * 10)
WR         <- yield_v2 %>% filter(Grain %in% c("Wheat", "Rye"))

fishprice  <- read.csv("processed data/fishprice_enso.csv") %>%
  mutate(Decade = floor(Year / 10) * 10,
         Year2  = Year^2)

grainprice <- read.csv("processed data/price_2023_enso.csv") %>%
  mutate(Decade = floor(Year / 10) * 10) %>%
  drop_na() %>% distinct() %>% arrange(Year)

# ── LP main specs ─────────────────────────────────────────────────────────────
# Grain harvest LP: run_lp_GT3() from fig3_harvest.R
irf_lp_wr <- lp_panel(
  data         = WR,
  outcome      = "logyield",
  main_var     = "nino34",
  controls     = "i(decade) + l(nino34, 1:3)",
  horizon      = HOR_LP_WR,
  fe           = "VarLocationGrain[Year]",
  panel_id     = c("VarLocationGrain", "Year"),
  vcov_formula = DK ~ Year
) %>% mutate(outcome = "Grain harvest")

# Grain price LP: run_lp_price() "All regions" spec from fig5_prices.R
irf_lp_gp <- lp_panel(
  data         = grainprice,
  outcome      = "logprice",
  main_var     = "nino34",
  controls     = "i(Decade) + l(nino34, 1:3)",
  horizon      = HOR_LP_P,
  fe           = "Location",
  panel_id     = c("Location", "Year"),
  vcov_formula = DK ~ Year
) %>% mutate(outcome = "Grain price")

# Fish price LP: deviation-coded Cod/Herring average from
# run_lp_fish_devcoded() in fig5_prices.R. lp_panel() does not support the
# nino34:cod_dev interaction natively, so implement inline.
irf_lp_fp <- {
  d <- fishprice %>%
    mutate(Decade  = as.factor(Decade),
           cod_dev = ifelse(Species == "Cod", 0.5, -0.5))
  bind_rows(lapply(0:HOR_LP_P, function(h) {
    fml <- as.formula(paste0(
      "f(logprice,", h, ") - l(logprice,1) ~ ",
      "nino34 + nino34:cod_dev + i(Decade) + l(nino34, 1:3) | LocationSpecies"
    ))
    m <- feols(fml, data = d,
               panel.id = c("LocationSpecies", "Year"),
               vcov     = DK ~ Year)
    b  <- as.numeric(coef(m)["nino34"])
    se <- sqrt(vcov(m)["nino34", "nino34"])
    data.frame(horizon = h, irf_mean = b, se = se,
               irf_down = b - 1.96 * se,
               irf_up   = b + 1.96 * se)
  }))
} %>% mutate(outcome = "Fish price")

irf_lp <- bind_rows(irf_lp_wr, irf_lp_gp, irf_lp_fp) %>%
  transmute(outcome, horizon,
            mean = irf_mean, lo = irf_down, hi = irf_up,
            estimator = "Local projection")

# ── DL cumulative coefficients using raw NINO3.4 ──────────────────────────────
# Fit  d(y) ~ l(nino34, 0:H) + f(nino34, 1:LEADS) + controls | FE
# with LEADS leads of the shock as pre-treatment controls.  Cumulative sums
# Σ_{k=0}^{h} γ_k of the lag coefficients recover the level impulse
# response β_h.
DL_SHOCK <- "nino34"
DL_LEADS <- 2

dl_cumulative <- function(data, outcome, main_var, hor, fe, panel_id,
                          extra_controls = NULL, extra_terms = NULL) {
  ctrl_parts <- c(paste0("f(", main_var, ", 1:", DL_LEADS, ")"),
                  extra_terms, extra_controls)
  ctrl_str <- paste("+", paste(ctrl_parts, collapse = " + "))
  fml <- as.formula(paste0(
    "d(", outcome, ") ~ l(", main_var, ", 0:", hor, ")", ctrl_str,
    " | ", fe
  ))
  m <- feols(fml, data = data, panel.id = panel_id, vcov = DK ~ Year)

  lag_names <- paste0("l(", main_var, ", ", 0:hor, ")")
  b <- coef(m)[lag_names]
  V <- vcov(m)[lag_names, lag_names]

  bind_rows(lapply(0:hor, function(h) {
    L <- c(rep(1, h + 1), rep(0, hor - h))
    est <- sum(L * b)
    se  <- sqrt(as.numeric(t(L) %*% V %*% L))
    data.frame(horizon = h, mean = est,
               lo = est - 1.96 * se, hi = est + 1.96 * se)
  }))
}

# One separate DL fit per lag count in the grid ─────────────────────────────
run_dl_wr <- function(lags) {
  dl_cumulative(
    WR, "logyield", DL_SHOCK, lags,
    fe = "VarLocationGrain[Year]",
    panel_id = c("VarLocationGrain", "Year"),
    extra_controls = "i(decade)"
  ) %>% mutate(outcome = "Grain harvest", dl_lags = lags)
}

run_dl_gp <- function(lags) {
  dl_cumulative(
    grainprice, "logprice", DL_SHOCK, lags,
    fe = "Location",
    panel_id = c("Location", "Year"),
    extra_controls = "i(Decade)"
  ) %>% mutate(outcome = "Grain price", dl_lags = lags)
}

run_dl_fp <- function(lags) {
  d <- fishprice %>%
    mutate(Decade  = as.factor(Decade),
           cod_dev = ifelse(Species == "Cod", 0.5, -0.5))
  dl_cumulative(
    d, "logprice", DL_SHOCK, lags,
    fe = "LocationSpecies",
    panel_id = c("LocationSpecies", "Year"),
    extra_controls = "i(Decade)",
    extra_terms = c(
      paste0("l(", DL_SHOCK, ", 0:", lags,     "):cod_dev"),
      paste0("f(", DL_SHOCK, ", 1:", DL_LEADS, "):cod_dev")
    )
  ) %>% mutate(outcome = "Fish price", dl_lags = lags)
}

dl_all <- bind_rows(
  bind_rows(lapply(LAGS_WR, run_dl_wr)),
  bind_rows(lapply(LAGS_P,  run_dl_gp)),
  bind_rows(lapply(LAGS_P,  run_dl_fp))
) %>%
  mutate(dl_lags_lbl = paste0("DL (", dl_lags, " lags)"))

# ── Plot ──────────────────────────────────────────────────────────────────────
lag_levels <- sort(unique(dl_all$dl_lags))
lag_labels <- paste0("DL (", lag_levels, " lags)")

curve_levels <- c("Local projection", lag_labels)

lp_curves <- irf_lp %>%
  transmute(outcome, horizon, mean, lo, hi,
            curve = factor("Local projection", levels = curve_levels))
dl_curves <- dl_all %>%
  transmute(outcome, horizon, mean, lo, hi,
            curve = factor(dl_lags_lbl, levels = curve_levels))

# Sequential palette for DL curves (light → dark orange); LP stays blue.
dl_palette <- c("#F16913", "#D94801", "#A63603", "#7F2704",
                "#54180C", "#2E0000")[seq_along(lag_levels)]
names(dl_palette) <- lag_labels
cols <- c("Local projection" = "#0072B2", dl_palette)

outcome_order <- c("Grain harvest", "Grain price", "Fish price")
lp_curves <- lp_curves %>% mutate(outcome = factor(outcome, levels = outcome_order))
dl_curves <- dl_curves %>% mutate(outcome = factor(outcome, levels = outcome_order))

p <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey40",
             linetype = "dashed", linewidth = 0.3) +
  geom_ribbon(data = lp_curves,
              aes(horizon, ymin = lo, ymax = hi),
              fill = cols[["Local projection"]], alpha = 0.18) +
  geom_line(data = lp_curves,
            aes(horizon, mean, colour = curve), linewidth = 0.9) +
  geom_line(data = dl_curves,
            aes(horizon, mean, colour = curve, group = curve),
            linetype = "dashed", linewidth = 0.75) +
  facet_wrap(~outcome, scales = "free", ncol = 3,
             labeller = as_labeller(c(
               "Grain harvest" = "a   Grain harvest",
               "Grain price"   = "b   Grain price",
               "Fish price"    = "c   Fish price"
             ))) +
  scale_x_continuous(breaks = 0:HOR_LP_P) +
  scale_colour_manual(values = cols, drop = FALSE) +
  labs(x = "Horizon (years)", y = "Response to +1 NINO3.4 shock",
       colour = NULL) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 12,
                                        hjust = 0),
        legend.position  = "bottom",
        legend.text      = element_text(size = 9)) +
  guides(colour = guide_legend(nrow = 1))

out_path <- file.path(OUT_FIG, "lp_vs_dl_irf.pdf")
ggsave(out_path, p, width = 10, height = 3.8)
cat("Saved", out_path, "\n")
