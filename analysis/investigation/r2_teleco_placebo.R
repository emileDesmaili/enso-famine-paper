###############################################################################
# Two-panel figure for the WR Teleconnected IRF:
#   Panel (a) — point + 95% CI for the main spec and three alt specs
#               (controls, 5 ENSO lags, quadratic trend). AMJJ partition
#               removed at the user's request.
#   Panel (b) — 1000-shuffle placebo IRF: 95% null envelope (2.5–97.5
#               percentiles of the null distribution) vs the main point.
#
# Output:
#   analysis/output/figures/investigation/fig3c_robust_and_placebo.pdf
###############################################################################
suppressPackageStartupMessages({
  library(dplyr); library(fixest); library(ggplot2); library(scales)
  library(patchwork)
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
SCRIPT_DIR <- if (exists("SCRIPT_PATH_HINT")) dirname(SCRIPT_PATH_HINT) else .script_dir()
ANALYSIS_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
ROOT_DIR     <- normalizePath(file.path(ANALYSIS_DIR, ".."), mustWork = FALSE)
setwd(ROOT_DIR)

OUT_FIG <- "analysis/output/figures/investigation"
OUT_DAT <- "analysis/output/data"
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

yield_v2 <- read.csv("processed data/yield_ljungqvist_v2.csv") |>
  mutate(decade = floor(Year / 10) * 10,
         Year2  = Year^2)
WR <- yield_v2 |> filter(Grain %in% c("Wheat", "Rye"))

HORIZONS <- 0:3

fit_teleco <- function(df, h, fe_str = "VarLocationGrain[Year]",
                       n_lags = 3, extra_controls = character()) {
  ctrl <- c("i(decade)", paste0("l(nino34, 1:", n_lags, ")"), extra_controls)
  fml <- as.formula(paste0(
    "f(logyield,", h, ") - l(logyield,1) ~ ",
    "i(teleco_PDSI_10, nino34) + ", paste(ctrl, collapse = " + "),
    " | ", fe_str))
  m <- tryCatch(
    feols(fml, data = df,
          panel.id = c("VarLocationGrain", "Year"),
          vcov     = DK ~ Year),
    error = function(e) NULL)
  if (is.null(m)) return(list(b = NA_real_, se = NA_real_))
  cn <- names(coef(m))
  nm <- grep("teleco_PDSI_10::1:nino34", cn, value = TRUE)
  if (length(nm) == 0) return(list(b = NA_real_, se = NA_real_))
  list(b  = as.numeric(unname(coef(m)[nm])),
       se = as.numeric(sqrt(vcov(m)[nm, nm])))
}

run_spec <- function(label, ...) {
  bind_rows(lapply(HORIZONS, function(h) {
    r <- fit_teleco(WR, h, ...)
    data.frame(spec = label, horizon = as.integer(h),
               irf_mean = r$b, se = r$se,
               irf_up   = r$b + 1.96 * r$se,
               irf_down = r$b - 1.96 * r$se)
  }))
}

# ── Panel (a): four specs ────────────────────────────────────────────────────
main_irf <- run_spec("Main")
ctrl_irf <- run_spec("+ controls (JSL, NAO, wars, deaths)",
                      extra_controls = c("JSL", "NAO_cal",
                                          "ongoing_wars", "Deaths"))
lag5_irf <- run_spec("5 ENSO lags", n_lags = 5)
quad_irf <- run_spec("Quadratic record trend",
                      fe_str = "VarLocationGrain[Year, Year2]")

all_irf <- bind_rows(main_irf, ctrl_irf, lag5_irf, quad_irf) |>
  mutate(spec = factor(spec, levels = c(
    "Main",
    "+ controls (JSL, NAO, wars, deaths)",
    "5 ENSO lags",
    "Quadratic record trend")))

n_specs <- length(levels(all_irf$spec))
offsets <- setNames(seq(-0.27, 0.27, length.out = n_specs),
                    levels(all_irf$spec))
all_irf$x <- all_irf$horizon + offsets[as.character(all_irf$spec)]

spec_cols <- c(
  "Main"                                = "black",
  "+ controls (JSL, NAO, wars, deaths)" = "firebrick",
  "5 ENSO lags"                          = "darkorange2",
  "Quadratic record trend"               = "purple3"
)

pA <- ggplot(all_irf, aes(x = x, y = irf_mean, colour = spec)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_errorbar(aes(ymin = irf_down, ymax = irf_up),
                width = 0.08, linewidth = 0.8) +
  geom_point(size = 2.6) +
  scale_x_continuous(breaks = HORIZONS) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_colour_manual(values = spec_cols, name = NULL) +
  labs(x = "Horizon (years)",
       y = "% Harvest response per +1 SD NINO3.4",
       title = "Robustness across specifications",
       tag   = "a") +
  theme_classic(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.box       = "vertical",
        panel.grid.minor = element_blank(),
        plot.tag         = element_text(face = "bold", size = 16)) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

# ── Panel (b): 1000-shuffle placebo ──────────────────────────────────────────
N_SHUFFLE <- 2000
SEED <- 20260626

year_enso <- WR |> distinct(Year, nino34) |> arrange(Year)
WR_no_nino <- WR |> dplyr::select(-nino34)

fit_b_only <- function(df, h) {
  fml <- as.formula(paste0(
    "f(logyield,", h, ") - l(logyield,1) ~ ",
    "i(teleco_PDSI_10, nino34) + i(decade) + l(nino34, 1:3) | ",
    "VarLocationGrain[Year]"))
  m <- tryCatch(
    feols(fml, data = df,
          panel.id = c("VarLocationGrain", "Year"),
          vcov     = "iid"),  # vcov not used; speed up
    error = function(e) NULL)
  if (is.null(m)) return(NA_real_)
  cn <- names(coef(m))
  nm <- grep("teleco_PDSI_10::1:nino34", cn, value = TRUE)
  if (length(nm) == 0) return(NA_real_)
  as.numeric(unname(coef(m)[nm]))
}

actual <- main_irf$irf_mean
cat("Actual Teleco coefs (h = 0..3):\n"); print(round(actual, 4))

set.seed(SEED)
null_mat <- matrix(NA_real_, nrow = N_SHUFFLE, ncol = length(HORIZONS),
                   dimnames = list(NULL, paste0("h", HORIZONS)))
t0 <- Sys.time()
for (b in seq_len(N_SHUFFLE)) {
  if (b %% 50 == 0) cat("shuffle ", b, "/", N_SHUFFLE,
                         " (",
                         sprintf("%.1f", as.numeric(difftime(Sys.time(), t0,
                                                              units = "mins"))),
                         " min)\n", sep = "")
  shuffled <- year_enso |> mutate(nino34 = sample(nino34))
  WR_b <- WR_no_nino |> left_join(shuffled, by = "Year")
  for (j in seq_along(HORIZONS)) {
    null_mat[b, j] <- fit_b_only(WR_b, HORIZONS[j])
  }
}

null_summary <- data.frame(
  horizon = HORIZONS,
  null_mean = colMeans(null_mat, na.rm = TRUE),
  null_lo95 = apply(null_mat, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
  null_hi95 = apply(null_mat, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
)
write.csv(null_mat, file.path(OUT_DAT, "placebo_shuffle_null_2000.csv"),
          row.names = FALSE)
write.csv(null_summary, file.path(OUT_DAT, "placebo_shuffle_null_summary.csv"),
          row.names = FALSE)

pvals <- sapply(seq_along(HORIZONS), function(j) {
  null_h <- null_mat[, j]
  mean(abs(null_h) >= abs(actual[j]), na.rm = TRUE)
})
cat("\nPlacebo two-sided p-values:\n")
print(data.frame(h = HORIZONS, b = round(actual, 4),
                 p = round(pvals, 4)))

actual_df <- data.frame(horizon = HORIZONS, b = actual,
                        p = pvals)

# Annotate h = 0 and h = 3 with permutation p-values
pval_df <- actual_df |> filter(horizon %in% c(0, 3)) |>
  mutate(label = sprintf("p = %.3f", p))

pB <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_errorbar(data = null_summary,
                aes(x = horizon, ymin = null_lo95, ymax = null_hi95),
                width = 0.12, linewidth = 0.9, colour = "cornflowerblue") +
  geom_point(data = null_summary,
             aes(x = horizon, y = null_mean),
             colour = "cornflowerblue", size = 2.4) +
  geom_point(data = actual_df,
             aes(x = horizon, y = b),
             colour = "firebrick", size = 3.0) +
  geom_text(data = pval_df,
            aes(x = horizon, y = b, label = label),
            hjust = -0.25, vjust = 0.5,
            colour = "firebrick", size = 4.0, fontface = "bold") +
  scale_x_continuous(breaks = HORIZONS,
                     limits = c(min(HORIZONS) - 0.25, max(HORIZONS) + 0.7)) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(x = "Horizon (years)",
       y = "Teleconnected slope per +1 SD NINO3.4",
       title = "Placebo shuffle (2000 reps)",
       tag   = "b") +
  theme_classic(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 16))

# ── Combine ──────────────────────────────────────────────────────────────────
combined <- pA + pB + plot_layout(widths = c(1.15, 1.0))
ggsave(file.path(OUT_FIG, "fig3c_robust_and_placebo.pdf"),
       combined, width = 12, height = 5.0, device = cairo_pdf)
cat("\nSaved fig3c_robust_and_placebo.pdf\n")
