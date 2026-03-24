###############################################################################
# Figure 2 – ENSO and famine onset
# This script:
#   1. Computes and exports tidy plot-data CSVs for assemble_figures.py
#   2. Saves LaTeX regression tables
#
# CSV exports → analysis/output/data/
#   fig2A_onset_box.csv        – NINO3.4 by group (famine onset / no famine)
#   fig2B_onset_coef.csv       – LPM coefficients by region × model
#
# Tables → analysis/output/tables/
#   famine_starts_reg.tex
#   figA_loo_lpm_onset.csv     (appendix LOO data)
#   figA_permutation_onset.csv (appendix permutation null distribution)
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(fixest)
  library(broom)
  library(purrr)
  library(modelsummary)
})

# ── paths ─────────────────────────────────────────────────────────────────────
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

SCRIPT_DIR <- .script_dir()
ROOT       <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
DATA_PROC  <- file.path(ROOT, "processed data")
out_data   <- file.path(SCRIPT_DIR, "output", "data")
fig_out    <- file.path(ROOT, "output", "figures", "appendix")
tab_out    <- file.path(SCRIPT_DIR, "output", "tables")
for (d in c(out_data, fig_out, tab_out)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

data <- read.csv(file.path(DATA_PROC, "famine_region_data.csv")) |>
  mutate(Decade = floor(Year / 10) * 10) |>
  arrange(Region, Year)


###############################################################################
# Fig 2A – NINO3.4 boxplot data
###############################################################################
export_fig2A <- function(data) {
  central <- data |>
    filter(Region == "Central Europe") |>
    arrange(Year) |>
    mutate(Famine_onset = (Famine_start == 1))

  bind_rows(
    central |> filter(Famine_onset) |>
      transmute(nino34, Group = "Famine Onset"),
    central |> filter(Famine == 0) |>
      transmute(nino34, Group = "No Famine")
  ) |>
    write.csv(file.path(out_data, "fig2A_onset_box.csv"), row.names = FALSE)
  message("Exported fig2A_onset_box.csv")
}


###############################################################################
# Fig 2B – LPM coefficients
###############################################################################
export_fig2B <- function(data) {
  est0 <- feols(Famine_start ~ nino34:Region,
                data = data, vcov = NW ~ Year,
                panel.id = c("Region", "Year"))

  est2 <- feols(Famine_start ~ nino34:Region + i(Decade) | Region,
                data = data, vcov = DK ~ Year,
                panel.id = c("Region", "Year"))

  est3 <- feols(
    Famine_start ~ nino34:Region + i(Decade) +
      PDSI:Region + temp_summer:Region + temp_winter:Region +
      precip_summer:Region + precip_winter:Region +
      ongoing_wars:Region + log(1 + Deaths):Region | Region,
    data = data, vcov = DK ~ Year,
    panel.id = c("Region", "Year")
  )

  bind_rows(
    broom::tidy(est0, conf.int = TRUE) |> mutate(model = "No FEs"),
    broom::tidy(est2, conf.int = TRUE) |> mutate(model = "FEs"),
    broom::tidy(est3, conf.int = TRUE) |> mutate(model = "FEs + Controls")
  ) |>
    filter(grepl("nino34", term)) |>
    mutate(Region = trimws(gsub("nino34:Region", "", term))) |>
    write.csv(file.path(out_data, "fig2B_onset_coef.csv"), row.names = FALSE)
  message("Exported fig2B_onset_coef.csv")
}


###############################################################################
# Regression table
###############################################################################
save_onset_table <- function(data) {
  est0 <- feols(Famine_start ~ nino34:Region,
                data = data, vcov = NW ~ Year, panel.id = c("Region", "Year"))
  est1 <- feols(Famine_start ~ nino34:Region + i(Decade),
                data = data, vcov = NW ~ Year, panel.id = c("Region", "Year"))
  est2 <- feols(Famine_start ~ nino34:Region + i(Decade) | Region,
                data = data, vcov = DK ~ Year, panel.id = c("Region", "Year"))
  est3 <- feols(
    Famine_start ~ nino34:Region + i(Decade) + PDSI:Region +
      temp_summer:Region + temp_winter:Region + precip_summer:Region +
      precip_winter:Region + ongoing_wars:Region + log(1 + Deaths):Region | Region,
    data = data, vcov = DK ~ Year, panel.id = c("Region", "Year")
  )
  est4 <- feglm(
    Famine_start ~ nino34:Region + i(Decade):Region + PDSI:Region +
      temp_summer:Region + temp_winter:Region + precip_summer:Region +
      precip_winter:Region + ongoing_wars:Region + log(1 + Deaths):Region | Region,
    vcov = "iid", data = data, panel.id = c("Region", "Year"), family = "logit"
  )
  est5 <- feols(
    Famine_start ~ nino34 + i(Decade) + PDSI + temp_summer + temp_winter +
      precip_summer + precip_winter + ongoing_wars + log(1 + Deaths) | Region,
    data = data, vcov = DK ~ Year, panel.id = c("Region", "Year")
  )
  est6 <- feglm(
    Famine_start ~ nino34 + i(Decade) + PDSI + temp_summer + temp_winter +
      precip_summer + precip_winter + ongoing_wars + log(1 + Deaths) | Region,
    vcov = "iid", data = data, panel.id = c("Region", "Year"), family = "logit"
  )

  etable(
    est0, est1, est2, est3, est4, est5, est6,
    keep         = c("Nino3.4", "(Intercept)"),
    drop.section = "fixef",
    fitstat      = c("AIC", "BIC", "N", "cor2"),
    dict         = c(nino34 = "Nino3.4", Famine_start = "Famine Start"),
    extralines   = list(
      "Decade FE"       = c("No",  "Yes", "Yes", "Yes", "Yes", "Yes", "Yes"),
      "Region FE"       = c("No",  "No",  "Yes", "Yes", "Yes", "Yes", "Yes"),
      "Controls"        = c("No",  "No",  "No",  "Yes", "Yes", "Yes", "Yes"),
      "Standard Errors" = c("NW",  "NW",  "DK",  "DK",  "IID", "DK",  "IID")
    ),
    file    = file.path(tab_out, "famine_starts_reg.tex"),
    label   = "tab:famine_starts",
    title   = "Effect of a 1\\textdegree C Anomaly in the Nino 3.4 Index on the Probability of a Famine Start",
    replace = TRUE
  )
  message("Saved famine_starts_reg.tex")
}


###############################################################################
# Appendix: LOO sensitivity (exported as CSV for Python if desired)
###############################################################################
export_loo_onset <- function(data) {
  est_full <- feols(Famine_start ~ nino34:Region + i(Decade),
                    data = data, vcov = NW ~ Year,
                    panel.id = c("Region", "Year"))
  full_coef        <- broom::tidy(est_full)
  baseline_central <- full_coef |>
    filter(term == "nino34:RegionCentral Europe") |>
    pull(estimate)

  famine_ids <- which(data$Famine_start == 1 & data$Region == "Central Europe")

  loo_results <- map_dfr(seq_along(famine_ids), function(i) {
    data_loo <- data[-famine_ids[i], ]
    mod <- feols(
      Famine_start ~ nino34:Region + i(Decade):Region,
      data = data_loo, vcov = NW ~ Year,
      panel.id = c("Region", "Year"), notes = FALSE
    )
    broom::tidy(mod) |> mutate(left_out = paste0("Famine_", i))
  }) |>
    filter(term == "nino34:RegionCentral Europe") |>
    mutate(baseline = baseline_central)

  write.csv(loo_results, file.path(out_data, "figA_loo_lpm_onset.csv"),
            row.names = FALSE)
  message("Exported figA_loo_lpm_onset.csv")
}


###############################################################################
# Appendix: permutation test (exported as CSV)
###############################################################################
export_permutation_onset <- function(data, n_iter = 10000) {
  set.seed(42)
  formula    <- Famine_start ~ nino34:Region + i(Decade)
  true_model <- feols(formula, data = data, warn = FALSE, notes = FALSE)
  true_est   <- broom::tidy(true_model) |>
    filter(term == "nino34:RegionCentral Europe") |>
    pull(estimate)

  coef_vec <- vapply(seq_len(n_iter), function(i) {
    df_sh <- data |> mutate(nino34 = sample(nino34))
    m     <- try(feols(formula, data = df_sh, warn = FALSE, notes = FALSE),
                 silent = TRUE)
    if (inherits(m, "try-error")) return(NA_real_)
    broom::tidy(m) |>
      filter(term == "nino34:RegionCentral Europe") |>
      pull(estimate)
  }, numeric(1))

  data.frame(estimate = coef_vec, true_estimate = rep(true_est, length(coef_vec))) |>
    write.csv(file.path(out_data, "figA_permutation_onset.csv"), row.names = FALSE)
  message("Exported figA_permutation_onset.csv")
}


###############################################################################
# Appendix plots
###############################################################################
fig_out <- file.path(SCRIPT_DIR, "output", "figures", "appendix")
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

plot_loo_onset <- function() {
  path <- file.path(out_data, "figA_loo_lpm_onset.csv")
  if (!file.exists(path)) { message("figA_loo_lpm_onset.csv not found"); return() }
  df <- read.csv(path) |>
    arrange(estimate) |>
    mutate(i     = seq_len(n()),
           lower = estimate - 1.96 * std.error,
           upper = estimate + 1.96 * std.error)
  p <- ggplot(df, aes(x = i, y = estimate, ymin = lower, ymax = upper)) +
    geom_hline(aes(yintercept = baseline), linetype = "dashed",
               color = "red", linewidth = 0.8) +
    geom_pointrange(color = "grey40", size = 0.4) +
    labs(x = "Famine onset left out",
         y = "LPM coefficient (NINO3.4)") +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank())
  ggsave(file.path(fig_out, "figA_loo_onset.pdf"), p,
         width = 8, height = 5, device = cairo_pdf)
  message("Saved figA_loo_onset.pdf")
}

plot_permutation_onset <- function() {
  path <- file.path(out_data, "figA_permutation_onset.csv")
  if (!file.exists(path)) { message("figA_permutation_onset.csv not found"); return() }
  df    <- read.csv(path) |> dplyr::filter(!is.na(estimate))
  true_v <- df$true_estimate[1]
  p_val  <- round(mean(df$estimate >= true_v), 3)
  p <- ggplot(df, aes(x = estimate)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50,
                   fill = "cornflowerblue", color = "white", alpha = 0.7) +
    geom_vline(xintercept = true_v, color = "orange", linetype = "solid", linewidth = 1) +
    annotate("text", x = true_v, y = Inf, hjust = -0.15, vjust = 1.5,
             label = paste0("True est.\np = ", p_val),
             color = "black", size = 4) +
    labs(x = "Permuted LPM coefficient",
         y = "Density") +
    theme_classic(base_size = 13) +
    theme(panel.grid = element_blank())
  ggsave(file.path(fig_out, "figA_permutation_onset.pdf"), p,
         width = 7, height = 5, device = cairo_pdf)
  message("Saved figA_permutation_onset.pdf")
}


###############################################################################
# Run
###############################################################################
if (!interactive()) {
  export_fig2A(data)
  export_fig2B(data)
  save_onset_table(data)
  export_loo_onset(data)
  export_permutation_onset(data, n_iter = 100)
  plot_loo_onset()
  plot_permutation_onset()
  message("Fig 2 data exports complete.")
} else {
  message("Call functions interactively or source this script.")
}
