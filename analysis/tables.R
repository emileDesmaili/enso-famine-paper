# Replication: all LaTeX regression tables and summary statistics
#   summary_price.tex              – grain price summary statistics
#   summary_yield_TI.tex           – grain tithes summary statistics
#   summary_yield_YR.tex           – grain yield/ratio summary statistics
#   summary_fishprice.tex          – fish price summary statistics
#   Ljungqvist_ENSO_price_main.tex – grain price FE-DL regression
#   FE_ENSO_yieldPDSI.tex          – grain harvest FE-DL regression (PDSI teleconnection)
#   ENSO_fishprice_main.tex        – fish price FE-DL regression

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(fixest)
  library(modelsummary)
  library(lmtest)
  library(haven)
})

# ── paths ──────────────────────────────────────────────────────────────────────
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
OUT_TAB    <- file.path(SCRIPT_DIR, "output", "tables")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# HELPER — cumulative DL table writer
# Fits are assumed to have LHS d(y) and RHS l(nino34, 0:3) + controls | FE.
# Under a level DL  y_t = Σ β_k x_{t-k},  the differenced coefficients satisfy
# γ_k = β_k − β_{k-1}, so cumulative sums Σ_{k=0}^{h} γ_k recover the level
# impulse response β_h.  Driscoll-Kraay SEs applied via vcov(m, vcov = DK ~ Year).
# ══════════════════════════════════════════════════════════════════════════════
save_cum_dl_table <- function(models, output_path,
                              title, label, dep_label,
                              fe_rows,
                              varying_slopes_rows = NULL,
                              lag_names = c("l(nino34, 0)", "l(nino34, 1)",
                                            "l(nino34, 2)", "l(nino34, 3)"),
                              vcov_formula = DK ~ Year) {

  cum_effect <- function(m, k) {
    V <- vcov(m, vcov = vcov_formula)
    b <- coef(m)[lag_names]
    V <- V[lag_names, lag_names]
    L <- c(rep(1, k + 1), rep(0, length(lag_names) - k - 1))
    est_val <- sum(L * b)
    se_val  <- sqrt(as.numeric(t(L) %*% V %*% L))
    z       <- est_val / se_val
    p       <- 2 * (1 - pnorm(abs(z)))
    list(est = est_val, se = se_val, p = p)
  }

  stars <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.01)      "$^{***}$"
    else if (p < 0.05) "$^{**}$"
    else if (p < 0.1)  "$^{*}$"
    else               ""
  }
  fmt_est <- function(r) paste0(sprintf("%.4f", r$est), stars(r$p))
  fmt_se  <- function(r) sprintf("(%.4f)", r$se)

  row_labels <- c(
    "NINO3.4 T (sum lag 0)",
    "NINO3.4 T to T-1 (sum lags 0-1)",
    "NINO3.4 T to T-2 (sum lags 0-2)",
    "NINO3.4 T to T-3 (sum lags 0-3)"
  )

  make_body_row <- function(k) {
    rs   <- lapply(models, cum_effect, k = k)
    ests <- vapply(rs, fmt_est, character(1))
    ses  <- vapply(rs, fmt_se,  character(1))
    c(
      paste0("      ", format(row_labels[k + 1], width = 35),
             " & ", paste(ests, collapse = " & "), "\\\\"),
      paste0("      ", strrep(" ", 35),
             " & ", paste(ses,  collapse = " & "), "\\\\")
    )
  }
  body_rows <- unlist(lapply(0:(length(lag_names) - 1), make_body_row))

  n_obs <- vapply(models, function(m) as.integer(nobs(m)), integer(1))
  r2    <- vapply(models, function(m) as.numeric(fitstat(m, "r2")$r2), numeric(1))

  fmt_row <- function(name, vals)
    paste0("      ", format(name, width = 35), " & ",
           paste(sprintf("%-13s", vals), collapse = " & "), "\\\\")

  fe_tex <- vapply(names(fe_rows),
                   function(nm) fmt_row(nm, fe_rows[[nm]]),
                   character(1))

  vs_tex <- if (!is.null(varying_slopes_rows)) {
    c("      \\midrule",
      "      \\emph{Varying Slopes}\\\\",
      vapply(names(varying_slopes_rows),
             function(nm) fmt_row(nm, varying_slopes_rows[[nm]]),
             character(1)))
  } else character(0)

  tex <- c(
    "",
    "\\begin{table}[htbp]",
    paste0("   \\caption{\\label{", label, "} ", title, "}"),
    "   \\centering",
    "   \\begin{tabular}{lcccc}",
    "      \\tabularnewline \\midrule \\midrule",
    paste0("      Dependent Variable: & \\multicolumn{4}{c}{", dep_label, "}\\\\"),
    "      Model:              & (1)           & (2)           & (3)           & (4)\\\\",
    "      \\midrule",
    "      \\emph{Cumulative NINO3.4 effects}\\\\",
    body_rows,
    "      \\midrule",
    fe_tex,
    vs_tex,
    "      \\midrule",
    "      \\emph{Fit statistics}\\\\",
    fmt_row("Observations", format(n_obs, big.mark = ",")),
    fmt_row("R$^2$",        sprintf("%.5f", r2)),
    "      \\midrule \\midrule",
    "      \\multicolumn{5}{l}{\\emph{Driscoll-Kraay (L=4) standard-errors in parentheses}}\\\\",
    "      \\multicolumn{5}{l}{\\emph{Signif. Codes: ***: 0.01, **: 0.05, *: 0.1}}\\\\",
    "   \\end{tabular}",
    "\\end{table}",
    "",
    ""
  )

  writeLines(tex, output_path)
  cat("Saved", basename(output_path), "\n")
}


# ══════════════════════════════════════════════════════════════════════════════
# 1. GRAIN PRICES – summary + FE-DL
# ══════════════════════════════════════════════════════════════════════════════
save_price_tables <- function() {
  data <- read.csv(file.path(DATA_PROC, "price_2023_enso.csv")) %>%
    mutate(Decade = floor(Year / 10) * 10) %>%
    dplyr::select(Location, Year, Decade, everything()) %>%
    drop_na() %>%
    distinct() %>%
    arrange(Year)

  ## 1a. Summary statistics
  data %>%
    rename(
      `temp AMJJ`   = temp_summer,
      `temp NDJF`   = temp_winter,
      `precip AMJJ` = precip_summer,
      `precip NDJF` = precip_winter,
      NINO3.4       = nino34,
      `ongoing wars` = ongoing_wars
    ) %>%
    datasummary(
      formula = Price + NINO3.4 + PDSI + JSL + NAO_cal +
        `temp AMJJ` + `temp NDJF` + `precip AMJJ` + `precip NDJF` +
        `ongoing wars` + Deaths ~ Mean + SD + Min + Max + N,
      output  = file.path(OUT_TAB, "summary_price.tex"),
      title   = "Summary Statistics for Grain Price Data"
    )
  cat("Saved summary_price.tex\n")

  ## 1b. FE-DL regressions on Δ log price (cumulative sums recover level β_h).
  ## Location[Year] adds a Location-specific linear trend on the differenced
  ## series (= Location-specific parabolic trend on level logprice).
  est1 <- feols(d(logprice) ~ l(nino34, 0:3) + i(Decade) | Location[Year],
                data = data, panel.id = c("Location", "Year"), vcov = DK ~ Year)
  est2 <- feols(d(logprice) ~ l(nino34, 0:3) + i(Decade) +
                  PDSI + temp_summer + temp_winter + precip_summer + precip_winter |
                  Location[Year],
                data = data, panel.id = c("Location", "Year"), vcov = DK ~ Year)
  est3 <- feols(d(logprice) ~ l(nino34, 0:3) + i(Decade) +
                  ongoing_wars + log(1 + Deaths) | Location[Year],
                data = data, panel.id = c("Location", "Year"), vcov = DK ~ Year)
  est4 <- feols(d(logprice) ~ l(nino34, 0:3) + i(Decade) +
                  PDSI + temp_summer + temp_winter + precip_summer + precip_winter +
                  ongoing_wars + log(1 + Deaths) | Location[Year],
                data = data, panel.id = c("Location", "Year"), vcov = DK ~ Year)

  save_cum_dl_table(
    list(est1, est2, est3, est4),
    output_path = file.path(OUT_TAB, "Ljungqvist_ENSO_price_main.tex"),
    title       = "Effect of a 1\\textdegree C Anomaly in the Nino 3.4 Index on Grain Prices --- Cumulative $\\Delta$ Log Response.",
    label       = "tab:grainprice_FE",
    dep_label   = "$\\Delta$ Log Grain Price",
    fe_rows = list(
      "Decade FE"                = rep("Yes", 4),
      "Location FE"              = rep("Yes", 4),
      "Location-specific trend"  = rep("Yes", 4),
      "Controls"                 = c("None", "Climate", "Conflict", "Climate + Conflict")
    )
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# 2. GRAIN YIELDS – summary + FE-DL
# ══════════════════════════════════════════════════════════════════════════════
save_yield_tables <- function() {
  yield_2023 <- read.csv(file.path(DATA_PROC, "yield_ljungqvist_v2.csv")) %>%
    mutate(decade = floor(Year / 10) * 10,
           Year2  = Year^2) %>%
    dplyr::select(VarLocationGrain, Year, decade, Year2, everything())

  TI_2023 <- yield_2023 %>%
    filter(Type == "TI") %>%
    drop_na() %>%
    distinct() %>%
    dplyr::select(LocationGrain, Year, decade, everything())

  YR_2023 <- yield_2023 %>%
    filter(Type %in% c("YR", "YI")) %>%
    drop_na() %>%
    distinct() %>%
    dplyr::select(VarLocationGrain, Year, decade, everything())

  ## 2a. Summary statistics – tithes
  TI_2023 %>%
    rename(
      `temp AMJJ`    = temp_summer,
      `temp NDJF`    = temp_winter,
      `precip AMJJ`  = precip_summer,
      `precip NDJF`  = precip_winter,
      NINO3.4        = nino34,
      `ongoing wars` = ongoing_wars
    ) %>%
    datasummary(
      formula = logyield + NINO3.4 + PDSI + JSL + NAO_cal +
        `temp AMJJ` + `temp NDJF` + `precip AMJJ` + `precip NDJF` +
        `ongoing wars` + Deaths ~ Mean + SD + Min + Max + N,
      output = file.path(OUT_TAB, "summary_yield_TI.tex"),
      title  = "Summary Statistics for Grain Tithes Data"
    )
  cat("Saved summary_yield_TI.tex\n")

  ## 2b. Summary statistics – yield ratios
  YR_2023 %>%
    rename(
      `temp AMJJ`    = temp_summer,
      `temp NDJF`    = temp_winter,
      `precip AMJJ`  = precip_summer,
      `precip NDJF`  = precip_winter,
      NINO3.4        = nino34,
      `ongoing wars` = ongoing_wars
    ) %>%
    datasummary(
      formula = logyield + NINO3.4 + PDSI + JSL + NAO_cal +
        `temp AMJJ` + `temp NDJF` + `precip AMJJ` + `precip NDJF` +
        `ongoing wars` + Deaths ~ Mean + SD + Min + Max + N,
      output = file.path(OUT_TAB, "summary_yield_YR.tex"),
      title  = "Summary Statistics for Grain Yields Data"
    )
  cat("Saved summary_yield_YR.tex\n")

  ## 2c. FE-DL on Δ log yield (cumulative sums recover level β_h).
  ## Per-unit linear trend in levels becomes a constant after differencing,
  ## absorbed by the VarLocationGrain FE — [Year] varying slope dropped.
  WR <- yield_2023 %>% filter(Grain %in% c("Wheat", "Rye"))
  est1 <- feols(
    d(logyield) ~ l(nino34, 0:3) + i(decade) | VarLocationGrain,
    data = WR, panel.id = c("VarLocationGrain", "Year"), vcov = DK ~ Year
  )
  est2 <- feols(
    d(logyield) ~ l(nino34, 0:3) + i(decade) +
      PDSI + temp_summer + temp_winter + precip_summer + precip_winter |
      VarLocationGrain,
    data = WR, panel.id = c("VarLocationGrain", "Year"), vcov = DK ~ Year
  )
  est3 <- feols(
    d(logyield) ~ l(nino34, 0:3) + i(decade) +
      ongoing_wars + log(1 + Deaths) | VarLocationGrain,
    data = WR, panel.id = c("VarLocationGrain", "Year"), vcov = DK ~ Year
  )
  est4 <- feols(
    d(logyield) ~ l(nino34, 0:3) + i(decade) +
      ongoing_wars + log(1 + Deaths) +
      temp_summer + temp_winter + precip_summer + precip_winter |
      VarLocationGrain,
    data = WR, panel.id = c("VarLocationGrain", "Year"), vcov = DK ~ Year
  )

  save_cum_dl_table(
    list(est1, est2, est3, est4),
    output_path = file.path(OUT_TAB, "FE_ENSO_yieldPDSI.tex"),
    title       = "Effect of a 1\\textdegree C Anomaly in the Nino 3.4 Index on Wheat and Rye Grain Harvests --- Cumulative $\\Delta$ Log Response.",
    label       = "tab:FE_ENSO_yieldPDSI",
    dep_label   = "$\\Delta$ Log Grain Harvest",
    fe_rows = list(
      "Decade FE"                = rep("Yes", 4),
      "Record-Location-Grain FE" = rep("Yes", 4),
      "Controls"                 = c("None", "Climate", "Conflict", "Climate + Conflict")
    )
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# 3. FISH PRICES – summary + FE-DL
# ══════════════════════════════════════════════════════════════════════════════
save_fishprice_tables <- function() {
  fishprice <- read.csv(file.path(DATA_PROC, "fishprice_enso.csv")) %>%
    mutate(Decade = floor(Year / 10) * 10) %>%
    dplyr::select(LocationSpecies, Year, Decade, everything())

  ## 3a. Summary statistics
  fishprice %>%
    rename(
      `temp AMJJ`    = temp_summer,
      `temp NDJF`    = temp_winter,
      `precip AMJJ`  = precip_summer,
      `precip NDJF`  = precip_winter,
      NINO3.4        = nino34,
      `ongoing wars` = ongoing_wars
    ) %>%
    datasummary(
      formula = Price + NINO3.4 + PDSI + JSL + NAO_cal +
        `temp AMJJ` + `temp NDJF` + `precip AMJJ` + `precip NDJF` +
        `ongoing wars` + Deaths ~ Mean + SD + Min + Max + N,
      output = file.path(OUT_TAB, "summary_fishprice.tex"),
      title  = "Summary Statistics for Fish Prices Data"
    )
  cat("Saved summary_fishprice.tex\n")

  ## 3b. FE-DL on Δ log price (cumulative sums recover level β_h).
  ## Location[Year] adds a Location-specific linear trend on the differenced
  ## series (= Location-specific parabolic trend on level logprice).
  est1 <- feols(d(logprice) ~ l(nino34, 0:3) + i(Decade) | Location[Year],
                data = fishprice, panel.id = c("LocationSpecies", "Year"),
                vcov = DK ~ Year)
  est2 <- feols(d(logprice) ~ l(nino34, 0:3) + i(Decade) +
                  PDSI + temp_summer + temp_winter + precip_summer + precip_winter |
                  Location[Year],
                data = fishprice, panel.id = c("LocationSpecies", "Year"),
                vcov = DK ~ Year)
  est3 <- feols(d(logprice) ~ l(nino34, 0:3) + i(Decade) +
                  ongoing_wars + log(1 + Deaths) | Location[Year],
                data = fishprice, panel.id = c("LocationSpecies", "Year"),
                vcov = DK ~ Year)
  est4 <- feols(d(logprice) ~ l(nino34, 0:3) + i(Decade) +
                  PDSI + temp_summer + temp_winter + precip_summer + precip_winter +
                  ongoing_wars + log(1 + Deaths) | Location[Year],
                data = fishprice, panel.id = c("LocationSpecies", "Year"),
                vcov = DK ~ Year)

  save_cum_dl_table(
    list(est1, est2, est3, est4),
    output_path = file.path(OUT_TAB, "ENSO_fishprice_main.tex"),
    title       = "Effect of a 1\\textdegree C Anomaly in the Nino 3.4 Index on Fish Prices --- Cumulative $\\Delta$ Log Response.",
    label       = "tab:fishprice_FE",
    dep_label   = "$\\Delta$ Log Fish Price",
    fe_rows = list(
      "Decade FE"                = rep("Yes", 4),
      "Location-Species FE"      = rep("Yes", 4),
      "Location-specific trend"  = rep("Yes", 4),
      "Controls"                 = c("None", "Climate", "Conflict", "Climate + Conflict")
    )
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════
if (!interactive()) {
  cat("=== Grain Price Tables ===\n")
  save_price_tables()

  cat("\n=== Grain Yield Tables ===\n")
  save_yield_tables()

  cat("\n=== Fish Price Tables ===\n")
  save_fishprice_tables()

  cat("\nAll tables saved to", OUT_TAB, "\n")
}
