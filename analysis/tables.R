# Replication: all LaTeX regression tables and summary statistics
#   summary_price.tex              – grain price summary statistics
#   summary_yield_TI.tex           – grain tithes summary statistics
#   summary_yield_YR.tex           – grain yield/ratio summary statistics
#   summary_fishprice.tex          – fish price summary statistics
#   Ljungqvist_ENSO_price_main.tex – grain price FE-DL regression
#   FE_ENSO_yieldPDSI.tex          – grain harvest FE-DL regression (PDSI teleconnection)
#   ENSO_fishprice_main.tex        – fish price FE-DL regression
#   granger_enso.tex               – Granger causality checks

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

  ## 1b. FE-DL regressions
  est1 <- feols(logprice ~ Year:Location + l(nino34, 0:3) + i(Decade) | Location,
                data = data, panel.id = c("Location", "Year"), vcov = DK ~ Year)
  est2 <- feols(logprice ~ Year:Location + l(nino34, 0:3) + i(Decade) +
                  PDSI + temp_summer + temp_winter + precip_summer + precip_winter |
                  Location,
                data = data, panel.id = c("Location", "Year"), vcov = DK ~ Year)
  est3 <- feols(logprice ~ Year:Location + l(nino34, 0:3) + i(Decade) +
                  ongoing_wars + log(1 + Deaths) | Location,
                data = data, panel.id = c("Location", "Year"), vcov = DK ~ Year)
  est4 <- feols(logprice ~ Year:Location + l(nino34, 0:3) + i(Decade) +
                  PDSI + temp_summer + temp_winter + precip_summer + precip_winter +
                  ongoing_wars + log(1 + Deaths) | Location,
                data = data, panel.id = c("Location", "Year"), vcov = DK ~ Year)

  dict_price <- c(
    "l(nino34, 0)"    = "NINO3.4 T",
    "nino34"          = "NINO3.4 T",
    "l(nino34, 1)"    = "NINO3.4 T-1",
    "l(nino34, 2)"    = "NINO3.4 T-2",
    "l(nino34, 3)"    = "NINO3.4 T-3",
    "ongoing_wars"    = "Ongoing Wars",
    "log(1+Deaths)"   = "Log(1+Deaths)",
    "PDSI"            = "PDSI",
    "temp_summer"     = "Summer Temp",
    "temp_winter"     = "Winter Temp",
    "precip_summer"   = "Summer Precip",
    "precip_winter"   = "Winter Precip",
    "logprice"        = "Log Grain Price"
  )

  etable(
    est1, est2, est3, est4,
    vcov         = DK ~ Year,
    keep         = "%nino34",
    dict         = dict_price,
    label        = "tab:grainprice_FE",
    title        = "Effect of a 1\\textdegree C Anomaly in the Nino 3.4 Index on Grain Prices",
    drop.section = "fixef",
    fitstat      = c("n", "r2"),
    extralines   = list(
      "Decade FE"   = c("Yes", "Yes", "Yes", "Yes"),
      "Location FE" = c("Yes", "Yes", "Yes", "Yes"),
      "Controls"    = c("None", "Climate", "Conflict", "Climate + Conflict")
    ),
    file    = file.path(OUT_TAB, "Ljungqvist_ENSO_price_main.tex"),
    replace = TRUE
  )
  cat("Saved Ljungqvist_ENSO_price_main.tex\n")
}


# ══════════════════════════════════════════════════════════════════════════════
# 2. GRAIN YIELDS – summary + FE-DL
# ══════════════════════════════════════════════════════════════════════════════
save_yield_tables <- function() {
  yield_2023 <- read.csv(file.path(DATA_PROC, "yield_ljungqvist_2025.csv")) %>%
    mutate(decade = floor(Year / 10) * 10) %>%
    dplyr::select(VarLocationGrain, Year, decade, everything())

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

  ## 2c. FE-DL regression: wheat/rye harvest × PDSI teleconnection
  est1 <- feols(
    logyield ~ Year:Location + l(nino34, 0:3):i(teleco_PDSI_10) + i(decade) |
      VarLocationGrain,
    data      = yield_2023 %>% filter(Grain %in% c("Wheat", "Rye")),
    panel.id  = ~ VarLocationGrain + Year
  )
  est2 <- feols(
    logyield ~ Year:Location + l(nino34, 0:3):i(teleco_PDSI_10) + i(decade) +
      PDSI + temp_summer + temp_winter + precip_summer + precip_winter |
      VarLocationGrain,
    data      = yield_2023 %>% filter(Grain %in% c("Wheat", "Rye")),
    panel.id  = c("VarLocationGrain", "Year")
  )
  est3 <- feols(
    logyield ~ Year:Location + l(nino34, 0:3):i(teleco_PDSI_10) + i(decade) +
      ongoing_wars + log(1 + Deaths) | VarLocationGrain,
    data      = yield_2023 %>% filter(Grain %in% c("Wheat", "Rye")),
    panel.id  = c("VarLocationGrain", "Year")
  )
  est4 <- feols(
    logyield ~ Year:LocationGrain + l(nino34, 0:3):i(teleco_PDSI_10) + i(decade) +
      ongoing_wars + log(1 + Deaths) +
      temp_summer_europe + temp_winter_europe +
      precip_summer_europe + precip_winter_europe | VarLocationGrain,
    data      = yield_2023 %>% filter(Grain %in% c("Wheat", "Rye")),
    panel.id  = c("VarLocationGrain", "Year")
  )

  dict_yield <- c(
    "l(nino34, 0)"        = "NINO3.4 T",
    "nino34"              = "NINO3.4 T",
    "l(nino34, 1)"        = "NINO3.4 T-1",
    "l(nino34, 2)"        = "NINO3.4 T-2",
    "l(nino34, 3)"        = "NINO3.4 T-3",
    "ongoing_wars"        = "Ongoing Wars",
    "log(1+Deaths)"       = "Log(1+Deaths)",
    "PDSI"                = "PDSI",
    "temp_summer"         = "Summer Temp",
    "temp_winter"         = "Winter Temp",
    "precip_summer"       = "Summer Precip",
    "precip_winter"       = "Winter Precip",
    "teleco_PDSI_10"      = "scPDSI Teleconnection",
    "logyield"            = "Log Grain Harvest"
  )

  etable(
    est1, est2, est3, est4,
    keep         = "%nino34",
    vcov         = DK ~ Year,
    dict         = dict_yield,
    label        = "tab:FE_ENSO_yieldPDSI",
    title        = "Effect of a 1\\textdegree C Anomaly in the Nino 3.4 Index on Wheat and Rye Grain Harvests.",
    drop.section = "fixef",
    fitstat      = c("n", "r2"),
    extralines   = list(
      "Decade FE"               = c("Yes", "Yes", "Yes", "Yes"),
      "Record-Location-Grain FE" = c("Yes", "Yes", "Yes", "Yes"),
      "Controls"                = c("None", "Climate", "Conflict", "Climate + Conflict")
    ),
    file    = file.path(OUT_TAB, "FE_ENSO_yieldPDSI.tex"),
    replace = TRUE
  )
  cat("Saved FE_ENSO_yieldPDSI.tex\n")
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

  ## 3b. FE-DL regressions
  est1 <- feols(logprice ~ Year:Location + l(nino34, 0:3) + i(Decade) | Location,
                data = fishprice, panel.id = c("LocationSpecies", "Year"),
                vcov = DK ~ Year)
  est2 <- feols(logprice ~ Year:Location + l(nino34, 0:3) + i(Decade) +
                  PDSI + temp_summer + temp_winter + precip_summer + precip_winter |
                  Location,
                data = fishprice, panel.id = c("LocationSpecies", "Year"),
                vcov = DK ~ Year)
  est3 <- feols(logprice ~ Year:Location + l(nino34, 0:3) + i(Decade) +
                  ongoing_wars + log(1 + Deaths) | Location,
                data = fishprice, panel.id = c("LocationSpecies", "Year"),
                vcov = DK ~ Year)
  est4 <- feols(logprice ~ Year:Location + l(nino34, 0:3) + i(Decade) +
                  PDSI + temp_summer + temp_winter + precip_summer + precip_winter +
                  ongoing_wars + log(1 + Deaths) | Location,
                data = fishprice, panel.id = c("LocationSpecies", "Year"),
                vcov = DK ~ Year)

  dict_fish <- c(
    "l(nino34, 0)"  = "NINO3.4 T",
    "nino34"        = "NINO3.4 T",
    "l(nino34, 1)"  = "NINO3.4 T-1",
    "l(nino34, 2)"  = "NINO3.4 T-2",
    "l(nino34, 3)"  = "NINO3.4 T-3",
    "ongoing_wars"  = "Ongoing Wars",
    "log(1+Deaths)" = "Log(1+Deaths)",
    "PDSI"          = "PDSI",
    "temp_summer"   = "Summer Temp",
    "temp_winter"   = "Winter Temp",
    "precip_summer" = "Summer Precip",
    "precip_winter" = "Winter Precip",
    "logprice"      = "Log Fish Price"
  )

  etable(
    est1, est2, est3, est4,
    vcov         = DK ~ Year,
    keep         = "%nino34",
    dict         = dict_fish,
    label        = "tab:fishprice_FE",
    title        = "Effect of a 1\\textdegree C Anomaly in the Nino 3.4 Index on Fish Prices",
    drop.section = "fixef",
    fitstat      = c("n", "r2"),
    extralines   = list(
      "Decade FE"              = c("Yes", "Yes", "Yes", "Yes"),
      "Location-Species FE"    = c("Yes", "Yes", "Yes", "Yes"),
      "Controls"               = c("None", "Climate", "Conflict", "Climate + Conflict")
    ),
    file    = file.path(OUT_TAB, "ENSO_fishprice_main.tex"),
    replace = TRUE
  )
  cat("Saved ENSO_fishprice_main.tex\n")
}


# ══════════════════════════════════════════════════════════════════════════════
# 4. GRANGER CAUSALITY CHECKS
# ══════════════════════════════════════════════════════════════════════════════
save_granger_table <- function() {
  data      <- read.csv(file.path(DATA_PROC, "price_2023_enso.csv")) %>%
    mutate(Decade = floor(Year / 10) * 10) %>% drop_na() %>% distinct()
  fishprice <- read.csv(file.path(DATA_PROC, "fishprice_enso.csv")) %>% drop_na()
  yield_2023 <- read.csv(file.path(DATA_PROC, "yield_ljungqvist_2025.csv")) %>%
    filter(Grain %in% c("Wheat", "Rye"))

  # Annual averages for panel Granger tests
  ljunqgvist_avg <- data %>%
    group_by(Year) %>%
    summarise(avg_price  = mean(logprice, na.rm = TRUE),
              avg_nino34 = mean(nino34,   na.rm = TRUE),
              .groups = "drop")

  ljunqgvist_avg_yield <- yield_2023 %>%
    filter(teleco_PDSI_10 == 1) %>%
    group_by(Year) %>%
    summarise(avg_yield  = mean(logyield, na.rm = TRUE),
              avg_nino34 = mean(nino34,   na.rm = TRUE),
              .groups = "drop")

  avg_fishprice <- fishprice %>%
    group_by(Location, Year) %>%
    summarise(avg_price  = mean(logprice, na.rm = TRUE),
              avg_nino34 = mean(nino34,   na.rm = TRUE),
              .groups = "drop")

  # Helper: safe p-value extraction
  granger_p <- function(formula, data, order) {
    tryCatch(
      grangertest(formula, data = data, order = order)$`Pr(>F)`[2],
      error = function(e) NA_real_
    )
  }

  results <- tibble::tibble(
    Variable  = c(
      "Grain Price (all cities)", "Grain Price (all cities)",
      "Harvests (Teleconnected)", "Harvests (Teleconnected)",
      "Fish Price (all cities)",  "Fish Price (all cities)"
    ),
    Direction = c(
      "Grain Price $\\to$ NINO3.4",  "NINO3.4 $\\to$ Grain Price",
      "Harvest $\\to$ NINO3.4",      "NINO3.4 $\\to$ Harvest",
      "Fish Price $\\to$ NINO3.4",   "NINO3.4 $\\to$ Fish Price"
    ),
    p_value = c(
      granger_p(avg_nino34 ~ avg_price, ljunqgvist_avg,       order = 10),
      granger_p(avg_price  ~ avg_nino34, ljunqgvist_avg,      order = 1),
      granger_p(avg_nino34 ~ avg_yield, ljunqgvist_avg_yield, order = 10),
      granger_p(avg_yield  ~ avg_nino34, ljunqgvist_avg_yield, order = 1),
      granger_p(avg_nino34 ~ avg_price, avg_fishprice,        order = 10),
      granger_p(avg_price  ~ avg_nino34, avg_fishprice,       order = 1)
    )
  ) %>%
    mutate(p_value = round(p_value, 3))

  latex_table <- knitr::kable(
    results,
    format    = "latex",
    booktabs  = TRUE,
    escape    = FALSE,
    col.names = c("Variable", "Direction", "\\textit{p}-value"),
    caption   = "Granger causality tests for the relationship between ENSO (NINO3.4), prices, and yields.",
    label     = "granger_enso"
  )

  writeLines(latex_table, file.path(OUT_TAB, "granger_enso.tex"))
  cat("Saved granger_enso.tex\n")
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

  cat("\n=== Granger Causality Table ===\n")
  save_granger_table()

  cat("\nAll tables saved to", OUT_TAB, "\n")
}
