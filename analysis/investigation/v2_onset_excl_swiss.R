###############################################################################
# 30_v2_onset_excl_swiss.R
#
# Re-estimate the main famine-onset model (matching fig 2B / Extended Data
# Table 1) after excluding "Swiss" famines from the Central Europe
# chronology.
#
# A Central-Europe famine-start year is flagged as Swiss-driven if the
# average within-record z-score of Swiss WR yields (Basel rye, Zurich rye,
# Winterthur rye) in that year is below −0.5 SD.
#
# We then set Famine_start = 0 in those CE rows and re-run:
#   est0 = LPM, no FEs
#   est2 = LPM + Decade + Region FE
#   est3 = LPM + Decade + Region FE + Region × {JSL, NAO, wars, deaths}
#
# Outputs:
#   analysis/output/data/figED_v2_onset_excl_swiss.csv
#   analysis/output/data/figED_v2_onset_swiss_years.csv (which years were flagged)
###############################################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(fixest)
})
setFixest_notes(FALSE)

.script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("--file=", "", file_arg)))
  else tryCatch(dirname(normalizePath(sys.frame(1)$ofile)),
                error = function(e) getwd())
}
SCRIPT_DIR <- if (exists("SCRIPT_PATH_HINT")) dirname(SCRIPT_PATH_HINT) else .script_dir()
ANALYSIS_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
ROOT         <- normalizePath(file.path(ANALYSIS_DIR, ".."), mustWork = FALSE)
OUT_DATA     <- file.path(ANALYSIS_DIR, "output", "data")

# ── 1. Identify Swiss-driven CE famine years from v2 yield panel ─────────────
yields <- read.csv(file.path(ROOT, "processed data",
                             "yield_ljungqvist_v2.csv"))
yields$decade <- floor(yields$Year / 10) * 10
swiss_WR <- yields |>
  filter(Country == "Switzerland", Grain %in% c("Wheat", "Rye"))
swiss_z <- swiss_WR |>
  group_by(VarLocationGrain) |>
  mutate(z = (logyield - mean(logyield, na.rm = TRUE)) /
              sd(logyield, na.rm = TRUE)) |>
  ungroup() |>
  group_by(Year) |>
  summarise(z_swiss = mean(z, na.rm = TRUE), n = n(), .groups = "drop")
# NOTE: z-score above is only for *flagging* Swiss-driven famine years;
# the famine onset model below regresses raw Famine_start (not logyield),
# so no outcome transformation is needed.

data <- read.csv(file.path(ROOT, "processed data",
                           "famine_region_data.csv"))
data$Decade <- floor(data$Year / 10) * 10

ce_starts <- data |>
  filter(Region == "Central Europe", Famine_start == 1) |>
  select(Year) |>
  left_join(swiss_z, by = "Year")

swiss_thresh <- -0.5
swiss_years  <- ce_starts |>
  filter(!is.na(z_swiss), z_swiss < swiss_thresh) |>
  pull(Year)

cat("CE famine-start years (n=", nrow(ce_starts), "):\n", sep = "")
print(ce_starts)
cat("\nSwiss-driven CE famine-start years (z_swiss < ", swiss_thresh, "):\n",
    sep = "")
print(swiss_years)

write.csv(ce_starts,
          file.path(OUT_DATA, "figED_v2_onset_swiss_years.csv"),
          row.names = FALSE)

# ── 2. Re-estimate main onset model after dropping Swiss CE famines ──────────
data_excl <- data |>
  mutate(Famine_start = ifelse(Region == "Central Europe" &
                                  Year %in% swiss_years,
                                0, Famine_start))

cat("\nDropped ", length(swiss_years), " Swiss-driven CE famine starts.\n",
    sep = "")
cat("CE Famine_start sum: original=",
    sum(data$Famine_start[data$Region == "Central Europe"]),
    "  excl Swiss=",
    sum(data_excl$Famine_start[data_excl$Region == "Central Europe"]),
    "\n", sep = "")

run_one <- function(df, label_tag) {
  est0 <- feols(Famine_start ~ nino34:Region,
                data = df, vcov = NW ~ Year,
                panel.id = c("Region", "Year"))
  est2 <- feols(Famine_start ~ nino34:Region + i(Decade) | Region,
                data = df, vcov = DK ~ Year,
                panel.id = c("Region", "Year"))
  est3 <- feols(
    Famine_start ~ nino34:Region + i(Decade) +
      JSL:Region + NAO_cal:Region +
      ongoing_wars:Region + log(1 + Deaths):Region | Region,
    data = df, vcov = DK ~ Year,
    panel.id = c("Region", "Year")
  )
  bind_rows(
    broom::tidy(est0, conf.int = TRUE) |> mutate(model = "No FEs"),
    broom::tidy(est2, conf.int = TRUE) |> mutate(model = "FEs"),
    broom::tidy(est3, conf.int = TRUE) |> mutate(model = "FEs + Controls")
  ) |>
    filter(grepl("nino34", term)) |>
    mutate(Region = trimws(gsub("nino34:Region", "", term)),
           sample = label_tag)
}

res_full <- run_one(data,      "Original")
res_excl <- run_one(data_excl, "Excl Swiss famines")
res_all  <- bind_rows(res_full, res_excl)
write.csv(res_all,
          file.path(OUT_DATA, "figED_v2_onset_excl_swiss.csv"),
          row.names = FALSE)
cat("\nSaved figED_v2_onset_excl_swiss.csv\n")
print(res_all |>
        filter(Region == "Central Europe") |>
        select(sample, model, estimate, conf.low, conf.high))
