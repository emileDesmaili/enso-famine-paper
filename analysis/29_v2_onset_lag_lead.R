###############################################################################
# 29_v2_onset_lag_lead.R
#
# Famine onset model with Region × {lag, contemp, lead} NINO3.4. Report
# coefficients on Central Europe for each timing.
#
# Spec (LPM):
#   Famine_start ~ i(Region, l(nino34, 1))     # one lag
#                + i(Region, nino34)            # contemporaneous
#                + i(Region, f(nino34, 1))     # one lead
#                + i(Decade) | Region
# DK SE, panel.id = (Region, Year).
#
# Output:
#   analysis/output/data/figED_v2_onset_lag_lead.csv
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
SCRIPT_DIR <- .script_dir()
ROOT       <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
OUT_DATA   <- file.path(SCRIPT_DIR, "output", "data")

data <- read.csv(file.path(ROOT, "processed data",
                           "famine_region_data.csv"))
data$Decade <- floor(data$Year / 10) * 10
data <- data |> arrange(Region, Year)

# Three models: (no FE) / (Decade+Region FE) / (FE + region-interacted controls)
fit_one <- function(mod_label, fml, vcov_arg) {
  mod <- feols(fml, data = data,
               panel.id = c("Region", "Year"),
               vcov     = vcov_arg)
  co  <- coef(mod)
  se  <- sqrt(diag(vcov(mod)))
  cn  <- names(co)
  pat <- grep("Region.*Central Europe.*nino34|nino34.*Central Europe",
              cn, value = TRUE)
  bind_rows(lapply(pat, function(p) {
    lag_type <- if (grepl("l\\(nino34", p)) "Lag (t-1)"
                else if (grepl("f\\(nino34", p)) "Lead (t+1)"
                else "Contemporaneous (t)"
    data.frame(model = mod_label, term = p, timing = lag_type,
               estimate = co[p], se = se[p],
               conf.low  = co[p] - 1.96 * se[p],
               conf.high = co[p] + 1.96 * se[p])
  }))
}

m1 <- fit_one("No FEs",
              Famine_start ~ i(Region, l(nino34, 1))
                           + i(Region, nino34)
                           + i(Region, f(nino34, 1)),
              NW ~ Year)
m2 <- fit_one("FEs",
              Famine_start ~ i(Region, l(nino34, 1))
                           + i(Region, nino34)
                           + i(Region, f(nino34, 1))
                           + i(Decade) | Region,
              DK ~ Year)
m3 <- fit_one("FEs + Controls",
              Famine_start ~ i(Region, l(nino34, 1))
                           + i(Region, nino34)
                           + i(Region, f(nino34, 1))
                           + i(Decade)
                           + JSL:Region + NAO_cal:Region
                           + ongoing_wars:Region
                           + log(1 + Deaths):Region | Region,
              DK ~ Year)

res <- bind_rows(m1, m2, m3)
res$timing <- factor(res$timing,
                     levels = c("Lag (t-1)", "Contemporaneous (t)",
                                "Lead (t+1)"))
write.csv(res, file.path(OUT_DATA, "figED_v2_onset_lag_lead.csv"),
          row.names = FALSE)
cat("Saved figED_v2_onset_lag_lead.csv\n")
print(res |> select(model, timing, estimate, se, conf.low, conf.high))
