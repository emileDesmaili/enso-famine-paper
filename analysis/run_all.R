# Master runner – sources all R analysis scripts in order.
# Run from the repo root with:
#   Rscript analysis/run_all.R
# or interactively:
#   source("analysis/run_all.R")
#
# Outputs land in:
#   analysis/output/figures/main/
#   analysis/output/figures/appendix/
#   analysis/output/figures/extended data/
#   analysis/output/tables/

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

ANALYSIS_DIR <- .script_dir()
ROOT_DIR     <- normalizePath(file.path(ANALYSIS_DIR, ".."), mustWork = FALSE)

scripts <- c(
  "02_fig2_onset.R",
  "03_fig3_harvest.R",
  "04_fig4_survival.R",
  "05_fig5_prices.R",
  "29_v2_onset_lag_lead.R",
  "30_v2_onset_excl_swiss.R",
  "34_ed_volcanic.R",         # figED_volcanic.pdf
  "35_ed_nao_jsl_enso.R",     # figED_NAO_JSL_ENSO.pdf
  "38_si_nonlinearity.R",
  "tables.R"
)

for (s in scripts) {
  path <- file.path(ANALYSIS_DIR, s)
  cat(rep("=", 60), "\n", sep = "")
  cat("Running:", s, "\n")
  cat(rep("=", 60), "\n", sep = "")
  source(path, local = new.env())
  cat("Done:", s, "\n\n")
}

cat("All R scripts complete.\n")
