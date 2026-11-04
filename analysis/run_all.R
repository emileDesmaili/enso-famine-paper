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

# Plots the scripts display on screen (not the saved figures) go to a null
# device, so a non-interactive run does not leave an Rplots.pdf behind.
if (!interactive()) pdf(NULL)

scripts <- c(
  "fig2_onset.R",             # fig 2 (LPM + ML) CSVs + famine_starts_reg.tex
  "fig3_harvest.R",           # fig 3 harvest LP CSVs + all SI harvest figures
  "fig4_survival.R",          # fig 4 (Cox + KM) CSVs + all SI survival figs
  "fig5_prices.R",            # fig 5 (grain/fish prices) CSVs + SI price figs
  "ed_volcanic.R",            # figED_volcanic.pdf
  "ed_nao_jsl_enso.R",        # figED_NAO_JSL_ENSO.pdf
  "si_nonlinearity.R",        # FigA_nonlin.pdf
  "tables.R"                   # SI LaTeX tables
)

run_one <- function(rel_path) {
  path <- file.path(ANALYSIS_DIR, rel_path)
  cat(rep("=", 60), "\n", sep = "")
  cat("Running:", rel_path, "\n")
  cat(rep("=", 60), "\n", sep = "")
  source(path, local = new.env())
  cat("Done:", rel_path, "\n\n")
}

for (s in scripts) run_one(s)

cat("All R scripts complete.\n")
