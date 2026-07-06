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
  "fig2_onset.R",             # fig 2 (LPM + ML) CSVs + famine_starts_reg.tex
  "fig3_harvest.R",           # fig 3 harvest LP CSVs + all SI harvest figures
  "fig4_survival.R",          # fig 4 (Cox + KM) CSVs + all SI survival figs
  "fig5_prices.R",            # fig 5 (grain/fish prices) CSVs + SI price figs
  "ed_volcanic.R",            # figED_volcanic.pdf
  "ed_nao_jsl_enso.R",        # figED_NAO_JSL_ENSO.pdf
  "si_nonlinearity.R",        # FigA_nonlin.pdf
  "tables.R"                   # SI LaTeX tables
)

# Investigation / rebuttal-only scripts (kept separate from the main pipeline
# so the paper build stays clean). Their outputs land in
# analysis/output/figures/investigation/.
investigation_scripts <- c(
  "investigation/v2_onset_lag_lead.R",   # CSV feed for v2 onset panel
  "investigation/v2_onset_excl_swiss.R", # CSV feed for v2 onset panel
  "investigation/r2_teleco_placebo.R",   # fig3c_robust_and_placebo.pdf (R2)
  "investigation/r2_geo_placebo.R"       # figA_geopartition_robust_and_placebo.pdf (R2)
)

run_one <- function(rel_path) {
  path <- file.path(ANALYSIS_DIR, rel_path)
  cat(rep("=", 60), "\n", sep = "")
  cat("Running:", rel_path, "\n")
  cat(rep("=", 60), "\n", sep = "")
  env <- new.env()
  # Hint each script about its own on-disk location, so scripts in
  # subfolders (analysis/investigation/) don't misresolve their paths.
  assign("SCRIPT_PATH_HINT", path, envir = env)
  source(path, local = env)
  cat("Done:", rel_path, "\n\n")
}

for (s in scripts)               run_one(s)
for (s in investigation_scripts) run_one(s)

cat("All R scripts complete.\n")
