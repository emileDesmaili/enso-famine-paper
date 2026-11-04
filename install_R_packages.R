# Installs the R packages used by analysis/*.R.
# The published results were produced with R 4.5.1 and the versions below.
#
#   Rscript install_R_packages.R          # install any missing package (current CRAN version)
#   Rscript install_R_packages.R --exact  # install the exact versions listed below

pkgs <- c(
  fixest       = "0.14.0",   # linear probability, logit and local-projection models
  survival     = "3.8-3",    # Cox proportional hazards models
  sandwich     = "3.1-1",
  dplyr        = "1.1.4",
  tidyr        = "1.3.1",
  tibble       = "3.3.0",
  purrr        = "1.1.0",
  furrr        = "0.3.1",
  future       = "1.67.0",
  data.table   = "1.18.2.1",
  stringr      = "1.5.1",
  glue         = "1.8.0",
  haven        = "2.5.5",
  readxl       = "1.4.5",
  broom        = "1.0.12",
  modelsummary = "2.4.0",
  stargazer    = "5.2.3",
  ggplot2      = "3.5.2",
  cowplot      = "1.2.0",
  patchwork    = "1.3.2",
  ggpubr       = "0.6.1",
  scales       = "1.4.0"
)

repos <- "https://cloud.r-project.org"
exact <- "--exact" %in% commandArgs(trailingOnly = TRUE)

if (exact) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = repos)
  for (p in names(pkgs)) {
    have <- tryCatch(packageVersion(p), error = function(e) NULL)
    if (is.null(have) || have != package_version(pkgs[[p]]))
      remotes::install_version(p, version = pkgs[[p]], repos = repos, upgrade = "never")
  }
} else {
  missing <- names(pkgs)[!vapply(names(pkgs), requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) install.packages(missing, repos = repos)
}

installed <- vapply(names(pkgs), function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING"), "")
print(data.frame(package = names(pkgs), used_in_paper = pkgs, installed = installed, row.names = NULL))
