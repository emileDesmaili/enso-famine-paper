###############################################################################
# Figure 4 – ENSO anomalies and famine duration
# This script exports tidy CSVs for assemble_figures.py
#
# CSV exports → analysis/output/data/
#   fig4A_cox_coefs.csv     – NINO3.4 log-HR and CIs across 6 Cox specs
#   fig4B_survcurves.csv    – Survival curves for 4 ENSO scenarios
#
# Appendix CSVs:
#   figA_surv_strata.csv
#   figA_permutation_survival.csv
#   figA_loo_duration.csv
#
# Tables → analysis/output/tables/
#   dyn_survival_models_reg.tex
#   survival_models_parametric.tex
#   famine_duration_fe_avgnino.tex
#   famine_duration_fe_maxnino.tex
#   summary_famine.tex
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(data.table)
  library(survival)
  library(survminer)
  library(coxme)
  library(fixest)
  library(modelsummary)
  library(stargazer)
  library(broom)
  library(purrr)
  library(scales)
  library(sandwich)
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

# ── data ──────────────────────────────────────────────────────────────────────
load_data <- function() {
  data <- read.csv(file.path(DATA_PROC, "famine_region_data.csv")) |>
    mutate(Decade = floor(Year / 10) * 10) |>
    arrange(Region, Year)

  famine_summary <- read.csv(file.path(DATA_PROC, "famine_survival.csv"))

  dt <- data |> as.data.table()
  dt[, famine_id := rleid(Famine) * Famine]
  dt <- dt[Famine == 1]

  famine_long <- dt |>
    group_by(famine_id, Region) |>
    arrange(Year) |>
    mutate(
      start = row_number() - 1,
      stop  = row_number(),
      event = if_else(stop == max(stop), 1L, 0L)
    ) |>
    ungroup()

  list(data = data, famine_summary = famine_summary, famine_long = famine_long)
}

fit_cox_models <- function(famine_long) {
  list(
    "Base"                      = coxph(Surv(start, stop, event) ~ nino34,
                                        data = famine_long),
    "Strata"                    = coxph(Surv(start, stop, event) ~ nino34 + strata(Region),
                                        data = famine_long),
    "Strata + Clustered S.E."   = coxph(Surv(start, stop, event) ~ nino34 + cluster(famine_id) + strata(Region),
                                        data = famine_long, robust = TRUE),
    "Frailty"                   = coxph(Surv(start, stop, event) ~ nino34 + frailty(Region),
                                        data = famine_long),
    "Controls + Clustered S.E." = coxph(
      Surv(start, stop, event) ~ nino34 + ongoing_wars + temp_winter +
        precip_winter + PDSI + temp_summer + log(1 + Deaths) + cluster(famine_id),
      data = famine_long, robust = TRUE
    ),
    "Controls + Strata"         = coxph(
      Surv(start, stop, event) ~ nino34 + ongoing_wars + temp_winter +
        precip_winter + PDSI + temp_summer + log(1 + Deaths) + strata(Region),
      data = famine_long
    )
  )
}


###############################################################################
# Fig 4A – Cox coefficient data
###############################################################################
export_fig4A <- function(models) {
  bind_rows(lapply(names(models), function(nm) {
    m <- models[[nm]]
    if (!inherits(m, "coxph")) return(NULL)
    coefs  <- summary(m)$coefficients
    if (!"nino34" %in% rownames(coefs)) return(NULL)
    se_col <- if ("robust se" %in% colnames(coefs)) "robust se" else "se(coef)"
    coef_val <- coefs["nino34", "coef"]
    se_val   <- coefs["nino34", se_col]
    data.frame(Model = nm, coef = coef_val, se = se_val,
               lower = coef_val - 1.96 * se_val,
               upper = coef_val + 1.96 * se_val)
  })) |>
    write.csv(file.path(out_data, "fig4A_cox_coefs.csv"), row.names = FALSE)
  message("Exported fig4A_cox_coefs.csv")
}


###############################################################################
# Fig 4B – Survival curve data
###############################################################################
export_fig4B <- function(models, famine_summary) {
  cox_base <- models[["Base"]]
  new_data <- expand.grid(
    nino34 = c(-1, 0, 1, 2),
    Region = unique(famine_summary$Region)
  )
  sf <- summary(survfit(cox_base, newdata = new_data))

  enso_labels <- c("La Niña (-1°)", "Neutral", "El Niño (+1°)", "Strong El Niño (+2°)")
  # survfit columns correspond to rows of new_data in order
  n_per_enso <- length(unique(famine_summary$Region))

  # Average across regions (newdata has 4 ENSO × n_regions rows,
  # survfit gives columns in the same order)
  n_enso <- 4
  surv_mat <- sf$surv          # rows = time points, cols = newdata rows
  times    <- sf$time

  # Mean over region columns for each ENSO level
  surv_avg <- do.call(cbind, lapply(seq_len(n_enso), function(e) {
    idx <- seq(e, ncol(surv_mat), by = n_enso)
    rowMeans(surv_mat[, idx, drop = FALSE])
  }))
  colnames(surv_avg) <- enso_labels

  as.data.frame(surv_avg) |>
    mutate(time = times) |>
    pivot_longer(-time, names_to = "ENSO", values_to = "surv") |>
    write.csv(file.path(out_data, "fig4B_survcurves.csv"), row.names = FALSE)
  message("Exported fig4B_survcurves.csv")
}


###############################################################################
# Appendix: stratified survival curves
###############################################################################
export_surv_strata <- function(models, famine_summary) {
  cox_strata <- models[["Strata"]]
  new_data   <- expand.grid(
    nino34 = c(-1, 0, 1),
    Region = unique(famine_summary$Region)
  )
  sf      <- summary(survfit(cox_strata, newdata = new_data))
  n_cols  <- ncol(sf$surv)
  lengths <- as.numeric(table(sf$strata))
  repeated <- new_data[rep(seq_len(nrow(new_data)), times = lengths), ]

  data.frame(time = sf$time, surv = sf$surv,
             lower = sf$lower, upper = sf$upper) |>
    bind_cols(repeated) |>
    mutate(ENSO = case_when(nino34 == -1 ~ "La Niña",
                            nino34 == 0  ~ "Neutral",
                            nino34 == 1  ~ "El Niño")) |>
    write.csv(file.path(out_data, "figA_surv_strata.csv"), row.names = FALSE)
  message("Exported figA_surv_strata.csv")
}


###############################################################################
# Appendix: permutation test (survival)
###############################################################################
export_permutation_survival <- function(famine_long, n_iter = 10000) {
  set.seed(42)
  formula    <- Surv(start, stop, event) ~ nino34 + strata(Region)
  true_model <- coxph(formula, data = famine_long)
  true_coef  <- coef(true_model)[["nino34"]]

  coef_vec <- vapply(seq_len(n_iter), function(i) {
    df_sh <- famine_long |> mutate(nino34 = sample(nino34))
    m     <- try(coxph(formula, data = df_sh), silent = TRUE)
    if (inherits(m, "try-error")) return(NA_real_)
    coef(m)[["nino34"]]
  }, numeric(1))

  data.frame(estimate = coef_vec, true_estimate = rep(true_coef, length(coef_vec))) |>
    write.csv(file.path(out_data, "figA_permutation_survival.csv"), row.names = FALSE)
  message("Exported figA_permutation_survival.csv")
}


###############################################################################
# Appendix: leave-one-out (duration)
###############################################################################
export_loo_duration <- function(famine_long) {
  cox_full  <- coxph(Surv(start, stop, event) ~ nino34 + strata(Region),
                     data = famine_long)
  full_coef <- coef(cox_full)[["nino34"]]
  full_se   <- sqrt(vcov(cox_full)["nino34","nino34"])

  map_dfr(unique(famine_long$famine_id), function(fam_id) {
    m <- coxph(Surv(start, stop, event) ~ nino34 + strata(Region),
               data = famine_long |> filter(famine_id != fam_id))
    broom::tidy(m) |>
      filter(term == "nino34") |>
      mutate(left_out = fam_id,
             baseline = full_coef, baseline_se = full_se)
  }) |>
    write.csv(file.path(out_data, "figA_loo_duration.csv"), row.names = FALSE)
  message("Exported figA_loo_duration.csv")
}


###############################################################################
# Tables
###############################################################################
save_survival_tables <- function(famine_long, famine_summary) {
  cox_models <- list(
    "(1)" = coxph(Surv(start, stop, event) ~ nino34, data = famine_long),
    "(2)" = coxph(Surv(start, stop, event) ~ nino34 + strata(Region),
                  data = famine_long),
    "(3)" = coxph(Surv(start, stop, event) ~ nino34 + cluster(famine_id) + strata(Region),
                  data = famine_long, robust = TRUE),
    "(4)" = coxph(Surv(start, stop, event) ~ nino34 + frailty(Region),
                  data = famine_long),
    "(5)" = coxph(
      Surv(start, stop, event) ~ nino34 + ongoing_wars + temp_winter +
        precip_winter + PDSI + temp_summer + log(1 + Deaths) + cluster(famine_id),
      data = famine_long, robust = TRUE
    ),
    "(6)" = coxph(
      Surv(start, stop, event) ~ nino34 + ongoing_wars + temp_winter +
        precip_winter + PDSI + temp_summer + log(1 + Deaths) + strata(Region),
      data = famine_long
    )
  )

  modelsummary(
    cox_models,
    title    = "\\label{tab:dynsurvENSO}Log Hazard Ratio of Nino 3.4 on Famine Duration",
    escape   = FALSE,
    coef_omit = "^(?!.*nino34).*",
    gof_omit  = ".*",
    stars    = c(`*`=0.1, `**`=0.05, `***`=0.01),
    output   = file.path(tab_out, "dyn_survival_models_reg.tex")
  )
  message("Saved dyn_survival_models_reg.tex")

  weib   <- survreg(Surv(duration) ~ avg_nino34 + strata(Region),
                    data = famine_summary, dist = "weibull")
  loglog <- survreg(Surv(duration) ~ avg_nino34 + strata(Region),
                    data = famine_summary, dist = "loglogistic")
  gauss  <- survreg(Surv(duration) ~ avg_nino34 + strata(Region),
                    data = famine_summary, dist = "gaussian")

  stargazer(weib, loglog, gauss, type = "latex",
            title = "Parametric Survival Models of Famine Duration",
            dep.var.labels = "Duration",
            covariate.labels = "Nino34",
            single.row = TRUE,
            star.char = c("*","**","***"),
            star.cutoffs = c(0.1, 0.05, 0.01),
            out = file.path(tab_out, "survival_models_parametric.tex"))
  message("Saved survival_models_parametric.tex")

  dict <- c(avg_nino34 = "Nino3.4", max_2y_nino = "Max Nino3.4 (2y)",
            avg_ongoing_wars = "Ongoing wars", duration = "Duration (years)")

  for (nino_var in c("avg_nino34", "max_2y_nino")) {
    ols1  <- feols(duration ~ get(nino_var) | Region,
                   data = famine_summary, cluster = ~Region)
    ols2  <- feols(duration ~ get(nino_var) + avg_ongoing_wars +
                     avg_temp_summer + avg_temp_winter + avg_PDSI +
                     avg_precip_winter | Region,
                   data = famine_summary, cluster = ~Region)
    pois1 <- fepois(duration ~ get(nino_var) | Region,
                    data = famine_summary, cluster = ~Region)
    pois2 <- fepois(duration ~ get(nino_var) + avg_ongoing_wars +
                      avg_temp_summer + avg_temp_winter + avg_PDSI +
                      avg_precip_winter | Region,
                    data = famine_summary, cluster = ~Region)
    nb1   <- fenegbin(duration ~ get(nino_var) | Region,
                      data = famine_summary, cluster = ~Region)
    nb2   <- fenegbin(duration ~ get(nino_var) + avg_ongoing_wars +
                        avg_temp_summer + avg_temp_winter + avg_PDSI +
                        avg_precip_winter | Region,
                      data = famine_summary, cluster = ~Region)

    fname <- if (nino_var == "avg_nino34") "famine_duration_fe_avgnino.tex"
             else "famine_duration_fe_maxnino.tex"
    etable(ols1, ols2, pois1, pois2, nb1, nb2,
           keep = "%nino", dict = dict,
           label = paste0("tab:duration_", nino_var),
           title = paste("Effect of a 1\\textdegree C Anomaly in the",
                         ifelse(nino_var == "avg_nino34", "Average", "Max 2-year"),
                         "Nino 3.4 Index on Famine Duration"),
           drop.section = "fixef", fitstat = c("AIC","BIC","N","cor2"),
           extralines = list(
             "Controls"             = c("No","Yes","No","Yes","No","Yes"),
             "Region Fixed Effects" = c("Yes","Yes","Yes","Yes","Yes","Yes")
           ),
           file = file.path(tab_out, fname), replace = TRUE)
    message("Saved ", fname)
  }

  datasummary(
    duration * Region ~ mean + sd + min + max + N,
    data   = famine_summary,
    output = file.path(tab_out, "summary_famine.tex"),
    title  = "Summary Statistics of Famine Episodes"
  )
  message("Saved summary_famine.tex")
}


###############################################################################
# Appendix plots
###############################################################################
fig_out <- file.path(SCRIPT_DIR, "output", "figures", "appendix")
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

plot_surv_strata <- function() {
  path <- file.path(out_data, "figA_surv_strata.csv")
  if (!file.exists(path)) { message("figA_surv_strata.csv not found"); return() }
  df <- read.csv(path) |>
    mutate(ENSO = factor(ENSO, levels = c("La Niña", "Neutral", "El Niño")))
  enso_colors <- c("La Niña" = "cornflowerblue", "Neutral" = "black", "El Niño" = "coral")
  enso_shapes <- c("La Niña" = 16, "Neutral" = 17, "El Niño" = 18)
  p <- ggplot(df, aes(x = time, y = surv, color = ENSO, shape = ENSO)) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2) +
    facet_wrap(~ Region, scales = "free", ncol = 3) +
    scale_color_manual(values = enso_colors) +
    scale_shape_manual(values = enso_shapes) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(breaks = 1:12) +
    labs(x = "Famine duration (years)", y = "Survival probability",
         color = "ENSO anomaly", shape = "ENSO anomaly") +
    theme_classic(base_size = 13) +
    theme(strip.background = element_blank(),
          strip.text       = element_text(face = "bold"),
          panel.grid       = element_blank(),
          legend.position  = "top")
  ggsave(file.path(fig_out, "figA_surv_strata.pdf"), p,
         width = 12, height = 8, device = cairo_pdf)
  message("Saved figA_surv_strata.pdf")
}

plot_permutation_survival <- function() {
  path <- file.path(out_data, "figA_permutation_survival.csv")
  if (!file.exists(path)) { message("figA_permutation_survival.csv not found"); return() }
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
    labs(x = "Permuted NINO3.4 log hazard ratio",
         y = "Density") +
    theme_classic(base_size = 13) +
    theme(panel.grid = element_blank())
  ggsave(file.path(fig_out, "figA_permutation_cox.pdf"), p,
         width = 7, height = 5, device = cairo_pdf)
  message("Saved figA_permutation_cox.pdf")
}

plot_loo_duration <- function() {
  path <- file.path(out_data, "figA_loo_duration.csv")
  if (!file.exists(path)) { message("figA_loo_duration.csv not found"); return() }
  df <- read.csv(path) |>
    arrange(estimate) |>
    mutate(i     = seq_len(n()),
           lower = estimate - 1.96 * std.error,
           upper = estimate + 1.96 * std.error)
  p <- ggplot(df, aes(x = i, y = estimate, ymin = lower, ymax = upper)) +
    geom_hline(aes(yintercept = baseline), linetype = "dashed",
               color = "red", linewidth = 0.8) +
    geom_pointrange(color = "grey", size = 0.4) +
    labs(x = "Famine episode left out",
         y = "NINO3.4 log hazard ratio") +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank())
  ggsave(file.path(fig_out, "figA_loo_duration.pdf"), p,
         width = 9, height = 5, device = cairo_pdf)
  message("Saved figA_loo_duration.pdf")
}


###############################################################################
# Run
###############################################################################
if (!interactive()) {
  lst            <- load_data()
  data           <- lst$data
  famine_summary <- lst$famine_summary
  famine_long    <- lst$famine_long

  models <- fit_cox_models(famine_long)

  export_fig4A(models)
  export_fig4B(models, famine_summary)
  export_surv_strata(models, famine_summary)
  export_permutation_survival(famine_long, n_iter=100)
  export_loo_duration(famine_long)
  save_survival_tables(famine_long, famine_summary)
  plot_surv_strata()
  plot_permutation_survival()
  plot_loo_duration()

  message("Fig 4 data exports complete.")
} else {
  message("Call functions interactively or source this script.")
}
