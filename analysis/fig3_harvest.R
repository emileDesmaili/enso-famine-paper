###############################################################################
# Figure 3 – ENSO teleconnections and grain harvest responses
# This script exports tidy CSVs for assemble_figures.py
#
# CSV exports → analysis/output/data/
#   fig3C_irf_yield_WR.csv    – Wheat/Rye harvest IRF
#   fig3D_irf_yield_noWR.csv  – Other grain harvest IRF
#
# Appendix CSVs:
#   figA_irf_yield63.csv
#   figA_irf_subsample_WR.csv
#   figA_irf_20y_WR.csv
#   figA_irf_WR_bootstrap_null.csv / _true.csv
#
# Tables → analysis/output/tables/
#   FE_ENSO_yieldPDSI.tex
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(scales)
  library(fixest)
  library(haven)
  library(purrr)
  library(furrr)
  library(future)
})

setFixest_notes(FALSE)

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

# ── load emileRegs helpers ────────────────────────────────────────────────────
source(file.path(ROOT, "emileRegs.R"))

# ── helpers ───────────────────────────────────────────────────────────────────
run_lp <- function(df, endog, shock = "nino34",
                   controls = "i(decade)", hor = 3) {
  lp_panel(
    data          = df,
    outcome       = endog,
    main_var      = shock,
    controls      = if (identical(controls, "decade")) "i(decade)"
                    else paste(c("i(decade)", setdiff(controls, "decade")),
                               collapse = " + "),
    horizon       = hor,
    fe            = "VarLocationGrain",
    panel_id      = c("VarLocationGrain", "Year"),
    vcov_formula  = DK ~ Year
  )
}

extract_irf <- function(res, label, n_loc = NA_integer_) {
  res |>
    transmute(label = label, n_loc = n_loc, horizon, irf_mean, irf_up, irf_down)
}

irf_teleco_loop <- function(df, endog, hor_max,
                             teleco_col     = "teleco_PDSI_10",
                             controls_full  = c("decade", "JSL", "NAO_cal",
                                                "ongoing_wars", "Deaths")) {
  conditions <- list(
    "All regions"        = NULL,
    "Weakly affected"    = paste0(teleco_col, " == 0"),
    "Teleconnected"      = paste0(teleco_col, " == 1"),
    "Teleco w. controls" = paste0(teleco_col, " == 1")
  )
  bind_rows(lapply(names(conditions), function(label) {
    cond   <- conditions[[label]]
    df_sub <- if (!is.null(cond)) df |> filter(!!rlang::parse_expr(cond)) else df
    ctrl   <- if (label == "Teleco w. controls") controls_full else "decade"
    n_loc  <- dplyr::n_distinct(df_sub$VarLocationGrain)
    extract_irf(run_lp(df_sub, endog, controls = ctrl, hor = hor_max), label, n_loc)
  }))
}

load_yields <- function() {
  # Uses the v2 Ljungqvist yield panel.
  read.csv(file.path(DATA_PROC, "yield_ljungqvist_v2.csv")) |>
    mutate(decade = floor(Year / 10) * 10) |>
    dplyr::select(VarLocationGrain, Year, decade, everything())
}

# ── v2 Ljungqvist yield panel loader (raw logyield, decade + Year features).
#    Used for appendix SI figures. β reads directly as % harvest.
load_yields_GT3 <- function() {
  read.csv(file.path(DATA_PROC, "yield_ljungqvist_v2.csv")) |>
    mutate(decade = floor(Year / 10) * 10,
           Year2  = Year^2,
           Year3  = Year^3) |>
    dplyr::select(VarLocationGrain, Year, decade, Year2, Year3, everything())
}

FE_GT3 <- "VarLocationGrain[Year]"   # linear per-record time trend

run_lp_GT3 <- function(df, endog = "logyield", shock = "nino34",
                       extra_controls = NULL, hor = 3) {
  ctrl <- c("i(decade)", paste0("l(", shock, ", 1:3)"))
  if (!is.null(extra_controls) && length(extra_controls) > 0)
    ctrl <- c(ctrl, extra_controls)
  lp_panel(
    data         = df,
    outcome      = endog,
    main_var     = shock,
    controls     = paste(ctrl, collapse = " + "),
    horizon      = hor,
    fe           = FE_GT3,
    panel_id     = c("VarLocationGrain", "Year"),
    vcov_formula = DK ~ Year
  )
}

# Interaction LP — both teleco categories estimated jointly on full sample.
run_lp_GT3_inter <- function(df, endog = "logyield", shock = "nino34",
                              extra_controls = NULL, hor = 3,
                              interact_var = "teleco_PDSI_10") {
  ctrl <- c("i(decade)", paste0("l(", shock, ", 1:3)"))
  if (!is.null(extra_controls) && length(extra_controls) > 0)
    ctrl <- c(ctrl, extra_controls)
  lp_panel_inter(
    data         = df,
    outcome      = endog,
    main_var     = shock,
    interact_var = interact_var,
    controls     = paste(ctrl, collapse = " + "),
    horizon      = hor,
    fe           = FE_GT3,
    panel_id     = c("VarLocationGrain", "Year"),
    vcov_formula = DK ~ Year
  )
}

# Teleconnection partition: pooled "All regions" + interaction split by teleco
# + "Teleco w. controls" (interaction + extra controls).
irf_teleco_loop_GT3 <- function(df, endog = "logyield", hor_max = 3,
                                 teleco_col = "teleco_PDSI_10",
                                 controls_full = c("JSL", "NAO_cal",
                                                   "ongoing_wars", "Deaths")) {
  res_all <- run_lp_GT3(df, endog, hor = hor_max)
  out_all <- extract_irf(res_all, "All regions",
                          dplyr::n_distinct(df$VarLocationGrain))

  res_int <- run_lp_GT3_inter(df, endog, hor = hor_max,
                                interact_var = teleco_col)
  n0 <- dplyr::n_distinct(df$VarLocationGrain[df[[teleco_col]] == 0])
  n1 <- dplyr::n_distinct(df$VarLocationGrain[df[[teleco_col]] == 1])
  out_w <- res_int |> filter(category == "0") |>
    transmute(label = "Weakly affected", n_loc = n0,
              horizon, irf_mean, irf_up, irf_down)
  out_t <- res_int |> filter(category == "1") |>
    transmute(label = "Teleconnected", n_loc = n1,
              horizon, irf_mean, irf_up, irf_down)

  res_int_c <- run_lp_GT3_inter(df, endog,
                                  extra_controls = controls_full,
                                  hor = hor_max,
                                  interact_var = teleco_col)
  out_tc <- res_int_c |> filter(category == "1") |>
    transmute(label = "Teleco w. controls", n_loc = n1,
              horizon, irf_mean, irf_up, irf_down)

  bind_rows(out_all, out_w, out_t, out_tc)
}


###############################################################################
# Fig 3C – Wheat/Rye harvest IRF (linear record trend, interaction with
# teleco_PDSI_10, 3 ENSO lags, v2 panel via load_yields_GT3)
###############################################################################
export_fig3C <- function(yield_2023) {
  WR <- load_yields_GT3() |> filter(Grain %in% c("Wheat", "Rye"))
  irf_teleco_loop_GT3(WR, "logyield", 3) |>
    write.csv(file.path(out_data, "fig3C_irf_yield_WR.csv"), row.names = FALSE)
  message("Exported fig3C_irf_yield_WR.csv")
}


###############################################################################
# Fig 3D – Other-grain harvest IRF (same spec as fig 3C)
###############################################################################
export_fig3D <- function(yield_2023) {
  noWR <- load_yields_GT3() |> filter(!Grain %in% c("Wheat", "Rye"))
  irf_teleco_loop_GT3(noWR, "logyield", 3) |>
    write.csv(file.path(out_data, "fig3D_irf_yield_noWR.csv"), row.names = FALSE)
  message("Exported fig3D_irf_yield_noWR.csv")
}


###############################################################################
# Appendix: SE robustness – Wheat/Rye IRF under 4 vcov types
###############################################################################
# SE types:
#   DK          – Driscoll-Kraay (baseline)
#   2way        – two-way cluster (VarLocationGrain × Year)
#   Conley      – spatial HAC via vcov_conley
#   Spatial HAC – Conley + NW - Hetero (combined spatial+serial+correct)
#
# Per-horizon model fitted with feols(), then vcov swapped per SE type.
run_lp_se_robust <- function(df, endog, shock = "nino34",
                              controls = "i(decade) + l(nino34, 1:3)",
                              hor = 3,
                              lat_col = "Latitude", lon_col = "Longitude") {
  fml_rhs <- paste0(shock, " + ", controls)
  fe_str  <- "VarLocationGrain[Year]"   # linear per-record trend
  pid     <- c("VarLocationGrain", "Year")

  bind_rows(lapply(0:hor, function(h) {
    fml <- as.formula(paste0(
      "f(", endog, ",", h, ") - l(", endog, ",1) ~ ",
      fml_rhs, " | ", fe_str
    ))
    mod <- feols(fml, data = df, panel.id = pid, vcov = DK ~ Year)

    # Spatial Conley vcov (lat/lon must be in data)
    v_conley <- tryCatch(
      vcov(mod, vcov = vcov_conley(lat = lat_col, lon = lon_col,
                                   cutoff = 500, distance = "spherical")),
      error = function(e) NULL
    )
    v_nw     <- vcov(mod, vcov = NW ~ Year)
    v_het    <- vcov(mod, vcov = "hetero")
    v_spat   <- if (!is.null(v_conley)) v_conley + v_nw - v_het else NULL

    se_list <- list(
      "DK"          = vcov(mod, vcov = DK ~ Year),
      "2-way"       = vcov(mod, vcov = cluster ~ VarLocationGrain + Year),
      "Conley"      = v_conley,
      "Spatial HAC" = v_spat
    )

    bind_rows(lapply(names(se_list), function(nm) {
      V <- se_list[[nm]]
      if (is.null(V)) return(NULL)
      b  <- coef(mod)[shock]
      se <- sqrt(V[shock, shock])
      data.frame(se_type  = nm, horizon = h,
                 irf_mean = b,
                 irf_up   = b + 1.96 * se,
                 irf_down = b - 1.96 * se)
    }))
  }))
}

export_se_robust_WR <- function(yield_2023) {
  # Appendix uses v2 Ljungqvist panel + main spec
  WR <- load_yields_GT3() |>
    filter(Grain %in% c("Wheat", "Rye"))
  run_lp_se_robust(WR, "logyield", hor = 3) |>
    write.csv(file.path(out_data, "figA_irf_WR_serobust.csv"), row.names = FALSE)
  message("Exported figA_irf_WR_serobust.csv")
}


###############################################################################
# Appendix: Yield63 IRF
###############################################################################
export_yield63_irf <- function() {
  data_region <- read.csv(file.path(DATA_PROC, "famine_region_data.csv"))
  dta_path    <- file.path(ROOT, "processed data", "yieldratio_dataset.dta")
  if (!file.exists(dta_path)) {
    message("  yieldratio_dataset.dta not found – skipping yield63")
    return(invisible(NULL))
  }
  region_year <- data_region |>
    dplyr::select(Year, nino34, PDSI_europe = PDSI,
                  Deaths, ongoing_wars, precip_summer_europe = precip_summer) |>
    dplyr::distinct(Year, .keep_all = TRUE)

  yield63 <- haven::read_dta(dta_path) |>
    left_join(region_year, by = c("year" = "Year")) |>
    mutate(decade = floor(year / 10) * 10, logyield = log(yield)) |>
    distinct() |> drop_na() |>
    filter(unique_sumtotal >= 10, year <= 1800) |>
    rename(Year = year)

  ctrl_dict <- list(
    "Base"               = "decade",
    "Climate"            = c("decade", "growseasontemp", "ngrowseasontempn1",
                             "growseasonrain", "ngrowseasonrainn1"),
    "Conflict"           = c("decade", "ongoing_wars", "Deaths"),
    "Climate + Conflict" = c("decade", "growseasontemp", "ngrowseasontempn1",
                             "growseasonrain", "ngrowseasonrainn1",
                             "ongoing_wars", "Deaths")
  )
  bind_rows(lapply(names(ctrl_dict), function(nm) {
    ctrl <- ctrl_dict[[nm]]
    ctrl_fml <- paste(c("i(decade)", setdiff(ctrl, "decade")), collapse = " + ")
    lp_panel(
      data         = yield63,
      outcome      = "logyield",
      main_var     = "nino34",
      controls     = ctrl_fml,
      horizon      = 4,
      fe           = "fid",
      panel_id     = c("fid", "Year"),
      vcov_formula = DK ~ Year
    ) |>
      transmute(controls = nm, horizon, irf_mean, irf_up, irf_down)
  })) |>
    write.csv(file.path(out_data, "figA_irf_yield63.csv"), row.names = FALSE)
  message("Exported figA_irf_yield63.csv")
}


###############################################################################
# Appendix: LOO subsampling (WR)
###############################################################################
export_loo_WR <- function(yield_2023) {
  WR <- load_yields_GT3() |>
    filter(Grain %in% c("Wheat", "Rye"))
  locations <- unique(WR$loc_id)
  combos    <- combn(locations, length(locations) - 1, simplify = FALSE)

  bind_rows(lapply(seq_along(combos), function(i) {
    left_out <- setdiff(locations, combos[[i]])
    run_lp_GT3(WR |> filter(loc_id %in% combos[[i]]),
               "logyield", hor = 3) |>
      transmute(subsample = paste("No", left_out), horizon, irf_mean, irf_up, irf_down)
  })) |>
    write.csv(file.path(out_data, "figA_irf_subsample_WR.csv"), row.names = FALSE)
  message("Exported figA_irf_subsample_WR.csv")
}


###############################################################################
# Appendix: 20-year block LOO (WR)
###############################################################################
export_20y_WR <- function(yield_2023) {
  WR <- load_yields_GT3() |>
    filter(Grain %in% c("Wheat", "Rye"))
  years        <- sort(unique(WR$Year))
  block_starts <- seq(min(years), max(years) - 20, by = 20)

  bind_rows(lapply(block_starts, function(bs) {
    block_yrs <- bs:(bs + 19)
    run_lp_GT3(WR |> filter(!Year %in% block_yrs),
               "logyield", hor = 3) |>
      transmute(block = paste(bs, bs + 19, sep = "-"), horizon, irf_mean, irf_up, irf_down)
  })) |>
    write.csv(file.path(out_data, "figA_irf_20y_WR.csv"), row.names = FALSE)
  message("Exported figA_irf_20y_WR.csv")
}


###############################################################################
# Appendix: Bootstrap permutation (WR)
###############################################################################
export_bootstrap_WR <- function(yield_2023, n_iter = 1000) {
  set.seed(42)
  plan(multisession)
  WR <- load_yields_GT3() |>
    filter(Grain %in% c("Wheat", "Rye"))
  true_res <- run_lp_GT3(WR, "logyield", hor = 3)

  null_list <- future_map(seq_len(n_iter), function(i) {
    df_sh <- WR |> mutate(nino34 = sample(nino34))
    res   <- try(run_lp_GT3(df_sh, "logyield", hor = 3), silent = TRUE)
    if (inherits(res, "try-error")) return(NULL)
    res |> transmute(iter = i, horizon, irf_mean)
  }, .progress = TRUE) |>
    dplyr::bind_rows()

  null_list |>
    write.csv(file.path(out_data, "figA_irf_WR_bootstrap_null.csv"),
              row.names = FALSE)
  true_res |>
    transmute(horizon, irf_mean) |>
    write.csv(file.path(out_data, "figA_irf_WR_bootstrap_true.csv"),
              row.names = FALSE)
  message("Exported figA_irf_WR_bootstrap CSVs")
}


###############################################################################
# Appendix: Yield IRF under different teleconnection partitions
# Reproduces: IRF_ENSO_yield_WR_teleco, IRF_ENSO_yield_noWR_teleco,
#             figA_irf_yield_WR_AMJJ_P, figA_irf_yield_NDJF_P,
#             figA_irf_yield_NDJF_T, figA_irf_yield_AMJJ_T
###############################################################################
export_yield_teleco_plots <- function(yield_2023) {
  yields_GT3 <- load_yields_GT3()
  WR   <- yields_GT3 |>
    filter(Grain %in% c("Wheat", "Rye"))
  noWR <- yields_GT3 |>
    filter(!Grain %in% c("Wheat", "Rye"))

  plot_teleco <- function(df, teleco_col, title_lab, file_name) {
    df_irf <- irf_teleco_loop_GT3(df, "logyield", hor_max = 3,
                                   teleco_col = teleco_col)
    n_ctrl  <- length(unique(df_irf$label))
    offsets <- seq(-0.12, 0.12, length.out = n_ctrl)
    df_irf  <- df_irf |>
      dplyr::group_by(label) |>
      dplyr::mutate(offset = offsets[dplyr::cur_group_id()]) |>
      dplyr::ungroup()

    p <- ggplot(df_irf, aes(x = horizon + offset, y = irf_mean, color = label)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
      geom_errorbar(aes(ymin = irf_down, ymax = irf_up), width = 0, linewidth = 0.9) +
      geom_point(size = 2.5) +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
      scale_color_manual(values = c(
        "All regions"        = "black",
        "Teleconnected"      = "firebrick",
        "Weakly affected"    = "cornflowerblue",
        "Teleco w. controls" = "orange"
      )) +
      labs(x = "Horizon (years)", y = title_lab, color = NULL) +
      theme_classic(base_size = 18) +
      theme(panel.grid = element_blank(), legend.position = "bottom") +
      guides(color = guide_legend(nrow = 2, byrow = TRUE))

    ggsave(file.path(fig_out, file_name), p,
           width = 6, height = 5, device = cairo_pdf)
    message("Saved ", file_name)
    invisible(df_irf)
  }

  # PDSI partition (Wheat/Rye and non-WR)
  irf_wr <- plot_teleco(WR, "teleco_PDSI_10",
                         "% Harvest response",
                         "IRF_ENSO_yield_WR_teleco.pdf")
  write.csv(irf_wr, file.path(out_data, "figA_irf_yield_WR_teleco.csv"), row.names = FALSE)

  # Summer precipitation (AMJJ) partition — used in the SI as an alternative
  # teleconnection criterion.
  if ("teleco_precip_summer_10" %in% names(yield_2023)) {
    plot_teleco(WR, "teleco_precip_summer_10",
                "% Harvest response\n(AMJJ Precip. partition)",
                "figA_irf_yield_WR_AMJJ_P.pdf")
  }
}


###############################################################################
# Regression table (appendix)
###############################################################################
save_yield_table <- function(yield_2023) {
  suppressPackageStartupMessages(library(modelsummary))
  # v2 yield panel, raw logyield outcome, pooled distributed-lag spec.
  # FE: VarLocationGrain[Year] (linear per-record trend) + i(decade).
  WR <- load_yields_GT3() |>
    filter(Grain %in% c("Wheat", "Rye"))
  pid_spec <- c("VarLocationGrain", "Year")

  est1 <- feols(logyield ~ l(nino34, 0:3) +
                  i(decade) | VarLocationGrain[Year],
                data = WR, panel.id = pid_spec)
  est2 <- feols(logyield ~ l(nino34, 0:3) +
                  i(decade) + PDSI + temp_summer + temp_winter +
                  precip_summer + precip_winter |
                  VarLocationGrain[Year],
                data = WR, panel.id = pid_spec)
  est3 <- feols(logyield ~ l(nino34, 0:3) +
                  i(decade) + ongoing_wars + log(1 + Deaths) |
                  VarLocationGrain[Year],
                data = WR, panel.id = pid_spec)
  est4 <- feols(logyield ~ l(nino34, 0:3) +
                  i(decade) + ongoing_wars + log(1 + Deaths) +
                  temp_summer + temp_winter + precip_summer + precip_winter |
                  VarLocationGrain[Year],
                data = WR, panel.id = pid_spec)

  dict_yield <- c(
    "nino34"         = "NINO3.4 T",
    "l(nino34, 0)"   = "NINO3.4 T",
    "l(nino34, 1)"   = "NINO3.4 T-1",
    "l(nino34, 2)"   = "NINO3.4 T-2",
    "l(nino34, 3)"   = "NINO3.4 T-3",
    "ongoing_wars"   = "Ongoing Wars",
    "log(1+Deaths)"  = "Log(1+Deaths)",
    "logyield"       = "Log Grain Harvest"
  )

  etable(
    est1, est2, est3, est4,
    keep         = "%nino34",
    vcov         = DK ~ Year,
    dict         = dict_yield,
    label        = "tab:FE_ENSO_yieldPDSI",
    title        = "Effect of a 1\\textdegree C Anomaly in the Nino 3.4 Index on Wheat and Rye Grain Harvests.",
    drop.section = c("fixef", "slopes"),
    fitstat      = c("n", "r2"),
    extralines   = list(
      "Decade FE"                = c("Yes","Yes","Yes","Yes"),
      "Record-Location-Grain FE" = c("Yes","Yes","Yes","Yes"),
      "Controls"                 = c("None","Climate","Conflict","Climate + Conflict")
    ),
    file    = file.path(tab_out, "FE_ENSO_yieldPDSI.tex"),
    replace = TRUE
  )
  message("Saved FE_ENSO_yieldPDSI.tex")
}


###############################################################################
# Appendix plots
###############################################################################
fig_out <- file.path(SCRIPT_DIR, "output", "figures", "appendix")
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

.loo_overlay_plot <- function(base_df, loo_df, group_col, y_label) {
  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_line(data = loo_df,
              aes(x = horizon, y = irf_mean, group = .data[[group_col]]),
              color = "cornflowerblue", linewidth = 0.9, linetype = "dashed") +
    geom_errorbar(data = base_df,
                  aes(x = horizon, ymin = irf_down, ymax = irf_up),
                  color = "black", width = 0.15, linewidth = 0.9) +
    geom_point(data = base_df, aes(x = horizon, y = irf_mean),
               color = "black", size = 2.5) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    labs(x = "Horizon (years)", y = y_label) +
    theme_classic(base_size = 18) +
    theme(panel.grid = element_blank())
}

plot_loo_WR <- function() {
  base_path <- file.path(out_data, "fig3C_irf_yield_WR.csv")
  loc_path  <- file.path(out_data, "figA_irf_subsample_WR.csv")
  y20_path  <- file.path(out_data, "figA_irf_20y_WR.csv")
  for (p in c(base_path, loc_path, y20_path)) {
    if (!file.exists(p)) { message(basename(p), " not found"); return() }
  }
  base <- read.csv(base_path) |> dplyr::filter(label == "All regions")
  loc  <- read.csv(loc_path)
  y20  <- read.csv(y20_path)
  y_lab <- "% Harvest response"

  ggsave(file.path(fig_out, "figA_irf_loo_loc_WR.pdf"),
         .loo_overlay_plot(base, loc, "subsample", y_lab),
         width = 6, height = 5, device = cairo_pdf)
  message("Saved figA_irf_loo_loc_WR.pdf")

  ggsave(file.path(fig_out, "figA_irf_loo_20y_WR.pdf"),
         .loo_overlay_plot(base, y20, "block", y_lab),
         width = 6, height = 5, device = cairo_pdf)
  message("Saved figA_irf_loo_20y_WR.pdf")
}


plot_yield63 <- function() {
  path <- file.path(out_data, "figA_irf_yield63.csv")
  if (!file.exists(path)) { message("figA_irf_yield63.csv not found"); return() }
  df <- read.csv(path)
  n_ctrl  <- length(unique(df$controls))
  offsets <- seq(-0.1, 0.1, length.out = n_ctrl)
  df      <- df |>
    dplyr::group_by(controls) |>
    dplyr::mutate(offset = offsets[dplyr::cur_group_id()]) |>
    dplyr::ungroup()
  p <- ggplot(df, aes(x = horizon + offset, y = irf_mean, color = controls)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_errorbar(aes(ymin = irf_down, ymax = irf_up), width = 0.1, linewidth = 0.9) +
    geom_point(size = 2.5) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    labs(x = "Horizon (years)",
         y = "% Harvest response",
         color = "Controls") +
    theme_classic(base_size = 18) +
    theme(panel.grid = element_blank(),
          legend.position  = "bottom",
          legend.box       = "horizontal",
          legend.text      = element_text(size = 11),
          legend.title     = element_text(size = 11),
          plot.margin      = margin(5, 5, 5, 5)) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE,
                                title.position = "top"))
  ggsave(file.path(fig_out, "figA_irf_yield63.pdf"), p,
         width = 8, height = 6, device = cairo_pdf)
  message("Saved figA_irf_yield63.pdf")
}

plot_bootstrap_WR <- function() {
  null_path <- file.path(out_data, "figA_irf_WR_bootstrap_null.csv")
  true_path <- file.path(out_data, "figA_irf_WR_bootstrap_true.csv")
  for (p in c(null_path, true_path)) {
    if (!file.exists(p)) { message(basename(p), " not found"); return() }
  }
  null_df  <- read.csv(null_path)
  true_df  <- read.csv(true_path)
  null_env <- null_df |>
    dplyr::group_by(horizon) |>
    dplyr::summarise(lo = quantile(irf_mean, 0.025),
                     hi = quantile(irf_mean, 0.975), .groups = "drop")
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_ribbon(data = null_env, aes(x = horizon, ymin = lo, ymax = hi),
                fill = "cornflowerblue", alpha = 0.25) +
    geom_line(data = true_df, aes(x = horizon, y = irf_mean),
              color = "black", linewidth = 1.1) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                       limits = c(-0.06, NA),
                       breaks = c(-0.05, 0, 0.05)) +
    labs(x = "Horizon (years)",
         y = "% Harvest response",
         caption = "Blue band: 95% permutation null interval. Black line: true IRF.") +
    theme_classic(base_size = 18) +
    theme(panel.grid = element_blank())
  ggsave(file.path(fig_out, "figA_irf_bootstrap_WR.pdf"), p,
         width = 6, height = 5, device = cairo_pdf)
  message("Saved figA_irf_bootstrap_WR.pdf")
}


plot_se_robust_WR <- function() {
  path <- file.path(out_data, "figA_irf_WR_serobust.csv")
  if (!file.exists(path)) { message("figA_irf_WR_serobust.csv not found"); return() }
  df <- read.csv(path)
  p <- ggplot(df, aes(x = horizon, y = irf_mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_errorbar(aes(ymin = irf_down, ymax = irf_up),
                  color = "black", width = 0.15, linewidth = 0.9) +
    geom_point(color = "black", size = 2.5) +
    facet_wrap(~ se_type, ncol = 4) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    labs(x = "Horizon (years)",
         y = "% Harvest response") +
    theme_classic(base_size = 18) +
    theme(strip.background = element_blank(),
          strip.text       = element_text(face = "bold"),
          panel.grid       = element_blank())
  ggsave(file.path(fig_out, "figA_irf_serobust_WR.pdf"), p,
         width = 10, height = 4, device = cairo_pdf)
  message("Saved figA_irf_serobust_WR.pdf")
}


###############################################################################
# Appendix: CE vs. rest robustness (geographic partition instead of teleconnection)
###############################################################################
export_irf_geo_partition <- function(yield_2023) {
  # Central Europe = Germany, Switzerland, Hungary, Austria
  # (mirrors the famine-onset "Central Europe" region)
  ce_countries <- c("Germany", "Switzerland", "Hungary", "Austria")

  # Interaction-based partition: fit one LP on the joint sample with
  # i(CE, nino34); pull CE = 1 / CE = 0 slopes from the same fit. The
  # pooled "All regions" line stays as a separate pooled LP for reference.
  run_geo_irf_inter <- function(df, endog, hor_max) {
    df <- df |> mutate(CE = as.integer(Country %in% ce_countries))
    if (nrow(df) < 50) return(NULL)

    out_all <- extract_irf(run_lp_GT3(df, endog, hor = hor_max),
                            "All regions",
                            dplyr::n_distinct(df$VarLocationGrain))

    res_int <- run_lp_GT3_inter(df, endog, hor = hor_max,
                                 interact_var = "CE")
    n_ce   <- dplyr::n_distinct(df$VarLocationGrain[df$CE == 1])
    n_rest <- dplyr::n_distinct(df$VarLocationGrain[df$CE == 0])
    out_ce   <- res_int |> dplyr::filter(category == "1") |>
      transmute(label = "Central Europe", n_loc = n_ce,
                horizon, irf_mean, irf_up, irf_down)
    out_rest <- res_int |> dplyr::filter(category == "0") |>
      transmute(label = "Rest of Europe", n_loc = n_rest,
                horizon, irf_mean, irf_up, irf_down)
    bind_rows(out_all, out_ce, out_rest)
  }

  yields_GT3 <- load_yields_GT3()
  WR   <- yields_GT3 |>
    filter(Grain %in% c("Wheat", "Rye"))
  noWR <- yields_GT3 |>
    filter(!Grain %in% c("Wheat", "Rye"))

  irf_wr   <- run_geo_irf_inter(WR,   "logyield", 3)
  irf_nowr <- run_geo_irf_inter(noWR, "logyield", 3)

  make_geo_panel <- function(df, title_lab) {
    groups  <- unique(df$label)
    n_g     <- length(groups)
    offsets <- setNames(seq(-0.12, 0.12, length.out = n_g), groups)
    df      <- df |> dplyr::mutate(offset = offsets[label])
    # build color mapping: match labels by prefix
    clr <- c("All regions"    = "black",
             "Central Europe" = "firebrick",
             "Rest of Europe" = "cornflowerblue")
    clr <- clr[names(clr) %in% groups]
    ggplot(df, aes(x = horizon + offset, y = irf_mean, color = label)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
      geom_errorbar(aes(ymin = irf_down, ymax = irf_up), width = 0, linewidth = 0.9) +
      geom_point(size = 2.5) +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
      scale_color_manual(values = clr) +
      labs(x = "Horizon (years)", y = title_lab, color = NULL) +
      theme_classic(base_size = 16) +
      theme(panel.grid = element_blank(), legend.position = "bottom")
  }

  # shared y-axis limits
  y_range <- range(c(irf_wr$irf_down, irf_wr$irf_up,
                     irf_nowr$irf_down, irf_nowr$irf_up), na.rm = TRUE)

  # extract legend from a standalone plot before suppressing it in panels
  legend <- cowplot::get_legend(
    make_geo_panel(irf_wr, "") + theme(legend.position = "bottom")
  )

  p_wr   <- make_geo_panel(irf_wr,   "% Harvest response") +
    coord_cartesian(ylim = y_range) +
    labs(title = "a  Wheat & Rye") +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))
  p_nowr <- make_geo_panel(irf_nowr, "% Harvest response") +
    coord_cartesian(ylim = y_range) +
    labs(title = "b  Other grains") +
    theme(legend.position = "none", axis.title.y = element_blank(),
          plot.title = element_text(face = "bold"))

  combined <- cowplot::plot_grid(p_wr, p_nowr, nrow = 1, align = "h",
                                  rel_widths = c(1, 1))
  final    <- cowplot::plot_grid(combined, legend, ncol = 1,
                                  rel_heights = c(1, 0.12))

  # Extended Data figure – save under figED_ prefix into the ED folder.
  ed_dir <- file.path(SCRIPT_DIR, "output", "figures", "extended data")
  dir.create(ed_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(ed_dir, "figED_yield_geopartition.pdf"), final,
         width = 10, height = 5, device = cairo_pdf)
  message("Saved figED_yield_geopartition.pdf")
}


###############################################################################
# Run
###############################################################################
if (!interactive()) {
  yield_2023 <- load_yields()
  export_fig3C(yield_2023)
  export_fig3D(yield_2023)
  save_yield_table(yield_2023)
  export_yield63_irf()
  export_se_robust_WR(yield_2023)
  export_loo_WR(yield_2023)
  export_20y_WR(yield_2023)
  export_bootstrap_WR(yield_2023, n_iter = 500)
  plot_loo_WR()
  plot_yield63()
  plot_bootstrap_WR()
  plot_se_robust_WR()
  export_yield_teleco_plots(yield_2023)
  export_irf_geo_partition(yield_2023)
  message("Fig 3 data exports complete.")
} else {
  message("Call functions interactively or source this script.")
}
