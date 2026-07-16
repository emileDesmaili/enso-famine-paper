###############################################################################
# Figure 5 – ENSO shocks broadly raised staple food prices
# This script exports tidy CSVs for assemble_figures.py
#
# CSV exports → analysis/output/data/
#   fig5A_irf_grain_price.csv   – Grain price IRF (teleconnected vs weakly)
#   fig5B_irf_fish_price.csv    – Fish price IRF by species (All/Herring/Cod)
#
# Appendix CSVs:
#   figA_irf_FSV.csv
#   figA_irf_subsample_price.csv
#   figA_irf_subsample_20y_price.csv
#   figA_irf_subsample_fishprice.csv
#   figA_irf_subsample_20y_fishprice.csv
#   figA_irf_price_bootstrap_null.csv / _true.csv
#   figA_irf_fishprice_bootstrap_null.csv / _true.csv
#   figA_irf_allen.csv
#
# Appendix PDFs:
#   figA_irf_allen.pdf  – Allen wheat price IRF under 4 control specs
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(scales)
  library(haven)
  library(purrr)
  library(furrr)
  library(future)
  library(fixest)
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
for (d in c(out_data, fig_out)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── load emileRegs helpers ────────────────────────────────────────────────────
source(file.path(ROOT, "emileRegs.R"))

# ── helpers ───────────────────────────────────────────────────────────────────
# Includes 3 lags of the shock in the controls (NINO3.4 persistence).
run_lp_price <- function(df, endog, shock = "nino34",
                          controls = "i(Decade)",
                          entity   = "Location",
                          hor      = 10) {
  lag_ctrl <- paste0("l(", shock, ", 1:3)")
  ctrl_str <- if (nzchar(controls)) paste(controls, lag_ctrl, sep = " + ")
              else lag_ctrl
  lp_panel(
    data         = df,
    outcome      = endog,
    main_var     = shock,
    controls     = ctrl_str,
    horizon      = hor,
    fe           = paste0(entity, "[Year]"),   # entity FE + entity-specific linear trend
    panel_id     = c(entity, "Year"),
    vcov_formula = DK ~ Year
  )
}

# Build a controls string from a character vector, replacing bare "Decade"
# with the fixest factor notation i(Decade)
.build_ctrl <- function(ctrl_vec) {
  parts <- ifelse(ctrl_vec == "Decade", "i(Decade)", ctrl_vec)
  paste(parts, collapse = " + ")
}

load_prices <- function() {
  read.csv(file.path(DATA_PROC, "price_2023_enso.csv")) |>
    mutate(Decade = floor(Year / 10) * 10) |>
    dplyr::select(Location, Year, Decade, everything()) |>
    drop_na() |> distinct() |> arrange(Year)
}

load_fishprice <- function() {
  read.csv(file.path(DATA_PROC, "fishprice_enso.csv")) |>
    mutate(Decade = floor(Year / 10) * 10) |>
    dplyr::select(LocationSpecies, Year, Decade, everything())
}

load_fishcatch <- function() {
  read.csv(file.path(DATA_PROC, "fishcatch_enso.csv")) |>
    mutate(Decade = floor(Year / 10) * 10) |>
    dplyr::select(TypeLocation, Year, Decade, everything())
}


###############################################################################
# Fig 5A – Grain price IRF
# Teleconnected / Weakly-affected partition uses an interacted LP
# (i(teleco_PDSI_10, nino34)) so both slopes and CIs come from a single fit.
# "All regions" and "All w. controls" remain pooled fits.
###############################################################################
export_fig5A <- function(data) {
  hor <- 10

  # Interacted LP for the teleconnection partition
  df_inter <- data |>
    mutate(Decade = as.factor(Decade),
           teleco_PDSI_10 = as.integer(teleco_PDSI_10))
  n_tele <- dplyr::n_distinct(filter(data, teleco_PDSI_10 == 1)$Location)
  n_weak <- dplyr::n_distinct(filter(data, teleco_PDSI_10 == 0)$Location)

  irf_inter <- lp_panel_inter(
    data         = df_inter,
    outcome      = "logprice",
    main_var     = "nino34",
    interact_var = "teleco_PDSI_10",
    controls     = "i(Decade) + l(nino34, 1:3)",
    horizon      = hor,
    fe           = "Location[Year]",
    panel_id     = c("Location", "Year"),
    vcov_formula = DK ~ Year
  ) |>
    mutate(label = case_when(
      category == "1" ~ "Teleconnected",
      category == "0" ~ "Weakly affected",
      TRUE ~ as.character(category)
    ),
    n_loc = if_else(label == "Teleconnected", n_tele, n_weak)) |>
    transmute(label, n_loc, horizon, irf_mean, irf_up, irf_down)

  # Pooled LP for full sample and full sample with controls
  pooled_specs <- list(
    "All regions"     = list(filter = NULL, ctrl = "Decade"),
    "All w. controls" = list(filter = NULL,
                              ctrl = c("Decade", "JSL", "NAO_cal",
                                       "ongoing_wars", "Deaths"))
  )
  irf_pooled <- bind_rows(lapply(names(pooled_specs), function(label) {
    spec  <- pooled_specs[[label]]
    n_loc <- dplyr::n_distinct(data$Location)
    run_lp_price(data, "logprice",
                 controls = .build_ctrl(spec$ctrl),
                 hor      = hor) |>
      transmute(label = label, n_loc = n_loc, horizon, irf_mean, irf_up, irf_down)
  }))

  bind_rows(irf_pooled, irf_inter) |>
    write.csv(file.path(out_data, "fig5A_irf_grain_price.csv"), row.names = FALSE)
  message("Exported fig5A_irf_grain_price.csv")
}


###############################################################################
# Fig 5B – Fish price IRF by species
###############################################################################
# Sum-to-zero (deviation) coding so the bare `shock` coefficient is the
# unweighted (Cod + Herring) / 2 average, with a native SE. The Cod and Herring
# species-specific slopes are pulled from `lp_panel_inter` (mathematically
# identical to shock ± 0.5 * shock:cod_dev in this same fit).
run_lp_fish_devcoded <- function(fishprice, shock = "nino34", hor = 10,
                                  controls = "i(Decade)") {
  d <- fishprice |>
    mutate(Decade  = as.factor(floor(Year / 10) * 10),
           cod_dev = ifelse(Species == "Cod", 0.5, -0.5))
  lag_ctrl <- paste0("l(", shock, ", 1:3)")
  ctrl_str <- if (nzchar(controls)) paste(controls, lag_ctrl, sep = " + ")
              else lag_ctrl
  bind_rows(lapply(0:hor, function(h) {
    fml <- as.formula(paste0(
      "f(logprice,", h, ") - l(logprice,1) ~ ",
      shock, " + ", shock, ":cod_dev + ", ctrl_str,
      " | LocationSpecies[Year]"))
    m  <- tryCatch(
      feols(fml, data = d,
            panel.id = c("LocationSpecies", "Year"),
            vcov     = DK ~ Year),
      error = function(e) NULL)
    if (is.null(m) || !(shock %in% rownames(vcov(m)))) {
      return(data.frame(horizon = as.integer(h),
                        irf_mean = NA_real_,
                        irf_up   = NA_real_,
                        irf_down = NA_real_))
    }
    V  <- vcov(m)
    b  <- as.numeric(unname(coef(m)[shock]))
    se <- as.numeric(sqrt(V[shock, shock]))
    data.frame(horizon  = as.integer(h),
               irf_mean = b,
               irf_up   = b + 1.96 * se,
               irf_down = b - 1.96 * se)
  }))
}

export_fig5B <- function(fishprice) {
  # Main IRF: deviation-coded "All" = unweighted (Cod + Herring) / 2 with
  # native SE from a single fit that also estimates species heterogeneity.
  irf_main <- run_lp_fish_devcoded(fishprice, hor = 10) |>
    transmute(species = "All", horizon, irf_mean, irf_up, irf_down)

  # Interaction IRFs: Cod and Herring via lp_panel_inter
  irf_inter <- lp_panel_inter(
    data         = fishprice |> mutate(Decade = as.factor(Decade)),
    outcome      = "logprice",
    main_var     = "nino34",
    interact_var = "Species",
    controls     = "i(Decade) + l(nino34, 1:3)",
    horizon      = 10,
    fe           = "LocationSpecies[Year]",
    panel_id     = c("LocationSpecies", "Year"),
    vcov_formula = DK ~ Year
  ) |>
    transmute(species = category, horizon, irf_mean, irf_up, irf_down)

  bind_rows(irf_main, irf_inter) |>
    write.csv(file.path(out_data, "fig5B_irf_fish_price.csv"), row.names = FALSE)
  message("Exported fig5B_irf_fish_price.csv")
}


###############################################################################
# Appendix: SE robustness – grain & fish price IRFs under 4 vcov types
###############################################################################
# SE types:
#   DK          – Driscoll-Kraay (baseline)
#   2-way       – two-way cluster (entity × Year)
#   Conley      – spatial HAC via vcov_conley
# Combined Spatial HAC = Conley + NW − Hetero is not guaranteed to be PSD;
# when a diagonal entry goes non-positive, fall back to the larger of the
# Conley and NW variances for that entry (still honest, avoids sqrt(NaN)).
.clip_spatial_hac <- function(v_spat, v_conley, v_nw) {
  d       <- diag(v_spat)
  bad     <- which(!is.finite(d) | d <= 0)
  if (length(bad) > 0) {
    d[bad] <- pmax(diag(v_conley)[bad], diag(v_nw)[bad])
    diag(v_spat) <- d
  }
  v_spat
}
#   Spatial HAC – Conley + NW - Hetero (combined spatial+serial+correct)
run_lp_se_robust_price <- function(df, endog, entity,
                                    lat_col = "Latitude", lon_col = "Longitude",
                                    hor = 10) {
  pid <- c(entity, "Year")
  bind_rows(lapply(0:hor, function(h) {
    fml <- as.formula(paste0(
      "f(", endog, ",", h, ") - l(", endog, ",1) ~ nino34 + l(nino34, 1:3) ",
      "+ i(Decade) | ", entity, "[Year]"
    ))
    mod <- feols(fml, data = df, panel.id = pid, vcov = DK ~ Year)

    v_conley <- tryCatch(
      vcov(mod, vcov = vcov_conley(lat = lat_col, lon = lon_col,
                                   cutoff = 300, distance = "triangular")),
      error = function(e) NULL
    )
    v_nw   <- vcov(mod, vcov = NW ~ Year)
    v_het  <- vcov(mod, vcov = "hetero")
    v_spat <- if (!is.null(v_conley))
                .clip_spatial_hac(v_conley + v_nw - v_het, v_conley, v_nw)
              else NULL

    se_list <- list(
      "DK"          = vcov(mod, vcov = DK ~ Year),
      "2-way"       = vcov(mod, vcov = as.formula(paste0("cluster ~ ", entity, " + Year"))),
      "Conley"      = v_conley,
      "Spatial HAC" = v_spat
    )

    bind_rows(lapply(names(se_list), function(nm) {
      V <- se_list[[nm]]
      if (is.null(V)) return(NULL)
      b  <- coef(mod)["nino34"]
      se <- sqrt(V["nino34", "nino34"])
      data.frame(se_type   = nm, horizon = h,
                 irf_mean  = b,
                 irf_up95  = b + 1.960 * se,
                 irf_down95 = b - 1.960 * se,
                 irf_up90  = b + 1.645 * se,
                 irf_down90 = b - 1.645 * se)
    }))
  }))
}

export_se_robust_price <- function(data) {
  run_lp_se_robust_price(data, "logprice", "Location", hor = 10) |>
    write.csv(file.path(out_data, "figA_irf_price_serobust.csv"), row.names = FALSE)
  message("Exported figA_irf_price_serobust.csv")
}

# Deviation-coded SE-robustness: bare `nino34` = (Cod+Herring)/2 average,
# four vcov variants applied to the same dev-coded LP fit.
run_lp_se_robust_fish_devcoded <- function(fishprice,
                                            lat_col = "Latitude",
                                            lon_col = "Longitude",
                                            hor = 10) {
  d <- fishprice |>
    mutate(Decade  = as.factor(floor(Year / 10) * 10),
           cod_dev = ifelse(Species == "Cod", 0.5, -0.5))
  pid <- c("LocationSpecies", "Year")
  bind_rows(lapply(0:hor, function(h) {
    fml <- as.formula(paste0(
      "f(logprice,", h, ") - l(logprice,1) ~ ",
      "nino34 + nino34:cod_dev + l(nino34, 1:3) + i(Decade) | LocationSpecies[Year]"))
    mod <- feols(fml, data = d, panel.id = pid, vcov = DK ~ Year)

    v_conley <- tryCatch(
      vcov(mod, vcov = vcov_conley(lat = lat_col, lon = lon_col,
                                   cutoff = 300, distance = "triangular")),
      error = function(e) NULL
    )
    v_nw   <- vcov(mod, vcov = NW ~ Year)
    v_het  <- vcov(mod, vcov = "hetero")
    v_spat <- if (!is.null(v_conley))
                .clip_spatial_hac(v_conley + v_nw - v_het, v_conley, v_nw)
              else NULL

    se_list <- list(
      "DK"          = vcov(mod, vcov = DK ~ Year),
      "2-way"       = vcov(mod, vcov = ~ LocationSpecies + Year),
      "Conley"      = v_conley,
      "Spatial HAC" = v_spat
    )

    bind_rows(lapply(names(se_list), function(nm) {
      V <- se_list[[nm]]
      if (is.null(V)) return(NULL)
      b  <- coef(mod)["nino34"]
      se <- sqrt(V["nino34", "nino34"])
      data.frame(se_type    = nm, horizon = h,
                 irf_mean   = b,
                 irf_up95   = b + 1.960 * se,
                 irf_down95 = b - 1.960 * se,
                 irf_up90   = b + 1.645 * se,
                 irf_down90 = b - 1.645 * se)
    }))
  }))
}

export_se_robust_fish <- function(fishprice) {
  run_lp_se_robust_fish_devcoded(fishprice, hor = 10) |>
    write.csv(file.path(out_data, "figA_irf_fishprice_serobust.csv"), row.names = FALSE)
  message("Exported figA_irf_fishprice_serobust.csv")
}


###############################################################################
# Appendix: Allen wheat price IRF
###############################################################################
export_allen_irf <- function() {
  allen_path <- file.path(ROOT, "processed data", "wheatprices_dataset.dta")
  if (!file.exists(allen_path)) {
    message("wheatprices_dataset.dta not found – skipping Allen IRF")
    return(invisible(NULL))
  }

  # Year-level ENSO + controls from famine_region_data (averaged across regions)
  region_data <- read.csv(file.path(DATA_PROC, "famine_region_data.csv")) |>
    dplyr::group_by(Year) |>
    dplyr::summarise(
      nino34          = mean(nino34,          na.rm = TRUE),
      PDSI_europe     = mean(PDSI_europe,     na.rm = TRUE),
      ongoing_wars    = mean(ongoing_wars,    na.rm = TRUE),
      Deaths          = mean(Deaths,          na.rm = TRUE),
      .groups = "drop"
    )

  data_allen <- haven::read_dta(allen_path) |>
    dplyr::rename(Year = year) |>
    dplyr::left_join(region_data, by = "Year") |>
    mutate(decade = floor(Year / 10) * 10) |>
    drop_na(citynum, Year, logwheatprice, nino34) |>
    distinct()

  ctrl_dict <- list(
    "Base"               = "decade",
    "Climate"            = c("decade", "PDSI_europe", "growseasontemp",
                             "meantemp", "growseasonrain"),
    "Conflict"           = c("decade", "ongoing_wars", "Deaths"),
    "Climate + Conflict" = c("decade", "PDSI_europe", "growseasontemp",
                             "meantemp", "growseasonrain", "ongoing_wars", "Deaths")
  )

  results <- bind_rows(lapply(names(ctrl_dict), function(nm) {
    avail <- names(data_allen)
    ctrl_vec <- ctrl_dict[[nm]]
    ctrl_vec <- ctrl_vec[ctrl_vec %in% c("decade", avail)]
    ctrl_str <- paste(
      c(ifelse(ctrl_vec == "decade", "i(decade)", ctrl_vec),
        "l(nino34, 1:3)"),
      collapse = " + "
    )
    tryCatch(
      lp_panel(
        data         = data_allen,
        outcome      = "logwheatprice",
        main_var     = "nino34",
        controls     = ctrl_str,
        horizon      = 10,
        fe           = "citynum",
        panel_id     = c("citynum", "Year"),
        vcov_formula = DK ~ Year
      ) |> transmute(controls = nm, horizon, irf_mean, irf_up, irf_down),
      error = function(e) {
        message("Allen IRF failed for controls '", nm, "': ", e$message)
        NULL
      }
    )
  }))

  if (is.null(results) || nrow(results) == 0) {
    message("Allen IRF: no results to export")
    return(invisible(NULL))
  }

  write.csv(results, file.path(out_data, "figA_irf_allen.csv"), row.names = FALSE)
  message("Exported figA_irf_allen.csv")
}

plot_allen_irf <- function() {
  path <- file.path(out_data, "figA_irf_allen.csv")
  if (!file.exists(path)) { message("figA_irf_allen.csv not found"); return() }
  .plot_ctrl_irf(read.csv(path), "Price variation", "figA_irf_allen.pdf")
}


###############################################################################
# Appendix: Federico (FSV) wheat price IRF
###############################################################################
export_FSV_irf <- function() {
  data_federico <- read.csv(file.path(DATA_PROC, "federico_prices_enso.csv")) |>
    mutate(Decade = floor(Year / 10) * 10) |>
    drop_na() |> distinct() |> arrange(Year)

  ctrl_dict <- list(
    "Base"               = "Decade",
    "Conflict"           = c("Decade","ongoing_wars","Deaths"),
    "Climate"            = c("Decade","JSL","NAO_cal"),
    "Climate + Conflict" = c("Decade","JSL","NAO_cal","ongoing_wars","Deaths")
  )

  bind_rows(lapply(names(ctrl_dict), function(nm) {
    run_lp_price(data_federico, "logprice",
                 controls = .build_ctrl(ctrl_dict[[nm]]),
                 hor      = 10) |>
      transmute(controls = nm, horizon, irf_mean, irf_up, irf_down)
  })) |>
    write.csv(file.path(out_data, "figA_irf_FSV.csv"), row.names = FALSE)
  message("Exported figA_irf_FSV.csv")
}


###############################################################################
# Appendix: LOO subsampling (grain prices)
###############################################################################
export_loo_prices <- function(data) {
  locations <- unique(data$Location)
  combos    <- combn(locations, length(locations) - 1, simplify = FALSE)

  bind_rows(lapply(seq_along(combos), function(i) {
    left_out <- setdiff(locations, combos[[i]])
    run_lp_price(data |> filter(Location %in% combos[[i]]),
                 "logprice", controls = "i(Decade)", hor = 10) |>
      transmute(subsample = paste("No", left_out), horizon, irf_mean, irf_up, irf_down)
  })) |>
    write.csv(file.path(out_data, "figA_irf_subsample_price.csv"), row.names = FALSE)
  message("Exported figA_irf_subsample_price.csv")
}


###############################################################################
# Appendix: 20-year block LOO (grain prices)
###############################################################################
export_20y_prices <- function(data) {
  years        <- sort(unique(data$Year))
  block_starts <- seq(min(years), max(years) - 20, by = 20)

  bind_rows(lapply(block_starts, function(s) {
    run_lp_price(data |> filter(!Year %in% s:(s + 19)),
                 "logprice", controls = "i(Decade)", hor = 10) |>
      transmute(block = paste0("Drop ", s, "-", s + 19), horizon, irf_mean, irf_up, irf_down)
  })) |>
    write.csv(file.path(out_data, "figA_irf_subsample_20y_price.csv"), row.names = FALSE)
  message("Exported figA_irf_subsample_20y_price.csv")
}


###############################################################################
# Appendix: LOO subsampling (fish prices)
# Deviation-coded LP: bare `nino34` = unweighted (Cod+Herring)/2.
# Dropping all-Cod or all-Herring subsamples falls back to a pooled fit
# (single species => cod_dev is constant, the interaction is collinear).
###############################################################################
export_loo_fish <- function(fishprice) {
  locations <- unique(fishprice$LocationSpecies)
  combos    <- combn(locations, length(locations) - 1, simplify = FALSE)

  bind_rows(lapply(seq_along(combos), function(i) {
    left_out <- setdiff(locations, combos[[i]])
    sub <- fishprice |> filter(LocationSpecies %in% combos[[i]])
    irf <- if (length(unique(sub$Species)) >= 2) {
      run_lp_fish_devcoded(sub, controls = "i(Decade)", hor = 10)
    } else {
      run_lp_price(sub, "logprice", controls = "i(Decade)",
                   entity = "LocationSpecies", hor = 10) |>
        transmute(horizon, irf_mean, irf_up, irf_down)
    }
    irf |> transmute(subsample = paste("No", left_out),
                     horizon, irf_mean, irf_up, irf_down)
  })) |>
    write.csv(file.path(out_data, "figA_irf_subsample_fishprice.csv"), row.names = FALSE)
  message("Exported figA_irf_subsample_fishprice.csv")
}


###############################################################################
# Appendix: 20-year block LOO (fish prices)
# Deviation-coded LP: bare `nino34` = unweighted (Cod+Herring)/2.
###############################################################################
export_20y_fish <- function(fishprice) {
  years        <- sort(unique(fishprice$Year))
  block_starts <- seq(min(years), max(years) - 20, by = 20)

  bind_rows(lapply(block_starts, function(s) {
    sub <- fishprice |> filter(!Year %in% s:(s + 19))
    irf <- if (length(unique(sub$Species)) >= 2) {
      run_lp_fish_devcoded(sub, controls = "i(Decade)", hor = 10)
    } else {
      run_lp_price(sub, "logprice", controls = "i(Decade)",
                   entity = "LocationSpecies", hor = 10) |>
        transmute(horizon, irf_mean, irf_up, irf_down)
    }
    irf |> transmute(block = paste0("Drop ", s, "-", s + 19),
                     horizon, irf_mean, irf_up, irf_down)
  })) |>
    write.csv(file.path(out_data, "figA_irf_subsample_20y_fishprice.csv"), row.names = FALSE)
  message("Exported figA_irf_subsample_20y_fishprice.csv")
}


###############################################################################
# Appendix: fish catch IRF – NAO and NINO3.4 shocks
###############################################################################
export_fishcatch_irf <- function(fishcatch) {
  shocks <- c("nino34" = "NINO3.4", "NAO_cal" = "NAO")
  bind_rows(lapply(names(shocks), function(shock) {
    run_lp_price(fishcatch, "logcatch",
                 shock    = shock,
                 controls = "i(Decade)",
                 entity   = "TypeLocation",
                 hor      = 10) |>
      transmute(shock = shocks[[shock]], horizon, irf_mean, irf_up, irf_down)
  })) |>
    write.csv(file.path(out_data, "figA_irf_catches.csv"), row.names = FALSE)
  message("Exported figA_irf_catches.csv")
}

plot_fishcatch_irf <- function() {
  path <- file.path(out_data, "figA_irf_catches.csv")
  if (!file.exists(path)) { message("figA_irf_catches.csv not found"); return() }
  df <- read.csv(path)
  shock_colors <- c("NINO3.4" = "#1b9e77", "NAO" = "#d95f02")
  shock_lty    <- c("NINO3.4" = "dashed",  "NAO" = "solid")
  p <- ggplot(df, aes(x = horizon, color = shock, fill = shock, linetype = shock)) +
    geom_ribbon(aes(ymin = irf_down, ymax = irf_up), alpha = 0.2, color = NA) +
    geom_line(aes(y = irf_mean), linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    scale_x_continuous(breaks = 0:10) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    scale_color_manual(values = shock_colors) +
    scale_fill_manual(values = shock_colors) +
    scale_linetype_manual(values = shock_lty) +
    labs(x = "Horizon (Years)",
         y = "% Fish catch response",
         color = "Shock", fill = "Shock", linetype = "Shock") +
    theme_classic(base_size = 18) +
    theme(legend.position = "top",
          panel.grid = element_blank())
  ggsave(file.path(fig_out, "figA_irf_catches.pdf"), p,
         width = 6, height = 5, device = cairo_pdf)
  message("Saved figA_irf_catches.pdf")
}


###############################################################################
# Appendix: bootstrap permutation
###############################################################################
export_bootstrap <- function(data, fishprice, n_iter = 500) {
  set.seed(42)
  plan(multisession)

  # Pooled bootstrap for grain prices
  run_boot_pooled <- function(df, endog, entity, tag) {
    true_res <- run_lp_price(df, endog, controls = "i(Decade)",
                              entity = entity, hor = 10)
    true_res |>
      transmute(horizon, irf_mean) |>
      write.csv(file.path(out_data, paste0(tag, "_true.csv")), row.names = FALSE)

    future_map_dfr(seq_len(n_iter), function(i) {
      sh <- df |> mutate(nino34 = sample(nino34))
      tryCatch({
        run_lp_price(sh, endog, controls = "i(Decade)",
                     entity = entity, hor = 10) |>
          transmute(iter = i, horizon, irf_mean)
      }, error = function(e) NULL)
    }, .options = furrr_options(seed = TRUE)) |>
      write.csv(file.path(out_data, paste0(tag, "_null.csv")), row.names = FALSE)

    message("Exported bootstrap CSVs for ", tag)
  }

  # Dev-coded bootstrap for fish prices: bare `nino34` = (Cod+Herring)/2
  run_boot_fish_dev <- function(df, tag) {
    true_res <- run_lp_fish_devcoded(df, controls = "i(Decade)", hor = 10)
    true_res |>
      transmute(horizon, irf_mean) |>
      write.csv(file.path(out_data, paste0(tag, "_true.csv")), row.names = FALSE)

    future_map_dfr(seq_len(n_iter), function(i) {
      sh <- df |> mutate(nino34 = sample(nino34))
      tryCatch({
        run_lp_fish_devcoded(sh, controls = "i(Decade)", hor = 10) |>
          transmute(iter = i, horizon, irf_mean)
      }, error = function(e) NULL)
    }, .options = furrr_options(seed = TRUE)) |>
      write.csv(file.path(out_data, paste0(tag, "_null.csv")), row.names = FALSE)

    message("Exported bootstrap CSVs for ", tag)
  }

  run_boot_pooled(data, "logprice", "Location", "figA_irf_price_bootstrap")
  run_boot_fish_dev(fishprice,                 "figA_irf_fishprice_bootstrap")
}


###############################################################################
# Appendix plots
###############################################################################
fig_out <- file.path(SCRIPT_DIR, "output", "figures", "appendix")
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

# Shared palette for 4-control IRF plots (Allen & FSV).
# Main "Base" spec in black with cornflowerblue ribbon (matches Fig 5 style);
# alternates overlaid as coloured dashed lines without ribbons.
.CTRL_COLORS <- c(
  "Base"               = "black",
  "Climate"            = "#D55E00",   # vermillion
  "Conflict"           = "#0072B2",   # blue
  "Climate + Conflict" = "#009E73"    # bluish green
)
.CTRL_LTY <- c(
  "Base"               = "solid",
  "Climate"            = "dashed",
  "Conflict"           = "dotdash",
  "Climate + Conflict" = "longdash"
)
.CTRL_RIBBON_COL <- "cornflowerblue"

.plot_ctrl_irf <- function(df, y_label, out_file) {
  df <- df |>
    dplyr::mutate(controls = factor(controls,
                                    levels = names(.CTRL_COLORS)))
  base_df <- dplyr::filter(df, controls == "Base")

  p <- ggplot(df, aes(x = horizon, y = irf_mean,
                      colour = controls, linetype = controls)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "gray40", linewidth = 0.6) +
    geom_ribbon(data = base_df,
                aes(x = horizon, ymin = irf_down, ymax = irf_up),
                fill = .CTRL_RIBBON_COL, colour = NA, alpha = 0.20,
                inherit.aes = FALSE, show.legend = FALSE) +
    geom_line(data = dplyr::filter(df, controls != "Base"), linewidth = 1.0) +
    geom_line(data = base_df, linewidth = 1.4) +
    scale_x_continuous(breaks = 0:10) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    scale_colour_manual(values = .CTRL_COLORS) +
    scale_linetype_manual(values = .CTRL_LTY) +
    labs(x = "Horizon (years)", y = y_label,
         colour = "Controls", linetype = "Controls") +
    theme_classic(base_size = 16) +
    theme(legend.position  = "bottom",
          legend.text      = element_text(size = 11),
          legend.title     = element_text(size = 11),
          panel.grid       = element_blank()) +
    guides(colour   = guide_legend(nrow = 1, byrow = TRUE,
                                   override.aes = list(linewidth = 1.2)),
           linetype = guide_legend(nrow = 1, byrow = TRUE))

  ggsave(file.path(fig_out, out_file), p,
         width = 7, height = 5, device = cairo_pdf)
  message("Saved ", out_file)
}

.loo_overlay_plot <- function(base_df, loo_df, group_col, y_label) {
  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_ribbon(data = base_df,
                aes(x = horizon, ymin = irf_down, ymax = irf_up),
                fill = "cornflowerblue", alpha = 0.25) +
    geom_line(data = loo_df,
              aes(x = horizon, y = irf_mean, group = .data[[group_col]]),
              color = "cornflowerblue", linewidth = 1.1, linetype = "dashed") +
    geom_line(data = base_df, aes(x = horizon, y = irf_mean),
              color = "black", linewidth = 1.1) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    labs(x = "Horizon (years)", y = y_label) +
    theme_classic(base_size = 18) +
    theme(panel.grid = element_blank())
}

.bootstrap_plot <- function(null_csv, true_csv, y_label, out_file) {
  np <- file.path(out_data, null_csv)
  tp <- file.path(out_data, true_csv)
  for (f in c(np, tp)) {
    if (!file.exists(f)) { message(basename(f), " not found"); return() }
  }
  null_df  <- read.csv(np)
  true_df  <- read.csv(tp)
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
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    labs(x = "Horizon (years)", y = y_label,
         caption = "Blue band: 95% permutation null interval. Black line: true IRF.") +
    theme_classic(base_size = 18) +
    theme(panel.grid = element_blank())
  ggsave(file.path(fig_out, out_file), p,
         width = 6, height = 5, device = cairo_pdf)
  message("Saved ", out_file)
}

plot_FSV_irf <- function() {
  path <- file.path(out_data, "figA_irf_FSV.csv")
  if (!file.exists(path)) { message("figA_irf_FSV.csv not found"); return() }
  .plot_ctrl_irf(read.csv(path), "Price variation", "figA_irf_FSV.pdf")
}

plot_bootstrap_price <- function() {
  y_price <- "Price variation"
  y_fish  <- "Price variation"
  .bootstrap_plot("figA_irf_price_bootstrap_null.csv",
                  "figA_irf_price_bootstrap_true.csv",
                  y_price, "figA_irf_bootstrap_price.pdf")
  .bootstrap_plot("figA_irf_fishprice_bootstrap_null.csv",
                  "figA_irf_fishprice_bootstrap_true.csv",
                  y_fish, "figA_irf_bootstrap_fishprice.pdf")
}

plot_loo_price <- function() {
  y_price <- "Price variation"
  y_fish  <- "Price variation"

  .make_loo <- function(base_csv, base_grp, loo_csv, grp_col, y_label, out_file) {
    for (f in c(file.path(out_data, base_csv), file.path(out_data, loo_csv))) {
      if (!file.exists(f)) { message(basename(f), " not found"); return() }
    }
    base_df <- read.csv(file.path(out_data, base_csv))
    base_df <- base_df |>
      dplyr::filter(.data[[names(base_df)[1]]] == base_grp)
    loo_df  <- read.csv(file.path(out_data, loo_csv))
    ggsave(file.path(fig_out, out_file),
           .loo_overlay_plot(base_df, loo_df, grp_col, y_label),
           width = 6, height = 5, device = cairo_pdf)
    message("Saved ", out_file)
  }

  .make_loo("fig5A_irf_grain_price.csv", "All regions",
            "figA_irf_subsample_price.csv",     "subsample",
            y_price, "figA_irf_loo_loc_price.pdf")
  .make_loo("fig5A_irf_grain_price.csv", "All regions",
            "figA_irf_subsample_20y_price.csv", "block",
            y_price, "figA_irf_loo_20y_price.pdf")
  .make_loo("fig5B_irf_fish_price.csv", "All",
            "figA_irf_subsample_fishprice.csv",     "subsample",
            y_fish, "figA_irf_loo_loc_fishprice.pdf")
  .make_loo("fig5B_irf_fish_price.csv", "All",
            "figA_irf_subsample_20y_fishprice.csv", "block",
            y_fish, "figA_irf_loo_20y_fishprice.pdf")
}


.se_robust_plot <- function(csv_file, y_label, out_file) {
  path <- file.path(out_data, csv_file)
  if (!file.exists(path)) { message(csv_file, " not found"); return() }
  df <- read.csv(path)
  p <- ggplot(df, aes(x = horizon, y = irf_mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_ribbon(aes(ymin = irf_down95, ymax = irf_up95), fill = "cornflowerblue", alpha = 0.20) +
    geom_ribbon(aes(ymin = irf_down90, ymax = irf_up90), fill = "cornflowerblue", alpha = 0.25) +
    geom_line(color = "black", linewidth = 1.1) +
    facet_wrap(~ se_type, ncol = 4) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    labs(x = "Horizon (years)", y = y_label) +
    theme_classic(base_size = 18) +
    theme(strip.background = element_blank(),
          strip.text       = element_text(face = "bold"),
          panel.grid       = element_blank())
  ggsave(file.path(fig_out, out_file), p,
         width = 10, height = 4, device = cairo_pdf)
  message("Saved ", out_file)
}

plot_se_robust_price <- function() {
  .se_robust_plot("figA_irf_price_serobust.csv",
                  "Price variation",
                  "figA_irf_serobust_price.pdf")
}

plot_se_robust_fish <- function() {
  .se_robust_plot("figA_irf_fishprice_serobust.csv",
                  "Price variation",
                  "figA_irf_serobust_fishprice.pdf")
}


###############################################################################
# Run
###############################################################################
if (!interactive()) {
  data      <- load_prices()
  fishprice <- load_fishprice()
  fishcatch <- load_fishcatch()

  export_fig5A(data)
  export_fig5B(fishprice)
  export_se_robust_price(data)
  export_se_robust_fish(fishprice)
  export_allen_irf()
  export_FSV_irf()
  export_loo_prices(data)
  export_20y_prices(data)
  export_loo_fish(fishprice)
  export_20y_fish(fishprice)
  export_fishcatch_irf(fishcatch)
  export_bootstrap(data, fishprice, n_iter = 500)
  plot_allen_irf()
  plot_FSV_irf()
  plot_bootstrap_price()
  plot_loo_price()
  plot_se_robust_price()
  plot_se_robust_fish()
  plot_fishcatch_irf()

  message("Fig 5 data exports complete.")
} else {
  message("Call functions interactively or source this script.")
}
