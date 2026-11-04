
library(fixest)
library(dplyr)
library(stringr)



lp_panel_inter <- function(
    data,
    outcome,
    main_var,
    interact_var,
    controls = NULL,
    horizon = 12,
    fe = "countrycode + year",
    panel_id = c("countrycode", "year"),
    vcov_formula = DK ~ year,
    endog = NULL,
    instruments = NULL
) {
  
  irf_list <- vector("list", horizon + 1)
  
  rhs_controls <- if (!is.null(controls)) paste0(" + ", controls) else ""
  
  for (h in 0:horizon) {
    
    # -----------------------------
    # Main LP formula
    # -----------------------------
    base_fml <- paste0(
      "f(", outcome, ", ", h, ") - l(", outcome, ", 1) ~ ",
      "i(", interact_var, ", ", main_var, ")",
      rhs_controls,
      " | ", fe
    )
    
    # -----------------------------
    # Add IV part if requested
    # -----------------------------
    if (!is.null(endog) && !is.null(instruments)) {
      fml <- as.formula(
        paste0(
          base_fml,
          " | ",
          endog, " ~ ", instruments
        )
      )
    } else {
      fml <- as.formula(base_fml)
    }
    
    fe_model <- feols(
      fml,
      data = data,
      panel.id = panel_id,
      vcov = vcov_formula
    )
    
    # -----------------------------
    # Extract interacted IRFs
    # -----------------------------
    coef_df <- data.frame(
      horizon  = h,
      term     = names(coef(fe_model)),
      irf_mean = coef(fe_model),
      se       = sqrt(diag(vcov(fe_model))),
      stringsAsFactors = FALSE
    ) %>%
      filter(str_detect(term, paste0("^", interact_var, "::"))) %>%
      mutate(
        category = str_extract(term, "(?<=::).*?(?=:)"),
        main_var = main_var,
        irf_down = irf_mean - 1.96 * se,
        irf_up   = irf_mean + 1.96 * se
      )
    
    irf_list[[h + 1]] <- coef_df
  }
  
  bind_rows(irf_list)
}

lp_panel <- function(
    data,
    outcome,
    main_var,
    controls = NULL,
    horizon = 12,
    fe = "countrycode + year",
    panel_id = c("countrycode", "year"),
    vcov_formula = DK ~ year
) {
  
  irf_list <- vector("list", horizon + 1)
  
  rhs_controls <- if (!is.null(controls)) paste0(" + ", controls) else ""
  
  for (h in 0:horizon) {
    
    fml <- as.formula(
      paste0(
        "f(", outcome, ", ", h, ") - l(", outcome, ", 1) ~ ",
        main_var,
        rhs_controls,
        " | ", fe
      )
    )
    
    fe_model <- feols(
      fml,
      data = data,
      panel.id = panel_id,
      vcov = vcov_formula
    )
    
    beta <- coef(fe_model)[main_var]
    se   <- sqrt(vcov(fe_model)[main_var, main_var])
    
    irf_list[[h + 1]] <- data.frame(
      horizon = h,
      term = main_var,
      irf_mean = beta,
      se = se,
      irf_down = beta - 1.96 * se,
      irf_up   = beta + 1.96 * se
    )
  }

  bind_rows(irf_list)
}
