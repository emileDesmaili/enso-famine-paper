
library(fixest)
library(dplyr)
library(stringr)



cumulative_dl <- function(data,
                          outcome,
                          main_var,
                          lags = 0:10,
                          fe = NULL,
                          controls = NULL,
                          panel_id = NULL,
                          vcov_formula = NULL,
                          shock_label = NULL) {
  
  # Build formula
  rhs_parts <- paste0("l(", main_var, ", ", min(lags), ":", max(lags), ")")
  
  if (!is.null(controls) && length(controls) > 0) {
    rhs_parts <- c(rhs_parts, controls)
  }
  
  rhs <- paste(rhs_parts, collapse = " + ")
  
  if (!is.null(fe)) {
    formula_str <- paste0(outcome, " ~ ", rhs, " | ", fe)
  } else {
    formula_str <- paste0(outcome, " ~ ", rhs)
  }
  
  # Estimate model
  model_args <- list(
    fml = as.formula(formula_str),
    data = data
  )
  
  if (!is.null(panel_id)) {
    model_args$panel.id <- panel_id
  }
  
  if (!is.null(vcov_formula)) {
    model_args$vcov <- vcov_formula
  }
  
  model <- do.call(feols, model_args)
  
  # Extract coefficients
  coefs <- coef(model)
  vcov_mat <- vcov(model)
  
  coef_df <- data.frame(
    term = names(coefs),
    estimate = coefs,
    stringsAsFactors = FALSE
  ) %>%
    filter(str_detect(term, main_var)) %>%
    mutate(
      lag = case_when(
        str_detect(term, "^l\\(") ~ 
          as.integer(str_extract(term, "\\d+(?=\\))")),
        str_detect(term, "^f\\(") ~ 
          -as.integer(str_extract(term, "\\d+(?=\\))")),
        term == main_var ~ 0L,
        TRUE ~ NA_integer_
      )
    ) %>%
    filter(!is.na(lag)) %>%
    arrange(lag)
  
  # Compute cumulative effects
  results <- compute_cumulative_effects(coef_df, vcov_mat, model)
  
  # Add shock label if provided
  if (!is.null(shock_label)) {
    results$shock <- shock_label
  }
  
  return(results)
}

#' Internal function to compute cumulative effects with proper variance
#'
#' @keywords internal
compute_cumulative_effects <- function(df, vcov_mat, model) {
  
  # Separate leads and lags
  df_lags <- df %>% filter(lag >= 0) %>% arrange(lag)
  df_leads <- df %>% filter(lag < 0) %>% arrange(lag)
  
  # Process lags: cumulative from lag 0
  if (nrow(df_lags) > 0) {
    n_lags <- nrow(df_lags)
    
    df_lags$cumulative_effect <- numeric(n_lags)
    df_lags$se <- numeric(n_lags)
    
    for (k in 1:n_lags) {
      # Get terms from lag 0 to lag k
      terms_to_sum <- df_lags$term[1:k]
      
      # Cumulative coefficient
      df_lags$cumulative_effect[k] <- sum(df_lags$estimate[1:k])
      
      # Cumulative SE using the sqrt(sum(vcov)) method
      df_lags$se[k] <- sqrt(sum(vcov_mat[terms_to_sum, terms_to_sum]))
    }
    
    df_lags$lower_ci <- df_lags$cumulative_effect - 1.96 * df_lags$se
    df_lags$upper_ci <- df_lags$cumulative_effect + 1.96 * df_lags$se
  }
  
  # Process leads: cumulative backwards from -1 to most negative
  if (nrow(df_leads) > 0) {
    # Reverse order: -1, -2, -3, ... becomes the cumulative order
    df_leads <- df_leads %>% arrange(desc(lag))  # -1, -2, -3, ...
    n_leads <- nrow(df_leads)
    
    df_leads$cumulative_effect <- numeric(n_leads)
    df_leads$se <- numeric(n_leads)
    
    for (k in 1:n_leads) {
      # Get terms from -1 back to current lead k
      terms_to_sum <- df_leads$term[1:k]
      
      # Cumulative coefficient
      df_leads$cumulative_effect[k] <- sum(df_leads$estimate[1:k])
      
      # Cumulative SE using the sqrt(sum(vcov)) method
      df_leads$se[k] <- sqrt(sum(vcov_mat[terms_to_sum, terms_to_sum]))
    }
    
    df_leads$lower_ci <- df_leads$cumulative_effect - 1.96 * df_leads$se
    df_leads$upper_ci <- df_leads$cumulative_effect + 1.96 * df_leads$se
    
    # Sort back to natural order for output
    df_leads <- df_leads %>% arrange(lag)
  }
  
  bind_rows(df_leads, df_lags) %>% arrange(lag)
}




cumulative_dl_interact <- function(data, outcome, main_var, interact_var,
                                   lags = 0:10, fe = NULL, controls = NULL,
                                   panel_id = NULL, vcov_formula = NULL) {
  
  
  # TODO FIX BUG WITH LEADS AND INTERACTION

  # Ensure it is sorted for lagging
  df_work <- data %>%
    arrange(!!sym(panel_id[1]), !!sym(panel_id[2]))
  
  # 2. Get unique levels of the interaction variable (e.g., "Pre-1980", "Post-1980")
  # We remove NA and ensure they are valid names
  levels <- unique(na.omit(df_work[[interact_var]]))
  
  interacted_names <- c()
  
  # 3. Manually create Lag x Category variables
  for (lev in levels) {
    # Create a clean label for the category
    lev_clean <- str_replace_all(lev, "[^[:alnum:]]", "_")
    
    for (l in lags) {
      # Handle lag vs lead naming
      lag_type <- if(l >= 0) "lag" else "lead"
      v_name <- paste0("INT_", lev_clean, "_", lag_type, "_", abs(l))
      
      # Generate the variable: (Category == Level) * Lagged_Variable
      df_work <- df_work %>%
        group_by(!!sym(panel_id[1])) %>%
        mutate(!!v_name := (get(interact_var) == lev) * lag(get(main_var), l)) %>%
        ungroup()
      
      interacted_names <- c(interacted_names, v_name)
    }
  }
  


  rhs <- paste(c(interacted_names, controls), collapse = " + ")
  fml_str <- if(!is.null(fe)) {
    paste0(outcome, " ~ ", rhs, " | ", fe)
  } else {
    paste0(outcome, " ~ ", rhs)
  }
  
  model <- feols(as.formula(fml_str), data = df_work, vcov = vcov_formula, panel.id = panel_id)
  

  coef_df <- data.frame(
    term = names(coef(model)),
    estimate = as.numeric(coef(model)),
    stringsAsFactors = FALSE
  ) %>%
    filter(str_detect(term, "^INT_")) %>%
    mutate(
      # Extract level name from our custom string "INT_LEVEL_lag_X"
      category = str_extract(term, "(?<=INT_).*(?=_(lag|lead))"),
      lag_num = as.integer(str_extract(term, "\\d+$")),
      type = str_extract(term, "(lag|lead)"),
      lag = ifelse(type == "lag", lag_num, -lag_num)
    )
  
  # Compute Cumulative Effects
  vcov_mat <- vcov(model)
  
  results <- coef_df %>%
    group_split(category) %>%
    map_dfr(function(sub_df) {
      # Sort to ensure cumsum makes sense
      sub_df <- sub_df %>% arrange(lag)
      
      # Split into leads and lags for independent accumulation
      df_lags <- sub_df %>% filter(lag >= 0) %>% arrange(lag)
      df_leads <- sub_df %>% filter(lag < 0) %>% arrange(desc(lag))
      
      calc_stats <- function(block) {
        if(nrow(block) == 0) return(NULL)
        n <- nrow(block)
        block$cumulative_effect <- cumsum(block$estimate)
        block$se <- sapply(1:n, function(i) {
          sqrt(sum(vcov_mat[block$term[1:i], block$term[1:i]]))
        })
        return(block)
      }
      
      res <- bind_rows(calc_stats(df_leads), calc_stats(df_lags)) %>%
        mutate(
          lower_ci = cumulative_effect - 1.96 * se,
          upper_ci = cumulative_effect + 1.96 * se
        )
      return(res)
    })
  
  return(results %>% arrange(category, lag))
}





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

# =============================================================================
# lp_state_dep: State-dependent LP (continuous interaction)
# =============================================================================
# Estimates feols(f(outcome,h) - l(outcome,1) ~ main_var + main_var:state_var
#                                               + controls | FE)
# and evaluates the marginal effect ME = β_main + β_inter × state_val
# at specified quantiles (or explicit eval_values) with delta-method SEs.
#
# Returns a tidy data frame with columns:
#   horizon, irf_mean, se, irf_down, irf_up, quantile, state_val
# =============================================================================
lp_state_dep <- function(data, outcome, main_var, state_var,
                         controls = NULL, horizon = 10,
                         fe = "countrycode + year",
                         panel_id = c("countrycode", "year"),
                         vcov_formula = DK ~ year,
                         eval_quantiles = c(0.10, 0.50, 0.90),
                         eval_labels = NULL,
                         eval_values = NULL) {
  rhs_controls <- if (!is.null(controls)) paste0(" + ", controls) else ""

  if (is.null(eval_values)) {
    state_vals <- quantile(data[[state_var]], eval_quantiles, na.rm = TRUE)
  } else {
    state_vals <- eval_values
  }
  if (is.null(eval_labels)) eval_labels <- paste0("p", eval_quantiles * 100)

  irf_list <- list()
  for (h in 0:horizon) {
    fml <- as.formula(paste0(
      "f(", outcome, ", ", h, ") - l(", outcome, ", 1) ~ ",
      main_var, " + ", main_var, ":", state_var,
      rhs_controls, " | ", fe
    ))
    mod <- tryCatch(
      feols(fml, data = data, panel.id = panel_id, vcov = vcov_formula),
      error = function(e) { message("lp_state_dep h=", h, ": ", e$message); NULL }
    )
    if (is.null(mod)) next

    b_main  <- coef(mod)[main_var]
    int_nm  <- grep(paste0(main_var, ":", state_var, "|", state_var, ":", main_var),
                    names(coef(mod)), value = TRUE)[1]
    if (is.na(int_nm)) next
    b_inter <- coef(mod)[int_nm]
    V       <- vcov(mod)[c(main_var, int_nm), c(main_var, int_nm)]

    for (j in seq_along(state_vals)) {
      me   <- b_main + b_inter * state_vals[j]
      grad <- c(1, state_vals[j])
      se   <- as.numeric(sqrt(t(grad) %*% V %*% grad))
      irf_list[[length(irf_list) + 1]] <- data.frame(
        horizon   = h,
        irf_mean  = me,
        se        = se,
        irf_down  = me - 1.96 * se,
        irf_up    = me + 1.96 * se,
        quantile  = eval_labels[j],
        state_val = state_vals[j]
      )
    }
  }
  bind_rows(irf_list)
}
