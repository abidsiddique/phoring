#' Plain-English Econometric Model Interpretation
#'
#' Automatically inspects regression functional forms (log-log, log-level, 
#' level-log, level-level) and provides human-readable, plain-English 
#' economic interpretations for every coefficient, including elasticities,
#' semi-elasticities, percentage effects, and quadratic turning points.
#'
#' @param model An \code{econ_reg} object produced by \code{\link{reg}} or a standard \code{lm} model.
#' @param digits Number of decimal places to format (default is 4).
#'
#' @return An object of class \code{econ_interpret} containing interpretation summaries.
#' @export
#'
#' @examples
#' data(wage_data, package = "phoring")
#' # Quadratic wage equation
#' mod <- reg(wage ~ educ + exper + I(exper^2), data = wage_data)
#' interpret(mod)
#'
#' # Log-log elasticity model
#' data(house_prices, package = "phoring")
#' mod_log <- reg(log(price) ~ log(rooms) + crime, data = house_prices)
#' interpret(mod_log)
interpret <- function(model, digits = 4) {
  if (inherits(model, "econ_reg")) {
    base_lm <- model$lm_model
    formula <- model$formula
    coef_tab <- model$coef_table
    se_type <- model$se_type
  } else if (inherits(model, "lm")) {
    base_lm <- model
    formula <- stats::formula(model)
    sm <- summary(model)
    coef_tab <- as.data.frame(sm$coefficients)
    colnames(coef_tab) <- c("Estimate", "Std_Error", "t_stat", "p_value")
    coef_tab$Signif <- sapply(coef_tab$p_value, get_stars)
    se_type <- "Classical (OLS)"
  } else {
    stop("`model` must be an object of class `econ_reg` or `lm`.", call. = FALSE)
  }
  
  # Inspect LHS
  lhs_term <- deparse(formula[[2]])
  is_log_lhs <- grepl("^log\\(", lhs_term)
  
  # Extract predictor info
  mf <- stats::model.frame(base_lm)
  terms_obj <- stats::terms(formula)
  pred_names <- rownames(coef_tab)
  
  interpretations <- list()
  quad_turning_points <- list()
  
  # Check for quadratic terms: e.g., exper and I(exper^2)
  for (i in seq_along(pred_names)) {
    vname <- pred_names[i]
    if (vname == "(Intercept)") next
    
    b_val <- coef_tab[vname, "Estimate"]
    p_val <- coef_tab[vname, "p_value"]
    sig_str <- if (p_val < 0.05) "statistically significant (p < 0.05)" else "not statistically significant (p >= 0.05)"
    
    # Check if this variable is part of a quadratic pair
    quad_match <- regexec("^I\\((.*)\\^2\\)$", vname)
    quad_matches <- regmatches(vname, quad_match)[[1]]
    
    if (length(quad_matches) == 2) {
      base_var <- quad_matches[2]
      if (base_var %in% pred_names) {
        b1 <- coef_tab[base_var, "Estimate"]
        b2 <- b_val
        turning_point <- -b1 / (2 * b2)
        shape <- if (b2 < 0) "Inverted U-shaped (Maximum)" else "U-shaped (Minimum)"
        
        # Check empirical range if available
        range_str <- ""
        if (base_var %in% names(mf)) {
          x_vals <- mf[[base_var]]
          min_x <- min(x_vals, na.rm = TRUE)
          max_x <- max(x_vals, na.rm = TRUE)
          in_range <- (turning_point >= min_x && turning_point <= max_x)
          range_str <- sprintf(" [Sample range: %.1f to %.1f; %s]", 
                               min_x, max_x, 
                               if (in_range) "turning point is within sample range" else "turning point is outside sample range")
        }
        
        quad_turning_points[[base_var]] <- list(
          linear_var = base_var,
          quad_var = vname,
          b1 = b1,
          b2 = b2,
          turning_point = turning_point,
          shape = shape,
          range_info = range_str
        )
      }
    }
    
    # Check RHS functional form
    is_log_rhs <- grepl("^log\\(", vname)
    is_quad_rhs <- grepl("\\^2", vname)
    
    # Check if dummy variable (only takes 0 and 1)
    is_dummy <- FALSE
    clean_vname <- gsub("^factor\\((.*?)\\)", "\\1", vname)
    clean_vname <- gsub("^as\\.factor\\((.*?)\\)", "\\1", clean_vname)
    if (vname %in% names(mf)) {
      is_dummy <- all(mf[[vname]] %in% c(0, 1, NA))
    }
    
    # Classify specification
    if (is_log_lhs && is_log_rhs) {
      type <- "Log-Log (Elasticity)"
      summary_text <- sprintf("A 1%% increase in %s is associated with a %.*f%% change in %s.",
                              vname, digits, b_val, lhs_term)
      detail <- sprintf("Constant elasticity: Each 1%% increase in predictor yields a %.*f%% expected change in outcome (%s).",
                        digits, b_val, sig_str)
    } else if (is_log_lhs && !is_log_rhs) {
      type <- "Log-Level (Semi-Elasticity)"
      pct_approx <- b_val * 100
      pct_exact  <- (exp(b_val) - 1) * 100
      if (is_dummy) {
        summary_text <- sprintf("Being %s (vs. base) is associated with a %.*f%% change in %s.",
                                vname, digits, pct_exact, lhs_term)
        detail <- sprintf("Halvorsen-Palmquist adjustment for dummy: [exp(%.*f) - 1]*100 = %.*f%% difference in outcome (%s).",
                          digits, b_val, digits, pct_exact, sig_str)
      } else {
        summary_text <- sprintf("A 1-unit increase in %s is associated with a %.*f%% change in %s.",
                                vname, digits, pct_approx, lhs_term)
        detail <- sprintf("Semi-elasticity: Approx %.*f%% change (exact: %.*f%%) per 1-unit increase in predictor (%s).",
                          digits, pct_approx, digits, pct_exact, sig_str)
      }
    } else if (!is_log_lhs && is_log_rhs) {
      type <- "Level-Log"
      unit_change <- b_val / 100
      summary_text <- sprintf("A 1%% increase in %s is associated with a %.*f unit change in %s.",
                              vname, digits, unit_change, lhs_term)
      detail <- sprintf("Level-Log effect: beta/100 = %.*f unit change in outcome per 1%% increase in predictor (%s).",
                        digits, unit_change, sig_str)
    } else {
      type <- if (is_quad_rhs) "Polynomial (Quadratic)" else "Level-Level (Marginal Effect)"
      summary_text <- sprintf("A 1-unit increase in %s is associated with a %.*f unit change in %s.",
                              vname, digits, b_val, lhs_term)
      detail <- sprintf("Linear marginal effect: A 1-unit change in predictor changes outcome by %.*f units (%s).",
                        digits, b_val, sig_str)
    }
    
    interpretations[[vname]] <- list(
      variable = vname,
      estimate = b_val,
      p_value = p_val,
      type = type,
      summary = summary_text,
      detail = detail
    )
  }
  
  res <- list(
    formula = formula,
    lhs = lhs_term,
    interpretations = interpretations,
    turning_points = quad_turning_points,
    digits = digits
  )
  class(res) <- "econ_interpret"
  return(res)
}

#' @export
print.econ_interpret <- function(x, ...) {
  print_line("=")
  cat("                PLAIN-ENGLISH ECONOMETRIC INTERPRETATION                \n")
  cat(paste0("  Model: ", deparse(x$formula), "\n"))
  cat(paste0("  Dependent Variable: ", x$lhs, "\n"))
  print_line("=")
  
  v_w <- max(16, max(nchar(names(x$interpretations))))
  line_w <- max(74, v_w + 58)
  
  cat(sprintf(paste0("  %-", v_w, "s %-22s %10s  %-20s\n"), 
              "Variable", "Specification Type", "Estimate", "Interpretation Summary"))
  print_line("-", width = line_w)
  
  for (item in x$interpretations) {
    cat(sprintf(paste0("  %-", v_w, "s %-22s %10.*f  %s\n"),
                item$variable,
                item$type,
                x$digits, item$estimate,
                item$summary))
  }
  print_line("-", width = line_w)
  
  # Quadratic Turning Points Card
  if (length(x$turning_points) > 0) {
    cat("\n")
    print_line("=")
    cat("                     NONLINEAR TURNING POINT ANALYSIS                   \n")
    print_line("=")
    for (tp in x$turning_points) {
      cat(sprintf("  Variable: %s (Linear: %.4f | Quadratic %s: %.4f)\n", 
                  tp$linear_var, tp$b1, tp$quad_var, tp$b2))
      cat(sprintf("  * Shape:          %s\n", tp$shape))
      cat(sprintf("  * Turning Point:  x* = -beta1 / (2 * beta2) = %.*f\n", x$digits, tp$turning_point))
      if (nchar(tp$range_info) > 0) {
        cat(sprintf("  * Empirical Check:%s\n", tp$range_info))
      }
      if (grepl("Maximum", tp$shape)) {
        cat(sprintf("  --> Economic Insight: %s increases %s up to %.*f, after which diminishing returns dominate.\n",
                    tp$linear_var, x$lhs, x$digits, tp$turning_point))
      } else {
        cat(sprintf("  --> Economic Insight: %s decreases %s down to %.*f, after which it begins rising.\n",
                    tp$linear_var, x$lhs, x$digits, tp$turning_point))
      }
    }
    print_line("=")
  }
  invisible(x)
}

#' Detect Outliers and Influential Observations in Regression Models
#'
#' Automatically computes Cook's Distance, studentized residuals, and leverage 
#' (hat values) for regression models, flagging influential data points that could 
#' distort econometric estimates.
#'
#' @param model An \code{econ_reg} or \code{lm} model object.
#' @param top_n Number of top influential observations to display (default is 5).
#' @param digits Number of decimal places to format (default is 4).
#'
#' @return An object of class \code{econ_outliers}.
#' @export
#'
#' @examples
#' data(wage_data, package = "phoring")
#' mod <- reg(wage ~ educ + exper, data = wage_data)
#' outliers(mod)
outliers <- function(model, top_n = 5, digits = 4) {
  if (inherits(model, "econ_reg")) {
    base_lm <- model$lm_model
    e <- model$residuals
    n <- model$n
    k <- model$k
  } else if (inherits(model, "lm")) {
    base_lm <- model
    e <- stats::residuals(model)
    n <- length(e)
    k <- length(stats::coef(model))
  } else {
    stop("`model` must be an object of class `econ_reg` or `lm`.", call. = FALSE)
  }
  
  # Calculate leverage (hat values)
  X <- stats::model.matrix(base_lm)
  bread <- safe_solve(crossprod(X))
  h <- rowSums((X %*% bread) * X)
  h <- pmin(h, 0.9999) # prevent divide by zero
  
  # Residual variance
  s2 <- sum(e^2) / (n - k)
  s <- sqrt(s2)
  
  # Studentized residuals (internally studentized)
  r <- e / (s * sqrt(1 - h))
  
  # Cook's Distance: D_i = (r_i^2 / k) * (h_i / (1 - h_i))
  cooks_d <- (r^2 / k) * (h / (1 - h))
  
  # Thresholds
  cooks_threshold <- 4 / n
  leverage_threshold <- 2 * k / n
  resid_threshold <- 2.5
  
  # Flags
  flag_cooks <- cooks_d > cooks_threshold
  flag_leverage <- h > leverage_threshold
  flag_resid <- abs(r) > resid_threshold
  
  flag_summary <- character(n)
  for (i in 1:n) {
    f_list <- character(0)
    if (cooks_d[i] > 1) {
      f_list <- c(f_list, "Extreme Cook's D (> 1)")
    } else if (flag_cooks[i]) {
      f_list <- c(f_list, "High Cook's D (> 4/n)")
    }
    if (flag_leverage[i]) f_list <- c(f_list, "High Leverage (> 2k/n)")
    if (flag_resid[i]) f_list <- c(f_list, "Large Residual (|r| > 2.5)")
    
    if (length(f_list) > 0) {
      flag_summary[i] <- paste(f_list, collapse = "; ")
    } else {
      flag_summary[i] <- "Normal"
    }
  }
  
  # Top N by Cook's distance
  order_idx <- order(cooks_d, decreasing = TRUE)
  top_indices <- order_idx[seq_len(min(top_n, n))]
  
  top_table <- data.frame(
    Obs_ID        = top_indices,
    Cooks_D       = cooks_d[top_indices],
    Studentized_R = r[top_indices],
    Leverage      = h[top_indices],
    Flag          = flag_summary[top_indices],
    stringsAsFactors = FALSE
  )
  
  res <- list(
    top_table = top_table,
    n = n,
    k = k,
    n_high_cooks = sum(flag_cooks),
    n_high_leverage = sum(flag_leverage),
    n_high_resid = sum(flag_resid),
    cooks_threshold = cooks_threshold,
    leverage_threshold = leverage_threshold,
    resid_threshold = resid_threshold,
    cooks_d = cooks_d,
    studentized_residuals = r,
    leverage = h,
    digits = digits
  )
  class(res) <- "econ_outliers"
  return(res)
}

#' @export
print.econ_outliers <- function(x, ...) {
  d <- x$digits
  print_line("=")
  cat("             OUTLIER & INFLUENTIAL OBSERVATION DIAGNOSTICS              \n")
  cat(sprintf("  Sample Size (n): %-10d Predictors (k): %d\n", x$n, x$k))
  cat(sprintf("  Cook's D Threshold (4/n):       %.4f (Flagged: %d obs)\n", 
              x$cooks_threshold, x$n_high_cooks))
  cat(sprintf("  Leverage Threshold (2k/n):      %.4f (Flagged: %d obs)\n", 
              x$leverage_threshold, x$n_high_leverage))
  cat(sprintf("  Studentized Resid (|r| > 2.5):  2.5000 (Flagged: %d obs)\n", 
              x$n_high_resid))
  print_line("=")
  
  cat(sprintf("  TOP %d MOST INFLUENTIAL OBSERVATIONS (BY COOK'S DISTANCE):\n", nrow(x$top_table)))
  print_line("-")
  cat(sprintf("  %-8s %10s %14s %10s   %-30s\n", 
              "Obs #", "Cook's D", "Studentized R", "Leverage", "Diagnosis"))
  print_line("-")
  
  for (i in 1:nrow(x$top_table)) {
    row <- x$top_table[i, ]
    cat(sprintf("  %-8d %10.*f %14.*f %10.*f   %-30s\n",
                row$Obs_ID,
                d, row$Cooks_D,
                d, row$Studentized_R,
                d, row$Leverage,
                row$Flag))
  }
  print_line("-")
  
  cat("  Econometric Guidance:\n")
  if (x$n_high_cooks == 0) {
    cat("  * No excessively influential points detected; model stability is good.\n")
  } else {
    cat("  * Noticeable influential points detected. To ensure robustness, compare estimates\n")
    cat("    with and without the top flagged observation(s) or use `se = 'robust'`.\n")
  }
  print_line("=")
  invisible(x)
}

#' Chow Test for Parameter Stability and Structural Breaks
#'
#' Tests whether regression coefficients are constant across two subsamples
#' (e.g. before/after a policy shock, male vs. female, or high vs. low income).
#'
#' @param formula Regression formula (or an existing \code{econ_reg} or \code{lm} model).
#' @param split_by A binary variable or condition in \code{data} defining the two subsamples.
#' @param data A data.frame containing the variables (required if formula is used).
#' @param digits Number of decimal places to format (default is 4).
#'
#' @return An object of class \code{econ_chow}.
#' @export
#'
#' @examples
#' data(wage_data, package = "phoring")
#' # Test if wage equation differs by gender
#' chow_test(wage ~ educ + exper, split_by = "female", data = wage_data)
chow_test <- function(formula, split_by, data = NULL, digits = 4) {
  if (inherits(formula, "econ_reg") || inherits(formula, "lm")) {
    base_mod <- formula
    f <- if (inherits(base_mod, "econ_reg")) base_mod$formula else stats::formula(base_mod)
    if (is.null(data)) {
      data <- stats::model.frame(base_mod)
    }
  } else {
    f <- formula
    if (is.null(data)) stop("`data` must be provided when `formula` is specified.", call. = FALSE)
  }
  
  # Resolve split_by
  if (is.character(split_by) && length(split_by) == 1) {
    if (!split_by %in% names(data)) {
      stop(sprintf("Column `%s` not found in data.", split_by), call. = FALSE)
    }
    split_vec <- data[[split_by]]
    split_name <- split_by
  } else {
    split_vec <- split_by
    split_name <- "Subsample Condition"
  }
  
  # Must evaluate to binary
  split_fac <- as.factor(split_vec)
  if (nlevels(split_fac) != 2) {
    stop("`split_by` must define exactly two distinct subsamples (e.g. 0 and 1, or 2 levels).", call. = FALSE)
  }
  
  levels_found <- levels(split_fac)
  idx1 <- which(split_fac == levels_found[1])
  idx2 <- which(split_fac == levels_found[2])
  
  # Fit pooled model
  mod_pooled <- stats::lm(f, data = data)
  ssr_pooled <- sum(stats::residuals(mod_pooled)^2)
  k <- length(stats::coef(mod_pooled))
  n <- nrow(data)
  
  # Fit group 1
  mod_g1 <- stats::lm(f, data = data[idx1, , drop = FALSE])
  ssr_g1 <- sum(stats::residuals(mod_g1)^2)
  n1 <- length(idx1)
  
  # Fit group 2
  mod_g2 <- stats::lm(f, data = data[idx2, , drop = FALSE])
  ssr_g2 <- sum(stats::residuals(mod_g2)^2)
  n2 <- length(idx2)
  
  # Chow F-statistic:
  # F = [(SSR_p - (SSR_1 + SSR_2)) / k] / [(SSR_1 + SSR_2) / (n1 + n2 - 2k)]
  ssr_unres <- ssr_g1 + ssr_g2
  df1 <- k
  df2 <- n1 + n2 - (2 * k)
  
  if (df2 <= 0) {
    stop("Insufficient degrees of freedom to perform Chow test.", call. = FALSE)
  }
  
  f_stat <- ((ssr_pooled - ssr_unres) / df1) / (ssr_unres / df2)
  p_val <- stats::pf(f_stat, df1 = df1, df2 = df2, lower.tail = FALSE)
  
  res <- list(
    formula = f,
    split_name = split_name,
    group1_label = levels_found[1],
    group2_label = levels_found[2],
    n1 = n1,
    n2 = n2,
    n_total = n,
    k = k,
    ssr_pooled = ssr_pooled,
    ssr_g1 = ssr_g1,
    ssr_g2 = ssr_g2,
    f_stat = f_stat,
    df1 = df1,
    df2 = df2,
    p_value = p_val,
    data = data,
    mod_pooled = mod_pooled,
    mod_g1 = mod_g1,
    mod_g2 = mod_g2,
    split_fac = split_fac,
    digits = digits
  )
  class(res) <- "econ_chow"
  return(res)
}

#' @export
print.econ_chow <- function(x, ...) {
  d <- x$digits
  print_line("=")
  cat("                CHOW TEST FOR STRUCTURAL BREAK / STABILITY              \n")
  cat(paste0("  Formula: ", deparse(x$formula), "\n"))
  cat(sprintf("  Split Variable: %s ('%s' vs. '%s')\n", 
              x$split_name, x$group1_label, x$group2_label))
  cat(sprintf("  Sample Sizes:   Group 1: %d  |  Group 2: %d  |  Pooled: %d\n", 
              x$n1, x$n2, x$n_total))
  print_line("=")
  
  cat(sprintf("  H0: Coefficients are identical across '%s' and '%s' (Parameter Stability)\n",
              x$group1_label, x$group2_label))
  cat(sprintf("  H1: Coefficients differ across groups (Structural Break)\n"))
  print_line("-")
  cat(sprintf("  Pooled SSR:    %14.*f\n", d, x$ssr_pooled))
  cat(sprintf("  Group 1 SSR:   %14.*f\n", d, x$ssr_g1))
  cat(sprintf("  Group 2 SSR:   %14.*f\n", d, x$ssr_g2))
  cat(sprintf("  Numerator df:  %14d\n", x$df1))
  cat(sprintf("  Denominator df:%14d\n", x$df2))
  cat(sprintf("  Chow F-stat:   %14.*f\n", d, x$f_stat))
  cat(sprintf("  p-value:       %14s  %s\n", format_pval(x$p_value), get_stars(x$p_value)))
  print_line("-")
  
  if (x$p_value < 0.05) {
    verdict <- "REJECT H0 (p < 0.05) -> Significant structural break detected."
    advice <- "Coefficients differ systematically between groups; pooled estimation is misspecified."
  } else {
    verdict <- "FAIL TO REJECT H0 (p >= 0.05) -> Parameter stability confirmed."
    advice <- "No evidence of structural break; pooling data across groups is appropriate."
  }
  cat(sprintf("  Verdict: %s\n", verdict))
  cat(sprintf("  Note:    %s\n", advice))
  print_line("=")
  invisible(x)
}
