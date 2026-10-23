#' Econometric Diagnostic Tests
#'
#' Computes a complete battery of econometric diagnostic tests for linear models,
#' including Breusch-Pagan, White, Durbin-Watson, Breusch-Godfrey, Ramsey RESET,
#' Multicollinearity (VIF), and Jarque-Bera normality tests.
#'
#' @param model A fitted model object from `reg()` or `stats::lm()`.
#' @param ... Additional arguments.
#'
#' @importFrom stats model.frame model.response model.matrix lm embed
#' @return An object of class `econ_diagnostics` containing test statistics,
#'   p-values, decisions, and individual VIF values.
#' @export
#'
#' @examples
#' data(mtcars)
#' m <- reg(mpg ~ wt + hp + qsec, data = mtcars)
#' diagnostics(m)
diagnostics <- function(model, ...) {
  UseMethod("diagnostics")
}

#' @export
diagnostics.default <- function(model, alpha = 0.05, ...) {
  if (!inherits(model, c("lm", "econ_reg"))) {
    stop("`diagnostics()` requires an `econ_reg` or `lm` model object.", call. = FALSE)
  }
  
  # Extract components
  if (inherits(model, "econ_reg")) {
    base_lm <- model$lm_model
  } else {
    base_lm <- model
  }
  
  e <- stats::residuals(base_lm)
  y_hat <- stats::fitted(base_lm)
  mf <- stats::model.frame(base_lm)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(base_lm)
  n <- length(e)
  k <- ncol(X)
  
  # Check for degrees of freedom
  if (n <= k + 2) {
    stop("Not enough degrees of freedom to run diagnostic tests.", call. = FALSE)
  }
  
  test_rows <- list()
  
  # -------------------------------------------------------------
  # 1. Breusch-Pagan Test (Heteroskedasticity)
  # -------------------------------------------------------------
  e2 <- e^2
  aux_bp <- stats::lm(e2 ~ X[, -1, drop = FALSE])
  r2_bp <- summary(aux_bp)$r.squared
  df_bp <- k - 1
  stat_bp <- n * r2_bp
  p_bp <- stats::pchisq(stat_bp, df = df_bp, lower.tail = FALSE)
  verdict_bp <- if (p_bp < alpha) "REJECT (Het present)" else "FAIL TO REJECT (Homoskedastic)"
  
  test_rows[["BP"]] <- data.frame(
    Test = "Breusch-Pagan",
    Null = "Homoskedasticity",
    Stat = stat_bp,
    df = df_bp,
    p_value = p_bp,
    Verdict = verdict_bp,
    stringsAsFactors = FALSE
  )
  
  # -------------------------------------------------------------
  # 2. White Test (Special Form: Fitted Values y_hat and y_hat^2)
  # -------------------------------------------------------------
  aux_white <- stats::lm(e2 ~ y_hat + I(y_hat^2))
  r2_white <- summary(aux_white)$r.squared
  stat_white <- n * r2_white
  df_white <- 2
  p_white <- stats::pchisq(stat_white, df = df_white, lower.tail = FALSE)
  verdict_white <- if (p_white < alpha) "REJECT (Het present)" else "FAIL TO REJECT (Homoskedastic)"
  
  test_rows[["White"]] <- data.frame(
    Test = "White Test (Special)",
    Null = "Homoskedasticity",
    Stat = stat_white,
    df = df_white,
    p_value = p_white,
    Verdict = verdict_white,
    stringsAsFactors = FALSE
  )
  
  # -------------------------------------------------------------
  # 3. Durbin-Watson Test (Autocorrelation)
  # -------------------------------------------------------------
  diff_e <- diff(e)
  stat_dw <- sum(diff_e^2) / sum(e^2)
  # Normal approximation for DW p-value around 2
  z_dw <- (stat_dw - 2) / (2 / sqrt(n))
  p_dw <- 2 * stats::pnorm(-abs(z_dw))
  verdict_dw <- if (p_dw < alpha) "REJECT (Autocorr present)" else "FAIL TO REJECT (No Autocorr)"
  
  test_rows[["DW"]] <- data.frame(
    Test = "Durbin-Watson",
    Null = "No Autocorrelation",
    Stat = stat_dw,
    df = 1,
    p_value = p_dw,
    Verdict = verdict_dw,
    stringsAsFactors = FALSE
  )
  
  # -------------------------------------------------------------
  # 4. Breusch-Godfrey Test (Serial Correlation, order p = 2)
  # -------------------------------------------------------------
  order_bg <- min(2, n - k - 1)
  if (order_bg >= 1) {
    lag_e <- embed(c(rep(0, order_bg), e), order_bg + 1)[, -1, drop = FALSE]
    colnames(lag_e) <- paste0("lag_e", 1:order_bg)
    aux_bg <- stats::lm(e ~ X[, -1, drop = FALSE] + lag_e)
    r2_bg <- summary(aux_bg)$r.squared
    stat_bg <- n * r2_bg
    p_bg <- stats::pchisq(stat_bg, df = order_bg, lower.tail = FALSE)
    verdict_bg <- if (p_bg < alpha) "REJECT (Autocorr present)" else "FAIL TO REJECT (No Autocorr)"
    
    test_rows[["BG"]] <- data.frame(
      Test = paste0("Breusch-Godfrey (p=", order_bg, ")"),
      Null = "No Autocorrelation",
      Stat = stat_bg,
      df = order_bg,
      p_value = p_bg,
      Verdict = verdict_bg,
      stringsAsFactors = FALSE
    )
  }
  
  # -------------------------------------------------------------
  # 5. Ramsey RESET Test (Functional Form)
  # -------------------------------------------------------------
  if (n > k + 3) {
    aux_reset <- stats::lm(y ~ X[, -1, drop = FALSE] + I(y_hat^2) + I(y_hat^3))
    r2_unres <- summary(aux_reset)$r.squared
    r2_res <- summary(base_lm)$r.squared
    df_num <- 2
    df_denom <- n - k - 2
    stat_reset <- ((r2_unres - r2_res) / df_num) / ((1 - r2_unres) / df_denom)
    stat_reset <- max(0, stat_reset)
    p_reset <- stats::pf(stat_reset, df1 = df_num, df2 = df_denom, lower.tail = FALSE)
    verdict_reset <- if (p_reset < alpha) "REJECT (Misspecified)" else "FAIL TO REJECT (Specification OK)"
    
    test_rows[["RESET"]] <- data.frame(
      Test = "Ramsey RESET",
      Null = "Correct Functional Form",
      Stat = stat_reset,
      df = df_num,
      p_value = p_reset,
      Verdict = verdict_reset,
      stringsAsFactors = FALSE
    )
  }
  
  # -------------------------------------------------------------
  # 6. Jarque-Bera Test (Normality of Residuals)
  # -------------------------------------------------------------
  m2_res <- mean(e^2)
  m3_res <- mean(e^3)
  m4_res <- mean(e^4)
  skew <- m3_res / (m2_res^(1.5))
  kurt <- m4_res / (m2_res^2)
  stat_jb <- (n / 6) * (skew^2 + ((kurt - 3)^2) / 4)
  p_jb <- stats::pchisq(stat_jb, df = 2, lower.tail = FALSE)
  verdict_jb <- if (p_jb < alpha) "REJECT (Non-normal)" else "FAIL TO REJECT (Normality OK)"
  
  test_rows[["JB"]] <- data.frame(
    Test = "Jarque-Bera",
    Null = "Normal Residuals",
    Stat = stat_jb,
    df = 2,
    p_value = p_jb,
    Verdict = verdict_jb,
    stringsAsFactors = FALSE
  )
  
  # -------------------------------------------------------------
  # 7. Variance Inflation Factors (VIF)
  # -------------------------------------------------------------
  vif_vals <- numeric()
  preds <- colnames(X)
  preds <- preds[preds != "(Intercept)"]
  
  if (length(preds) > 1) {
    for (p_name in preds) {
      y_p <- X[, p_name]
      x_others <- X[, setdiff(preds, p_name), drop = FALSE]
      mod_p <- lm(y_p ~ x_others)
      r2_p <- summary(mod_p)$r.squared
      vif_vals[p_name] <- 1 / max(1e-4, 1 - r2_p)
    }
  } else if (length(preds) == 1) {
    vif_vals[preds] <- 1.00
  }
  
  table_df <- do.call(rbind, test_rows)
  rownames(table_df) <- NULL
  
  res <- list(
    tests = table_df,
    vifs = vif_vals,
    n = n,
    k = k
  )
  class(res) <- "econ_diagnostics"
  return(res)
}

#' @export
print.econ_diagnostics <- function(x, ...) {
  print_line("=")
  cat("                   ALL-IN-ONE DIAGNOSTIC TESTS                          \n")
  print_line("=")
  cat(sprintf("  %-23s %-24s %9s %10s  %-20s\n", "Test Name", "Null Hypothesis (H0)", "Stat", "p-value", "Result"))
  print_line("-")
  
  for (i in 1:nrow(x$tests)) {
    row <- x$tests[i, ]
    cat(sprintf("  %-23s %-24s %9.3f %10s  %-20s\n",
                row$Test,
                row$Null,
                row$Stat,
                format_pval(row$p_value),
                row$Verdict))
  }
  print_line("-")
  
  # Print VIF summary
  if (length(x$vifs) > 0) {
    cat("  Multicollinearity (VIF):\n    ")
    vif_strs <- paste0(names(x$vifs), ": ", sprintf("%.2f", x$vifs))
    cat(paste(vif_strs, collapse = "  |  "))
    max_vif <- max(x$vifs, na.rm = TRUE)
    if (max_vif < 5) {
      cat("\n    --> Assessment: Safe (all VIF < 5)\n")
    } else if (max_vif < 10) {
      cat("\n    --> Assessment: Moderate multicollinearity detected (5 <= VIF < 10)\n")
    } else {
      cat("\n    --> Assessment: Severe multicollinearity detected (VIF >= 10)\n")
    }
  }
  print_line("=")
  invisible(x)
}

#' Goldfeld-Quandt Test for Heteroskedasticity
#'
#' Tests for heteroskedasticity by ordering the sample according to a regressor,
#' omitting central observations, and testing whether the error variance differs 
#' significantly between the lower and upper subsamples (Goldfeld and Quandt, 1965).
#'
#' @param formula Model formula (e.g., \code{wage ~ educ + exper}) or a fitted model object.
#' @param split_by Character name of the explanatory variable to sort by. If \code{NULL},
#'   uses the first continuous predictor.
#' @param data A data.frame containing the variables (required if formula is used).
#' @param omit_fraction Fraction of central observations to omit (default is 0.20, roughly 1/5th).
#' @param digits Number of decimal places to format (default is 4).
#'
#' @return An object of class \code{econ_gq}.
#' @export
#'
#' @examples
#' data(wage_data, package = "phoring")
#' gq_test(wage ~ educ + exper, split_by = "educ", data = wage_data)
gq_test <- function(formula, split_by = NULL, data = NULL, omit_fraction = 0.20, digits = 4) {
  if (inherits(formula, "econ_reg") || inherits(formula, "lm")) {
    base_mod <- formula
    f <- if (inherits(base_mod, "econ_reg")) base_mod$formula else stats::formula(base_mod)
    if (is.null(data)) data <- stats::model.frame(base_mod)
  } else {
    f <- formula
    if (is.null(data)) stop("`data` must be provided when `formula` is specified.", call. = FALSE)
  }
  
  mf <- stats::model.frame(f, data = data)
  n <- nrow(mf)
  
  # Determine sorting variable
  all_preds <- names(mf)[-1]
  if (is.null(split_by)) {
    split_var_name <- all_preds[1]
  } else {
    split_var_name <- split_by
    if (!split_var_name %in% names(data)) {
      stop(sprintf("Sorting variable `%s` not found in data.", split_var_name), call. = FALSE)
    }
  }
  
  sort_vec <- data[[split_var_name]]
  order_idx <- order(sort_vec)
  data_sorted <- data[order_idx, , drop = FALSE]
  
  # Number of observations to omit
  c_omit <- floor(n * omit_fraction)
  n_sub <- floor((n - c_omit) / 2)
  
  if (n_sub < 5) {
    stop("Sample size too small to omit central observations for Goldfeld-Quandt test.", call. = FALSE)
  }
  
  idx1 <- 1:n_sub
  idx2 <- (n - n_sub + 1):n
  
  # Fit lower subsample (Group 1)
  mod1 <- stats::lm(f, data = data_sorted[idx1, , drop = FALSE])
  ssr1 <- sum(stats::residuals(mod1)^2)
  k <- length(stats::coef(mod1))
  df1 <- n_sub - k
  
  # Fit upper subsample (Group 2)
  mod2 <- stats::lm(f, data = data_sorted[idx2, , drop = FALSE])
  ssr2 <- sum(stats::residuals(mod2)^2)
  df2 <- n_sub - k
  
  if (df1 <= 0 || df2 <= 0) {
    stop("Insufficient degrees of freedom in subsamples for Goldfeld-Quandt test.", call. = FALSE)
  }
  
  # GQ F-statistic: F = (SSR2 / df2) / (SSR1 / df1)
  f_stat <- (ssr2 / df2) / (ssr1 / df1)
  
  # Two-tailed p-value
  if (f_stat >= 1) {
    p_val <- 2 * stats::pf(f_stat, df1 = df2, df2 = df1, lower.tail = FALSE)
  } else {
    p_val <- 2 * stats::pf(f_stat, df1 = df2, df2 = df1, lower.tail = TRUE)
  }
  p_val <- min(1, p_val)
  
  res <- list(
    formula = f,
    sort_variable = split_var_name,
    n_total = n,
    n_sub = n_sub,
    n_omitted = n - (2 * n_sub),
    ssr1 = ssr1,
    ssr2 = ssr2,
    f_stat = f_stat,
    df1 = df2,
    df2 = df1,
    p_value = p_val,
    idx1 = idx1,
    idx2 = idx2,
    sort_vals = sort_vec[order_idx],
    e_all = stats::residuals(stats::lm(f, data = data_sorted)),
    digits = digits
  )
  class(res) <- "econ_gq"
  return(res)
}

#' @export
print.econ_gq <- function(x, ...) {
  d <- x$digits
  print_line("=")
  cat("                GOLDFELD-QUANDT HETEROSCEDASTICITY TEST                 \n")
  cat(paste0("  Formula: ", deparse(x$formula), "\n"))
  cat(sprintf("  Ordered by: %s | Omitted Middle Obs: %d | Subsample n: %d each\n", 
              x$sort_variable, x$n_omitted, x$n_sub))
  print_line("=")
  cat(sprintf("  H0: Homoskedasticity (Equal error variance across %s)\n", x$sort_variable))
  cat(sprintf("  H1: Heteroskedasticity (Error variance increases with %s)\n", x$sort_variable))
  print_line("-")
  cat(sprintf("  Lower Subsample SSR: %12.*f  (df = %d)\n", d, x$ssr1, x$df2))
  cat(sprintf("  Upper Subsample SSR: %12.*f  (df = %d)\n", d, x$ssr2, x$df1))
  cat(sprintf("  Goldfeld-Quandt F:   %12.*f\n", d, x$f_stat))
  cat(sprintf("  p-value:             %12s  %s\n", format_pval(x$p_value), get_stars(x$p_value)))
  print_line("-")
  if (x$p_value < 0.05) {
    verdict <- "REJECT H0 (p < 0.05) -> Significant heteroskedasticity detected."
    advice <- "Variance of residuals differs across values of the sorting variable; use `se = 'robust'`."
  } else {
    verdict <- "PASS H0 (p >= 0.05) -> No significant evidence of heteroskedasticity."
    advice <- "Error variances appear stable across the distribution."
  }
  cat(sprintf("  Verdict: %s\n", verdict))
  cat(sprintf("  Note:    %s\n", advice))
  print_line("=")
  invisible(x)
}

#' Davidson-MacKinnon J-Test for Non-Nested Models
#'
#' Tests competing non-nested econometric models (neither model is a subset of the other),
#' following the specification testing framework of Davidson and MacKinnon (1981).
#'
#' @param model1 Fitted model object or formula for Model 1.
#' @param model2 Fitted model object or formula for Model 2.
#' @param data A data.frame containing variables (required if formulas are passed).
#' @param digits Number of decimal places to format (default is 4).
#'
#' @return An object of class \code{econ_jtest}.
#' @export
#'
#' @examples
#' data(wage_data, package = "phoring")
#' # Model 1: Education & Experience vs. Model 2: Tenure & Mother's Education
#' m1 <- reg(wage ~ educ + exper, data = wage_data)
#' m2 <- reg(wage ~ tenure + mother_educ, data = wage_data)
#' j_test(m1, m2)
j_test <- function(model1, model2, data = NULL, digits = 4) {
  # Resolve formulas and data
  f1 <- if (inherits(model1, "econ_reg")) model1$formula else if (inherits(model1, "lm")) stats::formula(model1) else model1
  f2 <- if (inherits(model2, "econ_reg")) model2$formula else if (inherits(model2, "lm")) stats::formula(model2) else model2
  
  if (is.null(data)) {
    if ((inherits(model1, "econ_reg") || inherits(model1, "lm")) &&
        (inherits(model2, "econ_reg") || inherits(model2, "lm"))) {
      mf1 <- stats::model.frame(if (inherits(model1, "econ_reg")) model1$lm_model else model1)
      mf2 <- stats::model.frame(if (inherits(model2, "econ_reg")) model2$lm_model else model2)
      if (nrow(mf1) != nrow(mf2)) {
        stop("Sample sizes of Model 1 and Model 2 do not match.", call. = FALSE)
      }
      diff_cols <- setdiff(names(mf2), names(mf1))
      data <- cbind(mf1, mf2[, diff_cols, drop = FALSE])
    } else {
      stop("`data` must be supplied when formulas are passed to `j_test()`.", call. = FALSE)
    }
  }
  
  # Fit base models
  fit1 <- stats::lm(f1, data = data)
  fit2 <- stats::lm(f2, data = data)
  
  # Check same dependent variable
  y1 <- stats::model.response(stats::model.frame(fit1))
  y2 <- stats::model.response(stats::model.frame(fit2))
  if (length(y1) != length(y2) || max(abs(y1 - y2)) > 1e-6) {
    stop("Both models must have the exact same dependent variable and sample observations.", call. = FALSE)
  }
  
  yhat1 <- stats::fitted(fit1)
  yhat2 <- stats::fitted(fit2)
  
  # Direction 1: Test Model 1 by adding yhat2
  aux_data1 <- data
  aux_data1$yhat_m2 <- yhat2
  f1_aux <- stats::update(f1, . ~ . + yhat_m2)
  fit1_aux <- stats::lm(f1_aux, data = aux_data1)
  coefs1 <- summary(fit1_aux)$coefficients
  t1 <- coefs1["yhat_m2", "t value"]
  p1 <- coefs1["yhat_m2", "Pr(>|t|)"]
  
  # Direction 2: Test Model 2 by adding yhat1
  aux_data2 <- data
  aux_data2$yhat_m1 <- yhat1
  f2_aux <- stats::update(f2, . ~ . + yhat_m1)
  fit2_aux <- stats::lm(f2_aux, data = aux_data2)
  coefs2 <- summary(fit2_aux)$coefficients
  t2 <- coefs2["yhat_m1", "t value"]
  p2 <- coefs2["yhat_m1", "Pr(>|t|)"]
  
  # Verdict logic
  rej1 <- p1 < 0.05
  rej2 <- p2 < 0.05
  
  if (!rej1 && rej2) {
    decision <- "Model 1 is preferred (Model 2 is rejected, Model 1 is not rejected)."
  } else if (rej1 && !rej2) {
    decision <- "Model 2 is preferred (Model 1 is rejected, Model 2 is not rejected)."
  } else if (rej1 && rej2) {
    decision <- "Both models are rejected (neither specification is adequate; consider a broader model)."
  } else {
    decision <- "Neither model is rejected (the data cannot discriminate between the two specifications)."
  }
  
  res <- list(
    formula1 = f1,
    formula2 = f2,
    t1 = t1,
    p1 = p1,
    t2 = t2,
    p2 = p2,
    decision = decision,
    digits = digits
  )
  class(res) <- "econ_jtest"
  return(res)
}

#' @export
print.econ_jtest <- function(x, ...) {
  d <- x$digits
  print_line("=")
  cat("           DAVIDSON-MacKINNON J-TEST FOR NON-NESTED MODELS              \n")
  cat(paste0("  Model 1: ", deparse(x$formula1), "\n"))
  cat(paste0("  Model 2: ", deparse(x$formula2), "\n"))
  print_line("=")
  cat(sprintf("  %-25s %12s %12s %12s\n", "Hypothesis", "t-statistic", "p-value", "Verdict"))
  print_line("-")
  v1 <- if (x$p1 < 0.05) "REJECT M1" else "ACCEPT M1"
  v2 <- if (x$p2 < 0.05) "REJECT M2" else "ACCEPT M2"
  cat(sprintf("  %-25s %12.*f %12s %12s\n",
              "H0: Model 1 is correct", d, x$t1, format_pval(x$p1), v1))
  cat(sprintf("  %-25s %12.*f %12s %12s\n",
              "H0: Model 2 is correct", d, x$t2, format_pval(x$p2), v2))
  print_line("-")
  cat(sprintf("  Conclusion: %s\n", x$decision))
  print_line("=")
  invisible(x)
}

