#' Complete Time Series Econometric Toolkit
#'
#' Provides a comprehensive suite for time series econometrics:
#' 1. `unit_root()`: Multi-test stationarity battery (ADF, Phillips-Perron, KPSS).
#' 2. `coint()`: Engle-Granger two-step cointegration test.
#' 3. `granger()`: Bidirectional Granger causality test.
#' 4. `arch_test()`: Engle's ARCH LM test for volatility clustering.
#' 5. `ardl()`: Autoregressive Distributed Lag model with short-run and long-run multipliers.
#'
#' @name timeseries
NULL

#' Multi-Test Stationarity Battery (ADF, Phillips-Perron, KPSS)
#'
#' Executes Augmented Dickey-Fuller (ADF), Phillips-Perron (PP), and KPSS tests
#' simultaneously with MacKinnon response surface p-values.
#'
#' @param x A numeric vector or time series.
#' @param lags Number of lags for ADF test (optional; automatically selected by AIC if NULL).
#' @param type Model specification: \code{"trend"} (constant and linear trend, default),
#'   \code{"drift"} (constant only), or \code{"none"} (neither constant nor trend).
#' @param digits Number of decimals to display (default is 4).
#'
#' @return An object of class \code{econ_unit_root}.
#' @export
#'
#' @examples
#' data(macro_data, package = "phoring")
#' unit_root(macro_data$inflation, type = "trend")
#' unit_root(macro_data$unemployment, type = "drift")
unit_root <- function(x, lags = NULL, type = c("trend", "drift", "none"), digits = 4) {
  type <- match.arg(type)
  x <- stats::na.omit(as.numeric(x))
  n <- length(x)
  if (n < 15) stop("Time series must have at least 15 observations.", call. = FALSE)
  
  max_allowed_lags <- max(1, floor(4 * (n / 100)^(2 / 9)))
  
  dx <- diff(x)
  x_lag <- x[-n]
  n_adf <- length(dx)
  time_trend <- 1:n_adf
  
  # Automatic lag selection by AIC/BIC if lags is NULL or "aic"/"bic"/"auto"
  if (is.null(lags) || (is.character(lags) && tolower(lags[1]) %in% c("aic", "bic", "auto"))) {
    crit <- if (is.character(lags) && tolower(lags[1]) == "bic") "bic" else "aic"
    best_crit <- Inf
    best_lag <- 1
    for (cand_lag in 1:max_allowed_lags) {
      if (n_adf <= cand_lag + 4) break
      cand_ldx <- matrix(NA, nrow = n_adf, ncol = cand_lag)
      for (l in 1:cand_lag) cand_ldx[, l] <- c(rep(NA, l), dx[1:(n_adf - l)])
      c_cases <- (cand_lag + 1):n_adf
      
      if (type == "trend") {
        x_m <- cbind(1, time_trend[c_cases], x_lag[c_cases], cand_ldx[c_cases, , drop = FALSE])
      } else if (type == "drift") {
        x_m <- cbind(1, x_lag[c_cases], cand_ldx[c_cases, , drop = FALSE])
      } else {
        x_m <- cbind(x_lag[c_cases], cand_ldx[c_cases, , drop = FALSE])
      }
      fit_c <- stats::lm.fit(x = x_m, y = dx[c_cases])
      sse_c <- sum(fit_c$residuals^2)
      n_c <- length(c_cases)
      k_c <- ncol(x_m)
      val_c <- if (crit == "bic") n_c * log(sse_c / n_c) + k_c * log(n_c) else n_c * log(sse_c / n_c) + 2 * k_c
      if (!is.na(val_c) && val_c < best_crit) {
        best_crit <- val_c
        best_lag <- cand_lag
      }
    }
    lags <- best_lag
  } else {
    lags <- as.integer(lags)
  }
  
  # -------------------------------------------------------------
  # 1. Augmented Dickey-Fuller (ADF)
  # -------------------------------------------------------------
  if (lags > 0) {
    lag_dx <- matrix(NA, nrow = n_adf, ncol = lags)
    for (l in 1:lags) {
      lag_dx[, l] <- c(rep(NA, l), dx[1:(n_adf - l)])
    }
    complete_cases <- (lags + 1):n_adf
    y_reg <- dx[complete_cases]
    
    if (type == "trend") {
      x_reg <- cbind(1, time_trend[complete_cases], x_lag[complete_cases], lag_dx[complete_cases, , drop = FALSE])
      gamma_col <- 3
    } else if (type == "drift") {
      x_reg <- cbind(1, x_lag[complete_cases], lag_dx[complete_cases, , drop = FALSE])
      gamma_col <- 2
    } else {
      x_reg <- cbind(x_lag[complete_cases], lag_dx[complete_cases, , drop = FALSE])
      gamma_col <- 1
    }
  } else {
    y_reg <- dx
    if (type == "trend") {
      x_reg <- cbind(1, time_trend, x_lag)
      gamma_col <- 3
    } else if (type == "drift") {
      x_reg <- cbind(1, x_lag)
      gamma_col <- 2
    } else {
      x_reg <- matrix(x_lag, ncol = 1)
      gamma_col <- 1
    }
  }
  
  mod_adf <- stats::lm.fit(x = x_reg, y = y_reg)
  b_adf <- mod_adf$coefficients
  e_adf <- mod_adf$residuals
  s2_adf <- sum(e_adf^2) / mod_adf$df.residual
  se_gamma <- sqrt(s2_adf * safe_solve(crossprod(x_reg))[gamma_col, gamma_col])
  stat_adf <- b_adf[gamma_col] / se_gamma
  
  # Continuous MacKinnon response surface p-value
  p_adf <- calc_mackinnon_p(stat_adf, n = n, type = type, test = "unitroot")
  crit_5_adf <- if (type == "trend") -3.41 else if (type == "drift") -2.86 else -1.94
  verdict_adf <- if (p_adf < 0.05) "REJECT H0 (Stationary)" else "FAIL TO REJECT (Unit Root)"
  
  # -------------------------------------------------------------
  # 2. Phillips-Perron (PP test with Newey-West variance)
  # -------------------------------------------------------------
  if (type == "trend") {
    mod_pp <- stats::lm(dx ~ 1 + time_trend + x_lag)
    pp_gamma_col <- 3
  } else if (type == "drift") {
    mod_pp <- stats::lm(dx ~ 1 + x_lag)
    pp_gamma_col <- 2
  } else {
    mod_pp <- stats::lm(dx ~ 0 + x_lag)
    pp_gamma_col <- 1
  }
  
  e_pp <- stats::residuals(mod_pp)
  s2_e <- sum(e_pp^2) / length(e_pp)
  
  gamma_0 <- s2_e
  gamma_sum <- 0
  for (l in 1:lags) {
    w_l <- 1 - (l / (lags + 1))
    c_l <- sum(e_pp[(l + 1):length(e_pp)] * e_pp[1:(length(e_pp) - l)]) / length(e_pp)
    gamma_sum <- gamma_sum + 2 * w_l * c_l
  }
  lambda_sq <- max(1e-8, gamma_0 + gamma_sum)
  
  se_pp_gamma <- summary(mod_pp)$coefficients[pp_gamma_col, 2]
  t_pp <- summary(mod_pp)$coefficients[pp_gamma_col, 3]
  stat_pp <- sqrt(gamma_0 / lambda_sq) * t_pp - (0.5 * (lambda_sq - gamma_0) / sqrt(lambda_sq)) * (length(e_pp) * se_pp_gamma / sqrt(s2_e))
  p_pp <- calc_mackinnon_p(stat_pp, n = n, type = type, test = "unitroot")
  verdict_pp <- if (p_pp < 0.05) "REJECT H0 (Stationary)" else "FAIL TO REJECT (Unit Root)"
  
  # -------------------------------------------------------------
  # 3. KPSS Test
  # -------------------------------------------------------------
  if (type == "trend") {
    mod_kpss <- stats::lm(x ~ 1 + seq_along(x))
    kpss_null <- "Trend Stationary"
    crit_5_kpss <- 0.146
  } else {
    mod_kpss <- stats::lm(x ~ 1)
    kpss_null <- "Level Stationary"
    crit_5_kpss <- 0.463
  }
  
  e_kpss <- stats::residuals(mod_kpss)
  S_t <- cumsum(e_kpss)
  
  gamma_kpss_0 <- sum(e_kpss^2) / n
  gamma_kpss_sum <- 0
  for (l in 1:lags) {
    w_l <- 1 - (l / (lags + 1))
    c_l <- sum(e_kpss[(l + 1):n] * e_kpss[1:(n - l)]) / n
    gamma_kpss_sum <- gamma_kpss_sum + 2 * w_l * c_l
  }
  s_sq_kpss <- max(1e-8, gamma_kpss_0 + gamma_kpss_sum)
  stat_kpss <- sum(S_t^2) / (n^2 * s_sq_kpss)
  
  # Continuous smooth p-value for KPSS
  if (type == "trend") {
    knots_kpss <- c(0.05, 0.08, 0.119, 0.146, 0.176, 0.216, 0.35)
    knots_kpss_p <- c(0.90, 0.50, 0.10, 0.05, 0.025, 0.01, 0.001)
  } else {
    knots_kpss <- c(0.15, 0.25, 0.347, 0.463, 0.574, 0.739, 1.00)
    knots_kpss_p <- c(0.90, 0.50, 0.10, 0.05, 0.025, 0.01, 0.001)
  }
  p_kpss <- stats::approx(knots_kpss, knots_kpss_p, xout = stat_kpss, rule = 2)$y
  verdict_kpss <- if (stat_kpss < crit_5_kpss) "FAIL TO REJECT (Stationary)" else "REJECT H0 (Non-Stationary)"
  
  tests_df <- data.frame(
    Test = c(paste0("ADF (lags=", lags, ", ", type, ")"), paste0("Phillips-Perron (PP, ", type, ")"), paste0("KPSS (", type, ")")),
    Null = c("Unit Root (Non-stationary)", "Unit Root (Non-stationary)", kpss_null),
    Stat = c(stat_adf, stat_pp, stat_kpss),
    Crit_5pct = c(sprintf("%.3f", crit_5_adf), sprintf("%.3f", crit_5_adf), sprintf("%.3f", crit_5_kpss)),
    p_value = c(p_adf, p_pp, p_kpss),
    Result = c(verdict_adf, verdict_pp, verdict_kpss),
    stringsAsFactors = FALSE
  )
  
  is_stationary <- (p_adf < 0.05 && p_kpss >= 0.05)
  
  res <- list(
    series = x,
    tests = tests_df,
    is_stationary = is_stationary,
    n = n,
    lags = lags,
    type = type,
    digits = digits
  )
  class(res) <- "econ_unit_root"
  return(res)
}

#' @export
print.econ_unit_root <- function(x, digits = NULL, ...) {
  d <- if (is.null(digits)) x$digits else digits
  
  print_line("=")
  cat("           TIME SERIES STATIONARITY & UNIT ROOT BATTERY                 \n")
  cat(sprintf("  Observations (T): %-12d Selected Lags: %-4d Model: %s\n", x$n, x$lags, x$type))
  print_line("=")
  cat(sprintf("  %-24s %-26s %8s %10s %10s  %-20s\n", 
              "Test Name", "Null Hypothesis (H0)", "Stat", "5% Crit", "p-value", "Verdict"))
  print_line("-")
  
  for (i in 1:nrow(x$tests)) {
    row <- x$tests[i, ]
    cat(sprintf("  %-24s %-26s %8.*f %10s %10s  %-20s\n",
                row$Test,
                row$Null,
                d, row$Stat,
                row$Crit_5pct,
                format_pval(row$p_value),
                row$Result))
  }
  print_line("-")
  cat("  Econometric Synthesis:\n")
  if (x$is_stationary) {
    cat("  * The series appears STATIONARY [I(0)] at the 5% significance level.\n")
    cat("    --> Safe for standard OLS regression without differencing.\n")
  } else {
    cat("  * Evidence of NON-STATIONARITY [I(1)].\n")
    cat("    --> Running OLS in levels may lead to spurious regression.\n")
    cat("    --> Recommended: Take first difference `diff(y)` or test for Cointegration (`coint()`).\n")
  }
  print_line("=")
  invisible(x)
}

#' Engle-Granger Cointegration Test
#'
#' Tests whether two or more non-stationary time series share a stationary
#' long-run equilibrium cointegrating relationship.
#'
#' @param formula Model formula specifying cointegrating regression (e.g., \code{y ~ x} or \code{log(y) ~ log(x)}).
#' @param data A data.frame containing the time series variables.
#' @param lags Number of lags for residual ADF test (default is 1).
#' @param digits Number of decimal places to display (default is 4).
#'
#' @return An object of class \code{econ_coint}.
#' @export
#'
#' @examples
#' data(macro_data, package = "phoring")
#' coint(inflation ~ unemployment + fed_funds, data = macro_data)
coint <- function(formula, data, lags = 1, digits = 4) {
  if (missing(data) || is.null(data)) stop("Please provide data.", call. = FALSE)
  
  # Step 1: Static OLS with full formula handling
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(formula, data = mf)
  
  base_mod <- stats::lm.fit(x = X, y = y)
  b <- stats::coef(base_mod)
  e <- stats::residuals(base_mod)
  n <- length(e)
  k <- length(b)
  
  # Step 2: ADF test on residuals without constant (Engle-Granger)
  de <- diff(e)
  e_lag <- e[-n]
  n_de <- length(de)
  
  # Automatic lag selection if lags is NULL or "aic"/"bic"/"auto"
  if (is.null(lags) || (is.character(lags) && tolower(lags[1]) %in% c("aic", "bic", "auto"))) {
    crit <- if (is.character(lags) && tolower(lags[1]) == "bic") "bic" else "aic"
    max_coint_lags <- max(1, floor(4 * (n_de / 100)^(2 / 9)))
    best_crit <- Inf
    best_lag <- 1
    for (cand_l in 1:max_coint_lags) {
      if (n_de <= cand_l + 3) break
      cand_lde <- matrix(NA, nrow = n_de, ncol = cand_l)
      for (l in 1:cand_l) cand_lde[, l] <- c(rep(NA, l), de[1:(n_de - l)])
      c_cases <- (cand_l + 1):n_de
      fit_c <- stats::lm.fit(x = cbind(e_lag[c_cases], cand_lde[c_cases, , drop = FALSE]), y = de[c_cases])
      sse_c <- sum(fit_c$residuals^2)
      n_c <- length(c_cases)
      k_c <- 1 + cand_l
      val_c <- if (crit == "bic") n_c * log(sse_c / n_c) + k_c * log(n_c) else n_c * log(sse_c / n_c) + 2 * k_c
      if (!is.na(val_c) && val_c < best_crit) {
        best_crit <- val_c
        best_lag <- cand_l
      }
    }
    lags <- best_lag
  } else {
    lags <- as.integer(lags)
  }
  
  if (lags > 0) {
    lag_de <- matrix(NA, nrow = n_de, ncol = lags)
    for (l in 1:lags) lag_de[, l] <- c(rep(NA, l), de[1:(n_de - l)])
    idx <- (lags + 1):n_de
    mod_adf_e <- stats::lm(de[idx] ~ 0 + e_lag[idx] + lag_de[idx, , drop = FALSE])
  } else {
    mod_adf_e <- stats::lm(de ~ 0 + e_lag)
  }
  
  stat_tau <- summary(mod_adf_e)$coefficients[1, 3]
  
  # MacKinnon cointegration critical values and continuous p-value
  crit_5pct <- -3.34 - 0.38 * (k - 2)
  p_val <- calc_mackinnon_p(stat_tau, n = n, test = "coint", k = k)
  is_cointegrated <- (p_val < 0.05)
  
  res <- list(
    formula = formula,
    coefficients = b,
    static_model = base_mod,
    tau_stat = stat_tau,
    crit_5pct = crit_5pct,
    p_value = p_val,
    is_cointegrated = is_cointegrated,
    residuals = e,
    fitted.values = as.vector(X %*% b),
    y = y,
    X = X,
    n = n,
    k = k,
    digits = digits
  )
  class(res) <- "econ_coint"
  return(res)
}

#' @export
print.econ_coint <- function(x, digits = NULL, ...) {
  d <- if (is.null(digits)) x$digits else digits
  
  print_line("=")
  cat("           ENGLE-GRANGER TWO-STEP COINTEGRATION TEST                    \n")
  cat(paste0("  Formula: ", deparse(x$formula), "\n"))
  cat(sprintf("  Observations (T): %-12d Regressors: %d\n", x$n, x$k - 1))
  print_line("=")
  
  cat("  1. ESTIMATED LONG-RUN EQUILIBRIUM EQUATION:\n")
  print_line("-")
  b <- x$coefficients
  for (i in seq_along(b)) {
    cat(sprintf("    %-16s : %10.*f\n", names(b)[i], d, b[i]))
  }
  print_line("-")
  
  cat("  2. RESIDUAL UNIT ROOT TEST (COINTEGRATION TEST):\n")
  print_line("-")
  cat(sprintf("  Test Statistic (Tau):  %10.*f\n", d, x$tau_stat))
  cat(sprintf("  MacKinnon 5%% Critical: %10.*f\n", d, x$crit_5pct))
  cat(sprintf("  MacKinnon p-value:     %10s\n", format_pval(x$p_value)))
  print_line("-")
  
  cat("  Econometric Verdict:\n")
  if (x$is_cointegrated) {
    cat("  * Tau is more negative than 5% critical value (p < 0.05).\n")
    cat("    --> REJECT H0 (No cointegration) -> The series ARE COINTEGRATED.\n")
    cat("    --> A stable, non-spurious long-run relationship exists!\n")
    cat("    --> You can estimate an Error Correction Model (ECM) using these residuals.\n")
  } else {
    cat("  * Tau does not exceed 5% critical value (p >= 0.05).\n")
    cat("    --> FAIL TO REJECT H0 -> NO cointegration detected.\n")
    cat("    --> Running regression in levels would yield spurious results.\n")
  }
  print_line("=")
  invisible(x)
}

#' @export
coef.econ_coint <- function(object, ...) {
  object$coefficients
}

#' @export
residuals.econ_coint <- function(object, ...) {
  object$residuals
}

#' @export
fitted.econ_coint <- function(object, ...) {
  object$fitted.values
}

#' Bidirectional Granger Causality Test
#'
#' Tests whether past values of X forecast Y, and whether past values of Y forecast X.
#'
#' @param formula Model formula specifying bivariate relationship \code{y ~ x}.
#' @param data A data.frame containing the time series.
#' @param lags Number of autoregressive lags to include (default is 2).
#' @param digits Number of decimals to display (default is 4).
#'
#' @return An object of class \code{econ_granger}.
#' @export
#'
#' @examples
#' data(macro_data, package = "phoring")
#' granger(inflation ~ fed_funds, data = macro_data, lags = 2)
granger <- function(formula, data, lags = 2, digits = 4) {
  if (missing(data) || is.null(data)) stop("Please provide data.", call. = FALSE)
  
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  x_name <- names(mf)[2]
  x <- mf[[x_name]]
  y_name <- names(mf)[1]
  n <- length(y)
  
  # Function to run Granger regression
  run_granger_dir <- function(dep, indep) {
    n_t <- length(dep)
    lag_dep <- embed(dep, lags + 1)[, -1, drop = FALSE]
    lag_indep <- embed(indep, lags + 1)[, -1, drop = FALSE]
    y_target <- dep[(lags + 1):n_t]
    
    mod_unres <- stats::lm(y_target ~ lag_dep + lag_indep)
    mod_res <- stats::lm(y_target ~ lag_dep)
    
    an <- stats::anova(mod_res, mod_unres)
    f_stat <- an$F[2]
    p_val <- an$`Pr(>F)`[2]
    list(f_stat = f_stat, p_val = p_val)
  }
  
  res_xy <- run_granger_dir(y, x) # X -> Y
  res_yx <- run_granger_dir(x, y) # Y -> X
  
  res <- list(
    y_name = y_name,
    x_name = x_name,
    x_series = x,
    y_series = y,
    lags = lags,
    n = n - lags,
    f_xy = res_xy$f_stat,
    p_xy = res_xy$p_val,
    f_yx = res_yx$f_stat,
    p_yx = res_yx$p_val,
    digits = digits
  )
  class(res) <- "econ_granger"
  return(res)
}

#' @export
print.econ_granger <- function(x, digits = NULL, ...) {
  d <- if (is.null(digits)) x$digits else digits
  
  print_line("=")
  cat("                BIDIRECTIONAL GRANGER CAUSALITY TEST                    \n")
  cat(sprintf("  Variables: %s <---> %s (Lags = %d, T = %d)\n", 
              x$x_name, x$y_name, x$lags, x$n))
  print_line("=")
  cat(sprintf("  %-32s %10s %10s  %-20s\n", 
              "Null Hypothesis (H0)", "F-stat", "p-value", "Conclusion"))
  print_line("-")
  
  verdict_xy <- if (x$p_xy < 0.05) "REJECT (X causes Y)" else "FAIL TO REJECT (No causality)"
  cat(sprintf("  %-32s %10.*f %10s  %-20s\n",
              paste0(x$x_name, " does not Granger-cause ", x$y_name),
              d, x$f_xy,
              format_pval(x$p_xy),
              verdict_xy))
  
  verdict_yx <- if (x$p_yx < 0.05) "REJECT (Y causes X)" else "FAIL TO REJECT (No causality)"
  cat(sprintf("  %-32s %10.*f %10s  %-20s\n",
              paste0(x$y_name, " does not Granger-cause ", x$x_name),
              d, x$f_yx,
              format_pval(x$p_yx),
              verdict_yx))
  print_line("-")
  
  cat("  Causal Direction Summary:\n")
  if (x$p_xy < 0.05 && x$p_yx < 0.05) {
    cat(sprintf("  * BIDIRECTIONAL FEEDBACK detected between %s and %s.\n", x$x_name, x$y_name))
  } else if (x$p_xy < 0.05) {
    cat(sprintf("  * UNIDIRECTIONAL: %s Granger-causes %s.\n", x$x_name, x$y_name))
  } else if (x$p_yx < 0.05) {
    cat(sprintf("  * UNIDIRECTIONAL: %s Granger-causes %s.\n", x$y_name, x$x_name))
  } else {
    cat(sprintf("  * NO CAUSALITY detected between %s and %s.\n", x$x_name, x$y_name))
  }
  print_line("=")
  invisible(x)
}

#' Autoregressive Conditional Heteroskedasticity (ARCH) Test
#'
#' Engle's Lagrange Multiplier test for ARCH effects / volatility clustering in time series models.
#'
#' @param model A fitted model object from \code{\link{reg}} or \code{\link[stats]{lm}}.
#' @param lags Number of squared residual lags to test (default is 2).
#' @param digits Number of decimals to display (default is 4).
#'
#' @return An object of class \code{econ_arch}.
#' @export
#'
#' @examples
#' data(macro_data, package = "phoring")
#' m <- reg(inflation ~ fed_funds, data = macro_data, diagnostics = FALSE)
#' arch_test(m, lags = 2)
arch_test <- function(model, lags = 2, digits = 4) {
  e <- stats::residuals(model)
  n <- length(e)
  e2 <- e^2
  
  lag_e2 <- embed(e2, lags + 1)[, -1, drop = FALSE]
  y_e2 <- e2[(lags + 1):n]
  
  aux_mod <- stats::lm(y_e2 ~ lag_e2)
  r2 <- summary(aux_mod)$r.squared
  lm_stat <- (n - lags) * r2
  p_val <- stats::pchisq(lm_stat, df = lags, lower.tail = FALSE)
  
  res <- list(
    residuals = e,
    stat = lm_stat,
    df = lags,
    p_value = p_val,
    digits = digits
  )
  class(res) <- "econ_arch"
  return(res)
}

#' @export
print.econ_arch <- function(x, digits = NULL, ...) {
  d <- if (is.null(digits)) x$digits else digits
  
  print_line("=")
  cat("              ENGLE'S ARCH LM TEST FOR VOLATILITY CLUSTERING            \n")
  print_line("=")
  cat(sprintf("  H0: No autoregressive conditional heteroskedasticity (No ARCH)\n"))
  cat(sprintf("  LM Statistic:  Chi-Sq = %10.*f  (df = %d)\n", d, x$stat, x$df))
  cat(sprintf("  p-value:       %10s\n", format_pval(x$p_value)))
  print_line("-")
  verdict <- if (x$p_value < 0.05) {
    "REJECT H0 (p < 0.05) -> Significant ARCH effects / Volatility clustering detected"
  } else {
    "FAIL TO REJECT H0 (p >= 0.05) -> No ARCH effects detected"
  }
  cat(sprintf("  Conclusion:    %s\n", verdict))
  print_line("=")
  invisible(x)
}

#' Autoregressive Distributed Lag (ARDL) Model
#'
#' Estimates an ARDL(p, q) model and calculates both short-run impact multipliers
#' and long-run cumulative propensity multipliers with delta-method standard errors.
#'
#' @param formula Model formula \code{y ~ x}.
#' @param data A data.frame containing the time series.
#' @param p Number of autoregressive lags of y (default is 1).
#' @param q Number of distributed lags of x (default is 1).
#' @param digits Number of decimals to display (default is 4).
#'
#' @return An object of class \code{econ_ardl}.
#' @export
#'
#' @examples
#' data(macro_data, package = "phoring")
#' ardl(inflation ~ unemployment, data = macro_data, p = 1, q = 1)
ardl <- function(formula, data, p = 1, q = 1, digits = 4) {
  if (missing(data) || is.null(data)) stop("Please provide data.", call. = FALSE)
  
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  x_name <- names(mf)[2]
  x <- mf[[x_name]]
  y_name <- names(mf)[1]
  n <- length(y)
  
  max_lag <- max(p, q)
  lag_y <- matrix(NA, nrow = n, ncol = p)
  for (i in 1:p) lag_y[, i] <- c(rep(NA, i), y[1:(n - i)])
  colnames(lag_y) <- paste0(y_name, "_lag", 1:p)
  
  lag_x <- matrix(NA, nrow = n, ncol = q + 1)
  lag_x[, 1] <- x
  colnames(lag_x) <- c(x_name, paste0(x_name, "_lag", 1:q))
  if (q > 0) {
    for (j in 1:q) lag_x[, j + 1] <- c(rep(NA, j), x[1:(n - j)])
  }
  
  df_ardl <- stats::na.omit(data.frame(y = y, lag_y, lag_x))
  mod <- stats::lm(y ~ ., data = df_ardl)
  
  b <- stats::coef(mod)
  V <- stats::vcov(mod)
  
  impact_mult <- b[x_name]
  se_impact <- sqrt(V[x_name, x_name])
  
  x_cols <- colnames(lag_x)
  y_cols <- colnames(lag_y)
  sum_beta <- sum(b[x_cols])
  sum_alpha <- sum(b[y_cols])
  denom <- 1 - sum_alpha
  long_run_mult <- sum_beta / denom
  
  grad <- numeric(length(b))
  names(grad) <- names(b)
  grad[x_cols] <- 1 / denom
  grad[y_cols] <- sum_beta / (denom^2)
  
  var_lr <- as.numeric(t(grad) %*% V %*% grad)
  se_lr <- sqrt(max(0, var_lr))
  t_lr <- long_run_mult / se_lr
  p_lr <- 2 * stats::pt(-abs(t_lr), df = stats::df.residual(mod))
  
  # Pesaran, Shin & Smith (2001) Bounds F-statistic on lagged level variables
  R_mat <- matrix(0, nrow = 2, ncol = length(b))
  colnames(R_mat) <- names(b)
  R_mat[1, y_cols] <- 1
  R_mat[2, x_cols] <- 1
  q_vec <- c(1, 0)
  r_vec <- as.vector(R_mat %*% b - q_vec)
  V_r <- R_mat %*% V %*% t(R_mat)
  inv_V_r <- tryCatch(solve(V_r), error = function(e) NULL)
  f_bounds <- if (!is.null(inv_V_r)) {
    max(0, as.numeric(t(r_vec) %*% inv_V_r %*% r_vec) / 2)
  } else {
    as.numeric(t_lr^2)
  }

  res <- list(
    formula = formula,
    model = mod,
    coefficients = b,
    vcov = V,
    residuals = stats::residuals(mod),
    fitted.values = stats::fitted(mod),
    n = nrow(df_ardl),
    p = p,
    q = q,
    impact_mult = impact_mult,
    se_impact = se_impact,
    long_run_mult = long_run_mult,
    se_lr = se_lr,
    t_lr = t_lr,
    p_lr = p_lr,
    f_bounds = f_bounds,
    digits = digits
  )
  class(res) <- "econ_ardl"
  return(res)
}

#' @export
print.econ_ardl <- function(x, digits = NULL, ...) {
  d <- if (is.null(digits)) x$digits else digits
  
  print_line("=")
  cat(sprintf("        AUTOREGRESSIVE DISTRIBUTED LAG MODEL: ARDL(%d, %d)       \n", x$p, x$q))
  cat(paste0("  Formula: ", deparse(x$formula), "\n"))
  print_line("=")
  
  cat("  1. REGRESSION ESTIMATES:\n")
  print(summary(x$model)$coefficients)
  print_line("-")
  
  cat("  2. MULTIPLIER ANALYSIS (SHORT-RUN VS LONG-RUN DYNAMICS):\n")
  print_line("-")
  cat(sprintf("  %-25s : %10.*f  (SE = %.*f)\n", 
              "Short-Run Impact Multiplier", d, x$impact_mult, d, x$se_impact))
  cat(sprintf("  %-25s : %10.*f  (SE = %.*f, t = %.3f, p = %s) %s\n", 
              "Long-Run Cumulative Propensity", d, x$long_run_mult, d, x$se_lr, x$t_lr,
              format_pval(x$p_lr), get_stars(x$p_lr)))
  print_line("=")
  invisible(x)
}

#' @export
coef.econ_ardl <- function(object, ...) {
  stats::coef(object$model)
}

#' @export
vcov.econ_ardl <- function(object, ...) {
  stats::vcov(object$model)
}

#' @export
residuals.econ_ardl <- function(object, ...) {
  stats::residuals(object$model)
}

#' @export
fitted.econ_ardl <- function(object, ...) {
  stats::fitted(object$model)
}

#' @export
nobs.econ_ardl <- function(object, ...) {
  object$n
}

#' @export
confint.econ_ardl <- function(object, parm, level = 0.95, ...) {
  stats::confint(object$model, parm = parm, level = level, ...)
}

#' @export
predict.econ_ardl <- function(object, newdata = NULL, ...) {
  stats::predict(object$model, newdata = newdata, ...)
}
