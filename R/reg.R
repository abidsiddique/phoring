#' One-Command Econometric Regression with Automated Diagnostics
#'
#' Fits an Ordinary Least Squares (OLS) or Weighted Least Squares (WLS) regression model,
#' calculates robust (HC0-HC3), clustered, or HAC standard errors, and automatically
#' runs and prints a complete battery of econometric diagnostic tests.
#'
#' @param formula An object of class \code{\link[stats]{formula}} (e.g., \code{y ~ x1 + x2}).
#' @param data A data.frame containing the variables in the model.
#' @param se Type of standard errors to compute: \code{"robust"} (White HC1, default),
#'   \code{"standard"} (classical OLS), \code{"cluster"} (cluster-robust),
#'   \code{"hac"} (Newey-West autocorrelation-consistent), or specific types
#'   \code{"HC0"}, \code{"HC1"}, \code{"HC2"}, \code{"HC3"}.
#' @param cluster Vector or column name indicating cluster grouping for clustered standard errors.
#' @param lag Number of lags for Newey-West HAC standard errors (optional, defaults to automatic rule).
#' @param weights Optional weights vector or column name for Weighted Least Squares (WLS).
#' @param digits Number of decimal places to display in output tables (default is 4).
#' @param diagnostics Logical: whether to automatically compute and display diagnostic tests (default is TRUE).
#' @param ... Additional arguments passed to \code{\link[stats]{lm}}.
#'
#' @return An object of class \code{econ_reg} containing the regression results,
#'   standard errors, coefficient table, model statistics, and diagnostic tests.
#' @export
#'
#' @examples
#' data(mtcars)
#' reg(mpg ~ wt + hp + qsec, data = mtcars)
#' reg(mpg ~ wt + hp, data = mtcars, cluster = "cyl")
reg <- function(formula, data, se = c("robust", "standard", "cluster", "hac", "HC0", "HC1", "HC2", "HC3"),
                cluster = NULL, lag = NULL, weights = NULL, digits = 4, diagnostics = TRUE, ...) {
  
  # Handle se argument default
  if (!missing(cluster) && !is.null(cluster)) {
    se <- "cluster"
  } else if (missing(se)) {
    se <- "robust"
  } else {
    se <- match.arg(se)
  }
  
  if (missing(data) || is.null(data)) {
    stop("Please provide a dataset in `data`.", call. = FALSE)
  }
  
  # Check if cluster is passed as a column name
  cluster_vec <- NULL
  if (!is.null(cluster)) {
    if (is.character(cluster) && length(cluster) == 1 && cluster %in% names(data)) {
      cluster_vec <- data[[cluster]]
    } else if (length(cluster) == nrow(data)) {
      cluster_vec <- cluster
    } else {
      stop("`cluster` must be a column name in `data` or a vector matching `nrow(data)`.", call. = FALSE)
    }
  }
  
  # Check if weights is passed as a column name
  weights_vec <- NULL
  if (!is.null(weights)) {
    if (is.character(weights) && length(weights) == 1 && weights %in% names(data)) {
      weights_vec <- data[[weights]]
    } else if (length(weights) == nrow(data)) {
      weights_vec <- weights
    } else {
      stop("`weights` must be a column name in `data` or a vector matching `nrow(data)`.", call. = FALSE)
    }
  }
  
  # Fit base OLS or WLS model
  if (!is.null(weights_vec)) {
    w_num <- as.numeric(weights_vec)
    cl <- bquote(stats::lm(.(formula), data = data, weights = .(w_num), x = TRUE, y = TRUE, model = TRUE, ...))
    base_lm <- eval(cl)
  } else {
    base_lm <- stats::lm(formula, data = data, x = TRUE, y = TRUE, model = TRUE, ...)
  }
  
  # Extract components
  X <- base_lm$x
  y <- base_lm$y
  e <- stats::residuals(base_lm)
  n <- nrow(X)
  k <- ncol(X)
  df_resid <- n - k
  
  if (df_resid <= 0) {
    stop("Degrees of freedom <= 0. Not enough observations relative to number of predictors.", call. = FALSE)
  }
  
  # Check for multicollinearity / dummy variable trap
  check_collinearity(X)
  
  # Handle subsetted cluster / weights vector if missing rows were dropped by lm
  na_action <- stats::na.action(base_lm)
  if (!is.null(na_action)) {
    if (!is.null(cluster_vec)) cluster_vec <- cluster_vec[-na_action]
    if (!is.null(weights_vec)) weights_vec <- weights_vec[-na_action]
  }
  
  # Compute Variance-Covariance Matrix
  df_inference <- df_resid
  
  if (se %in% c("robust", "HC1")) {
    vcov_mat <- calc_vcov_hc(X, e, type = "HC1", weights = weights_vec)
    se_label <- "Robust (White HC1)"
  } else if (se == "HC0") {
    vcov_mat <- calc_vcov_hc(X, e, type = "HC0", weights = weights_vec)
    se_label <- "Robust (White HC0)"
  } else if (se == "HC2") {
    vcov_mat <- calc_vcov_hc(X, e, type = "HC2", weights = weights_vec)
    se_label <- "Robust (MacKinnon-White HC2)"
  } else if (se == "HC3") {
    vcov_mat <- calc_vcov_hc(X, e, type = "HC3", weights = weights_vec)
    se_label <- "Robust (Davidson-MacKinnon HC3)"
  } else if (se == "cluster") {
    if (is.null(cluster_vec)) {
      stop("Please specify the `cluster` argument when `se = 'cluster'`.", call. = FALSE)
    }
    cluster_fac <- as.factor(cluster_vec)
    G <- nlevels(cluster_fac)
    df_inference <- max(1, G - 1)
    vcov_mat <- calc_vcov_cluster(X, e, cluster = cluster_vec, weights = weights_vec)
    se_label <- paste0("Clustered (by ", if (is.character(cluster)) cluster else "group", ", G = ", G, ")")
  } else if (se == "hac") {
    vcov_mat <- calc_vcov_hac(X, e, lag = lag, weights = weights_vec)
    se_label <- "Newey-West HAC"
  } else {
    if (!is.null(weights_vec)) {
      s2 <- sum(weights_vec * e^2) / df_resid
      bread <- safe_solve(crossprod(X * sqrt(weights_vec)))
    } else {
      s2 <- sum(e^2) / df_resid
      bread <- safe_solve(crossprod(X))
    }
    vcov_mat <- s2 * bread
    dimnames(vcov_mat) <- list(colnames(X), colnames(X))
    se_label <- if (is.null(weights_vec)) "Standard (Classical OLS)" else "Standard (Classical WLS)"
  }
  
  # Coefficients & Inference Table
  b <- stats::coef(base_lm)
  se_b <- sqrt(pmax(0, diag(vcov_mat)))
  t_b <- b / se_b
  p_b <- 2 * stats::pt(-abs(t_b), df = df_inference)
  
  coef_table <- data.frame(
    Estimate   = b,
    Std_Error  = se_b,
    t_stat     = t_b,
    p_value    = p_b,
    Signif     = sapply(p_b, get_stars),
    stringsAsFactors = FALSE
  )
  
  # Model Summary Statistics
  r2 <- summary(base_lm)$r.squared
  adj_r2 <- summary(base_lm)$adj.r.squared
  
  # Overall F-statistic (Wald test matching SE type)
  if (se == "standard") {
    f_stat <- summary(base_lm)$fstatistic
    f_val <- if (!is.null(f_stat)) f_stat[1] else NA
    f_p <- if (!is.null(f_stat)) stats::pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE) else NA
  } else {
    rob_f <- calc_robust_f(b, vcov_mat, df_denom = df_inference, has_intercept = "(Intercept)" %in% names(b))
    f_val <- rob_f$f_stat
    f_p <- rob_f$p_value
  }
  
  sigma_e <- sqrt(sum(if (!is.null(weights_vec)) weights_vec * e^2 else e^2) / df_resid)
  
  # Run automated diagnostics battery
  diag_res <- NULL
  if (isTRUE(diagnostics)) {
    diag_res <- tryCatch({
      diagnostics(base_lm)
    }, error = function(err) {
      NULL
    })
  }
  
  res <- list(
    formula = formula,
    coefficients = b,
    coef_table = coef_table,
    vcov = vcov_mat,
    residuals = e,
    fitted.values = stats::fitted(base_lm),
    n = n,
    k = k,
    df_residual = df_resid,
    df_inference = df_inference,
    r_squared = r2,
    adj_r_squared = adj_r2,
    sigma = sigma_e,
    f_statistic = f_val,
    f_p_value = f_p,
    se_type = se_label,
    weights = weights_vec,
    digits = digits,
    diagnostics = diag_res,
    lm_model = base_lm
  )
  class(res) <- "econ_reg"
  return(res)
}

#' @export
print.econ_reg <- function(x, digits = NULL, ...) {
  d <- if (is.null(digits)) x$digits else digits
  
  print_line("=")
  cat("                      OLS REGRESSION RESULTS                            \n")
  cat(paste0("  Formula: ", deparse(x$formula), "\n"))
  cat(sprintf("  Sample Size: %-13d R-squared:    %-10.4f F-stat:  %.2f\n", 
              x$n, x$r_squared, x$f_statistic))
  cat(sprintf("  Residual SE: %-13.4f Adj. R-sq:    %-10.4f p-value: %s\n", 
              x$sigma, x$adj_r_squared, format_pval(x$f_p_value)))
  cat(paste0("  Standard Errors: ", x$se_type, "\n"))
  v_w <- max(16, max(nchar(rownames(x$coef_table))))
  line_w <- max(74, v_w + 58)
  print_line("=", width = line_w)
  
  cat(sprintf(paste0("  %-", v_w, "s %10s %12s %10s %10s   %s\n"), 
              "Variable", "Estimate", "Std. Error", "t-stat", "p-value", "Signif"))
  print_line("-", width = line_w)
  
  for (i in 1:nrow(x$coef_table)) {
    vname <- rownames(x$coef_table)[i]
    row <- x$coef_table[i, ]
    cat(sprintf(paste0("  %-", v_w, "s %10.*f %12.*f %10.3f %10s   %s\n"),
                vname,
                d, row$Estimate,
                d, row$Std_Error,
                row$t_stat,
                format_pval(row$p_value),
                row$Signif))
  }
  print_line("-", width = line_w)
  cat("  Signif. codes:  *** p < 0.001,  ** p < 0.01,  * p < 0.05,  . p < 0.1\n")
  
  # Print Diagnostic Card if available
  if (!is.null(x$diagnostics)) {
    cat("\n")
    print_line("=")
    cat("                   ALL-IN-ONE DIAGNOSTIC TESTS                          \n")
    print_line("=")
    cat(sprintf("  %-23s %-24s %9s %10s  %-20s\n", "Test Name", "Null Hypothesis (H0)", "Stat", "p-value", "Result"))
    print_line("-")
    
    tests_df <- x$diagnostics$tests
    for (i in 1:nrow(tests_df)) {
      row <- tests_df[i, ]
      cat(sprintf("  %-23s %-24s %9.3f %10s  %-20s\n",
                  row$Test,
                  row$Null,
                  row$Stat,
                  format_pval(row$p_value),
                  row$Verdict))
    }
    print_line("-")
    
    # Print VIFs
    if (length(x$diagnostics$vifs) > 0) {
      cat("  Multicollinearity (VIF):\n    ")
      vif_strs <- paste0(names(x$diagnostics$vifs), ": ", sprintf("%.2f", x$diagnostics$vifs))
      cat(paste(vif_strs, collapse = "  |  "))
      max_vif <- max(x$diagnostics$vifs, na.rm = TRUE)
      if (max_vif < 5) {
        cat("\n    --> Assessment: Safe (all VIF < 5)\n")
      } else if (max_vif < 10) {
        cat("\n    --> Assessment: Moderate multicollinearity detected (5 <= VIF < 10)\n")
      } else {
        cat("\n    --> Assessment: Severe multicollinearity detected (VIF >= 10)\n")
      }
    }
    
    # Econometric guidance notes
    print_line("-")
    cat("  Econometric Guidance Notes:\n")
    bp_row <- tests_df[tests_df$Test == "Breusch-Pagan", ]
    if (nrow(bp_row) > 0 && bp_row$p_value < 0.05) {
      cat("  * Heteroskedasticity detected (Breusch-Pagan p < 0.05).\n")
      if (grepl("Robust", x$se_type)) {
        cat("    --> Robust standard errors are properly applied above.\n")
      } else {
        cat("    --> WARNING: You are using classical SEs. Consider using default robust SEs.\n")
      }
    } else {
      cat("  * No significant evidence of heteroskedasticity (Breusch-Pagan p >= 0.05).\n")
    }
    
    dw_row <- tests_df[tests_df$Test == "Durbin-Watson", ]
    if (nrow(dw_row) > 0 && dw_row$p_value < 0.05) {
      cat("  * Serial correlation detected (Durbin-Watson p < 0.05). Note: assumes time series.\n")
    }
    
    reset_row <- tests_df[tests_df$Test == "Ramsey RESET", ]
    if (nrow(reset_row) > 0 && reset_row$p_value < 0.05) {
      cat("  * Functional form misspecification detected (Ramsey RESET p < 0.05).\n")
      cat("    --> Consider adding polynomial terms (e.g. exper^2) or log transformations.\n")
    }
    print_line("=")
  }
  invisible(x)
}

#' @importFrom stats coef vcov residuals fitted nobs confint predict logLik AIC BIC
#' @export
summary.econ_reg <- function(object, ...) {
  object
}

#' @export
coef.econ_reg <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.econ_reg <- function(object, ...) {
  object$vcov
}

#' @export
residuals.econ_reg <- function(object, ...) {
  object$residuals
}

#' @export
fitted.econ_reg <- function(object, ...) {
  object$fitted.values
}

#' @export
nobs.econ_reg <- function(object, ...) {
  object$n
}

#' @export
confint.econ_reg <- function(object, parm, level = 0.95, ...) {
  b <- stats::coef(object)
  se <- sqrt(pmax(0, diag(stats::vcov(object))))
  df_inf <- if (!is.null(object$df_inference)) object$df_inference else object$df_residual
  a <- (1 - level) / 2
  crit <- stats::qt(1 - a, df = df_inf)
  ci <- cbind(b - crit * se, b + crit * se)
  colnames(ci) <- c(sprintf("%.1f %%", 100 * a), sprintf("%.1f %%", 100 * (1 - a)))
  rownames(ci) <- names(b)
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

#' @export
predict.econ_reg <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    return(object$fitted.values)
  }
  stats::predict(object$lm_model, newdata = newdata, ...)
}

#' @export
logLik.econ_reg <- function(object, ...) {
  stats::logLik(object$lm_model, ...)
}

#' @export
AIC.econ_reg <- function(object, ..., k = 2) {
  stats::AIC(object$lm_model, ..., k = k)
}

#' @export
BIC.econ_reg <- function(object, ...) {
  stats::BIC(object$lm_model, ...)
}

