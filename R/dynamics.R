#' Dynamic Econometric Models: Koyck, Durbin's h, Prais-Winsten, ECM, and VAR
#'
#' Implements dynamic time series models, autoregressive distributed lags, 
#' serial correlation corrections, error correction models, and vector autoregressions.
#'
#' @name dynamics
NULL

#' Durbin's h Test for Autocorrelation in Autoregressive Models
#'
#' Tests for first-order autocorrelation in dynamic models containing a lagged 
#' dependent variable, where the standard Durbin-Watson d test is biased toward 2 (Durbin, 1970).
#'
#' @param model A fitted \code{econ_reg} or \code{lm} model object containing a lagged dependent variable.
#' @param lag_var Optional character name of the lagged dependent variable. If \code{NULL},
#'   auto-detected by looking for terms containing \code{"lag"} or matching the dependent variable name.
#' @param digits Number of decimal places to format (default is 4).
#'
#' @return An object of class \code{econ_durbin_h}.
#' @export
#'
#' @examples
#' data(macro_data, package = "phoring")
#' n <- nrow(macro_data)
#' macro_data$inf_lag <- c(NA, macro_data$inflation[-n])
#' m <- reg(inflation ~ inf_lag + unemployment, data = macro_data[-1, ])
#' durbin_h(m, lag_var = "inf_lag")
durbin_h <- function(model, lag_var = NULL, digits = 4) {
  if (inherits(model, "econ_reg")) {
    base_lm <- model$lm_model
    e <- model$residuals
    f <- model$formula
  } else if (inherits(model, "lm")) {
    base_lm <- model
    e <- stats::residuals(model)
    f <- stats::formula(model)
  } else {
    stop("`model` must be an `econ_reg` or `lm` object.", call. = FALSE)
  }
  
  n <- length(e)
  dep_var <- deparse(f[[2]])
  all_preds <- names(stats::coef(base_lm))
  
  # Auto-detect lagged dependent variable if not supplied
  if (is.null(lag_var)) {
    candidates <- all_preds[grepl("lag", all_preds, ignore.case = TRUE) | 
                            grepl(dep_var, all_preds, fixed = TRUE)]
    candidates <- candidates[candidates != dep_var]
    if (length(candidates) == 0) {
      stop("Could not auto-detect lagged dependent variable. Please specify `lag_var = 'var_name'`.", call. = FALSE)
    }
    lag_var <- candidates[1]
  }
  
  if (!lag_var %in% all_preds) {
    stop(sprintf("Lagged variable `%s` not found in model predictors.", lag_var), call. = FALSE)
  }
  
  # Variance of the coefficient on the lagged dependent variable
  V <- stats::vcov(base_lm)
  var_alpha <- V[lag_var, lag_var]
  
  # Estimate first-order autocorrelation rho from residuals: rho = sum(e_t * e_{t-1}) / sum(e_{t-1}^2)
  e_curr <- e[-1]
  e_lag  <- e[-n]
  rho_hat <- sum(e_curr * e_lag) / sum(e_lag^2)
  
  # Standard Durbin-Watson d
  d_stat <- sum((e_curr - e_lag)^2) / sum(e^2)
  
  # Check condition: 1 - n * Var(alpha) > 0
  denom_term <- 1 - (n * var_alpha)
  
  if (denom_term > 0) {
    h_stat <- rho_hat * sqrt(n / denom_term)
    p_val  <- 2 * stats::pnorm(-abs(h_stat))
    is_valid <- TRUE
  } else {
    # Durbin's alternative test: regress e_t on e_{t-1} and all model X's
    mf <- stats::model.frame(base_lm)
    aux_data <- as.data.frame(stats::model.matrix(base_lm))
    aux_data$e_curr <- c(NA, e[-1])
    aux_data$e_lag  <- c(NA, e[-n])
    aux_fit <- stats::lm(e_curr ~ . - 1 + e_lag, data = aux_data[-1, ])
    aux_coef <- summary(aux_fit)$coefficients
    t_alt <- aux_coef["e_lag", "t value"]
    p_alt <- aux_coef["e_lag", "Pr(>|t|)"]
    
    h_stat <- t_alt
    p_val  <- p_alt
    is_valid <- FALSE
  }
  
  res <- list(
    formula = f,
    lag_var = lag_var,
    d_stat = d_stat,
    rho_hat = rho_hat,
    h_stat = h_stat,
    p_value = p_val,
    n = n,
    var_alpha = var_alpha,
    is_valid = is_valid,
    digits = digits
  )
  class(res) <- "econ_durbin_h"
  return(res)
}

#' @export
print.econ_durbin_h <- function(x, ...) {
  d <- x$digits
  print_line("=")
  cat("             DURBIN'S h TEST FOR AUTOREGRESSIVE MODELS                  \n")
  cat(paste0("  Formula: ", deparse(x$formula), "\n"))
  cat(sprintf("  Lagged Dependent Variable: %s | Sample Size (n): %d\n", x$lag_var, x$n))
  print_line("=")
  cat(sprintf("  H0: No first-order serial correlation in errors (rho = 0)\n"))
  cat(sprintf("  H1: First-order serial correlation is present (rho != 0)\n"))
  print_line("-")
  cat(sprintf("  Residual Autocorrelation (rho): %10.*f\n", d, x$rho_hat))
  cat(sprintf("  Standard Durbin-Watson (d):     %10.*f  (biased toward 2 in dynamic models)\n", d, x$d_stat))
  
  if (x$is_valid) {
    cat(sprintf("  Durbin's h Statistic:           %10.*f  ~ N(0, 1)\n", d, x$h_stat))
    cat(sprintf("  p-value:                        %10s  %s\n", format_pval(x$p_value), get_stars(x$p_value)))
  } else {
    cat("  Note: 1 - n*Var(alpha) <= 0; standard h cannot be computed.\n")
    cat(sprintf("  Durbin Alternative LM t-stat:   %10.*f\n", d, x$h_stat))
    cat(sprintf("  p-value:                        %10s  %s\n", format_pval(x$p_value), get_stars(x$p_value)))
  }
  print_line("-")
  
  if (x$p_value < 0.05) {
    cat("  Verdict: REJECT H0 (p < 0.05) -> Significant serial correlation detected.\n")
    cat("  Warning: OLS estimators in this dynamic model are biased and inconsistent.\n")
  } else {
    cat("  Verdict: PASS H0 (p >= 0.05) -> No significant serial correlation detected.\n")
    cat("  Notice:  Dynamic estimates appear consistent and well-specified.\n")
  }
  print_line("=")
  invisible(x)
}

#' Koyck Distributed Lag Model
#'
#' Estimates geometric distributed lag models using the Koyck transformation
#' (Koyck, 1954), computing the Short-Run Multiplier, Long-Run 
#' Multiplier, Rate of Decay, Median Lag, and Mean Lag.
#'
#' @param formula Model formula with one explanatory variable (e.g., \code{y ~ x}).
#' @param data A data.frame containing the time series.
#' @param method Estimation method: \code{"ols"} (default) or \code{"iv"} (Liviatan IV using lagged X as instrument).
#' @param digits Number of decimal places to format (default is 4).
#'
#' @return An object of class \code{econ_koyck}.
#' @export
#'
#' @examples
#' data(macro_data, package = "phoring")
#' koyck(inflation ~ unemployment, data = macro_data, method = "ols")
#' koyck(inflation ~ unemployment, data = macro_data, method = "iv")
koyck <- function(formula, data, method = c("ols", "iv"), digits = 4) {
  method <- match.arg(method)
  if (missing(data) || is.null(data)) stop("`data` must be provided.", call. = FALSE)
  
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  x_name <- names(mf)[2]
  x <- mf[[x_name]]
  y_name <- names(mf)[1]
  n <- length(y)
  
  if (n < 15) stop("Time series must have at least 15 observations.", call. = FALSE)
  
  # Construct lagged dependent variable: Y_{t-1} and lagged predictor: X_{t-1}
  y_lag <- c(NA, y[-n])
  x_lag <- c(NA, x[-n])
  
  df_koyck <- data.frame(
    Y = y[-1],
    X = x[-1],
    Y_lag = y_lag[-1],
    X_lag = x_lag[-1]
  )
  
  if (method == "ols") {
    # Fit Koyck transformed regression: Y_t = alpha*(1-lambda) + beta0 * X_t + lambda * Y_{t-1} + v_t
    mod <- stats::lm(Y ~ X + Y_lag, data = df_koyck)
    b <- stats::coef(mod)
    V <- stats::vcov(mod)
    resids <- stats::residuals(mod)
    fits <- stats::fitted(mod)
    r2_val <- summary(mod)$r.squared
  } else {
    # Liviatan (1963) IV estimator: instrument Y_lag with X_lag to eliminate MA(1) correlation
    X_mat <- cbind(Intercept = 1, X = df_koyck$X, Y_lag = df_koyck$Y_lag)
    Z_mat <- cbind(Intercept = 1, X = df_koyck$X, X_lag = df_koyck$X_lag)
    qr_Z <- qr(Z_mat)
    X_hat <- qr.fitted(qr_Z, X_mat)
    bread <- safe_solve(crossprod(X_hat))
    b <- as.vector(bread %*% crossprod(X_hat, df_koyck$Y))
    names(b) <- c("(Intercept)", "X", "Y_lag")
    resids <- as.vector(df_koyck$Y - X_mat %*% b)
    s2 <- sum(resids^2) / (nrow(df_koyck) - 3)
    V <- s2 * bread
    dimnames(V) <- list(names(b), names(b))
    fits <- as.vector(X_mat %*% b)
    mod <- stats::lm(Y ~ X + Y_lag, data = df_koyck) # for aux diagnostics
    r2_val <- max(0, 1 - sum(resids^2) / sum((df_koyck$Y - mean(df_koyck$Y))^2))
  }
  
  se <- sqrt(pmax(0, diag(V)))
  
  beta0 <- b["X"]
  lambda <- b["Y_lag"]
  alpha_star <- b["(Intercept)"]
  
  # Rate of decay lambda must be strictly between 0 and 1 for stable lag distribution
  stable <- (!is.na(lambda) && lambda > 0 && lambda < 1)
  
  # Long-run multiplier: beta_LR = beta0 / (1 - lambda)
  beta_lr <- if (!is.na(lambda) && lambda < 1) beta0 / (1 - lambda) else NA
  
  # Delta method standard error for beta_LR
  if (stable) {
    d_b0 <- 1 / (1 - lambda)
    d_lam <- beta0 / ((1 - lambda)^2)
    grad <- c(d_b0, d_lam)
    V_sub <- V[c("X", "Y_lag"), c("X", "Y_lag")]
    se_lr <- sqrt(max(0, as.numeric(t(grad) %*% V_sub %*% grad)))
    t_lr <- beta_lr / se_lr
    p_lr <- 2 * stats::pt(-abs(t_lr), df = nrow(df_koyck) - 3)
    
    median_lag <- -log(2) / log(lambda)
    mean_lag <- lambda / (1 - lambda)
  } else {
    se_lr <- NA; t_lr <- NA; p_lr <- NA
    median_lag <- NA; mean_lag <- NA
  }
  
  # Extract standard coefficient table
  t_b <- b / se
  p_b <- 2 * stats::pt(-abs(t_b), df = nrow(df_koyck) - 3)
  coef_table <- data.frame(
    Estimate  = b,
    Std_Error = se,
    t_stat    = t_b,
    p_value   = p_b,
    Signif    = sapply(p_b, get_stars),
    stringsAsFactors = FALSE
  )
  rownames(coef_table) <- c("(Intercept)", x_name, paste0(y_name, "_{t-1}"))
  
  # Durbin's h diagnostic
  dh <- tryCatch(durbin_h(mod, lag_var = "Y_lag"), error = function(e) NULL)
  
  res <- list(
    formula = formula,
    outcome = y_name,
    predictor = x_name,
    method = if (method == "ols") "OLS" else "Liviatan IV",
    coefficients = b,
    coef_table = coef_table,
    vcov = V,
    residuals = resids,
    fitted.values = fits,
    short_run_multiplier = unname(beta0),
    long_run_multiplier = unname(beta_lr),
    se_lr = se_lr,
    t_lr = t_lr,
    p_lr = p_lr,
    lambda = unname(lambda),
    median_lag = median_lag,
    mean_lag = mean_lag,
    stable = stable,
    n = nrow(df_koyck),
    r_squared = r2_val,
    durbin_h = dh,
    digits = digits
  )
  class(res) <- "econ_koyck"
  return(res)
}

#' @export
print.econ_koyck <- function(x, ...) {
  d <- x$digits
  print_line("=")
  cat("                    KOYCK DISTRIBUTED LAG MODEL                         \n")
  cat(paste0("  Outcome (Y): ", x$outcome, " | Regressor (X): ", x$predictor, "\n"))
  cat(sprintf("  Sample Size (T): %-12d R-squared: %.4f\n", x$n, x$r_squared))
  print_line("=")
  
  cat("  1. ESTIMATED KOYCK REGRESSION EQUATION:\n")
  print_line("-")
  v_w <- max(16, max(nchar(rownames(x$coef_table))))
  line_w <- max(74, v_w + 58)
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
  cat("  Signif. codes:  *** p < 0.001,  ** p < 0.01,  * p < 0.05,  . p < 0.1\n\n")
  
  cat("  2. MULTIPLIERS & LAG STRUCTURE (KOYCK GEOMETRIC DECAY):\n")
  print_line("-")
  cat(sprintf("  Short-Run Multiplier (beta0):    %10.*f  (Immediate 1-period impact)\n", d, x$short_run_multiplier))
  if (x$stable) {
    cat(sprintf("  Rate of Decay (lambda):          %10.*f  (Speed at which effect fades)\n", d, x$lambda))
    cat(sprintf("  Long-Run Multiplier (beta_LR):   %10.*f  (Total cumulative impact; p = %s)\n", 
                d, x$long_run_multiplier, format_pval(x$p_lr)))
    cat(sprintf("  Median Lag (-ln2 / ln lambda):   %10.2f periods (50%% of total impact realized)\n", x$median_lag))
    cat(sprintf("  Mean Lag (lambda / [1-lambda]):  %10.2f periods (Weighted average delay)\n", x$mean_lag))
  } else {
    cat(sprintf("  Rate of Decay (lambda):          %10.*f  (WARNING: lambda >= 1; non-convergent lag)\n", d, x$lambda))
  }
  print_line("-")
  
  # Durbin's h test if available
  if (!is.null(x$durbin_h)) {
    cat(sprintf("  Durbin's h Autocorrelation Test: h = %8.3f (p = %s)\n", 
                x$durbin_h$h_stat, format_pval(x$durbin_h$p_value)))
    if (x$durbin_h$p_value < 0.05) {
      cat("  --> WARNING: Significant autocorrelation in errors; estimates may be inconsistent.\n")
    } else {
      cat("  --> PASS: No significant autocorrelation detected in Koyck residuals.\n")
    }
  }
  print_line("=")
  invisible(x)
}

#' Prais-Winsten & Cochrane-Orcutt FGLS Autocorrelation Correction
#'
#' Fits Feasible Generalized Least Squares (FGLS) to correct for first-order AR(1) 
#' serial correlation (Prais and Winsten, 1954; Cochrane and Orcutt, 1949), using Prais-Winsten (retaining 
#' the first observation) or Cochrane-Orcutt transformation until convergence.
#'
#' @param formula Model formula (e.g., \code{y ~ x1 + x2}).
#' @param data A data.frame containing time series variables.
#' @param method Either \code{"prais"} (default, Prais-Winsten) or \code{"cochrane"} (Cochrane-Orcutt).
#' @param max_iter Maximum iterations for rho convergence (default is 50).
#' @param tol Convergence tolerance (default is 1e-4).
#' @param digits Number of decimal places to format (default is 4).
#'
#' @return An object of class \code{econ_prais}.
#' @export
#'
#' @examples
#' data(macro_data, package = "phoring")
#' prais_winsten(inflation ~ fed_funds + unemployment, data = macro_data)
prais_winsten <- function(formula, data, method = c("prais", "cochrane"), max_iter = 50, tol = 1e-4, digits = 4) {
  method <- match.arg(method)
  if (missing(data) || is.null(data)) stop("`data` must be provided.", call. = FALSE)
  
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(formula, data = mf)
  n <- length(y)
  k <- ncol(X)
  var_names <- colnames(X)
  
  # Step 1: Initial OLS
  init_lm <- stats::lm.fit(X, y)
  b <- stats::coef(init_lm)
  e <- y - as.vector(X %*% b)
  
  # Initial DW
  dw_init <- sum((e[-1] - e[-n])^2) / sum(e^2)
  
  rho <- 0
  iter <- 0
  converged <- FALSE
  
  while (iter < max_iter) {
    iter <- iter + 1
    # Estimate rho
    rho_new <- sum(e[-1] * e[-n]) / sum(e[-n]^2)
    rho_new <- max(-0.999, min(0.999, rho_new)) # bound
    
    if (abs(rho_new - rho) < tol) {
      converged <- TRUE
      rho <- rho_new
      break
    }
    rho <- rho_new
    
    # Transform data
    if (method == "prais") {
      # Retain t = 1 weighted by sqrt(1 - rho^2)
      y_trans <- c(sqrt(1 - rho^2) * y[1], y[-1] - rho * y[-n])
      X_trans <- rbind(sqrt(1 - rho^2) * X[1, ], X[-1, ] - rho * X[-n, ])
    } else {
      # Drop t = 1
      y_trans <- y[-1] - rho * y[-n]
      X_trans <- X[-1, ] - rho * X[-n, ]
    }
    
    # Re-estimate beta
    fit <- stats::lm.fit(X_trans, y_trans)
    b <- stats::coef(fit)
    e <- y - as.vector(X %*% b)
  }
  
  # Final standard errors from transformed system
  bread <- safe_solve(crossprod(X_trans))
  s2 <- sum(fit$residuals^2) / (length(y_trans) - k)
  vcov_mat <- s2 * bread
  dimnames(vcov_mat) <- list(var_names, var_names)
  
  se_b <- sqrt(pmax(0, diag(vcov_mat)))
  t_b <- b / se_b
  df_resid <- length(y_trans) - k
  p_b <- 2 * stats::pt(-abs(t_b), df = df_resid)
  
  # Transformed DW
  e_trans <- fit$residuals
  dw_final <- sum((e_trans[-1] - e_trans[-length(e_trans)])^2) / sum(e_trans^2)
  
  coef_table <- data.frame(
    Estimate  = b,
    Std_Error = se_b,
    t_stat    = t_b,
    p_value   = p_b,
    Signif    = sapply(p_b, get_stars),
    stringsAsFactors = FALSE
  )
  
  res <- list(
    formula = formula,
    method = if (method == "prais") "Prais-Winsten FGLS" else "Cochrane-Orcutt FGLS",
    coefficients = b,
    coef_table = coef_table,
    vcov = vcov_mat,
    residuals = e,
    residuals_transformed = e_trans,
    fitted.values = as.vector(X %*% b),
    rho = rho,
    iterations = iter,
    converged = converged,
    dw_initial = dw_init,
    dw_transformed = dw_final,
    n = n,
    k = k,
    digits = digits
  )
  class(res) <- "econ_prais"
  return(res)
}

#' @export
print.econ_prais <- function(x, ...) {
  d <- x$digits
  print_line("=")
  cat(sprintf("        %s AR(1) AUTOCORRELATION CORRECTION\n", toupper(x$method)))
  cat(paste0("  Formula: ", deparse(x$formula), "\n"))
  cat(sprintf("  Sample Size: %-12d Iterations to Convergence: %d (rho = %.4f)\n", 
              x$n, x$iterations, x$rho))
  cat(sprintf("  Initial DW:  %-12.4f Transformed DW (Target ~ 2.0): %.4f\n", 
              x$dw_initial, x$dw_transformed))
  print_line("=")
  
  v_w <- max(16, max(nchar(rownames(x$coef_table))))
  line_w <- max(74, v_w + 58)
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
  cat("  Econometric Note:\n")
  cat("  * By estimating rho and quasi-differencing the variables, serial correlation has\n")
  cat("    been purged from the disturbances, restoring efficiency and valid standard errors.\n")
  print_line("=")
  invisible(x)
}

#' Error Correction Model (ECM)
#'
#' Estimates an Error Correction Model based on the Granger Representation Theorem 
#' (Engle and Granger, 1987), capturing both short-run dynamic responses and the 
#' speed of adjustment towards long-run cointegrating equilibrium.
#'
#' @param formula Model formula (e.g., \code{y ~ x}).
#' @param data A data.frame containing the cointegrated time series.
#' @param lags Number of short-run difference lags to include (default is 1).
#' @param digits Number of decimal places to format (default is 4).
#'
#' @return An object of class \code{econ_ecm}.
#' @export
#'
#' @examples
#' data(macro_data, package = "phoring")
#' ecm(inflation ~ unemployment + fed_funds, data = macro_data, lags = 1)
ecm <- function(formula, data, lags = 1, digits = 4) {
  if (missing(data) || is.null(data)) stop("`data` must be provided.", call. = FALSE)
  
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  all_vars <- names(mf)
  dep_name <- all_vars[1]
  indep_names <- all_vars[-1]
  n <- length(y)
  
  # Step 1: Long-run cointegrating regression: Y_t = beta0 + beta1 X_t + ... + e_t
  mod_longrun <- stats::lm(formula, data = data)
  e_longrun <- stats::residuals(mod_longrun)
  b_longrun <- stats::coef(mod_longrun)
  
  # Step 2: Form differences and lagged cointegrating residual
  dy <- diff(y)
  dx_mat <- apply(as.matrix(data[, indep_names, drop = FALSE]), 2, diff)
  
  # Lagged error correction term: e_{t-1}
  ect_lag <- e_longrun[-n]
  
  # Subsample for regression with lags
  start_idx <- lags + 1
  n_eff <- length(dy) - lags
  
  y_diff <- dy[start_idx:length(dy)]
  ect <- ect_lag[start_idx:length(dy)]
  
  # Construct difference matrix
  X_diff <- cbind(Intercept = 1, ECT_Lag = ect)
  
  # Contemporaneous differences dX
  for (cn in indep_names) {
    col_v <- dx_mat[start_idx:length(dy), cn]
    X_diff <- cbind(X_diff, col_v)
    colnames(X_diff)[ncol(X_diff)] <- paste0("d_", cn)
  }
  
  # Add lagged differences if lags > 0
  if (lags > 0) {
    for (l in 1:lags) {
      idx_l <- (start_idx - l):(length(dy) - l)
      col_dy_l <- dy[idx_l]
      X_diff <- cbind(X_diff, col_dy_l)
      colnames(X_diff)[ncol(X_diff)] <- paste0("d_", dep_name, "_lag", l)
    }
  }
  
  # Fit ECM regression
  mod_ecm <- stats::lm(y_diff ~ X_diff[, -1, drop = FALSE])
  b_ecm <- stats::coef(mod_ecm)
  names(b_ecm) <- colnames(X_diff)
  V_ecm <- stats::vcov(mod_ecm)
  se_ecm <- sqrt(pmax(0, diag(V_ecm)))
  names(se_ecm) <- colnames(X_diff)
  
  t_ecm <- b_ecm / se_ecm
  p_ecm <- 2 * stats::pt(-abs(t_ecm), df = length(y_diff) - ncol(X_diff))
  
  coef_table <- data.frame(
    Estimate  = b_ecm,
    Std_Error = se_ecm,
    t_stat    = t_ecm,
    p_value   = p_ecm,
    Signif    = sapply(p_ecm, get_stars),
    stringsAsFactors = FALSE
  )
  
  gamma <- unname(b_ecm["ECT_Lag"])
  gamma_se <- unname(se_ecm["ECT_Lag"])
  gamma_p <- unname(p_ecm["ECT_Lag"])
  
  # Half-life of deviation shock: ln(0.5) / ln(1 + gamma)
  half_life <- if (!is.na(gamma) && gamma < 0 && gamma > -1) log(0.5) / log(1 + gamma) else NA
  
  res <- list(
    formula = formula,
    outcome = dep_name,
    coefficients = b_ecm,
    coef_table = coef_table,
    vcov = V_ecm,
    residuals = stats::residuals(mod_ecm),
    fitted.values = stats::fitted(mod_ecm),
    longrun_coefficients = b_longrun,
    ect_lag = unname(X_diff[, "ECT_Lag"]),
    y_diff = unname(y_diff),
    speed_of_adjustment = gamma,
    gamma_se = gamma_se,
    gamma_p = gamma_p,
    half_life = half_life,
    r_squared = summary(mod_ecm)$r.squared,
    n = length(y_diff),
    digits = digits
  )
  class(res) <- "econ_ecm"
  return(res)
}

#' @export
print.econ_ecm <- function(x, ...) {
  d <- x$digits
  print_line("=")
  cat("                 ERROR CORRECTION MODEL (ECM)                           \n")
  cat(paste0("  Outcome: ", x$outcome, " (in differences: d_", x$outcome, ")\n"))
  cat(sprintf("  Effective Obs: %-12d R-squared: %.4f\n", x$n, x$r_squared))
  print_line("=")
  
  cat("  1. SHORT-RUN DYNAMICS & ERROR CORRECTION ESTIMATES:\n")
  print_line("-")
  v_w <- max(16, max(nchar(rownames(x$coef_table))))
  line_w <- max(74, v_w + 58)
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
  cat("  Signif. codes:  *** p < 0.001,  ** p < 0.01,  * p < 0.05,  . p < 0.1\n\n")
  
  cat("  2. COINTEGRATING EQUILIBRIUM & SPEED OF ADJUSTMENT:\n")
  print_line("-")
  cat(sprintf("  Speed of Adjustment (gamma):     %10.*f  (p = %s)\n", 
              d, x$speed_of_adjustment, format_pval(x$gamma_p)))
  if (!is.na(x$speed_of_adjustment) && x$speed_of_adjustment < 0) {
    pct_corrected <- abs(x$speed_of_adjustment) * 100
    cat(sprintf("  --> Economic Interpretation: %.1f%% of any disequilibrium shock is corrected\n", pct_corrected))
    cat("      per time period, pulling the series back towards long-run cointegrating parity.\n")
    if (!is.na(x$half_life)) {
      cat(sprintf("  --> Half-life of deviation: %.1f periods until 50%% of equilibrium is restored.\n", x$half_life))
    }
  } else {
    cat("  --> WARNING: Speed of adjustment gamma >= 0; system does not converge to equilibrium.\n")
  }
  print_line("=")
  invisible(x)
}

#' Vector Autoregression (VAR) Model
#'
#' Fits a Vector Autoregression system (Sims, 1980) for multiple 
#' time series, computing system equation estimates, information criteria (AIC, BIC),
#' and pairwise Granger causality summary matrices.
#'
#' @param data A data.frame containing numeric time series variables.
#' @param vars Character vector of column names to include in the VAR system.
#' @param p Lag order (default is 1).
#' @param digits Number of decimal places to format (default is 4).
#'
#' @return An object of class \code{econ_var}.
#' @export
#'
#' @examples
#' data(macro_data, package = "phoring")
#' var_model(macro_data, vars = c("inflation", "unemployment", "fed_funds"), p = 1)
var_model <- function(data, vars, p = 1, digits = 4) {
  if (missing(data) || is.null(data)) stop("`data` must be provided.", call. = FALSE)
  if (missing(vars) || length(vars) < 2) {
    stop("`vars` must specify at least two time series column names.", call. = FALSE)
  }
  
  sub_data <- as.matrix(data[, vars, drop = FALSE])
  n <- nrow(sub_data)
  m <- ncol(sub_data)
  
  if (n <= p * m + 5) {
    stop("Sample size too small for specified VAR lag order.", call. = FALSE)
  }
  
  # Construct lagged regressors matrix
  # Rows: t = (p+1) to n
  y_mat <- sub_data[(p + 1):n, , drop = FALSE]
  T_eff <- nrow(y_mat)
  
  X_lags <- matrix(1, nrow = T_eff, ncol = 1)
  colnames(X_lags) <- "(Intercept)"
  
  for (l in 1:p) {
    lag_m <- sub_data[(p + 1 - l):(n - l), , drop = FALSE]
    colnames(lag_m) <- paste0(vars, "_lag", l)
    X_lags <- cbind(X_lags, lag_m)
  }
  
  # Estimate each equation via OLS
  eq_results <- list()
  resid_matrix <- matrix(NA, nrow = T_eff, ncol = m)
  colnames(resid_matrix) <- vars
  
  k_per_eq <- ncol(X_lags)
  bread <- safe_solve(crossprod(X_lags))
  
  for (j in 1:m) {
    vname <- vars[j]
    y_j <- y_mat[, j]
    b_j <- bread %*% crossprod(X_lags, y_j)
    e_j <- y_j - as.vector(X_lags %*% b_j)
    resid_matrix[, j] <- e_j
    
    s2_j <- sum(e_j^2) / (T_eff - k_per_eq)
    se_j <- sqrt(pmax(0, diag(s2_j * bread)))
    t_j <- b_j / se_j
    p_j <- 2 * stats::pt(-abs(t_j), df = T_eff - k_per_eq)
    
    coef_tab <- data.frame(
      Estimate  = as.vector(b_j),
      Std_Error = se_j,
      t_stat    = as.vector(t_j),
      p_value   = as.vector(p_j),
      Signif    = sapply(p_j, get_stars),
      stringsAsFactors = FALSE
    )
    rownames(coef_tab) <- colnames(X_lags)
    
    r2_j <- 1 - (sum(e_j^2) / sum((y_j - mean(y_j))^2))
    
    eq_results[[vname]] <- list(
      coefficients = b_j,
      coef_table = coef_tab,
      residuals = e_j,
      fitted.values = as.vector(y_j - e_j),
      y = as.vector(y_j),
      r_squared = r2_j,
      sigma = sqrt(s2_j)
    )
  }
  
  # Residual covariance matrix: Sigma = (1 / T) * E'E
  sigma_u <- crossprod(resid_matrix) / T_eff
  log_det_sigma <- log(det(sigma_u))
  
  # Information criteria:
  # AIC = ln |Sigma| + 2 * m^2 * p / T
  # BIC = ln |Sigma| + m^2 * p * ln(T) / T
  n_params <- m * k_per_eq
  aic <- log_det_sigma + (2 * n_params) / T_eff
  bic <- log_det_sigma + (n_params * log(T_eff)) / T_eff
  
  # Pairwise Granger Causality Matrix
  # Element (j, i): does var i Granger-cause var j?
  granger_p_matrix <- matrix(NA, nrow = m, ncol = m)
  rownames(granger_p_matrix) <- vars
  colnames(granger_p_matrix) <- vars
  
  for (j in 1:m) {
    for (i in 1:m) {
      if (i != j) {
        # Test joint significance of lag terms of var i in equation j
        lag_cols_i <- paste0(vars[i], "_lag", 1:p)
        # Wald F-test
        b_sub <- eq_results[[vars[j]]]$coefficients[lag_cols_i, 1]
        V_sub <- (eq_results[[vars[j]]]$sigma^2) * bread[lag_cols_i, lag_cols_i, drop = FALSE]
        f_val <- as.numeric(t(b_sub) %*% safe_solve(V_sub) %*% b_sub) / length(lag_cols_i)
        p_val <- stats::pf(f_val, df1 = length(lag_cols_i), df2 = T_eff - k_per_eq, lower.tail = FALSE)
        granger_p_matrix[j, i] <- p_val
      }
    }
  }
  
  res <- list(
    vars = vars,
    p = p,
    equations = eq_results,
    sigma_u = sigma_u,
    aic = aic,
    bic = bic,
    granger_matrix = granger_p_matrix,
    n_effective = T_eff,
    digits = digits
  )
  class(res) <- "econ_var"
  return(res)
}

#' @export
print.econ_var <- function(x, ...) {
  d <- x$digits
  print_line("=")
  cat(sprintf("               VECTOR AUTOREGRESSION: VAR(%d) MODEL\n", x$p))
  cat(sprintf("  Endogenous Variables (%d): %s\n", length(x$vars), paste(x$vars, collapse = ", ")))
  cat(sprintf("  Effective Sample Size (T): %-12d AIC: %-10.4f BIC: %.4f\n", 
              x$n_effective, x$aic, x$bic))
  print_line("=")
  
  # Equation by equation summary
  for (vname in x$vars) {
    eq <- x$equations[[vname]]
    cat(sprintf("  EQUATION: %s (R-squared: %.4f, Residual SE: %.*f)\n", 
                vname, eq$r_squared, d, eq$sigma))
    print_line("-")
    v_w <- max(16, max(nchar(rownames(eq$coef_table))))
    line_w <- max(74, v_w + 58)
    cat(sprintf(paste0("  %-", v_w, "s %10s %12s %10s %10s   %s\n"), 
                "Variable", "Estimate", "Std. Error", "t-stat", "p-value", "Signif"))
    print_line("-", width = line_w)
    for (i in 1:nrow(eq$coef_table)) {
      rname <- rownames(eq$coef_table)[i]
      row <- eq$coef_table[i, ]
      cat(sprintf(paste0("  %-", v_w, "s %10.*f %12.*f %10.3f %10s   %s\n"),
                  rname,
                  d, row$Estimate,
                  d, row$Std_Error,
                  row$t_stat,
                  format_pval(row$p_value),
                  row$Signif))
    }
    print_line("-", width = line_w)
    cat("\n")
  }
  
  # Granger Causality Matrix
  cat("  PAIRWISE GRANGER CAUSALITY MATRIX:\n")
  cat("  (Rows = Dependent Variable | Columns = Explanatory Cause; p-values shown)\n")
  print_line("-")
  m <- length(x$vars)
  v_w <- max(14, max(nchar(x$vars)))
  col_w <- 14
  header <- sprintf(paste0("  %-", v_w, "s"), "Dep \\ Cause")
  for (v in x$vars) {
    header <- paste0(header, sprintf(paste0("%", col_w, "s"), v))
  }
  cat(paste0(header, "\n"))
  print_line("-", width = max(74, v_w + m * col_w + 4))
  
  for (i in 1:m) {
    r_str <- sprintf(paste0("  %-", v_w, "s"), x$vars[i])
    for (j in 1:m) {
      if (i == j) {
        r_str <- paste0(r_str, sprintf(paste0("%", col_w, "s"), "-"))
      } else {
        p_val <- x$granger_matrix[i, j]
        stars <- get_stars(p_val)
        val_str <- sprintf("%.*f%s", d, p_val, stars)
        r_str <- paste0(r_str, sprintf(paste0("%", col_w, "s"), val_str))
      }
    }
    cat(paste0(r_str, "\n"))
  }
  print_line("-", width = max(74, v_w + m * col_w + 4))
  cat("  Signif. codes:  *** p < 0.001,  ** p < 0.01,  * p < 0.05,  . p < 0.1\n")
  print_line("=")
  invisible(x)
}

# -------------------------------------------------------------
# S3 Methods for Dynamic Time Series Models
# -------------------------------------------------------------

#' @export
coef.econ_koyck <- function(object, ...) object$coefficients

#' @export
vcov.econ_koyck <- function(object, ...) object$vcov

#' @export
residuals.econ_koyck <- function(object, ...) object$residuals

#' @export
fitted.econ_koyck <- function(object, ...) object$fitted.values

#' @export
nobs.econ_koyck <- function(object, ...) object$n

#' @export
confint.econ_koyck <- function(object, parm, level = 0.95, ...) {
  b <- stats::coef(object)
  se <- sqrt(pmax(0, diag(stats::vcov(object))))
  a <- (1 - level) / 2
  crit <- stats::qt(1 - a, df = object$n - length(b))
  ci <- cbind(b - crit * se, b + crit * se)
  colnames(ci) <- c(sprintf("%.1f %%", 100 * a), sprintf("%.1f %%", 100 * (1 - a)))
  rownames(ci) <- names(b)
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

#' @export
predict.econ_koyck <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$fitted.values)
  stop("Out-of-sample prediction for Koyck model requires lag initialization.", call. = FALSE)
}

#' @export
coef.econ_prais <- function(object, ...) object$coefficients

#' @export
vcov.econ_prais <- function(object, ...) object$vcov

#' @export
residuals.econ_prais <- function(object, ...) object$residuals

#' @export
fitted.econ_prais <- function(object, ...) object$fitted.values

#' @export
nobs.econ_prais <- function(object, ...) object$n

#' @export
confint.econ_prais <- function(object, parm, level = 0.95, ...) {
  b <- stats::coef(object)
  se <- sqrt(pmax(0, diag(stats::vcov(object))))
  a <- (1 - level) / 2
  crit <- stats::qt(1 - a, df = object$n - length(b))
  ci <- cbind(b - crit * se, b + crit * se)
  colnames(ci) <- c(sprintf("%.1f %%", 100 * a), sprintf("%.1f %%", 100 * (1 - a)))
  rownames(ci) <- names(b)
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

#' @export
coef.econ_ecm <- function(object, ...) object$coefficients

#' @export
vcov.econ_ecm <- function(object, ...) object$vcov

#' @export
residuals.econ_ecm <- function(object, ...) object$residuals

#' @export
fitted.econ_ecm <- function(object, ...) object$fitted.values

#' @export
nobs.econ_ecm <- function(object, ...) object$n

#' @export
confint.econ_ecm <- function(object, parm, level = 0.95, ...) {
  b <- stats::coef(object)
  se <- sqrt(pmax(0, diag(stats::vcov(object))))
  a <- (1 - level) / 2
  crit <- stats::qt(1 - a, df = object$n - length(b))
  ci <- cbind(b - crit * se, b + crit * se)
  colnames(ci) <- c(sprintf("%.1f %%", 100 * a), sprintf("%.1f %%", 100 * (1 - a)))
  rownames(ci) <- names(b)
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

#' @export
coef.econ_var <- function(object, ...) {
  lapply(object$equations, function(eq) eq$coefficients[, 1])
}

#' @export
nobs.econ_var <- function(object, ...) object$n_effective

