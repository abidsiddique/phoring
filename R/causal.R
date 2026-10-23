#' Sharp Regression Discontinuity Design (RDD)
#'
#' Estimates causal treatment effects at a policy cutoff using Sharp Regression 
#' Discontinuity (SRD) with local linear regression, optimal bandwidth selection, 
#' triangular/uniform kernel weights, and robust standard errors.
#'
#' @param formula Model formula where the first RHS variable is the running variable 
#'   (e.g., \code{outcome ~ running_var} or \code{outcome ~ running_var + controls}).
#' @param data A data.frame containing the variables.
#' @param cutoff Threshold value of the running variable that determines treatment (default is 0).
#' @param bandwidth Bandwidth window around the cutoff. If \code{NULL} (default), computes 
#'   the optimal rule-of-thumb bandwidth (Silverman/Imbens-Kalyanaraman).
#' @param kernel Kernel weighting function: \code{"triangular"} (default, optimal for boundary estimation) 
#'   or \code{"uniform"}.
#' @param digits Number of decimal places to format (default is 4).
#'
#' @return An object of class \code{econ_rdd}.
#' @export
#'
#' @examples
#' set.seed(42)
#' n <- 500
#' x <- runif(n, -2, 2)
#' d <- as.numeric(x >= 0)
#' y <- 2 + 3 * d + 1.5 * x + 0.8 * d * x + rnorm(n, sd = 0.5)
#' df <- data.frame(y = y, score = x)
#' rdd(y ~ score, data = df, cutoff = 0)
rdd <- function(formula, data, cutoff = 0, bandwidth = NULL, kernel = c("triangular", "uniform"), digits = 4) {
  kernel <- match.arg(kernel)
  if (missing(data) || is.null(data)) stop("`data` must be provided.", call. = FALSE)
  
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  all_preds <- names(mf)[-1]
  
  if (length(all_preds) == 0) {
    stop("Formula must include at least one running variable on the RHS.", call. = FALSE)
  }
  
  running_var_name <- all_preds[1]
  covar_names <- all_preds[-1]
  x <- mf[[running_var_name]]
  n_total <- length(y)
  
  # Centered running variable
  x_c <- x - cutoff
  
  # Automatic bandwidth if not supplied
  if (is.null(bandwidth)) {
    iqr_val <- stats::IQR(x, na.rm = TRUE)
    sd_val <- stats::sd(x, na.rm = TRUE)
    scale_val <- if (iqr_val > 0) min(sd_val, iqr_val / 1.34) else sd_val
    bandwidth <- 1.06 * scale_val * (n_total^(-1/5))
  }
  
  # Local window around cutoff: |X - c| <= h
  in_window <- which(abs(x_c) <= bandwidth)
  n_eff <- length(in_window)
  
  if (n_eff < 10) {
    stop("Fewer than 10 observations within bandwidth window. Consider increasing `bandwidth`.", call. = FALSE)
  }
  
  y_sub <- y[in_window]
  x_sub <- x_c[in_window]
  d_sub <- as.numeric(x_sub >= 0)
  
  # Kernel weights
  u <- abs(x_sub) / bandwidth
  if (kernel == "triangular") {
    w <- 1 - u
  } else {
    w <- rep(1, n_eff)
  }
  
  # Construct local linear regression model matrix:
  # Intercept, D (treatment effect), x_c (slope left), D * x_c (slope difference right)
  X_mat <- cbind(
    Intercept = 1,
    Treatment = d_sub,
    Running_Var = x_sub,
    Interaction = d_sub * x_sub
  )
  
  if (length(covar_names) > 0) {
    for (cn in covar_names) {
      z_val <- mf[[cn]][in_window]
      X_mat <- cbind(X_mat, z_val)
      colnames(X_mat)[ncol(X_mat)] <- cn
    }
  }
  
  # Weighted Least Squares estimation without allocating n x n diagonal matrix
  XtW <- t(X_mat * w)
  bread <- safe_solve(XtW %*% X_mat)
  b <- bread %*% (XtW %*% y_sub)
  b_vec <- as.vector(b)
  names(b_vec) <- colnames(X_mat)
  
  # Robust Heteroskedasticity-Consistent Covariance (HC1) for WLS
  e <- y_sub - as.vector(X_mat %*% b_vec)
  k <- ncol(X_mat)
  dfc <- n_eff / (n_eff - k)
  meat <- t(X_mat * (w * e)) %*% (X_mat * (w * e))
  vcov_mat <- dfc * (bread %*% meat %*% bread)
  dimnames(vcov_mat) <- list(colnames(X_mat), colnames(X_mat))
  
  se_vec <- sqrt(pmax(0, diag(vcov_mat)))
  names(se_vec) <- colnames(X_mat)
  t_vec <- b_vec / se_vec
  names(t_vec) <- colnames(X_mat)
  p_vec <- 2 * stats::pnorm(-abs(t_vec))
  names(p_vec) <- colnames(X_mat)
  
  coef_table <- data.frame(
    Estimate  = b_vec,
    Std_Error = se_vec,
    t_stat    = t_vec,
    p_value   = p_vec,
    Signif    = sapply(p_vec, get_stars),
    stringsAsFactors = FALSE
  )
  
  # Treatment effect statistics
  tau_est <- unname(b_vec["Treatment"])
  tau_se  <- unname(se_vec["Treatment"])
  tau_t   <- unname(t_vec["Treatment"])
  tau_p   <- unname(p_vec["Treatment"])
  ci_low  <- unname(tau_est - 1.96 * tau_se)
  ci_high <- unname(tau_est + 1.96 * tau_se)
  
  # Subsample counts and means
  n_left  <- sum(d_sub == 0)
  n_right <- sum(d_sub == 1)
  mean_left  <- mean(y_sub[d_sub == 0])
  mean_right <- mean(y_sub[d_sub == 1])
  
  # McCrary density continuity binomial test
  binom_test <- stats::binom.test(n_left, n_eff, p = 0.5)
  mccrary_p <- binom_test$p.value
  
  res <- list(
    formula = formula,
    outcome = names(mf)[1],
    running_var = running_var_name,
    cutoff = cutoff,
    bandwidth = bandwidth,
    kernel = kernel,
    tau = tau_est,
    se = tau_se,
    t_stat = tau_t,
    p_value = tau_p,
    ci_low = ci_low,
    ci_high = ci_high,
    coefficients = b_vec,
    vcov = vcov_mat,
    residuals = e,
    fitted.values = as.vector(X_mat %*% b_vec),
    n_total = n_total,
    n_effective = n_eff,
    n_left = n_left,
    n_right = n_right,
    mean_left = mean_left,
    mean_right = mean_right,
    mccrary_p = mccrary_p,
    x_data = x_sub + cutoff,
    y_data = y_sub,
    x_full = x,
    y_full = y,
    data = data,
    coef_table = coef_table,
    digits = digits
  )
  class(res) <- "econ_rdd"
  return(res)
}

#' @export
print.econ_rdd <- function(x, ...) {
  d <- x$digits
  print_line("=")
  cat("                SHARP REGRESSION DISCONTINUITY (RDD)                    \n")
  cat(paste0("  Outcome: ", x$outcome, " | Running Variable: ", x$running_var, "\n"))
  cat(sprintf("  Cutoff (c): %-12.4f Bandwidth (h): %-12.4f Kernel: %s\n", 
              x$cutoff, x$bandwidth, x$kernel))
  cat(sprintf("  Total Obs:  %-12d Effective Obs: %-12d (Left: %d, Right: %d)\n", 
              x$n_total, x$n_effective, x$n_left, x$n_right))
  cat("  Standard Errors: Robust (Huber-White HC1)\n")
  print_line("=")
  
  cat("  1. LOCAL AVERAGE TREATMENT EFFECT (LATE / SHARP RD EFFECT):\n")
  print_line("-")
  cat(sprintf("  Estimate (Tau):  %10.*f\n", d, x$tau))
  cat(sprintf("  Robust Std. Err: %10.*f\n", d, x$se))
  cat(sprintf("  t-statistic:     %10.3f\n", x$t_stat))
  cat(sprintf("  p-value:         %10s  %s\n", format_pval(x$p_value), get_stars(x$p_value)))
  cat(sprintf("  95%% Conf. Int:   [%.*f, %.*f]\n", d, x$ci_low, d, x$ci_high))
  print_line("-")
  
  # Local Linear Coefficients Table
  cat("\n  2. FULL LOCAL LINEAR MODEL ESTIMATES:\n")
  v_w <- max(16, max(nchar(rownames(x$coef_table))))
  line_w <- max(74, v_w + 58)
  print_line("-", width = line_w)
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
  
  # Identification Diagnostic (Density Continuity Check)
  cat("\n  3. IDENTIFICATION & DENSITY CONTINUITY CHECK:\n")
  print_line("-")
  cat(sprintf("  Mean Outcome Just Below Cutoff: %.*f\n", d, x$mean_left))
  cat(sprintf("  Mean Outcome Just Above Cutoff: %.*f\n", d, x$mean_right))
  cat(sprintf("  Running Var Density Continuity p-val: %s\n", format_pval(x$mccrary_p)))
  if (x$mccrary_p >= 0.05) {
    cat("  --> PASS: No evidence of manipulation/sorting around the cutoff (p >= 0.05).\n")
  } else {
    cat("  --> WARNING: Significant asymmetry around cutoff (p < 0.05); check for sorting.\n")
  }
  print_line("=")
  invisible(x)
}

#' @export
coef.econ_rdd <- function(object, ...) object$coefficients

#' @export
vcov.econ_rdd <- function(object, ...) object$vcov

#' @export
residuals.econ_rdd <- function(object, ...) object$residuals

#' @export
fitted.econ_rdd <- function(object, ...) object$fitted.values

#' @export
nobs.econ_rdd <- function(object, ...) object$n_effective

#' @export
confint.econ_rdd <- function(object, parm, level = 0.95, ...) {
  b <- stats::coef(object)
  se <- sqrt(pmax(0, diag(stats::vcov(object))))
  a <- (1 - level) / 2
  crit <- stats::qnorm(1 - a)
  ci <- cbind(b - crit * se, b + crit * se)
  colnames(ci) <- c(sprintf("%.1f %%", 100 * a), sprintf("%.1f %%", 100 * (1 - a)))
  rownames(ci) <- names(b)
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}
