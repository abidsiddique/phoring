#' Quantile Regression Across the Distribution
#'
#' Estimates quantile regressions for single or multiple conditional quantiles 
#' (e.g., 10th, 25th, median 50th, 75th, 90th percentiles), allowing researchers 
#' to analyze how economic covariates impact different parts of the outcome distribution 
#' beyond the conditional mean.
#'
#' @param formula Model formula (e.g., \code{wage ~ educ + exper}).
#' @param data A data.frame containing the variables.
#' @param tau Quantile(s) to estimate, between 0 and 1 (default is \code{c(0.1, 0.25, 0.5, 0.75, 0.9)}).
#' @param se Standard error type: \code{"robust"} (Huber sandwich, default) or \code{"classical"}.
#' @param digits Number of decimal places to format (default is 4).
#'
#' @return An object of class \code{econ_qreg}.
#' @export
#'
#' @examples
#' data(wage_data, package = "phoring")
#' # Estimate median regression (tau = 0.5)
#' qreg(wage ~ educ + exper, data = wage_data, tau = 0.5)
#' 
#' # Estimate full distribution quantiles (10%, 25%, 50%, 75%, 90%)
#' q_all <- qreg(wage ~ educ + exper, data = wage_data, tau = c(0.25, 0.5, 0.75))
qreg <- function(formula, data, tau = c(0.1, 0.25, 0.5, 0.75, 0.9), se = c("robust", "classical"), digits = 4) {
  se <- match.arg(se)
  if (missing(data) || is.null(data)) stop("`data` must be provided.", call. = FALSE)
  if (any(tau <= 0 | tau >= 1)) stop("All values of `tau` must be strictly between 0 and 1.", call. = FALSE)
  
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(formula, data = mf)
  
  n <- nrow(X)
  k <- ncol(X)
  var_names <- colnames(X)
  
  # Starting values from OLS
  ols_fit <- stats::lm.fit(X, y)
  b_init <- stats::coef(ols_fit)
  b_init[is.na(b_init)] <- 0
  
  results_by_tau <- list()
  coef_matrix <- matrix(NA, nrow = k, ncol = length(tau))
  rownames(coef_matrix) <- var_names
  colnames(coef_matrix) <- paste0("tau_", tau)
  
  se_matrix <- matrix(NA, nrow = k, ncol = length(tau))
  rownames(se_matrix) <- var_names
  colnames(se_matrix) <- paste0("tau_", tau)
  
  # Smooth quantile loss function and gradient:
  # rho_tau^eps(u) = u * (tau - 0.5) + 0.5 * sqrt(u^2 + eps^2)
  eps <- 1e-4
  
  for (t_idx in seq_along(tau)) {
    curr_tau <- tau[t_idx]
    
    loss_fn <- function(b) {
      u <- y - as.vector(X %*% b)
      sum(u * (curr_tau - 0.5) + 0.5 * sqrt(u^2 + eps^2))
    }
    
    grad_fn <- function(b) {
      u <- y - as.vector(X %*% b)
      w <- (curr_tau - 0.5) + 0.5 * (u / sqrt(u^2 + eps^2))
      -as.vector(crossprod(X, w))
    }
    
    # BFGS optimization on smooth loss
    opt <- stats::optim(par = b_init, fn = loss_fn, gr = grad_fn, method = "BFGS")
    b_est <- opt$par
    
    # Exact check function for polishing: rho(u) = u * (tau - I(u < 0))
    exact_loss <- function(b) {
      u <- y - as.vector(X %*% b)
      sum(u * (curr_tau - (u < 0)))
    }
    
    opt_exact <- stats::optim(par = b_est, fn = exact_loss, method = "Nelder-Mead", 
                              control = list(maxit = 150))
    b_est <- opt_exact$par
    names(b_est) <- var_names
    
    # Residuals
    e <- y - as.vector(X %*% b_est)
    
    # Standard Errors: Koenker-Bassett / Powell kernel estimator
    # Bandwidth h for error density f_e(0): Bofinger / Silverman rule
    h_bw <- 1.06 * stats::sd(e) * (n^(-1/5))
    if (h_bw <= 0) h_bw <- 0.1
    # Epanechnikov kernel weights: K(u) = 0.75 * (1 - u^2) for |u| <= 1
    u_scaled <- e / h_bw
    k_weights <- ifelse(abs(u_scaled) <= 1, (0.75 / h_bw) * (1 - u_scaled^2), 0)
    
    bread <- safe_solve(crossprod(X))
    if (se == "robust") {
      # Huber sandwich: V = D^(-1) * Omega * D^(-1)
      # D = X' diag(f_i) X
      D_mat <- crossprod(X * sqrt(pmax(k_weights, 1e-6)))
      D_inv <- safe_solve(D_mat + diag(1e-6, k))
      Omega <- curr_tau * (1 - curr_tau) * crossprod(X)
      vcov_mat <- D_inv %*% Omega %*% D_inv
    } else {
      # Classical i.i.d. formula: V = [tau(1-tau) / f(0)^2] * (X'X)^(-1)
      f0 <- mean(k_weights)
      f0 <- max(f0, 1e-4)
      vcov_mat <- (curr_tau * (1 - curr_tau) / (f0^2)) * bread
    }
    dimnames(vcov_mat) <- list(var_names, var_names)
    
    se_b <- sqrt(pmax(0, diag(vcov_mat)))
    t_b <- b_est / se_b
    p_b <- 2 * stats::pnorm(-abs(t_b))
    
    coef_tab <- data.frame(
      Estimate  = b_est,
      Std_Error = se_b,
      t_stat    = t_b,
      p_value   = p_b,
      Signif    = sapply(p_b, get_stars),
      stringsAsFactors = FALSE
    )
    
    results_by_tau[[as.character(curr_tau)]] <- list(
      tau = curr_tau,
      coefficients = b_est,
      se = se_b,
      vcov = vcov_mat,
      coef_table = coef_tab,
      residuals = e,
      fitted.values = as.vector(X %*% b_est),
      loss = opt_exact$value
    )
    
    coef_matrix[, t_idx] <- b_est
    se_matrix[, t_idx]   <- se_b
  }
  
  # Default residuals and fitted values (tau = 0.5 or first tau)
  def_tau_str <- if ("0.5" %in% names(results_by_tau)) "0.5" else names(results_by_tau)[1]
  
  res <- list(
    formula = formula,
    taus = tau,
    coefficients = coef_matrix,
    se = se_matrix,
    results_by_tau = results_by_tau,
    residuals = results_by_tau[[def_tau_str]]$residuals,
    fitted.values = results_by_tau[[def_tau_str]]$fitted.values,
    losses = sapply(results_by_tau, function(r) r$loss),
    n = n,
    k = k,
    se_type = if (se == "robust") "Robust (Huber Sandwich)" else "Classical (I.I.D.)",
    digits = digits
  )
  class(res) <- "econ_qreg"
  return(res)
}

#' @export
print.econ_qreg <- function(x, ...) {
  d <- x$digits
  print_line("=")
  cat("                   QUANTILE REGRESSION RESULTS                          \n")
  cat(paste0("  Formula: ", deparse(x$formula), "\n"))
  cat(sprintf("  Sample Size: %-12d Number of Predictors: %d\n", x$n, x$k))
  cat(sprintf("  Estimated Quantiles: %s\n", paste(sprintf("%.2f", x$taus), collapse = ", ")))
  cat(paste0("  Standard Errors: ", x$se_type, "\n"))
  print_line("=")
  
  if (length(x$taus) == 1) {
    # Single Quantile Printout
    q_res <- x$results_by_tau[[1]]
    v_w <- max(16, max(nchar(rownames(q_res$coef_table))))
    line_w <- max(74, v_w + 58)
    cat(sprintf("  Quantile (tau = %.2f):\n", q_res$tau))
    print_line("-", width = line_w)
    cat(sprintf(paste0("  %-", v_w, "s %10s %12s %10s %10s   %s\n"), 
                "Variable", "Estimate", "Std. Error", "t-stat", "p-value", "Signif"))
    print_line("-", width = line_w)
    
    for (i in 1:nrow(q_res$coef_table)) {
      vname <- rownames(q_res$coef_table)[i]
      row <- q_res$coef_table[i, ]
      cat(sprintf(paste0("  %-", v_w, "s %10.*f %12.*f %10.3f %10s   %s\n"),
                  vname,
                  d, row$Estimate,
                  d, row$Std_Error,
                  row$t_stat,
                  format_pval(row$p_value),
                  row$Signif))
    }
    print_line("-", width = line_w)
  } else {
    # Multi-Quantile Comparative Distribution Table
    v_w <- max(16, max(nchar(rownames(x$coefficients))))
    col_w <- 14
    total_w <- max(74, v_w + length(x$taus) * col_w + 4)
    
    cat("  DISTRIBUTIONAL PROFILE: ESTIMATES & SE ACROSS QUANTILES:\n")
    print_line("-", width = total_w)
    header <- sprintf(paste0("  %-", v_w, "s"), "Variable")
    for (tu in x$taus) {
      header <- paste0(header, sprintf(paste0("%", col_w, "s"), paste0("tau = ", sprintf("%.2f", tu))))
    }
    cat(paste0(header, "\n"))
    print_line("-", width = total_w)
    
    for (i in 1:nrow(x$coefficients)) {
      vname <- rownames(x$coefficients)[i]
      b_line <- sprintf(paste0("  %-", v_w, "s"), vname)
      se_line <- sprintf(paste0("  %-", v_w, "s"), "")
      
      for (tu in x$taus) {
        q_item <- x$results_by_tau[[as.character(tu)]]
        b_val <- q_item$coef_table[vname, "Estimate"]
        se_val <- q_item$coef_table[vname, "Std_Error"]
        p_val <- q_item$coef_table[vname, "p_value"]
        stars <- get_stars(p_val)
        
        b_entry <- sprintf("%.*f%s", d, b_val, stars)
        se_entry <- sprintf("(%.*f)", d, se_val)
        
        b_line <- paste0(b_line, sprintf(paste0("%", col_w, "s"), b_entry))
        se_line <- paste0(se_line, sprintf(paste0("%", col_w, "s"), se_entry))
      }
      cat(paste0(b_line, "\n"))
      cat(paste0(se_line, "\n"))
    }
    print_line("-", width = total_w)
  }
  
  cat("  Signif. codes:  *** p < 0.001,  ** p < 0.01,  * p < 0.05,  . p < 0.1\n")
  cat("  Econometric Note:\n")
  cat("  * Coefficients show the marginal effect of each predictor on the tau-th quantile\n")
  cat("    of the outcome. Diverging estimates across tau indicate distributional heterogeneity.\n")
  print_line("=")
  invisible(x)
}

#' @export
coef.econ_qreg <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.econ_qreg <- function(object, tau = NULL, ...) {
  if (is.null(tau)) {
    lapply(object$results_by_tau, function(item) item$vcov)
  } else {
    object$results_by_tau[[as.character(tau)]]$vcov
  }
}

#' @export
residuals.econ_qreg <- function(object, ...) object$residuals

#' @export
fitted.econ_qreg <- function(object, ...) object$fitted.values

#' @export
nobs.econ_qreg <- function(object, ...) object$n

#' @export
confint.econ_qreg <- function(object, parm, level = 0.95, tau = 0.5, ...) {
  tau_str <- as.character(tau)
  if (!tau_str %in% names(object$results_by_tau)) tau_str <- names(object$results_by_tau)[1]
  item <- object$results_by_tau[[tau_str]]
  b <- item$coefficients
  se <- item$se
  a <- (1 - level) / 2
  crit <- stats::qnorm(1 - a)
  ci <- cbind(b - crit * se, b + crit * se)
  colnames(ci) <- c(sprintf("%.1f %%", 100 * a), sprintf("%.1f %%", 100 * (1 - a)))
  rownames(ci) <- names(b)
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

#' @export
predict.econ_qreg <- function(object, newdata = NULL, tau = 0.5, ...) {
  tau_str <- as.character(tau)
  if (!tau_str %in% names(object$results_by_tau)) tau_str <- names(object$results_by_tau)[1]
  b <- object$results_by_tau[[tau_str]]$coefficients
  if (is.null(newdata)) return(object$results_by_tau[[tau_str]]$fitted.values)
  form_no_y <- stats::delete.response(stats::terms(object$formula))
  X_new <- stats::model.matrix(form_no_y, data = newdata)
  as.vector(X_new[, names(b), drop = FALSE] %*% b)
}
