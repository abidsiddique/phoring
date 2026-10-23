#' Censored Regression (Tobit Model) with McDonald-Moffitt Decomposition
#'
#' Fits a Tobit regression model (left-censored at 0 or a user-specified threshold),
#' calculates robust or standard standard errors, and automatically produces the
#' McDonald-Moffitt decomposition of marginal effects.
#'
#' @param formula Model formula (e.g., \code{hours ~ educ + exper + age}).
#' @param data A data.frame containing the variables.
#' @param left Left-censoring cutoff point (default is 0).
#' @param se Type of standard errors: \code{"robust"} (sandwich, default) or \code{"standard"} (MLE Hessian).
#' @param digits Number of decimal places to display (default is 4).
#'
#' @return An object of class \code{econ_tobit}.
#' @export
#'
#' @examples
#' data(labor_supply, package = "phoring")
#' labor_supply$hours <- pmax(labor_supply$hours, 0)
#' tobit(hours ~ educ + exper + age, data = labor_supply, left = 0)
tobit <- function(formula, data, left = 0, se = c("robust", "standard"), digits = 4) {
  se <- match.arg(se)
  if (missing(data) || is.null(data)) stop("Please provide data.", call. = FALSE)
  
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(formula, data = data)
  
  n <- length(y)
  k <- ncol(X)
  c_idx <- (y <= left) # Censored indicator
  
  # Initial values from OLS
  ols_init <- stats::lm.fit(x = X, y = y)
  b_init <- stats::coef(ols_init)
  b_init[is.na(b_init)] <- 0
  sigma_init <- max(0.5, sqrt(sum(stats::residuals(ols_init)^2) / n))
  par_init <- c(b_init, log(sigma_init))
  
  # Log-likelihood function
  loglik_fn <- function(par) {
    b_v <- par[1:k]
    sig <- exp(par[k + 1])
    eta <- as.vector(X %*% b_v)
    z_c <- (left - eta[c_idx]) / sig
    z_u <- (y[!c_idx] - eta[!c_idx]) / sig
    
    ll_c <- sum(stats::pnorm(z_c, log.p = TRUE))
    ll_u <- sum(-log(sig) + stats::dnorm(z_u, log = TRUE))
    -(ll_c + ll_u)
  }
  
  opt <- stats::optim(par = par_init, fn = loglik_fn, method = "BFGS", hessian = TRUE)
  par_est <- opt$par
  b_est <- par_est[1:k]
  names(b_est) <- colnames(X)
  sigma_est <- exp(par_est[k + 1])
  
  # Covariance matrix via Hessian
  inv_hessian <- tryCatch(safe_solve(opt$hessian), error = function(e) diag(1e-4, length(par_init)))
  
  if (se == "robust") {
    scores <- matrix(0, nrow = n, ncol = length(par_est))
    eta <- as.vector(X %*% b_est)
    z_c <- (left - eta) / sigma_est
    z_u <- (y - eta) / sigma_est
    
    for (i in 1:n) {
      if (c_idx[i]) {
        mills <- stats::dnorm(z_c[i]) / pmax(1e-10, stats::pnorm(z_c[i]))
        scores[i, 1:k] <- - (1 / sigma_est) * mills * X[i, ]
        scores[i, k + 1] <- - (z_c[i]) * mills
      } else {
        scores[i, 1:k] <- (z_u[i] / sigma_est) * X[i, ]
        scores[i, k + 1] <- (z_u[i]^2 - 1)
      }
    }
    meat <- crossprod(scores)
    vcov_mat <- inv_hessian %*% meat %*% inv_hessian
  } else {
    vcov_mat <- inv_hessian
  }
  
  vcov_b <- vcov_mat[1:k, 1:k, drop = FALSE]
  dimnames(vcov_b) <- list(colnames(X), colnames(X))
  se_b <- sqrt(pmax(0, diag(vcov_b)))
  z_b <- b_est / se_b
  p_b <- 2 * stats::pnorm(-abs(z_b))
  
  coef_table <- data.frame(
    Estimate  = b_est,
    Std_Error = se_b,
    z_stat    = z_b,
    p_value   = p_b,
    Signif    = sapply(p_b, get_stars),
    stringsAsFactors = FALSE
  )
  
  # McDonald-Moffitt Decomposition
  mean_eta <- mean(as.vector(X %*% b_est))
  z_mean <- (mean_eta - left) / sigma_est
  phi_z <- stats::dnorm(z_mean)
  Phi_z <- stats::pnorm(z_mean)
  lambda_z <- phi_z / max(1e-10, Phi_z)
  
  preds <- colnames(X)[colnames(X) != "(Intercept)"]
  ame_uncensored_prob <- numeric(length(preds)) # d P(y > left) / d x
  ame_conditional     <- numeric(length(preds)) # d E(y | y > left) / d x
  ame_unconditional   <- numeric(length(preds)) # d E(y) / d x
  names(ame_uncensored_prob) <- preds
  names(ame_conditional) <- preds
  names(ame_unconditional) <- preds
  
  for (pn in preds) {
    b_val <- b_est[pn]
    ame_uncensored_prob[pn] <- (phi_z / sigma_est) * b_val
    ame_conditional[pn]     <- b_val * (1 - z_mean * lambda_z - lambda_z^2)
    ame_unconditional[pn]   <- Phi_z * b_val
  }
  
  mm_table <- data.frame(
    Variable       = preds,
    Prob_Uncens    = ame_uncensored_prob,
    Cond_Mean      = ame_conditional,
    Uncond_Mean    = ame_unconditional,
    stringsAsFactors = FALSE
  )
  
  fitted_val <- pmax(left, as.vector(X %*% b_est))
  residuals_val <- y - fitted_val
  
  res <- list(
    formula = formula,
    coefficients = b_est,
    sigma = sigma_est,
    coef_table = coef_table,
    mm_table = mm_table,
    vcov = vcov_b,
    residuals = residuals_val,
    fitted.values = fitted_val,
    n = n,
    k = k,
    n_censored = sum(c_idx),
    n_uncensored = sum(!c_idx),
    left = left,
    y = y,
    X = X,
    linear_predictor = as.vector(X %*% b_est),
    log_lik = -opt$value,
    se_type = if (se == "robust") "Robust (Sandwich)" else "Standard (MLE)",
    digits = digits
  )
  class(res) <- "econ_tobit"
  return(res)
}

#' @export
print.econ_tobit <- function(x, digits = NULL, ...) {
  d <- if (is.null(digits)) x$digits else digits
  
  print_line("=")
  cat("                TOBIT CENSORED REGRESSION (TOBIT MODEL)                 \n")
  cat(paste0("  Formula: ", deparse(x$formula), "\n"))
  cat(sprintf("  Sample Size: %-10d Uncensored: %-10d Censored (<= %g): %d\n", 
              x$n, x$n_uncensored, x$left, x$n_censored))
  cat(sprintf("  Log-Likelihood: %-10.2f Scale (Sigma): %.4f\n", x$log_lik, x$sigma))
  cat(paste0("  Standard Errors: ", x$se_type, "\n"))
  print_line("=")
  
  v_w <- max(16, max(nchar(rownames(x$coef_table))))
  line_w <- max(74, v_w + 58)
  cat("  1. LATENT VARIABLE COEFFICIENTS (INDEX MODEL):\n")
  print_line("-", width = line_w)
  cat(sprintf(paste0("  %-", v_w, "s %10s %12s %10s %10s   %s\n"), 
              "Variable", "Estimate", "Std. Error", "z-stat", "p-value", "Signif"))
  print_line("-", width = line_w)
  
  for (i in 1:nrow(x$coef_table)) {
    vname <- rownames(x$coef_table)[i]
    row <- x$coef_table[i, ]
    cat(sprintf(paste0("  %-", v_w, "s %10.*f %12.*f %10.3f %10s   %s\n"),
                vname,
                d, row$Estimate,
                d, row$Std_Error,
                row$z_stat,
                format_pval(row$p_value),
                row$Signif))
  }
  print_line("-", width = line_w)
  cat("  Signif. codes:  *** p < 0.001,  ** p < 0.01,  * p < 0.05,  . p < 0.1\n\n")
  
  # McDonald-Moffitt Decomposition Table
  v_w_mm <- max(16, max(nchar(x$mm_table$Variable)))
  line_w_mm <- max(74, v_w_mm + 56)
  print_line("=", width = line_w_mm)
  cat("  2. McDONALD-MOFFITT MARGINAL EFFECTS (DECOMPOSITION):\n")
  print_line("=", width = line_w_mm)
  cat(sprintf(paste0("  %-", v_w_mm, "s %18s %18s %18s\n"), 
              "Variable", "d P(y > 0) / dx", "d E(y|y>0) / dx", "Total d E(y)/dx"))
  print_line("-", width = line_w_mm)
  
  for (i in 1:nrow(x$mm_table)) {
    row <- x$mm_table[i, ]
    cat(sprintf(paste0("  %-", v_w_mm, "s %18.*f %18.*f %18.*f\n"),
                row$Variable,
                d, row$Prob_Uncens,
                d, row$Cond_Mean,
                d, row$Uncond_Mean))
  }
  print_line("-", width = line_w_mm)
  cat("  Economic Notes:\n")
  cat("  * Total dE(y)/dx is decomposed into changes on the extensive margin (d P(y>0))\n")
  cat("    and the intensive margin (d E(y|y>0)). Total = P(y>0)*dE(y|y>0) + E(y|y>0)*dP(y>0).\n")
  print_line("=")
  invisible(x)
}

#' @export
coef.econ_tobit <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.econ_tobit <- function(object, ...) {
  object$vcov
}

#' @export
residuals.econ_tobit <- function(object, ...) {
  object$residuals
}

#' @export
fitted.econ_tobit <- function(object, ...) {
  object$fitted.values
}

#' @export
nobs.econ_tobit <- function(object, ...) {
  object$n
}

#' @export
confint.econ_tobit <- function(object, parm, level = 0.95, ...) {
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

#' @export
logLik.econ_tobit <- function(object, ...) {
  val <- object$log_lik
  attr(val, "df") <- object$k + 1
  attr(val, "nobs") <- object$n
  class(val) <- "logLik"
  val
}

#' @export
predict.econ_tobit <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$fitted.values)
  form_no_y <- stats::delete.response(stats::terms(object$formula))
  X_new <- stats::model.matrix(form_no_y, data = newdata)
  b <- stats::coef(object)
  pmax(object$left, as.vector(X_new[, names(b), drop = FALSE] %*% b))
}

#' Count Data Models (Poisson and Negative Binomial Regression)
#'
#' Fits Poisson or Negative Binomial regression for count data,
#' computes Incidence Rate Ratios (IRR) with 95\% confidence intervals,
#' robust standard errors, and conducts Cameron & Trivedi's overdispersion test
#' or Likelihood Ratio test against Poisson.
#'
#' @param formula Model formula (e.g., \code{docvis ~ educ + age + female}).
#' @param data A data.frame containing the variables.
#' @param model Model type: \code{"poisson"} (default) or \code{"negbin"}.
#' @param digits Number of decimals to display (default is 4).
#'
#' @return An object of class \code{econ_count}.
#' @export
#'
#' @examples
#' set.seed(42)
#' n <- 300
#' x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.5)
#' y <- rpois(n, lambda = exp(0.5 + 0.3 * x1 - 0.4 * x2))
#' df <- data.frame(y, x1, x2)
#' count_reg(y ~ x1 + x2, data = df, model = "poisson")
#' count_reg(y ~ x1 + x2, data = df, model = "negbin")
count_reg <- function(formula, data, model = c("poisson", "negbin"), digits = 4) {
  model <- match.arg(model)
  if (missing(data) || is.null(data)) stop("Please provide data.", call. = FALSE)
  
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(formula, data = mf)
  n <- length(y)
  k <- ncol(X)
  
  # Base Poisson GLM
  fam <- stats::poisson(link = "log")
  base_glm <- stats::glm(formula, data = data, family = fam, x = TRUE, y = TRUE)
  b_pois <- stats::coef(base_glm)
  ll_pois <- stats::logLik(base_glm)[1]
  mu_pois <- stats::fitted(base_glm)
  
  # Cameron & Trivedi (1990) Overdispersion Test
  y_disp <- ((y - mu_pois)^2 - y) / mu_pois
  mod_disp <- stats::lm(y_disp ~ 0 + mu_pois)
  t_disp <- summary(mod_disp)$coefficients[1, 3]
  p_disp <- stats::pt(t_disp, df = n - 1, lower.tail = FALSE)
  
  if (model == "poisson") {
    b <- b_pois
    mu_hat <- mu_pois
    ll_model <- ll_pois
    
    # Robust Huber-White standard errors for Poisson
    bread <- safe_solve(crossprod(X * sqrt(mu_hat)))
    scores <- (y - mu_hat) * X
    vcov_mat <- bread %*% crossprod(scores) %*% bread
    dimnames(vcov_mat) <- list(names(b), names(b))
    
    se_b <- sqrt(pmax(0, diag(vcov_mat)))
    z_b <- b / se_b
    p_b <- 2 * stats::pnorm(-abs(z_b))
    
    theta_est <- NA_real_
    lr_nb_stat <- NA_real_
    lr_nb_p <- NA_real_
  } else {
    # -------------------------------------------------------------
    # True Negative Binomial Estimation via Maximum Likelihood
    # -------------------------------------------------------------
    s2_y <- stats::var(y)
    mean_y <- mean(y)
    theta_init <- if (s2_y > mean_y) max(0.05, mean_y^2 / (s2_y - mean_y)) else 1.0
    gamma_init <- log(theta_init)
    
    # Negative Binomial log-likelihood with gamma = log(theta)
    nb_nll <- function(par) {
      beta_v <- par[1:k]
      theta_v <- exp(par[k + 1])
      eta_v <- as.vector(X %*% beta_v)
      mu_v <- exp(pmin(50, eta_v))
      
      ll_terms <- lgamma(y + theta_v) - lgamma(theta_v) - lfactorial(y) + 
                  theta_v * (log(theta_v) - log(theta_v + mu_v)) + 
                  y * (log(pmax(1e-10, mu_v)) - log(theta_v + mu_v))
      -sum(ll_terms)
    }
    
    nb_grad <- function(par) {
      beta_v <- par[1:k]
      theta_v <- exp(par[k + 1])
      eta_v <- as.vector(X %*% beta_v)
      mu_v <- exp(pmin(50, eta_v))
      
      # Gradient wrt beta
      w_resid <- (y - mu_v) / (1 + mu_v / theta_v)
      grad_beta <- -as.vector(crossprod(X, w_resid))
      
      # Gradient wrt gamma
      d_theta <- digamma(y + theta_v) - digamma(theta_v) + log(theta_v / (theta_v + mu_v)) + 1 - (y + theta_v) / (theta_v + mu_v)
      grad_gamma <- -theta_v * sum(d_theta)
      c(grad_beta, grad_gamma)
    }
    
    opt_nb <- stats::optim(par = c(b_pois, gamma_init), fn = nb_nll, gr = nb_grad, method = "BFGS", hessian = TRUE)
    b <- opt_nb$par[1:k]
    names(b) <- colnames(X)
    theta_est <- exp(opt_nb$par[k + 1])
    ll_model <- -opt_nb$value
    mu_hat <- exp(as.vector(X %*% b))
    
    inv_hess <- tryCatch(safe_solve(opt_nb$hessian), error = function(e) diag(1e-4, k + 1))
    vcov_mat <- inv_hess[1:k, 1:k, drop = FALSE]
    dimnames(vcov_mat) <- list(colnames(X), colnames(X))
    
    se_b <- sqrt(pmax(0, diag(vcov_mat)))
    z_b <- b / se_b
    p_b <- 2 * stats::pnorm(-abs(z_b))
    
    # Likelihood Ratio Test vs Poisson
    lr_nb_stat <- max(0, 2 * (ll_model - ll_pois))
    lr_nb_p <- stats::pchisq(lr_nb_stat, df = 1, lower.tail = FALSE)
  }
  
  # Incidence Rate Ratios: IRR = exp(b)
  irr <- exp(b)
  irr_low <- exp(b - 1.96 * se_b)
  irr_high <- exp(b + 1.96 * se_b)
  
  coef_table <- data.frame(
    Estimate  = b,
    Std_Error = se_b,
    z_stat    = z_b,
    p_value   = p_b,
    Signif    = sapply(p_b, get_stars),
    IRR       = irr,
    IRR_Lower = irr_low,
    IRR_Upper = irr_high,
    stringsAsFactors = FALSE
  )
  
  res <- list(
    formula = formula,
    model = model,
    coefficients = b,
    coef_table = coef_table,
    vcov = vcov_mat,
    residuals = y - mu_hat,
    fitted.values = mu_hat,
    y = y,
    X = X,
    n = n,
    k = k,
    theta = theta_est,
    disp_stat = t_disp,
    disp_p = p_disp,
    lr_nb_stat = lr_nb_stat,
    lr_nb_p = lr_nb_p,
    log_lik = ll_model,
    digits = digits,
    glm_model = base_glm
  )
  class(res) <- "econ_count"
  return(res)
}

#' @export
print.econ_count <- function(x, digits = NULL, ...) {
  d <- if (is.null(digits)) x$digits else digits
  
  print_line("=")
  if (x$model == "poisson") {
    cat("                  COUNT DATA REGRESSION: POISSON                        \n")
  } else {
    cat("             COUNT DATA REGRESSION: NEGATIVE BINOMIAL                   \n")
  }
  cat(paste0("  Formula: ", deparse(x$formula), "\n"))
  cat(sprintf("  Sample Size: %-12d Log-Likelihood: %.2f\n", x$n, x$log_lik))
  if (x$model == "negbin") {
    cat(sprintf("  Dispersion Theta: %-11.4f Alpha (1/Theta): %.4f\n", x$theta, 1 / x$theta))
  }
  cat(paste0("  Standard Errors: ", if (x$model == "poisson") "Robust (Huber-White Sandwich)" else "Standard (MLE Hessian)", "\n"))
  print_line("=")
  
  v_w <- max(16, max(nchar(rownames(x$coef_table))))
  line_w <- max(74, v_w + 64)
  cat("  1. LOG-LINEAR COEFFICIENTS & INCIDENCE RATE RATIOS (IRR):\n")
  print_line("-", width = line_w)
  cat(sprintf(paste0("  %-", v_w, "s %10s %10s %8s %10s %10s %16s\n"), 
              "Variable", "Estimate", "Std. Error", "z-stat", "p-value", "IRR", "95% IRR CI"))
  print_line("-", width = line_w)
  
  for (i in 1:nrow(x$coef_table)) {
    vname <- rownames(x$coef_table)[i]
    row <- x$coef_table[i, ]
    ci_str <- sprintf("[%.*f, %.*f]", d, row$IRR_Lower, d, row$IRR_Upper)
    cat(sprintf(paste0("  %-", v_w, "s %10.*f %10.*f %8.2f %10s %10.*f %16s %s\n"),
                vname,
                d, row$Estimate,
                d, row$Std_Error,
                row$z_stat,
                format_pval(row$p_value),
                d, row$IRR,
                ci_str,
                row$Signif))
  }
  print_line("-", width = line_w)
  cat("  Signif. codes:  *** p < 0.001,  ** p < 0.01,  * p < 0.05,  . p < 0.1\n\n")
  
  # Diagnostic Card
  print_line("=")
  if (x$model == "poisson") {
    cat("                 CAMERON & TRIVEDI OVERDISPERSION TEST                  \n")
    print_line("=")
    cat(sprintf("  H0: Equidispersion (Var(y) = E(y)) -> Poisson is adequate\n"))
    cat(sprintf("  H1: Overdispersion (Var(y) > E(y)) -> Negative Binomial is required\n"))
    cat(sprintf("  Test Statistic:  t = %10.*f\n", d, x$disp_stat))
    cat(sprintf("  p-value:         %10s\n", format_pval(x$disp_p)))
    print_line("-")
    verdict <- if (x$disp_p < 0.05) {
      "REJECT H0 (p < 0.05) -> Significant Overdispersion detected (Use Negative Binomial)"
    } else {
      "FAIL TO REJECT H0 (p >= 0.05) -> Equidispersion holds; Poisson model is appropriate"
    }
    cat(sprintf("  Conclusion:      %s\n", verdict))
  } else {
    cat("             LIKELIHOOD RATIO TEST: NEGATIVE BINOMIAL VS POISSON        \n")
    print_line("=")
    cat(sprintf("  H0: Alpha = 0 (Poisson equidispersion holds)\n"))
    cat(sprintf("  H1: Alpha > 0 (Overdispersion requires Negative Binomial)\n"))
    cat(sprintf("  LR Chi-Sq (1):   %10.*f\n", d, x$lr_nb_stat))
    cat(sprintf("  p-value:         %10s\n", format_pval(x$lr_nb_p)))
    print_line("-")
    verdict <- if (!is.na(x$lr_nb_p) && x$lr_nb_p < 0.05) {
      "REJECT H0 (p < 0.05) -> Significant overdispersion; Negative Binomial preferred"
    } else {
      "FAIL TO REJECT H0 -> No significant overdispersion; Poisson is sufficient"
    }
    cat(sprintf("  Conclusion:      %s\n", verdict))
  }
  print_line("=")
  invisible(x)
}

#' @export
coef.econ_count <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.econ_count <- function(object, ...) {
  object$vcov
}

#' @export
residuals.econ_count <- function(object, ...) {
  object$residuals
}

#' @export
fitted.econ_count <- function(object, ...) {
  object$fitted.values
}

#' @export
nobs.econ_count <- function(object, ...) {
  object$n
}

#' @export
confint.econ_count <- function(object, parm, level = 0.95, ...) {
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

#' @export
logLik.econ_count <- function(object, ...) {
  val <- object$log_lik
  attr(val, "df") <- object$k + if (object$model == "negbin") 1 else 0
  attr(val, "nobs") <- object$n
  class(val) <- "logLik"
  val
}

#' @export
predict.econ_count <- function(object, newdata = NULL, type = c("response", "link"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    if (type == "response") return(object$fitted.values)
    return(log(object$fitted.values))
  }
  form_no_y <- stats::delete.response(stats::terms(object$formula))
  X_new <- stats::model.matrix(form_no_y, data = newdata)
  b <- stats::coef(object)
  eta <- as.vector(X_new[, names(b), drop = FALSE] %*% b)
  if (type == "response") exp(eta) else eta
}
