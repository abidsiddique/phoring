#' Binary Choice Models (Logit and Probit) with Marginal Effects
#'
#' Fits logistic regression (Logit) or probit regression (Probit) models,
#' computes robust Huber-White or clustered standard errors, reports odds ratios (for logit),
#' ROC-AUC, classification accuracy, and calculates Average Marginal Effects (AME)
#' with delta-method standard errors.
#'
#' @param formula Model formula (e.g., \code{admit ~ gre + gpa + rank}).
#' @param data A data.frame containing the variables.
#' @param se Type of standard errors: \code{"robust"} (Huber-White sandwich, default),
#'   \code{"standard"} (classical MLE), or \code{"cluster"} (cluster-robust).
#' @param cluster Vector or column name for cluster-robust standard errors.
#' @param digits Number of decimal places to display (default is 4).
#' @param ... Additional arguments passed to \code{\link[stats]{glm}}.
#'
#' @return An object of class \code{econ_logit} containing index estimates, odds ratios,
#'   average marginal effects (AME), and model diagnostics.
#' @export
#'
#' @examples
#' data(mtcars)
#' logit(am ~ mpg + hp, data = mtcars)
#' probit(am ~ mpg + hp, data = mtcars)
logit <- function(formula, data, se = c("robust", "standard", "cluster"), cluster = NULL, digits = 4, ...) {
  fit_binary(formula = formula, data = data, link = "logit", se = se, cluster = cluster, digits = digits, ...)
}

#' @rdname logit
#' @export
probit <- function(formula, data, se = c("robust", "standard", "cluster"), cluster = NULL, digits = 4, ...) {
  fit_binary(formula = formula, data = data, link = "probit", se = se, cluster = cluster, digits = digits, ...)
}

#' Internal fitting function for binary choice models
#' @noRd
fit_binary <- function(formula, data, link = c("logit", "probit"), 
                       se = c("robust", "standard", "cluster"), cluster = NULL, digits = 4, ...) {
  link <- match.arg(link)
  
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
  
  # Resolve cluster column if provided
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
  
  # Fit base GLM
  fam <- stats::binomial(link = link)
  base_glm <- stats::glm(formula, data = data, family = fam, x = TRUE, y = TRUE, ...)
  
  X <- base_glm$x
  y <- base_glm$y
  b <- stats::coef(base_glm)
  n <- nrow(X)
  k <- ncol(X)
  df_resid <- n - k
  
  if (df_resid <= 0) {
    stop("Degrees of freedom <= 0.", call. = FALSE)
  }
  
  # Handle subsetted cluster vector if missing rows were dropped
  na_action <- stats::na.action(base_glm)
  if (!is.null(na_action) && !is.null(cluster_vec)) {
    cluster_vec <- cluster_vec[-na_action]
  }
  
  eta <- as.vector(X %*% b)
  p_hat <- stats::fitted(base_glm)
  
  # Universal exact GLM working weights and bread
  working_w <- stats::weights(base_glm, type = "working")
  working_res <- stats::residuals(base_glm, type = "working")
  bread <- safe_solve(crossprod(X * sqrt(working_w)))
  
  # Universal exact GLM score matrix: (y - mu) * (d eta / d mu) * w * X = working_res * working_w * X
  scores <- (working_res * working_w) * X
  
  df_inference <- df_resid
  # Standard errors
  if (se == "robust") {
    meat <- crossprod(scores)
    dfc <- n / (n - k)
    vcov_mat <- dfc * (bread %*% meat %*% bread)
    se_label <- "Robust (Huber-White Sandwich)"
  } else if (se == "cluster") {
    if (is.null(cluster_vec)) {
      stop("Please specify the `cluster` argument when `se = 'cluster'`.", call. = FALSE)
    }
    cluster_fac <- as.factor(cluster_vec)
    G <- nlevels(cluster_fac)
    df_inference <- max(1, G - 1)
    clust_scores <- rowsum(scores, cluster_fac)
    meat <- crossprod(clust_scores)
    dfc <- (G / (G - 1)) * ((n - 1) / (n - k))
    vcov_mat <- dfc * (bread %*% meat %*% bread)
    se_label <- paste0("Clustered (by ", if (is.character(cluster)) cluster else "group", ", G = ", G, ")")
  } else {
    vcov_mat <- bread
    se_label <- "Standard (MLE)"
  }
  dimnames(vcov_mat) <- list(colnames(X), colnames(X))
  
  se_b <- sqrt(pmax(0, diag(vcov_mat)))
  z_b <- b / se_b
  p_b <- 2 * stats::pnorm(-abs(z_b))
  
  coef_table <- data.frame(
    Estimate  = b,
    Std_Error = se_b,
    z_stat    = z_b,
    p_value   = p_b,
    Signif    = sapply(p_b, get_stars),
    stringsAsFactors = FALSE
  )
  
  # Model fit statistics
  ll_full <- stats::logLik(base_glm)[1]
  null_glm <- stats::glm(y ~ 1, family = fam)
  ll_null <- stats::logLik(null_glm)[1]
  
  pseudo_r2 <- 1 - (ll_full / ll_null)
  lr_stat <- 2 * (ll_full - ll_null)
  lr_p <- stats::pchisq(lr_stat, df = k - 1, lower.tail = FALSE)
  accuracy <- mean((p_hat >= 0.5) == y) * 100
  
  # ROC-AUC (Area Under the Curve)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  auc <- if (n1 > 0 && n0 > 0) {
    (sum(rank(p_hat)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  } else {
    NA
  }
  
  # Odds Ratios (for Logit)
  odds_ratio <- if (link == "logit") exp(b) else NULL
  or_ci_low  <- if (link == "logit") exp(b - 1.96 * se_b) else NULL
  or_ci_high <- if (link == "logit") exp(b + 1.96 * se_b) else NULL
  
  # -------------------------------------------------------------
  # AVERAGE MARGINAL EFFECTS (AME) & DELTA METHOD
  # -------------------------------------------------------------
  preds <- colnames(X)
  preds <- preds[preds != "(Intercept)"]
  
  calc_ame <- function(beta_vec) {
    eta_v <- as.vector(X %*% beta_vec)
    ame_vec <- numeric(length(preds))
    names(ame_vec) <- preds
    
    for (i in seq_along(preds)) {
      vname <- preds[i]
      col_val <- X[, vname]
      is_dummy <- all(col_val %in% c(0, 1))
      
      if (is_dummy) {
        # Discrete change from 0 to 1
        X_0 <- X; X_0[, vname] <- 0
        X_1 <- X; X_1[, vname] <- 1
        eta_0 <- as.vector(X_0 %*% beta_vec)
        eta_1 <- as.vector(X_1 %*% beta_vec)
        if (link == "logit") {
          p0 <- stats::plogis(eta_0)
          p1 <- stats::plogis(eta_1)
        } else {
          p0 <- stats::pnorm(eta_0)
          p1 <- stats::pnorm(eta_1)
        }
        ame_vec[vname] <- mean(p1 - p0)
      } else {
        # Continuous derivative
        if (link == "logit") {
          p_v <- stats::plogis(eta_v)
          d_v <- p_v * (1 - p_v)
        } else {
          d_v <- stats::dnorm(eta_v)
        }
        ame_vec[vname] <- mean(d_v * beta_vec[vname])
      }
    }
    return(ame_vec)
  }
  
  ame_values <- calc_ame(b)
  
  # Numerical Jacobian for Delta Method SE: J = d AME / d beta
  eps <- 1e-6
  J <- matrix(0, nrow = length(preds), ncol = length(b))
  rownames(J) <- preds
  colnames(J) <- names(b)
  
  for (j in seq_along(b)) {
    b_step_pos <- b; b_step_pos[j] <- b[j] + eps
    b_step_neg <- b; b_step_neg[j] <- b[j] - eps
    J[, j] <- (calc_ame(b_step_pos) - calc_ame(b_step_neg)) / (2 * eps)
  }
  
  vcov_ame <- J %*% vcov_mat %*% t(J)
  se_ame <- sqrt(pmax(0, diag(vcov_ame)))
  z_ame <- ame_values / se_ame
  p_ame <- 2 * stats::pnorm(-abs(z_ame))
  
  ame_table <- data.frame(
    AME       = ame_values,
    Std_Error = se_ame,
    z_stat    = z_ame,
    p_value   = p_ame,
    Signif    = sapply(p_ame, get_stars),
    stringsAsFactors = FALSE
  )
  
  res <- list(
    formula = formula,
    link = link,
    coefficients = b,
    coef_table = coef_table,
    ame_table = ame_table,
    vcov = vcov_mat,
    vcov_ame = vcov_ame,
    residuals = stats::residuals(base_glm, type = "response"),
    fitted.values = p_hat,
    n = n,
    k = k,
    df_residual = df_resid,
    df_inference = df_inference,
    pseudo_r2 = pseudo_r2,
    log_lik = ll_full,
    lr_stat = lr_stat,
    lr_p = lr_p,
    accuracy = accuracy,
    auc = auc,
    odds_ratio = odds_ratio,
    or_ci_low = or_ci_low,
    or_ci_high = or_ci_high,
    se_type = se_label,
    digits = digits,
    glm_model = base_glm
  )
  class(res) <- "econ_logit"
  return(res)
}

#' @export
print.econ_logit <- function(x, digits = NULL, ...) {
  d <- if (is.null(digits)) x$digits else digits
  model_name <- if (x$link == "logit") "LOGISTIC REGRESSION (LOGIT)" else "PROBIT REGRESSION (PROBIT)"
  
  print_line("=")
  cat(sprintf("                 %s\n", model_name))
  cat(paste0("  Formula: ", deparse(x$formula), "\n"))
  cat(sprintf("  Sample Size: %-13d McFadden Pseudo R2: %-10.4f Accuracy: %.1f%%\n", 
              x$n, x$pseudo_r2, x$accuracy))
  cat(sprintf("  ROC-AUC:     %-13.4f LR Chi-Sq (%d):      %-10.2f p-value:  %s\n", 
              x$auc, x$k - 1, x$lr_stat, format_pval(x$lr_p)))
  cat(paste0("  Standard Errors: ", x$se_type, "\n"))
  print_line("=")
  
  v_w <- max(16, max(nchar(rownames(x$coef_table))))
  line_w <- max(74, v_w + 58)
  
  cat("  1. RAW INDEX COEFFICIENTS (LATENT UTILITY):\n")
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
  
  # Odds Ratios (for Logit)
  if (x$link == "logit" && !is.null(x$odds_ratio)) {
    v_w_or <- max(16, max(nchar(names(x$odds_ratio))))
    line_w_or <- max(74, v_w_or + 48)
    print_line("=", width = line_w_or)
    cat("  2. ODDS RATIOS (EXPONENTIATED COEFFICIENTS) & 95% CI:\n")
    print_line("=", width = line_w_or)
    cat(sprintf(paste0("  %-", v_w_or, "s %12s %16s %16s\n"), "Variable", "Odds Ratio", "95% CI Lower", "95% CI Upper"))
    print_line("-", width = line_w_or)
    for (i in 1:length(x$odds_ratio)) {
      vname <- names(x$odds_ratio)[i]
      if (vname != "(Intercept)") {
        cat(sprintf(paste0("  %-", v_w_or, "s %12.*f %16.*f %16.*f\n"),
                    vname,
                    d, x$odds_ratio[i],
                    d, x$or_ci_low[i],
                    d, x$or_ci_high[i]))
      }
    }
    print_line("-", width = line_w_or)
    cat("  Interpretation: An odds ratio > 1 indicates increased odds of Y = 1.\n\n")
  }
  
  # Marginal Effects Table
  sec_num <- if (x$link == "logit") "3" else "2"
  v_w_ame <- max(16, max(nchar(rownames(x$ame_table))))
  line_w_ame <- max(74, v_w_ame + 58)
  print_line("=", width = line_w_ame)
  cat(sprintf("  %s. AVERAGE MARGINAL EFFECTS (AME) - INTERPRETABLE PROBABILITIES:\n", sec_num))
  print_line("=", width = line_w_ame)
  cat(sprintf(paste0("  %-", v_w_ame, "s %10s %12s %10s %10s   %s\n"), 
              "Variable", "AME", "Std. Error", "z-stat", "p-value", "Signif"))
  print_line("-", width = line_w_ame)
  
  for (i in 1:nrow(x$ame_table)) {
    vname <- rownames(x$ame_table)[i]
    row <- x$ame_table[i, ]
    cat(sprintf(paste0("  %-", v_w_ame, "s %10.*f %12.*f %10.3f %10s   %s\n"),
                vname,
                d, row$AME,
                d, row$Std_Error,
                row$z_stat,
                format_pval(row$p_value),
                row$Signif))
  }
  print_line("-", width = line_w_ame)
  cat("  Econometric Note:\n")
  cat("  * An AME represents the average percentage-point change in Pr(Y = 1)\n")
  cat("    associated with a 1-unit increase in the predictor (or 0 -> 1 for dummy).\n")
  print_line("=")
  invisible(x)
}

#' @export
summary.econ_logit <- function(object, ...) {
  object
}

#' @export
coef.econ_logit <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.econ_logit <- function(object, ...) {
  object$vcov
}

#' @export
residuals.econ_logit <- function(object, ...) {
  object$residuals
}

#' @export
fitted.econ_logit <- function(object, ...) {
  object$fitted.values
}

#' @export
nobs.econ_logit <- function(object, ...) {
  object$n
}

#' @export
confint.econ_logit <- function(object, parm, level = 0.95, ...) {
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
predict.econ_logit <- function(object, newdata = NULL, type = c("response", "link"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    if (type == "response") return(object$fitted.values)
    return(stats::predict(object$glm_model, type = "link"))
  }
  stats::predict(object$glm_model, newdata = newdata, type = type, ...)
}

#' @export
logLik.econ_logit <- function(object, ...) {
  stats::logLik(object$glm_model, ...)
}

#' @export
AIC.econ_logit <- function(object, ..., k = 2) {
  stats::AIC(object$glm_model, ..., k = k)
}

#' @export
BIC.econ_logit <- function(object, ...) {
  stats::BIC(object$glm_model, ...)
}
