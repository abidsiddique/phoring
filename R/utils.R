#' Format p-values for clear display
#' @noRd
format_pval <- function(p, digits = 4) {
  if (is.na(p)) return("    NA")
  if (p < 0.0001) {
    return("< 0.0001")
  } else {
    return(sprintf(paste0("%", digits + 2, ".", digits, "f"), p))
  }
}

#' Significance stars for p-values
#' @noRd
get_stars <- function(p) {
  if (is.na(p)) return("   ")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("** ")
  if (p < 0.05)  return("*  ")
  if (p < 0.1)   return(".  ")
  return("   ")
}

#' Format numbers cleanly with fixed decimals
#' @noRd
format_num <- function(x, digits = 4, width = 10) {
  if (is.na(x)) return(sprintf(paste0("%", width, "s"), "NA"))
  sprintf(paste0("%", width, ".", digits, "f"), x)
}

#' Print a horizontal divider line
#' @noRd
print_line <- function(char = "=", width = 74) {
  cat(paste(rep(char, width), collapse = ""), "\n", sep = "")
}

#' Safe matrix solver with rank check
#' @noRd
safe_solve <- function(A) {
  tryCatch({
    solve(A)
  }, error = function(err) {
    qr_A <- qr(A)
    if (qr_A$rank < ncol(A)) {
      stop("Matrix inversion failed: system is singular or rank-deficient.", call. = FALSE)
    }
    solve(qr_A)
  })
}

#' Check for exact multicollinearity and dummy variable trap
#' @noRd
check_collinearity <- function(X) {
  qr_x <- qr(X)
  k <- ncol(X)
  if (qr_x$rank < k) {
    piv <- qr_x$pivot
    rnk <- qr_x$rank
    bad_idx <- piv[(rnk + 1):k]
    bad_cols <- colnames(X)[bad_idx]
    if (length(bad_cols) == 0) bad_cols <- paste0("Column_", bad_idx)
    stop(paste0("Exact multicollinearity (dummy variable trap) detected!\n",
                "Model matrix is rank-deficient (rank ", rnk, " < ", k, " columns).\n",
                "The following predictor(s) are perfectly collinear with other variables:\n  * ",
                paste(bad_cols, collapse = ", "), "\n",
                "Suggestion: Drop one of the redundant categories or collinear variables."),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Calculate Heteroskedasticity-Consistent Covariance Matrix (HC0, HC1, HC2, HC3) with optional weights
#' @noRd
calc_vcov_hc <- function(X, e, type = c("HC1", "HC0", "HC2", "HC3"), weights = NULL) {
  type <- match.arg(type)
  n <- nrow(X)
  k <- ncol(X)
  
  if (!is.null(weights)) {
    w_sqrt <- sqrt(weights)
    Xw <- X * w_sqrt
    bread <- safe_solve(crossprod(Xw))
    eff_e <- e * weights
  } else {
    bread <- safe_solve(crossprod(X))
    eff_e <- e
  }
  
  if (type == "HC0") {
    meat <- crossprod(X * eff_e)
    vcov <- bread %*% meat %*% bread
  } else if (type == "HC1") {
    dfc <- n / (n - k)
    meat <- crossprod(X * eff_e)
    vcov <- dfc * (bread %*% meat %*% bread)
  } else if (type %in% c("HC2", "HC3")) {
    if (!is.null(weights)) {
      h <- rowSums((Xw %*% bread) * Xw)
    } else {
      h <- rowSums((X %*% bread) * X)
    }
    h <- pmin(h, 0.9999) # prevent division by zero
    
    if (type == "HC2") {
      w_adj <- eff_e / sqrt(1 - h)
    } else {
      w_adj <- eff_e / (1 - h)
    }
    meat <- crossprod(X * w_adj)
    vcov <- bread %*% meat %*% bread
  }
  
  dimnames(vcov) <- list(colnames(X), colnames(X))
  return(vcov)
}

#' Calculate Cluster-Robust Covariance Matrix with optional weights
#' @noRd
calc_vcov_cluster <- function(X, e, cluster, weights = NULL) {
  if (length(cluster) != nrow(X)) {
    stop("Length of `cluster` variable must match number of observations.", call. = FALSE)
  }
  
  cluster <- as.factor(cluster)
  G <- nlevels(cluster)
  n <- nrow(X)
  k <- ncol(X)
  
  if (G < 2) {
    stop("Clustering variable must contain at least 2 distinct clusters.", call. = FALSE)
  }
  
  if (!is.null(weights)) {
    bread <- safe_solve(crossprod(X * sqrt(weights)))
    scores <- X * (e * weights)
  } else {
    bread <- safe_solve(crossprod(X))
    scores <- X * e
  }
  
  cluster_scores <- rowsum(scores, cluster)
  meat <- crossprod(cluster_scores)
  
  # Standard Stata / Cameron-Miller finite-sample degrees-of-freedom correction
  dfc <- (G / (G - 1)) * ((n - 1) / (n - k))
  vcov <- dfc * (bread %*% meat %*% bread)
  
  dimnames(vcov) <- list(colnames(X), colnames(X))
  return(vcov)
}

#' Calculate Newey-West HAC Covariance Matrix for Time Series with optional weights
#' @noRd
calc_vcov_hac <- function(X, e, lag = NULL, weights = NULL) {
  n <- nrow(X)
  k <- ncol(X)
  
  if (is.null(lag)) {
    # Standard rule of thumb for lag truncation (Stock & Watson / Greene)
    lag <- max(1, floor(4 * (n / 100)^(2 / 9)))
  }
  
  if (!is.null(weights)) {
    bread <- safe_solve(crossprod(X * sqrt(weights)))
    scores <- X * (e * weights)
  } else {
    bread <- safe_solve(crossprod(X))
    scores <- X * e
  }
  
  meat <- crossprod(scores)
  
  # Autocovariance terms with Bartlett kernel weights
  if (lag > 0) {
    for (l in 1:lag) {
      weight <- 1 - (l / (lag + 1))
      gamma_l <- crossprod(scores[(l + 1):n, , drop = FALSE], scores[1:(n - l), , drop = FALSE])
      meat <- meat + weight * (gamma_l + t(gamma_l))
    }
  }
  
  dfc <- n / (n - k)
  vcov <- dfc * (bread %*% meat %*% bread)
  
  dimnames(vcov) <- list(colnames(X), colnames(X))
  return(vcov)
}

#' Calculate robust Wald F-statistic for overall regression
#' @noRd
calc_robust_f <- function(b, vcov_mat, df_denom, has_intercept = TRUE) {
  k <- length(b)
  if (has_intercept && "(Intercept)" %in% names(b)) {
    test_idx <- which(names(b) != "(Intercept)")
  } else {
    test_idx <- seq_along(b)
  }
  
  q <- length(test_idx)
  if (q == 0) {
    return(list(f_stat = NA_real_, p_value = NA_real_, df1 = 0, df2 = df_denom))
  }
  
  b_sub <- b[test_idx]
  v_sub <- vcov_mat[test_idx, test_idx, drop = FALSE]
  
  inv_v <- tryCatch(safe_solve(v_sub), error = function(e) NULL)
  if (is.null(inv_v)) {
    return(list(f_stat = NA_real_, p_value = NA_real_, df1 = q, df2 = df_denom))
  }
  
  wald_stat <- as.numeric(t(b_sub) %*% inv_v %*% b_sub)
  f_stat <- max(0, wald_stat / q)
  p_val <- stats::pf(f_stat, df1 = q, df2 = df_denom, lower.tail = FALSE)
  
  list(f_stat = f_stat, p_value = p_val, df1 = q, df2 = df_denom)
}

#' MacKinnon Continuous Response Surface p-value for Unit Root & Cointegration Tests
#' @noRd
calc_mackinnon_p <- function(stat, n = 100, type = c("trend", "drift", "none"), test = c("unitroot", "coint"), k = 2) {
  type <- match.arg(type)
  test <- match.arg(test)
  
  # Asymptotic surface response interpolation based on MacKinnon (1994, 1996)
  if (test == "unitroot") {
    if (type == "trend") {
      # ADF with constant and trend (mu and t)
      # 1% = -3.96, 2.5% = -3.66, 5% = -3.41, 10% = -3.13, 50% = -2.18, 90% = -1.25, 95% = -0.94
      knots_stat <- c(-5.5, -4.38, -3.96, -3.66, -3.41, -3.13, -2.57, -2.18, -1.62, -1.25, -0.94, 0.5)
      knots_p    <- c(0.0001, 0.001, 0.01, 0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.9999)
    } else if (type == "drift") {
      # ADF with constant (mu)
      # 1% = -3.43, 2.5% = -3.12, 5% = -2.86, 10% = -2.57, 50% = -1.57, 90% = -0.44
      knots_stat <- c(-5.0, -3.96, -3.43, -3.12, -2.86, -2.57, -2.04, -1.57, -0.98, -0.44, -0.07, 1.5)
      knots_p    <- c(0.0001, 0.001, 0.01, 0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.9999)
    } else {
      # ADF with no constant (none)
      # 1% = -2.56, 2.5% = -2.23, 5% = -1.94, 10% = -1.62, 50% = -0.52, 90% = 0.89
      knots_stat <- c(-4.0, -3.14, -2.56, -2.23, -1.94, -1.62, -1.04, -0.52, 0.17, 0.89, 1.33, 3.0)
      knots_p    <- c(0.0001, 0.001, 0.01, 0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.9999)
    }
  } else {
    # Cointegration (Engle-Granger on residuals of k variables)
    # k = total variables in cointegrating equation
    crit_1 <- -3.90 - 0.42 * (k - 2)
    crit_5 <- -3.34 - 0.38 * (k - 2)
    crit_10 <- -3.04 - 0.35 * (k - 2)
    knots_stat <- c(crit_1 - 2.0, crit_1 - 0.6, crit_1, (crit_1 + crit_5)/2, crit_5, crit_10, crit_10 + 0.6, crit_10 + 1.2, 0.0)
    knots_p    <- c(0.0001, 0.001, 0.01, 0.025, 0.05, 0.10, 0.25, 0.50, 0.9999)
  }
  
  # Smooth monotonic spline / log-odds interpolation
  log_p <- log(knots_p / (1 - knots_p))
  est_log_p <- stats::approx(x = knots_stat, y = log_p, xout = stat, rule = 2)$y
  p_val <- 1 / (1 + exp(-est_log_p))
  pmin(pmax(p_val, 0.0001), 0.9999)
}
