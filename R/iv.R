#' Instrumental Variables (2SLS) Regression with Automated Diagnostics
#'
#' Fits a Two-Stage Least Squares (2SLS) instrumental variables regression model,
#' calculates robust or clustered standard errors, and automatically runs diagnostic tests
#' (Weak Instruments First-Stage F, Wu-Hausman Endogeneity, and Sargan Overidentifying Restrictions).
#'
#' @param formula Two-part formula separated by \code{|}. Supports standard forms:
#'   \code{y ~ exog + endo | exog + iv} or \code{y ~ exog | endo ~ iv}.
#'   Example: \code{wage ~ exper + educ | exper + mother_educ + father_educ}.
#' @param data A data.frame containing the variables in the model.
#' @param se Type of standard errors to compute: \code{"robust"} (White HC1, default),
#'   \code{"standard"} (classical 2SLS), or \code{"cluster"} (cluster-robust).
#' @param cluster Vector or column name indicating cluster grouping.
#' @param digits Number of decimal places to display in output tables (default is 4).
#' @param ... Additional arguments.
#'
#' @return An object of class \code{econ_iv} containing 2SLS estimates, diagnostics,
#'   and standard errors.
#' @export
#'
#' @examples
#' data(wage_data, package = "phoring")
#' iv(wage ~ exper + educ | exper + mother_educ, data = wage_data)
#' iv(wage ~ exper | educ ~ mother_educ + father_educ, data = wage_data)
iv <- function(formula, data, se = c("robust", "standard", "cluster"), cluster = NULL, digits = 4, ...) {
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
  
  # Parse two-part formula separated by |
  form_str <- deparse(formula)
  form_str <- paste(form_str, collapse = " ")
  parts <- strsplit(form_str, "\\|")[[1]]
  
  if (length(parts) != 2) {
    stop("`iv()` formula must contain '|' separating regressors and instruments.\nExample: wage ~ exper + educ | exper + mother_educ + father_educ", call. = FALSE)
  }
  
  part1 <- trimws(parts[1])
  part2 <- trimws(parts[2])
  
  if (grepl("~", part2)) {
    # Syntax style: y ~ exog | endo ~ iv
    f1 <- stats::as.formula(part1)
    f2 <- stats::as.formula(part2)
    
    y_var <- all.vars(f1)[1]
    exog_terms <- attr(stats::terms(f1), "term.labels")
    endo_terms <- all.vars(f2)[1]
    iv_terms <- attr(stats::terms(f2), "term.labels")
    
    # Construct unified regressor and instrument formulas
    x_rhs <- if (length(exog_terms) > 0) paste(c(exog_terms, endo_terms), collapse = " + ") else endo_terms
    z_rhs <- if (length(exog_terms) > 0) paste(c(exog_terms, iv_terms), collapse = " + ") else paste(iv_terms, collapse = " + ")
    
    form_x <- stats::as.formula(paste(y_var, "~", x_rhs))
    form_z <- stats::as.formula(paste("~", z_rhs))
    
    endo_vars <- endo_terms
    iv_vars <- iv_terms
    exog_vars <- exog_terms
  } else {
    # Syntax style: y ~ exog + endo | exog + iv
    form_x <- stats::as.formula(part1)
    form_z <- stats::as.formula(paste("~", part2))
    
    vars_x <- all.vars(form_x)
    y_var <- vars_x[1]
    all_x <- setdiff(vars_x, y_var)
    all_z <- all.vars(form_z)
    
    exog_vars <- intersect(all_x, all_z)
    endo_vars <- setdiff(all_x, all_z)
    iv_vars <- setdiff(all_z, all_x)
  }
  
  if (length(endo_vars) == 0) {
    stop("No endogenous variables identified in formula. Use `reg()` for standard OLS.", call. = FALSE)
  }
  if (length(iv_vars) == 0) {
    stop("No excluded instruments identified in formula.", call. = FALSE)
  }
  if (length(iv_vars) < length(endo_vars)) {
    stop(sprintf("Model is underidentified! %d instruments provided for %d endogenous variables.",
                 length(iv_vars), length(endo_vars)), call. = FALSE)
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
  
  # Construct complete data frame using model.frame
  all_needed <- unique(c(all.vars(form_x), all.vars(form_z)))
  missing_vars <- setdiff(all_needed, names(data))
  if (length(missing_vars) > 0) {
    stop(paste0("Variables not found in data: ", paste(missing_vars, collapse = ", ")), call. = FALSE)
  }
  
  if (!is.null(cluster_vec)) {
    data$._cluster_id_. <- cluster_vec
    sub_data <- stats::na.omit(data[, c(all_needed, "._cluster_id_."), drop = FALSE])
    cluster_vec <- sub_data$._cluster_id_.
    sub_data$._cluster_id_. <- NULL
  } else {
    sub_data <- stats::na.omit(data[, all_needed, drop = FALSE])
  }
  
  n <- nrow(sub_data)
  
  # Build design matrices using model.frame and model.matrix preserving transformations
  mf_x <- stats::model.frame(form_x, data = sub_data)
  y <- stats::model.response(mf_x)
  X <- stats::model.matrix(form_x, data = mf_x)
  
  mf_z <- stats::model.frame(form_z, data = sub_data)
  Z <- stats::model.matrix(form_z, data = mf_z)
  
  k <- ncol(X)
  l <- ncol(Z)
  df_resid <- n - k
  
  if (df_resid <= 0) {
    stop("Degrees of freedom <= 0. Not enough observations.", call. = FALSE)
  }
  
  # First Stage Projections: X_hat = Z (Z'Z)^(-1) Z' X
  qr_Z <- qr(Z)
  if (qr_Z$rank < l) {
    stop("Instrument matrix Z is rank-deficient (collinear instruments).", call. = FALSE)
  }
  X_hat <- qr.fitted(qr_Z, X)
  colnames(X_hat) <- colnames(X)
  
  # 2SLS Estimation: beta_2sls = (X_hat' X_hat)^(-1) X_hat' y
  bread <- safe_solve(crossprod(X_hat))
  b <- as.vector(bread %*% crossprod(X_hat, y))
  names(b) <- colnames(X)
  
  # Structural Residuals (Crucial: using actual X, NOT X_hat)
  e <- as.vector(y - X %*% b)
  s2 <- sum(e^2) / df_resid
  
  # Covariance Matrix
  df_inference <- df_resid
  if (se == "robust") {
    dfc <- n / (n - k)
    meat <- crossprod(X_hat * e)
    vcov_mat <- dfc * (bread %*% meat %*% bread)
    se_label <- "Robust (White HC1)"
  } else if (se == "cluster") {
    if (is.null(cluster_vec)) {
      stop("Please specify `cluster` argument when `se = 'cluster'`.", call. = FALSE)
    }
    cluster_fac <- as.factor(cluster_vec)
    G <- nlevels(cluster_fac)
    df_inference <- max(1, G - 1)
    clust_scores <- rowsum(X_hat * e, cluster_fac)
    meat <- crossprod(clust_scores)
    dfc <- (G / (G - 1)) * ((n - 1) / (n - k))
    vcov_mat <- dfc * (bread %*% meat %*% bread)
    se_label <- paste0("Clustered (by ", if (is.character(cluster)) cluster else "group", ", G = ", G, ")")
  } else {
    vcov_mat <- s2 * bread
    se_label <- "Standard 2SLS"
  }
  dimnames(vcov_mat) <- list(colnames(X), colnames(X))
  
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
  
  # R-squared for 2SLS (based on structural residuals)
  tss <- sum((y - mean(y))^2)
  rss <- sum(e^2)
  r2 <- max(0, 1 - (rss / tss))
  adj_r2 <- max(0, 1 - ((rss / df_resid) / (tss / (n - 1))))
  
  # -------------------------------------------------------------
  # Automated Diagnostic Battery for 2SLS
  # -------------------------------------------------------------
  first_stage_list <- list()
  v_hat_mat <- matrix(NA, nrow = n, ncol = length(endo_vars))
  colnames(v_hat_mat) <- paste0("vhat_", endo_vars)
  
  for (i in seq_along(endo_vars)) {
    en <- endo_vars[i]
    y_en <- sub_data[[en]]
    
    # First stage regression
    mod_unres <- stats::lm(y_en ~ Z[, -1, drop = FALSE])
    
    if (length(exog_vars) > 0) {
      exog_cols <- intersect(colnames(Z), exog_vars)
      if (length(exog_cols) > 0) {
        mod_res <- stats::lm(y_en ~ Z[, exog_cols, drop = FALSE])
      } else {
        mod_res <- stats::lm(y_en ~ 1)
      }
    } else {
      mod_res <- stats::lm(y_en ~ 1)
    }
    
    f_test <- stats::anova(mod_res, mod_unres)
    f_stat <- f_test$F[2]
    p_val <- f_test$`Pr(>F)`[2]
    
    first_stage_list[[en]] <- list(
      f_stat = f_stat,
      p_val = p_val,
      is_weak = (!is.na(f_stat) && f_stat < 10)
    )
    v_hat_mat[, i] <- stats::residuals(mod_unres)
  }
  
  # Wu-Hausman Test for Endogeneity
  mod_hausman_unres <- stats::lm(y ~ X[, -1, drop = FALSE] + v_hat_mat)
  mod_hausman_res <- stats::lm(y ~ X[, -1, drop = FALSE])
  h_test <- stats::anova(mod_hausman_res, mod_hausman_unres)
  hausman_stat <- h_test$F[2]
  hausman_p <- h_test$`Pr(>F)`[2]
  
  # Sargan Test for Overidentifying Restrictions (if l > k)
  overid_df <- length(iv_vars) - length(endo_vars)
  sargan_stat <- NA
  sargan_p <- NA
  if (overid_df > 0) {
    aux_sargan <- stats::lm(e ~ Z[, -1, drop = FALSE])
    r2_sargan <- summary(aux_sargan)$r.squared
    sargan_stat <- n * r2_sargan
    sargan_p <- stats::pchisq(sargan_stat, df = overid_df, lower.tail = FALSE)
  }
  
  res <- list(
    formula = formula,
    coefficients = b,
    coef_table = coef_table,
    vcov = vcov_mat,
    residuals = e,
    fitted.values = as.vector(X %*% b),
    n = n,
    k = k,
    df_residual = df_resid,
    df_inference = df_inference,
    r_squared = r2,
    adj_r_squared = adj_r2,
    se_type = se_label,
    digits = digits,
    endo_vars = endo_vars,
    iv_vars = iv_vars,
    exog_vars = exog_vars,
    first_stage = first_stage_list,
    y = y,
    X = X,
    Z = Z,
    hausman_stat = hausman_stat,
    hausman_p = hausman_p,
    sargan_stat = sargan_stat,
    sargan_p = sargan_p,
    overid_df = overid_df
  )
  class(res) <- "econ_iv"
  return(res)
}

#' @export
print.econ_iv <- function(x, digits = NULL, ...) {
  d <- if (is.null(digits)) x$digits else digits
  
  print_line("=")
  cat("             INSTRUMENTAL VARIABLES (2SLS) REGRESSION                   \n")
  cat(paste0("  Formula: ", deparse(x$formula), "\n"))
  cat(sprintf("  Sample Size: %-13d R-squared:    %-10.4f Adj. R-sq: %.4f\n", 
              x$n, x$r_squared, x$adj_r_squared))
  cat(paste0("  Endogenous: ", paste(x$endo_vars, collapse = ", "), 
             " | Instruments: ", paste(x$iv_vars, collapse = ", "), "\n"))
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
  cat("  Signif. codes:  *** p < 0.001,  ** p < 0.01,  * p < 0.05,  . p < 0.1\n\n")
  
  # Diagnostics Suite
  print_line("=")
  cat("                   2SLS DIAGNOSTIC TESTS & VALIDATION                   \n")
  print_line("=")
  cat(sprintf("  %-25s %-22s %9s %10s  %-18s\n", "Diagnostic Check", "Null Hypothesis (H0)", "Stat", "p-value", "Verdict"))
  print_line("-")
  
  # First Stage F
  for (en in names(x$first_stage)) {
    fs <- x$first_stage[[en]]
    fs_verdict <- if (!is.na(fs$f_stat) && fs$f_stat >= 10) "PASS (Strong IV, F > 10)" else "WEAK IV (F < 10)"
    cat(sprintf("  %-25s %-22s %9.2f %10s  %-18s\n",
                paste0("Weak IV: ", en),
                "Instruments are Weak",
                fs$f_stat,
                format_pval(fs$p_val),
                fs_verdict))
  }
  
  # Wu-Hausman
  h_verdict <- if (!is.na(x$hausman_p) && x$hausman_p < 0.05) "REJECT (Endogenous; 2SLS needed)" else "Fail to reject (OLS consistent)"
  cat(sprintf("  %-25s %-22s %9.2f %10s  %-18s\n",
              "Wu-Hausman Endogeneity",
              "Variables are Exogenous",
              x$hausman_stat,
              format_pval(x$hausman_p),
              h_verdict))
  
  # Sargan
  if (!is.na(x$sargan_stat)) {
    s_verdict <- if (x$sargan_p < 0.05) "REJECT (Invalid instruments)" else "Fail to reject (Valid IVs)"
    cat(sprintf("  %-25s %-22s %9.2f %10s  %-18s\n",
                paste0("Sargan Overid (df=", x$overid_df, ")"),
                "Instruments are Exogenous",
                x$sargan_stat,
                format_pval(x$sargan_p),
                s_verdict))
  } else {
    cat(sprintf("  %-25s %-22s %9s %10s  %-18s\n",
                "Sargan Overid",
                "Instruments are Exogenous",
                "-", "-", "Exactly Identified"))
  }
  print_line("=")
  invisible(x)
}

#' @export
summary.econ_iv <- function(object, ...) {
  object
}

#' @export
coef.econ_iv <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.econ_iv <- function(object, ...) {
  object$vcov
}

#' @export
residuals.econ_iv <- function(object, ...) {
  object$residuals
}

#' @export
fitted.econ_iv <- function(object, ...) {
  object$fitted.values
}

#' @export
nobs.econ_iv <- function(object, ...) {
  object$n
}

#' @export
confint.econ_iv <- function(object, parm, level = 0.95, ...) {
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
predict.econ_iv <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    return(object$fitted.values)
  }
  form_no_y <- stats::delete.response(stats::terms(object$formula))
  X_new <- stats::model.matrix(form_no_y, data = newdata)
  b <- stats::coef(object)
  common <- intersect(colnames(X_new), names(b))
  as.vector(X_new[, common, drop = FALSE] %*% b[common])
}
