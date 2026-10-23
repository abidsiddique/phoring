#' Comprehensive Panel Data Econometrics with Two-Way FE and Specification Tests
#'
#' Fits panel data models (Fixed Effects, Two-Way Fixed Effects, Random Effects,
#' First Differences, Between Estimator, or Pooled OLS) with customizable standard errors
#' (clustered by entity id, White robust, or classical standard) and automatically
#' executes the Hausman test, Breusch-Pagan LM test, and Wooldridge panel autocorrelation test.
#'
#' @param formula Model formula (e.g., \code{inv ~ value + capital} or \code{log(inv) ~ value + I(capital^2)}).
#' @param data A data.frame containing the panel data.
#' @param id Column name (or unquoted symbol) indicating the entity/group identifier.
#' @param time Column name (or unquoted symbol) indicating the time period (optional, required for TWFE and FD).
#' @param model Model type: \code{"fe"} (One-way Fixed Effects / Within, default),
#'   \code{"twfe"} (Two-Way Fixed Effects: entity + time),
#'   \code{"re"} (Random Effects / Swamy-Arora GLS),
#'   \code{"fd"} (First Differences),
#'   \code{"between"} (Between / Group Means Estimator),
#'   or \code{"pooling"} (Pooled OLS).
#' @param se Type of standard errors: \code{"cluster"} (clustered by entity id, default),
#'   \code{"robust"} (White HC1), or \code{"standard"} (classical homoskedastic).
#' @param digits Number of decimal places to display (default is 4).
#' @param ... Additional arguments.
#'
#' @importFrom stats ave df.residual model.frame model.matrix model.response na.omit
#' @return An object of class \code{econ_panel} containing estimates, standard errors,
#'   and panel diagnostics.
#' @export
#'
#' @examples
#' data(firm_investment, package = "phoring")
#' # 1. Entity Fixed Effects (Within)
#' panel(inv ~ value + capital, data = firm_investment, id = "firm_id", time = "year", model = "fe")
#'
#' # 2. Two-Way Fixed Effects (Entity + Time FE)
#' panel(inv ~ value + capital, data = firm_investment, id = "firm_id", time = "year", model = "twfe")
#'
#' # 3. Random Effects (GLS)
#' panel(inv ~ value + capital, data = firm_investment, id = "firm_id", time = "year", model = "re")
panel <- function(formula, data, id, time = NULL,
                  model = c("fe", "twfe", "re", "fd", "between", "pooling"),
                  se = c("cluster", "robust", "standard"),
                  digits = 4, ...) {
  model <- match.arg(model)
  se <- match.arg(se)
  
  if (missing(data) || is.null(data)) {
    stop("Please provide a dataset in `data`.", call. = FALSE)
  }
  
  # Resolve id column
  id_col <- deparse(substitute(id))
  id_col <- gsub("\"", "", id_col)
  if (!id_col %in% names(data)) {
    stop(paste0("Entity ID column '", id_col, "' not found in data."), call. = FALSE)
  }
  
  # Resolve time column
  time_col <- NULL
  if (!missing(time) && !is.null(time)) {
    time_col <- deparse(substitute(time))
    time_col <- gsub("\"", "", time_col)
    if (!time_col %in% names(data)) {
      stop(paste0("Time column '", time_col, "' not found in data."), call. = FALSE)
    }
  }
  
  if (model %in% c("twfe", "fd") && is.null(time_col)) {
    stop(paste0("Model '", model, "' requires the `time` variable to be specified."), call. = FALSE)
  }
  
  # Extract variables and drop missing
  vars_in_form <- all.vars(formula)
  all_cols <- unique(c(vars_in_form, id_col, time_col))
  missing_cols <- setdiff(all_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(paste0("Variables not found in data: ", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  
  sub_data <- stats::na.omit(data[, all_cols, drop = FALSE])
  
  # Sort data by ID and time if time is provided
  if (!is.null(time_col)) {
    sub_data <- sub_data[order(sub_data[[id_col]], sub_data[[time_col]]), ]
    time_vec <- as.factor(sub_data[[time_col]])
    T_periods <- nlevels(time_vec)
  } else {
    sub_data <- sub_data[order(sub_data[[id_col]]), ]
    time_vec <- NULL
    T_periods <- NA
  }
  
  id_vec <- as.factor(sub_data[[id_col]])
  G <- nlevels(id_vec)
  n <- nrow(sub_data)
  T_bar <- n / G
  
  # Build full model frame and design matrix preserving all functional forms (log, I, factor)
  mf <- stats::model.frame(formula, data = sub_data)
  y_raw <- stats::model.response(mf)
  X_raw <- stats::model.matrix(formula, data = mf)
  X_no_int <- X_raw[, colnames(X_raw) != "(Intercept)", drop = FALSE]
  k_vars <- ncol(X_no_int)
  
  calc_means <- function(vec, grp) {
    stats::ave(vec, grp, FUN = mean)
  }
  
  # Helper to compute panel standard errors based on chosen 'se'
  calc_panel_se <- function(X_mat, e_vec, bread_mat, df_resid, se_type, grp_vec) {
    k_m <- ncol(X_mat)
    if (se_type == "cluster") {
      clust_scores <- rowsum(X_mat * e_vec, grp_vec)
      meat <- crossprod(clust_scores)
      G_c <- nlevels(grp_vec)
      dfc <- (G_c / (G_c - 1)) * ((nrow(X_mat) - 1) / df_resid)
      vcov <- dfc * (bread_mat %*% meat %*% bread_mat)
      df_inf <- max(1, G_c - 1)
      label <- paste0("Clustered (by ", id_col, ", G = ", G_c, ")")
    } else if (se_type == "robust") {
      dfc <- nrow(X_mat) / df_resid
      meat <- crossprod(X_mat * e_vec)
      vcov <- dfc * (bread_mat %*% meat %*% bread_mat)
      df_inf <- df_resid
      label <- "Robust (White HC1)"
    } else {
      s2 <- sum(e_vec^2) / df_resid
      vcov <- s2 * bread_mat
      df_inf <- df_resid
      label <- "Standard (Homoskedastic)"
    }
    dimnames(vcov) <- list(colnames(X_mat), colnames(X_mat))
    list(vcov = vcov, df_inference = df_inf, label = label)
  }
  
  # -------------------------------------------------------------
  # 1. FIXED EFFECTS (WITHIN ESTIMATOR)
  # -------------------------------------------------------------
  y_mean_id <- calc_means(y_raw, id_vec)
  y_fe <- y_raw - y_mean_id
  
  X_fe <- matrix(NA, nrow = n, ncol = k_vars)
  colnames(X_fe) <- colnames(X_no_int)
  for (j in 1:k_vars) {
    x_col <- X_no_int[, j]
    X_fe[, j] <- x_col - calc_means(x_col, id_vec)
  }
  
  col_vars <- apply(X_fe, 2, stats::var)
  valid_cols_fe <- col_vars >= 1e-12
  if (!any(valid_cols_fe)) {
    stop("All regressors are time-invariant; Fixed Effects estimator cannot be identified.", call. = FALSE)
  }
  X_fe <- X_fe[, valid_cols_fe, drop = FALSE]
  k_fe <- ncol(X_fe)
  
  bread_fe <- safe_solve(crossprod(X_fe))
  b_fe <- as.vector(bread_fe %*% crossprod(X_fe, y_fe))
  names(b_fe) <- colnames(X_fe)
  
  e_fe <- as.vector(y_fe - X_fe %*% b_fe)
  df_fe <- n - G - k_fe
  s2_fe <- sum(e_fe^2) / df_fe
  
  se_fe_res <- calc_panel_se(X_fe, e_fe, bread_fe, df_fe, se, id_vec)
  vcov_fe <- se_fe_res$vcov
  
  # -------------------------------------------------------------
  # 2. TWO-WAY FIXED EFFECTS (TWFE: Entity + Time)
  # -------------------------------------------------------------
  b_twfe <- NULL
  vcov_twfe <- NULL
  e_twfe <- NULL
  df_twfe <- NULL
  
  if (!is.null(time_col)) {
    # Exact Frisch-Waugh-Lovell demeaning for balanced AND unbalanced panels
    time_dummies <- stats::model.matrix(~ time_vec)[, -1, drop = FALSE]
    time_fe <- matrix(NA, nrow = n, ncol = ncol(time_dummies))
    colnames(time_fe) <- colnames(time_dummies)
    for (m in 1:ncol(time_dummies)) {
      t_col <- time_dummies[, m]
      time_fe[, m] <- t_col - calc_means(t_col, id_vec)
    }
    
    # Combined design matrix of demeaned X and demeaned time dummies
    X_all_twfe <- cbind(X_fe, time_fe)
    # Check for collinear time dummies in unbalanced structures
    tw_rank <- qr(X_all_twfe)
    keep_cols <- tw_rank$pivot[seq_len(tw_rank$rank)]
    X_all_twfe <- X_all_twfe[, keep_cols, drop = FALSE]
    
    bread_all_twfe <- safe_solve(crossprod(X_all_twfe))
    b_all_twfe <- as.vector(bread_all_twfe %*% crossprod(X_all_twfe, y_fe))
    names(b_all_twfe) <- colnames(X_all_twfe)
    
    e_twfe <- as.vector(y_fe - X_all_twfe %*% b_all_twfe)
    df_twfe <- n - G - ncol(X_all_twfe)
    
    se_twfe_res <- calc_panel_se(X_all_twfe, e_twfe, bread_all_twfe, df_twfe, se, id_vec)
    vcov_all_twfe <- se_twfe_res$vcov
    
    # Extract structural coefficients on X
    x_names_twfe <- intersect(colnames(X_fe), colnames(X_all_twfe))
    b_twfe <- b_all_twfe[x_names_twfe]
    vcov_twfe <- vcov_all_twfe[x_names_twfe, x_names_twfe, drop = FALSE]
  }
  
  # -------------------------------------------------------------
  # 3. RANDOM EFFECTS (SWAMY-ARORA GLS)
  # -------------------------------------------------------------
  y_grp <- tapply(y_raw, id_vec, mean)
  X_grp <- apply(X_no_int, 2, function(col) tapply(col, id_vec, mean))
  mod_between <- stats::lm(y_grp ~ X_grp)
  s2_between <- sum(stats::residuals(mod_between)^2) / max(1, (G - k_vars - 1))
  
  s2_u <- max(0, (s2_between - s2_fe / T_bar))
  theta <- 1 - sqrt(s2_fe / (s2_fe + T_bar * s2_u))
  
  y_re <- y_raw - theta * y_mean_id
  X_re <- X_raw
  for (j in 1:ncol(X_raw)) {
    x_col <- X_raw[, j]
    X_re[, j] <- x_col - theta * calc_means(x_col, id_vec)
  }
  
  bread_re <- safe_solve(crossprod(X_re))
  b_re <- as.vector(bread_re %*% crossprod(X_re, y_re))
  names(b_re) <- colnames(X_re)
  e_re <- as.vector(y_re - X_re %*% b_re)
  df_re <- n - ncol(X_re)
  
  se_re_res <- calc_panel_se(X_re, e_re, bread_re, df_re, se, id_vec)
  vcov_re <- se_re_res$vcov
  
  # -------------------------------------------------------------
  # 4. FIRST DIFFERENCES (FD)
  # -------------------------------------------------------------
  b_fd <- NULL
  vcov_fd <- NULL
  e_fd <- NULL
  df_fd_resid <- NULL
  g_fd <- NULL
  
  if (!is.null(time_col)) {
    diff_data_list <- list()
    for (g_id in levels(id_vec)) {
      sub_g <- sub_data[sub_data[[id_col]] == g_id, , drop = FALSE]
      if (nrow(sub_g) > 1) {
        mf_g <- stats::model.frame(formula, data = sub_g)
        y_g <- stats::model.response(mf_g)
        X_g <- stats::model.matrix(formula, data = mf_g)
        X_g_noint <- X_g[, colnames(X_g) != "(Intercept)", drop = FALSE]
        
        dy <- diff(y_g)
        dX <- apply(X_g_noint, 2, diff)
        diff_data_list[[g_id]] <- list(dy = dy, dX = dX, g_id = rep(g_id, length(dy)))
      }
    }
    
    if (length(diff_data_list) > 0) {
      y_fd <- unlist(lapply(diff_data_list, function(d) d$dy))
      X_fd <- do.call(rbind, lapply(diff_data_list, function(d) d$dX))
      g_fd <- as.factor(unlist(lapply(diff_data_list, function(d) d$g_id)))
      
      # Drop time-invariant columns
      fd_vars <- apply(X_fd, 2, stats::var) >= 1e-12
      X_fd <- X_fd[, fd_vars, drop = FALSE]
      n_fd <- length(y_fd)
      
      bread_fd <- safe_solve(crossprod(X_fd))
      b_fd <- as.vector(bread_fd %*% crossprod(X_fd, y_fd))
      names(b_fd) <- colnames(X_fd)
      e_fd <- as.vector(y_fd - X_fd %*% b_fd)
      df_fd_resid <- n_fd - ncol(X_fd)
      
      se_fd_res <- calc_panel_se(X_fd, e_fd, bread_fd, df_fd_resid, se, g_fd)
      vcov_fd <- se_fd_res$vcov
    }
  }
  
  # -------------------------------------------------------------
  # 5. POOLED OLS
  # -------------------------------------------------------------
  base_pool <- stats::lm(formula, data = sub_data)
  b_pool <- stats::coef(base_pool)
  e_pool <- stats::residuals(base_pool)
  X_pool <- stats::model.matrix(base_pool)
  bread_pool <- safe_solve(crossprod(X_pool))
  df_pool <- n - length(b_pool)
  
  se_pool_res <- calc_panel_se(X_pool, e_pool, bread_pool, df_pool, se, id_vec)
  vcov_pool <- se_pool_res$vcov
  
  # -------------------------------------------------------------
  # SPECIFICATION TESTS
  # -------------------------------------------------------------
  # A. Hausman Test (Mundlak Approach)
  X_means <- matrix(NA, nrow = n, ncol = k_fe)
  colnames(X_means) <- paste0("mean_", colnames(X_fe))
  for (j in 1:k_fe) {
    X_means[, j] <- calc_means(X_no_int[, colnames(X_fe)[j]], id_vec)
  }
  mod_mund_unres <- stats::lm(y_raw ~ X_raw[, -1, drop = FALSE] + X_means)
  mod_mund_res <- stats::lm(y_raw ~ X_raw[, -1, drop = FALSE])
  h_anova <- stats::anova(mod_mund_res, mod_mund_unres)
  hausman_stat <- h_anova$F[2] * k_fe
  hausman_df <- k_fe
  hausman_p <- stats::pchisq(hausman_stat, df = hausman_df, lower.tail = FALSE)
  
  # B. Breusch-Pagan LM Test for Panel (Pooled OLS vs RE)
  e_i_sum <- tapply(e_pool, id_vec, sum)
  lm_bp_panel <- (n * T_bar / (2 * max(1, T_bar - 1))) * ((sum(e_i_sum^2) / sum(e_pool^2) - 1)^2)
  lm_bp_panel <- max(0, lm_bp_panel)
  lm_bp_p <- stats::pchisq(lm_bp_panel, df = 1, lower.tail = FALSE)
  
  # C. Wooldridge Serial Correlation Test for Panel
  wooldridge_ar_stat <- NA
  wooldridge_ar_p <- NA
  if (!is.null(time_col) && !is.null(e_fd)) {
    lag_efd <- unlist(tapply(e_fd, g_fd, function(v) c(NA, v[-length(v)])))
    complete_ar <- !is.na(lag_efd)
    if (sum(complete_ar) > 10) {
      mod_ar <- stats::lm(e_fd[complete_ar] ~ lag_efd[complete_ar])
      rho_hat <- stats::coef(mod_ar)[2]
      se_rho <- summary(mod_ar)$coefficients[2, 2]
      t_ar <- (rho_hat - (-0.5)) / se_rho
      wooldridge_ar_stat <- t_ar^2
      wooldridge_ar_p <- 2 * stats::pt(-abs(t_ar), df = df.residual(mod_ar))
    }
  }
  
  # Select results based on requested model
  if (model == "fe") {
    b_out <- b_fe
    vcov_out <- vcov_fe
    e_out <- e_fe
    fit_out <- y_raw - e_fe
    df_out <- se_fe_res$df_inference
    df_resid_out <- df_fe
    r2_within <- 1 - (sum(e_fe^2) / sum((y_fe)^2))
    se_label <- se_fe_res$label
  } else if (model == "twfe") {
    if (is.null(b_twfe)) stop("TWFE requires valid time variation.", call. = FALSE)
    b_out <- b_twfe
    vcov_out <- vcov_twfe
    e_out <- e_twfe
    fit_out <- y_raw - e_twfe
    df_out <- se_twfe_res$df_inference
    df_resid_out <- df_twfe
    r2_within <- 1 - (sum(e_twfe^2) / sum((y_fe)^2))
    se_label <- se_twfe_res$label
  } else if (model == "re") {
    b_out <- b_re
    vcov_out <- vcov_re
    e_out <- e_re
    fit_out <- y_raw - e_re
    df_out <- se_re_res$df_inference
    df_resid_out <- df_re
    r2_within <- 1 - (sum(e_re^2) / sum((y_re)^2))
    se_label <- se_re_res$label
  } else if (model == "fd") {
    if (is.null(b_fd)) stop("FD requires at least 2 time periods per entity.", call. = FALSE)
    b_out <- b_fd
    vcov_out <- vcov_fd
    e_out <- e_fd
    fit_out <- y_fd - e_fd
    df_out <- se_fd_res$df_inference
    df_resid_out <- df_fd_resid
    r2_within <- 1 - (sum(e_fd^2) / sum((y_fd - mean(y_fd))^2))
    se_label <- se_fd_res$label
  } else if (model == "between") {
    b_out <- stats::coef(mod_between)
    vcov_out <- stats::vcov(mod_between)
    e_out <- stats::residuals(mod_between)
    fit_out <- as.vector(stats::fitted(mod_between))
    df_resid_out <- G - length(b_out)
    df_out <- df_resid_out
    r2_within <- summary(mod_between)$r.squared
    se_label <- "Standard (Homoskedastic)"
  } else {
    b_out <- b_pool
    vcov_out <- vcov_pool
    e_out <- e_pool
    fit_out <- as.vector(stats::fitted(base_pool))
    df_out <- se_pool_res$df_inference
    df_resid_out <- df_pool
    r2_within <- summary(base_pool)$r.squared
    se_label <- se_pool_res$label
  }
  
  se_b <- sqrt(pmax(0, diag(vcov_out)))
  t_b <- b_out / se_b
  p_b <- 2 * stats::pt(-abs(t_b), df = df_out)
  
  coef_table <- data.frame(
    Estimate  = b_out,
    Std_Error = se_b,
    t_stat    = t_b,
    p_value   = p_b,
    Signif    = sapply(p_b, get_stars),
    stringsAsFactors = FALSE
  )
  
  res <- list(
    formula = formula,
    model = model,
    coefficients = b_out,
    coef_table = coef_table,
    vcov = vcov_out,
    residuals = e_out,
    fitted.values = fit_out,
    y_raw = y_raw,
    X_raw = X_raw,
    id_vec = id_vec,
    time_vec = time_vec,
    n = n,
    G = G,
    T_bar = T_bar,
    df_residual = df_resid_out,
    df_inference = df_out,
    id_col = id_col,
    time_col = time_col,
    r2_within = r2_within,
    se_type = se_label,
    digits = digits,
    hausman_stat = hausman_stat,
    hausman_df = hausman_df,
    hausman_p = hausman_p,
    lm_bp_panel = lm_bp_panel,
    lm_bp_p = lm_bp_p,
    wooldridge_ar_stat = wooldridge_ar_stat,
    wooldridge_ar_p = wooldridge_ar_p,
    theta = theta
  )
  class(res) <- "econ_panel"
  return(res)
}

#' @export
print.econ_panel <- function(x, digits = NULL, ...) {
  d <- if (is.null(digits)) x$digits else digits
  
  model_title <- switch(x$model,
                        fe = "PANEL DATA REGRESSION: ONE-WAY FIXED EFFECTS (WITHIN)",
                        twfe = "PANEL DATA REGRESSION: TWO-WAY FIXED EFFECTS (ENTITY + TIME)",
                        re = "PANEL DATA REGRESSION: RANDOM EFFECTS (SWAMY-ARORA GLS)",
                        fd = "PANEL DATA REGRESSION: FIRST DIFFERENCES (FD)",
                        between = "PANEL DATA REGRESSION: BETWEEN ESTIMATOR (GROUP MEANS)",
                        pooling = "PANEL DATA REGRESSION: POOLED OLS")
  
  print_line("=")
  cat(sprintf("               %s\n", model_title))
  cat(paste0("  Formula: ", deparse(x$formula), "\n"))
  cat(sprintf("  Entities (G): %-12d Avg Time (T): %-11.2f Total Obs (N): %d\n",
              x$G, x$T_bar, x$n))
  cat(sprintf("  Within R-sq:  %-12.4f Standard Errors: %s\n",
              x$r2_within, x$se_type))
  if (x$model == "re") {
    cat(sprintf("  RE Quasi-Demeaning Weight (theta): %.4f\n", x$theta))
  }
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
  
  # Specification Tests Card
  cat("\n")
  print_line("=")
  cat("                   PANEL SPECIFICATION TESTS                            \n")
  print_line("=")
  cat(sprintf("  %-25s %-25s %8s %10s  %-18s\n", "Test Name", "Null Hypothesis (H0)", "Stat", "p-value", "Result"))
  print_line("-")
  
  h_verdict <- if (x$hausman_p < 0.05) "REJECT (Use FE)" else "Fail to reject (RE OK)"
  cat(sprintf("  %-25s %-25s %8.3f %10s  %-18s\n",
              "Hausman Test (Specification)", "RE is consistent (cov=0)", x$hausman_stat, format_pval(x$hausman_p), h_verdict))
  
  bp_verdict <- if (x$lm_bp_p < 0.05) "REJECT (Use RE/FE)" else "Fail to reject (Pooled OK)"
  cat(sprintf("  %-25s %-25s %8.3f %10s  %-18s\n",
              "Breusch-Pagan LM", "No entity effects", x$lm_bp_panel, format_pval(x$lm_bp_p), bp_verdict))
  
  if (!is.na(x$wooldridge_ar_stat)) {
    ar_verdict <- if (x$wooldridge_ar_p < 0.05) "REJECT (Serial AR1)" else "Fail to reject (No AR1)"
    cat(sprintf("  %-25s %-25s %8.3f %10s  %-18s\n",
                "Wooldridge Panel AR(1)", "No first-order serial corr", x$wooldridge_ar_stat, format_pval(x$wooldridge_ar_p), ar_verdict))
  }
  print_line("-")
  cat("  Econometric Guidance Notes:\n")
  if (x$hausman_p < 0.05) {
    cat("  * Hausman p < 0.05: Entity unobservables correlate with regressors -> Use Fixed Effects.\n")
  } else {
    cat("  * Hausman p >= 0.05: Random Effects is consistent and more efficient than FE.\n")
  }
  if (!is.na(x$wooldridge_ar_p) && x$wooldridge_ar_p < 0.05) {
    cat("  * Serial correlation detected (Wooldridge p < 0.05) -> Clustered SEs are strongly recommended.\n")
  }
  print_line("=")
  invisible(x)
}

#' @export
summary.econ_panel <- function(object, ...) {
  object
}

#' @export
coef.econ_panel <- function(object, ...) {
  object$coefficients
}

#' @export
vcov.econ_panel <- function(object, ...) {
  object$vcov
}

#' @export
residuals.econ_panel <- function(object, ...) {
  object$residuals
}

#' @export
fitted.econ_panel <- function(object, ...) {
  object$fitted.values
}

#' @export
nobs.econ_panel <- function(object, ...) {
  object$n
}

#' @export
confint.econ_panel <- function(object, parm, level = 0.95, ...) {
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
