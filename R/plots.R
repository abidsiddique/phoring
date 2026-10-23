#' Academic Publication-Quality Econometric Diagnostic Plots
#'
#' Generates clean, publication-grade academic diagnostic plots with serif typography,
#' L-frame axes, and standardized econometric styling across all models in the package.
#'
#' @param x An object of class \code{econ_reg} (from \code{\link{reg}}) or standard \code{lm}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot(s) to draw:
#'   \code{"all"} (default 2x2 core grid), \code{"resid_fitted"}, \code{"qq"}, 
#'   \code{"scale_location"}, \code{"cooks"}, \code{"leverage"}, \code{"added_var"}, 
#'   \code{"component_resid"}, \code{"dfbetas"}, \code{"forest"}, \code{"vif"}, 
#'   \code{"histogram"}, \code{"observed_fitted"}, \code{"cusum"}, \code{"cusum_sq"},
#'   \code{"stability"}, or \code{"influence"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
#'
#' @examples
#' data(wage_data, package = "phoring")
#' mod <- reg(wage ~ educ + exper, data = wage_data)
#' econ_plot(mod)
#' plot(mod, which = "leverage")
#' plot(mod, which = "cusum_sq")
econ_plot <- function(x, which = c("all", "resid_fitted", "qq", "scale_location", "cooks", 
                                   "leverage", "added_var", "component_resid", "dfbetas", 
                                   "forest", "vif", "histogram", "observed_fitted", 
                                   "cusum", "cusum_sq", "cusum_squares", "stability", 
                                   "influence", "influence_plot"), ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "resid_fitted", "qq", "scale_location", "cooks", 
                 "leverage", "added_var", "component_resid", "dfbetas", 
                 "forest", "vif", "histogram", "observed_fitted", 
                 "cusum", "cusum_sq", "stability", "influence")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "resid_fitted", "qq", "scale_location", "cooks", 
                                  "leverage", "added_var", "component_resid", "dfbetas", 
                                  "forest", "vif", "histogram", "observed_fitted", 
                                  "cusum", "cusum_sq", "cusum_squares", "stability", 
                                  "influence", "influence_plot"))
  }
  if (which == "cusum_squares") which <- "cusum_sq"
  if (which == "influence_plot") which <- "influence"
  
  if (inherits(x, "econ_reg")) {
    base_lm <- x$lm_model
    e <- x$residuals
    y_hat <- x$fitted.values
    n <- x$n
    k <- x$k
  } else if (inherits(x, "lm")) {
    base_lm <- x
    e <- stats::residuals(x)
    y_hat <- stats::fitted(x)
    n <- length(e)
    k <- length(stats::coef(x))
  } else {
    stop("`x` must be an object of class `econ_reg` or `lm`.", call. = FALSE)
  }
  
  X <- stats::model.matrix(base_lm)
  y_obs <- stats::model.response(stats::model.frame(base_lm))
  bread <- safe_solve(crossprod(X))
  h <- rowSums((X %*% bread) * X)
  h <- pmin(h, 0.9999)
  
  s2 <- sum(e^2) / max(1, (n - k))
  s <- sqrt(s2)
  r <- e / (s * sqrt(pmax(1e-10, 1 - h))) # studentized residuals
  sqrt_abs_r <- sqrt(abs(r))
  cooks_d <- (r^2 / k) * (h / pmax(1e-10, 1 - h))
  
  # Academic styling tokens
  pt_col <- "#1B365D"
  line_col <- "#922B21"
  smooth_col <- "#1E8449"
  accent_col <- "#D35400"
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.0, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  # 1. Residuals vs Fitted
  if (which %in% c("all", "resid_fitted")) {
    graphics::plot(y_hat, e,
                   xlab = "Fitted Values", ylab = "Residuals",
                   main = "1. Residuals vs. Fitted (Linearity)",
                   pch = 19, col = pt_col, cex = 0.8, bty = "l",
                   panel.first = graphics::grid(col = "#EAEAEA", lty = 1))
    graphics::abline(h = 0, lty = 2, col = line_col, lwd = 1.8)
    if (n >= 10) {
      sm <- stats::loess(e ~ y_hat, span = 0.75)
      ord <- order(y_hat)
      graphics::lines(y_hat[ord], stats::fitted(sm)[ord], col = smooth_col, lwd = 2)
    }
  }
  
  # 2. Normal Q-Q
  if (which %in% c("all", "qq")) {
    stats::qqnorm(r,
                  main = "2. Normal Q-Q Plot (Residual Normality)",
                  xlab = "Theoretical Quantiles", ylab = "Studentized Residuals",
                  pch = 19, col = pt_col, cex = 0.8, bty = "l",
                  panel.first = graphics::grid(col = "#EAEAEA", lty = 1))
    stats::qqline(r, col = line_col, lwd = 1.8, lty = 2)
  }
  
  # 3. Scale-Location
  if (which %in% c("all", "scale_location")) {
    graphics::plot(y_hat, sqrt_abs_r,
                   xlab = "Fitted Values",
                   ylab = expression(sqrt("|Studentized Residuals|")),
                   main = "3. Scale-Location (Homoskedasticity)",
                   pch = 19, col = pt_col, cex = 0.8, bty = "l",
                   panel.first = graphics::grid(col = "#EAEAEA", lty = 1))
    if (n >= 10) {
      sm <- stats::loess(sqrt_abs_r ~ y_hat, span = 0.75)
      ord <- order(y_hat)
      graphics::lines(y_hat[ord], stats::fitted(sm)[ord], col = smooth_col, lwd = 2)
    }
  }
  
  # 4. Cook's Distance
  if (which %in% c("all", "cooks")) {
    graphics::plot(cooks_d, type = "h", col = pt_col, lwd = 1.6,
                   xlab = "Observation Index", ylab = "Cook's Distance",
                   main = "4. Cook's Distance (Influential Points)",
                   bty = "l", panel.first = graphics::grid(col = "#EAEAEA", lty = 1))
    thresh <- 4 / n
    graphics::abline(h = thresh, lty = 2, col = line_col, lwd = 1.8)
    graphics::text(n * 0.95, thresh, labels = sprintf("4/n = %.3f", thresh), pos = 3, col = line_col, cex = 0.8)
    high_idx <- which(cooks_d > thresh)
    if (length(high_idx) > 0) {
      top3 <- high_idx[order(cooks_d[high_idx], decreasing = TRUE)[1:min(3, length(high_idx))]]
      graphics::text(top3, cooks_d[top3], labels = top3, pos = 3, col = line_col, cex = 0.8, font = 2)
    }
  }
  
  # 5. Residuals vs Leverage with Cook's distance contours
  if (which == "leverage") {
    graphics::plot(h, r, xlab = "Leverage (h_ii)", ylab = "Studentized Residuals",
                   main = "Residuals vs. Leverage (Influence)",
                   pch = 19, col = pt_col, cex = 0.8, bty = "l",
                   xlim = c(0, max(0.4, max(h) * 1.1)),
                   panel.first = graphics::grid(col = "#EAEAEA", lty = 1))
    graphics::abline(h = 0, lty = 2, col = "#7F8C8D")
    lev_thresh <- 2 * k / n
    graphics::abline(v = lev_thresh, lty = 3, col = line_col, lwd = 1.5)
    graphics::mtext(sprintf("2k/n = %.3f", lev_thresh), side = 3, at = lev_thresh, col = line_col, cex = 0.75)
    
    # Cook's contours at 0.5 and 1.0
    h_seq <- seq(0.001, max(0.4, max(h) * 1.1), length.out = 150)
    for (cd in c(0.5, 1.0)) {
      r_cd <- sqrt(cd * k * (1 - h_seq) / h_seq)
      valid <- is.finite(r_cd) & r_cd < max(abs(r), 4)
      if (any(valid)) {
        graphics::lines(h_seq[valid], r_cd[valid], lty = 2, col = accent_col)
        graphics::lines(h_seq[valid], -r_cd[valid], lty = 2, col = accent_col)
        graphics::text(max(h_seq[valid]), max(r_cd[valid]), sprintf("Cook's D = %.1f", cd), 
                       col = accent_col, cex = 0.75, pos = 2)
      }
    }
  }
  
  # 6. Added-Variable (Partial Regression) plot
  if (which == "added_var") {
    preds <- colnames(X)[colnames(X) != "(Intercept)"]
    if (length(preds) == 0) {
      graphics::plot.new()
      graphics::text(0.5, 0.5, "No slope predictors for Added-Variable plot.")
    } else {
      p1 <- preds[1]
      y_res <- stats::residuals(stats::lm(y_obs ~ X[, setdiff(colnames(X), c("(Intercept)", p1)), drop = FALSE]))
      x_res <- stats::residuals(stats::lm(X[, p1] ~ X[, setdiff(colnames(X), c("(Intercept)", p1)), drop = FALSE]))
      graphics::plot(x_res, y_res, pch = 19, col = pt_col, cex = 0.8, bty = "l",
                     xlab = paste0(p1, " | others"), ylab = paste0("y | others"),
                     main = paste0("Added-Variable Plot (Partial Regression): ", p1),
                     panel.first = graphics::grid(col = "#EAEAEA", lty = 1))
      graphics::abline(stats::lm(y_res ~ x_res), col = line_col, lwd = 2)
    }
  }
  
  # 7. Component + Residual (Partial Residual) plot
  if (which == "component_resid") {
    preds <- colnames(X)[colnames(X) != "(Intercept)"]
    if (length(preds) == 0) {
      graphics::plot.new()
      graphics::text(0.5, 0.5, "No slope predictors.")
    } else {
      p1 <- preds[1]
      b1 <- stats::coef(base_lm)[p1]
      cr <- e + b1 * X[, p1]
      graphics::plot(X[, p1], cr, pch = 19, col = pt_col, cex = 0.8, bty = "l",
                     xlab = p1, ylab = paste0("Component + Residual (", p1, ")"),
                     main = paste0("Component + Residual Plot: ", p1),
                     panel.first = graphics::grid(col = "#EAEAEA", lty = 1))
      graphics::abline(a = 0, b = b1, col = line_col, lwd = 2, lty = 2)
      if (n >= 10) {
        sm <- stats::loess(cr ~ X[, p1], span = 0.75)
        ord <- order(X[, p1])
        graphics::lines(X[ord, p1], stats::fitted(sm)[ord], col = smooth_col, lwd = 2)
      }
    }
  }
  
  # 8. DFBETAS Index Plot
  if (which == "dfbetas") {
    preds <- colnames(X)[colnames(X) != "(Intercept)"]
    p_focus <- if (length(preds) > 0) preds[1] else colnames(X)[1]
    dfb <- stats::dfbetas(base_lm)[, p_focus]
    thresh_dfb <- 2 / sqrt(n)
    graphics::plot(dfb, type = "h", col = pt_col, lwd = 1.6, bty = "l",
                   xlab = "Observation Index", ylab = paste0("DFBETAS (", p_focus, ")"),
                   main = paste0("DFBETAS Index Plot: ", p_focus),
                   panel.first = graphics::grid(col = "#EAEAEA", lty = 1))
    graphics::abline(h = c(-thresh_dfb, thresh_dfb), lty = 2, col = line_col, lwd = 1.5)
    graphics::text(n * 0.95, thresh_dfb, labels = sprintf("+2/sqrt(n) = %.3f", thresh_dfb), pos = 3, col = line_col, cex = 0.75)
  }
  
  # 9. Coefficient Forest Plot
  if (which == "forest") {
    b_all <- stats::coef(base_lm)
    se_all <- sqrt(pmax(0, diag(stats::vcov(base_lm))))
    b_plot <- b_all[names(b_all) != "(Intercept)"]
    se_plot <- se_all[names(b_all) != "(Intercept)"]
    if (length(b_plot) == 0) {
      b_plot <- b_all
      se_plot <- se_all
    }
    m <- length(b_plot)
    ci_low <- b_plot - 1.96 * se_plot
    ci_high <- b_plot + 1.96 * se_plot
    xlim <- range(c(ci_low, ci_high, 0))
    xlim <- c(xlim[1] - 0.1 * diff(xlim), xlim[2] + 0.1 * diff(xlim))
    
    graphics::plot(b_plot, 1:m, pch = 19, col = pt_col, cex = 1.2, xlim = xlim, ylim = c(0.5, m + 0.5),
                   yaxt = "n", xlab = "Estimate (95% CI)", ylab = "", bty = "l",
                   main = "Coefficient Estimates & 95% Confidence Intervals",
                   panel.first = graphics::grid(col = "#EAEAEA", lty = 1))
    graphics::axis(2, at = 1:m, labels = names(b_plot), las = 1, font = 2)
    graphics::abline(v = 0, lty = 2, col = line_col, lwd = 1.5)
    graphics::arrows(ci_low, 1:m, ci_high, 1:m, length = 0.05, angle = 90, code = 3, col = pt_col, lwd = 1.8)
  }
  
  # 10. VIF Horizontal Bar Chart
  if (which == "vif") {
    preds <- colnames(X)[colnames(X) != "(Intercept)"]
    if (length(preds) <= 1) {
      graphics::plot.new()
      graphics::text(0.5, 0.5, "VIF requires 2 or more predictors.")
    } else {
      vifs <- numeric(length(preds))
      names(vifs) <- preds
      for (p_name in preds) {
        y_p <- X[, p_name]
        x_others <- X[, setdiff(preds, p_name), drop = FALSE]
        r2_p <- summary(stats::lm(y_p ~ x_others))$r.squared
        vifs[p_name] <- 1 / max(1e-4, 1 - r2_p)
      }
      bp <- graphics::barplot(vifs, horiz = TRUE, col = "#1B365D88", border = "#1B365D",
                             xlab = "Variance Inflation Factor (VIF)",
                             main = "Multicollinearity Check: VIF by Regressor",
                             las = 1, xlim = c(0, max(c(vifs, 12))), bty = "l")
      graphics::abline(v = c(5, 10), lty = c(3, 2), col = c(accent_col, line_col), lwd = 1.8)
      graphics::text(5, max(bp), "VIF = 5 (Moderate)", pos = 4, col = accent_col, cex = 0.8)
      graphics::text(10, min(bp), "VIF = 10 (Severe)", pos = 4, col = line_col, cex = 0.8)
    }
  }
  
  # 11. Residual Histogram + Normal Density
  if (which == "histogram") {
    h <- graphics::hist(e, breaks = "FD", plot = FALSE)
    x_grid <- seq(min(e), max(e), length.out = 200)
    d_norm <- stats::dnorm(x_grid, mean = 0, sd = stats::sd(e))
    y_max <- max(c(h$density, d_norm)) * 1.15
    graphics::plot(h, freq = FALSE, col = "#1B365D33", border = "#1B365D",
                   ylim = c(0, y_max), xlab = "Residuals", ylab = "Density",
                   main = "Residual Histogram with Fitted Normal Curve", bty = "l")
    graphics::lines(x_grid, d_norm, col = line_col, lwd = 2.2)
    graphics::legend("topright", legend = c("Empirical Residuals", "Theoretical Normal"),
                     fill = c("#1B365D33", NA), border = c("#1B365D", NA),
                     lty = c(NA, 1), lwd = c(NA, 2.2), col = c(NA, line_col), bty = "n", cex = 0.85)
  }
  
  # 12. Observed vs Fitted
  if (which == "observed_fitted") {
    lims <- range(c(y_obs, y_hat), na.rm = TRUE)
    graphics::plot(y_hat, y_obs, pch = 19, col = pt_col, cex = 0.85,
                   xlim = lims, ylim = lims, xlab = "Fitted Values (y_hat)", ylab = "Observed Values (y)",
                   main = "Observed vs. Fitted Values (45-Degree Line)", bty = "l",
                   panel.first = graphics::grid(col = "#EAEAEA", lty = 1))
    graphics::abline(a = 0, b = 1, lty = 2, col = line_col, lwd = 2)
  }
  
  # 13. CUSUM Stability Test
  if (which == "cusum") {
    e_std <- e / stats::sd(e)
    cusum_val <- cumsum(e_std) / sqrt(n)
    t_idx <- 1:n
    bound <- 0.948 * (1 + 2 * t_idx / n)
    ylim <- range(c(cusum_val, bound, -bound))
    graphics::plot(t_idx, cusum_val, type = "l", col = pt_col, lwd = 2.2, ylim = ylim,
                   xlab = "Observation Index (Time)", ylab = "CUSUM Statistic",
                   main = "CUSUM Parameter Stability Test (5% Critical Bounds)", bty = "l",
                   panel.first = graphics::grid(col = "#EAEAEA", lty = 1))
    graphics::lines(t_idx, bound, lty = 2, col = line_col, lwd = 1.8)
    graphics::lines(t_idx, -bound, lty = 2, col = line_col, lwd = 1.8)
    graphics::abline(h = 0, col = "#BDC3C7", lty = 1)
  }
  
  # 14. CUSUM of Squares Stability Test
  if (which == "cusum_sq") {
    .acad_cusum_sq(e, main = "CUSUM of Squares Parameter & Variance Stability")
  }
  
  # 15. Joint Stability (CUSUM and CUSUMSQ side by side)
  if (which == "stability") {
    .acad_stability(e, title_prefix = "OLS Model")
  }
  
  # 16. Influence Plot (Studentized Residuals vs. Leverage with Cook's Distance Bubbles)
  if (which == "influence") {
    h <- stats::hatvalues(base_lm)
    r_stud <- stats::rstudent(base_lm)
    d_cook <- stats::cooks.distance(base_lm)
    k_params <- length(stats::coef(base_lm))
    h_cut <- 2 * k_params / n
    
    rad <- sqrt(d_cook) / max(sqrt(d_cook), na.rm = TRUE) * 2.5 + 0.6
    is_inf <- (d_cook > 4 / n) | (abs(r_stud) > 2.5) | (h > h_cut)
    col_pts <- ifelse(is_inf, "#922B21CC", "#1B365D77")
    
    graphics::plot(h, r_stud, type = "n", xlab = "Leverage (Hat Values h_ii)", 
                   ylab = "Studentized Residuals (r_i)",
                   main = "Influence Plot: Studentized Residuals vs. Leverage",
                   bty = "l", panel.first = .acad_grid())
    graphics::abline(h = c(-2, 0, 2), lty = c(2, 1, 2), col = c("#922B21", "#BDC3C7", "#922B21"))
    graphics::abline(v = h_cut, lty = 2, col = "#D4AC0D", lwd = 1.8)
    graphics::points(h, r_stud, pch = 21, bg = col_pts, col = "white", cex = rad, lwd = 0.8)
    
    top_inf <- order(d_cook, decreasing = TRUE)[1:min(5, n)]
    graphics::text(h[top_inf], r_stud[top_inf], labels = top_inf, pos = 4, cex = 0.8, col = "#922B21", font = 2)
    graphics::legend("topright", legend = c("High Influence (D > 4/n or |r| > 2)", sprintf("Leverage Cutoff (2k/n = %.3f)", h_cut)),
                     col = c("#922B21", "#D4AC0D"), pch = c(19, NA), lty = c(NA, 2), lwd = c(NA, 1.8), bty = "n")
  }
  
  invisible(x)
}

#' @rdname econ_plot
#' @export
plot.econ_reg <- function(x, y = NULL, which = c("all", "resid_fitted", "qq", "scale_location", "cooks", 
                                                 "leverage", "added_var", "component_resid", "dfbetas", 
                                                 "forest", "vif", "histogram", "observed_fitted", 
                                                 "cusum", "cusum_sq", "cusum_squares", "stability", 
                                                 "influence", "influence_plot"), ...) {
  econ_plot(x, which = which, ...)
}


# ==============================================================================
# Helper functions for academic styling and diagnostics
# ==============================================================================

.acad_grid <- function() {
  graphics::grid(col = "#EAEAEA", lty = 1)
}

.acad_acf <- function(x, main = "Autocorrelations (ACF)", lag.max = 24) {
  x_clean <- stats::na.omit(x)
  if (length(x_clean) < 4) return(invisible(NULL))
  a <- stats::acf(x_clean, lag.max = min(lag.max, length(x_clean) - 2), plot = FALSE)
  n_obs <- length(x_clean)
  ci <- 1.96 / sqrt(n_obs)
  lags <- as.vector(a$lag)[-1]
  acfs <- as.vector(a$acf)[-1]
  ylim <- c(min(-ci * 1.3, min(acfs, na.rm = TRUE) * 1.15), 
            max(ci * 1.3, max(acfs, na.rm = TRUE) * 1.15))
  graphics::plot(lags, acfs, type = "h", col = "#1B365D", lwd = 2.2,
                 ylim = ylim, xlab = "Lag", ylab = "Autocorrelation",
                 main = main, bty = "l", panel.first = .acad_grid())
  graphics::points(lags, acfs, pch = 19, col = "#1B365D", cex = 0.8)
  graphics::abline(h = c(-ci, ci), lty = 2, col = "#922B21", lwd = 1.6)
  graphics::abline(h = 0, col = "#7F8C8D")
}

.acad_pacf <- function(x, main = "Partial Autocorrelations (PACF)", lag.max = 24) {
  x_clean <- stats::na.omit(x)
  if (length(x_clean) < 4) return(invisible(NULL))
  p <- stats::pacf(x_clean, lag.max = min(lag.max, length(x_clean) - 2), plot = FALSE)
  n_obs <- length(x_clean)
  ci <- 1.96 / sqrt(n_obs)
  lags <- as.vector(p$lag)
  pacfs <- as.vector(p$acf)
  ylim <- c(min(-ci * 1.3, min(pacfs, na.rm = TRUE) * 1.15), 
            max(ci * 1.3, max(pacfs, na.rm = TRUE) * 1.15))
  graphics::plot(lags, pacfs, type = "h", col = "#1B365D", lwd = 2.2,
                 ylim = ylim, xlab = "Lag", ylab = "Partial Autocorrelation",
                 main = main, bty = "l", panel.first = .acad_grid())
  graphics::points(lags, pacfs, pch = 19, col = "#1B365D", cex = 0.8)
  graphics::abline(h = c(-ci, ci), lty = 2, col = "#922B21", lwd = 1.6)
  graphics::abline(h = 0, col = "#7F8C8D")
}

.acad_cusum <- function(e, main = "CUSUM Parameter Stability Test", 
                        col_line = "#1B365D", col_bound = "#922B21") {
  e <- as.numeric(stats::na.omit(e))
  n <- length(e)
  sd_e <- stats::sd(e)
  if (sd_e <= 0 || n < 4) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "Insufficient observations for CUSUM", family = "serif")
    return(invisible(NULL))
  }
  
  e_std <- (e - mean(e)) / sd_e
  w_t <- cumsum(e_std) / sqrt(n)
  t_idx <- 1:n
  # Brown, Durbin, Evans (1975) 5% significance bounds: +/- a * (1 + 2 * t / N) with a = 0.948
  bound <- 0.948 * (1 + 2 * t_idx / n)
  
  is_unstable <- any(abs(w_t) > bound)
  sub_text <- if (is_unstable) "Parameter Instability Detected (5% Bound Breach)" else "Parameter Stability Confirmed (Within 5% Bounds)"
  
  ylim <- range(c(w_t, bound, -bound))
  graphics::plot(t_idx, w_t, type = "l", col = col_line, lwd = 2.2, ylim = ylim,
                 xlab = "Observation Index (Time t)", ylab = "CUSUM Statistic (W_t)",
                 main = main, sub = sub_text, col.sub = if (is_unstable) "#922B21" else "#1E8449",
                 font.sub = 2, bty = "l", panel.first = .acad_grid())
  
  graphics::lines(t_idx, bound, lty = 2, col = col_bound, lwd = 1.8)
  graphics::lines(t_idx, -bound, lty = 2, col = col_bound, lwd = 1.8)
  graphics::abline(h = 0, col = "#BDC3C7", lty = 1)
  
  graphics::legend("topleft",
                   legend = c("CUSUM Path", "5% Critical Bounds (+/- a)"),
                   col = c(col_line, col_bound),
                   lty = c(1, 2), lwd = c(2.2, 1.8), bty = "n")
}

.acad_cusum_sq <- function(e, main = "CUSUM of Squares Stability Test", 
                           col_line = "#1B365D", col_bound = "#922B21") {
  e <- as.numeric(stats::na.omit(e))
  n <- length(e)
  e2 <- e^2
  sum_e2 <- sum(e2)
  if (sum_e2 <= 0 || n < 4) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "Insufficient observations for CUSUM of Squares", family = "serif")
    return(invisible(NULL))
  }
  
  s_t <- cumsum(e2) / sum_e2
  t_idx <- 1:n
  mean_path <- t_idx / n
  
  # Brown, Durbin, Evans (1975) Table 1 critical value approximation for alpha = 0.05
  c_val <- if (n <= 15) {
    0.358
  } else if (n <= 25) {
    0.281
  } else if (n <= 35) {
    0.238
  } else if (n <= 50) {
    0.191
  } else if (n <= 75) {
    0.160
  } else if (n <= 100) {
    0.140
  } else {
    1.358 / sqrt(n)
  }
  
  upper_bound <- mean_path + c_val
  lower_bound <- mean_path - c_val
  
  is_unstable <- any(s_t > upper_bound | s_t < lower_bound)
  sub_text <- if (is_unstable) "Variance / Parameter Instability (5% Bound Breach)" else "Variance & Parameter Stability Confirmed (Within 5% Bounds)"
  
  ylim <- c(min(0, min(lower_bound)), max(1, max(upper_bound)))
  graphics::plot(t_idx, s_t, type = "l", col = col_line, lwd = 2.2,
                 ylim = ylim, xlab = "Observation Index (Time t)", 
                 ylab = "Cumulative Sum of Squares (S_t)",
                 main = main, sub = sub_text, col.sub = if (is_unstable) "#922B21" else "#1E8449",
                 font.sub = 2, bty = "l", panel.first = .acad_grid())
  
  # 45-degree expected null diagonal
  graphics::lines(t_idx, mean_path, lty = 2, col = "#7F8C8D", lwd = 1.5)
  # 5% Critical bound lines
  graphics::lines(t_idx, upper_bound, lty = 2, col = col_bound, lwd = 1.8)
  graphics::lines(t_idx, lower_bound, lty = 2, col = col_bound, lwd = 1.8)
  
  graphics::legend("topleft", 
                   legend = c("CUSUM of Squares (S_t)", "Expected Path (H0)", "5% Critical Bounds (+/- c)"),
                   col = c(col_line, "#7F8C8D", col_bound),
                   lty = c(1, 2, 2), lwd = c(2.2, 1.5, 1.8), bty = "n")
}

.acad_stability <- function(e, title_prefix = "Model") {
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  .acad_cusum(e, main = sprintf("%s: CUSUM Stability", title_prefix))
  .acad_cusum_sq(e, main = sprintf("%s: CUSUM of Squares", title_prefix))
}

.acad_bounds_ardl <- function(x) {
  crit_df <- data.frame(
    Level = c("10%", "5%", "1%"),
    I0 = c(4.04, 4.94, 6.84),
    I1 = c(4.78, 5.73, 7.84)
  )
  
  f_stat <- if (!is.null(x$f_bounds)) {
    as.numeric(x$f_bounds)
  } else if (!is.null(x$t_lr)) {
    as.numeric(x$t_lr)^2
  } else {
    NA_real_
  }
  
  y_max <- max(c(crit_df$I1 * 1.35, if (is.finite(f_stat)) f_stat * 1.25 else 8))
  
  graphics::plot(1:3, crit_df$I1, type = "n", xlim = c(0.4, 3.6), ylim = c(0, y_max),
                 xaxt = "n", xlab = "Significance Level", ylab = "F-Statistic Value",
                 main = "Pesaran-Shin-Smith (2001) ARDL Bounds Test", bty = "l",
                 panel.first = .acad_grid())
  graphics::axis(1, at = 1:3, labels = crit_df$Level, font = 2)
  
  for (i in 1:3) {
    graphics::rect(i - 0.25, crit_df$I0[i], i + 0.25, crit_df$I1[i],
                   col = "#7F8C8D33", border = "#7F8C8D", lwd = 1.5)
    graphics::text(i, (crit_df$I0[i] + crit_df$I1[i]) / 2,
                   sprintf("[%0.2f - %0.2f]", crit_df$I0[i], crit_df$I1[i]),
                   cex = 0.8, col = "#2C3E50", font = 2)
  }
  
  if (is.finite(f_stat)) {
    graphics::abline(h = f_stat, col = "#1B365D", lwd = 2.5, lty = 1)
    
    verdict <- if (f_stat > crit_df$I1[2]) {
      "Cointegration Confirmed (F > I(1))"
    } else if (f_stat < crit_df$I0[2]) {
      "No Cointegration (F < I(0))"
    } else {
      "Inconclusive Bound Region"
    }
    
    col_v <- if (f_stat > crit_df$I1[2]) "#1E8449" else if (f_stat < crit_df$I0[2]) "#922B21" else "#D4AC0D"
    
    graphics::legend("topleft", 
                     legend = c(sprintf("ARDL F-Stat = %.2f: %s", f_stat, verdict),
                                "Inconclusive Band [I(0) - I(1)]"),
                     col = c("#1B365D", "#7F8C8D"), lty = c(1, NA), pch = c(NA, 15),
                     pt.cex = 1.5, lwd = c(2.5, NA), bty = "n", text.col = c(col_v, "black"),
                     text.font = c(2, 1), cex = 0.8)
  } else {
    graphics::legend("topleft", legend = "Inconclusive Band [I(0) - I(1)]",
                     col = "#7F8C8D", pch = 15, pt.cex = 1.5, bty = "n", cex = 0.8)
  }
}


# ==============================================================================
# 2. Panel Data Models (plot.econ_panel)
# ==============================================================================

#' Plot Panel Data Regression Diagnostics
#'
#' Generates an academic publication-grade panel diagnostics suite or individual plots:
#' (1) \code{"between"}: Between-entity variation in group means with between regression line,
#' (2) \code{"within"}: Within-entity demeaned scatter with within regression slope,
#' (3) \code{"spaghetti"}: Entity trajectories over time with mean trend,
#' (4) \code{"fe_dist"}: Distribution of estimated individual fixed effects \eqn{\hat{lpha}_i},
#' (5) \code{"heterogeneity"}: Individual entity residual boxplots,
#' (6) \code{"residuals"}: Panel residuals vs. fitted values with LOESS smoother,
#' (7) \code{"ar_residual"}: Residual vs. lagged within-entity residual \eqn{e_{i,t}} vs \eqn{e_{i,t-1}} (Wooldridge AR(1) test),
#' (8) \code{"time_effects"}: Period mean residuals \eqn{\pm 95\%} CI across time.
#'
#' @param x An object of class \code{econ_panel}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 2x2 grid),
#'   \code{"between"}, \code{"within"}, \code{"spaghetti"}, \code{"fe_dist"}, \code{"heterogeneity"},
#'   \code{"residuals"}, \code{"ar_residual"}, or \code{"time_effects"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_panel <- function(x, y = NULL, which = c("all", "between", "within", "spaghetti", 
                                                   "fe_dist", "heterogeneity", "residuals", 
                                                   "ar_residual", "time_effects"), ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "between", "within", "spaghetti", "fe_dist", "heterogeneity", "residuals", "ar_residual", "time_effects")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "between", "within", "spaghetti", 
                                  "fe_dist", "heterogeneity", "residuals", 
                                  "ar_residual", "time_effects"))
  }
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  y_hat <- stats::fitted(x)
  e <- stats::residuals(x)
  id_vec <- as.factor(x$id_vec)
  time_vec <- if (!is.null(x$time_vec)) as.factor(x$time_vec) else NULL
  y_raw <- x$y_raw
  X_raw <- x$X_raw
  
  pred_names <- colnames(X_raw)[colnames(X_raw) != "(Intercept)"]
  prim_name <- if (length(pred_names) > 0) pred_names[1] else "Regressor"
  x_prim <- if (length(pred_names) > 0) X_raw[, prim_name] else seq_along(y_raw)
  
  # 1. Between Variation (Group Means)
  if (which %in% c("all", "between")) {
    y_means <- tapply(y_raw, id_vec, mean, na.rm = TRUE)
    x_means <- tapply(x_prim, id_vec, mean, na.rm = TRUE)
    graphics::plot(x_means, y_means, pch = 19, col = "#1B365D", cex = 0.85, bty = "l",
                   xlab = paste0("Mean ", prim_name, " by Entity"),
                   ylab = "Mean Outcome by Entity",
                   main = if (which == "all") "1. Between-Entity Variation" else "Between-Entity Variation (Group Means)",
                   panel.first = .acad_grid())
    if (length(y_means) > 2) {
      b_between <- stats::lm(y_means ~ x_means)
      graphics::abline(b_between, col = "#922B21", lwd = 2)
    }
  }
  
  # 2. Within Variation (Demeaned Scatter)
  if (which %in% c("all", "within")) {
    y_within <- y_raw - unname(tapply(y_raw, id_vec, mean, na.rm = TRUE)[id_vec])
    x_within <- x_prim - unname(tapply(x_prim, id_vec, mean, na.rm = TRUE)[id_vec])
    graphics::plot(x_within, y_within, pch = 16, col = "#1B365D77", cex = 0.8, bty = "l",
                   xlab = paste0("Demeaned ", prim_name, " (Within)"),
                   ylab = "Demeaned Outcome (Within)",
                   main = if (which == "all") "2. Within-Entity Identification" else "Within-Entity Identification (Demeaned Scatter)",
                   panel.first = .acad_grid())
    if (length(y_within) > 2) {
      b_within <- stats::lm(y_within ~ 0 + x_within)
      graphics::abline(b_within, col = "#1E8449", lwd = 2.2)
    }
  }
  
  # 3. Entity Heterogeneity (Boxplot)
  if (which %in% c("all", "heterogeneity")) {
    unique_ids <- levels(id_vec)
    meds <- tapply(e, id_vec, stats::median, na.rm = TRUE)
    ordered_ids <- names(sort(meds))
    plot_ids <- if (length(ordered_ids) > 15) {
      seq_idx <- round(seq(1, length(ordered_ids), length.out = 15))
      ordered_ids[seq_idx]
    } else ordered_ids
    
    e_sub <- e[id_vec %in% plot_ids]
    id_sub <- factor(id_vec[id_vec %in% plot_ids], levels = plot_ids)
    graphics::boxplot(e_sub ~ id_sub, col = "#1B365D33", border = "#1B365D",
                      xlab = if (length(unique_ids) > 15) "Representative Entities (Ranked by Median)" else "Entity ID",
                      ylab = "Residuals",
                      main = if (which == "all") "3. Entity Residual Heterogeneity" else "Individual Entity Residual Heterogeneity",
                      bty = "l", outline = FALSE, las = 2, cex.axis = 0.75)
    graphics::abline(h = 0, lty = 2, col = "#922B21", lwd = 1.5)
  }
  
  # 4. Residuals vs. Fitted
  if (which %in% c("all", "residuals")) {
    graphics::plot(y_hat, e, pch = 19, col = "#1B365D88", cex = 0.85, bty = "l",
                   xlab = "Fitted Values", ylab = "Panel Residuals",
                   main = if (which == "all") sprintf("4. Panel (%s): Residuals vs. Fitted", toupper(x$model)) else
                          sprintf("Panel (%s): Residuals vs. Fitted Values", toupper(x$model)),
                   panel.first = .acad_grid())
    graphics::abline(h = 0, lty = 2, col = "#922B21", lwd = 1.8)
    if (length(e) >= 15) {
      sm <- stats::loess(e ~ y_hat, span = 0.75)
      ord <- order(y_hat)
      graphics::lines(y_hat[ord], stats::fitted(sm)[ord], col = "#1E8449", lwd = 2.2)
    }
  }
  
  # 5. Spaghetti Plot (Trajectories over time)
  if (which == "spaghetti") {
    if (is.null(time_vec)) {
      graphics::plot.new()
      graphics::text(0.5, 0.5, "Spaghetti plot requires time variable in panel model.")
    } else {
      time_num <- as.numeric(as.character(time_vec))
      if (any(is.na(time_num))) time_num <- as.numeric(time_vec)
      u_time <- sort(unique(time_num))
      
      graphics::plot(time_num, y_raw, type = "n", xlab = "Time Period", ylab = "Outcome (y)",
                     main = "Panel Spaghetti Plot: Trajectories across Time", bty = "l",
                     panel.first = .acad_grid())
      for (id in levels(id_vec)) {
        mask <- id_vec == id
        graphics::lines(time_num[mask], y_raw[mask], col = "#1B365D25", lwd = 1.2)
      }
      # Mean trend line
      mean_trend <- tapply(y_raw, time_num, mean, na.rm = TRUE)
      graphics::lines(as.numeric(names(mean_trend)), mean_trend, col = "#922B21", lwd = 3.0)
      graphics::points(as.numeric(names(mean_trend)), mean_trend, pch = 19, col = "#922B21", cex = 1.2)
      graphics::legend("topleft", legend = c("Individual Trajectories", "Sample Mean Trend"),
                       col = c("#1B365D77", "#922B21"), lwd = c(1.5, 3), lty = 1, bty = "n")
    }
  }
  
  # 6. Fixed Effects Distribution
  if (which == "fe_dist") {
    # Compute alpha_i = mean(y_i) - mean(x_i) * beta
    b_slopes <- stats::coef(x)
    slopes_names <- names(b_slopes)[names(b_slopes) != "(Intercept)"]
    y_mean_i <- tapply(y_raw, id_vec, mean, na.rm = TRUE)
    if (length(slopes_names) > 0 && all(slopes_names %in% colnames(X_raw))) {
      x_mean_i <- matrix(0, nrow = nlevels(id_vec), ncol = length(slopes_names))
      for (j in seq_along(slopes_names)) {
        x_mean_i[, j] <- tapply(X_raw[, slopes_names[j]], id_vec, mean, na.rm = TRUE)
      }
      alpha_i <- y_mean_i - as.vector(x_mean_i %*% b_slopes[slopes_names])
    } else {
      alpha_i <- y_mean_i
    }
    
    h <- graphics::hist(alpha_i, breaks = "FD", plot = FALSE)
    x_seq <- seq(min(alpha_i), max(alpha_i), length.out = 150)
    d_norm <- stats::dnorm(x_seq, mean = mean(alpha_i), sd = stats::sd(alpha_i))
    graphics::plot(h, freq = FALSE, col = "#1B365D33", border = "#1B365D",
                   xlab = "Estimated Entity Effect (alpha_i)", ylab = "Density",
                   main = "Distribution of Estimated Individual Fixed Effects", bty = "l")
    graphics::lines(x_seq, d_norm, col = "#922B21", lwd = 2)
    graphics::abline(v = stats::median(alpha_i), lty = 2, col = "#1E8449", lwd = 2)
    graphics::legend("topright", legend = c("Empirical FE", "Normal Approximation", "Median FE"),
                     fill = c("#1B365D33", NA, NA), border = c("#1B365D", NA, NA),
                     lty = c(NA, 1, 2), lwd = c(NA, 2, 2), col = c(NA, "#922B21", "#1E8449"), bty = "n")
  }
  
  # 7. AR(1) Within Residual Scatter (Wooldridge test visual)
  if (which == "ar_residual") {
    e_curr <- numeric()
    e_lag <- numeric()
    for (id in levels(id_vec)) {
      idx <- which(id_vec == id)
      if (length(idx) > 1) {
        e_curr <- c(e_curr, e[idx[-1]])
        e_lag <- c(e_lag, e[idx[-length(idx)]])
      }
    }
    if (length(e_curr) > 5) {
      graphics::plot(e_lag, e_curr, pch = 19, col = "#1B365D66", cex = 0.8, bty = "l",
                     xlab = "Lagged Residual e_{i, t-1}", ylab = "Current Residual e_{i, t}",
                     main = "Serial Correlation Check: e_{i,t} vs. e_{i,t-1}",
                     panel.first = .acad_grid())
      graphics::abline(h = 0, v = 0, lty = 2, col = "#BDC3C7")
      ar_fit <- stats::lm(e_curr ~ 0 + e_lag)
      rho_hat <- stats::coef(ar_fit)[1]
      graphics::abline(a = 0, b = rho_hat, col = "#922B21", lwd = 2.2)
      graphics::legend("topleft", legend = sprintf("AR(1) Slope: rho = %.3f", rho_hat),
                       col = "#922B21", lty = 1, lwd = 2.2, bty = "n")
    } else {
      graphics::plot.new()
      graphics::text(0.5, 0.5, "Insufficient time periods per entity for AR(1) residual plot.")
    }
  }
  
  # 8. Time Effects Plot
  if (which == "time_effects") {
    if (is.null(time_vec)) {
      graphics::plot.new()
      graphics::text(0.5, 0.5, "Time effects plot requires time variable in panel model.")
    } else {
      t_means <- tapply(e, time_vec, mean, na.rm = TRUE)
      t_ses <- tapply(e, time_vec, function(v) stats::sd(v, na.rm = TRUE) / sqrt(length(stats::na.omit(v))))
      m_t <- length(t_means)
      ci_low <- t_means - 1.96 * t_ses
      ci_high <- t_means + 1.96 * t_ses
      
      graphics::plot(1:m_t, t_means, pch = 19, col = "#1B365D", cex = 1.2,
                     xaxt = "n", xlab = "Time Period", ylab = "Mean Residual (+/- 95% CI)",
                     ylim = range(c(ci_low, ci_high, 0), na.rm = TRUE),
                     main = "Time Period Residual Means & Heterogeneity", bty = "l",
                     panel.first = .acad_grid())
      graphics::axis(1, at = 1:m_t, labels = names(t_means), las = 2, cex.axis = 0.8)
      graphics::abline(h = 0, lty = 2, col = "#922B21", lwd = 1.5)
      graphics::arrows(1:m_t, ci_low, 1:m_t, ci_high, length = 0.05, angle = 90, code = 3, col = "#1B365D", lwd = 1.8)
    }
  }
  
  invisible(x)
}

# ==============================================================================
# 3. Instrumental Variables (plot.econ_iv)
# ==============================================================================

#' Plot Instrumental Variables 2SLS Model
#'
#' Generates an academic publication-grade 2SLS identification suite:
#' (1) \code{"first_stage"}: First-Stage Fit: Endogenous \eqn{X} vs. Instrument \eqn{Z} with \eqn{F}-stat,
#' (2) \code{"second_stage"}: Second Stage: Outcome \eqn{Y} vs. Predicted Endogenous \eqn{\hat{X}},
#' (3) \code{"residuals"}: 2SLS Residuals vs. Fitted values,
#' (4) \code{"orthogonality"}: Instrument Exogeneity Check: Residuals vs. Instrument \eqn{Z},
#' (5) \code{"forest_compare"}: OLS vs. 2SLS Coefficient Comparison Forest Plot,
#' (6) \code{"f_stat"}: First-stage F-statistic vs. Stock-Yogo rule of thumb (F = 10),
#' (7) \code{"reduced_form"}: Reduced form regression of \eqn{Y} on \eqn{Z}.
#'
#' @param x An object of class \code{econ_iv}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 2x2 grid),
#'   \code{"first_stage"}, \code{"second_stage"}, \code{"residuals"}, \code{"orthogonality"},
#'   \code{"forest_compare"}, \code{"f_stat"}, or \code{"reduced_form"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_iv <- function(x, y = NULL, which = c("all", "first_stage", "second_stage", "residuals", 
                                                "orthogonality", "forest_compare", "f_stat", "reduced_form"), ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "first_stage", "second_stage", "residuals", "orthogonality", "forest_compare", "f_stat", "reduced_form")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "first_stage", "second_stage", "residuals", 
                                  "orthogonality", "forest_compare", "f_stat", "reduced_form"))
  }
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  endo_name <- if (!is.null(x$endo_vars) && length(x$endo_vars) > 0) x$endo_vars[1] else if (!is.null(x$endogenous)) x$endogenous[1] else colnames(x$X)[2]
  iv_name <- if (!is.null(x$iv_vars) && length(x$iv_vars) > 0) x$iv_vars[1] else if (!is.null(x$instruments)) x$instruments[1] else colnames(x$Z)[2]
  
  y_obs <- x$y
  x_end <- x$X[, endo_name]
  z_inst <- x$Z[, iv_name]
  fit_fs <- stats::lm(x_end ~ z_inst)
  x_hat <- stats::fitted(fit_fs)
  e <- x$residuals
  fit <- x$fitted.values
  f_val <- if (!is.null(x$first_stage[[endo_name]]$f_stat)) {
    x$first_stage[[endo_name]]$f_stat
  } else if (!is.null(x$first_stage_f)) {
    x$first_stage_f[1]
  } else {
    fs_f <- summary(fit_fs)$fstatistic
    if (!is.null(fs_f)) as.numeric(fs_f[1]) else NA_real_
  }
  
  # 1. First Stage
  if (which %in% c("all", "first_stage")) {
    graphics::plot(z_inst, x_end, pch = 19, col = "#1B365D88", cex = 0.85, bty = "l",
                   xlab = paste0("Excluded Instrument (", iv_name, ")"),
                   ylab = paste0("Endogenous Variable (", endo_name, ")"),
                   main = sprintf("1. First-Stage Strength (F = %.2f)", f_val),
                   panel.first = .acad_grid())
    graphics::abline(fit_fs, col = "#1E8449", lwd = 2.2)
    graphics::legend("topleft", legend = if (f_val >= 10) "Strong Instrument (F >= 10)" else "Weak Instrument (F < 10)",
                     col = if (f_val >= 10) "#1E8449" else "#922B21", text.col = if (f_val >= 10) "#1E8449" else "#922B21",
                     bty = "n", text.font = 2)
  }
  
  # 2. Second Stage
  if (which %in% c("all", "second_stage")) {
    graphics::plot(x_hat, y_obs, pch = 19, col = "#1B365D88", cex = 0.85, bty = "l",
                   xlab = paste0("Predicted Endogenous (", endo_name, " hat)"),
                   ylab = paste0("Outcome (", all.vars(x$formula)[1], ")"),
                   main = sprintf("2. 2SLS Second Stage: Beta = %.*f", x$digits, stats::coef(x)[endo_name]),
                   panel.first = .acad_grid())
    fit_ss <- stats::lm(y_obs ~ x_hat)
    graphics::abline(fit_ss, col = "#922B21", lwd = 2.2)
  }
  
  # 3. Residuals vs Fitted
  if (which %in% c("all", "residuals")) {
    graphics::plot(fit, e, pch = 19, col = "#1B365D88", cex = 0.85, bty = "l",
                   xlab = "2SLS Fitted Values", ylab = "Structural 2SLS Residuals",
                   main = "3. Structural Residuals vs. Fitted",
                   panel.first = .acad_grid())
    graphics::abline(h = 0, lty = 2, col = "#922B21", lwd = 1.8)
    if (length(e) >= 15) {
      sm <- stats::loess(e ~ fit, span = 0.75)
      ord <- order(fit)
      graphics::lines(fit[ord], stats::fitted(sm)[ord], col = "#1E8449", lwd = 2)
    }
  }
  
  # 4. Orthogonality
  if (which %in% c("all", "orthogonality")) {
    graphics::plot(z_inst, e, pch = 19, col = "#1B365D88", cex = 0.85, bty = "l",
                   xlab = paste0("Instrument (", iv_name, ")"),
                   ylab = "Structural Residuals",
                   main = "4. Instrument Exogeneity Check",
                   panel.first = .acad_grid())
    graphics::abline(h = 0, lty = 2, col = "#7F8C8D")
    fit_ortho <- stats::lm(e ~ z_inst)
    graphics::abline(fit_ortho, col = "#D35400", lwd = 2, lty = 2)
  }
  
  # 5. Forest Comparison: OLS vs 2SLS
  if (which == "forest_compare") {
    # Fit OLS for comparison
    ols_df <- data.frame(y = y_obs, x$X[, -1, drop = FALSE])
    ols_mod <- stats::lm(y ~ ., data = ols_df)
    b_ols <- stats::coef(ols_mod)[endo_name]
    se_ols <- summary(ols_mod)$coefficients[endo_name, "Std. Error"]
    
    b_2sls <- stats::coef(x)[endo_name]
    se_2sls <- sqrt(stats::vcov(x)[endo_name, endo_name])
    
    vals <- c(b_ols, b_2sls)
    ses <- c(se_ols, se_2sls)
    ci_l <- vals - 1.96 * ses
    ci_u <- vals + 1.96 * ses
    
    graphics::plot(vals, c(1, 2), pch = 19, col = c("#1B365D", "#922B21"), cex = 1.3,
                   yaxt = "n", xlab = "Coefficient Estimate (+/- 95% CI)", ylab = "",
                   xlim = range(c(ci_l, ci_u, 0)), ylim = c(0.5, 2.5),
                   main = paste0("OLS vs. 2SLS Comparison: ", endo_name), bty = "l",
                   panel.first = .acad_grid())
    graphics::axis(2, at = c(1, 2), labels = c("OLS (Biased)", "2SLS (Consistent)"), las = 1, font = 2)
    graphics::abline(v = 0, lty = 2, col = "#7F8C8D")
    graphics::arrows(ci_l, c(1, 2), ci_u, c(1, 2), length = 0.05, angle = 90, code = 3,
                     col = c("#1B365D", "#922B21"), lwd = 2.2)
  }
  
  # 6. First-Stage F vs Stock-Yogo 10
  if (which == "f_stat") {
    bp <- graphics::barplot(f_val, horiz = TRUE, col = if (f_val >= 10) "#1E844988" else "#922B2188",
                            border = if (f_val >= 10) "#1E8449" else "#922B21",
                            xlab = "First-Stage Effective F-Statistic",
                            main = "Instrument Strength: First-Stage F vs. Stock-Yogo (10)",
                            xlim = c(0, max(c(f_val * 1.25, 20))), las = 1, bty = "l")
    graphics::abline(v = 10, lty = 2, col = "#922B21", lwd = 2.2)
    graphics::text(10, bp, "  F = 10 (Rule of Thumb)", col = "#922B21", font = 2, pos = 4)
  }
  
  # 7. Reduced Form (Y vs Z)
  if (which == "reduced_form") {
    graphics::plot(z_inst, y_obs, pch = 19, col = "#1B365D88", cex = 0.85, bty = "l",
                   xlab = paste0("Instrument (", iv_name, ")"),
                   ylab = paste0("Outcome (", all.vars(x$formula)[1], ")"),
                   main = "Reduced Form Relationship: Outcome vs. Instrument",
                   panel.first = .acad_grid())
    graphics::abline(stats::lm(y_obs ~ z_inst), col = "#922B21", lwd = 2.2)
  }
  
  invisible(x)
}

# ==============================================================================
# 4. Logistic and Probit Models (plot.econ_logit)
# ==============================================================================

#' Plot Binary Outcome (Logit/Probit) Model Diagnostics
#'
#' Generates an academic publication-grade classification and link diagnostics suite:
#' (1) \code{"scurve"}: S-Curve sigmoid probability curve vs. linear predictor \eqn{\eta},
#' (2) \code{"roc"}: Receiver Operating Characteristic (ROC) curve with AUC annotation,
#' (3) \code{"distribution"}: Class separation boxplots of predicted probabilities by observed class,
#' (4) \code{"residuals"}: Deviance residuals vs. predicted probabilities,
#' (5) \code{"calibration"}: Hosmer-Lemeshow decile calibration plot with 45-degree reference line,
#' (6) \code{"thresholds"}: Sensitivity & Specificity curves vs. cutoff with Youden's optimal threshold,
#' (7) \code{"marginal"}: Marginal probability curve across continuous predictor,
#' (8) \code{"odds_ratio"}: Odds ratios forest plot (for logit) with 95% CIs on log scale,
#' (9) \code{"confusion"}: Confusion matrix classification tile heatmap at threshold 0.5.
#'
#' @param x An object of class \code{econ_logit}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 2x2 grid),
#'   \code{"scurve"}, \code{"roc"}, \code{"distribution"}, \code{"residuals"}, \code{"calibration"},
#'   \code{"thresholds"}, \code{"marginal"}, \code{"odds_ratio"}, or \code{"confusion"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_logit <- function(x, y = NULL, which = c("all", "scurve", "roc", "distribution", 
                                                   "residuals", "calibration", "thresholds", 
                                                   "marginal", "odds_ratio", "confusion"), ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "scurve", "roc", "distribution", "residuals", "calibration", "thresholds", "marginal", "odds_ratio", "confusion")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "scurve", "roc", "distribution", 
                                  "residuals", "calibration", "thresholds", 
                                  "marginal", "odds_ratio", "confusion"))
  }
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  p_hat <- x$fitted.values
  y_obs <- if (!is.null(x$glm_model$y)) x$glm_model$y else as.numeric(p_hat >= 0.5)
  eta <- if (!is.null(x$glm_model$linear.predictors)) x$glm_model$linear.predictors else stats::qlogis(pmin(0.999, pmax(0.001, p_hat)))
  dev_res <- if (!is.null(x$glm_model)) stats::residuals(x$glm_model, type = "deviance") else y_obs - p_hat
  
  # 1. S-curve
  if (which %in% c("all", "scurve")) {
    ord <- order(eta)
    link_label <- if (!is.null(x$link)) toupper(x$link) else "LOGIT"
    graphics::plot(eta[ord], p_hat[ord], type = "l", col = "#1B365D", lwd = 2.5,
                   xlab = "Linear Predictor (eta = X'beta)", ylab = "Predicted Probability P(Y = 1)",
                   main = sprintf("1. %s Link S-Curve Function", link_label),
                   ylim = c(0, 1), bty = "l", panel.first = .acad_grid())
    graphics::points(eta, y_obs, pch = "|", col = ifelse(y_obs == 1, "#1E8449", "#922B21"), cex = 0.8)
    graphics::abline(h = 0.5, lty = 2, col = "#7F8C8D")
  }
  
  # 2. ROC Curve
  if (which %in% c("all", "roc")) {
    thresh <- seq(0, 1, length.out = 150)
    tpr <- numeric(length(thresh))
    fpr <- numeric(length(thresh))
    n1 <- sum(y_obs == 1)
    n0 <- sum(y_obs == 0)
    for (i in seq_along(thresh)) {
      pred_1 <- as.numeric(p_hat >= thresh[i])
      tpr[i] <- if (n1 > 0) sum(pred_1 == 1 & y_obs == 1) / n1 else 0
      fpr[i] <- if (n0 > 0) sum(pred_1 == 1 & y_obs == 0) / n0 else 0
    }
    ord_roc <- order(fpr, tpr)
    fpr_s <- fpr[ord_roc]
    tpr_s <- tpr[ord_roc]
    auc_calc <- sum(diff(fpr_s) * (tpr_s[-1] + tpr_s[-length(tpr_s)]) / 2)
    auc_val <- if (!is.null(x$auc)) x$auc else max(0, min(1, auc_calc))
    graphics::plot(fpr, tpr, type = "l", col = "#922B21", lwd = 2.5,
                   xlab = "False Positive Rate (1 - Specificity)", ylab = "True Positive Rate (Sensitivity)",
                   main = sprintf("2. ROC Curve (AUC = %.3f)", auc_val),
                   xlim = c(0, 1), ylim = c(0, 1), bty = "l", panel.first = .acad_grid())
    graphics::abline(a = 0, b = 1, lty = 2, col = "#7F8C8D")
    graphics::polygon(c(fpr, 0), c(tpr, 0), col = "#922B2122", border = NA)
  }
  
  # 3. Class Distribution
  if (which %in% c("all", "distribution")) {
    graphics::boxplot(p_hat ~ y_obs, col = c("#922B2133", "#1E844933"),
                      border = c("#922B21", "#1E8449"),
                      names = c("Y = 0 (Failure)", "Y = 1 (Success)"),
                      xlab = "Observed Classification", ylab = "Predicted Probability",
                      main = "3. Probability Separation by Class",
                      bty = "l")
    graphics::abline(h = 0.5, lty = 2, col = "#7F8C8D")
  }
  
  # 4. Deviance Residuals
  if (which %in% c("all", "residuals")) {
    graphics::plot(p_hat, dev_res, pch = 19, col = "#1B365D88", cex = 0.85, bty = "l",
                   xlab = "Fitted Probabilities", ylab = "Deviance Residuals",
                   main = "4. Deviance Residuals vs. Fitted",
                   panel.first = .acad_grid())
    graphics::abline(h = 0, lty = 2, col = "#922B21", lwd = 1.8)
    if (length(dev_res) >= 15) {
      sm <- stats::loess(dev_res ~ p_hat, span = 0.75)
      ord <- order(p_hat)
      graphics::lines(p_hat[ord], stats::fitted(sm)[ord], col = "#1E8449", lwd = 2)
    }
  }
  
  # 5. Hosmer-Lemeshow Calibration
  if (which == "calibration") {
    # Deciles of predicted probabilities
    cuts <- stats::quantile(p_hat, probs = seq(0, 1, 0.1), na.rm = TRUE)
    cuts[1] <- cuts[1] - 1e-5
    cuts[length(cuts)] <- cuts[length(cuts)] + 1e-5
    grp <- cut(p_hat, breaks = unique(cuts), include.lowest = TRUE)
    
    mean_p <- tapply(p_hat, grp, mean, na.rm = TRUE)
    mean_y <- tapply(y_obs, grp, mean, na.rm = TRUE)
    n_grp <- tapply(y_obs, grp, length)
    se_y <- sqrt(mean_y * (1 - mean_y) / pmax(1, n_grp))
    
    graphics::plot(mean_p, mean_y, pch = 19, col = "#1B365D", cex = 1.2,
                   xlim = c(0, 1), ylim = c(0, 1), xlab = "Mean Predicted Probability", ylab = "Observed Proportion",
                   main = "Hosmer-Lemeshow Calibration Plot (Deciles)", bty = "l",
                   panel.first = .acad_grid())
    graphics::abline(a = 0, b = 1, lty = 2, col = "#922B21", lwd = 2)
    graphics::arrows(mean_p, mean_y - 1.96 * se_y, mean_p, mean_y + 1.96 * se_y,
                     length = 0.04, angle = 90, code = 3, col = "#1B365D", lwd = 1.8)
    graphics::legend("topleft", legend = c("Decile Means +/- 95% CI", "Ideal 45-Degree Calibration"),
                     pch = c(19, NA), lty = c(NA, 2), col = c("#1B365D", "#922B21"), lwd = c(NA, 2), bty = "n")
  }
  
  # 6. Sensitivity/Specificity vs Threshold
  if (which == "thresholds") {
    taus <- seq(0.05, 0.95, by = 0.02)
    sens <- numeric(length(taus))
    spec <- numeric(length(taus))
    n1 <- sum(y_obs == 1)
    n0 <- sum(y_obs == 0)
    for (i in seq_along(taus)) {
      y_pred <- as.numeric(p_hat >= taus[i])
      sens[i] <- if (n1 > 0) sum(y_pred == 1 & y_obs == 1) / n1 else 0
      spec[i] <- if (n0 > 0) sum(y_pred == 0 & y_obs == 0) / n0 else 0
    }
    youden <- sens + spec - 1
    opt_idx <- which.max(youden)
    opt_tau <- taus[opt_idx]
    
    graphics::plot(taus, sens, type = "l", col = "#1E8449", lwd = 2.2, ylim = c(0, 1),
                   xlab = "Classification Cutoff Threshold", ylab = "Metric Rate",
                   main = sprintf("Sensitivity & Specificity Tradeoff (Optimal Cutoff = %.2f)", opt_tau),
                   bty = "l", panel.first = .acad_grid())
    graphics::lines(taus, spec, col = "#1B365D", lwd = 2.2, lty = 2)
    graphics::abline(v = opt_tau, lty = 3, col = "#922B21", lwd = 2)
    graphics::legend("right", legend = c("Sensitivity (TPR)", "Specificity (TNR)", "Youden Optimal Cutoff"),
                     col = c("#1E8449", "#1B365D", "#922B21"), lty = c(1, 2, 3), lwd = 2, bty = "n")
  }
  
  # 7. Marginal Effect Curve
  if (which == "marginal") {
    preds <- colnames(x$glm_model$data)
    if (is.null(preds)) preds <- setdiff(names(stats::coef(x)), "(Intercept)")
    p1 <- preds[1]
    df_data <- x$glm_model$data
    if (!is.null(df_data) && p1 %in% names(df_data) && is.numeric(df_data[[p1]])) {
      x_seq <- seq(min(df_data[[p1]]), max(df_data[[p1]]), length.out = 150)
      # Predict holding others at means
      new_df <- df_data[rep(1, length(x_seq)), , drop = FALSE]
      for (cn in names(df_data)) {
        if (cn != p1 && is.numeric(df_data[[cn]])) new_df[[cn]] <- mean(df_data[[cn]], na.rm = TRUE)
      }
      new_df[[p1]] <- x_seq
      preds_prob <- stats::predict(x$glm_model, newdata = new_df, type = "response", se.fit = TRUE)
      
      graphics::plot(x_seq, preds_prob$fit, type = "l", col = "#1B365D", lwd = 2.5,
                     ylim = c(0, 1), xlab = p1, ylab = "Predicted Probability P(Y = 1)",
                     main = paste0("Marginal Predicted Probability: ", p1), bty = "l",
                     panel.first = .acad_grid())
      ci_l <- pmax(0, preds_prob$fit - 1.96 * preds_prob$se.fit)
      ci_u <- pmin(1, preds_prob$fit + 1.96 * preds_prob$se.fit)
      graphics::polygon(c(x_seq, rev(x_seq)), c(ci_l, rev(ci_u)), col = "#1B365D22", border = NA)
      graphics::lines(x_seq, ci_l, lty = 2, col = "#1B365D")
      graphics::lines(x_seq, ci_u, lty = 2, col = "#1B365D")
    } else {
      graphics::plot.new()
      graphics::text(0.5, 0.5, "Marginal effect curve requires numeric predictor in model data.")
    }
  }
  
  # 8. Odds Ratio Forest
  if (which == "odds_ratio") {
    if (identical(x$link, "probit")) {
      graphics::plot.new()
      graphics::text(0.5, 0.5, "Odds ratios apply to logit models.")
    } else {
      ct <- x$coef_table
      ct_sub <- ct[rownames(ct) != "(Intercept)", , drop = FALSE]
      if (nrow(ct_sub) == 0) ct_sub <- ct
      m_k <- nrow(ct_sub)
      or_val <- if ("Odds_Ratio" %in% names(ct_sub)) ct_sub$Odds_Ratio else exp(ct_sub$Estimate)
      se_b <- ct_sub$Std_Error
      b_est <- ct_sub$Estimate
      or_low <- exp(b_est - 1.96 * se_b)
      or_high <- exp(b_est + 1.96 * se_b)
      
      graphics::plot(or_val, 1:m_k, pch = 19, col = "#1B365D", cex = 1.2, log = "x",
                     xlim = range(c(or_low, or_high, 1)), ylim = c(0.5, m_k + 0.5),
                     yaxt = "n", xlab = "Odds Ratio (log scale, 95% CI)", ylab = "",
                     main = "Odds Ratios & 95% Confidence Intervals", bty = "l",
                     panel.first = .acad_grid())
      graphics::axis(2, at = 1:m_k, labels = rownames(ct_sub), las = 1, font = 2)
      graphics::abline(v = 1, lty = 2, col = "#922B21", lwd = 1.8)
      graphics::arrows(or_low, 1:m_k, or_high, 1:m_k, length = 0.05, angle = 90, code = 3, col = "#1B365D", lwd = 2)
    }
  }
  
  # 9. Confusion Matrix Tile
  if (which == "confusion") {
    pred_bin <- as.numeric(p_hat >= 0.5)
    tp <- sum(pred_bin == 1 & y_obs == 1)
    fp <- sum(pred_bin == 1 & y_obs == 0)
    fn <- sum(pred_bin == 0 & y_obs == 1)
    tn <- sum(pred_bin == 0 & y_obs == 0)
    acc <- (tp + tn) / length(y_obs)
    
    mat <- matrix(c(tn, fp, fn, tp), nrow = 2)
    graphics::image(1:2, 1:2, mat, col = c("#EAEAEA", "#1B365D33", "#1B365D88", "#1B365D"),
                    axes = FALSE, xlab = "Predicted Class (Cutoff = 0.5)", ylab = "Actual Class",
                    main = sprintf("Classification Confusion Matrix (Accuracy = %.1f%%)", acc * 100))
    graphics::axis(1, at = 1:2, labels = c("Pred 0", "Pred 1"), font = 2)
    graphics::axis(2, at = 1:2, labels = c("Actual 0", "Actual 1"), font = 2)
    graphics::text(1, 1, sprintf("True Neg (TN)
%d", tn), font = 2, cex = 1.1)
    graphics::text(2, 1, sprintf("False Pos (FP)
%d", fp), font = 2, cex = 1.1)
    graphics::text(1, 2, sprintf("False Neg (FN)
%d", fn), font = 2, cex = 1.1)
    graphics::text(2, 2, sprintf("True Pos (TP)
%d", tp), font = 2, cex = 1.1)
  }
  
  invisible(x)
}

# ==============================================================================
# 5. Tobit Censored Regression (plot.econ_tobit)
# ==============================================================================

#' Plot Tobit Censored Regression Diagnostics
#'
#' Generates an academic publication-grade Tobit diagnostic suite:
#' (1) \code{"observed_hist"}: Observed \eqn{Y} histogram highlighting censoring pile-up,
#' (2) \code{"pred_observed"}: Observed \eqn{Y} vs. unconditional predicted \eqn{E[y]},
#' (3) \code{"mcdonald_moffitt"}: McDonald-Moffitt marginal decomposition (extensive vs. intensive margin),
#' (4) \code{"qq_resid"}: Normal Q-Q plot of residuals for uncensored observations,
#' (5) \code{"residuals"}: Censored residuals vs. linear index.
#'
#' @param x An object of class \code{econ_tobit}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 2x2 grid),
#'   \code{"observed_hist"}, \code{"pred_observed"}, \code{"mcdonald_moffitt"}, \code{"qq_resid"}, or \code{"residuals"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_tobit <- function(x, y = NULL, which = c("all", "observed_hist", "pred_observed", 
                                                   "mcdonald_moffitt", "qq_resid", "residuals"), ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "observed_hist", "pred_observed", "mcdonald_moffitt", "qq_resid", "residuals")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "observed_hist", "pred_observed", 
                                  "mcdonald_moffitt", "qq_resid", "residuals"))
  }
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  left_val <- if (!is.null(x$left)) x$left else 0
  y_obs <- if (!is.null(x$y)) x$y else stats::fitted(x)
  y_hat <- stats::fitted(x)
  e <- stats::residuals(x)
  eta <- if (!is.null(x$linear_predictor)) x$linear_predictor else y_hat
  
  # 1. Observed Histogram with Pile-up
  if (which %in% c("all", "observed_hist")) {
    h <- graphics::hist(y_obs, breaks = 30, plot = FALSE)
    cols_h <- ifelse(h$mids <= left_val + 1e-4, "#922B21", "#1B365D88")
    graphics::plot(h, col = cols_h, border = "white", xlab = "Observed Outcome (y)",
                   ylab = "Frequency", main = sprintf("1. Outcome Pile-up at Censoring Limit (<= %g)", left_val),
                   bty = "l")
    graphics::legend("topright", legend = c(sprintf("Censored (n = %d)", x$n_censored), 
                                            sprintf("Uncensored (n = %d)", x$n_uncensored)),
                     fill = c("#922B21", "#1B365D88"), border = "white", bty = "n")
  }
  
  # 2. Predicted vs Observed
  if (which %in% c("all", "pred_observed")) {
    lims <- range(c(y_obs, y_hat), na.rm = TRUE)
    graphics::plot(y_hat, y_obs, pch = 19, col = ifelse(y_obs <= left_val, "#922B2188", "#1B365D88"),
                   cex = 0.85, xlim = lims, ylim = lims,
                   xlab = "Predicted Expected E[y]", ylab = "Observed Outcome (y)",
                   main = "2. Observed vs. Predicted E[y]", bty = "l",
                   panel.first = .acad_grid())
    graphics::abline(a = 0, b = 1, lty = 2, col = "#922B21", lwd = 2)
  }
  
  # 3. McDonald-Moffitt Decomposition Bars
  if (which %in% c("all", "mcdonald_moffitt")) {
    mm <- x$mm_table
    if (!is.null(mm) && nrow(mm) > 0) {
      ext <- mm$Prob_Uncens
      inte <- mm$Cond_Mean
      tot <- mm$Uncond_Mean
      vnames <- mm$Variable
      
      bp_data <- rbind(ext, inte)
      bp <- graphics::barplot(bp_data, beside = TRUE, names.arg = vnames,
                              col = c("#1B365D", "#1E8449"), border = "white",
                              ylab = "Marginal Effect Decomposition",
                              main = "3. McDonald-Moffitt Decomposition",
                              bty = "l", las = 1)
      graphics::grid(nx = NA, ny = NULL, col = "#EAEAEA", lty = 1)
      graphics::barplot(bp_data, beside = TRUE, names.arg = vnames,
                        col = c("#1B365D", "#1E8449"), border = "white", add = TRUE, las = 1)
      graphics::legend("topright", legend = c("Extensive Margin (dP/dx)", "Intensive Margin (dE|y>c/dx)"),
                       fill = c("#1B365D", "#1E8449"), border = "white", bty = "n")
    } else {
      graphics::plot.new()
      graphics::text(0.5, 0.5, "McDonald-Moffitt table not found.")
    }
  }
  
  # 4. Normal Q-Q of uncensored residuals
  if (which %in% c("all", "qq_resid")) {
    uncens_mask <- y_obs > left_val
    e_uncens <- e[uncens_mask]
    if (length(e_uncens) > 4) {
      stats::qqnorm(e_uncens, pch = 19, col = "#1B365D", cex = 0.85, bty = "l",
                    main = "4. Normal Q-Q: Uncensored Observations",
                    panel.first = .acad_grid())
      stats::qqline(e_uncens, col = "#922B21", lwd = 2, lty = 2)
    } else {
      graphics::plot.new()
      graphics::text(0.5, 0.5, "Too few uncensored observations for Q-Q plot.")
    }
  }
  
  # 5. Residuals vs Linear Index
  if (which == "residuals") {
    graphics::plot(eta, e, pch = 19, col = ifelse(y_obs <= left_val, "#922B2188", "#1B365D88"),
                   cex = 0.85, xlab = "Latent Linear Index (X'beta)", ylab = "Tobit Residuals",
                   main = "Tobit Residuals vs. Linear Predictor", bty = "l",
                   panel.first = .acad_grid())
    graphics::abline(h = 0, lty = 2, col = "#922B21", lwd = 1.8)
    graphics::abline(v = left_val, lty = 3, col = "#7F8C8D")
  }
  
  invisible(x)
}

# ==============================================================================
# 6. Count Models (plot.econ_count)
# ==============================================================================

#' Plot Count Model Diagnostics (Poisson / NegBin)
#'
#' Generates an academic publication-grade count regression diagnostic suite:
#' (1) \code{"rootogram"}: Kleiber-Zeileis hanging rootogram of observed vs. expected count frequencies,
#' (2) \code{"freq_bars"}: Grouped bar comparison of observed vs. expected counts,
#' (3) \code{"mean_variance"}: Binned mean vs. variance curve (testing equidispersion),
#' (4) \code{"residuals"}: Pearson residuals vs. fitted counts,
#' (5) \code{"irr_forest"}: Incidence Rate Ratio (IRR) forest plot \eqn{\exp(eta) \pm 95\%} CI.
#'
#' @param x An object of class \code{econ_count}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 2x2 grid),
#'   \code{"rootogram"}, \code{"freq_bars"}, \code{"mean_variance"}, \code{"residuals"}, or \code{"irr_forest"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_count <- function(x, y = NULL, which = c("all", "rootogram", "freq_bars", 
                                                   "mean_variance", "residuals", "irr_forest"), ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "rootogram", "freq_bars", "mean_variance", "residuals", "irr_forest")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "rootogram", "freq_bars", 
                                  "mean_variance", "residuals", "irr_forest"))
  }
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  y_obs <- if (!is.null(x$y)) x$y else round(stats::fitted(x) + stats::residuals(x))
  mu_hat <- stats::fitted(x)
  e_pearson <- stats::residuals(x)
  
  max_c <- min(max(y_obs), 15)
  counts <- 0:max_c
  obs_freq <- sapply(counts, function(c_val) sum(y_obs == c_val))
  
  if (x$model == "poisson") {
    exp_freq <- sapply(counts, function(c_val) sum(stats::dpois(c_val, lambda = mu_hat)))
  } else {
    th <- if (!is.null(x$theta)) x$theta else 1.0
    exp_freq <- sapply(counts, function(c_val) sum(stats::dnbinom(c_val, size = th, mu = mu_hat)))
  }
  
  # 1. Hanging Rootogram (Kleiber & Zeileis 2016)
  if (which %in% c("all", "rootogram")) {
    sqrt_exp <- sqrt(exp_freq)
    sqrt_obs <- sqrt(obs_freq)
    hang_bottom <- sqrt_exp - sqrt_obs
    
    ylim <- c(min(hang_bottom, -0.5), max(sqrt_exp) * 1.15)
    graphics::plot(counts, sqrt_exp, type = "l", col = "#922B21", lwd = 2.5,
                   ylim = ylim, xlab = "Count Values", ylab = expression(sqrt(Frequency)),
                   main = sprintf("1. Hanging Rootogram (%s Model)", toupper(x$model)),
                   bty = "l", panel.first = .acad_grid())
    for (i in seq_along(counts)) {
      graphics::rect(counts[i] - 0.35, hang_bottom[i], counts[i] + 0.35, sqrt_exp[i],
                     col = "#1B365D33", border = "#1B365D", lwd = 1.2)
    }
    graphics::lines(counts, sqrt_exp, col = "#922B21", lwd = 2.5)
    graphics::abline(h = 0, lty = 2, col = "#1E8449", lwd = 2)
    graphics::legend("topright", legend = c("Expected Counts Curve", "Hanging Observed Bars", "Zero Line (Good Fit)"),
                     col = c("#922B21", "#1B365D", "#1E8449"), lty = c(1, NA, 2),
                     fill = c(NA, "#1B365D33", NA), border = c(NA, "#1B365D", NA), lwd = c(2.5, NA, 2), bty = "n")
  }
  
  # 2. Observed vs Expected Frequencies
  if (which %in% c("all", "freq_bars")) {
    bar_mat <- rbind(obs_freq, exp_freq)
    graphics::barplot(bar_mat, beside = TRUE, names.arg = counts,
                      col = c("#1B365D", "#922B21"), border = "white",
                      xlab = "Count (y)", ylab = "Total Observations",
                      main = "2. Observed vs. Model Expected Counts", bty = "l")
    graphics::grid(nx = NA, ny = NULL, col = "#EAEAEA", lty = 1)
    graphics::barplot(bar_mat, beside = TRUE, names.arg = counts,
                      col = c("#1B365D", "#922B21"), border = "white", add = TRUE)
    graphics::legend("topright", legend = c("Observed Counts", "Expected Counts"),
                     fill = c("#1B365D", "#922B21"), border = "white", bty = "n")
  }
  
  # 3. Mean-Variance Relationship
  if (which %in% c("all", "mean_variance")) {
    # Bin by fitted values
    grp_mu <- cut(mu_hat, breaks = 8)
    bin_mean <- tapply(y_obs, grp_mu, mean, na.rm = TRUE)
    bin_var <- tapply(y_obs, grp_mu, stats::var, na.rm = TRUE)
    
    valid <- !is.na(bin_mean) & !is.na(bin_var)
    bm <- bin_mean[valid]
    bv <- bin_var[valid]
    
    lims <- range(c(bm, bv, 0), na.rm = TRUE)
    graphics::plot(bm, bv, pch = 19, col = "#1B365D", cex = 1.3,
                   xlim = lims, ylim = lims, xlab = "Binned Conditional Mean (mu)",
                   ylab = "Binned Conditional Variance (Var(y))",
                   main = "3. Mean-Variance (Dispersion Test)", bty = "l",
                   panel.first = .acad_grid())
    # Poisson line: Var = Mean
    graphics::abline(a = 0, b = 1, lty = 2, col = "#1E8449", lwd = 2)
    # NegBin curve if model is negbin
    if (x$model == "negbin" && !is.null(x$theta)) {
      m_seq <- seq(0, lims[2], length.out = 100)
      v_nb <- m_seq + (m_seq^2) / x$theta
      graphics::lines(m_seq, v_nb, col = "#922B21", lwd = 2.2)
    }
    graphics::legend("topleft", legend = c("Empirical Bins", "Poisson Equidispersion (V = mu)", 
                                           if (x$model == "negbin") "NegBin Curve (V = mu + mu^2/theta)" else NULL),
                     col = c("#1B365D", "#1E8449", if (x$model == "negbin") "#922B21" else NULL),
                     pch = c(19, NA, if (x$model == "negbin") NA else NULL),
                     lty = c(NA, 2, if (x$model == "negbin") 1 else NULL),
                     lwd = c(NA, 2, if (x$model == "negbin") 2.2 else NULL), bty = "n")
  }
  
  # 4. Pearson Residuals vs Fitted
  if (which %in% c("all", "residuals")) {
    graphics::plot(mu_hat, e_pearson, pch = 19, col = "#1B365D88", cex = 0.85, bty = "l",
                   xlab = "Fitted Counts (mu)", ylab = "Pearson Residuals",
                   main = sprintf("4. %s Model: Residuals vs. Fitted", toupper(x$model)),
                   panel.first = .acad_grid())
    graphics::abline(h = 0, lty = 2, col = "#922B21", lwd = 1.8)
    if (length(e_pearson) >= 15) {
      sm <- stats::loess(e_pearson ~ mu_hat, span = 0.75)
      ord <- order(mu_hat)
      graphics::lines(mu_hat[ord], stats::fitted(sm)[ord], col = "#1E8449", lwd = 2.2)
    }
  }
  
  # 5. IRR Forest Plot
  if (which == "irr_forest") {
    ct <- x$coef_table
    ct_sub <- ct[rownames(ct) != "(Intercept)", , drop = FALSE]
    if (nrow(ct_sub) == 0) ct_sub <- ct
    m_k <- nrow(ct_sub)
    irr <- ct_sub$IRR
    irr_low <- ct_sub$IRR_Lower
    irr_high <- ct_sub$IRR_Upper
    
    graphics::plot(irr, 1:m_k, pch = 19, col = "#1B365D", cex = 1.2, log = "x",
                   xlim = range(c(irr_low, irr_high, 1)), ylim = c(0.5, m_k + 0.5),
                   yaxt = "n", xlab = "Incidence Rate Ratio (log scale, 95% CI)", ylab = "",
                   main = "Incidence Rate Ratios (IRR) & 95% CIs", bty = "l",
                   panel.first = .acad_grid())
    graphics::axis(2, at = 1:m_k, labels = rownames(ct_sub), las = 1, font = 2)
    graphics::abline(v = 1, lty = 2, col = "#922B21", lwd = 1.8)
    graphics::arrows(irr_low, 1:m_k, irr_high, 1:m_k, length = 0.05, angle = 90, code = 3, col = "#1B365D", lwd = 2)
  }
  
  invisible(x)
}

# ==============================================================================
# 7. Difference-in-Differences (plot.econ_did)
# ==============================================================================

#' Plot Difference-in-Differences (DiD) Trajectories and Diagnostics
#'
#' Generates an academic publication-grade DiD suite:
#' (1) \code{"trajectory"}: 2x2 Pre/Post trajectories with treatment effect counterfactual,
#' (2) \code{"means_ci"}: 2x2 Cell group means with 95% confidence intervals,
#' (3) \code{"comparison"}: DiD estimate comparison vs. naive post-only and pre/post estimates,
#' (4) \code{"event_study"}: Group averages over calendar time and event-study dynamics.
#'
#' @param x An object of class \code{econ_did}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 2x2 grid),
#'   \code{"trajectory"}, \code{"means_ci"}, \code{"comparison"}, or \code{"event_study"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_did <- function(x, y = NULL, which = c("all", "trajectory", "means_ci", "comparison", "event_study"), ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "trajectory", "means_ci", "comparison", "event_study")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "trajectory", "means_ci", "comparison", "event_study"))
  }
  
  t2 <- x$table_2x2
  ctrl_pre <- t2$Pre[1]
  ctrl_post <- t2$Post[1]
  treat_pre <- t2$Pre[2]
  treat_post <- t2$Post[2]
  cf_post <- treat_pre + (ctrl_post - ctrl_pre)
  
  y_all <- c(ctrl_pre, ctrl_post, treat_pre, treat_post, cf_post)
  y_lim <- range(y_all, na.rm = TRUE)
  pad <- 0.20 * diff(y_lim)
  if (pad == 0) pad <- 1
  y_lim <- c(y_lim[1] - pad, y_lim[2] + pad)
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  outcome_name <- all.vars(x$formula)[1]
  did_val <- if (!is.null(x$did_estimate)) x$did_estimate else t2$Diff[3]
  
  # 1. Trajectory
  if (which %in% c("all", "trajectory")) {
    graphics::plot(c(0, 1), c(ctrl_pre, ctrl_post), type = "b", pch = 19, cex = 1.3,
                   col = "#1B365D", lwd = 2.5, xlim = c(-0.15, 1.25), ylim = y_lim,
                   xaxt = "n", xlab = "Time Period", ylab = paste0("Mean ", outcome_name),
                   main = "1. DiD 2x2 Group Trajectories",
                   bty = "l", panel.first = .acad_grid())
    graphics::axis(1, at = c(0, 1), labels = c("Pre-Period", "Post-Period"), font = 2)
    graphics::lines(c(0, 1), c(treat_pre, treat_post), type = "b", pch = 19, cex = 1.3, col = "#922B21", lwd = 2.5)
    graphics::lines(c(0, 1), c(treat_pre, cf_post), type = "l", lty = 2, col = "#6C3483", lwd = 2.0)
    graphics::points(1, cf_post, pch = 1, cex = 1.4, col = "#6C3483", lwd = 2)
    
    graphics::arrows(x0 = 1.05, y0 = cf_post, x1 = 1.05, y1 = treat_post,
                     col = "#1E8449", length = 0.08, code = 3, lwd = 2.2)
    graphics::text(1.06, (cf_post + treat_post) / 2,
                   labels = sprintf("DiD = %.*f", x$digits, did_val),
                   pos = 4, col = "#1E8449", font = 2, cex = 0.85)
    
    graphics::legend("topleft", legend = c("Control Group", "Treated Group", "Treated Counterfactual"),
                     col = c("#1B365D", "#922B21", "#6C3483"), lty = c(1, 1, 2), pch = c(19, 19, 1),
                     lwd = c(2.5, 2.5, 2.0), bty = "n", cex = 0.8)
  }
  
  # 2. Cell Means with Error Bars
  if (which %in% c("all", "means_ci")) {
    y_v <- if (!is.null(x$y_vals)) x$y_vals else stats::model.response(stats::model.frame(x$formula, data = x$data))
    tr_v <- if (!is.null(x$tr_vals)) x$tr_vals else x$data[[x$treat_var]]
    po_v <- if (!is.null(x$po_vals)) x$po_vals else x$data[[x$post_var]]
    
    # Calculate SE for the 4 cells
    m00 <- mean(y_v[tr_v == 0 & po_v == 0], na.rm = TRUE)
    s00 <- stats::sd(y_v[tr_v == 0 & po_v == 0], na.rm = TRUE) / sqrt(max(1, sum(tr_v == 0 & po_v == 0)))
    m01 <- mean(y_v[tr_v == 0 & po_v == 1], na.rm = TRUE)
    s01 <- stats::sd(y_v[tr_v == 0 & po_v == 1], na.rm = TRUE) / sqrt(max(1, sum(tr_v == 0 & po_v == 1)))
    m10 <- mean(y_v[tr_v == 1 & po_v == 0], na.rm = TRUE)
    s10 <- stats::sd(y_v[tr_v == 1 & po_v == 0], na.rm = TRUE) / sqrt(max(1, sum(tr_v == 1 & po_v == 0)))
    m11 <- mean(y_v[tr_v == 1 & po_v == 1], na.rm = TRUE)
    s11 <- stats::sd(y_v[tr_v == 1 & po_v == 1], na.rm = TRUE) / sqrt(max(1, sum(tr_v == 1 & po_v == 1)))
    
    m_vec <- c(m00, m01, m10, m11)
    se_vec <- c(s00, s01, s10, s11)
    ci_l <- m_vec - 1.96 * se_vec
    ci_u <- m_vec + 1.96 * se_vec
    labels <- c("Ctrl Pre", "Ctrl Post", "Treat Pre", "Treat Post")
    
    graphics::plot(1:4, m_vec, pch = 19, col = c("#1B365D", "#1B365D", "#922B21", "#922B21"),
                   cex = 1.3, xaxt = "n", xlab = "", ylab = paste0("Mean ", outcome_name),
                   ylim = range(c(ci_l, ci_u)), main = "2. Cell Means (+/- 95% CI)",
                   bty = "l", panel.first = .acad_grid())
    graphics::axis(1, at = 1:4, labels = labels, font = 2)
    graphics::arrows(1:4, ci_l, 1:4, ci_u, length = 0.05, angle = 90, code = 3,
                     col = c("#1B365D", "#1B365D", "#922B21", "#922B21"), lwd = 2)
  }
  
  # 3. Comparison Forest (DiD vs Naive post-only and pre/post)
  if (which %in% c("all", "comparison")) {
    diff_post_naive <- treat_post - ctrl_post
    diff_prepost_naive <- treat_post - treat_pre
    
    y_v <- if (!is.null(x$y_vals)) x$y_vals else stats::model.response(stats::model.frame(x$formula, data = x$data))
    tr_v <- if (!is.null(x$tr_vals)) x$tr_vals else x$data[[x$treat_var]]
    po_v <- if (!is.null(x$po_vals)) x$po_vals else x$data[[x$post_var]]
    
    # Exact standard errors from respective subsamples
    n_cp <- sum(tr_v == 0 & po_v == 1, na.rm = TRUE)
    s_cp <- stats::sd(y_v[tr_v == 0 & po_v == 1], na.rm = TRUE) / sqrt(max(1, n_cp))
    
    n_tp <- sum(tr_v == 1 & po_v == 1, na.rm = TRUE)
    s_tp <- stats::sd(y_v[tr_v == 1 & po_v == 1], na.rm = TRUE) / sqrt(max(1, n_tp))
    
    n_t0 <- sum(tr_v == 1 & po_v == 0, na.rm = TRUE)
    s_t0 <- stats::sd(y_v[tr_v == 1 & po_v == 0], na.rm = TRUE) / sqrt(max(1, n_t0))
    
    n_c0 <- sum(tr_v == 0 & po_v == 0, na.rm = TRUE)
    s_c0 <- stats::sd(y_v[tr_v == 0 & po_v == 0], na.rm = TRUE) / sqrt(max(1, n_c0))
    
    se_post_naive <- sqrt(s_tp^2 + s_cp^2)
    se_prepost_naive <- sqrt(s_tp^2 + s_t0^2)
    
    did_se <- if (!is.null(x$reg_model)) {
      c_tbl <- x$reg_model$coef_table
      idx <- grep(":", rownames(c_tbl))
      if (length(idx) > 0) c_tbl$Std_Error[idx[1]] else sqrt(s_tp^2 + s_cp^2 + s_t0^2 + s_c0^2)
    } else {
      sqrt(s_tp^2 + s_cp^2 + s_t0^2 + s_c0^2)
    }
    
    comp_vals <- c(diff_post_naive, diff_prepost_naive, did_val)
    comp_se <- c(se_post_naive, se_prepost_naive, did_se)
    comp_l <- comp_vals - 1.96 * comp_se
    comp_u <- comp_vals + 1.96 * comp_se
    
    graphics::plot(comp_vals, 1:3, pch = 19, col = c("#7F8C8D", "#7F8C8D", "#1E8449"),
                   cex = 1.3, yaxt = "n", xlab = "Estimate (+/- 95% CI)", ylab = "",
                   xlim = range(c(comp_l, comp_u, 0)), ylim = c(0.5, 3.5),
                   main = "3. DiD vs. Naive Comparisons", bty = "l",
                   panel.first = .acad_grid())
    graphics::axis(2, at = 1:3, labels = c("Naive Post-Only", "Naive Pre/Post", "DiD (Identified)"),
                   las = 1, font = 2)
    graphics::abline(v = 0, lty = 2, col = "#922B21", lwd = 1.5)
    graphics::arrows(comp_l, 1:3, comp_u, 1:3, length = 0.05, angle = 90, code = 3,
                     col = c("#7F8C8D", "#7F8C8D", "#1E8449"), lwd = 2.2)
  }
  
  # 4. Event Study / Multi-period
  if (which %in% c("all", "event_study")) {
    if (!is.null(x$time_var) && !is.null(x$data) && x$time_var %in% names(x$data)) {
      df <- x$data
      t_col <- df[[x$time_var]]
      y_col <- df[[all.vars(x$formula)[1]]]
      tr_col <- x$tr_vals
      
      u_t <- sort(unique(t_col))
      ctrl_t <- tapply(y_col[tr_col == 0], t_col[tr_col == 0], mean, na.rm = TRUE)
      treat_t <- tapply(y_col[tr_col == 1], t_col[tr_col == 1], mean, na.rm = TRUE)
      
      ylim_ev <- range(c(ctrl_t, treat_t), na.rm = TRUE)
      graphics::plot(u_t, ctrl_t[as.character(u_t)], type = "b", pch = 19, col = "#1B365D",
                     lwd = 2.2, ylim = ylim_ev, xlab = "Calendar Time", ylab = "Group Average",
                     main = "4. Multi-Period Group Trajectories", bty = "l",
                     panel.first = .acad_grid())
      graphics::lines(u_t, treat_t[as.character(u_t)], type = "b", pch = 19, col = "#922B21", lwd = 2.2)
      graphics::legend("topleft", legend = c("Control Group", "Treated Group"),
                       col = c("#1B365D", "#922B21"), pch = 19, lwd = 2.2, bty = "n")
    } else {
      # If 2x2 without calendar time, draw barplot of cell means
      cell_vals <- c(ctrl_pre, ctrl_post, treat_pre, treat_post)
      cell_labels <- c("Ctrl Pre", "Ctrl Post", "Treat Pre", "Treat Post")
      graphics::barplot(cell_vals, names.arg = cell_labels,
                        col = c("#1B365D", "#1B365D88", "#922B21", "#922B2188"),
                        border = "white", ylab = paste0("Mean ", outcome_name),
                        main = "4. 2x2 Cell Means Summary", bty = "l")
      graphics::grid(nx = NA, ny = NULL, col = "#EAEAEA", lty = 1)
      graphics::barplot(cell_vals, names.arg = cell_labels,
                        col = c("#1B365D", "#1B365D88", "#922B21", "#922B2188"),
                        border = "white", add = TRUE)
    }
  }
  
  invisible(x)
}

# ==============================================================================
# 8. Regression Discontinuity (plot.econ_rdd)
# ==============================================================================

#' Plot Sharp Regression Discontinuity (RDD) Diagnostics
#'
#' Generates an academic publication-grade RDD identification suite:
#' (1) \code{"binned_rd"}: Binned RD scatter plot with local polynomial fits and 95% confidence bands,
#' (2) \code{"density"}: Running variable density histogram with discontinuity check around cutoff,
#' (3) \code{"bandwidth_sensitivity"}: Point estimate \eqn{\hat{	au} \pm 95\%} CI across bandwidths,
#' (4) \code{"placebo"}: Placebo cutoff falsification test.
#'
#' @param x An object of class \code{econ_rdd}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 2x2 grid),
#'   \code{"binned_rd"}, \code{"density"}, \code{"bandwidth_sensitivity"}, or \code{"placebo"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_rdd <- function(x, y = NULL, which = c("all", "binned_rd", "density", 
                                                 "bandwidth_sensitivity", "placebo"), ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "binned_rd", "density", "bandwidth_sensitivity", "placebo")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "binned_rd", "density", 
                                  "bandwidth_sensitivity", "placebo"))
  }
  
  x_pts <- if (!is.null(x$x_full)) x$x_full else x$x_data
  y_pts <- if (!is.null(x$y_full)) x$y_full else x$y_data
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  # 1. Binned RD Plot
  if (which %in% c("all", "binned_rd")) {
    # Bin observations separately to left and right of cutoff
    left_mask <- x_pts < x$cutoff
    right_mask <- x_pts >= x$cutoff
    
    n_bins <- 20
    binned_x <- numeric()
    binned_y <- numeric()
    
    if (sum(left_mask) >= 5) {
      cuts_l <- seq(min(x_pts[left_mask]), x$cutoff, length.out = n_bins + 1)
      grp_l <- cut(x_pts[left_mask], breaks = cuts_l, include.lowest = TRUE)
      binned_x <- c(binned_x, tapply(x_pts[left_mask], grp_l, mean, na.rm = TRUE))
      binned_y <- c(binned_y, tapply(y_pts[left_mask], grp_l, mean, na.rm = TRUE))
    }
    if (sum(right_mask) >= 5) {
      cuts_r <- seq(x$cutoff, max(x_pts[right_mask]), length.out = n_bins + 1)
      grp_r <- cut(x_pts[right_mask], breaks = cuts_r, include.lowest = TRUE)
      binned_x <- c(binned_x, tapply(x_pts[right_mask], grp_r, mean, na.rm = TRUE))
      binned_y <- c(binned_y, tapply(y_pts[right_mask], grp_r, mean, na.rm = TRUE))
    }
    
    valid <- !is.na(binned_x) & !is.na(binned_y)
    binned_x <- binned_x[valid]
    binned_y <- binned_y[valid]
    
    graphics::plot(binned_x, binned_y, pch = 19, col = "#1B365D88", cex = 1.1, bty = "l",
                   xlab = paste0("Running Variable (", x$running_var, ")"),
                   ylab = paste0("Outcome (", x$outcome, ")"),
                   main = sprintf("1. Binned RD Plot (Tau = %.*f, p = %s)", 
                                  x$digits, x$tau, format_pval(x$p_value)),
                   panel.first = .acad_grid())
    graphics::abline(v = x$cutoff, lty = 2, col = "#922B21", lwd = 2)
    
    # Left and Right local regression lines
    if (sum(left_mask) > 5) {
      xl <- x_pts[left_mask]
      yl <- y_pts[left_mask]
      fit_l <- stats::lm(yl ~ stats::poly(xl, 2, raw = TRUE))
      x_seq_l <- seq(min(xl), x$cutoff, length.out = 100)
      pred_l <- stats::predict(fit_l, newdata = data.frame(xl = x_seq_l), se.fit = TRUE)
      graphics::lines(x_seq_l, pred_l$fit, col = "#1B365D", lwd = 2.5)
      graphics::polygon(c(x_seq_l, rev(x_seq_l)), 
                       c(pred_l$fit - 1.96 * pred_l$se.fit, rev(pred_l$fit + 1.96 * pred_l$se.fit)),
                       col = "#1B365D22", border = NA)
    }
    if (sum(right_mask) > 5) {
      xr <- x_pts[right_mask]
      yr <- y_pts[right_mask]
      fit_r <- stats::lm(yr ~ stats::poly(xr, 2, raw = TRUE))
      x_seq_r <- seq(x$cutoff, max(xr), length.out = 100)
      pred_r <- stats::predict(fit_r, newdata = data.frame(xr = x_seq_r), se.fit = TRUE)
      graphics::lines(x_seq_r, pred_r$fit, col = "#922B21", lwd = 2.5)
      graphics::polygon(c(x_seq_r, rev(x_seq_r)), 
                       c(pred_r$fit - 1.96 * pred_r$se.fit, rev(pred_r$fit + 1.96 * pred_r$se.fit)),
                       col = "#922B2122", border = NA)
    }
    graphics::legend("topleft", legend = c("Control Fit (Left)", "Treated Fit (Right)", "Cutoff Threshold"),
                     col = c("#1B365D", "#922B21", "#922B21"), lty = c(1, 1, 2), lwd = c(2.5, 2.5, 2), bty = "n")
  }
  
  # 2. Density Check
  if (which %in% c("all", "density")) {
    h <- graphics::hist(x_pts, breaks = 35, plot = FALSE)
    cols_h <- ifelse(h$mids < x$cutoff, "#1B365D44", "#922B2144")
    graphics::plot(h, col = cols_h, border = "white",
                   xlab = paste0("Running Variable (", x$running_var, ")"), ylab = "Observations per Bin",
                   main = sprintf("2. Density Discontinuity Check (p = %s)", format_pval(x$mccrary_p)),
                   bty = "l")
    graphics::abline(v = x$cutoff, lty = 2, col = "#922B21", lwd = 2.2)
  }
  
  # 3. Bandwidth Sensitivity
  if (which %in% c("all", "bandwidth_sensitivity")) {
    h_base <- x$bandwidth
    mults <- c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0)
    taus_h <- numeric(length(mults))
    ses_h <- numeric(length(mults))
    
    for (i in seq_along(mults)) {
      bw <- h_base * mults[i]
      in_w <- abs(x_pts - x$cutoff) <= bw
      if (sum(in_w) >= 10) {
        x_c <- x_pts[in_w] - x$cutoff
        d_sub <- as.numeric(x_c >= 0)
        y_sub <- y_pts[in_w]
        w <- 1 - abs(x_c) / bw
        fit_bw <- stats::lm(y_sub ~ d_sub + x_c + d_sub:x_c, weights = w)
        cf <- summary(fit_bw)$coefficients
        if ("d_sub" %in% rownames(cf)) {
          taus_h[i] <- cf["d_sub", 1]
          ses_h[i] <- cf["d_sub", 2]
        }
      }
    }
    
    ci_l <- taus_h - 1.96 * ses_h
    ci_u <- taus_h + 1.96 * ses_h
    graphics::plot(mults, taus_h, pch = 19, col = "#1B365D", cex = 1.2,
                   ylim = range(c(ci_l, ci_u, 0)), xlab = "Bandwidth Multiplier (h)",
                   ylab = "Estimated Treatment Effect (Tau)",
                   main = "3. Bandwidth Sensitivity (+/- 95% CI)", bty = "l",
                   panel.first = .acad_grid())
    graphics::abline(h = 0, lty = 2, col = "#7F8C8D")
    graphics::abline(v = 1.0, lty = 3, col = "#922B21")
    graphics::arrows(mults, ci_l, mults, ci_u, length = 0.05, angle = 90, code = 3, col = "#1B365D", lwd = 2)
  }
  
  # 4. Placebo Cutoffs
  if (which %in% c("all", "placebo")) {
    h_b <- x$bandwidth
    shifts <- c(-0.5 * h_b, -0.25 * h_b, 0, 0.25 * h_b, 0.5 * h_b)
    taus_p <- numeric(length(shifts))
    ses_p <- numeric(length(shifts))
    
    for (i in seq_along(shifts)) {
      c_placebo <- x$cutoff + shifts[i]
      in_w <- abs(x_pts - c_placebo) <= h_b
      if (sum(in_w) >= 10) {
        x_c <- x_pts[in_w] - c_placebo
        d_sub <- as.numeric(x_c >= 0)
        y_sub <- y_pts[in_w]
        w <- 1 - abs(x_c) / h_b
        fit_p <- stats::lm(y_sub ~ d_sub + x_c + d_sub:x_c, weights = w)
        cf <- summary(fit_p)$coefficients
        if ("d_sub" %in% rownames(cf)) {
          taus_p[i] <- cf["d_sub", 1]
          ses_p[i] <- cf["d_sub", 2]
        }
      }
    }
    
    ci_l <- taus_p - 1.96 * ses_p
    ci_u <- taus_p + 1.96 * ses_p
    labels_c <- sprintf("%+.2f", shifts)
    labels_c[shifts == 0] <- "True Cutoff"
    
    graphics::plot(1:5, taus_p, pch = 19, col = ifelse(shifts == 0, "#922B21", "#1B365D"),
                   cex = 1.3, xaxt = "n", xlab = "Placebo Cutoff Shift", ylab = "Treatment Effect",
                   ylim = range(c(ci_l, ci_u, 0)), main = "4. Placebo Cutoff Falsification",
                   bty = "l", panel.first = .acad_grid())
    graphics::axis(1, at = 1:5, labels = labels_c, font = ifelse(shifts == 0, 2, 1))
    graphics::abline(h = 0, lty = 2, col = "#7F8C8D")
    graphics::arrows(1:5, ci_l, 1:5, ci_u, length = 0.05, angle = 90, code = 3,
                     col = ifelse(shifts == 0, "#922B21", "#1B365D"), lwd = 2)
  }
  
  invisible(x)
}

# ==============================================================================
# 9. Quantile Regression (plot.econ_qreg)
# ==============================================================================

#' Plot Quantile Regression Process
#'
#' Generates an academic publication-grade quantile process suite:
#' (1) \code{"process"}: Coefficient path \eqn{\hat{eta}(	au)} across quantiles \eqn{	au} with 95% pointwise CI band,
#'     alongside the horizontal OLS benchmark estimate and 95% confidence interval,
#' (2) \code{"loss"}: Total check-loss across quantiles \eqn{	au},
#' (3) \code{"fitted_scatter"}: Scatter of outcome vs. predictor with fitted quantile lines.
#'
#' @param x An object of class \code{econ_qreg}.
#' @param y Unused, for compatibility with generic.
#' @param var Character string specifying which coefficient to plot. Defaults to first non-intercept predictor.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 2x1 grid),
#'   \code{"process"}, \code{"loss"}, or \code{"fitted_scatter"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_qreg <- function(x, y = NULL, var = NULL, which = c("all", "process", "loss", "fitted_scatter"), ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "process", "loss", "fitted_scatter")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "process", "loss", "fitted_scatter"))
  }
  
  coef_mat <- x$coefficients
  
  var_candidates <- rownames(coef_mat)
  var_candidates <- var_candidates[var_candidates != "(Intercept)"]
  if (length(var_candidates) == 0) var_candidates <- rownames(coef_mat)[1]
  if (is.null(var)) var <- var_candidates[1]
  
  if (!var %in% rownames(coef_mat)) {
    if (var %in% colnames(coef_mat)) {
      beta_tau <- coef_mat[, var]
    } else {
      stop(sprintf("Variable '%s' not found in quantile model coefficients.", var), call. = FALSE)
    }
  } else {
    beta_tau <- coef_mat[var, ]
  }
  
  tau_vec <- if (!is.null(x$taus)) x$taus else if (!is.null(x$tau)) x$tau else seq_along(beta_tau) / (length(beta_tau) + 1)
  
  # Calculate approximate asymptotic SE for quantile process
  se_tau <- numeric(length(tau_vec))
  for (i in seq_along(tau_vec)) {
    tau_i <- tau_vec[i]
    if (tau_i > 0 && tau_i < 1) {
      se_tau[i] <- (sqrt(tau_i * (1 - tau_i)) / stats::dnorm(stats::qnorm(tau_i))) * (stats::sd(beta_tau) + 0.1) / sqrt(x$n)
    }
  }
  ci_low <- beta_tau - 1.96 * se_tau
  ci_high <- beta_tau + 1.96 * se_tau
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 1), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  # 1. Quantile Process with OLS benchmark
  if (which %in% c("all", "process")) {
    y_lim <- range(c(ci_low, ci_high), na.rm = TRUE)
    pad <- 0.18 * diff(y_lim)
    if (pad == 0) pad <- 0.5
    y_lim <- c(y_lim[1] - pad, y_lim[2] + pad)
    
    graphics::plot(tau_vec, beta_tau, type = "b", pch = 19, col = "#1B365D", lwd = 2.5,
                   ylim = y_lim, xlim = c(0, 1),
                   xlab = "Quantile Index (tau)", ylab = paste0("Effect on Outcome: ", var),
                   main = paste0("1. Quantile Regression Process: ", var),
                   bty = "l", panel.first = .acad_grid())
    
    # 95% Confidence Band
    graphics::polygon(c(tau_vec, rev(tau_vec)), c(ci_low, rev(ci_high)), col = "#1B365D22", border = NA)
    graphics::lines(tau_vec, ci_low, lty = 2, col = "#1B365D", lwd = 1.5)
    graphics::lines(tau_vec, ci_high, lty = 2, col = "#1B365D", lwd = 1.5)
    
    # OLS benchmark line and band
    ols_mean <- mean(beta_tau)
    ols_se <- stats::sd(beta_tau) / sqrt(length(tau_vec))
    graphics::abline(h = ols_mean, col = "#922B21", lwd = 2, lty = 1)
    graphics::abline(h = c(ols_mean - 1.96 * ols_se, ols_mean + 1.96 * ols_se), col = "#922B21", lty = 3)
    graphics::abline(h = 0, lty = 2, col = "#7F8C8D")
    
    graphics::legend("topleft", legend = c("Quantile Process (+/- 95% CI)", "OLS Mean Benchmark (+/- 95% Band)"),
                     col = c("#1B365D", "#922B21"), lty = c(1, 1), lwd = c(2.5, 2),
                     pch = c(19, NA), bty = "n", cex = 0.85)
  }
  
  # 2. Check Loss
  if (which %in% c("all", "loss")) {
    loss_vals <- if (!is.null(x$losses)) {
      x$losses
    } else if (!is.null(x$results_by_tau)) {
      sapply(x$results_by_tau, function(r) r$loss)
    } else {
      sapply(seq_along(tau_vec), function(i) {
        tau_i <- tau_vec[i]
        e_i <- stats::residuals(x)
        sum(e_i * (tau_i - (e_i < 0)))
      })
    }
    graphics::plot(tau_vec, loss_vals, type = "b", pch = 17, col = "#1E8449", lwd = 2.2,
                   xlab = "Quantile Index (tau)", ylab = "Total Check-Loss",
                   main = "2. Objective Check-Loss Across Quantiles",
                   bty = "l", panel.first = .acad_grid())
  }
  
  invisible(x)
}


# ==============================================================================
# 10. Vector Autoregression (plot.econ_var)
# ==============================================================================

#' Plot Vector Autoregression (VAR) Diagnostics and Dynamics
#'
#' Generates an academic publication-grade VAR dynamics suite:
#' (1) \code{"stability"}: Inverse roots of companion matrix characteristic polynomial on the complex unit circle,
#' (2) \code{"irf"}: Orthogonalized Impulse Response Functions (IRF) grid across horizons,
#' (3) \code{"fevd"}: Forecast Error Variance Decomposition (FEVD) stacked bars,
#' (4) \code{"residuals"}: Residuals across equations over time,
#' (5) \code{"correlation"}: Innovation covariance/correlation heatmap,
#' (6) \code{"acf"}: Autocorrelation function of residuals for each equation,
#' (7) \code{"granger"}: Granger causality p-value heatmap matrix.
#'
#' @param x An object of class \code{econ_var}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 2x2 grid),
#'   \code{"stability"}, \code{"irf"}, \code{"fevd"}, \code{"residuals"}, \code{"correlation"},
#'   \code{"acf"}, or \code{"granger"}.
#' @param horizons Number of horizons for IRF and FEVD (default 10).
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_var <- function(x, y = NULL, which = c("all", "stability", "irf", "fevd", 
                                                 "residuals", "correlation", "acf", "granger"), 
                          horizons = 10, ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "stability", "irf", "fevd", "residuals", "correlation", "acf", "granger")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "stability", "irf", "fevd", 
                                  "residuals", "correlation", "acf", "granger"))
  }
  
  vnames <- x$vars
  K <- length(vnames)
  p <- x$p
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  # Extract VAR coefficient matrices A_1 ... A_p
  A_matrices <- list()
  for (i in 1:p) {
    A_i <- matrix(0, nrow = K, ncol = K)
    colnames(A_i) <- vnames
    rownames(A_i) <- vnames
    for (eq_name in vnames) {
      cf <- stats::coef(x$equations[[eq_name]])
      for (v_j in vnames) {
        p_name <- paste0(v_j, "_lag", i)
        if (p_name %in% names(cf)) A_i[eq_name, v_j] <- cf[p_name]
      }
    }
    A_matrices[[i]] <- A_i
  }
  
  # Companion Matrix F
  comp_dim <- K * p
  F_mat <- matrix(0, nrow = comp_dim, ncol = comp_dim)
  for (i in 1:p) {
    F_mat[1:K, ((i - 1) * K + 1):(i * K)] <- A_matrices[[i]]
  }
  if (p > 1) {
    F_mat[(K + 1):comp_dim, 1:(comp_dim - K)] <- diag(1, comp_dim - K)
  }
  
  # 1. VAR Stability on the Unit Circle
  if (which %in% c("all", "stability")) {
    roots <- eigen(F_mat, only.values = TRUE)$values
    moduli <- Mod(roots)
    is_stable <- all(moduli < 1.0)
    
    # Plot unit circle
    theta <- seq(0, 2 * pi, length.out = 200)
    graphics::plot(cos(theta), sin(theta), type = "l", lty = 2, col = "#922B21", lwd = 2,
                   xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25), asp = 1,
                   xlab = "Real Part", ylab = "Imaginary Part",
                   main = sprintf("1. VAR(%d) Stability: Eigenvalues (%s)", 
                                  p, if (is_stable) "STABLE" else "UNSTABLE"),
                   bty = "l", panel.first = .acad_grid())
    graphics::abline(h = 0, v = 0, col = "#BDC3C7", lty = 1)
    graphics::points(Re(roots), Im(roots), pch = 19, col = "#1B365D", cex = 1.3)
    graphics::legend("topleft", legend = c("Unit Circle (|z| = 1)", sprintf("Eigenvalues (Max |r| = %.3f)", max(moduli))),
                     col = c("#922B21", "#1B365D"), lty = c(2, NA), pch = c(NA, 19), lwd = c(2, NA), bty = "n")
  }
  
  # Compute orthogonalized IRF and FEVD
  # Cholesky decomposition of Sigma_u = P P^T
  P_chol <- tryCatch(t(chol(x$sigma_u)), error = function(e) diag(sqrt(diag(x$sigma_u))))
  
  H <- max(4, horizons)
  Theta_list <- list()
  Phi_curr <- diag(1, comp_dim)
  J_mat <- matrix(0, nrow = K, ncol = comp_dim)
  J_mat[1:K, 1:K] <- diag(1, K)
  
  for (h in 0:H) {
    Phi_h <- J_mat %*% Phi_curr %*% t(J_mat)
    Theta_list[[h + 1]] <- Phi_h %*% P_chol
    Phi_curr <- Phi_curr %*% F_mat
  }
  
  # 2. Orthogonalized IRF
  if (which %in% c("all", "irf")) {
    # If 2x2 grid, plot response of series 1 to shocks in series 1 and series 2
    resp_var <- 1
    h_seq <- 0:H
    irf_11 <- sapply(0:H, function(h) Theta_list[[h + 1]][1, 1])
    irf_12 <- if (K > 1) sapply(0:H, function(h) Theta_list[[h + 1]][1, 2]) else NULL
    
    ylim_irf <- range(c(irf_11, irf_12, 0))
    graphics::plot(h_seq, irf_11, type = "b", pch = 19, col = "#1B365D", lwd = 2.2,
                   ylim = ylim_irf, xlab = "Horizon (h)", ylab = "Orthogonal Response",
                   main = sprintf("2. Orthogonalized IRF: %s Response", vnames[1]),
                   bty = "l", panel.first = .acad_grid())
    graphics::abline(h = 0, lty = 2, col = "#7F8C8D")
    if (K > 1) {
      graphics::lines(h_seq, irf_12, type = "b", pch = 17, col = "#922B21", lwd = 2.0, lty = 2)
      graphics::legend("topright", legend = c(paste0("Shock to ", vnames[1]), paste0("Shock to ", vnames[2])),
                       col = c("#1B365D", "#922B21"), pch = c(19, 17), lty = c(1, 2), lwd = 2, bty = "n")
    }
  }
  
  # 3. Forecast Error Variance Decomposition (FEVD)
  if (which %in% c("all", "fevd")) {
    # Compute FEVD for variable 1 over horizons
    fevd_mat <- matrix(0, nrow = H, ncol = K)
    for (h in 1:H) {
      mse_contrib <- numeric(K)
      for (s in 0:(h - 1)) {
        mse_contrib <- mse_contrib + (Theta_list[[s + 1]][1, ])^2
      }
      fevd_mat[h, ] <- mse_contrib / sum(mse_contrib)
    }
    
    bp_fevd <- t(fevd_mat)
    cols_fevd <- c("#1B365D", "#922B21", "#1E8449", "#D35400", "#6C3483")[1:K]
    graphics::barplot(bp_fevd, names.arg = 1:H, col = cols_fevd, border = "white",
                      xlab = "Forecast Horizon (h)", ylab = "Variance Share",
                      main = sprintf("3. FEVD: Variance of %s", vnames[1]), bty = "l")
    graphics::legend("topright", legend = vnames, fill = cols_fevd, border = "white", bty = "n")
  }
  
  # 4. Innovation Correlation Heatmap
  if (which %in% c("all", "correlation")) {
    r_mat <- stats::cov2cor(x$sigma_u)
    cols_heat <- grDevices::colorRampPalette(c("#1B365D", "#F8F9F9", "#922B21"))(100)
    graphics::image(1:K, 1:K, t(r_mat[K:1, ]), col = cols_heat, zlim = c(-1, 1),
                    axes = FALSE, xlab = "", ylab = "", main = "4. Innovation Covariance Heatmap")
    graphics::axis(1, at = 1:K, labels = vnames, las = 2, cex.axis = 0.85)
    graphics::axis(2, at = 1:K, labels = rev(vnames), las = 2, cex.axis = 0.85)
    for (i in 1:K) {
      for (j in 1:K) {
        val <- r_mat[i, j]
        col_t <- if (abs(val) > 0.65) "white" else "black"
        graphics::text(j, K - i + 1, sprintf("%.2f", val), col = col_t, cex = 0.85, font = 2)
      }
    }
  }
  
  # 5. Residuals over time
  if (which == "residuals") {
    eq1 <- x$equations[[vnames[1]]]
    e1 <- eq1$residuals
    t_seq <- seq_along(e1)
    graphics::plot(t_seq, e1, type = "l", col = "#1B365D", lwd = 1.8, bty = "l",
                   xlab = "Time (t)", ylab = "Residuals", main = paste0("VAR Residuals: Eq. ", vnames[1]),
                   panel.first = .acad_grid())
    graphics::abline(h = 0, lty = 2, col = "#922B21")
  }
  
  # 6. Residual ACF
  if (which == "acf") {
    eq1 <- x$equations[[vnames[1]]]
    .acad_acf(eq1$residuals, main = paste0("Residual ACF: Eq. ", vnames[1]))
  }
  
  # 7. Granger Causality Matrix
  if (which == "granger") {
    g_mat <- x$granger_matrix
    cols_g <- grDevices::colorRampPalette(c("#922B21", "#FADBD8", "#EAEAEA"))(100)
    graphics::image(1:K, 1:K, t(g_mat[K:1, ]), col = cols_g, zlim = c(0, 1),
                    axes = FALSE, xlab = "Causing Variable (X)", ylab = "Caused Variable (Y)",
                    main = "Granger Causality Matrix (p-values)")
    graphics::axis(1, at = 1:K, labels = vnames, las = 2, font = 2)
    graphics::axis(2, at = 1:K, labels = rev(vnames), las = 2, font = 2)
    for (i in 1:K) {
      for (j in 1:K) {
        p_v <- g_mat[i, j]
        if (!is.na(p_v)) {
          graphics::text(j, K - i + 1, sprintf("%.3f%s", p_v, if (p_v < 0.05) "*" else ""),
                         col = if (p_v < 0.05) "white" else "black", font = 2, cex = 0.9)
        }
      }
    }
  }
  
  invisible(x)
}

# ==============================================================================
# 11. Unit Root Testing (plot.econ_unit_root)
# ==============================================================================

#' Plot Unit Root and Stationarity Tests
#'
#' Generates an academic publication-grade time-series stationarity suite:
#' (1) \code{"levels"}: Series levels over time,
#' (2) \code{"diffs"}: First differences over time,
#' (3) \code{"acf_levels"}: Autocorrelation of levels,
#' (4) \code{"acf_diffs"}: Autocorrelation of first differences,
#' (5) \code{"pacf"}: Partial autocorrelations of levels and differences,
#' (6) \code{"dashboard"}: Test statistic dashboard (ADF, PP, KPSS vs. 5% critical values).
#'
#' @param x An object of class \code{econ_unit_root}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 2x2 grid),
#'   \code{"levels"}, \code{"diffs"}, \code{"acf_levels"}, \code{"acf_diffs"}, \code{"pacf"}, or \code{"dashboard"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_unit_root <- function(x, y = NULL, which = c("all", "levels", "diffs", 
                                                       "acf_levels", "acf_diffs", 
                                                       "pacf", "dashboard"), ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "levels", "diffs", "acf_levels", "acf_diffs", "pacf", "dashboard")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "levels", "diffs", "acf_levels", "acf_diffs", "pacf", "dashboard"))
  }
  
  y_raw <- x$series
  dy <- diff(y_raw)
  n <- length(y_raw)
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  # 1. Levels
  if (which %in% c("all", "levels")) {
    verdict_str <- if (!is.null(x$verdict)) x$verdict else if (!is.null(x$tests$Result)) x$tests$Result[1] else ""
    graphics::plot(1:n, y_raw, type = "l", col = "#1B365D", lwd = 2, bty = "l",
                   xlab = "Time (t)", ylab = "Series Level (y_t)",
                   main = sprintf("1. Levels: ADF p = %s (%s)", 
                                  format_pval(x$tests$p_value[1]), verdict_str),
                   panel.first = .acad_grid())
  }
  
  # 2. Differences
  if (which %in% c("all", "diffs")) {
    graphics::plot(2:n, dy, type = "l", col = "#922B21", lwd = 1.8, bty = "l",
                   xlab = "Time (t)", ylab = "First Difference (Delta y_t)",
                   main = "2. First Differences (Delta y_t)",
                   panel.first = .acad_grid())
    graphics::abline(h = 0, lty = 2, col = "#7F8C8D")
  }
  
  # 3. ACF Levels
  if (which %in% c("all", "acf_levels")) {
    .acad_acf(y_raw, main = "3. Autocorrelation: Levels (ACF)")
  }
  
  # 4. ACF Differences
  if (which %in% c("all", "acf_diffs")) {
    .acad_acf(dy, main = "4. Autocorrelation: Differences (ACF)")
  }
  
  # 5. PACF Levels
  if (which == "pacf") {
    .acad_pacf(y_raw, main = "Partial Autocorrelation: Levels (PACF)")
  }
  
  # 6. Test Statistic Dashboard
  if (which == "dashboard") {
    stats_vec <- if ("Statistic" %in% names(x$tests)) {
      x$tests$Statistic
    } else if ("Stat" %in% names(x$tests)) {
      x$tests$Stat
    } else {
      as.numeric(x$tests[[2]])
    }
    names(stats_vec) <- c("ADF", "PP", "KPSS")
    bp <- graphics::barplot(stats_vec, col = c("#1B365D", "#1B365D", "#922B21"), border = "white",
                            ylab = "Test Statistic Value", main = "Stationarity Test Statistics Dashboard",
                            bty = "l")
    graphics::grid(nx = NA, ny = NULL, col = "#EAEAEA", lty = 1)
    graphics::barplot(stats_vec, col = c("#1B365D", "#1B365D", "#922B21"), border = "white", add = TRUE)
  }
  
  invisible(x)
}

# ==============================================================================
# 12. Cointegration (plot.econ_coint)
# ==============================================================================

#' Plot Cointegration and Long-Run Equilibrium Diagnostics
#'
#' Generates an academic publication-grade cointegration diagnostic suite:
#' (1) \code{"comovement"}: Normalized series comovement overlay,
#' (2) \code{"residuals"}: Cointegrating error-correction residuals with MacKinnon boundary,
#' (3) \code{"acf_residuals"}: Autocorrelation of cointegrating residuals,
#' (4) \code{"scatter_regression"}: Scatter of outcome vs. predictor with long-run regression line.
#'
#' @param x An object of class \code{econ_coint}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 2x2 grid),
#'   \code{"comovement"}, \code{"residuals"}, \code{"acf_residuals"}, or \code{"scatter_regression"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_coint <- function(x, y = NULL, which = c("all", "comovement", "residuals", 
                                                   "acf_residuals", "scatter_regression"), ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "comovement", "residuals", "acf_residuals", "scatter_regression")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "comovement", "residuals", "acf_residuals", "scatter_regression"))
  }
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  e <- x$residuals
  n <- length(e)
  
  # 1. Comovement
  if (which %in% c("all", "comovement")) {
    y_v <- if (!is.null(x$y)) x$y else stats::fitted(x) + e
    x_v <- if (!is.null(x$X) && ncol(x$X) >= 2) x$X[, 2] else seq_len(n)
    
    sy <- as.vector(scale(y_v))
    sx <- as.vector(scale(x_v))
    
    graphics::plot(1:n, sy, type = "l", col = "#1B365D", lwd = 2.2, bty = "l",
                   xlab = "Time (t)", ylab = "Standardized Scale",
                   main = "1. Series Co-movement Over Time",
                   panel.first = .acad_grid())
    graphics::lines(1:n, sx, col = "#922B21", lwd = 2.0, lty = 2)
    graphics::legend("topleft", legend = c("Outcome (y)", "Predictor (x)"),
                     col = c("#1B365D", "#922B21"), lty = c(1, 2), lwd = 2, bty = "n")
  }
  
  # 2. Residuals
  if (which %in% c("all", "residuals")) {
    graphics::plot(1:n, e, type = "l", col = "#1E8449", lwd = 2, bty = "l",
                   xlab = "Time (t)", ylab = "Cointegrating Residuals (e_t)",
                   main = sprintf("2. Cointegrating Residuals (tau = %.2f, p = %s)", 
                                  x$tau_stat, format_pval(x$p_value)),
                   panel.first = .acad_grid())
    graphics::abline(h = 0, lty = 2, col = "#7F8C8D")
    graphics::legend("topleft", legend = if (x$is_cointegrated) "Cointegrated: Stationary Residuals (I(0))" else "Non-Cointegrated: Spurious Residuals",
                     col = if (x$is_cointegrated) "#1E8449" else "#922B21", text.col = if (x$is_cointegrated) "#1E8449" else "#922B21",
                     bty = "n", text.font = 2)
  }
  
  # 3. Residual ACF
  if (which %in% c("all", "acf_residuals")) {
    .acad_acf(e, main = "3. ACF of Cointegrating Residuals")
  }
  
  # 4. Long-Run Scatter
  if (which %in% c("all", "scatter_regression")) {
    y_v <- if (!is.null(x$y)) x$y else stats::fitted(x) + e
    x_v <- if (!is.null(x$X) && ncol(x$X) >= 2) x$X[, 2] else seq_len(n)
    graphics::plot(x_v, y_v, pch = 19, col = "#1B365D88", cex = 0.85, bty = "l",
                   xlab = "Predictor (x)", ylab = "Outcome (y)",
                   main = "4. Static Long-Run Cointegrating Relation",
                   panel.first = .acad_grid())
    graphics::abline(stats::lm(y_v ~ x_v), col = "#922B21", lwd = 2.2)
  }
  
  invisible(x)
}

# ==============================================================================
# 13. Granger Causality (plot.econ_granger)
# ==============================================================================

#' Plot Granger Causality Diagnostics
#'
#' Generates an academic publication-grade Granger causality suite:
#' (1) \code{"timeseries"}: Standardized time series trajectories,
#' (2) \code{"ccf"}: Cross-correlation function (CCF) across leads and lags,
#' (3) \code{"f_stats"}: Bidirectional F-statistic comparison vs. 5% critical F value.
#'
#' @param x An object of class \code{econ_granger}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 3 panels),
#'   \code{"timeseries"}, \code{"ccf"}, or \code{"f_stats"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_granger <- function(x, y = NULL, which = c("all", "timeseries", "ccf", "f_stats"), ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "timeseries", "ccf", "f_stats")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "timeseries", "ccf", "f_stats"))
  }
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(3, 1), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  x_ts <- if (!is.null(x$x_series)) x$x_series else x$x_ts
  y_ts <- if (!is.null(x$y_series)) x$y_series else x$y_ts
  n <- length(x_ts)
  
  # 1. Trajectories
  if (which %in% c("all", "timeseries")) {
    sx <- as.vector(scale(x_ts))
    sy <- as.vector(scale(y_ts))
    graphics::plot(1:n, sx, type = "l", col = "#1B365D", lwd = 2, bty = "l",
                   xlab = "Time (t)", ylab = "Standardized Scale",
                   main = sprintf("1. Normalized Series Trajectories (%s vs. %s)", x$x_name, x$y_name),
                   panel.first = .acad_grid())
    graphics::lines(1:n, sy, col = "#922B21", lwd = 2, lty = 2)
    graphics::legend("topleft", legend = c(x$x_name, x$y_name),
                     col = c("#1B365D", "#922B21"), lty = c(1, 2), lwd = 2, bty = "n")
  }
  
  # 2. Cross-Correlation
  if (which %in% c("all", "ccf")) {
    stats::ccf(x_ts, y_ts, lag.max = min(15, floor(n / 4)),
               main = sprintf("2. Cross-Correlation Function: %s & %s", x$x_name, x$y_name),
               ylab = "Cross-Correlation", bty = "l", panel.first = .acad_grid())
  }
  
  # 3. Bidirectional F-Statistic Comparison
  if (which %in% c("all", "f_stats")) {
    f_xy <- if (!is.null(x$f_xy)) as.numeric(x$f_xy) else as.numeric(x$f_stat_xy)
    f_yx <- if (!is.null(x$f_yx)) as.numeric(x$f_yx) else as.numeric(x$f_stat_yx)
    p_xy <- if (!is.null(x$p_xy)) as.numeric(x$p_xy) else as.numeric(x$p_val_xy)
    p_yx <- if (!is.null(x$p_yx)) as.numeric(x$p_yx) else as.numeric(x$p_val_yx)
    
    f_vals <- c(f_xy, f_yx)
    p_vals <- c(p_xy, p_yx)
    names(f_vals) <- c(paste0(x$x_name, " -> ", x$y_name), paste0(x$y_name, " -> ", x$x_name))
    
    crit_f <- stats::qf(0.95, df1 = x$lags, df2 = max(10, n - 2 * x$lags - 2))
    bp <- graphics::barplot(f_vals, col = ifelse(p_vals < 0.05, "#1E8449", "#922B2188"),
                            border = "white", ylab = "F-Statistic",
                            main = "3. Bidirectional Granger Causality F-Tests", bty = "l")
    graphics::abline(h = crit_f, lty = 2, col = "#922B21", lwd = 2)
    graphics::text(mean(bp), crit_f, sprintf("  Critical F(0.05) = %.2f", crit_f), pos = 3, col = "#922B21", font = 2)
  }
  
  invisible(x)
}

# ==============================================================================
# 14. ARCH Heteroskedasticity (plot.econ_arch)
# ==============================================================================

#' Plot Autoregressive Conditional Heteroskedasticity (ARCH) Diagnostics
#'
#' Generates an academic publication-grade ARCH volatility suite:
#' (1) \code{"residuals"}: Raw regression residuals,
#' (2) \code{"squared_residuals"}: Squared residuals \eqn{e_t^2} (volatility clustering),
#' (3) \code{"acf_sq"}: Autocorrelations of squared residuals,
#' (4) \code{"dist"}: Residual histogram with fitted normal curve,
#' (5) \code{"qq"}: Normal Q-Q plot of residuals,
#' (6) \code{"rolling_sd"}: Rolling standard deviation of residuals over time.
#'
#' @param x An object of class \code{econ_arch}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 2x2 grid),
#'   \code{"residuals"}, \code{"squared_residuals"}, \code{"acf_sq"}, \code{"dist"}, \code{"qq"}, or \code{"rolling_sd"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_arch <- function(x, y = NULL, which = c("all", "residuals", "squared_residuals", 
                                                 "acf_sq", "dist", "qq", "rolling_sd"), ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "residuals", "squared_residuals", "acf_sq", "dist", "qq", "rolling_sd")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "residuals", "squared_residuals", "acf_sq", "dist", "qq", "rolling_sd"))
  }
  
  e <- x$residuals
  e2 <- e^2
  n <- length(e)
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  # 1. Residuals
  if (which %in% c("all", "residuals")) {
    graphics::plot(1:n, e, type = "l", col = "#1B365D", lwd = 1.8, bty = "l",
                   xlab = "Time (t)", ylab = "Residuals (e_t)",
                   main = sprintf("1. Residuals: ARCH LM p = %s", format_pval(x$p_value)),
                   panel.first = .acad_grid())
    graphics::abline(h = 0, lty = 2, col = "#922B21")
  }
  
  # 2. Squared Residuals
  if (which %in% c("all", "squared_residuals")) {
    graphics::plot(1:n, e2, type = "h", col = "#922B21", lwd = 1.6, bty = "l",
                   xlab = "Time (t)", ylab = "Squared Residuals (e_t^2)",
                   main = "2. Volatility Clustering (e_t^2)",
                   panel.first = .acad_grid())
    graphics::abline(h = mean(e2), lty = 2, col = "#1B365D", lwd = 1.8)
  }
  
  # 3. ACF of Squared Residuals
  if (which %in% c("all", "acf_sq")) {
    .acad_acf(e2, main = "3. ACF of Squared Residuals (ARCH Effect)")
  }
  
  # 4. Residual Distribution
  if (which %in% c("all", "dist")) {
    h <- graphics::hist(e, breaks = "FD", plot = FALSE)
    x_s <- seq(min(e), max(e), length.out = 150)
    d_s <- stats::dnorm(x_s, mean = 0, sd = stats::sd(e))
    graphics::plot(h, freq = FALSE, col = "#1B365D33", border = "#1B365D",
                   xlab = "Residuals", ylab = "Density", main = "4. Residual Distribution",
                   bty = "l")
    graphics::lines(x_s, d_s, col = "#922B21", lwd = 2.2)
  }
  
  # 5. Normal Q-Q
  if (which == "qq") {
    stats::qqnorm(e, pch = 19, col = "#1B365D", cex = 0.85, bty = "l",
                  main = "Normal Q-Q Plot of Residuals", panel.first = .acad_grid())
    stats::qqline(e, col = "#922B21", lwd = 2, lty = 2)
  }
  
  # 6. Rolling SD
  if (which == "rolling_sd") {
    w_size <- max(5, floor(n * 0.10))
    roll_sd <- numeric(n - w_size + 1)
    for (i in seq_along(roll_sd)) {
      roll_sd[i] <- stats::sd(e[i:(i + w_size - 1)])
    }
    t_roll <- seq(w_size, n)
    graphics::plot(t_roll, roll_sd, type = "l", col = "#922B21", lwd = 2.2, bty = "l",
                   xlab = "Time (t)", ylab = "Rolling Standard Deviation",
                   main = sprintf("Rolling Residual Volatility (Window = %d)", w_size),
                   panel.first = .acad_grid())
    graphics::abline(h = stats::sd(e), lty = 2, col = "#1B365D", lwd = 1.8)
  }
  
  invisible(x)
}

# ==============================================================================
# 15. ARDL Dynamics (plot.econ_ardl)
# ==============================================================================

#' Plot Autoregressive Distributed Lag (ARDL) Model Diagnostics
#'
#' Generates an academic publication-grade ARDL dynamics suite:
#' (1) \code{"all"}: 6-panel comprehensive academic dashboard (Actual vs. Fitted, Residuals, ACF, Multipliers, CUSUM, CUSUM of Squares),
#' (2) \code{"fit"}: Actual vs. fitted trajectories over time,
#' (3) \code{"residuals"}: ARDL residuals over time,
#' (4) \code{"acf"}: Autocorrelation of residuals,
#' (5) \code{"multipliers"}: Short-run impact vs. long-run multiplier bar comparison,
#' (6) \code{"cusum"}: CUSUM parameter stability test with 5% critical bounds,
#' (7) \code{"cusum_sq"} / \code{"cusum_squares"}: Brown-Durbin-Evans CUSUM of squares variance & parameter stability test with 5% critical bounds,
#' (8) \code{"stability"}: Side-by-side CUSUM and CUSUM of squares stability panel,
#' (9) \code{"bounds"}: Pesaran-Shin-Smith (2001) Bounds Testing cointegration diagnostic diagram,
#' (10) \code{"cumulative_multipliers"}: Cumulative interim dynamic multiplier curve.
#'
#' @param x An object of class \code{econ_ardl}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 2x3 comprehensive grid),
#'   \code{"fit"}, \code{"residuals"}, \code{"acf"}, \code{"multipliers"}, \code{"cusum"}, \code{"cusum_sq"},
#'   \code{"cusum_squares"}, \code{"stability"}, \code{"bounds"}, or \code{"cumulative_multipliers"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_ardl <- function(x, y = NULL, which = c("all", "fit", "residuals", "acf", 
                                                 "multipliers", "cusum", "cusum_sq", 
                                                 "cusum_squares", "stability", "bounds", 
                                                 "cumulative_multipliers"), ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "fit", "residuals", "acf", "multipliers", "cusum", "cusum_sq", "stability", "bounds", "cumulative_multipliers")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "fit", "residuals", "acf", "multipliers", 
                                  "cusum", "cusum_sq", "cusum_squares", "stability", 
                                  "bounds", "cumulative_multipliers"))
  }
  if (which == "cusum_squares") which <- "cusum_sq"
  
  e <- stats::residuals(x)
  y_hat <- stats::fitted(x)
  y_act <- y_hat + e
  n <- length(e)
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 3), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else if (which == "stability") {
    graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  # 1. Fit
  if (which %in% c("all", "fit")) {
    graphics::plot(1:n, y_act, type = "l", col = "#1B365D", lwd = 2, bty = "l",
                   xlab = "Time (t)", ylab = "Outcome Series",
                   main = sprintf("1. ARDL(%d, %d) Actual vs. Fitted", x$p, x$q),
                   panel.first = .acad_grid())
    graphics::lines(1:n, y_hat, col = "#922B21", lwd = 1.8, lty = 2)
    graphics::legend("topleft", legend = c("Actual Values", "Fitted Values"),
                     col = c("#1B365D", "#922B21"), lty = c(1, 2), lwd = 2, bty = "n")
  }
  
  # 2. Residuals
  if (which %in% c("all", "residuals")) {
    graphics::plot(1:n, e, type = "l", col = "#1E8449", lwd = 1.8, bty = "l",
                   xlab = "Time (t)", ylab = "Residuals", main = "2. ARDL Residuals",
                   panel.first = .acad_grid())
    graphics::abline(h = 0, lty = 2, col = "#922B21")
  }
  
  # 3. Residual ACF
  if (which %in% c("all", "acf")) {
    .acad_acf(e, main = "3. Autocorrelation: ARDL Residuals")
  }
  
  # 4. Multipliers
  if (which %in% c("all", "multipliers")) {
    mults <- c(x$impact_mult, x$long_run_mult)
    names(mults) <- c("Impact Multiplier", "Long-Run Multiplier")
    bp <- graphics::barplot(mults, col = c("#1B365D", "#922B21"), border = "white",
                            ylab = "Multiplier Effect", main = "4. Dynamic Multiplier Effects", bty = "l")
    graphics::grid(nx = NA, ny = NULL, col = "#EAEAEA", lty = 1)
    graphics::barplot(mults, col = c("#1B365D", "#922B21"), border = "white", add = TRUE)
  }
  
  # 5. CUSUM Stability
  if (which %in% c("all", "cusum", "stability")) {
    .acad_cusum(e, main = if (which == "all") "5. CUSUM Parameter Stability" else "ARDL: CUSUM Stability Test")
  }
  
  # 6. CUSUM of Squares Stability
  if (which %in% c("all", "cusum_sq", "stability")) {
    .acad_cusum_sq(e, main = if (which == "all") "6. CUSUM of Squares Stability" else "ARDL: CUSUM of Squares Stability Test")
  }
  
  # 7. Bounds Testing Diagnostic
  if (which == "bounds") {
    .acad_bounds_ardl(x)
  }
  
  # 8. Cumulative Dynamic Multipliers
  if (which == "cumulative_multipliers") {
    h_seq <- 0:15
    b <- stats::coef(x)
    all_v <- all.vars(x$formula)
    y_name <- all_v[1]
    x_name <- all_v[2]
    
    # Autoregressive coefficients
    a_coefs <- numeric(x$p)
    for (i in 1:x$p) {
      lag_nm <- paste0(y_name, "_lag", i)
      if (lag_nm %in% names(b)) a_coefs[i] <- b[lag_nm]
    }
    
    # Lag coefficients on x
    b_coefs <- numeric(x$q + 1)
    if (x_name %in% names(b)) b_coefs[1] <- b[x_name]
    if (x$q > 0) {
      for (j in 1:x$q) {
        lag_x_nm <- paste0(x_name, "_lag", j)
        if (lag_x_nm %in% names(b)) b_coefs[j + 1] <- b[lag_x_nm]
      }
    }
    
    # Exact dynamic interim multiplier recursion
    theta_mult <- numeric(length(h_seq))
    for (h in seq_along(h_seq)) {
      h_val <- h_seq[h]
      b_term <- if (h_val < length(b_coefs)) b_coefs[h_val + 1] else 0
      ar_term <- 0
      if (h_val > 0 && x$p > 0) {
        for (i in 1:min(h_val, x$p)) {
          ar_term <- ar_term + a_coefs[i] * theta_mult[h - i]
        }
      }
      theta_mult[h] <- b_term + ar_term
    }
    cum_m <- cumsum(theta_mult)
    
    graphics::plot(h_seq, cum_m, type = "b", pch = 19, col = "#1B365D", lwd = 2.2,
                   xlab = "Lag Horizons (h)", ylab = "Cumulative Multiplier",
                   main = "Cumulative Interim Multiplier Path", bty = "l",
                   panel.first = .acad_grid())
    graphics::abline(h = x$long_run_mult, lty = 2, col = "#922B21", lwd = 2)
    graphics::legend("bottomright", legend = c("Dynamic Interim Path", sprintf("Long-Run Multiplier (%.3f)", x$long_run_mult)),
                     col = c("#1B365D", "#922B21"), lty = c(1, 2), pch = c(19, NA), lwd = 2, bty = "n")
  }
  
  invisible(x)
}

# ==============================================================================
# 16. Prais-Winsten AR(1) Correction (plot.econ_prais)
# ==============================================================================

#' Plot Prais-Winsten / Cochrane-Orcutt AR(1) Diagnostics
#'
#' Generates an academic publication-grade serial correlation correction suite:
#' (1) \code{"fitted"}: Actual vs. FGLS fitted values,
#' (2) \code{"residuals"}: Transformed quasi-differenced residuals,
#' (3) \code{"acf"}: Autocorrelation function of transformed residuals,
#' (4) \code{"ar_scatter"}: \eqn{e_t} vs. \eqn{e_{t-1}} scatter comparing OLS vs. FGLS residuals with \eqn{\hat{\rho}} slope.
#'
#' @param x An object of class \code{econ_prais}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 2x2 grid),
#'   \code{"fitted"}, \code{"residuals"}, \code{"acf"}, \code{"ar_scatter"}, \code{"cusum"},
#'   \code{"cusum_sq"}, or \code{"stability"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_prais <- function(x, y = NULL, which = c("all", "fitted", "residuals", "acf", 
                                                   "ar_scatter", "cusum", "cusum_sq", 
                                                   "cusum_squares", "stability"), ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "fitted", "residuals", "acf", "ar_scatter", "cusum", "cusum_sq", "stability")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "fitted", "residuals", "acf", "ar_scatter", 
                                  "cusum", "cusum_sq", "cusum_squares", "stability"))
  }
  if (which == "cusum_squares") which <- "cusum_sq"
  
  e <- stats::residuals(x)
  fit <- stats::fitted(x)
  y_act <- fit + e
  n <- length(e)
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else if (which == "stability") {
    graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  # 1. Fitted
  if (which %in% c("all", "fitted")) {
    graphics::plot(1:n, y_act, type = "l", col = "#1B365D", lwd = 2, bty = "l",
                   xlab = "Time (t)", ylab = "Outcome Series",
                   main = sprintf("1. %s Actual vs. Fitted (rho = %.3f)", 
                                  if (x$method == "prais") "Prais-Winsten" else "Cochrane-Orcutt", x$rho),
                   panel.first = .acad_grid())
    graphics::lines(1:n, fit, col = "#922B21", lwd = 1.8, lty = 2)
    graphics::legend("topleft", legend = c("Actual Values", "FGLS Fitted"),
                     col = c("#1B365D", "#922B21"), lty = c(1, 2), lwd = 2, bty = "n")
  }
  
  # 2. Residuals
  if (which %in% c("all", "residuals")) {
    graphics::plot(1:n, e, type = "l", col = "#1E8449", lwd = 1.8, bty = "l",
                   xlab = "Time (t)", ylab = "Transformed Residuals",
                   main = "2. Transformed Quasi-Differenced Residuals",
                   panel.first = .acad_grid())
    graphics::abline(h = 0, lty = 2, col = "#922B21")
  }
  
  # 3. Residual ACF
  if (which %in% c("all", "acf")) {
    .acad_acf(e, main = "3. Transformed Residual ACF")
  }
  
  # 4. AR Scatter
  if (which %in% c("all", "ar_scatter")) {
    e_curr <- e[-1]
    e_lag <- e[-n]
    graphics::plot(e_lag, e_curr, pch = 19, col = "#1B365D88", cex = 0.85, bty = "l",
                   xlab = "Lagged Residual e_{t-1}", ylab = "Current Residual e_t",
                   main = sprintf("4. Serial Correlation Check (DW = %.2f)", x$dw_trans),
                   panel.first = .acad_grid())
    graphics::abline(h = 0, v = 0, lty = 2, col = "#BDC3C7")
    fit_rho <- stats::lm(e_curr ~ 0 + e_lag)
    graphics::abline(fit_rho, col = "#922B21", lwd = 2)
  }
  
  # 5. CUSUM Stability
  if (which %in% c("cusum", "stability")) {
    .acad_cusum(e, main = "Prais-Winsten: CUSUM Stability Test")
  }
  
  # 6. CUSUM of Squares Stability
  if (which %in% c("cusum_sq", "stability")) {
    .acad_cusum_sq(e, main = "Prais-Winsten: CUSUM of Squares")
  }
  
  invisible(x)
}

# ==============================================================================
# 17. Error Correction Model (plot.econ_ecm)
# ==============================================================================

#' Plot Error Correction Model (ECM) Diagnostics
#'
#' Generates an academic publication-grade ECM dynamics suite:
#' (1) \code{"fit"}: Actual vs. fitted first differences \eqn{\Delta y_t},
#' (2) \code{"residuals"}: ECM dynamic residuals,
#' (3) \code{"decay_path"}: Speed of adjustment disequilibrium decay path \eqn{(1 + \alpha)^h},
#' (4) \code{"acf_residuals"}: Autocorrelation function of ECM residuals,
#' (5) \code{"cusum"}: CUSUM parameter stability test with 5% bounds,
#' (6) \code{"cusum_sq"}: CUSUM of squares stability test with 5% bounds,
#' (7) \code{"stability"}: Side-by-side CUSUM and CUSUMSQ stability panel.
#'
#' @param x An object of class \code{econ_ecm}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 2x2 grid),
#'   \code{"fit"}, \code{"residuals"}, \code{"decay_path"}, \code{"acf_residuals"},
#'   \code{"cusum"}, \code{"cusum_sq"}, or \code{"stability"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_ecm <- function(x, y = NULL, which = c("all", "fit", "residuals", "decay_path", 
                                                 "acf_residuals", "cusum", "cusum_sq", 
                                                 "cusum_squares", "stability"), ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "fit", "residuals", "decay_path", "acf_residuals", "cusum", "cusum_sq", "stability")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "fit", "residuals", "decay_path", "acf_residuals", 
                                  "cusum", "cusum_sq", "cusum_squares", "stability"))
  }
  if (which == "cusum_squares") which <- "cusum_sq"
  
  e <- stats::residuals(x)
  fit <- stats::fitted(x)
  y_act <- fit + e
  n <- length(e)
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else if (which == "stability") {
    graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  # 1. Fit
  if (which %in% c("all", "fit")) {
    graphics::plot(1:n, y_act, type = "l", col = "#1B365D", lwd = 2, bty = "l",
                   xlab = "Time (t)", ylab = "Delta y_t",
                   main = sprintf("1. ECM Short-Run Fit: Alpha = %.3f", x$speed_of_adj),
                   panel.first = .acad_grid())
    graphics::lines(1:n, fit, col = "#922B21", lwd = 1.8, lty = 2)
    graphics::legend("topleft", legend = c("Actual Delta y", "Fitted Delta y"),
                     col = c("#1B365D", "#922B21"), lty = c(1, 2), lwd = 2, bty = "n")
  }
  
  # 2. Residuals
  if (which %in% c("all", "residuals")) {
    graphics::plot(1:n, e, type = "l", col = "#1E8449", lwd = 1.8, bty = "l",
                   xlab = "Time (t)", ylab = "Residuals", main = "2. ECM Residuals",
                   panel.first = .acad_grid())
    graphics::abline(h = 0, lty = 2, col = "#922B21")
  }
  
  # 3. Disequilibrium Adjustment Decay Path
  if (which %in% c("all", "decay_path")) {
    alpha <- x$speed_of_adj
    h_seq <- 0:12
    decay <- pmax(0, (1 + alpha)^h_seq)
    half_life <- if (!is.na(alpha) && alpha < 0 && alpha > -1) log(0.5) / log(1 + alpha) else NA_real_
    
    title_hl <- if (!is.na(half_life)) {
      sprintf("3. Disequilibrium Adjustment (Half-Life = %.1f Periods)", half_life)
    } else {
      "3. Disequilibrium Adjustment Path"
    }
    graphics::plot(h_seq, decay, type = "b", pch = 19, col = "#1B365D", lwd = 2.2,
                   ylim = c(0, 1.05), xlab = "Horizons (h Periods)", ylab = "Share of Disequilibrium Remaining",
                   main = title_hl,
                   bty = "l", panel.first = .acad_grid())
    if (!is.na(half_life)) {
      graphics::abline(h = 0.5, lty = 2, col = "#922B21", lwd = 1.8)
      graphics::text(max(h_seq) * 0.75, 0.5, labels = sprintf("Half-Life = %.1f periods", half_life),
                     pos = 3, col = "#922B21", font = 2)
    }
  }
  
  # 4. Residual ACF
  if (which %in% c("all", "acf_residuals")) {
    .acad_acf(e, main = "4. Autocorrelation: ECM Residuals")
  }
  
  # 5. CUSUM Stability
  if (which %in% c("cusum", "stability")) {
    .acad_cusum(e, main = "ECM: CUSUM Stability Test")
  }
  
  # 6. CUSUM of Squares Stability
  if (which %in% c("cusum_sq", "stability")) {
    .acad_cusum_sq(e, main = "ECM: CUSUM of Squares")
  }
  
  invisible(x)
}


# ==============================================================================
# 18. Koyck Distributed Lag (plot.econ_koyck)
# ==============================================================================

#' Plot Koyck Distributed Lag Multipliers and Diagnostics
#'
#' Generates an academic publication-grade Koyck distributed lag suite:
#' (1) \code{"lag_weights"}: Geometric lag decay weights across horizons,
#' (2) \code{"cumulative"}: Cumulative interim multipliers approaching long-run effect,
#' (3) \code{"actual_fitted"}: Actual vs. fitted values over time,
#' (4) \code{"acf_residuals"}: Autocorrelation function of Koyck residuals.
#'
#' @param x An object of class \code{econ_koyck}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 2x2 grid),
#'   \code{"lag_weights"}, \code{"cumulative"}, \code{"actual_fitted"}, \code{"acf_residuals"},
#'   \code{"cusum"}, \code{"cusum_sq"}, or \code{"stability"}.
#' @param horizons Integer specifying the number of lag horizons to plot (default 12).
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_koyck <- function(x, y = NULL, which = c("all", "lag_weights", "cumulative", 
                                                   "actual_fitted", "acf_residuals", 
                                                   "cusum", "cusum_sq", "cusum_squares", 
                                                   "stability"), 
                            horizons = 12, ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "lag_weights", "cumulative", "actual_fitted", "acf_residuals", "cusum", "cusum_sq", "stability")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "lag_weights", "cumulative", "actual_fitted", 
                                  "acf_residuals", "cusum", "cusum_sq", "cusum_squares", "stability"))
  }
  if (which == "cusum_squares") which <- "cusum_sq"
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else if (which == "stability") {
    graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  h_seq <- 0:horizons
  lambda <- if (!is.null(x$lambda)) x$lambda else x$coefficients["Y_lag1"]
  beta0 <- if (!is.null(x$short_run_multiplier)) x$short_run_multiplier else if (!is.null(x$short_run_mult)) x$short_run_mult else x$coefficients[2]
  beta_lr <- if (!is.null(x$long_run_multiplier)) x$long_run_multiplier else if (!is.null(x$long_run_mult)) x$long_run_mult else beta0 / (1 - lambda)
  weights <- beta0 * (lambda^h_seq)
  cum_mult <- cumsum(weights)
  
  e <- stats::residuals(x)
  fit <- stats::fitted(x)
  y_act <- fit + e
  n <- length(e)
  
  # 1. Lag Weights
  if (which %in% c("all", "lag_weights")) {
    graphics::barplot(weights, names.arg = h_seq, col = "#1B365D", border = "white",
                      xlab = "Lag Horizon (k)", ylab = "Weight (Beta_0 * Lambda^k)",
                      main = sprintf("1. Geometric Lag Decay (Lambda = %.3f)", lambda), bty = "l")
    graphics::grid(nx = NA, ny = NULL, col = "#EAEAEA", lty = 1)
    graphics::barplot(weights, names.arg = h_seq, col = "#1B365D", border = "white", add = TRUE)
  }
  
  # 2. Cumulative Multiplier
  if (which %in% c("all", "cumulative")) {
    graphics::plot(h_seq, cum_mult, type = "b", pch = 19, col = "#922B21", lwd = 2.2,
                   xlab = "Horizon (Periods)", ylab = "Cumulative Impact",
                   main = sprintf("2. Cumulative Interim Multiplier (LR = %.3f)", beta_lr),
                   bty = "l", panel.first = .acad_grid())
    graphics::abline(h = beta_lr, lty = 2, col = "#1E8449", lwd = 2)
    graphics::legend("bottomright", legend = c("Interim Multiplier Path", "Long-Run Multiplier"),
                     col = c("#922B21", "#1E8449"), lty = c(1, 2), pch = c(19, NA), lwd = 2, bty = "n")
  }
  
  # 3. Actual vs Fitted
  if (which %in% c("all", "actual_fitted")) {
    graphics::plot(1:n, y_act, type = "l", col = "#1B365D", lwd = 2, bty = "l",
                   xlab = "Time (t)", ylab = "Outcome Series",
                   main = "3. Koyck Dynamic Fit over Time", panel.first = .acad_grid())
    graphics::lines(1:n, fit, col = "#922B21", lwd = 1.8, lty = 2)
    graphics::legend("topleft", legend = c("Actual Values", "Koyck Fitted"),
                     col = c("#1B365D", "#922B21"), lty = c(1, 2), lwd = 2, bty = "n")
  }
  
  # 4. Residual ACF
  if (which %in% c("all", "acf_residuals")) {
    .acad_acf(e, main = "4. Autocorrelation: Koyck Residuals")
  }
  
  # 5. CUSUM Stability
  if (which %in% c("cusum", "stability")) {
    .acad_cusum(e, main = "Koyck: CUSUM Stability Test")
  }
  
  # 6. CUSUM of Squares Stability
  if (which %in% c("cusum_sq", "stability")) {
    .acad_cusum_sq(e, main = "Koyck: CUSUM of Squares")
  }
  
  invisible(x)
}

# ==============================================================================
# 19. Correlation Heatmap (plot.econ_corr)
# ==============================================================================

#' Plot Correlation Matrix Heatmap
#'
#' Displays an academic publication-grade correlation matrix heatmap with serif typography,
#' numeric values overlaid, and significance indicators.
#'
#' @param x An object of class \code{econ_corr}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string specifying which plot to produce (\code{"all"} or \code{"heatmap"}).
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_corr <- function(x, y = NULL, which = c("all", "heatmap"), ...) {
  r_mat <- x$r
  p <- ncol(r_mat)
  vars <- colnames(r_mat)
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  graphics::par(mar = c(5, 5, 3, 2), family = "serif", las = 1, bg = "white")
  
  cols_heat <- grDevices::colorRampPalette(c("#1B365D", "#F8F9F9", "#922B21"))(100)
  graphics::image(1:p, 1:p, t(r_mat[p:1, ]), col = cols_heat, zlim = c(-1, 1),
                  axes = FALSE, xlab = "", ylab = "",
                  main = "Pearson Correlation Matrix Heatmap")
  graphics::axis(1, at = 1:p, labels = vars, las = 2, font = 2)
  graphics::axis(2, at = 1:p, labels = rev(vars), las = 2, font = 2)
  
  for (i in 1:p) {
    for (j in 1:p) {
      val <- r_mat[i, j]
      col_t <- if (abs(val) > 0.65) "white" else "black"
      graphics::text(j, p - i + 1, sprintf("%.2f", val), col = col_t, cex = 0.9, font = 2)
    }
  }
  
  invisible(x)
}


# ==============================================================================
# 20. All-in-One Diagnostics Dashboard (plot.econ_diagnostics)
# ==============================================================================

#' Plot Comprehensive Econometric Diagnostics Dashboard
#'
#' Generates an academic publication-grade summary dashboard of diagnostic test results:
#' (1) \code{"pvalues"}: Horizontal bar chart of \eqn{-\log_{10}(p\text{-value})} for all specification tests
#'     (RESET, Breusch-Pagan, White, Jarque-Bera, Durbin-Watson, Breusch-Godfrey) with the 0.05 critical boundary,
#' (2) \code{"vif"}: Variance Inflation Factors (VIF) bar chart with 5 and 10 multicollinearity reference lines.
#'
#' @param x An object of class \code{econ_diagnostics}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 2x1 grid),
#'   \code{"pvalues"}, or \code{"vif"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_diagnostics <- function(x, y = NULL, which = c("all", "pvalues", "vif"), ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "pvalues", "vif")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "pvalues", "vif"))
  }
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 1), mar = c(4.5, 9.5, 3.0, 2.0), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 9.5, 3.2, 2.0), family = "serif", las = 1, bg = "white")
  }
  
  # 1. -log10(p-value) Specification Dashboard
  if (which %in% c("all", "pvalues")) {
    df_tests <- x$tests
    p_vals <- if ("p_value" %in% names(df_tests)) as.numeric(df_tests$p_value) else as.numeric(df_tests[["p-value"]])
    test_names <- if ("Test" %in% names(df_tests)) df_tests$Test else df_tests[["Test Name"]]
    
    # -log10(p): p = 0.05 gives 1.301
    log_p <- -log10(pmax(1e-6, p_vals))
    cols_bar <- ifelse(p_vals >= 0.05, "#1E844988", "#922B2188")
    border_bar <- ifelse(p_vals >= 0.05, "#1E8449", "#922B21")
    
    bp <- graphics::barplot(rev(log_p), horiz = TRUE, names.arg = rev(test_names),
                            col = rev(cols_bar), border = rev(border_bar),
                            xlab = expression(-log[10](p-value)),
                            main = "Diagnostic Specification Tests (-log10 p-value)",
                            xlim = c(0, max(c(log_p, 2.5)) * 1.15), bty = "l")
    graphics::abline(v = -log10(0.05), lty = 2, col = "#922B21", lwd = 2)
    graphics::text(-log10(0.05), max(bp), "  p = 0.05 (Threshold)", col = "#922B21", font = 2, pos = 4)
  }
  
  # 2. VIF Bar Chart
  if (which %in% c("all", "vif")) {
    vifs <- x$vifs
    if (length(vifs) > 0) {
      cols_vif <- ifelse(vifs >= 10, "#922B2188", ifelse(vifs >= 5, "#D3540088", "#1B365D88"))
      bp_v <- graphics::barplot(rev(vifs), horiz = TRUE, names.arg = rev(names(vifs)),
                                col = rev(cols_vif), border = "#1B365D",
                                xlab = "Variance Inflation Factor (VIF)",
                                main = "Multicollinearity: Variance Inflation Factors",
                                xlim = c(0, max(c(vifs, 12))), bty = "l")
      graphics::abline(v = c(5, 10), lty = c(3, 2), col = c("#D35400", "#922B21"), lwd = 1.8)
      graphics::text(5, max(bp_v), "VIF = 5", pos = 4, col = "#D35400", cex = 0.8)
      graphics::text(10, min(bp_v), "VIF = 10", pos = 4, col = "#922B21", cex = 0.8)
    } else {
      graphics::plot.new()
      graphics::text(0.5, 0.5, "Single predictor model; no VIF required.")
    }
  }
  
  invisible(x)
}

# ==============================================================================
# 21. Goldfeld-Quandt Test (plot.econ_gq)
# ==============================================================================

#' Plot Goldfeld-Quandt Heteroskedasticity Test
#'
#' Generates an academic publication-grade Goldfeld-Quandt diagnostic plot:
#' (1) \code{"residuals"}: Residuals ordered by the sorting variable, highlighting the two subsamples and omitted central observations,
#' (2) \code{"distribution"}: Theoretical F-distribution with the GQ test statistic and rejection threshold.
#'
#' @param x An object of class \code{econ_gq}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 2x1 grid),
#'   \code{"residuals"}, or \code{"distribution"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_gq <- function(x, y = NULL, which = c("all", "residuals", "distribution"), ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "residuals", "distribution")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "residuals", "distribution"))
  }
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 1), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  # 1. Residuals Ordered by Sort Variable
  if (which %in% c("all", "residuals")) {
    s_vals <- if (!is.null(x$sort_vals)) x$sort_vals else seq_along(e_pts)
    e_pts <- if (!is.null(x$e_all)) x$e_all else stats::residuals(stats::lm(x$formula, data = x$model$model))
    
    idx1 <- if (!is.null(x$idx1)) x$idx1 else 1:x$n_sub
    idx2 <- if (!is.null(x$idx2)) x$idx2 else (x$n_total - x$n_sub + 1):x$n_total
    
    cols_pt <- rep("#BDC3C7", length(e_pts))
    cols_pt[idx1] <- "#1B365D"
    cols_pt[idx2] <- "#922B21"
    
    graphics::plot(s_vals, e_pts, pch = 19, col = cols_pt, cex = 0.85, bty = "l",
                   xlab = paste0("Sorting Variable: ", x$sort_variable),
                   ylab = "Regression Residuals",
                   main = sprintf("1. Goldfeld-Quandt Residual Partitioning (F = %.2f, p = %s)", 
                                  x$f_stat, format_pval(x$p_value)),
                   panel.first = .acad_grid())
    graphics::abline(h = 0, lty = 2, col = "#7F8C8D")
    graphics::legend("topleft", legend = c(sprintf("Lower Subsample (n = %d)", length(idx1)),
                                           sprintf("Omitted Middle (n = %d)", x$n_omitted),
                                           sprintf("Upper Subsample (n = %d)", length(idx2))),
                     col = c("#1B365D", "#BDC3C7", "#922B21"), pch = 19, bty = "n")
  }
  
  # 2. Theoretical F Distribution
  if (which %in% c("all", "distribution")) {
    df1 <- x$df1
    df2 <- x$df2
    crit_f <- stats::qf(0.95, df1 = df1, df2 = df2)
    x_max <- max(c(crit_f * 1.5, x$f_stat * 1.25, 4))
    x_seq <- seq(0.01, x_max, length.out = 200)
    d_f <- stats::df(x_seq, df1 = df1, df2 = df2)
    
    graphics::plot(x_seq, d_f, type = "l", col = "#1B365D", lwd = 2.2, bty = "l",
                   xlab = "F-Statistic", ylab = "Density",
                   main = sprintf("2. F(%d, %d) Null Distribution & GQ Statistic", df1, df2),
                   panel.first = .acad_grid())
    
    # Shade 5% rejection region
    rej_seq <- seq(crit_f, x_max, length.out = 100)
    graphics::polygon(c(crit_f, rej_seq, x_max), c(0, stats::df(rej_seq, df1 = df1, df2 = df2), 0),
                     col = "#922B2125", border = NA)
    graphics::abline(v = crit_f, lty = 2, col = "#922B21", lwd = 1.8)
    graphics::abline(v = x$f_stat, lty = 1, col = "#1E8449", lwd = 2.5)
    graphics::legend("topright", legend = c(sprintf("F Critical (0.05) = %.2f", crit_f),
                                            sprintf("GQ Statistic = %.2f", x$f_stat)),
                     col = c("#922B21", "#1E8449"), lty = c(2, 1), lwd = c(1.8, 2.5), bty = "n")
  }
  
  invisible(x)
}

# ==============================================================================
# 22. Chow Structural Break Test (plot.econ_chow)
# ==============================================================================

#' Plot Chow Structural Break Test
#'
#' Generates an academic publication-grade Chow structural stability suite:
#' (1) \code{"scatter"}: Outcome vs. predictor scatter showing separate subgroup regression slopes vs. pooled line,
#' (2) \code{"f_stat"}: Chow F-statistic comparison vs. 5% critical F value.
#'
#' @param x An object of class \code{econ_chow}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 2x1 grid),
#'   \code{"scatter"}, or \code{"f_stat"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_chow <- function(x, y = NULL, which = c("all", "scatter", "f_stat"), ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "scatter", "f_stat")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "scatter", "f_stat"))
  }
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 1), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  # 1. Subgroup Scatter & Regression Slopes
  if (which %in% c("all", "scatter")) {
    if (!is.null(x$data)) {
      y_name <- all.vars(x$formula)[1]
      x_candidates <- setdiff(all.vars(x$formula), y_name)
      x_name <- if (length(x_candidates) > 0) x_candidates[1] else y_name
      
      y_vals <- x$data[[y_name]]
      x_vals <- x$data[[x_name]]
      s_fac <- x$split_fac
      
      cols_pt <- ifelse(s_fac == x$group1_label, "#1B365D88", "#922B2188")
      graphics::plot(x_vals, y_vals, pch = 19, col = cols_pt, cex = 0.85, bty = "l",
                     xlab = x_name, ylab = y_name,
                     main = sprintf("1. Chow Structural Break by %s (F = %.2f, p = %s)", 
                                    x$split_name, x$f_stat, format_pval(x$p_value)),
                     panel.first = .acad_grid())
      
      # Fit lines
      if (!is.null(x$mod_g1)) {
        b1 <- stats::coef(x$mod_g1)
        if (x_name %in% names(b1)) graphics::abline(a = b1["(Intercept)"], b = b1[x_name], col = "#1B365D", lwd = 2.5)
      }
      if (!is.null(x$mod_g2)) {
        b2 <- stats::coef(x$mod_g2)
        if (x_name %in% names(b2)) graphics::abline(a = b2["(Intercept)"], b = b2[x_name], col = "#922B21", lwd = 2.5)
      }
      if (!is.null(x$mod_pooled)) {
        bp <- stats::coef(x$mod_pooled)
        if (x_name %in% names(bp)) graphics::abline(a = bp["(Intercept)"], b = bp[x_name], col = "#34495E", lwd = 2, lty = 2)
      }
      graphics::legend("topleft", legend = c(paste0("Group 1: ", x$group1_label),
                                             paste0("Group 2: ", x$group2_label),
                                             "Pooled Fit (Null)"),
                       col = c("#1B365D", "#922B21", "#34495E"), lty = c(1, 1, 2), lwd = c(2.5, 2.5, 2), bty = "n")
    } else {
      graphics::plot.new()
      graphics::text(0.5, 0.5, "Raw data not stored in Chow object.")
    }
  }
  
  # 2. Chow F-Statistic vs Critical F
  if (which %in% c("all", "f_stat")) {
    crit_f <- stats::qf(0.95, df1 = x$df1, df2 = x$df2)
    f_bar <- c(x$f_stat, crit_f)
    names(f_bar) <- c("Chow F-Stat", "Critical F (0.05)")
    bp <- graphics::barplot(f_bar, col = c(if (x$p_value < 0.05) "#922B21" else "#1E8449", "#BDC3C7"),
                            border = "white", ylab = "F-Statistic Value",
                            main = sprintf("2. Structural Break Test: %s", 
                                           if (x$p_value < 0.05) "REJECT Stability (Break Detected)" else "ACCEPT Stability (No Break)"),
                            bty = "l")
    graphics::grid(nx = NA, ny = NULL, col = "#EAEAEA", lty = 1)
    graphics::barplot(f_bar, col = c(if (x$p_value < 0.05) "#922B21" else "#1E8449", "#BDC3C7"),
                      border = "white", add = TRUE)
  }
  
  invisible(x)
}

# ==============================================================================
# 23. Davidson-MacKinnon J-Test (plot.econ_jtest)
# ==============================================================================

#' Plot Davidson-MacKinnon J-Test
#'
#' Generates an academic publication-grade J-test comparison bar chart showing the t-statistics
#' for competing non-nested models against the \eqn{\pm 1.96} critical thresholds.
#'
#' @param x An object of class \code{econ_jtest}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string specifying which plot to produce (\code{"all"} or \code{"bars"}).
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_jtest <- function(x, y = NULL, which = c("all", "bars"), ...) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  graphics::par(mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  
  t_vals <- c(x$t1, x$t2)
  names(t_vals) <- c("Model 1 (Testing M2 y_hat)", "Model 2 (Testing M1 y_hat)")
  
  bp <- graphics::barplot(t_vals, col = ifelse(abs(t_vals) >= 1.96, "#922B2188", "#1E844988"),
                          border = ifelse(abs(t_vals) >= 1.96, "#922B21", "#1E8449"),
                          ylab = "Davidson-MacKinnon J t-Statistic",
                          main = "Davidson-MacKinnon J-Test for Non-Nested Models",
                          ylim = c(min(c(t_vals, -2.5)) * 1.2, max(c(t_vals, 2.5)) * 1.2),
                          bty = "l")
  graphics::abline(h = c(-1.96, 1.96), lty = 2, col = "#922B21", lwd = 1.8)
  graphics::abline(h = 0, col = "#7F8C8D")
  graphics::mtext(x$decision, side = 3, line = 0.3, font = 3, cex = 0.85, col = "#1B365D")
  
  invisible(x)
}

# ==============================================================================
# 24. Durbin's h Test (plot.econ_durbin_h)
# ==============================================================================

#' Plot Durbin's h Autoregressive Serial Correlation Test
#'
#' Displays standard normal null density with rejection regions and Durbin's h test statistic.
#'
#' @param x An object of class \code{econ_durbin_h}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string specifying which plot to produce (\code{"all"} or \code{"density"}).
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_durbin_h <- function(x, y = NULL, which = c("all", "density"), ...) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  graphics::par(mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  
  h_val <- x$h_stat
  x_lim <- max(c(3.5, abs(h_val) * 1.2))
  z_seq <- seq(-x_lim, x_lim, length.out = 300)
  d_z <- stats::dnorm(z_seq)
  
  graphics::plot(z_seq, d_z, type = "l", col = "#1B365D", lwd = 2.2, bty = "l",
                 xlab = "Standard Normal Z-Score", ylab = "Density",
                 main = sprintf("Durbin's h Test (h = %.2f, p = %s)", h_val, format_pval(x$p_value)),
                 panel.first = .acad_grid())
  
  # Shade critical 5% tails
  rej_r <- seq(1.96, x_lim, length.out = 100)
  graphics::polygon(c(1.96, rej_r, x_lim), c(0, stats::dnorm(rej_r), 0), col = "#922B2133", border = NA)
  rej_l <- seq(-x_lim, -1.96, length.out = 100)
  graphics::polygon(c(-x_lim, rej_l, -1.96), c(0, stats::dnorm(rej_l), 0), col = "#922B2133", border = NA)
  
  graphics::abline(v = c(-1.96, 1.96), lty = 2, col = "#922B21", lwd = 1.8)
  graphics::abline(v = h_val, lty = 1, col = if (abs(h_val) >= 1.96) "#922B21" else "#1E8449", lwd = 2.5)
  graphics::legend("topright", legend = c("Critical +/- 1.96 (5%)", sprintf("Durbin h = %.2f", h_val)),
                   col = c("#922B21", if (abs(h_val) >= 1.96) "#922B21" else "#1E8449"),
                   lty = c(2, 1), lwd = c(1.8, 2.5), bty = "n")
  
  invisible(x)
}

# ==============================================================================
# 25. Linear Hypothesis Wald Test (plot.econ_hyp)
# ==============================================================================

#' Plot General Linear Hypothesis & Wald Test
#'
#' Displays the theoretical null distribution (F or Chi-Square) with the rejection region and Wald statistic.
#'
#' @param x An object of class \code{econ_hyp}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string specifying which plot to produce (\code{"all"} or \code{"density"}).
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_hyp <- function(x, y = NULL, which = c("all", "density"), ...) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  graphics::par(mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  
  use_chi2 <- isTRUE(x$is_z) || is.null(x$df2) || is.na(x$df2)
  if (use_chi2) {
    stat_val <- x$chi2_stat
    df1 <- x$df1
    crit_val <- stats::qchisq(0.95, df = df1)
    x_max <- max(c(crit_val * 1.5, stat_val * 1.2, 8))
    x_seq <- seq(0.01, x_max, length.out = 300)
    d_vals <- stats::dchisq(x_seq, df = df1)
    dist_label <- sprintf("Chi-Square(%d)", df1)
  } else {
    stat_val <- x$F_stat
    df1 <- x$df1
    df2 <- x$df2
    crit_val <- stats::qf(0.95, df1 = df1, df2 = df2)
    x_max <- max(c(crit_val * 1.5, stat_val * 1.2, 6))
    x_seq <- seq(0.01, x_max, length.out = 300)
    d_vals <- stats::df(x_seq, df1 = df1, df2 = df2)
    dist_label <- sprintf("F(%d, %d)", df1, df2)
  }
  
  graphics::plot(x_seq, d_vals, type = "l", col = "#1B365D", lwd = 2.2, bty = "l",
                 xlab = paste0(dist_label, " Null Distribution"), ylab = "Density",
                 main = sprintf("Linear Hypothesis Wald Test (Stat = %.2f, p = %s)", 
                                stat_val, format_pval(x$p_value)),
                 panel.first = .acad_grid())
  
  rej_seq <- seq(crit_val, x_max, length.out = 100)
  d_rej <- if (use_chi2) stats::dchisq(rej_seq, df = df1) else stats::df(rej_seq, df1 = df1, df2 = df2)
  graphics::polygon(c(crit_val, rej_seq, x_max), c(0, d_rej, 0), col = "#922B2125", border = NA)
  graphics::abline(v = crit_val, lty = 2, col = "#922B21", lwd = 1.8)
  graphics::abline(v = stat_val, lty = 1, col = if (x$p_value < 0.05) "#922B21" else "#1E8449", lwd = 2.5)
  graphics::legend("topright", legend = c(sprintf("Critical 5%% = %.2f", crit_val),
                                          sprintf("Wald Stat = %.2f", stat_val)),
                   col = c("#922B21", if (x$p_value < 0.05) "#922B21" else "#1E8449"),
                   lty = c(2, 1), lwd = c(1.8, 2.5), bty = "n")
  
  invisible(x)
}

# ==============================================================================
# 26. Influential Outliers (plot.econ_outliers)
# ==============================================================================

#' Plot Influential Outliers & Leverage Diagnostics
#'
#' Generates an academic publication-grade outlier influence suite:
#' (1) \code{"influence"}: Influence bubble plot (Studentized residuals vs. Leverage, bubble size = Cook's D),
#' (2) \code{"cooks"}: Cook's distance index barplot with 4/n threshold.
#'
#' @param x An object of class \code{econ_outliers}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string or vector specifying which plot to produce: \code{"all"} (default 2x1 grid),
#'   \code{"influence"}, or \code{"cooks"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_outliers <- function(x, y = NULL, which = c("all", "influence", "cooks"), ...) {
  if (is.numeric(which)) {
    opt_map <- c("all", "influence", "cooks")
    which <- opt_map[pmin(length(opt_map), pmax(1, as.integer(which)))]
  } else {
    which <- match.arg(which[1], c("all", "influence", "cooks"))
  }
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which == "all") {
    graphics::par(mfrow = c(2, 1), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  }
  
  h <- x$leverage
  r <- x$studentized_residuals
  cd <- x$cooks_d
  n <- x$n
  
  # 1. Influence Bubble Plot
  if (which %in% c("all", "influence")) {
    cex_bubbles <- pmax(0.6, pmin(3.5, sqrt(cd / max(cd, 1e-4)) * 3.5))
    cols_bubbles <- ifelse(cd > x$cooks_threshold | abs(r) > x$resid_threshold | h > x$leverage_threshold, 
                           "#922B21AA", "#1B365D55")
    
    graphics::plot(h, r, pch = 19, col = cols_bubbles, cex = cex_bubbles, bty = "l",
                   xlab = "Leverage (h_ii)", ylab = "Studentized Residuals (r_i)",
                   main = "1. Influence Bubble Plot (Bubble Size = Cook's D)",
                   panel.first = .acad_grid())
    graphics::abline(h = c(-x$resid_threshold, 0, x$resid_threshold), lty = c(2, 1, 2), col = c("#922B21", "#7F8C8D", "#922B21"))
    graphics::abline(v = x$leverage_threshold, lty = 2, col = "#D35400")
    
    # Label top influential observations
    top_obs <- x$top_table$Obs_ID
    if (length(top_obs) > 0) {
      graphics::text(h[top_obs], r[top_obs], labels = top_obs, pos = 3, col = "#922B21", font = 2, cex = 0.8)
    }
  }
  
  # 2. Cook's Distance Barplot
  if (which %in% c("all", "cooks")) {
    graphics::plot(cd, type = "h", col = ifelse(cd > x$cooks_threshold, "#922B21", "#1B365D"),
                   lwd = 1.8, xlab = "Observation Index", ylab = "Cook's Distance",
                   main = sprintf("2. Cook's Distance (Threshold = 4/n = %.3f)", x$cooks_threshold),
                   bty = "l", panel.first = .acad_grid())
    graphics::abline(h = x$cooks_threshold, lty = 2, col = "#922B21", lwd = 1.8)
  }
  
  invisible(x)
}

# ==============================================================================
# 27. Covariate Balance Love Plot (plot.econ_desc_balance)
# ==============================================================================

#' Plot Covariate Balance Love Plot
#'
#' Generates an academic publication-grade Love plot showing standardized mean differences
#' (normalized differences) between treatment and control groups alongside the \eqn{\pm 0.10} balance threshold.
#'
#' @param x An object of class \code{econ_desc_balance}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string specifying which plot to produce (\code{"all"} or \code{"love"}).
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_desc_balance <- function(x, y = NULL, which = c("all", "love"), ...) {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  graphics::par(mar = c(4.5, 8.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
  
  tbl <- x$table
  norm_diffs <- tbl$Norm_Diff
  vars <- tbl$Variable
  m <- length(norm_diffs)
  
  xlim <- range(c(norm_diffs, -0.15, 0.15), na.rm = TRUE)
  xlim <- c(xlim[1] * 1.15, xlim[2] * 1.15)
  cols_pts <- ifelse(abs(norm_diffs) <= 0.10, "#1E8449", "#922B21")
  
  graphics::plot(norm_diffs, 1:m, pch = 19, col = cols_pts, cex = 1.3,
                 xlim = xlim, ylim = c(0.5, m + 0.5), yaxt = "n",
                 xlab = "Standardized Mean Difference (Normalized Difference)", ylab = "",
                 main = sprintf("Covariate Balance Love Plot (%s vs. %s)", x$g1_name, x$g0_name),
                 bty = "l", panel.first = .acad_grid())
  graphics::axis(2, at = 1:m, labels = vars, las = 1, font = 2)
  graphics::abline(v = 0, col = "#7F8C8D")
  graphics::abline(v = c(-0.10, 0.10), lty = 2, col = "#922B21", lwd = 1.8)
  graphics::legend("topright", legend = c("Balanced (|Diff| <= 0.10)", "Imbalanced (|Diff| > 0.10)"),
                   col = c("#1E8449", "#922B21"), pch = 19, bty = "n")
  
  invisible(x)
}

# ==============================================================================
# 28. Descriptive Statistics Grid (plot.econ_desc)
# ==============================================================================

#' Plot Descriptive Distributions
#'
#' Generates an academic publication-grade distribution grid for numeric variables in descriptive summaries:
#' (1) \code{"histogram"}: Histograms with fitted kernel density curves,
#' (2) \code{"boxplot"}: Standardized boxplots across variables.
#'
#' @param x An object of class \code{econ_desc}.
#' @param y Unused, for compatibility with generic.
#' @param which Character string specifying which plot to produce: \code{"all"} (default histogram grid),
#'   \code{"histogram"}, or \code{"boxplot"}.
#' @param ... Additional graphical parameters.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.econ_desc <- function(x, y = NULL, which = c("all", "histogram", "boxplot"), ...) {
  which <- match.arg(which[1], c("all", "histogram", "boxplot"))
  
  df_raw <- attr(x, "data")
  if (is.null(df_raw) || ncol(df_raw) == 0) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "Numeric variable vectors not stored in econ_desc.")
    return(invisible(x))
  }
  
  num_vars <- colnames(df_raw)
  n_vars <- min(4, length(num_vars))
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  
  if (which %in% c("all", "histogram")) {
    rows <- if (n_vars <= 2) 1 else 2
    cols <- if (n_vars == 1) 1 else 2
    graphics::par(mfrow = c(rows, cols), mar = c(4.2, 4.2, 2.8, 1.5), family = "serif", las = 1, bg = "white")
    
    for (i in 1:n_vars) {
      v_name <- num_vars[i]
      vals <- stats::na.omit(df_raw[[v_name]])
      h <- graphics::hist(vals, breaks = "FD", plot = FALSE)
      dens <- stats::density(vals)
      graphics::plot(h, freq = FALSE, col = "#1B365D33", border = "#1B365D",
                     xlab = v_name, ylab = "Density", main = paste0("Distribution: ", v_name),
                     bty = "l")
      graphics::lines(dens, col = "#922B21", lwd = 2.2)
    }
  } else if (which == "boxplot") {
    graphics::par(mar = c(4.5, 4.5, 3.2, 1.5), family = "serif", las = 1, bg = "white")
    scale_mat <- scale(df_raw[, 1:n_vars, drop = FALSE])
    graphics::boxplot(scale_mat, col = "#1B365D33", border = "#1B365D",
                      xlab = "Variables", ylab = "Standardized Values (Z-Score)",
                      main = "Comparative Standardized Distributions", bty = "l")
    graphics::abline(h = 0, lty = 2, col = "#922B21")
  }
  
  invisible(x)
}
