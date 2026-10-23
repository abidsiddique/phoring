test_that("All econometric visual plot methods execute cleanly", {
  pdf(NULL)
  on.exit(grDevices::dev.off())
  
  # 1. econ_plot and plot.econ_reg
  data(wage_data, package = "phoring")
  m_reg <- reg(wage ~ educ + exper, data = wage_data)
  expect_silent(econ_plot(m_reg, which = "resid_fitted"))
  expect_silent(plot(m_reg, which = "qq"))
  expect_silent(plot(m_reg, which = "all"))
  
  # 2. plot.econ_panel
  data(crime_panel, package = "phoring")
  m_panel <- panel(crime_rate ~ prbarr + police, id = "county_id", time = "year", 
                   data = crime_panel, model = "fe")
  expect_silent(plot(m_panel))
  
  # 3. plot.econ_iv
  set.seed(123)
  z <- rnorm(100)
  x_end <- 0.8 * z + rnorm(100, 0, 0.5)
  y_iv <- 2 * x_end + rnorm(100)
  df_iv <- data.frame(y = y_iv, x = x_end, z = z)
  m_iv <- iv(y ~ x | z, data = df_iv)
  expect_silent(plot(m_iv))
  
  # 4. plot.econ_unit_root
  data(macro_data, package = "phoring")
  m_ur <- unit_root(macro_data$unemployment, type = "drift")
  expect_silent(plot(m_ur))
  
  # 5. plot.econ_coint
  m_coint <- coint(inflation ~ unemployment, data = macro_data)
  expect_silent(plot(m_coint))
  
  # 6. plot.econ_granger
  m_granger <- granger(inflation ~ unemployment, data = macro_data, lags = 2)
  expect_silent(plot(m_granger))
  
  # 7. plot.econ_arch
  m_arch <- arch_test(m_reg, lags = 2)
  expect_silent(plot(m_arch))
  
  # 8. plot.econ_ardl
  m_ardl <- ardl(inflation ~ unemployment, data = macro_data, p = 1, q = 1)
  expect_silent(plot(m_ardl))
  
  # 9. plot.econ_prais
  m_prais <- prais_winsten(inflation ~ unemployment, data = macro_data)
  expect_silent(plot(m_prais))
  
  # 10. plot.econ_ecm
  m_ecm <- ecm(inflation ~ unemployment, data = macro_data)
  expect_silent(plot(m_ecm))
  
  # 11. plot.econ_var
  m_var <- var_model(macro_data, vars = c("gdp_growth", "inflation"), p = 1)
  expect_silent(plot(m_var))
  
  # 12. plot.econ_koyck
  set.seed(123)
  x_k <- rnorm(100)
  y_k <- numeric(100)
  for(t in 2:100) y_k[t] <- 0.5 * y_k[t-1] + 1.2 * x_k[t] + rnorm(1, 0, 0.2)
  df_k <- data.frame(y = y_k, x = x_k)
  m_koyck <- koyck(y ~ x, data = df_k)
  expect_silent(plot(m_koyck, horizons = 8))
  
  # 13. plot.econ_logit
  set.seed(123)
  df_bin <- data.frame(y = rbinom(100, 1, 0.5), x1 = rnorm(100), x2 = rnorm(100))
  m_logit <- logit(y ~ x1 + x2, data = df_bin)
  expect_silent(plot(m_logit))
  
  # 14. plot.econ_tobit
  data(labor_supply, package = "phoring")
  m_tobit <- tobit(hours ~ educ + exper + age, data = labor_supply, left = 0)
  expect_silent(plot(m_tobit))
  
  # 15. plot.econ_count
  set.seed(123)
  df_cnt <- data.frame(y = rpois(100, lambda = 3), x = rnorm(100))
  m_cnt <- count_reg(y ~ x, data = df_cnt, model = "poisson")
  expect_silent(plot(m_cnt))
  
  # 16. plot.econ_did
  set.seed(123)
  df_did <- data.frame(
    y = c(rnorm(50, 10), rnorm(50, 11), rnorm(50, 12), rnorm(50, 16)),
    treat = rep(c(0, 1), each = 100),
    post = rep(rep(c(0, 1), each = 50), 2)
  )
  m_did <- did(y ~ treat * post, data = df_did)
  expect_silent(plot(m_did))
  
  # 17. plot.econ_rdd
  set.seed(123)
  x_rdd <- runif(200, -5, 5)
  y_rdd <- 2 + 1.5 * x_rdd + 3 * (x_rdd >= 0) + rnorm(200)
  df_rdd <- data.frame(y = y_rdd, score = x_rdd)
  m_rdd <- rdd(y ~ score, data = df_rdd, cutoff = 0)
  expect_silent(plot(m_rdd))
  
  # 18. plot.econ_qreg
  m_qreg <- qreg(wage ~ educ + exper, data = wage_data, tau = c(0.2, 0.5, 0.8))
  expect_silent(plot(m_qreg, var = "educ"))
  
  # 19. plot.econ_corr
  m_corr <- corr(wage_data[, c("wage", "educ", "exper", "tenure")])
  expect_silent(plot(m_corr))
  
  # 20. plot.econ_diagnostics
  m_diag <- diagnostics(m_reg)
  expect_silent(plot(m_diag))
  expect_silent(plot(m_diag, which = "vif"))
  
  # 21. plot.econ_gq
  m_gq <- gq_test(wage ~ educ + exper, data = wage_data, split_by = "educ")
  expect_silent(plot(m_gq))
  expect_silent(plot(m_gq, which = "distribution"))
  
  # 22. plot.econ_chow
  m_chow <- chow_test(wage ~ educ + exper, split_by = "female", data = wage_data)
  expect_silent(plot(m_chow))
  expect_silent(plot(m_chow, which = "f_stat"))
  
  # 23. plot.econ_jtest
  m_j <- j_test(m_reg, reg(wage ~ tenure, data = wage_data))
  expect_silent(plot(m_j))
  
  # 24. plot.econ_durbin_h
  m_dh <- durbin_h(m_reg, lag_var = "exper")
  expect_silent(plot(m_dh))
  
  # 25. plot.econ_hyp
  m_hyp <- test_hyp(m_reg, "educ = 0")
  expect_silent(plot(m_hyp))
  
  # 26. plot.econ_outliers
  m_out <- outliers(m_reg)
  expect_silent(plot(m_out))
  expect_silent(plot(m_out, which = "cooks"))
  
  # 27. plot.econ_desc_balance (Love plot)
  m_bal <- desc(wage_data[, c("wage", "educ", "female")], by = "female")
  expect_silent(plot(m_bal))
  
  # 28. plot.econ_desc
  m_desc <- desc(wage_data[, c("wage", "educ", "exper")])
  expect_silent(plot(m_desc))
  expect_silent(plot(m_desc, which = "boxplot"))
  
  # Test specialized which options across econometric models
  expect_silent(plot(m_reg, which = "leverage"))
  expect_silent(plot(m_reg, which = "added_var"))
  expect_silent(plot(m_reg, which = "component_resid"))
  expect_silent(plot(m_reg, which = "dfbetas"))
  expect_silent(plot(m_reg, which = "forest"))
  expect_silent(plot(m_reg, which = "vif"))
  expect_silent(plot(m_reg, which = "histogram"))
  expect_silent(plot(m_reg, which = "observed_fitted"))
  expect_silent(plot(m_reg, which = "cusum"))
  
  expect_silent(plot(m_panel, which = "between"))
  expect_silent(plot(m_panel, which = "within"))
  expect_silent(plot(m_panel, which = "spaghetti"))
  expect_silent(plot(m_panel, which = "fe_dist"))
  expect_silent(plot(m_panel, which = "ar_residual"))
  expect_silent(plot(m_panel, which = "time_effects"))
  
  expect_silent(plot(m_iv, which = "forest_compare"))
  expect_silent(plot(m_iv, which = "f_stat"))
  expect_silent(plot(m_iv, which = "reduced_form"))
  
  expect_silent(plot(m_logit, which = "calibration"))
  expect_silent(plot(m_logit, which = "thresholds"))
  expect_silent(plot(m_logit, which = "marginal"))
  expect_silent(plot(m_logit, which = "odds_ratio"))
  expect_silent(plot(m_logit, which = "confusion"))
  
  expect_silent(plot(m_tobit, which = "observed_hist"))
  expect_silent(plot(m_tobit, which = "pred_observed"))
  expect_silent(plot(m_tobit, which = "mcdonald_moffitt"))
  expect_silent(plot(m_tobit, which = "qq_resid"))
  
  expect_silent(plot(m_cnt, which = "rootogram"))
  expect_silent(plot(m_cnt, which = "freq_bars"))
  expect_silent(plot(m_cnt, which = "mean_variance"))
  expect_silent(plot(m_cnt, which = "irr_forest"))
  
  expect_silent(plot(m_did, which = "trajectory"))
  expect_silent(plot(m_did, which = "means_ci"))
  expect_silent(plot(m_did, which = "comparison"))
  
  expect_silent(plot(m_rdd, which = "binned_rd"))
  expect_silent(plot(m_rdd, which = "density"))
  expect_silent(plot(m_rdd, which = "bandwidth_sensitivity"))
  expect_silent(plot(m_rdd, which = "placebo"))
  
  expect_silent(plot(m_qreg, which = "process"))
  expect_silent(plot(m_qreg, which = "loss"))
  
  expect_silent(plot(m_var, which = "stability"))
  expect_silent(plot(m_var, which = "irf"))
  expect_silent(plot(m_var, which = "fevd"))
  expect_silent(plot(m_var, which = "residuals"))
  expect_silent(plot(m_var, which = "correlation"))
  expect_silent(plot(m_var, which = "acf"))
  expect_silent(plot(m_var, which = "granger"))
  
  expect_silent(plot(m_ur, which = "pacf"))
  expect_silent(plot(m_ur, which = "dashboard"))
  
  expect_silent(plot(m_coint, which = "acf_residuals"))
  expect_silent(plot(m_coint, which = "scatter_regression"))
  
  expect_silent(plot(m_granger, which = "f_stats"))
  
  expect_silent(plot(m_arch, which = "qq"))
  expect_silent(plot(m_arch, which = "rolling_sd"))
  
  expect_silent(plot(m_ardl, which = "cusum"))
  expect_silent(plot(m_ardl, which = "cusum_sq"))
  expect_silent(plot(m_ardl, which = "cusum_squares"))
  expect_silent(plot(m_ardl, which = "stability"))
  expect_silent(plot(m_ardl, which = "bounds"))
  expect_silent(plot(m_ardl, which = "cumulative_multipliers"))
  
  expect_silent(plot(m_reg, which = "cusum_sq"))
  expect_silent(plot(m_reg, which = "stability"))
  expect_silent(plot(m_reg, which = "influence"))
  
  expect_silent(plot(m_prais, which = "ar_scatter"))
  expect_silent(plot(m_prais, which = "cusum_sq"))
  expect_silent(plot(m_prais, which = "stability"))
  
  expect_silent(plot(m_ecm, which = "decay_path"))
  expect_silent(plot(m_ecm, which = "cusum_sq"))
  expect_silent(plot(m_ecm, which = "stability"))
  
  expect_silent(plot(m_koyck, which = "actual_fitted"))
  expect_silent(plot(m_koyck, which = "cusum_sq"))
  expect_silent(plot(m_koyck, which = "stability"))
})

