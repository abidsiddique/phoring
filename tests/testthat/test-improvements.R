test_that("Dynamic column width formatting avoids variable name truncation", {
  data(wage_data)
  # Create a variable with name longer than 16 characters
  wage_data$very_long_variable_name_experience <- wage_data$exper
  mod <- reg(wage ~ educ + very_long_variable_name_experience, data = wage_data)
  
  # Capture output
  out <- capture.output(print(mod))
  # Must contain full variable name, not truncated
  expect_true(any(grepl("very_long_variable_name_experience", out)))
  
  # Also test logit with long variable name
  wage_data$very_long_binary_female <- wage_data$female
  mod_l <- logit(very_long_binary_female ~ educ + exper, data = wage_data)
  out_l <- capture.output(print(mod_l))
  expect_true(any(grepl("very_long_binary_female", out_l)))
})

test_that("interpret() correctly decomposes models and computes turning points", {
  data(wage_data)
  # Quadratic model
  mod_quad <- reg(wage ~ educ + exper + I(exper^2), data = wage_data)
  interp_quad <- interpret(mod_quad)
  
  expect_s3_class(interp_quad, "econ_interpret")
  expect_true("exper" %in% names(interp_quad$turning_points))
  
  tp <- interp_quad$turning_points$exper
  expect_equal(tp$shape, "Inverted U-shaped (Maximum)")
  # Turning point = -b1 / (2 * b2)
  b1 <- mod_quad$coefficients["exper"]
  b2 <- mod_quad$coefficients["I(exper^2)"]
  expect_equal(tp$turning_point, as.numeric(-b1 / (2 * b2)))
  
  out <- capture.output(print(interp_quad))
  expect_true(any(grepl("NONLINEAR TURNING POINT ANALYSIS", out)))
  expect_true(any(grepl("Turning Point", out)))
  
  # Log-log model
  data(house_prices)
  mod_log <- reg(log(price) ~ log(rooms) + crime, data = house_prices)
  interp_log <- interpret(mod_log)
  expect_equal(interp_log$interpretations[["log(rooms)"]]$type, "Log-Log (Elasticity)")
  expect_equal(interp_log$interpretations[["crime"]]$type, "Log-Level (Semi-Elasticity)")
  
  out_log <- capture.output(print(interp_log))
  expect_true(any(grepl("Elasticity", out_log)))
  
  # Works on standard lm
  base_lm <- lm(wage ~ educ + exper, data = wage_data)
  expect_s3_class(interpret(base_lm), "econ_interpret")
})

test_that("outliers() calculates Cook's D, studentized residuals, and leverage", {
  data(wage_data)
  mod <- reg(wage ~ educ + exper, data = wage_data)
  res <- outliers(mod, top_n = 5)
  
  expect_s3_class(res, "econ_outliers")
  expect_equal(nrow(res$top_table), 5)
  expect_true(res$cooks_threshold > 0)
  expect_true(all(res$cooks_d >= 0))
  expect_true(all(res$leverage >= 0 & res$leverage <= 1))
  
  out <- capture.output(print(res))
  expect_true(any(grepl("OUTLIER & INFLUENTIAL OBSERVATION DIAGNOSTICS", out)))
  expect_true(any(grepl("Cook's D", out)))
  
  # Works on standard lm
  base_lm <- lm(wage ~ educ + exper, data = wage_data)
  expect_s3_class(outliers(base_lm), "econ_outliers")
})

test_that("chow_test() evaluates parameter stability across subsamples", {
  data(wage_data)
  ch <- chow_test(wage ~ educ + exper, split_by = "female", data = wage_data)
  
  expect_s3_class(ch, "econ_chow")
  expect_true(ch$f_stat > 0)
  expect_true(ch$p_value >= 0 && ch$p_value <= 1)
  expect_equal(ch$df1, 3) # intercept + educ + exper
  expect_equal(ch$n1 + ch$n2, nrow(wage_data))
  
  out <- capture.output(print(ch))
  expect_true(any(grepl("CHOW TEST FOR STRUCTURAL BREAK", out)))
  
  # Test invalid split
  expect_error(chow_test(wage ~ educ, split_by = "nonexistent_col", data = wage_data))
})

test_that("desc() handles Table 1 covariate balance when by is specified", {
  data(wage_data)
  bal <- desc(wage_data, vars = c("wage", "educ", "exper"), by = "female")
  
  expect_s3_class(bal, "econ_desc_balance")
  expect_equal(nrow(bal$table), 3)
  expect_true(all(c("Mean_G0", "Mean_G1", "Difference", "t_stat", "p_value", "Norm_Diff", "Balance") %in% names(bal$table)))
  
  out <- capture.output(print(bal))
  expect_true(any(grepl("COVARIATE BALANCE TABLE", out)))
  expect_true(any(grepl("Difference", out)))
  expect_true(any(grepl("Norm Diff", out)))
  
  # Standard desc still works
  d_std <- desc(wage_data, vars = c("wage", "educ"))
  expect_s3_class(d_std, "econ_desc")
})

test_that("word_table() supports custom labels", {
  data(wage_data)
  m1 <- reg(wage ~ educ, data = wage_data)
  m2 <- reg(wage ~ educ + exper, data = wage_data)
  
  # Test custom labels
  out <- capture.output(word_table(m1, m2, labels = c("Bivariate", "Multivariate")))
  expect_true(any(grepl("Bivariate", out)))
  expect_true(any(grepl("Multivariate", out)))
  
  # Length mismatch raises error
  expect_error(word_table(m1, m2, labels = c("Only One")))
})

test_that("econ_plot() creates 2x2 diagnostic plots without error", {
  data(wage_data)
  mod <- reg(wage ~ educ + exper, data = wage_data)
  
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  expect_silent(econ_plot(mod, which = "all"))
  expect_silent(econ_plot(mod, which = "resid_fitted"))
  expect_silent(econ_plot(mod, which = "qq"))
  expect_silent(econ_plot(mod, which = "scale_location"))
  expect_silent(econ_plot(mod, which = "cooks"))
  grDevices::dev.off()
  if (file.exists(tf)) unlink(tf)
})

test_that("rdd() estimates Sharp Regression Discontinuity causal treatment effects", {
  set.seed(123)
  n <- 300
  x <- runif(n, -2, 2)
  d <- as.numeric(x >= 0)
  y <- 3 + 2.5 * d + 1.2 * x + rnorm(n, sd = 0.3)
  df <- data.frame(y = y, score = x)
  
  res <- rdd(y ~ score, data = df, cutoff = 0)
  expect_s3_class(res, "econ_rdd")
  expect_true(isTRUE(as.numeric(abs(res$tau - 2.5)) < 0.5)) # Estimate is close to true effect of 2.5
  expect_true(isTRUE(res$se > 0))
  expect_true(isTRUE(res$p_value < 0.001))
  expect_equal(res$cutoff, 0)
  expect_true(res$n_effective <= n && res$n_effective > 20)
  
  out <- capture.output(print(res))
  expect_true(any(grepl("SHARP REGRESSION DISCONTINUITY", out)))
  expect_true(any(grepl("LOCAL AVERAGE TREATMENT EFFECT", out)))
  
  # Uniform kernel and custom bandwidth
  res_unif <- rdd(y ~ score, data = df, cutoff = 0, bandwidth = 1.0, kernel = "uniform")
  expect_equal(res_unif$kernel, "uniform")
  expect_equal(res_unif$bandwidth, 1.0)
})

test_that("qreg() estimates single and multiple quantile regressions", {
  data(wage_data)
  
  # Single median regression
  q_med <- qreg(wage ~ educ + exper, data = wage_data, tau = 0.5)
  expect_s3_class(q_med, "econ_qreg")
  expect_equal(length(q_med$taus), 1)
  expect_true("educ" %in% rownames(coef(q_med)))
  
  out_single <- capture.output(print(q_med))
  expect_true(any(grepl("QUANTILE REGRESSION RESULTS", out_single)))
  expect_true(any(grepl("Quantile \\(tau = 0.50\\)", out_single)))
  
  # Multiple quantiles
  q_multi <- qreg(wage ~ educ + exper, data = wage_data, tau = c(0.25, 0.5, 0.75))
  expect_s3_class(q_multi, "econ_qreg")
  expect_equal(ncol(coef(q_multi)), 3)
  
  out_multi <- capture.output(print(q_multi))
  expect_true(any(grepl("DISTRIBUTIONAL PROFILE", out_multi)))
  expect_true(any(grepl("tau = 0.25", out_multi)))
  expect_true(any(grepl("tau = 0.75", out_multi)))
  
  # Classical SEs
  q_class <- qreg(wage ~ educ + exper, data = wage_data, tau = 0.5, se = "classical")
  expect_true(q_class$results_by_tau[[1]]$coef_table["educ", "Std_Error"] > 0)
  
  # Test invalid tau
  expect_error(qreg(wage ~ educ, data = wage_data, tau = 1.5))
})
