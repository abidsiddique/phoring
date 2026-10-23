test_that("gq_test() computes Goldfeld-Quandt heteroskedasticity test", {
  data(wage_data)
  res <- gq_test(wage ~ educ + exper, split_by = "educ", data = wage_data)
  
  expect_s3_class(res, "econ_gq")
  expect_true(res$f_stat > 0)
  expect_true(res$p_value >= 0 && res$p_value <= 1)
  expect_equal(res$sort_variable, "educ")
  expect_true(res$n_sub > 50)
  
  out <- capture.output(print(res))
  expect_true(any(grepl("GOLDFELD-QUANDT HETEROSCEDASTICITY TEST", out)))
  expect_true(any(grepl("Goldfeld-Quandt F", out)))
  
  # Error if sorting variable does not exist
  expect_error(gq_test(wage ~ educ, split_by = "nonexistent", data = wage_data))
})

test_that("j_test() evaluates non-nested competing model specifications", {
  data(wage_data)
  m1 <- reg(wage ~ educ + exper, data = wage_data)
  m2 <- reg(wage ~ tenure + mother_educ, data = wage_data)
  
  res <- j_test(m1, m2)
  expect_s3_class(res, "econ_jtest")
  expect_true(is.numeric(res$t1))
  expect_true(is.numeric(res$t2))
  expect_true(res$p1 >= 0 && res$p1 <= 1)
  expect_true(res$p2 >= 0 && res$p2 <= 1)
  expect_true(nchar(res$decision) > 10)
  
  out <- capture.output(print(res))
  expect_true(any(grepl("DAVIDSON-MacKINNON J-TEST", out)))
  expect_true(any(grepl("H0: Model 1 is correct", out)))
  expect_true(any(grepl("Conclusion", out)))
  
  # Also works when passing formula directly with data
  res_f <- j_test(wage ~ educ, wage ~ tenure, data = wage_data)
  expect_s3_class(res_f, "econ_jtest")
})

test_that("durbin_h() tests autocorrelation in dynamic autoregressive models", {
  data(macro_data)
  n <- nrow(macro_data)
  macro_data$inf_lag <- c(NA, macro_data$inflation[-n])
  m <- reg(inflation ~ inf_lag + unemployment, data = macro_data[-1, ])
  
  res <- durbin_h(m, lag_var = "inf_lag")
  expect_s3_class(res, "econ_durbin_h")
  expect_true(is.numeric(res$h_stat))
  expect_true(res$p_value >= 0 && res$p_value <= 1)
  expect_equal(res$lag_var, "inf_lag")
  
  out <- capture.output(print(res))
  expect_true(any(grepl("DURBIN'S h TEST", out)))
  expect_true(any(grepl("Residual Autocorrelation", out)))
})

test_that("koyck() estimates geometric distributed lag model and multipliers", {
  data(macro_data)
  res <- koyck(inflation ~ unemployment, data = macro_data)
  
  expect_s3_class(res, "econ_koyck")
  expect_true(is.numeric(res$short_run_multiplier))
  expect_true(is.numeric(res$lambda))
  expect_equal(res$outcome, "inflation")
  expect_equal(res$predictor, "unemployment")
  expect_equal(nrow(res$coef_table), 3)
  
  out <- capture.output(print(res))
  expect_true(any(grepl("KOYCK DISTRIBUTED LAG MODEL", out)))
  expect_true(any(grepl("Short-Run Multiplier", out)))
  expect_true(any(grepl("Rate of Decay", out)))
})

test_that("prais_winsten() and cochrane_orcutt() correct AR(1) serial correlation", {
  data(macro_data)
  
  # Prais-Winsten
  pw <- prais_winsten(inflation ~ fed_funds + unemployment, data = macro_data, method = "prais")
  expect_s3_class(pw, "econ_prais")
  expect_true(pw$converged)
  expect_true(pw$rho >= -1 && pw$rho <= 1)
  expect_true(pw$dw_transformed > pw$dw_initial) # DW should move closer towards 2
  
  out_pw <- capture.output(print(pw))
  expect_true(any(grepl("PRAIS-WINSTEN FGLS", out_pw)))
  expect_true(any(grepl("Target ~ 2.0", out_pw)))
  
  # Cochrane-Orcutt
  co <- prais_winsten(inflation ~ fed_funds + unemployment, data = macro_data, method = "cochrane")
  expect_s3_class(co, "econ_prais")
  expect_true(co$converged)
})

test_that("ecm() estimates Error Correction Model with speed of adjustment", {
  data(macro_data)
  res <- ecm(inflation ~ unemployment + fed_funds, data = macro_data, lags = 1)
  
  expect_s3_class(res, "econ_ecm")
  expect_true(is.numeric(res$speed_of_adjustment))
  expect_true(res$gamma_p >= 0 && res$gamma_p <= 1)
  expect_true("ECT_Lag" %in% rownames(res$coef_table))
  expect_equal(res$outcome, "inflation")
  
  out <- capture.output(print(res))
  expect_true(any(grepl("ERROR CORRECTION MODEL \\(ECM\\)", out)))
  expect_true(any(grepl("Speed of Adjustment", out)))
})

test_that("var_model() estimates Vector Autoregression and Granger causality matrix", {
  data(macro_data)
  v <- var_model(macro_data, vars = c("inflation", "unemployment", "fed_funds"), p = 1)
  
  expect_s3_class(v, "econ_var")
  expect_equal(length(v$vars), 3)
  expect_true(all(v$vars %in% names(v$equations)))
  expect_true(is.numeric(v$aic))
  expect_true(is.numeric(v$bic))
  expect_equal(dim(v$granger_matrix), c(3, 3))
  
  out <- capture.output(print(v))
  expect_true(any(grepl("VECTOR AUTOREGRESSION: VAR\\(1\\) MODEL", out)))
  expect_true(any(grepl("PAIRWISE GRANGER CAUSALITY MATRIX", out)))
})
