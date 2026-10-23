test_that("panel() TWFE, First Differences, and panel diagnostics work", {
  data(firm_investment, package = "phoring")
  
  # Two-Way Fixed Effects
  m_twfe <- panel(inv ~ value + capital, data = firm_investment, id = "firm_id", time = "year", model = "twfe")
  expect_s3_class(m_twfe, "econ_panel")
  expect_equal(m_twfe$model, "twfe")
  expect_true(!is.na(m_twfe$hausman_stat))
  expect_true(!is.na(m_twfe$lm_bp_panel))
  expect_output(print(m_twfe), "TWO-WAY FIXED EFFECTS")
  expect_output(print(m_twfe), "Breusch-Pagan LM")
  
  # First Differences
  m_fd <- panel(inv ~ value + capital, data = firm_investment, id = "firm_id", time = "year", model = "fd")
  expect_s3_class(m_fd, "econ_panel")
  expect_equal(m_fd$model, "fd")
  expect_output(print(m_fd), "FIRST DIFFERENCES")
  
  # Between Estimator
  m_bet <- panel(inv ~ value + capital, data = firm_investment, id = "firm_id", time = "year", model = "between")
  expect_s3_class(m_bet, "econ_panel")
  expect_equal(m_bet$model, "between")
})

test_that("unit_root, coint, granger, arch_test, and ardl work for time series", {
  data(macro_data, package = "phoring")
  
  # 1. unit_root
  ur <- unit_root(macro_data$inflation)
  expect_s3_class(ur, "econ_unit_root")
  expect_equal(nrow(ur$tests), 3)
  expect_output(print(ur), "STATIONARITY & UNIT ROOT BATTERY")
  expect_output(print(ur), "Phillips-Perron")
  expect_output(print(ur), "KPSS")
  
  # 2. coint
  co <- coint(inflation ~ unemployment + fed_funds, data = macro_data)
  expect_s3_class(co, "econ_coint")
  expect_true(!is.na(co$tau_stat))
  expect_output(print(co), "ENGLE-GRANGER TWO-STEP COINTEGRATION")
  
  # 3. granger
  gr <- granger(inflation ~ fed_funds, data = macro_data, lags = 2)
  expect_s3_class(gr, "econ_granger")
  expect_true(!is.na(gr$f_xy))
  expect_true(!is.na(gr$f_yx))
  expect_output(print(gr), "GRANGER CAUSALITY TEST")
  
  # 4. arch_test
  m_ts <- reg(inflation ~ fed_funds, data = macro_data, diagnostics = FALSE)
  ar <- arch_test(m_ts, lags = 2)
  expect_s3_class(ar, "econ_arch")
  expect_true(!is.na(ar$stat))
  expect_output(print(ar), "ENGLE'S ARCH LM TEST")
  
  # 5. ardl
  m_ardl <- ardl(inflation ~ unemployment, data = macro_data, p = 1, q = 1)
  expect_s3_class(m_ardl, "econ_ardl")
  expect_true(is.numeric(m_ardl$impact_mult))
  expect_true(is.numeric(m_ardl$long_run_mult))
  expect_output(print(m_ardl), "ARDL(1, 1)", fixed = TRUE)
  expect_output(print(m_ardl), "Short-Run Impact Multiplier")
  expect_output(print(m_ardl), "Long-Run Cumulative Propensity")
})

test_that("tobit and count_reg work for advanced microeconometrics", {
  data(labor_supply, package = "phoring")
  
  # 1. tobit
  m_tobit <- tobit(hours ~ educ + exper + age, data = labor_supply, left = 0)
  expect_s3_class(m_tobit, "econ_tobit")
  expect_true(m_tobit$n_censored > 0)
  expect_equal(nrow(m_tobit$mm_table), 3)
  expect_output(print(m_tobit), "TOBIT CENSORED REGRESSION")
  expect_output(print(m_tobit), "McDONALD-MOFFITT MARGINAL EFFECTS")
  
  # 2. count_reg
  m_count <- count_reg(docvis ~ educ + age + kidslt6, data = labor_supply)
  expect_s3_class(m_count, "econ_count")
  expect_true(!is.na(m_count$disp_stat))
  expect_true("IRR" %in% names(m_count$coef_table))
  expect_output(print(m_count), "COUNT DATA REGRESSION: POISSON")
  expect_output(print(m_count), "CAMERON & TRIVEDI OVERDISPERSION TEST")
  
  # 3. logit Odds Ratios and AUC
  m_logit <- logit(inlf ~ educ + exper + age, data = labor_supply)
  expect_false(is.null(m_logit$odds_ratio))
  expect_false(is.na(m_logit$auc))
  expect_output(print(m_logit), "ODDS RATIOS")
  expect_output(print(m_logit), "ROC-AUC")
})
