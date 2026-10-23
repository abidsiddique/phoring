test_that("panel() computes fixed effects, random effects, and Hausman test", {
  set.seed(42)
  n_firms <- 30
  n_years <- 4
  id <- rep(1:n_firms, each = n_years)
  time <- rep(1:n_years, times = n_firms)
  firm_fe <- rep(rnorm(n_firms, mean = 2), each = n_years)
  x1 <- rnorm(n_firms * n_years) + 0.6 * firm_fe
  x2 <- rnorm(n_firms * n_years)
  y <- 2 + 1.2 * x1 + 0.5 * x2 + firm_fe + rnorm(n_firms * n_years)
  df <- data.frame(id, time, y, x1, x2)
  
  # Fixed Effects
  m_fe <- panel(y ~ x1 + x2, data = df, id = "id", time = "time", model = "fe")
  expect_s3_class(m_fe, "econ_panel")
  expect_equal(m_fe$G, 30)
  expect_equal(m_fe$n, 120)
  expect_true(!is.na(m_fe$hausman_stat))
  
  # Random Effects
  m_re <- panel(y ~ x1 + x2, data = df, id = "id", time = "time", model = "re")
  expect_s3_class(m_re, "econ_panel")
  expect_true(!is.null(m_re$theta))
  
  # Output printing
  expect_output(print(m_fe), "FIXED EFFECTS")
  expect_output(print(m_fe), "Hausman Test")
})
