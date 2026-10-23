test_that("iv() computes 2SLS, first stage F-stat, Hausman, and Sargan tests", {
  set.seed(123)
  n <- 300
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  u <- rnorm(n)
  educ <- 2 * z1 + 1.8 * z2 + u + rnorm(n, sd = 0.5)
  exper <- rnorm(n)
  wage <- 4 + 1.2 * educ + 0.5 * exper + 0.8 * u + rnorm(n)
  df <- data.frame(wage, educ, exper, z1, z2)
  
  # Run 2SLS
  m_iv <- iv(wage ~ exper | educ ~ z1 + z2, data = df)
  expect_s3_class(m_iv, "econ_iv")
  expect_equal(m_iv$n, 300)
  expect_true(m_iv$first_stage$educ$f_stat > 10) # Strong instruments
  expect_true(!is.na(m_iv$hausman_stat))
  expect_true(!is.na(m_iv$sargan_stat))
  expect_equal(m_iv$overid_df, 1) # 2 instruments - 1 endogenous = 1
  
  # Output printing
  expect_output(print(m_iv), "INSTRUMENTAL VARIABLES")
  expect_output(print(m_iv), "Weak IV")
  expect_output(print(m_iv), "Wu-Hausman Endogeneity")
  expect_output(print(m_iv), "Sargan Overid")
})
