test_that("reg() works with standard, robust, clustered, and hac standard errors", {
  data(mtcars)
  
  # Standard OLS
  m_std <- reg(mpg ~ wt + hp, data = mtcars, se = "standard", diagnostics = FALSE)
  expect_s3_class(m_std, "econ_reg")
  expect_equal(m_std$n, 32)
  expect_true(m_std$r_squared > 0.8)
  
  # Robust OLS (White HC1)
  m_rob <- reg(mpg ~ wt + hp, data = mtcars, se = "robust", diagnostics = FALSE)
  expect_s3_class(m_rob, "econ_reg")
  expect_true(grepl("Robust", m_rob$se_type))
  # Point estimates must match standard OLS
  expect_equal(unname(m_std$coefficients), unname(m_rob$coefficients))
  # Robust SEs should differ slightly from standard
  expect_false(isTRUE(all.equal(m_std$coef_table$Std_Error, m_rob$coef_table$Std_Error)))
  
  # Clustered OLS
  m_clust <- reg(mpg ~ wt + hp, data = mtcars, cluster = "cyl", diagnostics = FALSE)
  expect_s3_class(m_clust, "econ_reg")
  expect_true(grepl("Clustered", m_clust$se_type))
  
  # HAC OLS
  m_hac <- reg(mpg ~ wt + hp, data = mtcars, se = "hac", lag = 2, diagnostics = FALSE)
  expect_s3_class(m_hac, "econ_reg")
  expect_true(grepl("HAC", m_hac$se_type))
})

test_that("reg() prints all-in-one output including diagnostics", {
  data(mtcars)
  m <- reg(mpg ~ wt + hp, data = mtcars)
  expect_false(is.null(m$diagnostics))
  expect_true(is.data.frame(m$diagnostics$tests))
  expect_output(print(m), "OLS REGRESSION RESULTS")
  expect_output(print(m), "ALL-IN-ONE DIAGNOSTIC TESTS")
  expect_output(print(m), "Breusch-Pagan")
  expect_output(print(m), "Durbin-Watson")
  expect_output(print(m), "Ramsey RESET")
})
