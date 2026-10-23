test_that("diagnostics() computes all tests and VIFs correctly", {
  data(mtcars)
  m <- lm(mpg ~ wt + hp + qsec, data = mtcars)
  d <- diagnostics(m)
  
  expect_s3_class(d, "econ_diagnostics")
  expect_true(is.data.frame(d$tests))
  
  test_names <- d$tests$Test
  expect_true(any(grepl("Breusch-Pagan", test_names)))
  expect_true(any(grepl("White", test_names)))
  expect_true(any(grepl("Durbin-Watson", test_names)))
  expect_true(any(grepl("Breusch-Godfrey", test_names)))
  expect_true(any(grepl("Ramsey RESET", test_names)))
  expect_true(any(grepl("Jarque-Bera", test_names)))
  
  # Check VIFs
  expect_equal(length(d$vifs), 3)
  expect_true(all(d$vifs >= 1))
})
