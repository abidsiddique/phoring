test_that("logit() and probit() compute AME and Pseudo R2", {
  data(mtcars)
  
  # Logit
  m_l <- logit(am ~ mpg + hp, data = mtcars)
  expect_s3_class(m_l, "econ_logit")
  expect_equal(m_l$link, "logit")
  expect_true(is.data.frame(m_l$ame_table))
  expect_equal(nrow(m_l$ame_table), 2)
  expect_output(print(m_l), "LOGISTIC REGRESSION")
  expect_output(print(m_l), "AVERAGE MARGINAL EFFECTS")
  
  # Probit
  m_p <- probit(am ~ mpg + hp, data = mtcars)
  expect_s3_class(m_p, "econ_logit")
  expect_equal(m_p$link, "probit")
  expect_output(print(m_p), "PROBIT REGRESSION")
})

test_that("did() computes 2x2 table and treatment interaction effect", {
  set.seed(42)
  n <- 200
  state <- sample(1:10, n, replace = TRUE)
  treat <- ifelse(state <= 5, 1, 0)
  post <- sample(0:1, n, replace = TRUE)
  y <- 5 + 1.5 * treat + 2 * post + 4 * (treat * post) + rnorm(n)
  df <- data.frame(y, treat, post, state)
  
  m_did <- did(y ~ treat * post, data = df, cluster = "state")
  expect_s3_class(m_did, "econ_did")
  expect_equal(nrow(m_did$table_2x2), 3)
  expect_output(print(m_did), "DIFFERENCE-IN-DIFFERENCES")
  expect_output(print(m_did), "2x2 TREATMENT & PERIOD MEANS")
})

test_that("desc(), corr(), test_hyp(), and word_table() work properly", {
  data(mtcars)
  
  # desc
  d <- desc(mtcars, vars = c("mpg", "wt"))
  expect_s3_class(d, "econ_desc")
  expect_equal(nrow(d), 2)
  
  # corr
  c_mat <- corr(mtcars, vars = c("mpg", "wt"))
  expect_s3_class(c_mat, "econ_corr")
  
  # test_hyp
  m <- reg(mpg ~ wt + hp, data = mtcars, diagnostics = FALSE)
  h <- test_hyp(m, "wt = 0")
  expect_s3_class(h, "econ_hyp")
  expect_true(h$F_stat > 0)
  
  # word_table
  tmp_docx <- tempfile(fileext = ".docx")
  expect_output(word_table(m, file = tmp_docx), "Regression Results")
  expect_true(file.exists(tmp_docx))
})
