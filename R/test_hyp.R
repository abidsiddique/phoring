#' General Linear Hypothesis & Wald F-Test
#'
#' Tests linear hypotheses of model parameters using the robust or clustered
#' covariance matrix from the estimated model. Supports interaction terms,
#' factor levels, polynomial terms, and multi-parameter restrictions.
#'
#' @param model A fitted model object from \code{\link{reg}}, \code{\link{iv}},
#'   \code{\link{panel}}, \code{\link{logit}}, \code{\link{tobit}}, or \code{\link[stats]{lm}}.
#' @param hypothesis Character vector of linear restrictions (e.g., \code{"educ + exper = 1"},
#'   \code{"treat:post = 0"}, or \code{c("wt = 0", "hp = 0")}).
#' @param ... Additional arguments.
#'
#' @return An object of class \code{econ_hyp} containing the Wald test statistic,
#'   degrees of freedom, and p-value.
#' @export
#'
#' @examples
#' data(mtcars)
#' m <- reg(mpg ~ wt + hp + qsec, data = mtcars)
#' test_hyp(m, "wt + hp = 0")
#' test_hyp(m, c("wt = 0", "hp = 0"))
test_hyp <- function(model, hypothesis, ...) {
  if (missing(hypothesis) || length(hypothesis) == 0) {
    stop("Please specify a hypothesis string (e.g. 'x1 + x2 = 1' or 'treat:post = 0').", call. = FALSE)
  }
  
  b <- stats::coef(model)
  if (is.list(b)) b <- unlist(b)
  V <- stats::vcov(model)
  
  # Determine residual degrees of freedom
  df_r <- if (!is.null(model$df_inference)) {
    model$df_inference
  } else if (!is.null(model$df_residual)) {
    model$df_residual
  } else {
    length(stats::residuals(model)) - length(b)
  }
  
  is_z_model <- inherits(model, c("econ_logit", "econ_tobit", "econ_count"))
  
  p_names <- names(b)
  k <- length(b)
  J <- length(hypothesis)
  
  R_mat <- matrix(0, nrow = J, ncol = k)
  colnames(R_mat) <- p_names
  q_vec <- numeric(J)
  
  # Sort names by descending length so complex/interaction terms match first
  order_p <- order(nchar(p_names), decreasing = TRUE)
  sorted_names <- p_names[order_p]
  sorted_idx <- seq_along(p_names)[order_p]
  
  for (i in 1:J) {
    hyp_str <- hypothesis[i]
    if (!grepl("=", hyp_str)) {
      lhs <- hyp_str
      rhs <- "0"
    } else {
      parts <- strsplit(hyp_str, "=")[[1]]
      lhs <- parts[1]
      rhs <- parts[2]
    }
    
    q_val <- as.numeric(trimws(rhs))
    if (is.na(q_val)) {
      lhs <- paste0(lhs, " - (", rhs, ")")
      q_val <- 0
    }
    q_vec[i] <- q_val
    
    # Tokenize LHS expression: replace each parameter name with ..p_j..
    lhs_tokenized <- lhs
    for (m_idx in seq_along(sorted_names)) {
      orig_name <- sorted_names[m_idx]
      j <- sorted_idx[m_idx]
      token <- paste0("..p_", j, "..")
      
      # Escape regex special characters in variable name
      esc_name <- gsub("([.|()\\^\\$+*?\\[\\]\\-])", "\\\\\\1", orig_name)
      
      # Match with or without backticks
      pattern_backtick <- paste0("`", esc_name, "`")
      lhs_tokenized <- gsub(pattern_backtick, token, lhs_tokenized, fixed = TRUE)
      
      pattern_word <- paste0("(?<![a-zA-Z0-9_.])", esc_name, "(?![a-zA-Z0-9_.])")
      lhs_tokenized <- gsub(pattern_word, token, lhs_tokenized, perl = TRUE)
    }
    
    # Parse tokenized expression
    expr <- tryCatch(parse(text = lhs_tokenized)[[1]], error = function(e) {
      stop(paste0("Could not parse linear hypothesis: '", hypothesis[i], "'.\nCheck variable names against: ",
                  paste(p_names, collapse = ", ")), call. = FALSE)
    })
    
    # Linear evaluation via finite differences on tokens
    env0 <- as.list(rep(0, k))
    names(env0) <- paste0("..p_", 1:k, "..")
    val0 <- eval(expr, envir = env0)
    
    for (j in 1:k) {
      env_j <- env0
      env_j[[paste0("..p_", j, "..")]] <- 1
      val_j <- eval(expr, envir = env_j)
      R_mat[i, j] <- as.numeric(val_j - val0)
    }
  }
  
  # Wald Statistic: W = (R b - q)' [R V R']^(-1) (R b - q)
  diff_vec <- as.vector(R_mat %*% b - q_vec)
  middle <- R_mat %*% V %*% t(R_mat)
  middle_inv <- safe_solve(middle)
  
  chi2_stat <- as.numeric(t(diff_vec) %*% middle_inv %*% diff_vec)
  chi2_stat <- max(0, chi2_stat)
  
  if (is_z_model || is.null(df_r) || is.na(df_r) || df_r <= 0) {
    p_val <- stats::pchisq(chi2_stat, df = J, lower.tail = FALSE)
    f_stat <- chi2_stat / J
    df_denom <- NA_integer_
  } else {
    f_stat <- chi2_stat / J
    p_val <- stats::pf(f_stat, df1 = J, df2 = df_r, lower.tail = FALSE)
    df_denom <- df_r
  }
  
  res <- list(
    hypothesis = hypothesis,
    R = R_mat,
    q = q_vec,
    F_stat = f_stat,
    chi2_stat = chi2_stat,
    df1 = J,
    df2 = df_denom,
    p_value = p_val,
    is_z = is_z_model
  )
  class(res) <- "econ_hyp"
  return(res)
}

#' @export
print.econ_hyp <- function(x, ...) {
  print_line("=")
  if (isTRUE(x$is_z)) {
    cat("              LINEAR HYPOTHESIS WALD CHI-SQUARE TEST                    \n")
  } else {
    cat("                 LINEAR HYPOTHESIS WALD F-TEST                          \n")
  }
  print_line("=")
  cat("  Hypothesis Tested:\n")
  for (i in seq_along(x$hypothesis)) {
    cat(sprintf("   [%d] %s\n", i, x$hypothesis[i]))
  }
  print_line("-")
  if (isTRUE(x$is_z)) {
    cat(sprintf("  Chi-Square:   %10.4f  (df = %d)\n", x$chi2_stat, x$df1))
  } else {
    cat(sprintf("  F-statistic:  %10.4f  (df1 = %d, df2 = %d)\n", x$F_stat, x$df1, x$df2))
  }
  cat(sprintf("  p-value:      %10s   %s\n", format_pval(x$p_value), get_stars(x$p_value)))
  print_line("-")
  verdict <- if (x$p_value < 0.05) {
    "REJECT H0 (p < 0.05) -> The restriction does NOT hold"
  } else {
    "FAIL TO REJECT H0 (p >= 0.05) -> The restriction is consistent with data"
  }
  cat(sprintf("  Conclusion:   %s\n", verdict))
  print_line("=")
  invisible(x)
}
