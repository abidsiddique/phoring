#' Publication-Ready Summary Statistics and Covariate Balance for Economics
#'
#' Generates standard summary statistics (Obs, Mean, Std Dev, Min, P25, Median, P75, Max, Missing)
#' for numeric variables in a dataset, or computes a full "Table 1" covariate balance test across two 
#' treatment/subsample groups with normalized differences and t-tests.
#'
#' @param data A data.frame.
#' @param vars Optional character vector of column names to summarize.
#' @param by Optional character string specifying a binary grouping column (e.g. treatment, gender)
#'   for two-group covariate balance testing ("Table 1").
#' @param digits Number of decimal places to display (default is 3).
#'
#' @return A data.frame of class \code{econ_desc} or \code{econ_desc_balance}.
#' @export
#'
#' @examples
#' data(wage_data, package = "phoring")
#' # Standard descriptive statistics
#' desc(wage_data, vars = c("wage", "educ", "exper"))
#' 
#' # Covariate balance table ("Table 1") by gender
#' desc(wage_data, vars = c("wage", "educ", "exper"), by = "female")
desc <- function(data, vars = NULL, by = NULL, digits = 3) {
  if (missing(data) || is.null(data)) {
    stop("Please provide a dataset in `data`.", call. = FALSE)
  }
  
  if (is.null(vars)) {
    num_cols <- names(data)[sapply(data, is.numeric)]
    if (!is.null(by) && by %in% num_cols) {
      num_cols <- setdiff(num_cols, by)
    }
  } else {
    num_cols <- vars
  }
  
  if (length(num_cols) == 0) {
    stop("No numeric variables found in dataset.", call. = FALSE)
  }
  
  # If `by` is provided, compute Covariate Balance ("Table 1")
  if (!is.null(by)) {
    if (!by %in% names(data)) {
      stop(sprintf("Grouping column `%s` not found in data.", by), call. = FALSE)
    }
    grp_vec <- as.factor(data[[by]])
    lvls <- levels(grp_vec)
    if (length(lvls) != 2) {
      stop("`by` grouping variable must have exactly two distinct levels.", call. = FALSE)
    }
    
    grp0_idx <- which(grp_vec == lvls[1])
    grp1_idx <- which(grp_vec == lvls[2])
    
    rows <- list()
    for (v in num_cols) {
      y0 <- stats::na.omit(data[[v]][grp0_idx])
      y1 <- stats::na.omit(data[[v]][grp1_idx])
      
      n0 <- length(y0)
      n1 <- length(y1)
      
      if (n0 > 1 && n1 > 1) {
        m0 <- mean(y0)
        s0 <- stats::sd(y0)
        m1 <- mean(y1)
        s1 <- stats::sd(y1)
        
        diff_val <- m1 - m0
        tt <- stats::t.test(y1, y0)
        t_stat <- tt$statistic
        p_val <- tt$p.value
        
        # Imbens & Rubin Normalized Difference: (m1 - m0) / sqrt((s1^2 + s2^2)/2)
        norm_diff <- diff_val / sqrt((s0^2 + s1^2) / 2)
        balance_check <- if (abs(norm_diff) < 0.25) "Balanced" else "Imbalanced (|D| >= 0.25)"
        
        rows[[v]] <- data.frame(
          Variable    = v,
          Mean_G0     = m0,
          SD_G0       = s0,
          Mean_G1     = m1,
          SD_G1       = s1,
          Difference  = diff_val,
          t_stat      = as.numeric(t_stat),
          p_value     = p_val,
          Signif      = get_stars(p_val),
          Norm_Diff   = norm_diff,
          Balance     = balance_check,
          stringsAsFactors = FALSE
        )
      }
    }
    
    res_df <- do.call(rbind, rows)
    rownames(res_df) <- NULL
    res <- list(
      table = res_df,
      by = by,
      g0_name = lvls[1],
      g1_name = lvls[2],
      n0 = length(grp0_idx),
      n1 = length(grp1_idx),
      data = data,
      digits = digits
    )
    class(res) <- "econ_desc_balance"
    return(res)
  }
  
  # Standard Descriptive Statistics
  rows <- list()
  for (v in num_cols) {
    vec <- data[[v]]
    n_total <- length(vec)
    n_miss <- sum(is.na(vec))
    pct_miss <- (n_miss / n_total) * 100
    v_clean <- stats::na.omit(vec)
    n_obs <- length(v_clean)
    
    if (n_obs > 0) {
      m <- mean(v_clean)
      s <- stats::sd(v_clean)
      qs <- stats::quantile(v_clean, probs = c(0, 0.25, 0.5, 0.75, 1))
      
      rows[[v]] <- data.frame(
        Variable = v,
        Obs      = n_obs,
        Mean     = m,
        Std_Dev  = s,
        Min      = qs[1],
        P25      = qs[2],
        Median   = qs[3],
        P75      = qs[4],
        Max      = qs[5],
        Miss_Pct = pct_miss,
        stringsAsFactors = FALSE
      )
    }
  }
  
  res_df <- do.call(rbind, rows)
  rownames(res_df) <- NULL
  attr(res_df, "digits") <- digits
  attr(res_df, "data") <- data[, num_cols, drop = FALSE]
  class(res_df) <- c("econ_desc", "data.frame")
  return(res_df)
}

#' @export
print.econ_desc <- function(x, digits = NULL, ...) {
  d <- if (is.null(digits)) attr(x, "digits") else digits
  if (is.null(d)) d <- 3
  
  v_w <- max(16, max(nchar(x$Variable)))
  line_w <- max(74, v_w + 64)
  
  print_line("=", width = line_w)
  cat("                      DESCRIPTIVE STATISTICS                            \n")
  print_line("=", width = line_w)
  cat(sprintf(paste0("  %-", v_w, "s %6s %10s %10s %10s %10s %10s %10s\n"),
              "Variable", "Obs", "Mean", "Std. Dev.", "Min", "Median", "Max", "% Miss"))
  print_line("-", width = line_w)
  
  for (i in 1:nrow(x)) {
    row <- x[i, ]
    cat(sprintf(paste0("  %-", v_w, "s %6d %10.*f %10.*f %10.*f %10.*f %10.*f %9.1f%%\n"),
                row$Variable,
                row$Obs,
                d, row$Mean,
                d, row$Std_Dev,
                d, row$Min,
                d, row$Median,
                d, row$Max,
                row$Miss_Pct))
  }
  print_line("=", width = line_w)
  invisible(x)
}

#' @export
print.econ_desc_balance <- function(x, ...) {
  d <- x$digits
  tab <- x$table
  v_w <- max(16, max(nchar(tab$Variable)))
  line_w <- max(82, v_w + 72)
  
  print_line("=", width = line_w)
  cat("                   COVARIATE BALANCE TABLE ('TABLE 1')                  \n")
  cat(sprintf("  Group Variable: %s | Group 0 (%s): n = %d | Group 1 (%s): n = %d\n",
              x$by, x$g0_name, x$n0, x$g1_name, x$n1))
  print_line("=", width = line_w)
  
  cat(sprintf(paste0("  %-", v_w, "s %14s %14s %12s %10s %11s  %-12s\n"),
              "Variable", 
              paste0("G0 (", x$g0_name, ")"), 
              paste0("G1 (", x$g1_name, ")"), 
              "Difference", "p-value", "Norm Diff", "Balance Status"))
  print_line("-", width = line_w)
  
  for (i in 1:nrow(tab)) {
    row <- tab[i, ]
    g0_str <- sprintf("%.*f (%.*f)", d, row$Mean_G0, d, row$SD_G0)
    g1_str <- sprintf("%.*f (%.*f)", d, row$Mean_G1, d, row$SD_G1)
    diff_str <- sprintf("%.*f%s", d, row$Difference, row$Signif)
    
    cat(sprintf(paste0("  %-", v_w, "s %14s %14s %12s %10s %11.*f  %-12s\n"),
                row$Variable,
                g0_str,
                g1_str,
                diff_str,
                format_pval(row$p_value),
                d, row$Norm_Diff,
                row$Balance))
  }
  print_line("-", width = line_w)
  cat("  Signif. codes:  *** p < 0.001,  ** p < 0.01,  * p < 0.05,  . p < 0.1\n")
  cat("  Econometric Guidance:\n")
  cat("  * Imbens & Rubin (2015) rule of thumb: |Norm Diff| < 0.25 indicates acceptable balance.\n")
  cat("  * Statistically significant differences (p < 0.05) indicate potential selection bias.\n")
  print_line("=", width = line_w)
  invisible(x)
}

#' Correlation Matrix with Significance Stars
#'
#' Computes pairwise Pearson correlation coefficients with significance stars.
#'
#' @param data A data.frame.
#' @param vars Optional character vector of column names.
#' @param digits Number of decimal places (default is 3).
#'
#' @return An object of class \code{econ_corr}.
#' @export
#'
#' @examples
#' data(mtcars)
#' corr(mtcars, vars = c("mpg", "hp", "wt", "qsec"))
corr <- function(data, vars = NULL, digits = 3) {
  if (missing(data) || is.null(data)) {
    stop("Please provide a dataset in `data`.", call. = FALSE)
  }
  
  if (is.null(vars)) {
    num_cols <- names(data)[sapply(data, is.numeric)]
  } else {
    num_cols <- vars
  }
  
  if (length(num_cols) < 2) {
    stop("Need at least 2 numeric variables to compute correlations.", call. = FALSE)
  }
  
  sub_data <- data[, num_cols, drop = FALSE]
  p <- length(num_cols)
  r_mat <- matrix(1, nrow = p, ncol = p)
  p_mat <- matrix(0, nrow = p, ncol = p)
  colnames(r_mat) <- rownames(r_mat) <- num_cols
  colnames(p_mat) <- rownames(p_mat) <- num_cols
  
  for (i in 1:p) {
    for (j in 1:p) {
      if (i != j) {
        ct <- stats::cor.test(sub_data[[i]], sub_data[[j]], use = "complete.obs")
        r_mat[i, j] <- ct$estimate
        p_mat[i, j] <- ct$p.value
      }
    }
  }
  
  res <- list(
    r = r_mat,
    p = p_mat,
    vars = num_cols,
    digits = digits
  )
  class(res) <- "econ_corr"
  return(res)
}

#' @export
print.econ_corr <- function(x, digits = NULL, ...) {
  d <- if (is.null(digits)) x$digits else digits
  
  print_line("=")
  cat("                PEARSON CORRELATION MATRIX WITH STARS                   \n")
  print_line("=")
  
  p <- length(x$vars)
  v_w <- max(14, max(nchar(x$vars)))
  header <- sprintf(paste0("%-", v_w, "s"), "Variable")
  for (j in 1:p) {
    header <- paste0(header, sprintf("%12s", x$vars[j]))
  }
  cat(paste0("  ", header, "\n"))
  print_line("-", width = max(74, v_w + p * 12 + 4))
  
  for (i in 1:p) {
    row_str <- sprintf(paste0("  %-", v_w, "s"), x$vars[i])
    for (j in 1:p) {
      if (j > i) {
        row_str <- paste0(row_str, sprintf("%12s", ""))
      } else if (i == j) {
        row_str <- paste0(row_str, sprintf("%12s", "1.000   "))
      } else {
        r_val <- x$r[i, j]
        p_val <- x$p[i, j]
        stars <- get_stars(p_val)
        val_str <- sprintf("%.*f%s", d, r_val, stars)
        row_str <- paste0(row_str, sprintf("%12s", val_str))
      }
    }
    cat(paste0(row_str, "\n"))
  }
  print_line("-", width = max(74, v_w + p * 12 + 4))
  cat("  Signif. codes:  *** p < 0.001,  ** p < 0.01,  * p < 0.05\n")
  print_line("=", width = max(74, v_w + p * 12 + 4))
  invisible(x)
}
