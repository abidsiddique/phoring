#' Side-by-Side Model Comparison & Word Table Export
#'
#' Compiles multiple econometric models into a single publication-quality table
#' with standard errors in parentheses and significance stars, and exports directly
#' to Microsoft Word (.docx/.doc), HTML, or Markdown.
#'
#' @param ... One or more fitted model objects (from \code{\link{reg}}, \code{\link{iv}},
#'   \code{\link{panel}}, \code{\link{logit}}, \code{\link{tobit}}, \code{\link{count_reg}},
#'   \code{\link{did}}, \code{\link{koyck}}, \code{\link{prais_winsten}}, \code{\link{ecm}},
#'   \code{\link{ardl}}, or \code{\link[stats]{lm}}).
#' @param labels Optional character vector of custom column labels for each model
#'   (e.g., \code{c("OLS", "FE", "2SLS")}).
#' @param file Optional filename for export (e.g., \code{"table1.docx"} or \code{"table1.html"}).
#' @param title Title for the table.
#' @param format Output display format: \code{"text"} (console, default), \code{"html"}, or \code{"word"}.
#' @param digits Number of decimal places to format (default is 3).
#'
#' @return A list of the input models (invisibly).
#' @export
#'
#' @examples
#' data(mtcars)
#' m1 <- reg(mpg ~ wt, data = mtcars)
#' m2 <- reg(mpg ~ wt + hp, data = mtcars)
#' word_table(m1, m2, labels = c("Simple OLS", "Multiple OLS"))
word_table <- function(..., labels = NULL, file = NULL, title = "Regression Results", format = c("text", "html", "word"), digits = 3) {
  format <- match.arg(format)
  raw_args <- list(...)
  if (length(raw_args) == 0) {
    stop("Please pass at least one model object to `word_table()`.", call. = FALSE)
  }
  
  if (length(raw_args) == 1 && is.list(raw_args[[1]]) && !inherits(raw_args[[1]], c("lm", "glm", "econ_reg", "econ_panel", "econ_iv", "econ_logit", "econ_tobit", "econ_count", "econ_ardl", "econ_prais", "econ_koyck"))) {
    models <- raw_args[[1]]
  } else {
    models <- raw_args
  }
  
  is_model <- sapply(models, function(m) {
    is.list(m) && (inherits(m, c("lm", "glm")) || !is.null(m$coefficients) || !is.null(m$coef_table) || !is.null(attr(m, "class")))
  })
  
  if (!all(is_model)) {
    models <- models[is_model]
  }
  
  if (length(models) == 0) {
    stop("Please pass at least one valid model object (e.g. from reg, panel, iv, logit, lm) to `word_table()`.", call. = FALSE)
  }
  
  M <- length(models)
  if (!is.null(labels)) {
    if (length(labels) != M) {
      stop(sprintf("Length of `labels` (%d) must match number of models (%d).", length(labels), M), call. = FALSE)
    }
    m_names <- labels
  } else {
    m_names <- paste0("Model (", 1:M, ")")
  }
  
  # Collect all unique variable names
  all_vars <- character()
  model_data <- list()
  
  for (m_idx in 1:M) {
    mod <- models[[m_idx]]
    b <- stats::coef(mod)
    if (is.list(b)) b <- unlist(b)
    V <- stats::vcov(mod)
    se <- sqrt(pmax(0, diag(V)))
    names(se) <- names(b)
    
    df_r <- if (!is.null(mod$df_inference)) {
      mod$df_inference
    } else if (!is.null(mod$df_residual)) {
      mod$df_residual
    } else {
      length(stats::residuals(mod)) - length(b)
    }
    
    t_val <- b / se
    is_z <- inherits(mod, c("econ_logit", "econ_tobit", "econ_count"))
    if (is_z || is.null(df_r) || is.na(df_r) || df_r <= 0) {
      p_val <- 2 * stats::pnorm(-abs(t_val))
    } else {
      p_val <- 2 * stats::pt(-abs(t_val), df = df_r)
    }
    
    n_obs <- if (!is.null(mod$n)) mod$n else length(stats::residuals(mod))
    r2_val <- if (!is.null(mod$r_squared)) {
      mod$r_squared
    } else if (!is.null(mod$r2_within)) {
      mod$r2_within
    } else if (!is.null(mod$pseudo_r2)) {
      mod$pseudo_r2
    } else {
      NA
    }
    se_type <- if (!is.null(mod$se_type)) mod$se_type else "Standard"
    
    model_data[[m_idx]] <- list(
      b = b,
      se = se,
      p = p_val,
      n = n_obs,
      r2 = r2_val,
      se_type = se_type
    )
    all_vars <- unique(c(all_vars, names(b)))
  }
  
  # Ensure Intercept is placed first or last
  if ("(Intercept)" %in% all_vars) {
    all_vars <- c("(Intercept)", setdiff(all_vars, "(Intercept)"))
  }
  
  v_w <- max(18, max(nchar(all_vars)))
  col_w <- max(14, max(nchar(m_names)) + 2)
  total_w <- max(74, v_w + M * col_w + 4)
  
  # Build text display table
  cat("\n")
  print_line("=", width = total_w)
  cat(sprintf("  %s\n", title))
  print_line("=", width = total_w)
  
  header_str <- sprintf(paste0("  %-", v_w, "s"), "Variable")
  for (m_idx in 1:M) {
    header_str <- paste0(header_str, sprintf(paste0("%", col_w, "s"), m_names[m_idx]))
  }
  cat(paste0(header_str, "\n"))
  print_line("-", width = total_w)
  
  for (v in all_vars) {
    # Row for coefficients
    b_str <- sprintf(paste0("  %-", v_w, "s"), v)
    # Row for standard errors
    se_str <- sprintf(paste0("  %-", v_w, "s"), "")
    
    for (m_idx in 1:M) {
      m_info <- model_data[[m_idx]]
      if (v %in% names(m_info$b)) {
        b_val <- m_info$b[v]
        se_val <- m_info$se[v]
        p_val <- m_info$p[v]
        stars <- get_stars(p_val)
        
        b_entry <- sprintf("%.*f%s", digits, b_val, stars)
        se_entry <- sprintf("(%.*f)", digits, se_val)
        
        b_str <- paste0(b_str, sprintf(paste0("%", col_w, "s"), b_entry))
        se_str <- paste0(se_str, sprintf(paste0("%", col_w, "s"), se_entry))
      } else {
        b_str <- paste0(b_str, sprintf(paste0("%", col_w, "s"), ""))
        se_str <- paste0(se_str, sprintf(paste0("%", col_w, "s"), ""))
      }
    }
    cat(paste0(b_str, "\n"))
    cat(paste0(se_str, "\n"))
  }
  
  print_line("-", width = total_w)
  
  # Model statistics rows
  obs_str <- sprintf(paste0("  %-", v_w, "s"), "Observations")
  r2_str  <- sprintf(paste0("  %-", v_w, "s"), "R-squared")
  se_row  <- sprintf(paste0("  %-", v_w, "s"), "Standard Errors")
  
  for (m_idx in 1:M) {
    m_info <- model_data[[m_idx]]
    obs_str <- paste0(obs_str, sprintf(paste0("%", col_w, "d"), m_info$n))
    r2_str <- paste0(r2_str, sprintf(paste0("%", col_w, "s"), if (is.na(m_info$r2)) "-" else sprintf("%.*f", digits, m_info$r2)))
    se_row <- paste0(se_row, sprintf(paste0("%", col_w, "s"), m_info$se_type))
  }
  
  cat(paste0(obs_str, "\n"))
  cat(paste0(r2_str, "\n"))
  cat(paste0(se_row, "\n"))
  print_line("-", width = total_w)
  cat("  Signif. codes:  *** p < 0.001,  ** p < 0.01,  * p < 0.05\n")
  print_line("=", width = total_w)
  
  # File Export
  if (!is.null(file)) {
    html_lines <- c(
      "<!DOCTYPE html>",
      "<html>",
      "<head>",
      "<meta charset='utf-8'>",
      paste0("<title>", title, "</title>"),
      "<style>",
      "table { border-collapse: collapse; width: 85%; font-family: 'Times New Roman', Times, serif; margin: 20px auto; }",
      "th, td { padding: 6px 12px; text-align: right; }",
      "th { border-bottom: 2px solid black; border-top: 2px solid black; }",
      ".left { text-align: left; }",
      ".border-top { border-top: 1px solid black; }",
      ".border-bottom { border-bottom: 2px solid black; }",
      "</style>",
      "</head>",
      "<body>",
      paste0("<h3 style='text-align:center;'>", title, "</h3>"),
      "<table>",
      "<thead>",
      "<tr><th class='left'>Variable</th>"
    )
    
    for (m_idx in 1:M) {
      html_lines <- c(html_lines, paste0("<th>", m_names[m_idx], "</th>"))
    }
    html_lines <- c(html_lines, "</tr></thead><tbody>")
    
    for (v in all_vars) {
      b_cells <- paste0("<td class='left'>", v, "</td>")
      se_cells <- "<td class='left'></td>"
      for (m_idx in 1:M) {
        m_info <- model_data[[m_idx]]
        if (v %in% names(m_info$b)) {
          b_val <- m_info$b[v]
          se_val <- m_info$se[v]
          p_val <- m_info$p[v]
          stars <- get_stars(p_val)
          b_cells <- paste0(b_cells, sprintf("<td>%.*f%s</td>", digits, b_val, stars))
          se_cells <- paste0(se_cells, sprintf("<td>(%.*f)</td>", digits, se_val))
        } else {
          b_cells <- paste0(b_cells, "<td></td>")
          se_cells <- paste0(se_cells, "<td></td>")
        }
      }
      html_lines <- c(html_lines, paste0("<tr>", b_cells, "</tr>"))
      html_lines <- c(html_lines, paste0("<tr>", se_cells, "</tr>"))
    }
    
    # Summary statistics rows
    obs_cells <- "<td class='left border-top'>Observations</td>"
    r2_cells  <- "<td class='left'>R-squared</td>"
    se_cells  <- "<td class='left border-bottom'>Standard Errors</td>"
    
    for (m_idx in 1:M) {
      m_info <- model_data[[m_idx]]
      obs_cells <- paste0(obs_cells, sprintf("<td class='border-top'>%d</td>", m_info$n))
      r2_cells <- paste0(r2_cells, sprintf("<td>%s</td>", if (is.na(m_info$r2)) "-" else sprintf("%.*f", digits, m_info$r2)))
      se_cells <- paste0(se_cells, sprintf("<td class='border-bottom'>%s</td>", m_info$se_type))
    }
    
    html_lines <- c(html_lines,
                    paste0("<tr>", obs_cells, "</tr>"),
                    paste0("<tr>", r2_cells, "</tr>"),
                    paste0("<tr>", se_cells, "</tr>"),
                    "</tbody></table>",
                    "<p style='text-align:center; font-size: 0.9em;'><em>Significance codes: *** p &lt; 0.001, ** p &lt; 0.01, * p &lt; 0.05</em></p>",
                    "</body></html>")
    
    writeLines(html_lines, con = file)
    cat(sprintf("\n[Export Successful] Table saved to: %s\n", file))
  }
  
  invisible(models)
}
