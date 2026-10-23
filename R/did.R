#' Difference-in-Differences (DiD) Policy Evaluation
#'
#' Fits a Difference-in-Differences (DiD) model with optional control covariates,
#' computes cluster-robust standard errors, and outputs the classic 2x2 pre/post
#' treatment table alongside full regression estimates.
#'
#' @param formula Model formula specifying outcome, treatment, and post periods,
#'   optionally with control covariates (e.g., \code{y ~ treat * post} or
#'   \code{y ~ treat * post + x1 + x2}).
#' @param data A data.frame containing the variables.
#' @param treat Column name for treatment indicator (optional, automatically inferred from interaction).
#' @param post Column name for post-period indicator (optional, automatically inferred from interaction).
#' @param cluster Optional column name for cluster-robust standard errors (recommended: state/unit level).
#' @param digits Number of decimals to display (default is 4).
#' @param time Optional name of time period variable (e.g. \code{"year"}) for multi-period
#'   trajectories and event-study plotting.
#' @param ... Additional arguments passed to \code{\link{reg}}.
#'
#' @return An object of class \code{econ_did} containing the DiD effect, 2x2 table, and estimates.
#' @export
#'
#' @examples
#' set.seed(42)
#' n <- 400
#' state <- sample(1:20, n, replace = TRUE)
#' treat <- ifelse(state <= 10, 1, 0)
#' post <- sample(0:1, n, replace = TRUE)
#' x1 <- rnorm(n)
#' y <- 10 + 2 * treat + 1.5 * post + 3.2 * (treat * post) + 0.5 * x1 + rnorm(n, sd = 2)
#' df <- data.frame(y, treat, post, state, x1)
#' did(y ~ treat * post + x1, data = df, cluster = "state")
did <- function(formula, data, treat = NULL, post = NULL, time = NULL, cluster = NULL, digits = 4, ...) {
  if (missing(data) || is.null(data)) {
    stop("Please provide a dataset in `data`.", call. = FALSE)
  }
  
  vars <- all.vars(formula)
  y_var <- vars[1]
  
  # Identify treat and post variables
  term_labels <- attr(stats::terms(formula), "term.labels")
  inter_terms <- grep(":", term_labels, value = TRUE)
  
  treat_var <- treat
  post_var <- post
  
  if (is.null(treat_var) || is.null(post_var)) {
    if (length(inter_terms) > 0) {
      inter_vars <- strsplit(inter_terms[1], ":")[[1]]
      # If one is named post, time, after, year, period, treat it as post
      post_synonyms <- c("post", "after", "time", "period", "t", "year")
      is_post_match <- tolower(inter_vars) %in% post_synonyms
      if (any(is_post_match)) {
        post_var <- inter_vars[is_post_match][1]
        treat_var <- setdiff(inter_vars, post_var)[1]
      } else {
        treat_var <- inter_vars[1]
        post_var <- inter_vars[2]
      }
    } else {
      if (length(vars) < 3) {
        stop("`did()` requires treatment and post indicators either in the formula (e.g. y ~ treat * post) or via `treat` and `post` arguments.", call. = FALSE)
      }
      treat_var <- vars[2]
      post_var <- vars[3]
    }
  }
  
  # Ensure treat_var * post_var is included in regression formula
  has_interaction <- length(inter_terms) > 0 && any(grepl(paste0("(", treat_var, ":", post_var, "|", post_var, ":", treat_var, ")"), inter_terms))
  reg_formula <- formula
  if (!has_interaction) {
    rhs_terms <- attr(stats::terms(formula), "term.labels")
    covs <- setdiff(rhs_terms, c(treat_var, post_var))
    did_spec <- paste0(treat_var, " * ", post_var)
    new_rhs <- if (length(covs) > 0) paste(c(did_spec, covs), collapse = " + ") else did_spec
    reg_formula <- stats::as.formula(paste(y_var, "~", new_rhs))
  }
  
  # Keep all required variables (formula + treat + post + cluster + time)
  cluster_name <- if (is.character(cluster) && length(cluster) == 1) cluster else NULL
  time_name <- if (is.character(time) && length(time) == 1) time else NULL
  all_needed <- unique(c(all.vars(reg_formula), treat_var, post_var, cluster_name, time_name))
  missing_vars <- setdiff(all_needed, names(data))
  if (length(missing_vars) > 0) {
    stop(paste0("Variables not found in data: ", paste(missing_vars, collapse = ", ")), call. = FALSE)
  }
  
  sub_data <- stats::na.omit(data[, all_needed, drop = FALSE])
  
  y_vals <- sub_data[[y_var]]
  tr_raw <- sub_data[[treat_var]]
  po_raw <- sub_data[[post_var]]
  
  # Standardize 0/1 binary indicators
  if (is.logical(tr_raw)) {
    tr_vals <- as.numeric(tr_raw)
  } else if (is.numeric(tr_raw)) {
    u_tr <- sort(unique(tr_raw))
    if (length(u_tr) == 2 && all(u_tr == c(0, 1))) tr_vals <- tr_raw else tr_vals <- as.numeric(as.factor(tr_raw)) - 1
  } else {
    tr_vals <- as.numeric(as.factor(tr_raw)) - 1
  }
  
  if (is.logical(po_raw)) {
    po_vals <- as.numeric(po_raw)
  } else if (is.numeric(po_raw)) {
    u_po <- sort(unique(po_raw))
    if (length(u_po) == 2 && all(u_po == c(0, 1))) po_vals <- po_raw else po_vals <- as.numeric(as.factor(po_raw)) - 1
  } else {
    po_vals <- as.numeric(as.factor(po_raw)) - 1
  }
  
  # 2x2 Means Breakdown (Unadjusted)
  ctrl_pre  <- mean(y_vals[tr_vals == 0 & po_vals == 0], na.rm = TRUE)
  ctrl_post <- mean(y_vals[tr_vals == 0 & po_vals == 1], na.rm = TRUE)
  ctrl_diff <- ctrl_post - ctrl_pre
  
  treat_pre  <- mean(y_vals[tr_vals == 1 & po_vals == 0], na.rm = TRUE)
  treat_post <- mean(y_vals[tr_vals == 1 & po_vals == 1], na.rm = TRUE)
  treat_diff <- treat_post - treat_pre
  
  did_unadj <- treat_diff - ctrl_diff
  
  # Run regression with clustered or robust SE (handles covariates)
  se_mode <- if (!is.null(cluster)) "cluster" else "robust"
  reg_res <- reg(reg_formula, data = sub_data, se = se_mode, cluster = cluster, digits = digits, diagnostics = FALSE, ...)
  
  # Identify interaction row in regression table
  c_names <- rownames(reg_res$coef_table)
  inter_idx <- grep(paste0("(", treat_var, ":", post_var, "|", post_var, ":", treat_var, ")"), c_names)
  if (length(inter_idx) == 0) {
    inter_idx <- grep(":", c_names)
  }
  
  did_row_reg <- if (length(inter_idx) > 0) reg_res$coef_table[inter_idx[1], ] else NULL
  did_est <- if (!is.null(did_row_reg)) did_row_reg$Estimate else did_unadj
  
  table_2x2 <- data.frame(
    Group = c("Control (Untreated)", "Treated"),
    Pre   = c(ctrl_pre, treat_pre),
    Post  = c(ctrl_post, treat_post),
    Diff  = c(ctrl_diff, treat_diff),
    stringsAsFactors = FALSE
  )
  
  did_summary_row <- data.frame(
    Group = "Difference-in-Differences (Unadjusted)",
    Pre = NA,
    Post = NA,
    Diff = did_unadj,
    stringsAsFactors = FALSE
  )
  table_2x2 <- rbind(table_2x2, did_summary_row)
  
  res <- list(
    formula = formula,
    did_estimate = did_est,
    table_2x2 = table_2x2,
    reg_model = reg_res,
    coefficients = stats::coef(reg_res),
    vcov = stats::vcov(reg_res),
    residuals = stats::residuals(reg_res),
    fitted.values = stats::fitted(reg_res),
    cluster = cluster,
    treat_var = treat_var,
    post_var = post_var,
    time_var = time_name,
    y_vals = y_vals,
    tr_vals = tr_vals,
    po_vals = po_vals,
    data = sub_data,
    n = nrow(sub_data),
    digits = digits
  )
  class(res) <- "econ_did"
  return(res)
}

#' @export
print.econ_did <- function(x, digits = NULL, ...) {
  d <- if (is.null(digits)) x$digits else digits
  
  print_line("=")
  cat("                 DIFFERENCE-IN-DIFFERENCES (DiD)                        \n")
  cat(paste0("  Formula: ", deparse(x$formula), "\n"))
  cat(sprintf("  Sample Size: %-15d Clustered By: %s\n", 
              x$n, if (is.null(x$cluster)) "None (Robust HC1)" else x$cluster))
  cat(sprintf("  Treatment:   %-15s Post Period:  %s\n", x$treat_var, x$post_var))
  print_line("=")
  
  cat("  1. 2x2 TREATMENT & PERIOD MEANS BREAKDOWN (UNADJUSTED):\n")
  print_line("-")
  cat(sprintf("  %-38s %11s %11s %11s\n", "Group", "Pre-Period", "Post-Period", "Difference"))
  print_line("-")
  
  for (i in 1:nrow(x$table_2x2)) {
    row <- x$table_2x2[i, ]
    if (i == 3) print_line("-")
    cat(sprintf("  %-38s %11s %11s %11.*f\n",
                row$Group,
                if (is.na(row$Pre)) "-" else sprintf("%.*f", d, row$Pre),
                if (is.na(row$Post)) "-" else sprintf("%.*f", d, row$Post),
                d, row$Diff))
  }
  print_line("-")
  
  cat("\n  2. REGRESSION ESTIMATION (CAUSAL INTERACTION & COVARIATES):\n")
  print(x$reg_model, digits = d)
  invisible(x)
}

#' @export
coef.econ_did <- function(object, ...) {
  stats::coef(object$reg_model)
}

#' @export
vcov.econ_did <- function(object, ...) {
  stats::vcov(object$reg_model)
}

#' @export
residuals.econ_did <- function(object, ...) {
  stats::residuals(object$reg_model)
}

#' @export
fitted.econ_did <- function(object, ...) {
  stats::fitted(object$reg_model)
}

#' @export
nobs.econ_did <- function(object, ...) {
  object$n
}

#' @export
confint.econ_did <- function(object, parm, level = 0.95, ...) {
  stats::confint(object$reg_model, parm = parm, level = level, ...)
}

#' @export
predict.econ_did <- function(object, newdata = NULL, ...) {
  stats::predict(object$reg_model, newdata = newdata, ...)
}
