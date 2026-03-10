#' Input Parser Functions for bayesianatoR
#'
#' @description
#' Functions for converting between common reporting formats used in clinical
#' papers (confidence intervals, standard errors, standard deviations,
#' p-values) and the sufficient statistics needed for Bayesian updating.
#' Log-transform handling is included for ratio measures (OR, RR, HR).
#'
#' @name input_parsers
NULL

# ── CI → SE ──────────────────────────────────────────────────────────────────

#' Convert a confidence interval to a standard error
#'
#' @param lower Numeric. Lower bound of the confidence interval.
#' @param upper Numeric. Upper bound of the confidence interval.
#' @param conf_level Numeric. Confidence level (0–1). Default `0.95`.
#'
#' @return Numeric scalar: the standard error.
#' @export
#'
#' @examples
#' ci_to_se(lower = 1.2, upper = 3.8, conf_level = 0.95)
ci_to_se <- function(lower, upper, conf_level = 0.95) {
  rlang::check_required(lower)
  rlang::check_required(upper)

  if (!is.numeric(lower) || !is.numeric(upper)) {
    rlang::abort("`lower` and `upper` must be numeric.", .call = rlang::caller_env())
  }
  if (lower >= upper) {
    rlang::abort("`lower` must be less than `upper`.", .call = rlang::caller_env())
  }
  if (!is.numeric(conf_level) || conf_level <= 0 || conf_level >= 1) {
    rlang::abort(
      "`conf_level` must be a number strictly between 0 and 1.",
      .call = rlang::caller_env()
    )
  }

  z <- qnorm(1 - (1 - conf_level) / 2)
  (upper - lower) / (2 * z)
}

# ── SE → CI ──────────────────────────────────────────────────────────────────

#' Convert an estimate and standard error to a confidence interval
#'
#' @param estimate Numeric. Point estimate (mean, difference, log ratio, etc.).
#' @param se Numeric. Standard error of the estimate.
#' @param conf_level Numeric. Confidence level (0–1). Default `0.95`.
#'
#' @return Named numeric vector with elements `lower` and `upper`.
#' @export
#'
#' @examples
#' se_to_ci(estimate = 2.5, se = 0.66, conf_level = 0.95)
se_to_ci <- function(estimate, se, conf_level = 0.95) {
  rlang::check_required(estimate)
  rlang::check_required(se)

  if (!is.numeric(estimate) || !is.numeric(se)) {
    rlang::abort("`estimate` and `se` must be numeric.", .call = rlang::caller_env())
  }
  if (se <= 0) {
    rlang::abort("`se` must be positive.", .call = rlang::caller_env())
  }
  if (!is.numeric(conf_level) || conf_level <= 0 || conf_level >= 1) {
    rlang::abort(
      "`conf_level` must be a number strictly between 0 and 1.",
      .call = rlang::caller_env()
    )
  }

  z <- qnorm(1 - (1 - conf_level) / 2)
  c(lower = estimate - z * se, upper = estimate + z * se)
}

# ── SD → SE (two-group) ───────────────────────────────────────────────────────

#' Convert standard deviations to standard error for a two-group comparison
#'
#' Computes the SE of the mean difference under the assumption of equal
#' variances (pooled SD) or unequal variances (Welch-style).
#'
#' @param sd1 Numeric. Standard deviation in group 1.
#' @param sd2 Numeric. Standard deviation in group 2.
#' @param n1 Integer. Sample size in group 1.
#' @param n2 Integer. Sample size in group 2.
#' @param pooled Logical. Use pooled SD? Default `TRUE`.
#'
#' @return Numeric scalar: the standard error of the mean difference.
#' @export
#'
#' @examples
#' sd_to_se(sd1 = 2.1, sd2 = 1.9, n1 = 50, n2 = 48)
sd_to_se <- function(sd1, sd2, n1, n2, pooled = TRUE) {
  for (nm in c("sd1", "sd2")) {
    val <- get(nm)
    if (!is.numeric(val) || val <= 0) {
      rlang::abort(
        sprintf("`%s` must be a positive number.", nm),
        .call = rlang::caller_env()
      )
    }
  }
  for (nm in c("n1", "n2")) {
    val <- get(nm)
    if (!is.numeric(val) || val < 2 || val != round(val)) {
      rlang::abort(
        sprintf("`%s` must be an integer >= 2.", nm),
        .call = rlang::caller_env()
      )
    }
  }

  if (pooled) {
    sp <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
    sp * sqrt(1 / n1 + 1 / n2)
  } else {
    sqrt(sd1^2 / n1 + sd2^2 / n2)
  }
}

# ── p-value → SE ─────────────────────────────────────────────────────────────

#' Approximate standard error from a p-value and point estimate
#'
#' This is an approximation and should be used with caution. A warning is
#' always issued when this function is called. Prefer CI- or SD-based inputs
#' wherever available.
#'
#' @param estimate Numeric. Point estimate.
#' @param pvalue Numeric. Two-sided p-value (strictly between 0 and 1).
#' @param df Numeric or `NULL`. Degrees of freedom. If `NULL`, the
#'   z-distribution is used (appropriate for large samples).
#'
#' @return Numeric scalar: approximate standard error.
#' @export
#'
#' @examples
#' pvalue_to_se(estimate = 2.5, pvalue = 0.03)
pvalue_to_se <- function(estimate, pvalue, df = NULL) {
  rlang::check_required(estimate)
  rlang::check_required(pvalue)

  if (!is.numeric(pvalue) || pvalue <= 0 || pvalue >= 1) {
    rlang::abort(
      "`pvalue` must be strictly between 0 and 1.",
      .call = rlang::caller_env()
    )
  }
  if (!is.numeric(estimate)) {
    rlang::abort("`estimate` must be numeric.", .call = rlang::caller_env())
  }

  rlang::warn(
    paste0(
      "SE derived from p-value is approximate. ",
      "Prefer CI- or SD-based inputs where available."
    ),
    .frequency    = "regularly",
    .frequency_id = "pvalue_se_approx"
  )

  if (is.null(df)) {
    z <- qnorm(1 - pvalue / 2)
    abs(estimate) / z
  } else {
    t_stat <- qt(1 - pvalue / 2, df = df)
    abs(estimate) / t_stat
  }
}

# ── Log-transform helpers ─────────────────────────────────────────────────────

#' Log-transform a ratio estimate and its confidence interval
#'
#' Converts OR, RR, or HR (on the natural scale) to the log scale, returning
#' the log estimate and SE derived from the CI.
#'
#' @param estimate Numeric. Ratio estimate on the natural scale (must be > 0).
#' @param lower Numeric. Lower CI bound on the natural scale (must be > 0).
#' @param upper Numeric. Upper CI bound on the natural scale (must be > 0).
#' @param conf_level Numeric. Confidence level (0–1). Default `0.95`.
#'
#' @return Named list with elements:
#'   * `log_estimate` — log-transformed point estimate
#'   * `log_se`       — SE on the log scale
#'   * `log_lower`    — lower CI on the log scale
#'   * `log_upper`    — upper CI on the log scale
#' @export
#'
#' @examples
#' log_transform_ratio(estimate = 1.5, lower = 1.1, upper = 2.1)
log_transform_ratio <- function(estimate, lower, upper, conf_level = 0.95) {
  for (nm in c("estimate", "lower", "upper")) {
    val <- get(nm)
    if (!is.numeric(val) || any(val <= 0)) {
      rlang::abort(
        paste0(
          "All ratio inputs must be positive (> 0). ",
          "Did you forget to exponentiate?"
        ),
        .call = rlang::caller_env()
      )
    }
  }
  if (lower >= upper) {
    rlang::abort("`lower` must be less than `upper`.", .call = rlang::caller_env())
  }

  log_est <- log(estimate)
  log_lo  <- log(lower)
  log_hi  <- log(upper)
  log_se  <- ci_to_se(log_lo, log_hi, conf_level = conf_level)

  list(
    log_estimate = log_est,
    log_se       = log_se,
    log_lower    = log_lo,
    log_upper    = log_hi
  )
}

#' Back-transform log-scale results to the natural scale
#'
#' @param log_estimate Numeric. Estimate on the log scale.
#' @param log_lower Numeric. Lower CI on the log scale.
#' @param log_upper Numeric. Upper CI on the log scale.
#'
#' @return Named list with `estimate`, `lower`, `upper` on the natural scale.
#' @export
#'
#' @examples
#' exp_transform_result(log_estimate = 0.405, log_lower = 0.095, log_upper = 0.742)
exp_transform_result <- function(log_estimate, log_lower, log_upper) {
  list(
    estimate = exp(log_estimate),
    lower    = exp(log_lower),
    upper    = exp(log_upper)
  )
}

# ── Soft validation ───────────────────────────────────────────────────────────

#' Validate and cross-check reported statistics
#'
#' Performs soft (non-blocking) consistency checks and returns a character
#' vector of warning messages. An empty character vector means no issues were
#' detected. All messages are safe to display as UI callouts without
#' interrupting the analysis workflow.
#'
#' Checks performed:
#' * SE > |estimate| (implausible precision)
#' * CI-derived SE vs reported SE (>20% discrepancy)
#' * p-value implied by estimate+SE vs reported p-value (>20% discrepancy)
#' * SE implausibly small for reported sample size
#'
#' @param estimate Numeric. Point estimate.
#' @param se Numeric. Standard error.
#' @param lower Numeric or `NULL`. Lower CI bound.
#' @param upper Numeric or `NULL`. Upper CI bound.
#' @param pvalue Numeric or `NULL`. Reported p-value.
#' @param n Integer or `NULL`. Total sample size.
#' @param conf_level Numeric. Confidence level for CI checks. Default `0.95`.
#'
#' @return Character vector of warning messages (length 0 if no issues).
#' @export
#'
#' @examples
#' validate_inputs(estimate = 2.5, se = 0.66, lower = 1.2, upper = 3.8, pvalue = 0.03)
validate_inputs <- function(estimate, se, lower = NULL, upper = NULL,
                            pvalue = NULL, n = NULL, conf_level = 0.95) {
  warnings <- character(0)

  # SE > |estimate|
  if (is.numeric(se) && is.numeric(estimate) && se > abs(estimate)) {
    warnings <- c(
      warnings,
      "SE exceeds |estimate|: reported precision may be implausibly low."
    )
  }

  # CI consistency check
  if (!is.null(lower) && !is.null(upper) && is.numeric(lower) && is.numeric(upper)) {
    tryCatch(
      {
        se_from_ci    <- ci_to_se(lower, upper, conf_level = conf_level)
        discrepancy   <- abs(se_from_ci - se) / se
        if (discrepancy > 0.20) {
          warnings <- c(
            warnings,
            sprintf(
              paste0(
                "SE derived from CI (%.4f) differs from reported SE (%.4f) ",
                "by %.0f%%. Check CI level or rounding."
              ),
              se_from_ci, se, 100 * discrepancy
            )
          )
        }
      },
      error = function(e) NULL
    )
  }

  # p-value consistency
  if (!is.null(pvalue) && is.numeric(pvalue) && pvalue > 0 && pvalue < 1 &&
      is.numeric(estimate) && is.numeric(se) && se > 0) {
    z_reported  <- abs(estimate) / se
    p_implied   <- 2 * pnorm(-abs(z_reported))
    p_discrepancy <- abs(p_implied - pvalue) / pvalue
    if (p_discrepancy > 0.20) {
      warnings <- c(
        warnings,
        sprintf(
          paste0(
            "Reported p-value (%.3f) is inconsistent with estimate and SE ",
            "(implied p = %.3f, %.0f%% discrepancy)."
          ),
          pvalue, p_implied, 100 * p_discrepancy
        )
      )
    }
  }

  # Implausibly small n for SE
  if (!is.null(n) && is.numeric(n) && n >= 2 && is.numeric(se) && se > 0) {
    # Rough lower bound: SE >= 1/sqrt(n) only if SD ~ 1 and all in one group
    # Use a very generous threshold (SE < 1/sqrt(n) / 10)
    min_plausible_se <- 1 / sqrt(n) / 10
    if (se < min_plausible_se) {
      warnings <- c(
        warnings,
        sprintf(
          "SE (%.4f) seems implausibly small for n = %d. Check units and scale.",
          se, as.integer(n)
        )
      )
    }
  }

  warnings
}
