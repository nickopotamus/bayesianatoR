#' Posterior Summary Methods
#'
#' @description
#' Functions for computing and formatting posterior summaries from `nv_result`
#' objects. These are called automatically by the conjugate update functions
#' via `.finalise_result()` and generally do not need to be called directly.
#'
#' @name posteriors
NULL

#' Compute posterior summary statistics
#'
#' Fills in the `credible_interval` and `probability_of_benefit` slots of an
#' `nv_result` object using the posterior `mean` and `sd` (Normal approximation).
#'
#' @param result An `nv_result` object with `posterior$mean`, `posterior$sd`,
#'   `meta$ci_level`, and `meta$mid` populated.
#'
#' @return The `nv_result` with `credible_interval` and `probability_of_benefit`
#'   populated.
#' @export
#'
#' @examples
#' prior  <- make_prior("sceptical", effect_size = 0.5)
#' result <- update_two_group_continuous(estimate = 0.4, se = 0.15, prior = prior)
#' result$credible_interval
#' result$probability_of_benefit
compute_posterior_summary <- function(result) {
  if (!inherits(result, "nv_result")) {
    rlang::abort(
      "`result` must be an `nv_result` object.",
      .call = rlang::caller_env()
    )
  }

  post  <- result$posterior
  alpha <- 1 - result$meta$ci_level
  mid   <- result$meta$mid

  # Credible interval (exact for Normal posterior)
  z <- qnorm(1 - alpha / 2)
  result$credible_interval <- list(
    lower = post$mean - z * post$sd,
    upper = post$mean + z * post$sd,
    level = result$meta$ci_level
  )

  # P(theta > MID)
  result$probability_of_benefit <- pnorm(
    mid,
    mean       = post$mean,
    sd         = post$sd,
    lower.tail = FALSE
  )

  result
}

#' Print method for nv_result
#'
#' @param x An `nv_result` object.
#' @param ... Ignored.
#' @export
print.nv_result <- function(x, ...) {
  cat("bayesianatoR result\n")
  cat(strrep("=", 45), "\n")
  cat(sprintf("Analysis  : %s\n\n", x$analysis_type))

  cat("Prior\n")
  cat(sprintf("  %s\n\n", x$prior$label))

  cat("Posterior\n")
  cat(sprintf("  Mean     : %.4f\n", x$posterior$mean))
  cat(sprintf("  SD       : %.4f\n", x$posterior$sd))

  if (!is.null(x$credible_interval)) {
    ci <- x$credible_interval
    cat(sprintf(
      "  %d%% CrI  : [%.4f, %.4f]\n",
      round(100 * ci$level), ci$lower, ci$upper
    ))
  }

  if (!is.null(x$probability_of_benefit)) {
    cat(sprintf(
      "\nP(benefit) P(\u03b8 > MID = %.3g) : %.3f\n",
      x$meta$mid, x$probability_of_benefit
    ))
  }

  if (!is.null(x$bayes_factor)) {
    bf <- x$bayes_factor
    cat(sprintf("\nBayes factor vs H\u2080 : %.3f  [%s]\n",
                bf$bf_null, interpret_bf(bf$bf_null)))
    if (!is.null(bf$bf_rope)) {
      cat(sprintf("Bayes factor vs ROPE : %.3f  [%s]\n",
                  bf$bf_rope, interpret_bf(bf$bf_rope)))
    }
  }

  invisible(x)
}

#' Summarise an nv_result as a one-row tibble
#'
#' @param result An `nv_result` object.
#'
#' @return A one-row [tibble::tibble] with key posterior summary statistics.
#' @export
#'
#' @examples
#' prior  <- make_prior("sceptical", effect_size = 0.5)
#' result <- update_two_group_continuous(estimate = 0.4, se = 0.15, prior = prior)
#' result_summary(result)
result_summary <- function(result) {
  if (!inherits(result, "nv_result")) {
    rlang::abort(
      "`result` must be an `nv_result` object.",
      .call = rlang::caller_env()
    )
  }

  tibble::tibble(
    analysis_type          = result$analysis_type,
    prior_type             = result$prior$type,
    prior_mean             = result$prior$mean,
    prior_sd               = result$prior$sd,
    likelihood_estimate    = result$likelihood$estimate %||% NA_real_,
    likelihood_se          = result$likelihood$se %||% NA_real_,
    posterior_mean         = result$posterior$mean,
    posterior_sd           = result$posterior$sd,
    cri_lower              = result$credible_interval$lower %||% NA_real_,
    cri_upper              = result$credible_interval$upper %||% NA_real_,
    ci_level               = result$meta$ci_level,
    mid                    = result$meta$mid,
    probability_of_benefit = result$probability_of_benefit %||% NA_real_,
    bf_null                = result$bayes_factor$bf_null %||% NA_real_,
    bf_rope                = result$bayes_factor$bf_rope %||% NA_real_
  )
}

# NULL-coalescing operator (re-exported from rlang)
`%||%` <- rlang::`%||%`
