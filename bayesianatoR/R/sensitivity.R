#' Prior Sensitivity Analysis
#'
#' @description
#' Functions for sweeping across a grid of prior specifications and examining
#' how posterior summaries change. Supports robustness assessments across both
#' prior mean and prior SD dimensions.
#'
#' The sensitivity sweep uses the Normal–Normal update and is therefore only
#' available for analysis types whose likelihood is summarised as a single
#' estimate + SE (i.e., all types except beta-binomial single-proportion).
#' For non-standard likelihoods an empty tibble is returned with a message.
#'
#' @name sensitivity
NULL

#' Compute a prior sensitivity sweep
#'
#' Evaluates the posterior mean, SD, credible interval, and P(benefit) across
#' a user-supplied (or auto-generated) grid of prior means and SDs. The result
#' is stored in the `prior_sensitivity` slot of the `nv_result`.
#'
#' @param result An `nv_result` object. The `likelihood$estimate` and
#'   `likelihood$se` slots are used to run the sweep.
#' @param prior_means Numeric vector. Prior means to evaluate. If `NULL`,
#'   a symmetric 7-point grid spanning \eqn{\pm 2 \times |\hat\theta|} is used.
#' @param prior_sds Numeric vector. Prior SDs to evaluate. If `NULL`, a
#'   5-point grid from \eqn{0.5 \times se} to \eqn{3 \times se} is used.
#' @param ci_level Numeric. Credible interval level for the sensitivity output.
#'   Default `0.95`.
#'
#' @return The `nv_result` with `prior_sensitivity` populated: a
#'   [tibble::tibble] with columns:
#'   * `prior_mean`, `prior_sd`
#'   * `posterior_mean`, `posterior_sd`
#'   * `cri_lower`, `cri_upper`
#'   * `probability_of_benefit`
#' @export
#'
#' @examples
#' prior  <- make_prior("sceptical", effect_size = 0.5)
#' result <- update_two_group_continuous(estimate = 0.4, se = 0.15, prior = prior)
#' nrow(result$prior_sensitivity)
compute_sensitivity <- function(result, prior_means = NULL, prior_sds = NULL,
                                ci_level = 0.95) {
  if (!inherits(result, "nv_result")) {
    rlang::abort(
      "`result` must be an `nv_result` object.",
      .call = rlang::caller_env()
    )
  }

  est <- result$likelihood$estimate
  se  <- result$likelihood$se

  # Non-standard likelihoods (e.g. raw beta-binomial without normal approximation)
  if (is.null(est) || is.null(se) || is.na(est) || is.na(se)) {
    rlang::warn(
      paste0(
        "Prior sensitivity sweep requires `likelihood$estimate` and `likelihood$se`. ",
        "Returning empty tibble."
      ),
      .call = rlang::caller_env()
    )
    result$prior_sensitivity <- tibble::tibble(
      prior_mean             = numeric(0),
      prior_sd               = numeric(0),
      posterior_mean         = numeric(0),
      posterior_sd           = numeric(0),
      cri_lower              = numeric(0),
      cri_upper              = numeric(0),
      probability_of_benefit = numeric(0)
    )
    return(result)
  }

  # Default grids
  if (is.null(prior_means)) {
    span <- max(abs(est), se) * 2
    prior_means <- seq(-span, span, length.out = 7)
  }
  if (is.null(prior_sds)) {
    prior_sds <- seq(0.5 * se, 3 * se, length.out = 5)
  }

  grid <- tidyr::expand_grid(prior_mean = prior_means, prior_sd = prior_sds)
  mid  <- result$meta$mid
  z    <- qnorm(1 - (1 - ci_level) / 2)

  sensitivity_tbl <- purrr::pmap_dfr(grid, function(prior_mean, prior_sd) {
    tmp_prior <- structure(
      list(
        type = "sweep", mean = prior_mean, sd = prior_sd,
        log_scale = FALSE, label = "sweep"
      ),
      class = "nv_prior"
    )
    post <- normal_normal_update(est, se, tmp_prior)

    tibble::tibble(
      prior_mean             = prior_mean,
      prior_sd               = prior_sd,
      posterior_mean         = post$mean,
      posterior_sd           = post$sd,
      cri_lower              = post$mean - z * post$sd,
      cri_upper              = post$mean + z * post$sd,
      probability_of_benefit = pnorm(mid, mean = post$mean, sd = post$sd,
                                     lower.tail = FALSE)
    )
  })

  result$prior_sensitivity <- sensitivity_tbl
  result
}

#' Extract sensitivity data from an nv_result
#'
#' Convenience accessor that returns the `prior_sensitivity` tibble or aborts
#' with an informative message if it is missing.
#'
#' @param result An `nv_result` object.
#'
#' @return A [tibble::tibble].
#' @export
#'
#' @examples
#' prior  <- make_prior("sceptical", effect_size = 0.5)
#' result <- update_two_group_continuous(estimate = 0.4, se = 0.15, prior = prior)
#' get_sensitivity(result)
get_sensitivity <- function(result) {
  if (!inherits(result, "nv_result")) {
    rlang::abort(
      "`result` must be an `nv_result` object.",
      .call = rlang::caller_env()
    )
  }
  sens <- result$prior_sensitivity
  if (is.null(sens)) {
    rlang::abort(
      paste0(
        "No sensitivity data found. ",
        "Run `compute_sensitivity()` first or ensure the update function ",
        "supports sensitivity analysis."
      ),
      .call = rlang::caller_env()
    )
  }
  sens
}
