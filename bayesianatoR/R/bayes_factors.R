#' Bayes Factor Computation
#'
#' @description
#' Functions for computing Bayes factors against a point null hypothesis and
#' against a Region of Practical Equivalence (ROPE). Uses the Savage–Dickey
#' density ratio for point null tests.
#'
#' @section Savage–Dickey density ratio:
#' For nested models where H\eqn{_0}: \eqn{\theta = 0} is a special case of
#' H\eqn{_1}: \eqn{\theta \sim \text{Normal}(\mu_0, \tau_0^2)}, the Bayes
#' factor in favour of H\eqn{_1} is:
#' \deqn{
#'   BF_{10} = \frac{p(\theta = 0 \mid \text{prior})}
#'                  {p(\theta = 0 \mid \text{posterior})}
#' }
#'
#' @section ROPE Bayes factor:
#' The ROPE BF is the ratio of prior probability mass inside the ROPE to
#' posterior probability mass inside the ROPE:
#' \deqn{
#'   BF_{\text{ROPE}} = \frac{P(\theta \in \text{ROPE} \mid \text{prior})}
#'                           {P(\theta \in \text{ROPE} \mid \text{posterior})}
#' }
#' Values > 1 indicate that data moved mass away from the ROPE (evidence
#' of a meaningful effect).
#'
#' @references
#' Wagenmakers, E.-J., Lodewyckx, T., Kuriyal, H., & Grasman, R. (2010).
#' Bayesian hypothesis testing for psychologists: A tutorial on the
#' Savage–Dickey method. *Cognitive Psychology*, 60(3), 158–189.
#' \doi{10.1016/j.cogpsych.2009.12.001}
#'
#' @name bayes_factors
NULL

#' Compute Bayes factors for an nv_result
#'
#' Populates the `bayes_factor` slot of the result using the Savage–Dickey
#' density ratio (point null) and optionally a ROPE-based BF.
#'
#' @param result An `nv_result` object with `prior` and `posterior` slots.
#'
#' @return The `nv_result` with `bayes_factor` populated: a named list with:
#'   * `bf_null` — BF\eqn{_{10}} vs H\eqn{_0}: \eqn{\theta = 0}
#'   * `bf_rope` — BF vs ROPE (or `NULL` if no ROPE was specified)
#' @export
#'
#' @examples
#' prior  <- make_prior("sceptical", effect_size = 0.5)
#' result <- update_two_group_continuous(estimate = 0.4, se = 0.15, prior = prior)
#' result$bayes_factor
compute_bayes_factor <- function(result) {
  if (!inherits(result, "nv_result")) {
    rlang::abort(
      "`result` must be an `nv_result` object.",
      .call = rlang::caller_env()
    )
  }

  prior <- result$prior
  post  <- result$posterior
  null  <- 0L  # point null (0 on analysis scale; log(1) = 0 for ratio measures)

  # Guard: if prior or posterior SD is NA (e.g. beta-binomial priors without
  # a Normal representation), skip BF computation gracefully.
  if (is.na(prior$sd) || is.na(post$sd)) {
    result$bayes_factor <- list(bf_null = NA_real_, bf_rope = NULL)
    return(result)
  }

  # Savage–Dickey ratio: BF10 = prior_density(null) / posterior_density(null)
  prior_at_null <- dnorm(null, mean = prior$mean, sd = prior$sd)
  post_at_null  <- dnorm(null, mean = post$mean,  sd = post$sd)

  bf_null <- if (post_at_null > 0) prior_at_null / post_at_null else Inf

  # ROPE BF
  bf_rope <- NULL
  rope    <- result$meta$rope
  if (!is.null(rope) && length(rope) == 2 && !any(is.na(rope))) {
    p_prior_rope <- diff(pnorm(rope, mean = prior$mean, sd = prior$sd))
    p_post_rope  <- diff(pnorm(rope, mean = post$mean,  sd = post$sd))
    bf_rope <- if (p_post_rope > 0) p_prior_rope / p_post_rope else Inf
  }

  result$bayes_factor <- list(bf_null = bf_null, bf_rope = bf_rope)
  result
}

#' Interpret a Bayes factor using the Jeffreys–Lee–Wagenmakers scale
#'
#' Returns a verbal interpretation of a BF\eqn{_{10}} value. If BF < 1, the
#' interpretation is given in terms of evidence for H\eqn{_0}.
#'
#' @param bf Numeric. Single Bayes factor value (BF\eqn{_{10}}).
#'
#' @return Character string giving a verbal interpretation and direction.
#' @export
#'
#' @examples
#' interpret_bf(3.2)
#' interpret_bf(0.08)
#' interpret_bf(125)
interpret_bf <- function(bf) {
  if (!is.numeric(bf) || length(bf) != 1) {
    rlang::abort("`bf` must be a single numeric value.", .call = rlang::caller_env())
  }
  if (is.na(bf)) return("Not computable")
  if (is.infinite(bf) && bf > 0) return("Extreme evidence in favour of H\u2081")
  if (is.infinite(bf) && bf < 0) return("Extreme evidence in favour of H\u2080")

  # Work with |BF|; determine direction
  if (bf >= 1) {
    direction <- "in favour of H\u2081"
    val       <- bf
  } else {
    direction <- "in favour of H\u2080"
    val       <- 1 / bf
  }

  label <- dplyr::case_when(
    val < 1   ~ "No evidence",
    val < 3   ~ "Anecdotal evidence",
    val < 10  ~ "Moderate evidence",
    val < 30  ~ "Strong evidence",
    val < 100 ~ "Very strong evidence",
    TRUE      ~ "Extreme evidence"
  )

  sprintf("%s %s (BF = %.2f)", label, direction, bf)
}
