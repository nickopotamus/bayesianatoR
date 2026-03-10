#' Core Bayesian Updating Functions
#'
#' @description
#' Conjugate (analytical) Bayesian updating for eight common clinical outcome
#' types. Each function accepts a likelihood summary (from the paper) and an
#' `nv_prior` object (from [make_prior()]), and returns a fully populated
#' `nv_result` S3 object.
#'
#' @section Normal–Normal conjugate update:
#' For a Normal likelihood with known variance, the conjugate prior is Normal.
#' Given prior N(\eqn{\mu_0}, \eqn{\tau_0^2}) and data summary
#' \eqn{\hat{\theta}} with standard error \eqn{\sigma}:
#' \deqn{
#'   \frac{1}{\tau_1^2} = \frac{1}{\tau_0^2} + \frac{1}{\sigma^2}
#'   \qquad
#'   \mu_1 = \tau_1^2 \left(\frac{\mu_0}{\tau_0^2} + \frac{\hat\theta}{\sigma^2}\right)
#' }
#'
#' @section Ratio measures (OR, RR, HR):
#' For ratio measures the Normal–Normal update is applied on the *log* scale.
#' Posterior summaries are back-transformed for display.
#'
#' @section Beta–Binomial update:
#' For proportion outcomes, the conjugate prior is Beta(\eqn{\alpha_0},
#' \eqn{\beta_0}). Given \eqn{x} events in \eqn{n} trials:
#' \deqn{
#'   \text{Posterior} = \text{Beta}(\alpha_0 + x,\; \beta_0 + n - x)
#' }
#'
#' @name conjugate_updates
NULL

# ── Internal constructor ──────────────────────────────────────────────────────

#' Build an nv_result skeleton
#'
#' @param analysis_type Character. Human-readable analysis label.
#' @param prior An `nv_prior` object.
#' @param likelihood Named list of sufficient statistics from the paper.
#' @param posterior Named list with at least `mean` and `sd`.
#' @param ci_level Numeric. Posterior CrI level (0–1). Default `0.95`.
#' @param mid Numeric. Minimum important difference for P(benefit). Default `0`.
#' @param rope Numeric vector of length 2 or `NULL`. ROPE bounds.
#' @param call The caller environment (for error attribution).
#'
#' @return An incomplete `nv_result` object (slots filled by `.finalise_result()`).
#' @keywords internal
new_nv_result <- function(analysis_type, prior, likelihood, posterior,
                          ci_level = 0.95, mid = 0, rope = NULL,
                          call = rlang::caller_env()) {
  structure(
    list(
      analysis_type          = analysis_type,
      prior                  = prior,
      likelihood             = likelihood,
      posterior              = posterior,
      credible_interval      = NULL,
      probability_of_benefit = NULL,
      bayes_factor           = NULL,
      prior_sensitivity      = NULL,
      plots                  = NULL,
      meta = list(ci_level = ci_level, mid = mid, rope = rope)
    ),
    class = c(
      paste0("nv_", gsub("[ /]", "_", analysis_type)),
      "nv_result"
    )
  )
}

# ── Normal–Normal update (internal workhorse) ─────────────────────────────────

#' Normal–Normal conjugate update
#'
#' @param estimate Numeric. Observed effect estimate.
#' @param se Numeric. Standard error of the estimate.
#' @param prior An `nv_prior` object.
#'
#' @return Named list with `mean` and `sd` of the posterior distribution.
#' @keywords internal
normal_normal_update <- function(estimate, se, prior) {
  tau0_sq  <- prior$sd^2
  sigma_sq <- se^2

  # Posterior precision = sum of precisions
  tau1_sq <- 1 / (1 / tau0_sq + 1 / sigma_sq)

  # Posterior mean = precision-weighted average
  mu1 <- tau1_sq * (prior$mean / tau0_sq + estimate / sigma_sq)

  list(mean = mu1, sd = sqrt(tau1_sq))
}

# ── Finaliser ─────────────────────────────────────────────────────────────────

#' Finalise an nv_result by computing all summary slots
#'
#' Calls [compute_posterior_summary()], [compute_bayes_factor()],
#' [compute_sensitivity()], and [make_plots()] in sequence, attaching results
#' to the supplied (incomplete) `nv_result` object.
#'
#' @param result An incomplete `nv_result` object.
#'
#' @return A fully populated `nv_result` object.
#' @keywords internal
.finalise_result <- function(result) {
  result <- compute_posterior_summary(result)
  result <- compute_bayes_factor(result)
  result <- compute_sensitivity(result)
  result$plots <- make_plots(result)
  result
}

# ── Input check helpers ───────────────────────────────────────────────────────

#' @keywords internal
.check_estimate_se <- function(estimate, se, call = rlang::caller_env()) {
  if (!is.numeric(estimate) || length(estimate) != 1) {
    rlang::abort("`estimate` must be a single numeric value.", .call = call)
  }
  if (!is.numeric(se) || length(se) != 1 || se <= 0) {
    rlang::abort("`se` must be a single positive number.", .call = call)
  }
}

#' @keywords internal
.check_prior <- function(prior, call = rlang::caller_env()) {
  if (!inherits(prior, "nv_prior")) {
    rlang::abort(
      "`prior` must be an `nv_prior` object created by `make_prior()`.",
      .call = call
    )
  }
}

# ── Two-group continuous (Normal–Normal) ──────────────────────────────────────

#' Bayesian update for a two-group continuous outcome
#'
#' Uses a Normal–Normal conjugate model. The likelihood is summarised as
#' N(\eqn{\hat\delta}, \eqn{\sigma^2}) where \eqn{\hat\delta} is the observed
#' mean difference and \eqn{\sigma} is its standard error.
#'
#' @param estimate Numeric. Observed mean difference (group 1 \eqn{-} group 2).
#' @param se Numeric. Standard error of the mean difference. Obtain from
#'   [sd_to_se()] or [ci_to_se()] if not directly reported.
#' @param prior An `nv_prior` object created by [make_prior()].
#' @param ci_level Numeric. Credible interval level (0–1). Default `0.95`.
#' @param mid Numeric. Minimum important difference for P(benefit). Default `0`.
#' @param rope Numeric vector of length 2 or `NULL`. ROPE bounds for the
#'   ROPE-based Bayes factor. Default `NULL`.
#'
#' @return An object of class `c("nv_two_group_continuous", "nv_result")`.
#' @export
#'
#' @examples
#' prior  <- make_prior("sceptical", effect_size = 0.5)
#' result <- update_two_group_continuous(estimate = 0.4, se = 0.15, prior = prior)
#' print(result)
update_two_group_continuous <- function(estimate, se, prior,
                                        ci_level = 0.95, mid = 0,
                                        rope = NULL) {
  .check_estimate_se(estimate, se)
  .check_prior(prior)

  likelihood <- list(estimate = estimate, se = se)
  posterior  <- normal_normal_update(estimate, se, prior)

  result <- new_nv_result(
    analysis_type = "two-group continuous",
    prior         = prior,
    likelihood    = likelihood,
    posterior     = posterior,
    ci_level      = ci_level,
    mid           = mid,
    rope          = rope
  )

  .finalise_result(result)
}

# ── Two-group proportions (Beta–Binomial) ─────────────────────────────────────

#' Bayesian update for a two-group proportion outcome
#'
#' Uses a Beta–Binomial conjugate model independently for each group. The
#' posterior risk difference is summarised via a normal approximation (delta
#' method) for credible interval and Bayes factor calculations.
#'
#' @param x1 Integer. Event count in group 1.
#' @param n1 Integer. Total in group 1.
#' @param x2 Integer. Event count in group 2.
#' @param n2 Integer. Total in group 2.
#' @param prior_alpha Numeric. Alpha (shape1) of the Beta prior for each group.
#'   Default `1` (uniform / Haldane-ish prior with `prior_beta = 1`).
#' @param prior_beta Numeric. Beta (shape2) of the Beta prior. Default `1`.
#' @param ci_level Numeric. Credible interval level. Default `0.95`.
#' @param mid Numeric. MID on the risk-difference scale. Default `0`.
#' @param rope Numeric vector of length 2 or `NULL`.
#'
#' @return An object of class `c("nv_two_group_proportions", "nv_result")`.
#' @export
#'
#' @examples
#' update_two_group_proportions(x1 = 30, n1 = 100, x2 = 20, n2 = 100)
update_two_group_proportions <- function(x1, n1, x2, n2,
                                         prior_alpha = 1, prior_beta = 1,
                                         ci_level = 0.95, mid = 0,
                                         rope = NULL) {
  for (nm in c("x1", "n1", "x2", "n2")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1 || val < 0 || val != round(val)) {
      rlang::abort(
        sprintf("`%s` must be a non-negative integer.", nm),
        .call = rlang::caller_env()
      )
    }
  }
  if (x1 > n1 || x2 > n2) {
    rlang::abort("Event counts cannot exceed group totals.", .call = rlang::caller_env())
  }

  # Posterior Beta parameters
  post_alpha1 <- prior_alpha + x1
  post_beta1  <- prior_beta  + (n1 - x1)
  post_alpha2 <- prior_alpha + x2
  post_beta2  <- prior_beta  + (n2 - x2)

  # Posterior means (event rates)
  mu1  <- post_alpha1 / (post_alpha1 + post_beta1)
  mu2  <- post_alpha2 / (post_alpha2 + post_beta2)

  # Posterior variances (exact Beta formula)
  var1 <- (post_alpha1 * post_beta1) /
    ((post_alpha1 + post_beta1)^2 * (post_alpha1 + post_beta1 + 1))
  var2 <- (post_alpha2 * post_beta2) /
    ((post_alpha2 + post_beta2)^2 * (post_alpha2 + post_beta2 + 1))

  prior_obj <- structure(
    list(
      type = "beta-binomial", alpha = prior_alpha, beta = prior_beta,
      mean = prior_alpha / (prior_alpha + prior_beta), sd = NA_real_,
      log_scale = FALSE,
      label = sprintf("Beta(\u03b1 = %.3g, \u03b2 = %.3g)", prior_alpha, prior_beta)
    ),
    class = "nv_prior"
  )

  likelihood <- list(
    x1 = x1, n1 = n1, x2 = x2, n2 = n2,
    p1_hat    = x1 / n1,
    p2_hat    = x2 / n2,
    risk_diff = x1 / n1 - x2 / n2,
    # Store estimate/se for downstream normal approximation
    estimate  = x1 / n1 - x2 / n2,
    se        = sqrt((x1 / n1) * (1 - x1 / n1) / n1 + (x2 / n2) * (1 - x2 / n2) / n2)
  )

  posterior <- list(
    alpha1   = post_alpha1, beta1 = post_beta1,
    alpha2   = post_alpha2, beta2 = post_beta2,
    mean     = mu1 - mu2,                  # posterior risk difference
    sd       = sqrt(var1 + var2)           # delta method / independence approx
  )

  result <- new_nv_result(
    analysis_type = "two-group proportions",
    prior         = prior_obj,
    likelihood    = likelihood,
    posterior     = posterior,
    ci_level      = ci_level,
    mid           = mid,
    rope          = rope
  )

  .finalise_result(result)
}

# ── Odds ratio / Risk ratio (log-Normal) ──────────────────────────────────────

#' Bayesian update for an odds ratio or risk ratio
#'
#' Applies a Normal–Normal conjugate update on the log scale. Input should be
#' on the natural scale; log-transformation is handled internally via
#' [log_transform_ratio()].
#'
#' @param estimate Numeric. OR or RR on the *natural* scale (must be > 0).
#' @param lower Numeric. Lower CI bound on the natural scale.
#' @param upper Numeric. Upper CI bound on the natural scale.
#' @param prior An `nv_prior` object specified on the *log* scale.
#'   Use `make_prior(..., log_scale = TRUE)`.
#' @param conf_level Numeric. Confidence level of the input CI. Default `0.95`.
#' @param ci_level Numeric. Posterior CrI level. Default `0.95`.
#' @param mid Numeric. MID on the log scale (e.g., `log(1.1)`). Default `0`
#'   (i.e., ratio = 1, no effect).
#' @param rope Numeric vector of length 2 or `NULL` (on log scale).
#' @param measure Character. `"OR"` or `"RR"`. Used for labelling only.
#'
#' @return An object of class `c("nv_odds_ratio", "nv_result")`.
#' @export
#'
#' @examples
#' prior <- make_prior("sceptical", effect_size = log(1.5), log_scale = TRUE)
#' update_odds_ratio(estimate = 1.5, lower = 1.1, upper = 2.1, prior = prior)
update_odds_ratio <- function(estimate, lower, upper, prior,
                               conf_level = 0.95, ci_level = 0.95,
                               mid = 0, rope = NULL,
                               measure = c("OR", "RR")) {
  measure <- rlang::arg_match(measure)
  .check_prior(prior)

  log_inputs <- log_transform_ratio(estimate, lower, upper, conf_level = conf_level)

  likelihood <- c(
    list(estimate = estimate, lower = lower, upper = upper, measure = measure),
    log_inputs
  )

  posterior_log <- normal_normal_update(
    log_inputs$log_estimate, log_inputs$log_se, prior
  )

  z_ci <- qnorm(1 - (1 - ci_level) / 2)
  posterior <- c(
    posterior_log,
    list(
      estimate_natural = exp(posterior_log$mean),
      lower_natural    = exp(posterior_log$mean - z_ci * posterior_log$sd),
      upper_natural    = exp(posterior_log$mean + z_ci * posterior_log$sd)
    )
  )

  result <- new_nv_result(
    analysis_type = tolower(measure),
    prior         = prior,
    likelihood    = likelihood,
    posterior     = posterior,
    ci_level      = ci_level,
    mid           = mid,
    rope          = rope
  )

  .finalise_result(result)
}

# ── Paired continuous ─────────────────────────────────────────────────────────

#' Bayesian update for a paired continuous outcome
#'
#' Mathematically identical to [update_two_group_continuous()], but accepts
#' the mean of within-pair differences and its standard error. Labelled
#' separately for clarity in output.
#'
#' @param estimate Numeric. Mean of within-pair differences.
#' @param se Numeric. Standard error of the mean difference.
#' @param prior An `nv_prior` object.
#' @param ci_level Numeric. Default `0.95`.
#' @param mid Numeric. MID for P(benefit). Default `0`.
#' @param rope Numeric vector of length 2 or `NULL`.
#'
#' @return An object of class `c("nv_paired_continuous", "nv_result")`.
#' @export
#'
#' @examples
#' prior <- make_prior("noninformative")
#' update_paired_continuous(estimate = -2.1, se = 0.8, prior = prior)
update_paired_continuous <- function(estimate, se, prior,
                                     ci_level = 0.95, mid = 0, rope = NULL) {
  .check_estimate_se(estimate, se)
  .check_prior(prior)

  likelihood <- list(estimate = estimate, se = se)
  posterior  <- normal_normal_update(estimate, se, prior)

  result <- new_nv_result(
    analysis_type = "paired continuous",
    prior         = prior,
    likelihood    = likelihood,
    posterior     = posterior,
    ci_level      = ci_level,
    mid           = mid,
    rope          = rope
  )

  .finalise_result(result)
}

# ── Single proportion vs null ─────────────────────────────────────────────────

#' Bayesian update for a single proportion
#'
#' Uses a Beta–Binomial conjugate model. The null proportion for Bayes factor
#' computation is set via `null_p` (default 0.5 for a two-sided null).
#'
#' @param x Integer. Number of events.
#' @param n Integer. Number of trials.
#' @param prior_alpha Numeric. Alpha parameter of Beta prior. Default `1`.
#' @param prior_beta Numeric. Beta parameter of Beta prior. Default `1`.
#' @param null_p Numeric. Null hypothesis proportion (0–1). Default `0.5`.
#' @param ci_level Numeric. Default `0.95`.
#' @param mid Numeric. MID probability threshold for P(benefit). Default `0.5`.
#'
#' @return An object of class `c("nv_single_proportion", "nv_result")`.
#' @export
#'
#' @examples
#' update_single_proportion(x = 35, n = 50, null_p = 0.5)
update_single_proportion <- function(x, n, prior_alpha = 1, prior_beta = 1,
                                     null_p = 0.5, ci_level = 0.95,
                                     mid = 0.5) {
  for (nm in c("x", "n")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1 || val < 0 || val != round(val)) {
      rlang::abort(
        sprintf("`%s` must be a non-negative integer.", nm),
        .call = rlang::caller_env()
      )
    }
  }
  if (x > n) {
    rlang::abort("`x` cannot exceed `n`.", .call = rlang::caller_env())
  }
  if (!is.numeric(null_p) || null_p <= 0 || null_p >= 1) {
    rlang::abort("`null_p` must be strictly between 0 and 1.", .call = rlang::caller_env())
  }

  post_alpha <- prior_alpha + x
  post_beta  <- prior_beta  + (n - x)
  post_mean  <- post_alpha / (post_alpha + post_beta)
  post_var   <- (post_alpha * post_beta) /
    ((post_alpha + post_beta)^2 * (post_alpha + post_beta + 1))

  prior_obj <- structure(
    list(
      type = "beta-binomial", alpha = prior_alpha, beta = prior_beta,
      mean = prior_alpha / (prior_alpha + prior_beta),
      sd   = NA_real_, log_scale = FALSE,
      label = sprintf("Beta(\u03b1 = %.3g, \u03b2 = %.3g)", prior_alpha, prior_beta)
    ),
    class = "nv_prior"
  )

  likelihood <- list(x = x, n = n, p_hat = x / n, null_p = null_p,
                     estimate = x / n, se = sqrt(x / n * (1 - x / n) / n))

  posterior  <- list(
    alpha = post_alpha, beta = post_beta,
    mean  = post_mean,  sd  = sqrt(post_var)
  )

  result <- new_nv_result(
    analysis_type = "single proportion",
    prior         = prior_obj,
    likelihood    = likelihood,
    posterior     = posterior,
    ci_level      = ci_level,
    mid           = mid,
    rope          = NULL
  )

  .finalise_result(result)
}

# ── Correlation coefficient (Fisher z) ────────────────────────────────────────

#' Bayesian update for a Pearson correlation coefficient
#'
#' Applies the Fisher z-transform (`atanh(r)`) to achieve approximate normality,
#' performs a Normal–Normal update in z-space (SE = 1 / sqrt(n - 3)), and
#' back-transforms the posterior mean to the correlation scale.
#' The credible interval is computed on the z-scale and transformed back.
#'
#' @param r Numeric. Observed Pearson correlation (\eqn{-1 < r < 1}).
#' @param n Integer. Sample size (must be >= 4).
#' @param prior An `nv_prior` object on the *z-transformed* scale. For a
#'   non-informative prior, `make_prior("noninformative")` is appropriate.
#' @param ci_level Numeric. Default `0.95`.
#' @param mid Numeric. MID on the z scale. Default `0` (i.e., \eqn{\rho = 0}).
#' @param rope Numeric vector of length 2 or `NULL` (on z scale).
#'
#' @return An object of class `c("nv_correlation", "nv_result")`. The
#'   `posterior$r_mean` element gives the back-transformed posterior mean.
#' @export
#'
#' @examples
#' prior <- make_prior("noninformative")
#' update_correlation(r = 0.45, n = 80, prior = prior)
update_correlation <- function(r, n, prior, ci_level = 0.95, mid = 0,
                                rope = NULL) {
  if (!is.numeric(r) || length(r) != 1 || r <= -1 || r >= 1) {
    rlang::abort(
      "`r` must be a single numeric value strictly between -1 and 1.",
      .call = rlang::caller_env()
    )
  }
  if (!is.numeric(n) || length(n) != 1 || n < 4 || n != round(n)) {
    rlang::abort(
      "`n` must be an integer >= 4 (Fisher z requires at least 4 observations).",
      .call = rlang::caller_env()
    )
  }
  .check_prior(prior)

  z_obs <- atanh(r)
  se_z  <- 1 / sqrt(n - 3)

  posterior_z <- normal_normal_update(z_obs, se_z, prior)

  # CrI on z scale, then back-transform
  z_ci <- qnorm(1 - (1 - ci_level) / 2)
  r_lower <- tanh(posterior_z$mean - z_ci * posterior_z$sd)
  r_upper <- tanh(posterior_z$mean + z_ci * posterior_z$sd)

  likelihood <- list(
    r = r, z = z_obs, se_z = se_z, n = n,
    estimate = z_obs, se = se_z     # for sensitivity sweep (z scale)
  )

  posterior <- c(
    posterior_z,
    list(
      r_mean  = tanh(posterior_z$mean),
      r_lower = r_lower,
      r_upper = r_upper
    )
  )

  result <- new_nv_result(
    analysis_type = "correlation",
    prior         = prior,
    likelihood    = likelihood,
    posterior     = posterior,
    ci_level      = ci_level,
    mid           = mid,
    rope          = rope
  )

  .finalise_result(result)
}

# ── Hazard ratio (log-Normal) ─────────────────────────────────────────────────

#' Bayesian update for a hazard ratio
#'
#' Uses a Normal–Normal conjugate model on the log scale; structurally
#' identical to [update_odds_ratio()] but labelled as a hazard ratio.
#'
#' @param estimate Numeric. HR on the natural scale (must be > 0).
#' @param lower Numeric. Lower CI bound on the natural scale.
#' @param upper Numeric. Upper CI bound on the natural scale.
#' @param prior An `nv_prior` object on the *log* scale.
#' @param conf_level Numeric. Confidence level of the input CI. Default `0.95`.
#' @param ci_level Numeric. Posterior CrI level. Default `0.95`.
#' @param mid Numeric. MID on the log scale. Default `0` (HR = 1).
#' @param rope Numeric vector of length 2 or `NULL` (log scale).
#'
#' @return An object of class `c("nv_hazard_ratio", "nv_result")`.
#' @export
#'
#' @examples
#' prior <- make_prior("sceptical", effect_size = log(0.7), log_scale = TRUE)
#' update_hazard_ratio(estimate = 0.72, lower = 0.55, upper = 0.94, prior = prior)
update_hazard_ratio <- function(estimate, lower, upper, prior,
                                conf_level = 0.95, ci_level = 0.95,
                                mid = 0, rope = NULL) {
  .check_prior(prior)

  log_inputs <- log_transform_ratio(estimate, lower, upper, conf_level = conf_level)

  likelihood <- c(
    list(estimate = estimate, lower = lower, upper = upper, measure = "HR"),
    log_inputs
  )

  posterior_log <- normal_normal_update(
    log_inputs$log_estimate, log_inputs$log_se, prior
  )

  z_ci <- qnorm(1 - (1 - ci_level) / 2)
  posterior <- c(
    posterior_log,
    list(
      estimate_natural = exp(posterior_log$mean),
      lower_natural    = exp(posterior_log$mean - z_ci * posterior_log$sd),
      upper_natural    = exp(posterior_log$mean + z_ci * posterior_log$sd)
    )
  )

  result <- new_nv_result(
    analysis_type = "hazard ratio",
    prior         = prior,
    likelihood    = likelihood,
    posterior     = posterior,
    ci_level      = ci_level,
    mid           = mid,
    rope          = rope
  )

  .finalise_result(result)
}

# ── Simple linear regression coefficient ──────────────────────────────────────

#' Bayesian update for a simple linear regression coefficient
#'
#' Uses a Normal–Normal conjugate model. Input should be the OLS estimate and
#' its standard error (as typically reported in regression output tables).
#'
#' @param estimate Numeric. OLS regression coefficient.
#' @param se Numeric. Standard error of the coefficient.
#' @param prior An `nv_prior` object.
#' @param ci_level Numeric. Default `0.95`.
#' @param mid Numeric. MID for P(benefit). Default `0`.
#' @param rope Numeric vector of length 2 or `NULL`.
#'
#' @return An object of class `c("nv_regression_coefficient", "nv_result")`.
#' @export
#'
#' @examples
#' prior <- make_prior("sceptical", effect_size = 0.3)
#' update_regression_coef(estimate = 0.25, se = 0.10, prior = prior)
update_regression_coef <- function(estimate, se, prior,
                                   ci_level = 0.95, mid = 0, rope = NULL) {
  .check_estimate_se(estimate, se)
  .check_prior(prior)

  likelihood <- list(estimate = estimate, se = se)
  posterior  <- normal_normal_update(estimate, se, prior)

  result <- new_nv_result(
    analysis_type = "regression coefficient",
    prior         = prior,
    likelihood    = likelihood,
    posterior     = posterior,
    ci_level      = ci_level,
    mid           = mid,
    rope          = rope
  )

  .finalise_result(result)
}
