# tests/testthat/setup.R
# Test helpers and shared fixtures loaded before every test file.

# ── Shared priors ──────────────────────────────────────────────────────────────

# A sceptical prior centred at the null
.prior_sceptical <- bayesianatoR::make_prior("sceptical", effect_size = 0.5)

# An enthusiastic prior centred on the expected effect
.prior_enthusiastic <- bayesianatoR::make_prior(
  "enthusiastic", effect_size = 0.5, se = 0.2
)

# A non-informative prior
.prior_noninformative <- bayesianatoR::make_prior("noninformative")

# A custom prior
.prior_custom <- bayesianatoR::make_prior("custom", mean = 0.3, sd = 0.15)

# A sceptical prior on the log scale (for ratio measures)
.prior_log_sceptical <- bayesianatoR::make_prior(
  "sceptical", effect_size = log(1.5), log_scale = TRUE
)

# ── Helper: exact Normal-Normal posterior ─────────────────────────────────────

#' Compute exact Normal-Normal posterior for test validation
exact_nn_posterior <- function(estimate, se, prior_mean, prior_sd) {
  tau0_sq  <- prior_sd^2
  sigma_sq <- se^2
  tau1_sq  <- 1 / (1 / tau0_sq + 1 / sigma_sq)
  mu1      <- tau1_sq * (prior_mean / tau0_sq + estimate / sigma_sq)
  list(mean = mu1, sd = sqrt(tau1_sq))
}
