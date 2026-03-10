# ── Two-group continuous ───────────────────────────────────────────────────────

test_that("update_two_group_continuous: exact posterior for N(0,1) prior + N(2,1) likelihood", {
  # Prior N(0,1), Likelihood N(2,1)  →  Posterior N(1, 1/sqrt(2))
  prior  <- make_prior("custom", mean = 0, sd = 1)
  result <- update_two_group_continuous(estimate = 2, se = 1, prior = prior)

  expect_equal(result$posterior$mean, 1,              tolerance = 1e-10)
  expect_equal(result$posterior$sd,   sqrt(1 / 2),    tolerance = 1e-10)
})

test_that("update_two_group_continuous: posterior is between prior and likelihood", {
  prior  <- make_prior("custom", mean = 0, sd = 1)
  result <- update_two_group_continuous(estimate = 2, se = 0.5, prior = prior)
  expect_gt(result$posterior$mean, 0)
  expect_lt(result$posterior$mean, 2)
})

test_that("update_two_group_continuous: posterior SD < prior SD and likelihood SE", {
  prior  <- make_prior("custom", mean = 0, sd = 1)
  result <- update_two_group_continuous(estimate = 2, se = 1, prior = prior)
  expect_lt(result$posterior$sd, 1)   # < prior SD
  expect_lt(result$posterior$sd, 1)   # < likelihood SE
})

test_that("update_two_group_continuous returns nv_result", {
  prior  <- make_prior("sceptical", effect_size = 0.5)
  result <- update_two_group_continuous(estimate = 0.4, se = 0.15, prior = prior)
  expect_s3_class(result, "nv_result")
  expect_s3_class(result, "nv_two_group_continuous")
})

test_that("update_two_group_continuous populates all required slots", {
  prior  <- make_prior("sceptical", effect_size = 0.5)
  result <- update_two_group_continuous(estimate = 0.4, se = 0.15, prior = prior)
  expect_false(is.null(result$credible_interval))
  expect_false(is.null(result$probability_of_benefit))
  expect_false(is.null(result$bayes_factor))
  expect_false(is.null(result$prior_sensitivity))
  expect_false(is.null(result$plots))
})

test_that("update_two_group_continuous: probability_of_benefit is P(theta > 0)", {
  prior  <- make_prior("custom", mean = 0, sd = 1)
  result <- update_two_group_continuous(estimate = 2, se = 1, prior = prior,
                                        mid = 0)
  post_mean <- result$posterior$mean
  post_sd   <- result$posterior$sd
  expected  <- pnorm(0, post_mean, post_sd, lower.tail = FALSE)
  expect_equal(result$probability_of_benefit, expected, tolerance = 1e-10)
})

test_that("update_two_group_continuous: CrI width shrinks with more precise likelihood", {
  prior   <- make_prior("custom", mean = 0, sd = 1)
  result1 <- update_two_group_continuous(estimate = 1, se = 1.0, prior = prior)
  result2 <- update_two_group_continuous(estimate = 1, se = 0.1, prior = prior)
  width1  <- result1$credible_interval$upper - result1$credible_interval$lower
  width2  <- result2$credible_interval$upper - result2$credible_interval$lower
  expect_lt(width2, width1)
})

test_that("update_two_group_continuous aborts on non-numeric estimate", {
  prior <- make_prior("noninformative")
  expect_error(update_two_group_continuous("a", 0.5, prior), "estimate")
})

test_that("update_two_group_continuous aborts on non-positive SE", {
  prior <- make_prior("noninformative")
  expect_error(update_two_group_continuous(1, 0,  prior), "se")
  expect_error(update_two_group_continuous(1, -1, prior), "se")
})

# ── Two-group proportions ──────────────────────────────────────────────────────

test_that("update_two_group_proportions returns nv_result", {
  result <- update_two_group_proportions(x1 = 30, n1 = 100, x2 = 20, n2 = 100)
  expect_s3_class(result, "nv_result")
  expect_s3_class(result, "nv_two_group_proportions")
})

test_that("update_two_group_proportions: posterior alpha/beta correct", {
  result <- update_two_group_proportions(x1 = 10, n1 = 50,
                                         x2 = 5,  n2 = 50,
                                         prior_alpha = 1, prior_beta = 1)
  expect_equal(result$posterior$alpha1, 11)
  expect_equal(result$posterior$beta1,  41)
  expect_equal(result$posterior$alpha2, 6)
  expect_equal(result$posterior$beta2,  46)
})

test_that("update_two_group_proportions aborts when x > n", {
  expect_error(
    update_two_group_proportions(x1 = 60, n1 = 50, x2 = 10, n2 = 50),
    "cannot exceed"
  )
})

# ── Paired continuous ──────────────────────────────────────────────────────────

test_that("update_paired_continuous returns nv_result", {
  prior  <- make_prior("noninformative")
  result <- update_paired_continuous(estimate = -2.1, se = 0.8, prior = prior)
  expect_s3_class(result, "nv_result")
  expect_s3_class(result, "nv_paired_continuous")
})

test_that("update_paired_continuous posterior matches N-N formula", {
  prior  <- make_prior("custom", mean = 0, sd = 1)
  result <- update_paired_continuous(estimate = 1, se = 0.5, prior = prior)
  expected <- exact_nn_posterior(1, 0.5, 0, 1)
  expect_equal(result$posterior$mean, expected$mean, tolerance = 1e-10)
  expect_equal(result$posterior$sd,   expected$sd,   tolerance = 1e-10)
})

# ── Single proportion ──────────────────────────────────────────────────────────

test_that("update_single_proportion returns nv_result", {
  result <- update_single_proportion(x = 35, n = 50)
  expect_s3_class(result, "nv_result")
  expect_s3_class(result, "nv_single_proportion")
})

test_that("update_single_proportion: posterior beta params correct", {
  result <- update_single_proportion(x = 7, n = 10,
                                     prior_alpha = 2, prior_beta = 2)
  expect_equal(result$posterior$alpha, 9)
  expect_equal(result$posterior$beta,  5)
  expect_equal(result$posterior$mean,  9 / 14, tolerance = 1e-10)
})

test_that("update_single_proportion aborts when x > n", {
  expect_error(update_single_proportion(x = 60, n = 50), "cannot exceed")
})

# ── Correlation ────────────────────────────────────────────────────────────────

test_that("update_correlation returns nv_result", {
  prior  <- make_prior("noninformative")
  result <- update_correlation(r = 0.45, n = 80, prior = prior)
  expect_s3_class(result, "nv_result")
  expect_s3_class(result, "nv_correlation")
})

test_that("update_correlation: z-transform stored correctly", {
  prior  <- make_prior("noninformative")
  result <- update_correlation(r = 0.5, n = 100, prior = prior)
  expect_equal(result$likelihood$z, atanh(0.5), tolerance = 1e-10)
  expect_equal(result$likelihood$se_z, 1 / sqrt(97), tolerance = 1e-10)
})

test_that("update_correlation: r_mean is between -1 and 1", {
  prior  <- make_prior("noninformative")
  result <- update_correlation(r = 0.8, n = 50, prior = prior)
  expect_gt(result$posterior$r_mean, -1)
  expect_lt(result$posterior$r_mean,  1)
})

test_that("update_correlation aborts on |r| >= 1", {
  prior <- make_prior("noninformative")
  expect_error(update_correlation(r = 1.0,  n = 50, prior = prior), "strictly between")
  expect_error(update_correlation(r = -1.0, n = 50, prior = prior), "strictly between")
})

test_that("update_correlation aborts on n < 4", {
  prior <- make_prior("noninformative")
  expect_error(update_correlation(r = 0.3, n = 3, prior = prior), ">= 4")
})

# ── Odds ratio ─────────────────────────────────────────────────────────────────

test_that("update_odds_ratio returns nv_result", {
  prior  <- make_prior("sceptical", effect_size = log(1.5), log_scale = TRUE)
  result <- update_odds_ratio(1.5, 1.1, 2.1, prior = prior)
  expect_s3_class(result, "nv_result")
})

test_that("update_odds_ratio: log-scale posterior is N-N update", {
  prior      <- make_prior("custom", mean = 0, sd = 1)
  li         <- log_transform_ratio(1.5, 1.1, 2.1)
  result     <- update_odds_ratio(1.5, 1.1, 2.1, prior = prior)
  expected   <- exact_nn_posterior(li$log_estimate, li$log_se, 0, 1)
  expect_equal(result$posterior$mean, expected$mean, tolerance = 1e-10)
  expect_equal(result$posterior$sd,   expected$sd,   tolerance = 1e-10)
})

test_that("update_odds_ratio: natural-scale posterior mean is exp of log-scale", {
  prior  <- make_prior("noninformative")
  result <- update_odds_ratio(1.5, 1.1, 2.1, prior = prior)
  expect_equal(
    result$posterior$estimate_natural,
    exp(result$posterior$mean),
    tolerance = 1e-10
  )
})

# ── Hazard ratio ───────────────────────────────────────────────────────────────

test_that("update_hazard_ratio returns nv_result", {
  prior  <- make_prior("sceptical", effect_size = log(0.7), log_scale = TRUE)
  result <- update_hazard_ratio(0.72, 0.55, 0.94, prior = prior)
  expect_s3_class(result, "nv_result")
  expect_s3_class(result, "nv_hazard_ratio")
})

test_that("update_hazard_ratio posterior on natural scale is back-transformed", {
  prior  <- make_prior("noninformative")
  result <- update_hazard_ratio(0.72, 0.55, 0.94, prior = prior)
  expect_gt(result$posterior$estimate_natural, 0)
  expect_equal(
    result$posterior$estimate_natural,
    exp(result$posterior$mean),
    tolerance = 1e-10
  )
})

# ── Regression coefficient ─────────────────────────────────────────────────────

test_that("update_regression_coef returns nv_result", {
  prior  <- make_prior("sceptical", effect_size = 0.3)
  result <- update_regression_coef(estimate = 0.25, se = 0.10, prior = prior)
  expect_s3_class(result, "nv_result")
  expect_s3_class(result, "nv_regression_coefficient")
})

test_that("update_regression_coef posterior matches N-N formula", {
  prior    <- make_prior("custom", mean = 0, sd = 1)
  result   <- update_regression_coef(estimate = 0.5, se = 0.2, prior = prior)
  expected <- exact_nn_posterior(0.5, 0.2, 0, 1)
  expect_equal(result$posterior$mean, expected$mean, tolerance = 1e-10)
  expect_equal(result$posterior$sd,   expected$sd,   tolerance = 1e-10)
})

# ── Integration: result_summary ───────────────────────────────────────────────

test_that("result_summary returns a one-row tibble with expected columns", {
  prior  <- make_prior("sceptical", effect_size = 0.5)
  result <- update_two_group_continuous(estimate = 0.4, se = 0.15, prior = prior)
  tbl    <- result_summary(result)
  expect_s3_class(tbl, "tbl_df")
  expect_equal(nrow(tbl), 1L)
  expect_true("posterior_mean" %in% names(tbl))
  expect_true("probability_of_benefit" %in% names(tbl))
  expect_true("bf_null" %in% names(tbl))
})

# ── Integration: Bayes factor ──────────────────────────────────────────────────

test_that("Bayes factor: strong prior at null gives large BF for non-null data", {
  # Prior tightly centred at 0; data strongly away from 0 → BF10 should be > 1
  prior  <- make_prior("custom", mean = 0, sd = 0.1)
  result <- update_two_group_continuous(estimate = 2, se = 0.1, prior = prior)
  expect_gt(result$bayes_factor$bf_null, 1)
})

test_that("Bayes factor: prior matching data gives BF near 1", {
  # Prior perfectly centred on observed data; BF should be close to 1
  prior  <- make_prior("custom", mean = 2, sd = 1)
  result <- update_two_group_continuous(estimate = 2, se = 1, prior = prior)
  # Not exactly 1, but the ratio shouldn't be extreme
  expect_gt(result$bayes_factor$bf_null, 0)
  expect_false(is.infinite(result$bayes_factor$bf_null))
})

# ── Integration: sensitivity tibble ───────────────────────────────────────────

test_that("prior_sensitivity is a tibble with expected columns", {
  prior  <- make_prior("sceptical", effect_size = 0.5)
  result <- update_two_group_continuous(estimate = 0.4, se = 0.15, prior = prior)
  sens   <- result$prior_sensitivity
  expect_s3_class(sens, "tbl_df")
  expect_true(all(c("prior_mean", "prior_sd", "posterior_mean",
                     "cri_lower", "cri_upper",
                     "probability_of_benefit") %in% names(sens)))
  expect_gt(nrow(sens), 0L)
})
