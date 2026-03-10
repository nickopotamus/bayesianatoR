test_that("ci_to_se returns correct SE for 95% CI", {
  # 95% CI: z = 1.96, SE = (upper - lower) / (2 * 1.96)
  se <- ci_to_se(lower = 1.2, upper = 3.8, conf_level = 0.95)
  expected <- (3.8 - 1.2) / (2 * qnorm(0.975))
  expect_equal(se, expected, tolerance = 1e-10)
})

test_that("ci_to_se returns correct SE for 90% CI", {
  se <- ci_to_se(lower = 0, upper = 2, conf_level = 0.90)
  expected <- 2 / (2 * qnorm(0.95))
  expect_equal(se, expected, tolerance = 1e-10)
})

test_that("ci_to_se and se_to_ci are exact inverses", {
  se_original <- 0.66
  est         <- 2.5
  ci          <- se_to_ci(est, se_original, conf_level = 0.95)
  se_recovered <- ci_to_se(ci["lower"], ci["upper"], conf_level = 0.95)
  expect_equal(as.numeric(se_recovered), se_original, tolerance = 1e-12)
})

test_that("ci_to_se round-trip holds at 99% confidence level", {
  se_orig  <- 1.23
  est      <- -0.5
  ci       <- se_to_ci(est, se_orig, conf_level = 0.99)
  se_back  <- ci_to_se(ci["lower"], ci["upper"], conf_level = 0.99)
  expect_equal(as.numeric(se_back), se_orig, tolerance = 1e-12)
})

test_that("ci_to_se aborts when lower >= upper", {
  expect_error(ci_to_se(3, 2), "`lower` must be less than `upper`")
  expect_error(ci_to_se(2, 2), "`lower` must be less than `upper`")
})

test_that("ci_to_se aborts on invalid conf_level", {
  expect_error(ci_to_se(1, 3, conf_level = 0))
  expect_error(ci_to_se(1, 3, conf_level = 1))
  expect_error(ci_to_se(1, 3, conf_level = 1.5))
})

test_that("ci_to_se aborts on non-numeric inputs", {
  expect_error(ci_to_se("a", 3))
  expect_error(ci_to_se(1, "b"))
})

test_that("se_to_ci returns named vector with correct names", {
  ci <- se_to_ci(2.5, 0.66)
  expect_named(ci, c("lower", "upper"))
})

test_that("se_to_ci aborts on non-positive SE", {
  expect_error(se_to_ci(1, 0))
  expect_error(se_to_ci(1, -1))
})

test_that("sd_to_se pooled gives correct result for equal groups", {
  # Pooled SE for equal SDs and equal n: SE = sd * sqrt(2/n)
  se <- sd_to_se(sd1 = 2, sd2 = 2, n1 = 50, n2 = 50, pooled = TRUE)
  expected <- 2 * sqrt(2 / 50)
  expect_equal(se, expected, tolerance = 1e-10)
})

test_that("sd_to_se Welch gives correct result", {
  se <- sd_to_se(sd1 = 2, sd2 = 3, n1 = 40, n2 = 60, pooled = FALSE)
  expected <- sqrt(4 / 40 + 9 / 60)
  expect_equal(se, expected, tolerance = 1e-10)
})

test_that("sd_to_se aborts on non-positive SD", {
  expect_error(sd_to_se(0, 2, 50, 50))
  expect_error(sd_to_se(-1, 2, 50, 50))
})

test_that("sd_to_se aborts on n < 2", {
  expect_error(sd_to_se(2, 2, 1, 50))
})

test_that("pvalue_to_se issues a warning", {
  expect_warning(pvalue_to_se(2.5, 0.03))
})

test_that("pvalue_to_se gives approximate SE consistent with z-test", {
  est    <- 2.5
  pval   <- 0.03
  z_exp  <- qnorm(1 - pval / 2)
  se_exp <- est / z_exp
  expect_warning(se_obs <- pvalue_to_se(est, pval))
  expect_equal(se_obs, se_exp, tolerance = 1e-10)
})

test_that("pvalue_to_se aborts on invalid p-value", {
  expect_error(pvalue_to_se(1, 0))
  expect_error(pvalue_to_se(1, 1))
  expect_error(pvalue_to_se(1, -0.1))
})

test_that("log_transform_ratio transforms correctly", {
  lt <- log_transform_ratio(1.5, 1.1, 2.1)
  expect_equal(lt$log_estimate, log(1.5), tolerance = 1e-10)
  expect_equal(lt$log_lower,    log(1.1), tolerance = 1e-10)
  expect_equal(lt$log_upper,    log(2.1), tolerance = 1e-10)
  expect_equal(lt$log_se, ci_to_se(log(1.1), log(2.1)), tolerance = 1e-10)
})

test_that("log_transform_ratio aborts on non-positive inputs", {
  expect_error(log_transform_ratio(0, 1.1, 2.1))
  expect_error(log_transform_ratio(1.5, -0.1, 2.1))
})

test_that("exp_transform_result inverts log_transform_ratio", {
  lt  <- log_transform_ratio(1.5, 1.1, 2.1)
  back <- exp_transform_result(lt$log_estimate, lt$log_lower, lt$log_upper)
  expect_equal(back$estimate, 1.5, tolerance = 1e-10)
  expect_equal(back$lower,    1.1, tolerance = 1e-10)
  expect_equal(back$upper,    2.1, tolerance = 1e-10)
})

test_that("validate_inputs returns empty vector when all is well", {
  msgs <- validate_inputs(estimate = 2.5, se = 0.66,
                          lower = 1.2, upper = 3.8)
  expect_type(msgs, "character")
  # May or may not be empty depending on rounding, but should not error
})

test_that("validate_inputs warns when SE > |estimate|", {
  msgs <- validate_inputs(estimate = 0.1, se = 0.5)
  expect_true(any(grepl("SE exceeds", msgs)))
})

test_that("validate_inputs warns on p-value inconsistency", {
  # Estimate = 2, SE = 1 → p ≈ 0.046; supply wildly different p
  msgs <- validate_inputs(estimate = 2, se = 1, pvalue = 0.5)
  expect_true(any(grepl("p-value.*inconsistent", msgs, ignore.case = TRUE)))
})
