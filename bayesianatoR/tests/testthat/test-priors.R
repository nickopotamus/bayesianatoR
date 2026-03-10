test_that("make_prior('enthusiastic') has correct parameters", {
  p <- make_prior("enthusiastic", effect_size = 0.5, se = 0.2)
  expect_s3_class(p, "nv_prior")
  expect_equal(p$mean, 0.5)
  expect_equal(p$sd,   0.2)
  expect_equal(p$type, "enthusiastic")
  expect_false(p$log_scale)
})

test_that("make_prior('enthusiastic') errors without effect_size", {
  expect_error(make_prior("enthusiastic", se = 0.2), "effect_size")
})

test_that("make_prior('enthusiastic') errors without se", {
  expect_error(make_prior("enthusiastic", effect_size = 0.5), "`se`")
})

test_that("make_prior('sceptical') has correct parameters", {
  p <- make_prior("sceptical", effect_size = 0.5)
  expect_s3_class(p, "nv_prior")
  expect_equal(p$mean, 0)
  expect_equal(p$sd,   0.25)   # |0.5| / 2
  expect_equal(p$type, "sceptical")
})

test_that("make_prior('sceptical') works with negative effect_size", {
  p <- make_prior("sceptical", effect_size = -0.8)
  expect_equal(p$mean, 0)
  expect_equal(p$sd,   0.4)   # |-0.8| / 2
})

test_that("make_prior('sceptical') errors without effect_size", {
  expect_error(make_prior("sceptical"), "effect_size")
})

test_that("make_prior('noninformative') has wide SD", {
  p <- make_prior("noninformative")
  expect_s3_class(p, "nv_prior")
  expect_equal(p$mean, 0)
  expect_equal(p$sd,   10)
  expect_equal(p$type, "noninformative")
})

test_that("make_prior('custom') has correct parameters", {
  p <- make_prior("custom", mean = 0.3, sd = 0.15)
  expect_s3_class(p, "nv_prior")
  expect_equal(p$mean, 0.3)
  expect_equal(p$sd,   0.15)
  expect_equal(p$type, "custom")
})

test_that("make_prior('custom') errors on missing mean or sd", {
  expect_error(make_prior("custom", mean = 0.3),      "`sd`")
  expect_error(make_prior("custom", sd   = 0.15),     "`mean`")
  expect_error(make_prior("custom"),                   "`mean`")
})

test_that("make_prior('custom') errors on non-positive sd", {
  expect_error(make_prior("custom", mean = 0, sd = 0),  "`sd` must be a positive")
  expect_error(make_prior("custom", mean = 0, sd = -1), "`sd` must be a positive")
})

test_that("make_prior log_scale flag is stored correctly", {
  p <- make_prior("sceptical", effect_size = log(1.5), log_scale = TRUE)
  expect_true(p$log_scale)
})

test_that("make_prior rejects invalid type", {
  expect_error(make_prior("magic"))
})

test_that("prior_summary returns a tibble with expected columns", {
  p   <- make_prior("custom", mean = 0.3, sd = 0.15)
  tbl <- prior_summary(p)
  expect_s3_class(tbl, "tbl_df")
  expect_named(tbl, c("type", "mean", "sd", "lower_95", "upper_95",
                       "log_scale", "label"))
  expect_equal(nrow(tbl), 1L)
})

test_that("prior_summary 95% range is correct", {
  p   <- make_prior("custom", mean = 1, sd = 2)
  tbl <- prior_summary(p)
  expect_equal(tbl$lower_95, qnorm(0.025, 1, 2), tolerance = 1e-10)
  expect_equal(tbl$upper_95, qnorm(0.975, 1, 2), tolerance = 1e-10)
})

test_that("prior_summary aborts on non-prior object", {
  expect_error(prior_summary(list(mean = 0, sd = 1)), "nv_prior")
})

test_that("print.nv_prior produces output without error", {
  p <- make_prior("sceptical", effect_size = 0.5)
  expect_output(print(p), "sceptical")
  expect_output(print(p), "SD")
})
