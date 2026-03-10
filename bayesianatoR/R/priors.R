#' Prior Construction Helpers
#'
#' @description
#' Functions for constructing named prior distributions for use in Bayesian
#' updating. Supports enthusiastic, sceptical, non-informative, and custom
#' priors, all represented as Normal distributions on the analysis scale.
#'
#' For ratio measures (OR, RR, HR), priors should be specified on the *log*
#' scale by setting `log_scale = TRUE`.
#'
#' @name priors
NULL

#' Construct a prior distribution
#'
#' @description
#' Creates a prior specification as a named list of class `"nv_prior"`. Four
#' canonical clinical priors are supported:
#'
#' * `"enthusiastic"` — centred on the expected effect size (e.g., from the
#'   power calculation), with SD equal to `se`. Expresses belief that the
#'   true effect is close to the anticipated value.
#' * `"sceptical"` — centred on the null (0), with SD = `effect_size / 2`.
#'   Represents mild prior scepticism that any effect exists.
#' * `"noninformative"` — wide Normal(0, 10) on the standardised scale.
#'   Effectively non-committal; lets the data dominate.
#' * `"custom"` — user-specified `mean` and `sd`.
#'
#' For ratio measures, work on the log scale (e.g., `effect_size = log(1.5)`)
#' and set `log_scale = TRUE`.
#'
#' @param type Character. One of `"enthusiastic"`, `"sceptical"`,
#'   `"noninformative"`, or `"custom"`.
#' @param effect_size Numeric. Expected effect size. Required for
#'   `"enthusiastic"` and `"sceptical"`. On the log scale for ratio measures.
#' @param se Numeric. Standard error of the expected effect. Used as the prior
#'   SD for `"enthusiastic"`. Ignored for other types.
#' @param mean Numeric. Prior mean. Required for `type = "custom"`.
#' @param sd Numeric. Prior SD (must be > 0). Required for `type = "custom"`.
#' @param log_scale Logical. Is the prior specified on the log scale?
#'   Default `FALSE`. Set to `TRUE` for OR, RR, and HR analyses.
#'
#' @return A list of class `"nv_prior"` with elements:
#'   * `type`       — character label
#'   * `mean`       — prior mean
#'   * `sd`         — prior SD
#'   * `log_scale`  — logical
#'   * `label`      — human-readable description (used in plots)
#' @export
#'
#' @examples
#' make_prior("enthusiastic", effect_size = 0.5, se = 0.2)
#' make_prior("sceptical", effect_size = 0.5)
#' make_prior("noninformative")
#' make_prior("custom", mean = 0.3, sd = 0.15)
#' # Ratio measure (prior on log scale)
#' make_prior("sceptical", effect_size = log(1.5), log_scale = TRUE)
make_prior <- function(
    type       = c("enthusiastic", "sceptical", "noninformative", "custom"),
    effect_size = NULL,
    se         = NULL,
    mean       = NULL,
    sd         = NULL,
    log_scale  = FALSE
) {
  type <- rlang::arg_match(type)

  prior <- switch(
    type,

    enthusiastic = {
      if (is.null(effect_size)) {
        rlang::abort(
          "`effect_size` is required for `type = 'enthusiastic'`.",
          .call = rlang::caller_env()
        )
      }
      if (is.null(se)) {
        rlang::abort(
          "`se` is required for `type = 'enthusiastic'`.",
          .call = rlang::caller_env()
        )
      }
      if (!is.numeric(se) || se <= 0) {
        rlang::abort("`se` must be a positive number.", .call = rlang::caller_env())
      }
      list(
        type      = "enthusiastic",
        mean      = effect_size,
        sd        = se,
        log_scale = log_scale,
        label     = sprintf(
          "Enthusiastic: Normal(\u03bc = %.3g, \u03c3 = %.3g)",
          effect_size, se
        )
      )
    },

    sceptical = {
      if (is.null(effect_size)) {
        rlang::abort(
          "`effect_size` is required for `type = 'sceptical'`.",
          .call = rlang::caller_env()
        )
      }
      prior_sd <- abs(effect_size) / 2
      list(
        type      = "sceptical",
        mean      = 0,
        sd        = prior_sd,
        log_scale = log_scale,
        label     = sprintf(
          "Sceptical: Normal(\u03bc = 0, \u03c3 = %.3g)",
          prior_sd
        )
      )
    },

    noninformative = {
      list(
        type      = "noninformative",
        mean      = 0,
        sd        = 10,
        log_scale = log_scale,
        label     = "Non-informative: Normal(\u03bc = 0, \u03c3 = 10)"
      )
    },

    custom = {
      if (is.null(mean) || is.null(sd)) {
        rlang::abort(
          "`mean` and `sd` are both required for `type = 'custom'`.",
          .call = rlang::caller_env()
        )
      }
      if (!is.numeric(mean)) {
        rlang::abort("`mean` must be numeric.", .call = rlang::caller_env())
      }
      if (!is.numeric(sd) || sd <= 0) {
        rlang::abort("`sd` must be a positive number.", .call = rlang::caller_env())
      }
      list(
        type      = "custom",
        mean      = mean,
        sd        = sd,
        log_scale = log_scale,
        label     = sprintf(
          "Custom: Normal(\u03bc = %.3g, \u03c3 = %.3g)",
          mean, sd
        )
      )
    }
  )

  structure(prior, class = "nv_prior")
}

#' Print method for nv_prior
#'
#' @param x An `nv_prior` object.
#' @param ... Ignored.
#' @export
print.nv_prior <- function(x, ...) {
  cat("bayesianatoR prior\n")
  cat(strrep("-", 40), "\n")
  cat(sprintf("Type      : %s\n", x$type))
  cat(sprintf("Mean      : %g\n", x$mean))
  cat(sprintf("SD        : %g\n", x$sd))
  cat(sprintf(
    "95%% range : [%g, %g]\n",
    round(qnorm(0.025, x$mean, x$sd), 4),
    round(qnorm(0.975, x$mean, x$sd), 4)
  ))
  if (x$log_scale) cat("Scale     : log (ratio measure)\n")
  cat(sprintf("Label     : %s\n", x$label))
  invisible(x)
}

#' Summarise a prior as a one-row tibble
#'
#' @param prior An `nv_prior` object created by [make_prior()].
#'
#' @return A one-row [tibble::tibble] with columns `type`, `mean`, `sd`,
#'   `lower_95`, `upper_95`, `log_scale`, and `label`.
#' @export
#'
#' @examples
#' prior_summary(make_prior("sceptical", effect_size = 0.5))
prior_summary <- function(prior) {
  if (!inherits(prior, "nv_prior")) {
    rlang::abort(
      "`prior` must be an `nv_prior` object created by `make_prior()`.",
      .call = rlang::caller_env()
    )
  }

  tibble::tibble(
    type      = prior$type,
    mean      = prior$mean,
    sd        = prior$sd,
    lower_95  = qnorm(0.025, prior$mean, prior$sd),
    upper_95  = qnorm(0.975, prior$mean, prior$sd),
    log_scale = prior$log_scale,
    label     = prior$label
  )
}
