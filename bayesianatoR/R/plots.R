#' ggplot2 Plot Methods for bayesianatoR
#'
#' @description
#' Functions for creating publication-ready visualisations from `nv_result`
#' objects. All plots use a consistent clinical theme (high-contrast blues and
#' whites). Three standard plots are produced for every analysis:
#'
#' 1. **Prior–posterior overlay** (`plot_prior_posterior`): density curves for
#'    the prior, likelihood, and posterior on a common x-axis.
#' 2. **Sensitivity plot** (`plot_sensitivity`): posterior mean and 95% CrI
#'    across the prior mean spectrum (median over the prior SD grid).
#' 3. **Probability of benefit** (`plot_probability`): posterior density shaded
#'    by benefit / no-benefit regions relative to the MID.
#'
#' @name plots
NULL

# ── Clinical theme ─────────────────────────────────────────────────────────────

#' bayesianatoR ggplot2 theme
#'
#' A minimal, high-contrast clinical theme suitable for embedding in reports
#' and Shiny applications.
#'
#' @return A [ggplot2::theme] object.
#' @export
#'
#' @examples
#' ggplot2::ggplot() + theme_bayesian()
theme_bayesian <- function() {
  ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 14,
                                               colour = "#1a3a5c"),
      plot.subtitle    = ggplot2::element_text(size = 10, colour = "#4a6a8c"),
      plot.caption     = ggplot2::element_text(size = 8,  colour = "#7a9abc"),
      axis.title       = ggplot2::element_text(colour = "#1a3a5c", size = 11),
      axis.text        = ggplot2::element_text(colour = "#2c4a6e"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#e8f0f7", linewidth = 0.4),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.position  = "bottom",
      legend.title     = ggplot2::element_blank(),
      legend.text      = ggplot2::element_text(colour = "#1a3a5c")
    )
}

# ── Internal colour palette ────────────────────────────────────────────────────

.bnr_colours <- list(
  prior      = "#a8c4e0",
  posterior  = "#1a6faf",
  likelihood = "#e07b39",
  null       = "#c0392b",
  rope       = "#27ae60",
  benefit    = "#2980b9",
  harm       = "#c0392b",
  neutral    = "#bdc3c7"
)

# ── Main dispatcher ────────────────────────────────────────────────────────────

#' Create all standard plots for an nv_result
#'
#' Convenience wrapper that calls [plot_prior_posterior()],
#' [plot_sensitivity()], and [plot_probability()] and returns them as a named
#' list.
#'
#' @param result An `nv_result` object.
#'
#' @return Named list with three [ggplot2::ggplot] objects:
#'   `prior_posterior`, `sensitivity`, `probability`.
#' @export
#'
#' @examples
#' prior  <- make_prior("sceptical", effect_size = 0.5)
#' result <- update_two_group_continuous(estimate = 0.4, se = 0.15, prior = prior)
#' plots  <- make_plots(result)
#' plots$prior_posterior
make_plots <- function(result) {
  if (!inherits(result, "nv_result")) {
    rlang::abort(
      "`result` must be an `nv_result` object.",
      .call = rlang::caller_env()
    )
  }

  list(
    prior_posterior = plot_prior_posterior(result),
    sensitivity     = plot_sensitivity(result),
    probability     = plot_probability(result)
  )
}

# ── Prior–posterior overlay ────────────────────────────────────────────────────

#' Plot prior, likelihood, and posterior distributions
#'
#' Overlays density curves for the prior, likelihood (if available), and
#' posterior. The posterior CrI is shaded. A dashed vertical line marks the
#' null value (0).
#'
#' @param result An `nv_result` object.
#' @param n_points Integer. Number of points for density evaluation. Default `500`.
#'
#' @return A [ggplot2::ggplot] object.
#' @export
#'
#' @examples
#' prior  <- make_prior("sceptical", effect_size = 0.5)
#' result <- update_two_group_continuous(estimate = 0.4, se = 0.15, prior = prior)
#' plot_prior_posterior(result)
plot_prior_posterior <- function(result, n_points = 500L) {
  if (!inherits(result, "nv_result")) {
    rlang::abort(
      "`result` must be an `nv_result` object.",
      .call = rlang::caller_env()
    )
  }

  prior <- result$prior
  post  <- result$posterior
  lik   <- result$likelihood

  # Defensive: if prior or posterior SD is NA, return empty plot
  if (is.na(prior$sd) || is.na(post$sd)) {
    return(.empty_plot("Density plot not available for this prior type."))
  }

  # x range: cover ±4 SD of prior and posterior
  x_lo <- min(prior$mean - 4 * prior$sd, post$mean - 4 * post$sd)
  x_hi <- max(prior$mean + 4 * prior$sd, post$mean + 4 * post$sd)

  # Extend to include likelihood if available
  if (!is.null(lik$estimate) && !is.null(lik$se) &&
      !is.na(lik$estimate) && !is.na(lik$se)) {
    x_lo <- min(x_lo, lik$estimate - 4 * lik$se)
    x_hi <- max(x_hi, lik$estimate + 4 * lik$se)
  }

  x <- seq(x_lo, x_hi, length.out = n_points)

  dist_df <- dplyr::bind_rows(
    tibble::tibble(x = x,
                   density      = dnorm(x, prior$mean, prior$sd),
                   Distribution = "Prior"),
    tibble::tibble(x = x,
                   density      = dnorm(x, post$mean, post$sd),
                   Distribution = "Posterior")
  )

  if (!is.null(lik$estimate) && !is.null(lik$se) &&
      !is.na(lik$estimate) && !is.na(lik$se)) {
    dist_df <- dplyr::bind_rows(
      dist_df,
      tibble::tibble(x = x,
                     density      = dnorm(x, lik$estimate, lik$se),
                     Distribution = "Likelihood")
    )
  }

  colours <- c(
    Prior       = .bnr_colours$prior,
    Likelihood  = .bnr_colours$likelihood,
    Posterior   = .bnr_colours$posterior
  )

  ci       <- result$credible_interval
  subtitle <- sprintf(
    "Posterior: N(\u03bc = %.3g, \u03c3 = %.3g)",
    post$mean, post$sd
  )
  if (!is.null(ci)) {
    subtitle <- paste0(
      subtitle,
      sprintf("  |  %d%% CrI [%.3g, %.3g]",
              round(100 * ci$level), ci$lower, ci$upper)
    )
  }

  p <- ggplot2::ggplot(
    dist_df,
    ggplot2::aes(x = x, y = density,
                 colour = Distribution, fill = Distribution)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_area(
      data  = dplyr::filter(dist_df, Distribution == "Posterior"),
      alpha = 0.20
    ) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed",
      colour = .bnr_colours$null, linewidth = 0.7
    ) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::scale_fill_manual(values   = colours) +
    ggplot2::labs(
      title    = paste("Bayesian update:", result$analysis_type),
      subtitle = subtitle,
      x        = "Effect",
      y        = "Density"
    ) +
    theme_bayesian()

  # Shade CrI region
  if (!is.null(ci)) {
    p <- p + ggplot2::annotate(
      "rect",
      xmin = ci$lower, xmax = ci$upper,
      ymin = -Inf,      ymax = Inf,
      fill = .bnr_colours$posterior, alpha = 0.08
    )
  }

  p
}

# ── Sensitivity plot ────────────────────────────────────────────────────────────

#' Plot posterior sensitivity to prior choice
#'
#' Shows how the posterior mean and 95% CrI shift across the prior mean
#' spectrum (averaging over the prior SD grid). Useful for communicating
#' robustness of conclusions.
#'
#' @param result An `nv_result` object with `prior_sensitivity` populated.
#'
#' @return A [ggplot2::ggplot] object, or an empty placeholder plot if
#'   sensitivity data are unavailable.
#' @export
#'
#' @examples
#' prior  <- make_prior("sceptical", effect_size = 0.5)
#' result <- update_two_group_continuous(estimate = 0.4, se = 0.15, prior = prior)
#' plot_sensitivity(result)
plot_sensitivity <- function(result) {
  if (!inherits(result, "nv_result")) {
    rlang::abort(
      "`result` must be an `nv_result` object.",
      .call = rlang::caller_env()
    )
  }

  sens <- result$prior_sensitivity

  if (is.null(sens) || nrow(sens) == 0) {
    return(.empty_plot("Sensitivity analysis not available for this outcome type."))
  }

  # Median CrI/mean across the prior SD dimension
  sens_summary <- sens |>
    dplyr::group_by(prior_mean) |>
    dplyr::summarise(
      post_mean = stats::median(posterior_mean),
      cri_lower = stats::median(cri_lower),
      cri_upper = stats::median(cri_upper),
      .groups   = "drop"
    )

  # Mark the actual prior used
  actual_prior_mean <- result$prior$mean

  ggplot2::ggplot(sens_summary, ggplot2::aes(x = prior_mean)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = cri_lower, ymax = cri_upper),
      fill = .bnr_colours$posterior, alpha = 0.20
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = post_mean),
      colour = .bnr_colours$posterior, linewidth = 1
    ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed",
      colour = .bnr_colours$null, linewidth = 0.7
    ) +
    ggplot2::geom_vline(
      xintercept = actual_prior_mean, linetype = "dotted",
      colour = "#4a8ec2", linewidth = 0.8
    ) +
    ggplot2::annotate(
      "text",
      x      = actual_prior_mean,
      y      = max(sens_summary$cri_upper) * 0.97,
      label  = "Selected\nprior",
      colour = "#4a8ec2", size = 3.2, hjust = -0.1
    ) +
    ggplot2::labs(
      title    = "Prior sensitivity analysis",
      subtitle = "Posterior mean (\u00b1 95% CrI) across prior mean spectrum",
      caption  = "Shaded band = range across prior SD grid. Dotted line = selected prior.",
      x        = "Prior mean",
      y        = "Posterior estimate"
    ) +
    theme_bayesian()
}

# ── Probability of benefit plot ────────────────────────────────────────────────

#' Plot the posterior probability of benefit
#'
#' Displays the posterior distribution shaded into two regions: above and
#' below the minimum important difference (MID). The P(benefit) annotation
#' is placed in the upper-right quadrant.
#'
#' @param result An `nv_result` object.
#' @param n_points Integer. Number of density evaluation points. Default `500`.
#'
#' @return A [ggplot2::ggplot] object.
#' @export
#'
#' @examples
#' prior  <- make_prior("sceptical", effect_size = 0.5)
#' result <- update_two_group_continuous(estimate = 0.4, se = 0.15, prior = prior)
#' plot_probability(result)
plot_probability <- function(result, n_points = 500L) {
  if (!inherits(result, "nv_result")) {
    rlang::abort(
      "`result` must be an `nv_result` object.",
      .call = rlang::caller_env()
    )
  }

  post <- result$posterior
  mid  <- result$meta$mid
  p_b  <- result$probability_of_benefit %||% NA_real_

  if (is.na(post$sd)) {
    return(.empty_plot("Probability plot not available for this posterior type."))
  }

  x_lo <- post$mean - 4 * post$sd
  x_hi <- post$mean + 4 * post$sd
  x    <- seq(x_lo, x_hi, length.out = n_points)
  dens <- dnorm(x, post$mean, post$sd)

  df <- tibble::tibble(
    x      = x,
    density = dens,
    region  = dplyr::if_else(x > mid, "Benefit (> MID)", "No benefit / harm (\u2264 MID)")
  )

  annotation_x <- x_hi - (x_hi - x_lo) * 0.05
  annotation_y <- max(dens) * 0.85

  ggplot2::ggplot(df, ggplot2::aes(x = x, y = density)) +
    ggplot2::geom_area(ggplot2::aes(fill = region), alpha = 0.55) +
    ggplot2::geom_line(linewidth = 1, colour = .bnr_colours$posterior) +
    ggplot2::geom_vline(
      xintercept = mid, linetype = "dashed",
      colour = "#2c4a6e", linewidth = 0.8
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Benefit (> MID)"             = .bnr_colours$benefit,
        "No benefit / harm (\u2264 MID)" = .bnr_colours$harm
      )
    ) +
    {
      if (!is.na(p_b)) {
        ggplot2::annotate(
          "text",
          x      = annotation_x,
          y      = annotation_y,
          label  = sprintf("P(benefit) =\n%.1f%%", 100 * p_b),
          colour = "#1a3a5c", fontface = "bold", size = 4.5, hjust = 1
        )
      }
    } +
    ggplot2::labs(
      title    = "Posterior probability of benefit",
      subtitle = sprintf("MID = %.3g", mid),
      caption  = "MID = minimum important difference",
      x        = "Effect",
      y        = "Density"
    ) +
    theme_bayesian()
}

# ── Internal helper ────────────────────────────────────────────────────────────

#' Return an empty placeholder ggplot with a centred message
#' @keywords internal
.empty_plot <- function(message) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0.5, y = 0.5, label = message,
      colour = "#4a6a8c", size = 5, hjust = 0.5
    ) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::theme_void() +
    theme_bayesian()
}
