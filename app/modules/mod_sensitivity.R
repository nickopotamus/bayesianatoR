# app/modules/mod_sensitivity.R
# Prior sensitivity analysis module.
#
# Displays how the posterior mean, CrI, and P(benefit) change across a
# spectrum of prior means and SDs. Users can customise the grid bounds.

# ── UI ────────────────────────────────────────────────────────────────────────

#' Sensitivity analysis UI
#'
#' @param id Module namespace ID.
#' @export
mod_sensitivity_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(

    bslib::layout_columns(
      col_widths = c(3, 9),

      # ── Controls ──────────────────────────────────────────────────────
      bslib::card(
        bslib::card_header(
          shiny::icon("sliders"), " Sweep settings"
        ),
        bslib::card_body(
          shiny::p(
            class = "text-muted small",
            "Adjust the grid over which the prior mean and SD are swept."
          ),
          shiny::numericInput(
            ns("mean_lo"),
            label = shiny::span(
              "Prior mean: lower bound",
              bslib::tooltip(shiny::icon("circle-question"),
                             "Minimum prior mean in the sweep grid.")
            ),
            value = -1, step = 0.1
          ),
          shiny::numericInput(
            ns("mean_hi"),
            label = "Prior mean: upper bound",
            value = 1, step = 0.1
          ),
          shiny::numericInput(
            ns("mean_n"),
            label = "Grid points (mean)",
            value = 15L, min = 5L, max = 50L, step = 1L
          ),
          shiny::hr(),
          shiny::numericInput(
            ns("sd_lo"),
            label = "Prior SD: lower bound",
            value = 0.1, min = 0.01, step = 0.05
          ),
          shiny::numericInput(
            ns("sd_hi"),
            label = "Prior SD: upper bound",
            value = 2, min = 0.05, step = 0.1
          ),
          shiny::numericInput(
            ns("sd_n"),
            label = "Grid points (SD)",
            value = 5L, min = 2L, max = 20L, step = 1L
          ),
          shiny::hr(),
          shiny::actionButton(
            ns("refresh"),
            label = "Rerun sweep",
            icon  = shiny::icon("rotate"),
            class = "btn-primary w-100"
          )
        )
      ),

      # ── Plots ──────────────────────────────────────────────────────────
      shiny::tagList(
        bslib::card(
          bslib::card_header(
            shiny::icon("chart-line"),
            " Posterior mean across prior mean spectrum"
          ),
          bslib::card_body(
            shiny::plotOutput(ns("plot_sensitivity_mean"), height = "320px")
          )
        ),
        bslib::layout_columns(
          col_widths = c(6, 6),
          bslib::card(
            bslib::card_header(
              shiny::icon("arrows-left-right"), " 95% CrI width"
            ),
            bslib::card_body(
              shiny::plotOutput(ns("plot_sensitivity_cri"), height = "250px")
            )
          ),
          bslib::card(
            bslib::card_header(
              shiny::icon("percent"), " P(benefit)"
            ),
            bslib::card_body(
              shiny::plotOutput(ns("plot_sensitivity_pbenefit"), height = "250px")
            )
          )
        ),
        bslib::card(
          bslib::card_header(
            shiny::icon("table"), " Sensitivity data"
          ),
          bslib::card_body(
            shiny::div(
              class = "text-end mb-2",
              shiny::downloadButton(ns("download_sens"), "Download CSV",
                                    class = "btn-sm btn-outline-secondary")
            ),
            shiny::tableOutput(ns("sensitivity_table"))
          )
        )
      )
    )
  )
}

# ── Server ────────────────────────────────────────────────────────────────────

#' Sensitivity analysis server
#'
#' @param id Module namespace ID.
#' @param analysis_out Reactive `nv_result` (or `NULL`) from `app.R`.
#' @export
mod_sensitivity_server <- function(id, analysis_out) {
  shiny::moduleServer(id, function(input, output, session) {

    # ── Recompute sensitivity on refresh or new result ────────────────────
    sens_data <- shiny::eventReactive(
      list(analysis_out(), input$refresh),
      {
        r <- analysis_out()
        shiny::req(r)

        prior_means <- seq(input$mean_lo %||% -1, input$mean_hi %||% 1,
                           length.out = max(5L, as.integer(input$mean_n %||% 15L)))
        prior_sds   <- seq(input$sd_lo   %||% 0.1, input$sd_hi  %||% 2,
                           length.out = max(2L, as.integer(input$sd_n   %||% 5L)))

        tryCatch(
          bayesianatoR::compute_sensitivity(
            r,
            prior_means = prior_means,
            prior_sds   = prior_sds,
            ci_level    = r$meta$ci_level
          )$prior_sensitivity,
          error = function(e) {
            shiny::showNotification(
              paste("Sensitivity error:", conditionMessage(e)),
              type = "warning"
            )
            NULL
          }
        )
      },
      ignoreNULL = TRUE
    )

    # ── Posterior mean sensitivity plot ───────────────────────────────────
    output$plot_sensitivity_mean <- shiny::renderPlot({
      sens <- sens_data()
      shiny::req(sens, nrow(sens) > 0)

      sens_summary <- sens |>
        dplyr::group_by(prior_mean) |>
        dplyr::summarise(
          post_mean = stats::median(posterior_mean),
          cri_lower = stats::median(cri_lower),
          cri_upper = stats::median(cri_upper),
          .groups   = "drop"
        )

      ggplot2::ggplot(sens_summary, ggplot2::aes(x = prior_mean)) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = cri_lower, ymax = cri_upper),
          fill = "#a8c4e0", alpha = 0.30
        ) +
        ggplot2::geom_line(
          ggplot2::aes(y = post_mean),
          colour = "#1a6faf", linewidth = 1
        ) +
        ggplot2::geom_hline(
          yintercept = 0, linetype = "dashed",
          colour = "#c0392b", linewidth = 0.7
        ) +
        ggplot2::labs(
          title    = "Posterior mean (median over prior SD grid)",
          subtitle = "Ribbon = 95% CrI (median over prior SD grid)",
          x        = "Prior mean",
          y        = "Posterior estimate"
        ) +
        bayesianatoR::theme_bayesian()
    })

    # ── CrI width sensitivity plot ────────────────────────────────────────
    output$plot_sensitivity_cri <- shiny::renderPlot({
      sens <- sens_data()
      shiny::req(sens, nrow(sens) > 0)

      sens$cri_width <- sens$cri_upper - sens$cri_lower

      sens_summary <- sens |>
        dplyr::group_by(prior_mean) |>
        dplyr::summarise(
          cri_width = stats::median(cri_width),
          .groups   = "drop"
        )

      ggplot2::ggplot(sens_summary, ggplot2::aes(x = prior_mean, y = cri_width)) +
        ggplot2::geom_line(colour = "#4a8ec2", linewidth = 1) +
        ggplot2::labs(
          title = "Posterior 95% CrI width",
          x     = "Prior mean",
          y     = "CrI width"
        ) +
        bayesianatoR::theme_bayesian()
    })

    # ── P(benefit) sensitivity plot ───────────────────────────────────────
    output$plot_sensitivity_pbenefit <- shiny::renderPlot({
      sens <- sens_data()
      shiny::req(sens, nrow(sens) > 0)

      sens_summary <- sens |>
        dplyr::group_by(prior_mean) |>
        dplyr::summarise(
          p_benefit = stats::median(probability_of_benefit),
          .groups   = "drop"
        )

      ggplot2::ggplot(sens_summary, ggplot2::aes(x = prior_mean, y = p_benefit)) +
        ggplot2::geom_line(colour = "#27ae60", linewidth = 1) +
        ggplot2::geom_hline(
          yintercept = 0.95, linetype = "dashed",
          colour = "#e67e22", linewidth = 0.6
        ) +
        ggplot2::scale_y_continuous(
          limits = c(0, 1),
          labels = scales::percent_format(accuracy = 1)
        ) +
        ggplot2::labs(
          title    = "P(benefit) across prior spectrum",
          subtitle = "Dashed line: 95% threshold",
          x        = "Prior mean",
          y        = "P(benefit)"
        ) +
        bayesianatoR::theme_bayesian()
    })

    # ── Sensitivity table ─────────────────────────────────────────────────
    output$sensitivity_table <- shiny::renderTable({
      sens <- sens_data()
      shiny::req(sens, nrow(sens) > 0)
      num_cols <- sapply(sens, is.numeric)
      sens[num_cols] <- lapply(sens[num_cols], round, 4)
      head(sens, 50)
    })

    output$download_sens <- shiny::downloadHandler(
      filename = function() {
        paste0("bayesianatoR_sensitivity_", Sys.Date(), ".csv")
      },
      content = function(file) {
        utils::write.csv(sens_data(), file, row.names = FALSE)
      }
    )
  })
}

`%||%` <- rlang::`%||%`
