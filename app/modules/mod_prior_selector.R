# app/modules/mod_prior_selector.R
# Prior selection module.
#
# Provides radio buttons for the four canonical prior types plus a live
# ggplot preview and numeric summary of the selected prior.

# ── UI ────────────────────────────────────────────────────────────────────────

#' Prior selector UI
#'
#' @param id Module namespace ID.
#' @export
mod_prior_selector_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::layout_columns(
    col_widths = c(5, 7),

    # ── Left: prior type selection ─────────────────────────────────────────
    bslib::card(
      bslib::card_header(
        shiny::icon("sliders"), " Select prior"
      ),
      bslib::card_body(
        shiny::radioButtons(
          inputId  = ns("prior_type"),
          label    = "Prior type",
          choices  = c(
            "Enthusiastic (centred on expected effect)" = "enthusiastic",
            "Sceptical (centred on null)"               = "sceptical",
            "Non-informative (wide Normal)"             = "noninformative",
            "Custom"                                    = "custom"
          ),
          selected = "sceptical"
        ),
        shiny::hr(),

        # Enthusiastic / sceptical settings
        shiny::conditionalPanel(
          condition = sprintf(
            "input['%s'] == 'enthusiastic' || input['%s'] == 'sceptical'",
            ns("prior_type"), ns("prior_type")
          ),
          shiny::numericInput(
            ns("effect_size"),
            label = shiny::span(
              "Expected effect size",
              bslib::tooltip(shiny::icon("circle-question"),
                             "Typically from the power calculation or literature.")
            ),
            value = 0.5, step = 0.01
          ),
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'enthusiastic'", ns("prior_type")),
            shiny::numericInput(
              ns("effect_se"),
              label = shiny::span(
                "SE of expected effect",
                bslib::tooltip(shiny::icon("circle-question"),
                               "Uncertainty in the expected effect (used as prior SD).")
              ),
              value = 0.2, min = 0.001, step = 0.01
            )
          )
        ),

        # Custom settings
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] == 'custom'", ns("prior_type")),
          shiny::numericInput(
            ns("custom_mean"),
            label = "Prior mean",
            value = 0, step = 0.01
          ),
          shiny::numericInput(
            ns("custom_sd"),
            label = "Prior SD",
            value = 0.5, min = 0.001, step = 0.01
          )
        ),

        # Log-scale toggle (for ratio measures)
        shiny::checkboxInput(
          ns("log_scale"),
          label = shiny::span(
            "Prior on log scale (for OR, RR, HR)",
            bslib::tooltip(shiny::icon("circle-question"),
                           "Check this when analysing ratio measures.")
          ),
          value = FALSE
        ),

        shiny::hr(),
        # Prior numeric summary
        shiny::uiOutput(ns("prior_summary_text"))
      )
    ),

    # ── Right: prior density preview ──────────────────────────────────────
    bslib::card(
      bslib::card_header(
        shiny::icon("chart-line"), " Prior distribution preview"
      ),
      bslib::card_body(
        shiny::plotOutput(ns("prior_plot"), height = "380px"),
        shiny::p(
          class = "text-muted small mt-2",
          "Density shown on the analysis scale. For ratio measures with log-scale",
          " prior, the x-axis is log(ratio)."
        )
      )
    )
  )
}

# ── Server ────────────────────────────────────────────────────────────────────

#' Prior selector server
#'
#' @param id Module namespace ID.
#' @param wizard_out Reactive list from `mod_input_wizard_server` (used to
#'   infer a sensible default effect size from the data).
#'
#' @return A [shiny::reactive()] returning an `nv_prior` object, or `NULL`
#'   if inputs are incomplete.
#' @export
mod_prior_selector_server <- function(id, wizard_out = shiny::reactive(NULL)) {
  shiny::moduleServer(id, function(input, output, session) {

    # ── Construct prior reactively ────────────────────────────────────────
    prior_rv <- shiny::reactive({
      type <- input$prior_type
      shiny::req(type)

      tryCatch(
        switch(
          type,
          enthusiastic = {
            shiny::req(input$effect_size, input$effect_se)
            bayesianatoR::make_prior(
              "enthusiastic",
              effect_size = input$effect_size,
              se          = input$effect_se,
              log_scale   = isTRUE(input$log_scale)
            )
          },
          sceptical = {
            shiny::req(input$effect_size)
            bayesianatoR::make_prior(
              "sceptical",
              effect_size = input$effect_size,
              log_scale   = isTRUE(input$log_scale)
            )
          },
          noninformative = {
            bayesianatoR::make_prior(
              "noninformative",
              log_scale = isTRUE(input$log_scale)
            )
          },
          custom = {
            shiny::req(input$custom_mean, input$custom_sd)
            if (is.na(input$custom_sd) || input$custom_sd <= 0) return(NULL)
            bayesianatoR::make_prior(
              "custom",
              mean      = input$custom_mean,
              sd        = input$custom_sd,
              log_scale = isTRUE(input$log_scale)
            )
          }
        ),
        error = function(e) {
          shiny::showNotification(
            paste("Prior error:", conditionMessage(e)),
            type = "warning"
          )
          NULL
        }
      )
    })

    # ── Prior numeric summary ─────────────────────────────────────────────
    output$prior_summary_text <- shiny::renderUI({
      p <- prior_rv()
      if (is.null(p)) return(NULL)

      tbl <- bayesianatoR::prior_summary(p)

      bslib::card(
        class = "bg-light border-0",
        bslib::card_body(
          class = "py-2",
          shiny::tags$dl(
            class = "row mb-0 small",
            shiny::tags$dt(class = "col-5", "Mean"),
            shiny::tags$dd(class = "col-7", sprintf("%.4g", tbl$mean)),
            shiny::tags$dt(class = "col-5", "SD"),
            shiny::tags$dd(class = "col-7", sprintf("%.4g", tbl$sd)),
            shiny::tags$dt(class = "col-5", "95% range"),
            shiny::tags$dd(
              class = "col-7",
              sprintf("[%.4g, %.4g]", tbl$lower_95, tbl$upper_95)
            )
          )
        )
      )
    })

    # ── Prior density plot ────────────────────────────────────────────────
    output$prior_plot <- shiny::renderPlot({
      p <- prior_rv()
      if (is.null(p)) return(NULL)

      x_lo <- p$mean - 4 * p$sd
      x_hi <- p$mean + 4 * p$sd
      x    <- seq(x_lo, x_hi, length.out = 400)
      df   <- data.frame(x = x, density = dnorm(x, p$mean, p$sd))

      ggplot2::ggplot(df, ggplot2::aes(x = x, y = density)) +
        ggplot2::geom_area(fill = "#a8c4e0", alpha = 0.5) +
        ggplot2::geom_line(colour = "#1a6faf", linewidth = 1) +
        ggplot2::geom_vline(
          xintercept = p$mean, linetype = "dashed",
          colour = "#1a6faf", linewidth = 0.7
        ) +
        ggplot2::labs(
          title    = p$label,
          subtitle = sprintf(
            "95%% range: [%.4g, %.4g]",
            qnorm(0.025, p$mean, p$sd),
            qnorm(0.975, p$mean, p$sd)
          ),
          x        = if (isTRUE(p$log_scale)) "log(ratio)" else "Effect",
          y        = "Density"
        ) +
        bayesianatoR::theme_bayesian()
    })

    return(prior_rv)
  })
}
