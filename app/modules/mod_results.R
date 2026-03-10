# app/modules/mod_results.R
# Results display module.
#
# Presents the posterior distribution (prior/likelihood/posterior overlay),
# numeric summaries (posterior mean, CrI, P(benefit), Bayes factor), and
# a downloadable result table. All outputs are reactive to prior changes.

# ── UI ────────────────────────────────────────────────────────────────────────

#' Results UI
#'
#' @param id Module namespace ID.
#' @export
mod_results_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(

    # ── Key numeric summary row ───────────────────────────────────────────
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),

      bslib::value_box(
        title    = "Posterior mean",
        value    = shiny::textOutput(ns("post_mean")),
        showcase = shiny::icon("bullseye"),
        theme    = "primary"
      ),
      bslib::value_box(
        title    = shiny::textOutput(ns("cri_label")),
        value    = shiny::textOutput(ns("post_cri")),
        showcase = shiny::icon("arrows-left-right"),
        theme    = "info"
      ),
      bslib::value_box(
        title    = "P(benefit)",
        value    = shiny::textOutput(ns("prob_benefit")),
        showcase = shiny::icon("chart-simple"),
        theme    = "success"
      ),
      bslib::value_box(
        title    = "Bayes factor (vs null)",
        value    = shiny::textOutput(ns("bf_null")),
        showcase = shiny::icon("scale-balanced"),
        theme    = "secondary"
      )
    ),

    # ── Plots ─────────────────────────────────────────────────────────────
    bslib::layout_columns(
      col_widths = c(8, 4),

      bslib::card(
        bslib::card_header(
          shiny::icon("chart-area"), " Prior \u2192 Posterior update"
        ),
        bslib::card_body(
          shiny::plotOutput(ns("plot_prior_posterior"), height = "400px")
        )
      ),

      bslib::card(
        bslib::card_header(
          shiny::icon("chart-column"), " Probability of benefit"
        ),
        bslib::card_body(
          shiny::plotOutput(ns("plot_probability"), height = "400px")
        )
      )
    ),

    # ── Full numeric table ────────────────────────────────────────────────
    bslib::card(
      bslib::card_header(
        bslib::layout_columns(
          col_widths = c(8, 4),
          shiny::div(shiny::icon("table"), " Full result summary"),
          shiny::div(
            class = "text-end",
            shiny::downloadButton(ns("download_csv"), "Download CSV",
                                  class = "btn-sm btn-outline-primary")
          )
        )
      ),
      bslib::card_body(
        shiny::tableOutput(ns("result_table"))
      )
    ),

    # ── Interpretation ────────────────────────────────────────────────────
    bslib::card(
      bslib::card_header(shiny::icon("comment-dots"), " Interpretation"),
      bslib::card_body(
        shiny::uiOutput(ns("interpretation"))
      )
    )
  )
}

# ── Server ────────────────────────────────────────────────────────────────────

#' Results server
#'
#' @param id Module namespace ID.
#' @param analysis_out Reactive `nv_result` object (or `NULL`) from the
#'   analysis dispatcher in `app.R`.
#' @export
mod_results_server <- function(id, analysis_out) {
  shiny::moduleServer(id, function(input, output, session) {

    result <- shiny::reactive({
      shiny::req(analysis_out())
      analysis_out()
    })

    # ── Value boxes ───────────────────────────────────────────────────────
    output$post_mean <- shiny::renderText({
      r <- result()
      sprintf("%.4f", r$posterior$mean)
    })

    output$cri_label <- shiny::renderText({
      r <- result()
      sprintf("%d%% CrI", round(100 * r$meta$ci_level))
    })

    output$post_cri <- shiny::renderText({
      r  <- result()
      ci <- r$credible_interval
      if (is.null(ci)) return("—")
      sprintf("[%.4f, %.4f]", ci$lower, ci$upper)
    })

    output$prob_benefit <- shiny::renderText({
      r <- result()
      if (is.null(r$probability_of_benefit)) return("—")
      sprintf("%.1f%%", 100 * r$probability_of_benefit)
    })

    output$bf_null <- shiny::renderText({
      r  <- result()
      bf <- r$bayes_factor$bf_null
      if (is.null(bf) || is.na(bf)) return("—")
      if (is.infinite(bf)) return("> 1000")
      sprintf("%.2f", bf)
    })

    # ── Plots ─────────────────────────────────────────────────────────────
    output$plot_prior_posterior <- shiny::renderPlot({
      result()$plots$prior_posterior
    })

    output$plot_probability <- shiny::renderPlot({
      result()$plots$probability
    })

    # ── Result table ──────────────────────────────────────────────────────
    output$result_table <- shiny::renderTable({
      tbl <- bayesianatoR::result_summary(result())
      # Round numerics for display
      num_cols <- sapply(tbl, is.numeric)
      tbl[num_cols] <- lapply(tbl[num_cols], round, 4)
      as.data.frame(t(tbl))
    }, rownames = TRUE, colnames = FALSE)

    output$download_csv <- shiny::downloadHandler(
      filename = function() {
        paste0("bayesianatoR_result_", Sys.Date(), ".csv")
      },
      content = function(file) {
        utils::write.csv(bayesianatoR::result_summary(result()), file,
                         row.names = FALSE)
      }
    )

    # ── Interpretation ────────────────────────────────────────────────────
    output$interpretation <- shiny::renderUI({
      r  <- result()
      bf <- r$bayes_factor$bf_null
      pb <- r$probability_of_benefit
      ci <- r$credible_interval

      bf_text <- if (!is.null(bf) && !is.na(bf)) {
        bayesianatoR::interpret_bf(bf)
      } else {
        "Bayes factor not computable for this prior type."
      }

      shiny::tagList(
        shiny::tags$ul(
          shiny::tags$li(
            sprintf(
              "The posterior mean is %.4f (95%% CrI: %.4f to %.4f), ",
              r$posterior$mean,
              if (!is.null(ci)) ci$lower else NA,
              if (!is.null(ci)) ci$upper else NA
            ),
            "based on combining the prior with the reported data."
          ),
          shiny::tags$li(
            sprintf(
              "The probability that the true effect exceeds the MID (%.3g) is %.1f%%.",
              r$meta$mid, 100 * (pb %||% NA_real_)
            )
          ),
          shiny::tags$li(
            sprintf("The Bayes factor vs H\u2080: %s.", bf_text)
          )
        ),
        shiny::p(
          class = "text-muted small",
          "This analysis assumes conjugate priors. Results are conditional on the ",
          "selected prior and should be interpreted alongside the original study ",
          "context. See the Sensitivity tab for robustness."
        )
      )
    })
  })
}

`%||%` <- rlang::`%||%`
