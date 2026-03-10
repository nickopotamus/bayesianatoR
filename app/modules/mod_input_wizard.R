# app/modules/mod_input_wizard.R
# Decision-tree input wizard module.
#
# The user first selects their *outcome type* (continuous, proportion, ratio,
# correlation, regression), then their *reporting format* (CI, SE, SD, etc.),
# and is routed to the appropriate input form. Validation warnings are
# displayed in a bslib callout block — never in a blocking modal.

# ── UI ────────────────────────────────────────────────────────────────────────

#' Input wizard UI
#'
#' @param id Module namespace ID.
#' @export
mod_input_wizard_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::layout_columns(
    col_widths = c(4, 8),

    # ── Left panel: outcome type selector ──────────────────────────────────
    bslib::card(
      bslib::card_header(
        shiny::icon("list-check"), " Step 1: Outcome type"
      ),
      bslib::card_body(
        shiny::radioButtons(
          inputId  = ns("outcome_type"),
          label    = "What type of outcome did the paper report?",
          choices  = c(
            "Two-group continuous (mean diff)" = "two_group_continuous",
            "Two-group proportions"            = "two_group_proportions",
            "Odds ratio (OR)"                  = "odds_ratio",
            "Risk ratio (RR)"                  = "risk_ratio",
            "Paired continuous (pre/post)"     = "paired_continuous",
            "Single proportion vs null"        = "single_proportion",
            "Correlation coefficient"          = "correlation",
            "Hazard ratio (HR)"                = "hazard_ratio",
            "Regression coefficient"           = "regression_coef"
          ),
          selected = "two_group_continuous"
        )
      )
    ),

    # ── Right panel: dynamic input form ────────────────────────────────────
    bslib::card(
      bslib::card_header(
        shiny::icon("keyboard"), " Step 2: Enter reported statistics"
      ),
      bslib::card_body(
        # Reporting format selector
        shiny::uiOutput(ns("format_selector")),
        shiny::hr(),
        # Dynamic form
        shiny::uiOutput(ns("input_form")),
        shiny::hr(),
        # Optional settings
        bslib::accordion(
          bslib::accordion_panel(
            title = "Additional settings",
            open  = FALSE,
            shiny::numericInput(
              ns("ci_level"),
              label = shiny::span(
                "Credible interval level",
                bslib::tooltip(shiny::icon("circle-question"),
                               "Probability mass inside the posterior CrI.")
              ),
              value = 0.95, min = 0.80, max = 0.99, step = 0.01
            ),
            shiny::numericInput(
              ns("mid"),
              label = shiny::span(
                "Minimum important difference (MID)",
                bslib::tooltip(shiny::icon("circle-question"),
                               "P(benefit) = P(\u03b8 > MID). Set to 0 for superiority.")
              ),
              value = 0, step = 0.01
            ),
            shiny::numericInput(
              ns("sample_size"),
              label = shiny::span(
                "Total sample size (optional, for validation)",
                bslib::tooltip(shiny::icon("circle-question"),
                               "Used for soft plausibility checks only.")
              ),
              value = NA
            )
          )
        ),
        shiny::hr(),
        # Validation warnings
        shiny::uiOutput(ns("validation_warnings")),
        # Submit button
        shiny::actionButton(
          ns("submit"),
          label = "Run analysis",
          icon  = shiny::icon("play"),
          class = "btn-primary mt-2"
        )
      )
    )
  )
}

# ── Server ────────────────────────────────────────────────────────────────────

#' Input wizard server
#'
#' @param id Module namespace ID.
#'
#' @return A [shiny::reactive()] returning a named list of parsed inputs,
#'   or `NULL` before first submission. The list always contains
#'   `$analysis_type`, `$ci_level`, and `$mid`.
#' @export
mod_input_wizard_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Format selector (depends on outcome type) ─────────────────────────
    output$format_selector <- shiny::renderUI({
      choices <- .format_choices(input$outcome_type)
      shiny::radioButtons(
        inputId  = ns("reporting_format"),
        label    = "What did the paper report?",
        choices  = choices,
        selected = choices[[1]]
      )
    })

    # ── Dynamic input form ────────────────────────────────────────────────
    output$input_form <- shiny::renderUI({
      shiny::req(input$outcome_type, input$reporting_format)
      .build_input_form(input$outcome_type, input$reporting_format, ns)
    })

    # ── Validation warnings ───────────────────────────────────────────────
    parsed_inputs <- shiny::reactive({
      .parse_form_inputs(input, input$outcome_type, input$reporting_format)
    })

    output$validation_warnings <- shiny::renderUI({
      pi <- parsed_inputs()
      if (is.null(pi) || is.null(pi$estimate) || is.null(pi$se)) return(NULL)

      msgs <- bayesianatoR::validate_inputs(
        estimate   = pi$estimate,
        se         = pi$se,
        lower      = pi$lower,
        upper      = pi$upper,
        pvalue     = pi$pvalue,
        n          = input$sample_size,
        conf_level = input$ci_level %||% 0.95
      )

      if (length(msgs) == 0) return(NULL)

      bslib::card(
        class = "border-warning",
        bslib::card_header(
          class = "bg-warning text-dark",
          shiny::icon("triangle-exclamation"), " Input validation warnings"
        ),
        bslib::card_body(
          shiny::tags$ul(
            lapply(msgs, function(m) shiny::tags$li(m))
          )
        )
      )
    })

    # ── Submitted result ──────────────────────────────────────────────────
    submitted <- shiny::reactiveVal(NULL)

    shiny::observeEvent(input$submit, {
      pi <- parsed_inputs()
      if (is.null(pi)) {
        shiny::showNotification(
          "Please complete the input form before running the analysis.",
          type = "warning"
        )
        return()
      }
      # Attach meta-settings
      pi$analysis_type <- input$outcome_type
      pi$ci_level      <- input$ci_level %||% 0.95
      pi$mid           <- input$mid      %||% 0
      submitted(pi)
    })

    return(submitted)
  })
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' Reporting format choices by outcome type
#' @keywords internal
.format_choices <- function(outcome_type) {
  ci_formats <- c(
    "Mean difference + 95% CI"  = "estimate_ci95",
    "Mean difference + 90% CI"  = "estimate_ci90",
    "Mean difference + SE"      = "estimate_se",
    "Mean difference + SD + n"  = "estimate_sd_n",
    "Mean difference + p-value" = "estimate_pval"
  )
  ratio_formats <- c(
    "Ratio + 95% CI"  = "ratio_ci95",
    "Ratio + 90% CI"  = "ratio_ci90",
    "Ratio + p-value" = "ratio_pval"
  )

  switch(
    outcome_type,
    two_group_continuous  = ci_formats,
    paired_continuous     = ci_formats,
    regression_coef       = ci_formats,
    odds_ratio            = ratio_formats,
    risk_ratio            = ratio_formats,
    hazard_ratio          = ratio_formats,
    two_group_proportions = c("Event counts (x/n)"         = "counts"),
    single_proportion     = c("Event count + total"        = "count_n"),
    correlation           = c("r + sample size"            = "r_n"),
    c("Standard inputs" = "default")
  )
}

#' Build the dynamic input form
#' @keywords internal
.build_input_form <- function(outcome_type, reporting_format, ns) {
  # Helper: numeric input with tooltip
  ni <- function(id, label, tooltip_text, value = NA, min = NA, step = 0.01) {
    shiny::div(
      class = "mb-2",
      shiny::numericInput(
        ns(id),
        label = shiny::span(
          label,
          bslib::tooltip(shiny::icon("circle-question"), tooltip_text)
        ),
        value = value, min = min, step = step
      )
    )
  }

  if (outcome_type %in% c("two_group_continuous", "paired_continuous",
                           "regression_coef")) {
    switch(
      reporting_format,
      estimate_ci95 = shiny::tagList(
        ni("estimate", "Mean difference",
           "Observed mean difference (group 1 \u2212 group 2)"),
        ni("lower", "Lower 95% CI", "Lower bound of the 95% confidence interval"),
        ni("upper", "Upper 95% CI", "Upper bound of the 95% confidence interval")
      ),
      estimate_ci90 = shiny::tagList(
        ni("estimate", "Mean difference", "Observed mean difference"),
        ni("lower", "Lower 90% CI", "Lower bound of the 90% CI"),
        ni("upper", "Upper 90% CI", "Upper bound of the 90% CI")
      ),
      estimate_se = shiny::tagList(
        ni("estimate", "Mean difference", "Observed mean difference"),
        ni("se_direct", "Standard error (SE)", "SE of the mean difference")
      ),
      estimate_sd_n = shiny::tagList(
        ni("estimate", "Mean difference", "Observed mean difference"),
        ni("sd1", "SD (group 1)", "Standard deviation in group 1"),
        ni("sd2", "SD (group 2)", "Standard deviation in group 2"),
        ni("n1", "n (group 1)", "Sample size in group 1", step = 1, min = 2),
        ni("n2", "n (group 2)", "Sample size in group 2", step = 1, min = 2)
      ),
      estimate_pval = shiny::tagList(
        ni("estimate", "Mean difference", "Observed mean difference"),
        ni("pvalue", "p-value", "Two-sided p-value", min = 0, step = 0.001)
      ),
      shiny::p("Form not yet implemented for this format.", class = "text-muted")
    )

  } else if (outcome_type %in% c("odds_ratio", "risk_ratio", "hazard_ratio")) {
    switch(
      reporting_format,
      ratio_ci95 = shiny::tagList(
        ni("estimate", "Ratio (natural scale)",
           "OR, RR or HR — must be > 0", value = 1),
        ni("lower", "Lower 95% CI", "Lower bound; must be > 0"),
        ni("upper", "Upper 95% CI", "Upper bound; must be > 0")
      ),
      ratio_ci90 = shiny::tagList(
        ni("estimate", "Ratio (natural scale)", "Must be > 0", value = 1),
        ni("lower", "Lower 90% CI", "Lower bound; must be > 0"),
        ni("upper", "Upper 90% CI", "Upper bound; must be > 0")
      ),
      ratio_pval = shiny::tagList(
        ni("estimate", "Ratio (natural scale)", "Must be > 0", value = 1),
        ni("pvalue", "p-value", "Two-sided p-value", min = 0, step = 0.001)
      )
    )

  } else if (outcome_type == "two_group_proportions") {
    shiny::tagList(
      ni("x1", "Events (group 1)", "Number of events in group 1",
         step = 1, min = 0),
      ni("n1", "Total (group 1)", "Total participants in group 1",
         step = 1, min = 1),
      ni("x2", "Events (group 2)", "Number of events in group 2",
         step = 1, min = 0),
      ni("n2", "Total (group 2)", "Total participants in group 2",
         step = 1, min = 1)
    )

  } else if (outcome_type == "single_proportion") {
    shiny::tagList(
      ni("x", "Events", "Number of events", step = 1, min = 0),
      ni("n", "Total", "Total participants", step = 1, min = 1),
      ni("null_p", "Null proportion", "H\u2080 value for Bayes factor",
         value = 0.5, min = 0.01)
    )

  } else if (outcome_type == "correlation") {
    shiny::tagList(
      ni("r", "Pearson r", "Correlation coefficient (-1 < r < 1)"),
      ni("n", "Sample size", "Total sample size", step = 1, min = 4)
    )

  } else {
    shiny::p("Input form not yet implemented for this outcome type.",
             class = "text-muted")
  }
}

#' Parse form inputs into a standardised list
#' @keywords internal
.parse_form_inputs <- function(input, outcome_type, reporting_format) {
  out <- list()

  tryCatch({
    if (outcome_type %in% c("two_group_continuous", "paired_continuous",
                             "regression_coef")) {
      out$estimate <- input$estimate
      out$lower    <- input$lower
      out$upper    <- input$upper

      out$se <- switch(
        reporting_format,
        estimate_ci95  = {
          shiny::req(input$estimate, input$lower, input$upper)
          bayesianatoR::ci_to_se(input$lower, input$upper, conf_level = 0.95)
        },
        estimate_ci90  = {
          shiny::req(input$estimate, input$lower, input$upper)
          bayesianatoR::ci_to_se(input$lower, input$upper, conf_level = 0.90)
        },
        estimate_se    = {
          shiny::req(input$se_direct)
          input$se_direct
        },
        estimate_sd_n  = {
          shiny::req(input$sd1, input$sd2, input$n1, input$n2)
          bayesianatoR::sd_to_se(input$sd1, input$sd2, input$n1, input$n2)
        },
        estimate_pval  = {
          shiny::req(input$estimate, input$pvalue)
          suppressWarnings(
            bayesianatoR::pvalue_to_se(input$estimate, input$pvalue)
          )
        },
        NA
      )

    } else if (outcome_type %in% c("odds_ratio", "risk_ratio", "hazard_ratio")) {
      shiny::req(input$estimate)
      out$estimate    <- input$estimate
      out$lower       <- input$lower
      out$upper       <- input$upper
      out$conf_level  <- if (reporting_format == "ratio_ci90") 0.90 else 0.95

      # Derive SE on log scale for validation display
      if (!is.null(input$lower) && !is.null(input$upper) &&
          !is.na(input$lower) && !is.na(input$upper) &&
          input$lower > 0 && input$upper > 0 && input$estimate > 0) {
        lt      <- bayesianatoR::log_transform_ratio(
          input$estimate, input$lower, input$upper,
          conf_level = out$conf_level
        )
        out$se    <- lt$log_se
        out$lower <- input$lower
        out$upper <- input$upper
      }

    } else if (outcome_type == "two_group_proportions") {
      shiny::req(input$x1, input$n1, input$x2, input$n2)
      out$x1 <- as.integer(input$x1)
      out$n1 <- as.integer(input$n1)
      out$x2 <- as.integer(input$x2)
      out$n2 <- as.integer(input$n2)
      # Approximate risk-diff for validation
      p1 <- out$x1 / out$n1
      p2 <- out$x2 / out$n2
      out$estimate <- p1 - p2
      out$se <- sqrt(p1 * (1 - p1) / out$n1 + p2 * (1 - p2) / out$n2)

    } else if (outcome_type == "single_proportion") {
      shiny::req(input$x, input$n)
      out$x      <- as.integer(input$x)
      out$n      <- as.integer(input$n)
      out$null_p <- input$null_p %||% 0.5
      out$estimate <- out$x / out$n
      out$se <- sqrt(out$estimate * (1 - out$estimate) / out$n)

    } else if (outcome_type == "correlation") {
      shiny::req(input$r, input$n)
      out$r  <- input$r
      out$n  <- as.integer(input$n)
      out$estimate <- atanh(input$r)
      out$se <- 1 / sqrt(out$n - 3)
    }
  }, error = function(e) NULL)

  if (length(out) == 0) return(NULL)
  out
}

`%||%` <- rlang::`%||%`
