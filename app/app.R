# app/app.R
# bayesianatoR Shiny application entry point.
#
# Launch with: shiny::runApp("app")
# Or from project root: shiny::runApp()

library(shiny)
library(bslib)
library(bayesianatoR)

# Load modules
source(file.path("modules", "mod_input_wizard.R"))
source(file.path("modules", "mod_prior_selector.R"))
source(file.path("modules", "mod_results.R"))
source(file.path("modules", "mod_sensitivity.R"))

# ── Clinical bslib theme ──────────────────────────────────────────────────────

bnr_theme <- bslib::bs_theme(
  version    = 5,
  bg         = "#ffffff",
  fg         = "#1a3a5c",
  primary    = "#1a6faf",
  secondary  = "#4a8ec2",
  success    = "#27ae60",
  warning    = "#e67e22",
  danger     = "#c0392b",
  base_font  = bslib::font_google("Inter", wght = c(400, 600)),
  code_font  = bslib::font_google("Fira Code"),
  bootswatch = "flatly"
)

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- bslib::page_navbar(
  title = shiny::div(
    shiny::img(src = "logo.png", height = "30px",
               style = "margin-right: 8px; vertical-align: middle;",
               onerror = "this.style.display='none'"),
    "bayesianatoR"
  ),
  theme    = bnr_theme,
  fillable = FALSE,
  id       = "main_navbar",

  # ── Tab 1: Input wizard ───────────────────────────────────────────────────
  bslib::nav_panel(
    title = "1  Input",
    icon  = shiny::icon("file-medical"),
    value = "tab_input",
    mod_input_wizard_ui("wizard")
  ),

  # ── Tab 2: Prior selection ────────────────────────────────────────────────
  bslib::nav_panel(
    title = "2  Prior",
    icon  = shiny::icon("sliders"),
    value = "tab_prior",
    mod_prior_selector_ui("prior")
  ),

  # ── Tab 3: Results ────────────────────────────────────────────────────────
  bslib::nav_panel(
    title = "3  Results",
    icon  = shiny::icon("chart-area"),
    value = "tab_results",
    mod_results_ui("results")
  ),

  # ── Tab 4: Sensitivity analysis ───────────────────────────────────────────
  bslib::nav_panel(
    title = "4  Sensitivity",
    icon  = shiny::icon("rotate"),
    value = "tab_sensitivity",
    mod_sensitivity_ui("sensitivity")
  ),

  bslib::nav_spacer(),

  bslib::nav_item(
    shiny::tags$a(
      shiny::icon("circle-info"),
      "About",
      href   = "#",
      id     = "about_link",
      style  = "color: #4a8ec2;"
    )
  ),

  bslib::nav_item(
    shiny::tags$a(
      shiny::icon("github"),
      "GitHub",
      href   = "https://github.com/user/bayesianatoR",
      target = "_blank",
      style  = "color: #4a8ec2;"
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # ── Module servers ──────────────────────────────────────────────────────
  # wizard_out: reactive list with $analysis_type, $estimate, $se (and
  #   optionally $lower, $upper, $x1/$n1/$x2/$n2, etc.)
  wizard_out <- mod_input_wizard_server("wizard")

  # prior_out: reactive nv_prior object
  prior_out  <- mod_prior_selector_server("prior", wizard_out)

  # analysis_out: reactive nv_result (computed from wizard + prior)
  analysis_out <- shiny::reactive({
    shiny::req(wizard_out(), prior_out())
    inputs <- wizard_out()
    prior  <- prior_out()

    tryCatch(
      .run_analysis(inputs, prior),
      error = function(e) {
        shiny::showNotification(
          paste("Analysis error:", conditionMessage(e)),
          type     = "error",
          duration = 8
        )
        NULL
      }
    )
  })

  mod_results_server("results", analysis_out)
  mod_sensitivity_server("sensitivity", analysis_out)

  # ── Auto-advance tabs on wizard completion ──────────────────────────────
  shiny::observeEvent(wizard_out(), {
    if (!is.null(wizard_out())) {
      bslib::nav_select("main_navbar", "tab_prior")
    }
  })

  shiny::observeEvent(analysis_out(), {
    if (!is.null(analysis_out())) {
      bslib::nav_select("main_navbar", "tab_results")
    }
  })
}

# ── Analysis dispatcher ───────────────────────────────────────────────────────

#' Dispatch to the appropriate conjugate update function
#'
#' @param inputs Named list from `mod_input_wizard_server`.
#' @param prior An `nv_prior` object from `mod_prior_selector_server`.
#' @return An `nv_result` or `NULL`.
#' @keywords internal
.run_analysis <- function(inputs, prior) {
  type <- inputs$analysis_type

  switch(
    type,
    "two_group_continuous" = update_two_group_continuous(
      estimate = inputs$estimate,
      se       = inputs$se,
      prior    = prior,
      ci_level = inputs$ci_level %||% 0.95,
      mid      = inputs$mid      %||% 0
    ),
    "two_group_proportions" = update_two_group_proportions(
      x1          = inputs$x1,
      n1          = inputs$n1,
      x2          = inputs$x2,
      n2          = inputs$n2,
      prior_alpha = inputs$prior_alpha %||% 1,
      prior_beta  = inputs$prior_beta  %||% 1,
      ci_level    = inputs$ci_level    %||% 0.95,
      mid         = inputs$mid         %||% 0
    ),
    "odds_ratio" = update_odds_ratio(
      estimate   = inputs$estimate,
      lower      = inputs$lower,
      upper      = inputs$upper,
      prior      = prior,
      conf_level = inputs$conf_level %||% 0.95,
      ci_level   = inputs$ci_level   %||% 0.95,
      mid        = inputs$mid        %||% 0,
      measure    = "OR"
    ),
    "risk_ratio" = update_odds_ratio(
      estimate   = inputs$estimate,
      lower      = inputs$lower,
      upper      = inputs$upper,
      prior      = prior,
      conf_level = inputs$conf_level %||% 0.95,
      ci_level   = inputs$ci_level   %||% 0.95,
      mid        = inputs$mid        %||% 0,
      measure    = "RR"
    ),
    "paired_continuous" = update_paired_continuous(
      estimate = inputs$estimate,
      se       = inputs$se,
      prior    = prior,
      ci_level = inputs$ci_level %||% 0.95,
      mid      = inputs$mid      %||% 0
    ),
    "single_proportion" = update_single_proportion(
      x           = inputs$x,
      n           = inputs$n,
      prior_alpha = inputs$prior_alpha %||% 1,
      prior_beta  = inputs$prior_beta  %||% 1,
      null_p      = inputs$null_p      %||% 0.5,
      ci_level    = inputs$ci_level    %||% 0.95,
      mid         = inputs$mid         %||% 0.5
    ),
    "correlation" = update_correlation(
      r        = inputs$r,
      n        = inputs$n,
      prior    = prior,
      ci_level = inputs$ci_level %||% 0.95,
      mid      = inputs$mid      %||% 0
    ),
    "hazard_ratio" = update_hazard_ratio(
      estimate   = inputs$estimate,
      lower      = inputs$lower,
      upper      = inputs$upper,
      prior      = prior,
      conf_level = inputs$conf_level %||% 0.95,
      ci_level   = inputs$ci_level   %||% 0.95,
      mid        = inputs$mid        %||% 0
    ),
    "regression_coef" = update_regression_coef(
      estimate = inputs$estimate,
      se       = inputs$se,
      prior    = prior,
      ci_level = inputs$ci_level %||% 0.95,
      mid      = inputs$mid      %||% 0
    ),
    rlang::abort(paste("Unknown analysis type:", type))
  )
}

# ── Launch ────────────────────────────────────────────────────────────────────

shiny::shinyApp(ui = ui, server = server)
