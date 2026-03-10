# bayesianatoR NEWS

## bayesianatoR 0.0.1 (development)

### New features

* Initial scaffold of the `bayesianatoR` R package.
* Eight conjugate Bayesian update functions:
  - `update_two_group_continuous()` — Normal–Normal conjugate
  - `update_two_group_proportions()` — Beta–Binomial conjugate
  - `update_odds_ratio()` — log-Normal (OR and RR)
  - `update_paired_continuous()` — Normal–Normal conjugate
  - `update_single_proportion()` — Beta–Binomial conjugate
  - `update_correlation()` — Fisher z + Normal–Normal
  - `update_hazard_ratio()` — log-Normal
  - `update_regression_coef()` — Normal–Normal conjugate
* Prior construction via `make_prior()` (enthusiastic / sceptical /
  noninformative / custom).
* Input parsers: `ci_to_se()`, `se_to_ci()`, `sd_to_se()`,
  `pvalue_to_se()`, `log_transform_ratio()`, `exp_transform_result()`.
* Soft input validation via `validate_inputs()` (non-blocking warnings).
* Posterior summaries: `compute_posterior_summary()`, `result_summary()`,
  `print.nv_result()`.
* Bayes factors: `compute_bayes_factor()`, `interpret_bf()` using the
  Savage–Dickey density ratio.
* Prior sensitivity sweep: `compute_sensitivity()`, `get_sensitivity()`.
* ggplot2 plots: `plot_prior_posterior()`, `plot_sensitivity()`,
  `plot_probability()`, `theme_bayesian()`.
* Shiny application scaffolded in `app/`:
  - `mod_input_wizard` — decision-tree input UI
  - `mod_prior_selector` — prior type selector with live preview
  - `mod_results` — posterior plots and numeric outputs
  - `mod_sensitivity` — sensitivity analysis tab
* Unit tests for input parsers, prior construction, and all conjugate
  update functions.

### Breaking changes

None (initial release).

### Bug fixes

None (initial release).

### Internal changes

* `rlang` used throughout for error handling.
* All functions return consistently structured `nv_result` S3 objects.
