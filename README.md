# bayesianatoR

> Bayesian reanalysis tool for clinical papers — R package + Shiny app

bayesianatoR lets clinicians and researchers re-examine the results of published clinical trials through a Bayesian lens, without requiring specialist statistical software. It uses conjugate (analytical) priors so updates are instant — no MCMC required.

---

## Project structure

```
bayesianatoR/
├── bayesianatoR/          # R package
│   ├── DESCRIPTION
│   ├── NAMESPACE
│   ├── R/
│   │   ├── input_parsers.R      # CI/SE/SD conversions, input validation
│   │   ├── priors.R             # Prior construction helpers
│   │   ├── conjugate_updates.R  # Core Bayesian updating functions
│   │   ├── posteriors.R         # Posterior summary methods
│   │   ├── bayes_factors.R      # BF vs null and ROPE
│   │   ├── sensitivity.R        # Prior sensitivity sweep
│   │   └── plots.R              # ggplot2 plot methods
│   └── tests/testthat/
├── app/                   # Shiny application
│   ├── app.R
│   ├── modules/
│   │   ├── mod_input_wizard.R   # Outcome type → reporting format router
│   │   ├── mod_prior_selector.R # Prior type + live preview
│   │   ├── mod_results.R        # Posterior plots + numeric outputs
│   │   └── mod_sensitivity.R    # Sensitivity analysis tab
│   └── www/
├── README.md
├── NEWS.md
└── CONTRIBUTING.md
```

---

## Supported analysis types

| Outcome type              | Model                        |
|---------------------------|------------------------------|
| Two-group continuous      | Normal–Normal conjugate      |
| Two-group proportions     | Beta–Binomial conjugate      |
| Odds ratio / Risk ratio   | Log-normal (log scale)       |
| Paired continuous         | Normal–Normal conjugate      |
| Single proportion vs null | Beta–Binomial conjugate      |
| Correlation coefficient   | Fisher z + Normal–Normal     |
| Hazard ratio              | Log-normal (log scale)       |
| Regression coefficient    | Normal–Normal conjugate      |

---

## Input formats accepted

- Mean difference + 95%/90%/99% CI
- Mean difference + SE
- Mean difference + SD + n (pooled or Welch SE)
- Mean difference + p-value (with warning — approximate)
- Ratio + CI (OR, RR, HR — log-transformed internally)

---

## Prior types

| Label           | Mean        | SD               |
|-----------------|-------------|------------------|
| Enthusiastic    | effect_size | se               |
| Sceptical       | 0           | \|effect_size\|/2|
| Non-informative | 0           | 10               |
| Custom          | user-specified | user-specified|

---

## Installation

```r
# Install from local source
pak::pkg_install("path/to/bayesianatoR/bayesianatoR")

# Or with devtools
devtools::install("bayesianatoR/bayesianatoR")
```

## Run the Shiny app

```r
shiny::runApp("app")
```

---

## Quick start (package only)

```r
library(bayesianatoR)

# 1. Specify a prior
prior <- make_prior("sceptical", effect_size = 0.5)

# 2. Run the Bayesian update
result <- update_two_group_continuous(
  estimate = 0.4,   # observed mean difference
  se       = 0.15,  # standard error
  prior    = prior,
  mid      = 0.2    # minimum important difference
)

# 3. Inspect results
print(result)
result_summary(result)

# 4. View plots
result$plots$prior_posterior
result$plots$sensitivity
result$plots$probability
```

---

## Running tests

```r
devtools::test("bayesianatoR/bayesianatoR")
```

---

## Tech stack

- **R** ≥ 4.1, **tidyverse**, **ggplot2**
- Conjugate/analytical Bayesian updating (no Stan/MCMC for v1)
- **Shiny** + **bslib** v5 UI
- **testthat** v3 for package testing
- **roxygen2** for documentation

---

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development guidelines, branch conventions, and the pull request checklist.

---

## Licence

MIT — see [bayesianatoR/LICENSE.md](bayesianatoR/LICENSE.md).
