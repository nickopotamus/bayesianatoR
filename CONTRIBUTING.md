# Contributing to bayesianatoR

Thank you for contributing! This document is aimed at developers working on
the R package (`bayesianatoR/`) or the Shiny application (`app/`).

---

## Project layout

```
bayesianatoR/
├── bayesianatoR/          # R package source
│   ├── R/                 # Function source files
│   └── tests/testthat/    # Unit tests
├── app/                   # Shiny app (separate from package)
│   └── modules/
├── README.md
├── NEWS.md
├── CONTRIBUTING.md
└── CLAUDE.md              # Claude Code persistent instructions
```

---

## Development setup

```r
# Install package dependencies
pak::pak(desc::desc_get_deps("bayesianatoR/bayesianatoR")$package)

# Load the package in development mode
devtools::load_all("bayesianatoR/bayesianatoR")

# Build documentation
devtools::document("bayesianatoR/bayesianatoR")

# Run tests
devtools::test("bayesianatoR/bayesianatoR")

# Run the Shiny app
shiny::runApp("app")
```

---

## Code style

* Follow the [tidyverse style guide](https://style.tidyverse.org/).
* Use `rlang::abort()` for all errors, `rlang::warn()` for warnings.
* Never use `stop()`, `warning()`, or `message()` directly.
* Use `rlang::check_required()` for mandatory arguments.
* All public functions must have complete roxygen2 documentation with at
  least one `@examples` entry.
* Internal helpers should be prefixed with `.` and tagged `@keywords internal`.

---

## Statistical conventions

* All Bayesian updates use conjugate priors — no MCMC or numerical integration.
* For ratio measures (OR, RR, HR), computations are on the **log scale**.
  Input/output conversion is handled in `input_parsers.R`.
* The Bayes factor against H₀ uses the **Savage–Dickey density ratio**.
  See `bayes_factors.R` for the reference.
* Credible intervals are always symmetric Normal intervals around the
  posterior mean (exact for Normal posteriors).

---

## Adding a new analysis type

1. Add a `update_<type>()` function in `bayesianatoR/R/conjugate_updates.R`
   following the established pattern:
   - Validate inputs with `.check_*` helpers.
   - Compute posterior via `normal_normal_update()` (or custom update).
   - Wrap in `new_nv_result()`.
   - Pass to `.finalise_result()`.
   - Return the result.

2. Export the function in `NAMESPACE` (or run `devtools::document()`).

3. Add unit tests in `tests/testthat/test-conjugate_updates.R`:
   - At least one test verifying the exact posterior against analytical values.
   - At least one integration test (all slots populated, `result_summary()` works).
   - At least one error test (bad inputs abort with informative messages).

4. Add the analysis type to the Shiny wizard:
   - `mod_input_wizard.R`: add to `radioButtons`, `.format_choices()`, and
     `.build_input_form()`.
   - `app.R`: add a case to `.run_analysis()`.

5. Update `NEWS.md` with the new function.

---

## Shiny module conventions

* All modules follow the `mod_<name>_ui()` / `mod_<name>_server()` naming
  convention.
* Module servers return a **reactive** (not a reactiveVal) where possible.
* Use `shiny::req()` defensively to avoid errors when upstream reactives
  are `NULL`.
* All validation warnings are shown as bslib cards with class
  `"border-warning"` — never as modal dialogs.
* Use `bslib::tooltip()` for every non-obvious input field.

---

## Testing

```r
# Run all tests
devtools::test("bayesianatoR/bayesianatoR")

# Run a single test file
testthat::test_file("bayesianatoR/tests/testthat/test-conjugate_updates.R")

# Check test coverage
covr::package_coverage("bayesianatoR/bayesianatoR")
```

### Testing checklist

- [ ] Unit tests for exact analytical results (tolerance 1e-10)
- [ ] Round-trip tests for input parsers (CI → SE → CI)
- [ ] Error tests for all input validation
- [ ] Integration test: all `nv_result` slots are populated
- [ ] Integration test: `result_summary()` returns a one-row tibble
- [ ] Integration test: `make_plots()` returns named list of ggplot objects

---

## Pull request checklist

- [ ] Tests pass: `devtools::test("bayesianatoR/bayesianatoR")`
- [ ] Documentation builds: `devtools::document("bayesianatoR/bayesianatoR")`
- [ ] `R CMD CHECK` passes with 0 errors, 0 warnings: `devtools::check("bayesianatoR/bayesianatoR")`
- [ ] `NEWS.md` updated with the change
- [ ] No new `stop()` / `warning()` / `message()` calls (use `rlang`)
- [ ] New public functions have `@examples` in roxygen docs

---

## Branch conventions

| Branch prefix | Purpose                              |
|---------------|--------------------------------------|
| `feat/`       | New features                         |
| `fix/`        | Bug fixes                            |
| `docs/`       | Documentation-only changes           |
| `test/`       | Test additions without code changes  |
| `refactor/`   | Code reorganisation without new features |

---

## Reporting issues

Please open a GitHub issue with:
1. The R session info (`sessionInfo()`)
2. A minimal reproducible example
3. The expected vs actual output
