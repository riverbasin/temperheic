# R Package Conventions (Always Active)

## Code style rules

- Use `%>%` pipe for data transformation chains
- Use purrr (`map`, `map_dfr`, `safely`, `possibly`) instead of for loops or apply family
- Use tidy evaluation (`{{ }}`) for column name arguments in user-facing functions
- Use `match.arg()` for method/type selection parameters
- Function factories for multiple implementations of similar functionality
- `group_by() %>% nest() %>% mutate(result = map(...))` for split-apply-combine

## roxygen2 documentation

- Every exported function gets `@description`, `@param`, `@return`, `@examples`
- Use `@details` for methodological decisions and equation references
- Cross-reference source papers by equation number: "Luce 2013 Eq. 5" or "observedAmpPhase.tex Eq. 1-3"
- Include working examples that run in < 5 seconds

## Testing conventions

- testthat 3rd edition
- Test file naming: `test-{function_name}.R`
- Round-trip validation pattern: generate synthetic data → fit → verify parameter recovery
- Tolerance: 1e-6 for analytical tests, 0.05 (5%) for real-data tests
- Edge cases to always test: NA handling, single-column input, short records, mismatched timestamps

## Error handling

- `cli::cli_abort()` for user-facing errors with informative messages
- `cli::cli_warn()` for recoverable issues (e.g., gaps in data, weak signal)
- Validate units early: check that K is in plausible m/s range, not m/d

## Commit messages

- Imperative mood: "Add sigma_eta uncertainty function" not "Added..."
- Reference phase/issue: "Phase 1: Add sigma_eta() — Luce 2013 Eq. 13"
