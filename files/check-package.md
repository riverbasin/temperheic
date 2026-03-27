# Check Package

Run the full temperheic quality assurance pipeline:

1. Run `devtools::document()` to rebuild roxygen docs
2. Run `devtools::test()` and report results (should be 0 failures)
3. Run `devtools::check()` and report errors/warnings/notes
4. If any tests fail, investigate the failure and suggest a fix
5. Summarize: total tests, pass/fail count, any check issues

Report results concisely. If everything passes, just say so.
