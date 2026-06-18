#!/usr/bin/env Rscript
# Run the test suite.
#   Rscript tests/run_tests.R                  # tiers 1 + 2 (offline)
#   RUN_LIVE_TESTS=true Rscript tests/run_tests.R   # also tier 3 (live network)
suppressMessages(library(testthat))
here::i_am("tests/run_tests.R")
testthat::test_dir(here::here("tests", "testthat"), stop_on_failure = TRUE)
