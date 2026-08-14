# Tier 1: the console reporting helpers in code/progress.R (offline).
#
# Two things matter here: with reporting OFF nothing is printed and no value is
# altered (the helpers are dropped inline throughout the pipeline, so they must be
# transparent), and with it ON the step timings are still recorded correctly.

on_prog <- function(code) withr::with_options(list(brapi.progress = TRUE), code)
off_prog <- function(code) withr::with_options(list(brapi.progress = FALSE), code)

test_that("reporting off: every helper is silent", {
  off_prog({
    expect_silent(say("hello"))
    expect_silent(note("detail"))
    expect_silent({
      h <- step_start("Step")
      step_done(h, "summary")
    })
    expect_silent({
      pb <- pb_start(10, "things")
      pb_tick(pb); pb_done(pb)
    })
    expect_silent(time_it("work", 1 + 1))
  })
})

test_that("reporting off: helpers pass values through unchanged", {
  off_prog({
    expect_equal(time_it("work", 6 * 7), 42)
    expect_null(pb_start(10, "things"))       # no bar to tick
    expect_false(pb_wrap("things"))           # purrr .progress = FALSE
  })
})

test_that("pb_start declines a bar for a degenerate total", {
  on_prog({
    expect_null(pb_start(0, "nothing"))
    expect_null(pb_start(NA_real_, "unknown"))
  })
})

test_that("pb_wrap returns a cli bar spec when reporting is on", {
  on_prog({
    spec <- pb_wrap("studies")
    expect_type(spec, "list")
    expect_match(spec$format, "studies", fixed = TRUE)
  })
})

test_that("step timings accumulate one row per step, in order", {
  timings_reset()
  expect_equal(nrow(pipeline_timings()), 0L)

  h <- step_start("1 First"); step_done(h)
  h <- step_start("2 Second"); step_done(h, cached = TRUE)
  tm <- pipeline_timings()

  expect_equal(tm$step, c("1 First", "2 Second"))
  expect_equal(tm$cached, c(FALSE, TRUE))
  expect_true(all(tm$seconds >= 0))
  timings_reset()
  expect_equal(nrow(pipeline_timings()), 0L)
})

test_that("timings are recorded even when reporting is off", {
  timings_reset()
  off_prog({ h <- step_start("Quiet step"); step_done(h) })
  expect_equal(pipeline_timings()$step, "Quiet step")
  timings_reset()
})

test_that("durations and sizes are formatted for humans", {
  expect_equal(.fmt_dur(0.94), "0.9s")
  expect_equal(.fmt_dur(84.2), "1m 24.2s")
  expect_equal(.fmt_dur(3720), "1h 02m")
  expect_equal(.fmt_bytes(2048), "2.0 KB")
  expect_match(.fmt_bytes(5.5 * 1024^3), "^5\\.5 GB$")
})

test_that("braces in a caller-supplied label cannot break cli's format string", {
  # cli interpolates {} in bar formats, so a label carrying a brace (a trait name,
  # a file name) must be escaped rather than evaluated.
  expect_equal(.cli_escape("a{b}c"), "a{{b}}c")
  on_prog({
    pb <- expect_silent(pb_start(2, "trait {x}"))
    pb_done(pb)
  })
})
