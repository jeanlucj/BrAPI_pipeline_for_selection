# Tier 1: the request-keyed cache in code/cache.R (offline).
#
# The behaviour under test is the one whose absence caused a step to analyse the
# wrong trials: a cached result must be served ONLY when it was built for the same
# request, and an unrecorded provenance must count as "rebuild", not "assume fine".

cache_file <- function() file.path(withr::local_tempdir(.local_envir = parent.frame()),
                                   "thing.rds")

test_that("a cache written with a key is served back for the same key", {
  p <- cache_file()
  k <- cache_key(trials = c("B", "A"), radius = 500)
  cache_write(p, mtcars, k)
  expect_true(file.exists(cache_key_file(p)))
  expect_equal(suppressMessages(cache_read(p, k)), mtcars)
})

test_that("order and duplicates in a character vector are not a different request", {
  p <- cache_file()
  cache_write(p, "v", cache_key(trials = c("A", "B", "B")))
  expect_equal(suppressMessages(cache_read(p, cache_key(trials = c("B", "A")))), "v")
})

test_that("a changed request returns NULL and says what changed", {
  p <- cache_file()
  cache_write(p, "v", cache_key(trials = c("A", "B", "C"), radius = 500))
  expect_message(
    res <- cache_read(p, cache_key(trials = c("A", "B", "C", "D"), radius = 500)),
    "different request")
  expect_null(res)
  msgs <- capture_messages(cache_read(p, cache_key(trials = "A", radius = 300)))
  expect_true(any(grepl("trials: 3 -> 1 names", msgs)))
  expect_true(any(grepl("radius: 500 -> 300", msgs)))
})

test_that("a cache with no recorded request is treated as stale", {
  # Every cache built before this mechanism existed looks like this. Unknown
  # provenance is not evidence of a match.
  p <- cache_file()
  saveRDS("v", p)                                   # payload only, no sidecar
  expect_message(res <- cache_read(p, cache_key(x = 1)), "no recorded request")
  expect_null(res)
})

test_that("refresh and a missing payload both force a rebuild", {
  p <- cache_file()
  k <- cache_key(x = 1)
  cache_write(p, "v", k)
  expect_null(suppressMessages(cache_read(p, k, refresh = TRUE)))
  unlink(p)
  expect_null(suppressMessages(cache_read(p, k)))   # sidecar alone is not a cache
})

test_that("the key diff reports added, removed, unset and scalar changes", {
  d <- .key_diff(cache_key(a = c("x", "y"), b = 1),
                 cache_key(a = c("y", "z"), b = 2, c = "new"))
  expect_true(any(grepl("\\+z", d)))
  expect_true(any(grepl("-x", d)))
  expect_true(any(grepl("b: 1 -> 2", d)))
  expect_true(any(grepl("c: unset -> new", d)))
  expect_length(.key_diff(cache_key(a = 1), cache_key(a = 1)), 0)
})

test_that("the payload file keeps its plain format", {
  # The diagnostics in EVALUATION.md §9 read data/genotypes.rds and data/gebv.rds
  # directly, so the value must not be wrapped in a key envelope.
  p <- cache_file()
  cache_write(p, list(G = 1:3), cache_key(x = 1))
  expect_equal(readRDS(p), list(G = 1:3))
})
