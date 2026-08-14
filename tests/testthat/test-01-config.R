# Tier 1: the data/config/*.txt readers in config.R (offline).

# Write a fixture file with NO trailing newline -- the real config files are
# exported that way, and readLines() would warn on them.
write_cfg <- function(lines) {
  p <- tempfile(fileext = ".txt")
  cat(paste(lines, collapse = "\n"), file = p)
  p
}

test_that("config_lines trims, drops blanks and comments, keeps values verbatim", {
  p <- write_cfg(c("# which accessions to predict", "", "  APPLER|CIAV2680  ",
                   "3|PI344799", "   ", "ARSNY25W001001"))
  expect_equal(config_lines(p),
               c("APPLER|CIAV2680", "3|PI344799", "ARSNY25W001001"))
})

test_that("config_lines drops duplicates and returns NULL for an empty list", {
  expect_message(x <- config_lines(write_cfg(c("x", "y", "x"))), "duplicate")
  expect_equal(x, c("x", "y"))
  expect_null(config_lines(write_cfg(c("", "# nothing here"))))
})

test_that("config_lines errors on a missing file rather than returning NULL", {
  # A silent NULL would mean "no training trials" -> fall back to the radius
  # search, i.e. a typo would quietly run a completely different analysis.
  expect_error(config_lines("no_such_file_xyz.txt"), "config file not found")
})

test_that("config_traits with one column gives names only", {
  p <- write_cfg(c("Grain yield - g/m2|CO_350:0000260",
                   "Freeze damage severity - 0-9 Rating|CO_350:0005001"))
  tr <- config_traits(p)
  expect_equal(tr$names, c("Grain yield - g/m2|CO_350:0000260",
                           "Freeze damage severity - 0-9 Rating|CO_350:0005001"))
  expect_null(tr$weights)
  expect_null(tr$short)
})

test_that("config_traits reads weights and short names keyed by trait name", {
  p <- write_cfg(c("Grain yield - g/m2|CO_350:0000260\t0.5\tYield g/m2",
                   "Lodging severity - 0-9 Rating|CO_350:0005007\t-1\tLodging 0-9"))
  tr <- config_traits(p)
  expect_equal(unname(tr$weights), c(0.5, -1))
  expect_equal(names(tr$weights), tr$names)
  expect_equal(names(tr$short), tr$names)
  expect_equal(unname(tr$short), c("Yield g/m2", "Lodging 0-9"))
})

test_that("config_traits rejects a half-filled or non-numeric weight column", {
  expect_error(config_traits(write_cfg(c("A\t1", "B"))), "some lines but not on: B")
  expect_error(config_traits(write_cfg(c("A\thigh", "B\t1"))), "non-numeric weight")
})

test_that("config_traits handles a single trait", {
  tr <- config_traits(write_cfg("Grain yield - g/m2|CO_350:0000260\t2"))
  expect_equal(tr$names, "Grain yield - g/m2|CO_350:0000260")
  expect_equal(unname(tr$weights), 2)
})
