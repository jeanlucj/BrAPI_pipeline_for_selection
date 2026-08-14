# Tier 1: the evaluation tooling in code/evaluation.R (offline).
#
# The point of the first test is to keep EVALUATION.md honest: if a function is
# renamed or removed, the group that names it must be updated too, or arming that
# group silently reaches fewer functions than the doc claims.

source(here::here("code", "evaluation.R"))   # not sourced by config.R, on purpose

test_that("every name in EVAL_GROUPS resolves to a real function", {
  # canonicalize_to_primary comes from T3_brapi_helpers via synonyms.R, and the
  # sommer helper needs the package; allow only those to be absent.
  optional <- c("canonicalize_to_primary")
  missing <- EVAL_GROUPS |>
    unlist(use.names = FALSE) |>
    unique() |>
    setdiff(optional) |>
    keep(~ !exists(.x, mode = "function"))
  expect_equal(missing, character(0))
})

test_that("EVAL_ORDER and EVAL_ONLINE refer to real groups", {
  expect_true(all(EVAL_ORDER %in% names(EVAL_GROUPS)))
  expect_true(all(EVAL_ONLINE %in% names(EVAL_GROUPS)))
  # Every group except the step4 aggregate should be walkable in order.
  expect_setequal(setdiff(names(EVAL_GROUPS), "step4"), EVAL_ORDER)
})

test_that("arm/disarm toggle isdebugged() and are idempotent", {
  withr::defer(disarm_evaluation())
  expect_false(isdebugged(get("stage1_blues")))

  expect_message(n1 <- arm_evaluation("stage1"), "armed")
  expect_true(isdebugged(get("stage1_blues")))
  expect_true(isdebugged(get(".fit_one")))

  expect_message(n2 <- arm_evaluation("stage1"), "armed 0")   # already armed
  expect_equal(n2, 0L)
  expect_gt(n1, 0L)

  expect_message(disarm_evaluation("stage1"), "disarmed")
  expect_false(isdebugged(get("stage1_blues")))
  expect_message(off2 <- disarm_evaluation("stage1"), "disarmed 0")
  expect_equal(off2, 0L)
})

test_that("armed_functions reports what is armed, and disarm clears everything", {
  withr::defer(disarm_evaluation())
  expect_equal(armed_functions(), character(0))
  arm_evaluation("select") |> suppressMessages()
  expect_true("select_parents" %in% armed_functions())
  disarm_evaluation() |> suppressMessages()
  expect_equal(armed_functions(), character(0))
})

test_that("an unknown group is an error, not a silent no-op", {
  expect_error(arm_evaluation("nope"), "unknown evaluation group")
  expect_error(disarm_evaluation("nope"), "unknown evaluation group")
})

test_that("eval_groups lists every group with its loaded count", {
  out <- suppressMessages(eval_groups())
  expect_setequal(out$group, names(EVAL_GROUPS))
  expect_true(all(out$loaded <= out$n_fns))
  expect_true(all(out$layer %in% c("offline", "online")))
})

test_that("peek returns its input unchanged for each shape it handles", {
  G <- std_grm(.Gmatrix(make_dosage(8, 30)))
  dos <- make_dosage(8, 30)
  ph  <- tibble(studyDbId = "s1", germplasmName = "A", trait = "t", value = 1)
  bl  <- tibble(trait = "t", studyDbId = "s1", genotype = "A", BLUE = 1,
                SE = 0.5, weight = 4)
  gb  <- tibble(trait = "t", genotype = c("A", "B"), GEBV = c(0.1, -0.2),
                phenotyped = c(TRUE, FALSE))
  for (obj in list(G, dos, ph, bl, gb, c(1, 2, 3))) {
    res <- NULL
    invisible(capture.output(res <- peek(obj)))   # peek prints; the value passes through
    expect_identical(res, obj)
  }
})

test_that("peek flags the signatures it exists to catch", {
  dos <- make_dosage(6, 20)
  out <- capture.output(peek(dos, accessions = c("ZZZ1", "ZZZ2")))
  expect_true(any(grepl("ZERO overlap", out)))              # synonym mismatch

  bad <- make_dosage(6, 20) - 1                             # {-1,0,1} coding
  expect_true(any(grepl("outside \\{0,1,2\\}", capture.output(peek(bad)))))

  ph <- tibble(studyDbId = "s1", germplasmName = "A", trait = "t", value = 1,
               rep = NA_character_, block = NA_character_)
  expect_true(any(grepl("ALL NA", capture.output(peek(ph)))))

  # An injected training accession: diagonal at the mean, off-diagonals all 0.
  G <- std_grm(.Gmatrix(make_dosage(6, 30)))
  G <- rbind(cbind(G, 0), 0); dimnames(G) <- list(c(rownames(G)[1:6], "GHOST"),
                                                  c(colnames(G)[1:6], "GHOST"))
  G["GHOST", "GHOST"] <- mean(diag(G)[1:6])
  expect_true(any(grepl("off-diagonals exactly 0", capture.output(peek(G)))))
})
