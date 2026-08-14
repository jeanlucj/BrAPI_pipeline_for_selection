# Tier 1: t3_login() credential handling (offline; no server contact).

# A stand-in for the BrAPI R6 connection: an environment so $login() can mutate
# $auth_token the way the real one does. `token` is what the "server" issues --
# "" reproduces a REJECTED password, which BrAPI signals by leaving auth_token
# empty rather than by erroring.
stub_conn <- function(token = "tok123") {
  e <- new.env()
  e$auth_token <- NULL
  e$calls <- list()
  e$login <- function(username = NULL, password = NULL, token_ = NULL) {
    e$calls <- c(e$calls, list(list(username = username, password = password)))
    e$auth_token <- token
    invisible(NULL)
  }
  e
}

# t3_login() reads the project .Renviron itself, which would drag a developer's
# real credentials into these tests. Neutralize that for the duration.
local_no_renviron <- function(env = parent.frame()) {
  old <- .load_project_renviron
  .load_project_renviron <<- function() invisible(FALSE)
  withr::defer(.load_project_renviron <<- old, envir = env)
}

test_that("t3_login errors when the credentials are unset, without calling login", {
  local_no_renviron()
  conn <- stub_conn()
  withr::with_envvar(c(T3_USERNAME = "", T3_PASSWORD = ""), {
    expect_error(t3_login(conn), class = "t3_missing_credentials")
  })
  withr::with_envvar(c(T3_USERNAME = "someone", T3_PASSWORD = ""), {
    expect_error(t3_login(conn), class = "t3_missing_credentials")
  })
  # Never reached conn$login() -- which would prompt interactively on empty args.
  expect_length(conn$calls, 0)
})

test_that("t3_login errors when the server issues no token", {
  local_no_renviron()
  conn <- stub_conn(token = "")          # rejected password: login() only warns
  withr::with_envvar(c(T3_USERNAME = "someone", T3_PASSWORD = "wrong"), {
    expect_error(t3_login(conn), class = "t3_bad_credentials")
  })
  expect_length(conn$calls, 1)
})

test_that("t3_login passes the credentials through and returns the connection", {
  local_no_renviron()
  conn <- stub_conn(token = "tok123")
  withr::with_envvar(c(T3_USERNAME = "someone", T3_PASSWORD = "secret"), {
    expect_silent(res <- t3_login(conn))
    expect_identical(res, conn)
  })
  expect_equal(conn$calls[[1]], list(username = "someone", password = "secret"))
  expect_equal(conn$auth_token, "tok123")
})

test_that("connect_t3(login = FALSE) skips authentication entirely", {
  local_no_renviron()
  # No env vars, no server: only reachable because login is skipped. The BrAPI
  # connection constructor itself does not touch the network.
  withr::with_envvar(c(T3_USERNAME = "", T3_PASSWORD = ""), {
    conn <- connect_t3(login = FALSE)
    expect_true(inherits(conn, "R6"))
    expect_null(conn$auth_token)
  })
})
