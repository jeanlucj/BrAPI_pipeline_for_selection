# Step 1: connect to the T3/Oat BrAPI server.

library(tidyverse)
# NOTE: BrAPI.R calls timeout() (from httr) unqualified, so httr must be
# attached to the search path for any conn$get()/conn$search() call to work.
library(httr)
here::i_am("code/01_connect.R")

source(here::here("code", "config.R"))

#' Read the project `.Renviron` (repo root), if present.
#'
#' R only auto-reads `.Renviron` from the directory it was STARTED in, and never
#' walks up to parents -- so credentials would be invisible when R is launched
#' from `analysis/`, from a `workflowr::wflow_build()` subprocess, or from
#' anywhere but the repo root. Reading it explicitly here also means an edit takes
#' effect on the next `connect_t3()` rather than only after restarting R.
#'
#' @return TRUE if the file existed and was read, FALSE otherwise (invisibly).
.load_project_renviron <- function() {
  f <- here::here(".Renviron")
  if (file.exists(f)) readRenviron(f)
  invisible(file.exists(f))
}

#' Log a BrAPI connection in with the T3_USERNAME / T3_PASSWORD env vars.
#'
#' Never prompts: `conn$login()` falls back to `readline()`/`askpass()` when its
#' arguments are empty, which would hang a non-interactive run (e.g. wflow_build).
#' Fails fast with a typed condition instead -- `t3_missing_credentials` when the
#' env vars are unset, `t3_bad_credentials` when the server rejects them.
#'
#' @param conn A BrAPI connection object (mutated in place: it stores the token).
#' @return `conn`, invisibly.
t3_login <- function(conn) {
  .load_project_renviron()
  user <- Sys.getenv("T3_USERNAME")
  pass <- Sys.getenv("T3_PASSWORD")

  if (!nzchar(user) || !nzchar(pass)) {
    stop(structure(
      class = c("t3_missing_credentials", "error", "condition"),
      list(message = paste(
        "T3 login needs T3_USERNAME and T3_PASSWORD in the environment.",
        "Copy .Renviron.example to .Renviron in the repo root and fill it in.",
        "Check with Sys.getenv(\"T3_USERNAME\")."),
        call = NULL)))
  }

  conn$login(username = user, password = pass)

  # conn$login() returns NORMALLY on a rejected password (it only warns), leaving
  # auth_token NULL -- so the token, not the return value, is what says it worked.
  if (!isTRUE(nzchar(conn$auth_token))) {
    stop(structure(
      class = c("t3_bad_credentials", "error", "condition"),
      list(message = paste0(
        "T3 login was REJECTED for user '", user, "' -- the server issued no ",
        "token (look for an 'Incorrect Password' warning above). Fix ",
        "T3_USERNAME / T3_PASSWORD in .Renviron."),
        call = NULL)))
  }

  invisible(conn)
}

#' Connect to a T3/Breedbase BrAPI server.
#'
#' Authenticates by default: T3/Oat no longer serves this pipeline's data
#' anonymously, so an un-logged-in connection returns empty results rather than an
#' obvious error. Credentials come from T3_USERNAME / T3_PASSWORD (see
#' `.Renviron.example`). Pass `login = FALSE` only for offline/mock use.
#'
#' @param db_name Pre-configured BrAPI connection name (see BrAPI::listBrAPIConnections()).
#' @param login   Whether to authenticate with T3_USERNAME / T3_PASSWORD.
#' @return A BrAPI connection object.
connect_t3 <- function(db_name = DB_NAME, login = TRUE) {
  # Status lines live here rather than in t3_login(), which stays silent (it is
  # the unit under test in tests/testthat/test-01-connect.R).
  say("Connecting to ", db_name, if (login) " (logging in)" else " (no login)", " ...")
  conn <- BrAPI::getBrAPIConnection(db_name)
  if (login) {
    t3_login(conn)
    note("logged in as ", Sys.getenv("T3_USERNAME"), " at ", conn$host)
  }
  conn
}
