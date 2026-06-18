# Step 1: connect to the T3/Oat BrAPI server.

library(tidyverse)
# NOTE: BrAPI.R calls timeout() (from httr) unqualified, so httr must be
# attached to the search path for any conn$get()/conn$search() call to work.
library(httr)
here::i_am("code/01_connect.R")

source(here::here("code", "config.R"))

#' Connect to a T3/Breedbase BrAPI server.
#'
#' Uses anonymous (public) access by default. Set `login = TRUE` to authenticate
#' with the T3_USER / T3_PASS environment variables -- needed only for restricted
#' data or if a public call unexpectedly returns HTTP 401.
#'
#' @param db_name Pre-configured BrAPI connection name (see BrAPI::listBrAPIConnections()).
#' @param login   Whether to authenticate using T3_USER / T3_PASS env vars.
#' @return A BrAPI connection object.
connect_t3 <- function(db_name = DB_NAME, login = FALSE) {
  conn <- BrAPI::getBrAPIConnection(db_name)
  if (login) {
    user <- Sys.getenv("T3_USER")
    pass <- Sys.getenv("T3_PASS")
    if (nzchar(user) && nzchar(pass)) {
      conn$login(user, pass)
    } else {
      warning("login = TRUE but T3_USER / T3_PASS are not set; continuing anonymously.")
    }
  }
  conn
}
