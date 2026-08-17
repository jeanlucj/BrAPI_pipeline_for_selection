# Console progress, status and timing reporting for the pipeline.
#
# Sourced by config.R, which every step sources, so these helpers are available
# everywhere without touching any step's source list.
#
# Everything here is GATED: it draws only in a live interactive console. Rscript,
# the test suite, and the non-interactive callr subprocess that wflow_build() uses
# all get complete silence, so nothing here can pollute the knitted report or a
# batch run. Override either way with options(brapi.progress = TRUE/FALSE).
#
# Bars are cli's (already a tidyverse/purrr dependency, and already what purrr's
# .progress uses). Two cli facts shaped this file:
#   - A progress bar is auto-terminated when the FRAME THAT CREATED IT returns, so
#     pb_start() must create it in the caller's frame (.envir = parent.frame());
#     creating it in the wrapper would kill the bar the moment pb_start() returned.
#   - cli displays one bar at a time (the innermost). "Nested" therefore means the
#     inner bar takes over the display and the outer resumes afterwards -- hence
#     inner labels carry their outer context ("protocol 66: imputing markers").

library(tidyverse)
here::i_am("code/progress.R")

# --- gate --------------------------------------------------------------------

# SHOW_PROGRESS is defined in config.R, which sources this file, so it does not
# exist yet at source time -- only when these functions are actually called.
.prog_on <- function() {
  isTRUE(getOption("brapi.progress",
                   if (exists("SHOW_PROGRESS")) SHOW_PROGRESS else interactive()))
}

# --- formatting helpers ------------------------------------------------------

.fmt_dur <- function(secs) {
  if (!is.finite(secs)) return("?")
  if (secs < 60)  return(sprintf("%.1fs", secs))
  mins <- floor(secs / 60); rem <- secs - 60 * mins
  if (mins < 60) return(sprintf("%dm %04.1fs", mins, rem))
  hrs <- floor(mins / 60)
  sprintf("%dh %02dm", hrs, mins - 60 * hrs)
}

.fmt_bytes <- function(n) {
  if (!is.finite(n)) return("? B")
  units <- c("B", "KB", "MB", "GB", "TB")
  k <- min(floor(log(max(n, 1)) / log(1024)), length(units) - 1)
  sprintf("%.1f %s", n / 1024^k, units[k + 1])
}

# Show project-relative paths -- an absolute here::here() path is noise.
.rel_path <- function(path) {
  root <- here::here()
  if (startsWith(path, root)) sub(paste0("^", root, "/?"), "", path) else path
}

# cli interpolates {} in its format strings, so any caller-supplied text must be
# passed as a VALUE (cli_alert_info("{txt}")), never pasted into the format.
# Literal text destined for a format string is escaped here instead.
.cli_escape <- function(x) gsub("\\}", "}}", gsub("\\{", "{{", x))

.txt <- function(...) paste0(...)

# --- status lines ------------------------------------------------------------

#' A status line: what the pipeline is about to do / just did.
say <- function(...) {
  if (.prog_on()) { txt <- .txt(...); cli::cli_alert_info("{txt}") }
  invisible(NULL)
}

#' A de-emphasised detail line (indented, grey) -- subordinate to the last say().
#' cli_verbatim, not cli_text: it neither interpolates {} (so arbitrary trait and
#' file names are safe) nor collapses the leading indent.
note <- function(...) {
  if (.prog_on()) cli::cli_verbatim(paste0("  ", cli::col_grey(.txt(...))))
  invisible(NULL)
}

#' Announce a cache hit. Steps return their cached result silently otherwise,
#' which makes a fast re-run indistinguishable from a hung one.
note_cache <- function(path) {
  if (.prog_on()) {
    info <- file.info(path)
    txt <- sprintf("using cached %s (%s, written %s) -- pass refresh = TRUE to recompute",
                   .rel_path(path), .fmt_bytes(info$size),
                   format(info$mtime, "%Y-%m-%d %H:%M"))
    cli::cli_alert_info("{txt}")
  }
  invisible(NULL)
}

# --- step banners + the run timing log ---------------------------------------

.prog_env <- new.env(parent = emptyenv())
.prog_env$timings <- list()

#' Start a pipeline step: print a banner and start its clock.
#' @return A handle to pass to step_done().
step_start <- function(title) {
  if (.prog_on()) cli::cli_rule(left = "{cli::style_bold(title)}")
  invisible(list(title = title, t0 = Sys.time()))
}

#' Finish a step: record its elapsed time and print it with an optional summary.
#' Timings are recorded even when reporting is off, so print_timings() still works.
step_done <- function(h, ..., cached = FALSE) {
  secs <- as.numeric(difftime(Sys.time(), h$t0, units = "secs"))
  .prog_env$timings <- c(.prog_env$timings,
                         list(tibble(step = h$title, seconds = secs, cached = cached)))
  if (.prog_on()) {
    # paste0() with no arguments is character(0), not "" -- so a step_done(h) with no
    # summary must not reach nzchar() unguarded (it would yield logical(0) and error).
    summary <- .txt(...)
    has_summary <- length(summary) == 1L && nzchar(summary)
    txt <- sprintf("%s (%s)%s", h$title, .fmt_dur(secs),
                   if (has_summary) paste0(" -- ", summary) else "")
    cli::cli_alert_success("{txt}")
  }
  invisible(secs)
}

#' Per-step elapsed times recorded so far, newest last.
pipeline_timings <- function() {
  if (!length(.prog_env$timings)) {
    return(tibble(step = character(), seconds = double(), cached = logical()))
  }
  list_rbind(.prog_env$timings)
}

#' Forget the recorded timings (a fresh run starts a fresh log).
timings_reset <- function() {
  .prog_env$timings <- list()
  invisible(NULL)
}

#' End-of-run summary table: every step, its elapsed time, and the total.
print_timings <- function(timings = pipeline_timings()) {
  if (!nrow(timings)) return(invisible(timings))
  total <- sum(timings$seconds)
  cli::cli_rule(left = "{cli::style_bold(paste0('Run complete (', .fmt_dur(total), ')'))}")
  width <- max(nchar(timings$step))
  pwalk(timings, function(step, seconds, cached) {
    txt <- sprintf("  %-*s %10s%s", width, step, .fmt_dur(seconds),
                   if (isTRUE(cached)) "   (cached)" else "")
    cli::cli_verbatim(txt)
  })
  invisible(timings)
}

# --- progress bars -----------------------------------------------------------

#' Start a progress bar in the CALLER's frame (see the .envir note at the top).
#'
#' @param total  Number of units of work.
#' @param name   Label; include the outer context when this bar nests inside another.
#' @return A bar id for pb_tick()/pb_done(), or NULL when reporting is off (both
#'   accept NULL, so callers need no `if`).
pb_start <- function(total, name, .envir = parent.frame()) {
  if (!.prog_on() || !is.finite(total) || total <= 0) return(NULL)
  cli::cli_progress_bar(
    format = paste0("{cli::pb_spin} ", .cli_escape(name),
                    " {cli::pb_current}/{cli::pb_total} {cli::pb_bar} ",
                    "{cli::pb_percent} | ETA {cli::pb_eta}"),
    total = total, .envir = .envir, clear = FALSE)
}

#' Advance a bar. `inc` batches ticks for very tight loops (redrawing every
#' iteration of a 17k-marker loop would itself cost time).
pb_tick <- function(id, inc = 1L) {
  if (!is.null(id)) cli::cli_progress_update(id = id, inc = inc)
  invisible(NULL)
}

pb_done <- function(id) {
  if (!is.null(id)) cli::cli_progress_done(id = id)
  invisible(NULL)
}

#' The `.progress =` argument for a purrr call: a cli bar spec, or FALSE when off.
pb_wrap <- function(name) {
  if (!.prog_on()) return(FALSE)
  list(format = paste0("{cli::pb_spin} ", .cli_escape(name),
                       " {cli::pb_current}/{cli::pb_total} {cli::pb_bar} ",
                       "{cli::pb_percent} | ETA {cli::pb_eta}"),
       clear = FALSE)
}

# --- phases: a timing breakdown WITHIN one step -------------------------------
# step_start/step_done time the seven pipeline steps; these time the phases inside
# one of them (step 4 spans downloads, marker parsing, imputation, per-protocol
# GRMs, the pedigree stitch and the EM combine, and "42 minutes" says nothing about
# which of those to attack).

.prog_env$phases <- list()

phases_reset <- function() {
  .prog_env$phases <- list()
  invisible(NULL)
}

#' Run `expr`, recording how long it took under `label`.
#'
#' Repeated labels accumulate, so a phase inside a per-protocol loop can be wrapped
#' at each iteration and still report one total. `cached = TRUE` marks the phase as
#' having been served from a cache rather than computed.
#' @return `expr`'s value, unchanged.
with_phase <- function(label, expr, cached = FALSE) {
  t0  <- Sys.time()
  res <- expr
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  prev <- .prog_env$phases[[label]]
  .prog_env$phases[[label]] <- list(
    seconds = (prev$seconds %||% 0) + secs,
    # a phase counts as cached only if EVERY contribution to it was
    cached  = (prev$cached %||% TRUE) && isTRUE(cached))
  res
}

#' Record phases that were skipped because a cache supplied their result, so they
#' still appear in the breakdown (marked cached) rather than vanishing from it.
phase_cached <- function(labels) {
  walk(labels, ~ with_phase(.x, NULL, cached = TRUE))
  invisible(NULL)
}

#' The phase breakdown recorded since the last phases_reset().
pipeline_phases <- function() {
  if (!length(.prog_env$phases)) {
    return(tibble(phase = character(), seconds = double(), cached = logical()))
  }
  tibble(phase   = names(.prog_env$phases),
         seconds = map_dbl(.prog_env$phases, "seconds"),
         cached  = map_lgl(.prog_env$phases, "cached"))
}

#' Print the phase breakdown, indented to sit under a step's own line.
print_phases <- function(phases = pipeline_phases()) {
  if (!.prog_on() || !nrow(phases)) return(invisible(phases))
  width <- max(nchar(phases$phase))
  pwalk(phases, function(phase, seconds, cached) {
    cli::cli_verbatim(sprintf("    %s %s %10s%s", phase,
                              strrep(".", width - nchar(phase) + 1),
                              .fmt_dur(seconds),
                              if (isTRUE(cached)) "  (cached)" else ""))
  })
  invisible(phases)
}

#' Time a single long call that has no loop to hang a bar on, e.g. the EM combine
#' or an eigendecomposition. Returns the expression's value unchanged.
time_it <- function(label, expr) {
  t0 <- Sys.time()
  say(label, " ...")
  res <- expr
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  note(label, ": ", .fmt_dur(secs))
  res
}
