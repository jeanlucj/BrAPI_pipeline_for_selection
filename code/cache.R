# Request-keyed caching: a cached result is reused only when it was produced by the
# SAME request.
#
# Why this exists: every step caches to data/*.rds, and until this file the FILENAME
# was the entire cache key -- so changing TRAINING_TRIALS, TRAIT_NAMES, a genotyping
# protocol or an EM weight left the old result in place and the run silently analysed
# the wrong data. The only defence was a line in README.md telling you to delete the
# caches yourself. Now each step records what it was asked for; a step whose request
# has changed says what changed and rebuilds.
#
#   key <- cache_key(training_trials = tt, test_trials = te)
#   hit <- cache_read(cache_path("ny_trials.rds"), key, refresh)
#   if (!is.null(hit)) return(hit)
#   ... compute ...
#   cache_write(cache_path("ny_trials.rds"), out, key)
#
# The request is stored in a SIDECAR (data/<name>.key.rds), not wrapped around the
# payload: the payload file keeps its exact format, so the diagnostics in
# EVALUATION.md §9 can go on reading data/genotypes.rds and data/gebv.rds directly.

library(tidyverse)
here::i_am("code/cache.R")

#' Assemble a request key from the things that determined a result.
#'
#' Pass only what actually shaped the output -- adding parameters the step ignored
#' causes rebuilds that change nothing. Character vectors are sorted and deduplicated
#' so that a reordered config file is not a different request.
cache_key <- function(...) {
  key <- list(...)
  map(key, function(v) {
    if (is.null(v)) return(NULL)
    if (is.character(v) || is.factor(v)) return(sort(unique(as.character(v))))
    if (is.numeric(v) && length(v) > 1) return(sort(unique(v)))
    v
  })
}

cache_key_file <- function(path) paste0(tools::file_path_sans_ext(path), ".key.rds")

# Human-readable differences between two keys, most useful first. Empty = same request.
.key_diff <- function(old, new) {
  fields <- union(names(old), names(new))
  diffs <- map_chr(fields, function(f) {
    a <- old[[f]]; b <- new[[f]]
    if (isTRUE(all.equal(a, b))) return(NA_character_)
    if (is.null(a)) return(sprintf("%s: unset -> %s", f, .key_show(b)))
    if (is.null(b)) return(sprintf("%s: %s -> unset", f, .key_show(a)))
    if (is.character(a) && is.character(b)) {
      added <- setdiff(b, a); dropped <- setdiff(a, b)
      return(sprintf("%s: %d -> %d names%s%s", f, length(a), length(b),
                     .key_examples("+", added), .key_examples("-", dropped)))
    }
    sprintf("%s: %s -> %s", f, .key_show(a), .key_show(b))
  })
  diffs[!is.na(diffs)]
}

.key_show <- function(v, n = 6) {
  if (is.null(v)) return("unset")
  if (length(v) == 1) return(as.character(v))
  # Show short vectors in full: "df: 60 60 30 -> 120 120 30" tells you what changed;
  # "2 values -> 2 values" does not.
  if (length(v) <= n && (is.numeric(v) || is.logical(v))) {
    return(paste(format(v, trim = TRUE), collapse = " "))
  }
  sprintf("%d values", length(v))
}

.key_examples <- function(sign, v, n = 3) {
  if (!length(v)) return("")
  sprintf(" (%s%s%s)", sign, paste(head(v, n), collapse = paste0(", ", sign)),
          if (length(v) > n) ", ..." else "")
}

#' Read a cached value, but only if it was built for `key`.
#'
#' @return The cached value, or NULL meaning "rebuild" -- which happens when
#'   `refresh` is set, the payload is missing, the request was never recorded
#'   (a cache built by hand or before this mechanism existed: unknown provenance
#'   is not evidence of a match), or the recorded request differs.
#' Note on visibility: a rebuild message is NOT routed through say()/note(), which
#' are gated by SHOW_PROGRESS. It explains why an expensive step is suddenly re-running
#' and which setting caused it, so it must survive a batch run. The cache-HIT line
#' stays gated -- that one is ordinary progress narration.
cache_read <- function(path, key, refresh = FALSE) {
  if (refresh || !file.exists(path)) return(NULL)
  kf <- cache_key_file(path)
  if (!file.exists(kf)) {
    message(.rel_path(path), ": no recorded request (built before request-keyed ",
            "caching, or by hand) -- rebuilding")
    return(NULL)
  }
  diffs <- .key_diff(read_rds(kf), key)
  if (length(diffs)) {
    message(.rel_path(path), " was built for a different request -- rebuilding\n",
            paste0("    ", diffs, collapse = "\n"))
    return(NULL)
  }
  note_cache(path)
  read_rds(path)
}

#' Write a value and record the request that produced it.
cache_write <- function(path, value, key) {
  write_rds(value, path)
  write_rds(key, cache_key_file(path))
  invisible(value)
}
