## This makes sure that R loads the workflowr package
## automatically, everytime the project is loaded
if (requireNamespace("workflowr", quietly = TRUE)) {
  message("Loading .Rprofile for the current workflowr project")
  library("workflowr")
} else {
  message("workflowr package not installed, please run install.packages(\"workflowr\") to use the workflowr functions")
}

## Silence macOS RStudio's benign "MallocStackLogging: can't turn off ..."
## diagnostic that surfaces during memory-heavy steps (Stage 1 lme4, Stage 2
## BGLR). It is environmental (RStudio injects MallocStackLogging), not from the
## pipeline, and does not affect results. Restart the R session for it to take
## full effect.
if (Sys.info()[["sysname"]] == "Darwin") {
  Sys.unsetenv(c("MallocStackLogging", "MallocStackLoggingNoCompact",
                 "MallocScribble", "MallocGuardEdges"))
}
