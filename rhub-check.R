### * Description

# This script can be used to check the package on various platforms using the
# rhub package and the winbuilder platform with the devtools package.

# See: https://cran.r-project.org/web/packages/rhub/vignettes/rhub.html

### * Setup

library(here)
library(crayon)
library(rhub)

# The first time rhub is used on your machine, you should validate your email
# with:
#   validate_email()
EMAIL <- "matthieu.bruneaux@gmail.com"

### * Functions

### ** start_report()

#' Open a file for report writing
#' @param name File name to which a timestamp and the ".md" extension are appended.
#' @return A file descriptor to the file `name`, located in the project root directory.
start_report <- function(name = "rhub-report") {
    name <- paste0(name, "_", format(Sys.time(), "%Y-%m-%d-%H%M%S"), ".md")
    file <- file(here::here(name), open = "w")
    cat("# Rhub check report, started on", format(Sys.time()), "\n", file = file)
    return(file)
}

### ** check_and_report()

#' Run a rhub check for one platform and save the report to an open file
check_and_report <- function(platform, file) {
    run <- check_and_wait(platform = platform)
    cat("\n## Check for", platform, "\n", file = file)
    write_report(run = run, file = file)
}

### ** close_report()

#' Close report file
close_report <- function(file) {
    cat("\n# End of report,", format(Sys.time()), "\n", file = file)
    close(file)
}

### ** check_and_wait()

#' Submit a check to rhub and wait for the run to complete
#'
#' @return Output from `check_for_cran()`
check_and_wait <- function(platform) {
    env_vars <- c("_R_CHECK_FORCE_SUGGESTS_" = "true",
                  "_R_CHECK_CRAN_INCOMING_USE_ASPELL_" = "true")
    if (grepl("macos", platform)) {
        warning("Spell check turned off for r-hub check on MacOS.",
                "See https://github.com/r-hub/rhub/issues/472")
        env_vars["_R_CHECK_CRAN_INCOMING_USE_ASPELL_"] <- "false"
    }
    run <- check_for_cran(platform = platform, email = EMAIL, env_vars = env_vars,
                          show_status = FALSE)
    completed <- FALSE
    Sys.sleep(15)
    while (!completed) {
        completed <- tryCatch(run$cran_summary(),
                              error = function(e) {
                                  if (grepl("At least one of the builds has not finished yet", e$message)) {
                                      message("Results not ready yet, pausing for 60 seconds.")
                                      Sys.sleep(60)
                                      return(FALSE)
                                  } else if (grepl("Timeout was reached", e$message)) {
                                      message("Timeout issue, pausing for 60 seconds.")
                                      Sys.sleep(60)
                                      return(FALSE)
                                  } else {
                                      return(TRUE)
                                  }})
        if (typeof(completed) == "environment") {
            completed <- TRUE
        }
    }
    return(run)
}

### ** write_report()

#' Write report for one check
write_report <- function(run, file) {
    z <- run$cran_summary()
    cat("\n```", crayon::strip_style(capture.output(z$print())), "```",
        file = file, sep = "\n")    
}

### * Check

# Checks with winbuilder
devtools::check_win_release(quiet = TRUE)
devtools::check_win_devel(quiet = TRUE)

# List rhub platforms
platforms()

# Platforms selected for check
targets <- c("macos-highsierra-release",
             "ubuntu-gcc-release")

report <- start_report()
for (x in targets) {
    check_and_report(platform = x, file = report)
}
close_report(report)

### * Alternative check on Windows

# For an alternative check on windows, one can run:
# devtools::check_win_devel(email = EMAIL)
