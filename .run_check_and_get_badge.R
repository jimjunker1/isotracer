### * Description

# This script runs devtools::check() and parse the returned value to download a
# badge from https://shields.io/
#
# Two files are produced by this script:
# - R-CMD-check_badge.svg (badge)
# - R-CMD-check_output.txt (errors, warnings and notes from the check() call)

### * Run

# Run check
output = devtools::check(error_on = "never")
n = c("errors" = length(output[["errors"]]),
      "warnings" = length(output[["warnings"]]),
      "notes" = length(output[["notes"]]))

# Get the badge
SITE = "https://img.shields.io/badge/"
MAIN = "R CMD check"
STATUS = ""
COLOR = "lightgrey"
## Parse the check results
if (all(n == 0)) {
    STATUS = "passed"
    COLOR = "brightgreen"
}
if (n["errors"] > 0) {
    COLOR = "red"
    STATUS = paste(n["errors"], "error")
    if (n["errors"] > 1) { STATUS = paste0(STATUS, "s") }
}
if (n["warnings"] > 0) {
    if (COLOR == "lightgrey") { COLOR = "orange" }
    if (STATUS != "") { STATUS = paste0(STATUS, " ") }
    STATUS = paste0(STATUS, n["warnings"], " warning")
    if (n["warnings"] > 1) { STATUS = paste0(STATUS, "s") }
}
if (n["notes"] > 0) {
    if (COLOR == "lightgrey") { COLOR = "yellow" }
    if (STATUS != "") { STATUS = paste0(STATUS, " ") }
    STATUS = paste0(STATUS, n["notes"], " note")
    if (n["notes"] > 1) { STATUS = paste0(STATUS, "s") }
}
## Download the badge
target = paste(c(MAIN, STATUS, COLOR), collapse = "-")
target = gsub(" ", "%20", target)
link = paste0(c(SITE, target, ".svg"), collapse = "")
print("Badge for R CMD CHECK downloaded from:")
print(link)
command = paste0("wget \"", link, "\" -O R-CMD-check_badge.svg")
system(command)

# Save check() output
outFile = file("R-CMD-check_output.txt", open = "w")
cat("* ERRORS\n\n", file = outFile)
lapply(output$errors, function(e) {
    cat("** ", e, "\n\n", file = outFile)
})
cat("* WARNINGS\n\n", file = outFile)
lapply(output$warnings, function(w) {
    cat("** ", w, "\n\n", file = outFile)
})
cat("* NOTES\n\n", file = outFile)
lapply(output$notes, function(n) {
    cat("** ", n, "\n\n", file = outFile)
})
close(outFile)
