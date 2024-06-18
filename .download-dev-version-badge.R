### * Description

# Download a CRAN badge with the current stable version

### * Setup

SITE = "https://img.shields.io/badge/"
MAIN = "In Dev"
STATUS = "1.1.5.9000"
COLOR = "orange"

### * Run

# Get the badge
target = paste(c(MAIN, STATUS, COLOR), collapse = "-")
target = gsub(" ", "%20", target)
link = paste0(c(SITE, target, ".svg"), collapse = "")
print("Badge for dev version downloaded from:")
print(link)
command = paste0("wget \"", link, "\" -O dev-version_badge.svg")
system(command)
