### * Description

# Script to prepare the dataset for the case study based on the article:
#
# Collins, Sarah M., Steven A. Thomas, Thomas Heatherly, Keeley L. MacNeill,
# Antoine O.H.C. Leduc, Andrés López-Sepulcre, Bradley A. Lamphere, et
# al. 2016. “Fish Introductions and Light Modulate Food Web Fluxes in Tropical
# Streams: A Whole-Ecosystem Experimental Approach.” Ecology, July,
# n/a-n/a. https://doi.org/10.1002/ecy.1530.

# Note: This script depends on the file `./data/trini_mod.rda`, which is
# created by the script `prep-large-model.R`.

### * Setup

suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(devtools::load_all())
options("dplyr.summarise.inform" = FALSE)

### * Main

# Load data
z <- load(here::here("data", "trini_mod.rda"))
stopifnot(z == "trini_mod")

# Rebuild "raw" data from trini_mod
z <- trini_mod
for (i in seq_len(nrow(z))) {
  z$observations[[i]][["stream"]] <- z$group[[i]][["stream"]]
  z$observations[[i]][["transect"]] <- z$group[[i]][["transect"]]
}
z <- dplyr::bind_rows(z$observations)
z$stream <- ifelse(z$stream == "LL", "LOL", "UPL")
z <- z[, c("stream", "transect", "compartment", "time", "size", "proportion")]
names(z) <- c("stream", "transect", "compartment", "time.days", "mgN.per.m2", "prop15N")

d <- z %>% 
  filter(compartment %in% c("NH4", "NO3", "epi", "FBOM", "petro", "tricor", "arg", "pseph")) %>%
  mutate(stream = ifelse(stream == "LOL", "LL", "UL")) %>%
  filter(stream == "UL")
# Get initial conditions for NH4 and NO3 (during drip)
inits <- d %>% filter(compartment %in% c("NH4", "NO3") & time.days <= 10) %>%
  group_by(stream, transect, compartment) %>%
  summarize(mgN.per.m2 = mean(mgN.per.m2, na.rm = TRUE),
            prop15N = mean(prop15N, na.rm = TRUE))
postdrip <- d %>% filter(compartment %in% c("NH4", "NO3") & time.days > 10) %>%
  group_by(stream, transect, compartment) %>%
  summarize(mgN.per.m2 = mean(mgN.per.m2, na.rm = TRUE),
            prop15N = mean(prop15N, na.rm = TRUE))
pulses <- (postdrip$prop15N - inits$prop15N) * inits[["mgN.per.m2"]]
inits2 <- d %>% filter(!compartment %in% c("NH4", "NO3")) %>%
  group_by(stream, compartment) %>%
  summarize(mgN.per.m2 = mean(mgN.per.m2, na.rm = TRUE)) %>%
  mutate(prop15N = delta2prop(0, "d15N")) %>%
  crossing(tibble(transect = c("transect.1", "transect.2", "transect.3")))
inits <- bind_rows(inits, inits2)

inits$time.days <- 0

lalaja <- dplyr::bind_rows(inits, d)

# Save data
save(lalaja, file = here::here("data", "lalaja.rda"), version = 3, compress = "xz")
