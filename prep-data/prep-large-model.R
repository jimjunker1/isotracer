### * Description

# Prepare a large network model (based on Collins 2016)

# This script prepares a network model similar to the one used in
# López-Sepulcre et al. 2020 (Trinidadian mountain streams based on Collins
# 2016).
#
# It is slightly different from the actual model used in López-Sepulcre
# et al., as it estimates a common zeta for all compartments in a given stream
# instead of using compartment-specific coefficients of variation for sizes.

### * Setup

suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(devtools::load_all())

### * Main

### ** Prepare data

# Load data
col_types <- c(
    stream = col_character(),
    transect = col_character(),
    compartment = col_character(),
    time.days = col_double(),
    mgN.per.m2 = col_double(),
    d15N = col_double(),
    prop15N = col_double()
)
data <- read_tsv("collins-dump.tsv", col_types = col_types)  %>%
  mutate(stream = ifelse(stream == "LOL", "LL", "UL")) %>%
  select(-d15N)

# Get initial conditions for NH4 and NO3 (during drip)
inits1 <- data %>% filter(compartment %in% c("NH4", "NO3") & time.days <= 10) %>%
  group_by(stream, transect, compartment) %>%
  summarize(mgN.per.m2 = mean(mgN.per.m2, na.rm = TRUE),
            prop15N = mean(prop15N, na.rm = TRUE),
            .groups = "drop")
inits2 <- data %>% filter(!compartment %in% c("NH4", "NO3")) %>%
  group_by(stream, compartment) %>%
    summarize(mgN.per.m2 = mean(mgN.per.m2, na.rm = TRUE),
              .groups = "drop") %>%
  mutate(prop15N = delta2prop(0, "d15N")) %>%
  crossing(tibble(transect = c("transect.1", "transect.2", "transect.3")))
inits <- bind_rows(inits1, inits2)
inits[["mgN.per.m2"]][inits[["stream"]] == "LL" & inits[["compartment"]] == "FBOM"] <- 5267 # Value from UL

# Get pulse info
postdrip <- data %>% filter(compartment %in% c("NH4", "NO3") & time.days > 10) %>%
  group_by(stream, transect, compartment) %>%
  summarize(mgN.per.m2 = mean(mgN.per.m2, na.rm = TRUE),
            prop15N = mean(prop15N, na.rm = TRUE), .groups = "drop")
pulses <- (postdrip$prop15N - inits1$prop15N) * inits1[["mgN.per.m2"]]

### ** Prepare model

links <- c("NH4, NO3 -> epi, seston, FBOM, CBOM",
           "epi -> petro, pseph", "seston -> lepto",
           "FBOM -> tricor", "CBOM -> eudan, phylo",
           "tricor -> arg, euthy", "petro -> arg")

m <- new_networkModel(quiet = TRUE) %>%
    set_topo(links) %>%
    set_steady(c("NH4", "NO3")) %>%
    set_split(c("epi", "CBOM", "FBOM", "seston"))

suppressMessages({
m <- m %>%
    set_init(inits, comp = "compartment", size = "mgN.per.m2", prop = "prop15N",
             group_by = c("stream", "transect"))

m <- m %>%
    set_obs(data %>% filter(time.days > 0), comp = "compartment",
            size = "mgN.per.m2", prop = "prop15N", time = "time.days",
            group_by = c("stream", "transect"))
})

m %<>% add_pulse_event(time = 11, comp = "NH4", unmarked = 0, marked = -0.00538, which = 1) %>%
  add_pulse_event(time = 11, comp = "NO3", unmarked = 0, marked = -0.0168, which = 1) %>%
  add_pulse_event(time = 11, comp = "NH4", unmarked = 0, marked = -0.00219, which = 2) %>%
  add_pulse_event(time = 11, comp = "NO3", unmarked = 0, marked = -0.021, which = 2) %>%
  add_pulse_event(time = 11, comp = "NH4", unmarked = 0, marked = -0.000752, which = 3) %>%
  add_pulse_event(time = 11, comp = "NO3", unmarked = 0, marked = -0.0216, which = 3) %>%
  add_pulse_event(time = 11, comp = "NH4", unmarked = 0, marked = -0.00569, which = 4) %>%
  add_pulse_event(time = 11, comp = "NO3", unmarked = 0, marked = -0.00852, which = 4) %>%
  add_pulse_event(time = 11, comp = "NH4", unmarked = 0, marked = -0.00264, which = 5) %>%
  add_pulse_event(time = 11, comp = "NO3", unmarked = 0, marked = -0.0112, which = 5) %>%
  add_pulse_event(time = 11, comp = "NH4", unmarked = 0, marked = -0.000752, which = 6) %>%
  add_pulse_event(time = 11, comp = "NO3", unmarked = 0, marked = -0.0125, which = 6)  

# Covariates
m <- m %>%
    add_covariates(. ~ stream,
                   eta + zeta ~ 1)

# Priors
m <- m %>%
    set_prior(normal_p(0, 5), "lambda", quiet = TRUE) %>%
    set_prior(normal_p(0, 5), "upsilon", quiet = TRUE) %>%
    set_prior(constant_p(0), "lambda_N", quiet = TRUE) %>%
    set_prior(hcauchy_p(150), "upsilon_N", quiet = TRUE) %>%
    set_prior(uniform_p(0, 1), "portion", quiet = TRUE) %>%
    set_prior(normal_p(0, 2), "eta", quiet = TRUE)

### ** Save the model

trini_mod <- m

save(trini_mod, file = file.path(here::here(), "data", "trini_mod.rda"),
     version = 3, compress = "xz")
