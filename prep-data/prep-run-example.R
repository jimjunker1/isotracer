### * Description

# Prepare a network model that can be run at once

# This script builds the network model that is described in
# vignettes/tutorial-010-quick-start.Rmd. The model is stored in the object
# "aquarium_mod".

### * Setup

suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(devtools::load_all())
set.seed(8)
options(mc.cores = parallel::detectCores())

### * Aquarium model

exp <- tibble::tribble(
  ~time.day,  ~species, ~biomass, ~prop15N,
          0,   "algae",     1.02,  0.00384,
          1,   "algae",       NA,   0.0534,
        1.5,   "algae",    0.951,       NA,
          2,   "algae",    0.889,   0.0849,
        2.5,   "algae",       NA,   0.0869,
          3,   "algae",    0.837,   0.0816,
          0, "daphnia",     1.74,  0.00464,
          1, "daphnia",       NA,  0.00493,
        1.5, "daphnia",     2.48,       NA,
          2, "daphnia",       NA,  0.00831,
        2.5, "daphnia",     2.25,       NA,
          3, "daphnia",     2.15,   0.0101,
          0,     "NH4",    0.208,     0.79,
          1,     "NH4",    0.227,       NA,
        1.5,     "NH4",       NA,    0.482,
          2,     "NH4",    0.256,    0.351,
        2.5,     "NH4",       NA,    0.295,
          3,     "NH4",     0.27,        NA
  )
inits <- exp %>% filter(time.day == 0)
obs <- exp %>% filter(time.day > 0)

m <- new_networkModel(quiet = TRUE) %>%
    set_topo("NH4 -> algae -> daphnia -> NH4") %>%
    set_init(inits, comp = "species", size = "biomass",
             prop = "prop15N") %>%
    set_obs(obs, comp = "species", size = "biomass",
            prop = "prop15N", time = "time.day") %>%
    set_priors(normal_p(0, 0.5), "", quiet = TRUE) %>%
    set_priors(normal_p(0, 0.2), "lambda", quiet = TRUE)

aquarium_mod <- m

save(aquarium_mod, file = file.path(here::here(), "data", "aquarium_mod.rda"),
     version = 3, compress = "xz")

### * Run model

cache_file <- file.path(here::here(), "prep-data", "z-cache-prep-run-example.rds")
if (!file.exists(cache_file)) {
    r <- run_mcmc(m, seed = 40, control = list(adapt_delta = 0.9), thin = 4)
    saveRDS(r, cache_file)
} else {
    r <- readRDS(cache_file)
}

aquarium_run <- r

save(aquarium_run, file = file.path(here::here(), "data", "aquarium_run.rda"),
     version = 3, compress = "xz")
