### * Description

# Implementing and testing Stan model with matrix exponential (in reply to
# reviews from submission to MEE, [2021-10-08 Fri])

# Solving ODE in R:
# - See the deSolve package
# - https://kinglab.eeb.lsa.umich.edu/480/nls/de.html

# Solving ODE in Stan:
# - https://mc-stan.org/users/documentation/case-studies/lotka-volterra-predator-prey.html

# Matrix exponential in R:
# - package expm, function expm
# - package Matrix, function expm
# (Some quick benchmarking suggests that expm::expm() is faster than Matrix::expm().)

# Profiling a Stan model:
# - https://mc-stan.org/cmdstanr/articles/profiling.html

### * Setup

library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
devtools::load_all()

source("helpers.R")

### * Functions

parse_stan_fit <- function(fit, stan.data) {
    # Get mcpars
    start <- fit@sim[["warmup"]] + 1
    end <- fit@sim[["iter"]]
    thin <- 1
    n_kept <- end - start + 1
    mcpars <- c(start, end, thin)
    # Prepare the mcmc.list object
    out <- rstan::As.mcmc.list(fit)
    for (i in seq_along(out)) {
        stopifnot(nrow(out[[i]]) == n_kept)
    }
    attr(out, "mcpar") <- mcpars
    rawNames <- coda::varnames(out)
    nonConstantParamNames <- rawNames[startsWith(rawNames, "nonConstantParams[")]
    loglikNames <- rawNames[startsWith(rawNames, "log_lik[")]
    ll <- out[, loglikNames]
    out <- out[, nonConstantParamNames]
    coda::varnames(out) <- stan.data[["allParams"]][stan.data[["mappingParamPriorType"]] != 0]
    llTrace <- lapply(ll, function(x) {
        out <- coda::as.mcmc(apply(as.matrix(x), 1, sum))
        attr(out, "mcpar") <- attr(x, "mcpar")
        return(out)
    })
    llTrace <- coda::as.mcmc.list(llTrace)
    attr(out, "loglik") <- llTrace
    attr(out, "mcpar") <- mcpars
    # Add values of constant parameters (if any)
    n_constant_params <- sum(stan.data[["mappingParamPriorType"]] == 0)
    if (n_constant_params > 0) {
        constant_params <- stan.data[["allParams"]][stan.data[["mappingParamPriorType"]] == 0]
        constant_values <- stan.data[["constantParams"]][stan.data[["mappingParamPriorType"]] == 0]
        constant_params <- setNames(constant_values, nm = constant_params)
        attr(out, "constant_params") <- constant_params
    }
    # Return the mcmc.list object
    return(out)
}

### * Main

### ** aquarium_mod

m <- aquarium_mod

# Reference run
system.time({
    r <- run_mcmc(m, iter = 1000)
})
##   user  system elapsed 
## 72.811   0.703  24.766
plot(r)

# Run with matrix exp model
d <- prep_stan_data_expm(m)

model <- stan_model(here::here("sandbox", "implement-matrix-exp",
                               "networkModel-matrix-exp.stan"))

## CAUTION: marked and unmarked arrays in stan model always starts at t0, which
## might be different from the first obs time, so there might be some mismatch
## between time indices of observations and of those arrays.

system.time({
    rm <- sampling(model, data = d, iter = 1000, chains = 4,
                   pars = c("nonConstantParams", "log_lik"),
                   control = list(adapt_delta = 0.95))
})
##   user  system elapsed 
## 14.304   0.397   4.473

z <- parse_stan_fit(rm, d)
class(z) <- unique(c("networkModelStanfit", class(z)))
plot(z)

ks.test(as.matrix(r)[, "zeta"], as.matrix(z)[, "zeta"])

### ** trini_mod

d_trini <- prep_stan_data_expm(trini_mod)
system.time({
    r_trini <- sampling(model, data = d_trini, iter = 50, chains = 1,
                        pars = c("nonConstantParams", "log_lik"),
                        control = list(adapt_delta = 0.95))
})

z <- parse_stan_fit(r_trini, d_trini)
class(z) <- unique(c("networkModelStanfit", class(z)))
plot(z)

### * Models from tutorials

### ** https://matthieu-bruneaux.gitlab.io/isotracer/articles/tutorial-040-pulse-drip-events.html

exp <- tibble::tribble(
  ~time.day, ~species, ~biomass, ~prop15N, ~aquariumID,
          0,    "NH4",   0.2054,   0.0045,      "aq01",
          2,    "NH4",       NA,   0.6815,      "aq01",
          4,    "NH4",   0.3883,       NA,      "aq01",
          6,    "NH4",   0.2804,       NA,      "aq01",
          8,    "NH4",   0.1604,   0.5817,      "aq01",
         10,    "NH4",       NA,   0.5216,      "aq01",
          0,  "algae",   0.8691,   0.0043,      "aq01",
          2,  "algae",   1.1449,   0.0029,      "aq01",
          4,  "algae",       NA,   0.1293,      "aq01",
          6,  "algae",   1.2182,       NA,      "aq01",
          8,  "algae",       NA,   0.1507,      "aq01",
         10,  "algae",   1.3169,   0.2549,      "aq01",
          0,    "NH4",   0.5637,   0.0035,      "aq02",
          2,    "NH4",   0.7241,   0.4896,      "aq02",
          6,    "NH4",       NA,   0.5623,      "aq02",
          8,    "NH4",       NA,   0.4578,      "aq02",
         10,    "NH4",   0.1515,   0.5202,      "aq02",
          0,  "algae",   1.4809,   0.0039,      "aq02",
          2,  "algae",   1.1758,   0.0041,      "aq02",
          4,  "algae",       NA,   0.0833,      "aq02",
          6,  "algae",   1.7349,    0.105,      "aq02",
          8,  "algae",       NA,   0.1899,      "aq02",
         10,  "algae",   2.2429,   0.1281,      "aq02"
  )
m <- new_networkModel() %>% set_topo("NH4 -> algae")
inits <- exp %>% filter(time.day == 0)
m <- m %>% set_init(inits, comp = "species", size = "biomass", prop = "prop15N",
                    group_by = "aquariumID")
obs <- exp %>% filter(time.day > 0)
m <- m %>% set_obs(obs, time = "time.day")
m <- m %>% add_pulse_event(time = 2, comp = "NH4", unmarked = 0, marked = 0.4)

d <- prep_stan_data_expm(m)
system.time({
    r <- sampling(model, data = d, iter = 2000, chains = 4,
                  pars = c("nonConstantParams", "log_lik"),
                  control = list(adapt_delta = 0.95))
})

z <- parse_stan_fit(r, d)
class(z) <- unique(c("networkModelStanfit", class(z)))
plot(z)

# Drip event

exp <- tibble::tribble(
  ~time.day,    ~species, ~biomass, ~prop15N,    ~transect,
          0,       "NH4",    0.313,   0.0038, "transect_1",
          4,       "NH4",   0.2746,       NA, "transect_1",
          8,       "NH4",   0.3629,   0.0295, "transect_1",
         12,       "NH4",       NA,   0.0032, "transect_1",
         16,       "NH4",       NA,   0.0036, "transect_1",
         20,       "NH4",   0.3414,   0.0038, "transect_1",
          0, "epilithon",  89.2501,   0.0022, "transect_1",
          8, "epilithon", 123.1212,    0.024, "transect_1",
         12, "epilithon",       NA,     0.02, "transect_1",
         16, "epilithon",  90.5919,   0.0107, "transect_1",
         20, "epilithon",  80.3261,   0.0062, "transect_1",
          0,       "NH4",   0.3525,   0.0035, "transect_2",
          4,       "NH4",   0.2958,   0.0362, "transect_2",
          8,       "NH4",       NA,     0.03, "transect_2",
         12,       "NH4",   0.3392,   0.0044, "transect_2",
         16,       "NH4",    0.212,   0.0026, "transect_2",
         20,       "NH4",   0.3818,   0.0046, "transect_2",
          0, "epilithon",  127.873,   0.0029, "transect_2",
          4, "epilithon",       NA,   0.0096, "transect_2",
          8, "epilithon", 123.3216,       NA, "transect_2",
         12, "epilithon",  89.8053,   0.0144, "transect_2",
         16, "epilithon",  74.9105,   0.0098, "transect_2",
         20, "epilithon",  88.0108,   0.0067, "transect_2"
  )
m <- new_networkModel() %>% set_topo("NH4 -> epilithon")
inits <- exp %>% filter(time.day == 0)
m <- m %>% set_init(inits, comp = "species", size = "biomass", prop = "prop15N",
                    group_by = "transect")
obs <- exp %>% filter(time.day > 0)
m <- m %>% set_obs(obs, time = "time.day")
m <- m %>% set_steady(comps = "NH4")
m <- m %>%
    add_pulse_event(time = 2, comp = "NH4", unmarked = 0, marked = 0.008)
m <- m %>%
    add_pulse_event(time = 10, comp = "NH4", unmarked = 0, marked = -0.008)

d <- prep_stan_data_expm(m)
system.time({
    r <- sampling(model, data = d, iter = 2000, chains = 4,
                  pars = c("nonConstantParams", "log_lik"),
                  control = list(adapt_delta = 0.95))
})

z <- parse_stan_fit(r, d)
class(z) <- unique(c("networkModelStanfit", class(z)))
plot(z)

### ** https://matthieu-bruneaux.gitlab.io/isotracer/articles/tutorial-050-fixed-effects.html

exp <- tibble::tribble(
  ~time.day,  ~species, ~biomass, ~prop15N, ~treatment,
          0,     "NH4",    0.205,    0.739,    "light",
          2,     "NH4",    0.232,    0.403,    "light",
          4,     "NH4",       NA,    0.199,    "light",
          6,     "NH4",       NA,    0.136,    "light",
          8,     "NH4",    0.306,       NA,    "light",
         10,     "NH4",    0.323,   0.0506,    "light",
          0,   "algae",    0.869,  0.00305,    "light",
          2,   "algae",       NA,   0.0875,    "light",
          4,   "algae",     0.83,    0.131,    "light",
          6,   "algae",    0.706,       NA,    "light",
         10,   "algae",    0.666,   0.0991,    "light",
          0, "daphnia",     2.13,  0.00415,    "light",
          2, "daphnia",     1.99,       NA,    "light",
          4, "daphnia",     1.97,   0.0122,    "light",
          6, "daphnia",       NA,   0.0284,    "light",
          8, "daphnia",       NA,   0.0439,    "light",
         10, "daphnia",      1.9,   0.0368,    "light",
          0,     "NH4",    0.474,     0.98,     "dark",
          2,     "NH4",    0.455,     0.67,     "dark",
          4,     "NH4",    0.595,    0.405,     "dark",
          6,     "NH4",       NA,    0.422,     "dark",
         10,     "NH4",    0.682,    0.252,     "dark",
          0,   "algae",     1.06,  0.00455,     "dark",
          2,   "algae",        1,   0.0637,     "dark",
          4,   "algae",    0.862,   0.0964,     "dark",
          6,   "algae",       NA,    0.222,     "dark",
          8,   "algae",       NA,    0.171,     "dark",
         10,   "algae",    0.705,    0.182,     "dark",
          0, "daphnia",     1.19,  0.00315,     "dark",
          4, "daphnia",     1.73,   0.0204,     "dark",
          6, "daphnia",     1.75,       NA,     "dark",
          8, "daphnia",     1.54,   0.0598,     "dark",
         10, "daphnia",     1.65,   0.0824,     "dark"
  )
inits <- exp %>% filter(time.day == 0)
obs <- exp %>% filter(time.day > 0)
m <- new_networkModel() %>%
  set_topo("NH4 -> algae -> daphnia -> NH4") %>%
  set_init(inits, comp = "species", size = "biomass", prop = "prop15N",
           group_by = "treatment") %>%
  set_obs(obs, time = "time.day")
m <- m %>% add_covariates(upsilon_NH4_to_algae ~ treatment)
m <- m %>% add_covariates(upsilon_algae_to_daphnia ~ treatment)
m <- m %>% add_covariates(lambda ~ treatment)
m <- m %>% add_covariates(. ~ treatment)
m <- m %>% add_covariates(zeta ~ 1)
m <- m %>% add_covariates(. ~ 1)
m <- m %>% add_covariates(upsilon ~ treatment)

d <- prep_stan_data_expm(m)
system.time({
    r <- sampling(model, data = d, iter = 200, chains = 4,
                  pars = c("nonConstantParams", "log_lik"),
                  control = list(adapt_delta = 0.95))
})

z <- parse_stan_fit(r, d)
class(z) <- unique(c("networkModelStanfit", class(z)))
plot(z)
