### * Setup

library(tidyverse)
library(magrittr)
library(isotracer)

devtools::load_all()
set.seed(4)

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

### * Test simulations

### ** Loop

x <- new_NetworkModel() %>%
    set_topo(c("a -> b -> c", "c -> a")) %>%
    set_init(tibble(comp = c("a", "b", "c"),
                   size = c(100, 50, 200),
                   prop = c(1, 0, 0.5)),
            comp = "comp", size = "size", prop = "prop")
params(x)

y <- x %>% set_params(c("eta" = 1, "lambda_a" = 0, "lambda_b" = 0, "lambda_c" = 0,
                       "upsilon_a_to_b" = 0.1, "upsilon_b_to_c" = 0.2,
                       "upsilon_c_to_a" = 0.2, "zeta" = 1))
y <- y %>% project(gridSize = 512, end = 30)

plotTraj(y)

### ** Branched with pulses

x <- new_NetworkModel() %>%
    set_topo(c("a -> b -> c", "a -> c -> d", "d -> a")) %>%
    set_init(tibble(comp = c("a", "b", "c", "d"),
                   size = c(100, 50, 200, 10),
                   prop = c(0.1, 0.1, 0.1, 0.1)),
            comp = "comp", size = "size", prop = "prop")
params(x)

y <- x %>% set_params(c("eta" = 1, "lambda_a" = 0, "lambda_b" = 0, "lambda_c" = 0,
                       "lambda_d" = 0,
                       "upsilon_a_to_b" = 0.1, "upsilon_b_to_c" = 0.15,
                       "upsilon_c_to_d" = 0.1, "upsilon_d_to_a" = 0.3,
                       "upsilon_a_to_c" = 0.2, "zeta" = 1)) %>%
    set_steady("a") %>%
    add_pulse_event(5, "a", -50, 50) %>%
    add_pulse_event(15, "a", 50, -50)

y <- y %>% project(gridSize = 512, end = 30)

plotTraj(y)

### ** Linear with drip

x <- new_NetworkModel() %>%
    set_topo("a -> b -> c -> d -> e") %>%
    set_init(tibble(comps = c("a", "b", "c", "d", "e"),
                   sizes = rep(100, 5),
                   props = rep(0, 5)),
            comp = "comps", size = "sizes", prop = "props") %>%
    set_steady("a")
params(x)
x <- x %>% set_params(c("eta" = 1, "zeta" = 1,
                       "lambda_a" = 0, "lambda_b" = 0, "lambda_c" = 0,
                       "lambda_d" = 0, "lambda_e" = 0.1,
                       "upsilon_a_to_b" = 0.1, "upsilon_b_to_c" = 0.1,
                       "upsilon_c_to_d" = 0.1, "upsilon_d_to_e" = 0.1)) %>%
    add_pulse_event(10, "a", -50, 50) %>%
    add_pulse_event(20, "a", 50, -50)

y <- x %>% project(end = 80)

plotTraj(y)

### * Test ABC

d <- quickABC_01()
x <- new_NetworkModel() %>%
    set_topo("A -> B,C") %>%
    set_init(d %>% filter(time == 0),
            comp = "compartment", size = "biomass",
            prop = "proportion", group_by = "loc") %>%
    add_obs(d %>% filter(time > 0),
            comp = "compartment", size = "biomass",
           prop = "proportion", time = "time",
           group_by = "loc")
#x$topology[[3]]["C", "B"] <- 1
x %<>% set_prior(uniform_p(0, 1), ".") %>%
    set_prior(hcauchy_p(1), "eta", FALSE) %>%
    set_steady("A") %>%
    add_pulse_event(2, "A", unmarked = 0, marked = 200) %>%
    add_pulse_event(5, "A", unmarked = 0, marked = -200) %>%
    add_covariates(upsilon ~ loc)
params(x)

z <- mugen_stan(x, iter = 1000)
plotTraces(z, drawHist = TRUE, nBars = 32)

p <- summary(z)$statistics
p <- setNames(p[, "Mean"], rownames(p))

y <- x %>%
    set_params(p) %>%
    project()

## d <- prep_stan_data(x)
## fit <- stan(file = file.path(here::here(), "inst", "stan",
##                                   "networkModelMugen.stan"),
##             data = d, chains = 4, iter = 500,
##             pars = c("params"))
## z <- As.mcmc.list(fit)
## rawNames <- coda::varnames(z)
## paramNames <- rawNames[startsWith(rawNames, "params[")]
## z <- z[, paramNames]
## coda::varnames(z) <- d[["allParams"]]
## plotTraces(z)

### * Test projections

y <- new_NetworkModel() %>%
    set_topo("A -> B,C") %>%
    set_init(tibble(comp = c("A", "B", "C"), size = c(100, 50, 80),
                   prop = c(1, 0, 0)),
            comp = "comp", size = "size", prop = "prop") %>%
    set_params(c(eta = 0.1, lambda_A = 0.1, lambda_B = 0.05, lambda_C = 0.02,
                upsilon_A_to_B = 0.02, upsilon_A_to_C = 0.04, zeta = 0.2)) %>%
    project(end = 10)

### * Tests

### ** Data

links <- c("grass -> cow -> wolf")
inits <- tribble(~ species, ~biomass, ~proportion, ~location,
                 "grass", 10, 0.1, "north", 
                 "cow", 20, 0.07, "north", 
                 "wolf", 5, 0.05, "north", 
                 "grass", 17, 0.1, "south",
                 "cow", 23, 0.06, "south", 
                 "wolf", 3, 0.03, "south")
obs <- tribble(~ species, ~day, ~biomass, ~proportion, ~location, ~year,
               "grass", 1, 10, 0.1, "north", 2014,
               "cow", 1, 20, 0.07, "north", 2014,
               "wolf", 1, 5, 0.05, "north", 2014,
               "grass", 1, 17, 0.1, "south", 2014,
               "cow", 1, 23, 0.06, "south", 2014,
               "wolf", 1, 3, 0.03, "south", 2014,
               "grass", 1, 11, 0.08, "north", 2015,
               "cow", 1, 22, 0.12, "north", 2015,
               "wolf", 5, 7, 0.1, "north", 2015,
               "grass", 5, 10, 0.1, "north", 2014,
               "cow", 5, 20, 0.07, "north", 2014,
               "wolf", 5, 5, 0.05, "north", 2014,
               "grass", 5, 17, 0.1, "south", 2014,
               "cow", 5, 23, 0.06, "south", 2014,
               "wolf", 5, 3, 0.03, "south", 2014,
               "grass", 5, 11, 0.08, "north", 2015,
               "cow", 5, 22, 0.12, "north", 2015,
               "wolf", 5, 7, 0.1, "north", 2015)

### ** Run

x <- new_NetworkModel() %>%
    set_topo(links) %>%
    set_init(data = inits,
            comp = "species",
            size = "biomass",
            prop = "proportion",
            group_by = c("location")) %>%
    add_obs(data =obs,
           comp = "species", size = "biomass",
           prop = "proportion", time = "day",
           group_by = c("location", "year"))

### ** Run with Collins2016

x <- newNetworkModel() %>%
    set_topo(links = "NH4 -> epi -> pseph")
x

# Initial conditions
d <- datasetObs("collins2016") %>%
    filter(stream == "UPL",
           compartment %in% c("NH4", "epi", "pseph")) %>%
    group_by(compartment) %>%
    summarize(init_size = mean(mgN.per.m2, na.rm = TRUE))
d$init_prop <- delta2prop(0, "d15N")

x %<>% set_init(data = d,
               comp = "compartment",
               size = "init_size",
               prop = "init_prop")

# Observations
d <- datasetObs("collins2016") %>%
    filter(stream == "UPL",
           compartment %in% c("NH4", "epi", "pseph")) %>%
    select(transect, compartment, time.days, d15N) %>%
    mutate(prop = delta2prop(d15N, "d15N")) %>%
    na.omit()
d$size <- NA

x %<>% add_obs(data = d,
              comp = "compartment",
              size = "size",
              prop = "prop",
              time = "time.days",
              group_by = "transect")

