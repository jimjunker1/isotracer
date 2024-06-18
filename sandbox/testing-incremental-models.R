### * Description

# Test models of incremental complexity to check that the R and Stan codes work
# properly

### * Setup

library(tidyverse)
library(magrittr)
library(isotracer)
devtools::load_all()
library(future)
plan(multiprocess)

### * Model: two compartments

m <- new_networkModel() %>%
    set_topo("A -> B") %>%
    set_init(tibble(comp = c("A", "B"),
                    size = c(100, 100),
                    prop = c(0.5, 0)),
             comp = "comp", size = "size", prop = "prop")
params(m)
m %<>% set_params(c("upsilon_A_to_B" = 0.1,
                    "eta" = 0.01, "zeta" = 0.01,
                    "lambda_A" = 0, "lambda_B" = 0)) %>%
    project(end = 50)
plot(m)

obs <- m %>% sample_from(at = 1:50)

m %<>% set_obs(obs, comp = "comp", size = "size", prop = "prop",
               time = "time")
plot(m)

f <- run_mcmc(m, iter = 500)
plot(f)

p <- predict(m, f)
plot(p)

### * Model: two compartments, one split

m <- new_networkModel() %>%
    set_topo("A -> B") %>%
    set_split("B") %>%
    set_init(tibble(comp = c("A", "B"),
                    size = c(100, 100),
                    prop = c(0.5, 0)),
             comp = "comp", size = "size", prop = "prop")
params(m)

m %<>% set_params(c("upsilon_A_to_B" = 0.1,
                    "eta" = 0.01, "zeta" = 0.01,
                    "lambda_A" = 0, "lambda_B" = 0,
                    "portion.act_B" = 0.5)) %>%
    project(end = 50)
plot(m)

obs <- m %>% sample_from(at = 1:50)

m %<>% set_obs(obs, comp = "comp", size = "size", prop = "prop",
               time = "time")
plot(m)

f <- run_mcmc(m, iter = 500)
plot(f)

p <- predict(m, f)
plot(p)

### * Model: three compartments

m <- new_networkModel() %>%
    set_topo("A -> B -> C") %>%
    set_init(tibble(comp = c("A", "B", "C"),
                    size = c(100, 100, 100),
                    prop = c(0.5, 0, 0)),
             comp = "comp", size = "size", prop = "prop")
params(m)
m %<>% set_params(c("upsilon_A_to_B" = 0.1,
                    "upsilon_B_to_C" = 0.1,
                    "eta" = 0.01, "zeta" = 0.01,
                    "lambda_A" = 0, "lambda_B" = 0,
                    "lambda_C" = 0)) %>%
    project(end = 50)
plot(m)

obs <- m %>% getSamples(at = 1:50)

m %<>% set_obs(obs, comp = "comp", size = "size", prop = "prop",
               time = "time")
plot(m)

f <- run_mcmc(m, iter = 500)
plot(f)

p <- predict(m, f)
plot(p)

### * Model: three compartments, one split

m <- new_networkModel() %>%
    set_topo("A -> B -> C") %>%
    set_split("B") %>%
    set_init(tibble(comp = c("A", "B", "C"),
                    size = c(100, 100, 100),
                    prop = c(0.5, 0, 0)),
             comp = "comp", size = "size", prop = "prop")
params(m)

m %<>% set_params(c("upsilon_A_to_B" = 0.1,
                    "upsilon_B_to_C" = 0.1,
                    "eta" = 0.01, "zeta" = 0.01,
                    "lambda_A" = 0, "lambda_B" = 0,
                    "lambda_C" = 0,
                    "portion.act_B" = 0.5)) %>%
    project(end = 50)
plot(m)

obs <- m %>% getSamples(at = (1:5)*10)

m %<>% set_obs(obs, comp = "comp", size = "size", prop = "prop",
               time = "time")
plot(m)

f <- run_mcmc(m, iter = 500)
plot(f)

p <- predict(m, f)
plot(p)

## Closed forms

q_a <- function(t, a0) {
    a0 * exp(- 0.1 * t)
}

q_b <- function(t, a0, b0) {
    (0.1 * a0 * t + b0) * exp(-0.1 * t)
}

q_c <- function(t, a0, b0, c0) {
    (-0.1 * a0 * t - a0 - b0) * exp(-0.1 * t) + c0 + a0 + b0
}

at <- seq(0, 50, length.out = 128)
plot(at, q_b(at, 100, 50) + 50, type = "l")
