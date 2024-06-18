### * Description

# Test modelling with split compartments

### * Setup

library(isotracer)
library(tidyverse)

### * Main

### ** Simulate data

x <- new_networkModel() %>%
    set_topo("NH4 -> algae -> daphnia -> NH4") %>%
    set_split("algae") %>%
    set_init(tibble(comps = c("NH4", "algae", "daphnia"),
                    sizes = c(0.2, 1, 2),
                    props = c(0.8, 0.004, 0.004)),
             comp = "comps", size = "sizes", prop = "props") %>%
    set_params(c("eta" = 0.2, "lambda_algae" = 0, "lambda_daphnia" = 0, 
                 "lambda_NH4" = 0, "upsilon_algae_to_daphnia" = 0.15, 
                 "upsilon_NH4_to_algae" = 0.25, "upsilon_daphnia_to_NH4" = 0.04,
                 "zeta" = 0.1, "portion.act_algae" = 0.1)) %>%
  project(end = 20)

y <- tibble(time = x$trajectory[[1]]$timepoints[[1]]) %>%
    bind_cols(as.data.frame(x$trajectory[[1]]$proportions[[1]])) %>%
    pivot_longer(-time, names_to = "compartment", values_to = "prop")

ggplot(y, aes(time, prop)) +
    geom_line(aes(col = compartment))

z <- sample_from(x, at = seq(0, 20, length.out = 7))

### ** Run model

init <- z %>% filter(time == 0)
obs <- z %>% filter(time > 0)

m <- new_networkModel() %>%
    set_topo("NH4 -> algae -> daphnia -> NH4") %>%
    set_split("algae") %>%
    set_init(init, comp = "comp", size = "size", prop = "prop") %>%
    set_obs(obs, comp = "comp", size = "size", prop = "prop",
            time = "time")

r <- run_mcmc(m, iter = 500)

plot(r)

p <- predict(m, r)

plot(p, facet_row = "type")
