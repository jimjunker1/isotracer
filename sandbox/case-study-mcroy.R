### * Description

# Running the case study for McRoy 1970 (eelgrasses)

### * Setup

library(isotracer)
library(tidyverse)
library(magrittr)

### * Run

d <- eelgrass

x <- new_NetworkModel() %>%
    set_topo("upperWater -> leavesAndStem, rootsAndRhizome") %>%
    set_init(d %>% filter(time.min == 1), comp = "compartment",
            size = "phosphorus.n", prop = "proportion") %>%
    add_obs(d %>% filter(time.min > 1), comp = "compartment",
           size = "phosphorus.n", prop = "proportion", time = "time.min") %>%
    addParamMapping() %>%
    set_prior(uniform_p(0, 0.01), ".") %>%
    set_prior(uniform_p(0, 10), "zeta") %>%
    set_prior(uniform_p(0, 10), "eta", FALSE)

priors(x)
priors(x)$prior

z <- mugen_stan(x, chains = 4, iter = 1000)
plotTraces(z, drawHist = TRUE, nBars = 32)
