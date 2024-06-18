### * Description

# Testing R's tools for ODE solving

### * Setup

library(tidyverse)
devtools::load_all()
library(deSolve)

w <- aquarium_mod
w <- set_params(w, sample_params(w)) %>%
    project(end = 10)
plot(w)

topo(z)

M <- matrix(c(- (0.195 + 0.125), 0, 0.171,
              0.125,  - (0.359 + 0.895), 0,
              0, 0.895, - (0.171 + 0.282)),
            ncol = 3, byrow = TRUE)

mod <- function(t, s, ...) {
    list(M %*% s)
}

y0 <- c(1.02, 1.74, 0.208)

w <- ode(y0, seq(0, 10, length.out = 261), mod)

library(microbenchmark)

microbenchmark(project(z, end = 10),
               ode(y0, seq(0, 10, length.out = 261), mod),
               times = 10)
