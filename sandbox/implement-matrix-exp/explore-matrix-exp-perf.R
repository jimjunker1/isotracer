### * Description

# Compare the performances of matrix_exp and euler methods.

### * Setup

library(tidyverse)
library(patchwork)
library(etalon)
devtools::load_all()
library(rstan)
options(mc.cores = parallel::detectCores())
library(future)
plan(multicore(workers = 6))

### * Functions

### ** random_topo()

#' @param n Integer, number of compartments.
random_topo <- function(n) {
    stopifnot(n <= 26)
    z <- igraph::as_edgelist(igraph::sample_pa(n = n)) %>%
        as.data.frame()
    all_comps <- unique(unlist(z))
    sources <- unique(z$V2[!z$V2 %in% z$V1])
    non_sources <- all_comps[!all_comps %in% sources]
    n_edges <- nrow(z)
    n_new_edges <- ceiling((n_edges)^(1/4))
    new_edges <- lapply(seq_len(n_new_edges), function(i) {
        c(V1 = sample(non_sources, 1), V2 = sample(non_sources, 1))
    })
    z <- rbind(z, bind_rows(new_edges))
    z <- z[, c("V2", "V1")]
    colnames(z) <- c("from", "to")
    z <- z %>% filter(from != to)
    z$from <- letters[z$from]
    z$to <- letters[z$to]
    z <- unique(z)
    z <- z[order(z$from, z$to), ]
    return(z)
}

### ** random_network()

#' @param n Integer, number of compartments.
#' @param p_steady For source compartments, probability of being set to steady
#'     state.
#' @param p_split For non-source compartments, probability of being set to
#'     split compartment.
random_network <- function(n, p_steady = 0.2, p_split = 0.1) {
    topo <- random_topo(n = n)
    comps <- unique(unlist(topo))
    sources <- unique(topo$from[!topo$from %in% topo$to])
    non_sources <- comps[!comps %in% sources]
    nm <- new_networkModel() %>%
        set_topo(topo, from = "from", to = "to")
    steady <- runif(length(sources)) < p_steady
    for (i in which(steady)) {
        nm <- set_steady(nm, comps = sources[i])
    }
    split <- runif(length(non_sources)) < p_split
    for (i in which(split)) {
        nm <- set_split(nm, comps = non_sources[i])
    }
    return(nm)
}

### ** random_init_obs()

#' Add init quantities and random observations to a model
init_random_obs <- function(nm, n_obs_times = 10) {
    comps <- comps(nm)[[1]]
    sources <- apply(topo(nm), 1, sum)
    sources <- names(sources)[which(sources == 0)]
    init <- tibble(comp = comps,
                   size = 10,
                   prop = 0.00366)
    init[["prop"]][init[["comp"]] %in% sources] <- 0.1
    nm <- set_init(nm, init, comp = "comp", size = "size", prop = "prop")
    stopifnot(length(sources)>0)
    nm <- add_pulse_event(nm, time = 5, comp = sources[1], unmarked = 0, marked = runif(1, 0, 2))
    p <- params(nm)
    p <- setNames(runif(length(p), 0, 0.1), nm = p)
    nm <- set_params(nm, p)
    z <- project(nm, end = 10, grid_size = 1024)
    obs <- sample_from(z, at = seq(0, 10, length.out = n_obs_times))
    nm <- set_obs(nm, obs, time = "time")
    nm
}

### ** random_model()

#' Generate a random network and add random observations
#' @return A network model ready to be run with \code{run_mcmc}.
random_model <- function(n_comps, n_obs_times) {
    nm <- random_network(n = n_comps)
    nm <- init_random_obs(nm, n_obs_times = n_obs_times)
    nm
}

### * Main

## ### ** Quick test

## m <- random_model(n_comps = 6, n_obs_times = 10)
## #plot(m)
## system.time({
##     r <- run_mcmc(m, iter = 1000, chains = 3)
## })
## #plot(r)

## ### ** Fitting functions

## run_time_euler <- function(dataset) {
##     nm <- dataset
##     system.time({run_mcmc(model = nm,
##                           iter = 1000,
##                           chains = 1,
##                           method = "euler")})[["elapsed"]]
## }

## run_time_matrix_exp <- function(dataset) {
##     nm <- dataset
##     system.time({run_mcmc(model = nm,
##                           iter = 1000,
##                           chains = 1,
##                           method = "matrix_exp")})[["elapsed"]]
## }

## ### ** Scaling with number of compartments

## z <- simtable() %>%
##     set_parameters(n_comps = c(3, 6, 9, 12),
##                    n_obs_times = c(3, 7, 13, 15)) %>%
##     replicate_sims(2) %>%
##     generate_datasets(random_model)
## z <- z[sample(seq_len(nrow(z))), ]
## z$id <- replicate(nrow(z), uuid())

## i <- which(etalon::params(z)[["n_comps"]] > 9 | etalon::params(z)[["n_obs_times"]] > 13)
## z <- z[i, ]

## z <- z %>%
##     fit_model(run_time_euler, method = "euler") %>%
##     fit_model(run_time_matrix_exp, method = "matrix_exp")
## z$run_time <- unlist(z$fit)
## w3 <- bind_cols(etalon::params(z), z)

## ggplot(w, aes(x = n_comps, y = run_time)) +
##     geom_jitter(width = 0.1, aes(col = method_fit)) 

## w <- bind_rows(w, w3)

## a <- w %>% filter(method_fit == "euler") %>% rename(euler = run_time)
## b <- w %>% filter(method_fit == "matrix_exp") %>% rename(matrix_exp = run_time)
## x <- left_join(a %>% select(euler, id), b %>% select(matrix_exp, id, n_comps, n_obs_times))
## x$ratio <- x$matrix_exp / x$euler

## m <- lm(log(ratio) ~ n_comps + n_obs_times, data = x)

## # Adding one compartment increases the run time ratio by about 21%.
## # Adding one obs time increases the run time ratio by about 16%.

## (p0 <- ggplot(x, aes(x = n_comps, y = ratio)) +
##     geom_point(aes(col = as.factor(n_obs_times))) +
##     scale_y_log10())

## (p1 <- ggplot(x, aes(x = n_obs_times, y = ratio)) +
##     geom_point(aes(col = as.factor(n_comps))) +
##     scale_y_log10() +
##     ylab("Ratio of run times (matrix exp)/(Euler)"))

## pdf("comparison-runtimes.pdf", width = 6, height = 4)
## print(p0)
## print(p1)
## dev.off()

## uuid <- function() {
##     paste0(sample(c(0:9, letters), 16, TRUE), collapse = "")
## }

## saveRDS(w, "toto.rds")
## w <- readRDS("toto.rds")

## #quit()

## ### ** Many tests

## # max n_obs_times 25?
## # max n_comps 10?
## z <- simtable() %>%
##     set_parameters(n_comps = c(3),
##                    n_obs_times = c(3, 5, 7)) %>%
##     replicate_sims(3)

## # There is a bug somewhere: sometimes the line below crashes?
## z <- z %>% generate_datasets(random_model)


## z2 <- z %>% fit_model(run_time_matrix_exp, method = "matrix_exp", .sequential = TRUE)
## z <- z %>% fit_model(run_time_euler, method = "euler", .sequential = TRUE)

## w <- z[, 1:2]
## w[["euler"]] <- unlist(z[["fit"]])
## w[["matrix_exp"]] <- unlist(z2[["fit"]])
## w <- bind_cols(etalon::params(w), w)

## plot(w$euler, w$matrix_exp, asp = 1, las = 1)
## abline(0, 1)

## ggplot(w, aes(x = euler, y = matrix_exp)) +
##     geom_point(aes(col = n_obs_times)) +
##     geom_abline(intercept = 0, slope = 1) +
##     coord_fixed()

## ggplot(w %>% select(-dataset, -parameters) %>%
##          pivot_longer(names_to = "method", - c("n_obs_times", "n_comps")),
##        aes(x = n_obs_times, y = value)) +
##     geom_jitter(aes(col = method), width = 0.1)

## ### * Comparing methods for trini data

## m <- trini_mod[1, ] %>%
##     add_covariates(. ~ 1)

## system.time({
##     re <- run_mcmc(m, chains = 1, iter = 500, refresh = 1,
##                    control = list(adapt_delta = 0.7, max_treedepth = 5))
## })

## system.time({
##     re <- run_mcmc(m, chains = 1, iter = 500, refresh = 1,
##                    method = "euler")
## })

### * Comparing run times with different Stan tunings

adapt_delta <- c(0.6, 0.7, 0.8, 0.9)
max_treedepth <- c(7, 8)
n_repl <- 3

### ** Prepare the model

m <- trini_mod[1, ] %>%
    add_covariates(. ~ 1) %>%
    set_size_family("normal_sd", by_compartment = TRUE) %>%
    set_prior(constant_p(0.1), "zeta_NH4") %>%  # Dummy values for the steady
    set_prior(constant_p(0.1), "zeta_NO3") %>%        # state dissolved nutrients
    set_prior(constant_p(115), "zeta_epi") %>%
    set_prior(constant_p(2436), "zeta_FBOM") %>%
    set_prior(constant_p(10.7), "zeta_pseph") %>%
    set_prior(constant_p(2.88), "zeta_arg") %>%
    set_prior(constant_p(13.3), "zeta_tricor") %>%
    set_prior(constant_p(2.48), "zeta_petro") %>%
    set_prior(hcauchy_p(5), "upsilon_N") %>%
    set_prior(constant_p(0), "lambda_N")

### ** Run the simulations

z <- simtable() %>%
    set_parameters(adapt_delta = adapt_delta,
                   max_treedepth = max_treedepth) %>%
    replicate_sims(n_repl)
d <- function(adapt_delta, max_treedepth) {
    list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
}
z <- z %>% generate_datasets(d)
z <- z[sample(seq_len(nrow(z))), ]
f <- function(dataset) {
    t <- system.time({ r <- run_mcmc(m, iter = 500, chains = 1,
                                     cores = 1, refresh = 10,
                                     control = dataset,
                                     stanfit = TRUE) })
    return(list(fit = r, time = t))
}
z <- z %>% fit_model(f)
old_z <- readRDS("toto-z.rds")
z <- bind_rows(z, old_z)
saveRDS(z, file = "toto-z.rds")

### ** Parse the simulations

z <- readRDS("toto-z.rds")

z <- z %>%
    mutate(time = sapply(fit, function(i) i$time[["elapsed"]]),
           n_divergent = sapply(fit, function(i) rstan::get_num_divergent(i$fit)),
           n_max_treedepth = sapply(fit, function(i) rstan::get_num_max_treedepth(i$fit)))

w <- bind_cols(etalon::params(z), z)
w <- arrange(w, adapt_delta, max_treedepth)

p_time1 <- ggplot(w, aes(x = adapt_delta, y = time, col = as.factor(max_treedepth))) +
    geom_jitter(width = 0.005) +
    scale_y_log10()
p_time2 <- ggplot(w, aes(x = max_treedepth, y = time, col = as.factor(adapt_delta))) +
    geom_jitter(width = 0.01) +
    scale_y_log10()
p_div1 <-ggplot(w, aes(x = adapt_delta, y = n_divergent, col = as.factor(max_treedepth))) +
    geom_jitter(width = 0.005)
p_div2 <- ggplot(w, aes(x = max_treedepth, y = n_divergent, col = as.factor(adapt_delta))) +
    geom_jitter(width = 0.01)
p_mtd1 <- ggplot(w, aes(x = adapt_delta, y = n_max_treedepth, col = as.factor(max_treedepth))) +
    geom_jitter(width = 0.005)
p_mtd2 <- ggplot(w, aes(x = max_treedepth, y = n_max_treedepth, col = as.factor(adapt_delta))) +
    geom_jitter(width = 0.01)

(p_time1 + p_time2) / (p_div1 + p_div2) / (p_mtd1 + p_mtd2)
