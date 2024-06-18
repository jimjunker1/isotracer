### * Setup

library(magrittr)
library(dplyr)
library(ggplot2)

ITERS <- 53
DRAWS <- 32
FORCE_STAN_COVR <- TRUE

n_cores <- min(2, parallel::detectCores())
n_chains <- max(n_cores, 2)

run_mcmc <- function(...) {
  isotracer:::run_mcmc(..., cores = n_cores, chains = n_chains,
                       seed = 43)
}

new_networkModel <- function() {
  isotracer::new_networkModel(quiet = TRUE)
}

### * McRoy & Barsdate 1970 - Phosphorus uptake in eelgrass

test_that("McRoy & Barsdate case study runs", {

  # Define the model
  
  eelgrass <- eelgrass %>% mutate(prop = 1)
  eelgrass <- eelgrass %>%
    mutate(time_day = time_min / (24 * 60),
           n_1e6_32P = n_32P / 1e6)
  eelgrass <- eelgrass %>%
    select(compartment, time_day, n_1e6_32P, prop, light_treatment, addition_site)
  expect_equal(dim(eelgrass), c(83, 6))
  m <- new_networkModel() %>%
    set_topo("upper_water -> leaves_stem -> roots_rhizome",
             "lower_water -> roots_rhizome -> leaves_stem")
  expect_equal(dim(m), c(1, 4))
  expect_error(ggtopo(m), NA)
  m <- m %>% set_half_life(14.268)
  init <- filter(eelgrass, time_day < 0.01)
  obs <- filter(eelgrass, time_day > 0.01)
  m <- m %>%
    set_init(init, comp = "compartment", size = "n_1e6_32P", prop = "prop",
             group_by = c("light_treatment", "addition_site")) %>%
    set_obs(obs, time = "time_day")
  m <- m %>% add_covariates(upsilon ~ light_treatment)
  expect_equal(dim(params(m)), c(14, 2))
  m <- m %>% set_priors(constant_p(0), "lambda")
  m <- m %>% set_priors(constant_p(1), "^eta")
  m <- m %>%
    set_priors(normal_p(0, 2), "upsilon") %>%
    set_priors(normal_p(0, 2), "zeta")
  expect_equal(dim(m), c(4, 5))
  expect_equal(dim(groups(m)), c(4, 2))
  
  # Run the model

  fit <- purrr::quietly(run_mcmc)(m, iter = ITERS)$result
  expect_error(plot(fit), NA)
  expect_error({
    ggplot(eelgrass, aes(x = time_day, y = n_1e6_32P)) +
      geom_line(aes(col = compartment)) +
      facet_wrap(~ light_treatment + addition_site) +
      coord_trans(y = "log10")
  }, NA)
  pred <- predict(m, fit, probs = 0.95, draws = DRAWS)
  expect_s3_class(pred, "networkModel")
  expect_equal(dim(pred), c(4, 7))
  expect_error({
    plot(pred, facet_row = "group", facet_column = "compartment", log = TRUE, type = "size")
  }, NA)
  expect_equal(varnames(fit),
               c("upsilon_leaves_stem_to_roots_rhizome|dark", "upsilon_leaves_stem_to_roots_rhizome|light", 
                 "upsilon_lower_water_to_roots_rhizome|dark", "upsilon_lower_water_to_roots_rhizome|light", 
                 "upsilon_roots_rhizome_to_leaves_stem|dark", "upsilon_roots_rhizome_to_leaves_stem|light", 
                 "upsilon_upper_water_to_leaves_stem|dark", "upsilon_upper_water_to_leaves_stem|light", 
                 "zeta"))

  # Downstream analyses
  
  z <- (fit[, "upsilon_lower_water_to_roots_rhizome|light"] /
        fit[, "upsilon_upper_water_to_leaves_stem|light"])
  expect_s3_class(z, "derived.mcmc.list")
  expect_length(summary(z)$quantiles, 5)
  z <- (fit[, "upsilon_upper_water_to_leaves_stem|light"] /
        fit[, "upsilon_upper_water_to_leaves_stem|dark"])
  expect_s3_class(z, "derived.mcmc.list")
  expect_length(summary(z)$quantiles, 5)
  z <- (fit[, "upsilon_lower_water_to_roots_rhizome|light"] /
        fit[, "upsilon_lower_water_to_roots_rhizome|dark"])
  expect_s3_class(z, "derived.mcmc.list")
  expect_length(summary(z)$quantiles, 5)
  init2 <- tibble(compartment = c("lower_water", "upper_water", "roots_rhizome",
                                  "leaves_stem"),
                  n_1e6_32P = c(1e5, 1e5, 0, 0),
                  prop = 1)
  # Create a `m2` model with the adjusted initial conditions
  m2 <- m %>%
    set_init(init2, comp = "compartment", size = "n_1e6_32P", prop = "prop")
  expect_s3_class(m2, "networkModel")
  expect_equal(dim(m2), c(4, 5))
  flows <- tidy_flows(m2, fit, n = DRAWS)
  expect_s3_class(flows, "tidy_flows")
  expect_equal(nrow(flows), 4 * DRAWS)
  expect_equal(dim(groups(m2)), c(4, 2))
  flows_light <- flows %>% filter_by_group(light_treatment == "light") 
  expect_equal(nrow(flows_light), 2 * DRAWS)
  expect_error(quick_sankey(flows_light, node_s = "roundsquare", edge_f = 0.25),
               NA)

})
