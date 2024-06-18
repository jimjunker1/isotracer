### * Setup

library(magrittr)
library(tibble)
library(dplyr)
library(tidyr)

ITERS <- 52
DRAWS <- 24
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

### * Collins et al. 2016 - Nitrogen in Trinidadian streams

test_that("Collins et al. case study runs", {

  # Define model
  
  m <- new_networkModel() %>%
    set_topo("NH4, NO3 -> epi, FBOM", "epi -> petro, pseph", "FBOM -> tricor",
             "petro, tricor -> arg")
  m <- m %>% set_steady(c("NH4", "NO3"))
  expect_error(ggtopo(m, "sugiyama"), NA)
  inits <- lalaja %>% filter(time.days == 0)
  m <- m %>% set_init(inits, comp = "compartment", size = "mgN.per.m2", prop = "prop15N",
                      group_by = c("transect"))
  expect_s3_class(m, "networkModel")
  expect_equal(dim(m), c(3, 5))
  expect_equal(dim(groups(m)), c(3, 1))
  obs <- lalaja %>% filter(time.days > 0)
  m <- m %>% set_obs(obs, time = "time.days")
  expect_error(plot(m, facet_row = "group", facet_col = "compartment", type = "prop",
                    scale = "all", log = TRUE,
                    comps = c("NH4", "NO3", "epi", "FBOM", "petro", "pseph", "tricor", "arg")),
               NA)
  m <- m %>% set_split(c("epi", "FBOM"))
  expect_equal(dim(m), c(3, 5))
  pulses <- tribble(
    ~ stream,    ~ transect, ~ comp, ~ time, ~ qty_14N, ~ qty_15N,
    "UL",  "transect.1",  "NH4",     11,         0,  -0.00569,
    "UL",  "transect.2",  "NH4",     11,         0,  -0.00264,
    "UL",  "transect.3",  "NH4",     11,         0, -0.000726,
    "UL",  "transect.1",  "NO3",     11,         0,  -0.00851,
    "UL",  "transect.2",  "NO3",     11,         0,  -0.01118,
    "UL",  "transect.3",  "NO3",     11,         0,  -0.01244,
    )
  m <- m %>% add_pulse_event(pulses = pulses, comp = "comp", time = "time",
                             unmarked = "qty_14N", marked = "qty_15N")
  obs %>% select(compartment, mgN.per.m2) %>%
    na.omit() %>%
    group_by(compartment) %>%
    summarize(mean = mean(mgN.per.m2),
              sd = sd(mgN.per.m2),
              cv = sd / mean) %>%
    arrange(cv)
  m <- m %>% set_size_family("normal_sd", by_compartment = TRUE)
  m <- set_prior(m, constant_p(0.1), "zeta_NH4") %>%  # Dummy values for the steady
    set_prior(constant_p(0.1), "zeta_NO3") %>%        # state dissolved nutrients
    set_prior(constant_p(115), "zeta_epi") %>%
    set_prior(constant_p(2436), "zeta_FBOM") %>%
    set_prior(constant_p(10.7), "zeta_pseph") %>%
    set_prior(constant_p(2.88), "zeta_arg") %>%
    set_prior(constant_p(13.3), "zeta_tricor") %>%
    set_prior(constant_p(2.48), "zeta_petro")
  expect_equal(dim(missing_priors(m)), c(20, 2))
  m <- set_priors(m, normal_p(0, 0.5), "upsilon|lambda")
  m <- m %>%
    set_prior(normal_p(0, 1000), "upsilon_N") %>%
    set_prior(constant_p(0), "lambda_N")
  expect_equal(dim(missing_priors(m)), c(3, 2))
  expect_equal(prop_family(m), "gamma_cv")
  m <- set_priors(m, normal_p(0, 3), "^eta")
  m <- set_priors(m, uniform_p(0, 1), "portion")

  # Run model
  
  run <- purrr::quietly(run_mcmc)(m, iter = ITERS)$result
  expect_error(plot(run), NA)
  pred <- predict(m, run, draws = DRAWS)
  expect_s3_class(pred, "networkModel")
  expect_equal(dim(pred), c(3, 7))
  expect_error(plot(pred, facet_row = "group", facet_col = "compartment", type = "prop",
                    scale = "all", log = TRUE,
                    comps = c("epi", "FBOM", "petro", "pseph", "tricor", "arg")),
               NA)
  m2 <- m %>% set_topo("NH4, NO3 -> epi, FBOM", "epi -> petro, pseph",
                       "FBOM -> tricor", "petro, pseph, tricor -> arg")
  m2 <- m2 %>%
    set_steady(c("NH4", "NO3")) %>%
    set_split(c("epi", "FBOM"))
  expect_error(ggtopo(m2, "sugiyama"), NA)
  m2 <- set_priors(m2, priors(m))
  expect_equal(dim(missing_priors(m2)), c(1, 2))
  m2 <- set_priors(m2, normal_p(0, 0.5), "upsilon_pseph_to_arg")
  expect_s3_class(m2, "networkModel")
  expect_equal(dim(m2), c(3, 6))
  run2 <- purrr::quietly(run_mcmc)(m2, iter = ITERS)$result
  expect_error(plot(run2), NA)
  pred2 <- predict(m2, run2, draws = DRAWS)
  expect_s3_class(pred2, "networkModel")
  expect_equal(dim(pred2), c(3, 7))
  expect_error(plot(pred2, facet_row = "group", facet_col = "compartment", type = "prop",
                    scale = "all", log = TRUE,
                    comps = c("epi", "FBOM", "petro", "pseph", "tricor", "arg")), NA)
  expect_equal(dim(dic(run, run2)), c(2, 6))
  flows <- tidy_flows(m, run, n = DRAWS, steady_state = TRUE)
  expect_equal(dim(flows), c(3 * DRAWS, 5))
  f <- flows %>%
    select(flows) %>%
    unnest(flows) %>%
    na.omit() %>% # We drop the rows corresponding to `lambda` losses
    group_by(from, to) %>%
    summarize(flow = mean(average_flow),
              sd = sd(average_flow),
              cv = sd / flow)
  expect_equal(dim(f), c(9, 5))
  nodes <- inits %>%
    ungroup() %>%
    filter(transect == "transect.1") %>%
    select(comp = compartment, size = mgN.per.m2) %>%
    mutate(label = comp)
  expect_equal(dim(nodes), c(8, 3))
  expect_error(sankey(topo(m), nodes, f %>% select(from, to, flow) %>% rename(width = flow),
                      edge_f = 0.2, layout = "left2right"),
               NA)
  expect_error(mcmc_heatmap(run), NA)
  
})
