### * Setup

library(magrittr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

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

### * Li et al. 2017 - Protein turnover in Arabidopsis

test_that("Li et al. case study runs", {

  # Define model
  
  m <- new_networkModel() %>%
    set_topo("medium -> RBCL, THI1, CPN60A, PSBA, PGK1")
  m <- m %>% set_steady("medium")
  m <- m %>% set_prop_family("beta_phi")
  expect_equal(dim(m), c(1, 4))
  expect_error(ggtopo(m), NA)
  prots <- c("ATCG00490.1" = "RBCL", "AT5G54770.1" = "THI1", "AT2G28000.1" = "CPN60A",
             "ATCG00020.1" = "PSBA", "AT3G12780.1" = "PGK1")
  data <- li2017 %>%
    filter(prot_id %in% names(prots)) %>%
    mutate(prot = prots[prot_id]) %>%
    rename(compartment = prot) %>%
    select(compartment, leaf_id, time_day, rel_abundance, labeled_fraction)
  expect_equal(dim(data), c(284, 5))
  inits <- data %>%
    filter(time_day == 0) %>%
    group_by(compartment, leaf_id, time_day) %>%
    summarize(rel_abundance = mean(rel_abundance), .groups = "drop")
  inits <- inits %>%
    mutate(labeled_fraction = 0.003663)
  # We use the same medium specifications for all proteins
  inits_medium <- tibble(compartment = "medium", rel_abundance = 1, labeled_fraction = 1) %>%
    crossing(inits %>% select(leaf_id) %>% unique())
  # Put it all together
  inits <- bind_rows(inits, inits_medium)
  expect_equal(dim(inits), c(18, 5))
  m <- m %>% set_init(inits, comp = "compartment", size = "rel_abundance",
                      prop = "labeled_fraction", group_by = "leaf_id")
  m <- m %>% set_obs(data, time = "time_day")
  m <- m %>% add_covariates(upsilon + lambda ~ leaf_id)
  expect_equal(dim(m), c(3, 5))
  expect_identical(groups(m),
                   structure(list(leaf_id = c("leaf_3", "leaf_5", "leaf_7")),
                             class = c("tbl_df", "tbl", "data.frame"),
                             row.names = c(NA, -3L)))
  
  # Add priors
  
  expect_equal(nrow(params(m)), 35)
  m <- set_priors(m, normal_p(0, 2), "lambda")
  m <- set_priors(m, normal_p(0, 8), "upsilon")
  m <- set_prior(m, constant_p(0), "lambda_medium")
  expect_equal(prop_family(m), "beta_phi")
  m <- set_prior(m, normal_p(mean = 50, sd = 50), "^eta")
  expect_equal(size_family(m), "normal_cv")
  m <- set_prior(m, normal_p(0, 3), "zeta")

  # Run MCMC
  
  run <- purrr::quietly(run_mcmc)(m, iter = ITERS)$result
  expect_error(plot(run), NA)
  pred <- predict(m, run, draws = DRAWS)
  expect_s3_class(pred, "networkModel")
  expect_equal(dim(pred), c(3, 7))
  expect_error({
    plot(pred, facet_row = "group", facet_col = "compartment", type = "prop",
         scale = "all", comps = c("CPN60A", "PGK1", "RBCL", "THI1", "PSBA"))
    plot(pred, facet_row = "group", facet_col = "compartment", type = "size",
         scale = "all", comps = c("CPN60A", "PGK1", "RBCL", "THI1", "PSBA"))
  }, NA)
  
  # Downstream analyses
  
  library(ggplot2)
  library(ggdist)
  lambdas <- run %>% select("lambda") %>%
    tidy_mcmc() %>%
    mutate(mcmc.parameters = map(mcmc.parameters, enframe)) %>%
    pull(mcmc.parameters) %>%
    bind_rows() %>%
    separate(name, sep = "[|]", into = c("param", "leaf")) %>%
    mutate(prot = substr(param, 8, nchar(param)))
  expect_equal(dim(lambdas), c(15 * ceiling(ITERS / 2) * n_chains, 4))
  expect_error({ lambdas %>% group_by(prot, leaf) %>%
                   median_qi(value, .width = c(0.8, 0.95)) %>%
                   ggplot(aes(x = leaf, y = value, ymin = .lower, ymax = .upper)) +
                   facet_grid(. ~ prot) +
                   geom_pointinterval() +
                   coord_trans(y = "log10") + 
                   labs(y = "lambda (day-1)")
  }, NA)
  totimes <- lambdas %>% mutate(value = 1/value)
  expect_error({ totimes %>% group_by(prot, leaf) %>%
                   median_qi(value, .width = c(0.8, 0.95)) %>%
                   ggplot(aes(x = leaf, y = value, ymin = .lower, ymax = .upper)) +
                   facet_grid(. ~ prot) +
                   geom_pointinterval() +
                   coord_trans(y = "log10") +
                   scale_y_continuous(breaks = c(1, 3, 5, 10, 20, 40, 60)) + 
                   labs(y = "turnover time (day)")
  }, NA)
  flows <- tidy_flows(m, run, n = DRAWS)
  expect_equal(dim(flows), c(3 * DRAWS, 5))
  library(grid)
  expect_error({
    grid.newpage()
    vp <- viewport(layout = grid.layout(ncol = 2))
    pushViewport(vp)
    pushViewport(viewport(layout.pos.col = 1))
    flows_3 <- flows %>% filter_by_group(leaf_id == "leaf_3")
    quick_sankey(flows_3, node_s = "roundsquare", edge_f = 0.25, new = FALSE)
    popViewport()
    pushViewport(viewport(layout.pos.col = 2))
    flows_7 <- flows %>% filter_by_group(leaf_id == "leaf_7")
    quick_sankey(flows_7, node_s = "roundsquare", edge_f = 0.25, new = FALSE)
    popViewport()
  }, NA)

})
