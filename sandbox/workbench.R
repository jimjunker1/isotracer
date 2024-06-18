### * Description

# A script to try and test code.

### * Setup

library(tidyverse)
devtools::load_all()

### * Prepare data

# Load data
col_types <- c(
    stream = col_character(),
    transect = col_character(),
    compartment = col_character(),
    time.days = col_double(),
    mgN.per.m2 = col_double(),
    d15N = col_double(),
    prop15N = col_double()
)
data <- read_tsv("collins-dump.tsv", col_types = col_types)  %>%
  mutate(stream = ifelse(stream == "LOL", "LL", "UL")) %>%
  select(-d15N)

# Get initial conditions for NH4 and NO3 (during drip)
inits1 <- data %>% filter(compartment %in% c("NH4", "NO3") & time.days <= 10) %>%
  group_by(stream, transect, compartment) %>%
  summarize(mgN.per.m2 = mean(mgN.per.m2, na.rm = TRUE),
            prop15N = mean(prop15N, na.rm = TRUE))
inits2 <- data %>% filter(!compartment %in% c("NH4", "NO3")) %>%
  group_by(stream, compartment) %>%
  summarize(mgN.per.m2 = mean(mgN.per.m2, na.rm = TRUE)) %>%
  mutate(prop15N = delta2prop(0, "d15N")) %>%
  crossing(tibble(transect = c("transect.1", "transect.2", "transect.3")))
inits <- bind_rows(inits1, inits2)
inits[["mgN.per.m2"]][inits[["stream"]] == "LL" & inits[["compartment"]] == "FBOM"] <- 5267 # Value from UL

# Get pulse info
postdrip <- data %>% filter(compartment %in% c("NH4", "NO3") & time.days > 10) %>%
  group_by(stream, transect, compartment) %>%
  summarize(mgN.per.m2 = mean(mgN.per.m2, na.rm = TRUE),
            prop15N = mean(prop15N, na.rm = TRUE))
pulses <- (postdrip$prop15N - inits1$prop15N) * inits1[["mgN.per.m2"]]

### ** Prepare model

## links <- c("NH4, NO3 -> epi, seston, FBOM, CBOM",
##            "epi -> petro, pseph", "seston -> lepto",
##            "FBOM -> tricor", "CBOM -> eudan, phylo",
##            "tricor -> arg, euthy", "petro -> arg")

links <- c("NH4, NO3 -> epi, CBOM", "epi -> pseph",
           "CBOM -> phylo")

m <- new_networkModel(quiet = TRUE) %>%
    set_topo(links) %>%
    set_steady(c("NH4", "NO3")) %>%
    set_split(c("epi", "CBOM"))

m <- m %>%
    set_init(inits, comp = "compartment", size = "mgN.per.m2", prop = "prop15N",
             group_by = c("stream", "transect"))

m <- m %>%
    set_obs(data %>% filter(time.days > 0), comp = "compartment",
            size = "mgN.per.m2", prop = "prop15N", time = "time.days",
            group_by = c("stream", "transect"))

m %<>% add_pulse_event(time = 11, comp = "NH4", unmarked = 0, marked = -0.00538, which = 1) %>%
  add_pulse_event(time = 11, comp = "NO3", unmarked = 0, marked = -0.0168, which = 1) %>%
  add_pulse_event(time = 11, comp = "NH4", unmarked = 0, marked = -0.00219, which = 2) %>%
  add_pulse_event(time = 11, comp = "NO3", unmarked = 0, marked = -0.021, which = 2) %>%
  add_pulse_event(time = 11, comp = "NH4", unmarked = 0, marked = -0.000752, which = 3) %>%
  add_pulse_event(time = 11, comp = "NO3", unmarked = 0, marked = -0.0216, which = 3) %>%
  add_pulse_event(time = 11, comp = "NH4", unmarked = 0, marked = -0.00569, which = 4) %>%
  add_pulse_event(time = 11, comp = "NO3", unmarked = 0, marked = -0.00852, which = 4) %>%
  add_pulse_event(time = 11, comp = "NH4", unmarked = 0, marked = -0.00264, which = 5) %>%
  add_pulse_event(time = 11, comp = "NO3", unmarked = 0, marked = -0.0112, which = 5) %>%
  add_pulse_event(time = 11, comp = "NH4", unmarked = 0, marked = -0.000752, which = 6) %>%
  add_pulse_event(time = 11, comp = "NO3", unmarked = 0, marked = -0.0125, which = 6)  

# Covariates
m <- m %>%
    add_covariates(. ~ stream,
                   eta + zeta ~ 1)

# Priors
m <- m %>%
    set_prior(constant_p(0), "lambda_N", quiet = TRUE) %>%
    set_prior(hcauchy_p(150), "upsilon_N", quiet = TRUE) #%>%
    #set_prior(uniform_p(0, 1), "portion", quiet = TRUE)

# Run
if (!file.exists("toto-r.rds")) {
    r <- run_mcmc(m, iter = 1000, cores = 4)
    saveRDS(r, "toto-r.rds")
} else {
    r <- readRDS("toto-r.rds")
}

plot(r)

### * Calculation of flows

tt <- tidy_trajectories(m, r, n = 20)
tf <- tidy_flows(m, r, n_per_chain = 5)
tfmc <- as.mcmc.list(tf)
pdf("toto.pdf", width = 10, height = 8)
plot(tfmc)
dev.off()

filter_by_group(tf, stream == "LL", transect == "transect.1")
filter_by_group(tt, stream == "LL", transect == "transect.1")

stf <- tidy_flows(m, r, n_per_chain = 3, steady_state = TRUE)

### * Sankey plotting

z <- filter_by_group(stf, stream == "UL", transect == "transect.1")
flows <- bind_rows(z$flows) %>%
    group_by(from, to) %>%
    summarize(width = mean(average_flow)) %>%
    na.omit()
flows$fill <- "grey"
flows$fill[flows$from == "NH4"] <- "darkgrey"
flows$fill[flows$from == "NO3"] <- "lightgrey"
nodes <- filter_by_group(m, stream == "LL", transect == "transect.1")$initial[[1]]
topo <- topo(m)[[1]]
nodes$comp <- nodes$compartment

sankey(topo, nodes, flows, edge_f = 0.3)

### * Playing with steady states

z <- m[1, ]
z <- set_params(z, summary(r)$statistics[, "Mean"])
proj <- project(z, end = 100)
plot(proj)

tail(proj$trajectory[[1]]$sizes[[1]])
calculate_steady_state_one_row(z)

m <- set_params(m, summary(r)$statistics[, "Mean"])
z <- calculate_steady_state(m)

z <- tidy_steady_states(m, r, n_per_chain = 5)
tsmcmc <- as.mcmc.list(z)
plot(tsmcmc %>% select("transect.1"))
summary(tsmcmc %>% select("transect.1"))

### * Debugging

x <- trini_mod
group_mapping <- x[["group"]]
names(group_mapping) <- sapply(group_mapping, group2string)
x[["group"]] <- names(group_mapping)
x <- x[, c("observations", "group")]

z <- tibble(a = 1:5,
            b = 6:10)
x$observations <- rep(list(z), 6)

# Not working
tbl_vars(x)

# Not working as expected
x$group <- letters[1:6]
tbl_vars(x)

class(x) <- c("tbl_df", "tbl", "data.frame")
tbl_vars(x)
# The problem was from the groups.networkModel method!

### * Radioactive decay

set.seed(4)

m <- new_networkModel() %>%
    set_topo("A -> B") %>%
    set_half_life(8)

inits <- tibble(comp = c("A", "B"),
                t = 0,
                size = 1,
                prop = c(0.5, 0))
m <- m %>% set_init(inits, "comp", "size", "prop")
params(m)
m <- set_params(m, c("eta" = 0.1, "zeta" = 0.1, "lambda_A" = 0,
                     "lambda_B" = 0, "upsilon_A_to_B" = 0.1))
m <- project(m, end = 20)
plot(m)

obs <- sample_from(m, at = 1:10 * 2)

m <- set_obs(m, obs, time = "time")

plot(m)

fit <- run_mcmc(m)
plot(fit)

m_no_hl <- set_half_life(m, 0)
fit_no_hl <- run_mcmc(m_no_hl)
plot(fit_no_hl)

### * Only tracking radioactive material (a la McRoy 1970)

set.seed(4)

m <- new_networkModel() %>%
    set_topo("A -> B -> C -> A") %>%
    set_half_life(8)

inits <- tibble(comp = c("A", "B", "C"),
                t = 0,
                size = c(1, 0, 0),
                prop = c(1, 0, 0))
m <- m %>% set_init(inits, "comp", "size", "prop")
params(m)

m <- set_prior(m, constant_p(1), "^eta") %>%
    set_prior(constant_p(0), "lambda")

m <- set_params(m, c("eta" = 0.1, "zeta" = 0.1, "lambda_A" = 0,
                     "lambda_B" = 0, "lambda_C" = 0,
                     "upsilon_A_to_B" = 0.1,
                     "upsilon_B_to_C" = 0.1,
                     "upsilon_C_to_A" = 0.05))
m <- project(m, end = 20)
plot(m)

obs <- sample_from(m, at = 1:10 * 2) %>%
    mutate(prop = NA)
m <- set_obs(m, obs, time = "time")

plot(m)

fit <- run_mcmc(m, iter = 500)
plot(fit)

pred <- predict(m, fit, draws = 100)
plot(pred)

flows <- tidy_flows(m, fit, n = 50)

flows_ss <- tidy_flows(m, fit, n = 50, steady_state = TRUE)

f <- flows %>% select(flows) %>%
    unnest(flows) %>%
    group_by(from, to) %>%
    summarize(avg = mean(average_flow),
              sd = sd(average_flow),
              cv = sd/avg) %>%
    mutate(width = avg)

f_ss <- flows_ss %>% select(flows) %>%
    unnest(flows) %>%
    group_by(from, to) %>%
    summarize(avg = mean(average_flow),
              sd = sd(average_flow),
              cv = sd/avg) %>%
    mutate(width = avg)

s <- nodes_from_topo(topo(m)) %>%
    mutate(size = 1)

sankey(topo(m), s, na.omit(f), debug = T)

sankey(topo(m), s, na.omit(f_ss), debug = T)
