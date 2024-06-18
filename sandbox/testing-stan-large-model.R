### * Description

# Testing the Stan implementation with the largest model from Collins 2016

### * Packages

library(isotracer)
library(tidyverse)

### * Run

targetStreams <- c("UPL", "LOL")
# Prepare data
myTopo <- makeTopo(c("NH4, NO3 -> epi, seston, FBOM, CBOM",
                     "epi -> petro, pseph", "seston -> lepto",
                     "FBOM -> tricor", "CBOM -> eudan, phylo",
                     "tricor -> arg, euthy"),
                   split = c("epi", "CBOM", "FBOM", "seston"))
bkg <- tibble(comp = c("NH4", "NO3", "seston", "epi", "CBOM", "FBOM", "eudan",
                       "phylo", "petro", "pseph", "arg", "lepto", "tricor", "euthy"),
              props = delta2prop(0, "d15N"))
myBkg <- bkg %>% makeInitProps(proportion = "props", compartment = "comp")
addition <- datasetAddition("collins2016")
myAddition <- addition %>%
    makeInput(compartment = "compartment", profile = "profile",
              group_by = c("stream", "transect"))
obs <- datasetObs("collins2016")
obs <- obs %>% mutate(proportion = delta2prop(obs$d15N, "d15N")) %>%
    filter(stream %in% targetStreams)
# Use FBOM biomass from UPL for LOL
lolFBOM <- obs %>% filter(compartment == "FBOM" & !is.na(mgN.per.m2)) %>%
    mutate(stream = "LOL")
obs <- bind_rows(obs, lolFBOM)
mySizes <- obs %>%
    makeCompSizes(size = "mgN.per.m2", time = "time.days", compartment = "compartment",
                  group_by = "stream")
myObs <- obs %>%
    makeTracer(prop = "proportion", time = "time.days", compartment = "compartment",
               group_by = c("stream", "transect"))

model <- networkModel(myTopo, mySizes, myBkg, myAddition, myObs)
model$parameterNames
model <- addFormula(model, . ~ stream, eta ~ 1)

# Run
r <- runStan(model, dt = 0.1, iter = 2000, refresh = 10,
             control = list(adapt_delta = 0.9))
plot(r)
saveRDS(r, "tmp.r.rds")



library(rstan)
readStanCsv <- function(basename) {
    stanFiles <- list.files(dirname(basename), pattern = basename)
    chains <- lapply(stanFiles, function(filename) {
        x <- rstan::read_stan_csv(file.path(dirname(basename), filename))
        s <- as.data.frame(x@sim$samples)
        s <- s[, grep("params", names(s))]
        z <- coda::as.mcmc(s)
        return(z)
    })
    return(chains)
}
# If call with
# r <- runStan(model, dt = 0.1, iter = 1000, sample_file = "toto")
# Then the chains can be loaded with:
# r <- readStanCsv("toto")

# Run with expm model
s <- networkModelToStanDataExpm(model)
stan.fit <- stan(file = file.path(here::here(), "sandbox", "stan-models", "networkModel-expm.stan"),
                 data = s, iter = 100, chains = 1,
                 pars = c("params"), control = list(adapt_delta = 0.95))
