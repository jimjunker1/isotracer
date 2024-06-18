### * Description

# Test the Stan model using adaptive dt

### * Setup

devtools::load_all()
library(tidyverse)
library(here)
options(scipen = 6)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(bayesplot)
library(coda)

### * Functions

### * myModel

myModel <- function() {
    myTopo <- makeTopo(c("NH4, NO3 -> epi, FBOM, seston", "epi -> pseph",
                         "FBOM -> seston", "seston -> FBOM"),
                       split = c("epi", "seston", "FBOM"))
    bkg <- tibble(comp = c("NH4", "NO3", "epi", "seston", "pseph", "FBOM"),
                  props = delta2prop(0, "d15N"))
    myBkg <- bkg %>% makeInitProps(proportion = "props", compartment = "comp")
    addition <- datasetAddition("collins2016")
    addition <- addition %>%
        filter(stream == "UPL" & compartment %in% c("NH4", "NO3"))
    myAddition <- addition %>% makeInput(compartment = "compartment",
                                         profile = "profile", group_by = "transect")
    obs <- datasetObs("collins2016")
    obs <- obs %>% filter(stream == "UPL")
    obs <- obs %>% mutate(proportion = delta2prop(obs$d15N, "d15N"))
    mySizes <- obs %>%
        makeCompSizes(size = "mgN.per.m2", time = "time.days", compartment = "compartment")
    myObs <- obs %>%
        makeTracer(prop = "proportion", time = "time.days", compartment = "compartment",
                   group_by = "transect")
    model <- networkModel(myTopo, mySizes, myBkg, myAddition, myObs)
    return(model)
}

### * Main

m <- myModel()

### ** Classic

r <- runStan(m, iter = 2000, gridSize = 1024,
             control = list(adapt_delta = 0.8))

pdf("toto.pdf", width = 9.5, height = 10, title = "Traces classic")
plot(r)
plotEstimatedProp(r, model = m)
plotEstimatedSize(r, model = m)
dev.off()

quit()

### ** Adaptive dt

stanModel <- stan_model(file = file.path(here::here(), "sandbox", "stan-models",
                                         "networkModel-adaptive-dt.stan"))
saveRDS(stanModel, file = "tmp.stan.model.rds")

stanModel <- readRDS("tmp.stan.model.rds")

s <- networkModelToStanData_adaptiveDt(m, minGrid = 128, f = 2^(0:8),
                                       relThreshold = 0.7)
stan.fit2 <- sampling(stanModel, data = s, chains = 1, iter = 200,
                 pars = c("params"),
                 control = list(adapt_delta = 0.8))

plot(stan.fit2)
plotTraces(As.mcmc.list(stan.fit2))

#plot(stan.fit)
z <- As.mcmc.list(stan.fit)
coda::varnames(z) <- c(s[["params_all"]], "lp__")
class(z) <- c("networkStanFit", "mcmc.list")
pdf("toto-solver.pdf", width = 9.5, height = 10, title = "Traces solver")
plotTraces(z)
plotEstimatedProp(z, model = m)
plotEstimatedSize(z, model = m)
dev.off()
