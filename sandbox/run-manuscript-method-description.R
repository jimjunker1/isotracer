### * Description

# This script runs the models used in the manuscript by Andrés López-Sepulcre
# entitled "A new method to reconstruct quantitative food webs and nutrient
# flows from isotope tracer addition experiments" and submitted to The American
# Naturalist.

### * Packages

library(isotracer)
library(coda)
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
options(scipen = 6)

### * Parameters

gridSize <- 256
iter <- 1000
chains <- 4

### * Functions

#' @param stream "UPL" or "LOL
#' @param psi A three-character string of 0 and 1 describing the foodweb
#'   topology used (see the manuscript for details)
#'
#' @return A networkModelR6 object ready for MCMC run

prepareModel <- function(stream, psi) {
    topoLinks <- c("NH4, NO3 -> epi, seston, FBOM, CBOM",
                   "epi -> petro, pseph", "seston -> lepto",
                   "FBOM -> tricor", "CBOM -> eudan, phylo",
                   "tricor -> arg, euthy")
    myStream <- stream
    # Parse psi code
    psi <- strsplit(psi, "")[[1]]
    if (psi[1] == "1") {
        topoLinks <- c(topoLinks, "petro -> arg")
    }
    if (psi[2] == "1") {
        topoLinks <- c(topoLinks, "pseph -> arg")
    }
    if (psi[3] == "1") {
        topoLinks <- c(topoLinks, "FBOM -> eudan")
    }
    # Make topo
    myTopo <- makeTopo(topoLinks,
                       split = c("epi", "CBOM", "FBOM", "seston"))
    # Make background
    bkg <- tibble(comp = c("NH4", "NO3", "seston", "epi", "CBOM", "FBOM", "eudan",
                           "phylo", "petro", "pseph", "arg", "lepto", "tricor", "euthy"),
                  props = delta2prop(0, "d15N"))
    myBkg <- bkg %>% makeInitProps(proportion = "props", compartment = "comp")
    # Make addition regime
    addition <- datasetAddition("collins2016")
    myAddition <- addition %>%
        filter(stream %in% myStream) %>%
        makeInput(compartment = "compartment", profile = "profile",
                  group_by = c("transect"))
    # Get observations
    obs <- datasetObs("collins2016")
    obs <- obs %>% mutate(proportion = delta2prop(obs$d15N, "d15N"))
    # Use FBOM biomass from UPL for LOL
    lolFBOM <- obs %>% filter(compartment == "FBOM" & !is.na(mgN.per.m2)) %>%
        mutate(stream = "LOL")
    obs <- bind_rows(obs, lolFBOM) %>%
        filter(stream %in% myStream)
    # Make sizes
    mySizes <- obs %>%
        makeCompSizes(size = "mgN.per.m2", time = "time.days", compartment = "compartment")
    # Make proportions
    myObs <- obs %>%
        makeTracer(prop = "proportion", time = "time.days", compartment = "compartment",
                   group_by = c("transect"))
    # Build model
    model <- networkModel(myTopo, mySizes, myBkg, myAddition, myObs)
    # Set prior boundaries
    ## Replace hcauchy_p(0.5) by scaled_beta_p(1, 20, 10) or scaled_beta(0.4, 12, 10)
    model <- set_prior(model, scaled_beta_p(1, 3, 1), param = "uptakeRate_from_", verbose = FALSE)
    model <- set_prior(model, hcauchy_p(250), param = "uptakeRate_from_NH4_", verbose = FALSE)
    model <- set_prior(model, hcauchy_p(250), param = "uptakeRate_from_NO3_", verbose = FALSE)
    model <- set_prior(model, scaled_beta_p(1, 3, 1), param = "lossRate_", verbose = FALSE)
    model <- set_prior(model, uniform_p(0, 1), param = "portion_", verbose = FALSE)
    model <- set_prior(model, hcauchy_p(1), param = "eta", verbose = FALSE)
    ## # Priors for seston
    ## model <- set_prior(model, scaled_beta_p(1, 20, 10), param = "uptakeRate_from_seston.act_to_lepto", verbose = FALSE)
    ## model <- set_prior(model, scaled_beta_p(1, 20, 10), param = "lossRate_seston.act", verbose = FALSE)
    ## # Priors for tricor
    ## model <- set_prior(model, scaled_beta_p(1, 3, 1), param = "lossRate_tricor", verbose = FALSE)
    # Priors for eudan and lepto
    model <- set_prior(model, scaled_beta_p(1, 3, 0.5), param = "uptakeRate_from_CBOM.act_to_eudan", verbose = FALSE)
    model <- set_prior(model, scaled_beta_p(1, 3, 0.5), param = "uptakeRate_from_seston.act_to_lepto", verbose = FALSE)
    # Return
    return(model)
}

prepareTestModel <- function(stream) {
    topoLinks <- c("NH4, NO3 -> seston -> lepto")
    myStream <- stream
    # Make topo
    myTopo <- makeTopo(topoLinks,
                       split = c("seston"))
    # Make background
    bkg <- tibble(comp = c("NH4", "NO3", "seston", "lepto"),
                  props = delta2prop(0, "d15N"))
    myBkg <- bkg %>% makeInitProps(proportion = "props", compartment = "comp")
    # Make addition regime
    addition <- datasetAddition("collins2016")
    myAddition <- addition %>%
        filter(stream %in% myStream) %>%
        makeInput(compartment = "compartment", profile = "profile",
                  group_by = c("transect"))
    # Get observations
    obs <- datasetObs("collins2016")
    obs <- obs %>% mutate(proportion = delta2prop(obs$d15N, "d15N"))
    # Use FBOM biomass from UPL for LOL
    lolFBOM <- obs %>% filter(compartment == "FBOM" & !is.na(mgN.per.m2)) %>%
        mutate(stream = "LOL")
    obs <- bind_rows(obs, lolFBOM) %>%
        filter(stream %in% myStream) %>%
        filter(compartment %in% c("NH4", "NO3", "seston", "lepto"))
    # Make sizes
    mySizes <- obs %>%
        makeCompSizes(size = "mgN.per.m2", time = "time.days", compartment = "compartment")
    # Make proportions
    myObs <- obs %>%
        makeTracer(prop = "proportion", time = "time.days", compartment = "compartment",
                   group_by = c("transect"))
    # Build model
    model <- networkModel(myTopo, mySizes, myBkg, myAddition, myObs)
    # Set prior boundaries
    model <- set_prior(model, hcauchy_p(0.5), param = "uptakeRate_from_", verbose = FALSE)
    model <- set_prior(model, hcauchy_p(5), param = "uptakeRate_from_NH4_", verbose = FALSE)
    model <- set_prior(model, hcauchy_p(5), param = "uptakeRate_from_NO3_", verbose = FALSE)
    model <- set_prior(model, hcauchy_p(0.5), param = "lossRate_", verbose = FALSE)
    model <- set_prior(model, uniform_p(0, 1), param = "portion_", verbose = FALSE)
    model <- set_prior(model, hcauchy_p(1), param = "eta", verbose = FALSE)
    # Return
    return(model)
}

### * Run

psi <- c("000", "001", "010", "100", "011", "101", "110", "111")
stream <- c("UPL", "LOL")
models <- crossing(stream, psi) %>%
    mutate(model = map2(stream, psi, prepareModel),
           date = Sys.time())

# Run one model at a time
if (FALSE) {
    models$run <- as.list(rep(NA, nrow(models)))
    for (i in seq_len(nrow(models))) {
        models$run[[i]] <- runStan(models$model[[i]], gridSize = gridSize,
                                   iter = iter, chains = chains,
                                   control = list(adapt_delta = 0.9))
        saveRDS(models$run[[i]], file = paste0("tmp.run.", models$stream[i], ".",
                                               models$psi[i], ".rds"))
    }
    saveRDS(models, file = "tmp.ms.models.rds")
}

# Rerun models with bad mixing
if (FALSE) {
    models <- readRDS("tmp.ms.models.rds")
    reruns <- c(1,2,3,5,7,14)
    for (i in reruns) {
        models$run[[i]] <- runStan(models$model[[i]], gridSize = gridSize,
                                   iter = iter, chains = chains,
                                   control = list(adapt_delta = 0.9))
        saveRDS(models$run[[i]], file = paste0("tmp.run.", models$stream[i], ".",
                                               models$psi[i], ".rds"))
    }
    saveRDS(models, file = "tmp.ms.models.rds")
}

models <- readRDS("tmp.ms.models.rds")
models$rhat <- sapply(models$run, function(x) gelman.diag(x)$mpsrf)
models$DIC <- sapply(models$run, dic)

pdf("toto.pdf", width = 19, height = 10)
invisible(sapply(models$run, function(x) plot(x, drawHist = TRUE, nBars = 32)))
dev.off()

### * Table 2: DIC table

### ** Function to save DIC table in tex format
tex_dicTable <- function(x, file = "tabular-model-comparison-DIC.tex") {
    nParams <- c("100" = 68, "000" = 66, "110" = 70, "010" = 68,
                 "011" = 70, "101" = 70, "001" = 68, "111" = 72)
    fo <- file(file, open = "w")
    cat("% Generated by script run-manuscript-method-description.R in folder projects_running/2018-05_isotracer/sandbox/

\\begin{tabular}{cp{2cm}p{2cm}p{2cm}cccc}
    \\toprule
    Model & \\multicolumn{3}{c}{Trophic link} & \\# param. & DIC & $\\Delta$DIC & $w_{DIC}$ \\\\
    \\midrule
    & \\textit{Petrophila} $\\rightarrow$ \\textit{Argia} & \\textit{Psephenus} $\\rightarrow$ \\textit{Argia} & FBOM $\\rightarrow$ \\textit{Eudaniella} & & & & \\\\
    \\midrule
    % Using the multicolumn hack to center Yes and No
", file = fo)
    for (i in seq_len(nrow(x))) {
        l <- paste0("$\\mathbf{\\Psi}_{", x$psi[i], "}$ & \\multicolumn{1}{c}{", c("0" = "No", "1" = "Yes")[x$psiCode[[i]][1]], "} & \\multicolumn{1}{c}{", c("0" = "No", "1" = "Yes")[x$psiCode[[i]][2]], "} & \\multicolumn{1}{c}{", c("0" = "No", "1" = "Yes")[x$psiCode[[i]][3]],"} & ", nParams[x$psi[i]], " & ", round(x$dic[i], 1), "  & ", round(-x$deltaDic[i], 1), " & ", round(x$weight[i], 3), " \\\\\n")
        cat(l, file = fo)
    }
    cat("    \\bottomrule
  \\end{tabular}
", file = fo)
    close(fo)
}

### ** Save DIC table
dicTable <- models %>% select(stream, psi, DIC) %>%
    spread(key = stream, value = DIC) %>%
    mutate(dic = LOL + UPL) %>%
    arrange(dic) %>%
    mutate(deltaDic = min(dic)-dic) %>%
    mutate(weight = exp(deltaDic/2)) %>%
    mutate(weight = weight / sum(weight)) %>%
    mutate(psiCode = strsplit(psi, ""))
tex_dicTable(dicTable)
models
dicTable

### ** Get best model
bestPsi <- dicTable$psi[1]
lolRow <- which(models$psi == bestPsi & models$stream == "LOL")
uplRow <- which(models$psi == bestPsi & models$stream == "UPL")

### * Figure S2: Primary traces

RES = 150
width = 8
height = 11
ratio = 16/9

png("supfig-02-traces-LL.png",
    width = width * RES, height = height * RES, res = RES,
    type = "cairo-png", bg = "transparent", family = "Century Schoolbook L")
plot(models$run[[lolRow]], ratio = ratio, drawHist = TRUE, nBars = 32)
dev.off()
cairo_pdf("supfig-02-traces-LL.pdf",
          width = width, height = height, family = "Century Schoolbook L")
plot(models$run[[lolRow]], ratio = ratio, drawHist = TRUE, nBars = 32)
dev.off()
png("supfig-02-traces-UL.png",
    width = width * RES, height = height * RES, res = RES,
    type = "cairo-png", bg = "transparent", family = "Century Schoolbook L")
plot(models$run[[uplRow]], ratio = ratio, drawHist = TRUE, nBars = 32)
dev.off()
cairo_pdf("supfig-02-traces-UL.pdf",
          width = width, height = height, family = "Century Schoolbook L")
plot(models$run[[uplRow]], ratio = ratio, drawHist = TRUE, nBars = 32)
dev.off()

### * Figure S1: Model fit (d15N)

comps <- c("epi", "CBOM", "FBOM", "seston", "petro", "pseph",
           "tricor", "lepto", "eudan", "phylo", "arg", "euthy")
compLabels <- c("epilithon", "CBOM", "FBOM", "seston", "Petrophila", "Psephenus",
                "Tricorythodes", "Leptonema", "Eudaniela", "Phylloicus", "Argia", "Euthyplocia")

pdf("supfig-01-posterior-predictive-fit-LL.pdf",
    width = 13, height = 10)
plotFit(models$run[[lolRow]], models$model[[lolRow]], n = 100,
        compartments = comps, compartmentLabels = compLabels,
        deltaHeavy = TRUE, predProb = 0.95, nGrid = 128)
dev.off()

pdf("supfig-01-posterior-predictive-fit-UL.pdf",
    width = 13, height = 10)
plotFit(models$run[[uplRow]], models$model[[uplRow]], n = 100,
        compartments = comps, compartmentLabels = compLabels,
        deltaHeavy = TRUE, predProb = 0.95, nGrid = 128)
dev.off()

### * Figure Sxxx: Model fit (biomass)

pdf("supfig-xx-biomass-trajectories-LL.pdf",
    width = 13, height = 10)
plotFit(models$run[[lolRow]], models$model[[lolRow]], n = 200,
        compartments = comps, compartmentLabels = compLabels,
        sizes = TRUE, yscale = "compZero")
dev.off()

pdf("supfig-xx-biomass-trajectories-UL.pdf",
    width = 13, height = 10)
plotFit(models$run[[uplRow]], models$model[[uplRow]], n = 200,
        compartments = comps, compartmentLabels = compLabels,
        sizes = TRUE, yscale = "compZero")
dev.off()

### * Table S1: Primary parameters estimates
tex_primary_estimates <- function(runLL, runUL, file = "tabular-S1.tex") {
    ### * Buid tibble
    params <- tibble(paramName = c("portion_epi.act", "portion_CBOM.act",
                                   "portion_FBOM.act", "portion_seston.act",
                                   "uptakeRate_from_NH4_to_epi.act",
                                   "uptakeRate_from_NH4_to_CBOM.act",
                                   "uptakeRate_from_NH4_to_FBOM.act",
                                   "uptakeRate_from_NH4_to_seston.act",
                                   "uptakeRate_from_NO3_to_epi.act",
                                   "uptakeRate_from_NO3_to_CBOM.act",
                                   "uptakeRate_from_NO3_to_FBOM.act",
                                   "uptakeRate_from_NO3_to_seston.act",
                                   "uptakeRate_from_epi.act_to_petro",
                                   "uptakeRate_from_epi.act_to_pseph",
                                   "uptakeRate_from_CBOM.act_to_eudan",
                                   "uptakeRate_from_CBOM.act_to_phylo",
                                   "uptakeRate_from_FBOM.act_to_tricor",
                                   "uptakeRate_from_seston.act_to_lepto",
                                   "uptakeRate_from_petro_to_arg",
                                   "uptakeRate_from_tricor_to_arg",
                                   "uptakeRate_from_tricor_to_euthy",
                                   "lossRate_epi.act", "lossRate_CBOM.act",
                                   "lossRate_FBOM.act", "lossRate_seston.act",
                                   "lossRate_petro", "lossRate_pseph",
                                   "lossRate_tricor", "lossRate_lepto",
                                   "lossRate_eudan", "lossRate_phylo",
                                   "lossRate_arg", "lossRate_euthy", "eta"),
                     paramOrder = 1:34,
                     paramLatex = c("$\\pi_\\textrm{epi}$", "$\\pi_\\textrm{CBOM}$",
                                    "$\\pi_\\textrm{FBOM}$", "$\\pi_\\textrm{ses}$",
                                    "$\\upsilon_{\\textrm{epi}, \\textrm{NH4}}$",
                                    "$\\upsilon_{\\textrm{CBOM}, \\textrm{NH4}}$",
                                    "$\\upsilon_{\\textrm{FBOM}, \\textrm{NH4}}$",
                                    " $\\upsilon_{\\textrm{ses}, \\textrm{NH4}}$",
                                    "$\\upsilon_{\\textrm{epi}, \\textrm{NO3}}$",
                                    "$\\upsilon_{\\textrm{CBOM}, \\textrm{NO3}}$",
                                    "$\\upsilon_{\\textrm{FBOM}, \\textrm{NO3}}$",
                                    "$\\upsilon_{\\textrm{ses}, \\textrm{NO3}}$",
                                    "$\\upsilon_{\\textrm{pet}, \\textrm{epi}}$",
                                    "$\\upsilon_{\\textrm{pse}, \\textrm{epi}}$",
                                    "$\\upsilon_{\\textrm{eud}, \\textrm{CBOM}}$",
                                    "$\\upsilon_{\\textrm{phy}, \\textrm{CBOM}}$",
                                    "$\\upsilon_{\\textrm{tri}, \\textrm{FBOM}}$",
                                    "$\\upsilon_{\\textrm{lep}, \\textrm{ses}}$",
                                    "$\\upsilon_{\\textrm{arg}, \\textrm{pet}}$",
                                    "$\\upsilon_{\\textrm{arg}, \\textrm{tri}}$",
                                    "$\\upsilon_{\\textrm{eut}, \\textrm{tri}}$",
                                    "$\\lambda_\\textrm{epi}$", "$\\lambda_\\textrm{CBOM}$",
                                    "$\\lambda_\\textrm{FBOM}$", "$\\lambda_\\textrm{ses}$",
                                    "$\\lambda_\\textrm{pet}$", "$\\lambda_\\textrm{pse}$",
                                    "$\\lambda_\\textrm{tri}$", "$\\lambda_\\textrm{lep}$",
                                    "$\\lambda_\\textrm{eud}$", "$\\lambda_\\textrm{phy}$",
                                    "$\\lambda_\\textrm{arg}$", "$\\lambda_\\textrm{eut}$", "$\\eta$"))
    params$LL <- unlist(lapply(1:nrow(params), function(i) {
        z <- signif(quantile(unlist(runLL[, params$paramName[i]]), probs = c(0.025, 0.5, 0.975)), 3)
        cell <- paste0("$", z[2], " \\textrm{~} [", z[1], "-", z[3], "]$")
        return(cell)
    }))
    params$UL <- unlist(lapply(1:nrow(params), function(i) {
        z <- signif(quantile(unlist(runUL[, params$paramName[i]]), probs = c(0.025, 0.5, 0.975)), 3)
        cell <- paste0("$", z[2], " \\textrm{~} [", z[1], "-", z[3], "]$")
        return(cell)
    }))
    ### * Build tex table
    fo <- file(file, open = "w")
    cat("% Generated by script run-manuscript-method-description.R in folder projects_running/2018-05_isotracer/sandbox/

\\resizebox{\\hsize}{!}{
\\begin{tabular}{lcc}
\\toprule
Parameter & LL estimate & UL estimate \\\\
\\midrule
", file = fo)
    for (i in 1:nrow(params)) {
        l <- paste0(params$paramLatex[i], " & ", params$LL[i], " & ", params$UL[i], " \\\\")
        cat(l, "\n", file = fo)
    }
    cat("\\bottomrule
\\end{tabular}
}
", file = fo)
    close(fo)
}

### ** Save the table
tex_primary_estimates(models$run[[lolRow]], models$run[[uplRow]], "tabular-S1.tex")

### * Table S2: Steady state compartment biomasses

sss <- list("ll" = tidySSSizes(models$run[[lolRow]], models$model[[lolRow]]),
            "ul" = tidySSSizes(models$run[[uplRow]], models$model[[uplRow]]))
sss <- lapply(sss, function(x) x %>% filter(transect == "transect.1"))
sss <- lapply(sss, function(x) x %>% pull(network.sssizes))
sss <- lapply(sss, function(x) {
    bind_rows(lapply(x, function(y) do.call(tibble, as.list(y))))
})
sss <- lapply(sss, function(d) {
    d$epi <- d$epi.act + d$epi.refr
    d$seston <- d$seston.act + d$seston.refr
    d$CBOM <- d$CBOM.act + d$CBOM.refr
    d$FBOM <- d$FBOM.act + d$FBOM.refr
    d
})
sss <- lapply(sss, function(d) {
    lapply(d, function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
})

### ** Function to buid the steady state biomasses table

tex_steady_biomasses <- function(sssList, file = "tabular-steady-state-biomasses.tex") {
    fo <- file(file, open = "w")
    wr <- function(x) {
        x <- signif(x, 3)
        paste0(" ", x[2], " [", x[1], "-", x[3], "] ")
    }
    sss <- sssList
    # Top
    cat("  \\begin{tabular}{lcc}
    \\toprule
    \\multirow{2}{*}{Compartment $i$} & \\multicolumn{2}{c}{Biomass $\\hat{X}_i$ (mgN/m$^2$)} \\\\
    \\cmidrule{2-3}
    & LL & UL \\\\
    \\midrule", "\n", file = fo)
    # Body
    cat("Epilithon &", wr(sss$ll$epi), "&", wr(sss$ul$epi), "\\\\", "\n", file = fo)
    cat("\\quad active &", wr(sss$ll$epi.act), "&", wr(sss$ul$epi.act), "\\\\", "\n", file = fo)
    cat("\\quad refractory &", wr(sss$ll$epi.refr), "&", wr(sss$ul$epi.refr), "\\\\", "\n", file = fo)
    cat("CBOM &", wr(sss$ll$CBOM), "&", wr(sss$ul$CBOM), "\\\\", "\n", file = fo)
    cat("\\quad active &", wr(sss$ll$CBOM.act), "&", wr(sss$ul$CBOM.act), "\\\\", "\n", file = fo)
    cat("\\quad refractory &", wr(sss$ll$CBOM.refr), "&", wr(sss$ul$CBOM.refr), "\\\\", "\n", file = fo)
    cat("FBOM &", wr(sss$ll$FBOM), "&", wr(sss$ul$FBOM), "\\\\", "\n", file = fo)
    cat("\\quad active &", wr(sss$ll$FBOM.act), "&", wr(sss$ul$FBOM.act), "\\\\", "\n", file = fo)
    cat("\\quad refractory &", wr(sss$ll$FBOM.refr), "&", wr(sss$ul$FBOM.refr), "\\\\", "\n", file = fo)
    cat("Seston &", wr(sss$ll$seston), "&", wr(sss$ul$seston), "\\\\", "\n", file = fo)
    cat("\\quad active &", wr(sss$ll$seston.act), "&", wr(sss$ul$seston.act), "\\\\", "\n", file = fo)
    cat("\\quad refractory &", wr(sss$ll$seston.refr), "&", wr(sss$ul$seston.refr), "\\\\", "\n", file = fo)
    cat("\\textit{Petrophila} &", wr(sss$ll$petro), "&", wr(sss$ul$petro), "\\\\", "\n", file = fo)
    cat("\\textit{Psephenus} &", wr(sss$ll$pseph), "&", wr(sss$ul$pseph), "\\\\", "\n", file = fo)
    cat("\\textit{Tricorythodes} &", wr(sss$ll$tricor), "&", wr(sss$ul$tricor), "\\\\", "\n", file = fo)
    cat("\\textit{Leptonema} &", wr(sss$ll$lepto), "&", wr(sss$ul$lepto), "\\\\", "\n", file = fo)
    cat("\\textit{Eudaniela} &", wr(sss$ll$eudan), "&", wr(sss$ul$eudan), "\\\\", "\n", file = fo)
    cat("\\textit{Phylloicus} &", wr(sss$ll$phylo), "&", wr(sss$ul$phylo), "\\\\", "\n", file = fo)
    cat("\\textit{Argia} &", wr(sss$ll$arg), "&", wr(sss$ul$arg), "\\\\", "\n", file = fo)
    cat("\\textit{Euthyplocia} &", wr(sss$ll$euthy), "&", wr(sss$ul$euthy), "\\\\", "\n", file = fo)
    # Bottom
    cat("    \\bottomrule 
  \\end{tabular}", "\n", file = fo)
    close(fo)
}

### ** Write the table
tex_steady_biomasses(sss, "tabular-steady-state-biomasses.tex")

### * Table S3: Comparison with Collins 2016 (steady state fluxes and turnover times)

# Fluxes
ssf <- list("ll" = tidySSFluxes(models$run[[lolRow]], models$model[[lolRow]]),
            "ul" = tidySSFluxes(models$run[[uplRow]], models$model[[uplRow]]))
ssf <- lapply(ssf, function(x) x %>% filter(transect == "transect.1"))
ssf <- lapply(ssf, function(x) x %>% pull(network.ssfluxes))
ssf <- lapply(ssf, function(x) {
    o <- x[[1]]
    fluxes <- list()
    for (i in seq_len(nrow(o))) {
        from <- sapply(x, function(z) z$source[i])
        to <- sapply(x, function(z) z$destination[i])
        stopifnot(all(from == o$source[i], na.rm = TRUE))
        stopifnot(all(to == o$destination[i], na.rm = TRUE))
        fluxes[[i]] <- sapply(x, function(z) z$flux[i])
    }
    o$flux <- fluxes
    as_tibble(o)
})
ssf <- lapply(ssf, function(x) {
    o <- list()
    for (comp in c("epi.act", "CBOM.act", "FBOM.act", "seston.act",
                   "petro", "pseph", "tricor", "lepto", "eudan",
                   "phylo", "arg", "euthy")) {
        indices <- which(x$destination == comp)
        traces <- x$flux[indices]
        trace <- apply(do.call(cbind, traces), 1, sum)
        o[[comp]] <- quantile(trace, probs = c(0.025, 0.5, 0.975))
    }
    return(o)
})
ssf_qt <- ssf

# Turnover times
ssp <- list("ll" = tidyOutput(models$run[[lolRow]]) %>% pull(mcmc.parameters),
            "ul" = tidyOutput(models$run[[uplRow]]) %>% pull(mcmc.parameters))
sstpr <- lapply(ssp, function(x) {
    params <- names(x[[1]])
    tpr <- list()
    for (comp in c("epi.act", "CBOM.act", "FBOM.act", "seston.act",
                   "petro", "pseph", "tricor", "lepto", "eudan",
                   "phylo", "arg", "euthy")) {
        indices <- c(grep(paste0("lossRate_", comp), params),
                     grep(paste0("_from_", comp), params))
        traces <- lapply(indices, function(i) {
            sapply(x, function(w) w[i])
        })
        trace <- apply(do.call(cbind, traces), 1, sum)
        # Check for split compartment
        if (paste0("portion_", comp) %in% params) {
            i <- grep(paste0("portion_", comp), params)
            portion <- sapply(x, function(w) w[i])
            trace <- trace * portion
        }
        tpr[[comp]] <- quantile(1/trace, probs = c(0.025, 0.5, 0.975))
    }
    return(tpr)
})
sstpr_qt <- sstpr

### ** Function to build the tex table
tex_steady_uptake_turnover <- function(ssf, sstpr, file = "tabular-steady-state-uptake-turnover.tex") {
    fo <- file(file, open = "w")
    wm <- function(comp) {
        xll <- signif(ssf$ll[[comp]], 3)
        xul <- signif(ssf$ul[[comp]], 3)
        yll <- signif(sstpr$ll[[comp]], 3)
        yul <- signif(sstpr$ul[[comp]], 3)
        paste0(" & ", xll[2], " & ", xul[2], " & ", yll[2], " & ", yul[2], "\\\\ \n")
    }
    wci <- function(comp) {
        xll <- signif(ssf$ll[[comp]], 3)
        xul <- signif(ssf$ul[[comp]], 3)
        yll <- signif(sstpr$ll[[comp]], 3)
        yul <- signif(sstpr$ul[[comp]], 3)
        paste0(" & ", xll[2], " & ", xul[2], " & ", yll[2], " & ", yul[2], "\\\\")
        paste0(" & [", xll[1], "-", xll[3], "] & [", xul[1], "-", xul[3], "] & [",
               yll[1], "-", yll[3], "] & [", yul[1], "-", yul[3], "] \\\\ \n")
    }
    # Top
    cat("  \\begin{tabular}{lcccc}
    \\toprule
    \\multirow{2}{*}{Compartment $i$} & \\multicolumn{2}{c}{Uptake $F_{i,.}$ (mgN m$^{-2}$ day$^{-1}$)} & \\multicolumn{2}{c}{Turnover time $T_i^\\prime$ (days)} \\\\
    \\cmidrule(r){2-3} \\cmidrule(r){4-5}
    & LL & UL & LL & UL \\\\
    \\midrule", "\n", file = fo)
    # Body
    cat("Epilithon", wm("epi.act"), wci("epi.act"), file = fo)
    cat("& \\textit{6.52} & \\textit{24.9} & \\textit{35.71} & \\textit{15.38} \\\\ \n", file = fo)
    cat("CBOM", wm("CBOM.act"), wci("CBOM.act"), file = fo)
    cat("& \\textit{6.52} & \\textit{24.9} & \\textit{35.71} & \\textit{15.38} \\\\ \n", file = fo)
    cat("FBOM", wm("FBOM.act"), wci("FBOM.act"), file = fo)
    cat("& \\textit{128} & \\textit{244} & \\textit{29.41} & \\textit{36.71} \\\\ \n", file = fo)
    cat("Seston", wm("seston.act"), wci("seston.act"), file = fo)
    cat("& \\textit{0.55} & \\textit{0.222} & \\textit{21.28} & \\textit{16.13} \\\\ \n", file = fo)
    cat("\\textit{Petrophila}", wm("petro"), wci("petro"), file = fo)
    cat("& \\textit{0.012} & \\textit{0.111} & \\textit{12.6} & \\textit{11.6} \\\\ \n", file = fo)
    cat("\\textit{Psephenus}", wm("pseph"), wci("pseph"), file = fo)
    cat("& \\textit{0.263} & \\textit{0.459} & \\textit{19.7} & \\textit{30} \\\\ \n", file = fo)
    cat("\\textit{Tricorythodes}", wm("tricor"), wci("tricor"), file = fo)
    cat("& \\textit{0.392} & \\textit{2.83} & \\textit{6.26} & \\textit{3.78} \\\\ \n", file = fo)
    cat("\\textit{Leptonema}", wm("lepto"), wci("lepto"), file = fo)
    cat("& \\textit{0.01} & \\textit{0.296} & \\textit{23.7} & \\textit{10.9} \\\\ \n", file = fo)
    cat("\\textit{Eudaniela}", wm("eudan"), wci("eudan"), file = fo)
    cat("& \\textit{4.3} & \\textit{2.7} & \\textit{173} & \\textit{149} \\\\ \n", file = fo)
    cat("\\textit{Phylloicus}", wm("phylo"), wci("phylo"), file = fo)
    cat("& \\textit{0.11} & \\textit{0.108} & \\textit{7.8} & \\textit{27.3} \\\\ \n",file = fo)
    cat("\\textit{Argia}", wm("arg"), wci("arg"), file = fo)
    cat("& \\textit{0.02} & \\textit{0.02} & \\textit{25.7} & \\textit{13.6} \\\\ \n", file = fo)
    cat("\\textit{Euthyplocia}", wm("euthy"), wci("euthy"), file = fo)
    cat("& \\textit{0.3} & \\textit{0.8} & \\textit{104.2} & \\textit{93.6} \\\\ \n", file = fo)
    # Bottom
    cat("    \\bottomrule    
  \\end{tabular}", "\n", file = fo)
    close(fo)
}

### ** Write the table
tex_steady_uptake_turnover(ssf, sstpr, file = "tabular-steady-state-uptake-turnover.tex")

### * Figure 6: Statistical comparisons of derived parameters between streams

### ** Calculation of traces

ssf <- list("ll" = tidySSFluxes(models$run[[lolRow]], models$model[[lolRow]]),
            "ul" = tidySSFluxes(models$run[[uplRow]], models$model[[uplRow]]))
ssf <- lapply(ssf, function(x) x %>% filter(transect == "transect.1"))
ssf <- lapply(ssf, function(x) x %>% pull(network.ssfluxes))

### ***  Total N flux

totalNflux <- lapply(ssf, function(x) sapply(x, function(z) {
    sum(z$flux[z$source %in% c("NH4", "NO3")])
}))
Nflux.yscale = c(0, max(unlist(totalNflux)))
Nflux.delta.yscale = range(totalNflux[["ul"]] - totalNflux[["ll"]])

### *** NO3/NH4 for epi

NO3_NH4.epi <- lapply(ssf, function(x) sapply(x, function(z) {
    z$destination[is.na(z$destination)] <- "none"
    (z$flux[z$source == "NO3" & z$destination == "epi.act"] /
     z$flux[z$source == "NH4" & z$destination == "epi.act"])
}))
log.NO3_NH4.epi.yscale = range(log(NO3_NH4.epi[["ll"]], 10), log(NO3_NH4.epi[["ul"]], 10))
log.NO3_NH4.epi.delta.yscale = range(log(NO3_NH4.epi[["ul"]], 10) - log(NO3_NH4.epi[["ll"]], 10))

### *** Proportion of CBOM flux going into Eudaniella

PeudanCBOM <- lapply(ssf, function(x) sapply(x, function(z) {
    z$destination[is.na(z$destination)] <- "none"
    (z$flux[z$source == "CBOM.act" & z$destination == "eudan"] /
     (z$flux[z$source == "CBOM.act" & z$destination == "eudan"] +
      z$flux[z$source == "CBOM.act" & z$destination == "phylo"] +
      z$flux[z$source == "CBOM.act" & z$destination == "none"]))
}))
prop.CBOM.to.eudan.yscale = c(0, range(PeudanCBOM[["ll"]], PeudanCBOM[["ul"]])[2])
prop.CBOM.to.eudan.delta.yscale = range(PeudanCBOM[["ul"]] - PeudanCBOM[["ll"]])
logit.PeudanCBOM <- lapply(PeudanCBOM, function(x) {
    log(x/(1-x))
})

### *** Proportion of Petrophila in Argia diet

PargPetro <- lapply(ssf, function(x) sapply(x, function(z) {
    z$destination[is.na(z$destination)] <- "none"
    (z$flux[z$source == "petro" & z$destination == "arg"] /
     (z$flux[z$source == "petro" & z$destination == "arg"] +
      z$flux[z$source == "tricor" & z$destination == "arg"]))
}))
prop.petro.for.argia.yscale = c(0, range(PargPetro[["ll"]], PargPetro[["ul"]])[2])
prop.petro.for.argia.delta.yscale = range(PargPetro[["ul"]] - PargPetro[["ll"]])
logit.PargPetro <- lapply(PargPetro, function(x) {
    log(x/(1-x))
})

### ** Functions for drawing

### *** draw_inGrid

draw_inGrid = function(expr, pos, margins = rep(0, 4)) {
    vpCell = grid::viewport(name = "vpCell", layout.pos.row = pos[1],
                            layout.pos.col = pos[2])
    grid::pushViewport(vpCell)
    vpMargedInCell = grid::viewport(name = "vpMargedInCell",
                                    width = 1 - margins[2] - margins[4], 
                                    height = 1 - margins[1] - margins[3],
                                    x = margins[2], y = margins[1],
                                    just = c("left", "bottom"))
    grid::pushViewport(vpMargedInCell)
    expr # As it appears in system.time, but replicate has a different way of doing it
    grid::upViewport(2)
}

### *** draw_grid

draw_grid = function(dim, rownames = NULL, colnames = NULL) {
    # Calculate widths and heights
    if (!is.null(colnames)) {
        widths = grid::unit(c(2, rep(1, dim[2])),
                            c("lines", rep("null", dim[2])))
        dim[2] = dim[2] + 1
    } else {
        widths = grid::unit(rep(1, dim[2]), rep("null", dim[2]))
    }
    if (!is.null(rownames)) {
        heights = grid::unit(c(2, rep(1, dim[1])),
                             c("lines", rep("null", dim[1])))
        dim[1] = dim[1] + 1
    } else {
        heights = grid::unit(rep(1, dim[1]), rep("null", dim[1]))
    }
    # Apply layout
    layout = grid::grid.layout(nrow = dim[1], ncol = dim[2],
                               widths = widths, heights = heights)
    vpGrid = grid::viewport(name = "vpGrid", layout = layout, clip = "on")
    grid::pushViewport(vpGrid)
    # Add row and column names
    if (!is.null(rownames)) {
        for (i in seq_along(rownames)) {
            vpRowname = grid::viewport(name = "vpRowname",
                                       layout.pos.row = i + ifelse(is.null(colnames), 0, 1),
                                       layout.pos.col = 1)
            grid::pushViewport(vpRowname)
            grid::grid.rect(gp = grid::gpar(fill = "lightgrey"))
            grid::grid.text(rownames[i], rot = 90)
            grid::upViewport(1)
        }
    }
    if (!is.null(colnames)) {
        for (i in seq_along(colnames)) {
            vpColname = grid::viewport(name = "vpColname",
                                       layout.pos.row = 1,
                                       layout.pos.col = i + ifelse(is.null(rownames), 0, 1))
            grid::pushViewport(vpColname)
            grid::grid.rect(gp = grid::gpar(fill = "lightgrey"))
            grid::grid.text(colnames[i])
            grid::upViewport(1)
        }
    }
}

### *** draw_comparison

draw_comparison = function(chains, xLabels = "", panel = "",
                           yTitle = "", deltaTitle = "", log = FALSE, yscale = NULL,
                           yaxisAt = NULL, yaxisLabel = TRUE,
                           xscale = NULL, xaxisAt = NULL, xaxisLabel = TRUE,
                           deltaChain = NULL, cexDeltaTitle = 0.95) {
    cex = 1.3 # For panel names and X top labels
    # Preparation
    if (log) {
        chains = lapply(chains, log)
    }
    prob = 0.95
    if (is.null(yscale)) {
        yscale = range(quantile(unlist(chains[[1]]), probs = c((1-prob)/2, 1 - ((1-prob)/2))),
                       quantile(unlist(chains[[2]]), probs = c((1-prob)/2, 1 - ((1-prob)/2))))
        yscale = extendrange(r = yscale, f = 0.1)
    }
    # Main vp
    layout = grid::grid.layout(nrow = 3, ncol = 1,
        heights = grid::unit(c(2, 2, 1.2), c("lines", "null", "null")))
    vpMain = grid::viewport(name = "vpMain", clip = "on", layout = layout)
    grid::pushViewport(vpMain)
    # Credible intervals vp
    plotAreaLayout = grid::grid.layout(nrow = 2, ncol = 2,
        widths = grid::unit(c(5, 1), c("lines", "null")),
        heights = grid::unit(c(1, 1), c("null", "lines")))
    vpCredibleIntervals = grid::viewport(name = "vpCredibleIntervals",
        layout = plotAreaLayout,
        layout.pos.row = 2)
    grid::pushViewport(vpCredibleIntervals)
    vpCredibleIntervalsGraph = grid::dataViewport(xscale = c(0.2, 2.8), yscale = yscale,
        name = "vpCredibleIntervalsGraph",
        layout.pos.row = 1, layout.pos.col = 2)
    grid::pushViewport(vpCredibleIntervalsGraph)
    for (i in c(1, 2)) {
        grid::grid.points(x = i, y = median(unlist(chains[[i]])),
                          default.units = "native", pch = 21,
                          gp = grid::gpar(col = "black", fill = "black"))
        grid::grid.lines(x = rep(i, 2), y = quantile(unlist(chains[[i]]), probs = c((1-prob)/2, 1 - ((1-prob)/2))),
                         default.units = "native", 
                         gp = grid::gpar(lwd = 3))
    }
    grid::grid.lines(y = grid::unit(c(0, 0), "npc"))
    grid::grid.yaxis(at = yaxisAt, label = yaxisLabel)
    grid::upViewport(1)
    vpYTitle = grid::viewport(name = "vpYTitle",
        layout.pos.row = 1, layout.pos.col = 1)
    grid::pushViewport(vpYTitle)
    grid::grid.text(yTitle, x = grid::unit(1.25, "line"), rot = 90)
    grid::upViewport(2)
    # Panel name vp
    panelLayout = grid::grid.layout(nrow = 2, ncol = 2,
        widths = grid::unit(c(2, 1), c("lines", "null")),
        heights = grid::unit(c(2, 1), c("lines", "null")))
    vpPanel = grid::viewport(name = "vpPanel",
        layout = panelLayout, layout.pos.row = 1)
    grid::pushViewport(vpPanel)
    vpPanelName = grid::viewport(name = "vpPanelName", layout.pos.row = 1,
                                 layout.pos.col = 1)
    grid::pushViewport(vpPanelName)
    grid::grid.text(panel, gp = grid::gpar(fontface = "bold", cex = cex))
    grid::upViewport(2)
    # X labels (on top)
    vpXlabelLayout = grid::grid.layout(nrow = 1, ncol = 2,
                                       widths = grid::unit(c(5, 1), c("lines", "null")))
    vpXlabel = grid::viewport(name = "vpXlabel", layout = vpXlabelLayout,
                              layout.pos.row = 1)
    grid::pushViewport(vpXlabel)
    vpXlabelLabels = grid::dataViewport(xscale = c(0.2, 2.8), yscale = yscale,
                                        name = "vpXlabelLabels", layout.pos.col = 2)
    grid::pushViewport(vpXlabelLabels)
    for (i in c(1, 2)) {
        grid::grid.text(x = grid::unit(i, "native"),
                        label = xLabels[i],
                        gp = grid::gpar(cex = cex))
    } 
    grid::upViewport(2)
    # Delta vp
    deltaAreaLayout = grid::grid.layout(nrow = 2, ncol = 2,
        widths = grid::unit(c(5, 1), c("lines", "null")),
        heights = grid::unit(c(1, 4), c("null", "lines")))
    vpDelta = grid::viewport(name = "vpDelta",
        layout = deltaAreaLayout, layout.pos.row = 3)
    grid::pushViewport(vpDelta)
    if (is.null(deltaChain)) {
        delta = unlist(chains[[2]] - chains[[1]])
    } else {
        delta = unlist(deltaChain)
    }
    deltaDensity = density(delta)
    if (is.null(xscale)) {
        xscale = c(-max(abs(deltaDensity$x)), max(abs(deltaDensity$x)))
    }
    yscale = c(0, max(deltaDensity$y))
    vpDeltaGraph = grid::dataViewport(xscale = xscale, yscale = yscale,
        name = "vpDeltaGraph",
        layout.pos.row = 1, layout.pos.col = 2)
    grid::pushViewport(vpDeltaGraph)
    grid::grid.polygon(x = deltaDensity$x, y = deltaDensity$y,
                       default.unit = "native",
                       gp = grid::gpar(col = "darkgrey",
                           fill = adjustcolor("lightgrey", alpha.f = 1)))
    grid::grid.lines(x = grid::unit(c(0, 0), "native"),
                     gp = grid::gpar(lty = 2))
    grid::grid.points(x = grid::unit(median(delta), "native"),
                      y = grid::unit(0.8, "lines"), pch = 21,
                      gp = grid::gpar(col = "black", fill = "black"))
    grid::grid.lines(x = quantile(delta, probs = c((1-prob)/2, 1 - ((1-prob)/2))),
                     y = grid::unit(rep(0.8, 2), "lines"),
                     default.units = "native", 
                     gp = grid::gpar(lwd = 3))
    grid::grid.xaxis(at = xaxisAt, label = xaxisLabel)
    grid::upViewport(1)
    vpXTitle = grid::viewport(name = "vpXTitle",
        layout.pos.row = 2, layout.pos.col = 2)
    grid::pushViewport(vpXTitle)
    grid::grid.text(deltaTitle,
                    y = grid::unit(1, "npc") - grid::unit(3, "lines"),
                    gp = grid::gpar(cex = cexDeltaTitle))
    grid::upViewport(3)
}

### *** drawFigure4

drawFigure4 = function() {
    cellMargins = rep(0, 4)
    grid::grid.newpage()
    draw_grid(c(1, 4))
    draw_inGrid(draw_comparison(totalNflux,
                                xLabels = c("LL", "UL"),
                                panel = "(A)", 
                                yTitle = latex2exp::TeX("$\\textit{F_T} (mgN/m^2/day)$"),
                                deltaTitle = latex2exp::TeX("$(\\textit{F_T})^{\\textit{UL}}-(\\textit{F_T})^{\\textit{LL}}$"),
                                yscale = c(0, 250)),
                margins = cellMargins,
                pos = c(1, 1))
    draw_inGrid(draw_comparison(NO3_NH4.epi,
                                xLabels = c("LL", "UL"),
                                panel = "(B)", 
                                yTitle = latex2exp::TeX("$\\textit{R}_{\\textit{epi},\\textit{NO3}}$"),
                                deltaTitle = latex2exp::TeX("$(\\textit{R}_{\\textit{epi},\\textit{NO3}})^{\\textit{UL}}/(\\textit{R}_{\\textit{epi},\\textit{NO3}})^{\\textit{LL}}$"),
                                yscale = c(log(10^(-4)), log(10)),
                                log = TRUE,
                                yaxisAt = log(c(10^-4, 10^-3, 0.01, 0.1, 1, 10)),
                                yaxisLabel = latex2exp::TeX(c("10^{-4}", "10^{-3}",
                                                              "10^{-2}", "0.1", "1", "10")),
                                xscale = log(c(10^-5, 10^5)),
                                xaxisAt = log(c(10^-4, 10^-2, 1, 10^2, 10^4)),
                                xaxisLabel = latex2exp::TeX(c("10^{-4}", "10^{-2}",
                                                              "1", "10^2", "10^4"))),
                margins = cellMargins,
                pos = c(1, 2))
    draw_inGrid(draw_comparison(PeudanCBOM,
                                xLabels = c("LL", "UL"),
                                panel = "(C)",
                                yTitle = latex2exp::TeX("$\\textit{P}_{\\textit{eud},\\textit{CBOM}}^{\\textit{K}}$"),
                                deltaChain = logit.PeudanCBOM[["ul"]] - logit.PeudanCBOM[["ll"]],
                                deltaTitle = latex2exp::TeX("$\\textit{logit}\\left((\\textit{P}^{\\textit{K}}_{\\textit{eud}, \\textit{CBOM}})^{\\textit{UL}}\\right) - \\textit{logit}\\left((\\textit{P}^{\\textit{K}}_{\\textit{eud}, \\textit{CBOM}})^{\\textit{LL}}\\right)$"),
                                cexDeltaTitle = 0.88, 
                                yscale = c(0, 1)),
                margins = cellMargins,
                pos = c(1, 3))
    draw_inGrid(draw_comparison(PargPetro,
                                xLabels = c("LL", "UL"),
                                panel = "(D)",
                                yTitle = latex2exp::TeX("$\\textit{P}_{\\textit{arg},\\textit{pet}}^{\\textit{U}}$"),
                                deltaChain = logit.PargPetro[["ul"]] - logit.PargPetro[["ll"]],
                                deltaTitle = latex2exp::TeX("$\\textit{logit}\\left((\\textit{P}^{\\textit{U}}_{\\textit{arg}, \\textit{pet}})^{\\textit{UL}}\\right) - \\textit{logit}\\left((\\textit{P}^{\\textit{U}}_{\\textit{arg}, \\textit{pet}})^{\\textit{LL}}\\right)$"),
                                yscale = c(0, 0.5)),
                margins = cellMargins,
                pos = c(1, 4))
}


### ** Save figure

cairo_pdf(file = "derived-parameters-comparison.pdf",
    width = 16, height = 5.1, family = "Century Schoolbook L")
drawFigure4()
dev.off()

### * Figure 5: Contour plot of N uptake

ssf <- list("ll" = tidySSFluxes(models$run[[lolRow]], models$model[[lolRow]]),
            "ul" = tidySSFluxes(models$run[[uplRow]], models$model[[uplRow]]))
ssf <- lapply(ssf, function(x) x %>% filter(transect == "transect.1"))
ssf <- lapply(ssf, function(x) x %>% pull(network.ssfluxes))
z <- lapply(ssf, bind_rows)
z[[1]]$stream <- "LOL"
z[[2]]$stream <- "UPL"
z <- bind_rows(z)
z$destination[is.na(z$destination)] <- "none"
ssf <- z

my.filled.contour<-function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
                                                       length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
          ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
          levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
          col = color.palette(length(levels) - 1), plot.title, plot.axes, 
          key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
          axes = TRUE, frame.plot = axes, ...) 
{
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  #on.exit(par(par.orig))
  #w <- (3 + mar.orig[2L]) * par("csi") * 2.54
  par(las = las)
  #mar <- mar.orig
  #mar[4L] <- mar[2L]
  #mar[2L] <- 1
  #par(mar = mar)
  #mar <- mar.orig
  #mar[4L] <- 1
  .filled.contour(x, y, z, levels, col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}

#' Plot contour plot of confidence bounds
#'
#' Given a 2D scatter of points, plots the 95% quantile contour
#' 
#' @param x, @param y Vectors of equal length indicating the coordinates of the points
#' @param interval the quantile to be plotted. Defaults to 0.95
#' @param ... Parameters passed to \code{\link{contour}}
#' @return A contour plot
#' @export
plot_2dCI<-function(x,y, interval=0.95, fill.col='grey', ...){
    kerneld <- MASS::kde2d(x, y, h=c(30,20),n = 150,lim=c(0,80,0,80))
    pp <- array()
    for (i in 1:length(x)){
        z.x <- max(which(kerneld$x < x[i]))
        z.y <- max(which(kerneld$y < y[i]))
        pp[i] <- kerneld$z[z.x, z.y]
    }
    confidencebound <- quantile(pp, 1-interval, na.rm = TRUE)
    #contour(kerneld, levels = confidencebound, drawlabels=F,...)
    my.filled.contour(kerneld,levels=c(0,confidencebound,1),col=c(NA,fill.col),xlim=c(0,80),ylim=c(0,80),plot.axes={contour(kerneld,levels=confidencebound,drawlabels=FALSE,col='white',add=TRUE,lwd=1.5)},...)
}

#' Plot 95% CI for basal flux MCMC
#' 
#' Plot a contour of the 95% confidence bounds for $NH_4$ and $NO_3$ flux to a basal compartment
#' @param site The stream, one of two: \code{LOL} or \code{UPL}
#' @param compartment The basal compartment to plot (\code{epi},\code{FBOM},\code{CBOM}, or \code{seston}) 
#' @export
plot_basal_contour<-function(site,compartment,...){
  basal <- ssf
  NH4<-subset(basal,source=='NH4' & stream==site & destination==compartment)
  NO3<-subset(basal,source=='NO3' & stream==site & destination==compartment)
  plot_2dCI(NH4$flux, NO3$flux,...)
}

plotContours <- function() {
    par(mar=c(5,5,2,2))
    plot(-10,-10,xlim=c(0,80),ylim=c(0,80),xlab=expression(paste(NH[4],' uptake ',(mgN/m^2/d))),ylab=expression(paste(NO[3],' uptake ',(mgN/m^2/d))),type='n', las=1, xaxs='i', yaxs='i',asp=1)
    streams<-c('LOL','UPL')
    color<-rbind(grey(rep(0.75,4)),grey(rep(0.9,4)))
    comps<-c('CBOM.act','FBOM.act','epi.act','seston.act')
    for(j in 1:4){
        for(i in 2:1){
            plot_basal_contour(site=streams[i],compartment=comps[j],fill.col=color[i,j],lty=i,add=TRUE)
        }
    }
    abline(0,1,lty=2)
    xy <- list(Epilithon = c(17, 11, 9.6, 1),
               CBOM = c(10.8, 35.7, 19.8, 64),
               FBOM = c(34, 9.4, 56, 11.7))
    for (i in seq_along(xy)) {
        text(xy[[i]][3], xy[[i]][4], names(xy)[i], cex = 0.8)
        text(xy[[i]][1], xy[[i]][2], names(xy)[i], cex = 0.8)
    }
    legend(54.5,79,fill=c(grey(0.9),grey(0.75)),legend=c('Open (UL)', 'Closed (LOL)'), title='Canopy Treatment')
}

### ** Save figure

cairo_pdf(file = "contour-N-uptake-basal-compartments.pdf",
    width = 6.5, height = 6.5, family = "Century Schoolbook L")
plotContours()
dev.off()

### * Figure 3: Comparison with Collins 2016 results

### ** Prep data

# This data was entered manually
data <- as_tibble(read.csv('collins2016_fluxComparison.csv'))
data <- data %>%
    mutate(compartment = as.character(compartment),
           stream = as.character(stream))
data$compartment[data$compartment == "Eudaniella"] <- "Eudaniela"
# Updating it with data from runs
ssf <- list("ll" = tidySSFluxes(models$run[[lolRow]], models$model[[lolRow]]),
            "ul" = tidySSFluxes(models$run[[uplRow]], models$model[[uplRow]]))
ssf <- lapply(ssf, function(x) x %>% filter(transect == "transect.1"))
ssf <- lapply(ssf, function(x) x %>% pull(network.ssfluxes))
comps <- ssf[["ll"]][[1]] %>% pull(destination) %>%
    na.omit() %>% unique()
uptakes <- lapply(ssf, function(x) {
    as_tibble(sapply(comps, function(y) {
        sapply(x, function(d) {
            sum(na.omit(d$flux[d$destination == y]))
        })
    }))
})
uptakes <- lapply(uptakes, function(u) {
    o <- tibble(compartment = colnames(u))
    o$uptake <- apply(u, 2, median)
    o$lo_uptake <- apply(u, 2, quantile, probs = 0.025)
    o$hi_uptake <- apply(u, 2, quantile, probs = 0.975)
    return(o)
})
uptakes[["ll"]]$stream <- "LL"
uptakes[["ul"]]$stream <- "UL"
uptakes <- bind_rows(uptakes)
turnovers <- sstpr_qt %>% lapply(function(x) {
    o <- tibble(compartment = names(x),
                quantiles = x)
    o$turnover <- sapply(o$quantiles, "[[", "50%")
    o$lo_turnover <- sapply(o$quantiles, "[[", "2.5%")
    o$hi_turnover <- sapply(o$quantiles, "[[", "97.5%")
    return(o)
})
turnovers[["ll"]]$stream <- "LL"
turnovers[["ul"]]$stream <- "UL"
turnovers <- bind_rows(turnovers) %>%
    select(-quantiles)
x <- left_join(uptakes, turnovers)
mapping <- c("CBOM.act" = "CBOM", "FBOM.act" = "FBOM", "arg" = "Argia",
             "epi.act" = "Epilithon", "eudan" = "Eudaniela",
             "euthy" = "Euthyplocia", "lepto" = "Leptonema",
             "petro" = "Petrophila", "phylo" = "Phylloicus",
             "pseph" = "Psephenus", "seston.act" = "Seston",
             "tricor" = "Tricorythodes")
x$compartment <- mapping[x$compartment]

data <- data %>%
    select(compartment, stream, collins_uptake, collins_turnover) %>%
    left_join(x)

### ** points_fluxcomp()

#' Draw estimates, credible intervals and Collins estimates on an existing plot
#'
#' @param data Data frame with columns "stream", "uptake", "lo_uptake",
#'     "hi_uptake", "collins_uptake", "turnover", "lo_turnover", "hi_turnover",
#'     "collins_turnover" and "y"
#' @param site String, stream ID
#' @param pch_collins pch value for Collins' estimates
#' @param pch0 pch value for this study's estimates
#' @param cex_collins cex for Collins' estimates (default 1)
#' @param cex0 cex for this study's estimates (default 1)
#' @param lty0 lty for credible intervals
#' @param dev Vertical adjustement for drawing (default 0)
#' @param variable String representing the variable to plot, either "flux" or
#'     "turnover"

points_fluxcomp <- function(data, site, pch_collins, pch0,
                            cex_collins = 1, cex0 = 1, lty0, dev = 0, variable){
    # Get the appropriate subset (flux or turnover)
    subdata <- subset(data, stream == site)
    if (variable == 'flux') {
        subdata$x <- subdata$uptake
        subdata$x_lo <- subdata$lo_uptake
        subdata$x_hi <- subdata$hi_uptake
        subdata$x_collins <- subdata$collins_uptake
    }
    else if (variable == 'turn') {
        subdata$x <- subdata$turnover
        subdata$x_lo <- subdata$lo_turnover
        subdata$x_hi <- subdata$hi_turnover
        subdata$x_collins <- subdata$collins_turnover
    }
    # Collins estimates
    points(I(y+dev) ~ log10(x_collins), subdata, pch = pch_collins, cex = cex_collins)
    # isotracerPrototype estimates
    points(I(y+dev) ~ log10(x), subdata, pch = pch0, cex = cex0)
    segments(log10(subdata$x_lo), subdata$y+dev, log10(subdata$x_hi), subdata$y+dev,
             lty = lty0)
}

### ** plot_fluxcomp()

#' Plot estimates, credible intervals and Collins estimates for fluxes and turn overs
#'
#' @param data Data frame with columns "stream", "uptake", "lo_uptake",
#'     "hi_uptake", "collins_uptake", "turnover", "lo_turnover", "hi_turnover",
#'     "collins_turnover", "compartment"
#' @param variable String representing the variable to plot, either "flux" or
#'     "turnover"

plot_fluxcomp <- function(data, variable = 'flux'){
    # Add y coordinates
    data$y <- seq(1,nrow(data)/2)
    # Calculate xlims
    if (variable == 'flux'){
        xlims = c(log10(0.001)-0.2, max(log10(c(data$hi_uptake,data$collins_uptake))+0.5,
                                        na.rm = T))
        titletxt = latex2exp::TeX('Uptake ($mgN/m^2/day$)')
    }
    else if (variable == 'turn'){
        xlims = c(log10(0.5), ceiling(max(log10(c(data$hi_turnover,data$collins_turnover)+100),
                                          na.rm = T)))
        titletxt = latex2exp::TeX('Turnover time (days)')
    }
    # Initialize plot and draw background
    plot(-100, 1000, xaxt = 'n', yaxt = 'n', ylab = '', xlab = '',
         xlim = xlims, ylim = c(0.5, 12.5), yaxs ='i', xaxs = 'i')
    rect(xlims[1]+0.01, data$y[seq(1,14,2)]-0.5, xlims[2]-0.01,
         data$y[seq(1,14,2)]+0.5, col = grey(0.95), border = NA)
    # Axes and Labels
    labels <- c(latex2exp::TeX(data$compartment[1:4]),
                latex2exp::TeX(paste('$\\textit{', data$compartment[5:12], '}$')))
    xaxis <- seq(floor(xlims[1]), round(xlims[2]))
    axis(1, at = xaxis, labels = 10^xaxis, cex.axis = 0.9)
    axis(3, at = xaxis, labels = 10^xaxis, cex.axis = 0.9)
    axis(2, at = seq(0,12), labels = c('', labels), las = 2, cex.axis = 0.9)
    mtext(titletxt, line = 3)
    # Data
    points_fluxcomp(data, site = 'LL', pch_collins = 17, pch0 = 16,
                    lty0 = 1, dev = +0.25, variable = variable)
    points_fluxcomp(data, site = 'UL', pch_collins = 2, pch0 = 1,
                    lty0 = 2, dev = -0.25, variable = variable)
}

### ** Plot

cairo_pdf("figure3-flux-turnover-comparison.pdf",
    width = 10, height = 5.5, family = "Century Schoolbook L")
par(fig = c(0, 0.49, 0.12, 1), mar = c(2.5, 4, 4, 2),
    oma = c(1, 3, 2, 0.5))
plot_fluxcomp(data, variable = 'flux')
par(fig = c(0.51, 1, 0.12, 1), new = TRUE)
plot_fluxcomp(data, variable = 'turn')
par(fig = c(0, 1, 0, 0.12), new = TRUE, mar = c(0, 0, 0, 0))
# Legend
xleg = 0.55
yleg = 0.5
plot(c(xleg-0.12, xleg-0.1), c(yleg-0.34, yleg-0.33),
     pch = c(17, 2), cex = 0.8, xlim = c(0, 1), ylim = c(0, 1),
     xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n')
legend(xleg, yleg , pch = c(16, 1, NA), lty =c (1, 2, NA),
       legend = c('Closed Canopy (LL)', 'Open Canopy (UL)',
                  'Collins (2016) estimates'),
       bty = 'n', inset = 10, xjust = 0.5, yjust = 0.5, cex = 0.9)
dev.off()

### * Figure 4: Sankey plot

### ** Functions (taken from manuscriptMaterial.Rmd in isotracerPrototype)

### *** ribbonSankey()

#' Draw a Sankey plot using constant-width ribbons
#'
#' This function uses the \code{sankey} package to calculate the position of
#' the nodes and edges, and then custom code for the actual drawing.
#'
#' @param edges A data frame with at least two columns: \code{source} and
#'     \code{destination}. Optionally, if a \code{width} column exists, it will
#'     be used as the width of the connecting ribbons. Default is \code{width =
#'     0.2} for all ribbons.
#' @param nodes An optional data frame with at least two columns: \code{id} and
#'     \code{width}. If provided, the width data is used for the boxes
#'     widths. Additionally, if a \code{width2} column is present, it is used
#'     to draw a second box next to the first.
#' @param ribbonTurn Numeric between 0 and 1, determining the steepness of the
#'     ribbon turns next to the boxes
#' @param col.edges Vector of colors used for ribbons fill color. If only one
#'     value is provided, it is recycled for all ribbons.
#' @param alpha.edges Vector of numeric values between 0 and 1, determining the
#'     transparency of the ribbons. If only one value is provided, it is
#'     recycled for all ribbons.
#' @param cex.y Scaling factor applied to the width of all ribbons
#' @param xlim,ylim Optional but useful to specify if one wants to produce
#'     several plots with comparable x and y scales.
#' @param debug Boolean, if \code{TRUE} axes are plotted
#'
#' @return (Invisibly) a list with \code{nodes} and \code{edges} elements
#'
#' @examples
#' set.seed(4)
#' # This "edges" data frame is taken from the help of sankey::make_sankey()
#' edges <- read.table(stringsAsFactors = FALSE, textConnection(
#'     "                get_deps          get_description
#'                      get_deps               parse_deps
#'                      get_deps                     %||%
#'                      get_deps            drop_internal
#'               get_description        pkg_from_filename
#'                    parse_deps                 str_trim
#'                     cran_file             get_pkg_type
#'                     cran_file          r_minor_version
#'                 download_urls split_pkg_names_versions
#'                 download_urls                cran_file
#'                  pkg_download               dir_exists
#'                  pkg_download            download_urls
#'                  pkg_download        filename_from_url
#'                  pkg_download             try_download
#'                       restore             pkg_download
#'                       restore        drop_missing_deps
#'                       restore            install_order
#'                       restore                 get_deps
#'      split_pkg_names_versions               data_frame
#'     "))
#' names(edges) = c("source", "destination")
#' edges$width = runif(nrow(edges), 0.05, 0.4)
#' nodes = data.frame(id = unique(c(edges$source, edges$destination)),
#'                    stringsAsFactors = FALSE)
#' nodes$width = runif(nrow(nodes), 0.2, 0.8)
#' df = data.frame(source = c("A", "B", "C"), destination = c("C", "C", "D"), stringsAsFactors = FALSE)
#' ribbonSankey(edges)
#' ribbonSankey(edges, nodes)
#' ribbonSankey(edges, nodes, cex.y = 3, debug = TRUE)
#' ribbonSankey(edges, nodes, cex.y = 3, debug = TRUE, xlim = c(-2, 3), ylim = c(-40, 10))
#' nodes2 = nodes
#' nodes2$width2 = NA
#' nodes2$width2[c(4, 8)] = c(0.3, 0.7)
#' nodes2$col = viridis::viridis(nrow(nodes2))
#' ribbonSankey(edges, nodes2, col.nodes = nodes2$col,
#'   col.edges = colorRampPalette(RColorBrewer::brewer.pal(6, "Set2"))(nrow(edges)),
#'   cex.y = 4, cex.x = 0.5, ribbonTurn = 0, col2.nodes = "magenta")
#'
#'
#' @export
ribbonSankey = function(edges, nodes = NULL,
                        ribbonTurn = 0.2,
                        col.edges = "lightgrey", alpha.edges = 0.95,
                        border.edges = "black",
                        border.nodes = "black",
                        col.nodes = "darkgrey",
                        col2.nodes = "black",
                        density2.nodes = NULL,
                        cex.x = 1, cex.y = 1,
                        cex.label = 1,
                        xlim = NULL, ylim = NULL,
                        debug = FALSE, main = NULL, sources.left = TRUE,
                        invert.y = FALSE,
                        spacers.x = NULL, y.magnification = 1) {
    # Get nodes and edges geometry
    geometry = prepareRibbonSankey(edges = edges, nodes = nodes,
                                   col.edges = col.edges, alpha.edges = alpha.edges,
                                   border.edges = border.edges, 
                                   col.nodes = col.nodes,
                                   col2.nodes = col2.nodes,
                                   cex.x = cex.x, cex.y = cex.y,
                                   spacers.x = spacers.x)
    sNodes = geometry[["nodes"]]
    sEdges = geometry[["edges"]]
    # Determine which are the pure sources (i.e. only output, no input)
    pureSources = unique(sEdges[["source"]][!sEdges[["source"]] %in% sEdges[["destination"]]])
    # Get labels
    if (is.null(nodes$label)) {
        sNodes$label = sNodes$id
    } else {
        sNodes$label = nodes$label[match(sNodes$id, nodes$id)]
    }
    # Get shading density
    sNodes$density2 = NA
    if (!is.null(density2.nodes)) {
        for (i in seq_len(nrow(sNodes))) {
            if (grepl(".sidebox$", sNodes$id[i])) {
                sNodes$density2[i] = density2.nodes
            }
        }
    }
    # Reorder nodes so that the sideboxes are drawn first
    sNodes = rbind(sNodes[grepl(".sidebox$", sNodes$id), ],
                   sNodes[!grepl(".sidebox$", sNodes$id), ])
    # Reorder edges
    sEdges = sEdges[order(sEdges$source, decreasing = TRUE), ]
    # Create plot
    if (is.null(xlim)) { xlim = range(sNodes$xleft, sNodes$xright) }
    if (is.null(ylim)) { ylim = range(sNodes$ybottom, sNodes$ytop) }
    if (!debug) { par(mar = rep(1, 4)) } else { par(mar = c(3, 3, 1, 1)) }
    plot(0, type = "n", xlim = xlim, ylim = ylim, xlab = "", ylab = "",
         axes = debug, main = main)
    # Invert y
    if (invert.y) {
        invert.y.fun = function(y) {
            y = y - mean(ylim)
            y = -y
            y = y + mean(ylim)
        }
        sNodes[["y"]] = invert.y.fun(sNodes[["y"]])
        sNodes[["ybottom"]] = invert.y.fun(sNodes[["ybottom"]])
        sNodes[["ytop"]] = invert.y.fun(sNodes[["ytop"]])
        tmp = sNodes[["ytop"]]
        sNodes[["ytop"]] = sNodes[["ybottom"]]
        sNodes[["ybottom"]] = tmp
        sEdges[["ystart"]] = invert.y.fun(sEdges[["ystart"]])
        sEdges[["yend"]] = invert.y.fun(sEdges[["yend"]])
    }
    # Draw edges
    for (i in seq_len(nrow(sEdges))) {
        ribbonTurn = ribbonTurn
        rEdge = sEdges[i, ]
        xStart = rEdge$xstart
        xEnd = rEdge$xend
        midX = (xStart + xEnd) / 2
        yStart = rEdge$ystart
        yEnd = rEdge$yend
        width = rEdge$width * y.magnification
        col = rEdge$col
        alpha = rEdge$alpha
        bordercol = rEdge$bordercol
        trajectory = bezierCurve(c(xStart, xStart + ribbonTurn * (midX - xStart), midX,
                                   midX, xEnd - ribbonTurn * (xEnd - midX), xEnd),
                                 c(yStart, yStart, yStart, yEnd, yEnd, yEnd), n = 128)
        polygon(ribbonFromTrajectory(trajectory, width = width),
                col = adjustcolor(col, alpha.f = alpha), border = bordercol)
    }
    # Draw nodes
    for (i in seq_len(nrow(sNodes))) {
        rNode = sNodes[i, ]
        yWidth = (rNode$ytop - rNode$ybottom) * y.magnification
        yMiddle = (rNode$ytop + rNode$ybottom) / 2
        yTop = yMiddle + yWidth/2
        yBottom = yMiddle - yWidth/2
        rect(rNode$xleft, yBottom, rNode$xright, yTop,
             col = rNode$col, density = rNode$density2,
             border = border.nodes)
    }
    # Draw labels
    for (i in seq_len(nrow(sNodes))) {
        rNode = sNodes[i, ]
        yWidth = (rNode$ytop - rNode$ybottom) * y.magnification
        yMiddle = (rNode$ytop + rNode$ybottom) / 2
        yTop = yMiddle + yWidth/2
        yBottom = yMiddle - yWidth/2
        if (!grepl(".sidebox$", rNode$id)) {
            xleft = rNode$xleft
            xright = rNode$xright
            if (paste0(rNode$id, ".sidebox") %in% sNodes$id) {
                xright = sNodes$xright[sNodes$id == paste0(rNode$id, ".sidebox")]
            }
            labelX = mean(c(xleft, xright))
            labelY = yBottom - strheight(rNode$label[[1]]) * cex.label
            if (rNode$id %in% pureSources) {
                labelX = xright - strwidth(rNode$label[[1]]) * cex.label
                labelY = yMiddle
            }
            text(labelX, labelY, labels = rNode$label[[1]], cex = cex.label)
        }
    }
    # Return
    return(invisible(geometry))
}

### *** prepareRibbonSankey()

#' Return a table of nodes and edges to draw a Sankey plot using constant-width ribbons
#'
#' This function uses the \code{sankey} package to calculate the basic
#' positions of the nodes and edges, and then custom code for customizing the
#' coordinates.
#'
#' @param edges A data frame with at least two columns: \code{source} and
#'     \code{destination}. Optionally, if a \code{width} column exists, it will
#'     be used as the width of the connecting ribbons. Default is \code{width =
#'     0.2} for all ribbons.
#' @param nodes An optional data frame with at least two columns: \code{id} and
#'     \code{width}. If provided, the width data is used for the boxes
#'     widths. Additionally, if a \code{width2} column is present, it is used
#'     to draw a second box next to the first.
#' @param col.edges Vector of colors used for ribbons fill color. If only one
#'     value is provided, it is recycled for all ribbons.
#' @param alpha.edges Vector of numeric values between 0 and 1, determining the
#'     transparency of the ribbons. If only one value is provided, it is
#'     recycled for all ribbons.
#' @param cex.y Scaling factor applied to the width of all ribbons
#'
#' @return A list with \code{nodes} and \code{edges} elements
#'
#' @examples
#' set.seed(4)
#' # This "edges" data frame is taken from the help of sankey::make_sankey()
#' edges <- read.table(stringsAsFactors = FALSE, textConnection(
#'     "                get_deps          get_description
#'                      get_deps               parse_deps
#'                      get_deps                     %||%
#'                      get_deps            drop_internal
#'               get_description        pkg_from_filename
#'                    parse_deps                 str_trim
#'                     cran_file             get_pkg_type
#'                     cran_file          r_minor_version
#'                 download_urls split_pkg_names_versions
#'                 download_urls                cran_file
#'                  pkg_download               dir_exists
#'                  pkg_download            download_urls
#'                  pkg_download        filename_from_url
#'                  pkg_download             try_download
#'                       restore             pkg_download
#'                       restore        drop_missing_deps
#'                       restore            install_order
#'                       restore                 get_deps
#'      split_pkg_names_versions               data_frame
#'     "))
#' names(edges) = c("source", "destination")
#' edges$width = runif(nrow(edges), 0.05, 0.4)
#' nodes = data.frame(id = unique(c(edges$source, edges$destination)),
#'                    stringsAsFactors = FALSE)
#' nodes$width = runif(nrow(nodes), 0.2, 0.8)
#' df = data.frame(source = c("A", "B", "C"), destination = c("C", "C", "D"), stringsAsFactors = FALSE)
#' prepareRibbonSankey(edges)
#' prepareRibbonSankey(edges, nodes)
#' nodes2 = nodes
#' nodes2$width2 = NA
#' nodes2$width2[c(4, 8)] = c(1.3, 0.7)
#' prepareRibbonSankey(edges, nodes2)
#'
#'
#' @export
prepareRibbonSankey = function(edges, nodes = NULL,
                               col.edges = "lightgrey", alpha.edges = 0.95,
                               border.edges = "black", 
                        col.nodes = "darkgrey",
                        col2.nodes = "black", cex.x = 1, cex.y = 1,
                        spacers.x = NULL) {
    # Default values
    edgesDefault = list("width" = 0.2, "col" = col.edges,
                        "alpha" = alpha.edges, "bordercol" = border.edges)
    nodesDefault = list("width" = 0.2, "width2" = NA, "col" = col.nodes,
                        "col2" = col2.nodes)
    # Process edges df
    stopifnot(all(c("source", "destination") %in% colnames(edges)))
    if (is.null(edges[["width"]])) {  edges[["width"]] = edgesDefault[["width"]] }
    for (column in c("col", "alpha", "bordercol")) {
            edges[[column]] = edgesDefault[[column]]
    }
    edges = edges[, c("source", "destination", "width", "col", "alpha", "bordercol")]
    edges$id = as.character(interaction(edges$source, edges$destination))
    # Update edge width
    edges$width = edges$width * cex.y
    # Process nodes df
    if (is.null(nodes)) {
        nodes = data.frame("id" = sort(unique(c(edges$source, edges$destination))),
                           stringsAsFactors = FALSE)
    }
    for (column in c("width", "width2")) {
        if (is.null(nodes[[column]])) {
            nodes[[column]] = nodesDefault[[column]]
        }
    }
    for (column in c("col", "col2")) {
            nodes[[column]] = nodesDefault[[column]]
    }
    nodes = nodes[, c("id", "width", "width2", "col", "col2")]
    # Get nodes and edges from sankey::make_sankey
    sankeyLayout = sankey::make_sankey(edges = as.data.frame(edges))
    sNodes = sankeyLayout$nodes
    sNodes = sNodes[, c("id", "x", "center")]
    colnames(sNodes) = c("id", "x", "y")
    sEdges = sankeyLayout$edges
    sEdges = sEdges[, c("source", "destination", "width", "col", "alpha", "bordercol")]
    sEdges$id = as.character(interaction(sEdges$source, sEdges$destination))
    # Check consistency between input edges and Sankey edges
    stopifnot(all(sort(edges$id) == sort(sEdges$id)))
    # Apply x spacers
    if (!is.null(spacers.x)) {
        stopifnot(length(spacers.x) == (length(unique(sNodes$x)) - 1))
        xRange = range(sNodes$x)
        sNodes$x = as.numeric(as.factor(rank(sNodes$x)))
        ext.spacers.x = c(0, spacers.x)
        sNodes$x = sapply(sNodes$x, function(pos) { sum(ext.spacers.x[seq_len(pos)]) })
        back.beta = diff(xRange) / diff(range(sNodes$x))
        back.alpha = xRange[1] - back.beta * range(sNodes$x)[1]
        sNodes$x = back.alpha + back.beta * sNodes$x
    }
    # Calculate the sizes of each node
    sumInputs = aggregate(width ~ destination, data = sEdges, FUN = sum)
    colnames(sumInputs) = c("id", "heightIn")
    sumOutputs = aggregate(width ~ source , data = sEdges, FUN = sum)
    colnames(sumOutputs) = c("id", "heightOut")
    nodeHeights = merge(sumInputs, sumOutputs, all = TRUE)
    nodeHeights$height = mapply(max, nodeHeights$heightIn, nodeHeights$heightOut,
                                na.rm = TRUE)
    sNodes = merge(sNodes, nodeHeights[, c("id", "height")], all = TRUE)
    # Merge node information
    sNodes = merge(sNodes, nodes, all.x = TRUE)
    # Get the plotting order (modified from sankey:::optimal_edge_order)
    edgePlottingOrder = order(-sNodes$x[match(sEdges$source, sNodes$id)],
                              sNodes$y[match(sEdges$destination, sNodes$id)],
                              sNodes$x[match(sEdges$destination, sNodes$id)])
    # Get the detailed coordinates of each box
    sNodes$ybottom = sNodes$y - sNodes$height/2
    sNodes$ytop = sNodes$y + sNodes$height/2
    sNodes$xleft = sNodes$x - sNodes$width/2 * cex.x
    sNodes$xright = sNodes$x + sNodes$width/2 * cex.x
    # Adjust widths
    sNodes$width = sNodes$width * cex.x
    sNodes$width2 = sNodes$width2 * cex.x
    # Add side-boxes if needed
    nBoxes = nrow(sNodes)
    for (i in seq_len(nBoxes)) {
        if (!is.na(sNodes$width2[i])) {
            # Add a side box
            sNodeInfo = sNodes[i, ]
            sNodeInfoLeftBox = sNodeInfo
            totalWidth = sNodeInfo$width + sNodeInfo$width2       
            sNodeInfoLeftBox$xleft = sNodeInfo$x - totalWidth / 2
            sNodeInfoLeftBox$xright = sNodeInfoLeftBox$xleft + sNodeInfo$width
            sNodeInfo$id = paste0(sNodeInfo$id, ".sidebox")
            sNodeInfo$xleft = sNodeInfoLeftBox$xright
            sNodeInfo$xright = sNodeInfo$xleft + sNodeInfo$width2
            sNodeInfo$col = sNodeInfo$col2
            sNodeInfo$col2 = NA
            sNodes[i, ] = sNodeInfoLeftBox
            sNodes = rbind(sNodes, sNodeInfo)
        }
    }
    # Calculate ribbons coordinates
    nodeSidesLeft = list()
    nodeSidesRight = list()
    edgesCoordinates = list()
    for (i in edgePlottingOrder) {
        mySource = sEdges[i, "source"]
        myDestination = sEdges[i, "destination"]
        if (is.null(nodeSidesRight[[mySource]])) {
            sourceLocZero = sNodes$ybottom[sNodes$id == mySource]
        } else {
            sourceLocZero = nodeSidesRight[[mySource]]
        }
        sourceLoc = sourceLocZero + sEdges$width[i] / 2
        nodeSidesRight[[mySource]] = sourceLocZero + sEdges$width[i]
        if (is.null(nodeSidesLeft[[myDestination]])) {
            destinationLocZero = sNodes$ybottom[sNodes$id == myDestination]
        } else {
            destinationLocZero = nodeSidesLeft[[myDestination]]
        }
        destinationLoc = destinationLocZero + sEdges$width[i] / 2
        nodeSidesLeft[[myDestination]] = destinationLocZero + sEdges$width[i]
        edgesCoordinates[[sEdges$id[i]]] = c(sourceLoc, destinationLoc)
    }
    # Center the ribbons coordinates on the box sides
    for (i in seq_len(nrow(sNodes))) {
        nodeId = sNodes$id[i]
        nodeHeight = sNodes$height[i]
        flushLeft = nodeSidesLeft[[nodeId]] - sNodes$ybottom[i]
        flushRight = nodeSidesRight[[nodeId]] - sNodes$ybottom[i]
        if (length(flushLeft)>0 && flushLeft < nodeHeight) {
            vShift = (nodeHeight - flushLeft)/2
            myEdges = sEdges$id[sEdges$destination == nodeId]
            for (e in myEdges) {
                edgesCoordinates[[e]][2] = edgesCoordinates[[e]][2] + vShift
            }
        }
        if (length(flushRight)>0 && flushRight < nodeHeight) {
            vShift = (nodeHeight - flushRight)/2
            myEdges = sEdges$id[sEdges$source == nodeId]
            for (e in myEdges) {
                edgesCoordinates[[e]][1] = edgesCoordinates[[e]][1] + vShift
            }
        }
    }
    # Order the edges
    sEdges = sEdges[edgePlottingOrder, ]
    # Add the edge coordinates to the edge table
    for (i in seq_len(nrow(sEdges))) {
        nodeStart = sEdges[["source"]][i]
        nodeEnd = sEdges[["destination"]][i]
        if (!is.na(sNodes$width2[sNodes$id == nodeStart])) {
            nodeStart = paste0(nodeStart, ".sidebox")
        }
        sEdges$xstart[i] = sNodes$xright[sNodes$id == nodeStart]
        sEdges$xend[i] = sNodes$xleft[sNodes$id == nodeEnd]
        sEdges$ystart[i] = edgesCoordinates[[sEdges$id[i]]][1]
        sEdges$yend[i] = edgesCoordinates[[sEdges$id[i]]][2]
    }
    # Return
    return(list(nodes = sNodes, edges = sEdges))
}

### ** Data prep

# Get steady-state biomasses and fluxes
sss <- list("ll" = tidySSSizes(models$run[[lolRow]], models$model[[lolRow]]),
            "ul" = tidySSSizes(models$run[[uplRow]], models$model[[uplRow]]))
sss <- lapply(sss, function(x) x %>% filter(transect == "transect.1"))
sss <- lapply(sss, function(x) x %>% pull(network.sssizes))
sss <- lapply(sss, function(x) {
    bind_rows(lapply(x, function(y) do.call(tibble, as.list(y))))
})

ssf <- list("ll" = tidySSFluxes(models$run[[lolRow]], models$model[[lolRow]]),
            "ul" = tidySSFluxes(models$run[[uplRow]], models$model[[uplRow]]))
ssf <- lapply(ssf, function(x) x %>% filter(transect == "transect.1"))
ssf <- lapply(ssf, function(x) x %>% pull(network.ssfluxes))
ssf <- lapply(ssf, function(x) {
    o <- x[[1]]
    fluxes <- list()
    for (i in seq_len(nrow(o))) {
        from <- sapply(x, function(z) z$source[i])
        to <- sapply(x, function(z) z$destination[i])
        stopifnot(all(from == o$source[i], na.rm = TRUE))
        stopifnot(all(to == o$destination[i], na.rm = TRUE))
        fluxes[[i]] <- sapply(x, function(z) z$flux[i])
    }
    o$flux <- fluxes
    as_tibble(o)
})

# Process data for Sankey plot
edges <- lapply(ssf, function(f) {
    o <- f %>% mutate(width = map(flux, median)) %>%
        mutate(width = unlist(width)) %>%
        select(-flux) %>%
        na.omit()
    o[["source"]][endsWith(o[["source"]], ".act")] <-
        substr(o[["source"]][endsWith(o[["source"]], ".act")],
               1, nchar(o[["source"]][endsWith(o[["source"]], ".act")])-4)
    o
    o[["destination"]][endsWith(o[["destination"]], ".act")] <-
        substr(o[["destination"]][endsWith(o[["destination"]], ".act")],
               1, nchar(o[["destination"]][endsWith(o[["destination"]], ".act")])-4)
    o %>% mutate(col.edges = ifelse(source == "NH4", grey(0.8),
                             ifelse(source == "NO3", grey(0.2),
                                    grey(0.5))))
})

nodes <- lapply(seq_along(sss), function(k) {
    s <- sss[[k]]
    o <- s %>% apply(2, median) %>%
        enframe() %>%
        rename(id = name, width = value)
    o$width2 <- NA
    for (i in seq_len(nrow(o))) {
        if (endsWith(o$id[i], ".act")) {
            comp <- substr(o$id[i], 1, nchar(o$id[i])-4)
            o$id[i] <- comp
            refr <- paste0(comp, ".refr")
            o$width2[i] <- o$width[which(o$id == refr)]
        }
    }
    o <- o[!endsWith(o$id, ".refr"),]
    o$label <- c("CBOM" = "CBOM", "FBOM" = "FBOM",
                 "NH4" = "NH4", "NO3" = "NO3",
                 "arg" = expression(italic("Argia")),
                 "epi" = "Epilithon",
                 "eudan" = expression(italic("Eudaniella")),
                 "euthy" = expression(italic("Euthyplocia")),
                 "lepto" = expression(italic("Leptonema")),
                 "petro" = expression(italic("Petrophila")),
                 "pseph" = expression(italic("Psephenus")),
                 "phylo" = expression(italic("Phylloicus")),
                 "seston" = "Seston",
                 "tricor" = expression(italic("Tricorythodes")))[o$id]
    # Divide biomass by total flux in so that width is turnover time
    # and area is biomass
    fluxes <- edges[[names(sss)[k]]]
    influxes <- fluxes %>% dplyr::select(id = destination, influx = width) %>%
        group_by(id) %>%
        summarize(influx = sum(influx))
    for (i in 1:nrow(o)) {
        if (o$id[i] %in% c("NH4", "NO3")) {
            # Set width for sources
            o$width[i] <- 1
        } else {
            factor <- influxes$influx[influxes$id == o$id[i]]
            o$width[i] <- o$width[i] / factor
            o$width2[i] <- o$width2[i] / factor
        }
    }
    return(o)
})
names(nodes) <- names(sss)

### ** Plot

### *** 1x

cairo_pdf("sankey-plot.pdf", width = 6, height = 8.5,
          family = "Century SchoolBook L")
par(mfrow = c(2, 1), mar = c(0, 0, 0, 0))
ylim = c(-20, 0)
xlim = c(-0.3, 3.4)
cex.y = 0.08
cex.x = 0.006
debugFlag = FALSE
col.act = "white"
border.nodes = grey(0.6)
col.refr = "darkgrey"
ribbonTurn = 0.4
shadingDensity = 30
cex.label = 0.8
spacers.x = c(50, 40, 25)
mains <- c("ll" = "Closed canopy (LL)",
           "ul" = "Open canopy (UL)")
geoms <- list()
for (stream in c("ll", "ul")) {
    geoms[[stream]] <- ribbonSankey(edges = na.omit(edges[[stream]]),
      nodes = nodes[[stream]],
      col.edges = na.omit(edges[[stream]])$col.edges,
      cex.x = cex.x, cex.y = cex.y, debug = debugFlag,
      xlim = xlim, ylim = ylim,
      border.edges = na.omit(edges[[stream]])$col.edges,
      col.nodes = col.act, density2.nodes = shadingDensity,
      col2.nodes = col.refr, ribbonTurn = ribbonTurn,
      cex.label = cex.label, main = mains[stream],
      invert.y = TRUE, spacers.x = spacers.x,
      border.nodes = border.nodes)
}
dev.off()

### *** 10x

edgesMag <- edges
nodesMag <- nodes
magnification <- 10
zoomed <- c("lepto", "phylo", "eudan", "petro", "pseph", "tricor", "arg",
            "euthy")
for (stream in c("ll", "ul")) {
    for (i in 1:nrow(edgesMag[[stream]])) {
        if (edgesMag[[stream]]$source[i] %in% zoomed |
            edgesMag[[stream]]$destination[i] %in% zoomed) {
            edgesMag[[stream]]$width[i] <- edgesMag[[stream]]$width[i] *
                magnification
        }
    }
}
cairo_pdf("sankey-plot-10x.pdf", width = 6, height = 8.5,
          family = "Century SchoolBook L")
par(mfrow = c(2, 1), mar = c(0, 0, 0, 0))
for (stream in c("ll", "ul")) {
    geoms[[stream]] <- ribbonSankey(edges = na.omit(edgesMag[[stream]]),
      nodes = nodesMag[[stream]],
      col.edges = na.omit(edgesMag[[stream]])$col.edges,
      cex.x = cex.x, cex.y = cex.y, debug = debugFlag,
      xlim = xlim, ylim = ylim,
      border.edges = na.omit(edgesMag[[stream]])$col.edges,
      col.nodes = col.act, density2.nodes = shadingDensity,
      col2.nodes = col.refr, ribbonTurn = ribbonTurn,
      cex.label = cex.label, main = mains[stream],
      invert.y = TRUE, spacers.x = spacers.x,
      border.nodes = border.nodes)
}
dev.off()

### * Move the output files to the manuscript folder

if (F) {
MS_FOLDER <- path.expand("~/work/papers/writing/2018-01-06_15N-andres/current-manuscript")
if (dir.exists(MS_FOLDER)) {
    # Tables
    file.copy(from = "tabular-model-comparison-DIC.tex",
              to = file.path(MS_FOLDER, "tables"), overwrite = TRUE)
    file.copy(from = "tabular-S1.tex",
              to = file.path(MS_FOLDER, "tables"), overwrite = TRUE)
    file.copy(from = "tabular-steady-state-biomasses.tex",
              to = file.path(MS_FOLDER, "tables"), overwrite = TRUE)
    file.copy(from = "tabular-steady-state-uptake-turnover.tex",
              to = file.path(MS_FOLDER, "tables"), overwrite = TRUE)
    # Figures
    file.copy(from = "figure3-flux-turnover-comparison.pdf",
              to = file.path(MS_FOLDER, "figures"), overwrite = TRUE)
    file.copy(from = "contour-N-uptake-basal-compartments.pdf",
              to = file.path(MS_FOLDER, "figures"), overwrite = TRUE)
    file.copy(from = "derived-parameters-comparison.pdf",
              to = file.path(MS_FOLDER, "figures"), overwrite = TRUE)
    file.copy(from = "supfig-01-posterior-predictive-fit-LL.pdf",
              to = file.path(MS_FOLDER, "figures"), overwrite = TRUE)
    file.copy(from = "supfig-01-posterior-predictive-fit-UL.pdf",
              to = file.path(MS_FOLDER, "figures"), overwrite = TRUE)
    file.copy(from = "supfig-02-traces-LL.pdf",
              to = file.path(MS_FOLDER, "figures"), overwrite = TRUE)
    file.copy(from = "supfig-02-traces-UL.pdf",
              to = file.path(MS_FOLDER, "figures"), overwrite = TRUE)
}
}
