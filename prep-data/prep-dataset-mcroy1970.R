### * Description

# Script to import into R the dataset from the article:
#
# McRoy, C. Peter, and Robert J. Barsdate. “Phosphate Absorption in Eelgrass1.”
# Limnology and Oceanography 15, no. 1 (January 1, 1970):
# 6–13. https://doi.org/10.4319/lo.1970.15.1.0006.

### * Setup

suppressMessages(library(tidyverse))
suppressMessages(library(here))

### * Run

### ** Tracer data

# Data is acquired using the digitize package on Figure 2 of the original
# publication.

x <- read.table(file.path(here::here(), "prep-data", "mcroy-1970-fig2-topleft.tsv"),
                header = TRUE, sep = "\t")
x$tracerAddition <- "leavesCompartment"
x$lightCondition <- "dark"
x2 <- read.table(file.path(here::here(), "prep-data", "mcroy-1970-fig2-topright.tsv"),
                 header = TRUE, sep = "\t")
x2$tracerAddition <- "leavesCompartment"
x2$lightCondition <- "light"
x3 <- read.table(file.path(here::here(), "prep-data", "mcroy-1970-fig2-bottomleft.tsv"),
                 header = TRUE, sep = "\t")
x3$tracerAddition <- "rootsCompartment"
x3$lightCondition <- "dark"
x4 <- read.table(file.path(here::here(), "prep-data", "mcroy-1970-fig2-bottomright.tsv"),
                 header = TRUE, sep = "\t")
x4$tracerAddition <- "rootsCompartment"
x4$lightCondition <- "light"
tracer = as_tibble(rbind(x, x2, x3, x4))

# Check that we can reproduce Figure 2
library(ggplot2)
ggplot(tracer, aes(x = x, y = log10(y), col = as.factor(z))) +
    geom_line() +
    facet_wrap(~ tracerAddition * lightCondition)

# Convert numeric codes to compartments
compartments = c("upperWater", "leafTip", "leafMiddle", "leafBase", "stem",
                 "rhizome", "root", "lowerWater")
comp_side <- c("leafTip" = "leavesAndStem", "leafMiddle" = "leavesAndStem",
               "leafBase" = "leavesAndStem",
               "stem" = "leavesAndStem", "rhizome" = "rootsAndRhizome",
               "root" = "rootsAndRhizome",
               "upperWater" = "upperWater", "lowerWater" = "lowerWater")
tracer$compartment = compartments[tracer$z]
tracer$main_comp <- comp_side[tracer$compartment]
tracer$z = NULL

# Rename x and y columns
tracer$cpm.per.mg = tracer$y
tracer$time.min = tracer$x
tracer$x = NULL
tracer$y = NULL

# Reorder columns
tracer = tracer[, c("lightCondition", "tracerAddition", "compartment",
                    "main_comp", "time.min", "cpm.per.mg")]

# Convert cpm/mg to n of phosphorus 32
hl = 14.268 * 24 * 3600 # 32P half-life in seconds
tracer$phosphorus32.n.per.mg = tracer$cpm.per.mg / 60 * hl / log(2)

# Summarized tracer data
summ_tracer <- tracer %>%
    group_by(lightCondition, tracerAddition, main_comp, time.min) %>%
    summarize(n.32P.per.mg = mean(phosphorus32.n.per.mg), .groups = "drop")

ggplot(summ_tracer, aes(x = time.min, y = log10(n.32P.per.mg),
                        col = as.factor(main_comp))) +
    geom_line() +
    facet_wrap(~ tracerAddition * lightCondition)

### ** Biomass data

# Data is taken from Table 1 in the original paper. Experimental containers had
# 160 cc of seawater in the upper compartment and a plant and 80 cc of seawater
# in the lower compartment. This data represents the distribution of phosphorus
# in initial conditions.

# Based on comparison with data from Risgaard-Petersen 1998, I assume that the
# biomasses for tissues are given in dry weight. I also assume that this is
# also the case for the cpm/mg data (i.e. cpm/mg of dry weight).

compartment = rep(rep(c("upperWater", "leavesAndStem", "rootsAndRhizome",
                        "lowerWater"), each = 2), 2)
lightCondition = rep(c("light", "dark"), 8)
tracerAddition = rep(c("leavesCompartment", "rootsCompartment"), each = 8)
mass.g = c(160, 160, 0.084, 0.046, 0.038, 0.034, 80, 80,
           160, 160, 0.071, 0.063, 0.068, 0.025, 80, 80)
phosphorus.ug = c(4, 4, 660, 587, 340, 124, 2, 2,
                  4, 4, 780, 460, 190, 170, 2, 2)
phosphorusMW = 31 # g/mol
avogadro = 6.022e23 # mol-1
phosphorus.n = phosphorus.ug * 1e-6 / phosphorusMW * avogadro

biomass = tibble(lightCondition = lightCondition,
                 tracerAddition = tracerAddition,
                 compartment = compartment,
                 mass.g = mass.g,
                 phosphorus.ug = phosphorus.ug,
                 phosphorus.n = phosphorus.n)

### ** Save the data
## mcroy1970 = list(biomass = biomass, tracer = tracer)
## save(mcroy1970, file = file.path(here::here(), "data", "mcroy1970.rda"),
##      version = 2)

### * Convert to 32P size data

# Take into account compartment sizes to convert 32P/mg to 32P/compartment
summ_biomass <- biomass %>%
    select(lightCondition, tracerAddition, compartment, mass.g) %>%
    rename(main_comp = compartment) %>%
    mutate(mass.mg = 1000 * mass.g) %>%
    select(-mass.g)

full_summ <- full_join(summ_tracer, summ_biomass,
                       by = c("lightCondition", "tracerAddition", "main_comp")) %>%
    mutate(n.32P = n.32P.per.mg * mass.mg)

ggplot(full_summ, aes(x = time.min, y = log10(n.32P),
                      col = as.factor(main_comp))) +
    geom_line() +
    facet_wrap(~ tracerAddition * lightCondition)

### * Save data

# Final polish
eelgrass <- full_summ %>%
    ungroup() %>%
    rename(compartment = main_comp) %>%
    rename(light_treatment = lightCondition) %>%
    rename(addition_site = tracerAddition) %>%
    mutate(addition_site = recode(addition_site,
                                  leavesCompartment = "upper",
                                  rootsCompartment = "lower")) %>%
    mutate(compartment = recode(compartment,
                                leavesAndStem = "leaves_stem",
                                rootsAndRhizome = "roots_rhizome",
                                lowerWater = "lower_water",
                                upperWater = "upper_water")) %>%
    rename(time_min = time.min,
           n_32P_per_mg = n.32P.per.mg,
           mass_mg = mass.mg,
           n_32P = n.32P)

# Save data
save(eelgrass, file = file.path(here::here(), "data", "eelgrass.rda"),
     version = 3, compress = "xz")
