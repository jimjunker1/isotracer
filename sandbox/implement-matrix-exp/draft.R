### * Description

# Functions used to prepare data for a Stan run from a networkModel object
#
# Note that names of vector elements, matrix columns and matrix rows are only
# provided for human-readability: once passed to Stan , the model does not rely
# on those names, but rather on the indices system.

### * prepComps()

#' Determine hidden, observed and source compartments
#'
#' @param model A networkModel object
#' @param x A list to be expanded with Stan data
#'
#' @return An updated x list
#'
#' @examples
#' m <- minimodRep()
#' s <- prepComps(model = m)
#'
#' @export

prepComps <- function(model, x = NULL) {
    d <- model$data
    if (is.null(x)) x <- list()
    x <- prepGroups(model, x)
    # Compartment census
    x[["comps_hidden"]] <- colnames(d$topology[[1]])
    x[["comps_merge"]] <- attr(d$topology[[1]], "merge")
    x[["comps_obs"]] <- c(x[["comps_hidden"]][!x[["comps_hidden"]] %in% unlist(x[["comps_merge"]])],
                         names(x[["comps_merge"]]))
    x[["comps_src"]] <- d$tracerInput[[1]]$compartment
    x <- x[c("comps_hidden", "comps_obs", "comps_merge", "comps_src")]
    x[["n_comps_hidden"]] <- length(x[["comps_hidden"]])
    x[["n_comps_obs"]] <- length(x[["comps_obs"]])
    x[["n_comps_src"]] <- length(x[["comps_src"]])
    # Check that the same compartments are sources in every group
    if (nrow(d) > 1) {
        for (j in 2:nrow(d)) {
            other_src <- d$tracerInput[[j]]$compartment
            if (!setequal(x[["comps_src"]], other_src)) {
                stop("Different source compartments between row 1 and ", j, " in model$data")
            }
        }
    }
    return(x)
}

### * prepGroups()

#' Process group information
#'
#' @param model A networkModel object
#' @param x A list to be expanded with Stan data
#'
#' @return An updated x list
#'
#' @examples
#' m <- minimodRep()
#' s <- prepGroups(model = m)
#'
#' @export

prepGroups <- function(model, x = NULL) {
    d <- model$data
    if (is.null(x)) x <- list()
    x[["n_groups"]] <- nrow(d)
    return(x)
}

### * prepParams()

#' Process parameters (including priors and parameter mapping from MCMC to singleNetwork instances)
#'
#' @param model A networkModel object
#' @param x A list to be expanded with Stan data
#'
#' @return An updated x list
#'
#' @examples
#' m <- minimodRep()
#' s <- prepParams(model = m)
#'
#' @export

prepParams <- function(model, x = NULL) {
    d <- model$data
    if (is.null(x)) x <- list()
    x <- prepGroups(model, x)
    nGroups <- x[["n_groups"]]
    tmp <- prepComps(model)
    x[["params_all"]] <- sort(unique(unlist(d$parameterMapping)))
    x[["n_params_all"]] <- length(x[["params_all"]])
    # Prior types (used to map the correct prior in the stan model)
    priorTypes <- c("uniform" = 1, "hcauchy" = 2)
    x[["priors_numcodes"]] <- priorTypes
    # Upper default bounds
    UPPER_ETA <- 5
    UPPER_LOSS_RATE <- 1
    UPPER_UPTAKE_RATE <- 1
    UPPER_UPTAKE_RATE_FROM_FLOW_INPUT <- 500
    UPPER_PORTION <- 1
    # Default parameter piors (uniform priors)
    x[["params_prior_type"]] <- rep(1, length(x[["params_all"]])) # All uniform priors
    x[["params_prior_uniform_lower"]] <- rep(0, length(x[["params_all"]]))
    x[["params_prior_uniform_upper"]] <- as.vector(sapply(x[["params_all"]], function(p) {
        if (startsWith(p, "eta")) {
            return(UPPER_ETA)
        } else if (startsWith(p, "lossRate_")) {
            return(UPPER_LOSS_RATE)
        } else if (startsWith(p, "uptakeRate_from_")) {
            from = strsplit(p, "_")[[1]][3]
            if (from %in% tmp[["comps_src"]]) {
                return(UPPER_UPTAKE_RATE_FROM_FLOW_INPUT)
            } else {
                return(UPPER_UPTAKE_RATE)
            }
        } else if (startsWith(p, "portion_")) {
            return(UPPER_PORTION)
        } else {
            stop("Unknown parameter type: ", p)
        }
    }))
    x[["params_prior_hcauchy_scale"]] <- rep(0.5, length(x[["params_all"]])) # Not used by default
    # Apply customized parameter priors
    for (i in 1:length(x[["params_all"]])) {
        focalParam <- x[["params_all"]][i]
        customPrior <- model$priors[[focalParam]]
        if (!is.null(customPrior) && customPrior[[1]] != "placeholder") {
            priorType <- customPrior[["type"]]
            x[["params_prior_type"]][i] <- priorTypes[priorType]
            if (priorType == "uniform") {
                x[["params_prior_uniform_lower"]][i] <- customPrior$parameters["min"]
                x[["params_prior_uniform_upper"]][i] <- customPrior$parameters["max"]
            } else if (priorType == "hcauchy") {
                x[["params_prior_hcauchy_scale"]][i] <- customPrior$parameters["scale"]
            } else {
                stop("Unknown prior type: ", priorType)
            }
        }
    }
    # Mapping between actual parameters and parameters by types
    x[["params_mapping_paramID"]] <- seq_along(x[["params_all"]])
    x[["params_mapping_priorType"]] <- x[["params_prior_type"]]
    x[["params_mapping_paramID_to_priorID"]] <- sapply(x[["params_mapping_paramID"]], function(i) {
        sum(x[["params_mapping_priorType"]][1:i] == x[["params_mapping_priorType"]][i])
    })
    x[["n_priors_uniform_numcode1"]] <- sum(x[["params_mapping_priorType"]] == 1)
    x[["n_priors_hcauchy_numcode2"]] <- sum(x[["params_mapping_priorType"]] == 2)
    # Mapping between MCMC parameters and singleNetwork parameters
    # Transitions
    transitions <- lapply(seq_len(nGroups), function(g) {
        out <- transition2indices(d$topology[[g]], d$parameterMapping[[g]], x[["params_all"]])
        return(out[, c("from", "to", "parameter")])
    })
    x[["n_params_transitions"]] <- c(sapply(transitions, nrow), 0) # Padding
    names(x[["n_params_transitions"]]) <- c(paste0("grp", 1:nGroups), "padding")
    x[["maxn_params_transitions"]] <- max(x[["n_params_transitions"]])
    x[["params_mapping_transitions"]] <- array(0,
        dim = c(x[["maxn_params_transitions"]], 3, nGroups))
    dimnames(x[["params_mapping_transitions"]])[[2]] <- c("j", "i", "paramID")
    dimnames(x[["params_mapping_transitions"]])[[3]] <- paste0("grp", 1:nGroups)
    for (g in seq_len(nGroups)) {
        x[["params_mapping_transitions"]][1:(x[["n_params_transitions"]][g]), , g] <- as.matrix(transitions[[g]])
    }
    # Losses
    losses <- lapply(seq_len(nrow(d)), function(i) {
        out <- loss2indices(d$topology[[i]], d$parameterMapping[[i]], x[["params_all"]])
        return(out[, c("from", "parameter")])
    })
    x[["n_params_losses"]] <- c(sapply(losses, nrow), 0) # Padding
    names(x[["n_params_losses"]]) <- c(paste0("grp", 1:nGroups), "padding")
    x[["maxn_params_losses"]] <- max(x[["n_params_losses"]])
    x[["params_mapping_losses"]] <- array(0, dim = c(x[["maxn_params_losses"]], 2, nGroups))
    dimnames(x[["params_mapping_losses"]])[[2]] <- c("j", "paramID")
    dimnames(x[["params_mapping_losses"]])[[3]] <- paste0("grp", 1:nGroups)
    for (g in seq_len(nGroups)) {
        x[["params_mapping_losses"]][1:(x[["n_params_losses"]][g]), , g] <- as.matrix(losses[[g]])
    }
    # Eta mapping
    x[["params_mapping_eta"]] <- c(sapply(seq_len(nGroups), function(g) {
        match(d$parameterMapping[[g]][["eta"]], x[["params_all"]])
    }), 0) # Padding
    names(x[["params_mapping_eta"]]) <- c(paste0("grp", 1:nGroups), "padding")
    # Return
    return(x)
}

### * prepCompBuildMats()

#' Prepare the matrices used to build compartments (hidden to obs and reverse)
#'
#' @param model A networkModel object
#' @param x A list to be expanded with Stan data
#'
#' @return An updated x list
#'
#' @examples
#' m <- minimodRep()
#' s <- prepComps(model = m, s)
#' s <- prepParams(model = m
#' s <- prepBuildMat(model = m, s)
#'
#' @export

prepCompBuildMats <- function(model, x) {
    d <- model$data
    if (is.null(x) || !"params_all" %in% names(x) || !"comps_hidden" %in% names(x)) {
        stop("prepComps() and prepParams() must be run before prepCompBuildMats()")
    }
    nGroups <- nrow(d)
    x <- prepGroups(model, x)
    # Matrix to build observed compartments (comps_obs = comp1 + comp2)
    # col1: comps_obs index
    # col2: comp1 index (e.g. pseph, epi.act)
    # col3: comp2 index, if any (e.g. epi.refr), 0 otherwise
    # col4 and following (col_j, j>=4): index for portion parameter in group j-3
    z <- matrix(0, ncol = 3 + nGroups, nrow = x[["n_comps_obs"]])
    for (i in seq_len(x[["n_comps_obs"]])) {
        obs <- x[["comps_obs"]][i]
        z[i, 1] <- i
        if (obs %in% names(x[["comps_merge"]])) {
            z[i, 2] <- match(x[["comps_merge"]][[obs]][1], x[["comps_hidden"]])
            z[i, 3] <- match(x[["comps_merge"]][[obs]][2], x[["comps_hidden"]])
            for (j in 4:(nGroups+4-1)) {
                portionName <- paste0("portion_", x[["comps_merge"]][[obs]][1])
                z[i, j] <- match(d$parameterMapping[[j-4+1]][[portionName]],
                                 x[["params_all"]])
            }
        } else {
            z[i, 2] <- match(obs, x[["comps_hidden"]])
        }
    }
    colnames(z) <- c("comp_obs", "comp_hidden_act", "comp_hidden_refr",
                     paste0("paramID_grp", 1:nGroups))
    x[["buildMatrix_hiddenComp_to_obsComp"]] <- z
    # Matrix to build hidden compartments
    # col1: comps_hidden index (e.g. pseph, epi.act)
    # col2: comps_obs index (e.g. pseph, epi)
    # col3: type (1=unsplit, 2=.act, 3=.refr)
    # col4 and following (col_j, j>=4): index for portion parameter in group j-3
    z <- matrix(0, ncol = 3 + nGroups, nrow = x[["n_comps_hidden"]])
    for (i in seq_len(x[["n_comps_hidden"]])) {
        comp <- x[["comps_hidden"]][i]
        z[i, 1] <- i
        if (comp %in% unlist(x[["comps_merge"]])) {
            if (endsWith(comp, "act")) {
                obsComp <- substr(comp, 1, nchar(comp)-4)
                z[i, 2] <- match(obsComp, x[["comps_obs"]])
                z[i, 3] <- 2
            } else {
                stopifnot(endsWith(comp, "refr"))
                obsComp <- substr(comp, 1, nchar(comp)-5)
                z[i, 2] <- match(obsComp, x[["comps_obs"]])
                z[i, 3] <- 3
            }
            for (j in 4:(nGroups+4-1)) {
                portionName <- paste0("portion_", obsComp, ".act")
                z[i, j] <- match(d$parameterMapping[[j-4+1]][[portionName]],
                                 x[["params_all"]])
            }
        } else {
            z[i, 2] <- match(comp, x[["comps_obs"]])
            z[i, 3] <- 1
        }
    }
    colnames(z) <- c("comp_hidden", "comp_obs", "comp_hidden_type",
                     paste0("paramID_grp", 1:nGroups))
    x[["buildMatrix_obsComp_to_hiddenComp"]] <- z
    # Return
    return(x)
}

### * prepStepSources()

#' Prepare the step information for input flows
#'
#' @param model A networkModel object
#' @param x A list to be expanded with Stan data
#'
#' @return An updated x list
#'
#' @examples
#' m <- minimodRep()
#' s <- prepComps(model = m)
#' s <- prepStepSources(model = m, s)
#'
#' @export

prepStepSources <- function(model, x) {
    d <- model$data
    if (is.null(x) || !"comps_src" %in% names(x)) {
        stop("prepComps() must be run before prepStepSources()")
    }
    nGroups <- nrow(d)
    x <- prepGroups(model, x)
    # Steps
    minTime <- 0
    maxTime <- max(sapply(1:nGroups, function(g) { max(d$tracerObservations[[g]]$time) }))
    # Gather all times at which observations occur
    allSamplingTimes <- sort(unique(unlist(sapply(1:nGroups, function(g) {
        d$tracerObservations[[g]]$time }))))
    # Find input flow step points
    bpsPerGroup <- list()
    for (i in seq_len(nGroups)) {
        inputs <- d$tracerInput[[i]]
        bps <- lapply(unlist(inputs$profile), findBreakpoints,
                      limits = c(0, maxTime), stopIfNull = TRUE)
        bps <- unique(sort(unlist(bps)))
        bpsPerGroup[[i]] <- bps
    }
    maxNBps <- max(sapply(bpsPerGroup, length))
    # Build breaking points table
    bpTable <- matrix(0, ncol = nGroups, nrow = maxNBps)
    colnames(bpTable) <- paste0("grp", 1:nGroups)
    rownames(bpTable) <- paste0("timeStpPt", 1:maxNBps)
    for (i in seq_len(nGroups)) {
        nBps <- length(bpsPerGroup[[i]])
        bpTable[1:nBps, i] <- bpsPerGroup[[i]]
    }
    x[["n_step_points"]] <- c(sapply(bpsPerGroup, length), 0) # Padded for safety
    names(x[["n_step_points"]]) <- c(paste0("grp", 1:nGroups), "padding")
    x[["maxn_step_points"]] <- maxNBps
    x[["step_points_table"]] <- bpTable
    # Process input source compartments
    x[["comps_src_indices"]] <- c(match(x[["comps_src"]], x[["comps_hidden"]]), 0) # Padded for safety
    names(x[["comps_src_indices"]]) <- c(x[["comps_src"]], "padding")
    x[["comps_src_marked_values"]] <- array(0,
        dim = c(x[["maxn_step_points"]], x[["n_comps_src"]], x[["n_groups"]]))
    dimnames(x[["comps_src_marked_values"]])[[1]] <- paste0("stpPt", 1:maxNBps)
    dimnames(x[["comps_src_marked_values"]])[[2]] <- x[["comps_src"]]
    dimnames(x[["comps_src_marked_values"]])[[3]] <- paste0("grp", 1:nGroups)
    x[["comps_src_unmarked_values"]] <- x[["comps_src_marked_values"]]
    epsBp <- maxTime/4096 # Epsilon value to calculate the input values on the right-hand side of a step point
    for (k in seq_len(nGroups)) {
        for (j in seq_len(x[["n_comps_src"]])) {
            focalChar <- x[["comps_src"]][j]
            focalNum <- match(focalChar, d$tracerInput[[k]]$compartment)
            nBps <- x[["n_step_points"]][k]
            bps <- x[["step_points_table"]][1:nBps, k]
            x[["comps_src_marked_values"]][1:nBps, j, k] <- d$tracerInput[[k]]$profile[[focalNum]]$marked(bps+epsBp)
            x[["comps_src_unmarked_values"]][1:nBps, j, k] <- d$tracerInput[[k]]$profile[[focalNum]]$unmarked(bps+epsBp)
        }
    }
    # Return
    return(x)
}

### * prepInitProps()

#' Prepare the initial proportions
#'
#' @param model A networkModel object
#' @param x A list to be expanded with Stan data
#'
#' @return An updated x list
#'
#' @examples
#' m <- minimodRep()
#' s <- prepComps(model = m)
#' s <- prepInitProps(model = m, s)
#'
#' @export

prepInitProps <- function(model, x) {
    d <- model$data
    if (is.null(x) || !"comps_src" %in% names(x)) {
        stop("prepComps() must be run before prepInitProps()")
    }
    nGroups <- nrow(d)
    x <- prepGroups(model, x)
    # Starting proportions
    x[["init_proportions"]] <- lapply(seq_len(nGroups), function(g) {
        d$initProportions[[g]] %>%
            mutate(group = g,
                   obsCompartment = match(compartment, x[["comps_obs"]]))
    }) %>% bind_rows %>%
        select(-compartment) %>%
        tidyr::spread(key = obsCompartment, value = proportion)
    stopifnot(all(x[["init_proportions"]]$group == 1:nGroups))
    x[["init_proportions"]]$group <- NULL
    stopifnot(all(colnames(x[["init_proportions"]]) == 1:x[["n_comps_obs"]]))
    x[["init_proportions"]] <- as.matrix(x[["init_proportions"]])
    colnames(x[["init_proportions"]]) <- x[["comps_obs"]]
    rownames(x[["init_proportions"]]) <- paste0("grp", 1:nGroups)
    x[["init_proportions"]] <- t(x[["init_proportions"]])
    return(x)
}

### * prepCompSizes()

#' Prepare the compartment size data (mean and sd)
#'
#' @param model A networkModel object
#' @param x A list to be expanded with Stan data
#'
#' @return An updated x list
#'
#' @examples
#' m <- minimodRep()
#' s <- prepComps(model = m)
#' s <- prepCompSizes(model = m, s)
#'
#' @export

prepCompSizes <- function(model, x) {
    d <- model$data
    if (is.null(x) || !"comps_src" %in% names(x)) {
        stop("prepComps() must be run before prepCompSizes()")
    }
    nGroups <- nrow(d)
    x <- prepGroups(model, x)
    # Biomass means and sd
    meanBM <- lapply(seq_len(nGroups), function(g) {
        d$compartmentSizes[[g]] %>% na.omit() %>%
            filter(compartment %in% x[["comps_obs"]]) %>%
            group_by(compartment) %>%
            summarize(meanBiomass = mean(size)) %>%
            mutate(group = g,
                   obsCompartment = match(compartment, x[["comps_obs"]])) %>%
            select(-compartment)
    }) %>% bind_rows %>%
    tidyr::spread(key = obsCompartment, value = meanBiomass)
    stopifnot(all(meanBM$group == 1:nGroups))
    meanBM$group <- NULL
    stopifnot(all(colnames(meanBM) == 1:x[["n_comps_obs"]]))
    meanBM <- as.matrix(meanBM)
    colnames(meanBM) <- x[["comps_obs"]]
    rownames(meanBM) <- paste0("grp", 1:nGroups)
    x[["sizes_mean"]] <- t(meanBM)
    sdBM <- lapply(1:nrow(d), function(g) {
        d$compartmentSizes[[g]] %>% na.omit() %>%
            filter(compartment %in% x[["comps_obs"]]) %>%
            group_by(compartment) %>%
            summarize(sdBiomass = sd(size)) %>%
            mutate(group = g,
                   obsCompartment = match(compartment, x[["comps_obs"]])) %>%
            select(-compartment)
    }) %>% bind_rows %>%
    tidyr::spread(key = obsCompartment, value = sdBiomass)
    stopifnot(all(sdBM$group == 1:nGroups))
    sdBM$group <- NULL
    stopifnot(all(colnames(sdBM) == 1:x[["n_comps_obs"]]))
    sdBM[sdBM == 0 | is.na(sdBM)] <- 0.1 # Hack to avoid 0 or NA sd for e.g. source compartments
    sdBM <- as.matrix(sdBM)
    colnames(sdBM) <- x[["comps_obs"]]
    rownames(sdBM) <- paste0("grp", 1:nGroups)
    x[["sizes_sd"]] <- t(sdBM)
    # Return
    return(x)
}

### * prepObservations()

#' Prepare the observed data (proportions)
#'
#' The NA in the time_index columns are filled by prepProjTimes().
#'
#' @param model A networkModel object
#' @param x A list to be expanded with Stan data
#'
#' @return An updated x list
#'
#' @examples
#' m <- minimodRep()
#' s <- prepComps(model = m)
#' s <- prepObservations(model = m, s)
#'
#' @export

prepObservations <- function(model, x) {
    d <- model$data
    if (is.null(x) || !"comps_src" %in% names(x)) {
        stop("prepComps() must be run before prepObservations()")
    }
    nGroups <- nrow(d)
    x <- prepGroups(model, x)
    # Observations
    propObs <- lapply(seq_len(nGroups), function(g) {
        out <- na.omit(d$tracerObservations[[g]] %>%
                       filter(compartment %in% x[["comps_obs"]] & !compartment %in% x[["comps_src"]]))
        out$compartment <- match(out$compartment, x[["comps_obs"]])
        return(out[, c("compartment", "time", "prop")])
    })
    x[["n_observations"]] <- c(sapply(propObs, nrow), 0) # Padding
    names(x[["n_observations"]]) <- c(paste0("grp", 1:nGroups), "padding")
    x[["maxn_observations"]] <- max(x[["n_observations"]])
    x[["props_obs_indep"]] <- array(0, dim = c(x[["maxn_observations"]], 3, nGroups))
    dimnames(x[["props_obs_indep"]])[[2]] <- c("compartment", "time", "time_index")
    dimnames(x[["props_obs_indep"]])[[3]] <- paste0("grp", 1:nGroups)
    x[["props_obs_dep"]] <- matrix(0, ncol = nGroups, nrow = x[["maxn_observations"]])
    colnames(x[["props_obs_dep"]]) <- paste0("grp", 1:nGroups)
    for (g in seq_len(nGroups)) {
        x[["props_obs_indep"]][1:(x[["n_observations"]][g]), 1:2 , g] <- as.matrix(propObs[[g]][, c("compartment", "time")])
        x[["props_obs_indep"]][, 3, g] <- 0
        x[["props_obs_dep"]][1:(x[["n_observations"]][g]), g] <- unlist(propObs[[g]][, c("prop")])
    }
    # Return
    return(x)
}

### * prepProjTimes()

#' Prepare the times at which projections must be calculated
#'
#' @param model A networkModel object
#' @param x A list to be expanded with Stan data
#'
#' @return An updated x list
#'
#' @examples
#' m <- minimodRep()
#' s <- prepComps(model = m)
#' s <- prepStepSources(model = m, s)
#' s <- prepObservations(model = m, s)
#' s <- prepProjTimes(model = m, s)
#'
#' @export

prepProjTimes <- function(model, x) {
    d <- model$data
    if (is.null(x) || !"comps_src" %in% names(x) ||
        !"n_step_points" %in% names(x) || !"n_observations" %in% names(x)) {
        stop("prepComps(), prepStepSources() and prepObservations() must be run ",
             "before prepProjTimes() (or prepDiscreteProjTimes()).")
    }
    nGroups <- nrow(d)
    x <- prepGroups(model, x)
    # Compile sampling and step times
    times <- lapply(seq_len(nGroups), function(g) {
        steps <- x[["step_points_table"]][,g]
        samples <- x[["props_obs_indep"]][, "time", g]
        o <- tibble(time = unique(sort(c(steps, samples))))
        stopifnot(nrow(o) > 1)
        o$step <- as.numeric(o$time %in% steps)
        o$projSteps <- 0
        if (all(o$step < 0.5)) {
            # No source to update
            o$projSteps[1] <- nrow(o) - 1
        } else {
            o$projSteps[1] <- min(which(o$step[2:nrow(o)] > 0.5), nrow(o)-1)
            currentStep <- 1
            nextStep <- currentStep + o$projSteps[currentStep]
            while (nextStep < nrow(o)) {
                currentStep <- nextStep
                nextStep <- min(which(o$step[(currentStep+1):nrow(o)] > 0.5), nrow(o)-currentStep)
                o$projSteps[currentStep] <- nextStep
                nextStep <- nextStep + currentStep
            }
        }
        return(o)
    })
    x[["n_proj_times"]] <- c(sapply(times, nrow), 0) # Padding
    names(x[["n_proj_times"]]) <- c(paste0("grp", seq_len(nGroups)), "padding")
    stopifnot(all(x[["n_proj_times"]][1:nGroups] > 1)) # To ensure that there is at least one point after the first one
    x[["maxn_proj_times"]] <- max(x[["n_proj_times"]])
    projTable <- array(0, dim = c(x[["maxn_proj_times"]], 3, nGroups))
    dimnames(projTable)[[2]] <- c("time", "updateInitState", "projectToFurtherSteps")
    dimnames(projTable)[[3]] <- paste0("grp", seq_len(nGroups))
    for (g in seq_len(nGroups)) {
        projTable[1:(x[["n_proj_times"]][g]), , g] <- as.matrix(times[[g]])
    }
    x[["proj_table"]] <- projTable
    # Split projTable into real and int parts
    x[["proj_table_time"]] <- x[["proj_table"]][,1,]
    x[["proj_table_steps"]] <- x[["proj_table"]][,2:3,]
    # Convert n_proj_times to the number of pieces of step functions
    x[["n_proj_steps"]] <- c(sapply(1:nGroups, function(g) {
        sum(x[["proj_table"]][,3,g] > 0)
    }), 0) # Padding
    names(x[["n_proj_steps"]]) <- c(paste0("grp", seq_len(nGroups)), "padding")
    stopifnot(all(x[["n_proj_steps"]][1:nGroups] > 0)) # To ensure that there is at least one piece
    x[["maxn_proj_steps"]] <- max(x[["n_proj_steps"]])
    # Update observation tables
    for (g in seq_len(nGroups)) {
        nObs <- x[["n_observations"]][g]
        x[["props_obs_indep"]][1:nObs, 3, g] <- match(x[["props_obs_indep"]][1:nObs, "time", g],
                                                      x[["proj_table"]][, "time", g])
    }
    # Return
    return(x)
}

### * prepDiscreteProjTimes()

#' Prepare the times at which projections must be calculated, using dt steps
#'
#' @param model A networkModel object
#' @param x A list to be expanded with Stan data
#' @param minGrid Minimum grid size (smallest grid size used)
#' @param f Multiplicative factor used to produce grids from minGrid (should be
#'     ordered in increasing order)
#'
#' @return An updated x list
#'
#' @examples
#' m <- minimodRep()
#' s <- prepComps(model = m)
#' s <- prepStepSources(model = m, s)
#' s <- prepObservations(model = m, s)
#' s <- prepDiscreteProjTimes(model = m, s, minGrid = 128, f = 2^(0:4))
#'
#' @export

prepDiscreteProjTimes <- function(model, x, minGrid, f) {
    d <- model$data
    if (is.null(x)) { x <- list() }
    x <- prepProjTimes(model, x)
    nGroups <- nrow(d)
    # Flag times to be returned
    dimnames(x$proj_table)[[2]][3] <- "returnProj"
    for (g in 1:nGroups) {
        x$proj_table[1:x$n_proj_times[g], 3, g] <- 1
    }
    # Discretize timeline
    minDt <- max(x$proj_table_time) / minGrid
    dt <- minDt / f
    x[["n_proj_discrete_dt_values"]] <- length(dt)
    nDiscreteTimesteps <- array(0, dim = c(nGroups, length(dt)))
    dimnames(nDiscreteTimesteps)[[1]] <- paste0("grp", seq_len(nGroups))
    dimnames(nDiscreteTimesteps)[[2]] <- paste0("dt", dt)
    o0 <-list()
    for (i in seq_along(dt)) {
        o1 <- list()
        for (g in 1:nGroups) {
            seed <- x$proj_table[1:x$n_proj_times[g],,g]
            stopifnot(nrow(seed) > 1)
            for (j in 1:(nrow(seed)-1)) {
                breaks <- seq(seed[j,1], seed[j+1,1], by = dt[i])
                if (breaks[length(breaks)] == seed[j+1,1]) {
                    stopifnot(length(breaks)>1)
                    breaks <- breaks[1:(length(breaks)-1)]
                }
                if (length(breaks)>1) {
                    breaks <- breaks[2:length(breaks)]
                }
                expansion <- array(0, dim = c(length(breaks), 3))
                expansion[,1] <- breaks
                seed <- rbind(seed, expansion)
            }
            seed <- seed[order(seed[,1]), ]
            seed <- cbind(seed, dt = c(diff(seed[,1]), 0))
            o1[[g]] <- list(discreteTimeline = seed,
                            nDiscreteTimesteps = nrow(seed)-1)
            nDiscreteTimesteps[g, i] <- nrow(seed)-1
        }
        o0[[i]] <- list(timelinePerGroup = o1,
                        maxnDiscreteTimesteps = max(sapply(o1, "[[",
                                                           "nDiscreteTimesteps")))
    }
    x[["n_proj_discrete_timesteps"]] <- nDiscreteTimesteps
    maxnDiscreteTimesteps <- max(sapply(o0, "[[", "maxnDiscreteTimesteps"))
    x[["maxn_proj_discrete_timesteps"]] <- maxnDiscreteTimesteps
    # Pool all timelines into a single array
    discreteTimeline <- array(0, dim = c(maxnDiscreteTimesteps+1, 4, nGroups,
                                         length(dt)))
    for (i in seq_along(dt)) {
        for (g in seq_len(nGroups)) {
            nTimepoints <- o0[[i]]$timelinePerGroup[[g]]$nDiscreteTimesteps + 1
            timeline <- o0[[i]]$timelinePerGroup[[g]]$discreteTimeline
            stopifnot(nTimepoints == nrow(timeline))
            discreteTimeline[1:nTimepoints,,g,i] <- timeline
        }
    }
    dimnames(discreteTimeline)[[2]] <- c("time", "updateInitState", "returnProj", "dtToNext")
    dimnames(discreteTimeline)[[3]] <- paste0("grp", seq_len(nGroups))
    dimnames(discreteTimeline)[[4]] <- paste0("dt", dt)
    x[["proj_discrete_table_time"]] <- discreteTimeline
    # Get projections to return
    n_proj_discrete_returned_proj <- array(0, dim = c(nGroups, length(dt)))
    dimnames(n_proj_discrete_returned_proj)[[1]] <- paste0("grp", seq_len(nGroups))
    dimnames(n_proj_discrete_returned_proj)[[2]] <- paste0("dt", dt)
    for (i in seq_along(dt)) {
        for (g in 1:nGroups) {
            n_proj_discrete_returned_proj[g,i] <- sum(discreteTimeline[,3,g,i])
        }
    }
    maxn_proj_discrete_returned_proj <- max(n_proj_discrete_returned_proj)
    proj_discrete_returned_proj <- array(0, dim = c(maxn_proj_discrete_returned_proj, nGroups, length(dt)))
    dimnames(proj_discrete_returned_proj)[[2]] <- paste0("grp", seq_len(nGroups))
    dimnames(proj_discrete_returned_proj)[[3]] <- paste0("dt", dt)
    for (i in seq_along(dt)) {
        for (g in 1:nGroups) {
            nProj <- n_proj_discrete_returned_proj[g,i]
            proj_discrete_returned_proj[1:nProj,g,i] <- which(discreteTimeline[,3,g,i] == 1)
        }
    }
    x[["n_proj_discrete_returned_proj"]] <- n_proj_discrete_returned_proj
    x[["maxn_proj_discrete_returned_proj"]] <- maxn_proj_discrete_returned_proj
    x[["proj_discrete_returned_proj"]] <- proj_discrete_returned_proj
    # Split table into integers and reals
    x[["proj_discrete_table_time_int"]] <- x[["proj_discrete_table_time"]][, c(2,3),,]
    x[["proj_discrete_table_time_real"]] <- x[["proj_discrete_table_time"]][, c(1,4),,]
    # Return
    return(x)
}

### * findBreakpoints

#' Find breakpoints in input flow profile
#'
#' This function is useful when applied to step functions with one or several
#' breakpoints separating the steps.
#'
#' @param f Function f(t) giving the profile as a function of time t
#'     (vectorized)
#' @param limits c(tmin, tmax) giving the time limits of the experiment
#' @param ngrid Number of grid points used for determining breakpoints (default
#'     to 4096)
#' @param maxBreakpoints Maximum number of breakpoints accepted before the
#'     function returns NULL (which would mean that the input flow is not
#'     suitable for a description with breakpoints separating constant
#'     plateaus)
#' @param stopIfNull Boolean, raise an error instead of returning NULL?
#'
#' @return A vecto containing the t values for breakpoints, or NULL if the
#'     profile is deemed inappropriate for this approach
#'
#' @examples
#' f <- function(t) ifelse(t < 10, 0.1, 0.02)
#' findBreakpoints(f, limits = c(0, 40))
#'
#' f <- function(t) ifelse(t < 10, 0.1, ifelse(t > 30, 0.01, 0.02))
#' findBreakpoints(f, limits = c(0, 40))
#'
#' f <- function(t) sin(t)
#' findBreakpoints(f, limits = c(0, 40))
#'
#' f <- function(t) sin(t)
#' findBreakpoints(f, limits = c(0, 40), stopIfNull = TRUE)
#'
#' @export

findBreakpoints <- function(f, limits, ngrid = 4096, maxBreakpoints = 5,
                            stopIfNull = FALSE) {
    # First pass
    grid <- seq(limits[1], limits[2], length.out = ngrid)
    values <- f(grid)
    steps <- which(diff(values) != 0)
    if (length(steps) > maxBreakpoints) {
        if (stopIfNull) { stop("findBreakpoints would return NULL") }
        return(NULL)
    }
    if (length(steps) == 0) { return(limits[1]) }
    intervalLeft <- grid[steps]
    intervalRight <- intervalLeft + (grid[2] - grid[1])
    # Second pass
    t <- c()
    for (i in seq_along(intervalLeft)) {
        grid <- seq(intervalLeft[i], intervalRight[i], length.out = ngrid)
        values <- f(grid)
        steps <- which(diff(values) != 0)
        if (length(steps) > maxBreakpoints) {
            if (stopIfNull) { stop("findBreakpoints would return NULL") }
            return(NULL)
        }
        tInterval <- grid[steps] + 1/2 * (grid[2] - grid[1])
        t <- c(t, tInterval)
    }
    out <- unique(sort(c(limits[1], t)))
    if (length(out) > maxBreakpoints) {
        if (stopIfNull) { stop("findBreakpoints would return NULL") }
        return(NULL)
    }
    return(out)
}

### * networkModelToStanDataExpm()

#' Convert a network model to stan data, to run with matrix exponential model
#'
#' For now this only works for sources of type input flows.
#'
#' @param model A network model built with \code{link{networkModel}}
#' 
#' @return A list that can be used as stan data to run with `stanmodels$networkModel`
#'
#' @examples
#' m <- minimod()
#' stan.data <- networkModelToStanDataExpm(m)
#' 
#' @export

networkModelToStanDataExpm <- function(model, rel_tol = 1e-6, abs_tol = 1e-6,
                                       max_steps = 1e6, solver = "rk45") {
    d <- model$data
    s <- list()
    # Compartments
    s <- prepComps(model, s)
    # Parameters (priors and mapping)
    s <- prepParams(model, s)
    # Matrices to build compartments
    s <- prepCompBuildMats(model, s)
    # Steps for input flow sources
    s <- prepStepSources(model, s)
    # Initial proportions
    s <- prepInitProps(model, s)
    # Compartment sizes
    s <- prepCompSizes(model, s)
    # Observations
    s <- prepObservations(model, s)
    # Projection times
    s <- prepProjTimes(model, s)
    # ODE tuning
    s[["rel_tol"]] <- rel_tol
    s[["abs_tol"]] <- abs_tol
    s[["max_steps"]] <- max_steps
    stopifnot(solver %in% c("rk45", "bdf"))
    s[["solver_selector"]] <- ifelse(solver == "rk45", 1, 2)
    # Return
    s <- s[sort(names(s))]
    return(s)
}
