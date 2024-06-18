### * Description

# Script to profile code

### * Setup

library(tidyverse)
library(profvis)
library(microbenchmark)
devtools::load_all()

library(furrr)
plan(multiprocess)

### * tidy_trajectories()

m <- aquarium_mod
r <- aquarium_run
plot(r)

system.time({tt <- tidy_trajectories(m, r, 200, use_cache = FALSE) })
system.time({tt <- tidy_trajectories(m, r, 100, use_cache = TRUE) })

n <- c(10, 20, 40, 100, 200)
t_no_cache <- c(1.07, 1.97, 3.57, 8.89, 14.01)
t_cache <- c(0.77, 1.49, 2.54, 6, 11.9)
t_cache_new_sp <- c(NA, NA, NA, 4.25, 8.45)
t_cache_quick_sp <- c(NA, NA, NA, 4.29, 8.15)
t_cache_parallel <- c(NA, NA, NA, 2.5, 4.35)

t_cache/t_no_cache
t_cache_new_sp/t_no_cache
t_cache_parallel/t_no_cache

plot(n, t_no_cache, type = "b", col = "red")
lines(n, t_cache, type = "b", col = "pink")
lines(n, t_cache_new_sp, type = "b", col = "magenta")

p <- profvis({
    tt <- tidy_trajectories(m, r, 400)
})

print(p)

# Points where speed can be increased:
#
# - cache the results of encode_events() and nm_row_get_time_scheme() when
#   running project() and project_row()
# - see if using df instead of tbl can accelerate things at the end of
#   project_row()
# - try to speed up set_params()

### * set_params()

set_params_new <- function(nm, params, force = TRUE, quick = FALSE) {
    if (is.vector(params)) {
        params <- data.frame(value = params, parameter = names(params))
    }
    prev_params <- attr(nm, "parameterValues")
    if (!is.null(prev_params) & !force) {
        kept_indices <- !(prev_params$parameter %in% params$parameter)
        kept_params <- prev_params[kept_indices, ]
        new_params <- dplyr::bind_rows(kept_params, params)
    } else {
        new_params <- params
    }
    if (!quick) {
        new_params <- new_params[order(new_params[["parameter"]]), ]
        attr(nm, "parameterValues") <- tibble::as_tibble(new_params)
    } else {
        attr(nm, "parameterValues") <- new_params
    }
    for (i in seq_len(nrow(nm))) {
        nm[["parameters"]][[i]]$value <- new_params[["value"]][match(nm[["parameters"]][[i]]$in_model,
                                                                     new_params[["parameter"]])]
    }
    return(nm)
}

set_params_ref <- function(nm, params, force = TRUE) {
    if (is.vector(params)) {
        params <- tibble::tibble(value = params, parameter = names(params))
    }
    prev_params <- attr(nm, "parameterValues")
    if (!is.null(prev_params) & !force) {
        kept_indices <- !(prev_params$parameter %in% params$parameter)
        kept_params <- prev_params[kept_indices, ]
        new_params <- dplyr::bind_rows(kept_params, params)
    } else {
        new_params <- params
    }
    new_params <- new_params[order(new_params[["parameter"]]), ]
    attr(nm, "parameterValues") <- new_params
    for (i in seq_len(nrow(nm))) {
        nm[["parameters"]][[i]]$value <- new_params[["value"]][match(nm[["parameters"]][[i]]$in_model,
                                                                     new_params[["parameter"]])]
    }
    return(nm)
}


# Test
m <- aquarium_mod
p <- sample_params(m)

microbenchmark(
    set_params_ref(m, p),
    set_params_new(m, p, quick = TRUE),
    set_params_new(m, p, quick = FALSE)
)

prof <- profvis(replicate(500, set_params_new(m, p)))
print(prof)

### * Microbenchmarks

library(data.table)
relig_income <- as.data.frame(relig_income)

microbenchmark(
    melt(as.data.table(relig_income), id.vars = "religion"),
    pivot_longer(relig_income, -religion, names_to = "income", values_to = "count")
)

z1 <- melt(as.data.table(relig_income), id.vars = "religion",
           variable.name = "income", value.name = "count")
z2 <- pivot_longer(relig_income, -religion, names_to = "income", values_to = "count")

z <- setNames(1:10, letters[1:10])
w <- data.frame(v = z, p = names(z))

microbenchmark(
    tibble(v = z, p = names(z)),
    data.frame(v = z, p = names(z)),
    as_tibble(data.frame(v = z, p = names(z))),
    as_tibble(w)
    )

microbenchmark(
    w[["p"]],
    w$p
)

f1 <- function() {
    a <- list(1:10)
    b <- list(letters)
    z <- data.frame(a = NA, b = NA)
    z$a <- a
    z$b <- b
    as_tibble(z)
}

f2 <- function() {
    a <- list(1:10)
    b <- list(letters)
    z <- tibble(a = a, b = b)
    z
}

microbenchmark(f1, f2)

### * tidy_flows()

m <- aquarium_mod
r <- aquarium_run
plot(r)

system.time({tf <- tidy_flows(m, r, 40, cores = 4)})

n <- c(20, 40)
t_ref <- c(8.11, 15.8)
t_new_1 <- c(5.2, 10.09)
t_new_4 <- c(3.06, 5.21)
t_newnew_4 <- c(1.46, 2.65)

t_newnew_4/t_ref

p <- profvis({
    tf <- tidy_flows(m, r, 100, cores = 1)
})

print(p)


### * predict()

m <- aquarium_mod
r <- aquarium_run
plot(r)

system.time({pr <- predict(m, r, draw = 200, cores = 4)})

n <- c(100, 200)
t_ref <- c(5.34, 10.24)
t_new_4 <- c(5.19, 9.15)
t_newnew_4 <- c(NA, 2.81)

p <- profvis({
    pr <- predict(m, r, draws = 200, cores = 1)
})

print(p)

### * tidy_mcmc_list()

p <- profvis({
    replicate(5, tidy_mcmc_list(aquarium_run, spread = FALSE))
})

p <- profvis({
    replicate(5, tidy_mcmc_list(aquarium_run, spread = TRUE))
})

print(p)

# New version

tidy_mcmc_list_new <- function(x, spread = FALSE) {
    z <- x
    paramNames <- colnames(z[[1]])
    if (!spread) {
        tables <- lapply(seq_along(z), function(i) {
            my_chain <- as.matrix(z[[i]])
            tibble::tibble(mcmc.chain = i,
                           mcmc.iteration = seq_len(nrow(my_chain)),
                           mcmc.parameters = lapply(seq_len(nrow(my_chain)),
                                                    function(k) {
                                                        setNames(my_chain[k,], paramNames)
                                                    })
                           )
        })
        out <- dplyr::bind_rows(tables)
        return(out)
    }
    mcmc_parameters <- lapply(seq_along(z), function(i) {
        my_chain <- tibble::as_tibble(as.matrix(z[[i]]))
        mcmc_pars <- tibble::tibble(mcmc.chain = i,
                                    mcmc.iteration = seq_len(nrow(my_chain)))
        dplyr::bind_cols(mcmc_pars, my_chain)
    })
    out <- dplyr::bind_rows(mcmc_parameters)
    return(out)
}

p <- profvis({
    replicate(5, tidy_mcmc_list_new(aquarium_run, spread = FALSE))
})

print(p)

z <- tidy_mcmc_list_new(aquarium_run)

microbenchmark(
    z <- tidy_mcmc_list(aquarium_run, spread = TRUE),
    z <- tidy_mcmc_list_new(aquarium_run, spread = TRUE),
    times = 5)
