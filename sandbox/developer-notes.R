# [2021-11-18 Thu]

# Working on the code for projections.

devtools::load_all()

### * Reference implementation before deSolve

w <- aquarium_mod
w <- set_params(w, sample_params(w)) %>%
    project(end = 10)

### * project_desolve()

z <- aquarium_mod

w <- project_desolve(nm = z, at = seq(0, 5, length.out = 128))

project_desolve <- function(nm, at = NULL) {
    if (is.null(at)) { at <- get_project_default_at(nm) }
    z <- prep_stan_data_expm(nm)
    for (g in seq_len(z[["nGroups"]])) {
        # Build the transfer matrix
        transfer <- build_transfer_matrix(topo = topo(nm[g,]),
                                          params = nm$parameters[[g]],
                                          dt = 1)
        transfer <- transfer - diag(z[["lambda_decay"]], ncol(transfer))
        # Build the model equation for deSolve integration
        model_equation <- function(t, state, ...) {
            list(transfer %*% state)
        }
        # Project in order to get the starting conditions of all intervals
        n_intervals <- z[["nTimeIntervals"]][g]
        intervals_init <- list()
        WIP
    }
}

### * get_project_default_at()

#' Get default projection time points for a network model
#'
#' @param nm A networkModel object.
#' @param n Length of the returned sequence of time points.
#'
#' @return A vector containing time points.
#'
#' @examples
#' get_project_default_at(aquarium_mod)
#' get_project_default_at(trini_mod)
#'
#' @keywords internal
#' @noRd

get_project_default_at <- function(nm, n = 128) {
    end_obs <- max(sapply(nm[["observations"]], function(x) {
        max(x[["time"]], na.rm = TRUE)
    }))
    if ("events" %in% colnames(nm)) {
        end_events <- sapply(nm[["events"]], function(x) {
            max(x[["time"]], na.rm = TRUE)
        })
    } else {
        end_events <- NA
    }
    end_time <- max(end_obs, end_events, na.rm = TRUE)
    stopifnot(!is.na(end_time))
    return(seq(0, end_time, length.out = 128))
}
