/// * Description

/* Stan model for isotracer

   - Works for any number of replication groups (nGroups >= 1)

   - Works for for split and non-split compartments

   It assumes all replication groups share the same (raw and observed)
   compartments, in the same order in the topology matrix (even if the connections
   between compartments can be different, i.e. the uptake mask content can
   change).

   It also assumes that the same compartments are input sources in all replicates,
   and only works with input_flow for now (but the input_flow profiles can differ
   across replicates).

   For varying bounds on a parameter vector, see:
   https://mc-stan.org/docs/2_18/stan-users-guide/vectors-with-varying-bounds.html

   In the comments below, "mapping" without further characterization is mostly to
   be understood as "mapping to one of all the estimated parameters".

   "Padded" refers to adding an extra element to an array, to ensure that it
   behaves as an array even if only one meaningful element is provided by the
   data.

   This version of the model uses the matrix exponential to solve the system of 
   differential equations governing the trajectories.

   For each step over which the input flows are constant, the trajectories for
   marker quantities are defined by:

   dy/dt = A * Y(t)

   with y(t) a vector of n_hidden_comps values, and A the transition matrix (without
   ones on the diagonal). If y0 is the condition at the left boundary, then y(t) can 
   be calculated for any time t with t>t0 as:

   y(t) = expm(A * (t-t0) ) * y0

*/

/// * Functions

functions {
  
  /// * buildTransitionMatrix()
  matrix buildTransitionMatrix(int n_comps,
                               int n_src,
                               int[] src_indices,
                               int n_tr,
                               int[,] mappingTr,
                               int n_ls,
                               int[,] mappingLs, 
                               real[] params) {
    matrix[n_comps, n_comps] transition = rep_matrix(0, n_comps, n_comps);
    vector[n_comps] lossRates = rep_vector(0, n_comps);
    for (k in 1:n_tr) {
      transition[mappingTr[k, 2], mappingTr[k, 1]] = params[mappingTr[k, 3]];
      lossRates[mappingTr[k, 1]] += params[mappingTr[k, 3]];
    }
    for (k in 1:n_ls) {
      lossRates[mappingLs[k, 1]] += params[mappingLs[k, 2]];
    }
    for (k in 1:n_comps) {
      transition[k,k] -= lossRates[k];
    }
    // Set row values to 0 for sources (constant, so their derivatives is zero)
    if (n_src > 0) {
      for (k in 1:n_src) {
        for (j in 1:n_comps) {
          transition[src_indices[k], j] = 0;
        }
      }
    }
    return(transition);
  }

  /// * projectTrajectories()

  // @param n_proj_times Number of rows to read in proj_table
  // @param proj_table First column is time value, second column is 1/0 for input step or not
  //
  // @return An array of vectors, each vector of length n_comps. First dimension is in sync
  //   with the time values of proj_table, second dimension is (marked, unmarked).
  
  vector[,] projectTrajectories(int n_comps,
                                matrix A,
                                vector marked_init,
                                vector unmarked_init,
                                int n_proj_times,
                                int maxn_proj_times,
                                real[,] proj_table,
                                int n_src,
                                int[] src_indices,
                                real[,] src_marked_values,
                                real[,] src_unmarked_values) {
    vector[n_comps] y[maxn_proj_times,2]; // col1 = marked, col2 = unmarked
    int marked = 1;
    int unmarked = 2;
    vector[n_comps] y_marked;
    vector[n_comps] y_unmarked;
    matrix[n_comps,n_comps] expAt;
    int src_step = 1;
    // Init
    y[1, marked] = marked_init;
    y[1, unmarked] = unmarked_init;
    // Update sources if needed
    if (proj_table[1, 2] == 1) {
      if (n_src > 0) {
        for (m in 1:n_src) {
          y[1, marked][src_indices[m]] = src_marked_values[src_step, m];
          y[1, unmarked][src_indices[m]] = src_unmarked_values[src_step, m];
        }
        src_step += 1;
      }
    }
    // Projection loop
    for (k in 2:n_proj_times) {
      expAt = matrix_exp((proj_table[k, 1] - proj_table[k-1, 1]) * A);
      y_marked = y[k-1, marked];
      y_unmarked = y[k-1, unmarked];
      y[k, marked] = expAt * y_marked;
      y[k, unmarked] = expAt * y_unmarked;
      // Update sources if needed
      if (proj_table[k, 2] == 1) {
        if (n_src > 0) {
          for (m in 1:n_src) {
            y[k, marked][src_indices[m]] = src_marked_values[src_step, m];
            y[k, unmarked][src_indices[m]] = src_unmarked_values[src_step, m];
          }
          src_step += 1;
        }
      }
    }
    // Return
    return(y);
  }
  
  /// * buildInitQuantities()
  vector buildInitQuantities(int nRawComps, int currentGroup, 
                             real[] initProps, real[] initBiomass,
                             int[,] obs2raw,
                             real[] params, int returnMarked) {
    vector[nRawComps] myQuantities;
    myQuantities = rep_vector(0, nRawComps);
    for (k in 1:nRawComps) {
      if (returnMarked) {
        myQuantities[k] = initBiomass[obs2raw[k,2]] * initProps[obs2raw[k,2]]; // No split (yet)
      } else {
        myQuantities[k] = initBiomass[obs2raw[k,2]] * (1-initProps[obs2raw[k,2]]); // No split (yet)
      }
      if (obs2raw[k,3] == 2) {
        myQuantities[k] *= params[obs2raw[k,3+currentGroup]]; // Active portion
      } else if (obs2raw[k,3] == 3) {
        myQuantities[k] *= (1-params[obs2raw[k,3+currentGroup]]); // Refractory portion
      }
    }
    return(myQuantities);
  }
  
  /// * buildPredictions()
  real[,] buildPredictions(int nObs, int currentGroup, int maxnObs,
                           vector[] unmarked, vector[] marked,
                           int[,] obs_x, int[,] raw2obs) {
    real pred[maxnObs,2]; // Cols: (biomass, prop)
    real unmarkedBM;
    real markedBM;
    for (k in 1:nObs) {
      if (raw2obs[obs_x[k,1],3] == 0) { // Unsplit compartment
        unmarkedBM = unmarked[obs_x[k,3]][raw2obs[obs_x[k,1],2]];
        markedBM = marked[obs_x[k,3]][raw2obs[obs_x[k,1],2]];
      } else { // Split compartment
        unmarkedBM = (unmarked[obs_x[k,3]][raw2obs[obs_x[k,1],2]] +
                      unmarked[obs_x[k,3]][raw2obs[obs_x[k,1],3]]);
        markedBM = (marked[obs_x[k,3]][raw2obs[obs_x[k,1],2]] +
                    marked[obs_x[k,3]][raw2obs[obs_x[k,1],3]]);
      }
      pred[k,1] = unmarkedBM + markedBM;
      pred[k,2] = markedBM / (unmarkedBM + markedBM);
    }
    return(pred);
  }
  
}

/// * Data

data {
  /// * Counts
  int<lower=0> n_groups; // Number of replication groups
  int<lower=1> n_comps_hidden; // Number of hidden compartments (including .act and .refr)
  int<lower=1> n_comps_obs; // Number of observed compartments (.act and .refr merged)
  int<lower=0> n_comps_src; // Number of input source compartments
  int<lower=0> n_params_all; // Total number of parameters to be estimated
  int<lower=0> n_priors_uniform_numcode1; // Number of uniform priors
  int<lower=0> n_priors_hcauchy_numcode2; // Number of half-Cauchy priors
  
  /// * Prior parameters
  real params_prior_uniform_lower[n_params_all]; // Lower boundaries (for uniform)
  real params_prior_uniform_upper[n_params_all]; // Upper boundaries (for uniform)
  real params_prior_hcauchy_scale[n_params_all]; // Scale parameter (for half-Cauchy)

  /// * Mapping between raw params and parameters
  int<lower=0> params_mapping_priorType[n_params_all]; // Prior type for each parameter
  int<lower=0> params_mapping_paramID_to_priorID[n_params_all]; // Mapping from parameters to prior ID within prior types
  
  /// * Recipe to build observed compartments from hidden compartments
  int<lower=0> buildMatrix_hiddenComp_to_obsComp[n_comps_obs, 3+n_groups];
  
  /// * Recipe to build hidden compartments from observed compartments
  int<lower=0> buildMatrix_obsComp_to_hiddenComp[n_comps_hidden, 3+n_groups];
  
  /// * Projection table and step points
  int<lower=0> n_proj_times[n_groups+1]; // Number of projection times per group (padded)
  int<lower=0> maxn_proj_times; // Maximum number of projection times
  real<lower=0> proj_table[maxn_proj_times, 3, n_groups]; // Projection table
  int<lower=0> n_step_points[n_groups+1]; // Number of sources step points per group (padded)
  int<lower=0> maxn_step_points; // Maximum number of step points
  real<lower=0> step_points_table[maxn_step_points, n_groups]; // Step points table
  
  /// * Input source compartments (values at step points)
  int<lower=0> comps_src_indices[n_comps_src+1]; // Indices of input compartments (padded)
  real<lower=0> comps_src_marked_values[maxn_step_points,n_comps_src,n_groups]; // Pre-calculated marked input
  real<lower=0> comps_src_unmarked_values[maxn_step_points,n_comps_src,n_groups]; // Pre-calculated unmarked input
  
  /// * Transition parameters
  int<lower=0> n_params_transitions[n_groups+1]; // Number of transition coefficients per group (padded)
  int<lower=0> maxn_params_transitions; // Maximum number of transition coefficients across groups
  int<lower=0> params_mapping_transitions[maxn_params_transitions, 3, n_groups]; // Cols: (from, to, paramID)

  /// * Loss parameters
  int<lower=0> n_params_losses[n_groups+1]; // Number of loss coefficients per group (padded)
  int<lower=0> maxn_params_losses; // Maximum number of loss coefficients across groups
  int<lower=0> params_mapping_losses[maxn_params_losses, 2, n_groups]; // Cols: (from, paramID)

  /// * Eta parameters
  int<lower=0> params_mapping_eta[n_groups+1]; // Indices for eta parameter for each group (padded)
  
  /// * Starting proportions
  real<lower=0> init_proportions[n_comps_obs, n_groups]; 

  /// * Biomasses (mean and sd)
  real<lower=0> sizes_mean[n_comps_obs, n_groups];
  real<lower=0> sizes_sd[n_comps_obs, n_groups];
  
  /// * Isotope proportion observations
  int<lower=0> n_observations[n_groups+1]; // Number of observations per group (padded)
  int<lower=0> maxn_observations; // Maximum number of observations across groups
  int<lower=0> props_obs_indep[maxn_observations, 3, n_groups]; // Cols: (obsComp, time, timeIndex)
  real<lower=0> props_obs_dep[maxn_observations, n_groups]; // Obs prop values
}

/// * Transformed data

transformed data {

  /// * Initialization
  real obs_sizes_mean[maxn_observations, n_groups];
  real obs_sizes_sd[maxn_observations, n_groups];
  int<lower=0> n_total; // Total number of isotope proportion observations
  n_total = 0;
  
  /// * Observed mean and sd biomass
  for (g in 1:n_groups) {
    for (k in 1:n_observations[g]) {
      obs_sizes_mean[k,g] = sizes_mean[props_obs_indep[k,1,g], g];
      obs_sizes_sd[k,g] = sizes_sd[props_obs_indep[k,1,g], g];
    }
  }

  /// * Count observations
  for (g in 1:n_groups) {
    n_total += n_observations[g];
  }
  
}

/// * Parameters

parameters {
  real<lower=0,upper=1> params_uniform_raw[n_priors_uniform_numcode1];
  real<lower=0> params_hcauchy_raw[n_priors_hcauchy_numcode2];
}

/// * Transformed parameters

transformed parameters {
  
  /// * Initialization

  // Parameters on usable scale (converted from raw parameters)
  real<lower=0> params[n_params_all];
  
  // Initialize the transition matrix (updated for each group)
  matrix[n_comps_hidden,n_comps_hidden] transitions;
  
  // Initialize the array for projected trajectories
  vector[n_comps_hidden] marked_init;
  vector[n_comps_hidden] unmarked_init;
  vector[n_comps_hidden] y[maxn_proj_times, 2, n_groups]; // Second dimension is for (marked, unmarked)

  // Initialize the arrays for predicted sizes and predicted proportions
  real pred_both[maxn_observations, 2, n_groups];
  real pred_sizes[maxn_observations, n_groups];
  real pred_sizes_sd[maxn_observations, n_groups];
  real pred_props[maxn_observations, n_groups];
  real pred_eta[maxn_observations, n_groups];
  real pred_alpha[maxn_observations, n_groups];
  real pred_beta[maxn_observations, n_groups];

  /// * Convert raw parameters to usable parameters
  for (i in 1:n_params_all) {
    if (params_mapping_priorType[i] == 1) {
      // Uniform
      params[i] = params_prior_uniform_lower[i] + (params_prior_uniform_upper[i] - params_prior_uniform_lower[i]) * params_uniform_raw[params_mapping_paramID_to_priorID[i]];
    }
    if (params_mapping_priorType[i] == 2) {
      // Half-Cauchy
      params[i] = params_prior_hcauchy_scale[i] * params_hcauchy_raw[params_mapping_paramID_to_priorID[i]];
    }
  }
  
  /// * Trajectory calculation per group
  
  for (g in 1:n_groups) {

    // Build the transition matrix
    transitions = buildTransitionMatrix(n_comps_hidden,
                                        n_comps_src,
                                        comps_src_indices,
                                        n_params_transitions[g],
                                        params_mapping_transitions[,,g],
                                        n_params_losses[g],
                                        params_mapping_losses[,,g],
                                        params);

    // Project the trajectories
    marked_init = buildInitQuantities(n_comps_hidden, g, 
                                      init_proportions[,g], sizes_mean[, g],
                                      buildMatrix_obsComp_to_hiddenComp, params, 1);
    unmarked_init = buildInitQuantities(n_comps_hidden, g, 
                                        init_proportions[,g], sizes_mean[, g],
                                        buildMatrix_obsComp_to_hiddenComp, params, 0);
    y[,,g] = projectTrajectories(n_comps_hidden,
                                 transitions,
                                 marked_init,
                                 unmarked_init,
                                 n_proj_times[g],
                                 maxn_proj_times,
                                 proj_table[,,g],
                                 n_comps_src,
                                 comps_src_indices,
                                 comps_src_marked_values[,,g],
                                 comps_src_unmarked_values[,,g]);

    // Predicted values and eta for each observation
    pred_both[,,g] = buildPredictions(n_observations[g], g, maxn_observations,
                                      y[,2,g], y[,1,g],
                                      props_obs_indep[,,g], buildMatrix_hiddenComp_to_obsComp);
    pred_sizes[,g] = pred_both[,1,g];
    pred_props[,g] = pred_both[,2,g];

    // Assign biomass sd (this can be moved outside the transformed parameters)
    for (k in 1:n_observations[g]) {
      pred_sizes_sd[k,g] = sizes_sd[props_obs_indep[k,1,g], g];
    }

    // Assign eta parameter and calculate alpha and beta
    for (k in 1:n_observations[g]) {
      pred_eta[k,g] = params[params_mapping_eta[g]];
      pred_alpha[k,g] = pow(pred_eta[k,g], -2);
      pred_beta[k,g] = pred_alpha[k,g] / pred_props[k,g];
    }

  } // End of groups loop

}

/// * Model

model {

  /// * Priors
  if (n_priors_uniform_numcode1 > 0) {
    for (i in 1:n_priors_uniform_numcode1) {
      params_uniform_raw[i] ~ uniform(0, 1);
    }
  }
  if (n_priors_hcauchy_numcode2 > 0) {
    for (i in 1:n_priors_hcauchy_numcode2) {
      params_hcauchy_raw[i] ~ cauchy(0, 1);
    }
  }

  /// * Likelihood for observed biomasses and proportions
  for (g in 1:n_groups) {
    obs_sizes_mean[1:n_observations[g],g] ~ normal(pred_sizes[1:n_observations[g],g],
                                                   pred_sizes_sd[1:n_observations[g],g]);
    props_obs_dep[1:n_observations[g],g] ~ gamma(pred_alpha[1:n_observations[g],g],
                                                 pred_beta[1:n_observations[g],g]);
  }
  
}

/// * Generated quantities

generated quantities {

  // Initialization
  vector[2*n_total] log_lik;
  int llIndexShift;
  //real<lower=0> ppc_props[maxn_observations, n_groups];
  //real<lower=0> ppc_sizes[maxn_observations, n_groups];
  llIndexShift = 0;

  // Biomasses
  for (g in 1:n_groups) {
    for (o in 1:n_observations[g]) {
      log_lik[o+llIndexShift] = normal_lpdf(obs_sizes_mean[o,g] | pred_sizes[o,g], pred_sizes_sd[o,g]);
    }
    llIndexShift += n_observations[g];
  }

  // Proportions
  for (g in 1:n_groups) {
    for (o in 1:n_observations[g]) {
      log_lik[o+llIndexShift] = gamma_lpdf(props_obs_dep[o,g] | pred_alpha[o,g], pred_beta[o,g]);
    }
    llIndexShift += n_observations[g];
  }

  // For ppc on truncated distributions:
  // https://discourse.mc-stan.org/t/rng-for-truncated-distributions/3122/12
  
  // // Posterior predictive check
  // for (g in 1:n_groups) {
  //   ppc_props[1:n_observations[g],g] = gamma_rng(pred_alpha[1:n_observations[g],g],
  //                                                pred_beta[1:n_observations[g],g]);
  //   ppc_sizes[1:n_observations[g],g] = normal_rng(pred_sizes[1:n_observations[g],g],
  //                                                 pred_sizes_sd[1:n_observations[g],g]);
  // }
  
}
