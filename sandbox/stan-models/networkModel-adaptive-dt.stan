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

   This version uses an adaptive dt. Pre-calculated dt tables are provided in
   the data, with increasing number of steps, and the projection function will
   try the first time plan, test for relative changes in compartment sizes
   larger than a threshold, and if relative changes are too large switch to the
   next (higher resolution) time plan.

   If all time plans are tried without success, an error message is produced.

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

  /// * encodeTransitionMatrix_indices()

  /* Encode a transition matrix produced by buildTransitionMatrix().

     Return the indices values, i.e. integers that can be used to rebuild the
     matrix if provided with the cell values.

     @return int[] of length (2 + 2* n_comps * n_comps). The first integer is
     the number of comparments n_comps and the second the number of non-zero
     cells in the transition matrix. The following integers works by pairs, and
     each pair describe the (i,j) indices of a non-zero cell, ordered by
     increasing i and j, with j increasing faster (i.e. scanning every column j
     for a given row i, and then going to the next row).

  */

  int[] encodeTransitionMatrix_indices(int n_comps,
                               int n_src,
                               int[] src_indices,
                               int n_tr,
                               int[,] mappingTr,
                               int n_ls,
                               int[,] mappingLs) {
    int mask[n_comps, n_comps] = rep_array(0, n_comps, n_comps);
    int encoding[2 + 2 * n_comps * n_comps]; // Max size of encoding
    int encoding_index = 2; // First position is n_comps, second is the number of non-zero cells
    int n_nonzero_cells = 0; // Counter
    // Build a mask array to be able to detect non-zero values
    for (k in 1:n_tr) {
      mask[mappingTr[k, 2], mappingTr[k, 1]] = 1;
      mask[mappingTr[k, 1], mappingTr[k, 1]] = 1;
    }
    for (k in 1:n_ls) {
      mask[mappingLs[k, 1], mappingLs[k, 1]] = 1;
    }
    // Set all mask values to zero in source rows
    if (n_src > 0) {
      for (k in 1:n_src) {
        for (j in 1:n_comps) {
          mask[src_indices[k], j] = 0;
        }
      }
    }
    // From the mask, determine the encoding of the indices of non-zero cells
    encoding = rep_array(0, 2 + 2 * n_comps * n_comps);
    // Scan the mask
    for (i in 1:n_comps) {
      for (j in 1:n_comps) {
        if (mask[i,j] > 0) {
          n_nonzero_cells += 1;
          encoding_index += 1;
          encoding[encoding_index] = i;
          encoding_index += 1;
          encoding[encoding_index] = j;
        }
      }
    }
    encoding[1] = n_comps;
    encoding[2] = n_nonzero_cells;
    // Return
    return(encoding);
  }

  /// * encodeTransitionMatrix_values()

  /* Encode a transition matrix produced by buildTransitionMatrix().

     Return the cell contents, i.e. values that can be used to rebuild the
     matrix if provided with the cell indices.

     @param A matrix, output from buildTransitionMatrix()

     @return real[] of length (n_comps * n_comps). The values are ordered so
     that they can be used with the output from
     encodeTransitionMatrix_indices() to rebuild the transition matrix.
  */

  real[] encodeTransitionMatrix_values(int n_comps,
                               int n_src,
                               int[] src_indices,
                               int n_tr,
                               int[,] mappingTr,
                               int n_ls,
                               int[,] mappingLs, 
                               real[] params,
                               matrix A) {
    int mask[n_comps, n_comps] = rep_array(0, n_comps, n_comps);
    real encoding[n_comps * n_comps]; // Max size of encoding
    int encoding_index = 0;
    // Build a mask array to be able to detect non-zero values
    for (k in 1:n_tr) {
      mask[mappingTr[k, 2], mappingTr[k, 1]] = 1;
      mask[mappingTr[k, 1], mappingTr[k, 1]] = 1;
    }
    for (k in 1:n_ls) {
      mask[mappingLs[k, 1], mappingLs[k, 1]] = 1;
    }
    // Set all mask values to zero in source rows
    if (n_src > 0) {
      for (k in 1:n_src) {
        for (j in 1:n_comps) {
          mask[src_indices[k], j] = 0;
        }
      }
    }
    // From the mask, pick the values of non-zero cells in A
    encoding = rep_array(0.0, n_comps * n_comps);
    // Scan the mask
    for (i in 1:n_comps) {
      for (j in 1:n_comps) {
        if (mask[i,j] > 0) {
          encoding_index += 1;
          encoding[encoding_index] = A[i,j];
        }
      }
    }
    // Return
    return(encoding);
  }
  
  /// * ode_step_input_flows()

  /* An ODE definition that can be used to describe a system with constant
     input flows. For use with one of Stan's ODE solvers.

     @param theta Transition matrix parameters as encoded by
       encodeTransitionMatrix_values() (parameters)
     @param x_r empty (data)
     @param x_i Transition matrix indices as encoded by 
       encodeTransitionMatrix_indices() (data)
  */

  real[] ode_step_input_flows(real time,
                              real[] state,
                              real[] theta,
                              real[] x_r,
                              int[] x_i) {
    // Init
    int n_comps;
    int n_nonzero_cells;
    int indices_index;
    int ind_i;
    int ind_j;
    real dy_dt[x_i[1]];
    int theta_index;
    n_comps = x_i[1];
    n_nonzero_cells = x_i[2];
    // Calculate dy_dt
    dy_dt = rep_array(0, n_comps);
    indices_index = 1;
    theta_index = 0;
    for (i in 1:n_nonzero_cells) {
      indices_index += 2;
      theta_index += 1;
      ind_i = x_i[indices_index];
      ind_j = x_i[indices_index+1];
      dy_dt[ind_i] += theta[i] * state[ind_j];
    }
    // Return
    return(dy_dt);
  }

  /// * projectTrajectories()

  // @param n_proj_times Number of rows to read in proj_table
  // @param n_proj_steps Number of pieces in the step functions describing
  //   sources input profiles
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
                                int n_proj_steps,
                                int maxn_proj_steps,
                                real[] proj_table_time,
                                int[,] proj_table_steps,
                                int n_src,
                                int[] src_indices,
                                real[,] src_marked_values,
                                real[,] src_unmarked_values,
                                real[] theta,
                                real[] x_r,
                                int[] x_i,
                                int solver_selector,
                                real rel_tol, real abs_tol, int max_steps) {
    vector[n_comps] y[maxn_proj_times,2]; // col1 = marked, col2 = unmarked
    int marked = 1; // Index for marked column in y[,]
    int unmarked = 2; // Index for unmarked column in y[,]
    int projStartIndex = 1;
    int projEndIndex;
    real y_marked[n_comps]; // Container for a state
    real y_unmarked[n_comps]; // Container for a state
    int src_step = 1; // Counter for pieces of sources step functions
    real solutions[max(proj_table_steps[,2])+1, n_comps]; // Container for ODE solutions
    // Init
    y[1, marked] = marked_init;
    y[1, unmarked] = unmarked_init;
    // Update sources if needed
    if (proj_table_steps[1, 1] == 1) {
      if (n_src > 0) {
        for (m in 1:n_src) {
          y[1, marked][src_indices[m]] = src_marked_values[src_step, m];
          y[1, unmarked][src_indices[m]] = src_unmarked_values[src_step, m];
        }
        src_step += 1;
      }
    }
    // Init containers
    y_marked = to_array_1d(y[1, marked]);
    y_unmarked = to_array_1d(y[1, unmarked]);
    // Projection loop
    for (k in 1:n_proj_steps) {
      projEndIndex = projStartIndex + proj_table_steps[projStartIndex, 2];
      // Marked
      if (solver_selector == 1) {
        solutions[2:(proj_table_steps[projStartIndex, 2]+1), ] =
          integrate_ode_rk45(ode_step_input_flows,
                             y_marked,
                             proj_table_time[projStartIndex],
                             proj_table_time[(projStartIndex+1):projEndIndex],
                             theta, x_r, x_i, rel_tol, abs_tol, max_steps);
      } else {
        solutions[2:(proj_table_steps[projStartIndex, 2]+1), ] =
          integrate_ode_bdf(ode_step_input_flows,
                             y_marked,
                             proj_table_time[projStartIndex],
                             proj_table_time[(projStartIndex+1):projEndIndex],
                             theta, x_r, x_i, rel_tol, abs_tol, max_steps);
      }
      for (i in 1:proj_table_steps[projStartIndex, 2]) {
        y[projStartIndex+i, marked] = to_vector(solutions[i+1, ]);
      }
      // Unmarked
      if (solver_selector == 1) {
        solutions[2:(proj_table_steps[projStartIndex, 2]+1), ] =
          integrate_ode_rk45(ode_step_input_flows,
                             y_unmarked,
                             proj_table_time[projStartIndex],
                             proj_table_time[(projStartIndex+1):projEndIndex],
                             theta, x_r, x_i, rel_tol, abs_tol, max_steps);
      } else {
        solutions[2:(proj_table_steps[projStartIndex, 2]+1), ] =
          integrate_ode_bdf(ode_step_input_flows,
                            y_unmarked,
                            proj_table_time[projStartIndex],
                            proj_table_time[(projStartIndex+1):projEndIndex],
                            theta, x_r, x_i, rel_tol, abs_tol, max_steps);
      }
      for (i in 1:proj_table_steps[projStartIndex, 2]) {
        y[projStartIndex+i, unmarked] = to_vector(solutions[i+1, ]);
      }
      // Prepare next loop
      projStartIndex += proj_table_steps[projStartIndex, 2];
      if (proj_table_steps[projStartIndex, 1] == 1) {
        if (n_src > 0) {
          for (m in 1:n_src) {
            y[projStartIndex, marked][src_indices[m]] = src_marked_values[src_step, m];
            y[projStartIndex, unmarked][src_indices[m]] = src_unmarked_values[src_step, m];
          }
          src_step += 1;
        }
      }
      y_marked = to_array_1d(y[projStartIndex, marked]);
      y_unmarked = to_array_1d(y[projStartIndex, unmarked]);
    }
    // Return
    return(y);
  }

  /// * projectTrajectoriesAdaptive()

  // @return An array of vectors, each vector of length n_comps. First dimension is in sync
  //   with the time values of proj_table, second dimension is (marked, unmarked).

  
  vector[,] projectTrajectoriesAdaptive(int n_comps,
                                        int[] A_ij, // Encoding produced by encodeTransitionMatrix_indices()
                                        real[] A_values, // Encoding produced by encodeTransitionMatrix_values()
                                        vector marked_init,
                                        vector unmarked_init,
                                        int maxn_proj_times, // Number of projection times returned
                                        int maxn_proj_steps, // Number of discrete steps
                                        int [,,] proj_discrete_table_time_int, // For one group, time and dtNext
                                        real [,,] proj_discrete_table_time_real, // For one group, updateInitState and returnProj
                                        int n_dt_sets, // Number of dt sets available to test
                                        int n_src,
                                        int[] src_indices,
                                        real[,] src_marked_values,
                                        real[,] src_unmarked_values,
                                        real relThreshold) {
    vector[n_comps] out[maxn_proj_times,2]; // col1 = marked, col2 = unmarked
    vector[n_comps] y[maxn_proj_steps+1, 2];
    int marked = 1; // Index for marked column in out[,]
    int unmarked = 2; // Index for unmarked column in out[,]
    int src_step; // Counter for pieces of sources step functions
    int return_step; // Counter for projected times to be returned
    int A_n_tr = A_ij[2]; // Number of non-zero values in the A matrix
    int tr_index; // Counter for transition index
    // Temp container to avoid aliasing of variables
    vector[n_comps] tmp_y;
    // Checking relative change
    vector[n_comps] relChange;
    // First dt set is attempted
    int dtSet = 1;
    // Flag for largeChange
    int largeChange = 0;
    // Flag for success;
    int success = 0;
    // Time loop
    while ((!success) && (dtSet <= n_dt_sets)) {
      //print("Trying dt set number: ", dtSet);
      return_step = 1;
      src_step = 1;
      // Init returned container
      out[1, marked] = marked_init;
      out[1, unmarked] = unmarked_init;
      return_step += 1;
      // Init storage container for all time steps
      y[1, marked] = marked_init;
      y[1, unmarked] = unmarked_init;
      // Update sources if needed
      if (proj_discrete_table_time_int[1, 1, 1] == 1) {
        if (n_src > 0) {
          for (m in 1:n_src) {
            y[1, marked][src_indices[m]] = src_marked_values[src_step, m];
            y[1, unmarked][src_indices[m]] = src_unmarked_values[src_step, m];
          }
          src_step += 1;
        }
      }
      // Discrete steps
      for (t in 1:maxn_proj_steps) {
        for (marking in 1:2) {
          for (tr in 1:A_n_tr) {
            tmp_y = y[t, marking];
            y[t+1, marking] = tmp_y;
            tr_index = 1;
            // dt is proj_discrete_table_time_real[t,2,dtSet]
            y[t+1, marking][A_ij[tr_index]] += (y[t, marking][A_ij[tr_index+1]] *
                                                A_values[tr] *
                                                proj_discrete_table_time_real[t,2,dtSet]);
            tr_index += 2;
          }
          // Check for large change
          relChange = (y[t+1, marking] - y[t, marking]) ./ y[t, marking];
          if (fabs(max(relChange)) > relThreshold || fabs(min(relChange)) > relThreshold) {
            //print("fabs(max(relChange)): ", fabs(max(relChange)));
            //print("fabs(min(relChange)): ", fabs(min(relChange)));
            largeChange = 1;
          }
        }
        if (largeChange) {
          break;
        }
        // Update source if needed
        if (proj_discrete_table_time_int[t+1,1,dtSet] == 1) {
          if (n_src > 0) {
            for (m in 1:n_src) {
              y[t+1, marked][src_indices[m]] = src_marked_values[src_step, m];
              y[t+1, unmarked][src_indices[m]] = src_unmarked_values[src_step, m];
            }
            src_step += 1;
          }
        }
        // Store projected time if needed
        if (proj_discrete_table_time_int[t+1,2,dtSet] == 1) {
          out[return_step, unmarked] = y[t+1, unmarked];
          out[return_step, marked] = y[t+1, marked];
          return_step += 1;
        }
      }
      if (largeChange == 0) {
        success = 1;
      } else {
        largeChange = 0;
        dtSet += 1;
        //print("Too large relative change, moving to next dt set: ", dtSet);
      }
    }
    // Check for success
    if (!success) {
      // Raise an exception
      reject("Relative change too large given the attempted dt values! Try increasing the grid size for smaller dt values.");
    }
    //print("Success at step: ", dtSet);
    // Return
    return(out);
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
  int<lower=0> n_proj_steps[n_groups+1]; // Number of projection steps per group (padded)
  int<lower=0> maxn_proj_steps; // Maximum number of projection steps
  real<lower=0> proj_table[maxn_proj_times, 3, n_groups]; // Projection table
  real<lower=0> proj_table_time[maxn_proj_times, n_groups]; // Projection table: time data
  int<lower=0> proj_table_steps[maxn_proj_times, 2, n_groups]; // Projection table: steps data
  int<lower=0> n_step_points[n_groups+1]; // Number of sources step points per group (padded)
  int<lower=0> maxn_step_points; // Maximum number of step points
  real<lower=0> step_points_table[maxn_step_points, n_groups]; // Step points table

  /// * Projection table using a discrete, adaptive approach
  int<lower=0> n_proj_discrete_dt_values; // Number of dt sets available
  int<lower=0> maxn_proj_discrete_timesteps; // Maximum number of time steps given the dt sets
  int<lower=0> proj_discrete_table_time_int[maxn_proj_discrete_timesteps+1, 2, n_groups, n_proj_discrete_dt_values]; // updateInitState and returnProj booleans
  real<lower=0> proj_discrete_table_time_real[maxn_proj_discrete_timesteps+1, 2, n_groups, n_proj_discrete_dt_values]; // time and dtToNext
  real<lower=0> relThreshold; // Threshold above which the next dt set is used

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
  int encoding_indices[2 + 2 * n_comps_hidden * n_comps_hidden, n_groups];
  int<lower=0> n_total; // Total number of isotope proportion observations
  real x_r[0]; // Empty real data for ODE
  
  n_total = 0;

  /// * Encoding container
  for (g in 1:n_groups) {
    encoding_indices[,g] = encodeTransitionMatrix_indices(n_comps_hidden,
                                                          n_comps_src,
                                                          comps_src_indices,
                                                          n_params_transitions[g],
                                                          params_mapping_transitions[,,g],
                                                          n_params_losses[g],
                                                          params_mapping_losses[,,g]);
  }

  
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

  // Encoding container
  real theta[n_comps_hidden * n_comps_hidden];
  
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

    // Encode the parameters of the transition matrix
    theta = encodeTransitionMatrix_values(n_comps_hidden,
                                          n_comps_src,
                                          comps_src_indices,
                                          n_params_transitions[g],
                                          params_mapping_transitions[,,g],
                                          n_params_losses[g],
                                          params_mapping_losses[,,g],
                                          params,
                                          transitions);
    
    // Project the trajectories
    marked_init = buildInitQuantities(n_comps_hidden, g, 
                                      init_proportions[,g], sizes_mean[, g],
                                      buildMatrix_obsComp_to_hiddenComp, params, 1);
    unmarked_init = buildInitQuantities(n_comps_hidden, g, 
                                        init_proportions[,g], sizes_mean[, g],
                                        buildMatrix_obsComp_to_hiddenComp, params, 0);
    y[,,g] = projectTrajectoriesAdaptive(n_comps_hidden,
                                         encoding_indices[,g],
                                         theta,
                                         marked_init,
                                         unmarked_init,
                                         maxn_proj_times,
                                         maxn_proj_discrete_timesteps,
                                         proj_discrete_table_time_int[,,g,],
                                         proj_discrete_table_time_real[,,g,],
                                         n_proj_discrete_dt_values,
                                         n_comps_src,
                                         comps_src_indices,
                                         comps_src_marked_values[,,g],
                                         comps_src_unmarked_values[,,g],
                                         relThreshold);

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
