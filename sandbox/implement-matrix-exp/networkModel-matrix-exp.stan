/// * Functions

functions {

  /// * buildTransferMatrix() (CHECKED)

  /// ** Doc
  
  /*
  Parameters:
  
  nComps Integer, number of compartments in the network. Used to define the
      size of the transfer matrix, which is a square matrix nComps x nComps.

  nSteady Integer, number of steady-state compartments.

  steadyIndices Array of integers, indices of the steady-state compartments.

  nUpsilons Integer, number of upsilon parameters to add to the transfer
      matrix. This is the number of non-zero uptake rate parameters. An uptake
      rate parameters is non-zero if the two corresponding compartments are
      connected.

  mappingU Array of integers, with 3 columns and nUpsilons rows. The k^th row
      defines the mapping between the location of the k^th upsilon parameter in
      the transfer matrix at [mappingU[k, 2], mappingU[k, 1]] and the
      parameter values in params[mappingU[k, 3]]. It is used to fill in the
      transfer matrix with the upsilon rates values.

  nLambdas Integer, number of lambda parameters to add to the transfer
      matrix.

  mappingL Array of integers, with 2 columns and nLambdas rows. The k^th row
      defines the mapping between the location of the k^th lambda parameter in
      the loss rate vector at [mappingL[k, 1]] and the parameter values in
      params[mappingL[k, 2]].

  params Array of real numbers, of length the number of model parameters. It
      contains the values of the parameters currently sampled by the MCMC. They
      are actually transformed from the raw parameters defined in the
      parameters{} block, so they are already on the "real-world" scale and
      ready to use for network calculations.

  Returned value:

  A nComps by nComps transfer matrix A that can be used to project a network
  state from (t) to (t+dt). It is part of the equation: dx/dt = A x(t).
  */

  /// ** Code
  
  matrix buildTransferMatrix(int nComps, int nSteady, int[] steadyIndices,
                             int nUpsilons, int[,] mappingU,
                             int nLambdas, int[,] mappingL, real[] params) {
    matrix[nComps, nComps] transfer = rep_matrix(0, nComps, nComps);
    vector[nComps] lossRates = rep_vector(0, nComps);
    for (k in 1:nUpsilons) {
      transfer[mappingU[k, 2], mappingU[k, 1]] = params[mappingU[k, 3]];
      lossRates[mappingU[k, 1]] += params[mappingU[k, 3]];
    }
    for (k in 1:nLambdas) {
      lossRates[mappingL[k, 1]] += params[mappingL[k, 2]];
    }
    for (k in 1:nComps) {
      transfer[k,k] -= lossRates[k];
    }
    // Apply steady states
    if (nSteady > 0) {
      for (k in 1:nSteady) {
        transfer[steadyIndices[k],steadyIndices[k]] = 0;
      }
    }
    return(transfer);
  }

  /// * buildTransferMatrixDecay() (CHECKED)

  /*
  Parameters:

  nComps Integer, number of compartments in the network. Used to define the
      size of the transfer matrix, which is a square matrix nComps x nComps.

  transferRef Matrix of size nComps x nComps. This is the transfer matrix
      to which a decay will be applied on the diagonal.

  lambda_decay Real, the decay constant to be applied on the diagonal of the
      transfer matrix.

  Returned value:

  A nComps by nComps transfer matrix similar to transferRef but to which a
  decay has been applied on the diagonal.
  */

  matrix buildTransferMatrixDecay(int nComps, matrix transferRef,
                                  real lambda_decay) {
    matrix[nComps, nComps] transfer = transferRef;
    for (k in 1:nComps) {
      transfer[k,k] -= lambda_decay;
    }
    return(transfer);
  }

  /// * buildSizePredictions() (CHECKED)

  /*
  Parameters:

  nObs Integer, the number of observations for the current group.

  currentGroup Integer, the indice corresponding to the current group.

  maxNobs Integer, the maximum number of observations across all groups in the
      network model.

  unmarked Array of vectors, of length maxNtimesteps+1. Each vector is of
      length nComps and contains the quantities of unmarked material for the
      network compartment.

  marked Same as unmarked, but for marked material.

  indices Array of integer, defined by
      sizesObsIndices[maxNsizesObs,3,nGroups][,,g]. It has three columns and
      maxNsizesObs rows. For each observation (k^th row), indices[k,1] is the
      compartment index and indices[k,2] is the timestep index.
  
  splitPresent Boolean, is there any split compartment?

  splitComps Array of booleans, of length nComps. Indicates which compartments
      are split compartments.

  initRefr Array of reals, with two columns and nComps rows. For each
      compartment, it contains the initial quantity of unmarked (column 1) and
      marked (column 2) material in the refractory sub-compartment.
      
  Returned value:

  An array of reals of length maxNobs, containing the predicted sizes for all
  compartments and all timesteps.
  
  */
  
  real[] buildSizePredictions(int nObs, int currentGroup, int maxNobs,
                              vector[] unmarked, vector[] marked,
                              int[,] indices, int splitPresent,
                              int[] splitComps, real[,] initRefr) {
    real pred[maxNobs] = rep_array(0.0, maxNobs); 
    real unmarkedQ;
    real markedQ;
    if (splitPresent > 0) {
      // Split compartments
      for (k in 1:nObs) {
        unmarkedQ = unmarked[indices[k,2]][indices[k,1]];
        markedQ = marked[indices[k,2]][indices[k,1]];
        if (splitComps[indices[k,1]] > 0) {
          unmarkedQ = unmarkedQ + initRefr[indices[k,1], 1];
          markedQ = markedQ + initRefr[indices[k,1], 2];
        }
        pred[k] = unmarkedQ + markedQ;
      }
      return(pred);
    } else {
      // No split compartments
      for (k in 1:nObs) {
        unmarkedQ = unmarked[indices[k,2]][indices[k,1]];
        markedQ = marked[indices[k,2]][indices[k,1]];
        pred[k] = unmarkedQ + markedQ;
      }
      return(pred);
    }
  }

  /// * buildPropPredictions() (CHECKED)

  /*
  Parameters:

  nObs Integer, the number of observations for the current group.

  currentGroup Integer, the indice corresponding to the current group.

  maxNobs Integer, the maximum number of observations across all groups in the
      network model.

  unmarked Array of vectors, of length maxNtimesteps+1. Each vector is of
      length nComps and contains the quantities of unmarked material for the
      network compartment.

  marked Same as unmarked, but for marked material.

  indices Array of integer, defined by
      sizesObsIndices[maxNsizesObs,3,nGroups][,,g]. It has three columns and
      maxNsizesObs rows. For each observation (k^th row), indices[k,1] is the
      compartment index and indices[k,2] is the timestep index.
  
  splitPresent Boolean, is there any split compartment?

  splitComps Array of booleans, of length nComps. Indicates which compartments
      are split compartments.

  initRefr Array of reals, with two columns and nComps rows. For each
      compartment, it contains the initial quantity of unmarked (column 1) and
      marked (column 2) material in the refractory sub-compartment.
      
  Returned value:

  An array of reals of length maxNobs, containing the predicted proportions for
  all compartments and all timesteps.
  
  */
  
  real[] buildPropPredictions(int nObs, int currentGroup, int maxNobs,
                              vector[] unmarked, vector[] marked,
                              int[,] indices, int splitPresent,
                              int[] splitComps, real[,] initRefr) {
    real pred[maxNobs] = rep_array(0.0, maxNobs); 
    real unmarkedQ;
    real markedQ;
    if (splitPresent > 0) {
      // Split compartments
      for (k in 1:nObs) {
        unmarkedQ = unmarked[indices[k,2]][indices[k,1]];
        markedQ = marked[indices[k,2]][indices[k,1]];
        if (splitComps[indices[k,1]] > 0) {
          unmarkedQ = unmarkedQ + initRefr[indices[k,1], 1];
          markedQ = markedQ + initRefr[indices[k,1], 2];
        }
        pred[k] = markedQ / (unmarkedQ + markedQ);
      }
    return(pred);
    } else {
      // No split compartments
      for (k in 1:nObs) {
        unmarkedQ = unmarked[indices[k,2]][indices[k,1]];
        markedQ = marked[indices[k,2]][indices[k,1]];
        pred[k] = markedQ / (unmarkedQ + markedQ);
      }
      return(pred);
    }
  }

}

/// * Data (CHECKED)

data {

  /// * Counts (CHECKED)
  int<lower=1> nComps; // Number of compartments in a topology
  int<lower=0> nGroups; // Number of replication groups (i.e. rows in the original networkModel object)
  int<lower=0> nParams; // Total number of parameters to be estimated
  // Priors
  int<lower=1> nNonConstantPriors; // Number of non-constant parameters
  int<lower=0> nPriorUniform_code1; // Number of uniform priors
  int<lower=0> nPriorHcauchy_code2; // Number of half-Cauchy priors
  int<lower=0> nPriorBeta_code3; // Number of beta priors
  int<lower=0> nPriorTrNormal_code4; // Number of truncated normal priors

  /// * Parameter priors (CHECKED)

  // Constant
  real constantParams[nParams]; // Fixed values (for constant parameters)
  
  // Uniform
  real lowerParams[nParams]; // Lower boundaries (for uniform)
  real upperParams[nParams]; // Upper boundaries (for uniform)

  // Half-Cauchy
  real hcauchyScaleParams[nParams]; // Scale parameter (for half-Cauchy)

  // Scaled Beta
  real<lower=0> rawBetaAlpha[nParams]; // Alpha parameter (for scaled beta)
  real<lower=0> rawBetaBeta[nParams]; // Beta parameter (for scaled beta)
  real<lower=0> betaScaleParams[nParams]; // Scale parameter (for scaled beta)

  // Truncated normal
  real trNormMeanParams[nParams]; // Mean of the untruncated normal
  real trNormSdParams[nParams]; // Sd of the untruncated normal

  /// * Mapping between unscaled and scaled parameters (CHECKED)
  int<lower=0> mappingParamPriorType[nParams]; // Prior type for each parameter (numeric code)
  int<lower=0> mappingParamPriorID[nParams]; // Param ID within each prior typeID (1:nParams)

  /// * Distribution family for observed proportions (CHECKED)
  int<lower=1,upper=4> propFamily;
  // Code is:
  // 1 for gamma parameterized on mean and cv (eta = cv)
  // 2 for normal parameterized on mean and cv (eta = cv)
  // 3 for normal parameterized on mean and sd (eta = sd)
  // 4 for beta parameterized on mean and precision phi (eta = precision)
  
  /// * Distribution family for observed sizes (CHECKED)
  int<lower=1,upper=2> sizeFamily;
  // Code is:
  // 1 for normal parameterized on mean and cv (zeta = cv)
  // 2 for normal parameterized on mean and sd (zeta = sd)
  
  /// * Initial conditions (CHECKED)
  real<lower=0> initialQuantities[nComps,2,nGroups]; // Columns are (unmarked, marked)

  /// * Steady-state compartments (CHECKED)
  int<lower=0> maxNsteady; // Maximum number of steady state comps across groups
  int<lower=0> nSteady[nGroups+1]; // Padded
  int<lower=0> steadyIndices[maxNsteady,nGroups];

  /// * Split compartments (CHECKED)
  int<lower=0,upper=1> splitPresent; // Boolean
  int<lower=0,upper=1> splitComps[nComps,nGroups]; // Booleans

  /// * Parameter mapping for pis (portions of active compartments) (CHECKED)
  int<lower=0> piMapping[nComps, nGroups]; // Mapping to params[x]

  /// * Lambda due to decay (CHECKED)
  real<lower=0> lambda_decay; // 0 is equivalent to stable isotopes (no decay)
  
  /// * Time intervals (CHECKED)
  int<lower=0> maxNtimeIntervals; // Maximum number of intervals across groups
  int<lower=0> nTimeIntervals[nGroups+1]; // Number of intervals separated by events, per group (padded)
  real<lower=0> intervalsLengths[maxNtimeIntervals,nGroups]; // Duration of time intervals
  
  /// * Pulse events (CHECKED)
  int<lower=0> maxNpulseEvents; // Maximum number of pulse events across groups
  int<lower=0> nPulseEvents[nGroups+1]; // Padded
  int<lower=0> pulseEventsIndices[maxNpulseEvents,2,nGroups]; // Indices mapping events to the interval at the beginning of which they happen
  real pulseEventsQuantities[maxNpulseEvents,2,nGroups];
  
  /// * Observations (CHECKED)

  // Unique observation times per group (CHECKED)
  int<lower=0> maxNobsTimes; // Maximum number of unique obs times per group, across groups
  int<lower=0> nObsTimes[nGroups+1]; // Number of unique obs times per group, padded
  real<lower=0> elapsedTimeSinceEvent[nGroups,maxNobsTimes]; // Elapsed duration between observation times and the previous event (or t0)
  int<lower=0> obsIntervalsIndices[nGroups,maxNobsTimes]; // Indices of the intervals corresponding to each observation time
  
  // Individual observations (CHECKED)
  int<lower=0> maxNsizesObs;
  int<lower=0> maxNpropsObs;
  int<lower=0> nSizesObs[nGroups+1]; // Padded
  int<lower=0> nPropsObs[nGroups+1]; // Padded
  int<lower=0> sizesObsIndices[maxNsizesObs,3,nGroups]; // Columns are compartment, timepoint, zeta param index
  int<lower=0> propsObsIndices[maxNpropsObs,3,nGroups]; // Columns are compartment, timepoint, eta param index
  real<lower=0> sizesObs[maxNsizesObs, nGroups];
  real<lower=0> propsObs[maxNpropsObs, nGroups];
  
  /// * Mapping for upsilons (CHECKED)
  int<lower=0> maxNupsilons;
  int<lower=0> nUpsilons[nGroups+1]; // Number of upsilon parameters per group, padded
  int<lower=0> upsilonMapping[maxNupsilons,3,nGroups];
  
  /// * Mapping for lambdas (CHECKED)
  int<lower=0> maxNlambdas;
  int<lower=0> nLambdas[nGroups+1]; // Number of lambda parameters per group, padded
  int<lower=0> lambdaMapping[maxNlambdas,2,nGroups];
    
}

/// * Transformed data (CHECKED)

transformed data {

  /// * Initialization
  int<lower=0> nTotal; // Total number of observations
  nTotal = 0;
  
  /// * Count observations
  for (g in 1:nGroups) {
    nTotal += nSizesObs[g] + nPropsObs[g];
  }

}

/// * Parameters (CHECKED)

parameters {
  real<lower=0,upper=1> rawUniformParams[nPriorUniform_code1];
  real<lower=0> rawHcauchyParams[nPriorHcauchy_code2];
  real<lower=0,upper=1> rawBetaParams[nPriorBeta_code3];
  real<lower=0> rawTrNormParams[nPriorTrNormal_code4];
}

/// * Transformed parameters

transformed parameters {

  /// * Initialization (CHECKED)

  // Parameters on usable scale (converted from raw parameters)
  real<lower=0> params[nParams];

  // Initialize the transfer matrices (updated/reused from group to group)
  matrix[nComps,nComps] transfer;
  matrix[nComps,nComps] transferDecay;
  matrix[nComps,nComps] transition_tmp; // Temporary variable used to store exp(t*A) and reuse it several times

  // Initialize the arrays containing the initial states of each timeline
  // interval (updated/reused from group to group)
  vector[nComps] intervals_init_states_marked[maxNtimeIntervals];
  vector[nComps] intervals_init_states_unmarked[maxNtimeIntervals];
  
  // Initialize the [nGroups x nSteps] array for unmarked tracer
  vector[nComps] unmarked[nGroups, maxNobsTimes+1]; // nObsTimes + t0
  // Initialize the [nGroups x nSteps] array for marked tracer
  vector[nComps] marked[nGroups, maxNobsTimes+1]; // nObsTimes + t0

  // Initialize the array to store initial quantities for refractory portions
  real<lower=0> initRefr[nComps, 2, nGroups] = rep_array(0.0, nComps, 2, nGroups); // Columns are unmarked, marked
  
  // Predicted values (to compare with observations)
  real<lower=0> sizesPred[maxNsizesObs, nGroups];
  real<lower=0> propsPred[maxNpropsObs, nGroups];

  // Variables for likelihood calculations
  real<lower=0> sizesPred_zeta[maxNsizesObs, nGroups] = rep_array(0.0, maxNsizesObs, nGroups);
  real<lower=0> sizesPred_alpha[maxNsizesObs, nGroups] = rep_array(0.0, maxNsizesObs, nGroups);
  real<lower=0> sizesPred_beta[maxNsizesObs, nGroups] = rep_array(0.0, maxNsizesObs, nGroups);
  real<lower=0> propsPred_eta[maxNpropsObs, nGroups] = rep_array(0.0, maxNpropsObs, nGroups);
  real<lower=0> propsPred_alpha[maxNpropsObs, nGroups] = rep_array(0.0, maxNpropsObs, nGroups);
  real<lower=0> propsPred_beta[maxNpropsObs, nGroups] = rep_array(0.0, maxNpropsObs, nGroups);

  // Start block to use int counter
  // Cf. https://discourse.mc-stan.org/t/integer-loop-index-in-transformed-parameters-block/9264/3
  {

    int pulseIndex;
  
  /// * Convert raw parameters to usable parameters (CHECKED)
  for (i in 1:nParams) {
    if (mappingParamPriorType[i] == 0) {
      params[i] = constantParams[i];
    }
    if (mappingParamPriorType[i] == 1) {
      params[i] = lowerParams[i] + (upperParams[i] - lowerParams[i]) * rawUniformParams[mappingParamPriorID[i]];
    }
    if (mappingParamPriorType[i] == 2) {
      params[i] = hcauchyScaleParams[i] * rawHcauchyParams[mappingParamPriorID[i]];
    }
    if (mappingParamPriorType[i] == 3) {
      params[i] = betaScaleParams[i] * rawBetaParams[mappingParamPriorID[i]];
    }
    if (mappingParamPriorType[i] == 4) {
      params[i] = rawTrNormParams[mappingParamPriorID[i]];
    }
  }

  /// * Trajectory calculation per group
  for (g in 1:nGroups) {
    
    /// ** Build the transfer matrix A for group g (CHECKED)
    transfer = buildTransferMatrix(nComps, nSteady[g], steadyIndices[,g],
                             nUpsilons[g], upsilonMapping[,,g],
                             nLambdas[g], lambdaMapping[,,g], params);
    transferDecay = buildTransferMatrixDecay(nComps, transfer, lambda_decay);

    /// ** Initialize event index (CHECKED)
    pulseIndex = 1;
    
    /// ** Initialize first rows of unmarked and marked tracer quantities (CHECKED)
    // Those are the quantities at t=0
    unmarked[g, 1] = to_vector(initialQuantities[,1,g]);
    marked[g, 1] = to_vector(initialQuantities[,2,g]);

    /// ** Adjust split compartments (CHECKED)
    if (splitPresent > 0) {
      for (j in 1:nComps) {
        if (splitComps[j,g] > 0) {
          // Store the quantities for the refractory part
          initRefr[j,1,g] = unmarked[g, 1][j] * (1 - params[piMapping[j,g]]);
          initRefr[j,2,g] = marked[g, 1][j] * (1 - params[piMapping[j,g]]);
          // Update the initial conditions for the active part
          unmarked[g, 1][j] = unmarked[g, 1][j] * params[piMapping[j,g]];
          marked[g, 1][j] = marked[g, 1][j] * params[piMapping[j,g]];
        }
      }
    }

    /// ** Apply pulses, if any at t=0 (CHECKED)
    if (nPulseEvents[g] > 0) {
      if (pulseIndex <= nPulseEvents[g]) {
        while(pulseEventsIndices[pulseIndex,1,g] == 1) {
          unmarked[g, 1][pulseEventsIndices[pulseIndex,2,g]] += pulseEventsQuantities[pulseIndex,1,g];
          marked[g, 1][pulseEventsIndices[pulseIndex,2,g]] += pulseEventsQuantities[pulseIndex,2,g];
          pulseIndex += 1;
          if (pulseIndex > nPulseEvents[g]) {
            break;
          }
        }
      }
    }

    /// ** Store the initial state of the first timeline interval (CHECKED)
    intervals_init_states_unmarked[1] = unmarked[g, 1];
    intervals_init_states_marked[1] = marked[g, 1];

    // unmarked[g,1] and marked[g,1] will be overwritten by the first
    // observation time in the 1:nObsTimes[g] loop below.
    
    /// ** Calculate initial state for each timeline interval after the first one (CHECKED)
    if (nTimeIntervals[g] > 1) {
      for (t in 2:nTimeIntervals[g]) {

        // Calculate the end-point of the previous interval (CHECKED)
        transition_tmp = matrix_exp(intervalsLengths[t-1,g] * transferDecay);
        intervals_init_states_unmarked[t] = transition_tmp * intervals_init_states_unmarked[t-1];
        intervals_init_states_marked[t] = transition_tmp * intervals_init_states_marked[t-1];

        // Apply pulses, if any (CHECKED)
        if (nPulseEvents[g] > 0) {
          if (pulseIndex <= nPulseEvents[g]) {
            while(pulseEventsIndices[pulseIndex,1,g] == t) {
              intervals_init_states_unmarked[t][pulseEventsIndices[pulseIndex,2,g]] += pulseEventsQuantities[pulseIndex,1,g];
              intervals_init_states_marked[t][pulseEventsIndices[pulseIndex,2,g]] += pulseEventsQuantities[pulseIndex,2,g];
              pulseIndex += 1;
              if (pulseIndex > nPulseEvents[g]) {
                break;
              }
            }
          }
        } // End of pulse block
      } // End of loop over time intervals
    } // End of calculation of intervals initial states

    /// ** Calculate projections for each observation time (CHECKED)
    for (k in 1:nObsTimes[g]) {
      transition_tmp = matrix_exp(elapsedTimeSinceEvent[g,k]*transferDecay);
      unmarked[g,k] = transition_tmp*intervals_init_states_unmarked[obsIntervalsIndices[g,k]];
      marked[g,k] = transition_tmp*intervals_init_states_marked[obsIntervalsIndices[g,k]];
    }
      
    /// ** Store predicted values (CHECKED)
    sizesPred[,g] = buildSizePredictions(nSizesObs[g], g, maxNsizesObs,
                                         unmarked[g,], marked[g,],
                                         sizesObsIndices[,,g],
                                         splitPresent, splitComps[,g],
                                         initRefr[,,g]);
    propsPred[,g] = buildPropPredictions(nPropsObs[g], g, maxNpropsObs,
                                         unmarked[g,], marked[g,],
                                         propsObsIndices[,,g],
                                         splitPresent, splitComps[,g],
                                         initRefr[,,g]);

    /// ** Prepare variables for likelihood calculations (CHECKED)
    for (k in 1:nSizesObs[g]) {
      sizesPred_zeta[k,g] = params[sizesObsIndices[k,3,g]];
      if (sizeFamily == 1) {
        // Normal(mean, zeta = cv)
        sizesPred_alpha[k,g] = sizesPred[k,g]; // mean
        sizesPred_beta[k,g] = sizesPred_zeta[k,g] * sizesPred[k,g]; // sd
      }
      if (sizeFamily == 2) {
        // Normal(mean, zeta = sd)
        sizesPred_alpha[k,g] = sizesPred[k,g]; // mean
        sizesPred_beta[k,g] = sizesPred_zeta[k,g]; // sd
      }
    }
    for (k in 1:nPropsObs[g]) {
      propsPred_eta[k,g] = params[propsObsIndices[k,3,g]];
      if (propFamily == 1) {
        // Gamma(mean, eta = cv)
        propsPred_alpha[k,g] = pow(propsPred_eta[k,g], -2);
        propsPred_beta[k,g] = propsPred_alpha[k,g] / propsPred[k,g];
      }
      if (propFamily == 2) {
        // Normal(mean, eta = cv)
        propsPred_alpha[k,g] = propsPred[k,g]; // Mean
        propsPred_beta[k,g] = propsPred_eta[k,g] * propsPred_alpha[k,g]; // Sd
      }
      if (propFamily == 3) {
        // Normal(mean, eta = sd)
        propsPred_alpha[k,g] = propsPred[k,g]; // Mean
        propsPred_beta[k,g] = propsPred_eta[k,g]; // Sd
      }
      if (propFamily == 4) {
        // Beta(mean, eta = precision phi)
        propsPred_alpha[k,g] = propsPred[k,g] * propsPred_eta[k,g]; // alpha = mean * phi
        propsPred_beta[k,g] = propsPred_eta[k,g] * (1 - propsPred[k,g]); // beta = phi * (1 - mu)
      }
    }
    
  } // End of groups loop

  } // End of block for pulseIndex counter
  
} // End of transformed parameters block

/// * Model (OK)

model {

  /// * Priors (OK)
  for (i in 1:nParams) {
    if (mappingParamPriorType[i] == 1) {
      rawUniformParams[mappingParamPriorID[i]] ~ uniform(0, 1);
    }
    if (mappingParamPriorType[i] == 2) {
      rawHcauchyParams[mappingParamPriorID[i]] ~ cauchy(0, 1);
    }
    if (mappingParamPriorType[i] == 3) {
      rawBetaParams[mappingParamPriorID[i]] ~ beta(rawBetaAlpha[i], rawBetaBeta[i]);
    }
    if (mappingParamPriorType[i] == 4) {
      rawTrNormParams[mappingParamPriorID[i]] ~ normal(trNormMeanParams[i], trNormSdParams[i]);
    }
  }

  /// * Likelihood (OK)
  for (g in 1:nGroups) {

    /// ** Sizes
    if (sizeFamily == 1) {
      // Normal(mean, zeta = cv)
      sizesObs[1:nSizesObs[g], g] ~ normal(sizesPred_alpha[1:nSizesObs[g],g],
                                           sizesPred_beta[1:nSizesObs[g],g]);
    }
    if (sizeFamily == 2) {
      // Normal(mean, zeta = sd)
      sizesObs[1:nSizesObs[g], g] ~ normal(sizesPred_alpha[1:nSizesObs[g],g],
                                           sizesPred_beta[1:nSizesObs[g],g]);
    }

    /// ** Proportions
    if (propFamily == 1) {
      // Gamma(mean, eta = cv)
      propsObs[1:nPropsObs[g], g] ~ gamma(propsPred_alpha[1:nPropsObs[g],g],
                                          propsPred_beta[1:nPropsObs[g],g]);
    }
    if (propFamily == 2) {
      // Normal(mean, eta = cv)
      propsObs[1:nPropsObs[g], g] ~ normal(propsPred_alpha[1:nPropsObs[g],g],
                                           propsPred_beta[1:nPropsObs[g],g]);
    }
    if (propFamily == 3) {
      // Normal(mean, eta = sd)
      propsObs[1:nPropsObs[g], g] ~ normal(propsPred_alpha[1:nPropsObs[g],g],
                                           propsPred_beta[1:nPropsObs[g],g]);
    }
    if (propFamily == 4) {
      // Beta(mean, eta = precision phi)
      propsObs[1:nPropsObs[g], g] ~ beta(propsPred_alpha[1:nPropsObs[g],g],
                                         propsPred_beta[1:nPropsObs[g],g]);
    }
    
  } // End of groups loop
  
}

/// * Generated quantities (OK)

generated quantities {

  // Initialization
  vector[nNonConstantPriors] nonConstantParams;
  int paramIndex;
  vector[nTotal] log_lik;
  int llIndexShift;
  llIndexShift = 0;

  // Non-constant parameters
  paramIndex = 1;
  for (i in 1:nParams) {
    if (mappingParamPriorType[i] != 0) {
      // This is a non-constant parameter, save it
      nonConstantParams[paramIndex] = params[i];
      paramIndex += 1;
    }
  }
  
  // Sizes
  for (g in 1:nGroups) {
    for (o in 1:nSizesObs[g]) {
      if (sizeFamily == 1) {
        log_lik[o+llIndexShift] = normal_lpdf(sizesObs[o,g] | sizesPred_alpha[o,g], sizesPred_beta[o,g]);
      }
      if (sizeFamily == 2) {
        log_lik[o+llIndexShift] = normal_lpdf(sizesObs[o,g] | sizesPred_alpha[o,g], sizesPred_beta[o,g]);
      }
    }
    llIndexShift += nSizesObs[g];
  }

  // Proportions
  for (g in 1:nGroups) {
    for (o in 1:nPropsObs[g]) {
      if (propFamily == 1) {
        log_lik[o+llIndexShift] = gamma_lpdf(propsObs[o,g] | propsPred_alpha[o,g], propsPred_beta[o,g]);
      }
      if (propFamily == 2) {
        log_lik[o+llIndexShift] = normal_lpdf(propsObs[o,g] | propsPred_alpha[o,g], propsPred_beta[o,g]);
      }
      if (propFamily == 3) {
        log_lik[o+llIndexShift] = normal_lpdf(propsObs[o,g] | propsPred_alpha[o,g], propsPred_beta[o,g]);
      }
      if (propFamily == 4) {
        log_lik[o+llIndexShift] = beta_lpdf(propsObs[o,g] | propsPred_alpha[o,g], propsPred_beta[o,g]);
      }
    }
    llIndexShift += nPropsObs[g];
  }
  
}
