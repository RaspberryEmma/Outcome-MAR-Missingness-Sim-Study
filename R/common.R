# ****************************************
# Missingness Simulation Study
#
# Common functions used in multiple places
#
# Emma Tarmey
#
# Started:          06/10/2025
# Most Recent Edit: 19/01/2026
# ****************************************



# ----- Helper functions -----

# logical inverse of %in% operator
'%!in%' <- function(x,y)!('%in%'(x,y))

# Helper function for recording coefficients fitted
# Accounts for the idea that some modelling methods will exclude variables
# i.e:  Variables correct order with NaNs filling in excluded values
fill_in_blanks <- function(coefs = NULL, labels = NULL) {
  for (label in labels) {
    # if a given variable doesn't exist, create as NaN
    if (!(label %in% names(coefs))) {
      coefs[label] <- NaN
    }
  }
  
  # assert variable ordering
  coefs <- coefs[order(factor(names(coefs), levels = labels))]
  
  return (coefs)
}

# Helper function for recording covariate selection
# i.e:  Binary indicators per variable, 1 = selected, 0 = not
fill_in_cov_selection <- function(coefs = NULL, labels = NULL) {
  for (label in labels) {
    if (!(label %in% names(coefs))) {
      coefs[label] <- 0
    } else {
      coefs[label] <- 1
    }
  }
  
  # assert variable ordering
  coefs <- coefs[order(factor(names(coefs), levels = labels))]
  
  return (coefs)
}

# Create a string corresponding to a regression of Y on a given set of covariates
make_model_formula <- function(vars_selected = NULL) {
  
  if (length(vars_selected) > 0) {
    # first item does not require a '+' sign
    formula_string <- paste("Y ~ ", vars_selected[1], sep = "")
    
    for (var in vars_selected[-1]) {
      formula_string <- paste(formula_string, " + ", var, sep = "")
    }
  }
  else {
    # if no variables are selected, use constant term only
    formula_string <- "Y ~ 0"
  }
  
  return(formula_string)
}


# Create a string corresponding to a regression of X on a given set of covariates
make_X_model_formula <- function(vars_selected = NULL) {
  # remove X
  if ('X' %in% vars_selected) {
    vars_selected <- vars_selected[! vars_selected%in% c('X')]
  }
  
  if (length(vars_selected) > 0) {
    # first item does not require a '+' sign
    formula_string <- paste("X ~ ", vars_selected[1], sep = "")
    
    for (var in vars_selected[-1]) {
      formula_string <- paste(formula_string, " + ", var, sep = "")
    }
  }
  else {
    # if no variables are selected, use constant term only
    formula_string <- "X ~ 0"
  }
  
  return(formula_string)
}

# Inverse-logit transform function
# Used for simulating binary data
inverse_logit <- function(real_values = NULL) {
  probabilities <- (1)/(1 + exp(-1 * real_values))
  return (probabilities)
}


# Estimate the variance in Y before error
# Idea is that we fix this for given values of R2 and b
# We then get a compatible error term later on
determine_subgroup_var_Y <- function(num_total_conf = NULL,
                                     beta_X         = NULL,
                                     causal         = NULL,
                                     Z_correlation  = NULL,
                                     target_r_sq_Y  = NULL) {
  
  subgroup_size <- num_total_conf / 4
  beta_X_high   <- beta_X
  beta_X_low    <- beta_X / 4
  
  A <- (causal * beta_X_low)  + beta_X_high
  B <- (causal * beta_X_high) + beta_X_high
  C <- (causal * beta_X_high) + beta_X_low
  D <- (causal * beta_X_low)  + beta_X_low
  
  pairwise_combinations <- (A*B) + (A*C) + (A*D) + (B*C) + (B*D) + (C*D)
  
  LHS <- (A^2 + B^2 + C^2 + D^2) * (subgroup_size + (Z_correlation * subgroup_size * (subgroup_size - 1)))
  MID <- 2 * Z_correlation * subgroup_size^2 * pairwise_combinations
  RHS <- causal^2
  
  value <- LHS + MID + RHS
  
  
  return (value)
}


# Estimate the variance in the error term for Y
# Takes a given variance of Y and inverts the R2 formula
determine_subgroup_var_error_Y <- function(var_Y          = NULL,
                                           target_r_sq_Y  = NULL) {
  
  LHS   <- 1 - target_r_sq_Y
  RHS   <- var_Y / target_r_sq_Y
  value <- sqrt(LHS * RHS)
  
  message("TODO: remove MAGIC NUMBER (Currently R2X = 0.3, R2Y = 0.3)")
  value <- 2.94
  
  return (value)
}


# The value for all beta coefficients used for generating X
beta_X_formula <- function(num_total_conf = NULL,
                           target_r_sq_X  = NULL,
                           Z_correlation  = NULL) {
  
  numerator   <- target_r_sq_X
  denominator <- num_total_conf * (1 - target_r_sq_X) * (1 + Z_correlation*(num_total_conf - 1))
  
  value <- sqrt(numerator / denominator)
  
  message("TODO: remove MAGIC NUMBER (Currently R2X = 0.3, R2Y = 0.3)")
  value <- 0.092
  return (value)
}


# The value for the beta coefficients used for generating X
# 4 different values corresponding to the 4 subgroups
beta_X_subgroups_formula <- function(beta_X = NULL) {
  beta_X_1 <- beta_X / 4 # low X, low Y
  beta_X_2 <- beta_X / 4 # low X, high Y
  beta_X_3 <- beta_X     # high X, low Y
  beta_X_4 <- beta_X     # high X, high Y
  
  beta_Xs <- c(beta_X_1, beta_X_2, beta_X_3, beta_X_4)
  
  return (beta_Xs)
}


# The value for the beta coefficients used for generating Y
# 4 different values corresponding to the 4 subgroups
beta_Y_subgroups_formula <- function(beta_X = NULL) {
  beta_Y_1 <- beta_X / 4 # low X, low Y
  beta_Y_2 <- beta_X     # low X, high Y
  beta_Y_3 <- beta_X / 4 # high X, low Y
  beta_Y_4 <- beta_X     # high X, high Y
  
  beta_Ys <- c(beta_Y_1, beta_Y_2, beta_Y_3, beta_Y_4)
  
  return (beta_Ys)
}

estimate_within_CI <- function(estimate       = NULL,
                               true_value     = NULL,
                               standard_error = NULL,
                               sample_size    = NULL) {
  within_CI   <- 0.0
  
  # inverse t score
  # 2 sided test for 95% confidence interval
  t_score <- qt(p = 0.975, df = sample_size - 2)
  
  # upper and lower bounds of the confidence interval
  upper_bound <- estimate + (t_score * standard_error)
  lower_bound <- estimate - (t_score * standard_error)
  
  if ((true_value <= upper_bound) && (true_value >= lower_bound)) {
    within_CI <- 1.0
  }
  
  return (within_CI)
}



# ----- Missingness mechanisms and missingness handling -----

apply_MCAR_missingness <- function(data = NULL, vars_to_censor = NULL) {
  # MCAR probability of censorship does not depend on data value
  psel <- 0.25
  
  for (var in vars_to_censor) {
    # MCAR selection vector
    MCAR_selection <- rbinom(n_obs, size=1, prob=psel)
    
    # Apply censorship to variable for every observation which is NOT selected
    data[, var][MCAR_selection == 0] <- NA
  }
  
  return (list(data, psel, MCAR_selection))
}


apply_MNAR_missingness <- function(data = NULL, vars_to_censor = NULL) {
  for (var in vars_to_censor) {
    # MNAR probability of censorship depends on data value
    psel <- rep(0.1, times = n_obs)
    
    for (i in 1:n_obs) {
      if (data[i, var] >= 1) {
        psel[i] <- 0.9
      }
      else if (data[i, var] <= -1.65) {
        psel[i] <- 0.7
      }
    }
    
    # MNAR selection vector
    MNAR_selection <- rbinom(n_obs, size=1, prob=psel)
    
    # Apply censorship to variable for every observation which is NOT selected
    data[, var][MNAR_selection == 0] <- NA
  }
  
  return (list(data, psel, MNAR_selection))
}


apply_MAR_missingness <- function(data = NULL, vars_to_censor = NULL) {
  # generate sub of all covariates in low, low subgroup (i.e Z_new = sum(Z1, ..., Z8))
  Z_sum_subgroup_LL <- rowSums(data[, c("Z2", "Z3", "Z4", "Z5", "Z6", "Z7", "Z8")])
  
  # find 75th percentile of this sum
  percentile_75 <- quantile(Z_sum_subgroup_LL, 0.75)
  psel <- 0.25
  
  # select into study population where value of sum is above 75th percentile
  MAR_selection <- (Z_sum_subgroup_LL >= percentile_75)
  
  # Apply censorship to variable for every observation which is NOT selected
  for (var in vars_to_censor) {
    data[, var][MAR_selection == 0] <- NA
  }
  
  return (list(data, psel, MAR_selection))
}


apply_outcome_MAR_missingness <- function(data = NULL, vars_to_censor = NULL) {
  # generate sub of all covariates in low, low subgroup (i.e Z_new = sum(Z1, ..., Z8))
  Z_sum_subgroup_LL <- rowSums(data[, c("Z2", "Z3", "Z4", "Z5", "Z6", "Z7", "Z8")])
  
  outcome <- data[, c("Y")]
  
  # find 75th percentile of this sum
  percentile_75 <- quantile(Z_sum_subgroup_LL, 0.75)
  psel <- 0.25
  
  # select into study population where value of sum is above 75th percentile
  MAR_selection <- ((Z_sum_subgroup_LL + outcome) >= percentile_75)
  
  # Apply censorship to variable for every observation which is NOT selected
  for (var in vars_to_censor) {
    data[, var][MAR_selection == 0] <- NA
  }
  
  return (list(data, psel, MAR_selection))
}


# ----- Complete data generation mechanism -----


# datasets representative of our DAG
generate_dataset <- function(n_obs      = NULL,
                             n_rep      = NULL,
                             
                             Z_correlation     = NULL,
                             Z_subgroups       = NULL,
                             target_r_sq_X     = NULL,
                             target_r_sq_Y     = NULL,
                             causal            = NULL,
                             
                             num_total_conf  = NULL,
                             num_meas_conf   = NULL,
                             num_unmeas_conf = NULL,
                             
                             vars_to_make_unmeasured = NULL,
                             vars_to_censor          = NULL,
                             var_names               = NULL) {
  
  # coefficient values for DAG
  alpha <- sqrt(Z_correlation)     # U on Z
  beta  <- sqrt(1 - Z_correlation) # error_Z on Z
  
  beta_X  <- beta_X_formula(num_total_conf = num_total_conf, # Z on X
                            target_r_sq_X  = target_r_sq_X,
                            Z_correlation  = Z_correlation)
  
  beta_Xs <- beta_X_subgroups_formula(beta_X = beta_X) # subgroup Z on X
  beta_Ys <- beta_Y_subgroups_formula(beta_X = beta_X) # subgroup Z on Y
  
  # variance of the error term for Y
  var_Y <- determine_subgroup_var_Y(num_total_conf = num_total_conf,
                                    beta_X         = beta_X,
                                    causal         = causal,
                                    Z_correlation  = Z_correlation,
                                    target_r_sq_Y  = target_r_sq_Y)
  
  var_error_Y <- determine_subgroup_var_error_Y(var_Y          = var_Y,
                                                target_r_sq_Y  = target_r_sq_Y)
  
  dataset <- data.frame(matrix(NaN, nrow = n_obs, ncol = length(var_names)))
  colnames(dataset) <- var_names
  
  # shared prior U for all Z_i
  prior_U <- rnorm(n = n_obs, mean = 0, sd = 1)
  
  # generate all confounders Z_i
  # low, low subgroup
  Z1  <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z2  <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z3  <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z4  <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z5  <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z6  <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z7  <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z8  <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  
  # low, high subgroup
  Z9  <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z10 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z11 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z12 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z13 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z14 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z15 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z16 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  
  # high, low subgroup
  Z17 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z18 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z19 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z20 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z21 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z22 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z23 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z24 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  
  # high, high subgroup
  Z25 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z26 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z27 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z28 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z29 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z30 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z31 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z32 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  
  # generate exposure X
  X <- (beta_Xs[1] * (Z1  + Z2  + Z3  + Z4  + Z5  + Z6  + Z7  + Z8))  + # low X, low Y
    (beta_Xs[2] * (Z9  + Z10 + Z11 + Z12 + Z13 + Z14 + Z15 + Z16)) +    # low X, high Y
    (beta_Xs[3] * (Z17 + Z18 + Z19 + Z20 + Z21 + Z22 + Z23 + Z24)) +    # high X, low Y
    (beta_Xs[4] * (Z25 + Z26 + Z27 + Z28 + Z29 + Z30 + Z31 + Z32))      # high X, high Y
  
  # generate outcome Y
  Y <- (beta_Ys[1] * (Z1  + Z2  + Z3  + Z4  + Z5  + Z6  + Z7  + Z8))  + # low X, low Y
    (beta_Ys[2] * (Z9  + Z10 + Z11 + Z12 + Z13 + Z14 + Z15 + Z16)) +    # low X, high Y
    (beta_Ys[3] * (Z17 + Z18 + Z19 + Z20 + Z21 + Z22 + Z23 + Z24)) +    # high X, low Y
    (beta_Ys[4] * (Z25 + Z26 + Z27 + Z28 + Z29 + Z30 + Z31 + Z32))      # high X, high Y
  
  # add error term to X if X
  error_X <- rnorm(n = n_obs, mean = 0, sd = 1)
  X <- X + error_X
  
  # add error term to Y if Y
  error_Y <- rnorm(n = n_obs, mean = 0, sd = sqrt(var_error_Y))
  Y <- Y + error_Y
  
  # add causal effect (X on Y)
  Y <- Y + (causal * X)
  
  # record X, Y and Z_i
  dataset[, 'Y']   <- Y
  dataset[, 'X']   <- X
  
  # Low, low subgroup
  dataset[, 'Z1']  <- Z1
  dataset[, 'Z2']  <- Z2
  dataset[, 'Z3']  <- Z3
  dataset[, 'Z4']  <- Z4
  dataset[, 'Z5']  <- Z5
  dataset[, 'Z6']  <- Z6
  dataset[, 'Z7']  <- Z7
  dataset[, 'Z8']  <- Z8
  
  # low, high subgroup
  dataset[, 'Z9']  <- Z9
  dataset[, 'Z10'] <- Z10
  dataset[, 'Z11'] <- Z11
  dataset[, 'Z12'] <- Z12
  dataset[, 'Z13'] <- Z13
  dataset[, 'Z14'] <- Z14
  dataset[, 'Z15'] <- Z15
  dataset[, 'Z16'] <- Z16
  
  # high, low subgroup
  dataset[, 'Z17'] <- Z17
  dataset[, 'Z18'] <- Z18
  dataset[, 'Z19'] <- Z19
  dataset[, 'Z20'] <- Z20
  dataset[, 'Z21'] <- Z21
  dataset[, 'Z22'] <- Z22
  dataset[, 'Z23'] <- Z23
  dataset[, 'Z24'] <- Z24
  
  # high, high subgroup
  dataset[, 'Z25'] <- Z25
  dataset[, 'Z26'] <- Z26
  dataset[, 'Z27'] <- Z27
  dataset[, 'Z28'] <- Z28
  dataset[, 'Z29'] <- Z29
  dataset[, 'Z30'] <- Z30
  dataset[, 'Z31'] <- Z31
  dataset[, 'Z32'] <- Z32
  
  # censor covariates as appropriate
  dataset <- dataset[, !(names(dataset) %in% vars_to_make_unmeasured)]
  
  return (dataset)
}

