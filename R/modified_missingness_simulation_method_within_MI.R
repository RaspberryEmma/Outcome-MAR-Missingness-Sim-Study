# ****************************************
# Missingness Simulation Study
#
# Simulation procedure for confounder-handling
# with missingness considerations
#
# Emma Tarmey
#
# Started:          30/01/2026
# Most Recent Edit: 01/02/2026
# ****************************************


# ----- Missingness mechanisms and missingness handling -----

# See: https://www.rdocumentation.org/packages/mice/versions/3.17.0/topics/mice
apply_within_MI <- function(data = NULL, num_datasets = NULL, repetitions = NULL, imp_method = NULL) {
  
  # # find all vars containing missingness to be imputed
  # if (binary_Z) {
  #   vars_with_missingness <- colnames(data)[ apply(data, 2, anyNA) ]
  #   for (var in vars_with_missingness) {
  #     data[, var] <- as.factor(data[, var])
  #   }
  # }
  
  capture.output(                      # suppress command line output
    imp <- mice(data,
                m      = num_datasets, # number of imputations
                maxit  = repetitions,  # number of iterations
                method = imp_method)
  )
  
  return (imp)
}

vars_selected_string_to_binary <- function(vars_selected = NULL, var_names = NULL) {
  vars_binary        <- rep(0, times = length(var_names))
  names(vars_binary) <- var_names
  
  for (var in vars_selected) {
    vars_binary[var] <- 1.0
  }
  
  return (vars_binary)
}

vars_selected_MI_means_to_binary <- function(MI_means = NULL) {
  vars_binary        <- rep(0, times = length(MI_means))
  
  for (i in c(1:length(MI_means))) {
    if (MI_means[i] >= 0.5) {
      vars_binary[i] <- 1.0
    }
  }
  
  # include intercept here
  vars_binary        <- c(1.0, vars_binary)
  names(vars_binary) <- c("(Intercept)", names(MI_means))
  
  return (vars_binary)
}

run_within_MI_simulation <- function(n_scenario = NULL,
                                      n_obs      = NULL,
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
                                      missingness_mechanism   = NULL) {

  # ----- Recording results -----
  
  var_sel_methods <- c("fully_adjusted", "unadjusted", "two_step_lasso", "two_step_lasso_X", "two_step_lasso_union")
  
  results_methods <- c("causal_true_value", "causal_estimate", "causal_bias", "causal_bias_proportion", "causal_coverage",
                       "open_paths", "blocked_paths", "proportion_paths",
                       "empirical_SE", "model_SE")
  
  var_names                         <- c("Y", "X", paste('Z', c(1:num_total_conf), sep=''))
  var_names_except_Y                <- var_names[ !var_names == 'Y']
  var_names_except_Y_with_intercept <- c("(Intercept)", var_names_except_Y)
  
  n_variables       <- length(var_names)
  n_var_sel_methods <- length(var_sel_methods)
  n_results         <- length(results_methods)
  
  missingness_coefs         <- array(data     = NaN,
                              dim      = c(n_var_sel_methods, (n_variables - 1 + 1), n_rep),
                              dimnames = list(var_sel_methods, var_names_except_Y_with_intercept, 1:n_rep) )
  missingness_cov_selection <- array(data     = NaN,
                              dim      = c(n_var_sel_methods, (n_variables - 1 + 1), n_rep),
                              dimnames = list(var_sel_methods, var_names_except_Y_with_intercept, 1:n_rep) )
  missingness_results       <- array(data     = NaN,
                              dim      = c(n_var_sel_methods, n_results, n_rep),
                              dimnames = list(var_sel_methods, results_methods, 1:n_rep))
  
  sample_size_table <- array(data     = NaN,
                             dim      = c(3, 2),
                             dimnames = list(c("None", "missingness", "MCAR"), c("complete_cases", "sample_size_after_handling")))
  
  
  # ----- Simulation procedure ------
  
  message("\n\n\n ***** Running simulation procedure, scenario ", n_scenario, " *****")
  
  # track time
  current_time  <- as.numeric(Sys.time())*1000
  remaining_rep <- n_rep + 1
  
  for (repetition in c(1:n_rep)) {
    # track time
    remaining_rep   <- remaining_rep - 1
    previous_time   <- current_time
    current_time    <- as.numeric(Sys.time())*1000
    remaining_time  <- remaining_rep * (current_time - previous_time)
    expected_finish <- signif(((current_time + remaining_time)/1000), digits = 13) %>% as.POSIXct()
    
    print(paste0("Running repetition ", repetition, "/", n_rep, ", Expected to finish running at: ", expected_finish))
  
    FULL_dataset <- generate_dataset(n_obs      = n_obs,
                                     n_rep      = n_rep,
                                     
                                     Z_correlation     = Z_correlation,
                                     Z_subgroups       = Z_subgroups,
                                     target_r_sq_X     = target_r_sq_X,
                                     target_r_sq_Y     = target_r_sq_Y,
                                     causal            = causal,
                                     
                                     num_total_conf  = num_total_conf,
                                     num_meas_conf   = num_meas_conf,
                                     num_unmeas_conf = num_unmeas_conf,
                                     
                                     vars_to_make_unmeasured = vars_to_make_unmeasured,
                                     vars_to_censor          = vars_to_censor,
                                     var_names               = var_names)
    
    if (missingness_mechanism == "MNAR") {
      missingness_data <- apply_MNAR_missingness(FULL_dataset, vars_to_censor = vars_to_censor)
    } else if (missingness_mechanism == "MCAR") {
      missingness_data <- apply_MCAR_missingness(FULL_dataset, vars_to_censor = vars_to_censor)
    } else if (missingness_mechanism == "MAR") {
      missingness_data <- apply_MAR_missingness(FULL_dataset, vars_to_censor = vars_to_censor)
    } else if (missingness_mechanism == "outcome_MAR") {
      missingness_data <- apply_outcome_MAR_missingness(FULL_dataset, vars_to_censor = vars_to_censor)
    } else {
      stop("TODO: full data here, structure needs slight work")
    }
    
    missingness_dataset    <- missingness_data[[1]]
    missingness_psel       <- missingness_data[[2]]
    missingness_censorship <- missingness_data[[3]]
    
    # apply within MI
    # imp_method = norm -> MI using Bayesian linear regression
    # see: section "Details" of https://www.rdocumentation.org/packages/mice/versions/3.17.0/topics/mice
    handled_missingness_imputation_object <- apply_within_MI(data         = missingness_dataset,
                                             num_datasets = 20,
                                             repetitions  = 20,
                                             imp_method   = "norm")
    
    # record sample sizes before and after missingness handling is applied
    # NB: for within MI, all imputed datasets will have the same sample size
    sample_size_table["missingness", "complete_cases"]             <- dim(missingness_dataset[complete.cases(missingness_dataset), ])[1]
    sample_size_table["missingness", "sample_size_after_handling"] <- dim(handled_missingness_imputation_object$data)[1]
    
    # dataframe to store variable selection results across all imputed datasets
    lasso_missingness_vars_selected_all_imputations <- data.frame()
    
    # fully-adjusted
    all_available_vars <- var_names_except_Y[var_names_except_Y %!in% vars_to_make_unmeasured]
    
    fully_adjusted_missingness_models   <- with(handled_missingness_imputation_object,
                                                eval(parse(text=paste0("glm(", make_model_formula(vars_selected = all_available_vars), ", family = 'gaussian')"))))
    fully_adjusted_missingness_estimate <- summary(pool(fully_adjusted_missingness_models))
    
    # unadjusted
    unadjusted_missingness_models   <- with(handled_missingness_imputation_object,
                                            eval(parse(text=paste0("glm(Y ~ X, family = 'gaussian')"))))
    unadjusted_missingness_estimate <- summary(pool(unadjusted_missingness_models))
    
    # MI LASSO
    for (i in c(1:20)) {
      handled_missingness_dataset <- complete(handled_missingness_imputation_object, action = i)
      
      # cut up versions of the data as needed
      X_handled_missingness_dataset <- subset(handled_missingness_dataset, select=-c(Y))
      Z_handled_missingness_dataset <- subset(handled_missingness_dataset, select=-c(Y, X))
      Y_handled_missingness_column  <- subset(handled_missingness_dataset, select=c(Y))
      X_handled_missingness_column  <- subset(handled_missingness_dataset, select=c(X))
      
      # MI LASSO
      cv_lasso_missingness_model <- cv.glmnet(x      = data.matrix(X_handled_missingness_dataset),
                                                y      = data.matrix(Y_handled_missingness_column),
                                                family = 'gaussian',
                                                alpha  = 1)
      
      lambda <- cv_lasso_missingness_model$lambda.min
      
      lasso_missingness_model <- glmnet(x      = data.matrix(X_handled_missingness_dataset),
                                          y      = data.matrix(Y_handled_missingness_column),
                                          family = 'gaussian',
                                          alpha  = 1,
                                          lambda = lambda)
      
      # Extract covariate selection
      lasso_missingness_coefs         <- as.vector(lasso_missingness_model$beta)
      names(lasso_missingness_coefs)  <- rownames(lasso_missingness_model$beta)
      lasso_missingness_vars_selected <- union(c('X'), names(lasso_missingness_coefs[lasso_missingness_coefs != 0.0])) # always select X
      lasso_missingness_vars_selected <- lasso_missingness_vars_selected[lasso_missingness_vars_selected != "(Intercept)"] # exclude intercept
      lasso_missingness_vars_selected <- vars_selected_string_to_binary(vars_selected = lasso_missingness_vars_selected, var_names = var_names_except_Y)
      
      # Store covariate selection
      lasso_missingness_vars_selected_all_imputations <- rbind(lasso_missingness_vars_selected_all_imputations, lasso_missingness_vars_selected)
    }
    
    # Aggregate covariate selection across all imputed datasets
    colnames(lasso_missingness_vars_selected_all_imputations) <- var_names_except_Y
    lasso_missingness_vars_selected_mean                      <- colMeans(lasso_missingness_vars_selected_all_imputations, na.rm = TRUE)
    
    lasso_missingness_vars_selected_more_than_half <- names(lasso_missingness_vars_selected_mean[lasso_missingness_vars_selected_mean >= 0.5])
    lasso_missingness_formula                      <- paste0("glm(", make_model_formula(vars_selected = lasso_missingness_vars_selected_more_than_half), ", family = 'gaussian')")
    
    two_step_LASSO_missingness_models   <- with(handled_missingness_imputation_object,
                                              eval(parse(text=lasso_missingness_formula)))
    two_step_LASSO_missingness_estimate <- summary(pool(two_step_LASSO_missingness_models))
    

    # dataframe to store variable selection results across all imputed datasets
    lasso_X_missingness_vars_selected_all_imputations <- data.frame()
    
    # MI LASSO_X
    for (i in c(1:20)) {
      handled_missingness_dataset <- complete(handled_missingness_imputation_object, action = i)
      
      # cut up versions of the data as needed
      X_handled_missingness_dataset <- subset(handled_missingness_dataset, select=-c(Y))
      Z_handled_missingness_dataset <- subset(handled_missingness_dataset, select=-c(Y, X))
      Y_handled_missingness_column  <- subset(handled_missingness_dataset, select=c(Y))
      X_handled_missingness_column  <- subset(handled_missingness_dataset, select=c(X))
      
      # MI LASSO_X
      cv_lasso_X_missingness_model <- cv.glmnet(x      = data.matrix(Z_handled_missingness_dataset),
                                                y      = data.matrix(X_handled_missingness_column),
                                                family = 'gaussian',
                                                alpha  = 1)
      
      lambda <- cv_lasso_X_missingness_model$lambda.min
      
      lasso_X_missingness_model <- glmnet(x      = data.matrix(Z_handled_missingness_dataset),
                                          y      = data.matrix(X_handled_missingness_column),
                                          family = 'gaussian',
                                          alpha  = 1,
                                          lambda = lambda)
      
      # Extract covariate selection
      lasso_X_missingness_coefs         <- as.vector(lasso_X_missingness_model$beta)
      names(lasso_X_missingness_coefs)  <- rownames(lasso_X_missingness_model$beta)
      lasso_X_missingness_vars_selected <- union(c('X'), names(lasso_X_missingness_coefs[lasso_X_missingness_coefs != 0.0])) # always select X
      lasso_X_missingness_vars_selected <- lasso_X_missingness_vars_selected[lasso_X_missingness_vars_selected != "(Intercept)"] # exclude intercept
      lasso_X_missingness_vars_selected <- vars_selected_string_to_binary(vars_selected = lasso_X_missingness_vars_selected, var_names = var_names_except_Y)
      
      # Store covariate selection
      lasso_X_missingness_vars_selected_all_imputations <- rbind(lasso_X_missingness_vars_selected_all_imputations, lasso_X_missingness_vars_selected)
    }
    
    # Aggregate covariate selection across all imputed datasets
    colnames(lasso_X_missingness_vars_selected_all_imputations) <- var_names_except_Y
    lasso_X_missingness_vars_selected_mean                      <- colMeans(lasso_X_missingness_vars_selected_all_imputations, na.rm = TRUE)
    
    lasso_X_missingness_vars_selected_more_than_half <- names(lasso_X_missingness_vars_selected_mean[lasso_X_missingness_vars_selected_mean >= 0.5])
    lasso_X_missingness_formula                      <- paste0("glm(", make_model_formula(vars_selected = lasso_X_missingness_vars_selected_more_than_half), ", family = 'gaussian')")
    
    two_step_LASSO_X_missingness_models   <- with(handled_missingness_imputation_object,
                                                  eval(parse(text=lasso_X_missingness_formula)))
    two_step_LASSO_X_missingness_estimate <- summary(pool(two_step_LASSO_X_missingness_models))
    

    # MI LASSO_UNION
    
    # add binary results of lasso and lasso-X together element-wise
    lasso_union_missingness_vars_selected_all_imputations <- Reduce("+", list( lasso_missingness_vars_selected_all_imputations, lasso_X_missingness_vars_selected_all_imputations))
    
    # if a variable is selected by both methods, we store this as binary 1 instead of 2
    lasso_union_missingness_vars_selected_all_imputations[lasso_union_missingness_vars_selected_all_imputations == 2] <- 1
    
    # Aggregate covariate selection across all imputed datasets
    lasso_union_missingness_vars_selected_mean           <- colMeans(lasso_union_missingness_vars_selected_all_imputations, na.rm = TRUE)
    lasso_union_missingness_vars_selected_more_than_half <- names(lasso_union_missingness_vars_selected_mean[lasso_union_missingness_vars_selected_mean >= 0.5])
    lasso_union_missingness_formula                      <- paste0("glm(", make_model_formula(vars_selected = lasso_union_missingness_vars_selected_more_than_half), ", family = 'gaussian')")
    
    two_step_LASSO_union_missingness_models   <- with(handled_missingness_imputation_object,
                                                      eval(parse(text=lasso_union_missingness_formula)))
    two_step_LASSO_union_missingness_estimate <- summary(pool(two_step_LASSO_union_missingness_models))
    
    
    # ----- Record covariate selection -----
    
    missingness_cov_selection["fully_adjusted", , repetition]       <- vars_selected_string_to_binary(vars_selected = var_names_except_Y,                                   var_names = var_names_except_Y_with_intercept)
    missingness_cov_selection["unadjusted", , repetition]           <- vars_selected_string_to_binary(vars_selected = c(),                                                  var_names = var_names_except_Y_with_intercept)
    missingness_cov_selection["two_step_lasso", , repetition]       <- vars_selected_string_to_binary(vars_selected = lasso_missingness_vars_selected_more_than_half,       var_names = var_names_except_Y_with_intercept)
    missingness_cov_selection["two_step_lasso_X", , repetition]     <- vars_selected_string_to_binary(vars_selected = lasso_X_missingness_vars_selected_more_than_half,     var_names = var_names_except_Y_with_intercept)
    missingness_cov_selection["two_step_lasso_union", , repetition] <- vars_selected_string_to_binary(vars_selected = lasso_union_missingness_vars_selected_more_than_half, var_names = var_names_except_Y_with_intercept)
    
    
    # ------ Record coefficients ------
    
    fully_adjusted_missingness_estimate_coefs        <- fully_adjusted_missingness_estimate$estimate
    names(fully_adjusted_missingness_estimate_coefs) <- fully_adjusted_missingness_estimate$term
    
    unadjusted_missingness_estimate_coefs        <- unadjusted_missingness_estimate$estimate
    names(unadjusted_missingness_estimate_coefs) <- unadjusted_missingness_estimate$term
    
    two_step_LASSO_missingness_estimate_coefs        <- two_step_LASSO_missingness_estimate$estimate
    names(two_step_LASSO_missingness_estimate_coefs) <- two_step_LASSO_missingness_estimate$term
    
    two_step_LASSO_X_missingness_estimate_coefs        <- two_step_LASSO_X_missingness_estimate$estimate
    names(two_step_LASSO_X_missingness_estimate_coefs) <- two_step_LASSO_X_missingness_estimate$term
    
    two_step_LASSO_union_missingness_estimate_coefs        <- two_step_LASSO_union_missingness_estimate$estimate
    names(two_step_LASSO_union_missingness_estimate_coefs) <- two_step_LASSO_union_missingness_estimate$term
    
    missingness_coefs["fully_adjusted", , repetition]       <- fill_in_blanks(fully_adjusted_missingness_estimate_coefs, var_names_except_Y_with_intercept)
    missingness_coefs["unadjusted", , repetition]           <- fill_in_blanks(unadjusted_missingness_estimate_coefs, var_names_except_Y_with_intercept)
    missingness_coefs["two_step_lasso", , repetition]       <- fill_in_blanks(two_step_LASSO_missingness_estimate_coefs, var_names_except_Y_with_intercept)
    missingness_coefs["two_step_lasso_X", , repetition]     <- fill_in_blanks(two_step_LASSO_X_missingness_estimate_coefs, var_names_except_Y_with_intercept)
    missingness_coefs["two_step_lasso_union", , repetition] <- fill_in_blanks(two_step_LASSO_union_missingness_estimate_coefs, var_names_except_Y_with_intercept)
    
    
    # ----- Record results -----
    
    missingness_results["fully_adjusted", "causal_true_value", repetition]      <- causal
    missingness_results["fully_adjusted", "causal_estimate", repetition]        <- fully_adjusted_missingness_estimate[which(fully_adjusted_missingness_estimate$term == "X")[[1]], "estimate"]
    missingness_results["fully_adjusted", "causal_bias", repetition]            <- (fully_adjusted_missingness_estimate[which(fully_adjusted_missingness_estimate$term == "X")[[1]], "estimate"] - causal)
    missingness_results["fully_adjusted", "causal_bias_proportion", repetition] <- ((fully_adjusted_missingness_estimate[which(fully_adjusted_missingness_estimate$term == "X")[[1]], "estimate"] - causal) / causal)
    missingness_results["fully_adjusted", "causal_coverage", repetition]        <- estimate_within_CI(estimate       = fully_adjusted_missingness_estimate[which(fully_adjusted_missingness_estimate$term == "X")[[1]], "estimate"],
                                                                                                      true_value     = causal,
                                                                                                      standard_error = fully_adjusted_missingness_estimate[which(fully_adjusted_missingness_estimate$term == "X")[[1]], "std.error"],
                                                                                                      sample_size    = dim(handled_missingness_imputation_object$data)[1])
    missingness_results["fully_adjusted", "open_paths", repetition]             <- num_total_conf
    missingness_results["fully_adjusted", "blocked_paths", repetition]          <- num_meas_conf
    missingness_results["fully_adjusted", "proportion_paths", repetition]       <- num_meas_conf / num_total_conf
    missingness_results["fully_adjusted", "empirical_SE", repetition]           <- NaN
    missingness_results["fully_adjusted", "model_SE", repetition]               <- fully_adjusted_missingness_estimate[which(fully_adjusted_missingness_estimate$term == "X")[[1]], "std.error"]
    
    missingness_results["unadjusted", "causal_true_value", repetition]      <- causal
    missingness_results["unadjusted", "causal_estimate", repetition]        <- unadjusted_missingness_estimate[which(unadjusted_missingness_estimate$term == "X")[[1]], "estimate"]
    missingness_results["unadjusted", "causal_bias", repetition]            <- (unadjusted_missingness_estimate[which(unadjusted_missingness_estimate$term == "X")[[1]], "estimate"] - causal)
    missingness_results["unadjusted", "causal_bias_proportion", repetition] <- ((unadjusted_missingness_estimate[which(unadjusted_missingness_estimate$term == "X")[[1]], "estimate"] - causal) / causal)
    missingness_results["unadjusted", "causal_coverage", repetition]        <- estimate_within_CI(estimate       = unadjusted_missingness_estimate[which(unadjusted_missingness_estimate$term == "X")[[1]], "estimate"],
                                                                                                  true_value     = causal,
                                                                                                  standard_error = unadjusted_missingness_estimate[which(unadjusted_missingness_estimate$term == "X")[[1]], "std.error"],
                                                                                                  sample_size    = dim(handled_missingness_imputation_object$data)[1])
    missingness_results["unadjusted", "open_paths", repetition]             <- num_total_conf
    missingness_results["unadjusted", "blocked_paths", repetition]          <- 0.0
    missingness_results["unadjusted", "proportion_paths", repetition]       <- 0.0 / num_total_conf
    missingness_results["unadjusted", "empirical_SE", repetition]           <- NaN
    missingness_results["unadjusted", "model_SE", repetition]               <- unadjusted_missingness_estimate[which(unadjusted_missingness_estimate$term == "X")[[1]], "std.error"]
    
    missingness_results["two_step_lasso", "causal_true_value", repetition]      <- causal
    missingness_results["two_step_lasso", "causal_estimate", repetition]        <- two_step_LASSO_missingness_estimate[which(two_step_LASSO_missingness_estimate$term == "X")[[1]], "estimate"]
    missingness_results["two_step_lasso", "causal_bias", repetition]            <- (two_step_LASSO_missingness_estimate[which(two_step_LASSO_missingness_estimate$term == "X")[[1]], "estimate"] - causal)
    missingness_results["two_step_lasso", "causal_bias_proportion", repetition] <- ((two_step_LASSO_missingness_estimate[which(two_step_LASSO_missingness_estimate$term == "X")[[1]], "estimate"] - causal) / causal)
    missingness_results["two_step_lasso", "causal_coverage", repetition]        <- estimate_within_CI(estimate       = two_step_LASSO_missingness_estimate[which(two_step_LASSO_missingness_estimate$term == "X")[[1]], "estimate"],
                                                                                                      true_value     = causal,
                                                                                                      standard_error = two_step_LASSO_missingness_estimate[which(two_step_LASSO_missingness_estimate$term == "X")[[1]], "std.error"],
                                                                                                      sample_size    = dim(handled_missingness_imputation_object$data)[1])
    missingness_results["two_step_lasso", "open_paths", repetition]             <- num_total_conf
    missingness_results["two_step_lasso", "blocked_paths", repetition]          <- length(lasso_missingness_vars_selected_more_than_half[lasso_missingness_vars_selected_more_than_half != 'X'])
    missingness_results["two_step_lasso", "proportion_paths", repetition]       <- length(lasso_missingness_vars_selected_more_than_half[lasso_missingness_vars_selected_more_than_half != 'X']) / num_total_conf
    missingness_results["two_step_lasso", "empirical_SE", repetition]           <- NaN
    missingness_results["two_step_lasso", "model_SE", repetition]               <- two_step_LASSO_missingness_estimate[which(two_step_LASSO_missingness_estimate$term == "X")[[1]], "std.error"]
    
    missingness_results["two_step_lasso_X", "causal_true_value", repetition]      <- causal
    missingness_results["two_step_lasso_X", "causal_estimate", repetition]        <- two_step_LASSO_X_missingness_estimate[which(two_step_LASSO_X_missingness_estimate$term == "X")[[1]], "estimate"]
    missingness_results["two_step_lasso_X", "causal_bias", repetition]            <- (two_step_LASSO_X_missingness_estimate[which(two_step_LASSO_X_missingness_estimate$term == "X")[[1]], "estimate"] - causal)
    missingness_results["two_step_lasso_X", "causal_bias_proportion", repetition] <- (two_step_LASSO_X_missingness_estimate[which(two_step_LASSO_X_missingness_estimate$term == "X")[[1]], "estimate"] - causal) / causal
    missingness_results["two_step_lasso_X", "causal_coverage", repetition]        <- estimate_within_CI(estimate       = two_step_LASSO_X_missingness_estimate[which(two_step_LASSO_X_missingness_estimate$term == "X")[[1]], "estimate"],
                                                                                                        true_value     = causal,
                                                                                                        standard_error = two_step_LASSO_X_missingness_estimate[which(two_step_LASSO_X_missingness_estimate$term == "X")[[1]], "std.error"],
                                                                                                        sample_size    = dim(handled_missingness_imputation_object$data)[1])
    missingness_results["two_step_lasso_X", "open_paths", repetition]             <- num_total_conf
    missingness_results["two_step_lasso_X", "blocked_paths", repetition]          <- length(lasso_X_missingness_vars_selected_more_than_half[lasso_X_missingness_vars_selected_more_than_half != 'X'])
    missingness_results["two_step_lasso_X", "proportion_paths", repetition]       <- length(lasso_X_missingness_vars_selected_more_than_half[lasso_X_missingness_vars_selected_more_than_half != 'X']) / num_total_conf
    missingness_results["two_step_lasso_X", "empirical_SE", repetition]           <- NaN
    missingness_results["two_step_lasso_X", "model_SE", repetition]               <- two_step_LASSO_X_missingness_estimate[which(two_step_LASSO_X_missingness_estimate$term == "X")[[1]], "std.error"]
    
    missingness_results["two_step_lasso_union", "causal_true_value", repetition]      <- causal
    missingness_results["two_step_lasso_union", "causal_estimate", repetition]        <- two_step_LASSO_union_missingness_estimate[which(two_step_LASSO_union_missingness_estimate$term == "X")[[1]], "estimate"]
    missingness_results["two_step_lasso_union", "causal_bias", repetition]            <- (two_step_LASSO_union_missingness_estimate[which(two_step_LASSO_union_missingness_estimate$term == "X")[[1]], "estimate"] - causal)
    missingness_results["two_step_lasso_union", "causal_bias_proportion", repetition] <- (two_step_LASSO_union_missingness_estimate[which(two_step_LASSO_union_missingness_estimate$term == "X")[[1]], "estimate"] - causal) / causal
    missingness_results["two_step_lasso_union", "causal_coverage", repetition]        <- estimate_within_CI(estimate       = two_step_LASSO_union_missingness_estimate[which(two_step_LASSO_union_missingness_estimate$term == "X")[[1]], "estimate"],
                                                                                                            true_value     = causal,
                                                                                                            standard_error = two_step_LASSO_union_missingness_estimate[which(two_step_LASSO_union_missingness_estimate$term == "X")[[1]], "std.error"],
                                                                                                            sample_size    = dim(handled_missingness_imputation_object$data)[1])
    missingness_results["two_step_lasso_union", "open_paths", repetition]             <- num_total_conf
    missingness_results["two_step_lasso_union", "blocked_paths", repetition]          <- length(lasso_union_missingness_vars_selected_more_than_half[lasso_union_missingness_vars_selected_more_than_half != 'X'])
    missingness_results["two_step_lasso_union", "proportion_paths", repetition]       <- length(lasso_union_missingness_vars_selected_more_than_half[lasso_union_missingness_vars_selected_more_than_half != 'X']) / num_total_conf
    missingness_results["two_step_lasso_union", "empirical_SE", repetition]           <- NaN
    missingness_results["two_step_lasso_union", "model_SE", repetition]               <- two_step_LASSO_union_missingness_estimate[which(two_step_LASSO_union_missingness_estimate$term == "X")[[1]], "std.error"]
    
  } # end repetitions loop
  
  # "level" coefficient for effect of Zs on X and Y before subgroups
  beta_X  <- beta_X_formula(num_total_conf = num_total_conf, # Z on X
                            target_r_sq_X  = target_r_sq_X,
                            Z_correlation  = Z_correlation)
  
  # variance of the error term for Y
  var_Y <- determine_subgroup_var_Y(num_total_conf = num_total_conf,
                                    beta_X         = beta_X,
                                    causal         = causal,
                                    Z_correlation  = Z_correlation,
                                    target_r_sq_Y  = target_r_sq_Y)
  
  TRUE_coefs        <- c(beta_X_subgroups_formula(beta_X = beta_X),
                         beta_Y_subgroups_formula(beta_X = beta_X),
                         causal,
                         determine_subgroup_var_error_Y(var_Y = var_Y, target_r_sq_Y  = target_r_sq_Y))
  
  names(TRUE_coefs) <- c("Z on X (Subgroup Low X, Low Y)",
                         "Z on X (Subgroup Low X, High Y)",
                         "Z on X (Subgroup High X, Low Y)",
                         "Z on X (Subgroup High X, High Y)",
                         
                         "Z on Y (Subgroup Low X, Low Y)",
                         "Z on Y (Subgroup Low X, High Y)",
                         "Z on Y (Subgroup High X, Low Y)",
                         "Z on Y (Subgroup High X, High Y)",
                         
                         "Causal effect X on Y",
                         "Var in error for Y")
  
  return (list(missingness_cov_selection,
               TRUE_coefs,
               missingness_coefs,
               missingness_results,
               sample_size_table,
               handled_missingness_dataset
  ))
} # function


