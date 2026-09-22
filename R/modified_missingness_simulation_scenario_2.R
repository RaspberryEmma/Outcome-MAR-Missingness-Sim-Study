# ****************************************
# Missingness Simulation Study
#
# Simulation procedure for confounder-handling
# with missingness considerations
# Scenario 2: within MI missingness handling
#
# Emma Tarmey
#
# Started:          06/10/2025
# Most Recent Edit: 16/02/2026
# ****************************************



# ----- Tech Setup ----- 

# clear R memory
rm(list=ls())

# check all external libraries
using<-function(...) {
  libs <- unlist(list(...))
  req  <- unlist(lapply(libs, require, character.only=TRUE))
  need <- libs[req==FALSE]
  if(length(need) > 0){ 
    install.packages(need)
    lapply(need, require, character.only=TRUE)
  }
}
using("dplyr", "glmnet", "mice", "speedglm", "tidyr")

# fix wd issue
# forces wd to be the location of this file
if (Sys.getenv("RSTUDIO") == "1") {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

source("common.R")



# ----- Parameters ------

n_scenario   <- 2

n_obs             <- 10000
n_rep             <- 2000
Z_correlation     <- 0.1
Z_subgroups       <- 4
target_r_sq_X     <- 0.3
target_r_sq_Y     <- 0.3
causal            <- 0.5

num_total_conf  <- 32
num_meas_conf   <- 28
num_unmeas_conf <- 4

# missingness handling mechanism
missingness_handling <- "within_MI"

# confounders to be unmeasured
#vars_to_make_unmeasured <- c()
vars_to_make_unmeasured <- c("Z1", "Z9", "Z17", "Z25")

# confounders to have missingness applied
vars_to_censor <- c("Z26")


# ----- RNG -----

# # test RNG
# set.seed(2025)

# fix RNG seed based on current scenario
seeds_df <- read.csv(file = "../data/precomputed_RNG_seeds.csv")
seed     <- seeds_df %>% filter(simulation_scenario == n_scenario)
set.seed(seed$seed)



# ------ Run simulation procedure ------

# within_MI
source("modified_missingness_simulation_method_within_MI.R")
simulation_results_MNAR <- run_within_MI_simulation(n_scenario = n_scenario,
                                                   n_obs      = n_obs,
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
                                                   missingness_mechanism   = "MNAR")

simulation_results_MCAR <- run_within_MI_simulation(n_scenario = n_scenario,
                                                   n_obs      = n_obs,
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
                                                   missingness_mechanism   = "MCAR")

simulation_results_MAR <- run_within_MI_simulation(n_scenario = n_scenario,
                                                  n_obs      = n_obs,
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
                                                  missingness_mechanism   = "outcome_MAR")



# ----- Extract results -----

MNAR_cov_selection <- simulation_results_MNAR[[1]]
MCAR_cov_selection <- simulation_results_MCAR[[1]]
MAR_cov_selection  <- simulation_results_MAR[[1]]

TRUE_coefs <- simulation_results_MNAR[[2]]

MNAR_coefs <- simulation_results_MNAR[[3]]
MCAR_coefs <- simulation_results_MCAR[[3]]
MAR_coefs  <- simulation_results_MAR[[3]]

MNAR_results <- simulation_results_MNAR[[4]]
MCAR_results <- simulation_results_MCAR[[4]]
MAR_results  <- simulation_results_MAR[[4]]

sample_size_table_MNAR <- simulation_results_MNAR[[5]]
sample_size_table_MCAR <- simulation_results_MCAR[[5]]
sample_size_table_MAR  <- simulation_results_MAR[[5]]

sample_handled_MNAR_dataset <- simulation_results_MNAR[[6]]
sample_handled_MCAR_dataset <- simulation_results_MCAR[[6]]
sample_handled_MAR_dataset  <- simulation_results_MCAR[[6]]


# ----- Save results -----

# take mean across all repetitions
final_MNAR_cov_selection <- as.data.frame(apply(MNAR_cov_selection, c(1,2), mean))
final_MCAR_cov_selection <- as.data.frame(apply(MCAR_cov_selection, c(1,2), mean))
final_MAR_cov_selection  <- as.data.frame(apply(MAR_cov_selection, c(1,2), mean))

final_MNAR_coefs <- as.data.frame(apply(MNAR_coefs, c(1,2), mean))
final_MCAR_coefs <- as.data.frame(apply(MCAR_coefs, c(1,2), mean))
final_MAR_coefs  <- as.data.frame(apply(MAR_coefs, c(1,2), mean))

final_MNAR_results <- as.data.frame(apply(MNAR_results, c(1,2), mean))
final_MCAR_results <- as.data.frame(apply(MCAR_results, c(1,2), mean))
final_MAR_results  <- as.data.frame(apply(MAR_results, c(1,2), mean))


# record empirical standard error separately
causal        <- 0.5
model_methods <- c("fully_adjusted", "unadjusted", "two_step_lasso", "two_step_lasso_X", "two_step_lasso_union")
for (method in model_methods) {
  
  # MNAR
  MNAR_causal_effect_estimates                <- c(MNAR_results[method, "causal_estimate", ])
  final_MNAR_results[ method, "empirical_SE"] <- sd(MNAR_causal_effect_estimates)
  
  # MCAR
  MCAR_causal_effect_estimates                <- c(MCAR_results[method, "causal_estimate", ])
  final_MCAR_results[ method, "empirical_SE"] <- sd(MCAR_causal_effect_estimates)
  
  # MAR
  MAR_causal_effect_estimates                <- c(MAR_results[method, "causal_estimate", ])
  final_MAR_results[ method, "empirical_SE"] <- sd(MAR_causal_effect_estimates)
}

# save simulation setup
sim_setup <- c(n_obs,
               n_rep,
               Z_correlation,
               Z_subgroups,
               target_r_sq_X,
               target_r_sq_Y,
               causal,
               num_total_conf,
               num_meas_conf,
               num_unmeas_conf,
               missingness_handling)

names(sim_setup) <- c("n_obs",
                      "n_rep",
                      "Z_correlation",
                      "Z_subgroups",
                      "target_r_sq_X",
                      "target_r_sq_Y",
                      "causal",
                      "num_total_conf",
                      "num_meas_conf",
                      "num_unmeas_conf",
                      "missingness_handling")


# Save to file
id_string <- paste("method_", missingness_handling, sep='')

write.csv(final_MNAR_cov_selection,   paste("../data/scenario_", n_scenario, "/", id_string, "_MNAR_cov_selection.csv", sep=''))
write.csv(final_MCAR_cov_selection,   paste("../data/scenario_", n_scenario, "/", id_string, "_MCAR_cov_selection.csv", sep=''))
write.csv(final_MAR_cov_selection,    paste("../data/scenario_", n_scenario, "/", id_string, "_MAR_cov_selection.csv", sep=''))

write.csv(TRUE_coefs,         paste("../data/scenario_", n_scenario, "/", id_string, "_TRUE_coefs.csv", sep=''))
write.csv(final_MNAR_coefs,   paste("../data/scenario_", n_scenario, "/", id_string, "_MNAR_coefs.csv", sep=''))
write.csv(final_MCAR_coefs,   paste("../data/scenario_", n_scenario, "/", id_string, "_MCAR_coefs.csv", sep=''))
write.csv(final_MAR_coefs,    paste("../data/scenario_", n_scenario, "/", id_string, "_MAR_coefs.csv", sep=''))

write.csv(final_MNAR_results,   paste("../data/scenario_", n_scenario, "/", id_string, "_MNAR_results.csv", sep=''))
write.csv(final_MCAR_results,   paste("../data/scenario_", n_scenario, "/", id_string, "_MCAR_results.csv", sep=''))
write.csv(final_MAR_results,    paste("../data/scenario_", n_scenario, "/", id_string, "_MAR_results.csv", sep=''))


