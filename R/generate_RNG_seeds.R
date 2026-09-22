# ****************************************
# Confounder Handling Simulation Study
#
# Random Number Generation
# Generates a CSV of random-number-generator seeds
# This ensures non-overlapping sequences between parallel scenarios
#
# Emma Tarmey
#
# Started:          06/11/2025
# Most Recent Edit: 16/12/2025
# ****************************************

# clear R memory
rm(list=ls())

# fix wd issue
# forces wd to be the location of this file
if (Sys.getenv("RSTUDIO") == "1") {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

# RNG State
set.seed(2025) # given from God
RNGkind(kind = "Mersenne-Twister")
RNGkind()
R.version

# Simulation parameters (upper limits where applicable)
n_scen        <- 6
n_rep         <- 2000
n_obs         <- 10000
n_base_var    <- (32 + 3) # number of variables with no prior cause (all confounders + unseen prior + 2 error terms)

# Required samples for each individual scenario
# Multiply by 1.16 (1 + 5/32) to account for MI
n_samples_per_scenario <- as.integer( n_rep * n_obs * n_base_var * 1.16)
print(n_samples_per_scenario)

# Required seeds
n_seeds <- n_scen
print(n_seeds)

# Check periodicity
RNG_period     <- (2^(19937) - 1)                  # Mersenne-Twister period
min_req_period <- n_samples_per_scenario * n_seeds # Minimum period value needed to
                                                   # guarantee zero overlap in simulation
message("\n\nChecking whether period of RNG algorithm is sufficient to produce non-overlapping sequences for all scenarios")
print(RNG_period)
print(min_req_period)
message(ifelse((RNG_period > min_req_period), "YES", "NO"))

# Generate all seeds
seeds_df           <- data.frame(matrix(NA, nrow = n_seeds, ncol = 2)) 
colnames(seeds_df) <- c("simulation_scenario", "seed")
random_numbers     <- runif(n_samples_per_scenario, min = 0, max = 1) # burn-in

for (i in 1:n_seeds) {
  # Extract next seed
  seed_sequence <- .Random.seed      # get current RNG state
  new_seed      <- seed_sequence[3]  # extract valid seed value
  
  message(paste("\n\nSeed", i))
  print(seed_sequence[c(1:10)])
  print(new_seed)

  seeds_df[i, "simulation_scenario"] <- i
  seeds_df[i, "seed"]                <- new_seed
  
  # generate and throw away 'n_samples_per_scenario'-many samples
  # sufficiently modifies RNG state
  set.seed(new_seed)
  random_numbers <- runif(n_samples_per_scenario, min = 0, max = 1)
}

print(seeds_df)
write.csv(seeds_df, file = "../data/precomputed_RNG_seeds.csv", row.names = FALSE)


