#!/bin/bash
#
# ****************************************
# Outcome-MAR-Missingness Simulation Study
#
# BluePebble Launch Bash Script
# Defines and runs the R simulation procedure on the BluePebble HPC
# 
# Emma Tarmey
#
# Started:          06/10/2025
# Most Recent Edit: 21/09/2026
# ****************************************
#
#SBATCH --partition=compute
#SBATCH --job-name=scenario_1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --account=MATH033344
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aa22294@bristol.ac.uk


# Change into working directory
cd ${SLURM_SUBMIT_DIR}
cd ..
cd R

# Record info
echo ""
echo "***** START *****"
echo "***** Outcome-MAR-Missingness Simulation Study - Simulation 1 *****"
echo Start Time:        $(date)
echo Working Directory: $(pwd)
echo JOB ID:            ${SLURM_JOBID}
echo SLURM ARRAY ID:    ${SLURM_ARRAY_TASK_ID}
echo ""

# Import R
module load languages/R/4.4.1

# Execute code
Rscript modified_missingness_simulation_scenario_1.R

echo ""
echo End Time: $(date)
echo "***** END *****"
echo ""

