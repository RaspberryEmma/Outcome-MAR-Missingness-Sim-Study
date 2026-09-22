#!/bin/bash
#
# ****************************************
# Modified-Missingness Simulation Study
#
# BluePebble Automation Script
# This bash scripts automates updating to the most recent
# code version from GitHub and submitting the job to BP
#
# Emma Tarmey
#
# Started:          06/10/2025
# Most Recent Edit: 21/09/2026
# ****************************************

echo ""

# delete older version
rm -f -r Modified-Missingness-Sim-Study

# clone most recent version
git clone https://github.com/RaspberryEmma/Modified-Missingness-Sim-Study

# change wd
cd Modified-Missingness-Sim-Study
cd bash

# import python
module load languages/python/3.12.3



# submit simulation to BP HPC
for i in 1 2 3 4 5;
do
	echo   "Submitting job: launch_BP_run_"$i".sh"
	sbatch "launch_BP_scenario_"$i".sh"
done


# check jobs submitted correctly
sleep 5.0
echo ""
sacct -X
echo ""

