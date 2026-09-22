#!/bin/bash
#
# ****************************************
# Outcome-MAR-Missingness Simulation Study
#
# BluePebble Automation Script
# Local launcher for running sim on BP
# 
# Emma Tarmey
#
# Started:          06/10/2025
# Most Recent Edit: 21/09/2026
# ****************************************

# send BP script to BP server home directory
scp update_and_launch_outcome_modified.sh aa22294@bp1-login.acrc.bris.ac.uk:~


# log-in to BP
ssh -X aa22294@bp1-login.acrc.bris.ac.uk

